"""The project report: one summary of a project, built from its registry and run log.

``build_project_report`` reads the registry (and, unless told not to, runs the environment checks
of ``doctor`` without network access) and returns a JSON-ready dict with one key per section in
``SECTIONS``: the project, the cross-stage funnel of ``status``, per-genome counts, the extraction
rows of ``results_table``, download verdicts, failed downloads, timing, the environment, recorded
outputs and the most recent runs of the run log. Long tables are cut to ``max_rows`` rows and say
how many there were. Nothing here opens a FASTQ file or writes anything; the Markdown and HTML
renderers are ``prepare.processing.project_report_markdown`` and
``prepare.visualization.project_report``, and the ``project_report`` command writes the files.
"""

import math
import statistics
from collections import Counter
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Sequence

from prepare import __version__
from prepare.data import registry_blocks as rb
from prepare.data import run_log
from prepare.data.registry import Registry, project_root, stage_members
from prepare.data.registry_timing import TIMED_DOWNLOAD_STATES, timing_summary
from prepare.data.sra.run_report import failure_reason
from prepare.processing.doctor_report import OK, overall_status, run_checks, summary_counts
from prepare.processing.project_funnel import funnel
from prepare.processing.results import results_rows
from prepare.processing.status_report import genome_counts

REPORT_VERSION = 1
DEFAULT_MAX_ROWS = 200
DEFAULT_RUNS_LIMIT = 10

SECTION_TITLES: Dict[str, str] = {
    "project": "Project",
    "funnel": "Funnel",
    "genomes": "Genomes",
    "extractions": "Extractions",
    "downloads": "Downloads",
    "failures": "Failed downloads",
    "timing": "Timing",
    "environment": "Environment",
    "outputs": "Outputs",
    "runs": "Recent runs",
}
SECTIONS = tuple(SECTION_TITLES)

# The results_table columns shown per extraction; results_table itself has every column.
EXTRACTION_COLUMNS = (
    "accession",
    "genome_id",
    "containment",
    "mapped_reads",
    "mapping_rate_to_reference",
    "breadth",
    "mean_depth",
    "contigs",
    "total_bp",
    "n50",
    "quality_grade",
    "quality_source",
    "extraction_seconds",
    "assembly_seconds",
)
FAILURE_COLUMNS = ("accession", "reason", "attempts_total", "date", "message")
RUN_FIELDS = ("run_id", "command", "started", "seconds", "exit_code", "summary")
DOWNLOAD_VERDICTS = ("complete", "truncated", "unverified")
TIMING_KINDS = ("download", "extraction", "assembly")
PERCENTILES = (("p25", 0.25), ("median", 0.5), ("p75", 0.75), ("p90", 0.9))


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _cut(items: Sequence[Any], max_rows: int) -> List[Any]:
    """The first ``max_rows`` of ``items``, or all of them when ``max_rows`` is 0."""
    return list(items[:max_rows]) if max_rows else list(items)


def _number(value: Any) -> Optional[float]:
    """``value`` as a finite float, or None for a missing, boolean or non-numeric value."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def percentile(ordered: Sequence[float], fraction: float) -> float:
    """The ``fraction`` quantile of the sorted, non-empty ``ordered``, by linear interpolation."""
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = min(lower + 1, len(ordered) - 1)
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def distribution(values: Sequence[float]) -> Dict[str, Any]:
    """Count, minimum, quartiles, 90th percentile and maximum of ``values`` (None when empty)."""
    ordered = sorted(values)
    result: Dict[str, Any] = {"count": len(ordered), "min": round(ordered[0], 3) if ordered else None}
    for name, fraction in PERCENTILES:
        result[name] = round(percentile(ordered, fraction), 3) if ordered else None
    result["max"] = round(ordered[-1], 3) if ordered else None
    return result


# --- sections -----------------------------------------------------------------------------------


def _project_section(registry: Registry) -> Dict[str, Any]:
    project = rb.project_block(registry)
    return {
        "registry": str(registry.path) if registry.path is not None else None,
        "root": str(project_root(registry)),
        "name": project.name,
        "id": project.id,
        "created": registry.created,
        "updated": registry.updated,
        "registry_version": registry.version,
        "prepare_version": __version__,
        "datasets": len(registry.datasets),
        "genomes": len(registry.genomes),
    }


def _median(values: List[float]) -> Optional[float]:
    return round(statistics.median(values), 4) if values else None


def _genomes_section(registry: Registry) -> Dict[str, Any]:
    """Per genome: extracted and assembled counts, zero-mapped count, median breadth and depth."""
    breadths: Dict[str, List[float]] = {}
    depths: Dict[str, List[float]] = {}
    for record in registry.datasets.values():
        for genome_id, extraction in (record.get("extractions") or {}).items():
            if not isinstance(extraction, dict) or not (extraction.get("mapped_reads") or 0):
                continue
            breadth, depth = _number(extraction.get("breadth")), _number(extraction.get("mean_depth"))
            if breadth is not None:
                breadths.setdefault(genome_id, []).append(breadth)
            if depth is not None:
                depths.setdefault(genome_id, []).append(depth)
    return {
        genome_id: {
            "extracted": counts["extracted"],
            "assembled": counts["assembled"],
            "zero_mapped": len(counts["zero_mapped"]),
            "median_breadth": _median(breadths.get(genome_id, [])),
            "median_depth": _median(depths.get(genome_id, [])),
        }
        for genome_id, counts in genome_counts(registry).items()
    }


def _extractions_section(registry: Registry, max_rows: int) -> Dict[str, Any]:
    """The results_table rows that have an extraction, by genome and then by decreasing mapped reads."""
    rows = [row for row in results_rows(registry) if row["mapped_reads"] is not None]
    rows.sort(key=lambda row: (row["genome_id"], -row["mapped_reads"], row["accession"]))
    shown = [{column: row[column] for column in EXTRACTION_COLUMNS} for row in _cut(rows, max_rows)]
    truncated = len(shown) < len(rows)
    note = None
    if truncated:
        note = (
            f"{len(shown)} of {len(rows)} extraction rows shown; "
            "`prepare results_table` writes every row with all columns"
        )
    return {
        "columns": list(EXTRACTION_COLUMNS),
        "rows": shown,
        "rows_total": len(rows),
        "truncated": truncated,
        "note": note,
    }


def _downloads_section(registry: Registry, max_rows: int) -> Dict[str, Any]:
    """Download states, completeness verdicts of the downloaded datasets, truncated and unverified lists.

    The verdict counts and the two lists cover the same datasets, those whose download state is
    ``downloaded`` (a verdict kept from an earlier download of a dataset now missing or failed is
    not counted or listed). Each list is cut to ``max_rows``; ``<name>_total`` is its full length
    and ``<name>_note`` says how many are shown when the cut removed any (None otherwise).
    """
    states: Counter = Counter()
    verdicts: Dict[str, int] = {**{name: 0 for name in DOWNLOAD_VERDICTS}, "none": 0}
    listed: Dict[str, List[str]] = {"truncated": [], "unverified": []}
    for accession in sorted(registry.datasets):
        state = rb.raw(registry, accession, "download", "state")
        if not state:
            continue
        states[str(state)] += 1
        if state != "downloaded":
            continue
        complete = rb.raw(registry, accession, "download", "complete")
        verdict = complete.get("verdict") if isinstance(complete, dict) else None
        verdicts[verdict if verdict in DOWNLOAD_VERDICTS else "none"] += 1
        if verdict in listed:
            listed[verdict].append(accession)
    section: Dict[str, Any] = {"states": dict(sorted(states.items())), "verdicts": verdicts}
    for name, accessions in listed.items():
        shown = _cut(accessions, max_rows)
        section[name] = shown
        section[f"{name}_total"] = len(accessions)
        section[f"{name}_note"] = (
            f"{len(shown)} of {len(accessions)} {name} datasets shown" if len(shown) < len(accessions) else None
        )
    return section


def _failures_section(registry: Registry, max_rows: int) -> Dict[str, Any]:
    """Every dataset whose download state is ``failed``, with the reason read from its message.

    ``attempts_total`` is the registry's attempt count, summed over all runs (``download_run.json``
    and the report CSV of ``download_sra`` count the attempts of one run).
    """
    rows = []
    for accession in sorted(registry.datasets):
        download = registry.datasets[accession].get("download")
        if not isinstance(download, dict) or download.get("state") != "failed":
            continue
        message = str(download.get("message") or "")
        rows.append(
            {
                "accession": accession,
                "reason": failure_reason(message),
                "attempts_total": download.get("attempts"),
                "date": download.get("date") or None,
                "message": message,
            }
        )
    return {
        "columns": list(FAILURE_COLUMNS),
        "rows": _cut(rows, max_rows),
        "rows_total": len(rows),
        "by_reason": dict(sorted(Counter(row["reason"] for row in rows).items())),
    }


def _timing_section(registry: Registry) -> Dict[str, Any]:
    """``timing_summary`` and, per kind, the spread of the recorded times (same selection rules)."""
    values: Dict[str, List[float]] = {kind: [] for kind in TIMING_KINDS}
    for record in registry.datasets.values():
        download = record.get("download")
        if isinstance(download, dict) and download.get("state") in TIMED_DOWNLOAD_STATES:
            _append(values["download"], download.get("seconds"))
        for extraction in (record.get("extractions") or {}).values():
            if not isinstance(extraction, dict):
                continue
            _append(values["extraction"], extraction.get("seconds"))
            assembly = extraction.get("assembly")
            if isinstance(assembly, dict):
                _append(values["assembly"], assembly.get("seconds"))
    return {
        "summary": timing_summary(registry),
        "distribution": {kind: distribution(values[kind]) for kind in TIMING_KINDS},
    }


def _append(values: List[float], raw: Any) -> None:
    number = _number(raw)
    if number is not None:
        values.append(number)


def _environment_section(registry: Registry, include: bool) -> Dict[str, Any]:
    """The ``doctor`` checks without network access: overall status, counts and the checks not ok."""
    if not include:
        return {"included": False}
    checks = run_checks(project=project_root(registry), network=False)
    return {
        "included": True,
        "status": overall_status(checks),
        "counts": summary_counts(checks),
        "problems": [
            {"name": check.name, "status": check.status, "detail": check.detail}
            for check in checks
            if check.status != OK
        ],
    }


def _outputs_section(registry: Registry, max_rows: int) -> Dict[str, Any]:
    """The recorded project exports, and per analysis name the datasets and distinct output paths."""
    exports = {name: entry.to_dict() for name, entry in sorted(rb.project_block(registry).exports.items())}
    analyses: Dict[str, Dict[str, Any]] = {}
    for record in registry.datasets.values():
        for name, entry in (record.get("analyses") or {}).items():
            if not isinstance(entry, dict):
                continue
            info = analyses.setdefault(name, {"datasets": 0, "latest": None, "outputs": set()})
            info["datasets"] += 1
            if entry.get("output"):
                info["outputs"].add(str(entry["output"]))
            date = str(entry.get("date") or "")
            if date and (info["latest"] is None or date > info["latest"]):
                info["latest"] = date
    return {
        "exports": exports,
        "analyses": {
            name: {
                "datasets": info["datasets"],
                "latest": info["latest"],
                "outputs": _cut(sorted(info["outputs"]), max_rows),
                "outputs_total": len(info["outputs"]),
            }
            for name, info in sorted(analyses.items())
        },
    }


def _runs_section(registry: Registry, runs_limit: int) -> Dict[str, Any]:
    """The last ``runs_limit`` lines of the project's run log, newest first."""
    records = run_log.read_runs(project_root(registry))
    newest = list(reversed(records))
    return {
        "total": len(records),
        "runs": [{name: record.to_dict()[name] for name in RUN_FIELDS} for record in _cut(newest, runs_limit)],
    }


def build_project_report(
    registry: Registry,
    *,
    max_rows: int = DEFAULT_MAX_ROWS,
    include_environment: bool = True,
    runs_limit: int = DEFAULT_RUNS_LIMIT,
) -> Dict[str, Any]:
    """The project report of ``registry`` as a JSON-ready dict, one key per name in ``SECTIONS``.

    ``max_rows`` cuts the extraction, failure and accession lists (0 keeps every row);
    ``include_environment`` runs the ``doctor`` checks without network access; ``runs_limit`` is
    the number of run-log lines kept, newest first (0 keeps all). Reads the registry, the run log
    and, for the environment, what ``doctor`` reads; never a FASTQ file.
    """
    members = stage_members(registry)
    return {
        "report_version": REPORT_VERSION,
        "generated": _now(),
        "max_rows": max_rows,
        "project": _project_section(registry),
        "funnel": funnel(registry, members),
        "genomes": _genomes_section(registry),
        "extractions": _extractions_section(registry, max_rows),
        "downloads": _downloads_section(registry, max_rows),
        "failures": _failures_section(registry, max_rows),
        "timing": _timing_section(registry),
        "environment": _environment_section(registry, include_environment),
        "outputs": _outputs_section(registry, max_rows),
        "runs": _runs_section(registry, runs_limit),
    }
