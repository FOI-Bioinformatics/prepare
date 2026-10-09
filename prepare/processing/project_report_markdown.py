"""The project report as Markdown: one heading per section, short text lines and pipe tables.

``section_blocks`` lays each section of a ``build_project_report`` dict out as a list of blocks
(a line of text, or a table with a header and rows); ``render_markdown`` writes those blocks as
ASCII Markdown, and the HTML renderer in ``prepare.visualization.project_report`` draws the same
blocks, so both formats show the same content. Table cells go through ``escape_cell``, so a ``|``
or a line break in a value (an exclusion reason, a download error message) cannot break a row.
"""

import json
import math
from typing import Any, Callable, Dict, List, NamedTuple, Sequence

from prepare.processing.project_report import SECTION_TITLES

FUNNEL_DETAIL = (
    ("excluded", "excluded"),
    ("failed", "failed"),
    ("pairs", "pairs"),
    ("bytes", "bytes"),
    ("total_bp", "total bp"),
    ("seconds", "seconds"),
    ("failed_seconds", "failed seconds"),
)
TIMING_FIELDS = ("count", "min", "p25", "median", "p75", "p90", "max")


class Block(NamedTuple):
    """One piece of a section: a line of ``text``, or a table when ``header`` is not empty."""

    text: str = ""
    header: Sequence[str] = ()
    rows: Sequence[Sequence[Any]] = ()


def format_cell(value: Any) -> str:
    """A value as table text: None empty, booleans yes/no, floats to at most 4 decimals, structures as JSON."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float) and math.isfinite(value):
        return str(int(value)) if value == int(value) else f"{value:.4f}".rstrip("0").rstrip(".")
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True)
    return str(value)


def escape_cell(value: Any) -> str:
    """One Markdown table cell: ``format_cell`` with ``|`` escaped, line breaks as ``<br>``, ASCII only."""
    text = format_cell(value).replace("\r\n", "\n").replace("\r", "\n")
    text = text.replace("|", "\\|").replace("\n", "<br>")
    return text.encode("ascii", "replace").decode("ascii")


def _text(text: str) -> Block:
    return Block(text=text)


def _table(header: Sequence[str], rows: Sequence[Sequence[Any]]) -> Block:
    return Block(header=tuple(header), rows=[list(row) for row in rows])


def _project(section: Dict[str, Any]) -> List[Block]:
    rows = [[key.replace("_", " "), value] for key, value in section.items() if value is not None]
    return [_table(["field", "value"], rows)]


def _funnel(section: Dict[str, Any]) -> List[Block]:
    rows = []
    for stage, values in section.items():
        detail = ", ".join(
            f"{label} {format_cell(values[key])}" for key, label in FUNNEL_DETAIL if values.get(key) is not None
        )
        rows.append([stage, values.get("accessions"), detail])
    return [_table(["stage", "accessions", "detail"], rows)]


def _genomes(section: Dict[str, Any]) -> List[Block]:
    if not section:
        return [_text("No genome is recorded in the registry.")]
    header = ["genome", "extracted", "assembled", "zero mapped", "median breadth", "median depth"]
    rows = [
        [genome, g["extracted"], g["assembled"], g["zero_mapped"], g["median_breadth"], g["median_depth"]]
        for genome, g in section.items()
    ]
    return [_table(header, rows)]


def _shown(shown: int, total: int) -> List[Block]:
    return [_text(f"{shown} of {total} rows shown.")] if shown < total else []


def _extractions(section: Dict[str, Any]) -> List[Block]:
    if not section["rows_total"]:
        return [_text("No extraction is recorded.")]
    columns = section["columns"]
    blocks = [_table(columns, [[row.get(column) for column in columns] for row in section["rows"]])]
    if section.get("note"):
        return blocks + [_text(section["note"])]
    return blocks + _shown(len(section["rows"]), section["rows_total"])


def _downloads(section: Dict[str, Any]) -> List[Block]:
    if not section["states"]:
        return [_text("No download is recorded.")]
    blocks = [
        _text("Download states:"),
        _table(["state", "datasets"], list(section["states"].items())),
        _text("Completeness verdicts of the downloaded datasets:"),
        _table(["verdict", "datasets"], list(section["verdicts"].items())),
    ]
    for name in ("truncated", "unverified"):
        listed = section[name]
        if listed:
            total = section[f"{name}_total"]
            more = f" (first {len(listed)} of {total})" if len(listed) < total else ""
            blocks.append(_text(f"{name.capitalize()}{more}: {', '.join(listed)}"))
    return blocks


def _failures(section: Dict[str, Any]) -> List[Block]:
    if not section["rows_total"]:
        return [_text("No download is recorded as failed.")]
    reasons = ", ".join(f"{reason} {count}" for reason, count in section["by_reason"].items())
    columns = section["columns"]
    blocks = [
        _text(f"Failed downloads by reason: {reasons}. attempts_total counts the attempts over all runs."),
        _table(columns, [[row.get(column) for column in columns] for row in section["rows"]]),
    ]
    return blocks + _shown(len(section["rows"]), section["rows_total"])


def _timing(section: Dict[str, Any]) -> List[Block]:
    rows = []
    for kind, spread in section["distribution"].items():
        total = section["summary"].get(f"{kind}_seconds_total")
        rows.append([kind, *(spread[field] for field in TIMING_FIELDS), total])
    return [
        _text("Seconds per step, from the times recorded in the registry:"),
        _table(["step", *TIMING_FIELDS, "total"], rows),
    ]


def _environment(section: Dict[str, Any]) -> List[Block]:
    if not section.get("included"):
        return [_text("Environment checks not included (--no-environment).")]
    counts = ", ".join(f"{status} {count}" for status, count in section["counts"].items())
    blocks = [_text(f"Overall status: {section['status']} ({counts}); network checks were not run.")]
    if not section["problems"]:
        return blocks + [_text("Every check passed.")]
    rows = [[p["name"], p["status"], p["detail"]] for p in section["problems"]]
    return blocks + [_table(["check", "status", "detail"], rows)]


def _outputs(section: Dict[str, Any]) -> List[Block]:
    blocks: List[Block] = []
    if section["exports"]:
        rows = [[name, e.get("date"), e.get("output"), e.get("summary")] for name, e in section["exports"].items()]
        blocks.append(_table(["export", "date", "output", "summary"], rows))
    if section["analyses"]:
        rows = [
            [name, a["datasets"], a["latest"], a["outputs_total"], a["outputs"][0] if a["outputs"] else None]
            for name, a in section["analyses"].items()
        ]
        blocks.append(_table(["analysis", "datasets", "latest", "outputs", "first output"], rows))
    return blocks or [_text("No export or analysis output is recorded.")]


def _runs(section: Dict[str, Any]) -> List[Block]:
    if not section["runs"]:
        return [_text("The run log holds no runs.")]
    rows = [
        [r["run_id"], r["command"], r["started"], r["seconds"], r["exit_code"], r["summary"]] for r in section["runs"]
    ]
    return [
        _table(["run", "command", "started (UTC)", "seconds", "exit", "summary"], rows),
        _text(f"{len(section['runs'])} of {section['total']} recorded runs, newest first."),
    ]


_LAYOUT: Dict[str, Callable[[Dict[str, Any]], List[Block]]] = {
    "project": _project,
    "funnel": _funnel,
    "genomes": _genomes,
    "extractions": _extractions,
    "downloads": _downloads,
    "failures": _failures,
    "timing": _timing,
    "environment": _environment,
    "outputs": _outputs,
    "runs": _runs,
}


def section_blocks(report: Dict[str, Any]) -> Dict[str, List[Block]]:
    """The blocks of every section of ``report``, keyed and ordered as ``SECTION_TITLES``."""
    return {name: _LAYOUT[name](report[name]) for name in SECTION_TITLES}


def subtitle(report: Dict[str, Any]) -> str:
    """The line under the title: when, by which PREPARE version and from which registry."""
    project = report.get("project", {})
    return (
        f"Generated {report.get('generated')} by PREPARE {project.get('prepare_version')} "
        f"from {project.get('registry')}."
    )


def _markdown(block: Block) -> List[str]:
    if not block.header:
        return [block.text.encode("ascii", "replace").decode("ascii")]
    lines = ["| " + " | ".join(escape_cell(cell) for cell in block.header) + " |"]
    lines.append("|" + "|".join("---" for _ in block.header) + "|")
    lines.extend("| " + " | ".join(escape_cell(cell) for cell in row) + " |" for row in block.rows)
    return lines


def render_markdown(report: Dict[str, Any]) -> str:
    """The report as one Markdown document: a title, then one ``##`` section per report section."""
    lines = ["# PREPARE project report", "", subtitle(report).encode("ascii", "replace").decode("ascii")]
    for name, blocks in section_blocks(report).items():
        lines += ["", f"## {SECTION_TITLES[name]}"]
        for block in blocks:
            lines += [""] + _markdown(block)
    return "\n".join(lines) + "\n"
