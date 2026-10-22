"""
Suggested next steps for `prepare status --next`.

Builds the `download_sra`, `select_datasets` and `extract_target_reads` commands that advance the
most accessions, from what the registry recorded. Recorded values are checked (the `_as_*`
coercers) and shell-quoted before they are pasted into a command.
"""

import logging
import math
import shlex
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from prepare.core.constants import DEFAULT_CONTAINMENT_THRESHOLD, DEFAULT_PARSED_CONTAINMENT_FILE, GENOME_FASTA_GLOBS
from prepare.data import registry_blocks as rb
from prepare.data.project_selection import selection_for_genome
from prepare.data.registry import (
    ProjectPaths,
    Registry,
    known_genome_ids,
    query,
    resolve_project_path,
)

logger = logging.getLogger(__name__)


def _as_float(value: Any) -> Optional[float]:
    """``value`` as a finite float, or None when it is not a number."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def as_positive_int(value: Any) -> Optional[int]:
    """``value`` as a positive int, or None when it is not one."""
    if isinstance(value, bool):
        return None
    try:
        number = int(value)
    except (TypeError, ValueError):
        return None
    return number if number > 0 else None


def as_non_negative_int(value: Any) -> Optional[int]:
    """``value`` as an int of zero or more, or None when it is not one."""
    if isinstance(value, bool):
        return None
    try:
        number = int(value)
    except (TypeError, ValueError):
        return None
    return number if number >= 0 else None


def _warn_malformed(name: str, value: Any, reason: str) -> None:
    """Warn that a recorded criteria value could not be used to build a reselect command."""
    logger.warning("Recorded %s %r %s", name, value, reason)


def _run_filter_flags(criteria: Dict[str, Any]) -> str:
    """The --max-run-size, --min-spots, --max-spots and --platform flags a selection recorded.

    The run size is a byte count above zero and the spot bounds are zero or more; a recorded value
    that fails that check is left out with a warning rather than pasted into the command. The
    platform is shell-quoted.
    """
    part = ""
    numeric = [
        ("max_run_size", "--max-run-size", as_positive_int),
        ("min_spots", "--min-spots", as_non_negative_int),
        ("max_spots", "--max-spots", as_non_negative_int),
    ]
    for key, flag, coerce in numeric:
        raw = criteria.get(key)
        if raw is None:
            continue
        number = coerce(raw)
        if number is None:
            _warn_malformed(key, raw, f"is not a valid count; the suggested command omits {flag}")
        else:
            part += f" {flag} {number}"
    platform = criteria.get("platform")
    if platform is not None:
        part += f" --platform {shlex.quote(str(platform))}"
    return part


def reselect_command(criteria: Dict[str, Any], output: str) -> str:
    """A runnable ``select_datasets`` command that redoes a selection with ``--skip-excluded``.

    Built from the criteria the original ``--no-skip-excluded`` run recorded, targeting the same
    ``--output`` so rerunning it corrects that selection's file in place. Reproduces every
    criterion ``record_selection`` stores that changes which accessions are chosen (metadata
    filter, top-N cap, run size, spot and platform filters, source table), not just the genome
    column and threshold, so the suggested command redoes the same selection rather than a
    looser one. Every value that came from the registry rather than this method's own literal
    flag text is passed through ``shlex.quote``, so a value containing a space or shell
    metacharacter (a metadata value like "New York", say) still produces a command that is safe
    to paste into a shell and run as-is. The threshold is coerced with ``float``, the top-N
    count with ``int`` (positive only) and ``require`` must be ``any`` or ``all``; a recorded
    value that fails that check (a hand-edited registry, say) leaves its flag out rather than
    being pasted into the command.
    """
    raw_threshold = criteria.get("threshold", DEFAULT_CONTAINMENT_THRESHOLD)
    threshold = _as_float(raw_threshold)
    if threshold is None and raw_threshold is not None:
        _warn_malformed("threshold", raw_threshold, "is not a number; the suggested command uses the default")
    threshold_part = f" --threshold {threshold}" if threshold is not None else ""
    genome_ids = criteria.get("genome_ids")
    if genome_ids:
        raw_require = criteria.get("require", "any")
        quoted_ids = " ".join(shlex.quote(str(g)) for g in genome_ids)
        genome_part = f"--genome-ids {quoted_ids}"
        if raw_require in ("any", "all"):
            genome_part += f" --require {raw_require}"
        elif raw_require is not None:
            _warn_malformed("require", raw_require, "is not 'any' or 'all'; the suggested command omits --require")
    else:
        column = criteria.get("column") or "max_containment"
        genome_part = f"--genome-id {shlex.quote(str(column))}"
    command = (
        f"prepare select_datasets {genome_part}{threshold_part} --skip-excluded --output {shlex.quote(str(output))}"
    )

    metadata_file = criteria.get("metadata_file")
    if metadata_file:
        command += f" --metadata-file {shlex.quote(str(metadata_file))}"
    metadata_column = criteria.get("metadata_column")
    metadata_value = criteria.get("metadata_value")
    if metadata_column and metadata_value is not None:
        command += (
            f" --metadata-column {shlex.quote(str(metadata_column))}"
            f" --metadata-value {shlex.quote(str(metadata_value))}"
        )
    raw_top_n = criteria.get("top_n")
    top_n = as_positive_int(raw_top_n)
    if top_n:
        command += f" --top-n {top_n}"
    elif raw_top_n is not None:
        _warn_malformed("top_n", raw_top_n, "is not a number; the suggested command uses the default")
    command += _run_filter_flags(criteria)
    table = criteria.get("table")
    if table and str(table) != DEFAULT_PARSED_CONTAINMENT_FILE:
        command += f" --parsed-containment {shlex.quote(str(table))}"
    return command


def download_next_steps(registry: Registry) -> List[Dict[str, Any]]:
    """Suggested commands for selected accessions that are neither excluded nor downloaded yet."""
    selected, excluded, downloaded = (set(query(registry, s)) for s in ("selected", "excluded", "downloaded"))
    wanted = selected - excluded - downloaded
    to_download = [acc for acc in registry.datasets if acc in wanted]
    if not to_download:
        return []
    # A selection recorded with --no-skip-excluded may still list an excluded
    # accession, so its output file is never suggested for direct download;
    # such accessions instead point at re-running select_datasets with
    # --skip-excluded so the excluded run is dropped before download.
    by_output: Dict[str, List[str]] = {}
    reselect_groups: Dict[str, Tuple[Dict[str, Any], List[str]]] = {}
    for acc in to_download:
        selection = rb.selection_block(registry, acc) or rb.SelectionBlock()
        criteria = selection.criteria or {}
        output = selection.output or "accessions.txt"
        if criteria.get("skip_excluded") is False:
            group = reselect_groups.setdefault(output, (criteria, []))
            group[1].append(acc)
            continue
        by_output.setdefault(output, []).append(acc)
    steps = [
        {"command": f"prepare download_sra --accessions-file {shlex.quote(str(output))}", "accessions": accs}
        for output, accs in by_output.items()
    ]
    for output, (criteria, accs) in reselect_groups.items():
        steps.append({"command": reselect_command(criteria, output), "accessions": accs})
    return steps


def _selection_table(registry: Registry) -> str:
    """The containment table a selection was recorded from, else the default name."""
    for acc in registry.datasets:
        selection = rb.selection_block(registry, acc) or rb.SelectionBlock()
        table = (selection.criteria or {}).get("table") if selection.selected else None
        if table:
            return str(table)
    return DEFAULT_PARSED_CONTAINMENT_FILE


def _genome_fasta(registry: Registry, paths: ProjectPaths, genome_id: str) -> Path:
    """The genome's FASTA: the recorded one, else a file on disk, else the conventional name."""
    recorded = (registry.genomes.get(genome_id) or {}).get("fasta")
    if recorded:
        return resolve_project_path(registry, recorded)
    for pattern in GENOME_FASTA_GLOBS:
        candidate = paths.genomes / pattern.replace("*", genome_id)
        if candidate.exists():
            return candidate
    return paths.genomes / f"{genome_id}.fna"


def _would_extract(registry: Registry, accession: str, genome_id: str, selection: set) -> bool:
    """Whether the suggested ``extract_target_reads`` command (default threshold) would take the accession.

    An accession with no screening block at all is kept, so a project built by ``status --init``
    from downloaded reads alone is not silenced. A screened one must name the genome at or above
    the default threshold, and be in the recorded selection when that selection covers the genome.
    """
    if selection and accession not in selection:
        return False
    genomes = rb.raw(registry, accession, "screening", "genomes")
    if not genomes:
        return True
    entry = genomes.get(genome_id)
    if not isinstance(entry, dict):
        return False
    containment = _as_float(entry.get("containment"))
    return containment is not None and containment >= DEFAULT_CONTAINMENT_THRESHOLD


def _extraction_next_steps(registry: Registry, paths: ProjectPaths) -> List[Dict[str, Any]]:
    steps: List[Dict[str, Any]] = []
    excluded = set(query(registry, "excluded"))
    downloaded = [acc for acc in query(registry, "downloaded") if acc not in excluded]
    table = _selection_table(registry)
    for genome_id in sorted(known_genome_ids(registry)):
        genome_fasta = _genome_fasta(registry, paths, genome_id)
        base = (
            f"prepare extract_target_reads --parsed-containment {shlex.quote(table)} "
            f"--genome-id {shlex.quote(genome_id)} --genome-fasta {shlex.quote(str(genome_fasta))}"
        )
        # Any record, including a zero-mapped one, means the sample has been tried.
        recorded = {acc for acc in registry.datasets if rb.extraction_block(registry, acc, genome_id) is not None}
        selection = selection_for_genome(registry, genome_id)
        to_extract = [
            acc for acc in downloaded if acc not in recorded and _would_extract(registry, acc, genome_id, selection)
        ]
        if to_extract:
            steps.append({"command": base, "accessions": to_extract})
        assembled = set(query(registry, "assembled", genome_id))
        to_assemble = [
            acc for acc in query(registry, "extracted", genome_id) if acc not in assembled and acc not in excluded
        ]
        if to_assemble:
            steps.append({"command": f"{base} --assemble", "accessions": to_assemble})
    return steps


def next_steps(registry: Registry, paths: ProjectPaths) -> List[Dict[str, Any]]:
    """Every suggested next command for the project: downloads first, then read extractions."""
    return download_next_steps(registry) + _extraction_next_steps(registry, paths)
