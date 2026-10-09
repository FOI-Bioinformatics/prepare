"""
Text rendering of the `status` report.

Each `_print_*` function writes one section of the report through the `emit` callable it is
given (the command's `BaseCommand.emit`), so stdout stays the command's single channel.
"""

import argparse
from typing import Any, Callable, Dict, List, Optional

from prepare.cli.commands.status.suggest import as_non_negative_int, as_positive_int
from prepare.data import registry_blocks as rb
from prepare.data.registry import Registry, STAGES, known_genome_ids, query
from prepare.data.registry_reconcile import NAMED_IN_WARNING
from prepare.processing.status_report import download_verdicts, stage_filter_accessions


def _format_bytes(count: int) -> str:
    """A byte count in decimal units, as the run size flags read them (e.g. 600 MB, 1.5 GB)."""
    for factor, unit in ((10**12, "TB"), (10**9, "GB"), (10**6, "MB"), (10**3, "KB")):
        if count >= factor:
            return f"{count / factor:g} {unit}"
    return f"{count} bytes"


def _format_hours(seconds: float) -> str:
    """``seconds`` as hours: one decimal place below 10 h, none at or above it."""
    hours = seconds / 3600
    return f"{hours:.0f} h" if hours >= 10 else f"{hours:.1f} h"


def _format_size(count: int) -> str:
    """A byte count with one decimal in the largest decimal unit it reaches (e.g. 4.5 GB, 1.2 TB)."""
    for factor, unit in ((10**12, "TB"), (10**9, "GB"), (10**6, "MB"), (10**3, "KB")):
        if count >= factor:
            return f"{count / factor:.1f} {unit}"
    return f"{count} bytes"


def funnel_line(funnel: Dict[str, Any]) -> str:
    """The one-line text form of a `status` report's `funnel` block.

    Always one line, in stage order, with the download stage's total size and, when any download
    was timed, its total time in parentheses (downloads recorded before 0.7.0 carry no time, so
    the time is left out rather than shown as 0). ``selected``'s excluded count,
    ``downloaded``'s failed count and time, and the ``extracted``/``assembled`` pair counts,
    bytes and seconds are detail this line leaves to the JSON report.
    """
    downloaded = funnel["downloaded"]
    parts = [_format_size(downloaded["bytes"])]
    if downloaded["seconds"] is not None:
        parts.append(_format_hours(downloaded["seconds"]))
    detail = f" ({', '.join(parts)})"
    return (
        f"funnel: {funnel['screened']['accessions']} screened, "
        f"{funnel['selected']['accessions']} selected, "
        f"{downloaded['accessions']} downloaded{detail}, "
        f"{funnel['extracted']['accessions']} extracted, "
        f"{funnel['assembled']['accessions']} assembled"
    )


def _run_filter_detail(criteria: Dict[str, Any]) -> List[str]:
    """The run filters a selection recorded, for the selected stage row.

    A malformed value is left out; ``status --next`` warns about it when it builds the
    reselect command.
    """
    parts = []
    size = as_positive_int(criteria.get("max_run_size"))
    if size is not None:
        parts.append(f"max size {_format_bytes(size)}")
    for key, symbol in (("min_spots", ">="), ("max_spots", "<=")):
        spots = as_non_negative_int(criteria.get(key))
        if spots is not None:
            parts.append(f"spots {symbol} {spots}")
    if criteria.get("platform") is not None:
        parts.append(f"platform {criteria['platform']}")
    return parts


def _print_store(store: Dict[str, Any], emit: Callable[[str], None]) -> None:
    emit("\nStore")
    emit("=====")
    emit(f"  Root : {store['root']}")
    if not store.get("available", True):
        emit("  Unavailable: the store could not be read from here")
        return
    for state, count in sorted(store["datasets"].items()):
        emit(f"  {state:<10s}: {count}")


def _print_inventory(report: Dict[str, Any], list_missing: bool, emit: Callable[[str], None]) -> None:
    od = report["on_disk"]
    emit("Local inventory")
    emit("===============")
    emit(f"  FASTQ accessions on disk : {od['fastq_accessions']}")
    emit(f"  Metadata XML on disk     : {od['metadata_xml']}")
    emit(f"  Genome FASTA on disk     : {od['genome_fasta']}")

    w = report.get("wanted")
    if w:
        incomplete_links = w.get("fastq_incomplete_store_links") or []
        emit(f"\nReconciled against {w['total']} wanted accession(s)")
        fastq_line = f"  FASTQ    : {w['fastq_present']} present, {len(w['fastq_missing'])} missing"
        if incomplete_links:
            fastq_line += f", {len(incomplete_links)} linked to a store dataset that is not complete"
        emit(fastq_line)
        emit(f"  Metadata : {w['metadata_present']} present, {len(w['metadata_missing'])} missing")
        if list_missing:
            if w["fastq_missing"]:
                emit("  Missing FASTQ    : " + ", ".join(w["fastq_missing"]))
            if incomplete_links:
                emit("  Incomplete store links : " + ", ".join(incomplete_links))
            if w["metadata_missing"]:
                emit("  Missing metadata : " + ", ".join(w["metadata_missing"]))


def selection_detail(registry: Registry) -> str:
    """The criteria and date of the most recent selection, for the selected stage row."""
    latest: Optional[rb.SelectionBlock] = None
    for acc in registry.datasets:
        selection = rb.selection_block(registry, acc)
        if selection is not None and selection.selected and str(selection.date) >= str(latest.date if latest else ""):
            latest = selection
    if latest is None:
        return ""
    criteria = latest.criteria or {}
    parts = []
    if criteria.get("column"):
        parts.append(f"column {criteria['column']}")
    if criteria.get("threshold") is not None:
        parts.append(f"threshold {criteria['threshold']}")
    if criteria.get("metadata_column"):
        parts.append(f"{criteria['metadata_column']} = {criteria.get('metadata_value')}")
    parts.extend(_run_filter_detail(criteria))
    parts.append(str(latest.date))
    return ", ".join(p for p in parts if p)


def _exclusion_detail(registry: Registry) -> str:
    """How many accessions carry each exclusion reason."""
    reasons: Dict[str, int] = {}
    for acc in registry.datasets:
        exclusion = rb.exclusion_block(registry, acc)
        if exclusion is not None and exclusion.excluded:
            reason = str(exclusion.reason or "no reason given")
            reasons[reason] = reasons.get(reason, 0) + 1
    return ", ".join(f"{reason}: {count}" for reason, count in sorted(reasons.items()))


def _print_stages(stages: Dict[str, Any], registry: Registry, emit: Callable[[str], None]) -> None:
    emit("\nStages")
    emit("======")
    details = {
        "selected": selection_detail(registry),
        "excluded": _exclusion_detail(registry),
    }
    for stage in STAGES:
        info = stages[stage]
        detail = details.get(stage)
        emit(f"  {stage:<10s} : {info['count']}" + (f"   {detail}" if detail else ""))
    truncated = download_verdicts(registry)["truncated"]
    if truncated:
        emit(f"  truncated downloads: {len(truncated)} (" + ", ".join(truncated) + ")")


def _print_funnel(report: Dict[str, Any], emit: Callable[[str], None]) -> None:
    """Write the report's one-line cross-stage funnel (``build_report`` always sets the key)."""
    f = report.get("funnel")
    if f:
        emit(funnel_line(f))


def _print_timing(timing: Optional[Dict[str, Any]], emit: Callable[[str], None]) -> None:
    """One line with how many downloads, extractions and assemblies were timed; none when nothing was."""
    kinds = (("downloads", "download"), ("extractions", "extraction"), ("assemblies", "assembly"))
    if not timing or not any(timing.get(f"{plural}_timed") for plural, _ in kinds):
        return
    parts = []
    for plural, prefix in kinds:
        count = timing.get(f"{plural}_timed") or 0
        part = f"{plural} {count}"
        if count:
            total = timing[f"{prefix}_seconds_total"]
            median = timing[f"{prefix}_seconds_median"]
            part += f" ({total:.1f} s in total, median {median:.1f} s)"
        parts.append(part)
    emit("  timing     : " + ", ".join(parts))


def _print_genomes(genomes: Dict[str, Any], emit: Callable[[str], None]) -> None:
    if not genomes:
        return
    emit("\nGenomes")
    emit("=======")
    for genome_id, info in genomes.items():
        zero_mapped = len(info["zero_mapped"])
        empty_dirs = len(info["empty_assembly_dirs"])
        emit(
            f"  {genome_id}: extracted {info['extracted']} ({zero_mapped} with 0 mapped reads), "
            f"assembled {info['assembled']}, {empty_dirs} empty assembly dir(s)"
        )


def _print_stage_filter(
    registry: Registry, stage: Optional[str], genomes: Optional[List[str]], emit: Callable[[str], None]
) -> None:
    if not stage:
        return
    accs = stage_filter_accessions(registry, stage, genomes)
    detail = f" (genome {', '.join(genomes)})" if genomes else ""
    emit(f"\nStage '{stage}'{detail}: " + (", ".join(accs) if accs else "(none)"))


def _print_gaps(registry: Registry, emit: Callable[[str], None]) -> None:
    selected = set(query(registry, "selected"))
    excluded = set(query(registry, "excluded"))
    downloaded_set = set(query(registry, "downloaded"))
    not_downloaded = sorted(selected - excluded - downloaded_set)
    emit("\nGaps")
    emit("====")
    emit("  Selected but not downloaded : " + (", ".join(not_downloaded) if not_downloaded else "(none)"))
    for genome_id in sorted(known_genome_ids(registry)):
        gap = sorted(downloaded_set - set(query(registry, "extracted", genome_id)))
        if gap:
            emit(f"  Downloaded but not extracted for {genome_id} : " + ", ".join(gap))


def _named(items: List[str]) -> str:
    """The first ``NAMED_IN_WARNING`` of ``items``, then how many more there are."""
    more = len(items) - NAMED_IN_WARNING
    return ", ".join(items[:NAMED_IN_WARNING]) + (f" and {more} more" if more > 0 else "")


def _print_reconcile_notes(drift: Dict[str, Any], emit: Callable[[str], None]) -> None:
    """Warn about datasets an unmounted store left alone, and report re-checks, fills and drops.

    Each line appears only when its count is non-zero, so a reconcile that finds none of these
    leaves the drift section's existing text exactly as before this was added.
    """
    unavailable = drift.get("store_unavailable") or []
    if unavailable:
        emit(
            f"WARNING: {len(unavailable)} dataset(s) link into a data store that is not mounted "
            "and were not marked missing: " + _named(unavailable)
        )
    rechecked = drift.get("verdicts_rechecked") or []
    if rechecked:
        emit(f"  Verdicts re-checked                     : {len(rechecked)}")
    filled = drift.get("metadata_filled") or []
    if filled:
        emit(f"  Metadata filled from XML                : {len(filled)}")
    dropped = drift.get("assemblies_dropped") or []
    if dropped:
        emit(f"  Assemblies older than their extraction  : {len(dropped)} record(s) dropped: " + _named(dropped))


def _print_drift(drift: Dict[str, Any], emit: Callable[[str], None]) -> None:
    if not drift:
        return
    emit("\nDrift against disk")
    emit("===================")
    emit(
        "  Recorded downloaded but missing on disk : "
        + (", ".join(drift["recorded_missing"]) if drift["recorded_missing"] else "(none)")
    )
    emit(
        "  On disk but not tracked as downloaded   : "
        + (", ".join(drift["untracked_fastq"]) if drift["untracked_fastq"] else "(none)")
    )
    if drift["untracked_extractions"]:
        pairs = ", ".join(f"{acc}/{genome_id}" for acc, genome_id in drift["untracked_extractions"])
        emit(f"  Untracked extractions                   : {pairs}")
    if drift["empty_assembly_dirs"]:
        pairs = ", ".join(f"{acc}/{genome_id}" for acc, genome_id in drift["empty_assembly_dirs"])
        emit(f"  Empty assembly directories               : {pairs}")
    if drift.get("dangling_links"):
        emit("  Store links with a missing target       : " + ", ".join(drift["dangling_links"]))
    _print_reconcile_notes(drift, emit)


def _print_next(steps: List[Dict[str, Any]], emit: Callable[[str], None]) -> None:
    if not steps:
        return
    emit("\nSuggested next steps")
    emit("=====================")
    for step in steps:
        emit(f"  {step['command']}")
        emit("    accessions: " + ", ".join(step["accessions"]))


def print_report(
    args: argparse.Namespace, report: Dict[str, Any], registry: Registry, emit: Callable[[str], None]
) -> None:
    """Write the text form of a ``status`` report through ``emit``, one section after another."""
    _print_inventory(report, args.list_missing, emit)
    if report.get("store"):
        _print_store(report["store"], emit)
    _print_stages(report["stages"], registry, emit)
    _print_funnel(report, emit)
    _print_timing(report.get("timing"), emit)
    _print_genomes(report["genomes"], emit)
    _print_stage_filter(registry, args.stage, args.genome, emit)
    if args.list_missing:
        _print_gaps(registry, emit)
    _print_drift(report["drift"], emit)
    if args.next:
        _print_next(report.get("next", []), emit)
