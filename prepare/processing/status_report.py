"""
Building the `status` report.

Pure report builders for `prepare status`: the local inventory reconciled against the wanted
accessions, the shared data store block, per-genome extraction and assembly counts, drift against
the disk and download verdicts. `build_report` assembles them into the dict that `status --json`
prints; the text rendering and the suggested next steps live in `prepare.cli.commands.status`.
"""

import argparse
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Tuple

from prepare.core.constants import GENOME_FASTA_GLOBS
from prepare.core.exceptions import DataAccessError, PrepareError
from prepare.data import registry_blocks as rb
from prepare.data.file_io import is_hidden_name, visible_files
from prepare.data.registry import (
    ProjectPaths,
    Registry,
    STAGES,
    empty_assembly_dirs,
    known_genome_ids,
    query,
    stage_members,
)
from prepare.data.registry_reconcile import StoreReconcileReport
from prepare.data.registry_timing import timing_summary
from prepare.data.sra import STORE_READY_STATES, accession_has_fastq, is_transient_folder
from prepare.processing.project_funnel import funnel
from prepare.store.catalog import Catalog
from prepare.store.layout import StorePaths, sidecar_path, store_paths
from prepare.store.link import is_store_link
from prepare.store.resolve import resolve_store_root
from prepare.store.sidecar import read_sidecar

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)


def _wanted_accessions(args: argparse.Namespace) -> List[str]:
    """Load the wanted accession list from a file and/or a parsed containment table."""
    wanted: List[str] = []
    if args.accessions_file:
        path = Path(args.accessions_file)
        if not path.exists():
            raise PrepareError(f"Accessions file not found: {path}")
        wanted += [ln.strip() for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    if args.parsed_containment:
        import pandas as pd

        path = Path(args.parsed_containment)
        if not path.exists():
            raise PrepareError(f"Parsed containment file not found: {path}")
        # Only the index (accession) column is needed, so the genome columns are never parsed.
        wanted += [str(i) for i in pd.read_csv(path, sep="\t", index_col=0, usecols=[0]).index]
    # de-duplicate, preserve order
    return list(dict.fromkeys(wanted))


def _resolve_wanted(args: argparse.Namespace, registry: Registry) -> List[str]:
    """Explicit accessions/containment file win; otherwise fall back to the registry's selection."""
    if args.accessions_file or args.parsed_containment:
        return _wanted_accessions(args)
    return query(registry, "selected")


def _reconcile_present_missing(wanted: List[str], present_fn) -> Tuple[List[str], List[str]]:
    """Split a wanted list into (present, missing) using a predicate."""
    present = [a for a in wanted if present_fn(a)]
    present_set = set(present)
    missing = [a for a in wanted if a not in present_set]
    return present, missing


def _listed_or_probed(
    listed: List[str], unlisted: Callable[[str], bool], probe: Callable[[str], bool]
) -> Callable[[str], bool]:
    """A presence test that looks ``name`` up in a folder listing already made.

    ``unlisted(name)`` is True for a name the listing leaves out by design (a hidden or a
    transient folder name); only such a name is checked with ``probe(name)`` on disk, so a
    wanted list costs no filesystem call per accession. Names are compared exactly, so a name
    that differs from a listed one only in case reads as missing, on every filesystem.
    """
    names = set(listed)
    return lambda name: name in names or (unlisted(name) and probe(name))


def inventory_report(
    args: argparse.Namespace, registry: Registry, store: Optional[StorePaths] = None
) -> Dict[str, Any]:
    """The on-disk inventory part of the report: FASTQ, metadata and genome files against what is wanted."""
    fastq_dir = Path(args.fastq_folder)
    meta_dir = Path(args.metadata_folder)
    genomes_dir = Path(args.genomes_folder)

    on_disk_fastq = sorted(
        d.name
        for d in visible_files(fastq_dir, dirs=True)
        if not is_transient_folder(d.name) and accession_has_fastq(d)
    )
    on_disk_meta = sorted(p.name[: -len("_metadata.xml")] for p in visible_files(meta_dir, "*_metadata.xml"))
    on_disk_genomes = sorted(p.name for p in visible_files(genomes_dir, *GENOME_FASTA_GLOBS))

    report: Dict[str, Any] = {
        "on_disk": {
            "fastq_accessions": len(on_disk_fastq),
            "metadata_xml": len(on_disk_meta),
            "genome_fasta": len(on_disk_genomes),
        }
    }

    wanted = _resolve_wanted(args, registry)
    if wanted:
        # Presence is read from the listings above: on_disk_fastq holds exactly the visible,
        # non-transient folders for which accession_has_fastq is True.
        fastq_present, fastq_missing = _reconcile_present_missing(
            wanted,
            _listed_or_probed(
                on_disk_fastq,
                lambda a: is_hidden_name(a) or is_transient_folder(a),
                lambda a: accession_has_fastq(fastq_dir / a),
            ),
        )
        meta_present, meta_missing = _reconcile_present_missing(
            wanted,
            _listed_or_probed(on_disk_meta, is_hidden_name, lambda a: (meta_dir / f"{a}_metadata.xml").exists()),
        )
        report["wanted"] = {
            "total": len(wanted),
            "fastq_present": len(fastq_present),
            "fastq_missing": fastq_missing,
            "fastq_incomplete_store_links": _incomplete_store_links(fastq_dir, fastq_missing, store),
            "metadata_present": len(meta_present),
            "metadata_missing": meta_missing,
        }
    return report


def _incomplete_store_links(fastq_dir: Path, missing: List[str], store: Optional[StorePaths]) -> List[str]:
    """Missing accessions whose ``fastq/<ACC>`` is a symlink into the store, but the store's
    recorded sidecar state for that dataset falls outside ``STORE_READY_STATES``.

    Such a link is not simply absent: a download into the store was attempted (and left a
    ``failed`` or ``partial`` sidecar, or one is still ``downloading``), so the report says
    why the accession reads as missing rather than leaving that to be rediscovered by hand.
    """
    if store is None:
        return []
    incomplete = []
    for acc in missing:
        if not is_store_link(fastq_dir / acc, store):
            continue
        sidecar = read_sidecar(sidecar_path(store, acc))
        if sidecar is not None and sidecar.state not in STORE_READY_STATES:
            incomplete.append(acc)
    return sorted(incomplete)


def stage_filter_accessions(registry: Registry, stage: str, genomes: Optional[List[str]]) -> List[str]:
    """Accessions at ``stage``, restricted to the given genomes when any are named, in first-seen order."""
    if not genomes:
        return query(registry, stage)
    seen: Dict[str, None] = {}
    for genome_id in genomes:
        seen.update(dict.fromkeys(query(registry, stage, genome_id)))
    return list(seen)


def genome_counts(registry: Registry) -> Dict[str, Dict[str, Any]]:
    """The ``"genomes"`` part of ``stage_counts``, from one pass over the recorded extractions.

    ``stage_counts`` runs three passes over every dataset for each genome; this reads each
    dataset's extractions once and gives the same counts and the same ``zero_mapped`` order.
    """
    genomes: Dict[str, Dict[str, Any]] = {
        g: {"extracted": 0, "assembled": 0, "zero_mapped": []} for g in sorted(known_genome_ids(registry))
    }
    for acc, record in registry.datasets.items():
        for genome_id, extraction in (record.get("extractions") or {}).items():
            if not isinstance(extraction, dict):
                continue
            info = genomes[genome_id]
            mapped = extraction.get("mapped_reads") or 0
            if mapped == 0:
                info["zero_mapped"].append(acc)
            info["extracted"] += int(mapped > 0)
            assembly = extraction.get("assembly")
            info["assembled"] += int(isinstance(assembly, dict) and (assembly.get("contigs") or 0) > 0)
    return genomes


def _genome_report(
    registry: Registry,
    paths: ProjectPaths,
    genome_filter: Optional[List[str]],
    counts: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    genome_ids = sorted(known_genome_ids(registry))
    if genome_filter:
        wanted = set(genome_filter)
        genome_ids = [g for g in genome_ids if g in wanted]

    empty_by_genome: Dict[str, List[str]] = {}
    for acc, genome_id in empty_assembly_dirs(paths.targeted, genome_ids):
        empty_by_genome.setdefault(genome_id, []).append(acc)

    report: Dict[str, Any] = {}
    for genome_id in genome_ids:
        info = counts.get(genome_id, {"extracted": 0, "assembled": 0, "zero_mapped": []})
        report[genome_id] = {
            "extracted": info.get("extracted", 0),
            "assembled": info.get("assembled", 0),
            "zero_mapped": info.get("zero_mapped", []),
            "empty_assembly_dirs": sorted(empty_by_genome.get(genome_id, [])),
        }
    return report


def _drift_report(drift: StoreReconcileReport) -> Dict[str, Any]:
    return {
        "recorded_missing": list(drift.recorded_missing),
        "untracked_fastq": list(drift.untracked_fastq),
        "untracked_extractions": [[acc, genome_id] for acc, genome_id in drift.untracked_extractions],
        "empty_assembly_dirs": [[acc, genome_id] for acc, genome_id in drift.empty_assembly_dirs],
        "dangling_links": list(drift.dangling_links),
        # An accession here also appears in dangling_links; its download is left as recorded
        # rather than marked missing, because the store it links into is not mounted at all.
        "store_unavailable": list(drift.store_unavailable),
        # Accessions whose metadata block was filled in from their XML file, and whose
        # unverified/missing download verdict was recomputed, by this reconcile.
        "metadata_filled": list(drift.metadata_filled),
        "verdicts_rechecked": list(drift.verdicts_rechecked),
        # Assembly records ("ACC/GENOME") dated before their extraction, removed by this reconcile.
        "assemblies_dropped": list(drift.assemblies_dropped),
    }


def download_verdicts(registry: Registry) -> Dict[str, List[str]]:
    """Accessions whose recorded completeness verdict is "truncated" or "unverified"."""
    truncated = []
    unverified = []
    for acc in registry.datasets:
        complete = rb.raw(registry, acc, "download", "complete")
        verdict = complete.get("verdict") if isinstance(complete, dict) else None
        if verdict == "truncated":
            truncated.append(acc)
        elif verdict == "unverified":
            unverified.append(acc)
    return {"truncated": sorted(truncated), "unverified": sorted(unverified)}


def _store_report(root: Path) -> Dict[str, Any]:
    """Dataset counts by state from the shared data store's catalogue at ``root``."""
    paths = store_paths(root)
    with Catalog(paths) as catalog:
        rows = catalog.conn.execute(
            "SELECT state, COUNT(*) AS n FROM datasets WHERE state IS NOT 'unknown' GROUP BY state"
        ).fetchall()
    return {"root": str(root), "available": True, "datasets": {row["state"]: row["n"] for row in rows}}


def _resolve_store_root(args, registry) -> Tuple[Optional[Path], bool]:
    """``(root, available)``: the store this project points at, and whether it can be read.

    A status report is about the project, so an unmounted volume or a root that has moved
    must not stop it: the block still names the root, marked unavailable, and everything
    else in the report (dangling links above all, which is exactly what a missing store
    produces) is reported as usual.
    """
    try:
        return resolve_store_root(args.data_root, rb.store_block(registry).root), True
    except DataAccessError as e:
        logger.warning("store unavailable: %s; continuing without it", e)
        return resolve_store_root(args.data_root, rb.store_block(registry).root, require_marker=False), False


def _store_block(root: Path, available: bool) -> Dict[str, Any]:
    """The report's store block: dataset counts when readable, else root and a flag."""
    if not available:
        return {"root": str(root), "available": False, "datasets": {}}
    try:
        return _store_report(root)
    except DataAccessError as e:
        logger.warning("store unavailable: %s; continuing without it", e)
        return {"root": str(root), "available": False, "datasets": {}}


def build_report(
    registry: Registry,
    args: Any,
    paths: ProjectPaths,
    registry_file: Path,
    existed: bool,
    drift: Optional[StoreReconcileReport] = None,
) -> Dict[str, Any]:
    """The status report for ``registry``: local inventory, registry file, store, stages, the
    cross-stage funnel, download verdicts, genomes, drift and timing, as the dict ``status
    --json`` prints.

    ``args`` carries the status command's options (the folders, ``accessions_file``,
    ``parsed_containment``, ``data_root``, ``genome`` and ``init``); ``drift`` is the result of
    ``reconcile`` when ``--reconcile`` ran. Suggested next steps are added by the caller.
    """
    store_root, store_available = _resolve_store_root(args, registry)
    store: Optional[StorePaths] = None
    if store_root is not None and store_available:
        store = store_paths(store_root)
        logger.info("Using shared data store at %s", store_root)

    report = inventory_report(args, registry, store)
    report["registry"] = {
        "path": str(registry_file),
        "version": registry.version,
        "exists": existed or args.init,
        "updated": registry.updated,
    }
    if store_root is not None:
        report["store"] = _store_block(store_root, store_available)
    # One pass over the registry lists every stage; the counts are the lengths of those lists.
    members = stage_members(registry)
    report["stages"] = {s: {"count": len(members[s]), "accessions": members[s]} for s in STAGES}
    report["funnel"] = funnel(registry, members)
    report["downloads"] = download_verdicts(registry)
    report["genomes"] = _genome_report(registry, paths, args.genome, genome_counts(registry))
    report["drift"] = _drift_report(drift) if drift else {}
    report["timing"] = timing_summary(registry)
    return report


def to_dataframes(registry: Registry) -> Tuple["pd.DataFrame", "pd.DataFrame"]:
    """Flat views of the registry for ``status --export-tsv``: one row per accession, and one row
    per (accession, genome) extraction, as two pandas DataFrames."""
    import pandas as pd

    rows = []
    ext_rows = []
    for acc, record in registry.datasets.items():
        screening = rb.screening_block(registry, acc) or rb.ScreeningBlock()
        exclusion = rb.exclusion_block(registry, acc) or rb.ExclusionBlock()
        download = rb.download_block(registry, acc) or rb.DownloadBlock()
        rows.append(
            {
                "accession": acc,
                "screened_genomes": ",".join(sorted(screening.genomes)),
                "selected": bool((rb.selection_block(registry, acc) or rb.SelectionBlock()).selected),
                "excluded": bool(exclusion.excluded),
                "exclusion_reason": exclusion.reason,
                "download_state": download.state,
                "download_date": download.date,
                "bytes_total": download.bytes_total,
                "metadata": "metadata" in record,
                "analyses": ",".join(sorted(record.get("analyses", {}))),
                **rb.profile_summary(registry, acc),
                "download_seconds": download.seconds,
            }
        )
        for genome_id, ext in rb.extraction_blocks(registry, acc).items():
            asm = ext.assembly
            ext_rows.append(
                {
                    "accession": acc,
                    "genome_id": genome_id,
                    "mapped_reads": ext.mapped_reads,
                    "breadth": ext.breadth,
                    "mean_depth": ext.mean_depth,
                    "extraction_date": ext.date or None,
                    "contigs": asm.contigs if asm else None,
                    "total_bp": asm.total_bp if asm else None,
                    "n50": asm.n50 if asm else None,
                    "assembly_date": asm.date if asm else None,
                    "extraction_seconds": ext.seconds,
                    "assembly_seconds": asm.seconds if asm else None,
                    # Appended after the existing columns so earlier exports stay byte-identical;
                    # coverage_tsv is the extraction's own field, the other three come from the
                    # assembly's extra contig stats (n90, largest) and its folder (dir).
                    "coverage_tsv": ext.coverage_tsv,
                    "n90": asm.extra.get("n90") if asm else None,
                    "largest": asm.largest if asm else None,
                    "assembly_dir": asm.dir if asm else None,
                }
            )
    datasets = pd.DataFrame(rows).set_index("accession") if rows else pd.DataFrame()
    extractions = pd.DataFrame(ext_rows)
    return datasets, extractions
