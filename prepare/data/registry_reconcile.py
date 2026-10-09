"""Reconcile the project registry with the disk in two steps: a slow scan, then a fast apply.

``status --reconcile`` walks the FASTQ and targeted folders, counts the reads of every
extraction the registry does not know about and recomputes download completeness verdicts
from the reads on disk. On a large project that takes minutes, far longer than any other
writer should wait for the registry lock. ``scan_reconcile`` therefore does all of that
reading on a snapshot of the registry, without the lock, and returns a ``ReconcilePlan``;
``apply_reconcile`` then applies the plan to the registry loaded inside the lock, checking
each condition again against that registry, and reads no reads. ``reconcile`` runs both
steps on one registry and keeps the behaviour it had before the split.

Before the verdicts, a metadata block recorded without a spot count is filled in from the
project's ``<accession>_metadata.xml`` file, when there is one. A download whose verdict is
missing or ``unverified`` is then checked again once a spot count is known: a plain project
download by counting its mate-1 file and its file of unpaired spots, a link into the shared
store by comparing the read count the store's sidecar records, without reading any FASTQ.
A download whose spot count is still unknown keeps its ``unverified`` verdict. The first
reconcile after upgrading therefore reads the mate-1 file (and any file of unpaired spots)
of every ``unverified`` plain download that has a spot count, which on a large project takes
about as long as reading those files once; later runs find those verdicts known and read
nothing again.

A link into a store that is not mounted (the link's target and its ``sra`` folder are both
absent) is reported in ``StoreReconcileReport.store_unavailable`` and its download is not
marked missing; a dangling link into a mounted store, whose dataset folder is gone, is.

An assembly record dated earlier than its extraction record was built from reads that a later
extraction replaced (0.7.0 carried the old assembly block forward on re-extraction). Reconcile
removes such a record (``clear_assembly``), so ``status`` and ``results_table`` no longer report
it, and lists it in ``StoreReconcileReport.assemblies_dropped``; the assembly folder on disk is
left alone and the next ``extract_target_reads --assemble`` builds the assembly again.

The module sits next to ``prepare.data.registry``, which is held at its current size by
the module size ceiling; that module re-exports ``reconcile``.
"""

import copy
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from prepare.core.exceptions import DataAccessError
from prepare.data import registry as _registry
from prepare.data import registry_blocks as rb
from prepare.data.read_extraction import coverage_table_path, summarise_contigs
from prepare.data.registry import (
    ProjectPaths,
    ReconcileReport,
    Registry,
    _genome_ids_on_disk,
    _infer_coverage,
    empty_assembly_dirs,
    clear_assembly,
    query,
    record_assembly,
    record_download,
    record_extraction,
    scan_assemblies,
    scan_downloads,
    scan_extractions,
    set_download_verdict,
    update_linked,
)
from prepare.data.registry_assembly import assembly_predates_extraction
from prepare.data.sra import accession_has_fastq, count_fastq_reads, verify_download
from prepare.data.sra.spots import expected_spots, merged_verdict

logger = logging.getLogger(__name__)

_CONTIGS_NAME = "final.contigs.fa"
# Returned by ReconcilePlan.value for a figure the scan did not compute.
_NOT_SCANNED = object()
_KNOWN_VERDICTS = ("complete", "truncated")
# Accessions named in the warning about an unmounted store, here and in the status text; the rest
# are counted.
NAMED_IN_WARNING = 5


@dataclass
class StoreReconcileReport(ReconcileReport):
    """``ReconcileReport`` plus what reconcile found about the shared store, metadata and verdicts.

    ``store_unavailable`` lists accessions whose project folder links into a store that is not
    mounted (the link's target and the store's ``sra`` folder are both absent); their downloads
    are left as recorded rather than marked missing. They also appear in ``dangling_links``.
    ``metadata_filled`` lists accessions whose metadata block gained its fields, spot count
    included, from the project's metadata XML file. ``verdicts_rechecked`` lists accessions whose
    ``unverified`` verdict became ``complete`` or ``truncated``. ``assemblies_dropped`` lists, as
    ``"ACC/GENOME"``, the assembly records removed because they were dated earlier than their
    extraction record. Each list is sorted.
    """

    store_unavailable: List[str] = field(default_factory=list)
    metadata_filled: List[str] = field(default_factory=list)
    verdicts_rechecked: List[str] = field(default_factory=list)
    assemblies_dropped: List[str] = field(default_factory=list)


@dataclass
class ReconcilePlan:
    """What ``scan_reconcile`` found on disk, and the figures it computed there.

    The folder listings (``on_disk``, ``extractions``, ``assemblies``) and the report-only
    findings (``empty_assembly_dirs``, ``dangling_links``, and ``unavailable_links``, the dangling
    links into a store that is not mounted) are taken once, by the scan. ``computed`` holds every
    figure that required reading files: the read count and coverage of an untracked extraction,
    the contig statistics of its assembly, the metadata a metadata XML file holds, and the read and
    spot counts behind a download's completeness verdict. Each is keyed by what it depends on (the
    accession, genome id, recorded spot count, download date or store root), so ``apply_reconcile``
    uses a figure only when the registry it is given still asks the same question, and skips it
    otherwise.
    """

    paths: ProjectPaths = field(default_factory=ProjectPaths)
    on_disk: Dict[str, Tuple[int, int]] = field(default_factory=dict)
    genome_ids: List[str] = field(default_factory=list)
    extractions: Dict[str, Dict[str, List[Path]]] = field(default_factory=dict)
    assemblies: Dict[str, Dict[str, Path]] = field(default_factory=dict)
    empty_assembly_dirs: List[Tuple[str, str]] = field(default_factory=list)
    dangling_links: List[str] = field(default_factory=list)
    unavailable_links: List[str] = field(default_factory=list)
    computed: Dict[Tuple[str, ...], Any] = field(default_factory=dict)

    def value(self, key: Tuple[str, ...], compute: Optional[Callable[[], Any]]) -> Any:
        """The figure stored under ``key``; with ``compute`` given and none stored, compute and store it.

        Returns ``_NOT_SCANNED`` when nothing is stored and ``compute`` is None, which is how the
        apply step learns that a condition appeared only after the scan.
        """
        if key not in self.computed:
            if compute is None:
                return _NOT_SCANNED
            self.computed[key] = compute()
        return self.computed[key]


def scan_reconcile(snapshot: Registry, paths: ProjectPaths) -> ReconcilePlan:
    """List the project's folders and compute every figure a reconcile of ``snapshot`` needs.

    Runs without the registry lock and may take minutes (it counts reads). ``snapshot`` is not
    changed: the reconcile is rehearsed on a copy of it, and each figure that rehearsal needs is
    computed and stored in the plan.
    """
    # Imported here rather than at module level: the store package imports the data layer.
    from prepare.store.link import dangling_links

    genome_ids = sorted(_genome_ids_on_disk(paths, snapshot))
    dangling = dangling_links(paths.fastq)
    plan = ReconcilePlan(
        paths=paths,
        on_disk=scan_downloads(paths.fastq),
        genome_ids=genome_ids,
        extractions=scan_extractions(paths.targeted, genome_ids),
        assemblies=scan_assemblies(paths.targeted, genome_ids),
        empty_assembly_dirs=empty_assembly_dirs(paths.targeted, genome_ids),
        dangling_links=dangling,
        unavailable_links=_unmounted_store_links(paths.fastq, dangling),
    )
    _reconcile(copy.deepcopy(snapshot), plan, compute=True)
    return plan


def apply_reconcile(registry: Registry, plan: ReconcilePlan) -> StoreReconcileReport:
    """Apply ``plan`` to ``registry`` in place and return what differed; reads no FASTQ, sidecar or XML.

    Meant to run inside a ``registry_transaction`` on the freshly loaded registry. Every
    condition is checked again against ``registry``: a download recorded since the scan is not
    marked missing twice, an extraction recorded since the scan is not inferred over, metadata
    is filled in only for a block that still lacks a spot count, and a verdict is filled in or
    re-checked only when the registry still records none (or ``unverified``) and the scan
    computed it for the same spot count, download date and store root. Work the scan did not see
    is left for the next reconcile.
    """
    return _reconcile(registry, plan, compute=False)


def reconcile(registry: Registry, paths: ProjectPaths) -> ReconcileReport:
    """Compare the registry with the disk: mark missing downloads, register untracked work.

    Untracked FASTQ and extractions are recorded the way ``bootstrap_from_disk`` records
    them, with ``"inferred": true``, so a project worked on outside PREPARE lands in the
    journal instead of being reported as drift on every run. They stay in the report.

    A link into a shared store whose target has gone (an unmounted store, a dataset removed
    from it) is reported as well, and is not repaired here: removing the link would lose the
    record of which datasets this project uses.

    Equivalent to ``apply_reconcile(registry, scan_reconcile(registry, paths))``; a caller that
    saves the result should run the scan outside the registry lock and only the apply inside it.
    """
    return apply_reconcile(registry, scan_reconcile(registry, paths))


def _unmounted_store_links(project_fastq: Path, dangling: List[str]) -> List[str]:
    """The names in ``dangling`` whose link target's parent folder (the store's ``sra``) is absent too.

    A dangling link whose parent folder exists points into a mounted store that no longer holds
    the dataset; one whose parent folder is gone as well points into a store that is not mounted.
    """
    return [name for name in dangling if not Path(os.path.realpath(project_fastq / name)).parent.is_dir()]


def _reconcile(registry: Registry, plan: ReconcilePlan, compute: bool) -> StoreReconcileReport:
    """The reconcile itself; with ``compute`` False it only uses figures ``plan`` already holds."""
    report = StoreReconcileReport()
    _mark_missing(registry, plan, report, warn=not compute)
    _record_untracked(registry, plan, report)
    for acc, per_genome_files in plan.extractions.items():
        for genome_id, files in per_genome_files.items():
            if rb.extraction_block(registry, acc, genome_id) is None:
                if _infer_extraction(registry, plan, acc, genome_id, files, compute):
                    report.untracked_extractions.append((acc, genome_id))
    report.empty_assembly_dirs = list(plan.empty_assembly_dirs)
    report.dangling_links = list(plan.dangling_links)
    report.store_unavailable = list(plan.unavailable_links)
    report.metadata_filled = sorted(_fill_metadata(registry, plan, compute))
    report.verdicts_rechecked = sorted(_fill_verdicts(registry, plan, compute))
    report.assemblies_dropped = _drop_stale_assemblies(registry)
    return report


def _drop_stale_assemblies(registry: Registry) -> List[str]:
    """Remove each assembly record dated earlier than its extraction record; return them as sorted "ACC/GENOME".

    Uses only the recorded dates (``assembly_predates_extraction``), so it reads no file and the
    apply step needs nothing from the scan. A record without either date is kept.
    """
    dropped = []
    for acc, entry in registry.datasets.items():
        extractions = entry.get("extractions") if isinstance(entry, dict) else None
        for genome_id in list(extractions) if isinstance(extractions, dict) else []:
            if assembly_predates_extraction(registry, acc, genome_id):
                clear_assembly(registry, acc, genome_id)
                dropped.append(f"{acc}/{genome_id}")
    return sorted(dropped)


def _mark_missing(registry: Registry, plan: ReconcilePlan, report: StoreReconcileReport, warn: bool) -> None:
    """Mark downloaded accessions whose reads are gone as missing, except links into an unmounted store.

    ``warn`` logs one warning naming the downloads left alone because their store is not mounted;
    the scan's rehearsal passes False so the warning is logged once per reconcile.
    """
    unavailable = set(plan.unavailable_links)
    left_alone: List[str] = []
    for acc in registry.datasets:
        download = rb.download_block(registry, acc)
        if download is None or download.state != "downloaded" or acc in plan.on_disk:
            continue
        if acc in unavailable:
            left_alone.append(acc)
            continue
        # A download another process finished after the scan is not in the scan's listing;
        # one look at its folder (a few stat calls) keeps it from being marked missing.
        if accession_has_fastq(plan.paths.fastq / acc):
            continue
        download.state, download.date = "missing", _registry._now()
        rb.set_download_block(registry, acc, download)
        report.recorded_missing.append(acc)
    if warn and left_alone:
        named = ", ".join(sorted(left_alone)[:NAMED_IN_WARNING])
        more = len(left_alone) - NAMED_IN_WARNING
        logger.warning(
            "%d downloaded dataset(s) link into a data store that is not mounted and were not marked missing: %s%s",
            len(left_alone),
            named,
            f" and {more} more" if more > 0 else "",
        )


def _links_into_store(registry: Registry, entry: Path) -> bool:
    """True when the project entry ``entry`` is a symlink, into the registry's store when one is recorded."""
    if not entry.is_symlink():
        return False
    root = rb.store_block(registry).root
    if not root:
        return True
    # Imported here, not at module level: prepare.store imports the data registry.
    from prepare.store.layout import store_paths
    from prepare.store.link import is_store_link

    return is_store_link(entry, store_paths(Path(root)))


def _record_untracked(registry: Registry, plan: ReconcilePlan, report: StoreReconcileReport) -> None:
    """Record FASTQ found on disk that the registry does not record as downloaded, marked inferred.

    An accession whose previous download record says it came from the shared store (e.g. one an
    earlier reconcile marked missing while the store was unmounted) is recorded as linked from the
    store again, with its ``store_name``, and put back on the registry's list of linked datasets,
    provided its project entry is still a symlink into the store; a real folder that replaced the
    link is recorded as a plain download.
    """
    tracked = set(query(registry, "downloaded"))
    # Checked on disk again, the mirror of the check in _mark_missing: an accession unlinked or
    # removed after the scan (store_unlink records it missing) must not be written back as downloaded.
    untracked = (acc for acc in plan.on_disk if acc not in tracked)
    report.untracked_fastq = sorted(acc for acc in untracked if accession_has_fastq(plan.paths.fastq / acc))
    for acc in report.untracked_fastq:
        previous = rb.download_block(registry, acc)
        if previous is not None and previous.source == "store" and _links_into_store(registry, plan.paths.fastq / acc):
            store_name = previous.store_name or acc
            record_download(
                registry, acc, "downloaded", plan.paths.fastq, attempt=False, source="store", store_name=store_name
            )
            update_linked(registry, acc, add=True)
        else:
            record_download(registry, acc, "downloaded", plan.paths.fastq, attempt=False)
        rb.mark_inferred(registry, acc, "download", attempts=0)


def _count_extraction(files: List[Path], genome_id: str) -> Tuple[int, Optional[Dict[str, Any]]]:
    """Reads in every file of one extraction, and the coverage its table records (or None)."""
    reads = sum(count_fastq_reads(f) for f in files)
    coverage = _infer_coverage(coverage_table_path(files[0].parent, genome_id)) if files else None
    return reads, coverage


def _infer_extraction(
    registry: Registry, plan: ReconcilePlan, acc: str, genome_id: str, files: List[Path], compute: bool
) -> bool:
    """Record one untracked extraction (and its assembly), marked inferred; False when not scanned.

    The read count covers both mates, singles and unpaired reads, as ``bootstrap_from_disk``
    counts them; the assembly is recorded only when it holds contigs.
    """
    counted = plan.value(
        ("extraction", acc, genome_id), (lambda: _count_extraction(files, genome_id)) if compute else None
    )
    if counted is _NOT_SCANNED:
        return False
    reads, coverage = counted
    record_extraction(registry, acc, genome_id, files, reads, False, {}, coverage=coverage)
    rb.mark_inferred(registry, acc, "extractions", genome_id=genome_id)
    asm_dir = plan.assemblies.get(acc, {}).get(genome_id)
    if asm_dir is not None:
        contigs = asm_dir / _CONTIGS_NAME
        stats = plan.value(("assembly", acc, genome_id), (lambda: summarise_contigs(contigs)) if compute else None)
        if stats is not _NOT_SCANNED and stats["contigs"] > 0:
            record_assembly(registry, acc, genome_id, asm_dir, stats, "", {})
            rb.mark_inferred(registry, acc, "extractions", genome_id=genome_id)
    return True


def _fill_missing_download_verdicts(registry: Registry, paths: ProjectPaths) -> None:
    """Compute a completeness verdict for a downloaded accession that has none, or ``unverified``.

    A project downloaded before completeness verification existed (or with
    ``--no-verify-downloads``) has metadata recorded but no ``download.complete`` verdict, and a
    download made while NCBI's spot count was unknown records ``unverified``. Once the spot count
    is on file, this recomputes the verdict the same way a fresh download would have, against the
    accession's files on disk. A record whose reads came from the shared store is filled in from
    that dataset's sidecar instead: the store counted the reads when it downloaded them, and
    counting them again through a link would repeat work another project has done.
    """
    _fill_verdicts(registry, ReconcilePlan(paths=paths), compute=True)


def _fill_metadata(registry: Registry, plan: ReconcilePlan, compute: bool) -> List[str]:
    """Fill metadata blocks recorded without a spot count from the project's XML files; return the accessions.

    With ``compute`` (the scan) this is ``fill_metadata_from_xml`` on the project's metadata
    folder, and each block it fills is kept in the plan. Without it (the apply) no XML is read:
    a block the scan filled is copied in, only while the registry's block still lacks a spot
    count, keeping that block's ``inferred`` mark.
    """
    if compute:
        # Imported here: metadata_fields pulls in the XML parsing stack, which reconcile otherwise
        # does not need at import time.
        from prepare.data.metadata_fields import fill_metadata_from_xml

        filled = fill_metadata_from_xml(registry, plan.paths.metadata)
        for acc in filled:
            plan.computed[("metadata", acc)] = copy.deepcopy(registry.datasets[acc]["metadata"])
        return filled
    filled = []
    for acc in registry.datasets:
        block = rb.metadata_block(registry, acc)
        if block is None or block.run_total_spots is not None:
            continue
        stored = plan.value(("metadata", acc), None)
        if stored is _NOT_SCANNED:
            continue
        new_block = rb.MetadataBlock.from_dict(copy.deepcopy(stored))
        if block.inferred:
            new_block.inferred = True
        rb.set_metadata_block(registry, acc, new_block)
        filled.append(acc)
    return filled


def _fill_verdicts(registry: Registry, plan: ReconcilePlan, compute: bool) -> List[str]:
    """``_fill_missing_download_verdicts`` on a plan; returns the accessions whose ``unverified`` verdict became known.

    With ``compute`` False, only figures the scan computed are used. A previous ``unverified``
    verdict is replaced only by a known one (``complete`` or ``truncated``); a download whose
    store is not mounted is skipped.
    """
    rechecked: List[str] = []
    unavailable = set(plan.unavailable_links)
    for acc in registry.datasets:
        download = rb.download_block(registry, acc) or rb.DownloadBlock()
        # A verdict recorded as an empty dict counts as none, as it did before the typed blocks.
        previous = download.complete.to_dict() if download.complete is not None else {}
        unverified = bool(previous) and previous.get("verdict") == "unverified"
        if download.state != "downloaded" or (previous and not unverified) or acc in unavailable:
            continue
        if download.source == "store":
            new = _store_recheck(registry, plan, acc, download, compute)
        else:
            new = _plain_recheck(registry, plan, acc, download, compute)
        if new is None or new == previous:
            continue
        known = new.get("verdict") in _KNOWN_VERDICTS
        if unverified and not known:
            continue
        set_download_verdict(registry, acc, new)
        if unverified:
            rechecked.append(acc)
    return rechecked


def _plain_recheck(
    registry: Registry, plan: ReconcilePlan, acc: str, download: rb.DownloadBlock, compute: bool
) -> Optional[Dict[str, Any]]:
    """The verdict of a plain project download from its reads on disk, or None when it cannot be had.

    Counts the mate-1 (or single-end) file and the file of unpaired spots once per download; the
    count is keyed on the spot count and the download's date, so files downloaded again between
    the scan and the apply are never stamped with the verdict computed on the files they replaced.
    """
    spots = expected_spots(registry, acc)
    if spots is None:
        return None
    acc_dir = plan.paths.fastq / acc
    if not acc_dir.is_dir():
        return None
    reads = plan.value(
        ("spots", acc, str(spots), str(download.date)),
        (lambda: verify_download(acc, acc_dir, spots)["reads_r1"]) if compute else None,
    )
    if reads is _NOT_SCANNED:
        return None
    return merged_verdict(download.complete, None, reads, spots)


def _store_recheck(
    registry: Registry, plan: ReconcilePlan, acc: str, download: rb.DownloadBlock, compute: bool
) -> Optional[Dict[str, Any]]:
    """The verdict of a download linked from the store, from its sidecar's counts; reads no FASTQ.

    The sidecar's ``reads_per_mate`` is compared with the expected spot count (the registry's
    metadata first, then the sidecar's ``ncbi.spots``). Without both counts the sidecar's own
    verdict is taken, as before re-checking existed. Keyed on the store root, the registry's spot
    count and the download's date.
    """
    root = rb.store_block(registry).root or ""
    recorded_spots = (rb.metadata_block(registry, acc) or rb.MetadataBlock()).run_total_spots
    figures = plan.value(
        ("store", acc, root, str(recorded_spots), str(download.date)),
        (lambda: _store_figures(registry, acc)) if compute else None,
    )
    if figures is _NOT_SCANNED or figures is None:
        return None
    sidecar_verdict, reads, spots = figures
    return merged_verdict(download.complete, sidecar_verdict, reads, spots)


def _store_figures(registry: Registry, accession: str) -> Optional[Tuple[Dict[str, Any], Optional[int], Optional[int]]]:
    """The store sidecar's verdict for ``accession``, its read count per mate and the expected spot count.

    Reads the store root the registry itself recorded; a project whose store has moved or is
    not mounted simply gets no verdict this time round (None), exactly as before.
    """
    root = rb.store_block(registry).root
    if not root:
        return None
    # Imported here, not at module level: prepare.store imports the data registry.
    from prepare.store.layout import sidecar_path, store_paths
    from prepare.store.sidecar import sidecar_completeness

    try:
        store = store_paths(Path(root))
        verdict = sidecar_completeness(sidecar_path(store, accession))
        if verdict is None:
            return None
        return verdict, verdict.get("reads_r1"), expected_spots(registry, accession, store=store)
    except (OSError, DataAccessError) as e:
        logger.warning("Could not read the store sidecar for %s: %s", accession, e)
        return None
