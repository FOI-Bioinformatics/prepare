"""Per-project registry of dataset state, the project journal.

The registry records decisions and provenance for every SRA accession a
project touches: screening results per target genome, the selection criteria,
exclusions with reasons, download outcomes with file sizes and dates, analyses,
and per-genome extraction and assembly results. Existence on disk is never
taken from the registry alone; ``status`` re-checks the filesystem and the
scanners here rebuild or reconcile the journal from what is on disk.
"""

import json
import logging
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple, Union

from prepare.core.constants import DEFAULT_REGISTRY_MAX_SCREENED, GENOME_FASTA_GLOBS
from prepare.core.constants import SHORT_LOCK_HEARTBEAT_SECONDS, SHORT_LOCK_POLL_SECONDS
from prepare.core.exceptions import DataAccessError
from prepare.core.settings import settings_or
from prepare.data.file_io import visible_files, write_text_atomic
from prepare.data import registry_blocks as rb
from prepare.data.read_extraction import coverage_table_path, summarise_contigs, summarise_coverage_table
from prepare.data.sra import accession_has_fastq, count_fastq_reads, fastq_files, is_transient_folder
from prepare.utils.lockfile import LockPolicy, held_lock, verify_held

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)

REGISTRY_FILENAME = "prepare_registry.json"
# Version 1: version, created, updated, genomes, datasets.
# Version 2 adds project (this project's identity, recorded by `store_init`) and store
# (the shared data store this project is bound to). Both default to {} so a version-1
# file loads unchanged; the next save writes it back as version 2.
SCHEMA_VERSION = 2
STAGES = ("screened", "selected", "excluded", "downloaded", "analysed", "extracted", "assembled")
# Lock defaults (repeated in core.settings), read at call time; settings override; the wait fits SLURM KillWait.
LOCK_STALE_SECONDS = 120.0
LOCK_WAIT_SECONDS = 30.0
_MATE_SUFFIXES = ("_1", "_2", "_s", "_0")
_ASSEMBLY_SUFFIX = "_assembly"
_CONTIGS_NAME = "final.contigs.fa"


def _now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


@dataclass
class ProjectPaths:
    """Folders of one project, relative to the project root unless absolute."""

    fastq: Path = Path("fastq")
    metadata: Path = Path("metadata")
    genomes: Path = Path("genomes")
    targeted: Path = Path("targeted")
    matches: Path = Path("matches")


@dataclass
class Registry:
    """In-memory form of one project's registry file: schema version, timestamps, and state.

    Holds the per-genome screening thresholds (``genomes``), the per-accession dataset records
    (``datasets``, keyed by accession, each carrying its stage history and provenance), this
    project's own identity (``project``) once bound by ``store_init``, and the shared data store
    it is linked to, if any (``store``). ``path`` is set once the registry is bound to a file on
    disk (by ``load_registry`` or the first ``save_registry``); it is ``None`` for a registry
    built only in memory, e.g. mid-``bootstrap_from_disk``, and callers that need a project root
    from it should go through ``project_root`` rather than reading ``path`` directly.
    """

    version: int = SCHEMA_VERSION
    created: str = field(default_factory=_now)
    updated: str = field(default_factory=_now)
    genomes: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    datasets: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    # This project's identity ({"id", "name", "path", "created"}), recorded by `store_init`.
    project: Dict[str, Any] = field(default_factory=dict)
    # The shared data store this project is bound to ({"root", "mode", "linked"}).
    store: Dict[str, Any] = field(default_factory=dict)
    path: Optional[Path] = None


def project_root(registry: Registry) -> Path:
    """The project root a registry's paths are recorded relative to: its file's parent folder.

    A registry not yet bound to a file (before its first save, e.g. mid-``bootstrap_from_disk``)
    falls back to the current working directory, so a path recorded at that point resolves the
    same way once the registry is later bound and saved from the same directory.
    """
    if registry.path is not None:
        return registry.path.parent.resolve()
    return Path.cwd().resolve()


def _project_relative(path: Union[str, Path], root: Path) -> str:
    """The form to record for ``path``: relative to ``root`` when inside it, else absolute.

    Keeps a project's registry movable: renaming or relocating the project directory does not
    break paths recorded under it. A path outside the project root (for example a genome FASTA
    shared from elsewhere on disk) has no project-relative form, so it is kept absolute. A path
    under a symlinked folder resolves to the real (target) location, so a project folder linked
    into an external store records that store's absolute path here; the store-relative name used
    to look the dataset up in the store is recorded separately by the store-linking feature.
    """
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError:
        return str(resolved)


def resolve_project_path(registry: Registry, value: Union[str, Path]) -> Path:
    """Resolve one path recorded in the registry against its project root.

    An absolute recorded value is returned unchanged, so entries written by a PREPARE
    version before this change (always absolute) keep resolving correctly. A relative value
    is joined to ``project_root``. For a registry written by an older version, whose relative
    paths were relative to the working directory of that run, this still resolves correctly in
    the common case, since that working directory is normally where the registry file itself
    lived, which is exactly the project root computed here. Callers with a value that may be
    ``None`` (an optional recorded field) must check that themselves before calling.
    """
    path = Path(value)
    if path.is_absolute():
        return path
    return project_root(registry) / path


@dataclass
class ReconcileReport:
    """Differences ``reconcile`` found between the registry and the project's filesystem.

    ``recorded_missing`` lists accessions the registry marks downloaded whose FASTQ files are no
    longer present; ``untracked_fastq`` lists accessions with FASTQ on disk that the registry does
    not yet record as downloaded. ``untracked_extractions`` and ``empty_assembly_dirs`` are lists
    of ``(accession, genome_id)`` pairs found on disk but missing from, respectively, the
    registry's extraction and assembly records. ``dangling_links`` holds accessions whose project
    folder is a symlink into a shared data store that is unmounted or has lost its copy, so the
    link resolves to nothing.
    """

    recorded_missing: List[str] = field(default_factory=list)
    untracked_fastq: List[str] = field(default_factory=list)
    untracked_extractions: List[Tuple[str, str]] = field(default_factory=list)
    empty_assembly_dirs: List[Tuple[str, str]] = field(default_factory=list)
    # Accessions whose project folder is a symlink with nothing at the other end, i.e. a
    # dataset the project reads from a shared store that is unmounted or has lost the copy.
    dangling_links: List[str] = field(default_factory=list)


# ----------------------------------------------------------------- persistence


def registry_path(explicit: Optional[Union[str, Path]] = None, start: Union[str, Path] = ".") -> Path:
    """The registry file to use: an explicit path, else the nearest one walking up from ``start``."""
    if explicit:
        return Path(explicit)
    current = Path(start).resolve()
    for folder in (current, *current.parents):
        candidate = folder / REGISTRY_FILENAME
        if candidate.exists():
            if folder != current:
                logger.info("Using the project registry found above the working directory: %s", candidate)
            return candidate
    return current / REGISTRY_FILENAME


def load_registry(path: Optional[Union[str, Path]] = None) -> Registry:
    """Load the registry, or an empty one bound to the path when the file does not exist."""
    target = registry_path(path)
    if not target.exists():
        return Registry(path=target)
    try:
        text = target.read_text()
    except OSError as e:
        raise DataAccessError(f"Cannot read registry {target}: {e}") from e
    try:
        data = json.loads(text)
    except json.JSONDecodeError as e:
        raise DataAccessError(f"Registry is not valid JSON: {target} ({e})") from e
    return Registry(
        version=int(data.get("version", SCHEMA_VERSION)),
        created=str(data.get("created", _now())),
        updated=str(data.get("updated", _now())),
        genomes=dict(data.get("genomes", {})),
        datasets=dict(data.get("datasets", {})),
        project=dict(data.get("project", {})),
        store=dict(data.get("store", {})),
        path=target,
    )


@contextmanager
def _acquire_lock(lock: Path) -> Iterator[Path]:
    """Hold the registry lock file ``lock`` for a with-block (``prepare.utils.lockfile``)."""
    stale, wait = settings_or(registry_lock_stale=LOCK_STALE_SECONDS, registry_lock_wait=LOCK_WAIT_SECONDS)
    limits = (stale, wait, SHORT_LOCK_POLL_SECONDS, SHORT_LOCK_HEARTBEAT_SECONDS)
    with held_lock(lock, LockPolicy("Registry", *limits)) as held:
        yield held


def _write_registry(registry: Registry, target: Path, lock: Path) -> Path:
    """Write the registry to ``target`` atomically (temp file, compact JSON, rename) while ``lock`` is still held."""
    verify_held(lock)  # LockLost (a DataAccessError) if a stall let another holder reclaim it
    target.parent.mkdir(parents=True, exist_ok=True)
    registry.updated = _now()
    registry.path = target
    # A registry of an older schema is upgraded on save; there is no migration step, since
    # every field new schema versions add defaults to {} when missing.
    registry.version = SCHEMA_VERSION
    payload = {
        "version": registry.version,
        "created": registry.created,
        "updated": registry.updated,
        "genomes": registry.genomes,
        "datasets": registry.datasets,
        "project": registry.project,
        "store": registry.store,
    }
    # write_text_atomic removes its temporary file on any exit that did not replace the target,
    # KeyboardInterrupt included, and fsync keeps a power loss from leaving an empty registry.
    try:
        return write_text_atomic(target, json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n", fsync=True)
    except OSError as e:
        raise DataAccessError(f"Cannot write registry {target}: {e}") from e


def save_registry(registry: Registry, path: Optional[Union[str, Path]] = None) -> Path:
    """Write an already loaded registry atomically under a lock file."""
    target = Path(path) if path else (registry.path or registry_path())
    target.parent.mkdir(parents=True, exist_ok=True)
    with _acquire_lock(target.with_name(target.name + ".lock")) as lock:
        return _write_registry(registry, target, lock)


@contextmanager
def registry_transaction(path: Optional[Union[str, Path]] = None) -> Iterator[Registry]:
    """Load, mutate and save the registry under one lock, so a concurrent edit is never reverted.

    Not re-entrant: nesting a second call to this function (or to `save_registry`) for the
    same registry file inside this block's body raises ``LockReentry`` at once.

    Use this instead of holding one loaded registry across a long run: the file is
    read inside the lock and written back at the end of the block. If the block
    raises, the lock is released and nothing is written.
    """
    target = registry_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with _acquire_lock(target.with_name(target.name + ".lock")) as lock:
        registry = load_registry(target)
        yield registry
        _write_registry(registry, target, lock)


# --------------------------------------------------------------------- writers


def upsert_dataset(registry: Registry, accession: str) -> Dict[str, Any]:
    """Return ``accession``'s dataset record, creating an empty one in ``registry.datasets`` first if absent."""
    return registry.datasets.setdefault(accession, {})


def record_screening(
    registry: Registry,
    accession: str,
    genome_id: str,
    containment: float,
    cani: Optional[float],
    source: str,
    query_threshold: float,
    csv_path: Optional[Union[str, Path]],
) -> None:
    """Record one genome's containment (and, when available, cANI) screening result for one accession.

    Rounds ``containment`` and ``cani`` to four decimal places, clears any earlier "inferred"
    flag on the accession's screening block (a real screening result supersedes one carried over
    from another run), and registers ``genome_id`` in ``registry.genomes`` if it is not already
    known. Returns nothing; raises nothing of its own.
    """
    screening = rb.screening_block(registry, accession) or rb.ScreeningBlock()
    screening.date = _now()
    screening.discard("inferred")
    screening.genomes[genome_id] = rb.ScreeningEntry(
        containment=round(float(containment), 4),
        cani=round(float(cani), 4) if cani is not None else None,
        csv=str(csv_path) if csv_path else None,
        source=source,
        query_threshold=query_threshold,
    )
    rb.set_screening_block(registry, accession, screening)
    registry.genomes.setdefault(genome_id, {})


def record_genome(registry: Registry, genome_id: str, fasta: Union[str, Path], manifest: Union[str, Path]) -> None:
    """Record where a target genome's FASTA lives, and the manifest it came from.

    ``manifest`` is empty (falsy) for a genome fetched directly (e.g. by ``genome_download``,
    with no manifest CSV written); stored as "" rather than resolved, so it does not turn
    into a meaningless path built from wherever the process happened to run.
    """
    root = project_root(registry)
    registry.genomes[genome_id] = {
        "fasta": _project_relative(fasta, root),
        "manifest": _project_relative(manifest, root) if manifest else "",
        "date": _now(),
    }


def record_selection(
    registry: Registry,
    accessions: Sequence[str],
    criteria: Dict[str, Any],
    output: Union[str, Path],
    ranked: Optional[List[Dict[str, Any]]] = None,
) -> None:
    """Mark ``accessions`` selected with ``criteria``; anything selected earlier but absent now becomes unselected.

    ``ranked`` is the per-accession rank detail (``{"accession", "rank", "column", "value"}``)
    for the selected accessions; entries for accessions outside ``accessions`` are dropped so
    the recorded list is always capped to what was actually selected.
    """
    chosen = set(accessions)
    for accession in registry.datasets:
        earlier = rb.selection_block(registry, accession)
        if earlier and earlier.selected and accession not in chosen:
            earlier.selected = False
            earlier.date = _now()
            rb.set_selection_block(registry, accession, earlier)
    ranked_by_accession = {entry["accession"]: entry for entry in ranked or [] if entry.get("accession") in chosen}
    for accession in accessions:
        selection = rb.SelectionBlock(selected=True, date=_now(), criteria=dict(criteria), output=str(output))
        if accession in ranked_by_accession:
            selection.ranked = [ranked_by_accession[accession]]
        rb.set_selection_block(registry, accession, selection)


def cap_screening(registry: Registry, genome_id: str, max_screened: int = DEFAULT_REGISTRY_MAX_SCREENED) -> int:
    """Keep only the ``max_screened`` highest containments for one genome; return how many were dropped.

    A broad search can match tens of thousands of metagenomes. The CSVs remain the raw
    record; the registry keeps the best matches so it stays small enough to rewrite after
    every completed accession.
    """
    screened = {acc: block for acc in registry.datasets if (block := rb.screening_block(registry, acc)) is not None}
    entries = [
        (acc, float(block.genomes[genome_id].containment or 0.0))
        for acc, block in screened.items()
        if genome_id in block.genomes
    ]
    if len(entries) <= max_screened:
        return 0
    for accession, _ in sorted(entries, key=lambda item: item[1], reverse=True)[max_screened:]:
        record, block = registry.datasets[accession], screened[accession]
        del block.genomes[genome_id]
        if block.genomes:
            rb.set_screening_block(registry, accession, block)
        else:
            del record["screening"]
        if not record:
            del registry.datasets[accession]
    dropped = len(entries) - max_screened
    logger.warning(
        "Kept the %d highest containments for %s in the registry and dropped %d more "
        "(the limit is --registry-max-screened); the match CSV keeps them all",
        max_screened,
        genome_id,
        dropped,
    )
    return dropped


def record_screening_from_table(
    registry: Registry,
    table: Union["pd.DataFrame", str, Path],
    matches_folder: Union[str, Path],
    max_screened: int = DEFAULT_REGISTRY_MAX_SCREENED,
) -> int:
    """Record a screening entry for every positive containment in a parsed containment table.

    ``table`` is the table ``parse_containment_data`` built (a DataFrame, as its summary's
    ``table``) or the path it was written to. For every column except ``max_containment`` and
    ``max_containment_annotation`` (each a genome), records one screening entry per row with a
    value greater than 0, keeping at most ``max_screened`` accessions per genome. Cells that do
    not hold a number are skipped; every block written gets one timestamp. Returns the number of
    entries recorded. If a table path does not exist, logs at debug level and returns 0.
    """
    from prepare.data.screening_table import record_screening_table

    return record_screening_table(registry, table, matches_folder, max_screened)


def record_exclusion(registry: Registry, accession: str, reason: str, source: str = "user") -> None:
    """Mark ``accession`` excluded with ``reason`` and ``source``, timestamped; overwrites any earlier exclusion."""
    block = rb.ExclusionBlock(excluded=True, reason=reason, source=source, date=_now())
    rb.set_exclusion_block(registry, accession, block)


def clear_exclusion(registry: Registry, accession: str) -> None:
    """Reverse a recorded exclusion for ``accession``, if any; a no-op when none was recorded.

    Sets the exclusion block's ``excluded`` flag to ``False`` and blanks its reason rather than
    deleting the block, so the accession keeps a record of having once been excluded.
    """
    if rb.exclusion_block(registry, accession) is not None:
        cleared = rb.ExclusionBlock(excluded=False, reason="", source="user", date=_now())
        rb.set_exclusion_block(registry, accession, cleared)


def _file_entries(paths: Iterable[Path], root: Path) -> List[rb.FileEntry]:
    entries = []
    for path in paths:
        stat = path.stat()
        mtime = datetime.fromtimestamp(stat.st_mtime).astimezone().isoformat(timespec="seconds")
        entries.append(rb.FileEntry(path=_project_relative(path, root), bytes=stat.st_size, mtime=mtime))
    return entries


def record_download(
    registry: Registry,
    accession: str,
    state: str,
    fastq_dir: Union[str, Path],
    message: str = "",
    attempt: bool = True,
    complete: Optional[Dict[str, Any]] = None,
    source: Optional[str] = None,
    store_name: Optional[str] = None,
) -> None:
    """Record a download outcome; ``state`` is downloaded, failed, missing or skipped.

    ``attempt`` counts this call against ``attempts``; pass ``False`` when recording a
    state without an actual download attempt (e.g. a file found already present on disk).
    ``complete`` is the completeness verdict from ``prepare.data.sra.verify_download``
    (via ``parse_verdict_message``); when omitted, any verdict already on file is left as is.
    When given, it replaces the block on file, except that a read count it carries as ``None``
    (e.g. a store sidecar with no read count of its own) keeps the previous count if, and only
    if, the previous verdict equals the new one; see ``Verdict.carrying_counts_from``.
    ``source`` says where the reads came from (``"store"`` for a dataset the shared store
    holds and the project only links to) and ``store_name`` is the dataset's name inside
    that store. Both describe this outcome, so a call that names neither clears whatever
    an earlier outcome recorded rather than leaving a stale claim behind.
    """
    download = rb.download_block_for_write(registry, accession)
    if attempt and state in ("downloaded", "failed"):
        download.attempts = int(download.attempts or 0) + 1
    files: List[rb.FileEntry] = []
    if state == "downloaded":
        files = _file_entries(fastq_files(Path(fastq_dir) / accession), project_root(registry))
    download.state, download.date, download.files, download.message = state, _now(), files, message
    download.bytes_total = sum(int(f.bytes) for f in files)
    if complete is not None:
        download.complete = rb.Verdict.from_dict(complete).carrying_counts_from(download.complete)
    if source is None:
        download.discard("source")
        download.discard("store_name")
    else:
        download.source = source
        if store_name is not None:
            download.store_name = store_name
    # The mate read counts cached by extract_target_reads describe the files this entry
    # replaces, so they cannot survive a new attempt or state. Keeping them would let a pair
    # counted from a truncated download force single-end mapping of the complete one.
    for stale in ("inferred", "mate_reads", "mate_reads_signature"):
        download.discard(stale)
    rb.set_download_block(registry, accession, download)


def set_download_verdict(registry: Registry, accession: str, verdict: Dict[str, Any]) -> None:
    """Set only a downloaded accession's completeness verdict, touching nothing else.

    Unlike ``record_download``, this never resets ``date``, re-scans ``files``/``bytes_total``,
    or clears ``message``/``source``/``store_name``; use it when only the completeness verdict
    needs to change, e.g. computing one that a download recorded before verification existed
    never got.
    """
    download = rb.download_block_for_write(registry, accession)
    download.complete = rb.Verdict.from_dict(verdict)
    rb.set_download_block(registry, accession, download)


def update_linked(registry: Registry, accession: str, add: bool) -> None:
    """Add ``accession`` to, or remove it from, the registry's list of datasets this project links
    from the store (``registry.store["linked"]``).

    The list is kept sorted and free of duplicates, so calling this twice with the same arguments
    leaves the registry as one call did.
    """
    store = rb.store_block(registry)
    linked = set(store.linked or [])
    if add:
        linked.add(accession)
    else:
        linked.discard(accession)
    store.linked = sorted(linked)
    rb.set_store_block(registry, store)


def nan_to_none(value: Any) -> Any:
    """Return ``None`` for a pandas NaN/NA value, else ``value`` unchanged."""
    import pandas as pd

    return None if pd.isna(value) else value


def _to_int_or_none(value: Any) -> Optional[int]:
    """Convert a numeric value (int or numeric string) to ``int``; ``None`` otherwise."""
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


# Public name for callers outside this module (e.g. the results table).
to_int_or_none = _to_int_or_none


def record_metadata(
    registry: Registry, accession: str, xml_path: Union[str, Path], fields: Dict[str, Any], root: Optional[Path] = None
) -> None:
    """Record ``accession``'s downloaded NCBI metadata under its dataset entry, replacing any earlier record.

    Copies a fixed set of fields out of ``fields`` (run size and md5, assay type, organism,
    collection date, library layout, platform, library strategy, and the total spot/base counts,
    coerced to ``int`` or ``None``), alongside the project-relative path to the metadata XML and
    a timestamp. Fields absent from ``fields`` are recorded as ``None`` rather than omitted.
    ``root`` defaults to ``project_root(registry)``; a caller recording many accessions passes it.
    """
    root = project_root(registry) if root is None else root
    block = rb.MetadataBlock(xml=_project_relative(xml_path, root), date=_now())
    for key in (
        "run_size",
        "run_md5",
        "assay_type",
        "organism",
        "collection_date",
        "library_layout",
        "platform",
        "library_strategy",
    ):
        setattr(block, key, fields.get(key))
    block.run_total_spots = _to_int_or_none(fields.get("run_total_spots"))
    block.run_total_bases = _to_int_or_none(fields.get("run_total_bases"))
    rb.set_metadata_block(registry, accession, block)


def record_analysis(
    registry: Registry, accession: str, analysis: str, output: Union[str, Path], summary: Dict[str, Any]
) -> None:
    """Record one named analysis's output path, summary and timestamp for ``accession``.

    Stored under the dataset entry's ``analyses`` mapping, keyed by ``analysis``, so a second
    call with the same name overwrites that analysis's previous run rather than accumulating a
    history of runs.
    """
    output_path = _project_relative(output, project_root(registry))
    entry = rb.AnalysisEntry(date=_now(), output=output_path, summary=dict(summary))
    upsert_dataset(registry, accession).setdefault("analyses", {})[analysis] = entry.to_dict()


def record_export(registry: Registry, name: str, output: Union[str, Path], summary: Dict[str, Any]) -> None:
    """Record a project-level export (e.g. the results table) under ``project["exports"][name]``.

    Only the latest run of each export is kept: its date, the project-relative output path
    and a summary of what it contained.
    """
    project = rb.project_block(registry)
    output_path = _project_relative(output, project_root(registry))
    project.exports[name] = rb.ExportEntry(date=_now(), output=output_path, summary=dict(summary))
    rb.set_project_block(registry, project)


def record_extraction(
    registry: Registry,
    accession: str,
    genome_id: str,
    files: Sequence[Path],
    mapped_reads: int,
    unequal_mates: bool,
    params: Dict[str, Any],
    mapped_total: Optional[int] = None,
    coverage: Optional[Dict[str, Any]] = None,
) -> None:
    """Record one sample's extraction against one genome.

    ``coverage`` is ``ExtractionResult.coverage``: its ``breadth`` and ``mean_depth`` are
    recorded as given and its ``coverage_tsv`` project-relative; all three are None when
    ``coverage`` is None (nothing mapped, or the coverage step failed). A recorded assembly is
    dropped: it was built from the reads this call replaces.
    """
    root = project_root(registry)
    genome_fasta = params.get("genome_fasta")
    coverage = coverage or {}
    coverage_tsv = coverage.get("coverage_tsv")
    block = rb.ExtractionBlock(
        date=_now(),
        genome_fasta=_project_relative(genome_fasta, root) if genome_fasta is not None else None,
        preset=params.get("preset"),
        threshold=params.get("threshold"),
        filter_flags=params.get("filter_flags"),
        min_mapq=params.get("min_mapq"),
        index=params.get("index"),
        mapped_reads=int(mapped_reads),
        mapped_total=int(mapped_total) if mapped_total is not None else None,
        unequal_mates=bool(unequal_mates),
        files=[_project_relative(p, root) for p in files],
        breadth=coverage.get("breadth"),
        mean_depth=coverage.get("mean_depth"),
        coverage_tsv=_project_relative(coverage_tsv, root) if coverage_tsv is not None else None,
        assembly=None,
    )
    rb.set_extraction_block(registry, accession, genome_id, block)
    registry.genomes.setdefault(genome_id, {})


# The four assembly stats every assembly block is guaranteed to carry, forced to int so a
# caller can always rely on their type regardless of what ``stats`` provides.
_REQUIRED_ASSEMBLY_STATS = ("contigs", "total_bp", "n50", "largest")


def record_assembly(
    registry: Registry,
    accession: str,
    genome_id: str,
    assembly_dir: Union[str, Path],
    stats: Dict[str, Any],
    tool_version: str,
    params: Dict[str, Any],
) -> None:
    """Record one assembly's stats, provenance and parameters.

    ``stats`` is recorded in full (contig-level metrics such as N90, GC content and
    contigs_ge_1kb, and read-mapping metrics such as reads_mapped, mapping_rate,
    mean_depth_estimate and genome_fraction_estimate, when the caller computed them), so a
    caller need not enumerate every field this function knows about. The four keys every
    caller has always been able to rely on (``contigs``, ``total_bp``, ``n50``, ``largest``)
    are still guaranteed present as ints, defaulting to 0 when ``stats`` omits them.
    """
    entry = rb.extraction_block(registry, accession, genome_id)
    if entry is None:
        entry = rb.ExtractionBlock.from_dict({"files": [], "mapped_reads": None})
    entry.discard("inferred")
    entry.assembly = rb.AssemblyBlock(
        date=_now(),
        dir=_project_relative(assembly_dir, project_root(registry)),
        tool="megahit",
        version=tool_version,
        params=dict(params),
        extra={key: value for key, value in stats.items() if key not in _REQUIRED_ASSEMBLY_STATS},
        **{key: int(stats.get(key, 0)) for key in _REQUIRED_ASSEMBLY_STATS},
    )
    rb.set_extraction_block(registry, accession, genome_id, entry)


def clear_assembly(registry: Registry, accession: str, genome_id: str) -> None:
    """Remove one extraction's recorded assembly block, leaving the extraction itself alone.

    Called before an ``--assemble`` redo (forced, or of a stale assembly): if megahit then
    fails, the previous folder is kept on disk but unrecorded, so the registry never claims
    an assembly its inputs no longer match. A no-op when there is no extraction record (or no
    assembly block) for this accession/genome.
    """
    entry = rb.extraction_block(registry, accession, genome_id)
    if entry is not None:
        entry.assembly = None
        rb.set_extraction_block(registry, accession, genome_id, entry)


# --------------------------------------------------------------------- queries


def _extraction_stage(registry: Registry, acc: str, stage: str, genome_id: Optional[str]) -> bool:
    """Handle the "extracted"/"assembled" stages of ``_in_stage`` (kept separate to bound complexity)."""
    if genome_id is None:
        chosen = [e for e in (registry.datasets[acc].get("extractions") or {}).values() if isinstance(e, dict)]
    else:
        chosen = [e for e in (rb.raw(registry, acc, "extractions", genome_id),) if isinstance(e, dict)]
    if stage == "extracted":
        return any((e.get("mapped_reads") or 0) > 0 for e in chosen)
    if stage == "assembled":
        return any(isinstance(e.get("assembly"), dict) and (e["assembly"].get("contigs") or 0) > 0 for e in chosen)
    raise DataAccessError(f"Unknown stage '{stage}'. Choose one of: {', '.join(STAGES)}")


def _in_stage(registry: Registry, acc: str, stage: str, genome_id: Optional[str]) -> bool:
    """Whether ``acc`` is in ``stage``; reads the one field it needs, without building a block."""
    if stage == "screened":
        genomes = rb.raw(registry, acc, "screening", "genomes") or {}
        return bool(genomes) if genome_id is None else genome_id in genomes
    if stage == "selected":
        return bool(rb.raw(registry, acc, "selection", "selected"))
    if stage == "excluded":
        return bool(rb.raw(registry, acc, "exclusion", "excluded"))
    if stage == "downloaded":
        return rb.raw(registry, acc, "download", "state") == "downloaded"
    if stage == "analysed":
        return bool(registry.datasets[acc].get("analyses"))
    return _extraction_stage(registry, acc, stage, genome_id)


def query(registry: Registry, stage: str, genome_id: Optional[str] = None) -> List[str]:
    """Accessions in ``stage`` (insertion order), optionally for one target genome."""
    if stage not in STAGES:
        raise DataAccessError(f"Unknown stage '{stage}'. Choose one of: {', '.join(STAGES)}")
    return [acc for acc in registry.datasets if _in_stage(registry, acc, stage, genome_id)]


def stage_members(registry: Registry) -> Dict[str, List[str]]:
    """Accessions in each of ``STAGES`` (insertion order), as ``query`` returns them, in one pass."""
    members: Dict[str, List[str]] = {stage: [] for stage in STAGES}
    for acc in registry.datasets:
        for stage in STAGES:
            if _in_stage(registry, acc, stage, None):
                members[stage].append(acc)
    return members


def stage_counts(registry: Registry) -> Dict[str, Any]:
    """Return per-stage and per-genome accession counts for the registry.

    The result has a ``"stages"`` mapping of each entry in ``STAGES`` to the number of
    accessions currently in it, and a ``"genomes"`` mapping of each known genome id to its
    extracted and assembled accession counts plus the list of accessions extracted against it
    with zero mapped reads.

    ``status`` uses ``processing.status_report.stage_members`` and ``genome_counts`` instead (one
    pass); this function is kept as the reference those two are tested against.
    """
    stages = {stage: len(query(registry, stage)) for stage in STAGES}
    genomes: Dict[str, Dict[str, Any]] = {}
    for genome_id in sorted(known_genome_ids(registry)):
        zero = [
            acc
            for acc in registry.datasets
            if isinstance(extraction := rb.raw(registry, acc, "extractions", genome_id), dict)
            and (extraction.get("mapped_reads") or 0) == 0
        ]
        genomes[genome_id] = {
            "extracted": len(query(registry, "extracted", genome_id)),
            "assembled": len(query(registry, "assembled", genome_id)),
            "zero_mapped": zero,
        }
    return {"stages": stages, "genomes": genomes}


def known_genome_ids(registry: Registry) -> Set[str]:
    """Return every genome id the registry knows about: recorded genomes plus any seen only in a dataset entry."""
    ids: Set[str] = set(registry.genomes)
    for acc in registry.datasets:
        ids.update(rb.raw(registry, acc, "screening", "genomes") or {})
        ids.update(registry.datasets[acc].get("extractions") or {})
    return ids


# -------------------------------------------------------------------- scanners


def split_extract_filename(name: str, genome_ids: Sequence[str]) -> Optional[Tuple[str, str]]:
    """Split ``<genome><mate>.fastq(.gz)`` into (genome_id, mate suffix); known genome ids win, longest first."""
    stem = name
    for suffix in (".fastq.gz", ".fastq", ".fq.gz", ".fq"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    else:
        return None
    for genome_id in sorted(genome_ids, key=len, reverse=True):
        if stem == genome_id:
            return genome_id, ""
        for mate in _MATE_SUFFIXES:
            if stem == genome_id + mate:
                return genome_id, mate
    for mate in _MATE_SUFFIXES:
        if stem.endswith(mate):
            return stem[: -len(mate)], mate
    return stem, ""


def scan_downloads(fastq_folder: Path) -> Dict[str, Tuple[int, int]]:
    """Accession -> (number of FASTQ files, total bytes) for every per-accession folder with reads.

    Transient folders (``<acc>_temp``, ``.sra-cache``) are skipped even when they hold a
    partial FASTQ file, so an in-progress or failed download is never mistaken for a
    downloaded accession.
    """
    found: Dict[str, Tuple[int, int]] = {}
    for folder in visible_files(fastq_folder, dirs=True):
        if is_transient_folder(folder.name):
            continue
        if accession_has_fastq(folder):
            # fastq_files, not a raw glob: a zero-byte file or the hidden temporary leftover of
            # an interrupted compression is not a downloaded read file.
            files = fastq_files(folder)
            found[folder.name] = (len(files), sum(p.stat().st_size for p in files))
    return found


def scan_metadata(metadata_folder: Path) -> Set[str]:
    """Return the accessions with an ``<accession>_metadata.xml`` file directly under ``metadata_folder``."""
    return {p.name[: -len("_metadata.xml")] for p in visible_files(metadata_folder, "*_metadata.xml")}


def _genome_ids_on_disk(paths: ProjectPaths, registry: Optional[Registry]) -> Set[str]:
    ids: Set[str] = set(known_genome_ids(registry)) if registry else set()
    for pattern in GENOME_FASTA_GLOBS:
        for p in visible_files(paths.genomes, pattern):
            ids.add(p.name[: -len(pattern[1:])])
    ids.update(p.stem for p in visible_files(paths.matches, "*.csv"))
    for acc_dir in visible_files(paths.targeted, dirs=True):
        for asm in visible_files(acc_dir, f"*{_ASSEMBLY_SUFFIX}", dirs=True):
            ids.add(asm.name[: -len(_ASSEMBLY_SUFFIX)])
    return ids


def scan_extractions(targeted_folder: Path, genome_ids: Sequence[str]) -> Dict[str, Dict[str, List[Path]]]:
    """Return the extracted read files found on disk, as accession -> genome id -> list of file paths.

    Walks each per-accession folder under ``targeted_folder`` and classifies its files with
    ``split_extract_filename`` against ``genome_ids``; files that do not match a known genome id
    or mate suffix pattern are skipped. Returns an empty mapping if ``targeted_folder`` does not
    exist as a directory.
    """
    found: Dict[str, Dict[str, List[Path]]] = {}
    if not targeted_folder.is_dir():
        return found
    for acc_dir in visible_files(targeted_folder, dirs=True):
        for path in visible_files(acc_dir):
            split = split_extract_filename(path.name, genome_ids)
            if split:
                found.setdefault(acc_dir.name, {}).setdefault(split[0], []).append(path)
    return found


def scan_assemblies(targeted_folder: Path, genome_ids: Sequence[str]) -> Dict[str, Dict[str, Path]]:
    """Return the assembly directories found on disk, as accession -> genome id -> assembly directory path.

    Walks each per-accession folder under ``targeted_folder`` for subfolders named
    ``<genome_id>_assembly``; ``genome_ids`` is accepted for symmetry with ``scan_extractions``
    but not otherwise used, since the assembly suffix alone identifies the genome id. Returns an
    empty mapping if ``targeted_folder`` does not exist as a directory.
    """
    found: Dict[str, Dict[str, Path]] = {}
    if not targeted_folder.is_dir():
        return found
    for acc_dir in visible_files(targeted_folder, dirs=True):
        for asm in visible_files(acc_dir, f"*{_ASSEMBLY_SUFFIX}", dirs=True):
            found.setdefault(acc_dir.name, {})[asm.name[: -len(_ASSEMBLY_SUFFIX)]] = asm
    return found


def empty_assembly_dirs(targeted_folder: Path, genome_ids: Sequence[str]) -> List[Tuple[str, str]]:
    """(accession, genome_id) pairs whose assembly directory holds zero contigs, sorted."""
    pairs: List[Tuple[str, str]] = []
    for acc, per_genome_asm in scan_assemblies(targeted_folder, genome_ids).items():
        for genome_id, asm_dir in per_genome_asm.items():
            if summarise_contigs(asm_dir / _CONTIGS_NAME)["contigs"] == 0:
                pairs.append((acc, genome_id))
    return sorted(pairs)


def _read_accession_list(path: Optional[Union[str, Path]]) -> List[str]:
    if not path or not Path(path).exists():
        return []
    return [ln.strip() for ln in Path(path).read_text().splitlines() if ln.strip() and not ln.startswith("#")]


def _screening_from_matches(registry: Registry, matches_folder: Path) -> None:
    import csv

    for csv_path in visible_files(matches_folder, "*.csv"):
        with open(csv_path, newline="") as handle:
            for row in csv.DictReader(handle):
                acc = (row.get("acc") or row.get("SRA accession") or "").strip()
                if not acc:
                    continue
                try:
                    containment = float(row.get("containment", ""))
                except ValueError:
                    continue
                cani_raw = row.get("cANI") or ""
                try:
                    cani = float(cani_raw) if cani_raw.strip() else None
                except ValueError:
                    cani = None
                record_screening(registry, acc, csv_path.stem, containment, cani, "matches", 0.0, csv_path)
                rb.mark_inferred(registry, acc, "screening")


def _bootstrap_selection(
    registry: Registry, accessions_file: Optional[Union[str, Path]], parsed_containment: Optional[Union[str, Path]]
) -> None:
    wanted = _read_accession_list(accessions_file)
    if parsed_containment and Path(parsed_containment).exists():
        with open(parsed_containment) as handle:
            next(handle, None)
            wanted += [ln.split("\t", 1)[0].strip() for ln in handle if ln.strip()]
    if wanted:
        record_selection(
            registry, list(dict.fromkeys(wanted)), {"source": "bootstrap"}, accessions_file or parsed_containment or ""
        )
        for acc in wanted:
            rb.mark_inferred(registry, acc, "selection")


def _bootstrap_downloads_and_metadata(registry: Registry, paths: ProjectPaths) -> None:
    for acc in scan_downloads(paths.fastq):
        record_download(registry, acc, "downloaded", paths.fastq)
        rb.mark_inferred(registry, acc, "download", attempts=0)
    root = project_root(registry)
    for acc in sorted(scan_metadata(paths.metadata)):
        record_metadata(registry, acc, paths.metadata / f"{acc}_metadata.xml", {}, root=root)
        rb.mark_inferred(registry, acc, "metadata")


def _infer_extraction(registry: Registry, acc: str, genome_id: str, files: Sequence[Path]) -> None:
    """Record one (accession, genome) extraction found on disk, marked as inferred.

    Reads are counted in every file of the pair (both mates, singles and unpaired), so the
    inferred count approximates the number of BAM records a real extraction records. When the
    extraction's ``<genome>_coverage.tsv`` is on disk, breadth and mean depth are read from it;
    otherwise (or when it cannot be read) they are None.
    """
    reads = sum(count_fastq_reads(f) for f in files)
    coverage = _infer_coverage(coverage_table_path(files[0].parent, genome_id)) if files else None
    record_extraction(registry, acc, genome_id, files, reads, False, {}, coverage=coverage)
    rb.mark_inferred(registry, acc, "extractions", genome_id=genome_id)


def _infer_coverage(tsv: Path) -> Optional[Dict[str, Any]]:
    """Breadth, mean depth and path of an existing ``samtools coverage`` table, else None."""
    if not tsv.is_file():
        return None
    try:
        summary = summarise_coverage_table(tsv)
    except (OSError, ValueError, KeyError) as e:
        logger.warning("Cannot read coverage table %s (%s); recording no coverage for it", tsv, e)
        return None
    return {"breadth": summary["breadth"], "mean_depth": summary["mean_depth"], "coverage_tsv": tsv}


def _infer_assembly(registry: Registry, acc: str, genome_id: str, asm_dir: Path) -> None:
    """Record the assembly in ``asm_dir``, marked as inferred, when it holds contigs."""
    stats = summarise_contigs(asm_dir / _CONTIGS_NAME)
    if stats["contigs"] > 0:
        record_assembly(registry, acc, genome_id, asm_dir, stats, "", {})
        rb.mark_inferred(registry, acc, "extractions", genome_id=genome_id)


def _bootstrap_extractions(registry: Registry, paths: ProjectPaths, genome_ids: List[str]) -> None:
    for acc, per_genome_files in scan_extractions(paths.targeted, genome_ids).items():
        for genome_id, files in per_genome_files.items():
            _infer_extraction(registry, acc, genome_id, files)
    for acc, per_genome_asm in scan_assemblies(paths.targeted, genome_ids).items():
        for genome_id, asm_dir in per_genome_asm.items():
            _infer_assembly(registry, acc, genome_id, asm_dir)


def bootstrap_from_disk(
    paths: ProjectPaths,
    accessions_file: Optional[Union[str, Path]] = None,
    parsed_containment: Optional[Union[str, Path]] = None,
    target_path: Optional[Path] = None,
) -> Registry:
    """Rebuild a registry from what is on disk; every reconstructed block carries ``"inferred": true``.

    ``target_path`` is the file the caller intends to save this registry to (``status --init``
    passes ``registry_file``); it is bound to the registry before any writer below runs, so every
    path recorded during bootstrap is already relative to the right project root, even when
    ``--registry`` points somewhere other than the working directory. Left as ``None`` (a caller
    with no file in mind, e.g. the in-memory report built when no registry exists yet), the
    writers fall back to the working directory as the project root, same as an unbound
    ``Registry()``.
    """
    registry = Registry(path=target_path)
    if paths.matches.is_dir():
        _screening_from_matches(registry, paths.matches)
    _bootstrap_selection(registry, accessions_file, parsed_containment)
    _bootstrap_downloads_and_metadata(registry, paths)
    genome_ids = sorted(_genome_ids_on_disk(paths, registry))
    for genome_id in genome_ids:
        registry.genomes.setdefault(genome_id, {})
    _bootstrap_extractions(registry, paths, genome_ids)
    return registry


def reconcile(registry: Registry, paths: ProjectPaths) -> ReconcileReport:
    """Compare the registry with the disk; see ``prepare.data.registry_reconcile.reconcile``.

    Kept here so existing callers keep importing it from this module. The import is deferred
    because ``registry_reconcile`` imports this module.
    """
    from prepare.data.registry_reconcile import reconcile as _reconcile

    return _reconcile(registry, paths)
