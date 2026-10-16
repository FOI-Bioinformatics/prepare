"""
Adopting a project's own downloaded FASTQ folders into the shared store.

A project that downloaded reads before it adopted the shared store (or ran a
download with the store unconfigured) ends up with real, project-owned
``fastq/<ACC>`` folders instead of links. ``adopt`` folds each of those folders
into the store: staging a copy, compressing any plain FASTQ files, writing the
dataset's sidecar and catalogue entry, and then either replacing the project's
folder with a link to the store's copy (``--move``) or leaving the project's
folder untouched and unlinked, with the store now also holding a copy
(``--copy``).

Staging always copies (never moves) the project's original folder into
``<store>/tmp/<ACC>_adopt``, regardless of ``--move``/``--copy``: the project's
own folder is never touched until the store's copy and its sidecar both exist,
so a crash mid-compression (or anywhere before that point) leaves the project's
data exactly as it was, with nothing to recover beyond removing the stale
staging copy. ``--move`` only ever removes the project's folder, and only
after that point. Peak disk use for one accession is therefore up to three
copies at once: the project's original folder, its staging copy, and (briefly,
during the final move) the store's copy.

Staging one accession runs under that accession's own lock
(``<store>/locks/<ACC>.lock``, see ``prepare.store.locks``), same as a store
download: two adopts of different accessions never touch each other's staging
folder, and two adopts of the same accession serialise on this lock rather
than trampling one another's ``<ACC>_adopt`` copy.

The sidecar is written into the staging folder before it is moved into
``<store>/sra/<ACC>``, so a published folder always carries its sidecar.

Adoption only ever claims the project's own ``fastq/<ACC>`` folders. A dataset
whose files made it into ``<store>/sra/<ACC>`` but whose sidecar was never
written (an adoption by an earlier version, interrupted between those two
steps) is finished in place on the next call, without re-touching the files,
but only when this project also has that accession: a sidecar-less folder for
any other accession is another project's interrupted or still running work,
and is reported rather than blessed with a sidecar and linked here.
"""

import fnmatch
import logging
import re
import shutil
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Union

from prepare.core.constants import SRA_ACCESSION_PATTERN
from prepare.core.exceptions import DataAccessError
from prepare.data.file_io import is_hidden_name, visible_files
from prepare.data.sra import (
    STORE_READY_STATES,
    compress_fastq,
    fastq_digest,
    fastq_files,
    fastq_stem,
    is_transient_folder,
)
from prepare.store.catalog import catalog_write
from prepare.store.layout import StorePaths, lock_path, sidecar_path, sra_dir
from prepare.store.link import link_dataset, room_shortfall
from prepare.store.locks import LockHeld, dataset_lock, lock_holder, lock_is_held
from prepare.store.sidecar import (
    Sidecar,
    build_sidecar,
    detect_layout,
    md5_file,
    ncbi_from_metadata_xml,
    read_sidecar,
    write_sidecar,
)
from prepare.utils.lockfile import LockLost, verify_held

logger = logging.getLogger(__name__)

# Threads handed to compress_fastq while adopting; adoption processes one accession's files at a
# time rather than many downloads in parallel, so a modest, fixed thread count is enough.
_ADOPT_COMPRESS_THREADS = 4

# Suffix on a staging folder under paths.tmp, e.g. "SRR1_adopt".
_STAGING_SUFFIX = "_adopt"

# Name patterns shutil.copytree leaves behind when staging a project folder: the AppleDouble
# sidecar files (``._<name>``) macOS writes next to every file on a volume without native
# extended attributes, and the per-folder Finder metadata file ``.DS_Store``. Neither is
# project data. Kept as one pattern list so the copy (``ADOPT_COPY_IGNORE``) and the
# free-space estimate (``_folder_bytes``, via ``_copy_ignores``) can never drift apart: a name
# the copy skips is a name the byte count must skip too, and nothing else -- an ordinary
# hidden directory or file (``.hidden/``, a pipeline's ``.snakemake/`` work folder) matches
# neither pattern, so the copy does not skip it and the byte count must not either. The two
# temporary-file patterns are what an interrupted write leaves: ``.<name>.<host>.<pid>.<token>.tmp``
# (``unique_temp_path``) and the ``<name>.tmp.<pid>`` of versions before 0.6.0.
_COPY_IGNORE_PATTERNS = ("._*", ".DS_Store", ".*.tmp", "*.tmp.[0-9]*")
ADOPT_COPY_IGNORE = shutil.ignore_patterns(*_COPY_IGNORE_PATTERNS)


def _copy_ignores(name: str) -> bool:
    """True when ``name`` matches one of ``_COPY_IGNORE_PATTERNS``, the same test
    ``ADOPT_COPY_IGNORE`` applies during the staging copy."""
    return any(fnmatch.fnmatch(name, pattern) for pattern in _COPY_IGNORE_PATTERNS)


@dataclass
class AdoptReport:
    """Outcome of one ``adopt`` call, one accession name per relevant list."""

    adopted: List[str] = field(default_factory=list)
    # --copy: the store now holds a copy too, but the project's folder is untouched and unlinked.
    copied: List[str] = field(default_factory=list)
    deduplicated: List[str] = field(default_factory=list)
    conflicts: List[str] = field(default_factory=list)
    # Entries in project_fastq that were already store links, so nothing to adopt.
    skipped: List[str] = field(default_factory=list)
    # Real folders with no FASTQ files at all (an empty or interrupted download); not adopted,
    # the project's folder is left exactly as it was.
    empty: List[str] = field(default_factory=list)
    # Staged (or already-in-store) copy came back in a non-ready sidecar state (state not in
    # STORE_READY_STATES, e.g. "failed" or "partial"); the project's copy is kept untouched.
    failed: List[str] = field(default_factory=list)
    # Populated instead of acting, when dry_run is True.
    planned: List[str] = field(default_factory=list)
    # A stale <ACC>_adopt staging folder from an earlier, interrupted run of this same
    # accession, removed (under this accession's lock) right before staging; the project's own
    # folder was untouched by that interruption, so the accession is simply re-adopted here.
    resumed: List[str] = field(default_factory=list)
    # Sidecar-less store folders for accessions this project does not have: another project's
    # interrupted (or still running) work, reported and never touched.
    foreign: List[str] = field(default_factory=list)
    # Accessions whose store lock is held right now, so another run is publishing them.
    in_progress: List[str] = field(default_factory=list)
    # Accessions left alone because the store's filesystem has no room to stage them.
    refused: List[str] = field(default_factory=list)


def _notify(on_progress: Optional[Callable[[str, str], None]], accession: str, event: str) -> None:
    """Call ``on_progress`` in isolation so a failing callback never breaks adoption itself."""
    if on_progress is None:
        return
    try:
        on_progress(accession, event)
    except Exception as e:  # pragma: no cover - defensive, mirrors prepare.data.sra._notify_result
        logger.warning("Progress callback failed for %s: %s", accession, e)


def _content_matches(path: Path, store_path: Path, record: Dict[str, Any]) -> bool:
    """True when ``path`` and the store copy hold the same reads, whatever their raw bytes.

    Compares the md5 of the decompressed content of both sides, and ``path``'s read count with
    the sidecar's; ``fastq_digest`` yields both from one pass over each file, so the callers
    use this only when the cheap byte comparison cannot settle the question.
    """
    if not store_path.is_file():
        return False
    project_digest = fastq_digest(path, content_md5=True)
    if project_digest.content_md5 != fastq_digest(store_path, content_md5=True).content_md5:
        return False
    return project_digest.records == record.get("reads")


def _files_match(project_dir: Path, store_dir: Path, sidecar: Sidecar) -> bool:
    """True when ``project_dir`` holds the same reads the sidecar's store copy records.

    Files are paired by ``fastq_stem`` (e.g. ``SRR1_1``, not ``SRR1_1.fastq.gz``), so a plain
    project copy matches a gzipped store copy of the same reads and vice versa. Compares the
    sidecar's recorded bytes and md5 first, which settles an identical copy without
    decompressing anything.

    Raw bytes are not decisive for two gzip files: the gzip header carries the compression
    time and the original file name, and a different compressor or level changes the body,
    so the same reads compressed twice differ byte for byte. A gzip pair whose recorded md5
    does not match therefore falls back to the decompressed comparison, as a pair with
    different compression on the two sides already does. Only that fallback can tell a
    duplicate from a genuine conflict here.
    """
    project_files = {fastq_stem(p): p for p in fastq_files(project_dir)}
    sidecar_files = {fastq_stem(Path(str(entry.get("name") or ""))): entry for entry in sidecar.files}
    if set(project_files) != set(sidecar_files):
        return False
    for stem, path in project_files.items():
        record = sidecar_files[stem]
        record_name = str(record.get("name") or "")
        compressed = str(path).endswith(".gz")
        if compressed == record_name.endswith(".gz"):
            if path.stat().st_size == record.get("bytes") and md5_file(path) == record.get("md5"):
                continue
            if not compressed:
                # Two plain files: the bytes are the content, so a mismatch is a real one.
                return False
        if not _content_matches(path, store_dir / record_name, record):
            return False
    return True


def _detect_compression(acc_dir: Path) -> str:
    """``"gzip"`` when every FASTQ file in ``acc_dir`` is gzip-compressed, else ``"none"``."""
    files = fastq_files(acc_dir)
    if files and all(str(p).endswith(".gz") for p in files):
        return "gzip"
    return "none"


def _ncbi_from_folders(accession: str, metadata_folders: Sequence[Union[str, Path]]) -> Dict[str, Any]:
    """The ``ncbi`` dict for ``accession`` from the first metadata folder that has its XML."""
    for folder in metadata_folders:
        candidate = Path(folder) / f"{accession}_metadata.xml"
        if candidate.exists():
            return ncbi_from_metadata_xml(candidate) or {}
    return {}


def _newest_file_time(acc_dir: Path) -> Optional[str]:
    """When the newest FASTQ file in ``acc_dir`` was last written, as an ISO timestamp.

    An adopted dataset was downloaded at some earlier point that nothing recorded, so its
    files' own age is the best available answer; ``store_gc --older-than`` then counts from
    when the reads were obtained rather than from when they were adopted.
    """
    times = []
    for file_path in fastq_files(acc_dir):
        try:
            times.append(file_path.stat().st_mtime)
        except OSError:
            continue
    if not times:
        return None
    return datetime.fromtimestamp(max(times), tz=timezone.utc).isoformat()


def _adopted_sidecar(accession: str, folder: Path, metadata_folders: Sequence[Union[str, Path]]) -> Sidecar:
    """The sidecar describing an adopted dataset's files in ``folder``.

    Records ``tool="adopted"`` and the files' own age rather than this moment: nothing here
    downloaded them, and no tool version is known.
    """
    return build_sidecar(
        accession,
        folder,
        _ncbi_from_folders(accession, metadata_folders),
        "",
        _detect_compression(folder),
        tool="adopted",
        downloaded=_newest_file_time(folder),
    )


def _write_staged_sidecar(accession: str, staged: Path, metadata_folders: Sequence[Union[str, Path]]) -> Sidecar:
    """Write ``<ACC>.json`` into the staging folder ``staged``, before it is published.

    The sidecar then travels with the files in the one rename that publishes the dataset, so
    a run killed at any point leaves either no store folder or one that carries its sidecar.
    """
    sidecar = _adopted_sidecar(accession, staged, metadata_folders)
    write_sidecar(staged / f"{accession}.json", sidecar)
    return sidecar


def _finish_sidecar(
    accession: str,
    store_dir: Path,
    sc_path: Path,
    paths: StorePaths,
    metadata_folders: Sequence[Union[str, Path]],
) -> Sidecar:
    """Write the sidecar and catalogue entry for a dataset whose files already sit in the store.

    Finishes an adoption that an earlier version interrupted after publishing the files but
    before writing their sidecar (a fresh adoption writes the sidecar into the staging folder
    instead, see ``_write_staged_sidecar``).
    """
    sidecar = _adopted_sidecar(accession, store_dir, metadata_folders)
    write_sidecar(sc_path, sidecar)
    with catalog_write(paths) as catalog:
        catalog.upsert_dataset(sidecar)
    return sidecar


# Error a rebuilt sidecar carries, so every reader treats the dataset as not yet usable and
# says how to check it; ``{acc}`` is the accession.
_REBUILT_ERROR = (
    "sidecar rebuilt by store_reindex from the files on disk; the files were not verified "
    "(run 'store_verify --rescan --fix-state {acc}')"
)


def rebuild_missing_sidecar(accession: str, store_dir: Path, paths: StorePaths) -> Sidecar:
    """Write a sidecar for a store folder that has none, from its files' names and sizes.

    For ``store_reindex``, which then catalogues it; nothing is catalogued here. The caller
    holds ``accession``'s lock. No file is read, so nothing is vouched for: the sidecar's state
    is ``failed`` with an error saying so, its completeness ``unverified``, and every file's
    md5 and read count are None, until ``store_verify --rescan --fix-state`` reads the files
    through. A sidecar that appeared since the caller looked (another run finished the
    dataset) is returned as it is, unchanged.
    """
    sc_path = sidecar_path(paths, accession)
    if sc_path.is_file():
        existing = read_sidecar(sc_path)
        if existing is None:
            # Never overwrite a sidecar that is there but unreadable; that needs a person.
            raise DataAccessError(f"{sc_path}: sidecar exists but cannot be read; not rebuilding it")
        return existing
    files = fastq_files(store_dir)
    sidecar = Sidecar(
        accession=accession,
        state="failed",
        layout=detect_layout(files),
        downloaded=_newest_file_time(store_dir),
        tool="unknown",
        compression=_detect_compression(store_dir),
        files=[{"name": path.name, "bytes": path.stat().st_size, "md5": None, "reads": None} for path in files],
        completeness={"method": "unverified", "ratio": None, "verdict": "unverified"},
        error=_REBUILT_ERROR.format(acc=accession),
    )
    verify_held(lock_path(paths, accession))
    write_sidecar(sc_path, sidecar)
    return sidecar


def _folder_bytes(folder: Path) -> int:
    """Total bytes of every file under ``folder`` that the staging copy
    (``shutil.copytree(..., ignore=ADOPT_COPY_IGNORE)``, see ``_stage_into_store``) will
    actually copy, so the free-space check in ``_has_room_for`` estimates the same folder the
    copy produces rather than a different one.

    Skips a path when some segment of its path relative to ``folder`` matches
    ``_copy_ignores`` (an AppleDouble ``._*`` sidecar or ``.DS_Store``, at any depth) -- the
    same names the copy itself skips, wherever they appear. A plain-named file inside an
    otherwise hidden directory (``.hidden/keep.txt``, a pipeline's ``.snakemake/`` work
    folder) is *not* skipped: ``.hidden``/``.snakemake`` match neither ignore pattern, so
    ``shutil.copytree`` copies that file too, and this must count it. Anything that cannot be
    stat'ed is skipped."""
    total = 0
    for sub in folder.rglob("*"):
        if any(_copy_ignores(p) for p in sub.relative_to(folder).parts):
            continue
        try:
            if sub.is_file():
                total += sub.stat().st_size
        except OSError:
            continue
    return total


def _has_room_for(accession: str, entry: Path, paths: StorePaths) -> bool:
    """True when the store's filesystem has room to stage ``entry``, with a copy to spare.

    Adoption peaks at three copies of one accession (the project's, the staging copy and,
    briefly, the store's), so it needs at least twice the folder's size free before it starts.
    A filesystem that cannot be measured is assumed to have room: refusing on an unreadable
    ``disk_usage`` would be worse than trying.
    """
    needed = _folder_bytes(entry) * 2
    free = room_shortfall(paths.root, needed)
    if free is None:
        return True
    logger.warning(
        "%s: not adopting it, the store's filesystem has %d bytes free and staging needs about %d",
        accession,
        free,
        needed,
    )
    return False


def _stage_into_store(
    entry: Path,
    accession: str,
    paths: StorePaths,
    compress: bool,
    metadata_folders: Sequence[Union[str, Path]],
    report: AdoptReport,
) -> None:
    """Copy ``entry``'s files into ``<store>/sra/<accession>``, compressing plain files on the way.

    Always stages through a *copy* into ``paths.tmp`` first, regardless of ``--move``/``--copy``:
    the project's original folder is only ever removed (for ``--move``, by the caller, after this
    returns) once the store's copy and its sidecar both exist. A crash anywhere in this function
    leaves ``entry`` untouched; a stale ``<ACC>_adopt`` staging copy left by such a crash is
    removed here, right before staging starts again, and reported in ``report.resumed``.

    The caller holds this accession's lock (``<store>/locks/<ACC>.lock``, see ``_adopt_one``).
    """
    store_dir = sra_dir(paths, accession)
    staged = paths.tmp / f"{accession}{_STAGING_SUFFIX}"
    paths.tmp.mkdir(parents=True, exist_ok=True)

    if staged.exists():
        shutil.rmtree(staged)
        report.resumed.append(accession)

    shutil.copytree(entry, staged, ignore=ADOPT_COPY_IGNORE)

    if compress:
        for file_path in fastq_files(staged):
            if not str(file_path).endswith(".gz"):
                compress_fastq(file_path, _ADOPT_COMPRESS_THREADS)

    # The sidecar goes into the staged folder first, so the folder is published with it.
    sidecar = _write_staged_sidecar(accession, staged, metadata_folders)

    store_dir.parent.mkdir(parents=True, exist_ok=True)
    # A holder stalled past the stale window may have lost the lock (and the staging path) to
    # another run: publish nothing then, and leave the staged copy to the new holder.
    verify_held(lock_path(paths, accession))
    shutil.move(str(staged), str(store_dir))

    with catalog_write(paths) as catalog:
        catalog.upsert_dataset(sidecar)


def _scan_project_dir(project_dir: Path, report: AdoptReport) -> Dict[str, Path]:
    """Real (non-symlink, non-transient) accession directories directly in ``project_dir``.

    A symlink is recorded in ``report.skipped`` (already a store link, nothing to adopt) rather
    than returned. A folder not named like an SRA accession (notes kept beside the datasets, for
    example) is not a dataset: it is logged as skipped and left in place, as ``store_reindex`` does.
    """
    real_dirs: Dict[str, Path] = {}
    if not project_dir.is_dir():
        return real_dirs
    for candidate in sorted(project_dir.iterdir()):
        name = candidate.name
        if is_hidden_name(name):
            continue
        if is_transient_folder(name):
            continue
        if candidate.is_symlink():
            report.skipped.append(name)
            continue
        if not candidate.is_dir():
            continue
        if not re.fullmatch(SRA_ACCESSION_PATTERN, name):
            logger.warning("%s in %s is not named like an SRA accession; skipped", name, project_dir)
            continue
        real_dirs[name] = candidate
    return real_dirs


def _scan_foreign_incomplete(paths: StorePaths, real_dirs: Dict[str, Path]) -> List[str]:
    """Sidecar-less store folders for accessions this project is not adopting.

    Such a folder is another project's work: an interrupted run of theirs, or one still in
    flight. Adoption reports it so the store's owner can see it, and never blesses it with a
    sidecar or links it here. Only the project's own folders (``real_dirs``) are ever adopted,
    so a dataset another project is mid-download on cannot be claimed by this run. A folder
    not named like an SRA accession is not a dataset: it is logged as skipped, as
    ``store_reindex`` does, before any path is built from its name.
    """
    foreign: List[str] = []
    for store_dir in visible_files(paths.sra, dirs=True):
        accession = store_dir.name
        if not re.fullmatch(SRA_ACCESSION_PATTERN, accession):
            logger.warning("%s in %s is not named like an SRA accession; skipped", accession, paths.sra)
            continue
        if accession in real_dirs or sidecar_path(paths, accession).is_file():
            continue
        foreign.append(accession)
        logger.warning("foreign: %s, no sidecar; left untouched", accession)
    return foreign


def _dedup_or_conflict(
    accession: str,
    entry: Optional[Path],
    store_dir: Path,
    sc_path: Path,
    paths: StorePaths,
    project_dir: Path,
    dry_run: bool,
    move: bool,
    on_progress: Optional[Callable[[str, str], None]],
    report: AdoptReport,
) -> None:
    """Handle an accession whose store folder already carries a sidecar: dedup or conflict.

    ``entry`` is None only for an accession already fully migrated (store has it, sidecar
    exists, nothing local left to reconcile); nothing to do in that case.

    A dedup is reported the same way whether or not the project's folder is actually removed:
    ``--copy`` (``move=False``) leaves it in place, unlinked, exactly as ``--copy`` does for a
    freshly adopted accession (see ``_apply_move_or_copy``).

    A store sidecar whose state is outside ``STORE_READY_STATES`` (left behind by an earlier
    run whose staged copy could not be verified, see ``_adopt_one``'s own check) is never
    treated as a match, even when its files are byte-identical to the project's: that copy is
    known bad, so an identical-content project folder is the only good copy left, and a dedup
    would remove it and link the project to the bad one instead. Such an accession is reported
    in ``failed`` and both copies are left exactly as they are.
    """
    if entry is None:
        return
    existing = read_sidecar(sc_path)
    if existing is not None and existing.state not in STORE_READY_STATES:
        report.failed.append(accession)
        logger.error(
            "%s: store copy's sidecar state is %s, not a verified copy; project copy kept "
            "(run 'store_verify --spots --fix-state %s' or remove the store copy and re-adopt)",
            accession,
            existing.state,
            accession,
        )
        _notify(on_progress, accession, "failed")
        return
    if existing is not None and _files_match(entry, store_dir, existing):
        if dry_run:
            report.planned.append(accession)
            return
        if not move:
            report.deduplicated.append(accession)
            logger.info("%s: identical copy already in the store; project copy kept (--copy)", accession)
            _notify(on_progress, accession, "deduplicated")
            return
        if _replace_with_link(accession, entry, paths, project_dir, on_progress, report):
            report.deduplicated.append(accession)
            _notify(on_progress, accession, "deduplicated")
        return
    report.conflicts.append(accession)
    logger.warning(
        "%s: project and store copies differ; leaving both in place (store copy left untouched)",
        accession,
    )
    _notify(on_progress, accession, "conflict")


def _apply_move_or_copy(
    accession: str,
    entry: Path,
    paths: StorePaths,
    project_dir: Path,
    move: bool,
    on_progress: Optional[Callable[[str, str], None]],
    report: AdoptReport,
) -> None:
    """Finish one adopted accession: ``--move`` links it, ``--copy`` leaves the project alone.

    ``--copy`` leaves the project's folder exactly as it was (a second, independent copy) and
    does not link it, since a project that asked to keep its own copy should not have it
    silently swapped for a link on the same run.
    """
    if not move:
        report.copied.append(accession)
        _notify(on_progress, accession, "copied")
        return

    if _replace_with_link(accession, entry, paths, project_dir, on_progress, report):
        report.adopted.append(accession)
        _notify(on_progress, accession, "adopted")


def _replace_with_link(
    accession: str,
    entry: Path,
    paths: StorePaths,
    project_dir: Path,
    on_progress: Optional[Callable[[str, str], None]],
    report: AdoptReport,
) -> bool:
    """Replace the project's folder ``entry`` with a link to the store's copy; return whether it was.

    The removal and the link run under one hold of the dataset lock, taken without waiting, so
    no other run (``store_gc`` in particular) can act on the store's copy in between. When the
    lock is held elsewhere the accession is reported in ``failed`` and the project's folder is
    left in place.
    """
    try:
        with dataset_lock(paths, accession, blocking=False):
            if entry.exists():
                shutil.rmtree(entry)
            link_dataset(project_dir, accession, paths, locked=True)
            return True
    except LockHeld:
        logger.error(
            "%s: %s is being written by %s; retry later (project copy kept)",
            accession,
            accession,
            lock_holder(paths, accession),
        )
        report.failed.append(accession)
        _notify(on_progress, accession, "failed")
        return False


def _adopt_one(
    accession: str,
    entry: Path,
    paths: StorePaths,
    project_dir: Path,
    move: bool,
    compress: bool,
    metadata_folders: Sequence[Union[str, Path]],
    on_progress: Optional[Callable[[str, str], None]],
    report: AdoptReport,
    lock_wait: float,
) -> None:
    """Get one of this project's accessions into the store, then apply ``--move``/``--copy``.

    Everything that writes to the store happens under this accession's own lock (see
    ``prepare.store.locks.dataset_lock``), and the store's state is read again inside it:
    a download that finished while this run waited leaves a sidecar behind, and this run must
    then compare copies rather than move its own on top of the published one.

    A sidecar-less store folder whose lock is held right now belongs to a run publishing it at
    this moment; that accession is reported ``in progress elsewhere`` and left alone.
    """
    store_dir = sra_dir(paths, accession)
    sc_path = sidecar_path(paths, accession)

    if store_dir.is_dir() and not sc_path.is_file() and lock_is_held(paths, accession):
        logger.warning(
            "%s: in progress elsewhere (%s); leaving both copies alone", accession, lock_holder(paths, accession)
        )
        report.in_progress.append(accession)
        return

    if not store_dir.is_dir() and not _has_room_for(accession, entry, paths):
        report.refused.append(accession)
        return

    paths.locks.mkdir(parents=True, exist_ok=True)
    try:
        with dataset_lock(paths, accession, wait_seconds=lock_wait):
            published_elsewhere = sc_path.is_file()
            if not published_elsewhere:
                if store_dir.is_dir():
                    _finish_sidecar(accession, store_dir, sc_path, paths, metadata_folders)
                else:
                    _stage_into_store(entry, accession, paths, compress, metadata_folders, report)
    except LockLost as e:
        logger.error("%s: lock lost: %s; nothing published, project copy kept", accession, e)
        report.failed.append(accession)
        return

    if published_elsewhere:
        # Someone published this accession while we waited for the lock; the project's copy is
        # now a second copy to compare, not something to move on top of theirs.
        _dedup_or_conflict(accession, entry, store_dir, sc_path, paths, project_dir, False, move, on_progress, report)
        return

    published = read_sidecar(sc_path)
    if published is None or published.state not in STORE_READY_STATES:
        state = published.state if published else "missing sidecar"
        logger.error(
            "%s: store copy is %s (%s); project copy kept",
            accession,
            state,
            (published.error if published else ""),
        )
        report.failed.append(accession)
        return

    _apply_move_or_copy(accession, entry, paths, project_dir, move, on_progress, report)


def adopt(
    project_fastq: Union[str, Path],
    paths: StorePaths,
    move: bool,
    dry_run: bool,
    on_progress: Optional[Callable[[str, str], None]] = None,
    compress: bool = True,
    metadata_folders: Sequence[Union[str, Path]] = (),
    lock_wait: float = 0.0,
    should_stop: Optional[Callable[[], bool]] = None,
) -> AdoptReport:
    """Fold every real accession folder in ``project_fastq`` into the shared store.

    For each real (non-symlink, non-transient) directory ``project_fastq/<ACC>``: if the store
    already has ``<ACC>`` with a sidecar in a ready state (``STORE_READY_STATES``), an identical
    set of files is deduplicated (the project copy is dropped and replaced with a link) and a
    differing set is left as a conflict (both copies kept); otherwise the folder is staged into
    the store, compressed, sidecar'd and catalogued, then either linked back (``--move``,
    removing the project's folder) or left alone (``--copy``, the project keeps its own folder,
    unlinked). A symlink already in ``project_fastq`` is counted in ``skipped`` and left
    untouched. A dedup under ``--copy`` is still counted in ``deduplicated``, but the project's
    folder is left in place and unlinked, the same as a fresh ``--copy`` adoption.

    A store sidecar outside a ready state is never matched against, even byte-identical: that
    copy already failed verification once, so an identical project folder is the only good copy
    and must not be removed on top of it. Such an accession lands in ``failed`` on this call too
    (not just on the run that first staged it), and both copies are left exactly as they are.

    A real folder with no FASTQ files at all is never adopted (the store gains nothing from an
    empty accession) and is listed in ``empty``; the project's folder is left exactly as it was.
    A staged (or already-in-store) copy whose sidecar comes back outside
    ``prepare.data.sra.STORE_READY_STATES`` (``"failed"`` or ``"partial"``) is listed in
    ``failed`` instead of being linked or copied over the project's folder, which is left
    untouched; the store keeps its own copy for inspection.

    Only the project's own folders are ever adopted. A sidecar-less ``<store>/sra/<ACC>`` for an
    accession this project also has is finished in place (an earlier run of this project died
    after the files reached the store but before the sidecar was written), under that accession's
    lock; a sidecar-less folder for any other accession belongs to another project, and is listed
    in ``foreign`` without being touched or linked. An accession whose lock is held right now is
    listed in ``in_progress`` and left alone, and one the store's filesystem has no room to stage
    is listed in ``refused``. Any stale ``<store>/tmp/<ACC>_adopt`` copy from a run interrupted
    mid-staging is removed (under that accession's lock, right before re-staging) and the
    accession reported in ``resumed``. Peak disk use for one accession adopted this way is up to
    three copies at once: see the module docstring.

    ``dry_run`` performs no filesystem changes and no database changes (including no sweep of a
    stale staging copy); every accession that would otherwise be dedup'd, adopted or restarted is
    instead listed in ``planned``. Conflicts are still detected and reported in ``conflicts`` even
    during a dry run, since detecting one never writes anything.

    ``should_stop`` is checked before each accession; once it returns True the loop stops there,
    leaving every later accession (in sort order) untouched -- neither staged nor counted in any
    report list. An accession already staged or adopted when the stop is noticed is left exactly
    as it finished: staging itself only ever touches ``<store>/tmp/<ACC>_adopt`` and the project's
    folder is replaced with a link (``--move``) only after the store's copy is verified, so a
    caller that stops the run here never has to distinguish a finished accession from one that
    was merely in flight.
    """
    report = AdoptReport()
    project_dir = Path(project_fastq)

    real_dirs = _scan_project_dir(project_dir, report)
    report.foreign.extend(_scan_foreign_incomplete(paths, real_dirs))

    for accession in sorted(real_dirs):
        if should_stop is not None and should_stop():
            logger.warning("Stopping before %s: interrupted", accession)
            break
        entry = real_dirs[accession]

        if not fastq_files(entry):
            logger.warning("%s: no FASTQ files in %s; not adopted", accession, entry)
            report.empty.append(accession)
            continue

        store_dir = sra_dir(paths, accession)
        sc_path = sidecar_path(paths, accession)

        if store_dir.is_dir() and sc_path.is_file():
            _dedup_or_conflict(
                accession, entry, store_dir, sc_path, paths, project_dir, dry_run, move, on_progress, report
            )
            continue

        # Either a fresh accession (no store folder yet) or one this project left incomplete
        # in the store (folder there, no sidecar).
        if dry_run:
            report.planned.append(accession)
            continue

        _adopt_one(
            accession,
            entry,
            paths,
            project_dir,
            move,
            compress,
            metadata_folders,
            on_progress,
            report,
            lock_wait,
        )

    return report
