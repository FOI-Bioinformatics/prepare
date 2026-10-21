"""Download of one SRA accession with prefetch and fasterq-dump, and classification of its failures."""

import functools
import logging
import re
import shutil
import subprocess
import threading
import zlib
from pathlib import Path
from typing import AbstractSet, Any, Callable, Dict, List, Optional, Tuple, Union

from prepare.core.constants import DATASET_LOCK_STALE_SECONDS, LOCK_HEARTBEAT_SECONDS
from prepare.core import settings
from prepare.core.settings import settings_or
from prepare.core.exceptions import DataAccessError, SecurityError
from prepare.data.accessions import require_accession
from prepare.data.sra import cleanup as cleanup_mod
from prepare.data.sra import fastq as fastq_mod
from prepare.utils.lockfile import LockHeld, LockLost, LockPolicy, LockWaitStopped, held_lock, verify_held
from prepare.utils import tools
from prepare.utils.progress import active_item_level
from prepare.utils.security import SecureSubprocess

logger = logging.getLogger(__name__)

# Set once the missing-prefetch fallback was logged at INFO; later accessions log it at item level.
_PREFETCH_FALLBACK_LOGGED = threading.Event()


# Regexes used by classify_download_error to sort a failure message into a coarse class that
# retry logic can act on: retry network/unknown failures, never retry not-found, and abort the
# whole run on disk-full.
# "could not resolve host" rather than a bare "resolve": prefetch reports a missing accession
# as "failed to resolve accession ... no data ( 404 )", which is a not-found, not a network error.
_NETWORK_ERROR_RE = re.compile(r"timeout|timed out|connection|could not resolve host|network|curl|ssl", re.IGNORECASE)
# "storage exhausted" and "disk-limit exeeded" (sic) are fasterq-dump's own out-of-space messages.
_DISK_FULL_ERROR_RE = re.compile(r"no space left|enospc|disk[ -]full|storage exhausted|disk-limit", re.IGNORECASE)
_NOT_FOUND_ERROR_RE = re.compile(
    r"not[ -]found|invalid accession|cannot be found|\b403\b|\b404\b|does not exist", re.IGNORECASE
)
# A lock message names the holder's pid, start time and lock path, any of which can contain
# "403" or "404"; it is classified before the other patterns are tried.
_LOCK_MESSAGE_RE = re.compile(r"(?:^|: )(?:locked|lock lost): | is locked by pid ")


# The result message for an accession whose FASTQ files were already in place, so nothing ran.
ALREADY_EXISTS = "already exists"

# The retry pass prefixes a result message with "Retry <n>: " (or "Retry <n> error: " for a
# worker that raised); classification reads the message underneath.
_RETRY_PREFIX_RE = re.compile(r"^Retry \d+(?: error)?: ")


def result_body(message: str) -> str:
    """``message`` without the retry pass's ``Retry <n>: `` prefix."""
    return _RETRY_PREFIX_RE.sub("", message or "", count=1)


# Process-wide emergency stop. Each download run has its own stop token (``download_sra``'s
# ``stop``), set when that run is interrupted; download_accession checks both before each
# prefetch or fasterq-dump call, so a worker thread that has not yet started a tool returns
# without starting one, and the retry pass does not run. Setting STOP stops every run in the
# process; no run clears it (tests and callers that set it clear it themselves).
STOP = threading.Event()


# A plain project's per-accession lock is <fastq>/.locks/<ACC>.lock; a waiter re-checks a held
# one this often. Read at call time so tests can shorten it.
PROJECT_LOCK_FOLDER = ".locks"
PROJECT_LOCK_POLL_SECONDS = 1.0


class _DownloadInterrupted(Exception):
    """Raised inside download_accession when the run is stopping before a tool call."""


def stop_requested(stop: Optional[threading.Event] = None) -> bool:
    """Whether a download should stop: the run's token ``stop`` or the process-wide ``STOP`` is set."""
    return (stop is not None and stop.is_set()) or STOP.is_set()


def _run_download_tool(executable: str, args: List[str], stop: Optional[threading.Event] = None) -> None:
    """Run prefetch or fasterq-dump through ``run_secure`` unless the run is stopping.

    The child is recorded under ``stop`` so ``terminate_children(stop=stop)`` stops this run's
    tools only. A tool that exits non-zero once the run is stopping was stopped by the
    interrupt (or killed as it started), so its failure is reported as an interruption rather
    than an error.
    """
    if stop_requested(stop):
        raise _DownloadInterrupted()
    try:
        SecureSubprocess.run_secure(executable, args, stop=stop)
    except subprocess.CalledProcessError as e:
        if stop_requested(stop):
            raise _DownloadInterrupted() from e
        raise


def classify_download_error(text: str) -> str:
    """Classify a download failure message into a coarse error class.

    Returns one of ``"network"``, ``"disk-full"``, ``"not-found"``, or ``"unknown"``. Pure
    function of the message text; used both to prefix ``download_accession``'s failure
    messages and to decide, in ``_retry_failed_downloads``, which accessions are worth
    retrying.

    A lock message (``locked: ...``, ``lock lost: ...``, a store lock that could not be taken) is
    ``unknown`` whatever pid or path it names. ``not-found`` is tested before ``network``
    because a prefetch not-found message mentions resolving a query, and each class name
    classifies back to its own class so an already prefixed message
    (``"disk-full: not attempted"``) keeps its class on a second pass.
    """
    if not text or _LOCK_MESSAGE_RE.search(text):
        return "unknown"
    if _NOT_FOUND_ERROR_RE.search(text):
        return "not-found"
    if _DISK_FULL_ERROR_RE.search(text):
        return "disk-full"
    if _NETWORK_ERROR_RE.search(text):
        return "network"
    return "unknown"


def _notify_result(
    on_result: Optional[Callable[[str, bool, str], None]], accession: str, success: bool, message: str
) -> None:
    """Call ``on_result`` in isolation so a failing callback never corrupts the download tally.

    A registry write (or any other callback) can raise, e.g. a lock timeout. That must not be
    mistaken for the download itself failing, so this is never allowed to propagate.
    """
    if on_result is None:
        return
    try:
        on_result(accession, success, message)
    # The callback is supplied by the caller and may raise anything (a registry lock timeout,
    # a JSON error); none of it may be mistaken for the download failing.
    except Exception as e:  # noqa: B902 - callback of unknown type, see the docstring
        logger.warning(f"Recording the result for {accession} failed: {e}")


def _check_existing_download(output_path, force):
    """
    Check if the accession is already downloaded.

    Args:
        output_path: Path to output directory
        force: Whether to force redownload

    Returns:
        True if already downloaded, False otherwise
    """
    # A dangling symlink (e.g. left by an interrupted run) reports False from `.exists()`
    # even though the directory entry itself is still there; clear it so a fresh download can
    # create a real directory at this path instead of failing on FileExistsError.
    if output_path.is_symlink() and not output_path.exists():
        try:
            output_path.unlink()
        except OSError as e:
            logger.warning(f"Could not remove dangling symlink {output_path}: {e}")
        return False

    if force and output_path.exists():
        # Force redownload - remove existing directory
        cleanup_mod._remove_stale_entry(output_path)
        if not output_path.exists():
            logger.info(f"Removed existing directory for force redownload: {output_path}")
        return False

    if not force and output_path.exists():
        if fastq_mod.accession_has_fastq(output_path):
            return True

        # Found empty (or not-yet-complete) directory, will redownload.
        try:
            if output_path.is_symlink():
                output_path.unlink()
            else:
                output_path.rmdir()
        except OSError as e:
            logger.warning(f"Could not remove empty directory {output_path}: {e}")

    return False


def _handle_download_output(
    temp_path,
    output_path,
    expected_spots: Optional[int] = None,
    compress: bool = False,
    num_threads: int = 4,
    stop: Optional[threading.Event] = None,
    accept_partial: bool = False,
):
    """
    Move downloaded files from temp path to output path, then verify completeness.

    Args:
        temp_path: Path to temporary folder
        output_path: Path to output directory
        expected_spots: NCBI's recorded total_spots for this accession, if known
        compress: If True, gzip each downloaded FASTQ file (via ``compress_fastq``)
            after the completeness verdict has been computed on the plain files
        num_threads: Thread count passed to ``compress_fastq`` (for pigz's ``-p``)
        stop: The run's stop token; once it (or ``STOP``) is set before or during
            compression the result is ``(False, "interrupted")`` and nothing is published
        accept_partial: Report a truncated download as a success. Without it a
            truncated verdict is a failure: the files stay in ``output_path`` for
            inspection, nothing is compressed, and the caller retries or records a failure,
            the same outcome the shared store gives an incomplete copy

    Returns:
        Tuple of (success, message); the message carries the completeness verdict
        (``"... complete (n of m spots)"``, ``"... truncated (n of m spots)"``, or
        ``"... unverified"`` when ``expected_spots`` is unknown).
    """
    # Check if files were actually created
    found = fastq_mod.fastq_files(temp_path)
    if not found:
        logger.error("No FASTQ files created despite successful command execution")
        # Kept for inspection, consistent with download_accession's failure paths.
        message = "No FASTQ files created"
        return False, f"{classify_download_error(message)}: {message}"

    # Move files to the final location
    # First ensure the output directory exists
    output_path.mkdir(parents=True, exist_ok=True)

    moved = []
    for file in found:
        dest = output_path / file.name
        shutil.move(str(file), str(dest))
        moved.append(dest)

    # Remove the temporary directory
    try:
        shutil.rmtree(temp_path)
    except OSError as e:
        logger.warning(f"Could not remove temp directory {temp_path}: {e}")

    logger.log(active_item_level(), f"Successfully downloaded: {len(found)} files")

    # Compute the verdict on the plain files first: counting reads in an uncompressed
    # file is cheaper, and the verdict message format must stay stable either way.
    verdict = fastq_mod.verify_download(output_path.name, output_path, expected_spots=expected_spots)
    message = _verdict_message(len(found), verdict)
    if verdict["verdict"] == "truncated":
        logger.warning(
            f"Truncated download for {output_path.name}: {verdict['reads_r1']} of "
            f"{verdict['expected_spots']} spots (ratio {verdict['ratio']})"
        )
        if not accept_partial:
            return False, message

    if compress:
        compression_failures, interrupted, corrupted = _compress_files(moved, output_path.name, num_threads, stop)
        if interrupted:
            # Not published: a folder with some files compressed and some not would be a
            # second layout for every later reader. The next run redoes the dump from the
            # prefetch archive, which is kept.
            return False, "interrupted"
        if corrupted:
            return False, "compression corrupted " + "; ".join(corrupted)
        if compression_failures:
            # Appended, never prepended, so parse_verdict_message's regex/substring checks
            # on the leading verdict text keep working unchanged.
            message += "; compression failed for " + ", ".join(compression_failures)

    return True, message


def _verdict_message(n_files: int, verdict: Dict[str, Any]) -> str:
    """The result message carrying the completeness verdict, in the format ``parse_verdict_message`` reads."""
    if verdict["verdict"] == "unverified":
        return f"Downloaded {n_files} files, unverified"
    spots_note = f"({verdict['reads_r1']} of {verdict['expected_spots']} spots)"
    return f"Downloaded {n_files} files, {verdict['verdict']} {spots_note}"


def _compress_files(
    moved: List[Path], name: str, num_threads: int, stop: Optional[threading.Event]
) -> Tuple[List[str], bool, List[str]]:
    """Compress each file in ``moved`` and read each result back.

    Returns the names that failed to compress, whether the run was interrupted before every
    file was compressed, and a description of each compressed file whose record count
    differs from the plain file's (``"<name>.gz: n of m reads"``). The read-back costs one
    gzip pass per file, as the store's sidecar does, and is what makes the verdict hold for
    the files published.
    """
    failures: List[str] = []
    corrupted: List[str] = []
    for file in moved:
        if stop_requested(stop):
            logger.info(f"Compression of {name} not started: the run was interrupted")
            return failures, True, corrupted
        try:
            before = fastq_mod.count_fastq_reads(file)
            target = fastq_mod.compress_fastq(file, num_threads, stop=stop)
        except fastq_mod._TOOL_ERRORS as e:
            if stop_requested(stop):
                # pigz terminated with the run's tools, or the Python fallback stopping on
                # the token: an interruption, not a compression failure.
                logger.info(f"Compression of {name} stopped: the run was interrupted")
                return failures, True, corrupted
            logger.warning(f"Could not compress {file}: {e}")
            failures.append(file.name)
            continue
        problem = _read_back_problem(target, before)
        if problem is not None:
            logger.error(f"Compressed {target.name} is not the plain file: {problem}")
            corrupted.append(f"{target.name}: {problem}")
    return failures, False, corrupted


def _read_back_problem(target: Path, expected_reads: int) -> Optional[str]:
    """Why the compressed ``target`` does not hold ``expected_reads`` records, or None when it does."""
    try:
        after = fastq_mod.count_fastq_reads(target)
    except (OSError, EOFError, zlib.error) as e:
        # A gzip stream that cannot be read back is a corrupt file, whatever its count.
        return f"unreadable ({e})"
    if after != expected_reads:
        return f"{after} of {expected_reads} reads"
    return None


def _fasterq_dump_args(
    source: str,
    temp_path: Path,
    temp_folder_path: Optional[Path],
    num_threads: int,
    using_prefetch: bool,
) -> List[str]:
    """Build the fasterq-dump argument list.

    ``source`` is the prefetched ``.sra`` archive when prefetch ran, otherwise the accession
    itself, which fasterq-dump then resolves and downloads on its own.
    """
    if using_prefetch:
        args = ["--split-3", "--skip-technical", "--threads", str(num_threads), "-O", str(temp_path)]
        tail = [source]
    else:
        args = ["--threads", str(num_threads), "--progress", source, "-O", str(temp_path)]
        tail = ["--split-3", "--skip-technical"]
    if temp_folder_path:
        args.extend(["--temp", str(temp_folder_path.absolute())])
    return args + tail


def _cached_sra_archive(acc_cache_dir: Path, accession: str) -> Path:
    """The archive prefetch wrote for ``accession``, preferring ``.sra`` over ``.sralite``.

    NCBI serves some runs only in the smaller ``.sralite`` format, which fasterq-dump reads
    just as well. When the folder holds neither, the ``.sra`` path is returned so
    fasterq-dump reports the missing file itself.
    """
    expected = acc_cache_dir / f"{accession}.sra"
    if expected.is_file():
        return expected
    candidates = sorted(p for p in acc_cache_dir.glob(f"{accession}.sra*") if p.is_file())
    return candidates[0] if candidates else expected


def _discard_cached_archive(cache_path: Path, accession: str, message: str, keep: bool = False) -> None:
    """Remove the prefetched ``.sra`` archive once its FASTQ files are on disk.

    The verdict is the one ``_handle_download_output`` put in ``message``. A truncated download
    drops the archive even when ``keep`` (``--keep-sra``) asks for it to be kept: keeping it
    would make the next attempt dump the same short archive again rather than fetch the run
    afresh. With ``keep``, a complete or unverified archive stays.
    """
    verdict = fastq_mod.parse_verdict_message(message)
    verdict_name = verdict.get("verdict") if verdict else None
    if verdict_name not in ("complete", "unverified", "truncated"):
        return
    if keep and verdict_name != "truncated":
        return
    cleanup_mod._safe_rmtree(cache_path / accession)
    if verdict_name == "truncated":
        logger.info(f"{accession}: archive removed so the next attempt fetches it again")


def _staging_root(output_folder: Union[str, Path], staging_folder: Optional[Union[str, Path]]) -> Path:
    """Where a download's ``<accession>_temp`` build folder goes: ``staging_folder`` if given,
    else beside the output. Registering it as an allowed root lets fasterq-dump write there."""
    root = Path(staging_folder or output_folder)
    SecureSubprocess.add_allowed_root(root)
    return root


def _log_prefetch_fallback(accession: str) -> None:
    """Say that fasterq-dump runs without prefetch: at INFO for a run's first accession, then at item level.

    Where the data comes from changes, so the first time is worth a line a user sees.
    """
    level = active_item_level() if _PREFETCH_FALLBACK_LOGGED.is_set() else logging.INFO
    _PREFETCH_FALLBACK_LOGGED.set()
    logger.log(
        level,
        f"prefetch not found on PATH; running fasterq-dump directly against {accession}, "
        "which downloads and dumps in one step",
    )


def download_accession(
    accession: str,
    output_folder: Union[str, Path],
    num_threads: int = 4,
    force: bool = False,
    temp_folder: Optional[Union[str, Path]] = None,
    expected_spots: Optional[int] = None,
    redownload_truncated: bool = False,
    sra_cache: Optional[Union[str, Path]] = None,
    use_prefetch: bool = True,
    keep_sra: bool = False,
    compress: bool = True,
    staging_folder: Optional[Union[str, Path]] = None,
    stop: Optional[threading.Event] = None,
    accept_partial: bool = False,
) -> Tuple[bool, str]:
    """
    Download a single SRA accession using prefetch + fasterq-dump --split-3.

    Args:
        accession: SRA accession to download
        output_folder: Folder to save the downloaded files
        num_threads: Number of threads to use for download
        force: If True, redownload even if files exist
        temp_folder: Directory for temporary files
        expected_spots: NCBI's recorded total_spots for this accession, used to verify
            completeness once the download finishes
        redownload_truncated: If True, treat an existing on-disk copy the same as ``force``
            (i.e. wipe it and redownload) rather than skipping it as already present; used
            for an accession whose registry verdict was "truncated"
        sra_cache: Directory prefetch downloads the ``.sra`` archive into; defaults to
            ``<output_folder>/.sra-cache``
        use_prefetch: If True and ``prefetch`` is on PATH, download the ``.sra`` archive
            first and run fasterq-dump against it; otherwise fasterq-dump is run directly
            against the accession, as before this became configurable
        keep_sra: If True, keep the downloaded ``.sra`` archive after a complete or
            unverified download rather than deleting it; a truncated archive is removed
            either way, and a redownload never reuses a cached archive
        compress: If True, gzip each downloaded FASTQ file once the completeness verdict
            has been computed
        staging_folder: Folder the ``<accession>_temp`` build directory is created in;
            defaults to ``output_folder``. The shared store points it at the store's own
            ``tmp`` folder so a half-written download never sits among finished datasets
        stop: The run's stop token. Once it (or the process-wide ``STOP``) is set no further
            tool is started and ``(False, "interrupted")`` is returned; the tools that are
            started are recorded under it for ``SecureSubprocess.terminate_children``
        accept_partial: Report a truncated download as a success rather than a failure
            (see ``_handle_download_output``)

    Returns:
        Tuple of (success, message)
    """
    require_accession(accession)
    output_path = Path(output_folder) / accession
    SecureSubprocess.add_allowed_root(Path(output_folder))
    staging_path = _staging_root(output_folder, staging_folder)
    cache_path = Path(sra_cache) if sra_cache else Path(output_folder) / ".sra-cache"

    # Check if already downloaded
    redownload = force or redownload_truncated
    if _check_existing_download(output_path, redownload):
        logger.log(active_item_level(), f"Skipping {accession}, FASTQ files already exist")
        return True, ALREADY_EXISTS

    # A redownload must not reuse a cached archive: prefetch treats an existing <acc>.sra as
    # already fetched, so a truncated archive would be dumped again and stay truncated. This
    # holds with keep_sra too, which keeps an archive after a download, never across a redownload.
    if redownload and (cache_path / accession).exists():
        cleanup_mod._safe_rmtree(cache_path / accession)
        logger.info(f"Removed the cached archive for {accession} so the redownload fetches it again")

    # Create a fresh temporary folder for download
    temp_path = staging_path / f"{accession}_temp"
    cleanup_mod._safe_rmtree(temp_path)
    temp_path.mkdir(parents=True, exist_ok=True)

    temp_folder_path = None
    try:
        logger.log(active_item_level(), f"Downloading SRA for {accession}")

        # Handle temp folder for fasterq-dump
        temp_folder_path = cleanup_mod._prepare_temp_folder(temp_folder)

        using_prefetch = use_prefetch and shutil.which("prefetch") is not None
        if use_prefetch and not using_prefetch:
            _log_prefetch_fallback(accession)

        if using_prefetch:
            SecureSubprocess.add_allowed_root(cache_path)
            _run_download_tool(
                "prefetch",
                ["-O", str(cache_path), "--max-size", settings.active().prefetch_max_size, "--progress", accession],
                stop,
            )
            source = str(_cached_sra_archive(cache_path / accession, accession))
        else:
            # Direct call against the accession, without going through prefetch's
            # on-disk .sra archive: used when use_prefetch is False, or prefetch is
            # not installed.
            source = accession

        # Run fasterq-dump command securely
        args = _fasterq_dump_args(source, temp_path, temp_folder_path, num_threads, using_prefetch)
        _run_download_tool("fasterq-dump", args, stop)

        # Handle download output
        success, message = _handle_download_output(
            temp_path,
            output_path,
            expected_spots=expected_spots,
            compress=compress,
            num_threads=num_threads,
            stop=stop,
            accept_partial=accept_partial,
        )

        # A truncated failure drops the archive too: prefetch would otherwise treat the short
        # archive as fetched and the retry would dump the same reads again.
        if using_prefetch and (success or "truncated" in message):
            _discard_cached_archive(cache_path, accession, message, keep=keep_sra)

        return success, message

    except _DownloadInterrupted:
        logger.info(f"Download of {accession} not started or continued: the run was interrupted")
        return False, "interrupted"

    except subprocess.CalledProcessError as e:
        logger.error(f"Error downloading {accession}: {e.stderr}")
        # <acc>_temp is kept for inspection; a new attempt starts clean (the _safe_rmtree
        # above always wipes it before the next fasterq-dump run, so this is not a resume).
        message = f"Download failed: {e.stderr}"
        return False, f"{classify_download_error(message)}: {message}"

    except SecurityError as e:
        logger.error(f"Security error downloading {accession}: {e}")
        # <acc>_temp is kept for inspection; a new attempt starts clean.
        message = f"Security error: {e}"
        return False, f"{classify_download_error(message)}: {message}"

    except (OSError, EOFError, ValueError, zlib.error, subprocess.TimeoutExpired, DataAccessError) as e:
        logger.error(f"Error downloading {accession}: {e}")
        # <acc>_temp is kept for inspection; a new attempt starts clean.
        message = f"Download failed: {str(e)}"
        return False, f"{classify_download_error(message)}: {message}"

    finally:
        # Clean up auto-created temp directory (from tempfile.mkdtemp)
        if temp_folder_path and not temp_folder:
            cleanup_mod._safe_rmtree(temp_folder_path)


def project_lock_path(fastq_folder: Union[str, Path], accession: str) -> Path:
    """The per-accession lock of a project without a shared store: ``<fastq>/.locks/<ACC>.lock``."""
    return Path(fastq_folder) / PROJECT_LOCK_FOLDER / f"{accession}.lock"


def project_scratch_path(fastq_folder: Union[str, Path], accession: str) -> Path:
    """fasterq-dump's scratch folder in a project without a shared store, when no temp folder is given.

    ``<fastq>/.prepare-tmp/<ACC>_fqtmp``: on the project's own filesystem rather than the
    system temporary folder, which on a cluster node is often small. The ``_fqtmp`` suffix
    makes it a transient folder (``is_transient_folder``), as is the store's scratch folder.
    """
    return Path(fastq_folder) / cleanup_mod.PROJECT_STAGING_FOLDER / f"{accession}_fqtmp"


def _project_lock_policy(accession: str, wait_seconds: float) -> LockPolicy:
    """The dataset lock's policy (stale threshold and heartbeat), applied to a project's accession."""
    stale, heartbeat = settings_or(dataset_lock_stale=DATASET_LOCK_STALE_SECONDS, lock_heartbeat=LOCK_HEARTBEAT_SECONDS)
    return LockPolicy(
        what=f"accession {accession}",
        stale_seconds=stale,
        wait_seconds=wait_seconds,
        poll_seconds=PROJECT_LOCK_POLL_SECONDS,
        heartbeat_seconds=heartbeat,
    )


def _stage_and_publish(
    accession: str,
    fastq: Path,
    lock: Path,
    num_threads: int,
    force: bool,
    temp_folder,
    download_kwargs: dict,
    stop: Optional[threading.Event] = None,
) -> Tuple[bool, str]:
    """Download ``accession`` into the staging folder and publish it to ``<fastq>/<ACC>``.

    Runs with the project lock held. The finished, verified and compressed folder
    ``<fastq>/.prepare-tmp/<ACC>`` becomes ``<fastq>/<ACC>`` with one rename, after
    ``verify_held`` confirms the lock is still this run's. A lost lock raises ``LockLost``
    and leaves the staged folder alone, since the new holder may be writing it.
    """
    staging = fastq / cleanup_mod.PROJECT_STAGING_FOLDER
    staged = staging / accession
    # Whatever an earlier interrupted attempt left staged is not a resume point: the download
    # would otherwise be skipped as "already exists" and that partial copy published.
    cleanup_mod._safe_rmtree(staged)
    if staged.exists():
        return False, f"unknown: could not clear the staged folder {staged}"
    success, message = download_accession(
        accession, staging, num_threads, force, temp_folder, staging_folder=staging, stop=stop, **download_kwargs
    )
    if not success:
        cleanup_mod._safe_rmtree(staged)
        return False, message
    if not staged.is_dir():
        return False, f"unknown: the download of {accession} reported success but staged no folder at {staged}"
    verify_held(lock)
    cleanup_mod.publish_folder(staged, fastq / accession, staging)
    return True, message


def _project_download(
    accession: str,
    output_folder: Union[str, Path],
    num_threads: int = 4,
    force: bool = False,
    temp_folder: Optional[Union[str, Path]] = None,
    *,
    lock_wait: float = 0.0,
    stop: Optional[threading.Event] = None,
    truncated: Optional[AbstractSet[str]] = None,
    **download_kwargs,
) -> Tuple[bool, str]:
    """Download ``accession`` into a project without a shared store, one process at a time.

    The worker ``download_sra`` uses when no store is configured; it takes
    ``download_accession``'s arguments. The per-accession lock ``<fastq>/.locks/<ACC>.lock``
    (the dataset lock's policy: heartbeat, stale takeover, no overall timeout unless
    ``lock_wait`` is positive) is held for the whole download. Inside it the folder is checked
    again, since another process may have finished the accession while this one waited, and
    an accession found complete is reported as "already exists". Otherwise the download is
    built under ``<fastq>/.prepare-tmp`` and published with one rename, so a second process
    never sees a half-moved or half-compressed folder; an existing folder is replaced only
    once its replacement is complete. The check is skipped for a forced download and for an
    accession in ``truncated`` (the registry's truncated verdicts); without that set,
    ``redownload_truncated`` skips it for every accession, as ``download_accession`` does.

    Without ``temp_folder``, fasterq-dump's scratch is ``project_scratch_path`` (under the
    staging folder), wiped under the lock before the download and again once it ends.

    The prefetch cache defaults to ``<fastq>/.sra-cache`` as before. A ``sra_cache`` folder
    the caller supplies and shares with another project is not covered by this lock (only a
    shared store locks across projects). ``stop`` is the run's stop token: when it (or the
    process-wide ``STOP``) is set, a wait for the lock ends with the result
    ``(False, "interrupted")``; it is passed on to ``download_accession``.
    """
    fastq = Path(output_folder)
    lock = project_lock_path(fastq, accession)
    if download_kwargs.get("sra_cache") is None:
        download_kwargs["sra_cache"] = fastq / ".sra-cache"
    if truncated is not None:
        download_kwargs["redownload_truncated"] = accession in truncated
    redownload = force or bool(download_kwargs.get("redownload_truncated"))
    policy = _project_lock_policy(accession, lock_wait)
    try:
        lock.parent.mkdir(parents=True, exist_ok=True)
        (fastq / cleanup_mod.PROJECT_STAGING_FOLDER).mkdir(parents=True, exist_ok=True)
        with held_lock(lock, policy, should_stop=functools.partial(stop_requested, stop)):
            if not redownload and _check_existing_download(fastq / accession, False):
                logger.log(active_item_level(), f"Skipping {accession}, FASTQ files already exist")
                return True, ALREADY_EXISTS
            scratch = None if temp_folder else project_scratch_path(fastq, accession)
            if scratch is not None:
                # Left by an interrupted attempt; this lock's holder is the only writer.
                cleanup_mod._safe_rmtree(scratch)
            try:
                return _stage_and_publish(
                    accession, fastq, lock, num_threads, force, scratch or temp_folder, download_kwargs, stop=stop
                )
            except LockLost:
                # The process that took the lock over uses the same scratch name: leave it alone.
                scratch = None
                raise
            finally:
                if scratch is not None:
                    cleanup_mod._safe_rmtree(scratch)
    except LockWaitStopped:
        logger.info(f"Stopped waiting for the lock on {accession}: the run was interrupted")
        return False, "interrupted"
    except LockHeld as e:
        logger.error(f"Gave up waiting for {accession}: {e}")
        return False, f"locked: {e}"
    except LockLost as e:
        logger.error(f"Not publishing {accession}: {e}")
        return False, f"lock lost: {e}"
    except (OSError, DataAccessError) as e:
        logger.error(f"Error downloading {accession}: {e}")
        message = f"Download failed: {e}"
        return False, f"{classify_download_error(message)}: {message}"


def fasterq_dump_version() -> str:
    """The installed fasterq-dump's version line (``fasterq-dump : 3.1.1``), or an empty string.

    Recorded in a store dataset's sidecar so a later reader knows which tool produced the
    files. Empty when the tool is missing or prints no version (``prepare.utils.tools.probe_tool``).
    """
    return tools.probe_tool("fasterq-dump").version_text
