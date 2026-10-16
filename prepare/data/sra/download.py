"""Download of every SRA accession listed in a file into a project's FASTQ folder (``download_sra``)."""

import functools
import logging
import threading
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple, Union

from prepare.core import settings
from prepare.core.constants import MAX_CONCURRENT_DOWNLOADS
from prepare.core.exceptions import DataAccessError, PrepareError, ValidationError
from prepare.data.accessions import read_accessions_file
from prepare.data.file_io import ensure_directory
from prepare.data.sra import accession as accession_mod
from prepare.data.sra import fastq as fastq_mod
from prepare.data.sra import retry as retry_mod
from prepare.data.sra import space as space_mod
from prepare.data.sra import store_handoff as store_handoff_mod
from prepare.utils import resources

if TYPE_CHECKING:  # pragma: no cover - import cycle: prepare.store imports this package
    from prepare.store.layout import StorePaths

logger = logging.getLogger(__name__)


def default_max_workers(num_threads: int) -> int:
    """Size the download worker pool from the CPUs available to this job and per-download thread use.

    Each worker runs its own ``fasterq-dump`` using ``num_threads`` threads, so the pool is
    sized to roughly saturate the CPUs without oversubscribing them: divide the CPUs this
    process may use (``resources.available_cpus``: the affinity mask, a SLURM allocation, or
    the CPU count) by the per-download thread count, floor at 1 worker, and cap at the
    ``max_workers_cap`` setting (4 unless ``PREPARE_MAX_WORKERS_CAP`` or the config file
    changes it; downloads are limited by the network more than by CPUs) and at
    ``MAX_CONCURRENT_DOWNLOADS`` (a hard ceiling).
    """
    cpus = resources.available_cpus()
    cap = settings.active().max_workers_cap
    return min(MAX_CONCURRENT_DOWNLOADS, max(1, cpus // max(1, num_threads)), cap)


def _read_blacklist_files(blacklist_files):
    """
    Read accessions from blacklist files.

    Args:
        blacklist_files: List of blacklist file paths

    Returns:
        Set of blacklisted accessions
    """
    blacklisted_accessions: set = set()

    if not blacklist_files:
        return blacklisted_accessions

    for blacklist_file in blacklist_files:
        try:
            file_accessions = set()
            with open(blacklist_file, "r") as f:
                for line in f:
                    accession = line.split("#", 1)[0].strip()
                    if accession:
                        file_accessions.add(accession)
                        blacklisted_accessions.add(accession)
            logger.info(f"Read {len(file_accessions)} blacklisted accessions from {blacklist_file}")
        except (OSError, UnicodeDecodeError) as e:
            logger.warning(f"Error reading blacklist file {blacklist_file}: {e}")

    return blacklisted_accessions


def _check_existing_downloads(
    accessions: List[str],
    fastq_path: Path,
    force: bool,
    blacklisted_accessions: Optional[Set[str]] = None,
    truncated_accessions: Optional[Set[str]] = None,
) -> Tuple[List[str], List[str], List[str]]:
    """
    Check which accessions need downloading and which are already downloaded or blacklisted.

    This is an up-front hint that keeps finished accessions out of the worker pool; the
    decision that counts is taken again inside each accession's lock by the worker, since
    another process may finish an accession after this check.

    Args:
        accessions: List of accessions
        fastq_path: Path to FASTQ directory
        force: Whether to force redownload
        blacklisted_accessions: Set of blacklisted accessions
        truncated_accessions: Accessions whose registry verdict is "truncated"; treated like
            ``force`` for that one accession, so a partial download on disk is redownloaded
            rather than counted as already present

    Returns:
        Tuple of (already_downloaded, to_download, blacklisted)
    """
    already_downloaded = []
    to_download = []
    blacklisted = []

    if blacklisted_accessions is None:
        blacklisted_accessions = set()
    if truncated_accessions is None:
        truncated_accessions = set()

    for acc in accessions:
        if acc in blacklisted_accessions:
            blacklisted.append(acc)
            continue

        if not force and acc not in truncated_accessions and fastq_mod.accession_has_fastq(fastq_path / acc):
            already_downloaded.append(acc)
        else:
            to_download.append(acc)

    return already_downloaded, to_download, blacklisted


def _resolve_fastq_path(fastq_folder: Union[str, Path], dry_run: bool) -> Path:
    """Return the FASTQ output path, creating it unless in dry-run mode.

    In dry-run mode the folder is not created, but an existing non-directory
    at that location is still rejected.
    """
    fastq_path = Path(fastq_folder)
    if dry_run:
        if fastq_path.exists() and not fastq_path.is_dir():
            raise DataAccessError(f"{fastq_folder} exists but is not a directory")
        logger.info(f"Dry run mode: Would use {fastq_path} for downloads")
        return fastq_path
    return ensure_directory(fastq_folder)


def _log_download_run_summary(
    all_accessions: List[str],
    already_downloaded: List[str],
    blacklisted: List[str],
    successful_count: int,
    failed_count: int,
    download_results: Dict[str, str],
    abort_reason: Optional[str],
    failed_accessions: List[str],
    fastq_path: Path,
) -> None:
    """Log the summary for a completed (non-dry-run) ``download_sra`` call.

    A result the store served from a copy it already had (on the first pass or, behind a
    ``Retry <n>: `` prefix, on the retry pass) is counted in ``successful_count``,
    but it downloaded nothing this run; it is reported separately ("Linked from store") rather
    than folded into "Newly downloaded", which would overstate how much this run actually
    fetched. Kept out of ``download_sra`` itself to keep that function's branching down.
    """
    linked_count = sum(
        1
        for message in download_results.values()
        if accession_mod.result_body(message).startswith(store_handoff_mod.STORE_LINKED_PREFIX)
    )

    logger.info("Download summary:")
    logger.info(f"  Total accessions: {len(all_accessions)}")
    logger.info(f"  Already downloaded: {len(already_downloaded)}")
    logger.info(f"  Blacklisted: {len(blacklisted)}")
    logger.info(f"  Newly downloaded: {successful_count - linked_count}")
    if linked_count:
        logger.info(f"  Linked from store: {linked_count}")
    logger.info(f"  Failed downloads: {failed_count}")

    if abort_reason:
        logger.error(f"Download run aborted: {abort_reason}")
    if failed_count > 0:
        logger.warning("Some downloads failed. Use --force to retry or --max-retries to enable automatic retry.")
        retry_mod._handle_download_failure(fastq_path, failed_accessions)


def _space_guard(
    fastq_path: Path,
    temp_folder: Optional[Union[str, Path]],
    sra_cache: Optional[Union[str, Path]],
    store: Optional["StorePaths"],
    use_prefetch: bool,
    run_sizes: Optional[Mapping[str, Any]],
    min_free_gb: Optional[float],
    force: bool,
    accessions: List[str],
    expected_spots: Optional[Mapping[str, int]] = None,
) -> Optional[space_mod.SpaceGuard]:
    """The free-space guard for this run, or None when ``min_free_gb`` (or its setting) is 0.

    An accession the shared store already holds is linked, not downloaded, so it needs no
    space (unless ``force`` downloads it again); nor does an incomplete copy whose refetches
    stopped gaining reads, which is linked or refused without a download. A copy that is
    ``unverified`` but short against its ``expected_spots`` count is refetched, so it is
    counted. Logs the preflight warnings.
    """
    floor_gb = settings.active().min_free_gb if min_free_gb is None else min_free_gb
    if not floor_gb or floor_gb <= 0:
        return None
    counts = expected_spots or {}
    exempt = (
        {acc for acc in accessions if store_handoff_mod._store_settles(store, acc, counts.get(acc))}
        if store is not None and not force
        else set()
    )
    guard = space_mod.SpaceGuard(
        space_mod.download_locations(fastq_path, temp_folder, sra_cache, store),
        int(floor_gb * space_mod.GB),
        run_sizes or {},
        use_prefetch,
        exempt=exempt,
    )
    for warning in guard.preflight(accessions):
        logger.warning(warning)
    return guard


def download_sra(
    fastq_folder: Union[str, Path],
    accessions_file: Union[str, Path],
    max_downloads: Optional[int] = None,
    dry_run: bool = False,
    num_threads: int = 4,
    max_workers: int = 4,
    force: bool = False,
    max_retries: int = 1,
    temp_folder: Optional[Union[str, Path]] = None,
    blacklist: Optional[List[Union[str, Path]]] = None,
    blacklist_accessions: Optional[Set[str]] = None,
    on_result: Optional[Callable[[str, bool, str], None]] = None,
    expected_spots: Optional[Dict[str, int]] = None,
    redownload_truncated: bool = False,
    truncated_accessions: Optional[Set[str]] = None,
    sra_cache: Optional[Union[str, Path]] = None,
    use_prefetch: bool = True,
    keep_sra: bool = False,
    compress: bool = True,
    store: Optional["StorePaths"] = None,
    link_mode: str = "auto",
    accept_partial: bool = False,
    resume_partial: bool = True,
    store_metadata: Optional[Union[str, Path, Sequence[Union[str, Path]]]] = None,
    lock_wait: float = 0.0,
    stop: Optional[threading.Event] = None,
    run_sizes: Optional[Mapping[str, Any]] = None,
    min_free_gb: Optional[float] = None,
    timings: Optional[retry_mod.Timings] = None,
) -> Dict[str, Any]:
    """
    Download multiple SRA datasets.

    Args:
        fastq_folder: Folder to save downloaded FASTQ files
        accessions_file: File containing SRA accessions, one per line
        max_downloads: Maximum number of datasets to download
        dry_run: If True, only count accessions without downloading
        num_threads: Number of threads for each fasterq-dump
        max_workers: Number of parallel downloads
        force: If True, redownload even if files exist
        max_retries: Maximum number of retry attempts for failed downloads
        temp_folder: Directory to use for fasterq-dump temporary files
        blacklist: One or more files containing accessions to skip
        blacklist_accessions: An additional set of accessions to skip (e.g. registry exclusions),
            joined with any accessions read from ``blacklist``
        on_result: Optional callback invoked with (accession, success, message) on the main
            thread as each download (and each retry) completes
        expected_spots: Per-accession NCBI ``run_total_spots``, used to verify each download's
            completeness once it finishes; an accession missing from this mapping is reported
            as "unverified" rather than "complete"/"truncated"
        redownload_truncated: Forwarded to every ``download_accession`` call so a partial copy
            left on disk by a prior truncated attempt is wiped and redownloaded rather than
            reused
        truncated_accessions: Accessions whose registry verdict is "truncated"; excluded from
            ``already_downloaded`` so they are redownloaded even though files exist on disk
        sra_cache: Forwarded to every ``download_accession`` call; directory prefetch downloads
            the ``.sra`` archive into (defaults to ``<fastq_folder>/.sra-cache`` per accession)
        use_prefetch: Forwarded to every ``download_accession`` call; download via prefetch then
            fasterq-dump when True and prefetch is on PATH, else fasterq-dump directly
        keep_sra: Forwarded to every ``download_accession`` call; keep the ``.sra`` archive
            after a successful, verified download instead of deleting it
        compress: Forwarded to every ``download_accession`` call; gzip each downloaded FASTQ
            file once its completeness verdict has been computed
        store: Layout of a shared data store. When given, every dataset is downloaded once
            into ``<store>/sra/<ACC>`` and this project's ``fastq/<ACC>`` becomes a link to
            it; a dataset another project already downloaded is linked without any network
            call
        link_mode: How the project points at the store: ``auto``, ``relative``, ``absolute``
            or ``copy`` (see ``prepare.store.link.link_dataset``)
        accept_partial: Link a store copy whose download is incomplete instead of refusing
            it; only consulted when ``resume_partial`` is off
        resume_partial: Download an incomplete store copy again rather than refusing to use
            it
        store_metadata: One folder, or an ordered list of folders, searched for
            ``<ACC>_metadata.xml`` to record NCBI's spot count in the dataset's sidecar; the
            store's own metadata folder is always tried last
        lock_wait: Seconds to wait for another run's lock on an accession (the store's dataset
            lock, or without a store the project's ``<fastq>/.locks/<ACC>.lock``) before giving
            up on that accession; zero (the default) waits for as long as the other run keeps
            working, since a download legitimately takes hours
        stop: This run's stop token; a new one is made when None. It reaches every worker,
            ``download_accession`` and each prefetch, fasterq-dump or pigz child. Setting it
            prevents this run from starting further tools and ends its lock waits, but a tool
            already running keeps going until it ends; ``SecureSubprocess.terminate_children(
            stop=stop)`` also stops those. Either way only this run is affected; another run in
            the same process keeps going. The process-wide ``accession.STOP`` still stops every
            run from starting further tools
        run_sizes: Per-accession ``.sra`` size in bytes (NCBI's run size, from the registry),
            used by the free-space guard to estimate what each download needs
        min_free_gb: Free space, in GB, an accession without a known run size needs on every
            filesystem a download writes to; 0 turns the free-space guard off; None uses the
            ``min_free_gb`` setting (``--min-free-gb``, ``PREPARE_MIN_FREE_GB``, default 10)
        timings: When given, receives ``accession -> (started, seconds)`` for every download
            attempt that ran, filled before that attempt's ``on_result`` call

    Returns:
        Dictionary with download statistics

    Raises:
        DataAccessError: If the download fails
    """
    try:
        # Read and validate the accessions before anything touches the disk
        all_accessions = read_accessions_file(accessions_file)

        # Handle the output folder based on dry run status
        fastq_path = _resolve_fastq_path(fastq_folder, dry_run)

        logger.info(f"Found {len(all_accessions)} accessions in file")

        # Read blacklisted accessions (from files, plus any passed in directly, e.g. registry exclusions)
        blacklisted_accessions = _read_blacklist_files(blacklist)
        if blacklist_accessions:
            blacklisted_accessions |= set(blacklist_accessions)
        if blacklisted_accessions:
            logger.info(f"Found total of {len(blacklisted_accessions)} blacklisted accessions")

        # Check which accessions need downloading
        already_downloaded, accessions_to_download, blacklisted = _check_existing_downloads(
            all_accessions, fastq_path, force, blacklisted_accessions, truncated_accessions
        )

        logger.info(f"{len(already_downloaded)} accessions already downloaded")
        logger.info(f"{len(blacklisted)} accessions blacklisted")
        logger.info(f"{len(accessions_to_download)} accessions need downloading")

        # Accessions that --max-downloads would cut off, computed before any truncation so a
        # dry run can report them too.
        skipped_accessions: List[str] = []
        if max_downloads is not None and max_downloads < len(accessions_to_download):
            skipped_accessions = accessions_to_download[max_downloads:]

        if dry_run:
            logger.info(f"Dry run: would download {len(accessions_to_download)} accessions")
            return {
                "total": len(all_accessions),
                "already_downloaded": len(already_downloaded),
                "blacklisted": len(blacklisted),
                "to_download": len(accessions_to_download),
                "successful": 0,
                "failed": 0,
                "already_downloaded_accessions": sorted(str(a) for a in already_downloaded),
                "blacklisted_accessions": sorted(str(a) for a in blacklisted),
                "skipped_accessions": sorted(str(a) for a in skipped_accessions),
                # No download ran, so nothing could abort; the key is present either way so
                # callers can read it without knowing which mode produced the stats.
                "aborted": None,
            }

        # Limit number of downloads if specified
        if max_downloads is not None and max_downloads < len(accessions_to_download):
            logger.info(f"Limiting to {max_downloads} downloads")
            accessions_to_download = accessions_to_download[:max_downloads]

        # With a shared store, every download goes through it: the store keeps the only copy
        # and the project gets a link to it. Without one, each accession is downloaded under
        # the project's own per-accession lock and published with one rename.
        downloader = store_handoff_mod._store_downloader(
            store, link_mode, accept_partial, resume_partial, store_metadata, lock_wait
        ) or functools.partial(
            accession_mod._project_download,
            lock_wait=lock_wait,
            truncated=frozenset(truncated_accessions) if truncated_accessions is not None else None,
        )

        guard = _space_guard(
            fastq_path,
            temp_folder,
            sra_cache,
            store,
            use_prefetch,
            run_sizes,
            min_free_gb,
            force,
            accessions_to_download,
            expected_spots=expected_spots,
        )

        # Download accessions in parallel, with an optional retry pass
        successful_count, failed_count, failed_accessions, download_results, abort_reason = (
            retry_mod._download_with_retries(
                accessions_to_download,
                fastq_path,
                num_threads,
                max_workers,
                force,
                temp_folder,
                max_retries,
                on_result,
                expected_spots,
                redownload_truncated,
                sra_cache,
                use_prefetch,
                keep_sra,
                compress,
                downloader,
                stop=stop if stop is not None else threading.Event(),
                guard=guard,
                timings=timings,
            )
        )

        _log_download_run_summary(
            all_accessions,
            already_downloaded,
            blacklisted,
            successful_count,
            failed_count,
            download_results,
            abort_reason,
            failed_accessions,
            fastq_path,
        )

        download_stats = {
            "total": len(all_accessions),
            "already_downloaded": len(already_downloaded),
            "blacklisted": len(blacklisted),
            "successful": successful_count,
            "failed": failed_count,
            "failed_accessions": failed_accessions,
            "results": download_results,
            "already_downloaded_accessions": sorted(str(a) for a in already_downloaded),
            "blacklisted_accessions": sorted(str(a) for a in blacklisted),
            "skipped_accessions": sorted(str(a) for a in skipped_accessions),
            "aborted": abort_reason,
        }

        return download_stats

    except ValidationError:
        raise
    except (OSError, ValueError, PrepareError) as e:
        raise DataAccessError(f"Downloading SRA data: {e}") from e
