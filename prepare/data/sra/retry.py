"""Parallel download of many accessions, tallying of their results, and the retry pass for failures."""

import functools
import logging
import subprocess
import threading
import time
import zlib
from concurrent.futures import CancelledError, Future, ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Set, Tuple, Union

from prepare.core import settings
from prepare.core.constants import FAILED_ACCESSIONS_FILE
from prepare.core.exceptions import NetworkError, PrepareError
from prepare.data.file_io import write_text_atomic
from prepare.data.sra import accession as accession_mod
from prepare.data.sra.space import INSUFFICIENT_SPACE_PREFIX
from prepare.data.sra.store_handoff import SETTLED_PREFIXES
from prepare.utils.progress import ProgressReporter, item_level
from prepare.utils.security import SecureSubprocess

if TYPE_CHECKING:  # pragma: no cover - typing only
    from prepare.data.sra.space import SpaceGuard

logger = logging.getLogger(__name__)


# What a download worker (download_accession, or the store downloader) raises for a dataset
# that could not be fetched: store and catalogue failures arrive as DataAccessError, file and
# FASTQ reading failures as OSError, EOFError, ValueError or (a corrupt gzip stream)
# zlib.error. Anything else is a programming error and propagates.
_DOWNLOAD_ERRORS = (PrepareError, OSError, EOFError, ValueError, zlib.error, subprocess.SubprocessError)

# The first word of every progress line, so a shared log file can be searched for one command's progress.
PROGRESS_LABEL = "download_sra"


# The message recorded for an accession a disk-full abort kept from starting, in either pass.
DISK_FULL_NOT_ATTEMPTED = "disk-full: not attempted"

# The message a worker returns when the run's stop token kept it from starting or finishing.
INTERRUPTED_MESSAGE = "interrupted"

# Per-accession timing of the last download attempt: (start as ISO 8601 UTC, seconds taken).
Timings = Dict[str, Tuple[str, float]]


def _instrumented(
    worker: Callable[..., Tuple[bool, str]], guard: Optional["SpaceGuard"], timings: Optional[Timings]
) -> Callable[..., Tuple[bool, str]]:
    """``worker`` with a free-space check before it starts and its run time recorded.

    Both download passes call the returned function in place of ``worker``, with the same
    arguments. Before the worker starts, ``guard.reserve`` must find room for the accession,
    waiting for running downloads to release theirs when that would be enough; the wait ends
    when the run's stop token (the ``stop`` keyword) or the process-wide ``STOP`` is set. When
    it finds no room (``insufficient-space: ...``) or is interrupted, the worker is not called
    and ``(False, <message>)`` is returned; an ``insufficient-space`` failure concerns that
    accession alone and does not stop the run. The reservation is released when the worker
    returns or raises. Each call the worker makes
    is timed into ``timings[accession] = (started, seconds)``, written in the worker thread
    before the result reaches the main thread's ``on_result``, so a later attempt overwrites an
    earlier one; a refused attempt removes the entry, so a refusal in the retry pass is not
    recorded with the time of the first pass's attempt. ``guard`` or ``timings`` may be None to
    skip that part.
    """

    def _run(accession: str, *args: Any, **kwargs: Any) -> Tuple[bool, str]:
        if guard is not None:
            stop = kwargs.get("stop")
            refusal = guard.reserve(accession, should_stop=lambda: accession_mod.stop_requested(stop))
            if refusal is not None:
                logger.log(
                    logging.INFO if refusal == "interrupted" else logging.ERROR,
                    "Not starting %s: %s",
                    accession,
                    refusal,
                )
                if timings is not None:
                    timings.pop(accession, None)
                return False, refusal
        started = datetime.now(timezone.utc)
        clock = time.monotonic()
        try:
            return worker(accession, *args, **kwargs)
        finally:
            if timings is not None:
                timings[accession] = (started.isoformat(timespec="seconds"), round(time.monotonic() - clock, 3))
            if guard is not None:
                guard.release(accession)

    return _run


def _progress_reporter(total: int) -> ProgressReporter:
    """A reporter over ``total`` downloads, at the interval of the ``progress_every`` setting.

    ``--progress-every`` reaches this through the runtime settings ``main()`` activates.
    """
    return ProgressReporter(PROGRESS_LABEL, total, settings.active().progress_every, logger=logger)


def _process_download_results(futures_results, accessions_to_download, download_results, failed_accessions):
    """
    Process download results from completed futures.

    Args:
        futures_results: List of (accession, future_result) tuples
        accessions_to_download: List of accessions that were attempted
        download_results: Dictionary to store results
        failed_accessions: List to store failed accessions

    Returns:
        Tuple of (successful_count, failed_count)
    """
    successful_count = 0
    failed_count = 0

    for accession, (success, message) in futures_results:
        download_results[accession] = message
        if success:
            successful_count += 1
        else:
            failed_count += 1
            failed_accessions.append(accession)
            logger.warning(f"Failed to download {accession}: {message}")

    logger.info("Downloaded %d of %d (%d failed)", successful_count, len(accessions_to_download), failed_count)

    return successful_count, failed_count


def _split_not_found(failed_accessions: List[str], download_results: Dict[str, Any]) -> Tuple[List[str], List[str]]:
    """Split failed accessions into (worth retrying, not worth retrying).

    An accession classified as not-found from its last attempt's message will not succeed on
    retry (the run genuinely does not exist, or the ID is invalid), and one the free-space guard
    refused does not fit even with no other download running; each is kept in the failed list
    rather than burning a retry round (and, for a refusal, a wait for running downloads) on it.
    A store copy kept but not linked (a message starting with one of ``SETTLED_PREFIXES``) is a
    settled outcome too: fetching it again in the same run gains nothing and uses up a refetch.
    """
    retry_batch: List[str] = []
    not_found: List[str] = []
    for accession in failed_accessions:
        message = accession_mod.result_body(str(download_results.get(accession, "")))
        if (
            INSUFFICIENT_SPACE_PREFIX in message
            or message.startswith(SETTLED_PREFIXES)
            or accession_mod.classify_download_error(message) == "not-found"
        ):
            not_found.append(accession)
        else:
            retry_batch.append(accession)
    return retry_batch, not_found


def _retry_failed_downloads(
    failed_accessions,
    max_retries,
    fastq_path,
    num_threads,
    temp_folder,
    download_results,
    on_result: Optional[Callable[[str, bool, str], None]] = None,
    expected_spots: Optional[Dict[str, int]] = None,
    redownload_truncated: bool = False,
    sra_cache: Optional[Union[str, Path]] = None,
    use_prefetch: bool = True,
    keep_sra: bool = False,
    compress: bool = True,
    downloader: Optional[Callable[..., Tuple[bool, str]]] = None,
    stop: Optional[threading.Event] = None,
    force: bool = False,
):
    """
    Retry failed downloads.

    Args:
        failed_accessions: List of accessions that failed
        max_retries: Maximum number of retry attempts
        fastq_path: Path to FASTQ directory
        num_threads: Number of threads to use
        temp_folder: Temporary folder path
        download_results: Dictionary to store results
        on_result: Optional callback invoked with (accession, success, message) after each retry
        expected_spots: Per-accession NCBI total_spots, used to verify completeness
        redownload_truncated: Forwarded to ``download_accession`` for each retry
        sra_cache: Forwarded to ``download_accession`` for each retry
        use_prefetch: Forwarded to ``download_accession`` for each retry
        keep_sra: Forwarded to ``download_accession`` for each retry
        compress: Forwarded to ``download_accession`` for each retry
        downloader: Callable used in place of ``download_accession``; the shared store
            passes one that links the project to the store's copy instead of downloading
            into the project folder
        stop: The run's stop token, forwarded to every retry. Once it (or the process-wide
            ``STOP``) is set no further round starts and the pause between rounds is skipped;
            a retry already queued in the current round returns "interrupted" from the worker
            without starting a tool
        force: The caller's ``--force``, forwarded to every retry so a forced run does not
            report an existing copy as "already exists" on the second attempt

    Returns:
        Tuple of (retried_successful, failed_accessions, abort_reason). ``abort_reason`` is
        ``"disk-full"`` when a retry attempt hit a disk-full error: the accession that hit it
        is recorded as failed and notified like any other failure, every other accession still
        queued in that round is marked failed with the message ``"disk-full: not attempted"``
        (and notified too) without ever calling ``download_accession``, and no further retry
        round runs. ``abort_reason`` is ``None`` when every round ran to completion normally.
    """
    if max_retries <= 0 or not failed_accessions or accession_mod.stop_requested(stop):
        return 0, failed_accessions, None

    logger.info(f"Retrying {len(failed_accessions)} failed downloads")
    retry_line_level = item_level(settings.active().progress_every)
    retry_count = 0
    retried_successful = 0
    expected_spots = expected_spots or {}
    abort_reason: Optional[str] = None

    for retry in range(max_retries):
        # A run stopped during the previous round starts no further round; within a round the
        # worker itself returns "interrupted" without starting a tool.
        if not failed_accessions or accession_mod.stop_requested(stop):
            break

        retry_batch, failed_accessions = _split_not_found(failed_accessions, download_results)

        if not retry_batch:
            break

        logger.info(f"Retry attempt {retry + 1}/{max_retries}")

        for index, accession in enumerate(retry_batch):
            if accession_mod.stop_requested(stop):
                # Not attempted and not recorded: the caller reports these as not completed.
                failed_accessions.extend(retry_batch[index:])
                break
            retry_count += 1
            attempt = functools.partial(
                downloader or accession_mod.download_accession,
                accession,
                fastq_path,
                num_threads,
                force=force,
                temp_folder=temp_folder,
                expected_spots=expected_spots.get(accession),
                redownload_truncated=redownload_truncated,
                sra_cache=sra_cache,
                use_prefetch=use_prefetch,
                keep_sra=keep_sra,
                compress=compress,
                stop=stop,
            )
            outcome = _retry_once(attempt, accession, retry, download_results, failed_accessions, on_result)
            if outcome is None:
                continue
            success, message = outcome

            if message == "interrupted":
                # Stopped before or while the tool ran: not an attempt, so not recorded or
                # notified; this and the rest of the batch are reported as not completed.
                failed_accessions.extend(retry_batch[index:])
                break

            download_results[accession] = f"Retry {retry + 1}: {message}"

            if success:
                retried_successful += 1
                logger.log(retry_line_level, f"Successfully downloaded {accession} on retry {retry + 1}")
                accession_mod._notify_result(on_result, accession, success, download_results[accession])
                continue

            failed_accessions.append(accession)
            logger.warning(f"Failed to download {accession} on retry {retry + 1}: {message}")
            accession_mod._notify_result(on_result, accession, success, download_results[accession])

            if accession_mod.classify_download_error(message) == "disk-full":
                abort_reason = "disk-full"
                _abort_disk_full(accession, retry_batch[index + 1 :], download_results, failed_accessions, on_result)
                break

        if abort_reason:
            break

        if failed_accessions and retry < max_retries - 1 and not accession_mod.stop_requested(stop):
            _back_off(2**retry, stop)

    return retried_successful, failed_accessions, abort_reason


def _retry_once(
    attempt: Callable[[], Tuple[bool, str]],
    accession: str,
    retry: int,
    download_results: Dict[str, Any],
    failed_accessions: List[str],
    on_result: Optional[Callable[[str, bool, str], None]],
) -> Optional[Tuple[bool, str]]:
    """Run one retry attempt; None when it raised, which is recorded and notified as a failure here."""
    try:
        return attempt()
    except _DOWNLOAD_ERRORS as e:
        failed_accessions.append(accession)
        logger.error(f"Error retrying download for {accession}: {e}")
        download_results[accession] = f"Retry {retry + 1} error: {str(e)}"
        accession_mod._notify_result(on_result, accession, False, download_results[accession])
        return None


def _abort_disk_full(
    accession: str,
    not_attempted: List[str],
    download_results: Dict[str, Any],
    failed_accessions: List[str],
    on_result: Optional[Callable[[str, bool, str], None]],
) -> None:
    """Mark the rest of a retry round as not attempted after ``accession`` filled the disk."""
    logger.error(f"Disk full while downloading {accession}; aborting remaining retries")
    for other in not_attempted:
        download_results[other] = DISK_FULL_NOT_ATTEMPTED
        failed_accessions.append(other)
        accession_mod._notify_result(on_result, other, False, download_results[other])


def _back_off(seconds: float, stop: Optional[threading.Event]) -> None:
    """Wait ``seconds`` between retry rounds, ending early once ``stop`` is set."""
    if stop is None:
        time.sleep(seconds)
    else:
        stop.wait(seconds)


def _handle_download_failure(fastq_path, failed_accessions):
    """
    Handle failed downloads by writing failed accessions to a file.

    Args:
        fastq_path: Path to FASTQ directory
        failed_accessions: List of failed accessions
    """
    failed_file = Path(fastq_path) / FAILED_ACCESSIONS_FILE
    if not failed_accessions:
        # A clean run leaves no stale list from an earlier run behind.
        if failed_file.exists():
            failed_file.unlink()
        return

    # Write failed accessions to file for easier retry; the command logs the retry hint.
    write_text_atomic(failed_file, "".join(f"{acc}\n" for acc in failed_accessions))
    logger.info(f"Failed accessions written to {failed_file}")


def _execute_parallel_downloads(
    accessions,
    fastq_path,
    num_threads,
    max_workers,
    force,
    temp_folder,
    download_results,
    failed_accessions,
    on_result: Optional[Callable[[str, bool, str], None]] = None,
    expected_spots: Optional[Dict[str, int]] = None,
    redownload_truncated: bool = False,
    sra_cache: Optional[Union[str, Path]] = None,
    use_prefetch: bool = True,
    keep_sra: bool = False,
    compress: bool = True,
    downloader: Optional[Callable[..., Tuple[bool, str]]] = None,
    stop: Optional[threading.Event] = None,
    progress: Optional[ProgressReporter] = None,
):
    """Download accessions concurrently and tally results. Returns (successful, failed).

    ``progress`` counts each result; when None a reporter is made here (before the first
    download starts, so a setting that does not parse stops the run first) and finished at the
    end of this pass. A caller that passes one logs its closing line itself, after a retry pass.

    ``downloader`` replaces ``download_accession`` when the project reads through a shared
    store; it takes the same arguments so the tally, retries and callbacks are unchanged.

    ``stop`` is the run's stop token (a new one when None), passed to every worker. On
    ``KeyboardInterrupt`` it is set, downloads not yet started are cancelled, this run's
    running prefetch or fasterq-dump children are terminated, and the executor waits for the
    downloads already running. Every download that then has a result not yet collected is
    collected (so ``on_result`` sees it), except one that returned ``"interrupted"``, and the
    interrupt is re-raised. Other runs in the process are not affected. The process-wide ``STOP`` and
    ``SecureSubprocess``'s stopping flag are left as they are: a run neither sets nor clears
    them, so an emergency stop set before the run still applies to it.

    The first result classified as disk-full (a tool that ran out of space) cancels every
    download not yet started: each is recorded as
    failed with ``"disk-full: not attempted"`` and notified at once, while downloads already
    running finish. The stop token is not used for this, so running tools are not stopped.

    Returns (successful, failed, abort_reason), ``abort_reason`` being ``"disk-full"`` after
    such an abort, else None.
    """
    if stop is None:
        stop = threading.Event()
    expected_spots = expected_spots or {}
    futures_results: list = []
    abort_reason: Optional[str] = None
    not_attempted: Set[str] = set()
    worker = downloader or accession_mod.download_accession
    finish_here = progress is None
    every = settings.active().progress_every  # read before the first download starts
    futures: Dict[Future, str] = {}
    collected: Set[Future] = set()
    interrupt: Optional[KeyboardInterrupt] = None
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        try:
            futures = {
                executor.submit(
                    worker,
                    acc,
                    fastq_path,
                    num_threads,
                    force,
                    temp_folder,
                    expected_spots=expected_spots.get(acc),
                    redownload_truncated=redownload_truncated,
                    sra_cache=sra_cache,
                    use_prefetch=use_prefetch,
                    keep_sra=keep_sra,
                    compress=compress,
                    stop=stop,
                ): acc
                for acc in accessions
            }
            if progress is None:
                progress = ProgressReporter(PROGRESS_LABEL, len(futures), every, logger=logger)
            for future in as_completed(futures):
                acc = futures[future]
                collected.add(future)
                if acc in not_attempted:
                    continue
                success, message = _collect_result(future, acc, futures_results, on_result, progress)
                if (
                    abort_reason is None
                    and not success
                    and accession_mod.classify_download_error(message) == "disk-full"
                ):
                    abort_reason = "disk-full"
                    logger.error("Disk full while downloading %s; downloads not yet started are cancelled", acc)
                    not_attempted = _cancel_pending(futures, futures_results, on_result, progress)
        except KeyboardInterrupt as exc:
            # Not re-raised here: leaving the block waits for the running workers, and one that
            # had already published its dataset (and was still releasing its lock or cleaning
            # up) finishes during that wait. It is collected below, before the re-raise.
            interrupt = exc
            stop.set()
            logger.warning("Interrupted; cancelling pending downloads and stopping running tools")
            executor.shutdown(wait=False, cancel_futures=True)
            SecureSubprocess.terminate_children(stop=stop)

    if interrupt is not None:
        # Every worker has returned by now, so each uncollected, uncancelled future is done.
        if progress is not None:
            finished = [f for f in futures if f not in collected and f.done() and not f.cancelled()]
            _collect_finished(finished, futures, not_attempted, futures_results, on_result, progress)
        raise interrupt

    if finish_here and progress is not None:
        progress.finish()
    successful, failed = _process_download_results(futures_results, accessions, download_results, failed_accessions)
    return successful, failed, abort_reason


def _classified_text(error: BaseException) -> str:
    """The exception's text, prefixed with ``network:`` when its class says so but its text does not."""
    text = str(error)
    if isinstance(error, NetworkError) and accession_mod.classify_download_error(text) != "network":
        return f"network: {text}"
    return text


def _collect_result(
    future: Future,
    acc: str,
    futures_results: list,
    on_result: Optional[Callable[[str, bool, str], None]],
    progress: ProgressReporter,
) -> Tuple[bool, str]:
    """Record one finished download on the main thread: tally entry, callback and progress.

    A worker that raised is tallied as ``(False, "Error: <text>")``, so the exception's text
    reaches the report, the failure class and the retry split like any other failure. A
    ``NetworkError`` whose text names no network symptom is prefixed with its class, so
    ``classify_download_error`` still reads it as retryable.
    """
    try:
        result = future.result()
    except (*_DOWNLOAD_ERRORS, CancelledError) as e:
        text = _classified_text(e)
        logger.error(f"Download failed for {acc}: {text}")
        futures_results.append((acc, (False, f"Error: {text}")))
        accession_mod._notify_result(on_result, acc, False, text)
        progress.update(ok=False)
        return False, text

    futures_results.append((acc, result))
    success, message = result
    if success:
        logger.log(progress.item_level, "%s: %s", acc, message)
    accession_mod._notify_result(on_result, acc, success, message)
    progress.update(ok=bool(success))
    return bool(success), message


def _collect_finished(
    finished: List[Future],
    futures: Dict[Future, str],
    not_attempted: Set[str],
    futures_results: list,
    on_result: Optional[Callable[[str, bool, str], None]],
    progress: ProgressReporter,
) -> None:
    """Collect downloads that finished after an interrupt was caught but were not yet collected.

    Called before the interrupt is re-raised, so nothing here may raise: a future whose worker
    raised something other than a download error (a programming error, or an interrupt delivered
    to the worker) is left out rather than re-raised over the ``KeyboardInterrupt``. A worker that
    returned ``(False, "interrupted")`` was cut short by the stop token and is left out too, so an
    interrupted re-download is not recorded as a new failure; a success or any other failure is
    collected.
    """
    for future in finished:
        acc = futures[future]
        if acc in not_attempted:
            continue
        error = future.exception()
        if error is not None and not isinstance(error, _DOWNLOAD_ERRORS):
            continue
        if error is None and _was_interrupted(future.result()):
            continue
        _collect_result(future, acc, futures_results, on_result, progress)


def _was_interrupted(result: Any) -> bool:
    """Whether a worker's result is the ``(False, "interrupted")`` a stopped run produces."""
    return isinstance(result, tuple) and len(result) == 2 and not result[0] and result[1] == INTERRUPTED_MESSAGE


def _cancel_pending(
    futures: Dict[Future, str],
    futures_results: list,
    on_result: Optional[Callable[[str, bool, str], None]],
    progress: ProgressReporter,
) -> Set[str]:
    """Cancel every download not yet started; record and notify each as ``"disk-full: not attempted"``.

    Returns the accessions cancelled, whose futures ``as_completed`` still yields later and the
    caller skips. A download already running cannot be cancelled and reports its own result.
    """
    cancelled: Set[str] = set()
    for future, acc in futures.items():
        if future.cancel():
            cancelled.add(acc)
            futures_results.append((acc, (False, DISK_FULL_NOT_ATTEMPTED)))
            accession_mod._notify_result(on_result, acc, False, DISK_FULL_NOT_ATTEMPTED)
            progress.update(ok=False)
    return cancelled


def _download_with_retries(
    accessions_to_download,
    fastq_path,
    num_threads,
    max_workers,
    force,
    temp_folder,
    max_retries,
    on_result: Optional[Callable[[str, bool, str], None]] = None,
    expected_spots: Optional[Dict[str, int]] = None,
    redownload_truncated: bool = False,
    sra_cache: Optional[Union[str, Path]] = None,
    use_prefetch: bool = True,
    keep_sra: bool = False,
    compress: bool = True,
    downloader: Optional[Callable[..., Tuple[bool, str]]] = None,
    stop: Optional[threading.Event] = None,
    guard: Optional["SpaceGuard"] = None,
    timings: Optional[Timings] = None,
) -> Tuple[int, int, List[str], Dict[str, Any], Optional[str]]:
    """Run the parallel downloads and optional retry pass.

    ``stop`` is the run's stop token (a new one when None), shared by both passes. The worker
    (``downloader``, or ``download_accession``) is wrapped once by ``_instrumented`` and both
    passes call the wrapped one: ``guard`` (a ``SpaceGuard``, or None for no free-space check)
    is asked for room before every attempt, and ``timings`` (when given) receives each
    accession's start time and duration.

    Returns (successful_count, failed_count, failed_accessions, download_results, abort_reason).
    ``abort_reason`` is ``"disk-full"`` when either pass hit a disk-full error (see
    ``_execute_parallel_downloads`` and ``_retry_failed_downloads``), else ``None``; after a
    first-pass abort no retry pass runs.
    """
    if stop is None:
        stop = threading.Event()
    failed_accessions: list = []
    download_results: dict = {}
    worker = _instrumented(downloader or accession_mod.download_accession, guard, timings)
    # Read before any download starts; its closing line follows the retry pass.
    progress = _progress_reporter(len(accessions_to_download))
    successful_count, failed_count, abort_reason = _execute_parallel_downloads(
        accessions_to_download,
        fastq_path,
        num_threads,
        max_workers,
        force,
        temp_folder,
        download_results,
        failed_accessions,
        on_result,
        expected_spots,
        redownload_truncated,
        sra_cache,
        use_prefetch,
        keep_sra,
        compress,
        worker,
        stop=stop,
        progress=progress,
    )

    if max_retries > 0 and failed_accessions and abort_reason is None:
        retried_successful, failed_accessions, abort_reason = _retry_failed_downloads(
            failed_accessions,
            max_retries,
            fastq_path,
            num_threads,
            temp_folder,
            download_results,
            on_result,
            expected_spots,
            redownload_truncated,
            sra_cache,
            use_prefetch,
            keep_sra,
            compress,
            worker,
            stop=stop,
            force=force,
        )
        successful_count += retried_successful
        failed_count -= retried_successful
        if retried_successful > 0:
            logger.info(f"Successfully downloaded {retried_successful} accessions on retry")
    progress.finish(ok=successful_count, failed=failed_count)
    return successful_count, failed_count, failed_accessions, download_results, abort_reason
