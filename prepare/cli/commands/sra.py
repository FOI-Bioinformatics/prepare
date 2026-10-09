"""
SRA-related CLI commands.
"""

import argparse
import functools
import time
from datetime import datetime, timezone
from typing import Callable, Dict, List, Optional, Set, Tuple

from prepare.cli.base import BaseCommand
from prepare.cli.commands.sra_verdicts import present_verdicts, recorded_verdict, registry_inputs
from prepare.cli.commands.sra_verdicts import store_sidecar_verdict
from pathlib import Path

from prepare.core.constants import FAILED_ACCESSIONS_FILE
from prepare.core.exceptions import DataAccessError, ExitCode, PrepareError
from prepare.core.settings import setting_for
from prepare.data import registry_blocks as rb
from prepare.data import run_log
from prepare.data.registry import Registry, load_registry, project_root, record_download, update_linked
from prepare.data.registry_batch import RegistryBatch, registry_batch
from prepare.data.registry_timing import set_download_timing
from prepare.data.sra import (
    ALREADY_EXISTS,
    SETTLED_PREFIXES,
    STORE_LINKED_PREFIX,
    classify_download_error,
    default_max_workers,
    download_sra,
    parse_verdict_message,
    transient_bytes,
)
from prepare.data.sra import run_report as rr
from prepare.data.sra.run_report import RunOutcomes, Timings, attempt_timing, write_run_document
from prepare.store.layout import StorePaths, store_paths
from prepare.store.link import LINK_MODES, is_store_link
from prepare.store.resolve import resolve_store_root
from prepare.store.usage import ensure_project_identity, record_usage_many
from prepare.utils import resources
from prepare.utils.termination import graceful_termination
from prepare.utils.tools import require_tools

# Marker the data layer puts in a result message for a dataset this run downloaded and
# saved into the store (as opposed to STORE_LINKED_PREFIX, imported above, for one the
# store already held).
STORE_SAVED_SUFFIX = "; stored"

# Kept .sra-cache archives and <ACC>_temp build folders bigger than this, summed across a
# run's candidate folders, are worth a warning: they are easy to forget about and can
# quietly use up a lot of disk.
TRANSIENT_BYTES_WARN_THRESHOLD = 1024**3

# Seconds to wait before retrying a final registry flush that failed (a lock held by another
# process, for example); module-level so tests can set it to zero.
FINAL_FLUSH_RETRY_SECONDS = 1.0


def _max_downloads(value: str) -> int:
    """argparse type for --max-downloads: rejects zero and negative values with a clear message."""
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError(f"--max-downloads must be a positive integer, got {value!r}")
    return parsed


# The download command's own name for the context every command now runs under (BaseCommand.run).
# Nested inside it, it installs nothing and shares its Termination; called on its own (tests that
# call ``execute`` directly) it installs the handlers itself.
_termination_raises_interrupt = graceful_termination


class DownloadSraCommand(BaseCommand):
    """Command for downloading SRA datasets."""

    @property
    def name(self) -> str:
        return "download_sra"

    @property
    def help(self) -> str:
        return "Download SRA datasets"

    @property
    def group(self) -> str:
        return "Reads"

    def records_run(self, args: argparse.Namespace) -> bool:
        return not getattr(args, "dry_run", False)

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--fastq-folder",
            default="fastq",
            help="Folder to save downloaded FASTQ files",
        )
        parser.add_argument(
            "--accessions-file",
            required=True,
            help="File containing SRA accessions, one per line",
        )
        parser.add_argument(
            "--max-downloads",
            type=_max_downloads,
            default=None,
            help="Maximum number of datasets to download",
        )
        parser.add_argument(
            "--num-threads",
            type=int,
            default=4,
            help="Number of threads for each fasterq-dump",
        )
        parser.add_argument(
            "--max-workers",
            type=int,
            default=None,
            help=(
                "Number of parallel downloads (default: CPUs available to this job / --num-threads, "
                "at most 4; change the cap with PREPARE_MAX_WORKERS_CAP)"
            ),
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Calculate number of accessions without downloading",
        )
        parser.add_argument("--force", action="store_true", help="Force redownload even if files exist")
        parser.add_argument(
            "--max-retries",
            type=int,
            default=1,
            help="Maximum number of retry attempts for failed downloads",
        )
        parser.add_argument(
            "--temp-folder",
            help="Directory to use for fasterq-dump temporary files (must be writable; default: PREPARE_TEMP_FOLDER)",
        )
        parser.add_argument(
            "--blacklist",
            nargs="+",
            help="One or more files containing blacklisted accessions, one per line",
        )
        parser.add_argument(
            "--report-file",
            default=None,
            help=(
                "Write a CSV of accession,status,message,seconds,reason,attempts after the run, an "
                "interrupted one included (statuses: downloaded, failed, already_present, blacklisted, "
                "skipped; reason for a failed row only); a summary of every run that is not a dry run "
                "is also written to <fastq-folder>/download_run.json"
            ),
        )
        parser.add_argument(
            "--registry",
            default=None,
            help="Path to the project registry file (defaults to the nearest prepare_registry.json)",
        )
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")
        parser.add_argument(
            "--verify-downloads",
            dest="verify_downloads",
            action="store_true",
            default=True,
            help="Verify each download's read count against NCBI's recorded total spots (default: on)",
        )
        parser.add_argument(
            "--no-verify-downloads",
            dest="verify_downloads",
            action="store_false",
            help="Skip completeness verification against NCBI spot counts",
        )
        parser.add_argument(
            "--redownload-truncated",
            action="store_true",
            help="Redownload accessions whose registry verdict is 'truncated' rather than skipping them",
        )
        parser.add_argument(
            "--sra-cache",
            default=None,
            help=(
                "Directory for prefetch's downloaded .sra archives (default: <fastq-folder>/.sra-cache); "
                "without a shared store, a directory shared by two projects is not locked"
            ),
        )
        parser.add_argument(
            "--no-prefetch",
            dest="use_prefetch",
            action="store_false",
            default=True,
            help="Run fasterq-dump directly against the accession instead of prefetch then fasterq-dump",
        )
        parser.add_argument(
            "--keep-sra",
            action="store_true",
            help=(
                "Keep the downloaded .sra archive after a complete or unverified download; a truncated archive "
                "is removed, and a redownload (--force, --redownload-truncated) always fetches the archive again"
            ),
        )
        parser.add_argument(
            "--compress",
            dest="compress",
            action="store_true",
            default=True,
            help="Gzip each downloaded FASTQ file (default: on)",
        )
        parser.add_argument(
            "--no-compress",
            dest="compress",
            action="store_false",
            help="Leave downloaded FASTQ files uncompressed",
        )
        parser.add_argument(
            "--link-mode",
            choices=list(LINK_MODES),
            default="auto",
            help=(
                "How this project points at a dataset in the shared store: a relative or "
                "absolute symlink, a copy of the folder, or auto (relative when the store "
                "and the project share a parent folder)"
            ),
        )
        parser.add_argument(
            "--accept-partial",
            action="store_true",
            help="Use a store dataset whose download is incomplete instead of refusing it",
        )
        parser.add_argument(
            "--no-resume-partial",
            dest="resume_partial",
            action="store_false",
            default=True,
            help="Do not download an incomplete store dataset again",
        )
        parser.add_argument(
            "--lock-wait",
            dest="lock_wait",
            type=float,
            default=None,
            help=(
                "Seconds to wait for another run's download of the same accession before giving up "
                "on it, whether through a shared store or, without one, this project's own "
                "per-accession lock (default: PREPARE_LOCK_WAIT, config [runtime] lock_wait, or 0, "
                "wait for as long as the other run keeps working)"
            ),
        )
        parser.add_argument(
            "--timeout",
            dest="timeout",
            type=float,
            default=None,
            help=(
                "Seconds before prefetch or fasterq-dump is stopped; 0 (the default) means no limit "
                "(default: PREPARE_TIMEOUT, config [runtime] timeout, or 0)"
            ),
        )

        parser.add_argument(
            "--min-free-gb",
            dest="min_free_gb",
            type=float,
            default=None,
            help=(
                "Free space, in GB, a download of unknown size needs on each filesystem it writes to; "
                "one with a registry run size needs about 8 times that size for fasterq-dump's temporary "
                "files and 10 times for the FASTQ folder (uncompressed and gzip files). A download that "
                "does not fit while others run waits for them; one that would not fit even alone fails "
                "with insufficient-space and the others continue. Only a tool's own out-of-space error "
                "stops the run. 0 turns the check off "
                "(default: PREPARE_MIN_FREE_GB, config [runtime] min_free_gb, or 10)"
            ),
        )

    def _log_dry_run_summary(self, args: argparse.Namespace, stats: dict) -> None:
        """Log the summary for a dry run."""
        self.logger.info(f"Dry run: would download {stats['to_download']} of {stats['total']} datasets")
        self.logger.info(f"  {stats['already_downloaded']} datasets would be skipped (already downloaded)")
        if stats.get("blacklisted", 0) > 0:
            self.logger.info(f"  {stats['blacklisted']} datasets would be skipped (blacklisted)")
        if stats.get("to_download", 0) > 0:
            self.logger.info(f"  Output folder would be: {args.fastq_folder}")
            if args.max_downloads is not None:
                self.logger.info(f"  Limited to {args.max_downloads} downloads")

    def _report_failed_downloads(self, args: argparse.Namespace, stats: dict) -> None:
        """Point at the retry file the data layer already wrote.

        The data layer (prepare.data.sra's download_sra) already logs "Some downloads
        failed..." itself; logging it again here would print it twice on a failed run.
        """
        if not stats.get("failed_accessions"):
            return
        failed_file = Path(args.fastq_folder) / FAILED_ACCESSIONS_FILE
        self.logger.info(
            f"To retry only failed accessions: prepare download_sra "
            f"--accessions-file {failed_file} "
            f"--fastq-folder {args.fastq_folder}"
        )

    @staticmethod
    def _failure_exit_code(stats: dict) -> int:
        """4 (retryable) when every failed accession failed for a network reason, else 1.

        A not-found, disk-full, lock or interrupted failure among them means a rerun alone
        would not succeed, so the run is a plain failure.
        """
        results = stats.get("results", {})
        failed = stats.get("failed_accessions", [])
        if failed and all(classify_download_error(results.get(acc, "")) == "network" for acc in failed):
            return int(ExitCode.TRANSIENT)
        return int(ExitCode.FAILURE)

    def _write_run_summary(
        self, args: argparse.Namespace, outcomes: RunOutcomes, started: datetime, code: int, aborted: Optional[str]
    ) -> None:
        """Write the ``--report-file`` CSV and ``<fastq>/download_run.json`` for a run that is not a dry run.

        Called after every such run, an interrupted one or one whose final registry flush failed
        included; the JSON document is left out for a run that stopped before the download stage
        when the FASTQ folder does not exist. Statistics come from ``download_sra`` when it returned
        them, else from the outcomes observed. A summary that cannot be written is logged and leaves
        ``code`` as it is.
        """
        stats = outcomes.stats if outcomes.stats is not None else rr.stats_from_outcomes(outcomes)
        rows = rr.report_rows(stats, outcomes.timings, outcomes.attempts)
        try:
            if args.report_file:
                rr.write_report_csv(args.report_file, rows)
                self.logger.info("Download report written to %s", args.report_file)
            document = rr.run_document(
                stats=stats,
                outcomes=outcomes,
                started=started,
                finished=datetime.now(timezone.utc),
                exit_code=code,
                aborted=aborted or stats.get("aborted"),
                settings=rr.run_settings(args),
                paths=rr.run_paths(args, outcomes.stats),
            )
            run_log.note_run(args, *rr.run_log_entries(document, rows))
            # A run stopped by its preflight checks does not create the FASTQ folder for this file.
            if outcomes.reached_download or Path(args.fastq_folder).is_dir():
                write_run_document(args.fastq_folder, document)
        except (OSError, PrepareError) as e:
            self.logger.warning("Could not write the download run summary: %s", e)

    @staticmethod
    def _record_skip(reg: Registry, acc: str, message: str, fastq_dir: Path) -> None:
        """Record a skipped accession, unless its record already says the reads are downloaded.

        A blacklist entry or a `--max-downloads` cut-off must never erase the files and
        sizes of an accession that was downloaded on an earlier run.
        """
        if (rb.download_block(reg, acc) or rb.DownloadBlock()).state == "downloaded":
            return
        record_download(reg, acc, "skipped", fastq_dir, message)
        set_download_timing(reg, acc, None, None)

    def _record_run_outcomes(
        self,
        args: argparse.Namespace,
        stats: dict,
        fastq_dir: Path,
        store: Optional[StorePaths] = None,
        expected_spots: Optional[Dict[str, int]] = None,
    ) -> None:
        """Record the outcomes the download loop could not report, in one registry transaction.

        Already-downloaded, blacklisted and skipped accessions are queued on one
        ``registry_batch`` that writes once, when the loops are done: a transaction rewrites the
        whole registry, so one per accession took about a second each on a large project.
        Every already-downloaded accession the shared store backs is also recorded as
        ``"linked"`` usage in the store catalogue, in one batched write after the registry
        write (``record_usage_many``) rather than one lock per accession, using the registry
        that write produced. An already-downloaded accession in a folder of its own without a
        ``downloaded`` record (``present_verdicts``) is verified against ``expected_spots`` first.
        """
        usage_rows = []
        spots = expected_spots or {}
        # Store links and sidecar verdicts are read here, before the batch writes, so the one
        # registry lock is never held across thousands of reads on a slow filesystem.
        present = []
        for acc in stats.get("already_downloaded_accessions", []):
            from_store = store is not None and is_store_link(fastq_dir / acc, store)
            present.append((acc, from_store, self._sidecar_completeness(store, acc) if from_store else None))
            if from_store:
                usage_rows.append((acc, "", "linked", "already downloaded"))
        local = [acc for acc, from_store, _ in present if not from_store]
        verdicts = present_verdicts(fastq_dir, local, load_registry(args.registry), spots) if local else {}
        with registry_batch(args.registry, flush_every=None, flush_seconds=None) as batch:
            for acc, from_store, complete in present:
                mutation = functools.partial(
                    self._record_present,
                    acc=acc,
                    fastq_dir=fastq_dir,
                    from_store=from_store,
                    complete=verdicts.get(acc) if complete is None else complete,
                    spots=spots.get(acc),
                )
                batch.apply(mutation, label=acc)
            for acc in stats.get("blacklisted_accessions", []):
                mutation = functools.partial(self._record_skip, acc=acc, message="blacklisted", fastq_dir=fastq_dir)
                batch.apply(mutation, label=acc)
            for acc in stats.get("skipped_accessions", []):
                mutation = functools.partial(self._record_skip, acc=acc, message="--max-downloads", fastq_dir=fastq_dir)
                batch.apply(mutation, label=acc)
        if usage_rows:
            usage_registry = batch.registry if batch.registry is not None else load_registry(args.registry)
            record_usage_many(store, usage_registry, usage_rows)

    @staticmethod
    def _record_present(
        reg: Registry,
        acc: str,
        fastq_dir: Path,
        from_store: bool,
        complete: Optional[dict],
        spots: Optional[int] = None,
    ) -> None:
        """Record an accession found already downloaded, unless its record already says so.

        ``complete`` is the store sidecar's completeness verdict for a store-linked accession, or
        ``present_verdicts``'s for one in a folder of its own, merged with the one on file by
        ``recorded_verdict`` against ``spots``; the caller reads either before the lock is taken.
        """
        if from_store:
            ensure_project_identity(reg)
        if (rb.download_block(reg, acc) or rb.DownloadBlock()).state == "downloaded":
            return
        complete = recorded_verdict(reg, acc, complete, spots, linked=from_store)
        record_download(
            reg,
            acc,
            "downloaded",
            fastq_dir,
            attempt=False,
            complete=complete,
            source="store" if from_store else None,
            store_name=acc if from_store else None,
        )
        if from_store:
            update_linked(reg, acc, add=True)

    def _resolve_max_workers(self, args: argparse.Namespace) -> int:
        """Resolve --max-workers, falling back to a CPU-derived default.

        Oversubscription is only worth warning about when the user chose the worker count:
        the derived default is already floored at one worker, so a single many-threaded
        download on a small machine is nothing the user could have set differently.
        """
        explicit = args.max_workers is not None
        max_workers = args.max_workers if explicit else default_max_workers(args.num_threads)
        cpu_count = resources.available_cpus()
        if explicit and max_workers * args.num_threads > cpu_count:
            self.logger.warning(
                "--max-workers %d x --num-threads %d = %d threads requested, which exceeds "
                "the %d CPUs available to this job; downloads may be slower than expected",
                max_workers,
                args.num_threads,
                max_workers * args.num_threads,
                cpu_count,
            )
        return max_workers

    def _resolve_store(self, args: argparse.Namespace, project_registry: Registry) -> Optional[StorePaths]:
        """Resolve the shared data store (if any), log it, and return its on-disk layout."""
        store_root = resolve_store_root(args.data_root, rb.store_block(project_registry).root)
        if store_root is None:
            return None
        self.logger.info("Using shared data store at %s", store_root)
        return store_paths(store_root)

    def _result_recorder(
        self,
        args: argparse.Namespace,
        fastq_dir: Path,
        store: Optional[StorePaths],
        batch: RegistryBatch,
        timings: Optional[Timings] = None,
        expected_spots: Optional[Dict[str, int]] = None,
    ) -> Callable[[str, bool, str], None]:
        """The callback the download loop uses to record each accession's outcome.

        Each outcome is queued on ``batch`` (which ``_run`` opens around the download call and
        flushes every 50 results, every 30 seconds, and on exit including an interrupt) rather
        than written in a transaction of its own.

        A dataset linked from the shared store is recorded as downloaded without counting an
        attempt against it, since no download ran; one this run downloaded into the store
        counts as an attempt like any other. Both are added to the project's list of linked
        datasets, and recorded as store catalogue usage: ``"linked"`` for a dataset the store
        already held, ``"downloaded"`` for one this run saved into it. The usage rows are
        written after each registry flush, once its lock is released: the catalogue has its own
        lock, and waiting for it while holding the project's registry lock can time a
        concurrent writer's registry write out.

        ``timings`` is the dict ``download_sra`` fills before each ``on_result`` call; the
        attempt's start and seconds are recorded with the outcome, and a dataset linked from
        the store, which no download produced, has any earlier time removed. The verdict recorded
        is ``recorded_verdict`` of the message's (or, for a store dataset whose message carries
        none, the sidecar's) and the one on file, with ``expected_spots`` as the spot count.
        """
        timings = timings if timings is not None else {}
        spot_counts = expected_spots or {}
        usage_rows: List[Tuple[str, str, str, str]] = []

        def _record_usage(reg: Registry) -> None:
            if usage_rows:
                record_usage_many(store, reg, list(usage_rows))
                usage_rows.clear()

        batch.add_flush_hook(_record_usage)

        def _record_result(accession: str, success: bool, message: str) -> None:
            if success and message == ALREADY_EXISTS:
                # Found in place (typically after waiting for another run that downloaded it):
                # not an attempt, and an existing record keeps the other run's message.
                present = functools.partial(
                    self._record_present, acc=accession, fastq_dir=fastq_dir, from_store=False, complete=None
                )
                batch.apply(present, label=accession)
                return
            linked = bool(success) and message.startswith(STORE_LINKED_PREFIX)
            from_store = linked or (bool(success) and message.endswith(STORE_SAVED_SUFFIX))
            if from_store:
                usage_rows.append((accession, "", "linked" if linked else "downloaded", message))
            complete = parse_verdict_message(message) if success else None
            # A store link or a store copy kept unlinked (SETTLED_PREFIXES) records its sidecar's verdict
            # (sra_verdicts.store_verdict); the sidecar is read here, outside the flush's lock.
            use_sidecar = complete is None and (from_store or (not success and message.startswith(SETTLED_PREFIXES)))
            sidecar = self._sidecar_completeness(store, accession) if use_sidecar else None
            started, seconds = attempt_timing(accession, success, message, timings) or (None, None)

            def _mutation(reg: Registry) -> None:
                if from_store:
                    ensure_project_identity(reg)
                verdict = recorded_verdict(
                    reg, accession, sidecar or complete, spot_counts.get(accession), bool(sidecar)
                )
                record_download(
                    reg,
                    accession,
                    "downloaded" if success else "failed",
                    fastq_dir,
                    message,
                    attempt=not linked,
                    complete=verdict,
                    source="store" if from_store else None,
                    store_name=accession if from_store else None,
                )
                set_download_timing(reg, accession, started, seconds)
                if from_store:
                    update_linked(reg, accession, add=True)

            batch.apply(_mutation, label=accession)

        return _record_result

    @staticmethod
    def _sidecar_completeness(store: Optional[StorePaths], accession: str) -> Optional[dict]:
        """The completeness verdict recorded in the store's sidecar for ``accession``, or None.

        None when there is no store, or no sidecar (not yet catalogued, or unreadable);
        ``read_sidecar`` already logs a warning for the latter case.
        """
        return store_sidecar_verdict(store, accession)

    def _transient_folders(self, args: argparse.Namespace, fastq_dir: Path, store: Optional[StorePaths]) -> Set[Path]:
        """Folders where ``download_accession`` can leave ``.sra-cache`` archives or
        ``<ACC>_temp`` build directories behind: the FASTQ output folder (always, since
        that is where ``.prepare-tmp/<ACC>`` staging folders, old-style ``<ACC>_temp``
        folders and, without a store, ``.sra-cache`` land; ``transient_bytes`` counts all
        three), any
        explicit ``--temp-folder`` or ``--sra-cache``, and, with a shared store, the
        store's own ``tmp`` folder (the default home for ``.sra-cache`` when downloading
        through a store)."""
        folders = {fastq_dir}
        temp_folder = setting_for(args, "temp_folder")
        if temp_folder:
            folders.add(Path(temp_folder))
        sra_cache = getattr(args, "sra_cache", None)
        if sra_cache:
            folders.add(Path(sra_cache).parent)
        if store is not None:
            folders.add(store.tmp)
        return folders

    def _warn_if_transient_bytes_large(
        self, args: argparse.Namespace, fastq_dir: Path, store: Optional[StorePaths]
    ) -> None:
        """Warn, naming each folder and its size, when kept transient artifacts add up."""
        sized = [(folder, transient_bytes(folder)) for folder in self._transient_folders(args, fastq_dir, store)]
        sized = [(folder, size) for folder, size in sized if size > 0]
        total = sum(size for _, size in sized)
        if total <= TRANSIENT_BYTES_WARN_THRESHOLD:
            return
        detail = ", ".join(f"{folder} ({size} bytes)" for folder, size in sorted(sized, key=lambda item: str(item[0])))
        self.logger.warning(
            "Kept .sra-cache archives, .prepare-tmp staging and <ACC>_temp build folders total %d bytes, "
            "over the "
            "%d byte warning threshold: %s",
            total,
            TRANSIENT_BYTES_WARN_THRESHOLD,
            detail,
        )

    def _store_options(self, args: argparse.Namespace, store: Optional[StorePaths], registry: Registry) -> dict:
        """The store-related keyword arguments for ``download_sra``; only ``lock_wait`` without a store.

        NCBI's spot count for an accession is read from the project's own metadata folder
        when it has one, since that is where ``download_metadata`` writes; the store's
        metadata folder is the fallback.
        """
        if store is None:
            # Without a store, --lock-wait still applies to the project's per-accession lock.
            return {"lock_wait": setting_for(args, "lock_wait")}
        return {
            "store": store,
            "link_mode": getattr(args, "link_mode", "auto"),
            "accept_partial": getattr(args, "accept_partial", False),
            "resume_partial": getattr(args, "resume_partial", True),
            "store_metadata": [project_root(registry) / "metadata", store.metadata],
            "lock_wait": setting_for(args, "lock_wait"),
        }

    def _download_batched(
        self,
        args: argparse.Namespace,
        fastq_dir: Path,
        store: Optional[StorePaths],
        project_registry: Registry,
        max_workers: int,
        timings: Optional[Timings] = None,
        inputs: Optional[Tuple[set, dict, set, dict]] = None,
        outcomes: Optional[RunOutcomes] = None,
    ) -> Optional[dict]:
        """Run ``download_sra`` with its outcomes queued in a registry batch; return its statistics.

        The batch is flushed when the download call returns or raises, a ``KeyboardInterrupt``
        included (SIGTERM and SIGHUP are turned into one for the duration, and a repeated signal
        during the flush is logged rather than raised), before the error propagates. A dry run
        records nothing, so its batch never writes. A final flush that fails with
        ``DataAccessError`` is retried once; None means the retry failed too, after the outcomes
        still queued have been logged. ``timings`` receives each download attempt's start and
        seconds; ``inputs`` is ``registry_inputs``'s result, computed here when not given.
        ``outcomes`` observes every result before it is recorded and receives the statistics,
        so the run summary has them even when the flush fails or the run is interrupted.
        """
        timings = timings if timings is not None else {}
        verify_downloads = getattr(args, "verify_downloads", True)
        excluded, expected_spots, truncated, run_sizes = inputs or registry_inputs(args, project_registry, store)
        on_result = None
        batch = registry_batch(args.registry)
        body_done = False
        # The run's stop token: BaseCommand.run's Termination, whose stop is set by the first
        # signal and targeted by its terminate_children call. Called without run (execute
        # directly), the context below makes a token of its own.
        run_term = getattr(args, "_termination", None)
        try:
            with _termination_raises_interrupt(run_term.stop if run_term is not None else None) as term, batch:
                if not args.dry_run:
                    on_result = self._result_recorder(args, fastq_dir, store, batch, timings, expected_spots)
                    on_result = outcomes.wrap(on_result) if outcomes is not None else on_result
                if outcomes is not None:
                    outcomes.reached_download = True
                download_stats = download_sra(
                    fastq_folder=args.fastq_folder,
                    accessions_file=args.accessions_file,
                    max_downloads=args.max_downloads,
                    dry_run=args.dry_run,
                    num_threads=args.num_threads,
                    max_workers=max_workers,
                    force=args.force,
                    max_retries=args.max_retries,
                    temp_folder=setting_for(args, "temp_folder"),
                    blacklist=args.blacklist,
                    blacklist_accessions=excluded,
                    on_result=on_result,
                    expected_spots=expected_spots if verify_downloads else None,
                    redownload_truncated=getattr(args, "redownload_truncated", False),
                    truncated_accessions=truncated,
                    sra_cache=getattr(args, "sra_cache", None),
                    use_prefetch=getattr(args, "use_prefetch", True),
                    keep_sra=getattr(args, "keep_sra", False),
                    compress=getattr(args, "compress", True),
                    stop=term.stop,
                    run_sizes=run_sizes,
                    min_free_gb=setting_for(args, "min_free_gb"),
                    timings=timings,
                    **self._store_options(args, store, project_registry),
                )
                if outcomes is not None:
                    outcomes.stats = download_stats
                body_done = True
        except DataAccessError as e:
            if not body_done:
                raise
            if not self._retry_final_flush(batch, e):
                return None
        return download_stats

    def _retry_final_flush(self, batch: RegistryBatch, error: DataAccessError) -> bool:
        """Retry a failed final flush once, after a short wait; log what stays unwritten if it fails again."""
        self.logger.warning(
            "Could not write %d queued download outcome(s) to %s (%s); retrying in %.0f s",
            len(batch),
            batch.path,
            error,
            FINAL_FLUSH_RETRY_SECONDS,
        )
        time.sleep(FINAL_FLUSH_RETRY_SECONDS)
        try:
            batch.flush()
        except DataAccessError as e:
            labels = batch.pending_labels()
            self.logger.error(
                "Could not write %d queued download outcome(s) to %s after a retry (%s); "
                "not recorded in the registry: %s",
                len(labels),
                batch.path,
                e,
                ", ".join(labels),
            )
            return False
        return True

    def execute(self, args: argparse.Namespace) -> int:
        outcomes, started, aborted = RunOutcomes(), datetime.now(timezone.utc), None
        try:
            code = self._run(args, outcomes)
        except KeyboardInterrupt:
            # download_sra has already cancelled pending downloads and stopped running tools.
            self.logger.error("Download interrupted by the user")
            code, aborted = 130, "interrupted"
        if not args.dry_run:
            self._write_run_summary(args, outcomes, started, code, aborted)
        return code

    def _run(self, args: argparse.Namespace, outcomes: RunOutcomes) -> int:
        try:
            if not args.dry_run:
                require_tools(["fasterq-dump"])

            max_workers = self._resolve_max_workers(args)
            fastq_dir = Path(args.fastq_folder)
            project_registry = load_registry(args.registry)
            store = self._resolve_store(args, project_registry)
            inputs = registry_inputs(args, project_registry, store)
            download_stats = self._download_batched(
                args, fastq_dir, store, project_registry, max_workers, outcomes.timings, inputs, outcomes
            )
            if download_stats is None:
                return 1
            if args.dry_run:
                self._log_dry_run_summary(args, download_stats)
            else:
                self._record_run_outcomes(args, download_stats, fastq_dir, store, inputs[1])
                self._warn_if_transient_bytes_large(args, fastq_dir, store)

            if not args.dry_run and download_stats.get("aborted"):
                self.logger.error(
                    "Download run aborted (%s); accessions attempted before the abort were " "still recorded above",
                    download_stats["aborted"],
                )
                return 1

            if not args.dry_run and download_stats["failed"] > 0:
                self._report_failed_downloads(args, download_stats)
                return self._failure_exit_code(download_stats)

            return 0

        except PrepareError as e:
            return self.fail(e, "Error downloading SRA data")
