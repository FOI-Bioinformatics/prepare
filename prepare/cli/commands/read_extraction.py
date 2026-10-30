"""CLI command for targeted read extraction before assembly."""

import argparse
import contextlib
import logging
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Set, Tuple

from prepare.cli.base import finite_float, BaseCommand, anchored_folder
from prepare.cli.commands.extraction_assembly import assemble_samples
from prepare.core.constants import DEFAULT_CONTAINMENT_THRESHOLD
from prepare.core.exceptions import PrepareError, TransientError
from prepare.core import settings
from prepare.core.settings import SETTINGS, setting_for
from prepare.data import registry_blocks as rb
from prepare.data import run_log
from prepare.data.read_extraction import (
    MINIMAP2_PRESETS,
    ExtractionResult,
    _record_matches,
    _sample_reads,
    extract_target_reads,
    resolve_index_path,
    selected_samples,
)
from prepare.data.registry import (
    Registry,
    load_registry,
    query,
    record_extraction,
    registry_transaction,
    resolve_project_path,
    scan_downloads,
)
from prepare.data.registry_timing import Stopwatch, set_extraction_timing
from prepare.data.sra import STORE_READY_STATES, count_fastq_reads
from prepare.data.sra_metadata import _resolved_sidecar_path
from prepare.store.layout import StorePaths
from prepare.store.link import dangling_links
from prepare.store.resolve import resolve_optional_store
from prepare.store.sidecar import Sidecar, read_sidecar
from prepare.store.stats import cached_stats
from prepare.store.usage import record_usage_safe
from prepare.utils.progress import DemoteInfo, ProgressReporter
from prepare.utils.tools import require_tools

# The per-sample INFO lines of data/read_extraction.py (held at a line ceiling), logged at the
# item level while a run reports progress; tests/test_cli_read_extraction.py checks they exist.
SAMPLE_LINE_TEMPLATES = (
    "%s: kept %d of %d mapped records (secondary/supplementary%s removed: %d)",
    "Extracted %d mapped records for %s -> %s",
    "%s already extracted against %s (%d mapped reads); use --force to redo",
    "%s: extraction against %s is in progress elsewhere; skipped",
    "%s: redoing extraction, %s",
)


def _store_copy_reason(sidecar: Sidecar) -> str:
    """Why a store copy whose sidecar state is not ready is skipped, for the skip line.

    A ``partial`` copy gives its read count against the NCBI spot count when both are
    recorded; a ``failed`` copy gives the verification error it recorded. Any other state,
    or a missing figure, is named as the state alone.
    """
    reads, spots = sidecar.reads_per_mate, sidecar.ncbi.get("spots")
    if sidecar.state == "partial" and reads is not None and spots is not None:
        return f"store copy partial ({reads} of {spots} spots)"
    if sidecar.state == "failed" and sidecar.error:
        return f"store copy failed verification: {sidecar.error}"
    return f"store copy {sidecar.state}"


def _non_negative_int(value: str) -> int:
    """argparse type for --min-mapq: a mapping quality is never negative.

    Without this, ``--min-mapq -5`` reaches samtools as ``-q -5``, where the security layer
    rejects ``-5`` as an unknown flag only after minimap2 has already run.
    """
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError(f"--min-mapq must be zero or a positive integer, got {value!r}")
    return parsed


class ExtractTargetReadsCommand(BaseCommand):
    """Map each sample's reads to a target genome and keep only the mapped reads."""

    @property
    def name(self) -> str:
        return "extract_target_reads"

    @property
    def help(self) -> str:
        return "Filter reads that map to a target genome for a small, targeted assembly"

    @property
    def group(self) -> str:
        return "Reads"

    def records_run(self, args: argparse.Namespace) -> bool:
        return not getattr(args, "dry_run", False)

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--parsed-containment",
            required=True,
            help="Parsed containment table (samples x genomes) from parse_containment",
        )
        parser.add_argument("--genome-id", required=True, help="Target genome column to extract against")
        parser.add_argument("--genome-fasta", required=True, help="FASTA file for the target genome")
        parser.add_argument("--fastq-folder", default="fastq", help="Root folder of per-accession FASTQ files")
        parser.add_argument("--output-folder", default="targeted", help="Root folder for the extracted reads")
        parser.add_argument(
            "--threshold",
            type=finite_float,
            default=DEFAULT_CONTAINMENT_THRESHOLD,
            help="Minimum containment for a sample to be included",
        )
        parser.add_argument(
            "--ignore-selection",
            action="store_true",
            help="Extract every sample at --threshold, not only the selected ones; excluded samples are always skipped",
        )
        parser.add_argument(
            "--preset", choices=sorted(MINIMAP2_PRESETS), default="sr", help="minimap2 preset for the read type"
        )
        parser.add_argument("--threads", type=int, default=4, help="Threads for minimap2 and samtools")
        parser.add_argument(
            "--min-mapq",
            type=_non_negative_int,
            default=0,
            help=(
                "Minimum mapping quality kept, in addition to dropping secondary/supplementary "
                "alignments (default: 0, keep every mapped record). Try 20 for a close relative of "
                "the target genome; a divergent strain can genuinely map with a low MAPQ, so raising "
                "this can discard real matches."
            ),
        )
        parser.add_argument(
            "--temp-folder",
            default=None,
            help=(
                "Where the intermediate SAM alignment(s) are written "
                "(default: PREPARE_TEMP_FOLDER, else alongside each sample's output)"
            ),
        )
        parser.add_argument(
            "--allow-truncated",
            action="store_true",
            help=(
                "Extract even for a sample whose registry download verdict is 'truncated', "
                "or whose store copy is recorded as partial or failed"
            ),
        )
        parser.add_argument(
            "--debug-keep-sam",
            action="store_true",
            help="Keep the intermediate SAM alignment(s) instead of removing them once the BAM exists "
            "(mapped records only)",
        )
        parser.add_argument(
            "--assemble", action="store_true", help="Assemble each sample's extracted reads with megahit"
        )
        parser.add_argument(
            "--assembly-threads",
            type=int,
            default=None,
            help="Threads for the megahit assembly (defaults to 1 on macOS, --threads elsewhere)",
        )
        parser.add_argument(
            "--assembly-memory",
            dest=SETTINGS["assembly_memory"].cli_dest,
            default=None,
            help=(
                "Memory for megahit (--memory): 'auto' gives 90%% of the memory limit detected for this job "
                "(cgroup or SLURM_MEM_PER_NODE) and leaves megahit's default when none is found; a size such "
                "as 32G or 32000M; or a fraction such as 0.5, which megahit applies to the whole node's "
                "memory, not to a cgroup or SLURM limit, so give a size under a scheduler "
                "(default: PREPARE_ASSEMBLY_MEMORY, config [runtime] assembly_memory, or auto)"
            ),
        )
        parser.add_argument("--min-contig-len", type=int, default=None, help="megahit minimum contig length")
        parser.add_argument(
            "--assembly-preset",
            choices=["default", "meta-sensitive", "meta-large"],
            default="meta-sensitive",
            help="megahit --presets value ('default' omits the flag)",
        )
        parser.add_argument(
            "--keep-intermediate",
            action="store_true",
            help="Keep megahit's intermediate_contigs/ folder instead of removing it after a successful assembly",
        )
        parser.add_argument(
            "--no-coverage",
            action="store_true",
            help="Skip mapping the extracted reads back onto the assembled contigs for coverage stats",
        )
        parser.add_argument(
            "--dry-run", action="store_true", help="List the qualifying samples without running any tool"
        )
        parser.add_argument(
            "--force",
            action="store_true",
            help="Redo extraction and assembly even when the registry says they are done",
        )
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")
        parser.add_argument(
            "--timeout",
            dest="timeout",
            type=finite_float,
            default=None,
            help=(
                "Seconds before minimap2, samtools or megahit is stopped; 0 (the default) means no "
                "limit (default: PREPARE_TIMEOUT, config [runtime] timeout, or 0)"
            ),
        )

    @staticmethod
    def _resolve_store(args: argparse.Namespace, registry: Registry) -> Optional[StorePaths]:
        """Resolve the shared data store (if any); returns None without one.

        Extraction reads the project's own ``fastq/`` folder either way, so a store that
        cannot be reached is a warning, not a reason to stop."""
        return resolve_optional_store(getattr(args, "data_root", None), rb.store_block(registry).root)

    def _warn_dangling_links(self, args: argparse.Namespace) -> None:
        """Log one WARNING naming every ``fastq/<ACC>`` symlink whose target is missing.

        A dangling link usually means the shared store is unmounted; the accession would
        otherwise fail silently with "No FASTQ files found" later, with no hint why.
        """
        dangling = dangling_links(args.fastq_folder)
        if dangling:
            self.logger.warning(
                "%d dangling fastq/<ACC> link(s) under %s (target missing; is the store mounted?): %s",
                len(dangling),
                args.fastq_folder,
                ", ".join(dangling),
            )

    @staticmethod
    def _mate_signature(reads: List[Path]) -> List[List[Any]]:
        """``[[name, size bytes, mtime], ...]`` for the mate files a count was taken from."""
        signature = []
        for path in reads:
            stat = path.stat()
            signature.append([path.name, stat.st_size, stat.st_mtime])
        return signature

    def _counts_from_store(self, accession: str, reads: List[Path]) -> Optional[Tuple[int, int]]:
        """The mate counts from the dataset's shared statistics record, when it has one.

        The record is size/mtime invalidated against the files on disk, so unlike the
        registry's own copy it cannot describe a download that has since been replaced.
        """
        acc_dir = reads[0].parent
        cached = cached_stats(acc_dir, _resolved_sidecar_path(acc_dir))
        reads_per_file = (cached or {}).get("reads_per_file") or {}
        first, second = (reads_per_file.get(path.name) for path in reads)
        if first is None or second is None:
            return None
        return int(first), int(second)

    def _cached_registry_counts(
        self, registry: Registry, accession: str, reads: List[Path]
    ) -> Optional[Tuple[int, int]]:
        """The registry's recorded mate pair, but only while it still describes these files.

        The pair is stored with a signature of the files it was counted from. A sample
        re-downloaded since (``--resume-partial`` completing a truncated pair) no longer
        matches, and the stale pair is ignored: an unequal pair would otherwise force
        single-end mapping of a now-complete paired run.
        """
        download = rb.download_block(registry, accession) or rb.DownloadBlock()
        cached = download.mate_reads
        if cached is None or len(cached) != 2:
            return None
        if download.mate_reads_signature != self._mate_signature(reads):
            self.logger.debug("%s: the recorded mate counts no longer match the files on disk", accession)
            return None
        return int(cached[0]), int(cached[1])

    def _mate_counts(
        self, args: argparse.Namespace, registry: Registry, selected: List[str]
    ) -> Dict[str, Tuple[int, int]]:
        """Accession -> (mate 1 reads, mate 2 reads) for every selected sample with both mate
        files present on disk.

        The counts come from the dataset's shared statistics record when it has one, then
        from the registry's own recorded pair while its file signature still matches, and
        otherwise from counting the files now, which is recorded with a signature so a later
        run does not recount them. Counting is best-effort: a file that cannot be read (e.g.
        corrupted, or not really gzip despite its name) is logged and skipped rather than
        stopping the extraction.
        """
        counts: Dict[str, Tuple[int, int]] = {}
        for accession in selected:
            reads = _sample_reads(Path(args.fastq_folder), accession)
            if len(reads) != 2:
                continue
            try:
                pair = self._counts_from_store(accession, reads) or self._cached_registry_counts(
                    registry, accession, reads
                )
                if pair is not None:
                    counts[accession] = pair
                    continue
                signature = self._mate_signature(reads)
                pair = (count_fastq_reads(reads[0]), count_fastq_reads(reads[1]))
            except Exception as e:
                self.logger.warning(
                    "Could not count reads in %s's mate files (%s); skipping the pre-count", accession, e
                )
                continue
            counts[accession] = pair
            with registry_transaction(args.registry) as reg:
                rb.set_mate_reads(reg, accession, [pair[0], pair[1]], signature)
        return counts

    def _samples_needing_mate_counts(
        self, selected: List[str], already_done: Dict[str, Any], truncated: Dict[str, Any], args: argparse.Namespace
    ) -> List[str]:
        """The selected samples that will really be mapped, so only their mates are counted.

        Counting reads both mate files of a sample, which is a streaming pass over every
        byte; doing it for a sample that is about to be skipped as already extracted or as an
        unusable download (``_unusable_downloads``) is pure cost.
        """
        genome_path = Path(args.genome_fasta)

        def will_be_skipped(accession: str) -> bool:
            record = already_done.get(accession)
            done = (
                record is not None
                and not args.force
                and _record_matches(record, genome_path, args.preset, args.threshold, min_mapq=args.min_mapq)
            )
            return done or (accession in truncated and not args.allow_truncated)

        return [accession for accession in selected if not will_be_skipped(accession)]

    @staticmethod
    def _available_accessions(args: argparse.Namespace, registry: Registry) -> Set[str]:
        """Accessions the registry or the FASTQ folder itself says are downloaded.

        Passed to ``extract_target_reads`` as ``available`` so a run over many screened-
        but-not-downloaded samples reports one summary line instead of a warning per
        missing sample.
        """
        recorded = set(query(registry, "downloaded"))
        return recorded | set(scan_downloads(Path(args.fastq_folder)))

    @staticmethod
    def _unusable_downloads(
        registry: Registry, fastq_folder: Path, accessions: Iterable[str]
    ) -> Dict[str, Dict[str, Any]]:
        """Accession -> why its download is skipped unless ``--allow-truncated`` is given.

        Holds every accession whose registry verdict is ``"truncated"`` (its verdict, as
        before) and every store copy among ``accessions`` (the selected samples) in
        ``fastq_folder`` whose sidecar records a state outside ``STORE_READY_STATES``, with a
        ``reason``. The sidecar is consulted whatever the registry verdict says: a record written
        as ``"unverified"`` may point at a store copy found short or failed since. Only the
        selected samples' sidecars are read, never a FASTQ file, so the cost follows the
        selection rather than the folder, a malformed sidecar of an unselected sample cannot stop
        the run, and a plain project download recorded as ``"unverified"`` is not skipped.
        """
        unusable: Dict[str, Dict[str, Any]] = {}
        for accession in registry.datasets:
            verdict = rb.download_verdict(registry, accession)
            if verdict is not None and verdict.verdict == "truncated":
                unusable[accession] = verdict.to_dict()
        for accession in dict.fromkeys(accessions):
            acc_dir = fastq_folder / accession
            if not acc_dir.is_dir():
                continue
            sidecar_file = _resolved_sidecar_path(acc_dir) or acc_dir / f"{acc_dir.name}.json"
            sidecar = read_sidecar(sidecar_file) if sidecar_file.is_file() else None
            if sidecar is not None and sidecar.state not in STORE_READY_STATES:
                unusable[acc_dir.name] = {**unusable.get(acc_dir.name, {}), "reason": _store_copy_reason(sidecar)}
        return unusable

    def _record_result(
        self,
        args: argparse.Namespace,
        accession: str,
        outcome: ExtractionResult,
        store: Optional[StorePaths] = None,
        timing: Optional[Tuple[str, float]] = None,
        recorded: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> None:
        """Checkpoint one extraction result, with its ``(started, seconds)`` timing when given.

        Skipped samples are already recorded, and their recorded timing is left alone; the run
        log still gets a row for one, with the values recorded before this run (``recorded``, the
        run's snapshot of the recorded extractions) and ``skipped: true``, so a later
        ``runs --diff`` does not read the sample as removed.
        """
        key = f"{accession}/{args.genome_id}"
        if outcome.skipped:
            before = (recorded or {}).get(accession) or {}
            row = {field: before.get(field) for field in ("mapped_reads", "breadth", "mean_depth")}
            run_log.note_rows(args, {key: {**row, "skipped": True}}, "extractions")
            return
        coverage = outcome.coverage or {}
        row = {"mapped_reads": outcome.mapped_records, "breadth": coverage.get("breadth")}
        row["mean_depth"] = coverage.get("mean_depth")
        run_log.note_rows(args, {key: row}, "extractions")
        index_dir = Path(args.output_folder) / ".index"
        with registry_transaction(args.registry) as reg:
            record_extraction(
                reg,
                accession,
                args.genome_id,
                outcome.files,
                outcome.mapped_records,
                outcome.unequal_mates,
                {
                    "genome_fasta": str(Path(args.genome_fasta)),
                    "preset": args.preset,
                    "threshold": args.threshold,
                    "filter_flags": "0x904",
                    "min_mapq": args.min_mapq,
                    "index": str(resolve_index_path(args.genome_fasta, args.preset, index_dir)),
                },
                mapped_total=outcome.mapped_total,
                coverage=outcome.coverage,
            )
            started, seconds = timing or (None, None)
            set_extraction_timing(reg, accession, args.genome_id, started, seconds)
        detail = f"{outcome.mapped_records} mapped reads"
        breadth = (outcome.coverage or {}).get("breadth")
        if breadth is not None:
            detail += f", breadth {breadth:.3f}"
        # After the registry lock is released: the catalogue has its own lock, and waiting for it
        # while holding the registry lock would hold up every other registry writer.
        record_usage_safe(store, reg, accession, args.genome_id, "extracted", detail=detail)

    def _record_result_and_check_stop(
        self,
        args: argparse.Namespace,
        accession: str,
        outcome: ExtractionResult,
        store: Optional[StorePaths],
        clock: Optional[Stopwatch] = None,
        progress: Optional[ProgressReporter] = None,
        recorded: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> None:
        """Checkpoint one sample, count it on ``progress``, then stop the run if a signal arrived.

        ``recorded`` is the run's snapshot of the recorded extractions (see ``_record_result``).

        ``clock`` times each sample from the moment the previous one was checkpointed (or the
        extraction started) until its result arrives here: the time ``extract_target_reads``
        spent on it, a module held at a frozen line ceiling, so it is measured from this side.
        The first sample also includes reading the containment table and selecting the samples,
        and the first that needs mapping the building of the minimap2 index.

        Checked here -- the boundary between one sample finishing and the next starting --
        rather than inside ``extract_target_reads``'s loop, a module held at a frozen line
        ceiling. Raising ``KeyboardInterrupt`` from an ``on_result`` callback is not swallowed
        by its caller (``_notify_result`` only catches ``Exception``), so it reaches ``execute``
        and then ``BaseCommand.run``, which logs, stops any running tool and returns 130. The
        sample just checkpointed above is recorded either way; the next one in
        ``extract_target_reads``'s ``samples`` list is never started.
        """
        timing = clock.lap() if clock is not None else None
        try:
            self._record_result(args, accession, outcome, store, timing, recorded)
            if progress is not None:
                if outcome.skipped:
                    state = "skipped"
                else:
                    state = f"{outcome.mapped_records} mapped reads" if outcome.files else "no reads mapped"
                self.logger.log(progress.item_level, "%s: %s", accession, state)
                progress.update(True)
        finally:
            # Checked even when the record failed (the failure is reported after the run), so
            # a stop that arrived during the write is not one sample late.
            if clock is not None:
                # The registry write above is bookkeeping, not part of the next sample's extraction.
                clock.restart()
            term = getattr(args, "_termination", None)
            if term is not None and term.stop.is_set():
                raise KeyboardInterrupt(f"extract_target_reads stopped after {accession}")

    def _fail_for_recording(self, failures: List[Tuple[str, PrepareError]], exit_code: int) -> int:
        """The exit code once the run's recording failures are reported: ``exit_code`` when there
        are none; otherwise the first transient failure's (4), or the first failure's (1) if that
        is above ``exit_code``."""
        if not failures:
            return exit_code
        named = ", ".join(f"{acc} ({(str(e).splitlines() or [''])[0]})" for acc, e in failures)
        self.logger.error("Recording failed for %d sample(s): %s", len(failures), named)
        first = next((e for _, e in failures if isinstance(e, TransientError)), failures[0][1])
        return max(exit_code, self.fail(first, "Error recording extraction results"))

    @staticmethod
    def _progress_reporter(selected: List[str], available: Set[str], already_done: Dict[str, Any]) -> ProgressReporter:
        """A reporter over the selected samples that have reads or a recorded extraction."""
        total = sum(1 for acc in selected if acc in available or acc in already_done)
        return ProgressReporter(
            "extract_target_reads", total, settings.active().progress_every, logger=logging.getLogger(__name__)
        )

    @staticmethod
    @contextlib.contextmanager
    def _sample_lines_at_item_level(progress: Optional[ProgressReporter]) -> Iterator[None]:
        """Log the data module's per-sample INFO lines at ``progress``'s item level while in the block."""
        if progress is None or progress.item_level == logging.INFO:
            yield
            return
        data_logger = logging.getLogger("prepare.data.read_extraction")
        demote = DemoteInfo(SAMPLE_LINE_TEMPLATES, progress.item_level)
        data_logger.addFilter(demote)
        try:
            yield
        finally:
            data_logger.removeFilter(demote)

    @staticmethod
    def _resolved_extraction_record(registry: Registry, accession: str, genome_id: str) -> Optional[Dict[str, Any]]:
        """The recorded extraction for one sample, with its ``genome_fasta`` and ``files``
        resolved against the project root.

        The registry stores these paths relative to its own location so the project can be
        moved; ``extract_target_reads`` and its ``_record_matches`` helper know nothing about
        the registry or its project root, so the paths must already be absolute (or otherwise
        directly usable) by the time they reach it.
        """
        block = rb.extraction_block(registry, accession, genome_id)
        if block is None:
            return None
        if block.genome_fasta is not None:
            block.genome_fasta = str(resolve_project_path(registry, block.genome_fasta))
        block.files = [str(resolve_project_path(registry, p)) for p in block.files]
        return block.to_dict()

    def _report_dry_run(self, args: argparse.Namespace, results: Dict[str, ExtractionResult]) -> None:
        """List the samples a real run would extract, and those it would skip."""
        would_skip = [acc for acc, r in results.items() if r.skipped]
        self.logger.info(
            "Dry run: %d sample(s) would be extracted for %s", len(results) - len(would_skip), args.genome_id
        )
        if would_skip:
            self.logger.info("  %d already extracted, would be skipped: %s", len(would_skip), ", ".join(would_skip))
        for accession, outcome in results.items():
            if not outcome.skipped:
                self.logger.info("  %s", accession)

    @staticmethod
    def _selected(args: argparse.Namespace, registry) -> list:
        """The accessions to extract: the table at the threshold, bounded by the registry (P2)."""
        return selected_samples(
            args.parsed_containment,
            args.genome_id,
            args.threshold,
            registry=registry,
            ignore_selection=getattr(args, "ignore_selection", False),
        )

    def _report_no_reads(self, args: argparse.Namespace, results: Dict[str, ExtractionResult]) -> None:
        """Say which of the three reasons left the run without a single mapped read."""
        registry = load_registry(args.registry)
        selected = self._selected(args, registry)
        if not selected:
            in_table = selected_samples(args.parsed_containment, args.genome_id, args.threshold)
            if in_table:
                self.logger.error(
                    "%d sample(s) meet containment >= %s for %s but all are excluded or outside the recorded "
                    "selection (%d excluded, %d selected); see blacklist, select_datasets or --ignore-selection",
                    len(in_table),
                    args.threshold,
                    args.genome_id,
                    len(query(registry, "excluded")),
                    len(query(registry, "selected")),
                )
            else:
                self.logger.error("No sample meets containment >= %s for %s", args.threshold, args.genome_id)
        elif not results:
            self.logger.error(
                "No FASTQ files found for the %d selected sample(s) under %s", len(selected), args.fastq_folder
            )
        else:
            self.logger.error("No reads mapped to %s in any sample; check the FASTQ files and --preset", args.genome_id)

    def _check_required_tools(self, args: argparse.Namespace) -> None:
        """Refuse to start when minimap2/samtools (and megahit with --assemble) are missing or too old.

        Checked once before any work starts, rather than surfacing as a raw subprocess
        error partway through extraction; skipped entirely for --dry-run, which never runs
        an external tool. Raises ``ConfigurationError`` (exit code 3) listing every problem.
        """
        if not args.dry_run:
            require_tools(["minimap2", "samtools"] + (["megahit"] if args.assemble else []))

    def execute(self, args: argparse.Namespace) -> int:
        try:
            self._check_required_tools(args)
            args.fastq_folder = str(anchored_folder(args.fastq_folder, getattr(args, "registry", None)))
            args.output_folder = str(anchored_folder(args.output_folder, getattr(args, "registry", None)))
            self._warn_dangling_links(args)
            registry = load_registry(args.registry)
            store = self._resolve_store(args, registry)
            already_done = {
                acc: rec
                for acc in registry.datasets
                if (rec := self._resolved_extraction_record(registry, acc, args.genome_id)) is not None
            }

            # Computed once: the unusable-download check reads only these samples' sidecars.
            selected = self._selected(args, registry)
            truncated_downloads = self._unusable_downloads(registry, Path(args.fastq_folder), selected)
            available = self._available_accessions(args, registry)

            mate_counts: Dict[str, Any] = {}
            progress: Optional[ProgressReporter] = None
            if not args.dry_run:
                to_count = self._samples_needing_mate_counts(selected, already_done, truncated_downloads, args)
                mate_counts = self._mate_counts(args, registry, to_count)
                progress = self._progress_reporter(selected, available, already_done)

            clock = Stopwatch()
            record_failures: List[Tuple[str, PrepareError]] = []
            with self._sample_lines_at_item_level(progress):
                results = extract_target_reads(
                    parsed_containment=args.parsed_containment,
                    genome_id=args.genome_id,
                    genome_fasta=args.genome_fasta,
                    fastq_folder=args.fastq_folder,
                    output_folder=args.output_folder,
                    threshold=args.threshold,
                    preset=args.preset,
                    threads=args.threads,
                    dry_run=args.dry_run,
                    force=args.force,
                    already_done=already_done,
                    on_result=lambda accession, outcome: self._record_result_and_check_stop(
                        args, accession, outcome, store, clock, progress, already_done
                    ),
                    min_mapq=args.min_mapq,
                    temp_folder=setting_for(args, "temp_folder"),
                    allow_truncated=args.allow_truncated,
                    mate_counts=mate_counts,
                    truncated_downloads=truncated_downloads,
                    keep_sam=args.debug_keep_sam,
                    samples=selected,
                    available=available,
                    record_failures=record_failures,
                )
            if progress is not None:
                progress.finish()

            if args.dry_run:
                self._report_dry_run(args, results)
                return 0

            # A zero-mapped sample can still have leftover files on disk from an earlier run;
            # they are not reads that mapped, and must never reach the assembler.
            with_reads = {acc: r.files for acc, r in results.items() if r.files and r.mapped_records > 0}
            self.logger.info("Extracted reads for %d of %d sample(s)", len(with_reads), len(results))
            skipped = sum(1 for r in results.values() if r.skipped)
            summary = {"genome_id": args.genome_id, "samples": len(results), "extracted": len(results) - skipped}
            run_log.note_run(args, summary={**summary, "with_reads": len(with_reads), "skipped": skipped})
            if not with_reads:
                self._report_no_reads(args, results)
                return self._fail_for_recording(record_failures, 1)

            if args.assemble:
                assembly = assemble_samples(self, args, with_reads, results, store)
                counts = {"assembled": len(assembly.assembled), "reused": len(assembly.reused)}
                run_log.note_run(args, summary={**counts, "assembly_failed": len(assembly.failed)})
                term = getattr(args, "_termination", None)
                if term is not None and term.stop.is_set():
                    raise KeyboardInterrupt("extract_target_reads assembly stopped")
                if assembly.failed:
                    # First line of each reason only; the per-sample ERROR line has the full text.
                    lines = {acc: (why.splitlines() or [""])[0] for acc, why in assembly.failed.items()}
                    failed = ", ".join(f"{acc} ({line})" for acc, line in lines.items())
                    self.logger.error("Assembly failed for %d sample(s): %s", len(assembly.failed), failed)
                    return self._fail_for_recording(record_failures, 1)
            return self._fail_for_recording(record_failures, 0)
        except PrepareError as e:
            return self.fail(e, "Error extracting target reads")
