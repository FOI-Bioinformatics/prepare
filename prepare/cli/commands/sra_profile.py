"""The ``sra_profile`` command: statistics and a quality profile for each downloaded dataset.

Replaces ``sra_stats`` and ``sra_profile_quality`` (0.5.0). Each accession is profiled once:
read and base totals, mean read length and GC come from the dataset's shared statistics
record (``prepare.store.stats.compute_dataset_stats``, cached in the store sidecar), the
per-read quality, complexity and contamination figures from a sample of every mate file.
The command writes one statistics table, one profile JSON per accession, and records a
``"profile"`` analysis per accession in the project registry. GC is in percent throughout.
"""

import argparse
import json
from functools import partial
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from prepare.cli.base import BaseCommand, accessions_from_args, anchored_folder, resolve_command_store
from prepare.core.exceptions import PrepareError, ValidationError
from prepare.data import run_log
from prepare.data.file_io import visible_files, write_text_atomic
from prepare.data.registry import load_registry, record_analysis
from prepare.data.registry_batch import registry_batch
from prepare.data.sra import is_transient_folder
from prepare.data.sra_metadata import generate_statistics_report
from prepare.sra.analytics import QualityProfile, SRADatasetAnalyzer
from prepare.sra.dataset_stats import DATASET_READ_ERRORS
from prepare.sra.profiles import ProfiledDataset, profile_accession, statistics_row, write_profile_json
from prepare.store.stats import DEFAULT_SAMPLE_SIZE
from prepare.store.usage import record_usage_many

ANALYSIS_NAME = "profile"

# Thresholds above which a profile line is flagged for the reader.
_HIGH_N_CONTENT = 0.01
_HIGH_DUPLICATION = 0.20
_HIGH_ADAPTER = 0.05


def positive_int(value: str) -> int:
    """argparse type for --sample-size: rejects zero and negative values with a clear message."""
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError(f"--sample-size must be a positive integer, got {value!r}")
    return parsed


def add_sampling_arguments(parser: argparse.ArgumentParser) -> None:
    """The --sample-size and --sampler options shared by sra_profile and sra_report."""
    parser.add_argument(
        "--sample-size",
        type=positive_int,
        default=DEFAULT_SAMPLE_SIZE,
        help="Reads sampled per dataset for the per-read figures. GC content is computed from a sample of "
        "the first mate file; per-read quality, complexity and contamination from a sample of all mates. "
        "Read totals are exact counts from the dataset's statistics record",
    )
    parser.add_argument(
        "--sampler",
        choices=["uniform", "head"],
        default="uniform",
        help="'uniform' samples reads across the whole of each file; 'head' takes the first reads of each file",
    )


def _flag_lines(profile: QualityProfile) -> List[str]:
    """Warnings for a profile's high N content, duplicate rate or adapter content."""
    lines = []
    if profile.n_content > _HIGH_N_CONTENT:
        lines.append(f"WARNING: high N content: {profile.n_content:.1%}")
    if profile.duplication_rate is not None and profile.duplication_rate > _HIGH_DUPLICATION:
        lines.append(f"WARNING: high duplicate rate: {profile.duplication_rate:.1%}")
    adapter = profile.contamination_indicators.get("adapter_contamination", 0)
    if adapter > _HIGH_ADAPTER:
        lines.append(f"WARNING: adapter contamination: {adapter:.1%}")
    return lines


def _analysis_summary(profile: QualityProfile) -> Dict[str, Any]:
    """The figures recorded for one accession: the registry's "profile" analysis and its run-log row."""
    return {
        "total_reads": profile.total_reads,
        "total_bases": profile.total_bases,
        "reads_sampled": profile.reads_sampled,
        "sampled": profile.sampled,
        "gc_percent": profile.gc_percent,
        "avg_read_length": profile.avg_read_length,
        "quality_grade": profile.quality_grade,
    }


class SRAProfileCommand(BaseCommand):
    """Profile downloaded SRA datasets: one statistics table, one profile JSON per accession."""

    @property
    def name(self) -> str:
        """Command name."""
        return "sra_profile"

    @property
    def help(self) -> str:
        """One-line help."""
        return "Statistics and a quality profile for each downloaded SRA dataset"

    @property
    def group(self) -> str:
        """Pipeline step."""
        return "Reads"

    def records_run(self, args: argparse.Namespace) -> bool:
        """Every run is added to the project's run log, one row per profiled accession."""
        return True

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add the command's options."""
        parser.add_argument("--fastq-folder", default="fastq", help="Folder holding one folder per accession")
        parser.add_argument("--accessions-file", default=None, help="Profile the accessions listed here, one per line")
        parser.add_argument(
            "--accession",
            action="append",
            default=None,
            help="Profile this accession (repeatable). With neither this nor --accessions-file, every "
            "accession folder in --fastq-folder is profiled",
        )
        add_sampling_arguments(parser)
        parser.add_argument("--output-report", default="sra_statistics.csv", help="Statistics table (CSV)")
        parser.add_argument(
            "--output-dir",
            default="sra_quality_profiles",
            help="Folder for the per-accession profile JSONs and quality_summary.json; "
            "sra_report --quality-profiles reads the JSONs back",
        )
        parser.add_argument(
            "--detailed-reports",
            action="store_true",
            help="Print the path of each per-accession profile JSON as it is written (they are always written)",
        )
        parser.add_argument(
            "--summary-only", action="store_true", help="Print only the summary, not each accession's profile"
        )
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")

    def _resolve_accessions(self, args: argparse.Namespace) -> List[str]:
        """The accessions named on the command line, else every accession folder in --fastq-folder."""
        named = accessions_from_args(args.accessions_file, args.accession)
        if named:
            return named
        folder = Path(args.fastq_folder)
        found = [d.name for d in visible_files(folder, dirs=True) if not is_transient_folder(d.name)]
        if not found:
            raise ValidationError(f"No accession folders in {folder}; give --accession or --accessions-file")
        return found

    def _print_profile(self, profile: QualityProfile) -> None:
        """Print one accession's profile."""
        self.emit(f"\nQuality profile: {profile.accession}")
        self.emit("=" * 50)
        self.emit(f"Total reads (mates counted): {profile.total_reads:,} (sampled {profile.reads_sampled:,})")
        self.emit(f"Total bases: {profile.total_bases:,}")
        self.emit(f"Average read length: {profile.avg_read_length:.1f}")
        self.emit(f"GC content: {profile.gc_percent:.1f}%")
        self.emit(f"Mean base quality: {profile.mean_quality:.1f}")
        self.emit(f"Quality grade: {profile.quality_grade}")
        self.emit(f"Sequence complexity: {profile.complexity_score:.3f}")
        for line in _flag_lines(profile):
            self.emit(line)

    def _profile_one(
        self, analyzer: SRADatasetAnalyzer, args: argparse.Namespace, accession: str
    ) -> Optional[ProfiledDataset]:
        """Profile one accession; None (logged) when it has no FASTQ files or they cannot be read."""
        try:
            dataset = profile_accession(analyzer, accession, args.sample_size, args.sampler)
        except DATASET_READ_ERRORS as e:
            self.logger.warning("Failed to profile %s: %s", accession, e)
            return None
        if dataset is None:
            self.logger.warning("No FASTQ files found for %s under %s", accession, args.fastq_folder)
        return dataset

    def _profile_all(
        self, analyzer: SRADatasetAnalyzer, args: argparse.Namespace, accessions: List[str], output_dir: Path
    ) -> Tuple[List[ProfiledDataset], List[str]]:
        """Profile every accession, write its JSON, and return (profiled datasets, failed accessions).

        Checked before each accession: a signal (``args._termination.stop``) stops the loop
        there, leaving every later accession unprofiled and its JSON unwritten. ``_run`` records
        the accessions profiled so far and then raises, so ``BaseCommand.run`` returns 130
        without costing the run what it already has.
        """
        profiled: List[ProfiledDataset] = []
        failed: List[str] = []
        term = getattr(args, "_termination", None)
        for i, accession in enumerate(accessions, 1):
            if term is not None and term.stop.is_set():
                self.logger.warning("Stopping before %s: interrupted", accession)
                break
            self.emit(f"[{i}/{len(accessions)}] Profiling {accession}...")
            dataset = self._profile_one(analyzer, args, accession)
            if dataset is None:
                failed.append(accession)
                continue
            path = write_profile_json(dataset.profile, output_dir)
            if args.detailed_reports:
                self.emit(f"  Profile saved: {path}")
            if not args.summary_only:
                self._print_profile(dataset.profile)
            profiled.append(dataset)
        return profiled, failed

    @staticmethod
    def _summary_stats(profiles: List[QualityProfile]) -> Optional[Dict[str, Any]]:
        """Read and base totals and mean GC percent across ``profiles``, or None when empty."""
        if not profiles:
            return None
        return {
            "total_reads": sum(p.total_reads for p in profiles),
            "total_bases": sum(p.total_bases for p in profiles),
            "avg_gc_percent": sum(p.gc_percent for p in profiles) / len(profiles),
        }

    def _write_quality_summary(self, output_dir: Path, profiles: List[QualityProfile], failed: List[str]) -> Path:
        """Write quality_summary.json (counts, failed accessions, totals) and return its path."""
        path = output_dir / "quality_summary.json"
        summary = {
            "total_analyzed": len(profiles),
            "total_failed": len(failed),
            "failed_accessions": failed,
            "summary_stats": self._summary_stats(profiles),
        }
        write_text_atomic(path, json.dumps(summary, indent=2))
        return path

    def _record(self, args: argparse.Namespace, profiled: List[ProfiledDataset], output_dir: Path) -> None:
        """Record a "profile" analysis and an "analysed" store usage for every profiled accession."""
        # A snapshot, only to find the store; the records go through a batch that loads the
        # registry inside the lock, and the catalogue is written once that lock is released.
        store = resolve_command_store(args, load_registry(args.registry))
        with registry_batch(args.registry, flush_every=None, flush_seconds=None) as batch:
            for dataset in profiled:
                profile = dataset.profile
                json_path = output_dir / f"{profile.accession}_quality_profile.json"
                batch.apply(
                    partial(
                        record_analysis,
                        accession=profile.accession,
                        analysis=ANALYSIS_NAME,
                        output=json_path,
                        summary=_analysis_summary(profile),
                    ),
                    profile.accession,
                )
        if batch.registry is not None:
            rows = [(d.profile.accession, "", "analysed", ANALYSIS_NAME) for d in profiled]
            record_usage_many(store, batch.registry, rows)

    def _note_run(
        self, args: argparse.Namespace, accessions: List[str], profiles: List[QualityProfile], failed: List[str]
    ) -> None:
        """Counts and totals for the run-log summary; one row per profiled accession for its detail."""
        summary = {"accessions": len(accessions), "profiled": len(profiles), "failed": len(failed)}
        summary.update(self._summary_stats(profiles) or {})
        run_log.note_run(args, summary=summary)
        run_log.note_rows(args, {p.accession: _analysis_summary(p) for p in profiles}, "analyses", ANALYSIS_NAME)

    def _run(self, args: argparse.Namespace) -> int:
        folder = Path(args.fastq_folder)
        if not folder.is_dir():
            self.logger.error("FASTQ folder %s does not exist", folder)
            return 1
        accessions = self._resolve_accessions(args)
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        profiled, failed = self._profile_all(SRADatasetAnalyzer(fastq_dir=folder), args, accessions, output_dir)
        rows = [statistics_row(d.profile, d.files, d.stats) for d in profiled]
        for line in generate_statistics_report(rows, args.output_report):
            self.emit(line)
        profiles = [d.profile for d in profiled]
        self._note_run(args, accessions, profiles, failed)
        flagged = [p.accession for p in profiles if _flag_lines(p)]
        if flagged:
            self.emit(f"\n{len(flagged)} dataset(s) with quality warnings: {', '.join(flagged)}")
        summary_path = self._write_quality_summary(output_dir, profiles, failed)
        if rows:
            self.emit(f"\nStatistics report saved to: {args.output_report}")
            self._record(args, profiled, output_dir)
        self.emit(f"Profiles and summary saved to: {summary_path.parent}")

        term = getattr(args, "_termination", None)
        if term is not None and term.stop.is_set():
            raise KeyboardInterrupt("sra_profile stopped")
        if failed:
            self.logger.warning("%d of %d accession(s) could not be profiled", len(failed), len(accessions))
            return 1
        return 0

    def execute(self, args: argparse.Namespace) -> int:
        """Profile the datasets; 1 when any accession could not be profiled."""
        try:
            args.fastq_folder = str(anchored_folder(args.fastq_folder, getattr(args, "registry", None)))
            return self._run(args)
        except (PrepareError, OSError) as e:
            return self.fail(e, "Profiling failed")
        except Exception as e:  # noqa: B902 - top-level catch: keep the traceback, return 1
            self.logger.exception("Profiling failed: %s", e)
            return 1
