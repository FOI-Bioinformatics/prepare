"""The ``sra_report`` command: one HTML report on a set of SRA datasets.

Replaces ``sra_dashboard`` and ``sra_compare`` (0.5.0). Each accession is profiled once, by the
same path as ``sra_profile`` (or read back from its saved profile JSON), and that one set of
profiles feeds both the quality section and, when ``--groups-file`` is given, the comparison
of groups with its statistical tests. The report is written as ``sra_report.html`` with its
figures in ``sra_report.json``, and a ``"report"`` analysis is recorded per accession.
"""

import argparse
import json
from collections import Counter
from functools import partial
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from prepare.cli.base import BaseCommand, anchored_folder, read_accessions_file, resolve_command_store
from prepare.cli.commands.sra_profile import add_sampling_arguments
from prepare.core.exceptions import DataAccessError, PrepareError, ValidationError
from prepare.core.optional import require
from prepare.data import run_log
from prepare.data.file_io import write_text_atomic
from prepare.data.registry import load_registry, record_analysis
from prepare.data.registry_batch import registry_batch
from prepare.sra.analytics import AnomalyReport, ComparativeAnalysis, QualityProfile, json_safe
from prepare.sra.dataset_stats import DATASET_READ_ERRORS
from prepare.sra.profiles import load_quality_profiles, profile_accession
from prepare.sra.reporting import SRAReportGenerator
from prepare.store.usage import record_usage_many
from prepare.utils.browser import open_in_browser

ANALYSIS_NAME = "report"
REPORT_HTML = "sra_report.html"
REPORT_JSON = "sra_report.json"

_GROUPS_EXAMPLE = json.dumps({"Group_A": ["SRR123456", "SRR123457"], "Group_B": ["SRR789012", "SRR789013"]})


def load_groups(path: str) -> Dict[str, List[str]]:
    """Accession groups from a JSON file mapping each group name to a list of accessions.

    Raises ``DataAccessError`` for a missing or unreadable file, and ``ValidationError`` for
    invalid JSON or any other shape, with an example of the expected format.
    """
    try:
        groups = json.loads(Path(path).read_text())
    except OSError as e:
        raise DataAccessError(f"Cannot read groups file {path}: {e}") from e
    except json.JSONDecodeError as e:
        raise ValidationError(f"Invalid JSON in groups file {path}: {e}; expected e.g. {_GROUPS_EXAMPLE}") from e
    valid = isinstance(groups, dict) and all(
        isinstance(accs, list) and all(isinstance(a, str) for a in accs) for accs in groups.values()
    )
    if not valid or not groups:
        raise ValidationError(f"Groups file {path} must map group names to accession lists, e.g. {_GROUPS_EXAMPLE}")
    return groups


def _analysis_summaries(
    profiles: Dict[str, QualityProfile], anomalies: AnomalyReport, groups: Optional[Dict[str, List[str]]]
) -> Dict[str, Dict[str, Any]]:
    """Accession -> the figures recorded for it: the registry's "report" analysis and its run-log row."""
    group_of = {acc: name for name, accs in (groups or {}).items() for acc in accs}
    return {
        accession: {
            "quality_grade": profile.quality_grade,
            "gc_percent": profile.gc_percent,
            "total_reads": profile.total_reads,
            "total_bases": profile.total_bases,
            "group": group_of.get(accession),
            "anomalous": accession in anomalies.anomalous_datasets,
        }
        for accession, profile in profiles.items()
    }


class SRAReportCommand(BaseCommand):
    """Write one HTML report on the quality of SRA datasets and, with groups, their comparison."""

    @property
    def name(self) -> str:
        """Command name."""
        return "sra_report"

    @property
    def help(self) -> str:
        """One-line help."""
        return "One HTML report on SRA dataset quality and, with --groups-file, a comparison of groups"

    @property
    def group(self) -> str:
        """Pipeline step."""
        return "Reads"

    def records_run(self, args: argparse.Namespace) -> bool:
        """Every run is added to the project's run log, one row per reported accession."""
        return True

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add the command's options."""
        parser.add_argument("--fastq-folder", default="fastq", help="Folder holding one folder per accession")
        parser.add_argument("--accessions-file", default=None, help="Report on the accessions listed here")
        parser.add_argument(
            "--groups-file",
            default=None,
            help="JSON file mapping group names to accession lists; adds the comparison of groups and "
            "the statistical tests (t-test for two groups, ANOVA for more). Its accessions are reported on too",
        )
        parser.add_argument(
            "--quality-profiles",
            default=None,
            help="Folder of profile JSONs written by sra_profile --output-dir; an accession found there is "
            "reused instead of being profiled again. With no other accession source, every profile there is used",
        )
        parser.add_argument("--output-dir", default="sra_reports", help="Folder for sra_report.html and .json")
        parser.add_argument("--title", default="SRA Report", help="Report title")
        parser.add_argument("--no-open", action="store_true", help="Do not open the report in a browser")
        parser.add_argument(
            "--no-report",
            action="store_true",
            help="Print the results and write sra_report.json only, without the HTML (needs no plotting packages)",
        )
        add_sampling_arguments(parser)
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")

    def _saved_profiles(self, path: Optional[str]) -> Dict[str, QualityProfile]:
        """Profiles saved by sra_profile under ``path``; a path that is not a folder is logged."""
        if not path:
            return {}
        folder = Path(path)
        if folder.is_dir():
            profiles = load_quality_profiles(folder)
            if profiles:
                self.emit(f"Reusing {len(profiles)} saved quality profile(s) from {folder}")
            return profiles
        if folder.exists():
            self.logger.warning("Quality profiles path is not a directory: %s", folder)
        else:
            self.logger.warning("Quality profiles directory not found: %s", folder)
        return {}

    @staticmethod
    def _resolve_accessions(
        args: argparse.Namespace, groups: Optional[Dict[str, List[str]]], saved: Dict[str, QualityProfile]
    ) -> List[str]:
        """Accessions from --accessions-file and the groups file, else every saved profile."""
        named = read_accessions_file(args.accessions_file) if args.accessions_file else []
        named += [acc for accs in (groups or {}).values() for acc in accs]
        if named:
            return list(dict.fromkeys(named))
        if saved:
            return sorted(saved)
        raise ValidationError(
            "Give --accessions-file, --groups-file or a --quality-profiles folder with saved profiles"
        )

    def _collect_profiles(
        self,
        reporter: SRAReportGenerator,
        args: argparse.Namespace,
        accessions: List[str],
        saved: Dict[str, QualityProfile],
    ) -> Tuple[Dict[str, QualityProfile], List[str]]:
        """One profile per accession, saved or computed once; returns (profiles, accessions left out)."""
        profiles: Dict[str, QualityProfile] = {}
        failed: List[str] = []
        for accession in accessions:
            if accession in saved:
                profiles[accession] = saved[accession]
                continue
            try:
                dataset = profile_accession(reporter.analyzer, accession, args.sample_size, args.sampler)
            except DATASET_READ_ERRORS as e:
                self.logger.warning("Failed to profile %s: %s", accession, e)
                dataset = None
            if dataset is None:
                failed.append(accession)
            else:
                profiles[accession] = dataset.profile
        return profiles, failed

    def _print_quality(self, profiles: Dict[str, QualityProfile], anomalies: AnomalyReport) -> None:
        """Print dataset count, grade distribution and anomalous datasets."""
        grades = Counter(p.quality_grade for p in profiles.values())
        self.emit(f"\nQuality of {len(profiles)} dataset(s):")
        self.emit("  Grades: " + ", ".join(f"{grade} {count}" for grade, count in sorted(grades.items())))
        mean_gc = sum(p.gc_percent for p in profiles.values()) / len(profiles)
        self.emit(f"  Mean GC content: {mean_gc:.1f}%")
        if anomalies.anomalous_datasets:
            self.emit(f"  Datasets flagged for review: {', '.join(anomalies.anomalous_datasets)}")

    def _print_comparison(self, comparison: ComparativeAnalysis) -> None:
        """Print per-group means and the statistical tests, marking p < 0.05 with '*'."""
        self.emit("\nComparison of groups:")
        for group_name, col_stats in comparison.summary_statistics.items():
            self.emit(f"  {group_name}: {len(comparison.dataset_groups.get(group_name, []))} dataset(s)")
            if "gc_percent" in col_stats:
                self.emit(f"    Mean GC content: {col_stats['gc_percent']['mean']:.1f}%")
            if "avg_read_length" in col_stats:
                self.emit(f"    Mean read length: {col_stats['avg_read_length']['mean']:.1f}")
            if "total_reads" in col_stats:
                self.emit(f"    Mean total reads (mates counted): {col_stats['total_reads']['mean']:,.0f}")
        for test_name, result in comparison.statistical_tests.items():
            p_value = result.get("p_value", 1.0)
            self.emit(f"  {test_name}: {result.get('test', '')} p={p_value:.4f} {'*' if p_value < 0.05 else ''}")

    @staticmethod
    def _comparison_payload(comparison: Optional[ComparativeAnalysis]) -> Optional[Dict[str, Any]]:
        """The comparison as a JSON-ready dict, or None without groups."""
        if comparison is None:
            return None
        tests = comparison.statistical_tests
        return {
            "groups": comparison.dataset_groups,
            "summary_statistics": comparison.summary_statistics,
            "statistical_tests": tests,
            "significant_differences": [col for col, test in tests.items() if test.get("significant")],
            "outlier_datasets": comparison.outlier_datasets,
            "recommendations": comparison.recommendations,
        }

    def _write_json(
        self,
        output_dir: Path,
        profiles: Dict[str, QualityProfile],
        failed: List[str],
        anomalies: AnomalyReport,
        comparison: Optional[ComparativeAnalysis],
    ) -> Path:
        """Write sra_report.json; numpy scalars from the statistical tests go through json_safe."""
        payload = {
            "accessions": sorted(profiles),
            "failed_accessions": failed,
            "quality_grades": {acc: p.quality_grade for acc, p in profiles.items()},
            "gc_percent": {acc: p.gc_percent for acc, p in profiles.items()},
            "anomalous_datasets": anomalies.anomalous_datasets,
            "anomaly_explanations": anomalies.explanations,
            "comparison": self._comparison_payload(comparison),
        }
        path = output_dir / REPORT_JSON
        write_text_atomic(path, json.dumps(json_safe(payload), indent=2))
        return path

    def _record(
        self,
        args: argparse.Namespace,
        output: Path,
        profiles: Dict[str, QualityProfile],
        anomalies: AnomalyReport,
        groups: Optional[Dict[str, List[str]]],
    ) -> None:
        """Record a "report" analysis and an "analysed" store usage for every reported accession."""
        # A snapshot, only to find the store; the records go through a batch that loads the
        # registry inside the lock, and the catalogue is written once that lock is released.
        store = resolve_command_store(args, load_registry(args.registry))
        with registry_batch(args.registry, flush_every=None, flush_seconds=None) as batch:
            for accession, summary in _analysis_summaries(profiles, anomalies, groups).items():
                batch.apply(
                    partial(
                        record_analysis, accession=accession, analysis=ANALYSIS_NAME, output=output, summary=summary
                    ),
                    accession,
                )
        if batch.registry is not None:
            record_usage_many(store, batch.registry, [(acc, "", "analysed", ANALYSIS_NAME) for acc in profiles])

    def _open(self, html: Path, no_open: bool) -> None:
        """Print the report's link and open it in a browser unless told not to."""
        self.emit(f"\nReport written: {html.resolve().as_uri()}")
        if no_open:
            return
        if open_in_browser(html):
            self.emit("Report opened in browser")
        else:
            self.emit("Could not open a browser automatically; open the link above manually.")

    def _run(self, args: argparse.Namespace) -> int:
        if not args.no_report:
            # Fail on a missing extra before any profiling or output, not after the JSON is written.
            require("plotly.graph_objects", "interactive", "The sra_report HTML report")
            require("jinja2", "interactive", "The sra_report HTML report")
        if args.groups_file:
            # Needed for the comparison's statistical tests; checked here (whether or not
            # --no-report is given) so a missing scipy fails at once rather than after profiling.
            require("scipy.stats", "analysis", "Comparing dataset groups (--groups-file)")
        groups = load_groups(args.groups_file) if args.groups_file else None
        saved = self._saved_profiles(args.quality_profiles)
        accessions = self._resolve_accessions(args, groups, saved)
        output_dir = Path(args.output_dir)
        reporter = SRAReportGenerator(output_dir=output_dir, fastq_dir=args.fastq_folder)

        profiles, failed = self._collect_profiles(reporter, args, accessions, saved)
        if not profiles:
            self.logger.error(
                "No dataset could be profiled: none of %d accession(s) has readable FASTQ files under %s",
                len(accessions),
                args.fastq_folder,
            )
            return 1
        if failed:
            self.logger.warning("Left out of the report (no readable FASTQ files): %s", ", ".join(failed))

        anomalies = reporter.analyzer.detect_dataset_anomalies(list(profiles), profiles=profiles)
        self._print_quality(profiles, anomalies)
        comparison = None
        if groups is not None:
            profiled_groups = {name: [a for a in accs if a in profiles] for name, accs in groups.items()}
            comparison = reporter.analyzer.compare_datasets(profiled_groups, profiles=profiles)
            self._print_comparison(comparison)

        summary = {
            "accessions": len(accessions),
            "reported": len(profiles),
            "failed": len(failed),
            "anomalous": len(anomalies.anomalous_datasets),
            "groups": len(groups) if groups is not None else 0,
            "html": not args.no_report,
        }
        run_log.note_run(args, summary=summary)
        run_log.note_rows(args, _analysis_summaries(profiles, anomalies, groups), "analyses", ANALYSIS_NAME)
        output = self._write_json(output_dir, profiles, failed, anomalies, comparison)
        self.emit(f"\nReport figures saved to: {output}")
        if not args.no_report:
            output = reporter.generate_report(
                profiles, title=args.title, comparison=comparison, anomalies=anomalies, filename=REPORT_HTML
            )
            self._open(output, args.no_open)
        self._record(args, output, profiles, anomalies, groups)
        return 0

    def execute(self, args: argparse.Namespace) -> int:
        """Write the report; 1 when no dataset could be profiled."""
        try:
            args.fastq_folder = str(anchored_folder(args.fastq_folder, getattr(args, "registry", None)))
            return self._run(args)
        except (PrepareError, OSError) as e:
            return self.fail(e, "Report failed")
        except Exception as e:  # noqa: B902 - top-level catch: keep the traceback, return 1
            self.logger.exception("Report failed: %s", e)
            return 1
