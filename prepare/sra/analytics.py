"""
Advanced SRA Dataset Analytics and Statistical Reporting.

This module provides comprehensive analysis capabilities for SRA datasets including:
- Quality profiling and contamination detection
- Comparative analysis across multiple datasets
- Technology and temporal trend analysis
- Advanced statistical reporting with visualizations
- Processing parameter recommendations
"""

import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union, Any

import numpy as np
import pandas as pd

from prepare.core.exceptions import DataAccessError
from prepare.data.sra_metadata import SRADatasetInfo
from prepare.processing.statistics import compare_group_means, holm_adjusted
from prepare.sra.dataset_stats import accession_fastq_files, load_dataset_stats
from prepare.sra.quality import SequenceQualityAnalyzer, _gc_histogram  # noqa: F401 - re-exported

logger = logging.getLogger(__name__)

# Per-dataset columns the comparative analysis summarises, tests and screens for outliers.
NUMERIC_COLUMNS = ["avg_read_length", "gc_percent", "total_reads", "complexity_score"]


def json_safe(value: Any) -> Any:
    """Return ``value`` with numpy scalars, sets and non-finite floats replaced by JSON-safe
    Python values, recursively.

    ``json.dump``/``json.dumps`` reject a numpy ``bool_``/integer/floating scalar (e.g. the
    ``p_value < 0.05`` comparison in ``_perform_statistical_tests`` produces a numpy ``bool_``,
    not a Python ``bool``) and a bare ``set``, and turn a non-finite float (``nan``/``inf``)
    into invalid JSON tokens rather than raising. Call this on a payload before dumping it.
    """
    if isinstance(value, (pd.Series, pd.DataFrame)):
        return json_safe(value.to_dict())
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return None if not np.isfinite(value) else float(value)
    return value


@dataclass
class QualityProfile:
    """Comprehensive quality profile for SRA dataset."""

    accession: str
    total_reads: int
    total_bases: int
    avg_read_length: float
    read_length_distribution: Dict[str, int]  # length_range -> count
    gc_percent: float  # 0-100
    gc_histogram: Dict[str, int]  # 5-percent-wide GC buckets, e.g. {"40-45": 12, ...}
    quality_distribution: Dict[str, float]  # quality_range -> percentage
    n_content: float
    contamination_indicators: Dict[str, float]
    complexity_score: float
    duplication_rate: Optional[float]
    technology_confidence: float
    quality_grade: str  # 'excellent', 'good', 'fair', 'poor'
    recommendations: List[str]
    # Records actually read for the per-read metrics above. ``total_reads``/``total_bases``
    # are the dataset totals from the shared statistics record when it is available;
    # ``sampled`` is True when they are the sample's own figures instead, i.e. a lower bound.
    reads_sampled: int = 0
    sampled: bool = False
    # Mean per-base Phred quality over the sampled reads.
    mean_quality: float = 0.0


@dataclass
class ComparativeAnalysis:
    """Results from comparative analysis across datasets."""

    dataset_groups: Dict[str, List[str]]
    summary_statistics: Dict[str, Dict[str, Any]]
    statistical_tests: Dict[str, Dict[str, Any]]
    outlier_datasets: List[str]
    clustering_results: Optional[Dict[str, Any]]
    batch_effects: Dict[str, float]
    recommendations: List[str]
    visualization_data: Dict[str, Any]


@dataclass
class AnomalyReport:
    """Report of detected anomalies in datasets."""

    anomalous_datasets: List[str]
    anomaly_types: Dict[str, List[str]]  # anomaly_type -> affected_datasets
    severity_scores: Dict[str, float]  # dataset -> severity (0-1)
    explanations: Dict[str, str]  # dataset -> explanation
    recommended_actions: Dict[str, List[str]]


@dataclass
class ProcessingRecommendations:
    """Recommended processing parameters based on dataset characteristics."""

    accession: str
    recommended_pipeline: str
    quality_trimming: Dict[str, Any]
    adapter_removal: Dict[str, Any]
    contamination_filtering: Dict[str, Any]
    assembly_parameters: Dict[str, Any]
    expected_coverage: Optional[float]
    computational_requirements: Dict[str, Any]
    estimated_processing_time: str


def _dataset_totals(
    dataset_stats: Optional[Dict[str, Any]], reads_sampled: int, mean_read_length: float, sample_gc: float
) -> Tuple[int, int, float, bool]:
    """``(total_reads, total_bases, gc_percent, sampled)`` for a profile, from the shared record.

    A dataset's shared statistics record (``prepare.store.stats.compute_dataset_stats``)
    holds an exact ``reads_total``, an exact ``bases_total`` when seqkit ran, and the GC
    fraction; those are what a profile reports, GC converted to percent. Without a record, or
    with one that carries no read total, the profile can only report what it sampled
    (``sample_gc`` is already in percent), flagged with ``sampled`` so a reader does not
    mistake the sample size for the dataset size.
    """
    stats = dataset_stats or {}
    reads_total = stats.get("reads_total")
    if not reads_total:
        return reads_sampled, int(reads_sampled * mean_read_length), sample_gc, True
    bases_total = stats.get("bases_total") or int(round(mean_read_length * reads_total))
    gc_fraction = stats.get("gc_content")
    gc_percent = sample_gc if gc_fraction is None else gc_fraction * 100
    return int(reads_total), int(bases_total), gc_percent, False


# Remediation actions implied by each anomaly type (ordered for stable output).
_ANOMALY_ACTIONS = {
    "high_contamination": ["Perform adapter trimming", "Check for vector contamination"],
    "poor_quality": ["Apply quality filtering", "Consider excluding from analysis"],
    "low_complexity": ["Check for PCR artifacts", "Verify sample preparation"],
}


def _detect_profile_anomalies(profile) -> List[Tuple[str, str, float]]:
    """Return (anomaly_type, human-readable reason, severity weight) for each anomaly tripped."""
    flags: List[Tuple[str, str, float]] = []
    if profile.complexity_score < 0.1:
        flags.append(("low_complexity", "Very low sequence complexity", 0.3))
    adapter = profile.contamination_indicators.get("adapter_contamination", 0)
    if adapter > 0.1:
        flags.append(("high_contamination", f"High adapter contamination: {adapter:.1%}", 0.4))
    if profile.gc_percent < 20 or profile.gc_percent > 80:
        flags.append(("unusual_gc", f"Unusual GC content: {profile.gc_percent:.1f}%", 0.2))
    if profile.n_content > 0.05:
        flags.append(("high_n_content", f"High N content: {profile.n_content:.1%}", 0.3))
    if profile.quality_grade == "poor":
        flags.append(("poor_quality", "Poor overall quality", 0.5))
    return flags


def _actions_for_anomalies(tripped_types: set) -> List[str]:
    """Remediation actions for the anomaly types a dataset tripped, in a stable order."""
    actions: List[str] = []
    for atype, acts in _ANOMALY_ACTIONS.items():
        if atype in tripped_types:
            actions.extend(acts)
    return actions


def _record_total(dataset_stats: Optional[Dict[str, Any]], files: Sequence[Path]) -> Optional[int]:
    """Reads across ``files`` from the statistics record's per-file counts, or None when any is missing.

    Passed to the quality sampler so it does not count the files again; a record computed for
    a different file set (a file added or renamed since) is not used.
    """
    per_file = (dataset_stats or {}).get("reads_per_file")
    if not isinstance(per_file, dict) or not all(f.name in per_file for f in files):
        return None
    return sum(int(per_file[f.name]) for f in files)


class SRADatasetAnalyzer:
    """Main analyzer for comprehensive SRA dataset analysis."""

    def __init__(self, fastq_dir: Optional[Union[str, Path]] = None):
        self.quality_analyzer = SequenceQualityAnalyzer()
        self.fastq_dir = Path(fastq_dir) if fastq_dir else None

    def profile_dataset_quality(
        self,
        accession: str,
        fastq_path: Optional[Union[str, Path, Sequence[Union[str, Path]]]] = None,
        metadata: Optional[SRADatasetInfo] = None,
        sample_size: int = 10000,
        sampler: str = "uniform",
        dataset_stats: Optional[Dict[str, Any]] = None,
    ) -> QualityProfile:
        """
        Generate comprehensive quality profile for SRA dataset.

        This is the one statistics path of ``sra_profile`` and ``sra_report``: read and base
        totals, mean read length and GC come from the dataset's shared statistics record, the
        per-read quality, complexity and contamination figures from a sample of every file.

        Args:
            accession: SRA accession
            fastq_path: The dataset's FASTQ file, or all of its mate files (located under
                the analyzer's folder when omitted)
            metadata: Dataset metadata (optional)
            sample_size: Reads sampled for the quality/complexity metrics
            sampler: "uniform" (default) or "head"; see ``SequenceQualityAnalyzer.analyze_fastq_quality``
            dataset_stats: The dataset's shared statistics record
                (``prepare.store.stats.compute_dataset_stats``). Loaded (from the store
                sidecar's cache, or computed) when omitted; when it cannot be had at all,
                the totals are the sample's own figures and ``sampled`` is set.

        Returns:
            QualityProfile with comprehensive analysis; GC in percent
        """
        logger.info(f"Profiling dataset quality for {accession}")

        if fastq_path is None:
            files = self.find_fastq_files(accession)
        elif isinstance(fastq_path, (str, Path)):
            files = [Path(fastq_path)]
        else:
            files = [Path(p) for p in fastq_path]
        files = [f for f in files if f.exists()]

        if files:
            if dataset_stats is None:
                dataset_stats = load_dataset_stats(files, sample_size=sample_size)
            quality_metrics = self.quality_analyzer.analyze_fastq_quality(
                files, sample_size=sample_size, sampler=sampler, total_records=_record_total(dataset_stats, files)
            )
        else:
            logger.warning(f"FASTQ file not found for {accession}, using metadata only")
            quality_metrics = {}

        # Extract metrics
        reads_sampled = quality_metrics.get("total_reads_sampled", 0)
        read_length_stats = quality_metrics.get("read_length_stats", {})
        gc_stats = quality_metrics.get("gc_content_stats", {})
        quality_stats = quality_metrics.get("quality_stats", {})
        complexity_metrics = quality_metrics.get("complexity_metrics", {})
        contamination = quality_metrics.get("contamination_indicators", {})

        # Calculate quality grade
        quality_grade = self._calculate_quality_grade(quality_stats, contamination, complexity_metrics)

        # Generate recommendations
        recommendations = self._generate_quality_recommendations(
            quality_stats, contamination, complexity_metrics, metadata
        )

        total_reads, total_bases, gc_percent, sampled = _dataset_totals(
            dataset_stats, reads_sampled, read_length_stats.get("mean", 0), gc_stats.get("mean", 0)
        )
        record_length = (dataset_stats or {}).get("avg_read_length")
        avg_read_length = record_length if (record_length and not sampled) else read_length_stats.get("mean", 0)

        return QualityProfile(
            accession=accession,
            total_reads=total_reads,
            total_bases=total_bases,
            avg_read_length=avg_read_length,
            read_length_distribution=read_length_stats.get("distribution", {}),
            gc_percent=gc_percent,
            gc_histogram=gc_stats.get("histogram", {}),
            quality_distribution=quality_stats.get("distribution", {}),
            n_content=quality_metrics.get("n_content_stats", {}).get("mean", 0),
            contamination_indicators=contamination,
            complexity_score=complexity_metrics.get("complexity_score", 0),
            duplication_rate=quality_metrics.get("duplication_rate"),
            technology_confidence=self._estimate_technology_confidence(metadata),
            quality_grade=quality_grade,
            recommendations=recommendations,
            reads_sampled=reads_sampled,
            sampled=sampled,
            mean_quality=quality_stats.get("mean", 0.0),
        )

    def _collect_group_profiles(self, groups: Dict[str, List[str]]) -> Dict[str, QualityProfile]:
        """Profile every accession across all groups, skipping and logging failures."""
        all_profiles: Dict[str, QualityProfile] = {}
        for accessions in groups.values():
            for accession in accessions:
                try:
                    all_profiles[accession] = self.profile_dataset_quality(accession)
                except Exception as e:
                    logger.warning(f"Failed to profile {accession}: {e}")
        return all_profiles

    def _merge_supplied_profiles(
        self, groups: Dict[str, List[str]], profiles: Dict[str, QualityProfile]
    ) -> Tuple[Dict[str, QualityProfile], List[str]]:
        """Fill in any accession from ``groups`` missing from a supplied ``profiles`` dict.

        A profile already in ``profiles`` is reused as-is; a gap is profiled fresh via
        ``profile_dataset_quality`` so a partial ``profiles`` dict (e.g. only some
        accessions had a saved quality profile) does not silently drop the rest of the
        comparison. Returns ``(profiles, failures)``: a new dict holding every supplied
        profile plus any freshly profiled one, and a human-readable "accession: reason"
        entry for each accession that could not be profiled either way (logged here, the
        same way ``detect_dataset_anomalies`` handles a profiling failure).
        """
        merged = dict(profiles)
        failures: List[str] = []
        all_accessions = [acc for accessions in groups.values() for acc in accessions]
        for accession in all_accessions:
            if accession in merged:
                continue
            try:
                merged[accession] = self.profile_dataset_quality(accession)
            except Exception as e:
                logger.warning(f"Failed to profile {accession}: {e}")
                failures.append(f"{accession}: {e}")
        return merged, failures

    @staticmethod
    def _group_of(accession: str, groups: Dict[str, List[str]]) -> Optional[str]:
        """Return the first group name containing the accession, or None."""
        for group_name, acc_list in groups.items():
            if accession in acc_list:
                return group_name
        return None

    @classmethod
    def _build_comparison_dataframe(cls, all_profiles, groups) -> pd.DataFrame:
        """Assemble a per-accession comparison DataFrame tagged with each group."""
        rows = []
        for accession, profile in all_profiles.items():
            group = cls._group_of(accession, groups)
            if not group:
                continue
            rows.append(
                {
                    "accession": accession,
                    "group": group,
                    "avg_read_length": profile.avg_read_length,
                    "gc_percent": profile.gc_percent,
                    "total_reads": profile.total_reads,
                    "complexity_score": profile.complexity_score,
                    "quality_grade": profile.quality_grade,
                    "adapter_contamination": profile.contamination_indicators.get("adapter_contamination", 0),
                }
            )
        return pd.DataFrame(rows)

    @staticmethod
    def _group_summary_statistics(comparison_df, groups, numeric_cols) -> Dict[str, Dict[str, Any]]:
        """Compute mean/std/median/min/max per numeric column for each non-empty group."""
        summary_stats: dict = {}
        for group in groups.keys():
            group_data = comparison_df[comparison_df["group"] == group]
            if group_data.empty:
                continue
            summary_stats[group] = {}
            for col in numeric_cols:
                if col not in group_data.columns:
                    continue
                values = group_data[col].dropna()
                if not values.empty:
                    summary_stats[group][col] = {
                        "mean": float(values.mean()),
                        "std": float(values.std()),
                        "median": float(values.median()),
                        "min": float(values.min()),
                        "max": float(values.max()),
                    }
        return summary_stats

    def compare_datasets(
        self,
        groups: Dict[str, List[str]],
        metadata_df: Optional[pd.DataFrame] = None,
        profiles: Optional[Dict[str, QualityProfile]] = None,
    ) -> ComparativeAnalysis:
        """
        Perform comparative analysis across dataset groups.

        Args:
            groups: Dictionary mapping group names to lists of accessions
            metadata_df: DataFrame with metadata for all datasets
            profiles: Previously computed profiles keyed by accession. An accession in
                ``groups`` covered here is reused as-is; one missing from ``profiles`` is
                still profiled fresh via ``profile_dataset_quality`` (a partial ``profiles``
                dict never silently drops accessions from the comparison), and a profiling
                failure is logged and excluded, with a note in ``recommendations``.

        Returns:
            ComparativeAnalysis with statistical comparisons
        """
        logger.info(f"Comparing {len(groups)} dataset groups")

        failures: List[str] = []
        if profiles is not None:
            all_profiles, failures = self._merge_supplied_profiles(groups, profiles)
        else:
            all_profiles = self._collect_group_profiles(groups)
        comparison_df = self._build_comparison_dataframe(all_profiles, groups)

        if comparison_df.empty:
            logger.warning("No data available for comparison")
            return ComparativeAnalysis(
                dataset_groups=groups,
                summary_statistics={},
                statistical_tests={},
                outlier_datasets=[],
                clustering_results=None,
                batch_effects={},
                recommendations=["No data available for analysis"] + [f"Could not profile {f}" for f in failures],
                visualization_data={},
            )

        numeric_cols = NUMERIC_COLUMNS
        summary_stats = self._group_summary_statistics(comparison_df, groups, numeric_cols)
        statistical_tests = self._perform_statistical_tests(comparison_df, groups)
        outliers = self._detect_outliers(comparison_df)
        recommendations = self._generate_comparative_recommendations(summary_stats, statistical_tests, outliers)
        recommendations = recommendations + [f"Could not profile {f}" for f in failures]

        return ComparativeAnalysis(
            dataset_groups=groups,
            summary_statistics=summary_stats,
            statistical_tests=statistical_tests,
            outlier_datasets=outliers,
            clustering_results=None,  # Would implement clustering if needed
            batch_effects={},  # Would implement batch effect detection
            recommendations=recommendations,
            visualization_data=self._prepare_visualization_data(comparison_df),
        )

    def detect_dataset_anomalies(
        self,
        accessions: List[str],
        metadata_df: Optional[pd.DataFrame] = None,
        profiles: Optional[Dict[str, QualityProfile]] = None,
    ) -> AnomalyReport:
        """
        Detect anomalies in SRA datasets.

        Args:
            accessions: List of SRA accessions to analyze
            metadata_df: DataFrame with metadata
            profiles: Previously computed profiles keyed by accession. When given, these are
                reused as-is instead of calling ``profile_dataset_quality`` again; an
                accession missing from ``profiles`` is treated the same as a profiling
                failure below.

        Returns:
            AnomalyReport with detected anomalies
        """
        logger.info(f"Detecting anomalies in {len(accessions)} datasets")

        anomalous_datasets = []
        anomaly_types = defaultdict(list)
        severity_scores = {}
        explanations = {}
        recommended_actions = defaultdict(list)

        for accession in accessions:
            try:
                if profiles is not None:
                    profile = profiles.get(accession)
                    if profile is None:
                        raise DataAccessError(f"No supplied quality profile for {accession}")
                else:
                    profile = self.profile_dataset_quality(accession)
                flags = _detect_profile_anomalies(profile)
                for atype, _reason, _weight in flags:
                    anomaly_types[atype].append(accession)

                severity = sum(weight for _atype, _reason, weight in flags)
                if severity > 0.2:  # Threshold for anomaly
                    anomalous_datasets.append(accession)
                    severity_scores[accession] = min(severity, 1.0)
                    explanations[accession] = "; ".join(reason for _atype, reason, _weight in flags)
                    recommended_actions[accession] = _actions_for_anomalies({a for a, _r, _w in flags})

            except Exception as e:
                logger.error(f"Error analyzing {accession} for anomalies: {e}")
                anomaly_types["analysis_failed"].append(accession)
                severity_scores[accession] = 0.1
                explanations[accession] = f"Analysis failed: {str(e)}"

        return AnomalyReport(
            anomalous_datasets=anomalous_datasets,
            anomaly_types=dict(anomaly_types),
            severity_scores=severity_scores,
            explanations=explanations,
            recommended_actions=dict(recommended_actions),
        )

    def recommend_processing_params(
        self, accession: str, profile: Optional[QualityProfile] = None
    ) -> ProcessingRecommendations:
        """
        Generate processing parameter recommendations based on dataset characteristics.

        Args:
            accession: SRA accession
            profile: Quality profile (will generate if not provided)

        Returns:
            ProcessingRecommendations with optimized parameters
        """
        if profile is None:
            profile = self.profile_dataset_quality(accession)

        # Determine recommended pipeline
        if profile.avg_read_length > 1000:
            pipeline = "long_read"
        elif profile.avg_read_length > 250:
            pipeline = "paired_end"
        else:
            pipeline = "short_read"

        # Quality trimming parameters
        quality_trimming = {
            "enabled": profile.quality_grade in ["poor", "fair"],
            "quality_threshold": 20 if profile.quality_grade == "poor" else 15,
            "min_length": max(50, int(profile.avg_read_length * 0.7)),
        }

        # Adapter removal
        adapter_contamination = profile.contamination_indicators.get("adapter_contamination", 0)
        adapter_removal = {
            "enabled": adapter_contamination > 0.01,
            "stringency": "high" if adapter_contamination > 0.1 else "medium",
        }

        # Contamination filtering
        contamination_filtering = {
            "enabled": any(v > 0.05 for v in profile.contamination_indicators.values()),
            "check_vector": True,
            "check_adapters": adapter_contamination > 0.01,
            "complexity_filter": profile.complexity_score < 0.3,
        }

        # Assembly parameters (basic recommendations)
        assembly_parameters = {
            "kmer_size": "auto" if profile.avg_read_length < 150 else 31,
            "coverage_cutoff": 5 if profile.quality_grade in ["good", "excellent"] else 10,
            "error_correction": profile.quality_grade != "excellent",
        }

        # Computational requirements estimation
        computational_requirements = {
            "memory_gb": max(8, int(profile.total_reads / 10000000) * 2),
            "cpu_cores": 4 if profile.total_reads < 50000000 else 8,
            "storage_gb": max(10, int(profile.total_bases / 1000000000) * 5),
        }

        # Processing time estimation
        base_time_hours = profile.total_reads / 5000000  # Rough estimate
        if quality_trimming["enabled"]:
            base_time_hours *= 1.5
        if contamination_filtering["enabled"]:
            base_time_hours *= 1.3

        processing_time = f"{base_time_hours:.1f} hours"

        return ProcessingRecommendations(
            accession=accession,
            recommended_pipeline=pipeline,
            quality_trimming=quality_trimming,
            adapter_removal=adapter_removal,
            contamination_filtering=contamination_filtering,
            assembly_parameters=assembly_parameters,
            expected_coverage=None,  # Would need genome size estimate
            computational_requirements=computational_requirements,
            estimated_processing_time=processing_time,
        )

    def find_fastq_files(self, accession: str) -> List[Path]:
        """Every FASTQ file of an accession under the configured download folder, mates in order.

        See ``prepare.sra.dataset_stats.accession_fastq_files`` for the layouts accepted.
        When no folder was given, ``fastq`` and then ``sra_downloads`` in the working
        directory are searched, and the first that holds the accession wins.
        """
        roots = [self.fastq_dir] if self.fastq_dir else [Path("fastq"), Path("sra_downloads")]
        for root in roots:
            files = accession_fastq_files(root, accession)
            if files:
                return files
        return []

    def find_fastq(self, accession: str) -> Optional[Path]:
        """The first of ``find_fastq_files``' files (mate 1 of a pair), or None."""
        files = self.find_fastq_files(accession)
        return files[0] if files else None

    def _calculate_quality_grade(self, quality_stats: Dict, contamination: Dict, complexity: Dict) -> str:
        """Calculate overall quality grade."""
        score = 0.0

        # Quality score component
        if quality_stats.get("mean", 0) >= 30:
            score += 0.4
        elif quality_stats.get("mean", 0) >= 20:
            score += 0.2

        # Contamination component
        adapter_cont = contamination.get("adapter_contamination", 0)
        if adapter_cont < 0.01:
            score += 0.3
        elif adapter_cont < 0.05:
            score += 0.1

        # Complexity component
        complexity_score = complexity.get("complexity_score", 0)
        if complexity_score > 0.7:
            score += 0.3
        elif complexity_score > 0.4:
            score += 0.15

        if score >= 0.8:
            return "excellent"
        elif score >= 0.6:
            return "good"
        elif score >= 0.4:
            return "fair"
        else:
            return "poor"

    def _generate_quality_recommendations(
        self, quality_stats: Dict, contamination: Dict, complexity: Dict, metadata: Optional[SRADatasetInfo]
    ) -> List[str]:
        """Generate quality improvement recommendations."""
        recommendations = []

        if quality_stats.get("mean", 0) < 20:
            recommendations.append("Apply quality trimming with Q20 threshold")

        if contamination.get("adapter_contamination", 0) > 0.05:
            recommendations.append("Remove adapter sequences")

        if complexity.get("complexity_score", 0) < 0.3:
            recommendations.append("Check for PCR duplicates and low-complexity regions")

        if not recommendations:
            recommendations.append("Dataset appears to be high quality")

        return recommendations

    def _estimate_technology_confidence(self, metadata: Optional[SRADatasetInfo]) -> float:
        """Estimate confidence in technology detection."""
        if not metadata:
            return 0.5

        # Simple heuristic based on metadata completeness
        confidence = 0.5
        if hasattr(metadata, "platform") and metadata.platform:
            confidence += 0.3
        if hasattr(metadata, "instrument") and metadata.instrument:
            confidence += 0.2

        return min(confidence, 1.0)

    def _perform_statistical_tests(
        self, comparison_df: pd.DataFrame, groups: Dict[str, List[str]]
    ) -> Dict[str, Dict[str, Any]]:
        """Perform statistical tests between groups."""
        tests: dict = {}

        numeric_cols = NUMERIC_COLUMNS
        group_names = list(groups.keys())

        if len(group_names) < 2:
            return tests

        for col in numeric_cols:
            if col not in comparison_df.columns:
                continue

            tests[col] = {}

            # Get data for each group
            group_data = []
            for group in group_names:
                data = comparison_df[comparison_df["group"] == group][col].dropna()
                if not data.empty:
                    group_data.append(data.values)

            if len(group_data) >= 2 and min(len(values) for values in group_data) < 2:
                # A group of one has no variance of its own; the test would give NaN and warnings.
                logger.info("Not testing %s: a group has fewer than 2 values", col)
                tests[col] = {
                    "test": "skipped (a group has fewer than 2 values)",
                    "statistic": np.nan,
                    "p_value": np.nan,
                    "significant": False,
                    "note": "Each group needs at least 2 datasets for a statistical comparison",
                }
            elif len(group_data) >= 2:
                # Check if data is nearly identical (would cause precision loss)
                all_values = np.concatenate(group_data)
                variance = np.var(all_values)

                if variance < 1e-10:  # Extremely low variance, skip statistical tests
                    tests[col]["test"] = "skipped (identical values)"
                    tests[col]["statistic"] = np.nan
                    tests[col]["p_value"] = 1.0  # No significant difference for identical values
                    tests[col]["significant"] = False
                    tests[col]["note"] = "Values too similar for meaningful statistical comparison"
                else:
                    # Perform t-test if 2 groups, ANOVA if more
                    try:
                        tests[col]["test"], statistic, p_value = compare_group_means(group_data)
                        tests[col]["statistic"] = float(statistic)
                        tests[col]["p_value"] = float(p_value)
                        tests[col]["significant"] = bool(p_value < 0.05)
                    except RuntimeWarning:
                        # Handle precision loss warnings gracefully
                        tests[col]["test"] = "failed (precision loss)"
                        tests[col]["statistic"] = np.nan
                        tests[col]["p_value"] = np.nan
                        tests[col]["significant"] = False

        # One correction over the tests that ran: four tests at 0.05 each would otherwise
        # report a chance difference far more often than one in twenty.
        ran = [col for col, result in tests.items() if result.get("test") in ("t-test", "ANOVA")]
        for col, adjusted in zip(ran, holm_adjusted([tests[col]["p_value"] for col in ran])):
            tests[col]["p_adjusted"] = adjusted
            tests[col]["significant"] = bool(adjusted < 0.05)
        return tests

    def _detect_outliers(self, comparison_df: pd.DataFrame) -> List[str]:
        """Detect outlier datasets using statistical methods."""
        outliers = []

        numeric_cols = NUMERIC_COLUMNS

        for col in numeric_cols:
            if col not in comparison_df.columns:
                continue

            values = comparison_df[col].dropna()
            if len(values) < 4:  # Need sufficient data
                continue

            # Use IQR method
            q25 = values.quantile(0.25)
            q75 = values.quantile(0.75)
            iqr = q75 - q25
            lower_bound = q25 - 1.5 * iqr
            upper_bound = q75 + 1.5 * iqr

            outlier_mask = (values < lower_bound) | (values > upper_bound)
            outlier_accessions = comparison_df.loc[values[outlier_mask].index, "accession"].tolist()
            outliers.extend(outlier_accessions)

        return list(set(outliers))  # Remove duplicates

    def _generate_comparative_recommendations(
        self, summary_stats: Dict, statistical_tests: Dict, outliers: List[str]
    ) -> List[str]:
        """Generate recommendations based on comparative analysis."""
        recommendations = []

        if statistical_tests:
            significant_tests = [col for col, test in statistical_tests.items() if test.get("significant", False)]
            if significant_tests:
                recommendations.append(
                    f"Significant differences (after Holm correction) detected in: {', '.join(significant_tests)}"
                )
                recommendations.append("Consider including group as covariate in downstream analysis")

        if outliers:
            recommendations.append(
                f"Outlier datasets detected: {', '.join(outliers[:5])}{'...' if len(outliers) > 5 else ''}"
            )
            recommendations.append("Review outlier datasets for quality issues")

        if not recommendations:
            recommendations.append("No major issues detected in comparative analysis")

        return recommendations

    def _prepare_visualization_data(self, comparison_df: pd.DataFrame) -> Dict[str, Any]:
        """Prepare data for visualization."""
        return {
            "boxplot_data": comparison_df.to_dict("records"),
            "summary_table": comparison_df.groupby("group").describe().to_dict(),
            "correlation_matrix": comparison_df.select_dtypes(include=[np.number]).corr().to_dict(),
        }
