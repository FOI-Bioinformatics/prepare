"""
Advanced SRA Analytics and Reporting Package for PREPARE.

This package provides comprehensive dataset quality analysis, statistical
comparison, and interactive HTML reporting for SRA-derived sequencing data.

Main Components:
- analytics: Quality profiling and comparative dataset analysis
- dataset_stats: A dataset's FASTQ files and its shared statistics record
- profiles: Per-accession profile JSONs and the sra_profile statistics row
- reporting: Interactive HTML reports and dashboards
"""

from .analytics import (
    SRADatasetAnalyzer,
    SequenceQualityAnalyzer,
    QualityProfile,
    ComparativeAnalysis,
    AnomalyReport,
    ProcessingRecommendations,
    json_safe,
)
from .profiles import load_quality_profiles

from .reporting import SRAReportGenerator

__version__ = "1.0.0"
__author__ = "PREPARE Development Team"

__all__ = [
    # Analytics
    "SRADatasetAnalyzer",
    "SequenceQualityAnalyzer",
    "QualityProfile",
    "ComparativeAnalysis",
    "AnomalyReport",
    "ProcessingRecommendations",
    "load_quality_profiles",
    "json_safe",
    # Reporting
    "SRAReportGenerator",
]
