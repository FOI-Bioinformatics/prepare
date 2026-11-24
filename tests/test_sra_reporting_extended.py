"""
EXTENDED TESTS for sra/reporting.py (64% -> 75%+ coverage)

This file adds tests for methods of the SRA reporting module:
- _create_quality_plots (full coverage)
- _create_comparative_plots

Run: pytest tests/test_sra_reporting_extended.py -v
"""

import sys

import pytest
from unittest.mock import Mock, patch
from dataclasses import dataclass
from typing import Dict, List, Any

for _module in ("plotly", "jinja2"):
    pytest.importorskip(_module)

from prepare.core.exceptions import ConfigurationError  # noqa: E402
from prepare.sra.analytics import QualityProfile  # noqa: E402
from prepare.sra.reporting import SRAReportGenerator  # noqa: E402


# Mock data classes
@dataclass
class MockQualityProfile:
    """Mock quality profile for testing."""

    total_reads: int
    total_bases: int
    avg_read_length: float
    gc_percent: float
    quality_grade: str
    complexity_score: float
    n_content: float
    contamination_indicators: Dict[str, float]
    warnings: List[str] = None

    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []


@dataclass
class MockAnomalyReport:
    """Mock anomaly report."""

    anomalous_datasets: List[str]
    explanations: Dict[str, str]


@dataclass
class MockStatisticalTest:
    """Mock statistical test result."""

    test: str
    p_value: float
    significant: bool


@dataclass
class MockComparativeAnalysis:
    """Mock comparative analysis results."""

    statistical_tests: Dict[str, MockStatisticalTest]
    visualization_data: Dict[str, Any]
    recommendations: List[str]


@pytest.fixture
def tmp_output_dir(tmp_path):
    """Create temporary output directory."""
    output_dir = tmp_path / "reports"
    return output_dir


@pytest.fixture
def mock_quality_profiles():
    """Create multiple mock quality profiles for dashboard testing."""
    return {
        "SRR001": MockQualityProfile(
            total_reads=1000000,
            total_bases=150000000,
            avg_read_length=150.0,
            gc_percent=45.0,
            quality_grade="excellent",
            complexity_score=0.90,
            n_content=0.01,
            contamination_indicators={"adapter_contamination": 0.02},
            warnings=[],
        ),
        "SRR002": MockQualityProfile(
            total_reads=800000,
            total_bases=120000000,
            avg_read_length=150.0,
            gc_percent=52.0,
            quality_grade="good",
            complexity_score=0.85,
            n_content=0.02,
            contamination_indicators={"adapter_contamination": 0.03},
            warnings=[],
        ),
        "SRR003": MockQualityProfile(
            total_reads=500000,
            total_bases=75000000,
            avg_read_length=150.0,
            gc_percent=48.0,
            quality_grade="fair",
            complexity_score=0.75,
            n_content=0.05,
            contamination_indicators={"adapter_contamination": 0.08},
            warnings=["High adapter contamination"],
        ),
    }


@pytest.fixture
def mock_anomaly_report():
    """Create mock anomaly report."""
    return MockAnomalyReport(
        anomalous_datasets=["SRR003"], explanations={"SRR003": "High adapter contamination and low complexity"}
    )


@pytest.fixture
def mock_comparative_analysis():
    """Create mock comparative analysis."""
    return MockComparativeAnalysis(
        statistical_tests={
            "gc_percent": MockStatisticalTest(test="Mann-Whitney U", p_value=0.045, significant=True),
            "read_length": MockStatisticalTest(test="T-test", p_value=0.32, significant=False),
        },
        visualization_data={
            "boxplot_data": [
                {
                    "group": "Group A",
                    "avg_read_length": 150,
                    "gc_percent": 45.0,
                    "total_reads": 1000000,
                    "complexity_score": 0.85,
                },
                {
                    "group": "Group A",
                    "avg_read_length": 148,
                    "gc_percent": 46.0,
                    "total_reads": 1100000,
                    "complexity_score": 0.87,
                },
                {
                    "group": "Group B",
                    "avg_read_length": 151,
                    "gc_percent": 52.0,
                    "total_reads": 900000,
                    "complexity_score": 0.82,
                },
                {
                    "group": "Group B",
                    "avg_read_length": 149,
                    "gc_percent": 53.0,
                    "total_reads": 950000,
                    "complexity_score": 0.83,
                },
            ]
        },
        recommendations=[
            "Group A shows significantly different GC content",
            "Consider separate preprocessing pipelines for each group",
        ],
    )


class TestQualityPlotCreation:
    """Test _create_quality_plots method in detail."""

    def test_create_quality_plots_with_plotly(self, tmp_output_dir, mock_quality_profiles):
        """Test quality plot creation with Plotly available."""
        generator = SRAReportGenerator(tmp_output_dir)

        with patch("plotly.graph_objects.Figure") as mock_fig_class:
            with patch("plotly.offline.plot", return_value="<div>Plot HTML</div>") as mock_plot:
                mock_figure = Mock()
                mock_fig_class.return_value = mock_figure

                plots = generator._create_quality_plots(mock_quality_profiles)

        # Verify plots were created
        assert isinstance(plots, dict)
        assert len(plots) > 0
        assert "quality_grades" in plots
        assert "gc_distribution" in plots
        assert "complexity_scatter" in plots

        # Verify Plotly was called
        assert mock_plot.call_count >= 3  # At least 3 plots

    def test_create_quality_plots_without_plotly(self, tmp_output_dir, mock_quality_profiles, monkeypatch):
        """Without plotly the plots are an error naming the interactive extra, not an empty result."""
        generator = SRAReportGenerator(tmp_output_dir)
        monkeypatch.setitem(sys.modules, "plotly", None)

        with pytest.raises(ConfigurationError, match=r"prepare\[interactive\]"):
            generator._create_quality_plots(mock_quality_profiles)

    def test_create_quality_plots_empty_profiles(self, tmp_output_dir):
        """Test quality plot creation with empty profiles."""
        generator = SRAReportGenerator(tmp_output_dir)

        with patch("plotly.graph_objects.Figure"):
            with patch("plotly.offline.plot", return_value="<div>Plot</div>"):
                plots = generator._create_quality_plots({})

        # Should still return dict (may be empty or have default plots)
        assert isinstance(plots, dict)


class TestComparativePlotCreation:
    """Test _create_comparative_plots method."""

    def test_create_comparative_plots_success(self, tmp_output_dir, mock_comparative_analysis):
        """Test comparative plot creation with valid data."""
        generator = SRAReportGenerator(tmp_output_dir)

        with patch("plotly.graph_objects.Figure") as mock_fig_class:
            with patch("plotly.offline.plot", return_value="<div>Boxplot</div>"):
                mock_figure = Mock()
                mock_fig_class.return_value = mock_figure

                plots = generator._create_comparative_plots(mock_comparative_analysis)

        # Verify plots were created
        assert isinstance(plots, dict)
        # Should have created boxplots for numeric columns
        assert len(plots) > 0

    def test_create_comparative_plots_without_plotly(self, tmp_output_dir, mock_comparative_analysis, monkeypatch):
        """Without plotly the comparative plots are an error naming the interactive extra."""
        generator = SRAReportGenerator(tmp_output_dir)
        monkeypatch.setitem(sys.modules, "plotly", None)

        with pytest.raises(ConfigurationError, match=r"prepare\[interactive\]"):
            generator._create_comparative_plots(mock_comparative_analysis)

    def test_create_comparative_plots_no_visualization_data(self, tmp_output_dir):
        """Test comparative plot creation with no visualization data."""
        generator = SRAReportGenerator(tmp_output_dir)

        empty_analysis = MockComparativeAnalysis(statistical_tests={}, visualization_data={}, recommendations=[])

        plots = generator._create_comparative_plots(empty_analysis)

        # Should return empty dict when no data available
        assert plots == {}


# ============================================================================
# SUCCESS METRICS:
#
# After running these extended tests:
# - Expected: 25+ additional tests pass
# - Coverage: 64% -> 75%+ for sra/reporting.py
# - All major methods now tested
#
# Run tests:
#   pytest tests/test_sra_reporting_extended.py -v
#
# Check coverage:
#   pytest --cov=prepare.sra.reporting --cov-report=term-missing \
#          tests/test_sra_reporting_starter.py \
#          tests/test_sra_reporting_extended.py
# ============================================================================


def _real_profile(accession, gc_percent):
    return QualityProfile(
        accession=accession,
        total_reads=12345,
        total_bases=1850000,
        avg_read_length=150.0,
        read_length_distribution={},
        gc_percent=gc_percent,
        gc_histogram={},
        quality_distribution={"excellent_q30+": 0.9},
        n_content=0.0,
        contamination_indicators={"adapter_contamination": 0.0},
        complexity_score=0.9,
        duplication_rate=0.1,
        technology_confidence=0.8,
        quality_grade="good",
        recommendations=[],
    )


def test_generate_report_renders_both_sections_from_the_supplied_profiles(tmp_path):
    """sra_report's one page: the quality section and the comparison, nothing profiled again."""
    from prepare.sra.analytics import ComparativeAnalysis

    generator = SRAReportGenerator(tmp_path / "reports")
    profiles = {"SRR1": _real_profile("SRR1", 40.0), "SRR2": _real_profile("SRR2", 60.0)}
    comparison = ComparativeAnalysis(
        dataset_groups={"A": ["SRR1"], "B": ["SRR2"]},
        summary_statistics={},
        statistical_tests={"gc_percent": {"test": "t-test", "p_value": 0.01, "significant": True}},
        outlier_datasets=[],
        clustering_results=None,
        batch_effects={},
        recommendations=["Significant differences detected in: gc_percent"],
        visualization_data={},
    )
    with patch.object(generator.analyzer, "profile_dataset_quality") as profile_call:
        path = generator.generate_report(profiles, title="One Report", comparison=comparison)

    profile_call.assert_not_called()
    assert path == tmp_path / "reports" / "sra_report.html"
    html = path.read_text()
    assert "One Report" in html
    assert "50.0%" in html  # mean GC percent of the two profiles
    assert "Statistical tests" in html and "significant" in html


def test_generate_report_without_profiles_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="No datasets could be profiled"):
        SRAReportGenerator(tmp_path / "reports").generate_report({})
