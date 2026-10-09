"""
EXTENDED TESTS for sra/reporting.py (64% -> 75%+ coverage)

This file adds tests for untested methods in the SRA reporting module:
- generate_quality_dashboard
- create_comparative_analysis
- _create_quality_plots (full coverage)
- _create_comparative_plots
- Jinja2 template rendering paths

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


class TestQualityDashboardGeneration:
    """Test generate_quality_dashboard method."""

    def test_generate_quality_dashboard_success(self, tmp_output_dir, mock_quality_profiles, mock_anomaly_report):
        """Test successful quality dashboard generation."""
        generator = SRAReportGenerator(tmp_output_dir)

        # Mock analyzer methods
        def profile_side_effect(accession, *args, **kwargs):
            return mock_quality_profiles[accession]

        with patch.object(generator.analyzer, "profile_dataset_quality", side_effect=profile_side_effect):
            with patch.object(generator.analyzer, "detect_dataset_anomalies", return_value=mock_anomaly_report):
                result_path = generator.generate_quality_dashboard(
                    accessions=["SRR001", "SRR002", "SRR003"], title="Test Quality Dashboard"
                )

        # Verify dashboard was created
        assert result_path.exists()
        assert result_path.suffix == ".html"
        assert "quality_dashboard" in result_path.name

        # Verify content
        html_content = result_path.read_text()
        assert "Test Quality Dashboard" in html_content

    def test_generate_quality_dashboard_with_jinja2(self, tmp_output_dir, mock_quality_profiles, mock_anomaly_report):
        """Test quality dashboard generation with Jinja2 templating."""
        generator = SRAReportGenerator(tmp_output_dir)

        def profile_side_effect(accession, *args, **kwargs):
            return mock_quality_profiles[accession]

        with patch.object(generator.analyzer, "profile_dataset_quality", side_effect=profile_side_effect):
            with patch.object(generator.analyzer, "detect_dataset_anomalies", return_value=mock_anomaly_report):
                with patch("jinja2.Environment") as mock_env:
                    mock_template = Mock()
                    mock_template.render.return_value = "<html><h1>Quality Dashboard</h1></html>"
                    mock_env.return_value.from_string.return_value = mock_template

                    result_path = generator.generate_quality_dashboard(
                        accessions=["SRR001", "SRR002"], title="Jinja2 Dashboard"
                    )

        assert result_path.exists()
        # Verify Jinja2 template was used
        mock_template.render.assert_called_once()
        call_kwargs = mock_template.render.call_args[1]
        assert "title" in call_kwargs
        assert call_kwargs["title"] == "Jinja2 Dashboard"

    def test_generate_quality_dashboard_with_plots(self, tmp_output_dir, mock_quality_profiles, mock_anomaly_report):
        """Test quality dashboard generation with Plotly plots."""
        generator = SRAReportGenerator(tmp_output_dir)

        def profile_side_effect(accession, *args, **kwargs):
            return mock_quality_profiles[accession]

        with patch.object(generator.analyzer, "profile_dataset_quality", side_effect=profile_side_effect):
            with patch.object(generator.analyzer, "detect_dataset_anomalies", return_value=mock_anomaly_report):
                with patch("plotly.graph_objects.Figure"):
                    with patch("plotly.offline.plot", return_value="<div>Plot</div>"):
                        result_path = generator.generate_quality_dashboard(
                            accessions=["SRR001", "SRR002"], title="Dashboard with Plots"
                        )

        assert result_path.exists()

    def test_generate_quality_dashboard_partial_failures(self, tmp_output_dir, mock_quality_profiles):
        """Test quality dashboard when some profiles fail."""
        generator = SRAReportGenerator(tmp_output_dir)

        def profile_side_effect(accession, *args, **kwargs):
            if accession == "SRR001":
                return mock_quality_profiles["SRR001"]
            elif accession == "SRR002":
                return mock_quality_profiles["SRR002"]
            else:
                raise Exception("Failed to profile")

        mock_anomaly = MockAnomalyReport(anomalous_datasets=[], explanations={})

        with patch.object(generator.analyzer, "profile_dataset_quality", side_effect=profile_side_effect):
            with patch.object(generator.analyzer, "detect_dataset_anomalies", return_value=mock_anomaly):
                # Should succeed with partial data
                result_path = generator.generate_quality_dashboard(
                    accessions=["SRR001", "SRR002", "SRR003"], title="Partial Dashboard"
                )

        assert result_path.exists()

    def test_generate_quality_dashboard_all_failures(self, tmp_output_dir):
        """Test quality dashboard when all profiles fail."""
        generator = SRAReportGenerator(tmp_output_dir)

        with patch.object(generator.analyzer, "profile_dataset_quality", side_effect=Exception("API Error")):
            with pytest.raises(ValueError, match="No datasets could be profiled"):
                generator.generate_quality_dashboard(accessions=["SRR001", "SRR002"], title="Failed Dashboard")

    def test_generate_quality_dashboard_reuses_supplied_profile(self, tmp_output_dir):
        """An accession covered by the ``profiles`` argument is not reprofiled from FASTQ."""
        generator = SRAReportGenerator(tmp_output_dir)

        profile = QualityProfile(
            accession="SRRX",
            total_reads=12345,
            total_bases=1850000,
            avg_read_length=150.0,
            read_length_distribution={},
            gc_percent=55.0,
            gc_histogram={},
            quality_distribution={"excellent_q30+": 0.9},
            n_content=0.01,
            contamination_indicators={"adapter_contamination": 0.01},
            complexity_score=0.9,
            duplication_rate=0.1,
            technology_confidence=0.8,
            quality_grade="excellent",
            recommendations=[],
        )
        mock_anomaly = MockAnomalyReport(anomalous_datasets=[], explanations={})

        with patch.object(generator.analyzer, "profile_dataset_quality") as mock_profile_call:
            with patch.object(generator.analyzer, "detect_dataset_anomalies", return_value=mock_anomaly):
                result_path = generator.generate_quality_dashboard(
                    accessions=["SRRX"],
                    title="Reused Profile Dashboard",
                    profiles={"SRRX": profile},
                )

        mock_profile_call.assert_not_called()
        assert result_path.exists()
        # The supplied profile's own values (not some freshly computed default) reached the
        # rendered summary: total_reads=12,345 and GC 55.0%.
        html_content = result_path.read_text()
        assert "12,345" in html_content
        assert "55.0%" in html_content


class TestComparativeAnalysisReports:
    """Test create_comparative_analysis method."""

    def test_create_comparative_analysis_success(self, tmp_output_dir, mock_comparative_analysis):
        """Test successful comparative analysis report generation."""
        generator = SRAReportGenerator(tmp_output_dir)

        groups = {"Group A": ["SRR001", "SRR002"], "Group B": ["SRR003", "SRR004"]}

        with patch.object(generator.analyzer, "compare_datasets", return_value=mock_comparative_analysis):
            result_path = generator.create_comparative_analysis(groups=groups, title="Test Comparative Analysis")

        # Verify report was created
        assert result_path.exists()
        assert result_path.suffix == ".html"
        assert "comparative_analysis" in result_path.name

        # Verify content
        html_content = result_path.read_text()
        assert "Group A" in html_content
        assert "Group B" in html_content

    def test_create_comparative_analysis_forwards_supplied_profiles(self, tmp_output_dir, mock_comparative_analysis):
        """profiles is threaded through to compare_datasets: the HTML report path must not
        silently reprofile from FASTQ when the caller already has saved profiles."""
        generator = SRAReportGenerator(tmp_output_dir)

        groups = {"Group A": ["SRR001", "SRR002"]}
        profiles = {"SRR001": Mock(), "SRR002": Mock()}

        with patch.object(
            generator.analyzer, "compare_datasets", return_value=mock_comparative_analysis
        ) as mock_compare:
            with patch.object(generator.analyzer, "profile_dataset_quality") as mock_profile:
                generator.create_comparative_analysis(groups=groups, title="X", profiles=profiles)

        mock_compare.assert_called_once_with(groups, profiles=profiles)
        mock_profile.assert_not_called()

    def test_create_comparative_analysis_with_jinja2(self, tmp_output_dir, mock_comparative_analysis):
        """Test comparative analysis with Jinja2 templating."""
        generator = SRAReportGenerator(tmp_output_dir)

        groups = {"Treatment": ["SRR001"], "Control": ["SRR002"]}

        with patch.object(generator.analyzer, "compare_datasets", return_value=mock_comparative_analysis):
            with patch("jinja2.Environment") as mock_env:
                mock_template = Mock()
                mock_template.render.return_value = "<html><h1>Comparative Analysis</h1></html>"
                mock_env.return_value.from_string.return_value = mock_template

                result_path = generator.create_comparative_analysis(groups=groups, title="Jinja2 Comparative Analysis")

        assert result_path.exists()
        mock_template.render.assert_called_once()

    def test_create_comparative_analysis_with_plots(self, tmp_output_dir, mock_comparative_analysis):
        """Test comparative analysis with Plotly plots."""
        generator = SRAReportGenerator(tmp_output_dir)

        groups = {"Group A": ["SRR001", "SRR002"], "Group B": ["SRR003", "SRR004"]}

        with patch.object(generator.analyzer, "compare_datasets", return_value=mock_comparative_analysis):
            with patch("plotly.graph_objects.Figure"):
                with patch("plotly.offline.plot", return_value="<div>Boxplot</div>"):
                    result_path = generator.create_comparative_analysis(groups=groups, title="Analysis with Plots")

        assert result_path.exists()


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


class TestJinja2TemplateRendering:
    """Test Jinja2 template rendering paths."""

    def test_generate_quality_html_with_anomalies(self, tmp_output_dir):
        """Test quality HTML generation with anomaly data."""
        generator = SRAReportGenerator(tmp_output_dir)

        dashboard_data = {
            "title": "Test Dashboard",
            "timestamp": "2025-01-01 12:00:00",
            "total_datasets": 3,
            "summary_stats": {
                "total_reads": 5000000,
                "average_gc_percent": 45.0,
                "high_contamination_count": 1,
                "quality_grade_distribution": {"excellent": 2, "good": 1},
            },
            "anomaly_report": MockAnomalyReport(
                anomalous_datasets=["SRR003"], explanations={"SRR003": "High contamination"}
            ),
            "plots": {},
        }

        with patch("jinja2.Environment") as mock_env:
            mock_template = Mock()
            mock_template.render.return_value = "<html>Dashboard with Anomalies</html>"
            mock_env.return_value.from_string.return_value = mock_template

            html = generator._generate_quality_html(dashboard_data)

        assert "<html>" in html
        mock_template.render.assert_called_once()

    def test_generate_comparative_html_with_statistical_tests(self, tmp_output_dir):
        """Test comparative HTML generation with statistical test results."""
        generator = SRAReportGenerator(tmp_output_dir)

        report_data = {
            "title": "Statistical Comparison",
            "timestamp": "2025-01-01 12:00:00",
            "group_counts": {"Group A": 10, "Group B": 15},
            "comparison": MockComparativeAnalysis(
                statistical_tests={"gc_percent": MockStatisticalTest("Mann-Whitney U", 0.03, True)},
                visualization_data={},
                recommendations=["Use different pipelines"],
            ),
            "plots": {},
        }

        with patch("jinja2.Environment") as mock_env:
            mock_template = Mock()
            mock_template.render.return_value = "<html>Statistical Results</html>"
            mock_env.return_value.from_string.return_value = mock_template

            html = generator._generate_comparative_html(report_data)

        assert "<html>" in html
        mock_template.render.assert_called_once()


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
