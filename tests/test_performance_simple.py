"""
EDGE CASE TESTS for PREPARE

Tests edge cases for validated components:
- Boundary conditions
- Numerical edge cases
- String encoding edge cases

The timed and memory-bounded performance checks that used to live in this module (a
pytest-benchmark timing of technology detection, plot-data preparation and metadata-report
writing, none of which asserted a threshold, and a memory-efficiency class whose tests asserted
nothing about memory) moved to tests/test_performance_regressions.py, where every check has an
explicit bound; see that module for the timing and memory guards, run with ``pytest -m perf``.

Run: pytest tests/test_performance_simple.py -v
"""

import pandas as pd
import numpy as np

from prepare.data.sra_metadata import (
    detect_sequencing_technology,
    SRADatasetInfo,
)
from prepare.plugins.visualizers.bar import (
    BarChartPlugin,
    _prepare_plot_data,
)

# ============================================================================
# TEST CLASS: Numerical Edge Cases
# ============================================================================


class TestNumericalEdgeCases:
    """Test edge cases with numerical values."""

    def test_zero_containment_values(self):
        """Test handling of zero values."""
        df = pd.DataFrame(
            {
                "category": ["A", "B", "C"],
                "value": [0.0, 0.0, 0.0],
            }
        )

        result_df, result_y = _prepare_plot_data(df, "category", "value", None)

        assert len(result_df) == 3
        assert all(result_df["value"] == 0.0)

    def test_very_large_numbers(self):
        """Test handling of very large numbers."""
        df = pd.DataFrame(
            {
                "category": ["A", "B"],
                "value": [1e15, 2e15],
            }
        )

        result_df, result_y = _prepare_plot_data(df, "category", "value", None)

        assert len(result_df) == 2
        assert result_df["value"].max() == 2e15

    def test_very_small_numbers(self):
        """Test handling of very small numbers."""
        df = pd.DataFrame(
            {
                "category": ["A", "B"],
                "value": [1e-10, 1e-9],
            }
        )

        result_df, result_y = _prepare_plot_data(df, "category", "value", None)

        assert len(result_df) == 2
        assert result_df["value"].min() > 0

    def test_mixed_positive_negative(self):
        """Test handling of mixed positive and negative values."""
        df = pd.DataFrame(
            {
                "category": ["A", "B", "C", "D"],
                "value": [-10, 5, -3, 8],
            }
        )

        result_df, result_y = _prepare_plot_data(df, "category", "value", None)

        assert len(result_df) == 4

    def test_identical_values(self):
        """Test dataset with all identical values."""
        df = pd.DataFrame(
            {
                "category": [f"Cat_{i}" for i in range(100)],
                "value": [0.75] * 100,
            }
        )

        result_df, result_y = _prepare_plot_data(df, "category", "value", 10)

        # Should return top 10 (or first 10 since all equal)
        assert len(result_df) == 10


# ============================================================================
# TEST CLASS: Boundary Conditions
# ============================================================================


class TestBoundaryConditions:
    """Test boundary conditions and limits."""

    def test_empty_dataframe(self):
        """Test processing empty DataFrame."""
        empty_df = pd.DataFrame(columns=["category", "value"])

        result_df, result_y = _prepare_plot_data(empty_df, "category", "value", None)

        assert result_df.empty

    def test_single_row_dataframe(self):
        """Test processing single-row DataFrame."""
        single_df = pd.DataFrame(
            {
                "category": ["A"],
                "value": [10],
            }
        )

        result_df, result_y = _prepare_plot_data(single_df, "category", "value", None)

        assert len(result_df) == 1

    def test_duplicate_categories(self):
        """Test handling of duplicate categories."""
        df = pd.DataFrame(
            {
                "category": ["A", "A", "B", "B"],
                "value": [10, 20, 15, 25],
            }
        )

        result_df, result_y = _prepare_plot_data(df, "category", "value", None)

        # Should handle duplicates (behavior depends on implementation)
        assert len(result_df) >= 2

    def test_very_long_category_names(self):
        """Test handling of very long strings."""
        long_name = "A" * 1000
        df = pd.DataFrame(
            {
                "category": [long_name, "B"],
                "value": [10, 20],
            }
        )

        result_df, result_y = _prepare_plot_data(df, "category", "value", None)

        assert len(result_df) == 2


# ============================================================================
# TEST CLASS: String and Encoding Edge Cases
# ============================================================================


class TestStringEncodingEdgeCases:
    """Test edge cases with strings and encodings."""

    def test_unicode_characters(self):
        """Test handling of unicode characters."""
        df = pd.DataFrame(
            {
                "category": ["Français", "中文", "Русский"],  # ascii-ok: non-Latin labels
                "value": [10, 20, 30],
            }
        )

        result_df, result_y = _prepare_plot_data(df, "category", "value", None)

        assert len(result_df) == 3

    def test_emoji_in_strings(self):
        """Test handling of emoji characters."""
        df = pd.DataFrame(
            {
                "category": ["Test 🦠", "Sample 🧬", "Data 📊"],  # ascii-ok: emoji labels
                "value": [10, 20, 30],
            }
        )

        result_df, result_y = _prepare_plot_data(df, "category", "value", None)

        assert len(result_df) == 3

    def test_empty_strings(self):
        """Test handling of empty strings."""
        df = pd.DataFrame(
            {
                "category": ["", "A", ""],
                "value": [10, 20, 30],
            }
        )

        result_df, result_y = _prepare_plot_data(df, "category", "value", None)

        assert len(result_df) == 3

    def test_special_characters(self):
        """Test handling of special characters."""
        df = pd.DataFrame(
            {
                "category": ["A/B", "C\\D", "E|F", "G*H"],
                "value": [10, 20, 30, 40],
            }
        )

        result_df, result_y = _prepare_plot_data(df, "category", "value", None)

        assert len(result_df) == 4

    def test_whitespace_variations(self):
        """Test handling of various whitespace."""
        df = pd.DataFrame(
            {
                "category": ["  A  ", "\tB\t", "\nC\n", "D"],
                "value": [10, 20, 30, 40],
            }
        )

        result_df, result_y = _prepare_plot_data(df, "category", "value", None)

        assert len(result_df) == 4


# ============================================================================
# TEST CLASS: Technology Detection Edge Cases
# ============================================================================


class TestTechnologyDetectionEdgeCases:
    """Test edge cases in technology detection."""

    def test_unknown_platform(self):
        """Test detection with unknown platform."""
        dataset = SRADatasetInfo(
            accession="SRR001",
            title="Test",
            organism="Unknown",
            platform="UNKNOWN_PLATFORM",
            instrument="Unknown",
            strategy="WGS",
            layout="PAIRED",
            spots=1000,
            bases=150000,
            avg_length=150.0,
            size_mb=100.0,
            release_date="2023-01-01",
            bioproject="",
            biosample="",
            library_selection="RANDOM",
            library_source="GENOMIC",
        )

        tech = detect_sequencing_technology(dataset)
        assert tech == "unknown"

    def test_mixed_case_platforms(self):
        """Test platform detection with mixed case."""
        dataset = SRADatasetInfo(
            accession="SRR001",
            title="Test",
            organism="E. coli",
            platform="iLLuMiNa",  # Mixed case
            instrument="HiSeq",
            strategy="WGS",
            layout="PAIRED",
            spots=1000,
            bases=150000,
            avg_length=150.0,
            size_mb=100.0,
            release_date="2023-01-01",
            bioproject="",
            biosample="",
            library_selection="RANDOM",
            library_source="GENOMIC",
        )

        tech = detect_sequencing_technology(dataset)
        assert tech == "illumina"

    def test_platform_in_instrument_field(self):
        """Test detection when platform is in instrument field."""
        dataset = SRADatasetInfo(
            accession="SRR001",
            title="Test",
            organism="E. coli",
            platform="",
            instrument="Illumina MiSeq",  # Platform in instrument
            strategy="WGS",
            layout="PAIRED",
            spots=1000,
            bases=150000,
            avg_length=150.0,
            size_mb=100.0,
            release_date="2023-01-01",
            bioproject="",
            biosample="",
            library_selection="RANDOM",
            library_source="GENOMIC",
        )

        tech = detect_sequencing_technology(dataset)
        assert tech == "illumina"


# ============================================================================
# TEST CLASS: Visualization Edge Cases
# ============================================================================


class TestVisualizationEdgeCases:
    """Test edge cases in visualization."""

    def test_plot_single_data_point(self, tmp_path):
        """Test plotting single data point."""
        df = pd.DataFrame(
            {
                "category": ["A"],
                "count": [10],
            }
        )

        output_file = tmp_path / "single_point.png"
        BarChartPlugin.create_plot(data=df, x_column="category", y_column="count", output_file=str(output_file))

        assert output_file.exists()

    def test_plot_with_zero_values(self, tmp_path):
        """Test plotting with zero values."""
        df = pd.DataFrame(
            {
                "category": ["A", "B", "C"],
                "count": [0, 0, 0],
            }
        )

        output_file = tmp_path / "zero_values.png"
        BarChartPlugin.create_plot(data=df, x_column="category", y_column="count", output_file=str(output_file))

        assert output_file.exists()

    def test_plot_many_categories(self, tmp_path):
        """Test plotting with many categories."""
        df = pd.DataFrame(
            {
                "category": [f"Cat_{i}" for i in range(100)],
                "count": np.random.randint(1, 100, 100),
            }
        )

        output_file = tmp_path / "many_categories.png"
        BarChartPlugin.create_plot(
            data=df, x_column="category", y_column="count", limit=20, output_file=str(output_file)  # Use limit
        )

        assert output_file.exists()


# ============================================================================
# SUCCESS METRICS:
#
# After running these tests:
# - Numerical edge cases (zero, large, small, negative)
# - Boundary conditions (empty, single row, duplicates)
# - String encoding (unicode, emoji, special chars)
# - Technology detection edge cases
#
# Run tests:
#   pytest tests/test_performance_simple.py -v
# ============================================================================
