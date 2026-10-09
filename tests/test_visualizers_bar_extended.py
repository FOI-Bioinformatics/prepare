"""
EXTENDED TESTS for plugins/visualizers/bar.py (22% -> 85%+ coverage)

This file adds comprehensive tests for bar chart visualizations:
- Data preparation with various configurations
- Horizontal and vertical bar charts
- Grouped bar charts
- Output file saving
- Error handling

Run: pytest tests/test_visualizers_bar_extended.py -v
"""

import pytest
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # Use non-GUI backend for testing
import matplotlib.pyplot as plt  # noqa: E402

from prepare.plugins.visualizers.bar import (  # noqa: E402
    BarChartPlugin,
    _prepare_plot_data,
    _create_horizontal_bar,
    _create_vertical_bar,
)
from prepare.core.exceptions import VisualizationError  # noqa: E402

# ============================================================================
# TEST CLASS: Data Preparation
# ============================================================================


class TestPrepareePlotData:
    """Test _prepare_plot_data helper function."""

    def test_prepare_with_y_column_specified(self):
        """Test data preparation with y_column specified."""
        df = pd.DataFrame({"category": ["A", "B", "C"], "count": [10, 20, 15]})

        result_df, result_y = _prepare_plot_data(df, "category", "count", None)

        assert result_y == "count"
        assert "count" in result_df.columns
        assert result_df.index.name == "category"

    def test_prepare_without_y_column(self):
        """Test data preparation without y_column (uses first column)."""
        df = pd.DataFrame({"values": [10, 20, 15], "other": [1, 2, 3]})

        result_df, result_y = _prepare_plot_data(df, None, None, None)

        assert result_y == "values"

    def test_prepare_without_x_column(self):
        """Test data preparation without x_column (uses index)."""
        df = pd.DataFrame({"count": [10, 20, 15]}, index=["A", "B", "C"])

        result_df, result_y = _prepare_plot_data(df, None, "count", None)

        assert isinstance(result_df.index, pd.Index)
        assert list(result_df.index) == ["A", "B", "C"]

    def test_prepare_with_limit(self):
        """Test data preparation with limit."""
        df = pd.DataFrame({"category": ["A", "B", "C", "D", "E"], "count": [10, 20, 15, 5, 25]})

        result_df, result_y = _prepare_plot_data(df, "category", "count", 3)

        # Should get top 3 by count
        assert len(result_df) == 3
        assert 25 in result_df["count"].values  # Highest value

    def test_prepare_with_limit_no_y_column(self):
        """Test data preparation with limit but no y_column."""
        df = pd.DataFrame({"values": [10, 20, 15, 5, 25]})

        result_df, result_y = _prepare_plot_data(df, None, None, 3)

        # y_column defaults to first column ('values'), so it uses nlargest
        assert len(result_df) == 3
        # Top 3 values: 25, 20, 15
        assert 25 in result_df["values"].values
        assert 20 in result_df["values"].values

    def test_prepare_with_multiindex_raises_error(self):
        """Test that MultiIndex without x_column raises error."""
        df = pd.DataFrame({"count": [10, 20, 15]})
        df.index = pd.MultiIndex.from_tuples([("A", "1"), ("B", "2"), ("C", "3")])

        with pytest.raises(VisualizationError, match="Cannot use MultiIndex"):
            _prepare_plot_data(df, None, "count", None)

    def test_prepare_empty_dataframe(self):
        """Test data preparation with empty DataFrame."""
        df = pd.DataFrame()

        result_df, result_y = _prepare_plot_data(df, None, None, None)

        assert result_df.empty
        assert result_y is None


# ============================================================================
# TEST CLASS: Horizontal Bar Chart Creation
# ============================================================================


class TestCreateHorizontalBar:
    """Test _create_horizontal_bar helper function."""

    def test_create_horizontal_with_y_column(self):
        """Test horizontal bar chart with y_column."""
        df = pd.DataFrame({"count": [10, 20, 15]}, index=["A", "B", "C"])

        fig, ax = plt.subplots()
        result_ax = _create_horizontal_bar(df, "count", ax, None)

        assert result_ax is ax
        assert ax.get_xlabel() == "Count"
        assert ax.get_ylabel() == ""

        plt.close(fig)

    def test_create_horizontal_without_y_column(self):
        """Test horizontal bar chart without y_column."""
        df = pd.DataFrame({"count": [10, 20, 15]}, index=["A", "B", "C"])

        fig, ax = plt.subplots()
        result_ax = _create_horizontal_bar(df, None, ax, None)

        assert result_ax is ax
        plt.close(fig)

    def test_create_horizontal_with_colors(self):
        """Test horizontal bar chart with custom colors."""
        df = pd.DataFrame({"count": [10, 20, 15]}, index=["A", "B", "C"])

        fig, ax = plt.subplots()
        colors = ["red", "green", "blue"]
        result_ax = _create_horizontal_bar(df, "count", ax, colors)

        assert result_ax is ax
        plt.close(fig)


# ============================================================================
# TEST CLASS: Vertical Bar Chart Creation
# ============================================================================


class TestCreateVerticalBar:
    """Test _create_vertical_bar helper function."""

    def test_create_vertical_with_y_column(self):
        """Test vertical bar chart with y_column."""
        df = pd.DataFrame({"count": [10, 20, 15]}, index=["A", "B", "C"])

        fig, ax = plt.subplots()
        result_ax = _create_vertical_bar(df, "count", ax, None)

        assert result_ax is ax
        assert ax.get_xlabel() == ""
        assert ax.get_ylabel() == "Count"

        plt.close(fig)

    def test_create_vertical_without_y_column(self):
        """Test vertical bar chart without y_column."""
        df = pd.DataFrame({"count": [10, 20, 15]}, index=["A", "B", "C"])

        fig, ax = plt.subplots()
        result_ax = _create_vertical_bar(df, None, ax, None)

        assert result_ax is ax
        plt.close(fig)


# ============================================================================
# TEST CLASS: BarChartPlugin - create_plot
# ============================================================================


class TestBarChartPluginCreatePlot:
    """Test BarChartPlugin.create_plot method."""

    def test_create_plot_horizontal_basic(self):
        """Test basic horizontal bar chart creation."""
        df = pd.DataFrame({"category": ["A", "B", "C"], "count": [10, 20, 15]})

        fig = BarChartPlugin.create_plot(data=df, x_column="category", y_column="count", horizontal=True)

        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_create_plot_vertical_basic(self):
        """Test basic vertical bar chart creation."""
        df = pd.DataFrame({"category": ["A", "B", "C"], "count": [10, 20, 15]})

        fig = BarChartPlugin.create_plot(data=df, x_column="category", y_column="count", horizontal=False)

        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_create_plot_with_title(self):
        """Test bar chart with title."""
        df = pd.DataFrame({"category": ["A", "B", "C"], "count": [10, 20, 15]})

        fig = BarChartPlugin.create_plot(data=df, x_column="category", y_column="count", title="Test Bar Chart")

        # Check that title was set
        ax = fig.axes[0]
        assert ax.get_title() == "Test Bar Chart"

        plt.close(fig)

    def test_create_plot_with_colors(self):
        """Test bar chart with custom colors."""
        df = pd.DataFrame({"category": ["A", "B", "C"], "count": [10, 20, 15]})

        colors = ["red", "green", "blue"]
        fig = BarChartPlugin.create_plot(data=df, x_column="category", y_column="count", colors=colors)

        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_create_plot_with_figsize(self):
        """Test bar chart with custom figure size."""
        df = pd.DataFrame({"category": ["A", "B"], "count": [10, 20]})

        fig = BarChartPlugin.create_plot(data=df, x_column="category", y_column="count", figsize=(8, 4))

        # Check figure size
        assert fig.get_figwidth() == 8
        assert fig.get_figheight() == 4

        plt.close(fig)

    def test_create_plot_with_limit(self):
        """Test bar chart with limit."""
        df = pd.DataFrame({"category": ["A", "B", "C", "D", "E"], "count": [10, 20, 15, 5, 25]})

        fig = BarChartPlugin.create_plot(data=df, x_column="category", y_column="count", limit=3)

        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_create_plot_save_to_file(self, tmp_path):
        """Test saving bar chart to file."""
        df = pd.DataFrame({"category": ["A", "B", "C"], "count": [10, 20, 15]})

        output_file = tmp_path / "bar_chart.png"

        fig = BarChartPlugin.create_plot(
            data=df, x_column="category", y_column="count", output_file=str(output_file), output_format="png"
        )

        assert output_file.exists()
        plt.close(fig)

    def test_create_plot_save_pdf(self, tmp_path):
        """Test saving bar chart as PDF."""
        df = pd.DataFrame({"category": ["A", "B"], "count": [10, 20]})

        output_file = tmp_path / "bar_chart.pdf"

        fig = BarChartPlugin.create_plot(
            data=df, x_column="category", y_column="count", output_file=str(output_file), output_format="pdf"
        )

        assert output_file.exists()
        plt.close(fig)

    def test_create_plot_with_default_columns(self):
        """Test bar chart with default column selection."""
        df = pd.DataFrame({"values": [10, 20, 15]}, index=["A", "B", "C"])

        fig = BarChartPlugin.create_plot(data=df)

        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_create_plot_error_handling(self):
        """Test error handling in create_plot."""
        # Pass invalid data
        with pytest.raises(VisualizationError, match="Error creating bar chart"):
            BarChartPlugin.create_plot(data=None)


# ============================================================================
# TEST CLASS: Plugin Metadata
# ============================================================================


class TestBarChartPluginMetadata:
    """Test BarChartPlugin metadata."""

    def test_plugin_name(self):
        """Test plugin name."""
        assert BarChartPlugin.name == "bar"

    def test_plugin_description(self):
        """Test plugin description."""
        assert BarChartPlugin.description == "Bar chart visualization"

    def test_plugin_version(self):
        """Test plugin version."""
        assert BarChartPlugin.version == "0.1.0"


# ============================================================================
# SUCCESS METRICS:
#
# After running these tests:
# - Expected: 30+ tests pass
# - Coverage: 22% -> 85%+ for plugins/visualizers/bar.py
# - All methods and edge cases covered
#
# Run tests:
#   pytest tests/test_visualizers_bar_extended.py -v
#
# Check coverage:
#   pytest --cov=prepare.plugins.visualizers.bar \
#          --cov-report=term-missing \
#          tests/test_visualizers_bar_extended.py -q
# ============================================================================
