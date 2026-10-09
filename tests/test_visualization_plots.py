"""
Test visualization.plots module functionality.

Tests for plotting functions including containment plots, metadata counts,
heatmaps, and correlation matrices with proper mocking of matplotlib.
"""

from unittest.mock import Mock, patch

import pytest
import pandas as pd

from prepare.visualization.plots import (
    _load_and_validate_data,
    _create_rank_plot,
    _create_histogram_plot,
    _create_box_plot,
    _create_violin_plot,
    _create_plot_by_type,
    _save_plot_if_needed,
    plot_containment,
    plot_metadata_counts,
    plot_correlation_matrix,
)
from prepare.core.exceptions import VisualizationError


class TestLoadAndValidateData:
    """Test _load_and_validate_data helper function."""

    def test_load_and_validate_data_success(self, tmp_path):
        """Test successful data loading and validation."""
        test_file = tmp_path / "test_data.tsv"
        test_data = pd.DataFrame(
            {
                "containment": [0.95, 0.87, 0.65, 0.12],
                "genome1": [0.95, 0.23, 0.15, 0.05],
                "genome2": [0.12, 0.87, 0.65, 0.12],
            },
            index=["Sample1", "Sample2", "Sample3", "Sample4"],
        )
        test_data.to_csv(test_file, sep="\t")

        result = _load_and_validate_data(test_file, "containment", threshold=None)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 4
        assert "containment" in result.columns
        pd.testing.assert_frame_equal(result, test_data)

    def test_load_and_validate_data_with_threshold(self, tmp_path):
        """Test data loading with threshold filtering."""
        test_file = tmp_path / "test_data.tsv"
        test_data = pd.DataFrame(
            {"containment": [0.95, 0.87, 0.65, 0.12], "genome1": [0.95, 0.23, 0.15, 0.05]},
            index=["Sample1", "Sample2", "Sample3", "Sample4"],
        )
        test_data.to_csv(test_file, sep="\t")

        result = _load_and_validate_data(test_file, "containment", threshold=0.8)

        assert len(result) == 2  # Only samples above 0.8
        assert result["containment"].min() >= 0.8

    def test_load_and_validate_data_threshold_no_data(self, tmp_path):
        """Test when threshold filters out all data."""
        test_file = tmp_path / "test_data.tsv"
        test_data = pd.DataFrame(
            {"containment": [0.15, 0.12, 0.08, 0.05]}, index=["Sample1", "Sample2", "Sample3", "Sample4"]
        )
        test_data.to_csv(test_file, sep="\t")

        with patch("prepare.visualization.plots.logger") as mock_logger:
            result = _load_and_validate_data(test_file, "containment", threshold=0.8)

        assert result.empty
        mock_logger.warning.assert_called_once_with("No data above threshold 0.8")

    def test_load_and_validate_data_missing_column(self, tmp_path):
        """Test error when specified column doesn't exist."""
        test_file = tmp_path / "test_data.tsv"
        test_data = pd.DataFrame(
            {"genome1": [0.95, 0.87, 0.65], "genome2": [0.12, 0.23, 0.15]}, index=["Sample1", "Sample2", "Sample3"]
        )
        test_data.to_csv(test_file, sep="\t")

        with pytest.raises(VisualizationError, match="Column 'missing_column' not found"):
            _load_and_validate_data(test_file, "missing_column", threshold=None)

    def test_load_and_validate_data_file_not_found(self):
        """Test error when file doesn't exist."""
        with pytest.raises(FileNotFoundError):
            _load_and_validate_data("nonexistent_file.tsv", "containment", threshold=None)


class TestCreateRankPlot:
    """Test _create_rank_plot helper function."""

    def test_create_rank_plot_success(self):
        """Test successful rank plot creation."""
        mock_ax = Mock()
        test_data = pd.DataFrame(
            {"containment": [0.95, 0.87, 0.65, 0.12]}, index=["Sample1", "Sample2", "Sample3", "Sample4"]
        )

        _create_rank_plot(mock_ax, test_data, "containment", "blue")

        # Verify scatter plot was called
        mock_ax.scatter.assert_called_once()
        call_args = mock_ax.scatter.call_args

        # Verify x-axis is ranks (1, 2, 3, 4)
        x_values = call_args[0][0].tolist()
        assert x_values == [1, 2, 3, 4]

        # Verify y-axis is sorted containment values (descending)
        y_values = call_args[0][1].tolist()
        assert y_values == [0.95, 0.87, 0.65, 0.12]

        # Verify color parameter
        assert call_args[1]["color"] == "blue"

    def test_create_rank_plot_empty_data(self):
        """Test rank plot with empty data."""
        mock_ax = Mock()
        empty_data = pd.DataFrame()

        _create_rank_plot(mock_ax, empty_data, "containment", "red")

        # Should still call scatter with empty arrays
        mock_ax.scatter.assert_called_once()


class TestCreateHistogramPlot:
    """Test _create_histogram_plot helper function."""

    def test_create_histogram_plot_success(self):
        """Test successful histogram creation."""
        mock_ax = Mock()
        test_data = pd.DataFrame({"containment": [0.95, 0.87, 0.65, 0.12, 0.95, 0.87]})

        _create_histogram_plot(mock_ax, test_data, "containment", "green")

        mock_ax.hist.assert_called_once()
        call_args = mock_ax.hist.call_args

        # Verify data passed to histogram
        data_values = call_args[0][0].tolist()
        assert set(data_values) == {0.95, 0.87, 0.65, 0.12, 0.95, 0.87}

        # Verify color and bin parameters
        assert call_args[1]["color"] == "green"
        assert call_args[1]["bins"] == 5  # Adjusted for 6 data points
        assert call_args[1]["alpha"] == 0.7

    def test_create_histogram_plot_small_dataset(self):
        """Test histogram with small dataset (adjusts bin count)."""
        mock_ax = Mock()
        test_data = pd.DataFrame({"containment": [0.95, 0.87]})  # Only 2 samples

        _create_histogram_plot(mock_ax, test_data, "containment", "red")

        call_args = mock_ax.hist.call_args
        assert call_args[1]["bins"] == 5  # Should use min(20, 5) = 5 bins


class TestCreateBoxPlot:
    """Test _create_box_plot helper function."""

    def test_create_box_plot_success(self):
        """Test successful box plot creation."""
        mock_ax = Mock()
        test_data = pd.DataFrame({"containment": [0.95, 0.87, 0.65, 0.12, 0.23, 0.45]})

        _create_box_plot(mock_ax, test_data, "containment")

        mock_ax.boxplot.assert_called_once()
        call_args = mock_ax.boxplot.call_args

        # Verify data passed to boxplot
        data_values = call_args[0][0].tolist()
        assert set(data_values) == {0.95, 0.87, 0.65, 0.12, 0.23, 0.45}


class TestCreateViolinPlot:
    """Test _create_violin_plot helper function."""

    def test_create_violin_plot_success(self):
        """Test successful violin plot creation."""
        mock_ax = Mock()
        test_data = pd.DataFrame({"containment": [0.95, 0.87, 0.65, 0.12, 0.23, 0.45, 0.67, 0.89]})

        _create_violin_plot(mock_ax, test_data, "containment")

        mock_ax.violinplot.assert_called_once()
        call_args = mock_ax.violinplot.call_args

        # Verify data passed to violin plot
        data_values = call_args[0][0].tolist()
        assert set(data_values) == {0.95, 0.87, 0.65, 0.12, 0.23, 0.45, 0.67, 0.89}


class TestCreatePlotByType:
    """Test _create_plot_by_type dispatcher function."""

    def test_create_plot_by_type_rank(self):
        """Test plot creation with rank type."""
        mock_ax = Mock()
        test_data = pd.DataFrame({"containment": [0.95, 0.87, 0.65]})

        with patch("prepare.visualization.plots._create_rank_plot") as mock_rank:
            _create_plot_by_type(mock_ax, test_data, "containment", "blue", "rank")

        mock_rank.assert_called_once_with(mock_ax, test_data, "containment", "blue")

    def test_create_plot_by_type_histogram(self):
        """Test plot creation with histogram type."""
        mock_ax = Mock()
        test_data = pd.DataFrame({"containment": [0.95, 0.87, 0.65]})

        with patch("prepare.visualization.plots._create_histogram_plot") as mock_hist:
            _create_plot_by_type(mock_ax, test_data, "containment", "red", "histogram")

        mock_hist.assert_called_once_with(mock_ax, test_data, "containment", "red")

    def test_create_plot_by_type_box(self):
        """Test plot creation with box type."""
        mock_ax = Mock()
        test_data = pd.DataFrame({"containment": [0.95, 0.87, 0.65]})

        with patch("prepare.visualization.plots._create_box_plot") as mock_box:
            _create_plot_by_type(mock_ax, test_data, "containment", "green", "box")

        mock_box.assert_called_once_with(mock_ax, test_data, "containment")

    def test_create_plot_by_type_violin(self):
        """Test plot creation with violin type."""
        mock_ax = Mock()
        test_data = pd.DataFrame({"containment": [0.95, 0.87, 0.65, 0.12, 0.23, 0.45]})

        with patch("prepare.visualization.plots._create_violin_plot") as mock_violin:
            _create_plot_by_type(mock_ax, test_data, "containment", "purple", "violin")

        mock_violin.assert_called_once_with(mock_ax, test_data, "containment")

    def test_create_plot_by_type_unknown(self):
        """Test error with unknown plot type."""
        mock_ax = Mock()
        test_data = pd.DataFrame({"containment": [0.95, 0.87]})

        with pytest.raises(VisualizationError, match="Unknown plot type"):
            _create_plot_by_type(mock_ax, test_data, "containment", "blue", "unknown")


class TestSavePlotIfNeeded:
    """Test _save_plot_if_needed helper function."""

    def test_save_plot_if_needed_with_format(self, tmp_path):
        """Test plot saving when format is specified."""
        output_file = tmp_path / "test_plot"

        with patch("matplotlib.pyplot.savefig") as mock_savefig:
            result = _save_plot_if_needed(output_file, "containment", "max_containment", "png")

        expected_path = tmp_path / "test_plot_containment_max_containment.png"
        assert result == expected_path
        mock_savefig.assert_called_once_with(expected_path, dpi=300, bbox_inches="tight")

    def test_save_plot_default_png_and_stem_name(self, tmp_path):
        """No format given still saves, defaulting to png and using the file's stem (not its
        full name with extension) so a table named parsed_containment.txt does not produce
        parsed_containment.txt_rank_max_containment.png."""
        import matplotlib.pyplot as plt

        plt.figure()
        out = _save_plot_if_needed(tmp_path / "parsed_containment.txt", "rank", "max_containment", None)
        assert out == tmp_path / "parsed_containment_rank_max_containment.png"
        assert out.exists()

    def test_save_plot_if_needed_different_formats(self, tmp_path):
        """Test plot saving with different formats."""
        output_file = tmp_path / "test_plot"

        for format_type in ["png", "pdf", "svg"]:
            with patch("matplotlib.pyplot.savefig") as mock_savefig:
                result = _save_plot_if_needed(output_file, "containment", "max_containment", format_type)

            expected_path = tmp_path / f"test_plot_containment_max_containment.{format_type}"
            assert result == expected_path
            mock_savefig.assert_called_once_with(expected_path, dpi=300, bbox_inches="tight")


class TestPlotContainment:
    """Test plot_containment public API function."""

    def test_plot_containment_success(self, tmp_path):
        """Test successful containment plot creation."""
        test_file = tmp_path / "containment_data.tsv"
        test_data = pd.DataFrame(
            {
                "max_containment": [0.95, 0.87, 0.65, 0.05],  # 0.05 is below 0.1 threshold
                "sample_count": [100, 85, 65, 12],
            },
            index=["Sample1", "Sample2", "Sample3", "Sample4"],
        )
        test_data.to_csv(test_file, sep="\t")

        mock_fig = Mock()
        mock_ax = Mock()

        with patch("matplotlib.pyplot.subplots", return_value=(mock_fig, mock_ax)):
            with patch("prepare.visualization.plots._create_plot_by_type") as mock_create:
                with patch("matplotlib.pyplot.show"):
                    result = plot_containment(
                        file_path=str(test_file),
                        column="max_containment",
                        plot_type="rank",
                        threshold=0.1,
                        colors="blue",
                    )

        assert result == mock_fig
        mock_create.assert_called_once()

        # Verify correct parameters passed
        call_args = mock_create.call_args
        assert call_args[0][1].shape[0] == 3  # 3 samples above 0.1 threshold
        assert call_args[0][2] == "max_containment"
        assert call_args[0][3] == "blue"
        assert call_args[0][4] == "rank"

    def test_plot_containment_all_plot_types(self, tmp_path):
        """Test all supported plot types."""
        test_file = tmp_path / "containment_data.tsv"
        test_data = pd.DataFrame(
            {"max_containment": [0.95, 0.87, 0.65, 0.12, 0.23, 0.45]}, index=[f"Sample{i}" for i in range(1, 7)]
        )
        test_data.to_csv(test_file, sep="\t")

        plot_types = ["rank", "histogram", "box", "violin"]

        for plot_type in plot_types:
            with patch("matplotlib.pyplot.subplots", return_value=(Mock(), Mock())):
                with patch("prepare.visualization.plots._create_plot_by_type") as mock_create:
                    with patch("matplotlib.pyplot.show"):
                        plot_containment(file_path=str(test_file), column="max_containment", plot_type=plot_type)

                mock_create.assert_called_once()

    def test_plot_containment_with_save(self, tmp_path):
        """Test containment plot with file saving."""
        test_file = tmp_path / "containment_data.tsv"
        tmp_path / "output_plot"
        test_data = pd.DataFrame({"max_containment": [0.95, 0.87, 0.65]}, index=["Sample1", "Sample2", "Sample3"])
        test_data.to_csv(test_file, sep="\t")

        with patch("matplotlib.pyplot.subplots", return_value=(Mock(), Mock())):
            with patch("prepare.visualization.plots._save_plot_if_needed") as mock_save:
                with patch("matplotlib.pyplot.show"):
                    plot_containment(file_path=str(test_file), column="max_containment", save_format="png")

        mock_save.assert_called_once_with(str(test_file), "rank", "max_containment", "png")

    def test_plot_containment_invalid_plot_type(self, tmp_path):
        """Test error with invalid plot type."""
        test_file = tmp_path / "containment_data.tsv"
        test_data = pd.DataFrame({"max_containment": [0.95, 0.87]})
        test_data.to_csv(test_file, sep="\t")

        with pytest.raises(VisualizationError, match="Unknown plot type"):
            plot_containment(str(test_file), "max_containment", plot_type="invalid")

    def test_plot_containment_empty_data_after_threshold(self, tmp_path):
        """Test containment plot when threshold filters out all data."""
        test_file = tmp_path / "containment_data.tsv"
        test_data = pd.DataFrame({"max_containment": [0.05, 0.03, 0.01]}, index=["Sample1", "Sample2", "Sample3"])
        test_data.to_csv(test_file, sep="\t")

        with patch("prepare.visualization.plots.logger") as mock_logger:
            with patch("matplotlib.pyplot.subplots", return_value=(Mock(), Mock())):
                with patch("matplotlib.pyplot.show"):
                    plot_containment(str(test_file), "max_containment", threshold=0.8)

        mock_logger.warning.assert_called()


class TestPlotMetadataCounts:
    """Test plot_metadata_counts public API function."""

    def test_plot_metadata_counts_success(self, tmp_path):
        """Test successful metadata counts plotting."""
        test_file = tmp_path / "metadata_counts.tsv"
        test_data = pd.DataFrame(
            {"category": ["E. coli", "S. enterica", "B. subtilis", "P. aeruginosa"], "count": [150, 120, 85, 60]}
        )
        test_data.to_csv(test_file, sep="\t", index=False, header=False)

        mock_plugin = Mock()
        mock_plugin.create_plot.return_value = Mock()

        with patch("prepare.plugins.base.visualizer_registry.get") as mock_get_plugin:
            mock_get_plugin.return_value = mock_plugin

            result = plot_metadata_counts(file_path=str(test_file), plot_type="bar")

        assert result is not None
        mock_plugin.create_plot.assert_called_once()

    def test_plot_metadata_counts_different_chart_types(self, tmp_path):
        """Test metadata counts with different chart types."""
        test_file = tmp_path / "metadata_counts.tsv"
        test_data = pd.DataFrame({"organism": ["E. coli", "S. enterica"], "count": [100, 80]})
        test_data.to_csv(test_file, sep="\t", index=False, header=False)

        # Test bar chart (uses plugin)
        mock_plugin = Mock()
        mock_plugin.create_plot.return_value = Mock()

        with patch("prepare.plugins.base.visualizer_registry.get") as mock_get_plugin:
            mock_get_plugin.return_value = mock_plugin

            result = plot_metadata_counts(file_path=str(test_file), plot_type="bar")

        mock_plugin.create_plot.assert_called_once()

        # Test pie chart (created directly with matplotlib)
        result = plot_metadata_counts(file_path=str(test_file), plot_type="pie")
        assert result is not None

    def test_plot_metadata_counts_radar(self, tmp_path):
        """Radar chart renders directly for >= 3 categories."""
        test_file = tmp_path / "metadata_counts.tsv"
        pd.DataFrame({"organism": ["E. coli", "S. enterica", "K. pneumoniae"], "count": [100, 80, 60]}).to_csv(
            test_file, sep="\t", index=False, header=False
        )

        result = plot_metadata_counts(file_path=str(test_file), plot_type="radar", colors="green")
        assert result is not None

    def test_plot_metadata_counts_radar_falls_back_to_bar(self, tmp_path):
        """Fewer than 3 categories falls back to the bar plugin."""
        test_file = tmp_path / "metadata_counts.tsv"
        pd.DataFrame({"organism": ["E. coli", "S. enterica"], "count": [100, 80]}).to_csv(
            test_file, sep="\t", index=False, header=False
        )

        mock_plugin = Mock()
        mock_plugin.create_plot.return_value = Mock()
        with patch("prepare.plugins.base.visualizer_registry.get", return_value=mock_plugin):
            result = plot_metadata_counts(file_path=str(test_file), plot_type="radar")

        assert result is not None
        mock_plugin.create_plot.assert_called_once()

    def test_plot_metadata_counts_with_limit(self, tmp_path):
        """Test metadata counts with top N limit."""
        test_file = tmp_path / "metadata_counts.tsv"
        test_data = pd.DataFrame(
            {
                "organism": ["E. coli", "S. enterica", "B. subtilis", "P. aeruginosa", "K. pneumoniae"],
                "count": [150, 120, 85, 60, 40],
            }
        )
        test_data.to_csv(test_file, sep="\t", index=False, header=False)

        mock_plugin = Mock()
        mock_plugin.create_plot.return_value = Mock()

        with patch("prepare.plugins.base.visualizer_registry.get") as mock_get_plugin:
            mock_get_plugin.return_value = mock_plugin

            plot_metadata_counts(file_path=str(test_file), plot_type="bar", limit=3)

        # Verify only top 3 entries were passed to the plugin
        call_args = mock_plugin.create_plot.call_args
        data_passed = call_args.kwargs["data"]  # Data should be passed as keyword argument
        assert len(data_passed) == 3

    def test_plot_metadata_counts_missing_columns(self, tmp_path):
        """Test error when data contains no numeric counts at all.

        Every count cell is non-numeric text, so after the header-detection
        heuristic strips the first row there are still no numeric values left;
        this must raise rather than silently plotting fabricated zero counts.
        """
        test_file = tmp_path / "metadata_counts.tsv"
        test_data = pd.DataFrame(
            {
                "organism": ["E. coli", "S. enterica"],
                "wrong_data": ["text1", "text2"],  # Non-numeric data in count column
            }
        )
        test_data.to_csv(test_file, sep="\t", index=False, header=False)

        with pytest.raises(VisualizationError, match="No numeric count values"):
            plot_metadata_counts(file_path=str(test_file), plot_type="bar")

    def test_plot_metadata_counts_plugin_not_found(self, tmp_path):
        """Test error when chart plugin not found."""
        test_file = tmp_path / "metadata_counts.tsv"
        test_data = pd.DataFrame({"category": ["E. coli"], "count": [100]})
        test_data.to_csv(test_file, sep="\t", index=False, header=False)

        with patch("prepare.plugins.base.visualizer_registry.get", return_value=None):
            with pytest.raises(VisualizationError, match="Unknown plot type"):
                plot_metadata_counts(file_path=str(test_file), plot_type="unknown_type")


class TestPlotCorrelationMatrix:
    """Test plot_correlation_matrix public API function."""

    def test_plot_correlation_matrix_success(self, tmp_path):
        """Test successful correlation matrix creation."""
        test_data = pd.DataFrame(
            {
                "genome1": [0.95, 0.12, 0.23, 0.45],
                "genome2": [0.23, 0.87, 0.15, 0.67],
                "genome3": [0.15, 0.45, 0.67, 0.23],
                "sample_count": [100, 85, 65, 45],  # Non-genome column
            },
            index=["sample1", "sample2", "sample3", "sample4"],
        )

        mock_plugin = Mock()
        mock_fig = Mock()
        mock_plugin.create_correlation_heatmap.return_value = mock_fig

        with patch("prepare.plugins.base.visualizer_registry.get") as mock_get_plugin:
            mock_get_plugin.return_value = mock_plugin

            result = plot_correlation_matrix(test_data)

        assert result is not None
        mock_plugin.create_correlation_heatmap.assert_called_once()

    def test_plot_correlation_matrix_insufficient_genomes(self):
        """Test correlation matrix with mixed data types."""
        test_data = pd.DataFrame(
            {"genome1": [0.95, 0.12, 0.23], "sample_count": [100, 85, 65]},  # Non-genome column
            index=["sample1", "sample2", "sample3"],
        )

        # Should still work with mixed numeric data
        mock_plugin = Mock()
        mock_fig = Mock()
        mock_plugin.create_correlation_heatmap.return_value = mock_fig

        with patch("prepare.plugins.base.visualizer_registry.get") as mock_get_plugin:
            mock_get_plugin.return_value = mock_plugin

            result = plot_correlation_matrix(test_data)

        assert result is not None

    def test_plot_correlation_matrix_no_genome_columns(self):
        """Test error when non-numeric data prevents correlation."""
        test_data = pd.DataFrame(
            {"sample_count": [100, 85, 65], "metadata_field": ["A", "B", "C"]}, index=["sample1", "sample2", "sample3"]
        )

        # This should raise an error because correlation can't be computed on strings
        with pytest.raises(VisualizationError, match="could not convert string to float"):
            plot_correlation_matrix(test_data)

    def test_plot_correlation_matrix_plugin_not_found(self):
        """Test error when heatmap plugin not available."""
        test_data = pd.DataFrame({"GCF_001": [0.95, 0.12], "GCF_002": [0.23, 0.87]}, index=["sample1", "sample2"])

        with patch("prepare.plugins.base.visualizer_registry.get", return_value=None):
            with pytest.raises(VisualizationError, match="'NoneType' object has no attribute"):
                plot_correlation_matrix(test_data)


class TestPlotsIntegration:
    """Integration tests for plotting workflows."""

    def test_full_plotting_workflow(self, tmp_path):
        """Test complete plotting workflow with realistic data."""
        # Create realistic containment data
        containment_file = tmp_path / "containment_results.tsv"
        containment_data = pd.DataFrame(
            {
                "max_containment": [0.95, 0.87, 0.65, 0.45, 0.23, 0.12, 0.05],
                "max_containment_annotation": [
                    "Escherichia coli",
                    "Salmonella enterica",
                    "Klebsiella pneumoniae",
                    "Pseudomonas aeruginosa",
                    "Bacillus subtilis",
                    "Streptococcus pneumoniae",
                    "Enterococcus faecalis",
                ],
                "GCF_000005825.2": [0.95, 0.12, 0.05, 0.02, 0.01, 0.0, 0.0],  # E. coli
                "GCF_000006945.2": [0.23, 0.87, 0.15, 0.08, 0.03, 0.01, 0.0],  # Salmonella
                "GCF_000009605.1": [0.15, 0.25, 0.65, 0.12, 0.05, 0.02, 0.01],  # Klebsiella
            },
            index=[f"ERR{123456+i}" for i in range(7)],
        )
        containment_data.to_csv(containment_file, sep="\t")

        # Create metadata counts data
        counts_file = tmp_path / "organism_counts.tsv"
        counts_data = pd.DataFrame(
            {
                "organism": ["E. coli", "Salmonella enterica", "Klebsiella pneumoniae", "Other"],
                "count": [245, 189, 123, 87],
            }
        )
        counts_data.to_csv(counts_file, sep="\t", index=False)

        # Test containment plotting
        with patch("matplotlib.pyplot.subplots", return_value=(Mock(), Mock())):
            with patch("matplotlib.pyplot.show"):
                fig1 = plot_containment(
                    file_path=str(containment_file),
                    column="max_containment",
                    plot_type="rank",
                    threshold=0.1,
                    colors="viridis",
                )

        assert fig1 is not None

        # Test metadata counts plotting
        mock_bar_plugin = Mock()
        mock_bar_plugin.create_plot.return_value = Mock()

        with patch("prepare.plugins.base.visualizer_registry.get") as mock_get_plugin:
            mock_get_plugin.return_value = mock_bar_plugin

            fig2 = plot_metadata_counts(file_path=str(counts_file), plot_type="bar", limit=10)

        assert fig2 is not None

        # Test correlation matrix
        mock_heatmap_plugin = Mock()
        mock_heatmap_plugin.create_heatmap.return_value = Mock()

        with patch("prepare.plugins.base.visualizer_registry.get") as mock_get_plugin:
            mock_get_plugin.return_value = mock_heatmap_plugin

            fig3 = plot_correlation_matrix(containment_data)

        assert fig3 is not None

    def test_error_handling_comprehensive(self, tmp_path):
        """Test comprehensive error handling."""
        # Test file not found
        with pytest.raises(VisualizationError, match="No such file or directory"):
            plot_containment("nonexistent_file.tsv", "max_containment")

        # Test invalid data format
        invalid_file = tmp_path / "invalid_data.tsv"
        invalid_file.write_text("invalid,data\nformat,here")

        with pytest.raises(VisualizationError):
            plot_containment(str(invalid_file), "nonexistent_column")


class TestCountsLoader:
    def test_header_row_and_genome_columns_are_summed(self, tmp_path):
        from prepare.visualization.plots import _load_counts_df

        f = tmp_path / "metadata_counts.txt"
        f.write_text("Sample_Scientific_Name\tGCF_A\tGCF_B\nmetagenome\t3\t0\nfood metagenome\t0\t2\n")
        df = _load_counts_df(f, limit=20)
        assert list(df["category"]) == ["metagenome", "food metagenome"]
        assert list(df["count"]) == [3, 2]

    def test_headerless_two_columns(self, tmp_path):
        from prepare.visualization.plots import _load_counts_df

        f = tmp_path / "stats.txt"
        f.write_text("GCF_B\t7\nGCF_A\t9\n")
        df = _load_counts_df(f, limit=20)
        assert list(df["category"]) == ["GCF_A", "GCF_B"]
        assert list(df["count"]) == [9, 7]

    def test_bar_plot_saves_file(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        f = tmp_path / "counts.txt"
        f.write_text("a\t3\nb\t1\n")
        plot_metadata_counts(file_path=str(f), plot_type="bar", save_format="png")
        assert (tmp_path / "counts_bar.png").exists()


if __name__ == "__main__":
    pytest.main([__file__])
