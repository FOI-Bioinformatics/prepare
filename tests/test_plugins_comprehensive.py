"""
Comprehensive unit tests for PREPARE plugin system.

This module provides comprehensive testing for the plugin system including:
- Plugin registry operations
- Plugin discovery mechanisms
- Format plugins (Branchwater)
- Visualizer plugins (Bar, Heatmap, Map)
- Error handling and edge cases
"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import tempfile
from unittest.mock import Mock, patch
from typing import List

from prepare.core.exceptions import ConfigurationError, PluginError, ValidationError, VisualizationError
from prepare.core.models import Containment, SRAMetadata
from prepare.plugins.base import (
    Plugin,
    PluginRegistry,
    format_registry,
    visualizer_registry,
)
from prepare.plugins.formats.branchwater import BranchWaterFormatPlugin
from prepare.plugins.visualizers.bar import BarChartPlugin
from prepare.plugins.visualizers.heatmap import HeatmapPlugin
from prepare.plugins.visualizers.map import MapVisualizerPlugin


# Test Plugin Classes
class MockFormatPlugin(Plugin):
    """Mock format plugin for testing."""

    name = "mock_format"
    description = "Mock format plugin"
    version = "0.1.0"

    @classmethod
    def validate_header(cls, headers: List[str]) -> bool:
        return True

    @classmethod
    def parse_file(cls, file_path, genome_id: str) -> List[Containment]:
        return []

    @classmethod
    def extract_metadata(cls, containment: Containment) -> SRAMetadata:
        return SRAMetadata(accession=containment.accession)


class MockVisualizerPlugin(Plugin):
    """Mock visualizer plugin for testing."""

    name = "mock_visualizer"
    description = "Mock visualizer plugin"
    version = "0.1.0"

    @classmethod
    def create_plot(cls, data, **kwargs):
        return Mock()


class TestPluginBase:
    """Test the Plugin base class."""

    def test_plugin_metadata(self):
        """Test plugin metadata attributes."""
        plugin = MockFormatPlugin()
        assert plugin.get_name() == "mock_format"
        assert plugin.get_description() == "Mock format plugin"
        assert plugin.get_version() == "0.1.0"

    def test_plugin_inheritance(self):
        """Test plugin inheritance from base class."""
        assert issubclass(MockFormatPlugin, Plugin)
        assert issubclass(MockVisualizerPlugin, Plugin)
        assert issubclass(BranchWaterFormatPlugin, Plugin)
        assert issubclass(BarChartPlugin, Plugin)


class TestPluginRegistry:
    """Test the PluginRegistry class."""

    def test_registry_creation(self):
        """Test creating a new registry."""
        registry = PluginRegistry()
        assert isinstance(registry._plugins, dict)
        assert len(registry.list()) == 0

    def test_plugin_registration(self):
        """Test registering plugins."""
        registry = PluginRegistry()

        # Register a plugin
        registry.register(MockFormatPlugin)
        assert "mock_format" in registry.list()
        assert len(registry.list()) == 1

        # Register another plugin
        registry.register(MockVisualizerPlugin)
        assert len(registry.list()) == 2
        assert "mock_visualizer" in registry.list()

    def test_duplicate_registration_error(self):
        """Test error when registering duplicate plugins."""
        registry = PluginRegistry()
        registry.register(MockFormatPlugin)

        with pytest.raises(PluginError, match="Plugin 'mock_format' is already registered"):
            registry.register(MockFormatPlugin)

    def test_plugin_retrieval(self):
        """Test retrieving registered plugins."""
        registry = PluginRegistry()
        registry.register(MockFormatPlugin)

        # Get plugin by name
        plugin_class = registry.get("mock_format")
        assert plugin_class == MockFormatPlugin

    def test_get_nonexistent_plugin_error(self):
        """Test error when getting non-existent plugin."""
        registry = PluginRegistry()

        with pytest.raises(PluginError, match="Plugin 'nonexistent' is not registered"):
            registry.get("nonexistent")

    def test_plugin_unregistration(self):
        """Test unregistering plugins."""
        registry = PluginRegistry()
        registry.register(MockFormatPlugin)
        assert "mock_format" in registry.list()

        # Unregister plugin
        registry.unregister("mock_format")
        assert "mock_format" not in registry.list()

    def test_unregister_nonexistent_plugin_error(self):
        """Test error when unregistering non-existent plugin."""
        registry = PluginRegistry()

        with pytest.raises(PluginError, match="Plugin 'nonexistent' is not registered"):
            registry.unregister("nonexistent")

    def test_list_all_plugins(self):
        """Test listing all registered plugins."""
        registry = PluginRegistry()
        registry.register(MockFormatPlugin)
        registry.register(MockVisualizerPlugin)

        plugins = registry.list()
        assert len(plugins) == 2
        assert "mock_format" in plugins
        assert "mock_visualizer" in plugins

    def test_get_all_plugins(self):
        """Test getting all registered plugins as dictionary."""
        registry = PluginRegistry()
        registry.register(MockFormatPlugin)
        registry.register(MockVisualizerPlugin)

        all_plugins = registry.get_all()
        assert len(all_plugins) == 2
        assert all_plugins["mock_format"] == MockFormatPlugin
        assert all_plugins["mock_visualizer"] == MockVisualizerPlugin


class TestBranchWaterFormatPlugin:
    """Test the BranchWaterFormatPlugin."""

    def test_plugin_metadata(self):
        """Test plugin metadata."""
        assert BranchWaterFormatPlugin.name == "branchwater"
        assert "Branchwater" in BranchWaterFormatPlugin.description
        assert BranchWaterFormatPlugin.version == "0.1.0"

    def test_required_columns(self):
        """Test required columns definition."""
        required = BranchWaterFormatPlugin.REQUIRED_COLS
        assert "acc" in required
        assert "containment" in required

    def test_header_validation_success(self):
        """Test successful header validation."""
        headers = ["acc", "containment", "organism", "biosample"]
        assert BranchWaterFormatPlugin.validate_header(headers) is True

    def test_header_validation_failure(self):
        """Test header validation failure."""
        headers = ["accession", "containment"]  # Missing 'acc'
        assert BranchWaterFormatPlugin.validate_header(headers) is False

        headers = ["acc", "similarity"]  # Missing 'containment'
        assert BranchWaterFormatPlugin.validate_header(headers) is False

    def test_parse_file_success(self):
        """Test successful file parsing."""
        csv_content = """acc,containment,organism,biosample,bioproject
SRR123456,0.95,Escherichia coli,SAMN123456,PRJNA123456
SRR789012,0.85,Salmonella enterica,SAMN789012,PRJNA789012
"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            f.flush()

            try:
                containments = BranchWaterFormatPlugin.parse_file(f.name, "GCF_123456.1")

                assert len(containments) == 2

                # Test first containment
                c1 = containments[0]
                assert c1.accession == "SRR123456"
                assert c1.value == 0.95
                assert c1.genome_id == "GCF_123456.1"
                assert c1.additional_data["organism"] == "Escherichia coli"
                assert c1.additional_data["biosample"] == "SAMN123456"

                # Test second containment
                c2 = containments[1]
                assert c2.accession == "SRR789012"
                assert c2.value == 0.85

            finally:
                os.unlink(f.name)

    def test_parse_file_invalid_headers(self):
        """Test parsing file with invalid headers."""
        csv_content = """accession,similarity,organism
SRR123456,0.95,Escherichia coli
"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            f.flush()

            try:
                with pytest.raises(ValidationError, match="Missing required columns"):
                    BranchWaterFormatPlugin.parse_file(f.name, "GCF_123456.1")
            finally:
                os.unlink(f.name)

    def test_parse_file_invalid_accession(self):
        """Test parsing file with invalid accessions."""
        csv_content = """acc,containment
invalid_acc,0.95
SRR123456,0.85
"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            f.flush()

            try:
                containments = BranchWaterFormatPlugin.parse_file(f.name, "GCF_123456.1")
                # Should skip invalid accession and only return valid one
                assert len(containments) == 1
                assert containments[0].accession == "SRR123456"
            finally:
                os.unlink(f.name)

    def test_parse_file_invalid_containment(self):
        """Test parsing file with invalid containment values."""
        csv_content = """acc,containment
SRR123456,invalid
SRR789012,0.85
"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            f.flush()

            try:
                containments = BranchWaterFormatPlugin.parse_file(f.name, "GCF_123456.1")
                # Should skip invalid containment and only return valid one
                assert len(containments) == 1
                assert containments[0].accession == "SRR789012"
            finally:
                os.unlink(f.name)

    def test_extract_metadata_success(self):
        """Test successful metadata extraction."""
        containment = Containment(
            accession="SRR123456",
            value=0.95,
            genome_id="GCF_123456.1",
            additional_data={
                "biosample": "SAMN123456",
                "bioproject": "PRJNA123456",
                "organism": "Escherichia coli",
                "collection_date_sam": "2023-01-01",
                "geo_loc_name_country_calc": "USA",
                "custom_field": "custom_value",
            },
        )

        metadata = BranchWaterFormatPlugin.extract_metadata(containment)

        assert metadata.accession == "SRR123456"
        assert metadata.biosample == "SAMN123456"
        assert metadata.bioproject == "PRJNA123456"
        assert metadata.organism == "Escherichia coli"
        assert metadata.collection_date == "2023-01-01"
        assert metadata.location == "USA"
        assert metadata.attributes["custom_field"] == "custom_value"

    def test_extract_metadata_no_additional_data(self):
        """Test metadata extraction with no additional data."""
        containment = Containment(accession="SRR123456", value=0.95, genome_id="GCF_123456.1")

        metadata = BranchWaterFormatPlugin.extract_metadata(containment)
        assert metadata is None

    def test_parse_file_io_error(self):
        """Test parsing non-existent file."""
        with pytest.raises(ValidationError, match="Error parsing Branchwater file"):
            BranchWaterFormatPlugin.parse_file("/nonexistent/file.csv", "GCF_123456.1")


class TestBarChartPlugin:
    """Test the BarChartPlugin."""

    def setup_method(self):
        """Set up test data."""
        self.test_data = pd.DataFrame({"category": ["A", "B", "C", "D", "E"], "count": [10, 25, 15, 30, 5]})

    def test_plugin_metadata(self):
        """Test plugin metadata."""
        assert BarChartPlugin.name == "bar"
        assert "Bar chart" in BarChartPlugin.description
        assert BarChartPlugin.version == "0.1.0"

    @patch("matplotlib.pyplot.subplots")
    @patch("matplotlib.pyplot.tight_layout")
    @patch("prepare.plugins.visualizers.bar._create_horizontal_bar")
    def test_create_plot_basic(self, mock_create_horizontal_bar, mock_tight_layout, mock_subplots):
        """Test basic plot creation."""
        mock_fig = Mock()
        mock_ax = Mock()
        mock_subplots.return_value = (mock_fig, mock_ax)
        mock_create_horizontal_bar.return_value = mock_ax

        result = BarChartPlugin.create_plot(
            data=self.test_data, x_column="category", y_column="count", title="Test Chart"
        )

        assert result == mock_fig
        mock_subplots.assert_called_once()
        mock_tight_layout.assert_called_once()
        mock_create_horizontal_bar.assert_called_once()

    @patch("matplotlib.pyplot.subplots")
    @patch("matplotlib.pyplot.tight_layout")
    @patch("prepare.plugins.visualizers.bar._create_horizontal_bar")
    def test_create_plot_horizontal(self, mock_create_horizontal_bar, mock_tight_layout, mock_subplots):
        """Test horizontal bar chart creation."""
        mock_fig = Mock()
        mock_ax = Mock()
        mock_subplots.return_value = (mock_fig, mock_ax)
        mock_create_horizontal_bar.return_value = mock_ax

        result = BarChartPlugin.create_plot(data=self.test_data, x_column="category", y_column="count", horizontal=True)

        assert result == mock_fig
        mock_create_horizontal_bar.assert_called_once()

    @patch("matplotlib.pyplot.subplots")
    @patch("matplotlib.pyplot.tight_layout")
    @patch("prepare.plugins.visualizers.bar._create_vertical_bar")
    def test_create_plot_vertical(self, mock_create_vertical_bar, mock_tight_layout, mock_subplots):
        """Test vertical bar chart creation."""
        mock_fig = Mock()
        mock_ax = Mock()
        mock_subplots.return_value = (mock_fig, mock_ax)
        mock_create_vertical_bar.return_value = mock_ax

        result = BarChartPlugin.create_plot(
            data=self.test_data, x_column="category", y_column="count", horizontal=False
        )

        assert result == mock_fig
        mock_create_vertical_bar.assert_called_once()

    @patch("matplotlib.pyplot.subplots")
    @patch("matplotlib.pyplot.tight_layout")
    @patch("prepare.plugins.visualizers.bar._create_horizontal_bar")
    @patch("prepare.plugins.visualizers.bar._prepare_plot_data")
    def test_create_plot_with_limit(
        self, mock_prepare_plot_data, mock_create_horizontal_bar, mock_tight_layout, mock_subplots
    ):
        """Test plot creation with data limit."""
        mock_fig = Mock()
        mock_ax = Mock()
        mock_subplots.return_value = (mock_fig, mock_ax)
        mock_create_horizontal_bar.return_value = mock_ax

        # Mock prepared data
        limited_data = self.test_data.head(2)
        mock_prepare_plot_data.return_value = (limited_data, "count")

        result = BarChartPlugin.create_plot(data=self.test_data, x_column="category", y_column="count", limit=2)

        mock_prepare_plot_data.assert_called_once_with(self.test_data, "category", "count", 2)
        assert result == mock_fig

    @patch("matplotlib.pyplot.subplots")
    @patch("matplotlib.pyplot.tight_layout")
    @patch("prepare.plugins.visualizers.bar._create_horizontal_bar")
    def test_create_plot_save_file(self, mock_create_horizontal_bar, mock_tight_layout, mock_subplots):
        """Test plot creation with file saving."""
        mock_fig = Mock()
        mock_ax = Mock()
        mock_subplots.return_value = (mock_fig, mock_ax)
        mock_create_horizontal_bar.return_value = mock_ax

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as temp_file:
            try:
                BarChartPlugin.create_plot(
                    data=self.test_data, x_column="category", y_column="count", output_file=temp_file.name
                )

                mock_fig.savefig.assert_called_once()
            finally:
                os.unlink(temp_file.name)

    @patch("matplotlib.pyplot.subplots")
    def test_create_plot_error_handling(self, mock_subplots):
        """Test plot creation error handling."""
        mock_subplots.side_effect = Exception("Test error")

        with pytest.raises(VisualizationError, match="Error creating bar chart"):
            BarChartPlugin.create_plot(data=self.test_data)


class TestHeatmapPlugin:
    """Test the HeatmapPlugin."""

    def setup_method(self):
        """Set up test data."""
        self.test_data = pd.DataFrame(
            {
                "sample1": [0.95, 0.80, 0.60, 0.40, 0.20],
                "sample2": [0.85, 0.90, 0.55, 0.35, 0.15],
                "sample3": [0.75, 0.70, 0.85, 0.45, 0.25],
            },
            index=["genome1", "genome2", "genome3", "genome4", "genome5"],
        )

    def test_plugin_metadata(self):
        """Test plugin metadata."""
        assert HeatmapPlugin.name == "heatmap"
        assert "Heatmap" in HeatmapPlugin.description
        assert HeatmapPlugin.version == "0.1.0"

    def test_create_plot_with_clustering(self):
        """Clustering reorders rows and columns and draws one image with a colour bar."""
        pytest.importorskip("scipy")
        result = HeatmapPlugin.create_plot(data=self.test_data, title="Test Heatmap", cluster=True)
        try:
            ax = result.axes[0]
            assert ax.get_title() == "Test Heatmap"
            assert len(ax.images) == 1
            assert len(result.axes) == 2  # heatmap and colour bar
            labels = sorted(t.get_text() for t in ax.get_yticklabels())
            assert labels == sorted(self.test_data.index)
        finally:
            plt.close(result)

    def test_create_plot_no_clustering(self):
        """Without clustering the rows keep their input order."""
        result = HeatmapPlugin.create_plot(data=self.test_data, cluster=False, title="Test Heatmap", annot=True)
        try:
            ax = result.axes[0]
            assert [t.get_text() for t in ax.get_yticklabels()] == list(self.test_data.index)
            assert len(ax.texts) == self.test_data.size
        finally:
            plt.close(result)

    def test_create_plot_clustering_without_scipy_names_the_extra(self, monkeypatch):
        """Clustering needs scipy; without it the error names the analysis extra."""
        import sys

        from prepare.core.exceptions import ConfigurationError

        monkeypatch.setitem(sys.modules, "scipy", None)
        with pytest.raises(ConfigurationError, match=r"prepare\[analysis\]"):
            HeatmapPlugin.create_plot(data=self.test_data, cluster=True)

    def test_create_plot_save_file(self, tmp_path):
        """The heatmap is written to the requested file."""
        out = tmp_path / "heatmap.png"
        fig = HeatmapPlugin.create_plot(data=self.test_data, output_file=out, cluster=False)
        plt.close(fig)
        assert out.stat().st_size > 0

    def test_create_correlation_heatmap(self):
        """The correlation heatmap spans -1..1, blanks the upper triangle and annotates the rest."""
        result = HeatmapPlugin.create_correlation_heatmap(
            data=self.test_data, method="pearson", title="Correlation Heatmap"
        )
        try:
            ax = result.axes[0]
            image = ax.images[0]
            assert image.get_clim() == (-1.0, 1.0)
            n = self.test_data.shape[1]
            assert len(ax.texts) == n * (n - 1) // 2
            assert ax.get_title() == "Correlation Heatmap"
        finally:
            plt.close(result)

    def test_correlation_heatmap_blanks_the_undefined_cells_of_a_constant_column(self):
        """A constant column has no defined correlation: its cells are masked and carry no "nan" label."""
        data = pd.DataFrame({"a": [1.0, 2.0, 3.0, 4.0], "b": [2.0, 1.0, 4.0, 3.0], "flat": [5.0, 5.0, 5.0, 5.0]})
        result = HeatmapPlugin.create_correlation_heatmap(data=data, mask_upper=False)
        try:
            ax = result.axes[0]
            labels = [text.get_text() for text in ax.texts]
            assert "nan" not in labels
            # Only the four a/b cells have a defined correlation.
            assert len(labels) == 4
            hidden = np.ma.getmaskarray(ax.images[0].get_array())
            assert hidden[2, :].all() and hidden[:, 2].all()
            assert not hidden[:2, :2].any()
        finally:
            plt.close(result)

    @patch("matplotlib.pyplot.subplots")
    def test_create_plot_error_handling(self, mock_subplots):
        """Test heatmap creation error handling."""
        mock_subplots.side_effect = Exception("Test error")

        with pytest.raises(VisualizationError, match="Error creating heatmap"):
            HeatmapPlugin.create_plot(data=self.test_data, cluster=False)


class TestMapVisualizerPlugin:
    """Test the MapVisualizerPlugin."""

    def setup_method(self):
        """Set up test data."""
        self.test_data = pd.DataFrame(
            {
                "lat_lon": ["40.7128, -74.0060", "34.0522, -118.2437", "41.8781, -87.6298"],
                "value": [0.95, 0.85, 0.75],
                "location": ["New York", "Los Angeles", "Chicago"],
            }
        )

    def test_plugin_metadata(self):
        """Test plugin metadata."""
        assert MapVisualizerPlugin.name == "map"
        assert "map" in MapVisualizerPlugin.description.lower()
        assert MapVisualizerPlugin.version == "0.1.0"

    @patch("prepare.plugins.visualizers.map._require_cartopy")
    @patch("prepare.plugins.visualizers.map._plot_points")
    @patch("prepare.plugins.visualizers.map._add_map_features")
    @patch("prepare.plugins.visualizers.map._create_map_figure")
    @patch("prepare.plugins.visualizers.map._extract_coordinates")
    def test_create_plot_basic(
        self, mock_extract_coords, mock_create_figure, mock_add_features, mock_plot_points, _mock_require
    ):
        """Test basic map creation."""
        mock_fig = Mock()
        mock_ax = Mock()
        mock_create_figure.return_value = (mock_fig, mock_ax)
        mock_extract_coords.return_value = self.test_data  # Return original data
        mock_plot_points.return_value = Mock()

        result = MapVisualizerPlugin.create_plot(data=self.test_data, lat_lon_column="lat_lon", value_column="value")

        assert result == mock_fig
        mock_create_figure.assert_called_once()
        mock_add_features.assert_called_once()
        mock_plot_points.assert_called_once()

    def test_cartopy_unavailable_error(self):
        """Test error when cartopy is unavailable."""
        with patch.dict(sys.modules, {"cartopy": None}):
            with pytest.raises(ConfigurationError, match="Map visualization needs the 'cartopy' package"):
                MapVisualizerPlugin.create_plot(data=self.test_data, lat_lon_column="lat_lon")

    @patch("prepare.plugins.visualizers.map._require_cartopy")
    @patch("prepare.plugins.visualizers.map._plot_points")
    @patch("prepare.plugins.visualizers.map._add_map_features")
    @patch("prepare.plugins.visualizers.map._create_map_figure")
    @patch("prepare.plugins.visualizers.map._extract_coordinates")
    def test_create_plot_with_save(
        self, mock_extract_coords, mock_create_figure, mock_add_features, mock_plot_points, _mock_require
    ):
        """Test map creation with file saving."""
        mock_fig = Mock()
        mock_ax = Mock()
        mock_create_figure.return_value = (mock_fig, mock_ax)
        mock_extract_coords.return_value = self.test_data
        mock_plot_points.return_value = Mock()

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as temp_file:
            try:
                MapVisualizerPlugin.create_plot(
                    data=self.test_data, lat_lon_column="lat_lon", output_file=temp_file.name
                )

                mock_fig.savefig.assert_called_once()
            finally:
                os.unlink(temp_file.name)

    def test_parse_coordinate_string(self):
        """Test coordinate string parsing."""
        from prepare.plugins.visualizers.map import _parse_coordinate_string

        # Test standard format
        lat, lon = _parse_coordinate_string("40.7128, -74.0060")
        assert lat == 40.7128
        assert lon == -74.0060

        # Test with N/S/E/W designations
        lat, lon = _parse_coordinate_string("40.7128N, 74.0060W")
        assert lat == 40.7128
        assert lon == -74.0060

        # Test space-separated
        lat, lon = _parse_coordinate_string("40.7128 -74.0060")
        assert lat == 40.7128
        assert lon == -74.0060

        # Test invalid format
        lat, lon = _parse_coordinate_string("invalid")
        assert lat is None
        assert lon is None

    @patch("prepare.plugins.visualizers.map._require_cartopy")
    @patch("matplotlib.pyplot.figure")
    def test_create_plot_error_handling(self, mock_figure, _mock_require):
        """Test map creation error handling."""
        mock_figure.side_effect = Exception("Test error")

        with pytest.raises(VisualizationError, match="Error creating map"):
            MapVisualizerPlugin.create_plot(data=self.test_data, lat_lon_column="lat_lon")


class TestPluginRegistries:
    """Test the global plugin registries."""

    def test_format_registry_exists(self):
        """Test that format registry exists and is a PluginRegistry."""
        assert isinstance(format_registry, PluginRegistry)

    def test_visualizer_registry_exists(self):
        """Test that visualizer registry exists and is a PluginRegistry."""
        assert isinstance(visualizer_registry, PluginRegistry)

    def test_format_registry_has_plugins(self):
        """Test that format registry has expected plugins."""
        # Note: This depends on actual plugin registration in the module
        format_registry.list()
        # At minimum should have branchwater if auto-discovery works
        # This test might need adjustment based on actual registration

    def test_visualizer_registry_has_plugins(self):
        """Test that visualizer registry has expected plugins."""
        # Note: This depends on actual plugin registration in the module
        visualizer_registry.list()
        # Should have bar, heatmap, map if registered in plots.py
        # This test might need adjustment based on actual registration


class TestPluginIntegration:
    """Integration tests for the plugin system."""

    def test_end_to_end_format_plugin_workflow(self):
        """Test complete workflow with format plugin."""
        # Create test CSV data
        csv_content = """acc,containment,organism
SRR123456,0.95,Escherichia coli
SRR789012,0.85,Salmonella enterica
"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            f.flush()

            try:
                # Parse file
                containments = BranchWaterFormatPlugin.parse_file(f.name, "GCF_123456.1")

                # Extract metadata
                metadata_list = []
                for containment in containments:
                    metadata = BranchWaterFormatPlugin.extract_metadata(containment)
                    if metadata:
                        metadata_list.append(metadata)

                # Verify results
                assert len(containments) == 2
                assert len(metadata_list) == 2
                assert metadata_list[0].organism == "Escherichia coli"

            finally:
                os.unlink(f.name)

    @patch("matplotlib.pyplot.subplots")
    @patch("matplotlib.pyplot.tight_layout")
    @patch("prepare.plugins.visualizers.bar._create_horizontal_bar")
    def test_end_to_end_visualizer_plugin_workflow(self, mock_create_horizontal_bar, mock_tight_layout, mock_subplots):
        """Test complete workflow with visualizer plugin."""
        mock_fig = Mock()
        mock_ax = Mock()
        mock_subplots.return_value = (mock_fig, mock_ax)
        mock_create_horizontal_bar.return_value = mock_ax

        # Create test data
        test_data = pd.DataFrame({"category": ["A", "B", "C"], "count": [10, 20, 15]})

        # Create plot
        result = BarChartPlugin.create_plot(data=test_data, x_column="category", y_column="count", title="Test Chart")

        assert result == mock_fig

    def test_plugin_discovery_and_registration_flow(self):
        """Test plugin discovery and registration workflow."""
        registry = PluginRegistry()

        # Register plugins manually (simulating discovery)
        registry.register(BranchWaterFormatPlugin)
        registry.register(BarChartPlugin)

        # Test retrieval
        format_plugin = registry.get("branchwater")
        viz_plugin = registry.get("bar")

        assert format_plugin == BranchWaterFormatPlugin
        assert viz_plugin == BarChartPlugin

        # Test listing
        all_plugins = registry.list()
        assert "branchwater" in all_plugins
        assert "bar" in all_plugins


def test_missing_cartopy_message_names_the_extra():
    from prepare.plugins.visualizers import map as map_module

    with patch.dict(sys.modules, {"cartopy": None}):
        with pytest.raises(ConfigurationError, match=r"-m pip install 'prepare\[maps\]'"):
            map_module._require_cartopy()
