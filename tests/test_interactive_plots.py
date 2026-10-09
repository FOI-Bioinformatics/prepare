"""
Tests for interactive visualization functions.
"""

import numpy as np
import pandas as pd
import pytest
from unittest.mock import patch
from pathlib import Path
import tempfile

for _module in ("plotly", "sklearn", "scipy"):
    pytest.importorskip(_module)

from prepare.visualization.interactive import (  # noqa: E402
    create_interactive_pca,
    create_interactive_heatmap,
    create_diversity_comparison_plot,
)


class TestInteractivePlots:
    """Test interactive plotting functions."""

    def setup_method(self):
        """Set up test data."""
        # Create test data matrix
        np.random.seed(42)
        self.data_matrix = pd.DataFrame(
            np.random.randn(10, 5), index=[f"Sample_{i}" for i in range(10)], columns=[f"Feature_{i}" for i in range(5)]
        )

        self.metadata = pd.DataFrame(
            {"group": ["A"] * 5 + ["B"] * 5, "value": np.random.randn(10), "treatment": ["X", "Y"] * 5},
            index=self.data_matrix.index,
        )

        # Diversity data for testing
        self.alpha_diversity = pd.DataFrame(
            {
                "shannon": np.random.uniform(1, 3, 10),
                "simpson": np.random.uniform(0.5, 0.9, 10),
                "observed_species": np.random.randint(10, 50, 10),
            },
            index=self.data_matrix.index,
        )

    @patch("plotly.graph_objects.Figure.show")
    def test_create_interactive_pca_basic(self, mock_show):
        """Test basic PCA plot creation."""
        fig = create_interactive_pca(self.data_matrix, metadata=self.metadata, color_by="group", show_plot=False)

        # Check that figure was created
        assert fig is not None
        assert hasattr(fig, "data")
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_create_interactive_pca_3d(self, mock_show):
        """Test 3D PCA plot creation."""
        fig = create_interactive_pca(
            self.data_matrix,
            metadata=self.metadata,
            color_by="group",
            size_by="value",
            n_components=3,
            title="Test 3D PCA",
            show_plot=False,
        )

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_create_interactive_pca_numpy_input(self, mock_show):
        """Test PCA with numpy array input."""
        fig = create_interactive_pca(self.data_matrix.values, show_plot=False)

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_create_interactive_pca_save_file(self, mock_show):
        """Test PCA plot saving to file."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "test_pca.html"

            with patch("plotly.graph_objects.Figure.write_html") as mock_write:
                create_interactive_pca(self.data_matrix, output_file=output_file, show_plot=False)

                mock_write.assert_called_once_with(str(output_file))

    @patch("plotly.graph_objects.Figure.show")
    def test_create_interactive_heatmap_basic(self, mock_show):
        """Test basic heatmap creation."""
        fig = create_interactive_heatmap(self.data_matrix, cluster_samples=True, cluster_features=True, show_plot=False)

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_create_interactive_heatmap_no_clustering(self, mock_show):
        """Test heatmap without clustering."""
        fig = create_interactive_heatmap(
            self.data_matrix, cluster_samples=False, cluster_features=False, show_plot=False
        )

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_create_diversity_comparison_plot(self, mock_show):
        """Test diversity comparison plot."""
        fig = create_diversity_comparison_plot(
            self.alpha_diversity,
            self.metadata,
            group_by="group",
            diversity_metric="shannon",
            plot_type="box",
            show_plot=False,
        )

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_create_diversity_comparison_violin_plot(self, mock_show):
        """Test diversity comparison with violin plot."""
        fig = create_diversity_comparison_plot(
            self.alpha_diversity,
            self.metadata,
            group_by="group",
            diversity_metric="simpson",
            plot_type="violin",
            show_plot=False,
        )

        assert fig is not None
        mock_show.assert_not_called()

    def test_create_diversity_comparison_invalid_plot_type(self):
        """Test error handling for invalid plot type."""
        with pytest.raises(Exception):  # Should raise VisualizationError
            create_diversity_comparison_plot(
                self.alpha_diversity,
                self.metadata,
                group_by="group",
                diversity_metric="shannon",
                plot_type="invalid_type",
                show_plot=False,
            )

    @patch("plotly.graph_objects.Figure.show")
    def test_pca_missing_metadata_column(self, mock_show):
        """Test PCA with missing metadata column."""
        # Should not crash, just issue warning and proceed without coloring
        fig = create_interactive_pca(
            self.data_matrix, metadata=self.metadata, color_by="missing_column", show_plot=False
        )

        assert fig is not None
        mock_show.assert_not_called()


class TestInteractivePlotsEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_dataframe(self):
        """Test with empty DataFrame."""
        empty_df = pd.DataFrame()

        with pytest.raises(Exception):
            create_interactive_pca(empty_df, show_plot=False)

    def test_single_sample(self):
        """Test with single sample."""
        single_sample = pd.DataFrame([[1, 2, 3]], columns=["A", "B", "C"])

        with pytest.raises(Exception):  # PCA needs at least 2 samples
            create_interactive_pca(single_sample, show_plot=False)

    @patch("plotly.graph_objects.Figure.show")
    def test_small_dataset_pca(self, mock_show):
        """Test PCA with small dataset."""
        small_data = pd.DataFrame([[1, 2], [3, 4], [5, 6]], columns=["A", "B"])

        fig = create_interactive_pca(small_data, n_components=2, show_plot=False)

        assert fig is not None
        mock_show.assert_not_called()


if __name__ == "__main__":
    pytest.main([__file__])
