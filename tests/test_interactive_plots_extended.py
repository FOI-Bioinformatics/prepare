"""
EXTENDED TESTS for visualization/interactive.py (65% -> 80%+ coverage)

This file provides comprehensive testing for interactive visualization functions,
focusing on uncovered code paths, error handling, and edge cases.

Run: pytest tests/test_interactive_plots_extended.py -v
"""

import numpy as np
import pandas as pd
import pytest
from unittest.mock import patch

for _module in ("plotly", "sklearn", "scipy"):
    pytest.importorskip(_module)

from prepare.visualization.interactive import (  # noqa: E402
    create_interactive_pca,
    create_interactive_tsne,
    create_interactive_heatmap,
    create_diversity_comparison_plot,
)
from prepare.core.exceptions import VisualizationError  # noqa: E402


@pytest.fixture
def sample_data():
    """Create sample test data."""
    np.random.seed(42)
    data = pd.DataFrame(
        np.random.randn(10, 5), index=[f"Sample_{i}" for i in range(10)], columns=[f"Feature_{i}" for i in range(5)]
    )
    return data


@pytest.fixture
def sample_metadata():
    """Create sample metadata."""
    return pd.DataFrame(
        {
            "group": ["A"] * 5 + ["B"] * 5,
            "value": np.random.randn(10),
            "size_value": np.random.uniform(1, 10, 10),
            "treatment": ["X", "Y"] * 5,
        },
        index=[f"Sample_{i}" for i in range(10)],
    )


@pytest.fixture
def sample_alpha_diversity():
    """Create sample alpha diversity data."""
    np.random.seed(42)
    return pd.DataFrame(
        {
            "shannon": np.random.uniform(1, 3, 10),
            "simpson": np.random.uniform(0.5, 0.9, 10),
            "observed_species": np.random.randint(10, 50, 10),
        },
        index=[f"Sample_{i}" for i in range(10)],
    )


class TestInteractiveTSNE:
    """Test t-SNE plotting function."""

    @patch("plotly.graph_objects.Figure.show")
    def test_create_tsne_basic(self, mock_show, sample_data):
        """Test basic t-SNE plot creation."""
        fig = create_interactive_tsne(sample_data, perplexity=5.0, show_plot=False)  # Must be < n_samples (10)

        assert fig is not None
        assert hasattr(fig, "data")
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_create_tsne_with_metadata(self, mock_show, sample_data, sample_metadata):
        """Test t-SNE with metadata coloring."""
        fig = create_interactive_tsne(
            sample_data,
            metadata=sample_metadata,
            color_by="group",
            title="Test t-SNE",
            perplexity=5.0,  # Must be < n_samples (10)
            show_plot=False,
        )

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_create_tsne_with_numpy(self, mock_show, sample_data):
        """Test t-SNE with numpy array input."""
        fig = create_interactive_tsne(sample_data.values, perplexity=5.0, show_plot=False)  # Must be < n_samples (10)

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_create_tsne_custom_parameters(self, mock_show, sample_data):
        """Test t-SNE with custom perplexity and iterations."""
        fig = create_interactive_tsne(sample_data, perplexity=5.0, n_iter=500, show_plot=False)

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_create_tsne_save_file(self, mock_show, sample_data, tmp_path):
        """Test t-SNE plot saving to file."""
        output_file = tmp_path / "test_tsne.html"

        with patch("plotly.graph_objects.Figure.write_html") as mock_write:
            create_interactive_tsne(
                sample_data, output_file=output_file, perplexity=5.0, show_plot=False  # Must be < n_samples (10)
            )

            mock_write.assert_called_once_with(str(output_file))

    @patch("plotly.graph_objects.Figure.show")
    def test_create_tsne_missing_metadata_column(self, mock_show, sample_data, sample_metadata):
        """Test t-SNE with missing metadata column."""
        fig = create_interactive_tsne(
            sample_data,
            metadata=sample_metadata,
            color_by="nonexistent_column",
            perplexity=5.0,  # Must be < n_samples (10)
            show_plot=False,
        )

        assert fig is not None
        mock_show.assert_not_called()

    def test_create_tsne_error_handling(self, sample_data):
        """Test t-SNE error handling."""
        # Force an error by passing invalid data
        with pytest.raises(VisualizationError, match="Failed to create interactive t-SNE"):
            with patch("sklearn.manifold.TSNE") as mock_tsne:
                mock_tsne.return_value.fit_transform.side_effect = ValueError("Test error")
                create_interactive_tsne(sample_data, show_plot=False)


class TestInteractivePCAExtended:
    """Extended tests for PCA plotting."""

    @patch("plotly.graph_objects.Figure.show")
    def test_pca_with_nan_size_values(self, mock_show, sample_data, sample_metadata):
        """Test PCA with NaN values in size_by column."""
        # Add NaN to size values
        metadata_with_nan = sample_metadata.copy()
        metadata_with_nan.loc["Sample_0", "size_value"] = np.nan

        fig = create_interactive_pca(
            sample_data, metadata=metadata_with_nan, color_by="group", size_by="size_value", show_plot=False
        )

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_pca_with_missing_size_column(self, mock_show, sample_data, sample_metadata):
        """Test PCA with missing size_by column."""
        fig = create_interactive_pca(
            sample_data, metadata=sample_metadata, size_by="nonexistent_column", show_plot=False
        )

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_pca_2d_explicit(self, mock_show, sample_data):
        """Test PCA with explicit 2D components."""
        fig = create_interactive_pca(sample_data, n_components=2, show_plot=False)

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_pca_show_plot_true(self, mock_show, sample_data):
        """Test PCA with show_plot=True."""
        fig = create_interactive_pca(sample_data, show_plot=True)

        assert fig is not None
        mock_show.assert_called_once()

    def test_pca_insufficient_features(self):
        """Test PCA with no features."""
        # Test with no features
        data_no_features = pd.DataFrame(index=["A", "B"])
        with pytest.raises(VisualizationError, match="PCA requires at least 1 feature"):
            create_interactive_pca(data_no_features, show_plot=False)

    def test_pca_error_in_computation(self, sample_data):
        """Test PCA error handling during computation."""
        with pytest.raises(VisualizationError, match="Failed to create interactive PCA"):
            with patch("sklearn.decomposition.PCA") as mock_pca:
                mock_pca.return_value.fit_transform.side_effect = ValueError("PCA failed")
                create_interactive_pca(sample_data, show_plot=False)


class TestInteractiveHeatmapExtended:
    """Extended tests for heatmap plotting."""

    @patch("plotly.graph_objects.Figure.show")
    def test_heatmap_with_numpy_input(self, mock_show):
        """Test heatmap with numpy array input."""
        data = np.random.randn(5, 4)

        fig = create_interactive_heatmap(data, cluster_samples=False, cluster_features=False, show_plot=False)

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_heatmap_save_file(self, mock_show, sample_data, tmp_path):
        """Test heatmap saving to file."""
        output_file = tmp_path / "test_heatmap.html"

        with patch("plotly.graph_objects.Figure.write_html") as mock_write:
            create_interactive_heatmap(sample_data, output_file=output_file, show_plot=False)

            mock_write.assert_called_once_with(str(output_file))

    @patch("plotly.graph_objects.Figure.show")
    def test_heatmap_custom_clustering_parameters(self, mock_show, sample_data):
        """Test heatmap with custom clustering parameters."""
        fig = create_interactive_heatmap(
            sample_data,
            cluster_samples=True,
            cluster_features=True,
            distance_metric="cityblock",
            linkage_method="average",
            title="Custom Heatmap",
            show_plot=False,
        )

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_heatmap_show_plot_true(self, mock_show, sample_data):
        """Test heatmap with show_plot=True."""
        fig = create_interactive_heatmap(sample_data, show_plot=True)

        assert fig is not None
        mock_show.assert_called_once()

    def test_heatmap_error_handling(self, sample_data):
        """Test heatmap error handling."""
        with pytest.raises(VisualizationError, match="Failed to create interactive heatmap"):
            with patch("scipy.cluster.hierarchy.linkage") as mock_linkage:
                mock_linkage.side_effect = ValueError("Clustering failed")
                create_interactive_heatmap(sample_data, cluster_samples=True, show_plot=False)


class TestDiversityComparisonExtended:
    """Extended tests for diversity comparison plots."""

    @patch("plotly.graph_objects.Figure.show")
    def test_diversity_strip_plot(self, mock_show, sample_alpha_diversity, sample_metadata):
        """Test diversity comparison with strip plot."""
        fig = create_diversity_comparison_plot(
            sample_alpha_diversity,
            sample_metadata,
            group_by="group",
            diversity_metric="shannon",
            plot_type="strip",
            show_plot=False,
        )

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_diversity_save_file(self, mock_show, sample_alpha_diversity, sample_metadata, tmp_path):
        """Test diversity plot saving to file."""
        output_file = tmp_path / "test_diversity.html"

        with patch("plotly.graph_objects.Figure.write_html") as mock_write:
            create_diversity_comparison_plot(
                sample_alpha_diversity,
                sample_metadata,
                group_by="group",
                diversity_metric="shannon",
                output_file=output_file,
                show_plot=False,
            )

            mock_write.assert_called_once_with(str(output_file))

    @patch("plotly.graph_objects.Figure.show")
    def test_diversity_custom_title(self, mock_show, sample_alpha_diversity, sample_metadata):
        """Test diversity plot with custom title."""
        fig = create_diversity_comparison_plot(
            sample_alpha_diversity,
            sample_metadata,
            group_by="group",
            diversity_metric="simpson",
            title="Custom Diversity Plot",
            show_plot=False,
        )

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_diversity_show_plot_true(self, mock_show, sample_alpha_diversity, sample_metadata):
        """Test diversity plot with show_plot=True."""
        fig = create_diversity_comparison_plot(
            sample_alpha_diversity, sample_metadata, group_by="group", diversity_metric="shannon", show_plot=True
        )

        assert fig is not None
        mock_show.assert_called_once()

    def test_diversity_error_handling(self, sample_alpha_diversity, sample_metadata):
        """Test diversity plot error handling."""
        with pytest.raises(VisualizationError, match="Failed to create diversity comparison plot"):
            with patch("plotly.express.box") as mock_box:
                mock_box.side_effect = ValueError("Plotting failed")
                create_diversity_comparison_plot(
                    sample_alpha_diversity,
                    sample_metadata,
                    group_by="group",
                    diversity_metric="shannon",
                    show_plot=False,
                )


class TestEdgeCasesExtended:
    """Additional edge case tests."""

    def test_pca_single_feature(self):
        """Test PCA with single feature."""
        data = pd.DataFrame({"Feature_0": np.random.randn(5)})

        # Should work but only produce 1 component
        with patch("plotly.graph_objects.Figure.show"):
            fig = create_interactive_pca(data, n_components=2, show_plot=False)
            assert fig is not None

    @patch("plotly.graph_objects.Figure.show")
    def test_heatmap_small_dataset(self, mock_show):
        """Test heatmap with very small dataset."""
        small_data = pd.DataFrame([[1, 2], [3, 4]], columns=["A", "B"])

        fig = create_interactive_heatmap(small_data, cluster_samples=True, cluster_features=True, show_plot=False)

        assert fig is not None
        mock_show.assert_not_called()

    @patch("plotly.graph_objects.Figure.show")
    def test_pca_with_all_nan_size_values(self, mock_show, sample_data, sample_metadata):
        """Test PCA when all size values are NaN."""
        metadata_all_nan = sample_metadata.copy()
        metadata_all_nan["size_value"] = np.nan

        fig = create_interactive_pca(sample_data, metadata=metadata_all_nan, size_by="size_value", show_plot=False)

        assert fig is not None
        mock_show.assert_not_called()


# ============================================================================
# SUCCESS METRICS:
#
# After running these extended tests:
# - Expected: 37+ tests pass
# - Coverage: 65% -> 80%+ for visualization/interactive.py
# - All interactive functions comprehensively tested
#
# Run tests:
#   pytest tests/test_interactive_plots_extended.py -v
#
# Check coverage:
#   pytest --cov=prepare.visualization.interactive --cov-report=term-missing \
#          tests/test_interactive_plots.py tests/test_interactive_plots_extended.py
# ============================================================================
