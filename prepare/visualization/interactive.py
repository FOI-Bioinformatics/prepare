"""
Interactive visualization functions for PREPARE.

This module provides interactive plotting capabilities using Plotly. Plotly comes from the
``interactive`` extra, and PCA, t-SNE and clustering from the ``analysis`` extra (scikit-learn,
scipy); each function imports what it needs when called, so importing this module needs neither.
"""

from __future__ import annotations

import logging
import numpy as np
import pandas as pd
from types import ModuleType
from typing import TYPE_CHECKING, Optional, Union
from pathlib import Path

from prepare.core.exceptions import VisualizationError
from prepare.core.optional import require

if TYPE_CHECKING:
    import plotly.graph_objects as go

logger = logging.getLogger(__name__)


def _plotly_express(purpose: str) -> ModuleType:
    """plotly.express, or a ConfigurationError naming the interactive extra."""
    return require("plotly.express", "interactive", purpose)


def _pca_matrix(data):
    """Return (X, sample_names) from a DataFrame or ndarray, validating PCA shape."""
    if isinstance(data, pd.DataFrame):
        X = data.values
        sample_names = data.index
    else:
        X = data
        sample_names = [f"Sample_{i}" for i in range(X.shape[0])]

    if X.shape[0] < 2:
        raise VisualizationError("PCA requires at least 2 samples")
    if X.shape[1] < 1:
        raise VisualizationError("PCA requires at least 1 feature")
    return X, sample_names


def _add_color_column(plot_data, metadata, color_by, sample_names):
    """Attach the color metadata column to plot_data; return the resolved column or None."""
    if metadata is None or color_by is None:
        return color_by
    if color_by in metadata.columns:
        plot_data[color_by] = metadata.loc[sample_names, color_by]
        return color_by
    logger.warning(f"Column '{color_by}' not found in metadata")
    return None


def _add_size_column(plot_data, metadata, size_by, sample_names):
    """Attach the size metadata column to plot_data; drop/skip it on NaN. Return resolved column or None."""
    if metadata is None or size_by is None:
        return size_by
    if size_by not in metadata.columns:
        logger.warning(f"Column '{size_by}' not found in metadata")
        return None

    size_values = metadata.loc[sample_names, size_by]
    if size_values.isna().any():
        logger.warning(f"Column '{size_by}' contains NaN values, removing size mapping")
        return None

    plot_data[size_by] = size_values
    if plot_data[size_by].isna().any():
        logger.warning(f"NaN values detected in plot data for '{size_by}', removing size mapping")
        plot_data.drop(columns=[size_by], inplace=True)
        return None
    return size_by


def _build_pca_figure(plot_data, n_components, n_pcs, color_by, size_by, title, variance) -> go.Figure:
    """Build the 1D/2D/3D PCA scatter depending on the number of components available."""
    px = _plotly_express("An interactive PCA plot")
    common = dict(color=color_by, size=size_by, hover_name="Sample", title=title)
    if n_components >= 3 and n_pcs >= 3:
        return px.scatter_3d(
            plot_data,
            x="PC1",
            y="PC2",
            z="PC3",
            labels={
                "PC1": f"PC1 ({variance[0]:.1%} variance)",
                "PC2": f"PC2 ({variance[1]:.1%} variance)",
                "PC3": f"PC3 ({variance[2]:.1%} variance)",
            },
            **common,
        )
    if n_pcs >= 2:
        return px.scatter(
            plot_data,
            x="PC1",
            y="PC2",
            labels={
                "PC1": f"PC1 ({variance[0]:.1%} variance)",
                "PC2": f"PC2 ({variance[1]:.1%} variance)",
            },
            **common,
        )
    return px.scatter(
        plot_data,
        x="PC1",
        y=[0] * len(plot_data),  # Flat y-axis for 1D data
        labels={"PC1": f"PC1 ({variance[0]:.1%} variance)"},
        **common,
    )


def create_interactive_pca(
    data: Union[pd.DataFrame, np.ndarray],
    metadata: Optional[pd.DataFrame] = None,
    color_by: Optional[str] = None,
    size_by: Optional[str] = None,
    n_components: int = 3,
    title: str = "Interactive PCA Plot",
    output_file: Optional[Union[str, Path]] = None,
    show_plot: bool = True,
) -> go.Figure:
    """
    Create an interactive 3D PCA plot with metadata overlay.

    Args:
        data: Sample x Features matrix
        metadata: Sample metadata for coloring/sizing points
        color_by: Metadata column for point colors
        size_by: Metadata column for point sizes
        n_components: Number of PCA components (2 or 3)
        title: Plot title
        output_file: HTML file to save plot
        show_plot: Whether to display plot

    Returns:
        Plotly Figure object

    Raises:
        ConfigurationError: If plotly or scikit-learn is not installed
        VisualizationError: If PCA fails
    """
    _plotly_express("An interactive PCA plot")
    decomposition = require("sklearn.decomposition", "analysis", "An interactive PCA plot")
    try:
        X, sample_names = _pca_matrix(data)

        pca = decomposition.PCA(n_components=min(n_components, X.shape[1], X.shape[0]))
        X_pca = pca.fit_transform(X)
        n_pcs = X_pca.shape[1]

        # Assemble the plotting frame (PC columns + optional metadata mappings)
        plot_data = pd.DataFrame({"Sample": sample_names})
        plot_data["PC1"] = X_pca[:, 0]
        if n_pcs >= 2:
            plot_data["PC2"] = X_pca[:, 1]
        if n_components >= 3 and n_pcs >= 3:
            plot_data["PC3"] = X_pca[:, 2]

        color_by = _add_color_column(plot_data, metadata, color_by, sample_names)
        size_by = _add_size_column(plot_data, metadata, size_by, sample_names)

        fig = _build_pca_figure(plot_data, n_components, n_pcs, color_by, size_by, title, pca.explained_variance_ratio_)

        fig.update_layout(width=800, height=600, showlegend=True, hovermode="closest")
        total_variance = sum(pca.explained_variance_ratio_[:n_components])
        fig.update_layout(title=f"{title}<br><sub>Total variance explained: {total_variance:.1%}</sub>")

        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.write_html(str(output_path))
            logger.info(f"Interactive PCA plot saved to {output_path}")

        if show_plot:
            fig.show()

        return fig

    except Exception as e:
        raise VisualizationError(f"Failed to create interactive PCA plot: {e}")


def create_interactive_tsne(
    data: Union[pd.DataFrame, np.ndarray],
    metadata: Optional[pd.DataFrame] = None,
    color_by: Optional[str] = None,
    perplexity: float = 30.0,
    n_iter: int = 1000,
    title: str = "Interactive t-SNE Plot",
    output_file: Optional[Union[str, Path]] = None,
    show_plot: bool = True,
) -> go.Figure:
    """
    Create an interactive t-SNE plot.

    Args:
        data: Sample x Features matrix
        metadata: Sample metadata for coloring points
        color_by: Metadata column for point colors
        perplexity: t-SNE perplexity parameter
        n_iter: Number of iterations (max_iter parameter in sklearn)
        title: Plot title
        output_file: HTML file to save plot
        show_plot: Whether to display plot

    Returns:
        Plotly Figure object

    Raises:
        ConfigurationError: If plotly or scikit-learn is not installed
    """
    px = _plotly_express("An interactive t-SNE plot")
    manifold = require("sklearn.manifold", "analysis", "An interactive t-SNE plot")
    try:
        # Prepare data
        if isinstance(data, pd.DataFrame):
            X = data.values
            sample_names = data.index
        else:
            X = data
            sample_names = [f"Sample_{i}" for i in range(X.shape[0])]

        # Perform t-SNE (use max_iter parameter for sklearn compatibility)
        tsne = manifold.TSNE(n_components=2, perplexity=perplexity, max_iter=n_iter, random_state=42)
        X_tsne = tsne.fit_transform(X)

        # Create DataFrame for plotting
        plot_data = pd.DataFrame({"tSNE1": X_tsne[:, 0], "tSNE2": X_tsne[:, 1], "Sample": sample_names})

        # Add metadata if provided
        if metadata is not None and color_by is not None:
            if color_by in metadata.columns:
                plot_data[color_by] = metadata.loc[sample_names, color_by]
            else:
                logger.warning(f"Column '{color_by}' not found in metadata")
                color_by = None

        # Create plot
        fig = px.scatter(
            plot_data,
            x="tSNE1",
            y="tSNE2",
            color=color_by,
            hover_name="Sample",
            title=title,
        )

        # Customize layout
        fig.update_layout(width=800, height=600, showlegend=True, hovermode="closest")

        # Save and show plot
        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.write_html(str(output_path))
            logger.info(f"Interactive t-SNE plot saved to {output_path}")

        if show_plot:
            fig.show()

        return fig

    except Exception as e:
        raise VisualizationError(f"Failed to create interactive t-SNE plot: {e}")


def create_interactive_heatmap(
    data: Union[pd.DataFrame, np.ndarray],
    sample_metadata: Optional[pd.DataFrame] = None,
    feature_metadata: Optional[pd.DataFrame] = None,
    cluster_samples: bool = True,
    cluster_features: bool = True,
    distance_metric: str = "euclidean",
    linkage_method: str = "ward",
    title: str = "Interactive Heatmap",
    output_file: Optional[Union[str, Path]] = None,
    show_plot: bool = True,
) -> go.Figure:
    """
    Create an interactive clustered heatmap with dendrograms.

    Args:
        data: Data matrix (samples x features or features x samples)
        sample_metadata: Metadata for sample annotations
        feature_metadata: Metadata for feature annotations
        cluster_samples: Whether to cluster samples
        cluster_features: Whether to cluster features
        distance_metric: Distance metric for clustering
        linkage_method: Linkage method for clustering
        title: Plot title
        output_file: HTML file to save plot
        show_plot: Whether to display plot

    Returns:
        Plotly Figure object

    Raises:
        ConfigurationError: If plotly, or scipy when clustering, is not installed
    """
    go = require("plotly.graph_objects", "interactive", "An interactive heatmap")
    if cluster_samples or cluster_features:
        hierarchy = require("scipy.cluster.hierarchy", "analysis", "Clustering an interactive heatmap")
    try:
        # Prepare data
        if isinstance(data, pd.DataFrame):
            data_matrix = data.values
            sample_names = data.index
            feature_names = data.columns
        else:
            data_matrix = data
            sample_names = [f"Sample_{i}" for i in range(data_matrix.shape[0])]
            feature_names = [f"Feature_{i}" for i in range(data_matrix.shape[1])]

        # Clustering
        sample_order = list(range(len(sample_names)))
        feature_order = list(range(len(feature_names)))

        if cluster_samples:
            sample_linkage = hierarchy.linkage(data_matrix, method=linkage_method, metric=distance_metric)
            sample_dendro = hierarchy.dendrogram(sample_linkage, no_plot=True)
            sample_order = sample_dendro["leaves"]

        if cluster_features:
            feature_linkage = hierarchy.linkage(data_matrix.T, method=linkage_method, metric=distance_metric)
            feature_dendro = hierarchy.dendrogram(feature_linkage, no_plot=True)
            feature_order = feature_dendro["leaves"]

        # Reorder data
        clustered_data = data_matrix[np.ix_(sample_order, feature_order)]
        clustered_sample_names = [sample_names[i] for i in sample_order]
        clustered_feature_names = [feature_names[i] for i in feature_order]

        # Create heatmap
        fig = go.Figure(
            data=go.Heatmap(
                z=clustered_data,
                x=clustered_feature_names,
                y=clustered_sample_names,
                colorscale="RdBu_r",
                showscale=True,
                hoverongaps=False,
                hovertemplate="Sample: %{y}<br>Feature: %{x}<br>" "Value: %{z}<extra></extra>",
            )
        )

        # Update layout
        fig.update_layout(
            title=title,
            width=max(800, len(clustered_feature_names) * 20),
            height=max(600, len(clustered_sample_names) * 15),
            xaxis={"side": "top"},
            yaxis={"autorange": "reversed"},
        )

        # Save and show plot
        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.write_html(str(output_path))
            logger.info(f"Interactive heatmap saved to {output_path}")

        if show_plot:
            fig.show()

        return fig

    except Exception as e:
        raise VisualizationError(f"Failed to create interactive heatmap: {e}")


def create_diversity_comparison_plot(
    alpha_diversity: pd.DataFrame,
    metadata: pd.DataFrame,
    group_by: str,
    diversity_metric: str = "shannon",
    plot_type: str = "box",
    title: Optional[str] = None,
    output_file: Optional[Union[str, Path]] = None,
    show_plot: bool = True,
) -> go.Figure:
    """
    Create interactive diversity comparison plots.

    Args:
        alpha_diversity: DataFrame with diversity metrics
        metadata: Sample metadata
        group_by: Column to group samples by
        diversity_metric: Which diversity metric to plot
        plot_type: Type of plot ("box", "violin", "strip")
        title: Plot title
        output_file: HTML file to save plot
        show_plot: Whether to display plot

    Returns:
        Plotly Figure object

    Raises:
        ConfigurationError: If plotly is not installed
    """
    px = _plotly_express("A diversity comparison plot")
    try:
        # Prepare data
        plot_data = alpha_diversity.copy()
        plot_data[group_by] = metadata.loc[alpha_diversity.index, group_by]

        # Set default title
        if title is None:
            title = f"{diversity_metric.title()} Diversity by {group_by}"

        # Create plot based on type
        if plot_type == "box":
            fig = px.box(plot_data, x=group_by, y=diversity_metric, title=title, points="all")
        elif plot_type == "violin":
            fig = px.violin(
                plot_data,
                x=group_by,
                y=diversity_metric,
                title=title,
                box=True,
                points="all",
            )
        elif plot_type == "strip":
            fig = px.strip(plot_data, x=group_by, y=diversity_metric, title=title)
        else:
            raise VisualizationError(f"Unknown plot type: {plot_type}")

        # Update layout
        fig.update_layout(width=800, height=600, showlegend=True)

        # Save and show plot
        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.write_html(str(output_path))
            logger.info(f"Diversity comparison plot saved to {output_path}")

        if show_plot:
            fig.show()

        return fig

    except Exception as e:
        raise VisualizationError(f"Failed to create diversity comparison plot: {e}")
