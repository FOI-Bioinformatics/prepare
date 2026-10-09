"""
Heatmap visualization plugin for PREPARE.

Heatmaps are drawn with matplotlib's ``imshow`` and a colour bar. Clustering the rows and
columns of a plain heatmap uses scipy, from the ``analysis`` extra.
"""

import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, Union

from prepare.core.exceptions import ConfigurationError, VisualizationError
from prepare.core.optional import require
from prepare.plugins.base import Plugin

logger = logging.getLogger(__name__)


def _cluster_order(matrix: np.ndarray) -> np.ndarray:
    """Leaf order of an average-linkage clustering of the rows of ``matrix``."""
    hierarchy = require("scipy.cluster.hierarchy", "analysis", "Clustering a heatmap")
    if matrix.shape[0] < 2:
        return np.arange(matrix.shape[0])
    return hierarchy.leaves_list(hierarchy.linkage(matrix, method="average", metric="euclidean"))


def _draw_matrix(
    fig: plt.Figure,
    ax: plt.Axes,
    df: pd.DataFrame,
    cmap: str,
    annot: bool,
    fmt: str = ".2f",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    mask: Optional[np.ndarray] = None,
    linewidths: float = 0,
) -> None:
    """Draw ``df`` as an image with labelled axes and a colour bar; masked cells are left blank.

    Non-finite cells (for example the correlations of a constant column) are masked as well, so
    they are drawn blank and carry no "nan" annotation.
    """
    raw = df.to_numpy(dtype=float)
    hidden = ~np.isfinite(raw)
    if mask is not None:
        hidden |= np.asarray(mask, dtype=bool)
    values = np.ma.masked_array(raw, mask=hidden)
    image = ax.imshow(values, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto", interpolation="nearest")
    fig.colorbar(image, ax=ax)

    ax.set_xticks(range(df.shape[1]))
    ax.set_xticklabels([str(c) for c in df.columns], rotation=90)
    ax.set_yticks(range(df.shape[0]))
    ax.set_yticklabels([str(i) for i in df.index])

    if linewidths > 0:
        ax.set_xticks(np.arange(-0.5, df.shape[1], 1), minor=True)
        ax.set_yticks(np.arange(-0.5, df.shape[0], 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=linewidths)
        ax.tick_params(which="minor", length=0)

    if annot:
        hidden = np.ma.getmaskarray(values)
        for row in range(df.shape[0]):
            for col in range(df.shape[1]):
                if hidden[row, col]:
                    continue
                ax.text(col, row, format(values.data[row, col], fmt), ha="center", va="center", fontsize=8)


class HeatmapPlugin(Plugin):
    """Plugin for creating heatmap visualizations."""

    name = "heatmap"
    description = "Heatmap visualization"
    version = "0.1.0"

    @classmethod
    def create_plot(
        cls,
        data: pd.DataFrame,
        title: Optional[str] = None,
        cmap: str = "viridis",
        figsize: Tuple[int, int] = (12, 10),
        cluster: bool = True,
        annot: bool = False,
        linewidths: float = 0,
        output_file: Optional[Union[str, Path]] = None,
        output_format: str = "png",
        **kwargs,
    ) -> plt.Figure:
        """
        Create a heatmap visualization.

        Args:
            data: DataFrame containing data to plot
            title: Title for the plot
            cmap: Colormap to use
            figsize: Figure size (width, height) in inches
            cluster: If True, order rows and columns by hierarchical clustering (needs scipy)
            annot: If True, annotate cells with values
            linewidths: Width of lines separating cells
            output_file: Path to save the plot
            output_format: Format to save the plot (png, jpg, pdf, svg)
            **kwargs: Accepted for compatibility with other visualizer plugins and ignored

        Returns:
            Matplotlib Figure object

        Raises:
            ConfigurationError: If clustering is requested and scipy is not installed
            VisualizationError: If the plot cannot be created
        """
        try:
            df = data.copy()

            if cluster:
                matrix = df.to_numpy(dtype=float)
                df = df.iloc[_cluster_order(matrix), _cluster_order(matrix.T)]

            fig, ax = plt.subplots(figsize=figsize)
            _draw_matrix(fig, ax, df, cmap=cmap, annot=annot, linewidths=linewidths)

            if title:
                ax.set_title(title)

            plt.tight_layout()

            if output_file:
                fig.savefig(output_file, format=output_format, dpi=300, bbox_inches="tight")
                logger.info(f"Saved heatmap to {output_file}")

            return fig

        except ConfigurationError:
            raise
        except Exception as e:
            raise VisualizationError(f"Error creating heatmap: {e}")

    @classmethod
    def create_correlation_heatmap(
        cls,
        data: pd.DataFrame,
        method: str = "pearson",
        title: Optional[str] = None,
        cmap: str = "coolwarm",
        figsize: Tuple[int, int] = (12, 10),
        mask_upper: bool = True,
        annot: bool = True,
        output_file: Optional[Union[str, Path]] = None,
        output_format: str = "png",
        **kwargs,
    ) -> plt.Figure:
        """
        Create a correlation heatmap visualization.

        Args:
            data: DataFrame containing data to calculate correlations from
            method: Correlation method ('pearson', 'spearman', or 'kendall')
            title: Title for the plot
            cmap: Colormap to use
            figsize: Figure size (width, height) in inches
            mask_upper: If True, leave the upper triangle (and diagonal) blank
            annot: If True, annotate cells with values
            output_file: Path to save the plot
            output_format: Format to save the plot (png, jpg, pdf, svg)
            **kwargs: Accepted for compatibility with other visualizer plugins and ignored

        Returns:
            Matplotlib Figure object

        Raises:
            VisualizationError: If the plot cannot be created
        """
        try:
            corr_matrix = data.corr(method=method)
            mask = np.triu(np.ones(corr_matrix.shape, dtype=bool)) if mask_upper else None

            fig, ax = plt.subplots(figsize=figsize)
            _draw_matrix(fig, ax, corr_matrix, cmap=cmap, annot=annot, vmin=-1.0, vmax=1.0, mask=mask, linewidths=0.5)
            ax.set_aspect("equal")

            if title:
                ax.set_title(title)

            plt.tight_layout()

            if output_file:
                fig.savefig(output_file, format=output_format, dpi=300, bbox_inches="tight")
                logger.info(f"Saved correlation heatmap to {output_file}")

            return fig

        except Exception as e:
            raise VisualizationError(f"Error creating correlation heatmap: {e}")
