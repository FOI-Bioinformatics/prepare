"""
Bar chart visualization plugin for PREPARE.
"""

import logging
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from typing import List, Optional, Tuple, Union

from prepare.core.exceptions import VisualizationError
from prepare.plugins.base import Plugin

logger = logging.getLogger(__name__)


def _prepare_plot_data(data, x_column, y_column, limit):
    """
    Prepare data for plotting.

    Args:
        data: DataFrame containing data
        x_column: Column to use for x-axis
        y_column: Column to use for y-axis
        limit: Maximum number of items to plot

    Returns:
        DataFrame prepared for plotting
    """
    df = data.copy()

    # Use first column if y_column not specified
    if y_column is None and not df.empty and df.shape[1] > 0:
        y_column = df.columns[0]

    # Use index if x_column not specified
    if x_column is None:
        if isinstance(df.index, pd.MultiIndex):
            raise VisualizationError("Cannot use MultiIndex for x-axis. Please specify x_column.")
        plot_df = df.copy()
    else:
        plot_df = df.set_index(x_column)

    # Apply limit if specified
    if limit is not None and limit > 0:
        if y_column:
            plot_df = plot_df.nlargest(limit, y_column)
        else:
            plot_df = plot_df.head(limit)

    return plot_df, y_column


def _create_horizontal_bar(plot_df, y_column, ax, colors, **kwargs):
    """
    Create a horizontal bar chart.

    Args:
        plot_df: DataFrame containing data
        y_column: Column to plot
        ax: Matplotlib axes to plot on
        colors: Colors to use
        **kwargs: Additional arguments for plotting

    Returns:
        Matplotlib axes with plot
    """
    if y_column:
        plot_df[y_column].plot(kind="barh", ax=ax, color=colors, **kwargs)
    else:
        plot_df.plot(kind="barh", ax=ax, color=colors, **kwargs)

    # Add labels
    ax.set_xlabel("Count")
    ax.set_ylabel("")

    return ax


def _create_vertical_bar(plot_df, y_column, ax, colors, **kwargs):
    """
    Create a vertical bar chart.

    Args:
        plot_df: DataFrame containing data
        y_column: Column to plot
        ax: Matplotlib axes to plot on
        colors: Colors to use
        **kwargs: Additional arguments for plotting

    Returns:
        Matplotlib axes with plot
    """
    if y_column:
        plot_df[y_column].plot(kind="bar", ax=ax, color=colors, **kwargs)
    else:
        plot_df.plot(kind="bar", ax=ax, color=colors, **kwargs)

    # Add labels
    ax.set_xlabel("")
    ax.set_ylabel("Count")

    # Rotate x-axis labels for vertical bar chart
    plt.xticks(rotation=45, ha="right")

    return ax


class BarChartPlugin(Plugin):
    """Plugin for creating bar chart visualizations."""

    name = "bar"
    description = "Bar chart visualization"
    version = "0.1.0"

    @classmethod
    def create_plot(
        cls,
        data: pd.DataFrame,
        x_column: Optional[str] = None,
        y_column: Optional[str] = None,
        title: Optional[str] = None,
        colors: Optional[Union[str, List[str]]] = None,
        figsize: Tuple[int, int] = (10, 6),
        horizontal: bool = True,
        limit: Optional[int] = None,
        output_file: Optional[Union[str, Path]] = None,
        output_format: str = "png",
        **kwargs,
    ) -> plt.Figure:
        """
        Create a bar chart visualization.

        Args:
            data: DataFrame containing data to plot
            x_column: Column to use for x-axis (if None, use index)
            y_column: Column to use for y-axis (if None, use first column)
            title: Title for the plot
            colors: Colors for the bars
            figsize: Figure size (width, height) in inches
            horizontal: If True, create horizontal bar chart
            limit: Maximum number of items to plot
            output_file: Path to save the plot
            output_format: Format to save the plot (png, jpg, pdf, svg)
            **kwargs: Additional arguments to pass to plotting function

        Returns:
            Matplotlib Figure object

        Raises:
            VisualizationError: If the plot cannot be created
        """
        try:
            # Prepare data
            plot_df, y_column = _prepare_plot_data(data, x_column, y_column, limit)

            # Create figure
            fig, ax = plt.subplots(figsize=figsize)

            # Create plot based on orientation
            if horizontal:
                ax = _create_horizontal_bar(plot_df, y_column, ax, colors, **kwargs)
            else:
                ax = _create_vertical_bar(plot_df, y_column, ax, colors, **kwargs)

            # Add title if specified
            if title:
                ax.set_title(title)

            # Adjust layout
            plt.tight_layout()

            # Save plot if output_file specified
            if output_file:
                fig.savefig(output_file, format=output_format, dpi=300, bbox_inches="tight")
                logger.info(f"Saved bar chart to {output_file}")

            return fig

        except Exception as e:
            raise VisualizationError(f"Error creating bar chart: {e}")
