"""
Plotting functions for PREPARE.

This module provides functions for creating various types of plots.
"""

import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Optional, Type, Union, cast

from prepare.core.exceptions import VisualizationError
from prepare.plugins.base import visualizer_registry
from prepare.plugins.visualizers.heatmap import HeatmapPlugin

logger = logging.getLogger(__name__)


def _load_and_validate_data(
    file_path: Union[str, Path, pd.DataFrame], column: str, threshold: Optional[float]
) -> pd.DataFrame:
    """Load and validate containment data."""
    # Load data if file path, otherwise use the DataFrame directly
    if isinstance(file_path, pd.DataFrame):
        df = file_path.copy()
    else:
        df = pd.read_csv(file_path, sep="\t", index_col=0)

    # Check if column exists
    if column not in df.columns:
        raise VisualizationError(f"Column '{column}' not found in file. " f"Available columns: {', '.join(df.columns)}")

    # Apply threshold if specified
    if threshold is not None:
        df = df[df[column] >= threshold]
        if df.empty:
            logger.warning(f"No data above threshold {threshold}")
            return pd.DataFrame()  # Return empty DataFrame

    return df


def _create_rank_plot(ax: plt.Axes, df: pd.DataFrame, column: str, colors: Optional[Union[str, List[str]]]) -> None:
    """Create rank plot."""
    # Handle empty DataFrame or missing column
    if df.empty or column not in df.columns:
        ax.scatter([], [], color=colors)
        ax.set_xlabel("Rank")
        ax.set_ylabel(f"{column} Value")
        return

    df_sorted = df.sort_values(by=column, ascending=False)
    df_sorted["rank"] = np.arange(1, len(df_sorted) + 1)
    ax.scatter(df_sorted["rank"], df_sorted[column], color=colors)
    ax.set_xlabel("Rank")
    ax.set_ylabel(f"{column} Value")


def _create_histogram_plot(
    ax: plt.Axes, df: pd.DataFrame, column: str, colors: Optional[Union[str, List[str]]]
) -> None:
    """Create histogram plot."""
    # Handle empty DataFrame or missing column
    if df.empty or column not in df.columns:
        ax.hist([], bins=20, color=colors, alpha=0.7)
        ax.set_xlabel(column)
        ax.set_ylabel("Frequency")
        return

    # Adjust bin count for small datasets
    data_size = len(df)
    if data_size <= 2:
        bin_count = 5  # Minimum bins for very small datasets
    elif data_size <= 10:
        bin_count = min(20, max(5, data_size // 2 + 1))
    else:
        bin_count = 20  # Default for normal datasets

    ax.hist(df[column], bins=bin_count, color=colors, alpha=0.7)
    ax.set_xlabel(column)
    ax.set_ylabel("Frequency")


def _create_box_plot(ax: plt.Axes, df: pd.DataFrame, column: str) -> None:
    """Create box plot."""
    # Handle empty DataFrame or missing column
    if df.empty or column not in df.columns:
        ax.boxplot([[]])  # Empty boxplot
        ax.set_ylabel(column)
        ax.set_xticklabels([column])
        return

    ax.boxplot(df[column])
    ax.set_ylabel(column)
    ax.set_xticklabels([column])


def _create_violin_plot(ax: plt.Axes, df: pd.DataFrame, column: str) -> None:
    """Create violin plot."""
    # Handle empty DataFrame or missing column
    if df.empty or column not in df.columns:
        ax.violinplot([[]])  # Empty violin plot
        ax.set_ylabel(column)
        ax.set_xticklabels([column])
        return

    ax.violinplot(df[column])
    ax.set_ylabel(column)
    ax.set_xticklabels([column])


def _create_plot_by_type(
    ax: plt.Axes, df: pd.DataFrame, column: str, colors: Optional[Union[str, List[str]]], plot_type: str
) -> None:
    """Create plot based on plot type."""
    if plot_type == "rank":
        _create_rank_plot(ax, df, column, colors)
    elif plot_type == "histogram":
        _create_histogram_plot(ax, df, column, colors)
    elif plot_type == "box":
        _create_box_plot(ax, df, column)
    elif plot_type == "violin":
        _create_violin_plot(ax, df, column)
    else:
        raise VisualizationError(f"Unknown plot type: {plot_type}. " "Supported types: rank, histogram, box, violin")


def _beside_input(file_path: Union[str, Path], suffix: str, extension: str) -> Path:
    """``<stem><suffix>.<extension>`` in the folder of ``file_path``; the stem drops the input's own extension."""
    source = Path(file_path)
    return source.parent / f"{source.stem}{suffix}.{extension}"


def plot_output_path(file_path: Union[str, Path], plot_type: str, column: str, save_format: Optional[str]) -> Path:
    """Where a containment plot for ``file_path`` is saved; unset formats default to png.

    Named from the table's stem, not its full name, so a table named
    ``parsed_containment.txt`` produces ``parsed_containment_rank_max_containment.png``
    rather than carrying the ``.txt`` suffix into the middle of the image's name. The single
    source of truth for this naming, shared by ``_save_plot_if_needed`` and by
    ``PlotContainmentCommand`` (which logs the path but never receives it back from
    ``plot_containment``, since that keeps returning the ``Figure``).
    """
    fmt = save_format or "png"
    return _beside_input(file_path, f"_{plot_type}_{column}", fmt)


def _save_plot_if_needed(file_path: Union[str, Path], plot_type: str, column: str, save_format: Optional[str]) -> Path:
    """Save the current plot, always: unset formats default to png."""
    output_file = plot_output_path(file_path, plot_type, column, save_format)
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    logger.info(f"Plot saved to {output_file}")
    return output_file


def plot_containment(
    file_path: Union[str, Path, pd.DataFrame],
    column: str = "max_containment",
    title: Optional[str] = None,
    colors: Optional[Union[str, List[str]]] = None,
    show_title: bool = True,
    save_format: Optional[str] = None,
    threshold: Optional[float] = None,
    plot_type: str = "rank",
) -> Optional[plt.Figure]:
    """
    Plot containment data with various plot types.

    Args:
        file_path: Path to the containment file or DataFrame
        column: Column to plot
        title: Title for the plot
        colors: Colors to use in the plot
        show_title: Whether to display the title
        save_format: Format to save the plot as (png, jpg, pdf, svg); defaults to png when
            file_path is a table path. Skipped entirely when file_path is a DataFrame,
            since callers passing one already save the returned figure themselves.
        threshold: Minimum value to be included in the plot
        plot_type: Type of plot to generate ('rank', 'histogram', 'box', 'violin')

    Returns:
        Matplotlib Figure if successful, None otherwise

    Raises:
        VisualizationError: If the visualization fails
    """
    try:
        # Load and validate data
        df = _load_and_validate_data(file_path, column, threshold)
        if df.empty:
            return None

        # A title given explicitly is kept; the default one follows show_title.
        if title is None and show_title:
            title = f"{plot_type.capitalize()} Plot of {column}"

        # Set default colors if not specified; a single colour name is accepted as is
        if colors is None:
            colors = "blue"
        colors = [colors] if isinstance(colors, str) else colors

        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6))

        # Create plot based on plot_type
        _create_plot_by_type(ax, df, column, colors, plot_type)

        # Add title if specified
        if title:
            ax.set_title(title)

        # Adjust layout
        plt.tight_layout()

        # Save the plot, defaulting to png. A DataFrame has no path to name the file
        # from; callers passing one already save the returned figure themselves (e.g.
        # embedding it in a PDF page or report image), so saving is skipped for those.
        if not isinstance(file_path, pd.DataFrame):
            try:
                _save_plot_if_needed(file_path, plot_type, column, save_format)
            finally:
                plt.close(fig)

        return fig

    except Exception as e:
        if isinstance(e, VisualizationError):
            raise
        raise VisualizationError(f"Error plotting containment: {e}")


def _looks_numeric(value: object) -> bool:
    try:
        float(str(value))
    except ValueError:
        return False
    return True


def _load_counts_df(file_path: Union[str, Path, pd.DataFrame], limit: int) -> pd.DataFrame:
    """Load a counts table as [category, count] and keep the top `limit` rows.

    Two file shapes are accepted: the header-less ``<category>\\t<count>`` statistics
    file, and the ``count_metadata`` table whose first row is a header and whose
    remaining columns hold per-genome counts, which are summed per category.
    """
    if isinstance(file_path, pd.DataFrame):
        df = file_path.copy()
        if df.shape[1] < 2:
            raise VisualizationError("File must have at least two columns (category and count)")
        df.columns = ["category", "count"] + [f"col{i + 3}" for i in range(df.shape[1] - 2)]
    else:
        raw = pd.read_csv(file_path, sep="\t", header=None, dtype=str)
        if raw.shape[1] < 2:
            raise VisualizationError("File must have at least two columns (category and count)")
        if not _looks_numeric(raw.iat[0, 1]):
            raw = raw.iloc[1:].reset_index(drop=True)
        numeric = raw.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        if numeric.notna().to_numpy().sum() == 0:
            raise VisualizationError(f"No numeric count values found in {file_path}")
        counts = numeric.fillna(0).sum(axis=1)
        df = pd.DataFrame({"category": raw.iloc[:, 0].astype(str), "count": counts})

    return df.sort_values(by="count", ascending=False).head(limit)


def _metadata_pie_chart(df: pd.DataFrame, colors, title: Optional[str]) -> plt.Figure:
    """Render a pie chart of category counts."""
    fig, ax = plt.subplots(figsize=(10, 10))
    df["category"] = df["category"].astype(str)
    ax.pie(
        df["count"],
        labels=df["category"],
        autopct="%1.1f%%",
        startangle=90,
        colors=colors if colors else plt.cm.tab20.colors,  # type: ignore[attr-defined]
    )
    ax.axis("equal")  # Equal aspect ratio ensures that pie is drawn as a circle
    if title:
        ax.set_title(title)
    return fig


def _metadata_radar_chart(df: pd.DataFrame, colors, title: Optional[str]) -> plt.Figure:
    """Render a radar/polar chart of category counts (assumes >= 3 categories)."""
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "polar"})

    categories = df["category"].tolist()
    values = df["count"].tolist()
    angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()

    # Close the loop by repeating the first point.
    values += [values[0]]
    angles += [angles[0]]

    color = colors[0] if colors else "blue"
    ax.plot(angles, values, "o-", linewidth=2, color=color)
    ax.fill(angles, values, alpha=0.25, color=color)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories)
    if title:
        ax.set_title(title)
    return fig


def _metadata_plot_path(file_path: Union[str, Path], plot_type: str, save_format: str) -> Path:
    """Where a metadata counts plot is saved: beside its input, named from the input's stem."""
    return _beside_input(file_path, f"_{plot_type}", save_format)


def plot_metadata_counts(
    file_path: Union[str, Path, pd.DataFrame],
    title: Optional[str] = None,
    plot_type: str = "bar",
    colors: Optional[Union[str, List[str]]] = None,
    show_title: bool = True,
    save_format: Optional[str] = None,
    limit: int = 20,
) -> Optional[plt.Figure]:
    """
    Plot metadata counts with various plot types.

    Args:
        file_path: Path to the metadata counts file or DataFrame
        title: Title for the plot
        plot_type: Type of plot to generate ('bar', 'pie', 'radar')
        colors: Colors to use in the plot
        show_title: Whether to display the title
        save_format: Format to save the plot (png, jpg, pdf, svg)
        limit: Maximum number of items to include in the plot

    Returns:
        Matplotlib Figure if successful, None otherwise

    Raises:
        VisualizationError: If the visualization fails
    """
    try:
        df = _load_counts_df(file_path, limit)

        colors = [colors] if isinstance(colors, str) else colors

        # A title given explicitly is kept; the default one follows show_title.
        if title is None and show_title:
            title = f"Top {len(df)} Categories"

        # Create plot based on plot_type
        if plot_type == "bar":
            plugin = visualizer_registry.get("bar")
            bar_output = None
            if save_format and not isinstance(file_path, pd.DataFrame):
                bar_output = _metadata_plot_path(file_path, "bar", save_format)
            fig = plugin.create_plot(
                data=df,
                x_column="category",
                y_column="count",
                title=title,
                colors=colors,
                horizontal=True,
                output_file=bar_output,
                output_format=save_format or "png",
            )
        elif plot_type == "pie":
            fig = _metadata_pie_chart(df, colors, title)
        elif plot_type == "radar":
            # Radar needs at least 3 categories; otherwise fall back to a bar chart.
            if len(df) < 3:
                logger.warning("Radar chart requires at least 3 categories, falling back to bar chart")
                return plot_metadata_counts(file_path, title, "bar", colors, show_title, save_format, limit=limit)
            fig = _metadata_radar_chart(df, colors, title)
        else:
            raise VisualizationError(f"Unknown plot type: {plot_type}. Supported types: bar, pie, radar")

        # Save plot if format specified (the bar plugin handles its own saving)
        if save_format and plot_type != "bar" and not isinstance(file_path, pd.DataFrame):
            output_file = _metadata_plot_path(file_path, plot_type, save_format)
            try:
                fig.savefig(output_file, format=save_format, dpi=300, bbox_inches="tight")
                logger.info(f"Plot saved to {output_file}")
            finally:
                plt.close(fig)

        return fig

    except Exception as e:
        if isinstance(e, VisualizationError):
            raise
        raise VisualizationError(f"Error plotting metadata counts: {e}")


def plot_correlation_matrix(
    data: Union[str, Path, pd.DataFrame],
    title: Optional[str] = None,
    method: str = "pearson",
    output_file: Optional[Union[str, Path]] = None,
    output_format: str = "png",
) -> Optional[plt.Figure]:
    """
    Create a correlation matrix visualization.

    Args:
        data: DataFrame or path to data file
        title: Title for the plot
        method: Correlation method ('pearson', 'spearman', 'kendall')
        output_file: Path to save the plot
        output_format: Format to save the plot

    Returns:
        Matplotlib Figure if successful, None otherwise

    Raises:
        VisualizationError: If the visualization fails
    """
    try:
        # Load data if string or Path
        if isinstance(data, (str, Path)):
            df = pd.read_csv(data, sep="\t", index_col=0)
        else:
            df = data.copy()

        # Use heatmap plugin to create correlation matrix
        heatmap = cast(Type[HeatmapPlugin], visualizer_registry.get("heatmap"))
        fig = heatmap.create_correlation_heatmap(
            data=df,
            title=title,
            method=method,
            output_file=output_file,
            output_format=output_format,
        )

        return fig

    except Exception as e:
        if isinstance(e, VisualizationError):
            raise
        raise VisualizationError(f"Error creating correlation matrix: {e}")
