"""
Visualization modules for PREPARE.

This package contains modules for creating visualizations and reports.
"""

# Import and re-export for type checking
from prepare.visualization.plots import (
    plot_containment,
    plot_metadata_counts,
    plot_correlation_matrix,
)

__all__ = ["plot_containment", "plot_metadata_counts", "plot_correlation_matrix"]
