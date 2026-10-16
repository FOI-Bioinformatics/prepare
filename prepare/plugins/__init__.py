"""
Plugin system for PREPARE.

This package provides the plugin system for extending PREPARE functionality. The built-in
format and visualizer plugins are registered by load_builtin_plugins(), which the shared
registries in prepare.plugins.base call on their first lookup, so importing this package
(as the sourmash entry point does) does not import pandas or matplotlib. Third-party plugins
are found through the ``prepare_plugins`` entry-point group by discover_entry_points().
"""

import importlib.metadata
import logging
from typing import List

from prepare.core.constants import PLUGIN_REGISTRY_NAME
from prepare.core.exceptions import PluginError

logger = logging.getLogger(__name__)


def load_builtin_plugins() -> None:
    """Register the built-in plugins, then any found through entry points.

    Safe to call more than once: registering the same class again is a no-op.
    """
    from prepare.plugins.base import format_registry, visualizer_registry
    from prepare.plugins.formats.branchwater import BranchWaterFormatPlugin
    from prepare.plugins.visualizers.bar import BarChartPlugin
    from prepare.plugins.visualizers.heatmap import HeatmapPlugin
    from prepare.plugins.visualizers.map import MapVisualizerPlugin

    format_registry.register(BranchWaterFormatPlugin)
    for visualizer in (BarChartPlugin, HeatmapPlugin, MapVisualizerPlugin):
        visualizer_registry.register(visualizer)
    discover_entry_points()


def discover_entry_points(group: str = PLUGIN_REGISTRY_NAME) -> List[str]:
    """Load and register the plugin classes advertised under an entry-point group.

    A plugin that fails to import, is not a FormatPlugin or VisualizerPlugin subclass, or
    clashes with an already registered name is logged at WARNING and skipped.

    Returns:
        The names of the plugins registered from entry points.
    """
    from prepare.plugins.base import FormatPlugin, VisualizerPlugin, format_registry, visualizer_registry

    loaded: List[str] = []
    for entry_point in importlib.metadata.entry_points(group=group):
        try:
            plugin_class = entry_point.load()
            if isinstance(plugin_class, type) and issubclass(plugin_class, FormatPlugin):
                format_registry.register(plugin_class)
            elif isinstance(plugin_class, type) and issubclass(plugin_class, VisualizerPlugin):
                visualizer_registry.register(plugin_class)
            else:
                raise PluginError(f"{entry_point.value} is not a FormatPlugin or VisualizerPlugin subclass")
        except Exception as exc:  # noqa: B902 - a broken third-party plugin must not stop PREPARE
            logger.warning(f"Skipping plugin entry point '{entry_point.name}': {exc}")
            continue
        loaded.append(plugin_class.get_name())
    return loaded
