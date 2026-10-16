"""
Base plugin system for PREPARE.

This module provides the foundational classes and registries for PREPARE's plugin system.
Plugins allow extending functionality without modifying the core codebase.
"""

import abc
import inspect
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, Generic, List, Optional, Type, TypeVar, Union

from prepare.core.exceptions import PluginError
from prepare.core.models import Containment

if TYPE_CHECKING:
    import pandas as pd
    from matplotlib.figure import Figure

logger = logging.getLogger(__name__)

# Type variable for plugin classes
T = TypeVar("T", bound="Plugin")


class Plugin(abc.ABC):
    """Base class for all PREPARE plugins."""

    # Plugin metadata
    name: str = "base"
    description: str = "Base plugin class"
    version: str = "0.1.0"

    @classmethod
    def get_name(cls) -> str:
        """Get the name of the plugin."""
        return cls.name

    @classmethod
    def get_description(cls) -> str:
        """Get the description of the plugin."""
        return cls.description

    @classmethod
    def get_version(cls) -> str:
        """Get the version of the plugin."""
        return cls.version


class FormatPlugin(Plugin):
    """Base class for plugins that read a containment file format."""

    @classmethod
    @abc.abstractmethod
    def validate_header(cls, headers: List[str]) -> bool:
        """Return True if the header row carries the columns the format requires."""

    @classmethod
    @abc.abstractmethod
    def parse_file(cls, file_path: Union[str, Path], genome_id: str) -> List[Containment]:
        """Parse one file into Containment records for genome_id."""


class VisualizerPlugin(Plugin):
    """Base class for plugins that draw a figure from a DataFrame."""

    @classmethod
    @abc.abstractmethod
    def create_plot(cls, data: "pd.DataFrame", **kwargs: Any) -> "Figure":
        """Draw data and return the matplotlib Figure."""


class PluginRegistry(Generic[T]):
    """
    Registry for managing plugin classes.

    This class maintains a registry of plugin classes of a specific type.
    It provides methods for registering, retrieving, and discovering plugins.
    """

    def __init__(self, loader: Optional[Callable[[], None]] = None) -> None:
        """
        Args:
            loader: Called before every lookup; the shared registries pass the built-in
                loader here so the plugin modules (and pandas/matplotlib) are imported only
                when a plugin is first looked up, not when this package is imported.
        """
        self._plugins: Dict[str, Type[T]] = {}
        self._loader = loader

    def _load(self) -> None:
        if self._loader is not None:
            self._loader()

    def register(self, plugin_class: Type[T]) -> None:
        """
        Register a plugin class.

        Registering the class already held under its name again is a no-op, so repeated
        loading (a module reload, a second load_builtin_plugins call) is harmless.

        Args:
            plugin_class: The plugin class to register

        Raises:
            PluginError: If the class is abstract, or a different class is already
                registered under the same name
        """
        if inspect.isabstract(plugin_class):
            missing = ", ".join(sorted(getattr(plugin_class, "__abstractmethods__", ())))
            raise PluginError(f"Plugin class {plugin_class.__qualname__} is abstract; it does not implement: {missing}")
        name = plugin_class.get_name()
        existing = self._plugins.get(name)
        if existing is plugin_class:
            return
        if existing is not None:
            raise PluginError(f"Plugin '{name}' is already registered ({existing.__module__}.{existing.__qualname__})")

        self._plugins[name] = plugin_class
        logger.debug(f"Registered plugin: {name} ({plugin_class.__module__}.{plugin_class.__name__})")

    def unregister(self, name: str) -> None:
        """
        Unregister a plugin by name.

        Args:
            name: The name of the plugin to unregister

        Raises:
            PluginError: If the plugin is not registered
        """
        self._load()
        if name not in self._plugins:
            raise PluginError(f"Plugin '{name}' is not registered")

        del self._plugins[name]
        logger.debug(f"Unregistered plugin: {name}")

    def get(self, name: str) -> Type[T]:
        """
        Get a plugin class by name.

        Args:
            name: The name of the plugin to retrieve

        Returns:
            The plugin class

        Raises:
            PluginError: If the plugin is not registered
        """
        self._load()
        if name not in self._plugins:
            raise PluginError(f"Plugin '{name}' is not registered")

        return self._plugins[name]

    def list(self) -> List[str]:
        """
        List all registered plugin names.

        Returns:
            A list of registered plugin names
        """
        self._load()
        return list(self._plugins.keys())

    def get_all(self) -> Dict[str, Type[T]]:
        """
        Get all registered plugins.

        Returns:
            A dictionary mapping plugin names to plugin classes
        """
        self._load()
        return self._plugins.copy()


_builtins_requested = False


def _ensure_builtins() -> None:
    """Run prepare.plugins.load_builtin_plugins() once, on the first registry lookup.

    The flag is set before loading, so a lookup made while a built-in module is being
    imported returns at once instead of recursing; a failed load clears it for a retry.
    """
    global _builtins_requested
    if _builtins_requested:
        return
    _builtins_requested = True
    loaded = False
    try:
        from prepare.plugins import load_builtin_plugins

        load_builtin_plugins()
        loaded = True
    finally:
        if not loaded:
            _builtins_requested = False


format_registry: PluginRegistry[FormatPlugin] = PluginRegistry(loader=_ensure_builtins)
visualizer_registry: PluginRegistry[VisualizerPlugin] = PluginRegistry(loader=_ensure_builtins)
