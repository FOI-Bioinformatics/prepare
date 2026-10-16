"""
Unit tests for PREPARE plugin system.
"""

import importlib
import inspect
import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from prepare.core.exceptions import PluginError
from prepare.plugins.base import Plugin, PluginRegistry, VisualizerPlugin
from prepare.plugins.formats.branchwater import BranchWaterFormatPlugin


class TestPlugin(Plugin):
    """Test plugin for unit tests."""

    name = "test"
    description = "Test plugin"
    version = "0.1.0"


def test_plugin_base():
    """Test the Plugin base class."""
    plugin = TestPlugin()

    assert plugin.get_name() == "test"
    assert plugin.get_description() == "Test plugin"
    assert plugin.get_version() == "0.1.0"


def test_plugin_registry():
    """Test the PluginRegistry class."""
    registry = PluginRegistry()

    # Register a plugin
    registry.register(TestPlugin)
    assert "test" in registry.list()

    # Get a plugin
    plugin_class = registry.get("test")
    assert plugin_class == TestPlugin

    # Get all plugins
    all_plugins = registry.get_all()
    assert "test" in all_plugins
    assert all_plugins["test"] == TestPlugin

    # Registering the same class again is a no-op
    registry.register(TestPlugin)
    assert registry.list() == ["test"]

    # Unregister a plugin
    registry.unregister("test")
    assert "test" not in registry.list()

    # Try to get non-existent plugin
    with pytest.raises(PluginError):
        registry.get("non_existent")

    # Try to unregister non-existent plugin
    with pytest.raises(PluginError):
        registry.unregister("non_existent")


def test_branchwater_format_plugin():
    """Test the BranchWaterFormatPlugin."""
    # Test required columns
    assert "acc" in BranchWaterFormatPlugin.REQUIRED_COLS
    assert "containment" in BranchWaterFormatPlugin.REQUIRED_COLS

    # Create a temporary CSV file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("acc,containment,cANI,biosample,bioproject,organism\n")
        f.write("SRR123456,0.95,0.98,SAMN123456,PRJNA123456,Escherichia coli\n")
        f.write("SRR789012,0.85,0.92,SAMN789012,PRJNA123456,Escherichia coli\n")
        f.flush()

        # Test header validation
        assert BranchWaterFormatPlugin.validate_header(
            ["acc", "containment", "cANI", "biosample", "bioproject", "organism"]
        )

        # Test parsing
        containments = BranchWaterFormatPlugin.parse_file(f.name, "GCF_123456.1")

        # Verify results
        assert len(containments) == 2
        assert containments[0].accession == "SRR123456"
        assert containments[0].value == 0.95
        assert containments[0].genome_id == "GCF_123456.1"
        assert containments[0].additional_data["organism"] == "Escherichia coli"

        # Test metadata extraction
        metadata = BranchWaterFormatPlugin.extract_metadata(containments[0])
        assert metadata.accession == "SRR123456"
        assert metadata.biosample == "SAMN123456"
        assert metadata.bioproject == "PRJNA123456"
        assert metadata.organism == "Escherichia coli"

    # Clean up
    os.unlink(f.name)


def test_branchwater_parse_file_tolerates_a_bom(tmp_path):
    """A file written by a spreadsheet program starts with a byte order mark."""
    path = tmp_path / "bom.csv"
    path.write_bytes("\ufeffacc,containment\r\nSRR123456,0.95\r\n".encode("utf-8"))

    containments = BranchWaterFormatPlugin.parse_file(path, "GCF_1")

    assert [c.accession for c in containments] == ["SRR123456"]


class OtherTestPlugin(Plugin):
    """A different class that claims the same name as TestPlugin."""

    name = "test"


def test_register_same_class_twice_is_a_no_op():
    registry = PluginRegistry()
    registry.register(TestPlugin)
    registry.register(TestPlugin)
    assert registry.get("test") is TestPlugin
    assert registry.list() == ["test"]


def test_register_different_class_under_same_name_raises():
    registry = PluginRegistry()
    registry.register(TestPlugin)
    with pytest.raises(PluginError, match="already registered"):
        registry.register(OtherTestPlugin)
    assert registry.get("test") is TestPlugin


def test_builtin_plugins_implement_the_abstract_interfaces():
    from prepare.plugins.base import FormatPlugin, VisualizerPlugin
    from prepare.plugins.visualizers.bar import BarChartPlugin
    from prepare.plugins.visualizers.heatmap import HeatmapPlugin
    from prepare.plugins.visualizers.map import MapVisualizerPlugin

    assert issubclass(BranchWaterFormatPlugin, FormatPlugin)
    for cls in (BarChartPlugin, HeatmapPlugin, MapVisualizerPlugin):
        assert issubclass(cls, VisualizerPlugin)
        assert not inspect.isabstract(cls)
    assert not inspect.isabstract(BranchWaterFormatPlugin)

    class IncompleteFormat(FormatPlugin):
        name = "incomplete"

    assert inspect.isabstract(IncompleteFormat)
    assert inspect.isabstract(VisualizerPlugin)


def test_load_builtin_plugins_is_idempotent():
    import prepare.plugins
    from prepare.plugins.base import format_registry, visualizer_registry

    prepare.plugins.load_builtin_plugins()
    prepare.plugins.load_builtin_plugins()
    assert "branchwater" in format_registry.list()
    assert {"bar", "heatmap", "map"} <= set(visualizer_registry.list())


def test_reloading_plots_module_does_not_raise():
    import prepare.visualization.plots as plots

    importlib.reload(plots)
    importlib.reload(plots)


class _BrokenEntryPoint:
    name = "broken_plugin"
    value = "broken_package.plugin:Broken"

    def load(self):
        raise ImportError("no module named broken_package")


class _GoodEntryPoint:
    name = "good_plugin"
    value = "good_package.plugin:Good"

    def __init__(self, cls):
        self._cls = cls

    def load(self):
        return self._cls


def test_discover_entry_points_skips_a_plugin_that_fails_to_load(caplog):
    from prepare.core.constants import PLUGIN_REGISTRY_NAME
    from prepare.plugins import discover_entry_points
    from prepare.plugins.base import VisualizerPlugin, visualizer_registry

    class EntryPointVisualizer(VisualizerPlugin):
        name = "entry_point_visualizer_test"

        @classmethod
        def create_plot(cls, data, **kwargs):
            return None

    calls = []

    def fake_entry_points(*, group):
        calls.append(group)
        return [_BrokenEntryPoint(), _GoodEntryPoint(EntryPointVisualizer)]

    try:
        with patch("importlib.metadata.entry_points", fake_entry_points):
            with caplog.at_level(logging.WARNING, logger="prepare.plugins"):
                loaded = discover_entry_points()
        assert calls == [PLUGIN_REGISTRY_NAME]
        assert loaded == ["entry_point_visualizer_test"]
        assert visualizer_registry.get("entry_point_visualizer_test") is EntryPointVisualizer
        warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
        assert any("broken_plugin" in r.getMessage() for r in warnings)
    finally:
        if "entry_point_visualizer_test" in visualizer_registry.list():
            visualizer_registry.unregister("entry_point_visualizer_test")


REPO_ROOT = Path(__file__).resolve().parent.parent


def _run_python(code):
    # Run from the repository root so ``-c`` finds this checkout's package first.
    return subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, timeout=60, cwd=REPO_ROOT, check=False
    )


def test_importing_the_plugin_package_does_not_import_matplotlib():
    """The sourmash entry point imports prepare.plugins; it must not pay for plotting imports."""
    code = (
        "import sys\n"
        "import prepare.plugins\n"
        "import prepare.plugins.base\n"
        "import prepare.plugins.sourmash_plugin\n"
        "sys.exit(1 if 'matplotlib.pyplot' in sys.modules else 0)\n"
    )
    result = _run_python(code)
    assert result.returncode == 0, result.stdout + result.stderr


def test_registries_return_the_builtins_after_a_fresh_import():
    code = (
        "from prepare.plugins.base import format_registry, visualizer_registry\n"
        "assert format_registry.list() == ['branchwater'], format_registry.list()\n"
        "assert {'bar', 'heatmap', 'map'} <= set(visualizer_registry.get_all()), visualizer_registry.list()\n"
        "assert visualizer_registry.get('heatmap').__name__ == 'HeatmapPlugin'\n"
    )
    result = _run_python(code)
    assert result.returncode == 0, result.stdout + result.stderr


class IncompleteVisualizer(VisualizerPlugin):
    """A visualizer that forgets create_plot."""

    name = "incomplete_visualizer_test"


def test_register_rejects_an_abstract_plugin():
    registry = PluginRegistry()
    with pytest.raises(PluginError, match="IncompleteVisualizer.*create_plot"):
        registry.register(IncompleteVisualizer)
    assert registry.list() == []


def test_discover_entry_points_skips_an_abstract_plugin(caplog):
    from prepare.plugins import discover_entry_points
    from prepare.plugins.base import visualizer_registry

    def fake_entry_points(*, group):
        return [_GoodEntryPoint(IncompleteVisualizer)]

    with patch("importlib.metadata.entry_points", fake_entry_points):
        with caplog.at_level(logging.WARNING, logger="prepare.plugins"):
            loaded = discover_entry_points()
    assert loaded == []
    assert "incomplete_visualizer_test" not in visualizer_registry.list()
    assert any("good_plugin" in r.getMessage() and "create_plot" in r.getMessage() for r in caplog.records)
