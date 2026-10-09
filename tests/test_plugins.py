"""
Unit tests for PREPARE plugin system.
"""

import pytest
import tempfile
import os

from prepare.core.exceptions import PluginError
from prepare.plugins.base import Plugin, PluginRegistry
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

    # Try to register duplicate plugin
    with pytest.raises(PluginError):
        registry.register(TestPlugin)

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
