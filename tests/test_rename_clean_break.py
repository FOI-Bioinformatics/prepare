"""The 0.10.0 rename is a clean break: old MetaQuest names are not read."""

from importlib.metadata import entry_points

from prepare.core import constants
from prepare.core.settings import SETTINGS
from prepare.data.registry import REGISTRY_FILENAME, load_registry
from prepare.store.resolve import resolve_store_root


def test_registry_filename_is_prepare():
    assert REGISTRY_FILENAME == "prepare_registry.json"


def test_store_marker_is_prepare():
    assert constants.STORE_MARKER == "prepare_store.json"


def test_plugin_registry_name_is_prepare():
    assert constants.PLUGIN_REGISTRY_NAME == "prepare_plugins"


def test_every_setting_env_var_uses_prepare_prefix():
    for spec in SETTINGS.values():
        assert spec.env.startswith("PREPARE_"), spec.env
        # An alias may be a variable owned by another service (NCBI_API_KEY), never an old tool name.
        for name in spec.aliases:
            assert not name.upper().startswith("METAQUEST"), name


def test_old_registry_file_is_not_read(tmp_path, monkeypatch):
    (tmp_path / "metaquest_registry.json").write_text('{"datasets": {"SRR1": {}}}')
    monkeypatch.chdir(tmp_path)
    registry = load_registry()
    assert registry.path == tmp_path.resolve() / "prepare_registry.json"
    assert not (tmp_path / "prepare_registry.json").exists()
    assert registry.datasets == {}


def test_old_store_variable_is_ignored(tmp_path, monkeypatch):
    monkeypatch.setenv("METAQUEST_DATA", str(tmp_path / "old"))
    monkeypatch.delenv("PREPARE_DATA", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    assert resolve_store_root(explicit=None, registry_root=None) is None


def test_sourmash_entry_points_use_prepare_names():
    names = {ep.name for ep in entry_points(group="sourmash.cli_script")}
    assert "prepare_parse" in names
    assert not any(name.startswith("metaquest") for name in names)
