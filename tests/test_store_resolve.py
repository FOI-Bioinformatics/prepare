"""
Tests for prepare.store.resolve: store root discovery and user config.

Every test runs under tmp_path and monkeypatches HOME, XDG_CONFIG_HOME and
PREPARE_DATA so nothing here reads or writes the real user config.
"""

import os

import pytest

from prepare.core.constants import CONFIG_DIRNAME, CONFIG_FILENAME, STORE_ENV
from prepare.core.exceptions import ConfigurationError, DataAccessError
from prepare.store.layout import init_store
from prepare.store.resolve import (
    config_path,
    read_config,
    resolve_store_root,
    write_config_data_root,
)
from prepare.utils.security import SecureSubprocess


@pytest.fixture(autouse=True)
def isolated_env(tmp_path, monkeypatch):
    """Point HOME/XDG_CONFIG_HOME at tmp_path and clear PREPARE_DATA."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    monkeypatch.delenv(STORE_ENV, raising=False)
    yield


@pytest.fixture(autouse=True)
def reset_allowed_roots():
    """SecureSubprocess.allowed_roots keeps class-level state; reset around each test."""
    SecureSubprocess._extra_roots = []
    yield
    SecureSubprocess._extra_roots = []


# --- config_path -----------------------------------------------------------


def test_config_path_defaults_to_dot_config(tmp_path):
    path = config_path()
    assert path == tmp_path / ".config" / CONFIG_DIRNAME / CONFIG_FILENAME


def test_config_path_honors_xdg_config_home(tmp_path, monkeypatch):
    xdg = tmp_path / "xdg"
    monkeypatch.setenv("XDG_CONFIG_HOME", str(xdg))
    path = config_path()
    assert path == xdg / CONFIG_DIRNAME / CONFIG_FILENAME


# --- read_config / write_config_data_root -----------------------------------


def test_read_config_returns_empty_dict_when_absent():
    assert read_config() == {}


def test_write_config_data_root_round_trips(tmp_path):
    root = tmp_path / "mystore"
    write_config_data_root(root)

    config = read_config()
    assert config["store"]["data_root"] == root.as_posix()


def test_write_config_data_root_preserves_other_tables(tmp_path):
    path = config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("[other]\nx = 1\n")

    root = tmp_path / "mystore"
    write_config_data_root(root)

    config = read_config()
    assert config["other"]["x"] == 1
    assert config["store"]["data_root"] == root.as_posix()


def test_write_config_data_root_replaces_existing_store_table(tmp_path):
    first_root = tmp_path / "first"
    write_config_data_root(first_root)

    second_root = tmp_path / "second"
    write_config_data_root(second_root)

    config = read_config()
    assert config["store"]["data_root"] == second_root.as_posix()


def test_write_config_data_root_escapes_quotes_and_backslashes(tmp_path):
    # A path component containing a double quote and a backslash: both are
    # legal filesystem characters on macOS/Linux, and a naive f-string write
    # (`data_root = "{path}"`) would produce invalid, or worse silently
    # truncated, TOML for a path like this.
    tricky_name = 'weird"quote\\slash'
    root = tmp_path / tricky_name
    os.makedirs(root)

    write_config_data_root(root)

    config = read_config()
    assert config["store"]["data_root"] == root.as_posix()


# --- resolve_store_root precedence ------------------------------------------


def test_resolve_returns_none_when_nothing_set():
    assert resolve_store_root(explicit=None, registry_root=None) is None


def test_resolve_uses_explicit_over_everything(tmp_path, monkeypatch):
    explicit_root = tmp_path / "explicit"
    init_store(explicit_root)

    env_root = tmp_path / "env"
    init_store(env_root)
    monkeypatch.setenv(STORE_ENV, str(env_root))

    registry_root = tmp_path / "registry"
    init_store(registry_root)

    config_root = tmp_path / "config"
    init_store(config_root)
    write_config_data_root(config_root)

    resolved = resolve_store_root(explicit=str(explicit_root), registry_root=str(registry_root))
    assert resolved == explicit_root.resolve()


def test_resolve_uses_env_over_registry_and_config(tmp_path, monkeypatch):
    env_root = tmp_path / "env"
    init_store(env_root)
    monkeypatch.setenv(STORE_ENV, str(env_root))

    registry_root = tmp_path / "registry"
    init_store(registry_root)

    config_root = tmp_path / "config"
    init_store(config_root)
    write_config_data_root(config_root)

    resolved = resolve_store_root(explicit=None, registry_root=str(registry_root))
    assert resolved == env_root.resolve()


def test_resolve_treats_empty_env_as_unset(tmp_path, monkeypatch):
    # An empty PREPARE_DATA (e.g. exported but never assigned) must not
    # win over a real registry or config root, and must not itself raise.
    monkeypatch.setenv(STORE_ENV, "")

    registry_root = tmp_path / "registry"
    init_store(registry_root)

    resolved = resolve_store_root(explicit=None, registry_root=str(registry_root))
    assert resolved == registry_root.resolve()


def test_resolve_uses_registry_over_config(tmp_path):
    registry_root = tmp_path / "registry"
    init_store(registry_root)

    config_root = tmp_path / "config"
    init_store(config_root)
    write_config_data_root(config_root)

    resolved = resolve_store_root(explicit=None, registry_root=str(registry_root))
    assert resolved == registry_root.resolve()


def test_resolve_falls_back_to_config(tmp_path):
    config_root = tmp_path / "config"
    init_store(config_root)
    write_config_data_root(config_root)

    resolved = resolve_store_root(explicit=None, registry_root=None)
    assert resolved == config_root.resolve()


# --- resolve_store_root marker validation ------------------------------------


def test_resolve_raises_naming_explicit_rule_when_marker_missing(tmp_path):
    missing_root = tmp_path / "no-marker"
    missing_root.mkdir()

    with pytest.raises(DataAccessError, match="--data-root"):
        resolve_store_root(explicit=str(missing_root), registry_root=None)


def test_resolve_raises_naming_env_rule_when_marker_missing(tmp_path, monkeypatch):
    missing_root = tmp_path / "no-marker"
    missing_root.mkdir()
    monkeypatch.setenv(STORE_ENV, str(missing_root))

    with pytest.raises(DataAccessError, match=STORE_ENV):
        resolve_store_root(explicit=None, registry_root=None)


def test_resolve_raises_naming_registry_rule_when_marker_missing(tmp_path):
    missing_root = tmp_path / "no-marker"
    missing_root.mkdir()

    with pytest.raises(DataAccessError, match="registry store.root"):
        resolve_store_root(explicit=None, registry_root=str(missing_root))


def test_resolve_raises_naming_config_rule_when_marker_missing(tmp_path):
    missing_root = tmp_path / "no-marker"
    missing_root.mkdir()
    write_config_data_root(missing_root)

    with pytest.raises(DataAccessError, match="config data_root"):
        resolve_store_root(explicit=None, registry_root=None)


def test_resolve_skips_marker_check_when_require_marker_false(tmp_path):
    missing_root = tmp_path / "no-marker"
    missing_root.mkdir()

    resolved = resolve_store_root(explicit=str(missing_root), registry_root=None, require_marker=False)
    assert resolved == missing_root.resolve()


# --- add_allowed_root registration ------------------------------------------


def test_resolve_registers_root_with_secure_subprocess(tmp_path):
    root = tmp_path / "mystore"
    init_store(root)

    resolved = resolve_store_root(explicit=str(root), registry_root=None)

    assert resolved in SecureSubprocess.allowed_roots()


def test_resolve_refuses_a_store_entry_that_is_not_a_table(tmp_path):
    config_path().parent.mkdir(parents=True, exist_ok=True)
    config_path().write_text('store = "x"\n')
    with pytest.raises(ConfigurationError, match="store"):
        resolve_store_root(None, None)


def test_resolve_refuses_a_data_root_that_is_not_a_string(tmp_path):
    config_path().parent.mkdir(parents=True, exist_ok=True)
    config_path().write_text("[store]\ndata_root = 5\n")
    with pytest.raises(ConfigurationError, match="data_root"):
        resolve_store_root(None, None)


def test_resolve_ignores_a_broken_config_when_the_explicit_root_is_found(tmp_path):
    root = tmp_path / "store"
    init_store(root)
    config_path().parent.mkdir(parents=True, exist_ok=True)
    config_path().write_text("not toml = = =\n")
    assert resolve_store_root(str(root), None) == root.resolve()


def test_resolve_expands_a_home_relative_candidate(tmp_path):
    init_store(tmp_path / "store")
    assert resolve_store_root("~/store", None) == (tmp_path / "store").resolve()
