"""
Tests for prepare.core.settings: one resolution of every runtime setting.

The autouse fixtures in conftest.py point HOME and XDG_CONFIG_HOME at tmp_path, clear every
PREPARE_* variable and reset the cached settings, so each test starts from the defaults.
"""

import argparse
import dataclasses
import importlib
import logging
from unittest.mock import patch

import pytest

from prepare.core import settings
from prepare.core.constants import (
    CATALOG_LOCK_WAIT_SECONDS,
    DATASET_LOCK_STALE_SECONDS,
    LOCK_HEARTBEAT_SECONDS,
)
from prepare.core.exceptions import ConfigurationError


def _write_config(tmp_path, text):
    """Write ``text`` as the user config file under this test's XDG_CONFIG_HOME."""
    path = tmp_path / "config" / "prepare" / "config.toml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


# --- precedence, one test per level ----------------------------------------


def test_default_level():
    resolved = settings.resolve_setting("progress_every")
    assert resolved == settings.Resolved("progress_every", 50, "default")


def test_config_level(tmp_path):
    _write_config(tmp_path, "[runtime]\nprogress_every = 20\n")
    resolved = settings.resolve_setting("progress_every")
    assert resolved.value == 20
    assert resolved.source == "config [runtime] progress_every"


def test_environment_beats_config(tmp_path, monkeypatch):
    _write_config(tmp_path, "[runtime]\nprogress_every = 20\n")
    monkeypatch.setenv("PREPARE_PROGRESS_EVERY", "30")
    resolved = settings.resolve_setting("progress_every")
    assert resolved.value == 30
    assert resolved.source == "PREPARE_PROGRESS_EVERY"


def test_flag_beats_environment_and_config(tmp_path, monkeypatch):
    _write_config(tmp_path, "[runtime]\ntimeout = 20\n")
    monkeypatch.setenv("PREPARE_TIMEOUT", "30")
    resolved = settings.resolve_setting("subprocess_timeout", cli_value=5.0)
    assert resolved == settings.Resolved("subprocess_timeout", 5.0, "--timeout")


def test_timeout_uses_the_short_names(tmp_path, monkeypatch):
    _write_config(tmp_path, "[runtime]\ntimeout = 3600\n")
    assert settings.resolve_setting("subprocess_timeout").source == "config [runtime] timeout"
    monkeypatch.setenv("PREPARE_TIMEOUT", "0")
    resolved = settings.resolve_setting("subprocess_timeout")
    assert (resolved.value, resolved.source) == (0.0, "PREPARE_TIMEOUT")


def test_empty_environment_value_counts_as_unset(monkeypatch):
    monkeypatch.setenv("PREPARE_PROGRESS_EVERY", "")
    assert settings.resolve_setting("progress_every").source == "default"


def test_ncbi_api_key_keeps_the_ncbi_variable(monkeypatch):
    monkeypatch.setenv("NCBI_API_KEY", "plain")
    assert settings.resolve_setting("ncbi_api_key") == settings.Resolved("ncbi_api_key", "plain", "NCBI_API_KEY")
    monkeypatch.setenv("PREPARE_NCBI_API_KEY", "prefixed")
    resolved = settings.resolve_setting("ncbi_api_key")
    assert (resolved.value, resolved.source) == ("prefixed", "PREPARE_NCBI_API_KEY")


def test_boolean_from_environment_and_native_toml(tmp_path, monkeypatch):
    _write_config(tmp_path, "[runtime]\nlog_host = true\n")
    assert settings.resolve_setting("log_host").value is True
    monkeypatch.setenv("PREPARE_LOG_HOST", "no")
    assert settings.resolve_setting("log_host").value is False


def test_every_default_is_the_documented_one(monkeypatch):
    monkeypatch.delenv("PREPARE_MIN_FREE_GB")  # set to 0 by conftest for every other test
    values = {name: resolved.value for name, resolved in settings.resolve_all().items()}
    assert values["subprocess_timeout"] == 0.0
    assert values["max_workers_cap"] == 4
    assert values["lock_wait"] == 0.0
    assert values["progress_every"] == 50
    assert values["log_host"] is False
    assert values["min_free_gb"] == 10.0
    assert values["assembly_memory"] == "auto"
    assert values["log_level"] == "INFO"
    assert values["ncbi_email"] is None and values["ncbi_api_key"] is None
    assert values["dataset_lock_stale"] == DATASET_LOCK_STALE_SECONDS
    assert values["lock_heartbeat"] == LOCK_HEARTBEAT_SECONDS
    assert values["catalog_lock_wait"] == CATALOG_LOCK_WAIT_SECONDS


def test_registry_defaults_match_the_registry_module():
    from prepare.data import registry

    assert settings.resolve_setting("registry_lock_wait").value == registry.LOCK_WAIT_SECONDS
    assert settings.resolve_setting("registry_lock_stale").value == registry.LOCK_STALE_SECONDS


def test_runtime_settings_has_one_field_per_setting():
    fields = {f.name for f in dataclasses.fields(settings.RuntimeSettings)} - {"sources", "warnings"}
    assert fields == set(settings.SETTINGS)


def test_unknown_setting_name_is_a_programming_error():
    with pytest.raises(KeyError):
        settings.resolve_setting("no_such_setting")


# --- bad values name their source -------------------------------------------


@pytest.mark.parametrize(
    "variable, value",
    [
        ("PREPARE_PROGRESS_EVERY", "many"),
        ("PREPARE_PROGRESS_EVERY", "-1"),
        ("PREPARE_TIMEOUT", "-5"),
        ("PREPARE_TIMEOUT", "nan"),
        ("PREPARE_LOG_HOST", "perhaps"),
        ("PREPARE_LOG_LEVEL", "LOUD"),
        ("PREPARE_MAX_WORKERS_CAP", "0"),
        ("PREPARE_NCBI_EMAIL", "not-an-address"),
        ("PREPARE_ASSEMBLY_MEMORY", "lots"),
        ("PREPARE_REGISTRY_LOCK_WAIT", "0"),
    ],
)
def test_bad_environment_value_names_the_variable(monkeypatch, variable, value):
    monkeypatch.setenv(variable, value)
    with pytest.raises(ConfigurationError, match=variable):
        settings.resolve_all()


def test_bad_config_value_names_the_key_and_the_file(tmp_path):
    path = _write_config(tmp_path, '[runtime]\nprogress_every = "often"\n')
    with pytest.raises(ConfigurationError) as excinfo:
        settings.resolve_setting("progress_every")
    assert "config [runtime] progress_every" in str(excinfo.value)
    assert str(path) in str(excinfo.value)


def test_config_value_of_the_wrong_type_is_refused(tmp_path):
    _write_config(tmp_path, "[runtime]\nlog_host = [1, 2]\n")
    with pytest.raises(ConfigurationError, match="log_host"):
        settings.resolve_setting("log_host")


@pytest.mark.parametrize("value", ["auto", "0.5", "1", "32G", "32000M", "2000000", "16g"])
def test_assembly_memory_accepts_documented_forms(monkeypatch, value):
    monkeypatch.setenv("PREPARE_ASSEMBLY_MEMORY", value)
    assert settings.resolve_setting("assembly_memory").value == value


@pytest.mark.parametrize("value", ["2.5", "0", "lots", "-1G"])
def test_assembly_memory_rejects_other_forms(monkeypatch, value):
    # "2.5" is neither a fraction (at most 1) nor a whole number of bytes.
    monkeypatch.setenv("PREPARE_ASSEMBLY_MEMORY", value)
    with pytest.raises(ConfigurationError, match="PREPARE_ASSEMBLY_MEMORY"):
        settings.resolve_setting("assembly_memory")


@pytest.mark.parametrize("value", ["500", ".0", "0.0"])
def test_assembly_memory_agrees_with_the_assembly_parser(monkeypatch, value):
    # 500 bytes is below the 1M floor of resources.parse_memory, so startup refuses it too.
    monkeypatch.setenv("PREPARE_ASSEMBLY_MEMORY", value)
    with pytest.raises(ConfigurationError, match="PREPARE_ASSEMBLY_MEMORY"):
        settings.resolve_setting("assembly_memory")


def test_heartbeat_not_below_the_stale_limit_is_refused(monkeypatch):
    monkeypatch.setenv("PREPARE_LOCK_HEARTBEAT", "600")
    with pytest.raises(ConfigurationError, match="lock_heartbeat"):
        settings.activate(None)


def test_registry_stale_limit_must_exceed_its_heartbeat(monkeypatch):
    monkeypatch.setenv("PREPARE_REGISTRY_LOCK_STALE", "1")
    with pytest.raises(ConfigurationError, match="registry_lock_stale"):
        settings.activate(None)


# --- the config file ---------------------------------------------------------


def test_malformed_toml_names_the_file(tmp_path):
    path = _write_config(tmp_path, "[runtime\nprogress_every = 1\n")
    with pytest.raises(ConfigurationError, match=str(path)):
        settings.resolve_setting("progress_every")


def test_malformed_toml_through_the_store_reexport(tmp_path):
    from prepare.store.resolve import read_config

    path = _write_config(tmp_path, "not toml at all = = =\n")
    with pytest.raises(ConfigurationError, match=str(path)):
        read_config()


def test_unreadable_config_names_the_file(tmp_path):
    path = _write_config(tmp_path, "")
    path.unlink()
    path.mkdir()  # a folder where the file should be: opening it raises an OSError
    with pytest.raises(ConfigurationError, match=str(path)):
        settings.read_config()


def test_runtime_that_is_not_a_table_is_refused(tmp_path):
    _write_config(tmp_path, 'runtime = "fast"\n')
    with pytest.raises(ConfigurationError, match=r"\[runtime\]"):
        settings.resolve_all()


def test_unknown_runtime_key_is_kept_for_the_caller(tmp_path, caplog):
    _write_config(tmp_path, "[runtime]\nprogres_every = 5\n")
    with caplog.at_level(logging.WARNING, logger="prepare.core.settings"):
        runtime = settings.activate(None)
    assert runtime.source("progress_every") == "default"
    assert len(runtime.warnings) == 1 and "'progres_every'" in runtime.warnings[0]
    assert caplog.text == ""  # activate runs before logging is set up; main() logs the warnings


def test_unknown_runtime_key_is_logged_by_active_without_activate(tmp_path, caplog):
    _write_config(tmp_path, "[runtime]\nprogres_every = 5\n")
    with caplog.at_level(logging.WARNING, logger="prepare.core.settings"):
        settings.active()
    assert "progres_every" in caplog.text


def test_main_reports_an_unknown_config_key_on_stderr(tmp_path, monkeypatch, capsys):
    """The warning reaches the stderr handler setup_logging installs, not only a caplog handler."""
    from prepare.cli.main import main

    _write_config(tmp_path, "[runtime]\nprogres_every = 5\n")
    monkeypatch.chdir(tmp_path)
    root = logging.getLogger()
    level = root.level
    try:
        main(["blacklist", "--list", "--registry", str(tmp_path / "prepare_registry.json")])
    finally:
        for handler in root.handlers[:]:
            if getattr(handler, "_prepare", False):
                root.removeHandler(handler)
        root.setLevel(level)
    err = capsys.readouterr().err
    assert "WARNING" in err and "unknown key 'progres_every' in [runtime] is ignored" in err


def test_other_tables_are_ignored(tmp_path, caplog):
    _write_config(tmp_path, '[store]\ndata_root = "/nowhere"\n')
    with caplog.at_level(logging.WARNING, logger="prepare.core.settings"):
        settings.resolve_all()
    assert caplog.text == ""


def test_xdg_config_home_is_honoured(tmp_path, monkeypatch):
    home_config = tmp_path / "home" / ".config" / "prepare" / "config.toml"
    home_config.parent.mkdir(parents=True)
    home_config.write_text("[runtime]\nprogress_every = 7\n")
    _write_config(tmp_path, "[runtime]\nprogress_every = 9\n")
    assert settings.config_path() == tmp_path / "config" / "prepare" / "config.toml"
    assert settings.resolve_setting("progress_every").value == 9

    monkeypatch.delenv("XDG_CONFIG_HOME")
    assert settings.config_path() == home_config
    assert settings.resolve_setting("progress_every").value == 7


# --- activate, active, reset ------------------------------------------------


def test_activate_takes_flags_from_the_namespace(monkeypatch):
    monkeypatch.setenv("PREPARE_LOCK_WAIT", "5")
    runtime = settings.activate(argparse.Namespace(lock_wait=12.0, email=None))
    assert runtime.lock_wait == 12.0
    assert runtime.sources["lock_wait"] == "--lock-wait"
    assert runtime.ncbi_email is None
    assert settings.active() is runtime


def test_a_flag_left_at_none_falls_through(monkeypatch):
    monkeypatch.setenv("PREPARE_LOCK_WAIT", "5")
    runtime = settings.activate(argparse.Namespace(lock_wait=None))
    assert (runtime.lock_wait, runtime.sources["lock_wait"]) == (5.0, "PREPARE_LOCK_WAIT")


def test_active_without_activate_builds_from_environment(monkeypatch):
    monkeypatch.setenv("PREPARE_PROGRESS_EVERY", "3")
    assert settings.active().progress_every == 3
    monkeypatch.setenv("PREPARE_PROGRESS_EVERY", "4")
    assert settings.active().progress_every == 3  # cached
    settings.reset_for_tests()
    assert settings.active().progress_every == 4


def test_describe_hides_the_api_key(monkeypatch):
    monkeypatch.setenv("NCBI_API_KEY", "secret-value")
    lines = settings.active().describe()
    assert any(line.startswith("ncbi_api_key = ") and "NCBI_API_KEY" in line for line in lines)
    assert not any("secret-value" in line for line in lines)
    assert any(line == "progress_every = 50 (default)" for line in lines)


def test_setting_for_prefers_the_flag_then_the_active_settings(monkeypatch):
    monkeypatch.setenv("PREPARE_TEMP_FOLDER", "/scratch/tmp")
    assert settings.setting_for(argparse.Namespace(temp_folder="/given"), "temp_folder") == "/given"
    assert settings.setting_for(argparse.Namespace(temp_folder=None), "temp_folder") == "/scratch/tmp"
    assert settings.setting_for(argparse.Namespace(), "temp_folder") == "/scratch/tmp"


def test_require_email_raises_when_nothing_resolves():
    with pytest.raises(ConfigurationError, match="--email"):
        settings.require_email(argparse.Namespace(email=None))


def test_require_email_from_environment(monkeypatch):
    monkeypatch.setenv("PREPARE_NCBI_EMAIL", "someone@example.org")
    assert settings.require_email(argparse.Namespace(email=None)) == "someone@example.org"


# --- lock policies read the settings, module attributes stay the fallback ----


def test_setting_or_uses_the_fallback_for_a_default():
    assert settings.setting_or("dataset_lock_stale", 0.3) == 0.3


def test_setting_or_uses_a_value_the_user_set(monkeypatch):
    monkeypatch.setenv("PREPARE_DATASET_LOCK_STALE", "900")
    assert settings.setting_or("dataset_lock_stale", 0.3) == 900.0


def test_registry_lock_reads_the_settings(tmp_path, monkeypatch):
    from prepare.data import registry

    seen = []

    class _Held:
        def __init__(self, lock, policy):
            seen.append(policy)

        def __enter__(self):
            return tmp_path / "lock"

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(registry, "held_lock", _Held)
    monkeypatch.setattr(registry, "LOCK_STALE_SECONDS", 0.3)
    with registry._acquire_lock(tmp_path / "lock"):
        pass
    monkeypatch.setenv("PREPARE_REGISTRY_LOCK_WAIT", "45")
    settings.reset_for_tests()
    with registry._acquire_lock(tmp_path / "lock"):
        pass
    assert (seen[0].stale_seconds, seen[0].wait_seconds) == (0.3, registry.LOCK_WAIT_SECONDS)
    assert (seen[1].stale_seconds, seen[1].wait_seconds) == (0.3, 45.0)


def test_dataset_and_catalog_policies_read_the_settings(tmp_path, monkeypatch):
    from prepare.data.sra import accession
    from prepare.store import locks

    monkeypatch.setattr(locks, "DATASET_LOCK_STALE_SECONDS", 0.3)
    assert locks._dataset_policy("SRR1", 0.0).stale_seconds == 0.3
    monkeypatch.setenv("PREPARE_DATASET_LOCK_STALE", "900")
    monkeypatch.setenv("PREPARE_LOCK_HEARTBEAT", "20")
    settings.reset_for_tests()
    policy = locks._dataset_policy("SRR1", 0.0)
    assert (policy.stale_seconds, policy.heartbeat_seconds) == (900.0, 20.0)
    policy = accession._project_lock_policy("SRR1", 0.0)
    assert (policy.stale_seconds, policy.heartbeat_seconds) == (900.0, 20.0)


def test_catalog_write_reads_the_catalog_wait(tmp_path, monkeypatch):
    from prepare.store import catalog
    from prepare.store.layout import init_store, store_paths

    seen = []
    real_held_lock = catalog.held_lock

    def _spy(lock, policy, **kwargs):
        seen.append(policy)
        return real_held_lock(lock, policy, **kwargs)

    monkeypatch.setattr(catalog, "held_lock", _spy)
    monkeypatch.setenv("PREPARE_CATALOG_LOCK_WAIT", "7")
    init_store(tmp_path / "store")
    with catalog.catalog_write(store_paths(tmp_path / "store")):
        pass
    assert seen[-1].wait_seconds == 7.0


# --- the command line ----------------------------------------------------------


def test_main_resolves_the_log_level_from_the_environment(monkeypatch):
    from prepare.cli.main import main

    monkeypatch.setenv("PREPARE_LOG_LEVEL", "warning")
    with patch("prepare.cli.main.setup_logging") as setup:
        assert main(["store_status", "--data-root", "/nonexistent-store-root"]) == 1
    setup.assert_called_once()
    assert setup.call_args.kwargs["level"] == logging.WARNING
    assert settings.active().sources["log_level"] == "PREPARE_LOG_LEVEL"


def test_main_logs_the_settings_with_their_sources(monkeypatch, caplog):
    from prepare.cli.main import main

    monkeypatch.setenv("PREPARE_PROGRESS_EVERY", "25")
    with patch("prepare.cli.main.setup_logging"), caplog.at_level(logging.DEBUG):
        main(["--log-level", "DEBUG", "store_status", "--data-root", "/nonexistent-store-root"])
    assert "progress_every = 25 (PREPARE_PROGRESS_EVERY)" in caplog.text
    assert "log_level = 'DEBUG' (--log-level)" in caplog.text


def test_main_refuses_a_malformed_config(tmp_path, caplog):
    from prepare.cli.main import main

    path = _write_config(tmp_path, "[runtime\n")
    with patch("prepare.cli.main.setup_logging"), caplog.at_level(logging.ERROR):
        assert main(["store_status"]) == 3  # ConfigurationError
    assert str(path) in caplog.text


def test_download_metadata_without_an_email_is_refused(tmp_path, caplog):
    from prepare.cli.main import main

    with patch("prepare.cli.commands.metadata.download_metadata") as download, caplog.at_level(logging.ERROR):
        assert main(["download_metadata", "--matches-folder", str(tmp_path)]) == 3  # ConfigurationError
    download.assert_not_called()
    assert "PREPARE_NCBI_EMAIL" in caplog.text


def test_download_metadata_takes_email_and_key_from_settings(tmp_path, monkeypatch):
    from prepare.cli.commands.metadata import DownloadMetadataCommand

    monkeypatch.setenv("PREPARE_NCBI_EMAIL", "someone@example.org")
    monkeypatch.setenv("NCBI_API_KEY", "env-key")
    parser = argparse.ArgumentParser()
    DownloadMetadataCommand().configure_parser(parser)
    args = parser.parse_args(["--matches-folder", str(tmp_path), "--dry-run"])
    assert args.email is None and args.api_key is None
    with patch("prepare.cli.commands.metadata.download_metadata", return_value={}) as download:
        assert DownloadMetadataCommand().execute(args) == 0
    assert download.call_args.kwargs["email"] == "someone@example.org"
    assert download.call_args.kwargs["api_key"] == "env-key"


@pytest.mark.parametrize(
    "module, class_name, argv",
    [
        ("prepare.cli.commands.sra_enhanced", "SRAInfoCommand", ["--accessions-file", "a.txt"]),
        ("prepare.cli.commands.advanced_analysis", "TaxonomyValidationCommand", ["--species-file", "s.txt"]),
    ],
)
def test_other_ncbi_commands_accept_a_missing_email_flag(module, class_name, argv):
    command = getattr(importlib.import_module(module), class_name)()
    parser = argparse.ArgumentParser()
    command.configure_parser(parser)
    args = parser.parse_args(argv)
    assert args.email is None
    with pytest.raises(ConfigurationError, match="--email"):
        command.execute(args)


def test_lock_wait_flags_default_to_none_and_resolve_through_settings(monkeypatch):
    from prepare.cli.commands.sra import DownloadSraCommand
    from prepare.cli.commands.store.adopt import StoreAdoptCommand
    from prepare.data.registry import Registry

    for command in (DownloadSraCommand(), StoreAdoptCommand()):
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)
        assert parser.get_default("lock_wait") is None

    monkeypatch.setenv("PREPARE_LOCK_WAIT", "5")
    options = DownloadSraCommand()._store_options(argparse.Namespace(lock_wait=None), None, Registry())
    assert options == {"lock_wait": 5.0}
    options = DownloadSraCommand()._store_options(argparse.Namespace(lock_wait=2.0), None, Registry())
    assert options == {"lock_wait": 2.0}


# --- fix wave: flag values, secrets, TOML line, first active() call -------------


@pytest.mark.parametrize(
    "dest, value, flag",
    [("email", "nope", "--email"), ("min_free_gb", -1.0, "--min-free-gb"), ("log_level", "LOUD", "--log-level")],
)
def test_a_flag_value_is_checked_like_an_environment_value(dest, value, flag):
    with pytest.raises(ConfigurationError, match=f"from {flag}"):
        settings.activate(argparse.Namespace(**{dest: value}))


def test_a_flag_value_is_parsed_by_setting_for():
    assert settings.setting_for(argparse.Namespace(log_level="debug"), "log_level") == "DEBUG"
    with pytest.raises(ConfigurationError, match="--min-free-gb"):
        settings.setting_for(argparse.Namespace(min_free_gb=-2.0), "min_free_gb")


def test_a_bad_secret_value_is_not_echoed(tmp_path):
    _write_config(tmp_path, '[runtime]\nncbi_api_key = ["SECRETKEY"]\n')
    with pytest.raises(ConfigurationError) as excinfo:
        settings.activate(None)
    assert "SECRETKEY" not in str(excinfo.value)
    assert "(hidden)" in str(excinfo.value)


def test_malformed_toml_quotes_the_offending_line(tmp_path):
    _write_config(tmp_path, "# settings\n[runtime\nprogress_every = 1\n")
    with pytest.raises(ConfigurationError, match=r"line 2 is '\[runtime'"):
        settings.resolve_setting("progress_every")


def test_concurrent_first_calls_build_the_settings_once(monkeypatch):
    import threading
    import time

    built = []
    real_build = settings._build

    def slow_build(args):
        built.append(threading.get_ident())
        time.sleep(0.05)
        return real_build(args)

    monkeypatch.setattr(settings, "_build", slow_build)
    results = []
    threads = [threading.Thread(target=lambda: results.append(settings.active())) for _ in range(8)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert len(built) == 1
    assert all(result is results[0] for result in results)


# --- prefetch_max_size --------------------------------------------------------


def test_prefetch_max_size_defaults_to_the_constant():
    from prepare.core.constants import DEFAULT_PREFETCH_MAX_SIZE

    assert DEFAULT_PREFETCH_MAX_SIZE == "100G"
    assert settings.resolve_setting("prefetch_max_size").value == "100G"


@pytest.mark.parametrize("value", ["100G", "20g", "500M", "1T", "4096k", "1000000"])
def test_prefetch_max_size_accepts_a_size(monkeypatch, value):
    monkeypatch.setenv("PREPARE_PREFETCH_MAX_SIZE", value)
    assert settings.resolve_setting("prefetch_max_size").value == value


@pytest.mark.parametrize("value", ["lots", "1.5G", "-1G", "10GB", "G", "10 G"])
def test_prefetch_max_size_rejects_other_forms(monkeypatch, value):
    monkeypatch.setenv("PREPARE_PREFETCH_MAX_SIZE", value)
    with pytest.raises(ConfigurationError, match="PREPARE_PREFETCH_MAX_SIZE"):
        settings.resolve_setting("prefetch_max_size")


def test_prefetch_max_size_from_the_config_file(tmp_path):
    _write_config(tmp_path, '[runtime]\nprefetch_max_size = "50G"\n')
    resolved = settings.resolve_setting("prefetch_max_size")
    assert resolved.value == "50G"
    assert resolved.source == "config [runtime] prefetch_max_size"
