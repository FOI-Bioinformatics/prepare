"""
Shared test fixtures and configuration for PREPARE tests.
"""

import os

import matplotlib.pyplot as plt
import pytest

from prepare.core import settings
from prepare.core.constants import STORE_ENV
from _sra_fixtures import fast_project_lock_poll  # noqa: F401 - a fixture for tests/test_data_sra_download.py

# Read once when the perf module is imported, not a runtime setting; left for CI to set.
_KEPT_PREPARE_VARIABLES = {"PREPARE_PERF_SCALE"}


@pytest.fixture(autouse=True)
def isolate_runtime_settings(monkeypatch):
    """Start every test from the default runtime settings.

    ``prepare.core.settings`` reads ``PREPARE_*`` variables (and ``NCBI_API_KEY``) and caches
    the result; a developer's own shell values must not leak into a test, and a test's
    ``activate`` must not leak into the next one.
    """
    for variable in list(os.environ):
        if variable.startswith(settings.ENV_PREFIX) and variable not in _KEPT_PREPARE_VARIABLES:
            monkeypatch.delenv(variable, raising=False)
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    # The download free-space guard reads the host's real free space; a test that wants it sets
    # --min-free-gb (or deletes this variable) itself, so no result depends on the machine's disk.
    monkeypatch.setenv(settings.SETTINGS["min_free_gb"].env, "0")
    # No test writes a run log under its project unless it sets PREPARE_RUN_LOG itself.
    monkeypatch.setenv(settings.SETTINGS["run_log"].env, "false")
    settings.reset_for_tests()
    yield
    settings.reset_for_tests()


@pytest.fixture(autouse=True)
def isolate_store_discovery(tmp_path, monkeypatch):
    """Keep every test away from the maintainer's own shared data store.

    A store root is discovered from ``PREPARE_DATA``, a project registry, or the user
    config file at ``$XDG_CONFIG_HOME/prepare/config.toml`` (``~/.config`` without it).
    The documented migration tells a user to run ``store_init --set-default``, after which an
    unprotected test suite would open that real store's catalogue, write to it, and read
    whatever it holds. Pointing HOME and XDG_CONFIG_HOME at this test's own tmp_path and
    dropping PREPARE_DATA makes the discovery rules find nothing but what a test sets up.
    """
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    monkeypatch.delenv(STORE_ENV, raising=False)
    yield


@pytest.fixture(autouse=True)
def forget_probed_tool_versions():
    """The fasterq-dump version is probed once per process; each test starts without that memory."""
    from prepare.data.sra.accession import fasterq_dump_version

    fasterq_dump_version.cache_clear()
    yield
    fasterq_dump_version.cache_clear()


@pytest.fixture(autouse=True)
def close_matplotlib_figures():
    """Close all matplotlib figures after each test to prevent memory leaks."""
    yield
    plt.close("all")
