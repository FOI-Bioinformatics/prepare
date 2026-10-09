"""Tests for module-level mutable state that must survive two runs in one process.

Covers the Task 10 fixes: prepare no longer configures logging on import (a library host
or a second test run must not have the root logger reconfigured out from under it),
SecureSubprocess._extra_roots and the GTDB session are safe under concurrent access, and
Bio.Entrez's email/api_key module attributes are set only for the duration of one call.
"""

import importlib
import logging
import threading

import pytest
from Bio import Entrez

import prepare
from prepare.data import gtdb
from prepare.data.metadata import _entrez_credentials
from prepare.utils.logging import setup_logging
from prepare.utils.security import SecureSubprocess


def _prepare_tagged_handlers():
    """The handlers on the root logger that setup_logging installed."""
    return [h for h in logging.getLogger().handlers if getattr(h, "_prepare", False)]


class TestImportDoesNotConfigureLogging:
    """Importing prepare must not call setup_logging() or touch the root logger."""

    def test_import_leaves_a_preexisting_root_handler_in_place(self):
        root_logger = logging.getLogger()
        sentinel = logging.NullHandler()
        root_logger.addHandler(sentinel)
        try:
            importlib.reload(prepare)
            assert sentinel in root_logger.handlers
        finally:
            root_logger.removeHandler(sentinel)

    def test_prepare_logger_has_a_null_handler(self):
        importlib.reload(prepare)
        prepare_logger = logging.getLogger("prepare")
        assert any(isinstance(h, logging.NullHandler) for h in prepare_logger.handlers)


class TestSetupLoggingTagsItsOwnHandlers:
    """setup_logging only ever removes the handlers it installed itself."""

    def teardown_method(self):
        for handler in _prepare_tagged_handlers():
            logging.getLogger().removeHandler(handler)

    def test_called_twice_leaves_one_console_handler(self):
        setup_logging(level=logging.INFO)
        setup_logging(level=logging.INFO)
        assert len(_prepare_tagged_handlers()) == 1

    def test_caplog_still_captures_after_setup_logging(self, caplog):
        setup_logging(level=logging.INFO)
        with caplog.at_level(logging.INFO):
            logging.getLogger("prepare.test_globals").info("hello from setup_logging test")
        assert "hello from setup_logging test" in caplog.text

    def test_a_foreign_root_handler_survives_setup_logging(self):
        root_logger = logging.getLogger()
        foreign = logging.NullHandler()
        root_logger.addHandler(foreign)
        try:
            setup_logging(level=logging.INFO)
            assert foreign in root_logger.handlers
        finally:
            root_logger.removeHandler(foreign)


class TestAddAllowedRootIsThreadSafe:
    """20 threads registering the same root concurrently leave exactly one entry."""

    def test_concurrent_add_allowed_root_leaves_one_entry(self, tmp_path):
        target = tmp_path / "concurrent-root"
        target.mkdir()

        threads = [threading.Thread(target=SecureSubprocess.add_allowed_root, args=(target,)) for _ in range(20)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        resolved = target.resolve()
        assert SecureSubprocess._extra_roots.count(resolved) == 1
        assert SecureSubprocess.allowed_roots().count(resolved) == 1


class TestEntrezCredentialsContextManager:
    """_entrez_credentials sets Entrez.email/api_key for the call and always restores them."""

    def setup_method(self):
        self._previous_email = Entrez.email
        self._previous_api_key = Entrez.api_key

    def teardown_method(self):
        Entrez.email = self._previous_email
        Entrez.api_key = self._previous_api_key

    def test_sets_and_restores_after_normal_use(self):
        Entrez.email = "before@example.com"
        Entrez.api_key = "before-key"

        with _entrez_credentials("during@example.com", "during-key"):
            assert Entrez.email == "during@example.com"
            assert Entrez.api_key == "during-key"

        assert Entrez.email == "before@example.com"
        assert Entrez.api_key == "before-key"

    def test_restores_after_an_exception(self):
        Entrez.email = "before@example.com"
        Entrez.api_key = "before-key"

        with pytest.raises(ValueError):
            with _entrez_credentials("during@example.com", "during-key"):
                raise ValueError("simulated efetch failure")

        assert Entrez.email == "before@example.com"
        assert Entrez.api_key == "before-key"


class TestGetSessionIsThreadSafe:
    """20 threads calling get_session() concurrently all get the one process-wide session."""

    def setup_method(self):
        gtdb._session = None

    def teardown_method(self):
        gtdb._session = None

    def test_concurrent_get_session_returns_one_shared_object(self):
        results = []
        results_lock = threading.Lock()

        def worker():
            session = gtdb.get_session()
            with results_lock:
                results.append(session)

        threads = [threading.Thread(target=worker) for _ in range(20)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        assert len(results) == 20
        assert len({id(session) for session in results}) == 1
