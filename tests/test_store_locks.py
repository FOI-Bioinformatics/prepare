"""Tests for the long-lived per-accession dataset lock."""

import json
import logging
import os
import threading
import time

import pytest

from prepare.core.exceptions import DataAccessError
from prepare.store import locks as locks_module
from prepare.store.layout import init_store, lock_path
from prepare.store.locks import dataset_lock, lock_holder, lock_is_held


@pytest.fixture
def paths(tmp_path):
    return init_store(tmp_path / "store")


@pytest.fixture
def fast_lock(monkeypatch):
    """Shrink every lock interval so a contention test runs in about a second."""
    monkeypatch.setattr(locks_module, "LOCK_HEARTBEAT_SECONDS", 0.05)
    monkeypatch.setattr(locks_module, "DATASET_LOCK_STALE_SECONDS", 0.3)
    monkeypatch.setattr(locks_module, "LOCK_POLL_SECONDS", 0.02)
    monkeypatch.setattr(locks_module, "LOCK_WAIT_LOG_SECONDS", 0.2)


class TestDatasetLockContention:
    def test_waiter_does_not_take_over_a_heartbeating_lock(self, paths, fast_lock):
        """A holder that keeps its heartbeat running is never displaced, even when it holds
        the lock for far longer than the stale threshold."""
        events = []
        holder_inside = threading.Event()

        def holder():
            with dataset_lock(paths, "SRR1"):
                holder_inside.set()
                # Held for well over three stale windows.
                time.sleep(1.0)
                events.append("holder-release")

        def waiter():
            holder_inside.wait(timeout=5)
            with dataset_lock(paths, "SRR1"):
                events.append("waiter-acquire")

        threads = [threading.Thread(target=holder), threading.Thread(target=waiter)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=20)

        assert events == ["holder-release", "waiter-acquire"]

    def test_lock_of_a_dead_holder_is_taken_over_with_a_warning(self, paths, fast_lock, caplog):
        lock = lock_path(paths, "SRR2")
        lock.parent.mkdir(parents=True, exist_ok=True)
        lock.write_text(json.dumps({"pid": 999999, "host": "otherhost", "started": "2020-01-01T00:00:00+00:00"}))
        stale_time = time.time() - 60
        os.utime(lock, (stale_time, stale_time))

        with caplog.at_level(logging.WARNING):
            with dataset_lock(paths, "SRR2"):
                assert json.loads(lock.read_text())["pid"] == os.getpid()

        assert not lock.exists()
        assert any("stale" in record.message.lower() for record in caplog.records)

    def test_release_never_removes_a_foreign_lock(self, paths, fast_lock):
        lock = lock_path(paths, "SRR3")
        with dataset_lock(paths, "SRR3"):
            # Another holder took the lock over while we held it (a stale takeover elsewhere).
            lock.write_text(json.dumps({"pid": 424242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00"}))

        assert lock.exists()
        assert json.loads(lock.read_text())["pid"] == 424242
        lock.unlink()

    def test_lock_wait_gives_up_naming_the_holder(self, paths, fast_lock):
        lock = lock_path(paths, "SRR4")
        lock.parent.mkdir(parents=True, exist_ok=True)
        lock.write_text(json.dumps({"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00"}))
        os.utime(lock, None)

        start = time.monotonic()
        with pytest.raises(DataAccessError) as excinfo:
            with dataset_lock(paths, "SRR4", wait_seconds=0.2):
                pass
        waited = time.monotonic() - start

        message = str(excinfo.value)
        assert "SRR4" in message
        assert "4242" in message
        assert "otherhost" in message
        assert waited < 5


class TestLockInspection:
    def test_lock_is_held_only_while_the_lock_is_fresh(self, paths, fast_lock):
        assert lock_is_held(paths, "SRR5") is False
        with dataset_lock(paths, "SRR5"):
            assert lock_is_held(paths, "SRR5") is True
            assert "otherhost" not in lock_holder(paths, "SRR5")
        assert lock_is_held(paths, "SRR5") is False

    def test_a_stale_lock_does_not_count_as_held(self, paths, fast_lock):
        lock = lock_path(paths, "SRR6")
        lock.parent.mkdir(parents=True, exist_ok=True)
        lock.write_text(json.dumps({"pid": 1, "host": "h", "started": "2020-01-01T00:00:00+00:00"}))
        stale_time = time.time() - 60
        os.utime(lock, (stale_time, stale_time))
        assert lock_is_held(paths, "SRR6") is False


class TestDatasetLockStop:
    """A waiter gives up as soon as its stop predicate is true, so Ctrl-C never waits on a lock."""

    def test_wait_on_a_held_lock_ends_once_stop_is_requested(self, paths):
        from prepare.store.locks import LockWaitStopped

        with dataset_lock(paths, "SRR1"):
            started = time.monotonic()
            with pytest.raises(LockWaitStopped, match="SRR1"):
                with dataset_lock(paths, "SRR1", wait_seconds=0, should_stop=lambda: True):
                    pass
            assert time.monotonic() - started < 1.0

    def test_stop_is_a_data_access_error_for_existing_callers(self):
        from prepare.store.locks import LockWaitStopped

        assert issubclass(LockWaitStopped, DataAccessError)

    def test_a_free_lock_is_taken_whatever_the_predicate(self, paths):
        with dataset_lock(paths, "SRR2", should_stop=lambda: True) as lock:
            assert lock.exists()
