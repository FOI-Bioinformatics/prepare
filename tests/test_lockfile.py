"""Tests for the one lock-file mechanism every PREPARE lock is built on."""

import json
import logging
import multiprocessing
import os
import socket
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

from prepare.core.exceptions import DataAccessError, LockTimeoutError
from prepare.utils import lockfile
from prepare.utils.lockfile import (
    LockHeld,
    LockLost,
    LockPolicy,
    LockReentry,
    LockWaitStopped,
    describe_holder,
    held_lock,
    holder_is_dead,
    read_holder,
    verify_held,
)


def _policy(**overrides):
    """A policy with short intervals, so contention tests run in about a second."""
    values = {
        "what": "Test lock",
        "stale_seconds": 600.0,
        "wait_seconds": 5.0,
        "poll_seconds": 0.01,
        "heartbeat_seconds": 0.05,
    }
    values.update(overrides)
    return LockPolicy(**values)


def _finished_pid():
    """The pid of a subprocess that has already exited, so no process holds it right now."""
    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    proc.wait()
    return proc.pid


def _write_holder(lock, holder, age=0.0):
    lock.write_text(json.dumps(holder))
    stamp = time.time() - age
    os.utime(lock, (stamp, stamp))


def _local_holder(pid):
    holder = {"pid": pid, "host": socket.gethostname(), "started": "2026-01-01T00:00:00+00:00", "token": "t0"}
    pidns = lockfile._pid_namespace()
    if pidns is not None:
        holder["pidns"] = pidns
    return holder


class TestExceptions:
    def test_every_lock_error_is_a_data_access_error(self):
        for cls in (LockHeld, LockLost, LockReentry, LockWaitStopped):
            assert issubclass(cls, DataAccessError)

    def test_store_locks_re_exports_the_same_stop_error(self):
        from prepare.store import locks

        assert locks.LockWaitStopped is LockWaitStopped
        assert locks.read_holder is read_holder


class TestHolderRecord:
    def test_holder_carries_pid_host_start_and_token(self, tmp_path):
        lock = tmp_path / "a.lock"
        with held_lock(lock, _policy()):
            holder = read_holder(lock)
        assert holder["pid"] == os.getpid()
        assert holder["host"] == socket.gethostname()
        assert holder["started"]
        assert len(holder["token"]) == 16
        assert ("pidns" in holder) == (lockfile._pid_namespace() is not None)
        assert not lock.exists()

    def test_two_acquisitions_get_different_tokens(self, tmp_path):
        lock = tmp_path / "a.lock"
        with held_lock(lock, _policy()):
            first = read_holder(lock)["token"]
        with held_lock(lock, _policy()):
            second = read_holder(lock)["token"]
        assert first != second

    def test_bare_pid_reads_as_a_pid_only_record(self, tmp_path):
        lock = tmp_path / "a.lock"
        lock.write_text("12345")
        assert read_holder(lock) == {"pid": 12345}
        assert "12345" in describe_holder(read_holder(lock))

    def test_unreadable_holder_is_an_empty_record(self, tmp_path):
        lock = tmp_path / "a.lock"
        lock.write_text("not json")
        assert read_holder(lock) == {}
        assert read_holder(tmp_path / "missing.lock") == {}
        assert "unknown" in describe_holder({})

    def test_describe_names_pid_host_and_start(self):
        text = describe_holder({"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00"})
        assert text == "pid 4242 on host otherhost since 2026-01-01T00:00:00+00:00"


class TestHolderIsDead:
    def test_finished_local_process_is_dead(self):
        assert holder_is_dead(_local_holder(_finished_pid())) is True

    def test_running_local_process_is_alive(self):
        assert holder_is_dead(_local_holder(os.getpid())) is False

    def test_other_host_is_never_judged_dead(self):
        holder = _local_holder(_finished_pid())
        holder["host"] = "elsewhere.example"
        assert holder_is_dead(holder) is False

    def test_other_pid_namespace_is_never_judged_dead(self):
        holder = _local_holder(_finished_pid())
        holder["pidns"] = -1
        assert holder_is_dead(holder) is False

    def test_permission_error_means_alive(self, monkeypatch):
        def refuse(pid, sig):
            raise PermissionError("not yours")

        monkeypatch.setattr(lockfile.os, "kill", refuse)
        assert holder_is_dead(_local_holder(1)) is False

    def test_bare_pid_is_never_judged_dead(self):
        assert holder_is_dead({"pid": _finished_pid()}) is False
        assert holder_is_dead({}) is False


class TestReclaim:
    def test_reclaim_race_with_forced_stale_judgement_has_one_winner(self, tmp_path, monkeypatch):
        """Every contender judged the same old lock stale; exactly one removes it, and the new
        holder's lock survives the others' late reclaims."""
        monkeypatch.setattr(lockfile, "_judged_stale", lambda age, policy: True)
        policy = _policy()
        for round_number in range(10):
            lock = tmp_path / f"race{round_number}.lock"
            _write_holder(lock, {"pid": 1, "host": "elsewhere.example", "token": "old"}, age=3600)
            observed = lockfile._observe(lock)
            assert observed is not None
            barrier = threading.Barrier(12)
            winners = []

            def contender(index):
                barrier.wait()
                if lockfile._reclaim(lock, observed, policy):
                    holder = lockfile._create(lock)
                    winners.append((index, holder))

            threads = [threading.Thread(target=contender, args=(i,)) for i in range(12)]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join(timeout=10)

            assert len(winners) == 1
            assert winners[0][1] is not None
            assert read_holder(lock)["token"] == winners[0][1]["token"]
            assert not lock.with_name(lock.name + ".reclaim").exists()

    def test_dead_same_host_holder_is_taken_over_at_once(self, tmp_path, caplog):
        lock = tmp_path / "a.lock"
        _write_holder(lock, _local_holder(_finished_pid()))
        started = time.monotonic()
        with caplog.at_level(logging.WARNING):
            with held_lock(lock, _policy(wait_seconds=0.0)):
                assert read_holder(lock)["pid"] == os.getpid()
        assert time.monotonic() - started < 1.0
        assert "no longer running" in caplog.text

    def test_foreign_host_holder_is_respected_until_stale(self, tmp_path):
        lock = tmp_path / "a.lock"
        holder = _local_holder(_finished_pid())
        holder["host"] = "elsewhere.example"
        _write_holder(lock, holder)
        with pytest.raises(LockHeld, match="elsewhere.example"):
            with held_lock(lock, _policy(stale_seconds=0.5, wait_seconds=0.2)):
                pass
        started = time.monotonic()
        with held_lock(lock, _policy(stale_seconds=0.5, wait_seconds=5.0)):
            assert read_holder(lock)["pid"] == os.getpid()
        assert time.monotonic() - started > 0.1

    def test_bare_pid_lock_follows_the_age_rule(self, tmp_path):
        lock = tmp_path / "a.lock"
        lock.write_text(str(_finished_pid()))
        with pytest.raises(LockHeld):
            with held_lock(lock, _policy(), blocking=False):
                pass
        os.utime(lock, (0, 0))
        with held_lock(lock, _policy(), blocking=False):
            assert read_holder(lock)["pid"] == os.getpid()

    def test_old_reclaim_guard_is_removed(self, tmp_path):
        lock = tmp_path / "a.lock"
        guard = tmp_path / "a.lock.reclaim"
        _write_holder(lock, {"pid": 1, "host": "elsewhere.example"}, age=3600)
        guard.write_text("")
        os.utime(guard, (0, 0))
        with held_lock(lock, _policy(stale_seconds=60.0), blocking=False):
            assert read_holder(lock)["pid"] == os.getpid()
        assert not guard.exists()

    def test_fresh_reclaim_guard_defers_the_reclaim(self, tmp_path):
        lock = tmp_path / "a.lock"
        guard = tmp_path / "a.lock.reclaim"
        _write_holder(lock, {"pid": 1, "host": "elsewhere.example"}, age=3600)
        guard.write_text("")
        with pytest.raises(LockHeld):
            with held_lock(lock, _policy(stale_seconds=60.0), blocking=False):
                pass
        assert lock.exists() and guard.exists()


class TestHoldingAndRelease:
    def test_release_never_removes_a_foreign_lock(self, tmp_path):
        lock = tmp_path / "a.lock"
        with held_lock(lock, _policy()):
            lock.write_text(json.dumps({"pid": 424242, "host": "otherhost", "token": "theirs"}))
            with pytest.raises(LockLost, match="424242"):
                verify_held(lock)
        assert read_holder(lock)["token"] == "theirs"

    def test_heartbeat_keeps_a_holder_past_the_stale_window(self, tmp_path):
        lock = tmp_path / "a.lock"
        # A stale window twenty heartbeats wide, so a stall of the shared heartbeat thread on a
        # loaded runner cannot let the waiter in early.
        policy = _policy(stale_seconds=1.0, heartbeat_seconds=0.05, poll_seconds=0.02, wait_seconds=0.0)
        events = []
        inside = threading.Event()

        def holder():
            with held_lock(lock, policy):
                inside.set()
                time.sleep(2.5)
                verify_held(lock)
                events.append("holder-release")

        def waiter():
            inside.wait(timeout=5)
            with held_lock(lock, policy):
                events.append("waiter-acquire")

        threads = [threading.Thread(target=holder), threading.Thread(target=waiter)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=20)
        assert events == ["holder-release", "waiter-acquire"]

    def test_heartbeat_marks_a_removed_lock_lost(self, tmp_path):
        lock = tmp_path / "a.lock"
        with held_lock(lock, _policy(heartbeat_seconds=0.02)):
            verify_held(lock)
            lock.unlink()
            time.sleep(0.2)
            with pytest.raises(LockLost):
                verify_held(lock)

    def test_verify_held_rejects_a_lock_this_process_does_not_hold(self, tmp_path):
        with pytest.raises(LockLost):
            verify_held(tmp_path / "a.lock")

    def test_same_thread_reentry_raises_at_once(self, tmp_path):
        lock = tmp_path / "a.lock"
        with held_lock(lock, _policy(wait_seconds=30.0)):
            started = time.monotonic()
            with pytest.raises(LockReentry):
                with held_lock(lock, _policy(wait_seconds=30.0)):
                    pass
            assert time.monotonic() - started < 0.5
            verify_held(lock)
        assert not lock.exists()

    def test_reentry_with_a_stop_request_reports_the_stop(self, tmp_path):
        lock = tmp_path / "a.lock"
        with held_lock(lock, _policy()):
            with pytest.raises(LockWaitStopped, match="Test lock"):
                with held_lock(lock, _policy(), should_stop=lambda: True):
                    pass

    def test_non_blocking_raises_lock_held(self, tmp_path):
        lock = tmp_path / "a.lock"
        _write_holder(lock, {"pid": 4242, "host": "otherhost", "started": "then", "token": "x"})
        started = time.monotonic()
        with pytest.raises(LockHeld, match="Test lock is locked by pid 4242 on host otherhost since then"):
            with held_lock(lock, _policy(wait_seconds=30.0), blocking=False):
                pass
        assert time.monotonic() - started < 0.5
        assert read_holder(lock)["token"] == "x"

    def test_stop_request_ends_a_wait(self, tmp_path):
        lock = tmp_path / "a.lock"
        _write_holder(lock, {"pid": 4242, "host": "otherhost", "token": "x"})
        with pytest.raises(LockWaitStopped):
            with held_lock(lock, _policy(wait_seconds=0.0), should_stop=lambda: True):
                pass

    def test_lock_released_when_the_block_raises(self, tmp_path):
        lock = tmp_path / "a.lock"
        with pytest.raises(RuntimeError):
            with held_lock(lock, _policy()):
                raise RuntimeError("boom")
        assert not lock.exists()

    @pytest.mark.skipif(not hasattr(os, "fork"), reason="needs os.fork")
    @pytest.mark.filterwarnings("ignore::DeprecationWarning")
    def test_forked_child_neither_heartbeats_nor_releases_the_parent_lock(self, tmp_path):
        lock = tmp_path / "a.lock"
        context = held_lock(lock, _policy())
        context.__enter__()
        try:
            pid = os.fork()
            if pid == 0:  # pragma: no cover - runs in the child
                code = 0
                try:
                    if lockfile._HEARTBEAT.held_count() != 0:
                        code = 2
                    context.__exit__(None, None, None)
                    if not lock.exists():
                        code = 3
                finally:
                    os._exit(code)
            _, status = os.waitpid(pid, 0)
            assert os.waitstatus_to_exitcode(status) == 0
            assert lock.exists()
            verify_held(lock)
        finally:
            context.__exit__(None, None, None)
        assert not lock.exists()


class TestRoundOneFixes:
    """Same-process takeover, reclaim by rename, and a heartbeat that survives failures."""

    def test_same_process_takeover_is_reported_to_the_first_holder(self, tmp_path, monkeypatch):
        """A sibling thread takes over a lock whose heartbeat stalled: the first holder's
        verify_held raises, and the sibling's lock keeps its heartbeat after the first leaves."""
        lock = tmp_path / "a.lock"
        real_refresh = lockfile._refresh
        monkeypatch.setattr(lockfile, "_refresh", lambda entry: None)
        policy = _policy(stale_seconds=0.2, heartbeat_seconds=0.02, poll_seconds=0.01)
        state = {}
        thief_inside = threading.Event()
        holder_left = threading.Event()

        def thief():
            with held_lock(lock, policy):
                monkeypatch.setattr(lockfile, "_refresh", real_refresh)
                thief_inside.set()
                holder_left.wait(timeout=5)
                before = lock.stat().st_mtime
                time.sleep(0.3)
                state["refreshed"] = lock.stat().st_mtime > before
                state["entries"] = len(lockfile._HEARTBEAT.entries_for(os.path.realpath(lock)))
                verify_held(lock)
                state["thief_verify"] = "ok"

        with held_lock(lock, policy):
            time.sleep(0.4)
            thread = threading.Thread(target=thief)
            thread.start()
            assert thief_inside.wait(timeout=5)
            with pytest.raises(LockLost):
                verify_held(lock)
        holder_left.set()
        thread.join(timeout=10)

        assert state == {"refreshed": True, "entries": 1, "thief_verify": "ok"}
        assert not lock.exists()

    def test_contenders_admitted_together_still_have_one_winner(self, tmp_path, monkeypatch):
        """Even when the reclaim guard admits every contender (two waiters removed a crashed
        guard at once), the move-aside-and-compare step lets exactly one reclaim."""
        monkeypatch.setattr(lockfile, "_judged_stale", lambda age, policy: True)
        monkeypatch.setattr(lockfile, "_take_guard", lambda guard: True)
        policy = _policy()
        for round_number in range(20):
            lock = tmp_path / f"race{round_number}.lock"
            _write_holder(lock, {"pid": 1, "host": "elsewhere.example", "token": "old"}, age=3600)
            observed = lockfile._observe(lock)
            barrier = threading.Barrier(12)
            winners = []

            def contender(index):
                barrier.wait()
                if lockfile._reclaim(lock, observed, policy):
                    winners.append(lockfile._create(lock))

            threads = [threading.Thread(target=contender, args=(i,)) for i in range(12)]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join(timeout=10)

            assert len(winners) == 1 and winners[0] is not None
            assert read_holder(lock)["token"] == winners[0]["token"]
            assert sorted(p.name for p in tmp_path.glob(f"race{round_number}.lock*")) == [lock.name]

    def test_a_replaced_lock_moved_aside_is_put_back(self, tmp_path, monkeypatch):
        lock = tmp_path / "a.lock"
        _write_holder(lock, {"pid": 1, "host": "elsewhere.example", "token": "old"}, age=3600)
        observed = lockfile._observe(lock)
        real_observe = lockfile._observe
        calls = []

        def observe_then_replace(path):
            calls.append(path)
            if len(calls) == 1:
                result = real_observe(path)
                # Another waiter replaces the lock between the check and the move.
                lock.unlink()
                _write_holder(lock, {"pid": 2, "host": "elsewhere.example", "token": "new"})
                return result
            return real_observe(path)

        monkeypatch.setattr(lockfile, "_observe", observe_then_replace)
        assert lockfile._reclaim(lock, observed, _policy()) is False
        assert read_holder(lock)["token"] == "new"
        assert sorted(p.name for p in tmp_path.iterdir()) == ["a.lock"]

    def test_heartbeat_survives_an_unexpected_refresh_error(self, tmp_path, monkeypatch):
        broken, healthy = tmp_path / "broken.lock", tmp_path / "healthy.lock"
        real_refresh = lockfile._refresh

        def refresh(entry):
            if entry.path == broken:
                raise RuntimeError("unexpected")
            real_refresh(entry)

        monkeypatch.setattr(lockfile, "_refresh", refresh)
        policy = _policy(heartbeat_seconds=0.02)
        with held_lock(broken, policy), held_lock(healthy, policy):
            before = healthy.stat().st_mtime
            time.sleep(0.3)
            assert lockfile._HEARTBEAT.thread is not None and lockfile._HEARTBEAT.thread.is_alive()
            assert healthy.stat().st_mtime > before
            verify_held(healthy)
            with pytest.raises(LockLost):
                verify_held(broken)

    def test_verify_held_restarts_a_dead_heartbeat_thread(self, tmp_path):
        lock = tmp_path / "a.lock"
        with held_lock(lock, _policy(heartbeat_seconds=0.02)):
            dead = threading.Thread(target=lambda: None)
            dead.start()
            dead.join()
            lockfile._HEARTBEAT.thread = dead
            verify_held(lock)
            assert lockfile._HEARTBEAT.thread is not dead and lockfile._HEARTBEAT.thread.is_alive()

    def test_a_persistent_refresh_failure_is_logged_once(self, tmp_path, monkeypatch, caplog):
        lock = tmp_path / "a.lock"

        def refuse(*args, **kwargs):
            raise PermissionError("read-only file system")

        with caplog.at_level(logging.WARNING, logger="prepare.utils.lockfile"):
            with held_lock(lock, _policy(heartbeat_seconds=0.02)):
                monkeypatch.setattr(lockfile.os, "utime", refuse)
                time.sleep(0.3)
                monkeypatch.undo()
        assert caplog.text.count("Cannot refresh the lock") == 1

    def test_lock_is_held_is_false_when_the_lock_cannot_be_read(self, tmp_path):
        from unittest.mock import patch

        from prepare.store.layout import init_store
        from prepare.store.locks import lock_is_held

        paths = init_store(tmp_path / "store")
        with patch.object(Path, "stat", side_effect=PermissionError("denied")):
            assert lock_is_held(paths, "SRR1") is False

    def test_reentry_through_a_symlinked_folder_is_recognised(self, tmp_path):
        real = tmp_path / "real"
        real.mkdir()
        link = tmp_path / "link"
        link.symlink_to(real)
        with held_lock(real / "a.lock", _policy()):
            with pytest.raises(LockReentry):
                with held_lock(link / "a.lock", _policy(wait_seconds=0.0)):
                    pass

    def test_a_failed_registration_leaves_no_lock_behind(self, tmp_path, monkeypatch):
        lock = tmp_path / "a.lock"

        def refuse(entry):
            raise RuntimeError("cannot start a thread")

        monkeypatch.setattr(lockfile._HEARTBEAT, "add", refuse)
        with pytest.raises(RuntimeError):
            with held_lock(lock, _policy()):
                pass
        assert not lock.exists()
        monkeypatch.undo()
        with held_lock(lock, _policy(), blocking=False):
            assert lock.exists()


class TestStress:
    def test_fifty_threads_never_overlap(self, tmp_path):
        lock = tmp_path / "counter.lock"
        counter = tmp_path / "counter.txt"
        counter.write_text("0")
        markers = []
        policy = _policy(poll_seconds=0.001, wait_seconds=60.0)

        def worker(index):
            for _ in range(10):
                with held_lock(lock, policy):
                    markers.append(("enter", index))
                    value = int(counter.read_text())
                    counter.write_text(str(value + 1))
                    markers.append(("exit", index))

        threads = [threading.Thread(target=worker, args=(i,)) for i in range(50)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=60)

        assert int(counter.read_text()) == 500
        _assert_no_overlap(markers)
        assert not lock.exists()

    @pytest.mark.multiprocess
    def test_eight_processes_never_overlap(self, tmp_path):
        lock = tmp_path / "counter.lock"
        counter = tmp_path / "counter.txt"
        markers = tmp_path / "markers.txt"
        counter.write_text("0")
        context = multiprocessing.get_context("spawn")
        processes = [
            context.Process(target=_process_worker, args=(str(lock), str(counter), str(markers), 25)) for _ in range(8)
        ]
        started = time.monotonic()
        for process in processes:
            process.start()
        for process in processes:
            process.join(timeout=60)
        assert all(process.exitcode == 0 for process in processes)
        assert int(counter.read_text()) == 200
        pairs = [tuple(line.split()) for line in markers.read_text().splitlines()]
        _assert_no_overlap(pairs)
        assert not lock.exists()
        assert time.monotonic() - started < 30


def _assert_no_overlap(markers):
    assert len(markers) % 2 == 0
    for position in range(0, len(markers), 2):
        enter, leave = markers[position], markers[position + 1]
        assert enter[0] == "enter" and leave[0] == "exit" and enter[1] == leave[1], (position, enter, leave)


def _process_worker(lock, counter, markers, iterations):
    """Module-level so the spawn start method can pickle it."""
    policy = LockPolicy(
        "Stress lock", stale_seconds=600.0, wait_seconds=60.0, poll_seconds=0.002, heartbeat_seconds=1.0
    )
    counter_path = Path(counter)
    for _ in range(iterations):
        with held_lock(Path(lock), policy):
            _append(markers, f"enter {os.getpid()}\n")
            value = int(counter_path.read_text())
            counter_path.write_text(str(value + 1))
            _append(markers, f"exit {os.getpid()}\n")


def _append(path, line):
    handle = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
    try:
        os.write(handle, line.encode())
    finally:
        os.close(handle)


class TestCallers:
    def test_nested_registry_transaction_raises_reentry_at_once(self, tmp_path):
        from prepare.data import registry as reg

        target = tmp_path / "prepare_registry.json"
        started = time.monotonic()
        with reg.registry_transaction(target):
            with pytest.raises(LockReentry):
                with reg.registry_transaction(target):
                    pass
        assert time.monotonic() - started < 1.0
        assert not (tmp_path / "prepare_registry.json.lock").exists()

    def test_catalogue_lock_error_names_the_holder(self, tmp_path, monkeypatch):
        from prepare.store import catalog as catalog_module
        from prepare.store.layout import init_store

        paths = init_store(tmp_path / "store")
        monkeypatch.setattr(catalog_module, "CATALOG_LOCK_WAIT_SECONDS", 0.2)
        _write_holder(paths.catalog_lock, {"pid": 4242, "host": "otherhost", "started": "then", "token": "x"})
        expected = f"Store catalogue is locked by pid 4242 on host otherhost since then: {paths.catalog_lock}"
        with pytest.raises(LockHeld) as excinfo:
            with catalog_module.catalog_write(paths):
                pass
        assert str(excinfo.value) == expected
        assert isinstance(excinfo.value, LockTimeoutError)  # a wait that gave up: exit code 4

    def test_dataset_lock_of_a_dead_local_holder_is_taken_over(self, tmp_path):
        from prepare.store.layout import init_store, lock_path
        from prepare.store.locks import dataset_lock, lock_is_held

        paths = init_store(tmp_path / "store")
        lock = lock_path(paths, "SRR1")
        lock.parent.mkdir(parents=True, exist_ok=True)
        _write_holder(lock, _local_holder(_finished_pid()))
        assert lock_is_held(paths, "SRR1") is False
        with dataset_lock(paths, "SRR1", wait_seconds=0.5):
            assert lock_is_held(paths, "SRR1") is True


@pytest.mark.skipif(not hasattr(__import__("signal"), "pthread_sigmask"), reason="needs pthread_sigmask")
class TestSignalsDeferredAroundCreation:
    def test_a_signal_during_creation_is_delivered_after_the_lock_is_recorded_and_then_released(
        self, tmp_path, monkeypatch
    ):
        import signal

        lock = tmp_path / "x.lock"
        real_create = lockfile._create
        seen = {}

        def _create_and_signal(path):
            holder = real_create(path)
            # Without the mask this raises here, between creating the file and recording it.
            # Sent to the main thread itself, as the kernel may deliver a process signal to it.
            signal.pthread_kill(threading.main_thread().ident, signal.SIGINT)
            seen["after_kill"] = True
            return holder

        previous = signal.signal(signal.SIGINT, signal.default_int_handler)
        try:
            monkeypatch.setattr(lockfile, "_create", _create_and_signal)
            with pytest.raises(KeyboardInterrupt):
                with held_lock(lock, _policy()):
                    seen["body"] = True
        finally:
            signal.signal(signal.SIGINT, previous)
        # The signal is raised once the lock is recorded (as the block starts, or on the way in),
        # so the lock is released rather than left behind.
        assert seen["after_kill"] is True
        assert not lock.exists()
        assert signal.SIGINT not in signal.pthread_sigmask(signal.SIG_BLOCK, [])

    def test_no_mask_off_the_main_thread(self, tmp_path):
        import signal

        masks = []

        def _worker():
            with held_lock(tmp_path / "y.lock", _policy()):
                masks.append(signal.pthread_sigmask(signal.SIG_BLOCK, []))

        thread = threading.Thread(target=_worker)
        thread.start()
        thread.join()
        assert masks and signal.SIGINT not in masks[0]
