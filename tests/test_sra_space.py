"""Tests for the download free-space guard (prepare.data.sra.space) and the first-pass disk-full abort."""

import argparse
import threading
import time
from collections import Counter, namedtuple
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from prepare.cli.commands.sra import DownloadSraCommand
from prepare.data.registry import load_registry, record_metadata, save_registry
from prepare.data.sra import download as download_mod
from prepare.data.sra import retry as retry_mod
from prepare.data.sra import space as space_mod

GB = space_mod.GB
E = space_mod.FASTQ_EXPANSION
# The output folder also holds the gzip-compressed files while they are written.
OUT = space_mod.FASTQ_EXPANSION + space_mod.GZIP_EXPANSION
Usage = namedtuple("Usage", "total used free")


@pytest.fixture
def fake_disks(monkeypatch):
    """Map folder names to fake filesystems: ``disks.devices[path] = dev``, ``disks.free[dev] = bytes``.

    A path's device is looked up by its last component; an unknown name raises OSError, like a
    filesystem that cannot be read. ``disks.free[dev] = None`` makes ``disk_usage`` fail.
    """

    class Disks:
        devices: dict = {}
        free: dict = {}

    disks = Disks()
    disks.devices, disks.free = {}, {}

    def device(path):
        return disks.devices.get(Path(path).name)

    def disk_usage(path):
        free = disks.free[disks.devices[Path(path).name]]
        if free is None:
            raise OSError("cannot stat")
        return Usage(10 * free, 0, free)

    monkeypatch.setattr(space_mod, "_device", device)
    monkeypatch.setattr(space_mod, "_existing", lambda path: Path(path))
    monkeypatch.setattr(space_mod, "mount_point", lambda path: Path(path))
    monkeypatch.setattr(space_mod.shutil, "disk_usage", disk_usage)
    return disks


def _locations(tmp_path):
    return {space_mod.OUTPUT: tmp_path / "out", space_mod.TEMP: tmp_path / "tmp", space_mod.CACHE: tmp_path / "cache"}


class TestNeeds:
    def test_known_size_on_three_filesystems(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 2, "cache": 3})
        guard = space_mod.SpaceGuard(_locations(tmp_path), GB, {"SRR1": 100}, use_prefetch=True)
        assert guard.needs("SRR1") == {1: OUT * 100, 2: E * 100, 3: 100}

    def test_needs_of_one_filesystem_are_added(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 1, "cache": 1})
        guard = space_mod.SpaceGuard(_locations(tmp_path), GB, {"SRR1": "100"}, use_prefetch=True)
        assert guard.needs("SRR1") == {1: OUT * 100 + E * 100 + 100}

    def test_no_cache_without_prefetch(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 2, "cache": 3})
        guard = space_mod.SpaceGuard(_locations(tmp_path), GB, {"SRR1": 100}, use_prefetch=False)
        assert guard.needs("SRR1") == {1: OUT * 100, 2: E * 100}

    def test_unknown_size_needs_the_floor_on_each_filesystem(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 1, "cache": 2})
        guard = space_mod.SpaceGuard(_locations(tmp_path), 5 * GB, {"SRR1": None, "SRR2": "n/a"}, use_prefetch=True)
        assert guard.needs("SRR1") == {1: 5 * GB, 2: 5 * GB}
        assert guard.needs("SRR2") == {1: 5 * GB, 2: 5 * GB}

    def test_exempt_accession_needs_nothing(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 1, "cache": 1})
        guard = space_mod.SpaceGuard(_locations(tmp_path), GB, {}, use_prefetch=True, exempt={"SRR1"})
        assert guard.needs("SRR1") == {}
        assert guard.reserve("SRR1") is None


class TestReserve:
    def test_a_second_reserve_without_a_release_raises(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 2, "cache": 3})
        fake_disks.free.update({1: 100 * GB, 2: 100 * GB, 3: 100 * GB})
        guard = space_mod.SpaceGuard(_locations(tmp_path), GB, {"SRR1": GB // 2}, use_prefetch=True)
        assert guard.reserve("SRR1") is None
        reserved = dict(guard._reserved)
        with pytest.raises(RuntimeError, match="SRR1 is already reserved"):
            guard.reserve("SRR1")
        assert guard._reserved == reserved

    def test_reserve_release_twice_leaves_nothing_reserved(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 2, "cache": 3})
        fake_disks.free.update({1: 100 * GB, 2: 100 * GB, 3: 100 * GB})
        guard = space_mod.SpaceGuard(_locations(tmp_path), GB, {"SRR1": GB // 2}, use_prefetch=True)
        for _ in range(2):
            assert guard.reserve("SRR1") is None
            guard.release("SRR1")
        assert sum(guard._reserved.values()) == 0

    def test_refuses_when_one_filesystem_is_short(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 2, "cache": 3})
        fake_disks.free.update({1: 100 * GB, 2: 3 * GB, 3: 100 * GB})
        guard = space_mod.SpaceGuard(_locations(tmp_path), GB, {"SRR1": GB // 2}, use_prefetch=True)
        message = guard.reserve("SRR1")
        assert message == (
            f"insufficient-space: not enough free space on {tmp_path / 'tmp'}: 3.0 GB free, about 4.0 GB needed"
        )
        # Not a disk-full message: it fails this accession only and does not stop the run.
        assert retry_mod.accession_mod.classify_download_error(message) != "disk-full"

    def test_short_only_because_of_reservations_waits(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 1, "cache": 1})
        fake_disks.free[1] = 10 * GB
        guard = space_mod.SpaceGuard(_locations(tmp_path), 6 * GB, {}, use_prefetch=True)
        assert guard.reserve("SRR1") is None
        # SRR2 would fit once SRR1's reservation is returned, so it waits rather than refusing;
        # a stop request ends the wait without reserving anything.
        assert guard.reserve("SRR2", should_stop=lambda: True) == "interrupted"
        guard.release("SRR1")
        assert guard.reserve("SRR2") is None

    def test_larger_than_free_plus_reservations_is_refused_at_once(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 1, "cache": 1})
        fake_disks.free[1] = 10 * GB
        guard = space_mod.SpaceGuard(_locations(tmp_path), 6 * GB, {"SRR2": 2 * GB}, use_prefetch=True)
        assert guard.reserve("SRR1") is None
        stop = Mock(return_value=False)
        assert guard.reserve("SRR2", should_stop=stop).startswith("insufficient-space:")
        stop.assert_not_called()

    def test_two_workers_reserving_at_once_cannot_both_have_the_space(self, tmp_path, fake_disks, monkeypatch):
        monkeypatch.setattr(space_mod, "WAIT_POLL_SECONDS", 0.05)
        fake_disks.devices.update({"out": 1, "tmp": 1, "cache": 1})
        fake_disks.free[1] = 10 * GB
        guard = space_mod.SpaceGuard(_locations(tmp_path), 6 * GB, {}, use_prefetch=True)
        barrier = threading.Barrier(2)
        results = {}
        first = threading.Event()

        def worker(accession):
            barrier.wait()
            results[accession] = guard.reserve(accession)
            first.set()

        threads = [threading.Thread(target=worker, args=(acc,)) for acc in ("SRR1", "SRR2")]
        for thread in threads:
            thread.start()
        assert first.wait(5)
        time.sleep(0.2)
        assert len(results) == 1  # the other one is waiting, not refused
        guard.release(next(iter(results)))
        for thread in threads:
            thread.join(5)
        assert results == {"SRR1": None, "SRR2": None}

    def test_wait_ends_on_the_stop_token(self, tmp_path, fake_disks, monkeypatch):
        monkeypatch.setattr(space_mod, "WAIT_POLL_SECONDS", 0.05)
        fake_disks.devices.update({"out": 1, "tmp": 1, "cache": 1})
        fake_disks.free[1] = 10 * GB
        guard = space_mod.SpaceGuard(_locations(tmp_path), 6 * GB, {}, use_prefetch=True)
        assert guard.reserve("SRR1") is None
        stop = threading.Event()
        wrapped = retry_mod._instrumented(Mock(return_value=(True, "ok")), guard, None)
        result = {}
        thread = threading.Thread(target=lambda: result.setdefault("r", wrapped("SRR2", stop=stop)))
        thread.start()
        time.sleep(0.2)
        assert thread.is_alive()
        stop.set()
        thread.join(5)
        assert result["r"] == (False, "interrupted")

    def test_unknown_size_below_the_floor_is_refused(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 1, "cache": 1})
        fake_disks.free[1] = 4 * GB
        guard = space_mod.SpaceGuard(_locations(tmp_path), 5 * GB, {}, use_prefetch=True)
        assert "about 5.0 GB needed" in guard.reserve("SRR1")

    def test_known_size_is_checked_against_the_estimate_not_the_floor(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 1, "cache": 1})
        fake_disks.free[1] = 4 * GB
        guard = space_mod.SpaceGuard(_locations(tmp_path), 5 * GB, {"SRR1": 100}, use_prefetch=True)
        assert guard.reserve("SRR1") is None

    def test_guard_off_with_zero(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 1, "cache": 1})
        fake_disks.free[1] = 0
        guard = space_mod.SpaceGuard(_locations(tmp_path), 0, {"SRR1": GB}, use_prefetch=True)
        assert not guard.enabled
        assert guard.reserve("SRR1") is None
        assert guard.preflight(["SRR1"]) == []

    def test_unmeasurable_free_space_proceeds(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 1, "cache": 1})
        fake_disks.free[1] = None
        guard = space_mod.SpaceGuard(_locations(tmp_path), 5 * GB, {}, use_prefetch=True)
        assert guard.reserve("SRR1") is None

    def test_unreadable_filesystem_is_not_checked(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1})  # tmp and cache: device unknown
        fake_disks.free[1] = 100 * GB
        guard = space_mod.SpaceGuard(_locations(tmp_path), 5 * GB, {}, use_prefetch=True)
        assert guard.needs("SRR1") == {1: 5 * GB}
        assert guard.reserve("SRR1") is None


class TestPreflight:
    def test_warns_when_the_known_sizes_do_not_fit(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 2, "cache": 2})
        fake_disks.free.update({1: 100 * GB, 2: GB})
        guard = space_mod.SpaceGuard(_locations(tmp_path), GB, {"SRR1": GB // 8, "SRR2": GB // 8}, use_prefetch=True)
        warnings = guard.preflight(["SRR1", "SRR2", "SRR3"])
        assert len(warnings) == 1
        assert "The 2 downloads with a known size" in warnings[0]
        assert str(tmp_path / "tmp") in warnings[0]

    def test_quiet_when_everything_fits(self, tmp_path, fake_disks):
        fake_disks.devices.update({"out": 1, "tmp": 1, "cache": 1})
        fake_disks.free[1] = 100 * GB
        guard = space_mod.SpaceGuard(_locations(tmp_path), GB, {"SRR1": GB}, use_prefetch=True)
        assert guard.preflight(["SRR1"]) == []


class TestLocations:
    def test_plain_project_defaults(self, tmp_path):
        locations = space_mod.download_locations(tmp_path / "fastq")
        assert locations[space_mod.OUTPUT] == tmp_path / "fastq"
        assert locations[space_mod.CACHE] == tmp_path / "fastq" / ".sra-cache"

    def test_store_defaults_and_explicit_folders(self, tmp_path):
        store = Mock(tmp=tmp_path / "store" / "tmp")
        locations = space_mod.download_locations(tmp_path / "fastq", tmp_path / "scratch", None, store)
        assert locations == {
            space_mod.OUTPUT: tmp_path / "store" / "tmp",
            space_mod.TEMP: tmp_path / "scratch",
            space_mod.CACHE: tmp_path / "store" / "tmp" / ".sra-cache",
        }

    def test_real_filesystem_helpers(self, tmp_path):
        missing = tmp_path / "a" / "b"
        assert space_mod._existing(missing) == tmp_path
        assert space_mod._device(missing) == space_mod._device(tmp_path)
        assert space_mod.mount_point(tmp_path) in [tmp_path, *tmp_path.parents]


class TestInstrumented:
    def test_times_every_call(self):
        timings = {}
        wrapped = retry_mod._instrumented(Mock(return_value=(True, "ok")), None, timings)
        assert wrapped("SRR1", "fastq", 4, force=False) == (True, "ok")
        started, seconds = timings["SRR1"]
        assert started.endswith("+00:00")
        assert seconds >= 0

    def test_refusal_does_not_call_the_worker(self):
        guard = Mock()
        guard.reserve.return_value = "disk-full: insufficient free space on /x: 0.0 GB free, about 1.0 GB needed"
        worker = Mock()
        timings = {}
        result = retry_mod._instrumented(worker, guard, timings)("SRR1")
        assert result == (False, guard.reserve.return_value)
        worker.assert_not_called()
        guard.release.assert_not_called()
        assert timings == {}

    def test_release_even_when_the_worker_raises(self):
        guard = Mock()
        guard.reserve.return_value = None
        wrapped = retry_mod._instrumented(Mock(side_effect=OSError("boom")), guard, None)
        with pytest.raises(OSError):
            wrapped("SRR1")
        guard.release.assert_called_once_with("SRR1")


def _disk_full_worker(gate):
    """SRR1 fails with a full disk; any other accession waits for ``gate``, then succeeds."""
    calls = []

    def worker(accession, *args, **kwargs):
        calls.append(accession)
        if accession == "SRR1":
            return False, "disk-full: Download failed: no space left on device"
        gate.wait(5)
        return True, "Downloaded 2 files"

    return worker, calls


class TestFirstPassDiskFullAbort:
    def test_pending_downloads_are_not_attempted(self, tmp_path):
        accessions = ["SRR1", "SRR2", "SRR3", "SRR4", "SRR5"]
        gate = threading.Event()
        worker, calls = _disk_full_worker(gate)
        notified = []

        def on_result(accession, success, message):
            notified.append((accession, success, message))
            if message == retry_mod.DISK_FULL_NOT_ATTEMPTED:
                gate.set()

        results, failed = {}, []
        successful, failed_count, abort = retry_mod._execute_parallel_downloads(
            accessions, tmp_path, 1, 1, False, None, results, failed, on_result=on_result, downloader=worker
        )

        assert abort == "disk-full"
        # One worker: SRR2 may have started before the abort; SRR3 to SRR5 cannot have.
        assert set(calls) <= {"SRR1", "SRR2"}
        for acc in ("SRR3", "SRR4", "SRR5"):
            assert results[acc] == "disk-full: not attempted"
            assert acc in failed
        assert Counter(acc for acc, _, _ in notified) == Counter(accessions)
        assert successful + failed_count == len(accessions)

    def test_no_retry_pass_after_a_first_pass_abort(self, tmp_path):
        gate = threading.Event()
        gate.set()
        worker, calls = _disk_full_worker(gate)
        successful, failed_count, failed, results, abort = retry_mod._download_with_retries(
            ["SRR1"], tmp_path, 1, 1, False, None, 3, downloader=worker
        )
        assert abort == "disk-full"
        assert calls == ["SRR1"]
        assert failed == ["SRR1"]

    def test_guard_refusal_does_not_abort_the_run(self, tmp_path):
        guard = Mock()
        guard.reserve.side_effect = lambda acc, **kwargs: (
            "insufficient-space: not enough free space on /x" if acc == "SRR1" else None
        )
        worker = Mock(return_value=(True, "ok"))
        timings = {}
        successful, _, failed, results, abort = retry_mod._download_with_retries(
            ["SRR1", "SRR2", "SRR3"], tmp_path, 1, 1, False, None, 1, downloader=worker, guard=guard, timings=timings
        )
        assert abort is None
        assert successful == 2
        assert failed == ["SRR1"]
        assert "SRR1" not in timings
        assert results["SRR2"] == results["SRR3"] == "ok"
        # The retry pass skips a refusal, as it does a not-found accession: no second reservation.
        assert [c.args[0] for c in guard.reserve.call_args_list].count("SRR1") == 1
        assert results["SRR1"].startswith("insufficient-space:")

    def test_retry_pass_uses_the_instrumented_worker(self, tmp_path):
        attempts = Counter()

        def flaky(accession, *args, **kwargs):
            attempts[accession] += 1
            return (attempts[accession] > 1), "network: connection reset"

        guard = Mock()
        guard.reserve.return_value = None
        timings = {}
        with patch.object(retry_mod.time, "sleep"):
            successful, _, failed, _, abort = retry_mod._download_with_retries(
                ["SRR1"], tmp_path, 1, 1, False, None, 1, downloader=flaky, guard=guard, timings=timings
            )
        assert (successful, failed, abort) == (1, [], None)
        assert guard.reserve.call_count == 2
        assert guard.release.call_count == 2
        assert "SRR1" in timings


class TestGuardInARun:
    """The real SpaceGuard in the real download loops, on a fake 100 GB filesystem."""

    def _guard(self, tmp_path, fake_disks, monkeypatch, sizes):
        monkeypatch.setattr(space_mod, "WAIT_POLL_SECONDS", 0.05)
        fake_disks.devices.update({"out": 1})
        fake_disks.free[1] = 100 * GB  # constant: the fake downloads write nothing
        return space_mod.SpaceGuard({space_mod.OUTPUT: tmp_path / "out"}, 10 * GB, sizes, use_prefetch=False)

    def _counting_worker(self):
        state = {"running": 0, "peak": 0, "order": []}
        lock = threading.Lock()

        def worker(accession, *args, **kwargs):
            with lock:
                state["running"] += 1
                state["peak"] = max(state["peak"], state["running"])
            time.sleep(0.1)
            with lock:
                state["running"] -= 1
                state["order"].append(accession)
            return True, "Downloaded 2 files, complete"

        return worker, state

    def test_downloads_that_each_fit_run_one_after_another(self, tmp_path, fake_disks, monkeypatch):
        # Three accessions of about 60 GB each (8 x 7.5 GB) on a 100 GB disk with 2 workers.
        sizes = {acc: int(7.5 * GB) for acc in ("A", "B", "C")}
        guard = self._guard(tmp_path, fake_disks, monkeypatch, sizes)
        worker, state = self._counting_worker()
        successful, failed_count, failed, results, abort = retry_mod._download_with_retries(
            ["A", "B", "C"], tmp_path, 1, 2, False, None, 2, downloader=worker, guard=guard
        )
        assert (successful, failed_count, failed, abort) == (3, 0, [], None)
        assert state["peak"] == 1
        assert sorted(state["order"]) == ["A", "B", "C"]

    def test_an_accession_larger_than_the_disk_fails_alone(self, tmp_path, fake_disks, monkeypatch):
        sizes = {"A": GB, "HUGE": 20 * GB, "C": GB}  # HUGE needs 200 GB on the output folder
        guard = self._guard(tmp_path, fake_disks, monkeypatch, sizes)
        worker, state = self._counting_worker()
        successful, failed_count, failed, results, abort = retry_mod._download_with_retries(
            ["A", "HUGE", "C"], tmp_path, 1, 2, False, None, 1, downloader=worker, guard=guard
        )
        assert abort is None
        assert (successful, failed_count, failed) == (2, 1, ["HUGE"])
        assert "insufficient-space: not enough free space" in results["HUGE"]
        assert "about 200.0 GB needed" in results["HUGE"]
        assert sorted(state["order"]) == ["A", "C"]


class TestDownloadSraGuard:
    def _run(self, tmp_path, **kwargs):
        accessions = tmp_path / "acc.txt"
        accessions.write_text("SRR1\n")
        with patch.object(retry_mod, "_download_with_retries", return_value=(1, 0, [], {"SRR1": "ok"}, None)) as run:
            download_mod.download_sra(tmp_path / "fastq", accessions, **kwargs)
        return run.call_args.kwargs

    def test_zero_turns_the_guard_off(self, tmp_path):
        assert self._run(tmp_path, min_free_gb=0)["guard"] is None

    def test_setting_is_used_when_not_given(self, tmp_path, monkeypatch):
        monkeypatch.delenv("PREPARE_MIN_FREE_GB")
        guard = self._run(tmp_path, run_sizes={"SRR1": 123})["guard"]
        assert guard.floor_bytes == 10 * GB
        assert guard._sizes == {"SRR1": 123}

    def test_timings_are_passed_through(self, tmp_path):
        timings = {}
        assert self._run(tmp_path, min_free_gb=1, timings=timings)["timings"] is timings

    def test_store_ready_accessions_are_exempt(self, tmp_path):
        store = Mock(tmp=tmp_path / "store" / "tmp")
        with patch.object(download_mod.store_handoff_mod, "_store_state", return_value="ready"):
            guard = download_mod._space_guard(
                tmp_path / "fastq", None, None, store, True, {}, 1.0, False, ["SRR1", "SRR2"]
            )
        assert guard.needs("SRR1") == {} and guard.needs("SRR2") == {}

    def test_preflight_warnings_are_logged(self, tmp_path, caplog):
        with patch.object(space_mod.SpaceGuard, "preflight", return_value=["not enough room"]):
            with caplog.at_level("WARNING"):
                self._run(tmp_path, min_free_gb=1)
        assert "not enough room" in caplog.text


class TestCommand:
    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_run_sizes_and_min_free_gb_reach_download_sra(self, mock_download, _which, tmp_path):
        registry_path = tmp_path / "prepare_registry.json"
        seeded = load_registry(registry_path)
        record_metadata(seeded, "SRR1", tmp_path / "SRR1.xml", {"run_size": 5000, "run_total_spots": 10})
        record_metadata(seeded, "SRR2", tmp_path / "SRR2.xml", {"run_total_spots": 10})
        save_registry(seeded)
        mock_download.return_value = {"total": 1, "successful": 1, "failed": 0, "failed_accessions": []}
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\nSRR2\n")
        parser = argparse.ArgumentParser()
        command = DownloadSraCommand()
        command.configure_parser(parser)
        args = parser.parse_args(
            ["--accessions-file", str(acc), "--fastq-folder", str(tmp_path / "fastq"), "--registry", str(registry_path)]
            + ["--min-free-gb", "2.5", "--max-workers", "1"]
        )
        assert command.execute(args) == 0
        kwargs = mock_download.call_args.kwargs
        assert kwargs["run_sizes"] == {"SRR1": 5000}
        assert kwargs["min_free_gb"] == 2.5

    def test_help_names_the_cap(self):
        parser = argparse.ArgumentParser()
        DownloadSraCommand().configure_parser(parser)
        text = " ".join(parser.format_help().split())
        assert "CPUs available to this job" in text
        assert "PREPARE_MAX_WORKERS_CAP" in text
        assert "--min-free-gb" in text


def test_a_refused_retry_drops_the_first_attempt_time(tmp_path):
    """The first pass times a failed attempt; the guard refuses the retry; no time is left for it."""
    guard = Mock()
    calls = {"reserve": 0}

    def reserve(acc, **kwargs):
        calls["reserve"] += 1
        return None if calls["reserve"] == 1 else "interrupted"

    guard.reserve.side_effect = reserve
    worker = Mock(return_value=(False, "network: connection reset"))
    timings = {}
    retry_mod._download_with_retries(
        ["SRR1"], tmp_path, 1, 1, False, None, 1, downloader=worker, guard=guard, timings=timings
    )
    assert worker.call_count == 1
    assert "SRR1" not in timings
