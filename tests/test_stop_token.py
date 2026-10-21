"""Per-run stop token: two download runs in one process stop independently.

``download_sra`` carries one ``threading.Event`` per run through the worker pool, the workers,
``download_accession`` and ``SecureSubprocess.run_secure``. Setting one run's token, or calling
``terminate_children(stop=token)``, stops that run's workers and tools only. The module-level
``accession.STOP`` stays a process-wide emergency stop that no run clears.

Real ``sleep`` children stand in for prefetch and fasterq-dump, as in
``tests/test_interrupt_e2e.py``; ``ALLOWED_EXECUTABLES`` is patched on a copy of the set.
"""

import argparse
import shutil
import gzip
import subprocess
import threading
import time
from unittest.mock import Mock, patch

import pytest

import prepare.data.sra.accession as accession_mod
import prepare.data.sra.retry as retry_mod
import prepare.data.sra.store_handoff as store_handoff_mod
from prepare.data.sra.download import download_sra
from prepare.utils.security import SecureSubprocess

needs_sleep = pytest.mark.skipif(shutil.which("sleep") is None, reason="needs a 'sleep' executable on PATH")


@pytest.fixture(autouse=True)
def _clean_state():
    """Start and end every test with the process-wide stop cleared and no tracked children."""
    accession_mod.STOP.clear()
    SecureSubprocess.clear_stopping()
    yield
    SecureSubprocess.terminate_children(grace=1.0)
    accession_mod.STOP.clear()
    SecureSubprocess.clear_stopping()


@pytest.fixture
def allow_sleep(monkeypatch):
    """Let run_secure start ``sleep`` so a real child stands in for a download tool."""
    monkeypatch.setattr(SecureSubprocess, "ALLOWED_EXECUTABLES", SecureSubprocess.ALLOWED_EXECUTABLES | {"sleep"})


def _children_of(token):
    with SecureSubprocess._children_lock:
        return [proc for proc, owner in SecureSubprocess._children.items() if owner is token]


def _wait_for(condition, timeout=5.0):
    deadline = time.monotonic() + timeout
    while not condition() and time.monotonic() < deadline:
        time.sleep(0.02)
    return condition()


def _start_sleep(token, seconds):
    """Run ``sleep`` under ``token`` in a thread; return the thread and a dict holding its outcome."""
    outcome = {}

    def target():
        try:
            outcome["result"] = SecureSubprocess.run_secure("sleep", [str(seconds)], stop=token)
        except subprocess.CalledProcessError as e:
            outcome["error"] = e

    thread = threading.Thread(target=target)
    thread.start()
    return thread, outcome


class TestRunSecureToken:
    """run_secure records each child under its run's token; terminate_children can target one token."""

    def test_run_secure_with_a_set_token_starts_nothing(self):
        token = threading.Event()
        token.set()
        with patch("subprocess.Popen") as popen:
            with pytest.raises(subprocess.CalledProcessError):
                SecureSubprocess.run_secure("datasets", ["--version"], stop=token)
        popen.assert_not_called()

    def test_run_secure_with_a_set_token_and_no_check_returns_a_failed_result(self):
        token = threading.Event()
        token.set()
        with patch("subprocess.Popen") as popen:
            result = SecureSubprocess.run_secure("datasets", ["--version"], stop=token, check=False)
        popen.assert_not_called()
        assert result.returncode != 0

    @needs_sleep
    def test_terminate_children_with_a_token_kills_only_that_runs_child(self, allow_sleep):
        token_a, token_b = threading.Event(), threading.Event()
        thread_a, outcome_a = _start_sleep(token_a, 30)
        thread_b, outcome_b = _start_sleep(token_b, 30)
        try:
            assert _wait_for(lambda: len(_children_of(token_a)) == 1 and len(_children_of(token_b)) == 1)
            child_b = _children_of(token_b)[0]

            assert SecureSubprocess.terminate_children(grace=1.0, stop=token_a) == 1
            thread_a.join(5)
            assert not thread_a.is_alive()
            assert isinstance(outcome_a.get("error"), subprocess.CalledProcessError)
            assert token_a.is_set(), "terminating a run's children marks that run as stopping"

            assert child_b.poll() is None, "the other run's child must keep running"
            assert thread_b.is_alive()
            assert not token_b.is_set()
            assert SecureSubprocess._stopping is False

            assert SecureSubprocess.terminate_children(grace=1.0) == 1
            thread_b.join(5)
            assert not thread_b.is_alive()
            assert child_b.poll() is not None
            assert isinstance(outcome_b.get("error"), subprocess.CalledProcessError)
        finally:
            SecureSubprocess.terminate_children(grace=1.0)
            thread_a.join(5)
            thread_b.join(5)

    def test_a_child_that_starts_after_its_token_is_terminated_is_killed_at_once(self, monkeypatch):
        monkeypatch.setattr(SecureSubprocess, "_children", {})
        token = threading.Event()
        proc = Mock()
        proc.communicate.return_value = ("", "")
        proc.returncode = -9

        def popen(*args, **kwargs):
            # The token is set between run_secure's own check and the child's registration.
            SecureSubprocess.terminate_children(grace=0.0, stop=token)
            return proc

        with patch("subprocess.Popen", side_effect=popen):
            with pytest.raises(subprocess.CalledProcessError):
                SecureSubprocess.run_secure("datasets", ["--version"], stop=token)
        proc.kill.assert_called_once()
        assert not SecureSubprocess._children

    def test_a_token_that_is_not_set_does_not_stop_other_runs_children(self, monkeypatch):
        monkeypatch.setattr(SecureSubprocess, "_children", {})
        SecureSubprocess.terminate_children(grace=0.0, stop=threading.Event())
        proc = Mock()
        proc.communicate.return_value = ("", "")
        proc.returncode = 0
        with patch("subprocess.Popen", return_value=proc):
            SecureSubprocess.run_secure("datasets", ["--version"], stop=threading.Event())
        proc.kill.assert_not_called()


class TestTokenThroughTheDownloadChain:
    """download_accession, _project_download and _store_download honour the run's token and STOP."""

    def test_download_accession_runs_no_tool_once_its_token_is_set(self, tmp_path):
        token = threading.Event()
        token.set()
        with patch("prepare.data.sra.accession.SecureSubprocess.run_secure") as run:
            with patch("prepare.data.sra.accession.shutil.which", return_value="/usr/bin/prefetch"):
                result = accession_mod.download_accession("SRR2517620", tmp_path / "fastq", stop=token)
        assert result == (False, "interrupted")
        run.assert_not_called()

    def test_download_accession_passes_its_token_to_run_secure(self, tmp_path):
        token = threading.Event()
        seen = []

        def fake_run(executable, args, **kwargs):
            seen.append(kwargs.get("stop"))
            token.set()
            raise subprocess.CalledProcessError(-15, [executable], output="", stderr="")

        with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=fake_run):
            with patch("prepare.data.sra.accession.shutil.which", return_value=None):
                result = accession_mod.download_accession("SRR2517620", tmp_path / "fastq", stop=token)
        assert result == (False, "interrupted")
        assert seen == [token]

    def test_process_wide_stop_still_stops_a_run_with_its_own_token(self, tmp_path):
        accession_mod.STOP.set()
        with patch("prepare.data.sra.accession.SecureSubprocess.run_secure") as run:
            with patch("prepare.data.sra.accession.shutil.which", return_value=None):
                result = accession_mod.download_accession("SRR2517620", tmp_path / "fastq", stop=threading.Event())
        assert result == (False, "interrupted")
        run.assert_not_called()

    def test_compression_is_skipped_once_the_token_is_set(self, tmp_path):
        temp = tmp_path / "SRR1_temp"
        temp.mkdir()
        (temp / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        token = threading.Event()
        token.set()
        with patch("prepare.data.sra.fastq.compress_fastq") as compress:
            success, message = accession_mod._handle_download_output(temp, tmp_path / "SRR1", compress=True, stop=token)
        assert (success, message) == (False, "interrupted")
        compress.assert_not_called()

    def test_compression_runs_pigz_under_the_runs_token(self, tmp_path):
        temp = tmp_path / "SRR1_temp"
        temp.mkdir()
        (temp / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        token = threading.Event()
        seen = []

        def fake_run(executable, args, **kwargs):
            seen.append((executable, kwargs.get("stop")))
            plain = tmp_path / "SRR1" / "SRR1_1.fastq"
            with gzip.open(plain.with_suffix(".fastq.gz"), "wb") as sink:
                sink.write(plain.read_bytes())
            plain.unlink()

        with patch("prepare.data.sra.fastq.shutil.which", return_value="/usr/bin/pigz"):
            with patch("prepare.data.sra.fastq.SecureSubprocess.run_secure", side_effect=fake_run):
                success, _ = accession_mod._handle_download_output(temp, tmp_path / "SRR1", compress=True, stop=token)
        assert success is True
        assert seen == [("pigz", token)]

    def test_the_gzip_fallback_stops_between_blocks_once_the_token_is_set(self, tmp_path):
        from prepare.data.sra.fastq import compress_fastq

        source = tmp_path / "SRR1_1.fastq"
        body = b"@r\nACGT\n+\nIIII\n" * (3 * 1024 * 1024 // 16)  # three 1 MiB blocks
        source.write_bytes(body)

        class _SetAfterFirstBlock(threading.Event):
            checks = 0

            def is_set(self):
                self.checks += 1
                return self.checks > 1

        token = _SetAfterFirstBlock()
        with patch("prepare.data.sra.fastq.shutil.which", return_value=None):
            with pytest.raises(InterruptedError):
                compress_fastq(source, 1, stop=token)
        assert token.checks == 2
        assert source.read_bytes() == body
        assert sorted(p.name for p in tmp_path.iterdir()) == ["SRR1_1.fastq"]

    def test_pigz_stopped_by_the_runs_token_is_reported_as_interrupted(self, tmp_path):
        temp = tmp_path / "SRR1_temp"
        temp.mkdir()
        (temp / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        token = threading.Event()

        def fake_run(executable, args, **kwargs):
            # terminate_children(stop=token) kills pigz as the run stops.
            token.set()
            raise subprocess.CalledProcessError(-15, ["pigz", *args])

        with patch("prepare.data.sra.fastq.shutil.which", return_value="/usr/bin/pigz"):
            with patch("prepare.data.sra.fastq.SecureSubprocess.run_secure", side_effect=fake_run):
                success, message = accession_mod._handle_download_output(
                    temp, tmp_path / "SRR1", compress=True, stop=token
                )
        assert (success, message) == (False, "interrupted")

    def test_project_download_forwards_its_token_to_download_accession(self, tmp_path):
        token = threading.Event()
        seen = {}

        def fake_download(accession, output_folder, *args, **kwargs):
            seen["stop"] = kwargs.get("stop")
            return False, "interrupted"

        with patch.object(accession_mod, "download_accession", side_effect=fake_download):
            accession_mod._project_download("SRR1", tmp_path / "fastq", stop=token)
        assert seen["stop"] is token

    def test_project_download_waiting_on_a_lock_stops_on_the_process_wide_stop(self, tmp_path):
        fastq = tmp_path / "fastq"
        lock = accession_mod.project_lock_path(fastq, "SRR1")
        lock.parent.mkdir(parents=True)
        with patch.object(accession_mod, "PROJECT_LOCK_POLL_SECONDS", 0.05):
            policy = accession_mod._project_lock_policy("SRR1", 0.0)
            from prepare.utils.lockfile import held_lock

            with held_lock(lock, policy):
                accession_mod.STOP.set()
                started = time.monotonic()
                result = accession_mod._project_download("SRR1", fastq, stop=threading.Event())
        assert result == (False, "interrupted")
        assert time.monotonic() - started < 2.0

    def test_store_download_waiting_on_a_held_lock_stops_on_its_token(self, tmp_path):
        from prepare.store.layout import init_store
        from prepare.store.locks import dataset_lock

        store = init_store(tmp_path / "store")
        token = threading.Event()
        token.set()
        with dataset_lock(store, "SRR1"):
            with patch("prepare.data.sra.store_handoff._store_precheck", return_value=None):
                with patch("prepare.data.sra.store_handoff._store_fetch") as fetch:
                    result = store_handoff_mod._store_download(
                        "SRR1", tmp_path / "fastq", store, lock_wait=0.0, stop=token
                    )
        assert result == (False, "interrupted")
        fetch.assert_not_called()

    def test_store_fetch_receives_the_token(self, tmp_path):
        from prepare.store.layout import init_store

        store = init_store(tmp_path / "store")
        token = threading.Event()
        with patch("prepare.data.sra.store_handoff._store_precheck", return_value=None):
            with patch("prepare.data.sra.store_handoff._store_fetch", return_value=(True, "ok")) as fetch:
                store_handoff_mod._store_download("SRR1", tmp_path / "fastq", store, stop=token)
        assert fetch.call_args.kwargs["stop"] is token


class TestRunLevelStop:
    """The worker pool sets the run's token on an interrupt and never clears the process-wide stop."""

    def test_keyboard_interrupt_sets_the_runs_token_not_the_process_wide_stop(self, tmp_path):
        token = threading.Event()

        def worker(acc, *args, **kwargs):
            raise KeyboardInterrupt

        with patch.object(SecureSubprocess, "terminate_children", return_value=0) as terminate:
            with pytest.raises(KeyboardInterrupt):
                retry_mod._execute_parallel_downloads(
                    ["SRR1"], tmp_path, 1, 1, False, None, {}, [], downloader=worker, stop=token
                )
        assert token.is_set()
        assert not accession_mod.STOP.is_set()
        terminate.assert_called_once()
        assert terminate.call_args.kwargs.get("stop") is token

    def test_a_new_run_does_not_clear_the_process_wide_stop(self, tmp_path):
        accession_mod.STOP.set()
        worker = Mock(return_value=(False, "interrupted"))
        retry_mod._execute_parallel_downloads(["SRR1"], tmp_path, 1, 1, False, None, {}, [], downloader=worker)
        assert accession_mod.STOP.is_set()

    def test_workers_receive_the_runs_token(self, tmp_path):
        token = threading.Event()
        worker = Mock(return_value=(True, "ok"))
        retry_mod._execute_parallel_downloads(
            ["SRR1", "SRR2"], tmp_path, 1, 2, False, None, {}, [], downloader=worker, stop=token
        )
        assert [c.kwargs["stop"] for c in worker.call_args_list] == [token, token]

    def test_retry_pass_returns_early_once_the_runs_token_is_set(self, tmp_path):
        token = threading.Event()
        token.set()
        worker = Mock(return_value=(True, "ok"))
        result = retry_mod._retry_failed_downloads(
            ["SRR1"], 2, tmp_path, 1, None, {"SRR1": "Download failed: timeout"}, downloader=worker, stop=token
        )
        assert result == (0, ["SRR1"], None)
        worker.assert_not_called()

    def test_retry_pass_forwards_the_token(self, tmp_path):
        token = threading.Event()
        worker = Mock(return_value=(True, "ok"))
        retry_mod._retry_failed_downloads(
            ["SRR1"], 1, tmp_path, 1, None, {"SRR1": "Download failed: timeout"}, downloader=worker, stop=token
        )
        assert worker.call_args.kwargs["stop"] is token

    def test_a_stop_during_the_retry_pass_leaves_the_remaining_accessions_unrecorded(self, tmp_path):
        token = threading.Event()

        def worker(accession, *args, **kwargs):
            token.set()
            return True, "ok"

        results = {"SRR1": "Download failed: timeout", "SRR2": "Download failed: timeout"}
        seen = []
        retried, failed, abort = retry_mod._retry_failed_downloads(
            ["SRR1", "SRR2"],
            1,
            tmp_path,
            1,
            None,
            results,
            on_result=lambda *args: seen.append(args[0]),
            downloader=worker,
            stop=token,
        )
        assert retried == 1
        assert failed == ["SRR2"]
        assert abort is None
        assert seen == ["SRR1"]
        assert results["SRR2"] == "Download failed: timeout"


@needs_sleep
def test_cancelling_one_download_run_leaves_a_concurrent_run_running(tmp_path, allow_sleep):
    """Two download_sra runs in two threads: stopping A's token ends A's tool and worker only; B finishes."""
    token_a, token_b = threading.Event(), threading.Event()
    a_starting = threading.Event()
    seen = {}

    def fake_project_download(accession, output_folder, num_threads=4, force=False, temp_folder=None, **kwargs):
        stop = kwargs.get("stop")
        seen[accession] = stop
        # Run A's tool would run for 30 s; run B's finishes by itself after 3 s. B starts its
        # tool only once A is starting its own, so B's 3 s cannot run out before both are seen.
        if accession == "SRR1000001":
            a_starting.set()
        else:
            assert a_starting.wait(10), "run A never started its tool"
        seconds = "30" if accession == "SRR1000001" else "3"
        try:
            accession_mod._run_download_tool("sleep", [seconds], stop)
        except accession_mod._DownloadInterrupted:
            return False, "interrupted"
        return True, "Downloaded 1 files, unverified"

    stats = {}

    def run(name, accession, token):
        accessions = tmp_path / f"{name}.txt"
        accessions.write_text(f"{accession}\n")
        stats[name] = download_sra(tmp_path / f"fastq_{name}", accessions, max_retries=1, max_workers=1, stop=token)

    with patch.object(accession_mod, "_project_download", side_effect=fake_project_download):
        thread_a = threading.Thread(target=run, args=("a", "SRR1000001", token_a))
        thread_b = threading.Thread(target=run, args=("b", "SRR1000002", token_b))
        thread_a.start()
        thread_b.start()
        try:
            assert _wait_for(
                lambda: len(_children_of(token_a)) == 1 and len(_children_of(token_b)) == 1
            ), "both runs' children should be running"
            child_b = _children_of(token_b)[0]
            SecureSubprocess.terminate_children(grace=1.0, stop=token_a)
            thread_a.join(5)
            assert not thread_a.is_alive(), "run A did not stop"
            assert child_b.poll() is None, "run B's child was stopped with run A"
            assert not token_b.is_set()
            thread_b.join(10)
            assert not thread_b.is_alive()
        finally:
            SecureSubprocess.terminate_children(grace=1.0)
            thread_a.join(5)
            thread_b.join(5)

    assert seen == {"SRR1000001": token_a, "SRR1000002": token_b}
    assert stats["a"]["successful"] == 0 and stats["a"]["failed"] == 1
    assert stats["a"]["results"]["SRR1000001"] == "interrupted"
    assert stats["b"]["successful"] == 1 and stats["b"]["failed"] == 0
    assert not accession_mod.STOP.is_set()


def test_download_sra_without_a_token_makes_one_for_the_run(tmp_path):
    seen = []

    def fake_project_download(accession, *args, **kwargs):
        seen.append(kwargs.get("stop"))
        return True, "Downloaded 1 files, unverified"

    accessions = tmp_path / "acc.txt"
    accessions.write_text("SRR1\n")
    with patch.object(accession_mod, "_project_download", side_effect=fake_project_download):
        download_sra(tmp_path / "fastq", accessions, max_retries=0)
    assert len(seen) == 1
    assert isinstance(seen[0], threading.Event)
    assert seen[0] is not accession_mod.STOP


class TestCliWiring:
    """download_sra gets the command's Termination token; BaseCommand.run targets that token."""

    def test_download_sra_command_passes_its_termination_token(self, tmp_path):
        from prepare.cli.commands.sra import DownloadSraCommand

        accessions = tmp_path / "acc.txt"
        accessions.write_text("SRR1\n")
        parser = argparse.ArgumentParser()
        command = DownloadSraCommand()
        command.configure_parser(parser)
        args = parser.parse_args(
            [
                "--accessions-file",
                str(accessions),
                "--fastq-folder",
                str(tmp_path / "fastq"),
                "--registry",
                str(tmp_path / "prepare_registry.json"),
            ]
        )
        seen = {}

        def fake_download_sra(**kwargs):
            seen["stop"] = kwargs.get("stop")
            return {
                "total": 1,
                "already_downloaded": 0,
                "blacklisted": 0,
                "successful": 0,
                "failed": 0,
                "failed_accessions": [],
                "results": {},
                "already_downloaded_accessions": [],
                "blacklisted_accessions": [],
                "skipped_accessions": [],
                "aborted": None,
            }

        def execute(parsed, *rest):
            seen["term"] = parsed._termination
            return original_run(parsed, *rest)

        original_run = command._run
        with patch("prepare.cli.commands.sra.require_tools"):
            with patch("prepare.cli.commands.sra.download_sra", side_effect=fake_download_sra):
                with patch.object(command, "_run", side_effect=execute):
                    command.run(args)
        assert seen["stop"] is seen["term"].stop

    def test_base_command_run_terminates_only_its_own_runs_children(self):
        from prepare.cli.base import BaseCommand

        class _Interrupted(BaseCommand):
            @property
            def name(self):
                return "interrupted"

            @property
            def help(self):
                return "raises KeyboardInterrupt"

            def configure_parser(self, parser):
                pass

            def execute(self, args):
                raise KeyboardInterrupt

        args = argparse.Namespace()
        with patch("prepare.cli.base.SecureSubprocess.terminate_children", return_value=0) as terminate:
            assert _Interrupted().run(args) == 130
        terminate.assert_called_once_with(stop=args._termination.stop)
