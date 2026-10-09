"""End-to-end interrupt test.

A real Ctrl-C (SIGINT) during ``_execute_parallel_downloads`` must: re-raise
``KeyboardInterrupt`` promptly, release every dataset lock the interrupted workers held, and
leave no running child process behind. The unit tests around ``SecureSubprocess`` and
``dataset_lock`` cover each piece in isolation (mocked subprocesses, a fake ``should_stop``);
this test instead runs two workers concurrently, each holding a real ``dataset_lock`` while
running a real ``sleep 30`` child through ``SecureSubprocess.run_secure`` (standing in for
``prefetch``/``fasterq-dump``, which ``run_secure`` runs the same way), then sends this
process a real ``SIGINT`` once both children are alive.

``ALLOWED_EXECUTABLES`` is patched (on a copy of the set, matching the pattern in
``tests/test_security_comprehensive.py::TestChildProcessTracking``) to allow ``sleep`` so a
real, harmless, long-running child stands in for a real download tool without depending on
one being installed. Skipped when no ``sleep`` executable is on ``PATH``.
"""

import os
import shutil
import signal
import threading
import time

import pytest

import prepare.data.sra.accession as accession_mod
import prepare.data.sra.retry as retry_mod
from prepare.store.layout import init_store
from prepare.store.locks import dataset_lock
from prepare.utils.security import SecureSubprocess

pytestmark = pytest.mark.skipif(shutil.which("sleep") is None, reason="needs a 'sleep' executable on PATH")

ACCESSIONS = ["SRR1", "SRR2"]


@pytest.fixture(autouse=True)
def _clean_security_state():
    """Every test in this module starts and ends with STOP cleared and no tracked children,
    so a failure here cannot leave a later, unrelated test believing an interrupt is already
    under way or tracking a stale child process."""
    accession_mod.STOP.clear()
    SecureSubprocess.clear_stopping()
    yield
    SecureSubprocess.terminate_children(grace=1.0)
    accession_mod.STOP.clear()
    SecureSubprocess.clear_stopping()


def test_sigint_during_parallel_downloads_stops_children_and_releases_locks(tmp_path, monkeypatch):
    monkeypatch.setattr(SecureSubprocess, "ALLOWED_EXECUTABLES", SecureSubprocess.ALLOWED_EXECUTABLES | {"sleep"})

    paths = init_store(tmp_path / "store")
    fastq_dir = tmp_path / "fastq"

    both_running = threading.Event()
    running_count = {"n": 0}
    count_lock = threading.Lock()
    captured_children: list = []
    # Set in `finally` so the signaler thread never fires a real SIGINT once this test is
    # done, even if it failed (or was interrupted) before reaching the point where it would
    # normally join that thread: a daemon thread that still fires os.kill(SIGINT) 5-10s later
    # hits whatever test happens to be running at that moment, not this one.
    cancel = threading.Event()

    def worker(acc, *args, **kwargs):
        # As the real workers do: wait and run tools under the run's stop token, which the
        # interrupt sets and whose children it terminates.
        stop = kwargs["stop"]
        with dataset_lock(paths, acc, should_stop=lambda: accession_mod.stop_requested(stop)):
            with count_lock:
                running_count["n"] += 1
                if running_count["n"] == len(ACCESSIONS):
                    both_running.set()
            try:
                accession_mod._run_download_tool("sleep", ["30"], stop)
            except accession_mod._DownloadInterrupted:
                return False, "interrupted"
            return True, "ok"

    def send_sigint_once_both_children_are_alive():
        both_running.wait(timeout=5)
        deadline = time.monotonic() + 5
        while len(SecureSubprocess._children) < len(ACCESSIONS) and time.monotonic() < deadline:
            if cancel.is_set():
                return
            time.sleep(0.05)
        with SecureSubprocess._children_lock:
            captured_children.extend(SecureSubprocess._children)
        if cancel.is_set():
            return
        os.kill(os.getpid(), signal.SIGINT)

    signaler = threading.Thread(target=send_sigint_once_both_children_are_alive, daemon=True)
    signaler.start()

    try:
        t0 = time.monotonic()
        with pytest.raises(KeyboardInterrupt):
            retry_mod._execute_parallel_downloads(
                ACCESSIONS,
                fastq_dir,
                1,
                2,
                False,
                None,
                {},
                [],
                downloader=worker,
            )
        elapsed = time.monotonic() - t0

        assert elapsed < 5.0, f"KeyboardInterrupt took {elapsed:.1f}s to re-raise"
        assert len(captured_children) == len(ACCESSIONS), "expected one 'sleep' child per worker"

        # No dataset lock survives the interrupt: each worker's `dataset_lock` block exited
        # (via _DownloadInterrupted, caught inside `worker`) and released its lock.
        assert sorted(p.name for p in paths.locks.iterdir()) == []

        # No child process is still tracked as running.
        assert not SecureSubprocess._children

        # The real 'sleep' processes themselves are gone, not merely untracked: each
        # captured Popen handle reports that its OS process has actually exited.
        for proc in captured_children:
            assert proc.poll() is not None, "a 'sleep' child is still alive after the interrupt"
    finally:
        # Cleanup that must run even if an assertion above failed partway through (including
        # `pytest.raises` itself failing, e.g. no KeyboardInterrupt ever arrived), so this test
        # never leaves a 30s 'sleep' child, a lock file, or a live signaler thread behind for
        # later tests. `cancel` is set, and the signaler joined, before `terminate_children`:
        # once the signaler has stopped, no belated os.kill(SIGINT) can reach a later test.
        cancel.set()
        signaler.join(timeout=6)
        SecureSubprocess.terminate_children(grace=1.0)
        for lock_file in paths.locks.glob("*.lock"):
            lock_file.unlink(missing_ok=True)
