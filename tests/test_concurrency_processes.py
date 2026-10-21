"""Concurrency and signal handling of the prepare CLI, checked with real processes.

Every test starts ``python -m prepare.cli.main`` as separate processes against a project (and,
where named, a shared store) under ``tmp_path``, with fake ``fasterq-dump``, ``prefetch`` and
``pigz`` scripts first on ``PATH`` (``tests/helpers_processes.py``). A fake download tool waits
while ``<barrier>/hold-<ACC>`` exists, so a test keeps one download in progress for as long as it
needs and releases it by removing that file; every wait is on a file or a log line, never on a
fixed delay. The scenarios:

1. Two runs download the same accession into one project: one fetch, both succeed.
2. Two projects download the same accession through one store: one fetch, both linked.
3. ``blacklist`` and ``select_datasets`` run to completion while a download run is in progress,
   and the registry keeps all three commands' records.
4. SIGTERM, then a second SIGTERM, during a download: exit 130, the finished accession recorded,
   no tool left running, no lock file left, the unfinished accession not published. 4b: a tool
   that ignores SIGTERM is killed once the grace period ends, and the run still exits 130.
5. A lock holder killed with SIGKILL: the next run takes the lock over at once, naming the pid.
6. ``store_gc --yes`` while a download holds a dataset's lock: reported ``in_use``, kept.

Two related checks are not repeated here: that two download runs in one process do not share a
stop token is ``tests/test_stop_token.py::test_cancelling_one_download_run_leaves_a_concurrent_run_running``,
and the lock stress tests across threads and processes are ``tests/test_lockfile.py::TestStress``
(its eight-process test also carries the ``multiprocess`` marker).
"""

import json
import os
import shutil
import signal
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import pytest

from prepare.data.registry import REGISTRY_FILENAME
from prepare.data.registry_batch import registry_update
from prepare.data.registry_blocks import ProjectBlock, StoreBlock, set_project_block, set_store_block
from prepare.store.catalog import catalog_write
from prepare.store.layout import init_store
from prepare.utils.lockfile import read_holder
from tests.helpers_processes import (
    FAKE_READS,
    alive,
    cli_env,
    fake_pids,
    install_fake_tools,
    is_fake_tool,
    pid_of_started,
    run_cli,
    spawn_cli,
    started_files,
    stderr_of,
    stop_process,
    wait_for,
)

pytestmark = [
    pytest.mark.multiprocess,
    pytest.mark.skipif(sys.platform == "win32", reason="POSIX signals, sessions and shebang scripts"),
]

# Options every download in this module uses: no NCBI spot counts to verify against, and the
# fasterq-dump-only path, so one fake tool call stands for one fetch.
DOWNLOAD_OPTIONS = ["--no-prefetch", "--no-compress", "--no-verify-downloads", "--num-threads", "1"]
# Upper bound for one command that is not held by a barrier.
COMMAND_TIMEOUT = 20.0


class Harness:
    """One test's fake tools, barrier folder, optional store and the children it started."""

    def __init__(self, tmp_path: Path, with_store: bool = False) -> None:
        self.tmp = tmp_path
        self.bin = tmp_path / "bin"
        self.barrier = tmp_path / "barrier"
        install_fake_tools(self.bin, self.barrier)
        self.store: Optional[Path] = tmp_path / "store" if with_store else None
        self.env = cli_env(tmp_path, self.bin, store=self.store)
        self.children: List[subprocess.Popen] = []

    def project(self, name: str, accessions: List[str]) -> Path:
        """A project folder holding ``accessions.txt`` and a registry, joined to the store if any.

        Set up in this process through the library, as ``store_init`` (or ``status --init``)
        would: only the commands a scenario is about run as separate processes.
        """
        folder = self.tmp / name
        folder.mkdir()
        (folder / "accessions.txt").write_text("".join(f"{acc}\n" for acc in accessions))
        registry_file = folder / REGISTRY_FILENAME
        if self.store is None:
            registry_update(registry_file, lambda registry: None)
        else:
            paths = init_store(self.store)
            project = ProjectBlock(
                id=str(uuid.uuid4()),
                name=name,
                path=str(folder.resolve()),
                created=datetime.now(timezone.utc).isoformat(),
            )

            def _join(registry):
                set_project_block(registry, project)
                set_store_block(registry, StoreBlock(root=str(self.store.resolve())))

            registry_update(registry_file, _join)
            with catalog_write(paths) as catalog:
                catalog.upsert_project(project.id, project.name, project.path, str(registry_file))
        assert registry_file.is_file()
        return folder

    @staticmethod
    def check(result: subprocess.CompletedProcess, expected: int = 0) -> subprocess.CompletedProcess:
        """Assert a finished command's exit code, showing its output when it differs."""
        assert result.returncode == expected, (
            f"prepare {' '.join(result.args)} exited {result.returncode}, expected {expected}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
        return result

    def hold(self, accession: str) -> None:
        """Keep every fake download tool for ``accession`` waiting until ``release``."""
        (self.barrier / f"hold-{accession}").write_text("")

    def release(self, accession: str) -> None:
        """Let the fake download tools for ``accession`` finish."""
        (self.barrier / f"hold-{accession}").unlink(missing_ok=True)

    def download(self, project: Path, *extra: str) -> subprocess.Popen:
        """Start ``download_sra`` for the project's ``accessions.txt``."""
        args = ["download_sra", "--accessions-file", "accessions.txt", *DOWNLOAD_OPTIONS, *extra]
        proc = spawn_cli(args, project, self.env)
        self.children.append(proc)
        return proc

    def wait_started(self, accession: str, count: int = 1, proc: Optional[subprocess.Popen] = None) -> None:
        """Wait until ``count`` fake tool calls for ``accession`` have started."""
        wait_for(
            lambda: len(started_files(self.barrier, accession)) >= count,
            what=f"{count} fake tool start(s) for {accession}",
            detail=(lambda: stderr_of(proc)) if proc is not None else None,
        )

    def finish(self, proc: subprocess.Popen, expected: int = 0, timeout: float = COMMAND_TIMEOUT) -> str:
        """Wait for ``proc`` to exit with ``expected``; return its stderr."""
        try:
            code = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            raise AssertionError(f"process {proc.pid} still running after {timeout} s:\n{stderr_of(proc)}")
        text = stderr_of(proc)
        assert code == expected, f"process {proc.pid} exited {code}, expected {expected}:\n{text}"
        return text

    def cleanup(self) -> None:
        """Release every hold and end every child still running, so a failing test leaves nothing behind."""
        for hold in self.barrier.glob("hold-*"):
            hold.unlink(missing_ok=True)
        for proc in self.children:
            stop_process(proc)
        for pid in fake_pids(self.barrier):
            # A pid is killed only while it still names one of this test's fakes: a fake that
            # exited long ago may have had its pid reused by an unrelated process.
            if alive(pid) and is_fake_tool(pid, self.bin):
                try:
                    os.kill(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass


@pytest.fixture
def harness(tmp_path):
    """A harness for a plain project (no shared store)."""
    h = Harness(tmp_path)
    try:
        yield h
    finally:
        h.cleanup()


@pytest.fixture
def store_harness(tmp_path):
    """A harness whose projects share one store under ``tmp_path``."""
    h = Harness(tmp_path, with_store=True)
    try:
        yield h
    finally:
        h.cleanup()


def _datasets(project: Path) -> Dict[str, dict]:
    return json.loads((project / "prepare_registry.json").read_text()).get("datasets", {})


def _download_state(project: Path, accession: str) -> Optional[str]:
    return _datasets(project).get(accession, {}).get("download", {}).get("state")


def _fastq_names(folder: Path) -> List[str]:
    return sorted(path.name for path in folder.iterdir()) if folder.is_dir() else []


def _temp_names(root: Path) -> List[str]:
    """Paths under ``root`` that look like an unfinished write or build folder."""
    return sorted(
        str(path.relative_to(root))
        for path in root.rglob("*")
        if ".tmp" in path.name or path.name.endswith("_temp") or path.name.endswith("_fqtmp")
    )


def _lock_files(*folders: Path) -> List[str]:
    return sorted(str(path) for folder in folders if folder.is_dir() for path in folder.glob("*.lock"))


def _waiting_logged(proc: subprocess.Popen, accession: str) -> bool:
    """Whether ``proc`` has logged that it waits for ``accession``'s lock (project or store lock)."""
    text = stderr_of(proc)
    return f"waiting for accession {accession}:" in text or f"waiting for {accession}:" in text


# --------------------------------------------------------------------------------------- scenarios


def _with_spots(project: Path, accession: str, spots: int) -> None:
    """Record NCBI's spot count for ``accession`` in the project registry, as download_metadata would."""
    from prepare.data.registry_blocks import MetadataBlock, set_metadata_block

    registry_update(
        project / REGISTRY_FILENAME,
        lambda registry: set_metadata_block(registry, accession, MetadataBlock(run_total_spots=spots)),
    )


def _verified_download(harness: Harness, project: Path, *extra: str, **env: str) -> subprocess.Popen:
    """``download_sra`` with verification on and the fasterq-dump-only path, plus fake-tool env settings."""
    args = ["download_sra", "--accessions-file", "accessions.txt", "--no-compress", "--num-threads", "1", *extra]
    proc = spawn_cli(args, project, {**harness.env, **env})
    harness.children.append(proc)
    return proc


def test_truncated_download_fails_in_plain_project(harness):
    """A run short of NCBI's spot count is a failure, recorded with its verdict and listed for retry."""
    project = harness.project("project", ["SRR1"])
    _with_spots(project, "SRR1", 100)

    log = harness.finish(_verified_download(harness, project, "--no-prefetch", "--max-retries", "0"), expected=1)

    assert "truncated" in log
    dataset = _datasets(project)["SRR1"]["download"]
    assert dataset["state"] == "failed"
    assert dataset["complete"]["verdict"] == "truncated"
    assert dataset["complete"]["reads_r1"] == FAKE_READS
    assert (project / "fastq" / "failed_accessions.txt").read_text() == "SRR1\n"
    assert not (project / "fastq" / "SRR1").exists()
    assert _temp_names(project) == []


def test_truncated_download_accepted_with_accept_partial(harness):
    project = harness.project("project", ["SRR1"])
    _with_spots(project, "SRR1", 100)

    harness.finish(_verified_download(harness, project, "--no-prefetch", "--accept-partial"))

    dataset = _datasets(project)["SRR1"]["download"]
    assert dataset["state"] == "downloaded"
    assert dataset["complete"]["verdict"] == "truncated"
    assert _fastq_names(project / "fastq" / "SRR1") == ["SRR1_1.fastq", "SRR1_2.fastq"]
    assert not (project / "fastq" / "failed_accessions.txt").exists()


def test_tool_failure_then_retry_succeeds(harness):
    """The first fasterq-dump call fails; the retry pass fetches the run and both attempts are counted."""
    project = harness.project("project", ["SRR1"])

    log = harness.finish(_verified_download(harness, project, "--no-prefetch", FAKE_FASTERQ_FAIL_ONCE="1"))

    assert "connection reset by peer" in log
    assert len(started_files(harness.barrier, "SRR1")) == 2
    dataset = _datasets(project)["SRR1"]["download"]
    assert dataset["state"] == "downloaded"
    assert dataset["attempts"] == 2
    assert dataset["complete"]["verdict"] == "unverified"
    assert _temp_names(project) == []


def test_prefetch_without_archive_is_a_network_failure(harness):
    """prefetch exiting 0 with nothing written: a network failure (exit 4) before fasterq-dump runs."""
    project = harness.project("project", ["SRR1"])

    proc = _verified_download(harness, project, "--max-retries", "0", FAKE_PREFETCH_NOOP="1")
    log = harness.finish(proc, expected=4)

    assert "prefetch wrote no archive for SRR1" in log
    assert all(text.strip() == "prefetch" for text in _started_tools(harness.barrier, "SRR1"))
    dataset = _datasets(project)["SRR1"]["download"]
    assert dataset["state"] == "failed"
    assert "prefetch wrote no archive" in dataset["message"]
    assert _temp_names(project) == []


def test_unpaired_file_is_recorded_with_the_mates(harness):
    project = harness.project("project", ["SRR1"])

    harness.finish(_verified_download(harness, project, "--no-prefetch", FAKE_FASTERQ_UNPAIRED="1"))

    assert _fastq_names(project / "fastq" / "SRR1") == ["SRR1.fastq", "SRR1_1.fastq", "SRR1_2.fastq"]
    files = sorted(entry["path"] for entry in _datasets(project)["SRR1"]["download"]["files"])
    assert [Path(path).name for path in files] == ["SRR1.fastq", "SRR1_1.fastq", "SRR1_2.fastq"]


def _started_tools(barrier: Path, accession: str) -> List[str]:
    """The tool name each ``<ACC>.<pid>.started`` file records."""
    return [path.read_text() for path in started_files(barrier, accession)]


def test_two_runs_download_one_accession_into_one_project_once(harness):
    """Scenario 1: the second run waits on the project's accession lock and then finds the files."""
    project = harness.project("project", ["SRR1"])
    harness.hold("SRR1")
    first = harness.download(project)
    harness.wait_started("SRR1", proc=first)
    # --progress-every 0 logs the per-accession lines (here "Skipping ...") at INFO.
    second = harness.download(project, "--progress-every", "0")
    wait_for(
        lambda: _waiting_logged(second, "SRR1"),
        what="the second run to wait for SRR1",
        detail=lambda: stderr_of(second),
    )

    harness.release("SRR1")
    harness.finish(first)
    second_log = harness.finish(second)

    assert len(started_files(harness.barrier, "SRR1")) == 1, second_log
    assert "Skipping SRR1, FASTQ files already exist" in second_log
    assert _fastq_names(project / "fastq" / "SRR1") == ["SRR1_1.fastq", "SRR1_2.fastq"]
    for mate in (1, 2):
        lines = (project / "fastq" / "SRR1" / f"SRR1_{mate}.fastq").read_text().splitlines()
        assert len(lines) == 4 * FAKE_READS, second_log
    assert _download_state(project, "SRR1") == "downloaded"
    # One fetch, one attempt: the run that waited and found the files counts none.
    assert _datasets(project)["SRR1"]["download"]["attempts"] == 1, second_log
    assert _temp_names(project) == []
    assert _lock_files(project / "fastq" / ".locks") == []


def test_two_projects_download_one_accession_through_one_store_once(store_harness):
    """Scenario 2: the second project waits on the store's dataset lock and then links the copy."""
    h = store_harness
    first_project = h.project("first", ["SRR1"])
    second_project = h.project("second", ["SRR1"])
    h.hold("SRR1")
    first = h.download(first_project)
    h.wait_started("SRR1", proc=first)
    second = h.download(second_project)
    wait_for(
        lambda: _waiting_logged(second, "SRR1"),
        what="the second project to wait for SRR1",
        detail=lambda: stderr_of(second),
    )

    h.release("SRR1")
    h.finish(first)
    h.finish(second)

    assert len(started_files(h.barrier, "SRR1")) == 1
    store_copy = (h.store / "sra" / "SRR1").resolve()
    assert "SRR1_1.fastq" in _fastq_names(store_copy)
    for project in (first_project, second_project):
        link = project / "fastq" / "SRR1"
        assert link.is_symlink() and link.resolve() == store_copy
        assert _download_state(project, "SRR1") == "downloaded"
        assert _datasets(project)["SRR1"]["download"].get("source") == "store"
    assert _temp_names(h.store) == []
    assert _lock_files(h.store / "locks") == []


def test_registry_keeps_blacklist_and_selection_written_during_a_download(harness, tmp_path):
    """Scenario 3: two registry writers finish while a download run holds one of three accessions."""
    project = harness.project("project", ["SRR1", "SRR2", "SRR3"])
    (project / "parsed_containment.txt").write_text("\tGCF_A\tmax_containment\nSRR7\t0.9\t0.9\nSRR8\t0.05\t0.05\n")
    harness.hold("SRR2")
    download = harness.download(project, "--max-workers", "3")
    harness.wait_started("SRR2", proc=download)
    wait_for(
        lambda: all((project / "fastq" / acc / f"{acc}_1.fastq").is_file() for acc in ("SRR1", "SRR3")),
        what="SRR1 and SRR3 to be published",
        detail=lambda: stderr_of(download),
    )

    harness.check(run_cli(["blacklist", "--add", "SRR9", "--reason", "contaminated"], project, harness.env))
    harness.check(
        run_cli(
            ["select_datasets", "--parsed-containment", "parsed_containment.txt", "--output", "selected.txt"],
            project,
            harness.env,
        )
    )
    assert download.poll() is None, "the download run finished before the other writers ran"

    harness.release("SRR2")
    harness.finish(download)

    datasets = _datasets(project)
    assert datasets["SRR9"]["exclusion"]["reason"] == "contaminated"
    assert datasets["SRR7"]["selection"]["selected"] is True
    assert (project / "selected.txt").read_text().split() == ["SRR7"]
    for accession in ("SRR1", "SRR2", "SRR3"):
        assert datasets[accession]["download"]["state"] == "downloaded", accession


def _scratch(h, project, accession):
    """fasterq-dump's scratch folder for ``accession``: under the store's tmp, or the project's staging folder."""
    if h.store is not None:
        return h.store / "tmp" / f"{accession}_fqtmp"
    return project / "fastq" / ".prepare-tmp" / f"{accession}_fqtmp"


@pytest.mark.parametrize("with_store", [False, True], ids=["project", "store"])
def test_two_sigterms_during_a_download_stop_the_run_cleanly(tmp_path, with_store):
    """Scenario 4: the first SIGTERM stops the run, the second is only logged; exit 130."""
    h = Harness(tmp_path, with_store=with_store)
    try:
        project = h.project("project", ["SRR1", "SRR2"])
        h.hold("SRR2")
        download = h.download(project, "--max-workers", "2")
        h.wait_started("SRR2", proc=download)
        wait_for(
            lambda: (project / "fastq" / "SRR1" / "SRR1_1.fastq").is_file(),
            what="SRR1 to be published",
            detail=lambda: stderr_of(download),
        )
        assert _scratch(h, project, "SRR2").is_dir(), "the held tool's scratch folder should exist while it runs"

        # The held fake tool delays its exit on SIGTERM while linger-SRR2 exists, so the run is
        # still stopping when the second signal arrives, at least 0.2 s after the first.
        (h.barrier / "linger-SRR2").write_text("")
        stopped_at = time.monotonic()
        os.kill(download.pid, signal.SIGTERM)
        wait_for(
            lambda: any(h.barrier.glob("SRR2.*.terminating")),
            what="the run to pass SIGTERM on to its fake tool",
            detail=lambda: stderr_of(download),
        )
        time.sleep(max(0.0, 0.2 - (time.monotonic() - stopped_at)))
        os.kill(download.pid, signal.SIGTERM)
        wait_for(
            lambda: "Received SIGTERM again" in stderr_of(download),
            what="the second SIGTERM to be logged",
            detail=lambda: stderr_of(download),
        )
        (h.barrier / "linger-SRR2").unlink()
        log = h.finish(download, expected=130, timeout=10.0)
        assert time.monotonic() - stopped_at < 10.0

        pids = fake_pids(h.barrier)
        wait_for(
            lambda: not any(alive(pid) for pid in pids), timeout=5.0, what="every fake tool to exit", detail=lambda: log
        )
        assert _download_state(project, "SRR1") == "downloaded", log
        assert _download_state(project, "SRR2") != "downloaded", log
        assert not (project / "fastq" / "SRR2").exists()
        lock_folders = [project / "fastq" / ".locks"]
        if h.store is not None:
            lock_folders.append(h.store / "locks")
            assert not (h.store / "sra" / "SRR2").exists()
        assert _lock_files(*lock_folders) == [], log
        assert not _scratch(h, project, "SRR1").exists(), log
        assert not _scratch(h, project, "SRR2").exists(), log
    finally:
        h.cleanup()


def test_a_tool_that_ignores_sigterm_is_killed_after_the_grace_period(harness):
    """Scenario 4b: SIGTERM to a run whose tool ignores it; the tool is SIGKILLed, the run exits 130."""
    project = harness.project("project", ["SRR1"])
    harness.hold("SRR1")
    (harness.barrier / "ignore-term-SRR1").write_text("")
    download = harness.download(project)
    harness.wait_started("SRR1", proc=download)
    tool_pid = pid_of_started(started_files(harness.barrier, "SRR1")[0])
    assert _scratch(harness, project, "SRR1").is_dir(), "the held tool's scratch folder should exist while it runs"

    os.kill(download.pid, signal.SIGTERM)
    log = harness.finish(download, expected=130, timeout=20.0)
    assert not _scratch(harness, project, "SRR1").exists(), log

    wait_for(lambda: not alive(tool_pid), timeout=5.0, what="the fake tool to be killed", detail=lambda: log)
    assert any(harness.barrier.glob("SRR1.*.ignored")), log
    assert _download_state(project, "SRR1") != "downloaded", log
    assert not (project / "fastq" / "SRR1").exists()
    assert _lock_files(project / "fastq" / ".locks") == [], log


def test_a_killed_lock_holder_is_taken_over_at_once(store_harness):
    """Scenario 5: SIGKILL of the run holding a store lock; the next run takes it over, naming the pid."""
    h = store_harness
    first_project = h.project("first", ["SRR1"])
    second_project = h.project("second", ["SRR1"])
    h.hold("SRR1")
    first = h.download(first_project)
    h.wait_started("SRR1", proc=first)
    lock = h.store / "locks" / "SRR1.lock"
    assert read_holder(lock).get("pid") == first.pid, stderr_of(first)
    tool_pid = pid_of_started(started_files(h.barrier, "SRR1")[0])

    # run_secure starts each tool in a session of its own, so a kill of the run's process group
    # does not reach the tool; a scheduler ending the job (SLURM tracks it by cgroup) kills both.
    os.killpg(first.pid, signal.SIGKILL)
    os.killpg(tool_pid, signal.SIGKILL)
    first_log = h.finish(first, expected=-signal.SIGKILL, timeout=10.0)
    wait_for(
        lambda: not alive(tool_pid), timeout=5.0, what="the killed run's fake tool to exit", detail=lambda: first_log
    )
    assert lock.exists(), "the killed run's lock file should still be on disk"

    h.release("SRR1")
    second = h.download(second_project)
    log = h.finish(second, timeout=10.0)

    assert f"Took over the lock on SRR1: its holder is no longer running (pid {first.pid} " in log
    # At once: on the first look at the lock, before a waiter logs that it is waiting (after 1 s).
    assert not _waiting_logged(second, "SRR1"), log
    assert len(started_files(h.barrier, "SRR1")) == 2
    link = second_project / "fastq" / "SRR1"
    assert link.is_symlink() and "SRR1_1.fastq" in _fastq_names(link)
    assert _download_state(second_project, "SRR1") == "downloaded"
    assert _lock_files(h.store / "locks") == []


def test_store_gc_keeps_a_dataset_whose_download_is_in_progress(store_harness):
    """Scenario 6: a dataset that would be removed is reported in use and kept while its lock is held."""
    h = store_harness
    live = h.project("live", ["SRR1"])
    gone = h.project("gone", ["SRR1"])
    h.finish(h.download(gone))
    # The only project that used SRR1 disappears (folder, registry and symlink), and its last
    # use is two days old, so gc would remove the dataset with --include-stale. A symlink left
    # behind would keep it: a live link on this host protects a dataset whatever the registry.
    shutil.rmtree(gone)
    two_days_ago = time.time() - 2 * 86400
    os.utime(h.store / "locks" / "SRR1.used", (two_days_ago, two_days_ago))
    dry = h.check(run_cli(["store_gc", "--include-stale", "--json"], live, h.env))
    assert [entry["accession"] for entry in json.loads(dry.stdout)["datasets"]] == ["SRR1"], dry.stdout

    h.hold("SRR1")
    refetch = h.download(live, "--force")
    h.wait_started("SRR1", count=2, proc=refetch)
    gc = h.check(run_cli(["store_gc", "--yes", "--include-stale", "--json"], live, h.env))
    report = json.loads(gc.stdout)

    assert [entry["accession"] for entry in report["in_use"]] == ["SRR1"], gc.stdout
    assert report["datasets"] == [] or all(entry["accession"] != "SRR1" for entry in report["datasets"])
    assert "SRR1_1.fastq" in _fastq_names(h.store / "sra" / "SRR1")
    h.release("SRR1")
    h.finish(refetch)
    assert "SRR1_1.fastq" in _fastq_names(h.store / "sra" / "SRR1")
    assert (live / "fastq" / "SRR1").is_symlink()
