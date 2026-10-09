"""Tests for the index build lock and the per-sample extraction lock (prepare.data.extraction_locks)."""

import threading
import time
from unittest.mock import patch

import pytest

from prepare.data import extraction_locks
from prepare.data.extraction_locks import (
    SAMPLE_LOCKS_DIRNAME,
    LockHeld,
    index_build_lock,
    index_lock_path,
    sample_extraction_lock,
    sample_lock_path,
)
from prepare.data.file_io import visible_files
from prepare.data.read_extraction import build_index, extract_target_reads
from helpers_extraction import _fake_tools


def _make_tree(tmp_path):
    """One sample, one target genome, above the default threshold."""
    table = tmp_path / "parsed_containment.txt"
    table.write_text("\tGCF_1\nSRR1\t0.9\n")
    fastq = tmp_path / "fastq" / "SRR1"
    fastq.mkdir(parents=True)
    (fastq / "SRR1_1.fastq.gz").write_text("x")
    (fastq / "SRR1_2.fastq.gz").write_text("x")
    genome = tmp_path / "GCF_1.fna"
    genome.write_text(">s\nACGT\n")
    return table, genome


class TestLockPaths:
    def test_index_lock_path_appends_lock_suffix(self, tmp_path):
        index_path = tmp_path / ".index" / "g.sr.mmi"
        assert index_lock_path(index_path) == tmp_path / ".index" / "g.sr.mmi.lock"

    def test_sample_lock_path_is_under_a_hidden_locks_folder(self, tmp_path):
        path = sample_lock_path(tmp_path / "targeted", "SRR1", "GCF_1")
        assert path == tmp_path / "targeted" / SAMPLE_LOCKS_DIRNAME / "SRR1.GCF_1.lock"


class TestLocksFolderInvisible:
    def test_locks_folder_hidden_from_directory_listings(self, tmp_path):
        """``visible_files`` backs every folder listing in prepare (status included), so a
        lock file placed under the hidden ``.locks`` folder never shows up as a fake sample."""
        output_root = tmp_path / "targeted"
        lock = sample_lock_path(output_root, "SRR1", "GCF_1")
        lock.parent.mkdir(parents=True)
        lock.write_text("{}")
        (output_root / "SRR1").mkdir()

        listed = {p.name for p in visible_files(output_root, dirs=True)}
        assert listed == {"SRR1"}


class TestIndexBuildLock:
    def test_second_caller_waits_and_then_reuses_the_index(self, monkeypatch, tmp_path):
        """A blocking lock: the second caller waits for the first to finish, not raises."""
        monkeypatch.setattr(extraction_locks, "LOCK_POLL_SECONDS", 0.01)
        index_path = tmp_path / ".index" / "g.sr.mmi"
        index_path.parent.mkdir(parents=True)
        order = []
        entered = threading.Event()

        def holder():
            with index_build_lock(index_path):
                order.append("holder-in")
                entered.set()
                time.sleep(0.1)
                order.append("holder-out")

        thread = threading.Thread(target=holder)
        thread.start()
        assert entered.wait(timeout=5)
        with index_build_lock(index_path):
            order.append("waiter-in")
        thread.join(timeout=5)

        assert order == ["holder-in", "holder-out", "waiter-in"]

    def test_build_index_runs_minimap2_once_for_two_concurrent_callers(self, monkeypatch, tmp_path):
        """Two threads asking ``build_index`` for the same genome/preset build it once; the
        second reuses what the first built rather than running minimap2 a second time."""
        monkeypatch.setattr(extraction_locks, "LOCK_POLL_SECONDS", 0.01)
        genome = tmp_path / "g.fna"
        genome.write_text(">s\nACGT\n")
        index_dir = tmp_path / ".index"
        state = {}
        base_run = _fake_tools(state)
        holding = threading.Event()
        release = threading.Event()

        def run(executable, args, **kwargs):
            if executable == "minimap2" and "-d" in args:
                holding.set()
                release.wait(timeout=5)
            return base_run(executable, args, **kwargs)

        results = {}

        def worker(name):
            results[name] = build_index(genome, "sr", index_dir)

        with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=run):
            first = threading.Thread(target=worker, args=("first",))
            first.start()
            assert holding.wait(timeout=5)

            second = threading.Thread(target=worker, args=("second",))
            second.start()
            time.sleep(0.05)  # let the second thread reach and start waiting on the lock
            release.set()
            first.join(timeout=5)
            second.join(timeout=5)

        assert results["first"] == results["second"]
        assert results["first"].is_file()
        index_calls = [c for c in state["calls"] if c[0] == "minimap2" and "-d" in c[1]]
        assert len(index_calls) == 1


class TestSampleExtractionLock:
    def test_second_caller_gets_lock_held_at_once(self, tmp_path):
        """Non-blocking: a second caller never waits, it is told the lock is held."""
        output_root = tmp_path / "targeted"
        entered = threading.Event()
        release = threading.Event()

        def holder():
            with sample_extraction_lock(output_root, "SRR1", "GCF_1"):
                entered.set()
                release.wait(timeout=5)

        thread = threading.Thread(target=holder)
        thread.start()
        assert entered.wait(timeout=5)
        try:
            with pytest.raises(LockHeld):
                with sample_extraction_lock(output_root, "SRR1", "GCF_1"):
                    pass
        finally:
            release.set()
            thread.join(timeout=5)

    def test_concurrent_extraction_runs_minimap2_once_and_skips_the_loser(self, tmp_path):
        """Two threads extracting the same sample against the same genome: the mapping run
        happens once, and the thread that loses the lock is reported skipped, not failed."""
        table, genome = _make_tree(tmp_path)
        state = {}
        base_run = _fake_tools(state)
        holding = threading.Event()
        release = threading.Event()

        def run(executable, args, **kwargs):
            if executable == "minimap2" and "-o" in args:
                holding.set()
                release.wait(timeout=5)
            return base_run(executable, args, **kwargs)

        results = {}

        def worker(name):
            results[name] = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=tmp_path / "fastq",
                output_folder=tmp_path / "targeted",
                threshold=0.5,
            )

        with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=run):
            first = threading.Thread(target=worker, args=("first",))
            first.start()
            assert holding.wait(timeout=5)

            second = threading.Thread(target=worker, args=("second",))
            second.start()
            second.join(timeout=5)
            assert not second.is_alive()  # the lock is non-blocking, so this returns at once

            release.set()
            first.join(timeout=5)

        outcomes = [results["first"]["SRR1"], results["second"]["SRR1"]]
        winners = [r for r in outcomes if not r.skipped]
        losers = [r for r in outcomes if r.skipped]
        assert len(winners) == 1 and len(losers) == 1
        assert winners[0].files and winners[0].mapped_records > 0
        assert losers[0].files == [] and losers[0].mapped_records == 0

        map_calls = [c for c in state["calls"] if c[0] == "minimap2" and "-o" in c[1]]
        assert len(map_calls) == 1

        # The .locks folder is cleaned up once both runs finish; the lock file itself is
        # released, not left behind.
        assert not list((tmp_path / "targeted" / SAMPLE_LOCKS_DIRNAME).glob("*.lock"))
