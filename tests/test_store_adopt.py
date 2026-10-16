"""Tests for `prepare.store.adopt`: folding project-owned FASTQ folders into the store."""

import gzip
import shutil
from pathlib import Path
from unittest.mock import patch

import pytest

from prepare.core.exceptions import DataAccessError
from prepare.store.adopt import ADOPT_COPY_IGNORE, _folder_bytes, adopt
from prepare.store.catalog import Catalog, catalog_write
from prepare.store.layout import init_store, lock_path, sidecar_path, sra_dir
from prepare.store.locks import dataset_lock
from prepare.store.sidecar import build_sidecar, read_sidecar, write_sidecar


def _write_fastq(path: Path, text: str = "@r\nACGT\n+\nIIII\n") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _write_fastq_gz(path: Path, text: str = "@r\nACGT\n+\nIIII\n", mtime: int = 0) -> None:
    """Write ``text`` gzipped, with a fixed header timestamp by default.

    ``gzip.open`` stamps the current wall-clock second into the gzip header, so two calls
    that straddle a second boundary produce different bytes for identical reads. Pinning
    ``mtime`` keeps a test's two copies byte-identical unless it asks for otherwise.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(path, "wb", mtime=mtime) as handle:
        handle.write(text.encode())


class TestAdoptFresh:
    def test_move_removes_project_copy_writes_sidecar_and_links(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.adopted == ["SRR1"]
        assert report.copied == []
        assert report.deduplicated == []
        assert report.conflicts == []

        link = project_fastq / "SRR1"
        assert link.is_symlink()
        assert (sra_dir(paths, "SRR1")).is_dir()
        sidecar = read_sidecar(sidecar_path(paths, "SRR1"))
        assert sidecar is not None
        assert sidecar.tool == "adopted"
        assert sidecar.compression == "gzip"
        assert (sra_dir(paths, "SRR1") / "SRR1.fastq.gz").is_file()

        with Catalog(paths, create=True) as cat:
            cat.migrate()
            assert cat.get_dataset("SRR1") is not None

    def test_copy_leaves_project_folder_untouched_and_unlinked(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        entry = project_fastq / "SRR1"
        _write_fastq(entry / "SRR1.fastq")

        report = adopt(project_fastq, paths, move=False, dry_run=False)

        assert report.copied == ["SRR1"]
        assert report.adopted == []
        # The project's own folder is untouched: still a real directory, never linked.
        assert entry.is_dir() and not entry.is_symlink()
        assert (entry / "SRR1.fastq").is_file()
        # The store now holds its own copy too.
        assert (sra_dir(paths, "SRR1")).is_dir()
        sidecar = read_sidecar(sidecar_path(paths, "SRR1"))
        assert sidecar is not None

    def test_no_compress_leaves_files_plain(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")

        adopt(project_fastq, paths, move=True, dry_run=False, compress=False)

        assert (sra_dir(paths, "SRR1") / "SRR1.fastq").is_file()
        sidecar = read_sidecar(sidecar_path(paths, "SRR1"))
        assert sidecar.compression == "none"

    def test_transient_folder_is_ignored(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        (project_fastq / "SRR1_temp").mkdir(parents=True)
        (project_fastq / ".sra-cache").mkdir(parents=True)

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.adopted == []
        assert report.skipped == []
        assert (project_fastq / "SRR1_temp").is_dir()

    def test_already_linked_entry_is_skipped(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")
        adopt(project_fastq, paths, move=True, dry_run=False)

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.adopted == []
        assert report.skipped == ["SRR1"]

    def test_metadata_folder_supplies_ncbi_spot_count(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")
        metadata_dir = tmp_path / "project" / "metadata"
        metadata_dir.mkdir(parents=True)
        (metadata_dir / "SRR1_metadata.xml").write_text('<root><RUN total_spots="1" total_bases="4"></RUN></root>')

        adopt(project_fastq, paths, move=True, dry_run=False, metadata_folders=[metadata_dir])

        sidecar = read_sidecar(sidecar_path(paths, "SRR1"))
        assert sidecar.ncbi.get("spots") == 1


class TestAdoptDedupAndConflict:
    def _seed_store_dataset(self, paths, accession, text="@r\nACGT\n+\nIIII\n", mtime=0):
        acc_dir = sra_dir(paths, accession)
        _write_fastq_gz(acc_dir / f"{accession}.fastq.gz", text, mtime=mtime)
        sidecar = build_sidecar(accession, acc_dir, {}, "adopted", "gzip")
        write_sidecar(sidecar_path(paths, accession), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)
        return sidecar

    def test_identical_copy_is_deduplicated(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        self._seed_store_dataset(paths, "SRR1")

        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq_gz(project_fastq / "SRR1" / "SRR1.fastq.gz")

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.deduplicated == ["SRR1"]
        assert report.adopted == []
        assert (project_fastq / "SRR1").is_symlink()

    def test_same_reads_gzipped_at_different_times_are_deduplicated(self, tmp_path):
        """Two gzip copies of the same reads differ in raw bytes: the header carries the
        compression time and the original file name, and a different compressor or level
        changes the body too. Two projects that each compressed the same download are a
        duplicate, not a conflict.
        """
        paths = init_store(tmp_path / "store")
        self._seed_store_dataset(paths, "SRR1", mtime=1_000_000_000)

        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq_gz(project_fastq / "SRR1" / "SRR1.fastq.gz", mtime=1_000_000_001)

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.deduplicated == ["SRR1"]
        assert report.conflicts == []
        assert (project_fastq / "SRR1").is_symlink()

    def test_identical_content_dedups_across_differing_compression(self, tmp_path):
        """A plain project copy of the same reads as a gzipped store copy is a duplicate, not
        a conflict: files are compared by decompressed content, not raw bytes/name."""
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        self._seed_store_dataset(paths, "SRR1", text="@r\nACGT\n+\nIIII\n")

        project_fastq = tmp_path / "project" / "fastq"
        # Same reads, but plain (not gzipped) and named without the .gz suffix.
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq", text="@r\nACGT\n+\nIIII\n")

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.deduplicated == ["SRR1"]
        assert report.conflicts == []
        assert (project_fastq / "SRR1").is_symlink()

    def test_conflicting_copy_keeps_both(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        self._seed_store_dataset(paths, "SRR1", text="@r\nACGT\n+\nIIII\n")

        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq_gz(project_fastq / "SRR1" / "SRR1.fastq.gz", text="@r\nTTTT\n+\nIIII\n")

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.conflicts == ["SRR1"]
        assert report.adopted == []
        assert report.deduplicated == []
        # Both copies survive untouched.
        assert (project_fastq / "SRR1").is_dir() and not (project_fastq / "SRR1").is_symlink()
        assert (sra_dir(paths, "SRR1") / "SRR1.fastq.gz").is_file()

    def test_conflicting_read_count_across_differing_compression_is_not_deduplicated(self, tmp_path):
        """Different content (and so a different read count) across differing compression must
        still be flagged as a conflict, not silently treated as a match."""
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        self._seed_store_dataset(paths, "SRR1", text="@r\nACGT\n+\nIIII\n")

        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq", text="@r\nTTTT\n+\nIIII\n")

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.conflicts == ["SRR1"]
        assert report.deduplicated == []


class TestAdoptDryRun:
    def test_dry_run_writes_nothing(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")

        report = adopt(project_fastq, paths, move=True, dry_run=True)

        assert report.planned == ["SRR1"]
        assert report.adopted == []
        # Nothing moved, nothing created in the store.
        assert (project_fastq / "SRR1" / "SRR1.fastq").is_file()
        assert not sra_dir(paths, "SRR1").exists()
        # Not even the catalogue was created: a dry run opens nothing for writing.
        assert not paths.catalog.exists()

    def test_dry_run_still_reports_conflicts(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz", text="@r\nACGT\n+\nIIII\n")
        sidecar = build_sidecar("SRR1", acc_dir, {}, "adopted", "gzip")
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq_gz(project_fastq / "SRR1" / "SRR1.fastq.gz", text="@r\nTTTT\n+\nIIII\n")

        report = adopt(project_fastq, paths, move=True, dry_run=True)

        assert report.conflicts == ["SRR1"]
        assert report.planned == []
        # Original project copy untouched.
        assert (project_fastq / "SRR1" / "SRR1.fastq.gz").is_file()

    def test_dry_run_does_not_sweep_stale_staging(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        stale = paths.tmp / "SRR1_adopt"
        _write_fastq_gz(stale / "SRR1.fastq.gz")
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")

        report = adopt(project_fastq, paths, move=True, dry_run=True)

        assert report.resumed == []
        assert stale.is_dir()


class TestAdoptRestart:
    def test_interrupted_after_final_move_before_sidecar_move_mode(self, tmp_path):
        """The store's final move succeeded (files sit in sra/<ACC>) but the process died
        before the sidecar was written; the project's folder is still there too, since it is
        only ever removed after the sidecar exists. Rerunning with --move finishes the sidecar,
        removes the project's folder, and links it."""
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz")

        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq_gz(project_fastq / "SRR1" / "SRR1.fastq.gz")

        assert not sidecar_path(paths, "SRR1").exists()

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.adopted == ["SRR1"]
        assert sidecar_path(paths, "SRR1").exists()
        assert (project_fastq / "SRR1").is_symlink()
        with Catalog(paths, create=True) as cat:
            cat.migrate()
            assert cat.get_dataset("SRR1") is not None

    def test_interrupted_after_final_move_before_sidecar_copy_mode(self, tmp_path):
        """Same interruption, but this run asked for --copy: the sidecar is finished, and the
        project's folder is left exactly as it was (not removed, not linked)."""
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz")

        project_fastq = tmp_path / "project" / "fastq"
        entry = project_fastq / "SRR1"
        _write_fastq_gz(entry / "SRR1.fastq.gz")

        report = adopt(project_fastq, paths, move=False, dry_run=False)

        assert report.copied == ["SRR1"]
        assert sidecar_path(paths, "SRR1").exists()
        assert entry.is_dir() and not entry.is_symlink()

    def test_a_foreign_sidecar_less_folder_is_reported_and_left_alone(self, tmp_path):
        """A sidecar-less sra/<ACC> for an accession this project never had belongs to
        someone else's interrupted run, most likely one still in flight: adoption reports it
        and touches neither the folder nor the project."""
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        acc_dir = sra_dir(paths, "SRR9000001")
        _write_fastq_gz(acc_dir / "SRR9000001.fastq.gz")

        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.adopted == ["SRR1"]
        assert report.foreign == ["SRR9000001"]
        assert not (project_fastq / "SRR9000001").exists()
        assert not sidecar_path(paths, "SRR9000001").exists()
        assert (acc_dir / "SRR9000001.fastq.gz").is_file()

    def test_an_accession_locked_by_another_run_is_skipped_as_in_progress(self, tmp_path):
        """The store holds our accession without a sidecar and its lock is held: another run
        is publishing it right now, so this one leaves both copies alone."""
        import json

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz")

        project_fastq = tmp_path / "project" / "fastq"
        entry = project_fastq / "SRR1"
        _write_fastq_gz(entry / "SRR1.fastq.gz")

        lock = lock_path(paths, "SRR1")
        lock.parent.mkdir(parents=True, exist_ok=True)
        lock.write_text(json.dumps({"pid": 999999, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00"}))

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.in_progress == ["SRR1"]
        assert report.adopted == []
        assert not sidecar_path(paths, "SRR1").exists()
        assert entry.is_dir() and not entry.is_symlink()

    def test_adoption_refuses_an_accession_without_room_for_it(self, tmp_path):
        """Adoption peaks at three copies of one accession; without room for the staging copy
        it refuses that accession with a message rather than filling the store's filesystem."""
        import shutil as shutil_module
        from collections import namedtuple

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        entry = project_fastq / "SRR1"
        _write_fastq(entry / "SRR1.fastq")

        usage = namedtuple("usage", "total used free")
        with patch.object(shutil_module, "disk_usage", return_value=usage(1000, 1000, 1)):
            report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.refused == ["SRR1"]
        assert report.adopted == []
        assert not sra_dir(paths, "SRR1").exists()
        assert (entry / "SRR1.fastq").is_file()

    def test_an_adopted_sidecar_records_the_tool_and_the_files_own_age(self, tmp_path):
        """An adopted dataset was not downloaded now, and not by fasterq-dump: its sidecar
        says so, so --older-than counts from the reads' own age rather than from adoption."""
        import os
        import time
        from datetime import datetime, timezone

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        entry = project_fastq / "SRR1"
        _write_fastq_gz(entry / "SRR1.fastq.gz")
        old = time.time() - 90 * 86400
        os.utime(entry / "SRR1.fastq.gz", (old, old))

        adopt(project_fastq, paths, move=True, dry_run=False)

        sidecar = read_sidecar(sidecar_path(paths, "SRR1"))
        assert sidecar.tool == "adopted"
        age_days = (datetime.now(timezone.utc) - datetime.fromisoformat(sidecar.downloaded)).days
        assert age_days >= 89


class TestAdoptStagingCrashSafety:
    def test_compress_failure_leaves_project_folder_intact_and_no_store_folder(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        entry = project_fastq / "SRR1"
        _write_fastq(entry / "SRR1.fastq")

        with patch("prepare.store.adopt.compress_fastq", side_effect=RuntimeError("boom")):
            with pytest.raises(RuntimeError, match="boom"):
                adopt(project_fastq, paths, move=True, dry_run=False)

        # The project's original folder was never touched: staging always copies.
        assert entry.is_dir() and not entry.is_symlink()
        assert (entry / "SRR1.fastq").is_file()
        # The store never got a finished folder for this accession.
        assert not sra_dir(paths, "SRR1").exists()
        assert not sidecar_path(paths, "SRR1").exists()
        # Only the staging leftover exists, under tmp.
        assert (paths.tmp / "SRR1_adopt").is_dir()

    def test_rerun_after_compress_failure_adopts_cleanly_with_no_leftover(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        entry = project_fastq / "SRR1"
        _write_fastq(entry / "SRR1.fastq")

        with patch("prepare.store.adopt.compress_fastq", side_effect=RuntimeError("boom")):
            with pytest.raises(RuntimeError):
                adopt(project_fastq, paths, move=True, dry_run=False)

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.resumed == ["SRR1"]
        assert report.adopted == ["SRR1"]
        assert not (paths.tmp / "SRR1_adopt").exists()
        assert (project_fastq / "SRR1").is_symlink()
        sidecar = read_sidecar(sidecar_path(paths, "SRR1"))
        assert sidecar is not None


class TestAdoptPerAccessionLocking:
    """The sweep of a stale staging folder is scoped to one accession's own lock, not a
    store-wide sweep: two concurrent adopts must never touch each other's staging folder."""

    def test_stale_staging_for_the_adopted_accession_is_removed_and_reported(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        stale = paths.tmp / "SRR1_adopt"
        _write_fastq_gz(stale / "SRR1.fastq.gz")

        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.resumed == ["SRR1"]
        assert report.adopted == ["SRR1"]
        assert not stale.exists()
        assert (project_fastq / "SRR1").is_symlink()

    def test_stale_staging_for_a_different_accession_is_left_untouched(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        other_stale = paths.tmp / "SRR-OTHER_adopt"
        _write_fastq_gz(other_stale / "SRR-OTHER.fastq.gz")

        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.resumed == []
        assert report.adopted == ["SRR1"]
        # A stale staging folder for an accession this run never touches survives untouched.
        assert other_stale.is_dir()
        assert (other_stale / "SRR-OTHER.fastq.gz").is_file()

    def test_lock_file_is_created_and_released_on_success(self, tmp_path):
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.adopted == ["SRR1"]
        assert not lock_path(paths, "SRR1").exists()

    def test_a_lock_held_by_a_live_holder_makes_lock_wait_give_up(self, tmp_path):
        import json

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")

        # A live holder (a fresh lock file, so not stale): with no --lock-wait, adopt would
        # rightly wait for as long as that project keeps working, so this run gives up instead.
        lock = lock_path(paths, "SRR1")
        lock.parent.mkdir(parents=True, exist_ok=True)
        lock.write_text(json.dumps({"pid": 999999, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00"}))

        with pytest.raises(DataAccessError, match="SRR1"):
            adopt(project_fastq, paths, move=True, dry_run=False, lock_wait=0.1)

        # The lock was held by "another process": adopt() must not remove a lock it did not
        # create itself.
        assert lock.exists()
        assert json.loads(lock.read_text())["pid"] == 999999
        # Nothing was staged or moved, since the lock was never acquired.
        assert not sra_dir(paths, "SRR1").exists()
        assert (project_fastq / "SRR1" / "SRR1.fastq").is_file()


class TestAdoptMoveUnderTheLock:
    """``--move`` removes the project's folder and links it under one hold of the dataset lock."""

    @staticmethod
    def _hold_elsewhere(paths, accession="SRR1"):
        import json

        lock = lock_path(paths, accession)
        lock.parent.mkdir(parents=True, exist_ok=True)
        lock.write_text(json.dumps({"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00"}))
        return lock

    def test_a_dedup_whose_lock_is_held_keeps_the_project_folder(self, tmp_path, caplog):
        paths = init_store(tmp_path / "store")
        TestAdoptDedupAndConflict()._seed_store_dataset(paths, "SRR1")
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq_gz(project_fastq / "SRR1" / "SRR1.fastq.gz")
        self._hold_elsewhere(paths)

        report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.failed == ["SRR1"]
        assert report.deduplicated == []
        assert not (project_fastq / "SRR1").is_symlink()
        assert (project_fastq / "SRR1" / "SRR1.fastq.gz").is_file()
        assert "SRR1 is being written by" in caplog.text and "retry later" in caplog.text

    def test_a_fresh_adoption_whose_lock_is_taken_before_the_link_keeps_the_project_folder(self, tmp_path):
        import importlib

        adopt_mod = importlib.import_module("prepare.store.adopt")

        paths = init_store(tmp_path / "store")
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")
        original = adopt_mod._apply_move_or_copy

        def _taken_then_apply(*args, **kwargs):
            # Another run takes the lock between the staging and the move.
            self._hold_elsewhere(paths)
            return original(*args, **kwargs)

        with patch.object(adopt_mod, "_apply_move_or_copy", side_effect=_taken_then_apply):
            report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.failed == ["SRR1"]
        assert report.adopted == []
        assert (project_fastq / "SRR1" / "SRR1.fastq").is_file()
        assert (sra_dir(paths, "SRR1") / "SRR1.fastq.gz").is_file()


class TestAdoptIgnoresAppleDouble:
    def test_staged_copy_has_no_appledouble_and_sidecar_is_complete(self, tmp_path):
        paths = init_store(tmp_path / "store")
        project = tmp_path / "proj" / "fastq"
        acc = project / "SRR1"
        acc.mkdir(parents=True)
        _write_fastq_gz(acc / "SRR1_1.fastq.gz")
        (acc / "._SRR1_1.fastq.gz").write_bytes(b"\x00\x05\x16\x07")
        (project / "._SRR1").write_bytes(b"\x00\x05")
        (project / ".DS_Store").write_bytes(b"\x00")
        report = adopt(project, paths, move=True, dry_run=False, compress=True, metadata_folders=[], lock_wait=0)
        assert report.adopted == ["SRR1"]
        names = sorted(p.name for p in (paths.sra / "SRR1").iterdir())
        assert "._SRR1_1.fastq.gz" not in names
        sc = read_sidecar(sidecar_path(paths, "SRR1"))
        assert sc.state == "complete"
        assert [f["name"] for f in sc.files] == ["SRR1_1.fastq.gz"]

    def test_leftover_temporary_files_are_not_copied_into_the_store(self, tmp_path):
        paths = init_store(tmp_path / "store")
        project = tmp_path / "proj" / "fastq"
        acc = project / "SRR1"
        acc.mkdir(parents=True)
        _write_fastq_gz(acc / "SRR1_1.fastq.gz")
        # An interrupted compression of this version and of a version before 0.6.0.
        (acc / ".SRR1_2.fastq.gz.node1.4242.0badf00d.tmp").write_bytes(b"partial")
        (acc / "SRR1_2.fastq.gz.tmp.4242").write_bytes(b"partial")
        report = adopt(project, paths, move=False, dry_run=False, compress=True, metadata_folders=[], lock_wait=0)
        assert report.copied == ["SRR1"]
        assert sorted(p.name for p in (paths.sra / "SRR1").iterdir()) == ["SRR1.json", "SRR1_1.fastq.gz"]
        assert _folder_bytes(acc) == (acc / "SRR1_1.fastq.gz").stat().st_size


class TestFolderBytes:
    """``_folder_bytes`` feeds ``_has_room_for``'s free-space check, which must estimate the
    same folder the staging copy (``shutil.copytree(..., ignore=ADOPT_COPY_IGNORE)``) actually
    produces: disagreeing in either direction is a bug. Counting bytes the copy skips (e.g.
    the contents of a directory whose own name matches an ignore pattern, like ``._cache/``)
    over-counts and can make adoption refuse a folder that actually fits; skipping bytes the
    copy does include (e.g. a plain file inside an ordinary hidden directory like
    ``.snakemake/``, which the copy does *not* skip) under-counts and can let adoption start
    staging a folder that does not actually fit. Earlier code hit both directions at once: the
    original version only checked each file's own name, over-counting a nested ``._cache/``
    directory's contents; a first fix skipped every hidden directory's contents outright,
    under-counting an ordinary hidden directory like ``.hidden/`` whose files the copy does
    include. The current implementation skips a path only when ``_copy_ignores`` matches one
    of its path segments -- the same test ``ADOPT_COPY_IGNORE`` applies during the copy."""

    def test_a_plain_file_inside_an_otherwise_hidden_directory_is_counted(self, tmp_path):
        """A plainly-named file sitting inside a directory whose name matches neither ignore
        pattern (e.g. a `.snakemake/` work folder some pipelines leave under an accession's
        FASTQ directory, or any other ordinary dotted directory) is not skipped by the
        staging copy, so ``_folder_bytes`` must count it too: both the visible top-level
        file's 100 bytes and the 50 bytes inside `.snakemake/` should count."""
        folder = tmp_path / "SRR1"
        folder.mkdir()
        (folder / "visible.fastq").write_bytes(b"x" * 100)
        hidden_dir = folder / ".snakemake"
        hidden_dir.mkdir()
        (hidden_dir / "metadata.txt").write_bytes(b"y" * 50)

        assert _folder_bytes(folder) == 150

    def test_matches_the_bytes_the_staging_copy_actually_produces(self, tmp_path):
        """Pin ``_folder_bytes`` against the real staging copy rather than against a
        restatement of its own logic: build a folder with a plain file, a ``._plain``
        AppleDouble file, a ``.DS_Store``, a ``.hidden/keep.txt`` file, a ``sub/._junk``
        file and a ``._cache/x.bin`` ignored directory; run the same
        ``shutil.copytree(..., ignore=ADOPT_COPY_IGNORE)`` staging uses; and assert
        ``_folder_bytes`` on the source equals the real byte total of the staged copy.
        ``._plain``, ``.DS_Store``, ``sub/._junk`` and everything under ``._cache`` must not
        reach the copy (or the count); ``.hidden/keep.txt`` must reach both, since
        ``.hidden`` matches neither ignore pattern."""
        src = tmp_path / "SRR1"
        src.mkdir()
        (src / "plain.fastq").write_bytes(b"x" * 100)
        (src / "._plain.fastq").write_bytes(b"y" * 20)
        (src / ".DS_Store").write_bytes(b"z" * 5)
        (src / ".hidden").mkdir()
        (src / ".hidden" / "keep.txt").write_bytes(b"a" * 30)
        (src / "sub").mkdir()
        (src / "sub" / "._junk").write_bytes(b"b" * 15)
        (src / "sub" / "real.fastq").write_bytes(b"c" * 7)
        (src / "._cache").mkdir()
        (src / "._cache" / "x.bin").write_bytes(b"d" * 40)

        staged = tmp_path / "staged"
        shutil.copytree(src, staged, ignore=ADOPT_COPY_IGNORE)
        copied_bytes = sum(p.stat().st_size for p in staged.rglob("*") if p.is_file())

        assert copied_bytes == 100 + 30 + 7  # ._plain, .DS_Store, sub/._junk and ._cache excluded
        assert _folder_bytes(src) == copied_bytes


class TestAdoptSafety:
    def test_empty_folder_is_not_adopted(self, tmp_path):
        paths = init_store(tmp_path / "store")
        project = tmp_path / "proj" / "fastq"
        (project / "SRR9").mkdir(parents=True)
        report = adopt(project, paths, move=True, dry_run=False, compress=True, metadata_folders=[], lock_wait=0)
        assert report.adopted == []
        assert report.empty == ["SRR9"]
        assert not (paths.sra / "SRR9").exists()
        assert (project / "SRR9").is_dir() and not (project / "SRR9").is_symlink()

    def test_hidden_files_only_folder_reaches_empty(self, tmp_path):
        """An accession directory whose own name is not hidden, but which holds only hidden
        files (a stray .DS_Store, an AppleDouble sidecar), has no usable FASTQ file: it must
        be reported in ``empty`` exactly like a folder with nothing in it at all, not adopted
        as if the hidden files counted."""
        paths = init_store(tmp_path / "store")
        project = tmp_path / "proj" / "fastq"
        acc = project / "SRR9"
        acc.mkdir(parents=True)
        (acc / ".DS_Store").write_bytes(b"\x00\x01\x02")
        (acc / "._SRR9_1.fastq.gz").write_bytes(b"\x00\x05\x16\x07")
        report = adopt(project, paths, move=True, dry_run=False, compress=True, metadata_folders=[], lock_wait=0)
        assert report.adopted == []
        assert report.empty == ["SRR9"]
        assert not (paths.sra / "SRR9").exists()
        assert (acc / ".DS_Store").is_file()  # left exactly as it was

    def test_failed_staged_copy_keeps_project_folder(self, tmp_path, monkeypatch):
        paths = init_store(tmp_path / "store")
        project = tmp_path / "proj" / "fastq"
        acc = project / "SRR1"
        acc.mkdir(parents=True)
        _write_fastq_gz(acc / "SRR1_1.fastq.gz")
        # prepare.store's __init__ re-exports the `adopt` function under the same name as
        # this submodule, so `import prepare.store.adopt as x` (attribute lookup on the
        # package) would resolve to the function, not the module; importlib sidesteps that.
        import importlib

        adopt_mod = importlib.import_module("prepare.store.adopt")
        real_build = adopt_mod.build_sidecar

        def failing_build(*a, **kw):
            sc = real_build(*a, **kw)
            sc.state = "failed"
            sc.error = "simulated gzip error"
            return sc

        monkeypatch.setattr(adopt_mod, "build_sidecar", failing_build)
        report = adopt(project, paths, move=True, dry_run=False, compress=True, metadata_folders=[], lock_wait=0)
        assert report.failed == ["SRR1"]
        assert (project / "SRR1" / "SRR1_1.fastq.gz").exists()
        assert not (project / "SRR1").is_symlink()

    def test_retry_after_failed_stage_does_not_delete_project_copy(self, tmp_path, monkeypatch):
        """A second adopt() run must not dedup a project copy onto a store copy that a
        previous run already found unverifiable: _dedup_or_conflict's content match alone
        cannot tell a byte-identical bad copy from a good one, so it must also check the
        existing sidecar's state before treating anything as a dedup."""
        paths = init_store(tmp_path / "store")
        project = tmp_path / "proj" / "fastq"
        acc = project / "SRR1"
        acc.mkdir(parents=True)
        _write_fastq_gz(acc / "SRR1_1.fastq.gz")
        import importlib

        adopt_mod = importlib.import_module("prepare.store.adopt")
        real_build = adopt_mod.build_sidecar

        def failing_build(*a, **kw):
            sc = real_build(*a, **kw)
            sc.state = "failed"
            sc.error = "simulated gzip error"
            return sc

        monkeypatch.setattr(adopt_mod, "build_sidecar", failing_build)
        first = adopt(project, paths, move=True, dry_run=False, compress=True, metadata_folders=[], lock_wait=0)
        assert first.failed == ["SRR1"]

        monkeypatch.undo()
        second = adopt(project, paths, move=True, dry_run=False, compress=True, metadata_folders=[], lock_wait=0)
        assert second.failed == ["SRR1"]
        assert second.deduplicated == []
        assert (project / "SRR1" / "SRR1_1.fastq.gz").exists()
        assert not (project / "SRR1").is_symlink()

    def test_copy_dedup_leaves_project_folder(self, tmp_path):
        paths = init_store(tmp_path / "store")
        project = tmp_path / "proj" / "fastq"
        acc = project / "SRR1"
        acc.mkdir(parents=True)
        _write_fastq_gz(acc / "SRR1_1.fastq.gz")
        adopt(project, paths, move=False, dry_run=False, compress=True, metadata_folders=[], lock_wait=0)
        assert (project / "SRR1").is_dir() and not (project / "SRR1").is_symlink()
        report = adopt(project, paths, move=False, dry_run=False, compress=True, metadata_folders=[], lock_wait=0)
        assert report.deduplicated == ["SRR1"]
        assert not (project / "SRR1").is_symlink()
        assert (project / "SRR1" / "SRR1_1.fastq.gz").exists()


class TestContentComparisonPasses:
    """The decompressed comparison reads each side once and takes the read count from that pass."""

    def test_content_matches_reads_each_side_once(self, tmp_path, monkeypatch):
        import builtins
        from collections import Counter

        from prepare.store.adopt import _content_matches

        text = "@r\nACGT\n+\nIIII\n" * 50
        project = tmp_path / "project" / "SRR1.fastq.gz"
        store = tmp_path / "store" / "SRR1.fastq.gz"
        _write_fastq_gz(project, text=text, mtime=1)
        _write_fastq_gz(store, text=text, mtime=2)

        opened = Counter()
        real_open = builtins.open

        def spy(file, *args, **kwargs):
            opened[str(file)] += 1
            return real_open(file, *args, **kwargs)

        monkeypatch.setattr(builtins, "open", spy)
        assert _content_matches(project, store, {"reads": 50}) is True
        assert _content_matches(project, store, {"reads": 49}) is False
        monkeypatch.undo()

        assert opened[str(project)] == 2
        assert opened[str(store)] == 2


class TestAdoptLockLost:
    def test_a_lock_taken_over_during_staging_publishes_nothing(self, tmp_path):
        import json

        from prepare.data.sra import compress_fastq as real_compress

        paths = init_store(tmp_path / "store")
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")
        other = {"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00", "token": "feedbeef"}

        def _compress(path, threads):
            # Another run reclaimed the lock while this one was compressing.
            lock_path(paths, "SRR1").write_text(json.dumps(other))
            return real_compress(path, threads)

        with patch("prepare.store.adopt.compress_fastq", side_effect=_compress):
            report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.failed == ["SRR1"]
        assert report.adopted == []
        assert not sra_dir(paths, "SRR1").exists()
        # The project's folder is untouched and the lock still names the other holder.
        assert (project_fastq / "SRR1").is_dir() and not (project_fastq / "SRR1").is_symlink()
        assert (project_fastq / "SRR1" / "SRR1.fastq").is_file()
        assert json.loads(lock_path(paths, "SRR1").read_text())["token"] == "feedbeef"


class TestAdoptSidecarBeforePublish:
    """The sidecar is written into the staged folder, so the published folder always carries one."""

    def test_the_staged_folder_carries_its_sidecar_when_it_is_published(self, tmp_path):
        paths = init_store(tmp_path / "store")
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")
        seen = {}
        real_move = shutil.move

        def _spy_move(src, dst, *args, **kwargs):
            seen["sidecar_present"] = (Path(src) / "SRR1.json").is_file()
            return real_move(src, dst, *args, **kwargs)

        with patch("prepare.store.adopt.shutil.move", side_effect=_spy_move):
            report = adopt(project_fastq, paths, move=True, dry_run=False)

        assert report.adopted == ["SRR1"]
        assert seen == {"sidecar_present": True}
        sidecar = read_sidecar(sidecar_path(paths, "SRR1"))
        assert sidecar is not None and sidecar.tool == "adopted"
        assert [entry["name"] for entry in sidecar.files] == ["SRR1.fastq.gz"]

    def test_a_kill_right_after_publishing_leaves_a_sidecar_and_the_project_copy(self, tmp_path):
        paths = init_store(tmp_path / "store")
        project_fastq = tmp_path / "project" / "fastq"
        _write_fastq(project_fastq / "SRR1" / "SRR1.fastq")

        with patch("prepare.store.adopt.catalog_write", side_effect=KeyboardInterrupt):
            with pytest.raises(KeyboardInterrupt):
                adopt(project_fastq, paths, move=True, dry_run=False)

        assert sidecar_path(paths, "SRR1").is_file()
        assert read_sidecar(sidecar_path(paths, "SRR1")).state == "complete"
        # The project's own folder is still there, untouched and unlinked.
        assert (project_fastq / "SRR1" / "SRR1.fastq").is_file()
        assert not (project_fastq / "SRR1").is_symlink()


class TestRebuildMissingSidecar:
    """``rebuild_missing_sidecar`` describes a sidecar-less store folder without vouching for it."""

    def test_rebuilds_from_the_files_on_disk_without_cataloguing(self, tmp_path):
        from prepare.store.adopt import rebuild_missing_sidecar

        paths = init_store(tmp_path / "store")
        store_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(store_dir / "SRR1_1.fastq.gz")
        _write_fastq_gz(store_dir / "SRR1_2.fastq.gz")

        with dataset_lock(paths, "SRR1"):
            sidecar = rebuild_missing_sidecar("SRR1", store_dir, paths)

        assert sidecar_path(paths, "SRR1").is_file()
        assert read_sidecar(sidecar_path(paths, "SRR1")) == sidecar
        assert sidecar.accession == "SRR1"
        # Never claimed complete: nothing here read the files through.
        assert sidecar.state == "failed"
        assert "store_verify" in (sidecar.error or "")
        assert sidecar.completeness.get("verdict") == "unverified"
        assert sidecar.layout == "PAIRED"
        assert sidecar.compression == "gzip"
        assert sorted((f["name"], f["bytes"]) for f in sidecar.files) == [
            ("SRR1_1.fastq.gz", (store_dir / "SRR1_1.fastq.gz").stat().st_size),
            ("SRR1_2.fastq.gz", (store_dir / "SRR1_2.fastq.gz").stat().st_size),
        ]
        assert all(f["md5"] is None and f["reads"] is None for f in sidecar.files)
        with Catalog(paths, create=True) as cat:
            cat.migrate()
            assert cat.get_dataset("SRR1") is None

    def test_a_sidecar_written_meanwhile_is_kept(self, tmp_path):
        from prepare.store.adopt import rebuild_missing_sidecar

        paths = init_store(tmp_path / "store")
        store_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(store_dir / "SRR1.fastq.gz")
        existing = build_sidecar("SRR1", store_dir, {}, "3.0", "gzip")
        write_sidecar(sidecar_path(paths, "SRR1"), existing)

        with dataset_lock(paths, "SRR1"):
            sidecar = rebuild_missing_sidecar("SRR1", store_dir, paths)

        assert sidecar == existing
        assert read_sidecar(sidecar_path(paths, "SRR1")) == existing

    def test_an_unreadable_sidecar_is_never_overwritten(self, tmp_path):
        from prepare.store.adopt import rebuild_missing_sidecar

        paths = init_store(tmp_path / "store")
        store_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(store_dir / "SRR1.fastq.gz")
        sidecar_path(paths, "SRR1").write_text("{ not json")

        with dataset_lock(paths, "SRR1"):
            with pytest.raises(DataAccessError):
                rebuild_missing_sidecar("SRR1", store_dir, paths)

        assert sidecar_path(paths, "SRR1").read_text() == "{ not json"
