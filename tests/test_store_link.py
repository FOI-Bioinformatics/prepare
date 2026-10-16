"""
Tests for prepare.store.link: the symlink a project keeps into the shared store.
"""

import os
import shutil
from pathlib import Path

import pytest

from prepare.core.exceptions import DataAccessError
from prepare.data.registry import scan_downloads
from prepare.data.sra import accession_has_fastq
from prepare.store.layout import init_store
from prepare.store.link import dangling_links, is_store_link, link_dataset, unlink_dataset


@pytest.fixture
def paths(tmp_path):
    return init_store(tmp_path / "store")


def _store_dataset(paths, accession="SRR1", state="complete"):
    """Create a store dataset folder with one FASTQ file and a sidecar in ``state``."""
    acc_dir = paths.sra / accession
    acc_dir.mkdir(parents=True, exist_ok=True)
    (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
    (acc_dir / f"{accession}.json").write_text('{"accession": "%s", "state": "%s"}' % (accession, state))
    return acc_dir


# ------------------------------------------------------------------ link_dataset


def test_link_dataset_creates_a_relative_symlink_under_a_shared_parent(tmp_path, paths):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"

    link = link_dataset(project_fastq, "SRR1", paths)

    assert link == project_fastq / "SRR1"
    assert link.is_symlink()
    assert not os.path.isabs(os.readlink(link))
    assert link.resolve() == (paths.sra / "SRR1").resolve()


def test_link_dataset_absolute_mode_writes_an_absolute_target(tmp_path, paths):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"

    link = link_dataset(project_fastq, "SRR1", paths, mode="absolute")

    assert os.path.isabs(os.readlink(link))
    assert link.resolve() == (paths.sra / "SRR1").resolve()


def test_link_dataset_relative_mode_writes_a_relative_target(tmp_path, paths):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"

    link = link_dataset(project_fastq, "SRR1", paths, mode="relative")

    assert not os.path.isabs(os.readlink(link))


def test_link_dataset_auto_falls_back_to_absolute_without_a_shared_parent(tmp_path, paths, monkeypatch):
    """Two paths with no common parent above the filesystem root get an absolute link."""
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    monkeypatch.setattr("prepare.store.link.os.path.commonpath", lambda paths_: os.sep)

    link = link_dataset(project_fastq, "SRR1", paths)

    assert os.path.isabs(os.readlink(link))


def test_link_dataset_auto_falls_back_to_absolute_on_separate_drives(tmp_path, paths, monkeypatch):
    """commonpath raises for paths on different drives; the link is absolute rather than failing."""
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"

    def _raise(_paths):
        raise ValueError("paths do not share a drive")

    monkeypatch.setattr("prepare.store.link.os.path.commonpath", _raise)

    link = link_dataset(project_fastq, "SRR1", paths)

    assert os.path.isabs(os.readlink(link))


def test_link_dataset_copy_mode_copies_the_folder(tmp_path, paths):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"

    copied = link_dataset(project_fastq, "SRR1", paths, mode="copy")

    assert copied.is_dir() and not copied.is_symlink()
    assert (copied / "SRR1_1.fastq").read_text() == "@r\nACGT\n+\nIIII\n"
    # The store copy is untouched.
    assert (paths.sra / "SRR1" / "SRR1_1.fastq").exists()


def test_link_dataset_refuses_to_replace_a_real_directory(tmp_path, paths):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    real_dir = project_fastq / "SRR1"
    real_dir.mkdir(parents=True)
    (real_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")

    with pytest.raises(DataAccessError):
        link_dataset(project_fastq, "SRR1", paths)

    assert (real_dir / "SRR1_1.fastq").exists()


def test_link_dataset_replaces_an_existing_symlink(tmp_path, paths):
    _store_dataset(paths)
    _store_dataset(paths, accession="SRR2")
    project_fastq = tmp_path / "project" / "fastq"
    project_fastq.mkdir(parents=True)
    os.symlink(paths.sra / "SRR2", project_fastq / "SRR1")

    link = link_dataset(project_fastq, "SRR1", paths)

    assert link.resolve() == (paths.sra / "SRR1").resolve()


def test_link_dataset_raises_when_the_store_dataset_is_missing(tmp_path, paths):
    with pytest.raises(DataAccessError):
        link_dataset(tmp_path / "project" / "fastq", "SRR404", paths)


def test_link_dataset_rejects_an_unknown_mode(tmp_path, paths):
    _store_dataset(paths)
    with pytest.raises(DataAccessError):
        link_dataset(tmp_path / "project" / "fastq", "SRR1", paths, mode="hardlink")


# ---------------------------------------------------------------- unlink_dataset


def test_unlink_dataset_removes_only_symlinks(tmp_path, paths):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths)

    assert unlink_dataset(project_fastq, "SRR1") is True
    assert not (project_fastq / "SRR1").exists()
    # The store keeps its copy.
    assert (paths.sra / "SRR1" / "SRR1_1.fastq").exists()


def test_unlink_dataset_leaves_a_real_directory_alone(tmp_path, paths):
    project_fastq = tmp_path / "project" / "fastq"
    real_dir = project_fastq / "SRR1"
    real_dir.mkdir(parents=True)
    (real_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")

    assert unlink_dataset(project_fastq, "SRR1") is False
    assert (real_dir / "SRR1_1.fastq").exists()


def test_unlink_dataset_returns_false_when_nothing_is_there(tmp_path):
    assert unlink_dataset(tmp_path / "fastq", "SRR1") is False


def test_unlink_dataset_removes_a_dangling_symlink(tmp_path, paths):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths)
    shutil.rmtree(paths.sra / "SRR1")

    assert unlink_dataset(project_fastq, "SRR1") is True
    assert not (project_fastq / "SRR1").is_symlink()


# ----------------------------------------------------------------- is_store_link


def test_is_store_link_distinguishes_store_links_from_other_entries(tmp_path, paths):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths)

    other = tmp_path / "elsewhere" / "SRR2"
    other.mkdir(parents=True)
    os.symlink(other, project_fastq / "SRR2")

    real_dir = project_fastq / "SRR3"
    real_dir.mkdir()

    assert is_store_link(project_fastq / "SRR1", paths) is True
    assert is_store_link(project_fastq / "SRR2", paths) is False
    assert is_store_link(real_dir, paths) is False
    assert is_store_link(project_fastq / "SRR9", paths) is False


def test_is_store_link_true_for_a_dangling_store_link(tmp_path, paths):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths)
    shutil.rmtree(paths.sra / "SRR1")

    assert is_store_link(project_fastq / "SRR1", paths) is True


# --------------------------------------------------------------- dangling_links


def test_dangling_links_lists_links_whose_target_is_gone(tmp_path, paths):
    _store_dataset(paths)
    _store_dataset(paths, accession="SRR2")
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths)
    link_dataset(project_fastq, "SRR2", paths)
    (project_fastq / "SRR3").mkdir()
    shutil.rmtree(paths.sra / "SRR2")

    assert dangling_links(project_fastq) == ["SRR2"]


def test_dangling_links_empty_for_a_missing_folder(tmp_path):
    assert dangling_links(tmp_path / "nope") == []


# -------------------------------------------------- the rest of the tool sees it


def test_scan_downloads_counts_a_linked_dataset(tmp_path, paths):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths)

    found = scan_downloads(project_fastq)

    assert list(found) == ["SRR1"]
    assert found["SRR1"][0] == 1
    assert accession_has_fastq(project_fastq / "SRR1") is True


def test_scan_downloads_skips_a_link_to_a_partial_dataset(tmp_path, paths):
    _store_dataset(paths, state="partial")
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths)

    assert scan_downloads(project_fastq) == {}
    assert accession_has_fastq(project_fastq / "SRR1") is False


def test_scan_downloads_skips_a_dangling_link(tmp_path, paths):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths)
    shutil.rmtree(paths.sra / "SRR1")

    assert scan_downloads(project_fastq) == {}


def test_linked_dataset_is_readable_through_the_project_path(tmp_path, paths):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths)

    assert (Path(project_fastq) / "SRR1" / "SRR1_1.fastq").read_text().startswith("@r")


def test_link_dataset_relative_link_for_a_relative_project_folder(tmp_path, paths, monkeypatch):
    """A relative --fastq-folder under a shared parent still gets a working relative link."""
    _store_dataset(paths)
    monkeypatch.chdir(tmp_path)

    link = link_dataset(Path("project") / "fastq", "SRR1", paths)

    assert not os.path.isabs(os.readlink(link))
    assert link.resolve() == (paths.sra / "SRR1").resolve()
    assert (link / "SRR1_1.fastq").is_file()


class TestLinkModeAcrossVolumes:
    """Two mounted volumes are not one tree, whatever their paths look like."""

    def test_two_volumes_get_an_absolute_link(self):
        from prepare.store.link import _shares_a_parent

        # /Volumes/A and /Volumes/B mount and unmount independently: a relative link between
        # them (../../A/store/...) breaks the moment either moves.
        assert _shares_a_parent(Path("/Volumes/A/store"), Path("/Volumes/B/project/fastq")) is False
        assert _shares_a_parent(Path("/mnt/data/store"), Path("/mnt/scratch/project/fastq")) is False
        assert _shares_a_parent(Path("/media/alex/disk1/store"), Path("/media/alex/disk2/proj")) is False

    def test_one_volume_still_gets_a_relative_link(self):
        from prepare.store.link import _shares_a_parent

        assert _shares_a_parent(Path("/Volumes/lab/store"), Path("/Volumes/lab/project/fastq")) is True


def test_link_dataset_copy_mode_interrupted_leaves_no_visible_folder(tmp_path, paths, monkeypatch):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    real_copytree = shutil.copytree

    def _interrupted_copytree(src, dst, *args, **kwargs):
        real_copytree(src, dst, *args, **kwargs)
        raise KeyboardInterrupt

    monkeypatch.setattr("prepare.store.link.shutil.copytree", _interrupted_copytree)
    with pytest.raises(KeyboardInterrupt):
        link_dataset(project_fastq, "SRR1", paths, mode="copy")

    assert not (project_fastq / "SRR1").exists()
    assert list(project_fastq.iterdir()) == []
    assert not accession_has_fastq(project_fastq / "SRR1")


# ------------------------------------------------------------------ copy-mode relink


def test_link_dataset_copy_mode_replaces_an_earlier_store_copy(tmp_path, paths):
    """A second copy-mode link of the same accession refreshes the project's earlier copy."""
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths, mode="copy")
    (paths.sra / "SRR1" / "SRR1_1.fastq").write_text("@r\nGGGG\n+\nIIII\n")

    copied = link_dataset(project_fastq, "SRR1", paths, mode="copy")

    assert copied.is_dir() and not copied.is_symlink()
    assert (copied / "SRR1_1.fastq").read_text() == "@r\nGGGG\n+\nIIII\n"
    # Neither the staging copy nor the earlier copy is left behind.
    assert [entry.name for entry in project_fastq.iterdir()] == ["SRR1"]


def test_link_dataset_copy_mode_still_refuses_a_folder_without_a_sidecar(tmp_path, paths):
    """A real folder with no ``<ACC>.json`` holds the project's own reads and is never replaced."""
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    real_dir = project_fastq / "SRR1"
    real_dir.mkdir(parents=True)
    (real_dir / "SRR1_1.fastq").write_text("@own\nTTTT\n+\nIIII\n")

    with pytest.raises(DataAccessError):
        link_dataset(project_fastq, "SRR1", paths, mode="copy")

    assert (real_dir / "SRR1_1.fastq").read_text() == "@own\nTTTT\n+\nIIII\n"
    assert [entry.name for entry in project_fastq.iterdir()] == ["SRR1"]


def test_copy_relink_interrupted_during_the_copy_keeps_the_earlier_copy(tmp_path, paths, monkeypatch):
    """The earlier copy is moved aside only after the new copy is complete."""
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths, mode="copy")
    real_copytree = shutil.copytree

    def _interrupted_copytree(src, dst, *args, **kwargs):
        real_copytree(src, dst, *args, **kwargs)
        raise KeyboardInterrupt

    monkeypatch.setattr("prepare.store.link.shutil.copytree", _interrupted_copytree)
    with pytest.raises(KeyboardInterrupt):
        link_dataset(project_fastq, "SRR1", paths, mode="copy")

    assert (project_fastq / "SRR1" / "SRR1_1.fastq").read_text() == "@r\nACGT\n+\nIIII\n"
    assert (project_fastq / "SRR1" / "SRR1.json").is_file()
    assert [entry.name for entry in project_fastq.iterdir()] == ["SRR1"]


def test_copy_relink_whose_swap_fails_puts_the_earlier_copy_back(tmp_path, paths, monkeypatch):
    """A failed swap leaves the project with its earlier copy, not with nothing."""
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths, mode="copy")
    real_replace = os.replace
    calls = []

    def _failing_swap(src, dst):
        calls.append((src, dst))
        # The second rename is the staged copy moving into place.
        if len(calls) == 2:
            raise OSError("simulated rename failure")
        return real_replace(src, dst)

    monkeypatch.setattr("prepare.store.link.os.replace", _failing_swap)
    with pytest.raises(DataAccessError, match="simulated rename failure"):
        link_dataset(project_fastq, "SRR1", paths, mode="copy")

    assert len(calls) == 3  # aside, the failed swap, and the earlier copy renamed back
    assert (project_fastq / "SRR1" / "SRR1_1.fastq").read_text() == "@r\nACGT\n+\nIIII\n"
    assert [entry.name for entry in project_fastq.iterdir()] == ["SRR1"]


def test_copy_relink_interrupted_between_aside_and_swap_puts_the_earlier_copy_back(tmp_path, paths, monkeypatch):
    """An interrupt after the earlier copy was moved aside, before the swap, leaves it visible again."""
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths, mode="copy")
    real_replace = os.replace
    calls = []

    def _interrupted_swap(src, dst):
        calls.append((src, dst))
        if len(calls) == 2:
            raise KeyboardInterrupt
        return real_replace(src, dst)

    monkeypatch.setattr("prepare.store.link.os.replace", _interrupted_swap)
    with pytest.raises(KeyboardInterrupt):
        link_dataset(project_fastq, "SRR1", paths, mode="copy")

    assert len(calls) == 3  # aside, the interrupted swap, and the earlier copy renamed back
    assert (project_fastq / "SRR1" / "SRR1_1.fastq").read_text() == "@r\nACGT\n+\nIIII\n"
    assert (project_fastq / "SRR1" / "SRR1.json").is_file()
    assert [entry.name for entry in project_fastq.iterdir()] == ["SRR1"]


def test_copy_relink_without_room_names_the_accession_and_keeps_the_earlier_copy(tmp_path, paths, monkeypatch):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"
    link_dataset(project_fastq, "SRR1", paths, mode="copy")
    full = shutil._ntuple_diskusage(1000, 1000, 1)
    monkeypatch.setattr("prepare.store.link.shutil.disk_usage", lambda _path: full)

    with pytest.raises(DataAccessError, match=r"Cannot copy SRR1 into .*1 bytes free, the copy needs about \d+"):
        link_dataset(project_fastq, "SRR1", paths, mode="copy")

    assert (project_fastq / "SRR1" / "SRR1_1.fastq").read_text() == "@r\nACGT\n+\nIIII\n"
    assert [entry.name for entry in project_fastq.iterdir()] == ["SRR1"]


def test_copy_link_whose_copy_fails_raises_a_data_access_error_naming_the_accession(tmp_path, paths, monkeypatch):
    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"

    def _failing_copytree(src, dst, *args, **kwargs):
        raise shutil.Error([(str(src), str(dst), "No space left on device")])

    monkeypatch.setattr("prepare.store.link.shutil.copytree", _failing_copytree)
    with pytest.raises(DataAccessError, match="Cannot copy SRR1 into"):
        link_dataset(project_fastq, "SRR1", paths, mode="copy")

    assert list(project_fastq.iterdir()) == []


def test_room_shortfall_assumes_room_when_the_filesystem_cannot_be_measured(tmp_path, monkeypatch, caplog):
    from prepare.store.link import room_shortfall

    def _unreadable(_path):
        raise OSError("not mounted")

    monkeypatch.setattr("prepare.store.link.shutil.disk_usage", _unreadable)
    assert room_shortfall(tmp_path, 10**12) is None
    assert "Could not check free space" in caplog.text


# ------------------------------------------------------------------ dataset lock


def _hold_lock_elsewhere(paths, accession="SRR1"):
    """Write a live lock file for ``accession`` held by a process on another host."""
    import json

    from prepare.store.layout import lock_path

    lock = lock_path(paths, accession)
    lock.parent.mkdir(parents=True, exist_ok=True)
    lock.write_text(json.dumps({"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00"}))
    return lock


@pytest.mark.parametrize("mode", ["auto", "copy"])
def test_link_dataset_refuses_while_another_run_holds_the_dataset_lock(tmp_path, paths, mode):
    _store_dataset(paths)
    _hold_lock_elsewhere(paths)
    project_fastq = tmp_path / "project" / "fastq"

    with pytest.raises(DataAccessError, match="SRR1 is being written by .*otherhost.*; retry later"):
        link_dataset(project_fastq, "SRR1", paths, mode=mode)

    assert not os.path.lexists(project_fastq / "SRR1")


def test_link_dataset_with_locked_links_while_the_caller_holds_the_lock(tmp_path, paths):
    from prepare.store.locks import dataset_lock

    _store_dataset(paths)
    project_fastq = tmp_path / "project" / "fastq"

    with dataset_lock(paths, "SRR1"):
        link = link_dataset(project_fastq, "SRR1", paths, locked=True)

    assert link.is_symlink()
    assert link.resolve() == (paths.sra / "SRR1").resolve()


def test_link_dataset_releases_the_dataset_lock_it_took(tmp_path, paths):
    from prepare.store.layout import lock_path

    _store_dataset(paths)
    link_dataset(tmp_path / "project" / "fastq", "SRR1", paths)

    assert not lock_path(paths, "SRR1").exists()
