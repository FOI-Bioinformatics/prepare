"""Tests for the shared size walker and the call sites that use it."""

import os

import pytest

from prepare.core.exceptions import DataAccessError
from prepare.data.sra.cleanup import transient_bytes
from prepare.utils.sizes import tree_bytes


def test_tree_bytes_counts_file_and_directory(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "a" / "x.txt").write_bytes(b"12345")
    (tmp_path / "y.txt").write_bytes(b"123")
    assert tree_bytes(tmp_path) == 8
    assert tree_bytes(tmp_path / "y.txt") == 3
    assert tree_bytes(tmp_path / "missing") == 0


def test_tree_bytes_ignores_symlink_pointing_outside(tmp_path):
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "big.bin").write_bytes(b"x" * 1000)
    inside = tmp_path / "inside"
    inside.mkdir()
    (inside / "own.txt").write_bytes(b"abc")
    os.symlink(outside / "big.bin", inside / "file_link")
    os.symlink(outside, inside / "dir_link")
    assert tree_bytes(inside) == 3
    assert tree_bytes(inside / "file_link") == 0


def test_transient_bytes_sums_transient_folders_only(tmp_path):
    temp = tmp_path / "SRR1_temp"
    temp.mkdir()
    (temp / "part").write_bytes(b"1234")
    keep = tmp_path / "SRR1"
    keep.mkdir()
    (keep / "reads.fastq").write_bytes(b"x" * 50)
    outside = tmp_path / "elsewhere.bin"
    outside.write_bytes(b"y" * 99)
    os.symlink(outside, temp / "link")
    assert transient_bytes(tmp_path) == 4


def test_gc_and_link_use_the_shared_walker():
    from prepare.cli.commands.store import gc
    from prepare.store import link

    assert gc.tree_bytes is tree_bytes
    assert link.tree_bytes is tree_bytes
    assert not hasattr(gc, "_path_bytes")
    assert not hasattr(link, "_tree_bytes")


def test_copy_into_place_needs_size_of_dataset_tree(tmp_path, monkeypatch):
    from prepare.store import link

    store_dataset = tmp_path / "SRR1"
    store_dataset.mkdir()
    (store_dataset / "r.fastq").write_bytes(b"x" * 10)
    seen = []

    def fake_shortfall(location, needed):
        seen.append(needed)
        return 0

    monkeypatch.setattr(link, "room_shortfall", fake_shortfall)
    with pytest.raises(DataAccessError):
        link._copy_into_place(store_dataset, tmp_path / "proj" / "SRR1")
    assert seen == [10]
