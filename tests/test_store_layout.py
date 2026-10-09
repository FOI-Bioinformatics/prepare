"""
Tests for prepare.store.layout: on-disk layout of the shared data store.
"""

import json

from prepare.store.layout import (
    StorePaths,
    init_store,
    lock_path,
    read_marker,
    sidecar_path,
    sra_dir,
    store_paths,
)
from prepare.core.constants import STORE_LAYOUT, STORE_MARKER


def test_two_first_time_inits_keep_the_first_marker(tmp_path, monkeypatch):
    """Both runs find no marker; the one that publishes second must not overwrite the first."""
    import threading

    from prepare.store import layout as layout_module

    root = tmp_path / "store"
    real_write = layout_module.write_text_atomic
    both_checked = threading.Barrier(2, timeout=5)
    first_done = threading.Event()
    ids = {}

    def _write(path, text, **kwargs):
        # Both runs have found no marker; the second writes only once the first has finished.
        both_checked.wait()
        if threading.current_thread().name == "second":
            assert first_done.wait(timeout=5)
        return real_write(path, text, **kwargs)

    monkeypatch.setattr(layout_module, "write_text_atomic", _write)

    def _init(name):
        init_store(root)
        ids[name] = read_marker(root)["id"]
        if name == "first":
            first_done.set()

    threads = [threading.Thread(target=_init, args=(name,), name=name) for name in ("first", "second")]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=10)

    assert ids["first"] == ids["second"] == read_marker(root)["id"]
    assert sorted(p.name for p in root.iterdir() if p.name.startswith(".")) == []


def test_store_paths_layout(tmp_path):
    root = tmp_path / "store"
    paths = store_paths(root)

    assert isinstance(paths, StorePaths)
    assert paths.root == root
    assert paths.marker == root / STORE_MARKER
    assert paths.catalog == root / "catalog.sqlite"
    assert paths.catalog_lock == root / "catalog.sqlite.lock"
    assert paths.locks == root / "locks"
    assert paths.tmp == root / "tmp"
    assert paths.sra == root / "sra"
    assert paths.metadata == root / "metadata"


def test_init_store_creates_folders_and_marker(tmp_path):
    root = tmp_path / "store"
    paths = init_store(root)

    assert paths.root.is_dir()
    assert paths.locks.is_dir()
    assert paths.tmp.is_dir()
    assert paths.sra.is_dir()
    assert paths.metadata.is_dir()
    assert paths.marker.is_file()

    marker = json.loads(paths.marker.read_text())
    assert marker["version"] == 1
    assert marker["layout"] == STORE_LAYOUT
    assert "id" in marker
    assert "created" in marker


def test_init_store_is_idempotent_and_keeps_id(tmp_path):
    root = tmp_path / "store"
    first = init_store(root)
    first_marker = json.loads(first.marker.read_text())

    second = init_store(root)
    second_marker = json.loads(second.marker.read_text())

    assert first_marker["id"] == second_marker["id"]
    assert first_marker["created"] == second_marker["created"]


def test_read_marker_returns_none_when_absent(tmp_path):
    root = tmp_path / "store"
    assert read_marker(root) is None


def test_read_marker_returns_dict_when_present(tmp_path):
    root = tmp_path / "store"
    init_store(root)

    marker = read_marker(root)
    assert marker is not None
    assert marker["version"] == 1
    assert marker["layout"] == STORE_LAYOUT


def test_sra_dir_and_sidecar_and_lock_paths(tmp_path):
    root = tmp_path / "store"
    paths = store_paths(root)
    acc = "SRR11011981"

    assert sra_dir(paths, acc) == paths.sra / acc
    assert sidecar_path(paths, acc) == paths.sra / acc / f"{acc}.json"
    assert lock_path(paths, acc) == paths.locks / f"{acc}.lock"
