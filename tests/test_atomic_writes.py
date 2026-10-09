"""Tests for the atomic write helpers in prepare.data.file_io and the state files written through them."""

import gzip
import json
import os
import stat
import threading
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from prepare.data import file_io
from prepare.data.file_io import atomic_path, open_atomic, unique_temp_path, write_csv, write_text_atomic


def _no_temp_files(folder: Path) -> bool:
    return not [p for p in folder.iterdir() if p.name.endswith(".tmp")]


def test_unique_temp_path_is_hidden_and_distinct(tmp_path):
    target = tmp_path / "out.tsv"
    names = {unique_temp_path(target).name for _ in range(1000)}
    assert len(names) == 1000
    for name in names:
        assert name.startswith(".out.tsv.")
        assert name.endswith(".tmp")
        assert f".{os.getpid()}." in name
    assert unique_temp_path(target).parent == tmp_path


def test_unique_temp_path_uses_the_first_host_label(tmp_path):
    with patch.object(file_io.socket, "gethostname", return_value="node17.cluster.example.org"):
        name = unique_temp_path(tmp_path / "a.json").name
    assert name.startswith(".a.json.node17.")


def test_write_text_atomic_writes_and_returns_path(tmp_path):
    target = tmp_path / "sub" / "note.txt"
    assert write_text_atomic(target, "hello\n") == target
    assert target.read_text() == "hello\n"
    assert _no_temp_files(target.parent)


def test_write_text_atomic_interrupted_keeps_old_content(tmp_path):
    target = tmp_path / "note.txt"
    target.write_text("old\n")
    with patch.object(file_io.os, "replace", side_effect=KeyboardInterrupt):
        with pytest.raises(KeyboardInterrupt):
            write_text_atomic(target, "new\n")
    assert target.read_text() == "old\n"
    assert _no_temp_files(tmp_path)


def test_atomic_path_error_in_body_leaves_no_temp(tmp_path):
    target = tmp_path / "x.txt"
    with pytest.raises(RuntimeError):
        with atomic_path(target) as tmp:
            tmp.write_text("partial")
            raise RuntimeError("boom")
    assert not target.exists()
    assert _no_temp_files(tmp_path)


def test_open_atomic_streams_and_publishes_on_close(tmp_path):
    target = tmp_path / "rows.csv"
    target.write_text("old\n")
    with open_atomic(target, newline="") as handle:
        handle.write("a,b\r\n")
        assert target.read_text() == "old\n"
    assert target.read_bytes() == b"a,b\r\n"
    assert _no_temp_files(tmp_path)


def test_open_atomic_error_keeps_old_content(tmp_path):
    target = tmp_path / "rows.bin"
    target.write_bytes(b"old")
    with pytest.raises(ValueError):
        with open_atomic(target, "wb") as handle:
            handle.write(b"partial")
            raise ValueError("stop")
    assert target.read_bytes() == b"old"
    assert _no_temp_files(tmp_path)


def test_open_atomic_refuses_append(tmp_path):
    with pytest.raises(ValueError):
        with open_atomic(tmp_path / "x", "a"):
            pass


def test_write_csv_interrupted_keeps_old_content(tmp_path):
    target = tmp_path / "table.tsv"
    write_csv(pd.DataFrame({"a": [1]}), target, sep="\t", index=False)
    before = target.read_text()
    with patch.object(file_io.os, "replace", side_effect=KeyboardInterrupt):
        with pytest.raises(KeyboardInterrupt):
            write_csv(pd.DataFrame({"a": [2, 3]}), target, sep="\t", index=False)
    assert target.read_text() == before
    assert _no_temp_files(tmp_path)


def test_write_csv_gz_round_trips(tmp_path):
    target = tmp_path / "table.tsv.gz"
    df = pd.DataFrame({"acc": ["SRR1", "SRR2"], "value": [0.5, 0.25]})
    write_csv(df, target, sep="\t", index=False)
    with gzip.open(target, "rt") as handle:
        assert handle.readline().strip() == "acc\tvalue"
    pd.testing.assert_frame_equal(pd.read_csv(target, sep="\t"), df)


def test_write_csv_explicit_compression_is_kept(tmp_path):
    target = tmp_path / "table.tsv"
    write_csv(pd.DataFrame({"a": [1]}), target, sep="\t", index=False, compression="gzip")
    with gzip.open(target, "rt") as handle:
        assert handle.read().startswith("a")


def test_existing_mode_is_kept(tmp_path):
    target = tmp_path / "shared.json"
    target.write_text("{}")
    os.chmod(target, 0o640)
    write_text_atomic(target, '{"a": 1}')
    assert stat.S_IMODE(target.stat().st_mode) == 0o640


@pytest.mark.parametrize("fsync", [False, True])
def test_a_read_only_target_is_replaced_and_stays_read_only(tmp_path, fsync):
    target = tmp_path / "sidecar.json"
    target.write_text("{}")
    os.chmod(target, 0o444)
    write_text_atomic(target, '{"a": 1}', fsync=fsync)
    assert target.read_text() == '{"a": 1}'
    assert stat.S_IMODE(target.stat().st_mode) == 0o444
    assert _no_temp_files(tmp_path)


def test_new_file_follows_umask(tmp_path):
    old = os.umask(0o002)
    try:
        target = write_text_atomic(tmp_path / "new.txt", "x")
    finally:
        os.umask(old)
    assert stat.S_IMODE(target.stat().st_mode) == 0o664


def test_symlink_target_is_written_through(tmp_path):
    real = tmp_path / "real.txt"
    real.write_text("old")
    link = tmp_path / "link.txt"
    link.symlink_to(real)
    write_text_atomic(link, "new")
    assert link.is_symlink()
    assert real.read_text() == "new"


def test_replace_retries_on_permission_error(tmp_path):
    target = tmp_path / "x.txt"
    real_replace = os.replace
    calls = {"n": 0}

    def flaky(src, dst):
        calls["n"] += 1
        if calls["n"] < 3:
            raise PermissionError("in use")
        real_replace(src, dst)

    with patch.object(file_io.os, "replace", side_effect=flaky), patch.object(file_io.time, "sleep"):
        write_text_atomic(target, "done")
    assert target.read_text() == "done"
    assert calls["n"] == 3


def test_replace_gives_up_after_retries(tmp_path):
    target = tmp_path / "x.txt"
    with (
        patch.object(file_io.os, "replace", side_effect=PermissionError("in use")),
        patch.object(file_io.time, "sleep"),
    ):
        with pytest.raises(PermissionError):
            write_text_atomic(target, "done")
    assert not target.exists()
    assert _no_temp_files(tmp_path)


def test_reader_never_sees_a_partial_file(tmp_path):
    target = tmp_path / "data.json"
    payloads = [json.dumps({"i": i, "pad": "x" * (5000 + 97 * i)}) for i in range(200)]
    write_text_atomic(target, payloads[0])
    valid = set(payloads)
    bad = []
    stop = threading.Event()

    def reader():
        while not stop.is_set():
            try:
                text = target.read_text()
            except FileNotFoundError:
                bad.append("missing")
                continue
            if text not in valid:
                bad.append(text[:40])

    thread = threading.Thread(target=reader)
    thread.start()
    try:
        for text in payloads:
            write_text_atomic(target, text)
    finally:
        stop.set()
        thread.join()
    assert bad == []
    assert _no_temp_files(tmp_path)


def test_temporary_fastq_names_are_not_listed(tmp_path):
    from prepare.data.registry import scan_downloads
    from prepare.data.sra import fastq_files

    acc_dir = tmp_path / "SRR1"
    acc_dir.mkdir()
    (acc_dir / "SRR1_1.fastq.gz").write_bytes(b"reads")
    unique_temp_path(acc_dir / "SRR1_2.fastq.gz").write_bytes(b"partial")
    assert [p.name for p in fastq_files(acc_dir)] == ["SRR1_1.fastq.gz"]
    assert file_io.visible_files(acc_dir) == [acc_dir / "SRR1_1.fastq.gz"]
    assert scan_downloads(tmp_path) == {"SRR1": (1, 5)}


def test_compress_fastq_fallback_is_atomic(tmp_path):
    from prepare.data.sra import compress_fastq

    source = tmp_path / "SRR1_1.fastq"
    source.write_text("@r\nACGT\n+\nIIII\n")
    with patch("prepare.data.sra.fastq.shutil.which", return_value=None):
        with patch.object(file_io.os, "replace", side_effect=KeyboardInterrupt):
            with pytest.raises(KeyboardInterrupt):
                compress_fastq(source, 1)
        assert source.exists()
        assert sorted(p.name for p in tmp_path.iterdir()) == ["SRR1_1.fastq"]
        target = compress_fastq(source, 1)
    assert target.name == "SRR1_1.fastq.gz"
    assert not source.exists()
    with gzip.open(target, "rt") as handle:
        assert handle.read() == "@r\nACGT\n+\nIIII\n"


def test_sidecar_write_uses_fsync_and_leaves_no_temp(tmp_path):
    from prepare.store.sidecar import Sidecar, read_sidecar, write_sidecar

    target = tmp_path / "SRR1" / "SRR1.json"
    with patch.object(file_io.os, "fsync", wraps=os.fsync) as fsync:
        write_sidecar(target, Sidecar(accession="SRR1", state="complete"))
    assert fsync.called
    assert read_sidecar(target).state == "complete"
    assert _no_temp_files(target.parent)


def test_write_registry_interrupted_leaves_no_temp(tmp_path):
    from prepare.data import registry as reg

    target = tmp_path / reg.REGISTRY_FILENAME
    registry = reg.load_registry(target)
    reg.save_registry(registry, target)
    before = target.read_text()
    registry.datasets["SRR1"] = {}
    with reg._acquire_lock(target.with_name(target.name + ".lock")) as lock:
        with patch.object(file_io.os, "replace", side_effect=KeyboardInterrupt):
            with pytest.raises(KeyboardInterrupt):
                reg._write_registry(registry, target, lock)
    assert target.read_text() == before
    assert _no_temp_files(tmp_path)
    assert not [p for p in tmp_path.iterdir() if ".tmp" in p.name]


def test_write_registry_fsyncs(tmp_path):
    from prepare.data import registry as reg

    target = tmp_path / reg.REGISTRY_FILENAME
    with reg._acquire_lock(target.with_name(target.name + ".lock")) as lock:
        with patch.object(file_io.os, "fsync", wraps=os.fsync) as fsync:
            reg._write_registry(reg.load_registry(target), target, lock)
    assert fsync.called
