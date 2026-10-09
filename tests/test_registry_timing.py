"""Tests for the timing fields of the download, extraction and assembly blocks."""

import json
import re

import pytest

from prepare.data import registry as reg
from prepare.data import registry_blocks as rb
from prepare.data.registry_timing import (
    Stopwatch,
    set_assembly_timing,
    set_download_timing,
    set_extraction_timing,
    timing_summary,
)

ISO_SECONDS = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\+00:00$")


def _registry(tmp_path):
    r = reg.load_registry(tmp_path / "prepare_registry.json")
    (tmp_path / "fastq" / "SRR1").mkdir(parents=True)
    (tmp_path / "fastq" / "SRR1" / "SRR1_1.fastq").write_text("@r\nA\n+\nI\n")
    reg.record_download(r, "SRR1", "downloaded", tmp_path / "fastq", "ok")
    reg.record_extraction(r, "SRR1", "GCF_A", [], 10, False, {})
    reg.record_assembly(r, "SRR1", "GCF_A", tmp_path / "asm", {"contigs": 2}, "v1", {})
    return r


@pytest.mark.parametrize(
    "cls, data",
    [
        (rb.DownloadBlock, {"attempts": 1, "state": "downloaded", "date": "d", "files": [], "bytes_total": 0}),
        (rb.ExtractionBlock, {"date": "d", "mapped_reads": 3, "files": ["a"], "assembly": None}),
        (rb.AssemblyBlock, {"date": "d", "dir": "x", "contigs": 1, "N90": 5}),
    ],
)
def test_block_without_timing_round_trips_byte_identical(cls, data):
    before = json.dumps(data, separators=(",", ":"))
    block = cls.from_dict(json.loads(before))
    assert json.dumps(block.to_dict(), separators=(",", ":")) == before
    assert block.started is None and block.seconds is None


def test_fresh_blocks_leave_out_unset_timing():
    for cls in (rb.DownloadBlock, rb.ExtractionBlock, rb.AssemblyBlock):
        out = cls().to_dict()
        assert "started" not in out and "seconds" not in out


def test_block_with_timing_round_trips():
    data = {"attempts": 1, "state": "downloaded", "started": "2026-10-01T10:00:00+00:00", "seconds": 12.5}
    assert rb.DownloadBlock.from_dict(data).to_dict() == data


def test_set_download_timing_writes_and_clears(tmp_path):
    r = _registry(tmp_path)
    before = json.dumps(r.datasets["SRR1"]["download"])
    set_download_timing(r, "SRR1", "2026-10-01T10:00:00+00:00", 12.5)
    download = r.datasets["SRR1"]["download"]
    assert (download["started"], download["seconds"]) == ("2026-10-01T10:00:00+00:00", 12.5)
    set_download_timing(r, "SRR1", None, None)
    assert json.dumps(r.datasets["SRR1"]["download"]) == before


def test_set_extraction_and_assembly_timing(tmp_path):
    r = _registry(tmp_path)
    set_extraction_timing(r, "SRR1", "GCF_A", "2026-10-01T10:00:00+00:00", 3.25)
    set_assembly_timing(r, "SRR1", "GCF_A", "2026-10-01T11:00:00+00:00", 60.0)
    extraction = r.datasets["SRR1"]["extractions"]["GCF_A"]
    assert (extraction["started"], extraction["seconds"]) == ("2026-10-01T10:00:00+00:00", 3.25)
    assert (extraction["assembly"]["started"], extraction["assembly"]["seconds"]) == (
        "2026-10-01T11:00:00+00:00",
        60.0,
    )
    assert extraction["assembly"]["contigs"] == 2
    set_assembly_timing(r, "SRR1", "GCF_A", None, None)
    set_extraction_timing(r, "SRR1", "GCF_A", None, None)
    extraction = r.datasets["SRR1"]["extractions"]["GCF_A"]
    assert "seconds" not in extraction and "seconds" not in extraction["assembly"]
    assert "started" not in extraction and "started" not in extraction["assembly"]


def test_timing_setters_leave_unrecorded_blocks_alone(tmp_path):
    r = reg.load_registry(tmp_path / "prepare_registry.json")
    set_download_timing(r, "SRR9", "2026-10-01T10:00:00+00:00", 1.0)
    set_extraction_timing(r, "SRR9", "GCF_A", "2026-10-01T10:00:00+00:00", 1.0)
    set_assembly_timing(r, "SRR9", "GCF_A", "2026-10-01T10:00:00+00:00", 1.0)
    assert r.datasets == {}
    reg.record_extraction(r, "SRR9", "GCF_A", [], 0, False, {})
    set_assembly_timing(r, "SRR9", "GCF_A", "2026-10-01T10:00:00+00:00", 1.0)
    assert r.datasets["SRR9"]["extractions"]["GCF_A"]["assembly"] is None


def test_stopwatch_laps(monkeypatch):
    ticks = iter([100.0, 102.5, 102.5, 110.0, 110.0])
    monkeypatch.setattr("prepare.data.registry_timing.monotonic", lambda: next(ticks))
    watch = Stopwatch()
    started, seconds = watch.lap()
    assert ISO_SECONDS.match(started) and seconds == 2.5
    assert watch.lap()[1] == 7.5


def test_timing_summary_counts_only_downloaded_or_failed_blocks(tmp_path):
    r = _registry(tmp_path)
    set_download_timing(r, "SRR1", "t", 5.0)
    for acc, state in (("SRR2", "missing"), ("SRR3", "skipped")):
        reg.record_download(r, acc, state, tmp_path / "fastq", "x")
        set_download_timing(r, acc, "t", 42.0)
    summary = timing_summary(r)
    assert summary["downloads_timed"] == 1 and summary["download_seconds_total"] == 5.0


def test_empty_totals_are_floats(tmp_path):
    summary = timing_summary(_registry(tmp_path))
    totals = [value for key, value in summary.items() if key.endswith("_seconds_total")]
    assert totals == [0.0, 0.0, 0.0] and all(isinstance(value, float) for value in totals)
    assert json.dumps(summary["download_seconds_total"]) == "0.0"


def test_record_skip_clears_an_earlier_attempt_time(tmp_path):
    from prepare.cli.commands.sra import DownloadSraCommand

    r = _registry(tmp_path)
    reg.record_download(r, "SRR2", "failed", tmp_path / "fastq", "x")
    set_download_timing(r, "SRR2", "t", 3.0)
    DownloadSraCommand._record_skip(r, "SRR2", "blacklisted", tmp_path / "fastq")
    block = r.datasets["SRR2"]["download"]
    assert block["state"] == "skipped" and "seconds" not in block and "started" not in block


def test_timing_summary(tmp_path):
    r = _registry(tmp_path)
    assert timing_summary(r) == {
        "downloads_timed": 0,
        "download_seconds_total": 0.0,
        "download_seconds_median": None,
        "extractions_timed": 0,
        "extraction_seconds_total": 0.0,
        "extraction_seconds_median": None,
        "assemblies_timed": 0,
        "assembly_seconds_total": 0.0,
        "assembly_seconds_median": None,
    }
    for acc, seconds in (("SRR1", 10.0), ("SRR2", 20.0), ("SRR3", 40.0)):
        reg.record_download(r, acc, "failed", tmp_path / "fastq", "x")
        set_download_timing(r, acc, "2026-10-01T10:00:00+00:00", seconds)
    reg.record_extraction(r, "SRR2", "GCF_A", [], 0, False, {})
    set_extraction_timing(r, "SRR1", "GCF_A", "t", 1.5)
    set_extraction_timing(r, "SRR2", "GCF_A", "t", 2.5)
    set_assembly_timing(r, "SRR1", "GCF_A", "t", 30.0)
    summary = timing_summary(r)
    assert summary["downloads_timed"] == 3
    assert summary["download_seconds_total"] == 70.0 and summary["download_seconds_median"] == 20.0
    assert summary["extractions_timed"] == 2
    assert summary["extraction_seconds_total"] == 4.0 and summary["extraction_seconds_median"] == 2.0
    assert summary["assemblies_timed"] == 1 and summary["assembly_seconds_median"] == 30.0
