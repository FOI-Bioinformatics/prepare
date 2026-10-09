"""Tests for the assembly step of extract_target_reads --assemble (cli/commands/extraction_assembly.py)."""

import argparse
import itertools
import json
import logging
import subprocess
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path
from unittest.mock import patch

import pytest

import prepare.cli.commands.extraction_assembly as asm_mod
import prepare.data.registry as registry_mod
from helpers_extraction import _fake_tools, tools_present  # noqa: F401 (autouse fixture)
from prepare.cli.commands.extraction_assembly import assemble_samples
from prepare.cli.commands.read_extraction import ExtractTargetReadsCommand
from prepare.data.assembly_identity import MARKER_NAME, read_marker
from prepare.data.extraction_locks import sample_extraction_lock
from prepare.data.read_extraction import ExtractionResult
from prepare.data.registry import load_registry, save_registry

RUN_SECURE = "prepare.data.read_extraction.SecureSubprocess.run_secure"


@pytest.fixture(autouse=True)
def distinct_registry_dates(monkeypatch):
    """Give every registry date its own second, as a real re-extraction (minutes apart) would have."""
    start = datetime(2026, 10, 1, 10, 0, 0, tzinfo=timezone.utc)
    ticks = itertools.count()
    monkeypatch.setattr(
        registry_mod, "_now", lambda: (start + timedelta(seconds=next(ticks))).isoformat(timespec="seconds")
    )


def _project(root: Path, samples=("SRR1",)) -> Path:
    """Containment table, paired FASTQ files per sample and a genome under ``root``; returns the table."""
    table = root / "parsed_containment.txt"
    table.write_text("\tGCF_1\n" + "".join(f"{acc}\t0.9\n" for acc in samples))
    for acc in samples:
        folder = root / "fastq" / acc
        folder.mkdir(parents=True)
        (folder / f"{acc}_1.fastq.gz").write_text("x")
        (folder / f"{acc}_2.fastq.gz").write_text("x")
    (root / "GCF_1.fna").write_text(">s\nACGT\n")
    return table


def _args(root: Path, **overrides) -> argparse.Namespace:
    base = dict(
        parsed_containment=str(root / "parsed_containment.txt"),
        genome_id="GCF_1",
        genome_fasta=str(root / "GCF_1.fna"),
        fastq_folder=str(root / "fastq"),
        output_folder=str(root / "targeted"),
        threshold=0.5,
        preset="sr",
        threads=4,
        min_mapq=0,
        temp_folder=None,
        allow_truncated=False,
        debug_keep_sam=False,
        assemble=True,
        assembly_threads=None,
        min_contig_len=None,
        assembly_preset="meta-sensitive",
        keep_intermediate=False,
        no_coverage=True,
        dry_run=False,
        force=False,
        registry=str(root / "registry.json"),
        data_root=None,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _megahit_runs(state, accession=None):
    """The megahit assembly calls (not ``--version``) the fake recorded, optionally for one sample."""
    runs = [a for exe, a in state.get("calls", []) if exe == "megahit" and "-o" in a]
    if accession is not None:
        runs = [a for a in runs if f"/{accession}/" in a[a.index("-o") + 1]]
    return runs


def _assembly(root: Path, accession="SRR1"):
    data = json.loads((root / "registry.json").read_text())
    return data["datasets"][accession]["extractions"]["GCF_1"].get("assembly")


def _out_dir(root: Path, accession="SRR1") -> Path:
    return root / "targeted" / accession / "GCF_1_assembly"


@patch(RUN_SECURE)
def test_reextraction_with_a_new_min_mapq_redoes_the_assembly(mock_run, tmp_path):
    state = {}
    mock_run.side_effect = _fake_tools(state)
    _project(tmp_path)
    assert ExtractTargetReadsCommand().execute(_args(tmp_path)) == 0
    first = _assembly(tmp_path)
    assert first["inputs"]["extraction_date"] is not None

    assert ExtractTargetReadsCommand().execute(_args(tmp_path, min_mapq=20)) == 0

    second = _assembly(tmp_path)
    assert len(_megahit_runs(state)) == 2
    assert second["inputs"]["extraction_date"] != first["inputs"]["extraction_date"]
    marker = read_marker(_out_dir(tmp_path))
    assert marker["inputs"] == second["inputs"]


@patch(RUN_SECURE)
def test_legacy_record_predating_the_extraction_is_redone(mock_run, tmp_path):
    """An unmarked folder whose (pre-marker) record is older than the extraction is assembled again."""
    state = {}
    mock_run.side_effect = _fake_tools(state)
    _project(tmp_path)
    assert ExtractTargetReadsCommand().execute(_args(tmp_path)) == 0
    (_out_dir(tmp_path) / MARKER_NAME).unlink()
    registry = load_registry(tmp_path / "registry.json")
    block = registry.datasets["SRR1"]["extractions"]["GCF_1"]["assembly"]
    block.pop("inputs")
    block["date"] = "2000-01-01T00:00:00+00:00"
    save_registry(registry)

    assert ExtractTargetReadsCommand().execute(_args(tmp_path)) == 0

    assert len(_megahit_runs(state)) == 2
    assert (_out_dir(tmp_path) / MARKER_NAME).is_file()
    assert _assembly(tmp_path)["date"] != "2000-01-01T00:00:00+00:00"


@patch(RUN_SECURE)
def test_legacy_record_matching_the_run_is_accepted_with_an_unknown_version(mock_run, tmp_path):
    """An unmarked folder with a current legacy record is kept; its marker does not claim this megahit built it."""
    state = {}
    mock_run.side_effect = _fake_tools(state)
    _project(tmp_path)
    assert ExtractTargetReadsCommand().execute(_args(tmp_path)) == 0
    (_out_dir(tmp_path) / MARKER_NAME).unlink()
    registry = load_registry(tmp_path / "registry.json")
    registry.datasets["SRR1"]["extractions"]["GCF_1"]["assembly"].pop("inputs")
    save_registry(registry)

    with patch.object(asm_mod, "megahit_version", return_value="MEGAHIT v9.9.9"):
        assert ExtractTargetReadsCommand().execute(_args(tmp_path)) == 0

    assert len(_megahit_runs(state)) == 1
    marker = read_marker(_out_dir(tmp_path))
    assert marker["megahit_version"] == ""
    assert _assembly(tmp_path)["inputs"] == marker["inputs"]


@patch(RUN_SECURE)
def test_contig_less_folder_is_redone(mock_run, tmp_path, caplog):
    state = {}
    mock_run.side_effect = _fake_tools(state)
    _project(tmp_path)
    assert ExtractTargetReadsCommand().execute(_args(tmp_path, assemble=False)) == 0
    # An interrupted megahit (before staging existed) left the folder without contigs.
    _out_dir(tmp_path).mkdir(parents=True)
    (_out_dir(tmp_path) / "log").write_text("killed\n")

    with caplog.at_level(logging.WARNING):
        assert ExtractTargetReadsCommand().execute(_args(tmp_path)) == 0

    assert len(_megahit_runs(state)) == 1
    assert (_out_dir(tmp_path) / "final.contigs.fa").is_file()
    assert "interrupted run" in caplog.text
    assert _assembly(tmp_path)["contigs"] > 0


@pytest.mark.parametrize(
    "error",
    [subprocess.CalledProcessError(1, ["megahit"], stderr="out of memory"), OSError("Permission denied")],
    ids=["megahit-failed", "os-error"],
)
@patch(RUN_SECURE)
def test_one_failed_sample_does_not_stop_the_next_and_the_run_exits_1(mock_run, error, tmp_path, caplog):
    state = {}
    fake = _fake_tools(state)

    def run(executable, args, **kwargs):
        if executable == "megahit" and "-o" in args and "/SRR1/" in args[args.index("-o") + 1]:
            raise error
        return fake(executable, args, **kwargs)

    mock_run.side_effect = run
    _project(tmp_path, samples=("SRR1", "SRR2"))
    with caplog.at_level(logging.INFO):
        rc = ExtractTargetReadsCommand().execute(_args(tmp_path))

    assert rc == 1
    assert _assembly(tmp_path, "SRR1") is None
    assert _assembly(tmp_path, "SRR2")["contigs"] > 0
    errors = [r.getMessage() for r in caplog.records if r.levelname == "ERROR"]
    assert any("SRR1" in line and "SRR2" not in line for line in errors)
    summary = [line for line in errors if line.startswith("Assembly failed for 1 sample(s): SRR1 (")]
    assert len(summary) == 1 and "\n" not in summary[0]
    assert not list((tmp_path / "targeted" / "SRR1").glob(".GCF_1_assembly.*"))


@patch(RUN_SECURE)
def test_failed_coverage_mapping_does_not_stop_the_next_sample(mock_run, tmp_path, caplog):
    """minimap2 fails while mapping SRR1's reads onto its new contigs; SRR2 is still assembled."""
    state = {}
    fake = _fake_tools(state)

    def run(executable, args, **kwargs):
        if executable == "minimap2" and any("/SRR1/GCF_1_assembly/" in str(a) for a in args):
            raise subprocess.CalledProcessError(1, ["minimap2"], stderr="mapping failed\nsecond line")
        return fake(executable, args, **kwargs)

    mock_run.side_effect = run
    _project(tmp_path, samples=("SRR1", "SRR2"))
    with caplog.at_level(logging.INFO):
        rc = ExtractTargetReadsCommand().execute(_args(tmp_path, no_coverage=False))

    assert rc == 1
    assert len(_megahit_runs(state)) == 2
    assert _assembly(tmp_path, "SRR1") is None
    assert _assembly(tmp_path, "SRR2")["contigs"] > 0
    assert any(r.levelname == "ERROR" and r.getMessage().startswith("SRR1: assembly") for r in caplog.records)


@patch(RUN_SECURE)
def test_sample_locked_elsewhere_is_skipped_as_busy(mock_run, tmp_path, caplog):
    state = {}
    mock_run.side_effect = _fake_tools(state)
    _project(tmp_path, samples=("SRR1", "SRR2"))
    args = _args(tmp_path)
    cmd = ExtractTargetReadsCommand()
    assert cmd.execute(_args(tmp_path, assemble=False)) == 0
    reads = {acc: sorted((tmp_path / "targeted" / acc).glob("GCF_1_*.fastq.gz")) for acc in ("SRR1", "SRR2")}
    results = {acc: ExtractionResult(files, 10, False, skipped=True) for acc, files in reads.items()}

    held, release = threading.Event(), threading.Event()

    def hold_lock():
        # Another worker (in a real run, another process) is assembling SRR1.
        with sample_extraction_lock(tmp_path / "targeted", "SRR1", "GCF_1"):
            held.set()
            release.wait(10)

    holder = threading.Thread(target=hold_lock)
    holder.start()
    try:
        assert held.wait(10)
        with caplog.at_level(logging.INFO):
            outcome = assemble_samples(cmd, args, reads, results)
    finally:
        release.set()
        holder.join()

    assert outcome.busy == ["SRR1"]
    assert outcome.assembled == ["SRR2"]
    assert outcome.failed == {}
    assert _megahit_runs(state, "SRR1") == []
    assert not _out_dir(tmp_path, "SRR1").exists()
    assert any(
        r.levelname == "INFO" and "SRR1" in r.getMessage() and "elsewhere" in r.getMessage() for r in caplog.records
    )


@patch(RUN_SECURE)
def test_stale_staging_folder_is_swept(mock_run, tmp_path):
    state = {}
    mock_run.side_effect = _fake_tools(state)
    _project(tmp_path)
    assert ExtractTargetReadsCommand().execute(_args(tmp_path)) == 0
    leftover = tmp_path / "targeted" / "SRR1" / ".GCF_1_assembly.node1.4242.abcd.tmp"
    leftover.mkdir()
    (leftover / "final.contigs.fa").write_text(">half\nAC\n")

    assert ExtractTargetReadsCommand().execute(_args(tmp_path)) == 0

    assert not leftover.exists()
    assert len(_megahit_runs(state)) == 1


@patch(RUN_SECURE)
def test_unchanged_inputs_do_not_run_megahit(mock_run, tmp_path, caplog):
    state = {}
    mock_run.side_effect = _fake_tools(state)
    _project(tmp_path)
    assert ExtractTargetReadsCommand().execute(_args(tmp_path)) == 0
    first = _assembly(tmp_path)

    with caplog.at_level(logging.INFO):
        assert ExtractTargetReadsCommand().execute(_args(tmp_path)) == 0

    assert len(_megahit_runs(state)) == 1
    assert _assembly(tmp_path) == first
    assert "reused 1" in caplog.text


@patch(RUN_SECURE)
def test_force_redoes_a_current_assembly(mock_run, tmp_path):
    state = {}
    mock_run.side_effect = _fake_tools(state)
    _project(tmp_path)
    assert ExtractTargetReadsCommand().execute(_args(tmp_path)) == 0
    assert ExtractTargetReadsCommand().execute(_args(tmp_path, force=True)) == 0
    assert len(_megahit_runs(state)) == 2
    assert _assembly(tmp_path)["contigs"] > 0


@patch(RUN_SECURE)
def test_marked_folder_without_a_record_is_recorded_from_its_marker(mock_run, tmp_path):
    """The registry lost the assembly but the folder's marker matches: megahit is not rerun and the
    record takes its version and parameters from the marker, not from the megahit installed now."""
    state = {}
    mock_run.side_effect = _fake_tools(state)
    _project(tmp_path)
    with patch.object(asm_mod, "megahit_version", return_value="MEGAHIT v1.2.9"):
        assert ExtractTargetReadsCommand().execute(_args(tmp_path, min_contig_len=300)) == 0
    registry = load_registry(tmp_path / "registry.json")
    registry.datasets["SRR1"]["extractions"]["GCF_1"]["assembly"] = None
    save_registry(registry)

    with patch.object(asm_mod, "megahit_version", return_value="MEGAHIT v9.9.9"):
        assert ExtractTargetReadsCommand().execute(_args(tmp_path, min_contig_len=300, assembly_threads=2)) == 0

    assert len(_megahit_runs(state)) == 1
    record = _assembly(tmp_path)
    assert record["version"] == "MEGAHIT v1.2.9"
    assert record["params"]["min_contig_len"] == 300
    assert record["params"]["threads"] != 2
    assert record["inputs"] == read_marker(_out_dir(tmp_path))["inputs"]
    assert "seconds" not in record
