"""Tests for the assembly inputs record and the helpers that judge whether an assembly is stale."""

import json

import pytest

from prepare.data import registry as reg
from prepare.data import registry_blocks as rb
from prepare.data.registry_assembly import (
    assembly_predates_extraction,
    legacy_assembly_current,
    set_assembly_inputs,
)

INPUTS = {"reads": [{"name": "SRR1_1.fastq", "size": 10, "sha256": "ab"}], "mapped_reads": 10}


def _registry(tmp_path, assembly_params=None):
    r = reg.load_registry(tmp_path / "prepare_registry.json")
    reg.record_extraction(r, "SRR1", "GCF_A", [], 10, False, {"preset": "sr"})
    reg.record_assembly(r, "SRR1", "GCF_A", tmp_path / "asm", {"contigs": 2}, "v1", assembly_params or {})
    return r


def _set_dates(r, extraction_date, assembly_date):
    raw = r.datasets["SRR1"]["extractions"]["GCF_A"]
    if extraction_date is None:
        raw.pop("date", None)
    else:
        raw["date"] = extraction_date
    if raw.get("assembly") is not None:
        if assembly_date is None:
            raw["assembly"].pop("date", None)
        else:
            raw["assembly"]["date"] = assembly_date


# ------------------------------------------------------------------ set_assembly_inputs


def test_inputs_are_omitted_while_none(tmp_path):
    r = _registry(tmp_path)
    assert "inputs" not in r.datasets["SRR1"]["extractions"]["GCF_A"]["assembly"]
    assert rb.extraction_block(r, "SRR1", "GCF_A").assembly.inputs is None


def test_set_and_drop_inputs(tmp_path):
    r = _registry(tmp_path)
    before = json.dumps(r.datasets["SRR1"]["extractions"]["GCF_A"], separators=(",", ":"))
    set_assembly_inputs(r, "SRR1", "GCF_A", INPUTS)
    assembly = r.datasets["SRR1"]["extractions"]["GCF_A"]["assembly"]
    assert assembly["inputs"] == INPUTS
    assert rb.extraction_block(r, "SRR1", "GCF_A").assembly.inputs == INPUTS
    set_assembly_inputs(r, "SRR1", "GCF_A", None)
    assert "inputs" not in r.datasets["SRR1"]["extractions"]["GCF_A"]["assembly"]
    assert json.dumps(r.datasets["SRR1"]["extractions"]["GCF_A"], separators=(",", ":")) == before


def test_inputs_survive_a_save_and_load(tmp_path):
    r = _registry(tmp_path)
    set_assembly_inputs(r, "SRR1", "GCF_A", INPUTS)
    reg.save_registry(r)
    loaded = reg.load_registry(tmp_path / "prepare_registry.json")
    assert rb.extraction_block(loaded, "SRR1", "GCF_A").assembly.inputs == INPUTS


def test_set_inputs_is_a_noop_without_an_assembly(tmp_path):
    r = reg.load_registry(tmp_path / "prepare_registry.json")
    set_assembly_inputs(r, "SRR1", "GCF_A", INPUTS)
    assert r.datasets == {}
    reg.record_extraction(r, "SRR1", "GCF_A", [], 10, False, {})
    before = dict(r.datasets["SRR1"]["extractions"]["GCF_A"])
    set_assembly_inputs(r, "SRR1", "GCF_A", INPUTS)
    assert r.datasets["SRR1"]["extractions"]["GCF_A"] == before


def test_assembly_block_without_inputs_round_trips_byte_identical():
    data = {"date": "d", "dir": "x", "tool": "megahit", "params": {"preset": "sr"}, "contigs": 1, "N90": 5}
    before = json.dumps(data, separators=(",", ":"))
    assert json.dumps(rb.AssemblyBlock.from_dict(data).to_dict(), separators=(",", ":")) == before
    assert "inputs" not in rb.AssemblyBlock(date="d").to_dict()


def test_assembly_block_with_inputs_round_trips():
    data = {"date": "d", "dir": "x", "inputs": INPUTS, "contigs": 1}
    assert rb.AssemblyBlock.from_dict(data).to_dict() == data


# ------------------------------------------------------------------ assembly_predates_extraction


@pytest.mark.parametrize(
    "extraction_date, assembly_date, expected",
    [
        ("2026-10-01T10:00:00+00:00", "2026-10-01T09:00:00+00:00", True),
        ("2026-10-01T10:00:00+00:00", "2026-10-01T10:00:00+00:00", False),
        ("2026-10-01T10:00:00+00:00", "2026-10-01T11:00:00+00:00", False),
        # Same instant written with different UTC offsets (a daylight saving change): not earlier.
        ("2026-10-25T02:30:00+02:00", "2026-10-25T01:40:00+01:00", False),
        ("2026-10-25T02:30:00+02:00", "2026-10-25T01:20:00+01:00", True),
        (None, "2026-10-01T09:00:00+00:00", False),
        ("2026-10-01T10:00:00+00:00", None, False),
        ("", "2026-10-01T09:00:00+00:00", False),
        ("2026-10-01T10:00:00+00:00", "", False),
    ],
)
def test_assembly_predates_extraction(tmp_path, extraction_date, assembly_date, expected):
    r = _registry(tmp_path)
    _set_dates(r, extraction_date, assembly_date)
    assert assembly_predates_extraction(r, "SRR1", "GCF_A") is expected


def test_predates_is_false_without_blocks(tmp_path):
    r = reg.load_registry(tmp_path / "prepare_registry.json")
    assert assembly_predates_extraction(r, "SRR1", "GCF_A") is False
    reg.record_extraction(r, "SRR1", "GCF_A", [], 10, False, {})
    assert assembly_predates_extraction(r, "SRR1", "GCF_A") is False


def test_predates_compares_unparsable_dates_as_text(tmp_path):
    r = _registry(tmp_path)
    _set_dates(r, "b", "a")
    assert assembly_predates_extraction(r, "SRR1", "GCF_A") is True
    _set_dates(r, "a", "b")
    assert assembly_predates_extraction(r, "SRR1", "GCF_A") is False


# ------------------------------------------------------------------ legacy_assembly_current


@pytest.mark.parametrize(
    "recorded, preset, min_contig_len, expected",
    [
        ({"preset": "meta-sensitive", "min_contig_len": 500}, "meta-sensitive", 500, True),
        ({"preset": "meta-sensitive", "min_contig_len": 500}, "meta-large", 500, False),
        ({"preset": "meta-sensitive", "min_contig_len": 500}, "meta-sensitive", 1000, False),
        # A value missing on either side matches anything.
        ({}, "meta-sensitive", 1000, True),
        ({"preset": "meta-sensitive", "min_contig_len": 500}, None, None, True),
        ({"min_contig_len": 500}, "meta-large", 500, True),
        ({"preset": "meta-large"}, "meta-large", 2000, True),
        ({"preset": None, "min_contig_len": None}, "meta-large", 2000, True),
        # A recorded "default" is a preset like any other: it matches only "default" or a missing value.
        ({"preset": "default", "min_contig_len": 500}, "meta-large", 500, False),
        ({"preset": "default", "min_contig_len": 500}, "meta-sensitive", 500, False),
        ({"preset": "meta-large", "min_contig_len": 500}, "default", 500, False),
        ({"preset": "default", "min_contig_len": 500}, "default", 500, True),
        ({"preset": "default", "min_contig_len": 500}, None, 500, True),
    ],
)
def test_legacy_assembly_current_params(tmp_path, recorded, preset, min_contig_len, expected):
    r = _registry(tmp_path, assembly_params=recorded)
    _set_dates(r, "2026-10-01T10:00:00+00:00", "2026-10-01T11:00:00+00:00")
    assert legacy_assembly_current(r, "SRR1", "GCF_A", preset, min_contig_len) is expected


@pytest.mark.parametrize(
    "extraction_date, assembly_date, expected",
    [
        ("2026-10-01T10:00:00+00:00", "2026-10-01T11:00:00+00:00", True),
        ("2026-10-01T10:00:00+00:00", "2026-10-01T10:00:00+00:00", True),
        ("2026-10-01T10:00:00+00:00", "2026-10-01T09:00:00+00:00", False),
        (None, "2026-10-01T09:00:00+00:00", True),
        ("2026-10-01T10:00:00+00:00", None, True),
    ],
)
def test_legacy_assembly_current_dates(tmp_path, extraction_date, assembly_date, expected):
    r = _registry(tmp_path, assembly_params={"preset": "meta-large", "min_contig_len": 500})
    _set_dates(r, extraction_date, assembly_date)
    assert legacy_assembly_current(r, "SRR1", "GCF_A", "meta-large", 500) is expected


def test_legacy_assembly_current_needs_an_assembly(tmp_path):
    r = reg.load_registry(tmp_path / "prepare_registry.json")
    assert legacy_assembly_current(r, "SRR1", "GCF_A", None, None) is False
    reg.record_extraction(r, "SRR1", "GCF_A", [], 10, False, {})
    assert legacy_assembly_current(r, "SRR1", "GCF_A", None, None) is False


# ------------------------------------------------------------------ re-extraction and reconcile


def test_re_extraction_drops_the_recorded_assembly(tmp_path):
    r = _registry(tmp_path)
    set_assembly_inputs(r, "SRR1", "GCF_A", INPUTS)
    reg.record_extraction(r, "SRR1", "GCF_A", [], 12, False, {"preset": "sr"})
    assert rb.extraction_block(r, "SRR1", "GCF_A").assembly is None


def _project(tmp_path):
    return reg.ProjectPaths(
        fastq=tmp_path / "fastq",
        metadata=tmp_path / "metadata",
        genomes=tmp_path / "genomes",
        targeted=tmp_path / "targeted",
        matches=tmp_path / "matches",
    )


def _extraction_on_disk(paths, accession):
    paths.genomes.mkdir(parents=True, exist_ok=True)
    (paths.genomes / "GCF_1.fna").write_text(">c\nACGT\n")
    folder = paths.targeted / accession
    folder.mkdir(parents=True)
    (folder / "GCF_1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
    asm = folder / "GCF_1_assembly"
    asm.mkdir()
    (asm / "final.contigs.fa").write_text(">k len=8\nACGTACGT\n")


def test_inferred_extraction_then_assembly_still_records(tmp_path):
    """Reconcile records an untracked extraction and then its assembly; dropping the assembly on
    re-extraction must not lose the second record."""
    paths = _project(tmp_path)
    paths.fastq.mkdir()
    registry = reg.bootstrap_from_disk(paths, target_path=tmp_path / reg.REGISTRY_FILENAME)
    _extraction_on_disk(paths, "SRR1")
    report = reg.reconcile(registry, paths)
    assert report.untracked_extractions == [("SRR1", "GCF_1")]
    block = rb.extraction_block(registry, "SRR1", "GCF_1")
    assert block.mapped_reads == 1 and block.inferred is True
    assert block.assembly is not None and block.assembly.contigs == 1


def test_bootstrap_records_extraction_and_assembly(tmp_path):
    paths = _project(tmp_path)
    paths.fastq.mkdir()
    _extraction_on_disk(paths, "SRR1")
    registry = reg.bootstrap_from_disk(paths, target_path=tmp_path / reg.REGISTRY_FILENAME)
    assert rb.extraction_block(registry, "SRR1", "GCF_1").assembly.contigs == 1
