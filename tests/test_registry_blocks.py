"""Round-trip tests for the typed registry blocks in ``prepare.data.registry_blocks``.

The fixture ``tests/fixtures/registry_v0.4.0.json`` was written by the ``record_*`` writers as
they stood before the blocks existed. Every block must come back out of ``from_dict``/``to_dict``
exactly as it went in, including keys no writer knows about.
"""

import gzip
import json
import re
import shutil
from pathlib import Path

import pytest

from prepare.data import registry as R
from prepare.data import registry_blocks as B
from prepare.data.registry_blocks import (
    AnalysisEntry,
    AssemblyBlock,
    DownloadBlock,
    ExclusionBlock,
    ExportEntry,
    ExtractionBlock,
    FileEntry,
    MetadataBlock,
    ProjectBlock,
    ScreeningBlock,
    ScreeningEntry,
    SelectionBlock,
    StoreBlock,
    Verdict,
)

FIXTURE = Path(__file__).parent / "fixtures" / "registry_v0.4.0.json"

DATASET_BLOCKS = {
    "screening": ScreeningBlock,
    "selection": SelectionBlock,
    "exclusion": ExclusionBlock,
    "download": DownloadBlock,
    "metadata": MetadataBlock,
}


@pytest.fixture
def fixture_data():
    return json.loads(FIXTURE.read_text())


def test_fixture_covers_every_block(fixture_data):
    seen = set()
    for record in fixture_data["datasets"].values():
        seen.update(record)
    assert set(DATASET_BLOCKS) | {"extractions", "analyses"} <= seen
    assert len(fixture_data["datasets"]) == 5


def test_every_dataset_block_round_trips(fixture_data):
    for accession, record in fixture_data["datasets"].items():
        for key, cls in DATASET_BLOCKS.items():
            if key in record:
                assert cls.from_dict(record[key]).to_dict() == record[key], (accession, key)
        for genome_id, extraction in record.get("extractions", {}).items():
            assert ExtractionBlock.from_dict(extraction).to_dict() == extraction, (accession, genome_id)
        for name, analysis in record.get("analyses", {}).items():
            assert AnalysisEntry.from_dict(analysis).to_dict() == analysis, (accession, name)


def test_project_and_store_keep_unknown_keys(fixture_data):
    project = ProjectBlock.from_dict(fixture_data["project"])
    store = StoreBlock.from_dict(fixture_data["store"])
    assert project.extra == {"notes": "kept by round trip"}
    assert store.extra == {"comment": "kept by round trip"}
    assert project.to_dict() == fixture_data["project"]
    assert store.to_dict() == fixture_data["store"]


def test_nested_blocks_are_typed(fixture_data):
    first = fixture_data["datasets"]["SRR0000001"]
    download = DownloadBlock.from_dict(first["download"])
    assert isinstance(download.files[0], FileEntry)
    assert isinstance(download.complete, Verdict)
    assert download.complete.verdict == "complete"
    assert download.mate_reads == [2, 2]
    extraction = ExtractionBlock.from_dict(first["extractions"]["GCF_000001"])
    assert isinstance(extraction.assembly, AssemblyBlock)
    assert extraction.assembly.extra["genome_fraction_estimate"] == 0.8
    screening = ScreeningBlock.from_dict(first["screening"])
    assert isinstance(screening.genomes["GCF_000001"], ScreeningEntry)
    project = ProjectBlock.from_dict(fixture_data["project"])
    assert isinstance(project.exports["results_table"], ExportEntry)


def test_partial_blocks_keep_their_shape(fixture_data):
    """A block a writer only ever wrote in part must not gain keys on the way back out."""
    verdict_only = fixture_data["datasets"]["SRR0000003"]["download"]
    assert DownloadBlock.from_dict(verdict_only).to_dict() == {"attempts": 0, "complete": {"verdict": "unverified"}}
    assembly_only = fixture_data["datasets"]["SRR0000004"]["extractions"]["GCF_000002"]
    assert set(ExtractionBlock.from_dict(assembly_only).to_dict()) == {"assembly", "files", "inferred", "mapped_reads"}


def test_from_dict_tolerates_missing_keys():
    block = DownloadBlock.from_dict({})
    assert block.state == ""
    assert block.files == []
    assert block.complete is None
    assert block.to_dict() == {}
    assert ExtractionBlock.from_dict(None).to_dict() == {}


def test_null_nested_entries_stay_null():
    screening = {"date": "d", "genomes": {"G1": None}}
    assert ScreeningBlock.from_dict(screening).to_dict() == screening
    download = {"attempts": 1, "files": [None], "complete": None}
    assert DownloadBlock.from_dict(download).to_dict() == download
    assert ScreeningEntry.from_dict({}).containment is None


def test_assignment_after_load_writes_the_key():
    block = DownloadBlock.from_dict({"attempts": 0})
    block.message = ""
    assert block.to_dict() == {"attempts": 0, "message": ""}
    block.discard("message")
    assert block.to_dict() == {"attempts": 0}
    project = ProjectBlock.from_dict({"id": "p"})
    project.exports["t"] = ExportEntry(date="d", output="o", summary={})
    assert project.to_dict() == {"id": "p", "exports": {"t": {"date": "d", "output": "o", "summary": {}}}}


def test_fresh_blocks_omit_optional_keys_left_at_their_default():
    selection = SelectionBlock(selected=True, date="d", criteria={}, output="o")
    assert selection.to_dict() == {"selected": True, "date": "d", "criteria": {}, "output": "o"}
    store = StoreBlock(root="/r", linked=["A"])
    assert store.to_dict() == {"root": "/r", "mode": "symlink", "linked": ["A"]}
    assembly = AssemblyBlock(date="d", dir="x", tool="megahit", version="v1")
    assert "inputs" not in assembly.to_dict()
    assembly.inputs = {"reads": []}
    assert assembly.to_dict()["inputs"] == {"reads": []}


def test_set_mate_reads_starts_a_download_block_when_there_is_none():
    registry = R.Registry()
    B.set_mate_reads(registry, "SRR1", [3, 3], [["SRR1_1.fastq.gz", 1, 1.0]])
    assert registry.datasets["SRR1"]["download"] == {
        "attempts": 0,
        "mate_reads": [3, 3],
        "mate_reads_signature": [["SRR1_1.fastq.gz", 1, 1.0]],
    }
    B.set_mate_reads(registry, "SRR1", [4, 4], [])
    assert registry.datasets["SRR1"]["download"]["mate_reads"] == [4, 4]


def test_mark_inferred_leaves_an_unrecorded_block_alone():
    registry = R.Registry()
    B.mark_inferred(registry, "SRR1", "download")
    B.mark_inferred(registry, "SRR1", "extractions", genome_id="G1")
    assert registry.datasets == {}


# ------------------------------------------------------------ quality_summary


def _analyses(registry, accession, **named):
    """Set ``accession``'s "analyses" map directly, one ``name=(date, summary)`` pair per entry."""
    analyses = {name: {"date": date, "output": "", "summary": summary} for name, (date, summary) in named.items()}
    registry.datasets.setdefault(accession, {})["analyses"] = analyses


def test_quality_summary_only_report_recorded_fills_gc_and_grade():
    registry = R.Registry()
    _analyses(
        registry,
        "SRR1",
        report=("2026-02-01T00:00:00+00:00", {"total_reads": 500, "gc_percent": 41.5, "quality_grade": "good"}),
    )
    summary, source = B.quality_summary(registry, "SRR1")
    assert summary == {"total_reads": 500, "gc_percent": 41.5, "quality_grade": "good"}
    assert source == "report"
    assert B.profile_summary(registry, "SRR1") == summary


def test_quality_summary_newer_report_wins_and_fills_none_from_profile():
    registry = R.Registry()
    _analyses(
        registry,
        "SRR1",
        profile=("2026-01-01T00:00:00+00:00", {"total_reads": 100, "gc_percent": 30.0, "quality_grade": "fair"}),
        report=("2026-02-01T00:00:00+00:00", {"total_reads": None, "gc_percent": 45.0, "quality_grade": "good"}),
    )
    summary, source = B.quality_summary(registry, "SRR1")
    assert summary == {"total_reads": 100, "gc_percent": 45.0, "quality_grade": "good"}
    assert source == "report"


def test_quality_summary_equal_dates_favour_profile():
    registry = R.Registry()
    _analyses(
        registry,
        "SRR1",
        profile=("2026-03-01T00:00:00+00:00", {"total_reads": 200, "gc_percent": 50.0, "quality_grade": "excellent"}),
        report=("2026-03-01T00:00:00+00:00", {"total_reads": 999, "gc_percent": 10.0, "quality_grade": "poor"}),
    )
    summary, source = B.quality_summary(registry, "SRR1")
    assert summary == {"total_reads": 200, "gc_percent": 50.0, "quality_grade": "excellent"}
    assert source == "profile"


def test_quality_summary_unparsable_date_counts_as_oldest():
    registry = R.Registry()
    _analyses(
        registry,
        "SRR1",
        profile=("not-a-date", {"total_reads": 200, "gc_percent": 50.0, "quality_grade": "excellent"}),
        report=("2026-03-01T00:00:00+00:00", {"total_reads": 999, "gc_percent": 10.0, "quality_grade": "poor"}),
    )
    summary, source = B.quality_summary(registry, "SRR1")
    assert summary == {"total_reads": 999, "gc_percent": 10.0, "quality_grade": "poor"}
    assert source == "report"


def test_quality_summary_legacy_only_when_neither_profile_nor_report_exists():
    registry = R.Registry()
    R.record_analysis(registry, "SRR1", "sra_stats", "s.csv", {"total_reads": 900, "gc_content": 52.0})
    summary, source = B.quality_summary(registry, "SRR1")
    assert summary == {"total_reads": 900, "gc_percent": 52.0, "quality_grade": None}
    assert source == "legacy"


def test_quality_summary_fills_a_field_still_none_from_the_legacy_analyses():
    # A registry written before 0.5.0 (sra_stats with 100 reads) that sra_report later reported on:
    # a "report" without total_reads keeps its own GC and grade and takes the legacy read count.
    registry = R.Registry()
    _analyses(
        registry,
        "SRR1",
        sra_stats=("2024-05-01T00:00:00+00:00", {"total_reads": 100, "gc_content": 52.0}),
        report=("2026-02-01T00:00:00+00:00", {"gc_percent": 41.5, "quality_grade": "good"}),
    )
    summary, source = B.quality_summary(registry, "SRR1")
    assert summary == {"total_reads": 100, "gc_percent": 41.5, "quality_grade": "good"}
    assert source == "report"


def test_quality_summary_profile_and_report_values_take_precedence_over_legacy():
    registry = R.Registry()
    _analyses(
        registry,
        "SRR1",
        sra_stats=("2024-05-01T00:00:00+00:00", {"total_reads": 100, "gc_content": 52.0}),
        quality=("2024-05-01T00:00:00+00:00", {"grade": "poor", "gc_content": 0.52}),
        profile=("2026-01-01T00:00:00+00:00", {"total_reads": 300, "gc_percent": None, "quality_grade": None}),
        report=("2025-01-01T00:00:00+00:00", {"total_reads": 200, "gc_percent": 45.0}),
    )
    summary, source = B.quality_summary(registry, "SRR1")
    assert summary == {"total_reads": 300, "gc_percent": 45.0, "quality_grade": "poor"}
    assert source == "profile"


def test_quality_summary_none_when_nothing_was_ever_recorded():
    registry = R.Registry()
    summary, source = B.quality_summary(registry, "SRR1")
    assert summary == {"total_reads": None, "gc_percent": None, "quality_grade": None}
    assert source is None


def _mask(value):
    """Timestamps and compressed file sizes differ between runs; everything else must not."""
    if isinstance(value, dict):
        return {k: "<masked>" if k in ("bytes", "bytes_total") else _mask(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_mask(v) for v in value]
    if isinstance(value, str) and re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d[+-]\d\d:\d\d", value):
        return "<stamp>"
    return value


def test_writers_reproduce_the_fixture(tmp_path, monkeypatch, fixture_data):
    """The record_* writers, going through the typed blocks, write exactly the fixture's shapes and values."""
    monkeypatch.chdir(tmp_path)
    fastq = tmp_path / "fastq"
    for acc, mates in (("SRR0000001", ("_1", "_2")), ("SRR0000002", ("",)), ("SRR0000005", ("",))):
        (fastq / acc).mkdir(parents=True)
        for mate in mates:
            with gzip.open(fastq / acc / f"{acc}{mate}.fastq.gz", "wt") as handle:
                handle.write("@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n")
    registry = R.Registry(path=tmp_path / R.REGISTRY_FILENAME)
    R.record_genome(registry, "GCF_000001", "genomes/GCF_000001.fna", "genomes/manifest.csv")
    R.record_genome(registry, "GCF_000002", "genomes/GCF_000002.fna", "")

    R.record_screening(
        registry, "SRR0000001", "GCF_000001", 0.912345, 0.98765, "branchwater", 0.1, "matches/GCF_000001.csv"
    )
    R.record_screening(registry, "SRR0000001", "GCF_000002", 0.25, None, "matches", 0.0, None)
    R.record_screening(registry, "SRR0000002", "GCF_000001", 0.5, None, "matches", 0.0, "matches/GCF_000001.csv")
    B.mark_inferred(registry, "SRR0000002", "screening")
    R.record_selection(registry, ["SRR0000001", "SRR0000002"], {"min_containment": 0.2}, "selected.txt")
    ranked = [{"accession": "SRR0000001", "rank": 1, "column": "GCF_000001", "value": 0.9123}]
    R.record_selection(registry, ["SRR0000001"], {"min_containment": 0.3, "top": 1}, "selected.txt", ranked=ranked)
    R.record_download(
        registry,
        "SRR0000001",
        "downloaded",
        fastq,
        message="downloaded, complete (2 reads, 2 expected)",
        complete={"verdict": "complete", "reads_r1": 2, "expected_spots": 2, "ratio": 1.0},
        source="store",
        store_name="SRR0000001",
    )
    signature = [["SRR0000001_1.fastq.gz", 60, 1758780000.5], ["SRR0000001_2.fastq.gz", 60, 1758780000.5]]
    B.set_mate_reads(registry, "SRR0000001", [2, 2], signature)
    fields = {
        "run_size": "1234",
        "run_md5": "abc",
        "assay_type": "WGS",
        "organism": "vaginal metagenome",
        "collection_date": "2020-01-01",
        "library_layout": "PAIRED",
        "platform": "ILLUMINA",
        "library_strategy": "WGS",
        "run_total_spots": "2",
        "run_total_bases": 16,
    }
    R.record_metadata(registry, "SRR0000001", "metadata/SRR0000001_metadata.xml", fields)
    R.record_analysis(registry, "SRR0000001", "quality", "analysis/SRR0000001.json", {"mean_q": 38.5})
    targeted = tmp_path / "targeted" / "SRR0000001"
    params = {
        "genome_fasta": "genomes/GCF_000001.fna",
        "preset": "sr",
        "threshold": 0.9,
        "filter_flags": "-F 4",
        "min_mapq": 0,
        "index": "genomes/GCF_000001.mmi",
    }
    coverage = {"breadth": 0.75, "mean_depth": 3.2, "coverage_tsv": targeted / "GCF_000001_coverage.tsv"}
    files = [targeted / "GCF_000001_1.fastq.gz", targeted / "GCF_000001_2.fastq.gz"]
    R.record_extraction(registry, "SRR0000001", "GCF_000001", files, 40, False, params, 42, coverage)
    stats = {"contigs": 3, "total_bp": 1500, "n50": 600, "largest": 700, "n90": 300, "gc": 0.41}
    stats.update({"genome_fraction_estimate": 0.8, "mapping_rate": 0.95})
    assembly_params = {"min_contig_len": 200, "k_list": "21,29,39"}
    asm_dir = targeted / "GCF_000001_assembly"
    R.record_assembly(registry, "SRR0000001", "GCF_000001", asm_dir, stats, "MEGAHIT v1.2.9", assembly_params)

    R.record_exclusion(registry, "SRR0000002", "low quality", source="qc")
    R.record_download(registry, "SRR0000002", "failed", fastq, message="prefetch failed: timeout")

    R.set_download_verdict(registry, "SRR0000003", {"verdict": "unverified"})
    R.record_metadata(registry, "SRR0000003", "metadata/SRR0000003_metadata.xml", {})
    B.mark_inferred(registry, "SRR0000003", "metadata")

    asm_dir = tmp_path / "targeted/SRR0000004/GCF_000002_assembly"
    R.record_assembly(
        registry,
        "SRR0000004",
        "GCF_000002",
        asm_dir,
        {"contigs": 1, "total_bp": 300, "n50": 300, "largest": 300},
        "",
        {},
    )
    B.mark_inferred(registry, "SRR0000004", "extractions", genome_id="GCF_000002")
    R.record_exclusion(registry, "SRR0000004", "contaminated")
    R.clear_exclusion(registry, "SRR0000004")

    R.record_download(registry, "SRR0000005", "downloaded", fastq)
    B.mark_inferred(registry, "SRR0000005", "download", attempts=0)
    R.set_download_verdict(registry, "SRR0000005", {"method": "spots", "ratio": 0.5, "verdict": "truncated"})
    R.record_extraction(registry, "SRR0000005", "GCF_000002", [], 0, True, {"preset": "sr"})
    R.clear_assembly(registry, "SRR0000005", "GCF_000002")

    registry.project = {k: v for k, v in fixture_data["project"].items() if k != "exports"}
    R.record_export(registry, "results_table", "results/results.tsv", {"rows": 5})
    registry.store = dict(fixture_data["store"])

    assert _mask(registry.datasets) == _mask(fixture_data["datasets"])
    assert _mask(registry.genomes) == _mask(fixture_data["genomes"])
    assert _mask(registry.project) == _mask(fixture_data["project"])


def test_registry_accessors_and_setters_round_trip(tmp_path, fixture_data):
    """Reading every block through an accessor and writing it back through its setter leaves the file as it was."""
    target = tmp_path / R.REGISTRY_FILENAME
    shutil.copy(FIXTURE, target)
    registry = R.load_registry(target)
    for accession, record in list(registry.datasets.items()):
        for key, getter, setter in (
            ("screening", B.screening_block, B.set_screening_block),
            ("selection", B.selection_block, B.set_selection_block),
            ("exclusion", B.exclusion_block, B.set_exclusion_block),
            ("download", B.download_block, B.set_download_block),
            ("metadata", B.metadata_block, B.set_metadata_block),
        ):
            if key in record:
                setter(registry, accession, getter(registry, accession))
        for genome_id in record.get("extractions", {}):
            B.set_extraction_block(registry, accession, genome_id, B.extraction_block(registry, accession, genome_id))
    B.set_project_block(registry, B.project_block(registry))
    B.set_store_block(registry, B.store_block(registry))
    R.save_registry(registry)

    # Compared as text, so a value that changed type (1 against 1.0) is caught too. The fixture
    # is kept indented for reading; the registry is written compact, so the fixture is re-dumped
    # the same way (sorted keys, compact separators) before the comparison.
    def without_updated(text):
        return re.sub(r'"updated":"[^"]*"', '"updated":""', text)

    expected = json.dumps(json.loads(FIXTURE.read_text()), sort_keys=True, separators=(",", ":")) + "\n"
    assert without_updated(target.read_text()) == without_updated(expected)
