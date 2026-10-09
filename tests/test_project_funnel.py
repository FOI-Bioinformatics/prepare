"""Tests for `prepare.processing.project_funnel.funnel`."""

from prepare.data.registry import (
    Registry,
    record_assembly,
    record_download,
    record_exclusion,
    record_extraction,
    record_screening,
    record_selection,
    stage_members,
)
from prepare.data.registry_timing import set_assembly_timing, set_download_timing, set_extraction_timing
from prepare.processing.project_funnel import funnel


def _registry(tmp_path):
    """8 accessions moving through screening, selection, exclusion, download, extraction and assembly."""
    r = Registry(path=tmp_path / "prepare_registry.json")
    for n in range(1, 9):
        record_screening(r, f"SRR{n}", "G1", 0.1 * n, None, "matches", 0.0, None)
    record_selection(r, ["SRR1", "SRR2", "SRR3", "SRR5"], {"min_containment": 0.1}, tmp_path / "sel.txt")
    record_exclusion(r, "SRR3", "amplicon")
    for acc in ("SRR1", "SRR2", "SRR4"):
        record_download(r, acc, "downloaded", tmp_path / "fastq", attempt=False)
    record_download(r, "SRR6", "failed", tmp_path / "fastq")
    record_extraction(r, "SRR1", "G1", [], 50, False, {})
    record_extraction(r, "SRR2", "G1", [], 0, False, {})
    record_assembly(r, "SRR1", "G1", tmp_path / "asm", {"contigs": 4, "total_bp": 40000}, "v1", {})
    return r


def test_funnel_accession_counts_match_stage_members(tmp_path):
    r = _registry(tmp_path)
    members = stage_members(r)
    f = funnel(r)
    assert f["screened"]["accessions"] == len(members["screened"])
    assert f["selected"]["accessions"] == len(members["selected"])
    assert f["selected"]["excluded"] == len(members["excluded"])
    assert f["downloaded"]["accessions"] == len(members["downloaded"])
    assert f["analysed"]["accessions"] == len(members["analysed"])
    assert f["extracted"]["accessions"] == len(members["extracted"])
    assert f["assembled"]["accessions"] == len(members["assembled"])


def test_funnel_with_a_precomputed_members_dict_matches_funnel_without_one(tmp_path):
    r = _registry(tmp_path)
    members = stage_members(r)
    assert funnel(r, members) == funnel(r)


def test_funnel_downloaded_bytes_are_summed_from_bytes_total(tmp_path):
    root = tmp_path
    sizes = {}
    for acc in ("SRR1", "SRR2"):
        d = root / "fastq" / acc
        d.mkdir(parents=True)
        path = d / f"{acc}_1.fastq"
        path.write_text("@r\nACGT\n+\nIIII\n")
        sizes[acc] = path.stat().st_size
    r = Registry(path=root / "prepare_registry.json")
    record_download(r, "SRR1", "downloaded", root / "fastq")
    record_download(r, "SRR2", "downloaded", root / "fastq")
    f = funnel(r)
    assert f["downloaded"]["bytes"] == sum(sizes.values()) > 0


def test_funnel_counts_failed_downloads_separately_from_downloaded(tmp_path):
    r = Registry(path=tmp_path / "prepare_registry.json")
    record_download(r, "SRR1", "downloaded", tmp_path / "fastq", attempt=False)
    record_download(r, "SRR2", "failed", tmp_path / "fastq")
    f = funnel(r)
    assert f["downloaded"]["accessions"] == 1
    assert f["downloaded"]["failed"] == 1


def test_funnel_untimed_datasets_leave_seconds_none(tmp_path):
    r = _registry(tmp_path)
    f = funnel(r)
    assert f["downloaded"]["seconds"] is None
    assert f["downloaded"]["failed_seconds"] is None
    assert f["extracted"]["seconds"] is None
    assert f["assembled"]["seconds"] is None


def test_funnel_sums_seconds_only_over_timed_records(tmp_path):
    r = Registry(path=tmp_path / "prepare_registry.json")
    record_download(r, "SRR1", "downloaded", tmp_path / "fastq", attempt=False)
    record_download(r, "SRR2", "failed", tmp_path / "fastq")
    record_download(r, "SRR3", "downloaded", tmp_path / "fastq", attempt=False)  # left untimed
    set_download_timing(r, "SRR1", "2026-10-01T10:00:00+00:00", 10.0)
    set_download_timing(r, "SRR2", "2026-10-01T10:00:00+00:00", 5.0)
    record_extraction(r, "SRR1", "G1", [], 5, False, {})
    set_extraction_timing(r, "SRR1", "G1", "2026-10-01T10:00:00+00:00", 2.0)
    record_assembly(r, "SRR1", "G1", tmp_path / "asm", {"contigs": 2, "total_bp": 500}, "v1", {})
    set_assembly_timing(r, "SRR1", "G1", "2026-10-01T10:00:00+00:00", 3.0)
    f = funnel(r)
    # seconds covers the downloaded datasets only, like accessions and bytes; a failed
    # attempt's time is kept apart under failed_seconds.
    assert f["downloaded"]["seconds"] == 10.0
    assert f["downloaded"]["failed_seconds"] == 5.0
    assert f["extracted"]["seconds"] == 2.0
    assert f["assembled"]["seconds"] == 3.0
    assert f["assembled"]["total_bp"] == 500


def test_funnel_extraction_and_assembly_pairs_can_exceed_accessions(tmp_path):
    r = Registry(path=tmp_path / "prepare_registry.json")
    record_extraction(r, "SRR1", "G1", [], 10, False, {})
    record_extraction(r, "SRR1", "G2", [], 20, False, {})
    record_assembly(r, "SRR1", "G1", tmp_path / "asm1", {"contigs": 3, "total_bp": 100}, "v1", {})
    record_assembly(r, "SRR1", "G2", tmp_path / "asm2", {"contigs": 5, "total_bp": 200}, "v1", {})
    f = funnel(r)
    assert f["extracted"]["accessions"] == 1
    assert f["extracted"]["pairs"] == 2
    assert f["assembled"]["accessions"] == 1
    assert f["assembled"]["pairs"] == 2
    assert f["assembled"]["total_bp"] == 300


def test_funnel_extraction_with_zero_mapped_reads_does_not_count(tmp_path):
    r = Registry(path=tmp_path / "prepare_registry.json")
    record_extraction(r, "SRR1", "G1", [], 0, False, {})
    f = funnel(r)
    assert f["extracted"]["accessions"] == 0
    assert f["extracted"]["pairs"] == 0


def test_funnel_assembly_with_zero_contigs_does_not_count(tmp_path):
    r = Registry(path=tmp_path / "prepare_registry.json")
    record_extraction(r, "SRR1", "G1", [], 10, False, {})
    record_assembly(r, "SRR1", "G1", tmp_path / "asm", {"contigs": 0}, "v1", {})
    f = funnel(r)
    assert f["assembled"]["accessions"] == 0
    assert f["assembled"]["pairs"] == 0


def test_funnel_on_an_empty_registry_is_all_zero_and_untimed(tmp_path):
    r = Registry(path=tmp_path / "prepare_registry.json")
    f = funnel(r)
    assert f == {
        "screened": {"accessions": 0},
        "selected": {"accessions": 0, "excluded": 0},
        "downloaded": {"accessions": 0, "bytes": 0, "seconds": None, "failed": 0, "failed_seconds": None},
        "analysed": {"accessions": 0},
        "extracted": {"accessions": 0, "pairs": 0, "seconds": None},
        "assembled": {"accessions": 0, "pairs": 0, "total_bp": 0, "seconds": None},
    }
