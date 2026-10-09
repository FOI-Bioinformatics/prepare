"""Tests for the consolidated per-accession, per-genome results table."""

import pandas as pd
import pytest

from prepare.data import registry as reg
from prepare.processing.results import RESULTS_COLUMNS, results_dataframe, results_rows, screened_pairs


def _registry(tmp_path):
    r = reg.load_registry(tmp_path / "prepare_registry.json")
    reg.record_screening(r, "SRR1", "GCF_A", 0.9, 0.99, "matches", 0.0, None)
    reg.record_screening(r, "SRR1", "GCF_B", 0.2, None, "matches", 0.0, None)
    reg.record_screening(r, "SRR2", "GCF_A", 0.5, None, "matches", 0.0, None)
    reg.record_selection(r, ["SRR1"], {"column": "GCF_A", "threshold": 0.5}, tmp_path / "accessions.txt")
    reg.record_exclusion(r, "SRR2", "16S amplicon")
    reg.record_download(r, "SRR1", "downloaded", tmp_path / "fastq", attempt=False)
    reg.record_metadata(r, "SRR1", tmp_path / "metadata" / "SRR1.xml", {"run_size": "4000", "run_total_spots": "1000"})
    coverage = {"breadth": 0.8, "mean_depth": 12.5, "coverage_tsv": tmp_path / "targeted" / "cov.tsv"}
    reg.record_extraction(r, "SRR1", "GCF_A", [], 250, False, {}, coverage=coverage)
    stats = {"contigs": 12, "total_bp": 90000, "n50": 8000, "genome_fraction_estimate": 0.7, "mapping_rate": 0.95}
    reg.record_assembly(r, "SRR1", "GCF_A", tmp_path / "targeted" / "SRR1" / "GCF_A_assembly", stats, "v1", {})
    return r


def _by_pair(rows):
    return {(row["accession"], row["genome_id"]): row for row in rows}


def test_columns_in_stated_order():
    assert RESULTS_COLUMNS == [
        "accession",
        "genome_id",
        "containment",
        "selected",
        "excluded",
        "exclusion_reason",
        "download_state",
        "run_total_spots",
        "run_size",
        "total_reads",
        "gc_percent",
        "quality_grade",
        "mapped_reads",
        "mapping_rate_to_reference",
        "breadth",
        "mean_depth",
        "contigs",
        "total_bp",
        "n50",
        "genome_fraction_estimate",
        "assembly_mapping_rate",
        "download_seconds",
        "extraction_seconds",
        "assembly_seconds",
        "download_verdict",
        "quality_source",
        "assembly_largest",
        "assembly_n90",
        "assembly_gc_percent",
        "assembly_contigs_ge_1kb",
        "assembly_mean_depth_estimate",
        "assembly_dir",
        "coverage_tsv",
    ]


def test_timing_columns_come_from_registry(tmp_path):
    from prepare.data.registry_timing import set_assembly_timing, set_download_timing, set_extraction_timing

    r = _registry(tmp_path)
    set_download_timing(r, "SRR1", "2026-10-01T10:00:00+00:00", 12.5)
    set_extraction_timing(r, "SRR1", "GCF_A", "2026-10-01T10:01:00+00:00", 3.25)
    set_assembly_timing(r, "SRR1", "GCF_A", "2026-10-01T10:02:00+00:00", 60.0)
    rows = _by_pair(results_rows(r))
    row = rows[("SRR1", "GCF_A")]
    assert (row["download_seconds"], row["extraction_seconds"], row["assembly_seconds"]) == (12.5, 3.25, 60.0)
    untimed = rows[("SRR2", "GCF_A")]
    assert (untimed["download_seconds"], untimed["extraction_seconds"], untimed["assembly_seconds"]) == (
        None,
        None,
        None,
    )


def test_one_row_per_pair_sorted_by_containment(tmp_path):
    rows = results_rows(_registry(tmp_path))
    assert [(row["accession"], row["genome_id"]) for row in rows] == [
        ("SRR1", "GCF_A"),
        ("SRR2", "GCF_A"),
        ("SRR1", "GCF_B"),
    ]
    assert all(list(row) == RESULTS_COLUMNS for row in rows)


def test_row_fields_come_from_registry(tmp_path):
    rows = _by_pair(results_rows(_registry(tmp_path)))
    row = rows[("SRR1", "GCF_A")]
    assert row["containment"] == 0.9
    assert row["selected"] is True and row["excluded"] is False and row["exclusion_reason"] is None
    assert row["download_state"] == "downloaded"
    assert row["run_total_spots"] == 1000 and row["run_size"] == 4000
    assert row["mapped_reads"] == 250
    assert row["mapping_rate_to_reference"] == pytest.approx(0.25)
    assert row["breadth"] == 0.8 and row["mean_depth"] == 12.5
    assert (row["contigs"], row["total_bp"], row["n50"]) == (12, 90000, 8000)
    assert row["genome_fraction_estimate"] == 0.7 and row["assembly_mapping_rate"] == 0.95
    excluded = rows[("SRR2", "GCF_A")]
    assert excluded["selected"] is False and excluded["excluded"] is True
    assert excluded["exclusion_reason"] == "16S amplicon"
    assert excluded["download_state"] is None and excluded["mapped_reads"] is None


def test_download_verdict_comes_from_the_registry(tmp_path):
    r = _registry(tmp_path)
    reg.set_download_verdict(r, "SRR1", {"method": "spots", "ratio": 1.0, "verdict": "complete"})
    reg.record_download(r, "SRR2", "downloaded", tmp_path / "fastq", attempt=False)
    reg.set_download_verdict(r, "SRR2", {"method": "spots", "ratio": 0.4, "verdict": "truncated"})
    reg.record_extraction(r, "SRR9", "GCF_C", [], 7, False, {})
    reg.record_download(r, "SRR9", "downloaded", tmp_path / "fastq", attempt=False)
    reg.set_download_verdict(r, "SRR9", {"verdict": "unverified"})
    rows = _by_pair(results_rows(r))
    assert rows[("SRR1", "GCF_A")]["download_verdict"] == "complete"
    assert rows[("SRR2", "GCF_A")]["download_verdict"] == "truncated"
    assert rows[("SRR9", "GCF_C")]["download_verdict"] == "unverified"
    assert rows[("SRR1", "GCF_B")]["download_verdict"] == "complete"


def test_download_verdict_is_none_without_a_recorded_verdict(tmp_path):
    rows = _by_pair(results_rows(_registry(tmp_path)))
    assert rows[("SRR1", "GCF_A")]["download_verdict"] is None
    assert rows[("SRR2", "GCF_A")]["download_verdict"] is None


@pytest.mark.parametrize("value", ["complete", ["complete"], 1])
def test_download_verdict_ignores_a_hand_edited_non_mapping_value(tmp_path, value):
    """A ``complete`` that is not a mapping reads as no verdict, as in ``status``'s download verdicts."""
    r = _registry(tmp_path)
    r.datasets["SRR1"]["download"]["complete"] = value
    rows = _by_pair(results_rows(r))
    assert rows[("SRR1", "GCF_A")]["download_verdict"] is None


def test_min_containment_drops_low_pairs(tmp_path):
    rows = results_rows(_registry(tmp_path), min_containment=0.5)
    assert {(row["accession"], row["genome_id"]) for row in rows} == {("SRR1", "GCF_A"), ("SRR2", "GCF_A")}


def test_genome_id_restricts_rows(tmp_path):
    rows = results_rows(_registry(tmp_path), genome_id="GCF_B")
    assert [(row["accession"], row["genome_id"]) for row in rows] == [("SRR1", "GCF_B")]


def test_mapping_rate_none_without_spots(tmp_path):
    r = reg.load_registry(tmp_path / "prepare_registry.json")
    reg.record_screening(r, "SRR3", "GCF_A", 0.4, None, "matches", 0.0, None)
    reg.record_extraction(r, "SRR3", "GCF_A", [], 50, False, {})
    row = results_rows(r)[0]
    assert row["mapped_reads"] == 50 and row["run_total_spots"] is None
    assert row["mapping_rate_to_reference"] is None


def test_mapping_rate_is_rounded_to_four_decimals(tmp_path):
    r = reg.load_registry(tmp_path / "prepare_registry.json")
    reg.record_metadata(r, "SRR3", tmp_path / "SRR3.xml", {"run_total_spots": 3})
    reg.record_extraction(r, "SRR3", "GCF_A", [], 1, False, {})
    assert results_rows(r)[0]["mapping_rate_to_reference"] == 0.3333


def test_mapping_rate_none_with_zero_mapped(tmp_path):
    r = reg.load_registry(tmp_path / "prepare_registry.json")
    reg.record_metadata(r, "SRR3", tmp_path / "SRR3.xml", {"run_total_spots": 100})
    reg.record_extraction(r, "SRR3", "GCF_A", [], 0, False, {})
    assert results_rows(r)[0]["mapping_rate_to_reference"] is None


def test_older_extraction_without_coverage_keys_gives_empty_values(tmp_path):
    r = reg.load_registry(tmp_path / "prepare_registry.json")
    r.datasets["SRR4"] = {"extractions": {"GCF_A": {"mapped_reads": 5, "files": []}}}
    row = results_rows(r)[0]
    assert row["breadth"] is None and row["mean_depth"] is None and row["contigs"] is None
    assert row["assembly_largest"] is None and row["assembly_n90"] is None
    assert row["assembly_gc_percent"] is None and row["assembly_contigs_ge_1kb"] is None
    assert row["assembly_mean_depth_estimate"] is None
    assert row["assembly_dir"] is None and row["coverage_tsv"] is None


def test_assembly_extra_stats_give_the_new_columns(tmp_path):
    r = _registry(tmp_path)
    stats = {
        "contigs": 12,
        "total_bp": 90000,
        "n50": 8000,
        "largest": 25000,
        "n90": 1500,
        "gc": 0.4123,
        "contigs_ge_1kb": 7,
        "mean_depth_estimate": 33.5,
    }
    reg.record_assembly(r, "SRR1", "GCF_A", tmp_path / "targeted" / "SRR1" / "GCF_A_assembly", stats, "v1", {})
    row = _by_pair(results_rows(r))[("SRR1", "GCF_A")]
    assert row["assembly_largest"] == 25000
    assert row["assembly_n90"] == 1500
    assert row["assembly_gc_percent"] == 41.23
    assert row["assembly_contigs_ge_1kb"] == 7
    assert row["assembly_mean_depth_estimate"] == 33.5
    # Project-relative: tmp_path is the registry's own folder, so "targeted/..." not an absolute path.
    assert row["assembly_dir"] == "targeted/SRR1/GCF_A_assembly"
    assert row["coverage_tsv"] == "targeted/cov.tsv"


def test_quality_source_is_none_without_any_quality_analysis(tmp_path):
    row = _by_pair(results_rows(_registry(tmp_path)))[("SRR1", "GCF_A")]
    assert row["quality_source"] is None


def test_parsed_table_adds_pairs_the_registry_capped(tmp_path):
    r = _registry(tmp_path)
    table = pd.DataFrame(
        {
            "GCF_A": [0.91234, 0.0, 0.3],
            "GCF_B": [0.2, 0.0, 0.0],
            "max_containment": [0.91234, 0.0, 0.3],
            "max_containment_annotation": ["GCF_A", "", "GCF_A"],
        },
        index=["SRR1", "SRR5", "SRR6"],
    )
    pairs = screened_pairs(r, table)
    assert pairs[("SRR1", "GCF_A")] == 0.91234
    assert pairs[("SRR6", "GCF_A")] == 0.3
    assert pairs[("SRR2", "GCF_A")] == 0.5
    assert ("SRR5", "GCF_A") not in pairs
    assert not any(genome.startswith("max_containment") for _, genome in pairs)
    rows = _by_pair(results_rows(r, parsed_table=table))
    assert rows[("SRR6", "GCF_A")]["containment"] == 0.3
    assert rows[("SRR1", "GCF_A")]["containment"] == 0.91234


@pytest.mark.parametrize("cell", [0.0, float("nan")])
def test_parsed_table_zero_or_nan_keeps_the_registry_containment(tmp_path, cell):
    r = _registry(tmp_path)
    table = pd.DataFrame(
        {"GCF_A": [cell, 0.95], "max_containment": [cell, 0.95], "max_containment_annotation": ["", "GCF_A"]},
        index=["SRR2", "SRR1"],
    )
    pairs = screened_pairs(r, table)
    assert pairs[("SRR2", "GCF_A")] == 0.5
    assert _by_pair(results_rows(r, parsed_table=table))[("SRR2", "GCF_A")]["containment"] == 0.5


def test_extraction_without_screening_appears_with_empty_containment(tmp_path):
    r = _registry(tmp_path)
    reg.record_extraction(r, "SRR9", "GCF_C", [], 7, False, {})
    rows = results_rows(r)
    assert (rows[-1]["accession"], rows[-1]["genome_id"]) == ("SRR9", "GCF_C")
    assert rows[-1]["containment"] is None and rows[-1]["mapped_reads"] == 7


def test_min_containment_above_zero_drops_pairs_with_unknown_containment(tmp_path):
    r = _registry(tmp_path)
    reg.record_extraction(r, "SRR9", "GCF_C", [], 7, False, {})
    assert ("SRR9", "GCF_C") not in _by_pair(results_rows(r, min_containment=0.1))
    assert ("SRR9", "GCF_C") in _by_pair(results_rows(r))


def test_dataframe_has_exact_columns_and_is_header_only_when_empty(tmp_path):
    df = results_dataframe(results_rows(_registry(tmp_path)))
    assert list(df.columns) == RESULTS_COLUMNS and len(df) == 3
    empty = results_dataframe([])
    assert list(empty.columns) == RESULTS_COLUMNS and empty.empty


def test_dataframe_keeps_integers_as_integers(tmp_path):
    df = results_dataframe(results_rows(_registry(tmp_path)))
    text = df.to_csv(sep="\t", index=False)
    line = next(line for line in text.splitlines() if line.startswith("SRR1\tGCF_A"))
    fields = dict(zip(RESULTS_COLUMNS, line.split("\t")))
    assert fields["run_total_spots"] == "1000" and fields["contigs"] == "12"
    missing = next(line for line in text.splitlines() if line.startswith("SRR2\tGCF_A"))
    assert dict(zip(RESULTS_COLUMNS, missing.split("\t")))["mapped_reads"] == ""


def test_profile_columns_come_from_the_profile_analysis(tmp_path):
    r = _registry(tmp_path)
    summary = {"total_reads": 2000, "gc_percent": 41.5, "quality_grade": "good"}
    reg.record_analysis(r, "SRR1", "profile", tmp_path / "profiles" / "SRR1_quality_profile.json", summary)
    row = _by_pair(results_rows(r))[("SRR1", "GCF_A")]
    assert (row["total_reads"], row["gc_percent"], row["quality_grade"]) == (2000, 41.5, "good")
    assert row["quality_source"] == "profile"


def test_profile_columns_read_an_old_registrys_quality_summary(tmp_path):
    """Before 0.5.0 sra_profile_quality recorded "quality" with GC as a 0-1 fraction."""
    r = _registry(tmp_path)
    old = {"grade": "fair", "total_reads": 1500, "reads_sampled": 100, "gc_content": 0.38}
    reg.record_analysis(r, "SRR1", "quality", tmp_path / "SRR1_quality_profile.json", old)
    row = _by_pair(results_rows(r))[("SRR1", "GCF_A")]
    assert row["total_reads"] == 1500
    assert row["gc_percent"] == pytest.approx(38.0)
    assert row["quality_grade"] == "fair"


def test_profile_columns_read_an_old_registrys_sra_stats_summary(tmp_path):
    """Before 0.5.0 sra_stats recorded "sra_stats" with GC already in percent."""
    r = _registry(tmp_path)
    reg.record_analysis(r, "SRR1", "sra_stats", tmp_path / "s.csv", {"total_reads": 900, "gc_content": 52.0})
    row = _by_pair(results_rows(r))[("SRR1", "GCF_A")]
    assert (row["total_reads"], row["gc_percent"], row["quality_grade"]) == (900, 52.0, None)
    assert row["quality_source"] == "legacy"


def _reference_row(r, accession, genome_id, containment):
    """One row as the earlier implementation built it, every block converted afresh for each pair."""
    from prepare.data import registry_blocks as rb

    exclusion = rb.exclusion_block(r, accession) or rb.ExclusionBlock()
    metadata = rb.metadata_block(r, accession) or rb.MetadataBlock()
    download = rb.download_block(r, accession)
    extraction = rb.extraction_block(r, accession, genome_id) or rb.ExtractionBlock()
    assembly = extraction.assembly.to_dict() if extraction.assembly is not None else {}
    spots = reg.to_int_or_none(metadata.run_total_spots)
    mapped = reg.to_int_or_none(extraction.mapped_reads)
    profile, quality_source = rb.quality_summary(r, accession)
    excluded = bool(exclusion.excluded)
    gc = assembly.get("gc")
    return {
        "accession": accession,
        "genome_id": genome_id,
        "containment": containment,
        "selected": bool((rb.selection_block(r, accession) or rb.SelectionBlock()).selected),
        "excluded": excluded,
        "exclusion_reason": (exclusion.reason or None) if excluded else None,
        "download_state": (download.state or None) if download is not None else None,
        "run_total_spots": spots,
        "run_size": reg.to_int_or_none(metadata.run_size),
        "total_reads": reg.to_int_or_none(profile["total_reads"]),
        "gc_percent": profile["gc_percent"],
        "quality_grade": profile["quality_grade"],
        "mapped_reads": mapped,
        "mapping_rate_to_reference": (
            round(mapped / spots, 4) if mapped and spots and mapped > 0 and spots > 0 else None
        ),
        "breadth": extraction.breadth,
        "mean_depth": extraction.mean_depth,
        "contigs": assembly.get("contigs"),
        "total_bp": assembly.get("total_bp"),
        "n50": assembly.get("n50"),
        "genome_fraction_estimate": assembly.get("genome_fraction_estimate"),
        "assembly_mapping_rate": assembly.get("mapping_rate"),
        "download_seconds": download.seconds if download is not None else None,
        "extraction_seconds": extraction.seconds,
        "assembly_seconds": assembly.get("seconds"),
        "download_verdict": (
            download.complete.verdict if download is not None and download.complete is not None else None
        ),
        "quality_source": quality_source,
        "assembly_largest": assembly.get("largest"),
        "assembly_n90": assembly.get("n90"),
        "assembly_gc_percent": (round(gc * 100, 2) if gc is not None else None),
        "assembly_contigs_ge_1kb": assembly.get("contigs_ge_1kb"),
        "assembly_mean_depth_estimate": assembly.get("mean_depth_estimate"),
        "assembly_dir": assembly.get("dir"),
        "coverage_tsv": extraction.coverage_tsv,
    }


def _sort_key(row):
    return (row["containment"] is None, -(row["containment"] or 0.0), row["accession"], row["genome_id"])


def test_rows_match_rows_built_afresh_for_each_pair(tmp_path):
    """Per-accession columns built once and shared across genomes give the same rows as before."""
    r = _registry(tmp_path)
    reg.record_screening(r, "SRR1", "GCF_C", 0.4, None, "matches", 0.0, None)
    reg.record_extraction(r, "SRR1", "GCF_B", [], 40, False, {})
    reg.record_analysis(r, "SRR1", "profile", tmp_path / "p.json", {"total_reads": 10, "gc_percent": 40.0})
    reg.record_extraction(r, "SRR9", "GCF_C", [], 7, False, {})
    # SRR3 appears twice: a positive value in the later row replaces the earlier one, as before.
    table = pd.DataFrame(
        {
            "GCF_A": [0.91234567, 0.0, float("nan"), 0.7],
            "GCF_D": ["0.3", "x", 0.2, 0.0],
            "max_containment": [0.9, 0.1, 0.2, 0.7],
        },
        index=["SRR1", "SRR2", "SRR3", "SRR3"],
    )
    for parsed in (None, table):
        pairs = screened_pairs(r, parsed)
        expected = sorted(
            (_reference_row(r, acc, genome, value) for (acc, genome), value in pairs.items()), key=_sort_key
        )
        assert results_rows(r, parsed_table=parsed) == expected
    pairs = screened_pairs(r, table)
    assert pairs[("SRR1", "GCF_A")] == 0.91234567
    assert pairs[("SRR1", "GCF_D")] == 0.3
    assert ("SRR2", "GCF_D") not in pairs
    assert pairs[("SRR3", "GCF_A")] == 0.7
    assert pairs[("SRR3", "GCF_D")] == 0.2


def test_screened_pairs_reads_a_repeated_genome_column_once(tmp_path):
    """A hand-built table with a repeated column name uses the first column of that name."""
    r = _registry(tmp_path)
    table = pd.DataFrame([[0.5, 0.9, 0.2]], index=["SRR1"], columns=["GCF_A", "GCF_A", "GCF_B"])
    pairs = screened_pairs(r, table)
    assert pairs[("SRR1", "GCF_A")] == 0.5
    assert pairs[("SRR1", "GCF_B")] == 0.2
