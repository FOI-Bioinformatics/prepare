"""Tests for accession selection from containment and metadata tables."""

import pytest

from prepare.core.exceptions import ProcessingError
from prepare.processing.selection import select_accessions, select_accessions_ranked


@pytest.fixture
def tables(tmp_path):
    cont = tmp_path / "parsed_containment.txt"
    cont.write_text(
        "\tGCF_A\tGCF_B\tmax_containment\tmax_containment_annotation\n"
        "SRR1\t0.95\t0.10\t0.95\tGCF_A\n"
        "SRR2\t0.20\t0.80\t0.80\tGCF_B\n"
        "SRR3\t0.05\t0.02\t0.05\tGCF_A\n"
    )
    meta = tmp_path / "branchwater_metadata.txt"
    meta.write_text(
        "Run_ID\tSample_Scientific_Name\tgeo_loc_name_country_calc\n"
        "SRR1\tWolbachia pipientis\tAustria\n"
        "SRR2\tmosquito metagenome\tFrance\n"
        "SRR3\tmosquito metagenome\tFrance\n"
    )
    return cont, meta


def test_default_uses_max_containment(tables):
    cont, _ = tables
    assert select_accessions(cont, threshold=0.5) == ["SRR1", "SRR2"]


def test_genome_column_and_threshold(tables):
    cont, _ = tables
    assert select_accessions(cont, genome_id="GCF_B", threshold=0.5) == ["SRR2"]


def test_metadata_filter_is_case_insensitive(tables):
    cont, meta = tables
    result = select_accessions(
        cont, threshold=0.0, metadata_file=meta, metadata_column="geo_loc_name_country_calc", metadata_value="france"
    )
    assert result == ["SRR2", "SRR3"]


def test_unknown_genome_lists_columns(tables):
    cont, _ = tables
    with pytest.raises(ProcessingError, match="GCF_A"):
        select_accessions(cont, genome_id="GCF_Z")


def test_metadata_column_without_value_raises(tables):
    cont, meta = tables
    with pytest.raises(ProcessingError, match="metadata_value"):
        select_accessions(cont, metadata_file=meta, metadata_column="Sample_Scientific_Name")


def test_metadata_join_ignores_whitespace_in_run_ids(tmp_path):
    cont = tmp_path / "parsed_containment.txt"
    cont.write_text("\tGCF_A\tmax_containment\n" "SRR1 \t0.95\t0.95\n" "SRR2\t0.80\t0.80\n")
    meta = tmp_path / "branchwater_metadata.txt"
    meta.write_text("Run_ID\tgeo_loc_name_country_calc\n" " SRR1\tFrance\n" "SRR2\tFrance\n")
    result = select_accessions(
        cont, threshold=0.0, metadata_file=meta, metadata_column="geo_loc_name_country_calc", metadata_value="france"
    )
    assert result == ["SRR1", "SRR2"]


def test_top_n_applies_after_threshold(tables):
    cont, _ = tables
    assert select_accessions(cont, threshold=0.0, top_n=1) == ["SRR1"]


def test_top_n_applies_after_exclusions(tables):
    cont, _ = tables
    # SRR1 ranks first; excluding it should let SRR2 take the single top_n slot,
    # not truncate before exclusion is applied.
    result = select_accessions(cont, threshold=0.0, top_n=1, exclude={"SRR1"})
    assert result == ["SRR2"]


def test_top_n_applies_after_metadata_filter(tables):
    cont, meta = tables
    result = select_accessions(
        cont,
        threshold=0.0,
        top_n=1,
        metadata_file=meta,
        metadata_column="geo_loc_name_country_calc",
        metadata_value="france",
    )
    assert result == ["SRR2"]


def test_exclude_removes_accessions(tables):
    cont, _ = tables
    assert select_accessions(cont, threshold=0.0, exclude={"SRR2"}) == ["SRR1", "SRR3"]


def test_genome_id_and_genome_ids_are_mutually_exclusive(tables):
    cont, _ = tables
    with pytest.raises(ProcessingError, match="genome_id"):
        select_accessions(cont, genome_id="GCF_A", genome_ids=["GCF_A", "GCF_B"])


def test_genome_ids_any_ranks_on_row_wise_max(tables):
    cont, _ = tables
    # SRR1: max(0.95, 0.10)=0.95; SRR2: max(0.20, 0.80)=0.80; SRR3: max(0.05, 0.02)=0.05
    result = select_accessions(cont, genome_ids=["GCF_A", "GCF_B"], threshold=0.5, require="any")
    assert result == ["SRR1", "SRR2"]


def test_genome_ids_all_requires_every_column_at_threshold(tables):
    cont, _ = tables
    # SRR1 min(0.95, 0.10)=0.10, SRR2 min(0.20, 0.80)=0.20, SRR3 min(0.05, 0.02)=0.02 (excluded);
    # both remaining accessions rank on their min, best first.
    result = select_accessions(cont, genome_ids=["GCF_A", "GCF_B"], threshold=0.1, require="all")
    assert result == ["SRR2", "SRR1"]


def test_genome_ids_all_excludes_rows_below_threshold_on_any_column(tables):
    cont, _ = tables
    result = select_accessions(cont, genome_ids=["GCF_A", "GCF_B"], threshold=0.5, require="all")
    assert result == []


def test_require_invalid_value_raises(tables):
    cont, _ = tables
    with pytest.raises(ProcessingError, match="require"):
        select_accessions(cont, genome_ids=["GCF_A", "GCF_B"], threshold=0.1, require="both")


def test_top_n_zero_raises(tables):
    cont, _ = tables
    with pytest.raises(ProcessingError, match="--top-n must be a positive integer"):
        select_accessions(cont, threshold=0.0, top_n=0)


def test_top_n_negative_raises(tables):
    cont, _ = tables
    with pytest.raises(ProcessingError, match="--top-n must be a positive integer"):
        select_accessions(cont, threshold=0.0, top_n=-1)


def test_select_accessions_ranked_returns_accession_column_value(tables):
    cont, _ = tables
    ranked = select_accessions_ranked(cont, threshold=0.5)
    assert ranked == [("SRR1", "max_containment", 0.95), ("SRR2", "max_containment", 0.80)]


def test_select_accessions_ranked_respects_top_n_and_exclude(tables):
    cont, _ = tables
    ranked = select_accessions_ranked(cont, threshold=0.0, top_n=1, exclude={"SRR1"})
    assert ranked == [("SRR2", "max_containment", 0.80)]


def test_select_accessions_wraps_ranked(tables):
    cont, _ = tables
    assert select_accessions(cont, threshold=0.5) == [
        acc for acc, _, _ in select_accessions_ranked(cont, threshold=0.5)
    ]


# --- Run size, spot count and platform filters -------------------------------------------


@pytest.fixture
def run_tables(tmp_path):
    """A containment table and an NCBI-style metadata table with run size, spots and platform.

    SRR4 has an empty Run_Size; SRR5 is absent from the metadata table altogether.
    """
    cont = tmp_path / "parsed_containment.txt"
    cont.write_text(
        "\tGCF_A\tmax_containment\n"
        "SRR1\t0.95\t0.95\n"
        "SRR2\t0.80\t0.80\n"
        "SRR3\t0.05\t0.05\n"
        "SRR4\t0.60\t0.60\n"
        "SRR5\t0.30\t0.30\n"
    )
    meta = tmp_path / "metadata_table.txt"
    meta.write_text(
        "Run_ID\tRun_Size\tRun_Total_Spots\tPlatform\n"
        "SRR1\t2000000000\t10000000\tILLUMINA\n"
        "SRR2\t400000000\t2000000\tOXFORD_NANOPORE\n"
        "SRR3\t100000000\t500000\tILLUMINA\n"
        "SRR4\t\t3000000\tIllumina\n"
    )
    return cont, meta


@pytest.mark.parametrize(
    "text,expected",
    [
        ("1024", 1024),
        ("500M", 500_000_000),
        ("2G", 2_000_000_000),
        ("1.5G", 1_500_000_000),
        ("2GB", 2_000_000_000),
        ("750k", 750_000),
        (" 1T ", 1_000_000_000_000),
    ],
)
def test_parse_size_accepts_decimal_suffixes(text, expected):
    from prepare.processing.selection import parse_size

    assert parse_size(text) == expected


@pytest.mark.parametrize("text", ["abc", "", "G", "-1G", "1X", "1.5.2G", "nan"])
def test_parse_size_rejects_garbage(text):
    from prepare.processing.selection import parse_size

    with pytest.raises(ProcessingError):
        parse_size(text)


def test_run_filters_active():
    from prepare.processing.selection import RunFilters

    assert not RunFilters().active()
    assert RunFilters(min_spots=0).active()
    assert RunFilters(platform="ILLUMINA").active()


def test_max_run_size_drops_larger_and_unknown_runs(run_tables, caplog):
    from prepare.processing.selection import RunFilters

    cont, meta = run_tables
    with caplog.at_level("INFO"):
        result = select_accessions(
            cont, threshold=0.0, metadata_file=meta, run_filters=RunFilters(max_run_size=500_000_000)
        )
    # SRR1 is larger, SRR4 has no Run_Size value, SRR5 is absent from the metadata table.
    assert result == ["SRR2", "SRR3"]
    assert "3 accession(s) dropped by --max-run-size > 500000000 bytes (2 with no Run_Size value)" in caplog.text


def test_min_and_max_spots_bound_the_selection(run_tables, caplog):
    from prepare.processing.selection import RunFilters

    cont, meta = run_tables
    with caplog.at_level("INFO"):
        result = select_accessions(
            cont,
            threshold=0.0,
            metadata_file=meta,
            run_filters=RunFilters(min_spots=1_000_000, max_spots=5_000_000),
        )
    assert result == ["SRR2", "SRR4"]
    assert "dropped by --min-spots < 1000000 (1 with no Run_Total_Spots value)" in caplog.text
    assert "dropped by --max-spots > 5000000 (0 with no Run_Total_Spots value)" in caplog.text


def test_platform_filter_is_case_insensitive(run_tables, caplog):
    from prepare.processing.selection import RunFilters

    cont, meta = run_tables
    with caplog.at_level("INFO"):
        result = select_accessions(
            cont, threshold=0.0, metadata_file=meta, run_filters=RunFilters(platform=" illumina ")
        )
    assert result == ["SRR1", "SRR4", "SRR3"]
    assert "(1 with no Platform value)" in caplog.text


def test_run_filters_require_a_metadata_file(run_tables):
    from prepare.processing.selection import RunFilters

    cont, _ = run_tables
    with pytest.raises(ProcessingError, match="metadata_file"):
        select_accessions(cont, threshold=0.0, run_filters=RunFilters(max_run_size=1000))


def test_inactive_run_filters_do_not_require_a_metadata_file(run_tables):
    from prepare.processing.selection import RunFilters

    cont, _ = run_tables
    assert select_accessions(cont, threshold=0.5, run_filters=RunFilters()) == ["SRR1", "SRR2", "SRR4"]


@pytest.mark.parametrize(
    "filters,flag,column",
    [
        ({"max_run_size": 1000}, "--max-run-size", "Run_Size"),
        ({"min_spots": 5}, "--min-spots", "Run_Total_Spots"),
        ({"max_spots": 10}, "--max-spots", "Run_Total_Spots"),
        ({"platform": "ILLUMINA"}, "--platform", "Platform"),
    ],
)
def test_missing_run_column_is_an_error(tables, filters, flag, column):
    from prepare.processing.selection import RunFilters

    cont, meta = tables
    expected = (
        f"{flag} needs a {column} column; branchwater_metadata.txt has none. Run download_metadata for the "
        "candidate list, or use the NCBI metadata table"
    )
    with pytest.raises(ProcessingError) as excinfo:
        select_accessions(cont, threshold=0.0, metadata_file=meta, run_filters=RunFilters(**filters))
    assert str(excinfo.value) == expected


@pytest.mark.parametrize(
    "filters,flag,column",
    [
        ({"max_run_size": 1000}, "--max-run-size", "Run_Size"),
        ({"min_spots": 5}, "--min-spots", "Run_Total_Spots"),
        ({"platform": "ILLUMINA"}, "--platform", "Platform"),
    ],
)
def test_run_filter_warns_when_no_candidate_has_a_value(tables, caplog, filters, flag, column):
    from prepare.processing.selection import RunFilters

    cont, _ = tables
    meta = cont.parent / "metadata_table.txt"
    meta.write_text("Run_ID\tRun_Size\tRun_Total_Spots\tPlatform\nSRR1\t\t\t\nSRR2\t\t\t\n")
    with caplog.at_level("INFO"):
        result = select_accessions(cont, threshold=0.0, metadata_file=meta, run_filters=RunFilters(**filters))
    assert result == []
    warnings = [r.getMessage() for r in caplog.records if r.levelname == "WARNING"]
    assert warnings == [
        f"{flag} dropped all 3 remaining candidate(s); none has a {column} value. "
        "Run download_metadata for the candidate list first."
    ]
    assert f"3 accession(s) dropped by {flag}" in caplog.text


def test_run_filter_does_not_warn_when_some_candidates_have_a_value(run_tables, caplog):
    from prepare.processing.selection import RunFilters

    cont, meta = run_tables
    with caplog.at_level("INFO"):
        result = select_accessions(cont, threshold=0.0, metadata_file=meta, run_filters=RunFilters(max_run_size=1))
    assert result == []
    assert not [r for r in caplog.records if r.levelname == "WARNING"]


def test_run_filters_apply_before_top_n(run_tables):
    from prepare.processing.selection import RunFilters

    cont, meta = run_tables
    result = select_accessions(
        cont, threshold=0.0, top_n=1, metadata_file=meta, run_filters=RunFilters(max_run_size=500_000_000)
    )
    assert result == ["SRR2"]


def test_run_filters_apply_after_metadata_equality_filter(run_tables):
    from prepare.processing.selection import RunFilters

    cont, meta = run_tables
    result = select_accessions(
        cont,
        threshold=0.0,
        metadata_file=meta,
        metadata_column="Platform",
        metadata_value="illumina",
        run_filters=RunFilters(max_run_size=500_000_000),
    )
    assert result == ["SRR3"]


def test_min_spots_greater_than_max_spots_raises(run_tables):
    from prepare.processing.selection import RunFilters

    cont, meta = run_tables
    with pytest.raises(ProcessingError, match=r"--min-spots \(10\) is greater than --max-spots \(5\)"):
        select_accessions(cont, metadata_file=meta, run_filters=RunFilters(min_spots=10, max_spots=5))


def test_selected_volume_is_logged_after_top_n(run_tables, caplog):
    from prepare.processing.selection import RunFilters

    cont, meta = run_tables
    with caplog.at_level("INFO"):
        result = select_accessions(
            cont, threshold=0.0, top_n=2, metadata_file=meta, run_filters=RunFilters(platform="ILLUMINA")
        )
    assert result == ["SRR1", "SRR4"]
    assert "Selected volume: 2.00 GB across 2 run(s) (1 with unknown size)" in caplog.text


def test_repeated_run_id_rows_warn_and_use_the_first_row(tmp_path, caplog):
    from prepare.processing.selection import RunFilters

    cont = tmp_path / "parsed_containment.txt"
    cont.write_text("\tGCF_A\tmax_containment\nSRR1\t0.9\t0.9\nSRR2\t0.8\t0.8\n")
    meta = tmp_path / "metadata_table.txt"
    meta.write_text(
        "Run_ID\tRun_Size\tPlatform\n"
        "SRR1\t100\tILLUMINA\n"
        "SRR1\t900\tPACBIO_SMRT\n"
        " SRR2\t200\tILLUMINA\n"
        "SRR2\t300\tILLUMINA\n"
    )
    with caplog.at_level("INFO"):
        result = select_accessions(
            cont, threshold=0.0, metadata_file=meta, run_filters=RunFilters(max_run_size=150, platform="ILLUMINA")
        )
    # The first SRR1 row (100 bytes, ILLUMINA) is used, not the second (900 bytes, PACBIO_SMRT).
    assert result == ["SRR1"]
    warnings = [r.getMessage() for r in caplog.records if r.levelname == "WARNING"]
    assert warnings == ["2 repeated Run_ID row(s) in metadata_table.txt; the first row of each is used"]
    assert "Selected volume: 0.00 GB across 1 run(s) (0 with unknown size)" in caplog.text


def test_ties_resolve_by_run_name(tmp_path):
    cont = tmp_path / "parsed_containment.txt"
    cont.write_text(
        "\tGCF_A\tmax_containment\tmax_containment_annotation\n"
        "SRR9\t0.50\t0.50\tGCF_A\n"
        "SRR1\t0.50\t0.50\tGCF_A\n"
        "SRR5\t0.50\t0.50\tGCF_A\n"
        "SRR2\t0.90\t0.90\tGCF_A\n"
    )
    assert select_accessions(cont, threshold=0.1) == ["SRR2", "SRR1", "SRR5", "SRR9"]
