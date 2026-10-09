"""Tests for prepare/data/genome_taxonomy.py.

Covers GTDB taxonomy parsing, API-based enrichment, cache round-trips,
long-format annotation, taxonomy filtering, and aggregation.

Run: pytest tests/test_genome_taxonomy.py -v
"""

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock

import requests

from prepare.core.models import TaxonomyInfo
from prepare.core.exceptions import DataAccessError
from prepare.data.genome_taxonomy import (
    parse_gtdb_taxonomy_string,
    _lookup_genome_taxonomy_gtdb,
    _taxonomy_from_record,
    _extract_search_rows,
    enrich_genomes_with_taxonomy,
    load_taxonomy_cache,
    save_taxonomy_cache,
    annotate_containment_with_taxonomy,
    filter_by_taxonomy,
    summarize_by_taxonomy,
)

# ============================================================================
# parse_gtdb_taxonomy_string
# ============================================================================


class TestParseGtdbTaxonomyString:

    def test_valid_full_string(self):
        s = "d__Bacteria;p__Firmicutes;c__Bacilli;" "o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis"
        info = parse_gtdb_taxonomy_string(s, "GCF_000009045.1")
        assert info.genome_id == "GCF_000009045.1"
        assert info.phylum == "Firmicutes"
        assert info.class_name == "Bacilli"
        assert info.order == "Bacillales"
        assert info.family == "Bacillaceae"
        assert info.genus == "Bacillus"
        assert info.species == "Bacillus subtilis"

    def test_empty_string(self):
        info = parse_gtdb_taxonomy_string("", "GCF_000")
        assert info.genome_id == "GCF_000"
        assert info.species is None
        assert info.family is None

    def test_partial_string(self):
        info = parse_gtdb_taxonomy_string("d__Bacteria;p__Proteobacteria", "GCF_001")
        assert info.phylum == "Proteobacteria"
        assert info.family is None

    def test_malformed_tokens_ignored(self):
        info = parse_gtdb_taxonomy_string("bad_token;f__Fam1;also_bad", "g1")
        assert info.family == "Fam1"
        assert info.genus is None

    def test_empty_values_not_set(self):
        info = parse_gtdb_taxonomy_string("f__;g__Genus1", "g2")
        assert info.family is None
        assert info.genus == "Genus1"

    def test_whitespace_tokens(self):
        info = parse_gtdb_taxonomy_string(" f__Fam1 ; g__Gen1 ", "g3")
        assert info.family == "Fam1"
        assert info.genus == "Gen1"


# ============================================================================
# enrich_genomes_with_taxonomy (mocked API)
# ============================================================================


class TestEnrichGenomesWithTaxonomy:

    @patch("prepare.data.genome_taxonomy._lookup_genome_taxonomy_gtdb")
    def test_basic_enrichment(self, mock_lookup):
        mock_lookup.return_value = TaxonomyInfo(genome_id="GCF_001", family="Bacillaceae", genus="Bacillus")
        result = enrich_genomes_with_taxonomy(["GCF_001"])
        assert "GCF_001" in result
        assert result["GCF_001"].family == "Bacillaceae"
        mock_lookup.assert_called_once_with("GCF_001")

    @patch("prepare.data.genome_taxonomy._lookup_genome_taxonomy_gtdb")
    def test_uses_cache(self, mock_lookup, tmp_path):
        cache_file = tmp_path / "cache.tsv"
        tax = {"GCF_A": TaxonomyInfo(genome_id="GCF_A", family="FamA")}
        save_taxonomy_cache(tax, cache_file)

        result = enrich_genomes_with_taxonomy(["GCF_A"], cache_file=cache_file)
        assert result["GCF_A"].family == "FamA"
        mock_lookup.assert_not_called()

    @patch("prepare.data.genome_taxonomy._lookup_genome_taxonomy_gtdb")
    def test_api_failure_returns_empty_info(self, mock_lookup):
        mock_lookup.side_effect = DataAccessError("API down")
        result = enrich_genomes_with_taxonomy(["GCF_FAIL"])
        assert "GCF_FAIL" in result
        assert result["GCF_FAIL"].family is None

    @patch("prepare.data.genome_taxonomy._lookup_genome_taxonomy_gtdb")
    def test_none_result_gives_empty_info(self, mock_lookup):
        mock_lookup.return_value = None
        result = enrich_genomes_with_taxonomy(["GCF_NONE"])
        assert result["GCF_NONE"].family is None


# ============================================================================
# load_taxonomy_cache / save_taxonomy_cache
# ============================================================================


class TestTaxonomyCache:

    def test_roundtrip(self, tmp_path):
        cache_file = tmp_path / "tax.tsv"
        original = {
            "GCF_001": TaxonomyInfo(
                genome_id="GCF_001",
                species="Sp1",
                genus="Gen1",
                family="Fam1",
                order="Ord1",
                class_name="Cls1",
                phylum="Phy1",
                organism="Org1",
                tax_id="12345",
            ),
            "GCF_002": TaxonomyInfo(genome_id="GCF_002", family="Fam2"),
        }
        save_taxonomy_cache(original, cache_file)
        loaded = load_taxonomy_cache(cache_file)

        assert set(loaded.keys()) == {"GCF_001", "GCF_002"}
        assert loaded["GCF_001"].species == "Sp1"
        assert loaded["GCF_001"].tax_id == "12345"
        assert loaded["GCF_002"].family == "Fam2"
        assert loaded["GCF_002"].species is None

    def test_load_nonexistent_file(self, tmp_path):
        result = load_taxonomy_cache(tmp_path / "no_such_file.tsv")
        assert result == {}


# ============================================================================
# annotate_containment_with_taxonomy
# ============================================================================


class TestAnnotateContainment:

    def _make_df(self):
        return pd.DataFrame(
            {"genA": [0.5, 0.1], "genB": [0.8, 0.3]},
            index=["S1", "S2"],
        )

    def test_long_format_conversion(self):
        tax = {
            "genA": TaxonomyInfo(genome_id="genA", family="F1", genus="G1", species="Sp1"),
            "genB": TaxonomyInfo(genome_id="genB", family="F2", genus="G2", species="Sp2"),
        }
        result = annotate_containment_with_taxonomy(self._make_df(), tax)
        assert set(result.columns) == {"sample", "genome", "containment", "species", "genus", "family"}
        assert len(result) == 4  # 2 samples x 2 genomes

    def test_missing_taxonomy_gives_none(self):
        result = annotate_containment_with_taxonomy(self._make_df(), {})
        assert result.iloc[0]["family"] is None

    def test_zero_containment_rows_are_dropped(self):
        """Behaviour change from the previous row-by-row version: a zero (or missing) containment
        cell is excluded from the long-format table rather than kept as a zero-valued row."""
        df = pd.DataFrame({"genA": [0.5, 0.0], "genB": [0.0, 0.3]}, index=["S1", "S2"])

        result = annotate_containment_with_taxonomy(df, {})

        assert len(result) == 2
        assert set(zip(result["sample"], result["genome"])) == {("S1", "genA"), ("S2", "genB")}
        assert (result["containment"] > 0).all()

    def test_matches_the_old_row_by_row_output_when_there_are_no_zeros(self):
        """Pin: on a fixture without zeros, the melt-based table equals the previous
        iterrows-based implementation (kept here for comparison, not used by the module)."""

        def _old_annotate(containment_df, taxonomy):
            genome_cols = list(containment_df.columns)
            records = []
            for _, row in containment_df.iterrows():
                sample = row.name if isinstance(row.name, str) else str(row.name)
                for genome in genome_cols:
                    info = taxonomy.get(genome, TaxonomyInfo(genome_id=genome))
                    records.append(
                        {
                            "sample": sample,
                            "genome": genome,
                            "containment": float(row[genome]),
                            "species": info.species,
                            "genus": info.genus,
                            "family": info.family,
                        }
                    )
            return pd.DataFrame(records)

        df = self._make_df()
        tax = {
            "genA": TaxonomyInfo(genome_id="genA", family="F1", genus="G1", species="Sp1"),
            "genB": TaxonomyInfo(genome_id="genB", family="F2", genus="G2", species="Sp2"),
        }

        old = _old_annotate(df, tax)
        new = annotate_containment_with_taxonomy(df, tax)

        pd.testing.assert_frame_equal(old, new)

    def test_rows_are_sample_major_in_genome_column_order(self):
        """The 0.5.0 row order, written out by hand: all genomes of S2, then S1, in column order."""
        df = pd.DataFrame({"genB": [0.8, 0.3], "genA": [0.5, 0.1]}, index=["S2", "S1"])

        result = annotate_containment_with_taxonomy(df, {})

        assert list(zip(result["sample"], result["genome"], result["containment"])) == [
            ("S2", "genB", 0.8),
            ("S2", "genA", 0.5),
            ("S1", "genB", 0.3),
            ("S1", "genA", 0.1),
        ]
        assert list(result.index) == [0, 1, 2, 3]

    def test_order_is_kept_when_zero_cells_are_dropped(self):
        df = pd.DataFrame({"genA": [0.0, 0.2, 0.4], "genB": [0.6, 0.0, 0.7]}, index=["S1", "S2", "S3"])

        result = annotate_containment_with_taxonomy(df, {})

        assert list(zip(result["sample"], result["genome"])) == [
            ("S1", "genB"),
            ("S2", "genA"),
            ("S3", "genA"),
            ("S3", "genB"),
        ]
        assert list(result.index) == [0, 1, 2, 3]


# ============================================================================
# filter_by_taxonomy
# ============================================================================


class TestFilterByTaxonomy:

    @pytest.fixture()
    def annotated(self):
        return pd.DataFrame(
            {
                "sample": ["S1", "S1", "S2"],
                "genome": ["gA", "gB", "gA"],
                "containment": [0.5, 0.8, 0.1],
                "species": ["Sp1", "Sp2", "Sp1"],
                "genus": ["Gen1", "Gen2", "Gen1"],
                "family": ["Fam1", "Fam2", "Fam1"],
            }
        )

    def test_filter_by_family(self, annotated):
        result = filter_by_taxonomy(annotated, family="Fam1")
        assert len(result) == 2
        assert set(result["family"]) == {"Fam1"}

    def test_filter_by_genus(self, annotated):
        result = filter_by_taxonomy(annotated, genus="Gen2")
        assert len(result) == 1

    def test_filter_by_species(self, annotated):
        result = filter_by_taxonomy(annotated, species="sp1")
        assert len(result) == 2  # case-insensitive

    def test_filter_by_min_containment(self, annotated):
        result = filter_by_taxonomy(annotated, min_containment=0.5)
        assert len(result) == 2

    def test_combined_filters(self, annotated):
        result = filter_by_taxonomy(annotated, family="Fam1", min_containment=0.3)
        assert len(result) == 1
        assert result.iloc[0]["sample"] == "S1"


# ============================================================================
# summarize_by_taxonomy
# ============================================================================


class TestSummarizeByTaxonomy:

    @pytest.fixture()
    def annotated(self):
        return pd.DataFrame(
            {
                "sample": ["S1", "S1", "S2", "S2"],
                "genome": ["gA", "gB", "gA", "gB"],
                "containment": [0.5, 0.8, 0.1, 0.3],
                "family": ["Fam1", "Fam2", "Fam1", "Fam2"],
                "genus": ["G1", "G2", "G1", "G2"],
                "species": ["Sp1", "Sp2", "Sp1", "Sp2"],
            }
        )

    def test_summarize_family(self, annotated):
        result = summarize_by_taxonomy(annotated, level="family")
        assert "Fam1" in result.columns
        assert "Fam2" in result.columns
        assert result.loc["S1", "Fam1"] == 0.5

    def test_summarize_genus(self, annotated):
        result = summarize_by_taxonomy(annotated, level="genus")
        assert "G1" in result.columns

    def test_invalid_level_raises(self, annotated):
        with pytest.raises(ValueError, match="level must be"):
            summarize_by_taxonomy(annotated, level="kingdom")

    def test_empty_after_dropna(self):
        df = pd.DataFrame(
            {
                "sample": ["S1"],
                "genome": ["g1"],
                "containment": [0.5],
                "family": [None],
                "genus": [None],
                "species": [None],
            }
        )
        result = summarize_by_taxonomy(df, level="family")
        assert result.empty


# ============================================================================
# _lookup_genome_taxonomy_gtdb (mocked requests)
# ============================================================================

_GTDB_STRING = "d__Bacteria;p__Firmicutes;c__Bacilli;" "o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis"


def _resp(status_code, json_data):
    r = MagicMock()
    r.status_code = status_code
    r.json.return_value = json_data
    return r


class TestLookupGenomeTaxonomyGtdb:

    @patch("requests.Session.get")
    def test_genome_endpoint_success(self, mock_get):
        mock_get.return_value = _resp(
            200,
            {"gtdb_taxonomy": _GTDB_STRING, "organism_name": "Bacillus subtilis", "ncbi_taxid": 1423},
        )
        info = _lookup_genome_taxonomy_gtdb("GCF_000009045.1")
        assert info is not None
        assert info.species == "Bacillus subtilis"
        assert info.family == "Bacillaceae"
        assert info.organism == "Bacillus subtilis"
        assert info.tax_id == "1423"
        # genome endpoint hit returns immediately -- no fallback call
        assert mock_get.call_count == 1

    @patch("requests.Session.get")
    def test_genome_endpoint_empty_taxid_becomes_none(self, mock_get):
        mock_get.return_value = _resp(200, {"gtdb_taxonomy": _GTDB_STRING})
        info = _lookup_genome_taxonomy_gtdb("GCF_001")
        assert info is not None
        assert info.tax_id is None

    @patch("requests.Session.get")
    def test_falls_back_to_search_endpoint(self, mock_get):
        # Genome endpoint misses (404), search endpoint returns a result dict.
        mock_get.side_effect = [
            _resp(404, {}),
            _resp(200, {"results": [{"gtdb_taxonomy": _GTDB_STRING, "ncbi_taxid": 1423}]}),
        ]
        info = _lookup_genome_taxonomy_gtdb("GCF_002")
        assert info is not None
        assert info.genus == "Bacillus"
        assert info.tax_id == "1423"
        assert mock_get.call_count == 2

    @patch("requests.Session.get")
    def test_search_results_key_fallback(self, mock_get):
        mock_get.side_effect = [
            _resp(404, {}),
            _resp(200, {"results": [{"taxonomy": _GTDB_STRING}]}),
        ]
        info = _lookup_genome_taxonomy_gtdb("GCF_003")
        assert info is not None
        assert info.species == "Bacillus subtilis"

    @patch("requests.Session.get")
    def test_no_taxonomy_returns_none(self, mock_get):
        mock_get.side_effect = [
            _resp(404, {}),
            _resp(200, {"rows": []}),
        ]
        assert _lookup_genome_taxonomy_gtdb("GCF_004") is None

    @patch("requests.Session.get")
    def test_request_exception_raises_dataaccesserror(self, mock_get):
        mock_get.side_effect = requests.exceptions.ConnectionError("boom")
        with pytest.raises(DataAccessError, match="GTDB API error"):
            _lookup_genome_taxonomy_gtdb("GCF_005")

    @patch("requests.Session.get")
    def test_search_endpoint_bare_list(self, mock_get):
        # GTDB search may return a bare JSON list of records.
        mock_get.side_effect = [
            _resp(404, {}),
            _resp(200, [{"gtdb_taxonomy": _GTDB_STRING, "ncbi_taxid": 1423}]),
        ]
        info = _lookup_genome_taxonomy_gtdb("GCF_006")
        assert info is not None
        assert info.genus == "Bacillus"
        assert info.tax_id == "1423"

    @patch("requests.Session.get")
    def test_search_endpoint_rows_key(self, mock_get):
        # GTDB search may wrap records under a "rows" key.
        mock_get.side_effect = [
            _resp(404, {}),
            _resp(200, {"rows": [{"gtdb_taxonomy": _GTDB_STRING}]}),
        ]
        info = _lookup_genome_taxonomy_gtdb("GCF_007")
        assert info is not None
        assert info.species == "Bacillus subtilis"


class TestSearchResponseHelpers:

    def test_taxonomy_from_record_with_taxonomy(self):
        info = _taxonomy_from_record({"gtdb_taxonomy": _GTDB_STRING, "ncbi_taxid": 1423}, "GCF_X")
        assert info is not None
        assert info.species == "Bacillus subtilis"
        assert info.tax_id == "1423"

    def test_taxonomy_from_record_without_taxonomy_returns_none(self):
        assert _taxonomy_from_record({"organism_name": "Foo"}, "GCF_X") is None

    def test_extract_rows_from_list(self):
        assert _extract_search_rows([{"a": 1}]) == [{"a": 1}]

    def test_extract_rows_from_rows_key(self):
        assert _extract_search_rows({"rows": [{"a": 1}]}) == [{"a": 1}]

    def test_extract_rows_from_results_key(self):
        assert _extract_search_rows({"results": [{"b": 2}]}) == [{"b": 2}]

    def test_extract_rows_empty_dict(self):
        assert _extract_search_rows({}) == []

    def test_extract_rows_unexpected_type(self):
        assert _extract_search_rows("nonsense") == []


class TestCacheAndEnrichEdges:

    def test_load_cache_skips_rows_with_empty_genome_id(self, tmp_path):
        cache_file = tmp_path / "cache.tsv"
        header = "genome_id\tspecies\tgenus\tfamily\torder\tclass_name\tphylum\torganism\ttax_id\n"
        cache_file.write_text(header + "\tSp\t\t\t\t\t\t\t\n" + "GCF_OK\tSp2\t\t\t\t\t\t\t\n")
        cache = load_taxonomy_cache(cache_file)
        assert list(cache.keys()) == ["GCF_OK"]

    @patch("prepare.data.genome_taxonomy._lookup_genome_taxonomy_gtdb")
    def test_enrich_saves_cache_when_file_given(self, mock_lookup, tmp_path):
        mock_lookup.return_value = TaxonomyInfo(genome_id="GCF_NEW", family="Bacillaceae")
        cache_file = tmp_path / "out_cache.tsv"
        result = enrich_genomes_with_taxonomy(["GCF_NEW"], cache_file=cache_file)
        assert result["GCF_NEW"].family == "Bacillaceae"
        # missing genomes + cache_file => cache is written back
        assert cache_file.exists()
        assert "GCF_NEW" in cache_file.read_text()


class TestIncrementalCacheWrites:
    """The TSV cache is appended to after every successful lookup, not written once at the end."""

    @patch("prepare.data.genome_taxonomy._lookup_genome_taxonomy_gtdb")
    def test_a_later_failure_does_not_lose_earlier_successes(self, mock_lookup, tmp_path):
        cache_file = tmp_path / "cache.tsv"
        mock_lookup.side_effect = [
            TaxonomyInfo(genome_id="GCF_1", family="Fam1"),
            TaxonomyInfo(genome_id="GCF_2", family="Fam2"),
            DataAccessError("boom"),
        ]

        result = enrich_genomes_with_taxonomy(["GCF_1", "GCF_2", "GCF_3"], cache_file=cache_file)

        # the failed lookup still gets an empty in-memory placeholder for this run
        assert result["GCF_3"].family is None
        # but only the two successful lookups are persisted, so a later run retries GCF_3
        on_disk = load_taxonomy_cache(cache_file)
        assert set(on_disk.keys()) == {"GCF_1", "GCF_2"}

    @patch("prepare.data.genome_taxonomy._lookup_genome_taxonomy_gtdb")
    def test_a_genome_with_no_taxonomy_found_is_not_cached(self, mock_lookup, tmp_path):
        cache_file = tmp_path / "cache.tsv"
        mock_lookup.side_effect = [
            TaxonomyInfo(genome_id="GCF_1", family="Fam1"),
            None,
        ]

        enrich_genomes_with_taxonomy(["GCF_1", "GCF_2"], cache_file=cache_file)

        on_disk = load_taxonomy_cache(cache_file)
        assert set(on_disk.keys()) == {"GCF_1"}

    @patch("prepare.data.genome_taxonomy._lookup_genome_taxonomy_gtdb")
    def test_cache_file_grows_during_the_run_not_only_at_the_end(self, mock_lookup, tmp_path):
        cache_file = tmp_path / "cache.tsv"
        rows_seen_before_each_call = []

        def side_effect(accession):
            rows_seen_before_each_call.append(len(load_taxonomy_cache(cache_file)))
            return TaxonomyInfo(genome_id=accession, family=f"Fam_{accession}")

        mock_lookup.side_effect = side_effect

        enrich_genomes_with_taxonomy(["GCF_1", "GCF_2", "GCF_3"], cache_file=cache_file)

        assert rows_seen_before_each_call == [0, 1, 2]
