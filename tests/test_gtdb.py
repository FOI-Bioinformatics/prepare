"""Tests for the GTDB API client module."""

import logging
from unittest.mock import MagicMock, patch

import pytest
import requests

import prepare.data.gtdb as gtdb
from prepare.core.exceptions import DataAccessError
from prepare.data.gtdb import (
    get_accessions_for_genus,
    get_accessions_for_species,
    get_session,
    search_species,
    search_taxon,
)

# --- Realistic mock data ---

SPECIES_RESULT_LIST = [
    {
        "accession": "GCF_000005845.2",
        "species": "s__Escherichia coli",
        "isRep": True,
        "name": "Escherichia coli K-12",
    },
    {
        "accession": "GCF_000008865.2",
        "species": "s__Escherichia coli",
        "isRep": False,
        "name": "Escherichia coli O157:H7",
    },
]

SPECIES_RESULT_DICT = {
    "genomes": [
        {
            "gid": "GCF_000005845.2",
            "species": "s__Escherichia coli",
            "is_representative": True,
        },
    ]
}

TAXON_RESULT_LIST = [
    {
        "accession": "GCF_000005845.2",
        "species": "s__Escherichia coli",
        "isRep": True,
    },
    {
        "accession": "GCF_000240185.1",
        "species": "s__Escherichia fergusonii",
        "isRep": True,
    },
    {
        "accession": "GCF_001234567.1",
        "species": "s__Escherichia fergusonii",
        "isRep": False,
    },
]


# --- search_species ---


class TestSearchSpecies:
    @patch("requests.Session.get")
    def test_returns_list(self, mock_get):
        mock_response = MagicMock()
        mock_response.json.return_value = SPECIES_RESULT_LIST
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_species("Escherichia coli")

        assert len(result) == 2
        assert result[0]["accession"] == "GCF_000005845.2"

    @patch("requests.Session.get")
    def test_returns_dict_with_genomes(self, mock_get):
        mock_response = MagicMock()
        mock_response.json.return_value = SPECIES_RESULT_DICT
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_species("Escherichia coli")

        assert len(result) == 1
        assert result[0]["gid"] == "GCF_000005845.2"

    @patch("requests.Session.get")
    def test_returns_dict_without_genomes_key(self, mock_get):
        """A dict without 'genomes' wraps itself in a list."""
        single_record = {"accession": "GCF_000005845.2", "species": "E. coli"}
        mock_response = MagicMock()
        mock_response.json.return_value = single_record
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_species("Escherichia coli")

        assert len(result) == 1
        assert result[0]["accession"] == "GCF_000005845.2"

    @patch("requests.Session.get")
    def test_empty_results(self, mock_get):
        mock_response = MagicMock()
        mock_response.json.return_value = []
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_species("Nonexistent species")

        assert result == []

    @patch("requests.Session.get")
    def test_none_response(self, mock_get):
        mock_response = MagicMock()
        mock_response.json.return_value = None
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_species("Nonexistent species")

        assert result == []

    @patch("requests.Session.get")
    def test_api_error(self, mock_get):
        mock_get.side_effect = requests.exceptions.ConnectionError("Connection refused")

        with pytest.raises(DataAccessError, match="GTDB API error"):
            search_species("Escherichia coli")

    @patch("requests.Session.get")
    def test_http_error(self, mock_get):
        mock_response = MagicMock()
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError("500 Server Error")
        mock_get.return_value = mock_response

        with pytest.raises(DataAccessError, match="GTDB API error"):
            search_species("Escherichia coli")

    @patch("requests.Session.get")
    def test_http_error_400_hints_stale_name(self, mock_get):
        """A 400 usually means the name predates a GTDB genus reclassification, so say so."""
        mock_response = MagicMock()
        error_response = MagicMock(status_code=400)
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
            "400 Client Error", response=error_response
        )
        mock_get.return_value = mock_response

        with pytest.raises(DataAccessError, match="not a current GTDB name") as excinfo:
            search_species("Lactobacillus rhamnosus")

        assert "genus" in str(excinfo.value)

    @patch("requests.Session.get")
    def test_http_error_500_keeps_plain_message(self, mock_get):
        """Only a 400 gets the stale-name hint; other HTTP errors keep the existing message."""
        mock_response = MagicMock()
        error_response = MagicMock(status_code=500)
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
            "500 Server Error", response=error_response
        )
        mock_get.return_value = mock_response

        with pytest.raises(DataAccessError) as excinfo:
            search_species("Escherichia coli")

        assert "not a current GTDB name" not in str(excinfo.value)

    @patch("requests.Session.get")
    def test_timeout(self, mock_get):
        mock_get.side_effect = requests.exceptions.Timeout("Request timed out")

        with pytest.raises(DataAccessError, match="GTDB API error"):
            search_species("Escherichia coli")

    @patch("requests.Session.get")
    def test_unexpected_data_type(self, mock_get):
        """Non-list, non-dict response returns empty list."""
        mock_response = MagicMock()
        mock_response.json.return_value = "unexpected string"
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_species("Escherichia coli")

        assert result == []


# --- search_taxon ---


class TestSearchTaxon:
    @patch("requests.Session.get")
    def test_returns_list(self, mock_get):
        mock_response = MagicMock()
        mock_response.json.return_value = TAXON_RESULT_LIST
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_taxon("Escherichia")

        assert len(result) == 3

    @patch("requests.Session.get")
    def test_with_limit(self, mock_get):
        mock_response = MagicMock()
        mock_response.json.return_value = TAXON_RESULT_LIST[:1]
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_taxon("Escherichia", limit=1)

        mock_get.assert_called_once()
        call_kwargs = mock_get.call_args
        assert call_kwargs[1]["params"]["limit"] == 1
        assert len(result) == 1

    @patch("requests.Session.get")
    def test_returns_dict_with_results_key(self, mock_get):
        mock_response = MagicMock()
        mock_response.json.return_value = {"results": TAXON_RESULT_LIST}
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_taxon("Escherichia")

        assert len(result) == 3

    @patch("requests.Session.get")
    def test_returns_dict_without_results_key(self, mock_get):
        single_record = {"accession": "GCF_000005845.2", "name": "Escherichia"}
        mock_response = MagicMock()
        mock_response.json.return_value = single_record
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_taxon("Escherichia")

        assert len(result) == 1

    @patch("requests.Session.get")
    def test_empty_results(self, mock_get):
        mock_response = MagicMock()
        mock_response.json.return_value = []
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_taxon("Nonexistent")

        assert result == []

    @patch("requests.Session.get")
    def test_api_error(self, mock_get):
        mock_get.side_effect = requests.exceptions.ConnectionError("Failed")

        with pytest.raises(DataAccessError, match="GTDB API error"):
            search_taxon("Escherichia")

    @patch("requests.Session.get")
    def test_http_error_400_hints_stale_name(self, mock_get):
        mock_response = MagicMock()
        error_response = MagicMock(status_code=400)
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
            "400 Client Error", response=error_response
        )
        mock_get.return_value = mock_response

        with pytest.raises(DataAccessError, match="not a current GTDB name"):
            search_taxon("Lactobacillus")

    @patch("requests.Session.get")
    def test_unexpected_data_type(self, mock_get):
        mock_response = MagicMock()
        mock_response.json.return_value = 42
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_taxon("Escherichia")

        assert result == []


# --- get_accessions_for_species ---


class TestGetAccessionsForSpecies:
    @patch("prepare.data.gtdb.search_species")
    def test_representative_only(self, mock_search):
        mock_search.return_value = SPECIES_RESULT_LIST

        result = get_accessions_for_species("Escherichia coli", representative_only=True)

        assert result == ["GCF_000005845.2"]

    @patch("prepare.data.gtdb.search_species")
    def test_all_genomes(self, mock_search):
        mock_search.return_value = SPECIES_RESULT_LIST

        result = get_accessions_for_species("Escherichia coli", representative_only=False)

        assert "GCF_000005845.2" in result
        assert "GCF_000008865.2" in result
        assert len(result) == 2

    @patch("prepare.data.gtdb.search_species")
    def test_no_representative_falls_back_to_first(self, mock_search):
        """When representative_only but none flagged, returns the first accession."""
        mock_search.return_value = [
            {"accession": "GCF_111111111.1", "isRep": False},
            {"accession": "GCF_222222222.1", "isRep": False},
        ]

        result = get_accessions_for_species("Some species", representative_only=True)

        assert result == ["GCF_111111111.1"]

    @patch("prepare.data.gtdb.search_species")
    def test_empty_results(self, mock_search):
        mock_search.return_value = []

        result = get_accessions_for_species("Nonexistent species")

        assert result == []

    @patch("prepare.data.gtdb.search_species")
    def test_alternative_accession_keys(self, mock_search):
        """Tests gid and ncbi_accession fallback keys."""
        mock_search.return_value = [
            {"gid": "GCF_000005845.2", "isRep": True},
        ]

        result = get_accessions_for_species("E. coli", representative_only=True)

        assert result == ["GCF_000005845.2"]

    @patch("prepare.data.gtdb.search_species")
    def test_is_representative_key(self, mock_search):
        """Tests the is_representative key variant."""
        mock_search.return_value = [
            {"accession": "GCF_000005845.2", "is_representative": True},
        ]

        result = get_accessions_for_species("E. coli", representative_only=True)

        assert result == ["GCF_000005845.2"]

    @patch("prepare.data.gtdb.search_species")
    def test_records_without_accession(self, mock_search):
        mock_search.return_value = [
            {"species": "E. coli", "isRep": True},
        ]

        result = get_accessions_for_species("E. coli", representative_only=False)

        assert result == []


# --- get_accessions_for_genus ---


class TestGetAccessionsForGenus:
    @patch("prepare.data.gtdb.search_taxon")
    def test_representative_only(self, mock_search):
        mock_search.return_value = TAXON_RESULT_LIST

        result = get_accessions_for_genus("Escherichia", representative_only=True)

        # Should get one per species (only representatives)
        assert "GCF_000005845.2" in result
        assert "GCF_000240185.1" in result
        assert len(result) == 2

    @patch("prepare.data.gtdb.search_taxon")
    def test_deduplicates_species(self, mock_search):
        """Only one accession per unique species name."""
        mock_search.return_value = TAXON_RESULT_LIST

        result = get_accessions_for_genus("Escherichia", representative_only=True)

        # E. fergusonii appears twice but only the representative should be included
        assert result.count("GCF_000240185.1") == 1

    @patch("prepare.data.gtdb.search_taxon")
    def test_empty_results(self, mock_search):
        mock_search.return_value = []

        result = get_accessions_for_genus("Nonexistent")

        assert result == []

    @patch("prepare.data.gtdb.search_taxon")
    def test_no_representatives_falls_back(self, mock_search):
        """When no records are flagged as representative, collects all accessions."""
        mock_search.return_value = [
            {"accession": "GCF_111111111.1", "species": "sp1", "isRep": False},
            {"accession": "GCF_222222222.1", "species": "sp2", "isRep": False},
        ]

        result = get_accessions_for_genus("SomeGenus", representative_only=True)

        assert "GCF_111111111.1" in result
        assert "GCF_222222222.1" in result

    @patch("prepare.data.gtdb.search_taxon")
    def test_all_genomes(self, mock_search):
        """With representative_only=False, no filtering happens."""
        mock_search.return_value = [
            {"accession": "GCF_000005845.2", "species": "sp1", "isRep": True},
            {"accession": "GCF_000008865.2", "species": "sp1", "isRep": False},
        ]

        result = get_accessions_for_genus("Escherichia", representative_only=False)

        # First record included; second has same species so skipped by dedup
        assert result == ["GCF_000005845.2"]

    @patch("prepare.data.gtdb.search_taxon")
    def test_uses_name_key_for_species(self, mock_search):
        """Falls back to 'name' key when 'species' is missing."""
        mock_search.return_value = [
            {"accession": "GCF_000005845.2", "name": "sp1", "isRep": True},
            {"accession": "GCF_000240185.1", "name": "sp2", "isRep": True},
        ]

        result = get_accessions_for_genus("Escherichia", representative_only=True)

        assert len(result) == 2


TAXON_MATCHES_RESPONSE = {
    "matches": ["g__Wolbachia", "s__Wolbachia massiliensis", "s__Wolbachia pipientis", "s__Wolbachia pipientis_A"]
}

SPECIES_GENOMES_RESPONSE = {
    "name": "Wolbachia pipientis",
    "genomes": [
        {"accession": "GCA_000174095.1", "gtdb_species_rep": False, "ncbi_org_name": "Wolbachia sp."},
        {"accession": "GCA_021378375.1", "gtdb_species_rep": True, "ncbi_org_name": "Wolbachia pipientis"},
    ],
}


class TestLiveApiShapes:
    @patch("requests.Session.get")
    def test_search_taxon_matches_shape(self, mock_get):
        mock_response = MagicMock()
        mock_response.json.return_value = TAXON_MATCHES_RESPONSE
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = search_taxon("Wolbachia")

        assert [r["name"] for r in result] == TAXON_MATCHES_RESPONSE["matches"]

    @patch("prepare.data.gtdb.search_species")
    def test_species_rep_key_recognised(self, mock_search):
        mock_search.return_value = SPECIES_GENOMES_RESPONSE["genomes"]

        assert get_accessions_for_species("Wolbachia pipientis") == ["GCA_021378375.1"]

    @patch("prepare.data.gtdb.search_species")
    @patch("prepare.data.gtdb.search_taxon")
    def test_genus_fans_out_to_species(self, mock_taxon, mock_species):
        mock_taxon.return_value = [{"name": n} for n in TAXON_MATCHES_RESPONSE["matches"]]
        mock_species.side_effect = lambda name: [
            {"accession": f"GCA_{abs(hash(name)) % 10**9:09d}.1", "gtdb_species_rep": True}
        ]

        accessions = get_accessions_for_genus("Wolbachia")

        assert len(accessions) == 3
        assert mock_species.call_args_list[0][0][0] == "Wolbachia massiliensis"
        assert all(acc.startswith("GCA_") for acc in accessions)

    @patch("prepare.data.gtdb.search_species")
    @patch("prepare.data.gtdb.search_taxon")
    def test_genus_skips_failing_species(self, mock_taxon, mock_species, caplog):
        mock_taxon.return_value = [
            {"name": "s__Wolbachia massiliensis"},
            {"name": "s__Wolbachia pipientis"},
            {"name": "s__Wolbachia pipientis_A"},
        ]

        def side_effect(name):
            if name == "Wolbachia pipientis":
                raise DataAccessError("boom")
            return [{"accession": f"GCA_{abs(hash(name)) % 10**9:09d}.1", "gtdb_species_rep": True}]

        mock_species.side_effect = side_effect

        with caplog.at_level(logging.WARNING):
            accessions = get_accessions_for_genus("Wolbachia")

        assert len(accessions) == 2
        assert "Skipping species" in caplog.text

    @patch("prepare.data.gtdb.search_species")
    @patch("prepare.data.gtdb.search_taxon")
    def test_genus_raises_when_every_species_fails(self, mock_taxon, mock_species):
        mock_taxon.return_value = [
            {"name": "s__Wolbachia massiliensis"},
            {"name": "s__Wolbachia pipientis"},
            {"name": "s__Wolbachia pipientis_A"},
        ]
        mock_species.side_effect = DataAccessError("boom")

        with pytest.raises(DataAccessError, match="All 3 species lookups"):
            get_accessions_for_genus("Wolbachia")

    @patch("prepare.data.gtdb.search_species")
    @patch("prepare.data.gtdb.search_taxon")
    def test_genus_fanout_excludes_other_genera(self, mock_taxon, mock_species):
        """GTDB /taxon/search is a name-contains search; only the target genus's
        species should be resolved, not every genus whose species name contains it."""
        mock_taxon.return_value = [
            {"name": "s__Bacillus subtilis"},
            {"name": "s__Paenibacillus larvae"},
            {"name": "s__Lactobacillus casei"},
        ]
        mock_species.return_value = [{"accession": "GCF_000009045.1", "gtdb_species_rep": True}]

        accessions = get_accessions_for_genus("Bacillus")

        assert accessions == ["GCF_000009045.1"]
        mock_species.assert_called_once_with("Bacillus subtilis")

    @patch("prepare.data.gtdb.search_species")
    @patch("prepare.data.gtdb.search_taxon")
    def test_genus_dedupes_accessions_shared_across_species(self, mock_taxon, mock_species):
        """Two species resolving to the same accession should yield it only once."""
        mock_taxon.return_value = [
            {"name": "s__Bacillus subtilis"},
            {"name": "s__Bacillus subtilis_A"},
        ]
        mock_species.return_value = [{"accession": "GCF_000009045.1", "gtdb_species_rep": True}]

        accessions = get_accessions_for_genus("Bacillus")

        assert accessions == ["GCF_000009045.1"]


# --- get_session ---


class TestGetSession:
    def setup_method(self):
        gtdb._session = None

    def teardown_method(self):
        gtdb._session = None

    def test_returns_a_session(self):
        assert isinstance(get_session(), requests.Session)

    def test_returns_the_same_session_on_every_call(self):
        assert get_session() is get_session()

    def test_retry_adapter_is_mounted_on_both_schemes(self):
        session = get_session()
        for scheme in ("https://", "http://"):
            adapter = session.adapters[scheme]
            assert adapter.max_retries.total == 3
            assert adapter.max_retries.backoff_factor == 0.5
            assert set(adapter.max_retries.status_forcelist) == {429, 500, 502, 503, 504}
            assert "GET" in adapter.max_retries.allowed_methods

    @patch("requests.Session.get")
    def test_session_is_reused_across_lookups(self, mock_get):
        """search_species and search_taxon share the one process-wide session."""
        mock_response = MagicMock()
        mock_response.json.return_value = []
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        session_before = get_session()
        search_species("Escherichia coli")
        search_taxon("Escherichia")

        assert mock_get.call_count == 2
        assert get_session() is session_before


# --- final 5xx response after retries (raise_on_status=False) ---


class TestRetriesReturnTheFinalResponse:
    """With retries exhausted the last 5xx response reaches the caller, as it did before retries."""

    def setup_method(self):
        gtdb._session = None

    def teardown_method(self):
        gtdb._session = None

    def test_session_retry_does_not_raise_on_status(self):
        adapter = get_session().adapters["https://"]
        assert adapter.max_retries.raise_on_status is False

    def test_genome_endpoint_500_falls_back_to_search(self, monkeypatch):
        from prepare.data.genome_taxonomy import _lookup_genome_taxonomy_gtdb
        from tests.fake_http import serve

        body = (
            b'{"rows": [{"gtdbTaxonomy": "d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;'
            b'f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus crispatus"}]}'
        )

        def handler(path):
            return (500, b"{}") if path.startswith("/genome/") else (200, body)

        requested = serve(monkeypatch, handler)
        info = _lookup_genome_taxonomy_gtdb("GCF_000001.1")
        assert info is not None and info.species == "Lactobacillus crispatus"
        # One first attempt plus RETRY_TOTAL retries on the genome endpoint, then the search endpoint once.
        assert requested == ["/genome/GCF_000001.1"] * (gtdb.RETRY_TOTAL + 1) + ["/search/gtdb"]

    def test_species_search_500_keeps_the_http_error_message(self, monkeypatch):
        from tests.fake_http import serve

        serve(monkeypatch, lambda path: (500, b"{}"))
        with pytest.raises(DataAccessError, match="500 Server Error"):
            search_species("Escherichia coli")
