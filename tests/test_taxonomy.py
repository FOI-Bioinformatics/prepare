"""
Tests for taxonomy integration functions.
"""

import pytest
from unittest.mock import patch, MagicMock
import pandas as pd

from prepare.core.exceptions import ProcessingError
from prepare.data.taxonomy import (
    NCBITaxonomyClient,
    validate_taxonomic_assignments,
    create_taxonomic_summary,
    analyze_taxonomic_composition,
)


class TestNCBITaxonomyClient:
    """Test NCBI taxonomy client functionality."""

    def setup_method(self):
        """Set up test client."""
        self.client = NCBITaxonomyClient("test@example.com")

    def test_client_initialization(self):
        """Test client initialization."""
        assert self.client.email == "test@example.com"
        assert self.client.api_key is None
        assert self.client.request_delay > 0

    def test_client_with_api_key(self):
        """Test client with API key."""
        client = NCBITaxonomyClient("test@example.com", "api_key_123")
        assert client.api_key == "api_key_123"
        assert client.request_delay < 0.5  # Should be faster with API key

    def test_clean_species_name(self):
        """Test species name cleaning."""
        # Test various cleaning operations
        assert self.client._clean_species_name("candidate Bacteria") == "Bacteria"
        assert self.client._clean_species_name("uncultured bacterium") == "bacterium"
        assert self.client._clean_species_name("Escherichia coli strain K12") == "Escherichia coli"
        assert self.client._clean_species_name("Bacteria sp. ABC123") == "Bacteria"
        assert self.client._clean_species_name("Species [with brackets]") == "Species"
        assert self.client._clean_species_name("Species (with parentheses)") == "Species"

    @patch("requests.Session.get")
    def test_make_request_success(self, mock_get):
        """Test successful API request."""
        mock_response = MagicMock()
        mock_response.text = "<xml>test</xml>"
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        result = self.client._make_request("http://test.com", {"param": "value"})

        assert result == "<xml>test</xml>"
        mock_get.assert_called_once()

    @patch("requests.Session.get")
    def test_make_request_failure(self, mock_get):
        """Test failed API request."""
        mock_get.side_effect = Exception("Network error")

        with pytest.raises(Exception):
            self.client._make_request("http://test.com", {"param": "value"})

    @patch.object(NCBITaxonomyClient, "_make_request")
    def test_search_taxonomy(self, mock_request):
        """Test taxonomy search."""
        # Mock search response
        search_xml = """
        <eSearchResult>
            <IdList>
                <Id>12345</Id>
                <Id>67890</Id>
            </IdList>
        </eSearchResult>
        """

        # Mock fetch response
        fetch_xml = """
        <TaxaSet>
            <Taxon>
                <TaxId>12345</TaxId>
                <ScientificName>Escherichia coli</ScientificName>
                <Rank>species</Rank>
                <LineageEx>
                    <Taxon>
                        <ScientificName>Bacteria</ScientificName>
                        <Rank>superkingdom</Rank>
                    </Taxon>
                </LineageEx>
            </Taxon>
        </TaxaSet>
        """

        mock_request.side_effect = [search_xml, fetch_xml]

        results = self.client.search_taxonomy("Escherichia coli")

        assert len(results) == 1
        assert results[0]["tax_id"] == "12345"
        assert results[0]["scientific_name"] == "Escherichia coli"
        assert results[0]["rank"] == "species"

    def test_validate_species_name_mock(self):
        """Test species name validation with mocked data."""
        with patch.object(self.client, "search_taxonomy") as mock_search:
            mock_search.return_value = [
                {
                    "tax_id": "562",
                    "scientific_name": "Escherichia coli",
                    "rank": "species",
                    "lineage": "superkingdom:Bacteria;phylum:Proteobacteria",
                }
            ]

            result = self.client.validate_species_name("E. coli")

            assert result["is_valid"] is True
            assert result["tax_id"] == "562"
            assert result["validated_name"] == "Escherichia coli"
            assert "lineage" in result

    def test_validate_species_name_not_found(self):
        """Test species name validation when not found."""
        with patch.object(self.client, "search_taxonomy") as mock_search:
            mock_search.return_value = []

            result = self.client.validate_species_name("Nonexistent species")

            assert result["is_valid"] is False
            assert result["validated_name"] == ""
            assert result["tax_id"] == ""


class TestClientSession:
    """Each client owns one requests.Session, reused across its requests, with retries mounted."""

    def setup_method(self):
        self.client = NCBITaxonomyClient("test@example.com")

    def test_client_has_a_session(self):
        import requests

        assert isinstance(self.client.session, requests.Session)

    def test_session_is_reused_across_requests(self):
        with patch.object(self.client.session, "get") as mock_get:
            mock_response = MagicMock()
            mock_response.text = "<xml/>"
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            self.client._make_request("http://test.com", {"a": "1"})
            self.client._make_request("http://test.com", {"a": "2"})

            assert mock_get.call_count == 2

    def test_retry_adapter_is_mounted_on_both_schemes(self):
        for scheme in ("https://", "http://"):
            adapter = self.client.session.adapters[scheme]
            assert adapter.max_retries.total == 3
            assert adapter.max_retries.backoff_factor == 0.5
            assert set(adapter.max_retries.status_forcelist) == {429, 500, 502, 503, 504}
            assert "GET" in adapter.max_retries.allowed_methods

    def test_each_client_gets_its_own_session(self):
        other = NCBITaxonomyClient("other@example.com")
        assert other.session is not self.client.session


class TestTaxonomyValidation:
    """Test taxonomy validation functions."""

    def setup_method(self):
        """Set up test data."""
        self.species_list = ["Escherichia coli", "Bacillus subtilis", "Invalid species name"]

    @patch("prepare.data.taxonomy.NCBITaxonomyClient")
    def test_validate_taxonomic_assignments(self, mock_client_class):
        """Test batch taxonomy validation."""
        # Mock client
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client

        # Mock validation results
        mock_client.validate_species_name.side_effect = [
            {"original_name": "Escherichia coli", "is_valid": True, "tax_id": "562"},
            {"original_name": "Bacillus subtilis", "is_valid": True, "tax_id": "1423"},
            {"original_name": "Invalid species name", "is_valid": False, "tax_id": ""},
        ]

        result = validate_taxonomic_assignments(self.species_list, email="test@example.com")

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 3
        assert result["is_valid"].sum() == 2  # Two valid species

    def test_create_taxonomic_summary(self):
        """Test taxonomic summary creation."""
        # Create test data
        abundance_data = pd.DataFrame(
            {"Species A": [0.5, 0.3, 0.2], "Species B": [0.2, 0.4, 0.1], "Species C": [0.3, 0.3, 0.7]},
            index=["Sample1", "Sample2", "Sample3"],
        )

        taxonomy_data = pd.DataFrame(
            {
                "original_name": ["Species A", "Species B", "Species C"],
                "validated_name": ["Genus alpha sp1", "Genus beta sp1", "Genus alpha sp2"],
                "is_valid": [True, True, True],
                "lineage": [
                    "genus:Genus alpha;family:Family X",
                    "genus:Genus beta;family:Family Y",
                    "genus:Genus alpha;family:Family X",
                ],
            }
        )

        result = create_taxonomic_summary(abundance_data, taxonomy_data, level="genus")

        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == 3  # 3 samples
        assert "Genus alpha" in result.columns
        assert "Genus beta" in result.columns

        # Check that Species A and C are grouped under Genus alpha
        assert result.loc["Sample1", "Genus alpha"] == 0.8  # 0.5 + 0.3

    def test_analyze_taxonomic_composition(self):
        """Test multi-level taxonomic composition analysis."""
        # Create test data
        abundance_data = pd.DataFrame({"Species A": [1.0, 0.5], "Species B": [0.5, 1.0]}, index=["Sample1", "Sample2"])

        taxonomy_data = pd.DataFrame(
            {
                "original_name": ["Species A", "Species B"],
                "validated_name": ["Species A", "Species B"],
                "is_valid": [True, True],
                "lineage": ["phylum:Phylum1;genus:Genus1", "phylum:Phylum2;genus:Genus2"],
            }
        )

        result = analyze_taxonomic_composition(abundance_data, taxonomy_data, levels=["phylum", "genus"])

        assert isinstance(result, dict)
        assert "phylum" in result
        assert "genus" in result
        assert isinstance(result["phylum"], pd.DataFrame)
        assert isinstance(result["genus"], pd.DataFrame)


class TestTaxonomyEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_species_list(self):
        """Test with empty species list."""
        with patch("prepare.data.taxonomy.NCBITaxonomyClient"):
            result = validate_taxonomic_assignments([], email="test@example.com")

            assert isinstance(result, pd.DataFrame)
            assert len(result) == 0

    def test_create_summary_invalid_taxonomy(self):
        """Test taxonomic summary with invalid taxonomy data."""
        abundance_data = pd.DataFrame(
            {
                "Species A": [1.0],
            },
            index=["Sample1"],
        )

        taxonomy_data = pd.DataFrame(
            {
                "original_name": ["Species A"],
                "validated_name": ["Species A"],
                "is_valid": [False],  # Invalid taxonomy
                "lineage": [""],
            }
        )

        result = create_taxonomic_summary(abundance_data, taxonomy_data, level="genus")

        # Should handle invalid taxonomy gracefully
        assert isinstance(result, pd.DataFrame)
        assert result.shape[0] == 1  # 1 sample


class TestTaxonomyMapInput:
    def _abundance(self):
        return pd.DataFrame({"GCF_1": [0.9, 0.0], "GCF_2": [0.1, 0.7], "GCF_9": [0.0, 0.2]}, index=["S1", "S2"])

    def _taxonomy_map(self):
        return pd.DataFrame(
            {
                "genome_id": ["GCF_1", "GCF_2"],
                "species": ["Wolbachia pipientis", "Wolbachia sp947251865"],
                "genus": ["Wolbachia", "Wolbachia"],
                "family": ["Anaplasmataceae", "Anaplasmataceae"],
                "order": ["Rickettsiales", "Rickettsiales"],
                "class_name": ["Alphaproteobacteria", ""],
                "phylum": ["Pseudomonadota", "Pseudomonadota"],
                "organism": ["", ""],
                "tax_id": ["", ""],
            }
        )

    def test_genus_summary_from_map(self):
        result = create_taxonomic_summary(self._abundance(), self._taxonomy_map(), level="genus", min_abundance=0.0)
        assert result.loc["S1", "Wolbachia"] == pytest.approx(1.0)
        assert result.loc["S2", "Unclassified_genus"] == pytest.approx(0.2)

    def test_class_name_maps_to_class_rank(self):
        result = create_taxonomic_summary(self._abundance(), self._taxonomy_map(), level="class", min_abundance=0.0)
        assert result.loc["S1", "Alphaproteobacteria"] == pytest.approx(0.9)
        assert result.loc["S1", "Unclassified_class"] == pytest.approx(0.1)

    def test_unknown_table_shape_raises(self):
        with pytest.raises(ProcessingError, match="validate_taxonomy .* or enrich_taxonomy"):
            create_taxonomic_summary(self._abundance(), pd.DataFrame({"foo": [1]}), level="genus")


if __name__ == "__main__":
    pytest.main([__file__])


class TestNCBIRetriesReturnTheFinalResponse:
    """A persistent 5xx from NCBI ends as the same ProcessingError as before retries were added."""

    def test_retry_does_not_raise_on_status(self):
        client = NCBITaxonomyClient("test@example.com")
        assert client.session.adapters["https://"].max_retries.raise_on_status is False

    def test_persistent_500_is_a_processing_error_not_a_retry_error(self, monkeypatch):
        from tests.fake_http import serve

        requested = serve(monkeypatch, lambda path: (500, b"<error/>"))
        client = NCBITaxonomyClient("test@example.com")
        client.request_delay = 0
        with pytest.raises(ProcessingError, match="500 Server Error") as caught:
            client._make_request(client.base_url + "esearch.fcgi", {"db": "taxonomy"})
        assert "Max retries" not in str(caught.value)
        assert len(requested) == 4


class TestNCBIClientSessionIsClosed:
    """The client's HTTP session is closed once the client is done with."""

    def test_close_and_context_manager_close_the_session(self):
        client = NCBITaxonomyClient("test@example.com")
        with patch.object(client.session, "close") as closed:
            client.close()
        assert closed.call_count == 1
        spied = NCBITaxonomyClient("test@example.com")
        with patch.object(spied.session, "close") as closed:
            with spied:
                pass
        assert closed.call_count == 1

    @patch("prepare.data.taxonomy.NCBITaxonomyClient")
    def test_validate_taxonomic_assignments_closes_its_client(self, mock_client_class):
        client = mock_client_class.return_value
        client.validate_species_name.return_value = {
            "original_name": "E. coli",
            "validated_name": "Escherichia coli",
            "is_valid": True,
            "tax_id": "562",
            "rank": "species",
            "lineage": "",
            "confidence": "high",
        }
        validate_taxonomic_assignments(["E. coli"], "test@example.com")
        assert client.close.call_count == 1

    @patch("prepare.data.taxonomy.NCBITaxonomyClient")
    def test_client_is_closed_when_a_lookup_raises(self, mock_client_class):
        client = mock_client_class.return_value
        client.validate_species_name.side_effect = RuntimeError("boom")
        with pytest.raises(ProcessingError):
            validate_taxonomic_assignments(["E. coli"], "test@example.com")
        assert client.close.call_count == 1
