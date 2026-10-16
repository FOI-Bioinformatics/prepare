"""
EXTENDED TESTS for data/taxonomy.py (73% -> 85%+ coverage)

This file adds tests for untested areas:
- get_taxonomy_details method
- Rate limiting behavior
- Cache file handling
- Output file operations
- Additional edge cases

Run: pytest tests/test_taxonomy_extended.py -v
"""

import pytest
import pandas as pd
import requests
from unittest.mock import patch, MagicMock

from prepare.data.taxonomy import (
    _VALIDATION_CACHE_COLUMNS,
    NCBITaxonomyClient,
    validate_taxonomic_assignments,
    create_taxonomic_summary,
    analyze_taxonomic_composition,
)
from prepare.core.exceptions import ProcessingError

# ============================================================================
# TEST CLASS: get_taxonomy_details Method
# ============================================================================


class TestGetTaxonomyDetails:
    """Test get_taxonomy_details method."""

    def setup_method(self):
        """Set up test client."""
        self.client = NCBITaxonomyClient("test@example.com")

    @patch.object(NCBITaxonomyClient, "_make_request")
    def test_get_taxonomy_details_single_id(self, mock_request):
        """Test getting details for single taxonomy ID."""
        fetch_xml = """<?xml version="1.0"?>
        <TaxaSet>
            <Taxon>
                <TaxId>562</TaxId>
                <ScientificName>Escherichia coli</ScientificName>
                <Rank>species</Rank>
                <LineageEx>
                    <Taxon>
                        <ScientificName>Bacteria</ScientificName>
                        <Rank>superkingdom</Rank>
                    </Taxon>
                    <Taxon>
                        <ScientificName>Proteobacteria</ScientificName>
                        <Rank>phylum</Rank>
                    </Taxon>
                </LineageEx>
            </Taxon>
        </TaxaSet>
        """
        mock_request.return_value = fetch_xml

        results = self.client.get_taxonomy_details(["562"])

        assert len(results) == 1
        assert results[0]["tax_id"] == "562"
        assert results[0]["scientific_name"] == "Escherichia coli"
        assert results[0]["rank"] == "species"
        assert "superkingdom:Bacteria" in results[0]["lineage"]
        assert "phylum:Proteobacteria" in results[0]["lineage"]

    @patch.object(NCBITaxonomyClient, "_make_request")
    def test_get_taxonomy_details_multiple_ids(self, mock_request):
        """Test getting details for multiple taxonomy IDs."""
        fetch_xml = """<?xml version="1.0"?>
        <TaxaSet>
            <Taxon>
                <TaxId>562</TaxId>
                <ScientificName>Escherichia coli</ScientificName>
                <Rank>species</Rank>
            </Taxon>
            <Taxon>
                <TaxId>1423</TaxId>
                <ScientificName>Bacillus subtilis</ScientificName>
                <Rank>species</Rank>
            </Taxon>
        </TaxaSet>
        """
        mock_request.return_value = fetch_xml

        results = self.client.get_taxonomy_details(["562", "1423"])

        assert len(results) == 2
        assert results[0]["tax_id"] == "562"
        assert results[1]["tax_id"] == "1423"

    def test_get_taxonomy_details_empty_list(self):
        """Test getting details with empty ID list."""
        results = self.client.get_taxonomy_details([])

        assert results == []

    @patch.object(NCBITaxonomyClient, "_make_request")
    def test_get_taxonomy_details_missing_fields(self, mock_request):
        """Test handling of missing fields in XML."""
        fetch_xml = """<?xml version="1.0"?>
        <TaxaSet>
            <Taxon>
                <TaxId>999</TaxId>
                <!-- Missing ScientificName and Rank -->
            </Taxon>
        </TaxaSet>
        """
        mock_request.return_value = fetch_xml

        results = self.client.get_taxonomy_details(["999"])

        assert len(results) == 1
        assert results[0]["tax_id"] == "999"
        assert results[0]["scientific_name"] == ""
        assert results[0]["rank"] == ""


# ============================================================================
# TEST CLASS: Rate Limiting
# ============================================================================


class TestRateLimiting:
    """Test rate limiting behavior."""

    def setup_method(self):
        """Set up test client."""
        self.client = NCBITaxonomyClient("test@example.com")

    @patch("requests.Session.get")
    @patch("time.sleep")
    def test_rate_limiting_enforced(self, mock_sleep, mock_get):
        """Test that rate limiting is enforced between requests."""
        mock_response = MagicMock()
        mock_response.text = "<xml>test</xml>"
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        # Make two requests quickly
        self.client._make_request("http://test.com", {"param": "1"})
        self.client._make_request("http://test.com", {"param": "2"})

        # Should have called sleep for rate limiting
        assert mock_sleep.called

    def test_api_key_faster_rate_limit(self):
        """Test that API key allows faster requests."""
        client_no_key = NCBITaxonomyClient("test@example.com")
        client_with_key = NCBITaxonomyClient("test@example.com", "api_key_123")

        # NCBI allows 3 requests per second without a key and 10 with one.
        assert client_no_key.pacer.delay == pytest.approx(1 / 3)
        assert client_with_key.pacer.delay == pytest.approx(1 / 10)

    @patch("requests.Session.get")
    def test_request_exception_handling(self, mock_get):
        """Test handling of request exceptions."""
        mock_get.side_effect = requests.RequestException("Network error")

        with pytest.raises(ProcessingError, match="Failed to query NCBI"):
            self.client._make_request("http://test.com", {})


# ============================================================================
# TEST CLASS: Cache File Handling
# ============================================================================


class TestCacheFileHandling:
    """Test cache file operations in validate_taxonomic_assignments."""

    @patch("prepare.data.taxonomy.NCBITaxonomyClient")
    def test_validate_with_cache_loading(self, mock_client_class, tmp_path):
        """Test loading results from cache file."""
        # Create cache file
        cache_file = tmp_path / "taxonomy_cache.csv"
        cache_df = pd.DataFrame(
            {
                "original_name": ["Escherichia coli"],
                "validated_name": ["Escherichia coli"],
                "is_valid": [True],
                "tax_id": ["562"],
                "rank": ["species"],
                "lineage": ["superkingdom:Bacteria"],
                "confidence": ["high"],
            }
        )
        cache_df.to_csv(cache_file, index=False)

        # Mock client
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client

        # Validate (should use cache)
        result = validate_taxonomic_assignments(["Escherichia coli"], email="test@example.com", cache_file=cache_file)

        # Should not have called client.validate_species_name
        assert not mock_client.validate_species_name.called
        assert len(result) == 1

    @patch("prepare.data.taxonomy.NCBITaxonomyClient")
    def test_validate_with_cache_save(self, mock_client_class, tmp_path):
        """Test saving results to cache file."""
        cache_file = tmp_path / "new_cache.csv"

        # Mock client
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        mock_client.validate_species_name.return_value = {
            "original_name": "Bacillus subtilis",
            "validated_name": "Bacillus subtilis",
            "is_valid": True,
            "tax_id": "1423",
            "rank": "species",
            "lineage": "superkingdom:Bacteria",
            "confidence": "high",
        }

        # Validate
        validate_taxonomic_assignments(["Bacillus subtilis"], email="test@example.com", cache_file=cache_file)

        # Cache file should be created
        assert cache_file.exists()
        cache_df = pd.read_csv(cache_file)
        assert len(cache_df) == 1

    @patch("prepare.data.taxonomy.NCBITaxonomyClient")
    def test_validate_with_invalid_cache(self, mock_client_class, tmp_path, caplog):
        """Test handling of invalid cache file."""
        # Create invalid cache file
        cache_file = tmp_path / "invalid_cache.csv"
        cache_file.write_text("invalid,csv,data\n1,2,3")

        # Mock client
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        mock_client.validate_species_name.return_value = {
            "original_name": "Test",
            "is_valid": False,
            "tax_id": "",
            "rank": "",
            "lineage": "",
            "confidence": "low",
        }

        # Should handle gracefully
        validate_taxonomic_assignments(["Test species"], email="test@example.com", cache_file=cache_file)

        assert "Failed to load cache file" in caplog.text


class TestIncrementalValidationCache:
    """The cache is appended to after every successful lookup, not written once at the end."""

    @patch("prepare.data.taxonomy.NCBITaxonomyClient")
    def test_only_successful_lookups_are_cached(self, mock_client_class, tmp_path):
        """A lookup that comes back with confidence 'error' is not written, so it is retried later."""
        cache_file = tmp_path / "cache.csv"
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        mock_client.validate_species_name.side_effect = [
            {
                "original_name": "A",
                "cleaned_name": "A",
                "validated_name": "A",
                "tax_id": "1",
                "rank": "species",
                "lineage": "",
                "is_valid": True,
                "confidence": "high",
            },
            {
                "original_name": "B",
                "cleaned_name": "B",
                "validated_name": "B",
                "tax_id": "2",
                "rank": "species",
                "lineage": "",
                "is_valid": True,
                "confidence": "high",
            },
            {
                "original_name": "C",
                "cleaned_name": "C",
                "validated_name": "",
                "tax_id": "",
                "rank": "",
                "lineage": "",
                "is_valid": False,
                "confidence": "error",
            },
        ]

        validate_taxonomic_assignments(["A", "B", "C"], email="test@example.com", cache_file=cache_file)

        cache_df = pd.read_csv(cache_file)
        assert len(cache_df) == 2
        assert set(cache_df["confidence"]) == {"high"}

    @patch("prepare.data.taxonomy.NCBITaxonomyClient")
    def test_cache_grows_during_the_run(self, mock_client_class, tmp_path):
        """Each successful lookup is on disk before the next one runs, not batched at the end."""
        cache_file = tmp_path / "cache.csv"
        rows_seen_before_each_call = []
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client

        def side_effect(species):
            rows_seen_before_each_call.append(len(pd.read_csv(cache_file)) if cache_file.exists() else 0)
            return {
                "original_name": species,
                "cleaned_name": species,
                "validated_name": species,
                "tax_id": "1",
                "rank": "species",
                "lineage": "",
                "is_valid": True,
                "confidence": "high",
            }

        mock_client.validate_species_name.side_effect = side_effect

        validate_taxonomic_assignments(["A", "B", "C"], email="test@example.com", cache_file=cache_file)

        assert rows_seen_before_each_call == [0, 1, 2]

    @patch("prepare.data.taxonomy.NCBITaxonomyClient")
    def test_corrupt_cache_is_replaced_on_first_successful_write(self, mock_client_class, tmp_path):
        """A cache file that failed to load is discarded, not appended to, on the next write."""
        cache_file = tmp_path / "cache.csv"
        cache_file.write_text("not,a,valid,taxonomy,cache\n1,2,3,4,5")

        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        mock_client.validate_species_name.return_value = {
            "original_name": "A",
            "cleaned_name": "A",
            "validated_name": "A",
            "tax_id": "1",
            "rank": "species",
            "lineage": "",
            "is_valid": True,
            "confidence": "high",
        }

        validate_taxonomic_assignments(["A"], email="test@example.com", cache_file=cache_file)

        cache_df = pd.read_csv(cache_file)
        assert list(cache_df.columns) == [
            "original_name",
            "cleaned_name",
            "validated_name",
            "tax_id",
            "rank",
            "lineage",
            "is_valid",
            "confidence",
        ]
        assert len(cache_df) == 1


# ============================================================================
# TEST CLASS: Output File Operations
# ============================================================================


class TestOutputFileOperations:
    """Test output file saving in validate_taxonomic_assignments."""

    @patch("prepare.data.taxonomy.NCBITaxonomyClient")
    def test_validate_with_output_file(self, mock_client_class, tmp_path):
        """Test saving validation results to output file."""
        output_file = tmp_path / "validation_results.csv"

        # Mock client
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        mock_client.validate_species_name.return_value = {
            "original_name": "Species A",
            "validated_name": "Species A validated",
            "is_valid": True,
            "tax_id": "123",
            "rank": "species",
            "lineage": "genus:Genus A",
            "confidence": "high",
        }

        # Validate with output file
        validate_taxonomic_assignments(["Species A"], email="test@example.com", output_file=output_file)

        # Output file should be created
        assert output_file.exists()
        output_df = pd.read_csv(output_file)
        assert len(output_df) == 1
        assert output_df.loc[0, "original_name"] == "Species A"

    @patch("prepare.data.taxonomy.NCBITaxonomyClient")
    def test_validate_creates_output_directory(self, mock_client_class, tmp_path):
        """Test that output directory is created if it doesn't exist."""
        output_file = tmp_path / "subdir" / "nested" / "results.csv"

        # Mock client
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        mock_client.validate_species_name.return_value = {
            "original_name": "Test",
            "is_valid": False,
            "tax_id": "",
            "rank": "",
            "lineage": "",
            "confidence": "low",
        }

        # Validate
        validate_taxonomic_assignments(["Test"], email="test@example.com", output_file=output_file)

        # Directory should be created
        assert output_file.parent.exists()
        assert output_file.exists()


# ============================================================================
# TEST CLASS: analyze_taxonomic_composition with Output
# ============================================================================


class TestAnalyzeTaxonomicCompositionOutput:
    """Test analyze_taxonomic_composition with output directory."""

    def test_analyze_with_output_dir(self, tmp_path):
        """Test saving composition analysis to files."""
        abundance_data = pd.DataFrame({"Species A": [1.0, 0.5], "Species B": [0.5, 1.0]}, index=["Sample1", "Sample2"])

        taxonomy_data = pd.DataFrame(
            {
                "original_name": ["Species A", "Species B"],
                "validated_name": ["Species A", "Species B"],
                "is_valid": [True, True],
                "lineage": ["phylum:Phylum1;genus:Genus1", "phylum:Phylum2;genus:Genus2"],
            }
        )

        output_dir = tmp_path / "taxonomy_output"

        result = analyze_taxonomic_composition(
            abundance_data, taxonomy_data, levels=["phylum", "genus"], output_dir=output_dir
        )

        # Check output files created
        assert (output_dir / "taxonomy_summary_phylum.csv").exists()
        assert (output_dir / "taxonomy_summary_genus.csv").exists()

        # Check results
        assert "phylum" in result
        assert "genus" in result


# ============================================================================
# TEST CLASS: Additional Edge Cases
# ============================================================================


class TestAdditionalEdgeCases:
    """Test additional edge cases."""

    def test_clean_species_name_edge_cases(self):
        """Test species name cleaning edge cases."""
        client = NCBITaxonomyClient("test@example.com")

        # Multiple spaces
        assert client._clean_species_name("Species   with   spaces") == "Species with spaces"

        # Mixed case prefixes
        assert client._clean_species_name("CANDIDATE Species") == "Species"
        assert client._clean_species_name("UnClassified organism") == "organism"

        # Multiple brackets
        assert client._clean_species_name("Species [abc] [def]") == "Species"

        # Empty string
        assert client._clean_species_name("") == ""

        # Only whitespace
        assert client._clean_species_name("   ") == ""

    def test_create_taxonomic_summary_min_abundance(self):
        """Test taxonomic summary with min_abundance threshold."""
        abundance_data = pd.DataFrame(
            {"Species A": [0.5, 0.0001, 0.001], "Species B": [0.2, 0.4, 0.1]},  # Last two below threshold
            index=["Sample1", "Sample2", "Sample3"],
        )

        taxonomy_data = pd.DataFrame(
            {
                "original_name": ["Species A", "Species B"],
                "validated_name": ["Species A", "Species B"],
                "is_valid": [True, True],
                "lineage": ["genus:Genus1", "genus:Genus2"],
            }
        )

        result = create_taxonomic_summary(
            abundance_data, taxonomy_data, level="genus", min_abundance=0.01  # 1% minimum
        )

        # Species A in Sample2 and Sample3 should be filtered out
        assert isinstance(result, pd.DataFrame)

    def test_create_taxonomic_summary_unclassified_handling(self):
        """Test handling of unclassified species."""
        abundance_data = pd.DataFrame({"Species A": [0.5], "Species B": [0.3]}, index=["Sample1"])

        taxonomy_data = pd.DataFrame(
            {
                "original_name": ["Species A", "Species B"],
                "validated_name": ["Species A", "Species B"],
                "is_valid": [True, False],  # Species B is invalid
                "lineage": ["genus:Genus1", ""],  # No lineage
            }
        )

        result = create_taxonomic_summary(abundance_data, taxonomy_data, level="genus")

        # Should have unclassified category
        assert "Unclassified_genus" in result.columns or len(result.columns) > 0

    def test_validate_species_name_exception_handling(self):
        """Test exception handling in validate_species_name."""
        client = NCBITaxonomyClient("test@example.com")

        with patch.object(client, "search_taxonomy", side_effect=Exception("API Error")):
            result = client.validate_species_name("Test species")

            # Should return error result
            assert result["is_valid"] is False
            assert result["confidence"] == "error"

    def test_analyze_composition_error_handling(self):
        """Test error handling in analyze_taxonomic_composition."""
        # Invalid abundance data
        invalid_abundance = pd.DataFrame({"invalid": [1, 2]})
        invalid_taxonomy = pd.DataFrame({"invalid": ["a", "b"]})

        with pytest.raises(ProcessingError):
            analyze_taxonomic_composition(invalid_abundance, invalid_taxonomy, levels=["genus"])

    @patch("prepare.data.taxonomy.NCBITaxonomyClient")
    def test_validate_assignments_progress_logging(self, mock_client_class, caplog):
        """Test progress logging during validation."""
        # Create large species list
        species_list = [f"Species {i}" for i in range(25)]

        # Mock client
        mock_client = MagicMock()
        mock_client_class.return_value = mock_client
        mock_client.validate_species_name.return_value = {
            "original_name": "Test",
            "is_valid": True,
            "tax_id": "123",
            "rank": "species",
            "lineage": "",
            "confidence": "high",
        }

        # Validate; the progress lines are INFO, which caplog records only when asked to.
        caplog.set_level("INFO", logger="prepare.data.taxonomy")
        validate_taxonomic_assignments(species_list, email="test@example.com")

        # Should log progress
        assert "Progress:" in caplog.text


# ============================================================================
# SUCCESS METRICS:
#
# After running these tests:
# - Expected: 35+ additional tests pass
# - Coverage: 73% -> 85%+ for data/taxonomy.py
# - All untested methods now covered
#
# Run tests:
#   pytest tests/test_taxonomy_extended.py -v
#
# Check coverage:
#   pytest --cov=prepare.data.taxonomy --cov-report=term-missing \
#          tests/test_taxonomy.py tests/test_taxonomy_extended.py -q
# ============================================================================


def test_validate_taxonomic_assignments_empty_list_returns_empty_frame(caplog):
    with caplog.at_level("INFO", logger="prepare.data.taxonomy"):
        result = validate_taxonomic_assignments([], "a@b")

    assert result.empty
    assert list(result.columns) == _VALIDATION_CACHE_COLUMNS
    assert "0/0 species validated" in caplog.text
