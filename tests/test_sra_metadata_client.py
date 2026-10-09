"""
Tests for the SRA metadata client: NCBI queries and technology detection. Read statistics
are the shared statistics record (tests/test_store_stats.py) since 0.5.0.
"""

import pytest
from unittest.mock import patch, MagicMock
import sys


# Force fresh import of modules to avoid contamination from other tests
def _force_fresh_import():
    """Force fresh import of SRA modules to get real implementations."""
    modules_to_clear = [
        "prepare.data.sra_metadata",
    ]
    for module_name in modules_to_clear:
        if module_name in sys.modules:
            # Only remove if it's a Mock object (contaminated)
            if hasattr(sys.modules[module_name], "_mock_name"):
                del sys.modules[module_name]


# Apply fresh import at start of test session
_force_fresh_import()

from prepare.data.sra_metadata import (  # noqa: E402
    SRAMetadataClient,
    SRADatasetInfo,
    detect_sequencing_technology,
)


class TestSRAMetadataClient:
    """Test SRA metadata client functionality."""

    def setup_method(self):
        """Set up test client."""
        self.client = SRAMetadataClient("test@example.com")

    def test_client_initialization(self):
        """Test client initialization."""
        assert self.client.email == "test@example.com"
        assert self.client.api_key is None
        assert self.client.request_delay > 0

    def test_client_with_api_key(self):
        """Test client with API key."""
        client = SRAMetadataClient("test@example.com", "api_key_123")
        assert client.api_key == "api_key_123"
        assert client.request_delay < 0.5  # Should be faster with API key

    def test_make_request_success(self):
        """Test successful API request."""
        mock_response = MagicMock()
        mock_response.text = '{"test": "data"}'
        mock_response.raise_for_status.return_value = None
        with patch.object(self.client.session, "get", return_value=mock_response) as mock_get:
            result = self.client._make_request("http://test.com", {"param": "value"})

        assert result == '{"test": "data"}'
        mock_get.assert_called_once()

    def test_make_request_failure(self):
        """Test failed API request."""
        with patch.object(self.client.session, "get", side_effect=Exception("Network error")):
            with pytest.raises(Exception):
                self.client._make_request("http://test.com", {"param": "value"})

    def test_dataset_info_creation(self):
        """Test SRADatasetInfo creation."""
        info = SRADatasetInfo(
            accession="SRR123456",
            title="Test dataset",
            organism="Escherichia coli",
            platform="ILLUMINA",
            instrument="Illumina HiSeq 2500",
            strategy="WGS",
            layout="PAIRED",
            spots=1000000,
            bases=150000000,
            avg_length=150.0,
            size_mb=100.0,
            release_date="2023-01-01",
            bioproject="PRJNA123456",
            biosample="SAMN123456",
            library_selection="RANDOM",
            library_source="GENOMIC",
        )

        assert info.accession == "SRR123456"
        assert info.platform == "ILLUMINA"
        assert info.spots == 1000000


class TestSRAMetadataClientClose:
    def test_close_closes_the_session(self):
        client = SRAMetadataClient("a@b.c")
        with patch.object(client.session, "close") as close:
            client.close()
        close.assert_called_once()

    def test_context_manager_closes_the_session(self):
        with SRAMetadataClient("a@b.c") as client:
            session_close = patch.object(client.session, "close").start()
        session_close.assert_called_once()
        patch.stopall()


class TestTechnologyDetection:
    """Test sequencing technology detection."""

    def create_test_dataset_info(self, platform, instrument, strategy="WGS"):
        """Create test dataset info."""
        return SRADatasetInfo(
            accession="TEST123",
            title="Test",
            organism="Test organism",
            platform=platform,
            instrument=instrument,
            strategy=strategy,
            layout="PAIRED",
            spots=1000,
            bases=150000,
            avg_length=150.0,
            size_mb=100.0,
            release_date="2023-01-01",
            bioproject="",
            biosample="",
            library_selection="",
            library_source="",
        )

    def test_illumina_detection(self):
        """Test Illumina technology detection."""
        # Test platform detection
        info = self.create_test_dataset_info("ILLUMINA", "HiSeq 2500")
        assert detect_sequencing_technology(info) == "illumina"

        # Test instrument detection
        info = self.create_test_dataset_info("", "Illumina NovaSeq 6000")
        assert detect_sequencing_technology(info) == "illumina"

    def test_nanopore_detection(self):
        """Test Nanopore technology detection."""
        # Test platform detection
        info = self.create_test_dataset_info("OXFORD_NANOPORE", "MinION")
        assert detect_sequencing_technology(info) == "nanopore"

        # Test instrument detection
        info = self.create_test_dataset_info("", "GridION X5")
        assert detect_sequencing_technology(info) == "nanopore"

        # Test strategy detection
        info = self.create_test_dataset_info("", "", "NANOPORE")
        assert detect_sequencing_technology(info) == "nanopore"

    def test_pacbio_detection(self):
        """Test PacBio technology detection."""
        # Test platform detection
        info = self.create_test_dataset_info("PACBIO_SMRT", "PacBio RS")
        assert detect_sequencing_technology(info) == "pacbio"

        # Test instrument detection
        info = self.create_test_dataset_info("", "Sequel II")
        assert detect_sequencing_technology(info) == "pacbio"

    def test_unknown_technology(self):
        """Test unknown technology detection."""
        info = self.create_test_dataset_info("UNKNOWN", "Unknown Instrument")
        assert detect_sequencing_technology(info) == "unknown"


class TestSRAIntegration:
    """Integration tests for SRA functionality."""

    def test_metadata_to_dataframe_conversion(self):
        """Test converting metadata to DataFrame format."""
        dataset_info = SRADatasetInfo(
            accession="SRR123456",
            title="Test dataset",
            organism="Escherichia coli",
            platform="ILLUMINA",
            instrument="Illumina HiSeq 2500",
            strategy="WGS",
            layout="PAIRED",
            spots=1000000,
            bases=150000000,
            avg_length=150.0,
            size_mb=100.0,
            release_date="2023-01-01",
            bioproject="PRJNA123456",
            biosample="SAMN123456",
            library_selection="RANDOM",
            library_source="GENOMIC",
        )

        # Convert to dictionary (as would be done for CSV export)
        record = {
            "accession": dataset_info.accession,
            "platform": dataset_info.platform,
            "technology": detect_sequencing_technology(dataset_info),
            "spots": dataset_info.spots,
            "bases": dataset_info.bases,
        }

        assert record["accession"] == "SRR123456"
        assert record["technology"] == "illumina"
        assert record["spots"] == 1000000

    @patch("prepare.data.sra_metadata.save_metadata_report")
    def test_report_generation_integration(self, mock_save):
        """Test integration of metadata and report generation."""
        metadata = {
            "SRR123": SRADatasetInfo(
                accession="SRR123",
                title="Test",
                organism="E. coli",
                platform="ILLUMINA",
                instrument="HiSeq",
                strategy="WGS",
                layout="PAIRED",
                spots=1000,
                bases=150000,
                avg_length=150.0,
                size_mb=100.0,
                release_date="2023-01-01",
                bioproject="",
                biosample="",
                library_selection="",
                library_source="",
            )
        }

        # This should not raise an exception
        mock_save(metadata, "test_report.csv")
        mock_save.assert_called_once()


class TestSRAMetadataClientRetries:
    """``SRAMetadataClient`` retries a transient NCBI failure before raising ``NetworkError``.

    Uses ``tests/fake_http.py``, the same urllib3-connection-pool-level fake the taxonomy and
    GTDB retry tests use, so the client's real retrying session (not a mock of ``requests.get``)
    is exercised end to end.
    """

    def test_persistent_503_raises_network_error_after_retrying(self, monkeypatch):
        from prepare.core.exceptions import NetworkError
        from tests.fake_http import serve

        requested = serve(monkeypatch, lambda path: (503, b"unavailable"))
        client = SRAMetadataClient("test@example.com")
        client.request_delay = 0

        with pytest.raises(NetworkError, match="Failed to query NCBI"):
            client._make_request(client.base_url + "esearch.fcgi", {"db": "sra"})

        # One first attempt plus 3 retries (prepare.utils.http.RETRY_TOTAL).
        assert len(requested) == 4

    def test_success_after_one_503_returns_the_response_body(self, monkeypatch):
        from tests.fake_http import serve

        attempts = {"n": 0}

        def handler(path):
            attempts["n"] += 1
            if attempts["n"] < 2:
                return (503, b"unavailable")
            return (200, b'{"esearchresult": {"idlist": []}}')

        serve(monkeypatch, handler)
        client = SRAMetadataClient("test@example.com")
        client.request_delay = 0

        result = client._make_request(client.base_url + "esearch.fcgi", {"db": "sra"})

        assert result == '{"esearchresult": {"idlist": []}}'

    def test_client_owns_a_retrying_session(self):
        client = SRAMetadataClient("test@example.com")
        assert client.session.adapters["https://"].max_retries.raise_on_status is False


if __name__ == "__main__":
    pytest.main([__file__])
