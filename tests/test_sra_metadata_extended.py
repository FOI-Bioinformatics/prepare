"""
EXTENDED TESTS for data/sra_metadata.py (45% -> 75%+ coverage)

This file adds tests for untested methods:
- get_sra_metadata (batch processing)
- _fetch_batch_metadata (API integration)
- _parse_sra_xml (XML parsing)
- _extract_dataset_info (data extraction)
- _get_text (XML helper)
- generate_statistics_report (the sra_profile table writer)
- save_metadata_report
- create_download_preview

Run: pytest tests/test_sra_metadata_extended.py -v
"""

import pytest
import json
import requests
from unittest.mock import Mock, patch

from prepare.data.sra_metadata import (
    SRAMetadataClient,
    SRADatasetInfo,
    create_download_preview,
    save_metadata_report,
    generate_statistics_report,
)
from prepare.core.exceptions import DataAccessError

# Mock XML responses for testing (real efetch shape: RUN carries its own accession and
# numbers as attributes, not a nested Statistics child)
MOCK_SRA_XML = """<?xml version="1.0"?>
<EXPERIMENT_PACKAGE_SET>
    <EXPERIMENT_PACKAGE>
        <EXPERIMENT accession="SRX123456">
            <TITLE>Test Experiment</TITLE>
            <PLATFORM>
                <ILLUMINA>
                    <INSTRUMENT_MODEL>Illumina HiSeq 2500</INSTRUMENT_MODEL>
                </ILLUMINA>
            </PLATFORM>
            <DESIGN>
                <LIBRARY_DESCRIPTOR>
                    <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
                    <LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION>
                    <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
                    <LIBRARY_LAYOUT>
                        <PAIRED/>
                    </LIBRARY_LAYOUT>
                </LIBRARY_DESCRIPTOR>
            </DESIGN>
        </EXPERIMENT>
        <SAMPLE>
            <SCIENTIFIC_NAME>Escherichia coli</SCIENTIFIC_NAME>
        </SAMPLE>
        <STUDY>
            <EXTERNAL_ID namespace="BioProject">PRJNA123456</EXTERNAL_ID>
        </STUDY>
        <SUBMISSION accession="SRA100" received="2023-01-01"/>
        <RUN_SET>
            <RUN accession="SRR123456" total_spots="1000000" total_bases="150000000" size="100000000"
                 published="2023-01-02"/>
        </RUN_SET>
        <SAMPLE_ATTRIBUTE>
            <TAG>biosample</TAG>
            <VALUE>SAMN123456</VALUE>
        </SAMPLE_ATTRIBUTE>
    </EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>
"""


# ============================================================================
# TEST CLASS: SRAMetadataClient - API Integration
# ============================================================================


class TestSRAMetadataClientAPI:
    """Test SRA metadata client API methods."""

    def setup_method(self):
        """Set up test client."""
        self.client = SRAMetadataClient("test@example.com")

    def test_get_sra_metadata_empty_list(self):
        """Test handling of empty accession list."""
        result = self.client.get_sra_metadata([])
        assert result == {}

    def test_get_sra_metadata_single_batch(self):
        """Test fetching metadata for single batch."""
        accessions = ["SRR001", "SRR002"]

        with patch.object(self.client, "_fetch_batch_metadata") as mock_fetch:
            mock_fetch.return_value = {
                "SRR001": Mock(spec=SRADatasetInfo),
                "SRR002": Mock(spec=SRADatasetInfo),
            }

            result = self.client.get_sra_metadata(accessions)

            assert len(result) == 2
            assert "SRR001" in result
            assert "SRR002" in result
            mock_fetch.assert_called_once()

    def test_get_sra_metadata_multiple_batches(self):
        """Test fetching metadata in multiple batches."""
        # Create list that requires 2 batches (batch_size=200)
        accessions = [f"SRR{i:06d}" for i in range(250)]

        with patch.object(self.client, "_fetch_batch_metadata") as mock_fetch:
            mock_fetch.return_value = {}

            self.client.get_sra_metadata(accessions)

            # Should be called twice (200 + 50)
            assert mock_fetch.call_count == 2

    def test_get_sra_metadata_batch_failure_continues(self):
        """Test that batch failures don't stop processing."""
        accessions = [f"SRR{i:06d}" for i in range(250)]

        with patch.object(self.client, "_fetch_batch_metadata") as mock_fetch:
            # First batch fails, second succeeds
            mock_fetch.side_effect = [
                DataAccessError("API Error"),
                {"SRR000200": Mock(spec=SRADatasetInfo)},
            ]

            result = self.client.get_sra_metadata(accessions)

            # Should continue after first batch failure
            assert len(result) == 1
            assert mock_fetch.call_count == 2

    def test_fetch_batch_metadata_success(self):
        """Test successful batch metadata fetching."""
        accessions = ["SRR123456"]

        mock_search_response = json.dumps({"esearchresult": {"idlist": ["123456"]}})

        with patch.object(self.client, "_make_request") as mock_request:
            mock_request.side_effect = [
                mock_search_response,
                MOCK_SRA_XML,
            ]

            result = self.client._fetch_batch_metadata(accessions)

            assert len(result) > 0
            assert mock_request.call_count == 2

    def test_fetch_batch_metadata_no_results(self):
        """Test batch fetch when no results found."""
        accessions = ["NONEXISTENT"]

        mock_search_response = json.dumps({"esearchresult": {"idlist": []}})

        with patch.object(self.client, "_make_request") as mock_request:
            mock_request.return_value = mock_search_response

            result = self.client._fetch_batch_metadata(accessions)

            assert result == {}

    def test_fetch_batch_metadata_invalid_response(self):
        """Test batch fetch with invalid search response."""
        accessions = ["SRR001"]

        mock_search_response = json.dumps({"invalid": "response"})

        with patch.object(self.client, "_make_request") as mock_request:
            mock_request.return_value = mock_search_response

            result = self.client._fetch_batch_metadata(accessions)

            assert result == {}


# ============================================================================
# TEST CLASS: XML Parsing
# ============================================================================


class TestSRAXMLParsing:
    """Test SRA XML parsing methods."""

    def setup_method(self):
        """Set up test client."""
        self.client = SRAMetadataClient("test@example.com")

    def test_parse_sra_xml_success(self):
        """Test successful XML parsing."""
        result = self.client._parse_sra_xml(MOCK_SRA_XML)

        assert len(result) > 0
        assert "SRR123456" in result

        dataset = result["SRR123456"]
        assert dataset.accession == "SRR123456"
        assert dataset.platform == "ILLUMINA"
        assert dataset.spots == 1000000

    def test_parse_sra_xml_invalid_xml(self):
        """Test parsing invalid XML."""
        invalid_xml = "<invalid>xml</broken>"

        result = self.client._parse_sra_xml(invalid_xml)

        assert result == {}

    def test_parse_sra_xml_empty(self):
        """Test parsing empty XML."""
        empty_xml = "<?xml version='1.0'?><EXPERIMENT_PACKAGE_SET/>"

        result = self.client._parse_sra_xml(empty_xml)

        assert result == {}

    def test_parse_sra_xml_with_package_error(self):
        """Test XML parsing with malformed package."""
        malformed_xml = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE_SET>
            <EXPERIMENT_PACKAGE>
                <!-- Missing required fields -->
            </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>
        """

        result = self.client._parse_sra_xml(malformed_xml)

        # Should handle error gracefully
        assert isinstance(result, dict)

    def test_extract_dataset_info_complete(self):
        """Test extracting complete dataset info."""
        import xml.etree.ElementTree as ET

        root = ET.fromstring(MOCK_SRA_XML)
        package = root.find(".//EXPERIMENT_PACKAGE")

        results = self.client._extract_dataset_info(package)

        assert len(results) == 1
        result = results[0]
        assert result.accession == "SRR123456"
        assert result.platform == "ILLUMINA"
        assert result.instrument == "Illumina HiSeq 2500"
        assert result.layout == "PAIRED"
        assert result.organism == "Escherichia coli"
        assert result.bioproject == "PRJNA123456"
        assert result.biosample == "SAMN123456"
        assert result.spots == 1000000
        assert result.bases == 150000000
        assert result.release_date == "2023-01-02"

    def test_extract_dataset_info_no_experiment(self):
        """Test extraction when EXPERIMENT is missing."""
        import xml.etree.ElementTree as ET

        xml_no_experiment = "<EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE>"
        package = ET.fromstring(xml_no_experiment)

        result = self.client._extract_dataset_info(package)

        assert result == []

    def test_extract_dataset_info_single_layout(self):
        """Test extraction of SINGLE layout."""
        xml_single = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE>
            <EXPERIMENT accession="SRX999">
                <TITLE>Single End</TITLE>
                <PLATFORM>
                    <ILLUMINA>
                        <INSTRUMENT_MODEL>NextSeq</INSTRUMENT_MODEL>
                    </ILLUMINA>
                </PLATFORM>
                <DESIGN>
                    <LIBRARY_DESCRIPTOR>
                        <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
                        <LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>
                        <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
                        <LIBRARY_LAYOUT>
                            <SINGLE/>
                        </LIBRARY_LAYOUT>
                    </LIBRARY_DESCRIPTOR>
                </DESIGN>
            </EXPERIMENT>
            <RUN_SET>
                <RUN accession="SRR999" total_spots="500000" total_bases="75000000" size="50000000"
                     published="2023-02-01"/>
            </RUN_SET>
        </EXPERIMENT_PACKAGE>
        """

        import xml.etree.ElementTree as ET

        package = ET.fromstring(xml_single)

        results = self.client._extract_dataset_info(package)

        assert len(results) == 1
        result = results[0]
        assert result.accession == "SRR999"
        assert result.layout == "SINGLE"
        assert result.spots == 500000

    def test_extract_dataset_info_no_run_set_yields_one_record_keyed_by_experiment(self):
        """A package with no RUN_SET/RUN element at all (nothing has been submitted to SRA
        for this experiment yet, or the efetch response is trimmed) must still yield exactly
        one record, keyed by the EXPERIMENT accession, with zeroed run-level numbers rather
        than being silently dropped."""
        import xml.etree.ElementTree as ET

        xml_no_run_set = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE>
            <EXPERIMENT accession="SRX777">
                <TITLE>No runs yet</TITLE>
                <PLATFORM>
                    <ILLUMINA>
                        <INSTRUMENT_MODEL>NovaSeq</INSTRUMENT_MODEL>
                    </ILLUMINA>
                </PLATFORM>
                <DESIGN>
                    <LIBRARY_DESCRIPTOR>
                        <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
                        <LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION>
                        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
                        <LIBRARY_LAYOUT>
                            <SINGLE/>
                        </LIBRARY_LAYOUT>
                    </LIBRARY_DESCRIPTOR>
                </DESIGN>
            </EXPERIMENT>
        </EXPERIMENT_PACKAGE>
        """
        package = ET.fromstring(xml_no_run_set)

        results = self.client._extract_dataset_info(package)

        assert len(results) == 1
        result = results[0]
        assert result.accession == "SRX777"
        assert result.spots == 0
        assert result.bases == 0
        assert result.avg_length == 0.0

    def test_extract_dataset_info_run_without_accession_falls_back_to_experiment(self):
        """A RUN element present but missing its own ``accession`` attribute (a malformed or
        partial efetch record) must fall back to the EXPERIMENT accession rather than
        yielding a record keyed by an empty string."""
        import xml.etree.ElementTree as ET

        xml_run_no_accession = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE>
            <EXPERIMENT accession="SRX888">
                <TITLE>Run missing its own accession</TITLE>
                <PLATFORM>
                    <ILLUMINA>
                        <INSTRUMENT_MODEL>NovaSeq</INSTRUMENT_MODEL>
                    </ILLUMINA>
                </PLATFORM>
                <DESIGN>
                    <LIBRARY_DESCRIPTOR>
                        <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
                        <LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION>
                        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
                        <LIBRARY_LAYOUT>
                            <SINGLE/>
                        </LIBRARY_LAYOUT>
                    </LIBRARY_DESCRIPTOR>
                </DESIGN>
            </EXPERIMENT>
            <RUN_SET>
                <RUN total_spots="1000" total_bases="150000" size="100000" published="2023-02-01"/>
            </RUN_SET>
        </EXPERIMENT_PACKAGE>
        """
        package = ET.fromstring(xml_run_no_accession)

        results = self.client._extract_dataset_info(package)

        assert len(results) == 1
        result = results[0]
        assert result.accession == "SRX888"
        assert result.spots == 1000

    # Note: Removed test_extract_dataset_info_exception_handling because
    # xml.etree.ElementTree.Element.find is immutable and cannot be patched


# ============================================================================
# TEST CLASS: XML Helper Method
# ============================================================================


class TestGetTextHelper:
    """Test _get_text helper method."""

    def setup_method(self):
        """Set up test client."""
        self.client = SRAMetadataClient("test@example.com")

    def test_get_text_with_none_element(self):
        """Test get_text with None element."""
        result = self.client._get_text(None, ".//ANY", "default")
        assert result == "default"

    def test_get_text_attribute(self):
        """Test extracting attribute."""
        import xml.etree.ElementTree as ET

        xml = '<root attr="value"/>'
        elem = ET.fromstring(xml)

        result = self.client._get_text(elem, "./@attr", "default")
        assert result == "value"

    def test_get_text_nested_attribute(self):
        """Test extracting nested attribute."""
        import xml.etree.ElementTree as ET

        xml = '<root><child attr="nested_value"/></root>'
        elem = ET.fromstring(xml)

        result = self.client._get_text(elem, ".//child/@attr", "default")
        assert result == "nested_value"

    def test_get_text_element_text(self):
        """Test extracting element text."""
        import xml.etree.ElementTree as ET

        xml = "<root><child>text_value</child></root>"
        elem = ET.fromstring(xml)

        result = self.client._get_text(elem, ".//child", "default")
        assert result == "text_value"

    def test_get_text_missing_element(self):
        """Test get_text with missing element."""
        import xml.etree.ElementTree as ET

        xml = "<root/>"
        elem = ET.fromstring(xml)

        result = self.client._get_text(elem, ".//missing", "default")
        assert result == "default"

    def test_get_text_exception_handling(self):
        """Test exception handling in get_text."""
        import xml.etree.ElementTree as ET

        elem = ET.fromstring("<root/>")

        # Pass invalid xpath that might cause exception
        result = self.client._get_text(elem, "//invalid[xpath", "default")
        assert result == "default"


# ============================================================================
# TEST CLASS: Download Preview
# ============================================================================


class TestCreateDownloadPreview:
    """Test create_download_preview function."""

    def test_create_preview_success(self):
        """Test successful download preview creation."""
        mock_client = Mock(spec=SRAMetadataClient)

        mock_metadata = {
            "SRR001": SRADatasetInfo(
                accession="SRR001",
                title="Test 1",
                organism="E. coli",
                platform="ILLUMINA",
                instrument="HiSeq",
                strategy="WGS",
                layout="PAIRED",
                spots=1000000,
                bases=150000000,
                avg_length=150.0,
                size_mb=100.0,
                release_date="2023-01-01",
                bioproject="PRJNA001",
                biosample="SAMN001",
                library_selection="RANDOM",
                library_source="GENOMIC",
            ),
            "SRR002": SRADatasetInfo(
                accession="SRR002",
                title="Test 2",
                organism="S. aureus",
                platform="OXFORD_NANOPORE",
                instrument="MinION",
                strategy="WGS",
                layout="SINGLE",
                spots=500000,
                bases=500000000,
                avg_length=1000.0,
                size_mb=400.0,
                release_date="2023-01-02",
                bioproject="PRJNA002",
                biosample="SAMN002",
                library_selection="RANDOM",
                library_source="GENOMIC",
            ),
        }

        mock_client.get_sra_metadata.return_value = mock_metadata

        metadata, tech_counts, total_size_gb = create_download_preview(["SRR001", "SRR002"], mock_client)

        assert len(metadata) == 2
        assert tech_counts["illumina"] == 1
        assert tech_counts["nanopore"] == 1
        assert total_size_gb == pytest.approx(500.0 / 1024, rel=0.01)


# ============================================================================
# TEST CLASS: Save Metadata Report
# ============================================================================


class TestSaveMetadataReport:
    """Test save_metadata_report function."""

    def test_save_metadata_report_success(self, tmp_path):
        """Test successful metadata report saving."""
        metadata = {
            "SRR001": SRADatasetInfo(
                accession="SRR001",
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
                bioproject="PRJNA001",
                biosample="SAMN001",
                library_selection="RANDOM",
                library_source="GENOMIC",
            )
        }

        output_file = tmp_path / "metadata_report.csv"

        save_metadata_report(metadata, output_file)

        assert output_file.exists()

        # Check contents
        import pandas as pd

        df = pd.read_csv(output_file)
        assert len(df) == 1
        assert "SRR001" in df["accession"].values
        assert "illumina" in df["technology"].values

    def test_save_metadata_report_empty(self, tmp_path, caplog):
        """Test saving with empty metadata."""
        output_file = tmp_path / "empty_report.csv"

        save_metadata_report({}, output_file)

        assert "No metadata to save" in caplog.text


# ============================================================================
# TEST CLASS: Generate Statistics Report
# ============================================================================


class TestGenerateStatisticsReport:
    """generate_statistics_report writes the sra_profile table from rows it is given.

    Folder scanning, the statistics record and its cache moved to sra_profile in 0.5.0 and
    are tested in tests/test_cli_sra_profile.py.
    """

    @staticmethod
    def _row(accession, layout="PAIRED", total_reads=10, gc_percent=45.0, sampled=False):
        return {
            "accession": accession,
            "num_files": 2,
            "layout": layout,
            "total_reads": total_reads,
            "total_bases": total_reads * 100,
            "avg_read_length": 100.0,
            "gc_percent": gc_percent,
            "sampled": sampled,
        }

    def test_writes_one_row_per_accession_with_gc_in_percent(self, tmp_path):
        import pandas as pd

        output = tmp_path / "sra_statistics.csv"
        generate_statistics_report([self._row("SRR1"), self._row("SRR2", gc_percent=55.0)], output)

        df = pd.read_csv(output)
        assert list(df["accession"]) == ["SRR1", "SRR2"]
        assert list(df["gc_percent"]) == [45.0, 55.0]
        assert "gc_content" not in df.columns

    def test_returns_the_summary_lines_without_printing(self, tmp_path, capsys):
        lines = generate_statistics_report(
            [self._row("SRR1"), self._row("SRR2", layout="SINGLE", total_reads=20, gc_percent=55.0)],
            tmp_path / "report.csv",
        )
        # total_reads counts mates, not NCBI spots; the summary says so.
        assert "Total reads (mates counted): 30" in lines
        assert "Average GC content: 50.0%" in lines
        assert "  PAIRED: 1" in lines and "  SINGLE: 1" in lines
        assert capsys.readouterr().out == ""

    def test_says_when_some_totals_are_sample_counts(self, tmp_path):
        lines = generate_statistics_report([self._row("SRR1", sampled=True)], tmp_path / "report.csv")
        assert "Total reads (mates counted): 10 (lower bound: some totals are sample counts)" in lines

    def test_no_rows_writes_nothing(self, tmp_path, caplog):
        output = tmp_path / "report.csv"
        assert generate_statistics_report([], output) == []
        assert not output.exists()
        assert "No statistics to write" in caplog.text


# ============================================================================
# TEST CLASS: API Request Error Handling
# ============================================================================


class TestAPIRequestErrorHandling:
    """Test API request error handling."""

    def test_make_request_with_api_key(self):
        """Test request with API key."""
        client = SRAMetadataClient("test@example.com", "api_key_123")

        with patch.object(client.session, "get") as mock_get:
            mock_response = Mock()
            mock_response.text = "response"
            mock_response.raise_for_status = Mock()
            mock_get.return_value = mock_response

            client._make_request("http://test.com", {"param": "value"})

            # Verify API key was included
            call_args = mock_get.call_args
            assert "api_key" in call_args[1]["params"]
            assert call_args[1]["params"]["api_key"] == "api_key_123"

    def test_make_request_rate_limiting(self):
        """Test rate limiting between requests."""
        client = SRAMetadataClient("test@example.com")

        with patch.object(client.session, "get") as mock_get:
            with patch("time.sleep") as mock_sleep:
                mock_response = Mock()
                mock_response.text = "response"
                mock_response.raise_for_status = Mock()
                mock_get.return_value = mock_response

                # Make two requests quickly
                client._make_request("http://test.com", {})
                client._make_request("http://test.com", {})

                # Should have called sleep for rate limiting
                assert mock_sleep.called

    def test_make_request_exception(self):
        """Test handling of request exceptions."""
        client = SRAMetadataClient("test@example.com")

        with patch.object(client.session, "get") as mock_get:
            # Need to raise requests.RequestException (not plain Exception) for code to catch it
            mock_get.side_effect = requests.RequestException("Network error")

            with pytest.raises(DataAccessError, match="Failed to query NCBI"):
                client._make_request("http://test.com", {})


class TestEstimateDownloadTime:
    def test_scales_with_size_and_parallelism(self):
        from prepare.data.sra_metadata import estimate_download_time

        one = estimate_download_time(1.0, 100.0, 4)
        assert one == pytest.approx(1.0 * 1024 * 8 / (100.0 * 4 * 0.8) / 3600)
        assert estimate_download_time(2.0, 100.0, 4) == pytest.approx(2 * one)


# ============================================================================
# TEST CLASS: Real efetch XML shape (RUN attributes, not Statistics child)
# ============================================================================

REAL_EFETCH_XML = """<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE>
<EXPERIMENT accession="SRX100" alias="e"><TITLE>gut sample</TITLE>
<DESIGN><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY><LIBRARY_SOURCE>METAGENOMIC</LIBRARY_SOURCE>
<LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN>
<PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina NovaSeq 6000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM></EXPERIMENT>
<SUBMISSION accession="SRA100" received="2023-03-01"/>
<STUDY accession="SRP1"><IDENTIFIERS><EXTERNAL_ID namespace="BioProject">PRJNA1</EXTERNAL_ID></IDENTIFIERS></STUDY>
<SAMPLE accession="SRS1"><SAMPLE_NAME><SCIENTIFIC_NAME>gut metagenome</SCIENTIFIC_NAME></SAMPLE_NAME></SAMPLE>
<RUN_SET><RUN accession="SRR100" total_spots="4866463" total_bases="1459938900" size="482592813"
published="2023-03-23"/>
<RUN accession="SRR101" total_spots="10" total_bases="1500" size="1048576" published="2023-03-24"/></RUN_SET>
</EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>"""


def test_parse_real_efetch_shape_reports_runs():
    client = SRAMetadataClient(email="a@b.c")
    results = client._parse_sra_xml(REAL_EFETCH_XML)
    assert set(results) == {"SRR100", "SRR101"}
    run = results["SRR100"]
    assert run.spots == 4866463 and run.bases == 1459938900
    assert abs(run.size_mb - 482592813 / (1024 * 1024)) < 0.01
    assert run.release_date == "2023-03-23"
    assert run.strategy == "WGS" and run.layout == "PAIRED" and run.organism == "gut metagenome"
    assert run.bioproject == "PRJNA1"
    assert abs(run.avg_length - 300.0) < 0.01


def test_parse_sra_xml_requested_one_run_of_a_shared_package_keeps_its_sibling_run():
    """SRR100 and SRR101 sit in the same EXPERIMENT_PACKAGE (two lanes/runs of one
    experiment). Filtering happens at the package level, not per RUN: a request naming only
    SRR100 is a request for that package, so its sibling run SRR101 comes back too, rather
    than being dropped as "not requested". A request list may also legitimately hold an
    experiment, study or sample accession instead of a RUN accession (nothing about
    --accessions-file rules that out); per-RUN filtering would silently return nothing at all
    for such a request, which is the regression this rule avoids."""
    client = SRAMetadataClient(email="a@b.c")
    results = client._parse_sra_xml(REAL_EFETCH_XML, requested={"SRR100"})
    assert set(results) == {"SRR100", "SRR101"}


def test_parse_sra_xml_requested_experiment_accession_keeps_whole_package():
    """Requesting a package's EXPERIMENT (SRX) accession is as valid as naming one of its
    RUNs directly; every RUN in that package is returned."""
    client = SRAMetadataClient(email="a@b.c")
    results = client._parse_sra_xml(REAL_EFETCH_XML, requested={"SRX100"})
    assert set(results) == {"SRR100", "SRR101"}


def test_parse_sra_xml_requested_study_accession_keeps_whole_package():
    """Requesting a package's STUDY (SRP) accession likewise keeps every RUN in it."""
    client = SRAMetadataClient(email="a@b.c")
    results = client._parse_sra_xml(REAL_EFETCH_XML, requested={"SRP1"})
    assert set(results) == {"SRR100", "SRR101"}


def test_parse_sra_xml_requested_sample_accession_keeps_whole_package():
    """Requesting a package's SAMPLE (SRS) accession likewise keeps every RUN in it."""
    client = SRAMetadataClient(email="a@b.c")
    results = client._parse_sra_xml(REAL_EFETCH_XML, requested={"SRS1"})
    assert set(results) == {"SRR100", "SRR101"}


def test_parse_sra_xml_requested_matches_case_insensitively():
    """A lowercase (or any-case) requested accession still matches the XML's own casing."""
    client = SRAMetadataClient(email="a@b.c")
    results = client._parse_sra_xml(REAL_EFETCH_XML, requested={"srr100"})
    assert set(results) == {"SRR100", "SRR101"}


def test_parse_sra_xml_without_requested_keeps_every_run():
    """Called directly with no request set (a script, a REPL, or a test that hands it XML on
    its own), every RUN in the package is still returned; filtering only applies when a
    caller names the accessions it actually asked for."""
    client = SRAMetadataClient(email="a@b.c")
    results = client._parse_sra_xml(REAL_EFETCH_XML)
    assert set(results) == {"SRR100", "SRR101"}


XML_TWO_PACKAGES = """<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_PACKAGE_SET>
<EXPERIMENT_PACKAGE>
<EXPERIMENT accession="SRX100" alias="e"><TITLE>gut sample</TITLE>
<DESIGN><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY><LIBRARY_SOURCE>METAGENOMIC</LIBRARY_SOURCE>
<LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN>
<PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina NovaSeq 6000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM></EXPERIMENT>
<SUBMISSION accession="SRA100" received="2023-03-01"/>
<STUDY accession="SRP1"><IDENTIFIERS><EXTERNAL_ID namespace="BioProject">PRJNA1</EXTERNAL_ID></IDENTIFIERS></STUDY>
<SAMPLE accession="SRS1"><SAMPLE_NAME><SCIENTIFIC_NAME>gut metagenome</SCIENTIFIC_NAME></SAMPLE_NAME></SAMPLE>
<RUN_SET><RUN accession="SRR100" total_spots="4866463" total_bases="1459938900" size="482592813"
published="2023-03-23"/></RUN_SET>
</EXPERIMENT_PACKAGE>
<EXPERIMENT_PACKAGE>
<EXPERIMENT accession="SRX200" alias="e2"><TITLE>soil sample</TITLE>
<DESIGN><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY><LIBRARY_SOURCE>METAGENOMIC</LIBRARY_SOURCE>
<LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN>
<PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina NovaSeq 6000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM></EXPERIMENT>
<SUBMISSION accession="SRA200" received="2023-04-01"/>
<STUDY accession="SRP2"><IDENTIFIERS><EXTERNAL_ID namespace="BioProject">PRJNA2</EXTERNAL_ID></IDENTIFIERS></STUDY>
<SAMPLE accession="SRS2"><IDENTIFIERS><PRIMARY_ID>SRS2</PRIMARY_ID>
<EXTERNAL_ID namespace="BioSample">SAMN2</EXTERNAL_ID></IDENTIFIERS>
<SAMPLE_NAME><SCIENTIFIC_NAME>soil metagenome</SCIENTIFIC_NAME></SAMPLE_NAME></SAMPLE>
<RUN_SET><RUN accession="SRR200" total_spots="1000" total_bases="150000" size="100000"
published="2023-04-02"/></RUN_SET>
</EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>"""


def test_parse_sra_xml_requested_excludes_runs_from_a_different_package():
    """SRR100 and SRR200 sit in different EXPERIMENT_PACKAGEs (different experiments);
    requesting only SRR100 does not pull in SRR200's package, unlike the shared-package case
    above."""
    client = SRAMetadataClient(email="a@b.c")
    results = client._parse_sra_xml(XML_TWO_PACKAGES, requested={"SRR100"})
    assert set(results) == {"SRR100"}


def test_parse_sra_xml_requested_bioproject_keeps_its_package():
    """A BioProject accession is carried as an EXTERNAL_ID of the package's STUDY; a request
    for it keeps that package and not the other one."""
    client = SRAMetadataClient(email="a@b.c")
    results = client._parse_sra_xml(XML_TWO_PACKAGES, requested={"PRJNA1"})
    assert set(results) == {"SRR100"}


def test_parse_sra_xml_requested_biosample_keeps_its_package():
    """A BioSample accession is carried as an EXTERNAL_ID of the package's SAMPLE; matching is
    case-insensitive like the other levels."""
    client = SRAMetadataClient(email="a@b.c")
    results = client._parse_sra_xml(XML_TWO_PACKAGES, requested={"samn2"})
    assert set(results) == {"SRR200"}


def test_parse_sra_xml_requested_matching_nothing_keeps_the_batch(caplog):
    """When the requested accessions match no package of a non-empty reply, every run is
    listed with a WARNING instead of reporting a false failure."""
    client = SRAMetadataClient(email="a@b.c")
    with caplog.at_level("WARNING"):
        results = client._parse_sra_xml(REAL_EFETCH_XML, requested={"nomatch"})
    assert set(results) == {"SRR100", "SRR101"}
    assert any(
        "requested accessions matched no package in the reply; listing every run returned" in r.message
        for r in caplog.records
    )


def test_parse_sra_xml_package_inspection_failure_is_isolated_per_package(caplog):
    """One EXPERIMENT_PACKAGE's match check raising must not blank out every package's
    results: the failing package is logged and treated as not matched, while a package that
    matches normally still returns its runs. XML_TWO_PACKAGES holds two separate packages
    (SRX100/SRR100 and SRX200/SRR200); the check is made to raise only for the SRX200
    package."""
    client = SRAMetadataClient(email="a@b.c")
    real_matches = SRAMetadataClient._package_matches_requested

    def flaky_matches(package, requested_upper):
        if package.find(".//EXPERIMENT").get("accession") == "SRX200":
            raise ValueError("boom")
        return real_matches(package, requested_upper)

    with patch.object(SRAMetadataClient, "_package_matches_requested", side_effect=flaky_matches):
        with caplog.at_level("WARNING"):
            results = client._parse_sra_xml(XML_TWO_PACKAGES, requested={"SRX100"})

    assert set(results) == {"SRR100"}
    assert any("Failed to inspect dataset package" in r.message for r in caplog.records)


def test_parse_sra_xml_inspection_failure_with_no_match_stays_empty(caplog):
    """When a package's match check raises and no successfully-inspected package matches
    either, the result must stay empty rather than falling back to "keep everything": that
    fallback is only for a clean inspection that genuinely found no match, and applying it
    here would resurrect the package that could not even be checked. The warning must name
    the inspection failure, not claim the request "matched no package"."""
    client = SRAMetadataClient(email="a@b.c")
    real_matches = SRAMetadataClient._package_matches_requested

    def flaky_matches(package, requested_upper):
        if package.find(".//EXPERIMENT").get("accession") == "SRX100":
            raise ValueError("boom")
        return real_matches(package, requested_upper)

    with patch.object(SRAMetadataClient, "_package_matches_requested", side_effect=flaky_matches):
        with caplog.at_level("WARNING"):
            results = client._parse_sra_xml(XML_TWO_PACKAGES, requested={"nomatch"})

    assert results == {}
    assert any("could not be inspected" in r.message for r in caplog.records)
    assert not any("matched no package" in r.message for r in caplog.records)


XML_PACKAGE_WITHOUT_RUN_SET = """<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_PACKAGE_SET>
<EXPERIMENT_PACKAGE>
<EXPERIMENT accession="SRX9" alias="e9"><TITLE>study-level record</TITLE>
<DESIGN><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY><LIBRARY_SOURCE>METAGENOMIC</LIBRARY_SOURCE>
<LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN>
<PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina NovaSeq 6000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM></EXPERIMENT>
<SUBMISSION accession="SRA900" received="2023-05-01"/>
<STUDY accession="SRP9"><IDENTIFIERS><EXTERNAL_ID namespace="BioProject">PRJNA9</EXTERNAL_ID></IDENTIFIERS></STUDY>
<SAMPLE accession="SRS9"><SAMPLE_NAME><SCIENTIFIC_NAME>test metagenome</SCIENTIFIC_NAME></SAMPLE_NAME></SAMPLE>
</EXPERIMENT_PACKAGE>
<EXPERIMENT_PACKAGE>
<EXPERIMENT accession="SRX100" alias="e"><TITLE>gut sample</TITLE>
<DESIGN><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY><LIBRARY_SOURCE>METAGENOMIC</LIBRARY_SOURCE>
<LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN>
<PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina NovaSeq 6000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM></EXPERIMENT>
<SUBMISSION accession="SRA100" received="2023-03-01"/>
<STUDY accession="SRP1"><IDENTIFIERS><EXTERNAL_ID namespace="BioProject">PRJNA1</EXTERNAL_ID></IDENTIFIERS></STUDY>
<SAMPLE accession="SRS1"><SAMPLE_NAME><SCIENTIFIC_NAME>gut metagenome</SCIENTIFIC_NAME></SAMPLE_NAME></SAMPLE>
<RUN_SET><RUN accession="SRR100" total_spots="4866463" total_bases="1459938900" size="482592813"
published="2023-03-23"/></RUN_SET>
</EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>"""


def test_parse_sra_xml_requested_filter_with_package_without_run_set():
    """A package with no RUN_SET still yields one record keyed by its EXPERIMENT accession
    (see _extract_dataset_info's docstring); filtering by that accession keeps only that
    record and drops the other, normal package."""
    client = SRAMetadataClient(email="a@b.c")
    results = client._parse_sra_xml(XML_PACKAGE_WITHOUT_RUN_SET, requested={"SRX9"})
    assert set(results) == {"SRX9"}


def test_fetch_batch_metadata_filters_to_the_requested_batch():
    """sra_info's underlying batch fetch must not surface a RUN from a package the caller
    never asked for, even when NCBI's efetch response for the batch bundles another
    experiment's package alongside the requested one."""
    client = SRAMetadataClient(email="a@b.c")
    mock_search_response = json.dumps({"esearchresult": {"idlist": ["100"]}})

    with patch.object(client, "_make_request") as mock_request:
        mock_request.side_effect = [mock_search_response, XML_TWO_PACKAGES]
        result = client._fetch_batch_metadata(["SRR100"])

    assert set(result) == {"SRR100"}


def test_fetch_batch_metadata_uses_accn_terms_and_warns_about_missing_runs(caplog):
    """The esearch term restricts each accession to the accession field, retmax leaves
    headroom for several records per accession, and runs absent from the result are logged."""
    client = SRAMetadataClient(email="a@b.c")
    mock_search_response = json.dumps({"esearchresult": {"idlist": ["100"]}})

    with patch.object(client, "_make_request") as mock_request:
        mock_request.side_effect = [mock_search_response, XML_TWO_PACKAGES]
        with caplog.at_level("WARNING", logger="prepare.data.sra_metadata"):
            result = client._fetch_batch_metadata(["SRR100", "SRR555"])

    params = mock_request.call_args_list[0][0][1]
    assert params["term"] == "SRR100[ACCN] OR SRR555[ACCN]"
    assert params["retmax"] == "4"
    assert set(result) == {"SRR100"}
    assert "SRR555" in caplog.text


def test_fetch_batch_metadata_missing_run_check_ignores_case(caplog):
    client = SRAMetadataClient(email="a@b.c")
    mock_search_response = json.dumps({"esearchresult": {"idlist": ["100"]}})

    with patch.object(client, "_make_request") as mock_request:
        mock_request.side_effect = [mock_search_response, XML_TWO_PACKAGES]
        with caplog.at_level("WARNING", logger="prepare.data.sra_metadata"):
            client._fetch_batch_metadata(["srr100"])

    assert "No metadata returned" not in caplog.text


XML_RUN_MISSING_ATTRS_BUT_STATISTICS_CHILD = """<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE>
<EXPERIMENT accession="SRX200" alias="e"><TITLE>soil sample</TITLE>
<DESIGN><LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY><LIBRARY_SOURCE>METAGENOMIC</LIBRARY_SOURCE>
<LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN>
<PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina NovaSeq 6000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM></EXPERIMENT>
<SUBMISSION accession="SRA200" received="2023-04-01"/>
<STUDY accession="SRP2"><IDENTIFIERS><EXTERNAL_ID namespace="BioProject">PRJNA2</EXTERNAL_ID></IDENTIFIERS></STUDY>
<SAMPLE accession="SRS2"><SAMPLE_NAME><SCIENTIFIC_NAME>soil metagenome</SCIENTIFIC_NAME></SAMPLE_NAME></SAMPLE>
<RUN_SET><RUN accession="SRR200" published="2023-04-02">
<Statistics nspots="5000" nbases="750000"/>
</RUN></RUN_SET>
</EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>"""


def test_parse_falls_back_to_statistics_child_when_run_attributes_missing():
    """A RUN element without total_spots/total_bases attributes still reports real numbers,
    read from its Statistics child, rather than the zeroed defaults."""
    client = SRAMetadataClient(email="a@b.c")
    results = client._parse_sra_xml(XML_RUN_MISSING_ATTRS_BUT_STATISTICS_CHILD)
    run = results["SRR200"]
    assert run.spots == 5000
    assert run.bases == 750000


# ============================================================================
# SUCCESS METRICS:
#
# After running these tests:
# - Expected: 40+ additional tests pass
# - Coverage: 45% -> 75%+ for data/sra_metadata.py
# - All untested methods now covered
#
# Run tests:
#   pytest tests/test_sra_metadata_extended.py -v
#
# Check coverage:
#   pytest --cov=prepare.data.sra_metadata --cov-report=term-missing \
#          tests/test_sra_metadata_client.py tests/test_sra_metadata_extended.py
# ============================================================================


# ============================================================================
# Narrow exception handling (tech debt): a programming error must not be swallowed
# ============================================================================


def test_unexpected_error_in_batch_metadata_propagates(monkeypatch):
    """Kind (a): a bug in a batch fetch is no longer logged and skipped."""
    client = SRAMetadataClient(email="a@b.c")

    def buggy(batch):
        raise TypeError("bug")

    monkeypatch.setattr(client, "_fetch_batch_metadata", buggy)
    with pytest.raises(TypeError):
        client.get_sra_metadata(["SRR1"])


def test_ncbi_error_in_batch_metadata_is_still_skipped(monkeypatch):
    """Kind (a): an NCBI failure skips the batch and the call returns what it has."""
    client = SRAMetadataClient(email="a@b.c")

    def failing(batch):
        raise DataAccessError("NCBI down")

    monkeypatch.setattr(client, "_fetch_batch_metadata", failing)
    assert client.get_sra_metadata(["SRR1"]) == {}


def test_parse_sra_xml_returns_empty_for_malformed_xml_and_propagates_a_bug(monkeypatch):
    """Kind (b): malformed XML yields the default; a bug in extraction propagates."""
    client = SRAMetadataClient(email="a@b.c")
    assert client._parse_sra_xml("<not xml") == {}

    def buggy(package):
        raise TypeError("bug")

    monkeypatch.setattr(client, "_extract_dataset_info", buggy)
    with pytest.raises(TypeError):
        client._parse_sra_xml(MOCK_SRA_XML)
