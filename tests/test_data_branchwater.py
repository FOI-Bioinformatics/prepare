"""
Test data.branchwater module functionality.

Tests for Branchwater file processing, metadata extraction, and containment analysis.
"""

from pathlib import Path
from unittest.mock import Mock, patch
import numpy as np
import pandas as pd
import pytest
from collections import defaultdict

from prepare.data.branchwater import (
    process_branchwater_files,
    _process_branchwater_row,
    _validate_branchwater_file,
    extract_metadata_from_branchwater,
    _finalize_metadata_extraction,
    _process_genome_containments,
    _generate_containment_summary,
    parse_containment_data,
    _threshold_counts,
)
from prepare.core.exceptions import DataAccessError, ValidationError
from prepare.core.models import ContainmentSummary


class TestProcessBranchwaterFiles:
    """Test process_branchwater_files function."""

    def test_process_branchwater_files_success(self, tmp_path):
        """Test successful processing of Branchwater files."""
        source_dir = tmp_path / "source"
        target_dir = tmp_path / "target"
        source_dir.mkdir()

        # Create test CSV files with valid branchwater format
        csv_content = "acc,containment,ksize,scaled\nERR123,0.95,31,1000\nERR456,0.87,31,1000"
        (source_dir / "genome1.csv").write_text(csv_content)
        (source_dir / "genome2.csv").write_text(csv_content)

        with patch("prepare.data.branchwater.read_csv") as mock_read_csv:
            # Mock the header check for format validation
            mock_df = pd.DataFrame(columns=["acc", "containment", "ksize", "scaled"])
            mock_read_csv.return_value = mock_df

            result = process_branchwater_files(source_dir, target_dir)

        assert len(result) == 2
        assert "genome1" in result
        assert "genome2" in result
        assert target_dir.exists()

    def test_process_branchwater_files_skip_existing(self, tmp_path):
        """Test skipping files that already exist in target."""
        source_dir = tmp_path / "source"
        target_dir = tmp_path / "target"
        source_dir.mkdir()
        target_dir.mkdir()

        # Create source file
        csv_content = "acc,containment,ksize,scaled\nERR123,0.95,31,1000"
        (source_dir / "genome1.csv").write_text(csv_content)

        # An identical copy is already there: nothing to do.
        (target_dir / "genome1.csv").write_text(csv_content)

        result = process_branchwater_files(source_dir, target_dir)

        assert len(result) == 1
        assert "genome1" in result
        assert (target_dir / "genome1.csv").read_text() == csv_content

    def test_process_branchwater_files_refreshes_a_changed_source(self, tmp_path, caplog):
        """A matches file that differs from its Branchwater source is replaced, with a warning."""
        import logging

        source_dir = tmp_path / "source"
        target_dir = tmp_path / "target"
        source_dir.mkdir()
        target_dir.mkdir()
        (source_dir / "genome1.csv").write_text(
            "acc,containment,ksize,scaled\nERR123,0.95,31,1000\nERR124,0.5,31,1000\n"
        )
        (target_dir / "genome1.csv").write_text("acc,containment,ksize,scaled\nERR123,0.95,31,1000\n")

        with caplog.at_level(logging.WARNING):
            result = process_branchwater_files(source_dir, target_dir)

        assert "genome1" in result
        assert (target_dir / "genome1.csv").read_text() == (source_dir / "genome1.csv").read_text()
        assert "genome1.csv" in caplog.text and "Replacing" in caplog.text

    def test_process_branchwater_files_no_csv_files(self, tmp_path):
        """Test with no CSV files in source directory."""
        source_dir = tmp_path / "source"
        target_dir = tmp_path / "target"
        source_dir.mkdir()

        # Create non-CSV file
        (source_dir / "not_csv.txt").write_text("not a csv")

        result = process_branchwater_files(source_dir, target_dir)

        assert result == {}

    def test_process_branchwater_files_invalid_format(self, tmp_path):
        """Test with files that have invalid format."""
        source_dir = tmp_path / "source"
        target_dir = tmp_path / "target"
        source_dir.mkdir()

        # Create CSV with invalid format
        invalid_content = "invalid,header,format\ndata1,data2,data3"
        (source_dir / "invalid.csv").write_text(invalid_content)

        result = process_branchwater_files(source_dir, target_dir)

        # Should have empty result due to validation errors
        assert result == {}

    def test_process_branchwater_files_bad_header_recorded(self, tmp_path):
        """A file without the required columns is not copied and is reported."""
        source_dir = tmp_path / "source"
        target_dir = tmp_path / "target"
        source_dir.mkdir()
        (source_dir / "bad.csv").write_text("foo,bar\n1,2\n")

        errors = []
        result = process_branchwater_files(source_dir, target_dir, errors)

        assert result == {}
        assert errors == ["bad.csv"]
        assert not (target_dir / "bad.csv").exists()

    def test_process_branchwater_files_source_not_exist(self):
        """Test with non-existent source folder."""
        with pytest.raises(ValidationError):
            process_branchwater_files("/nonexistent/source", "/tmp/target")

    def test_unreadable_csv_is_reported_to_caller(self, tmp_path):
        source = tmp_path / "bw"
        target = tmp_path / "matches"
        source.mkdir()
        (source / "good.csv").write_text("acc,containment,cANI\nSRR1,0.9,0.99\n")
        (source / "bad.csv").write_bytes(b"\x00\x05\x16\x07\xb0")
        errors: list = []
        result = process_branchwater_files(source, target, errors=errors)
        assert "good" in result
        assert errors == ["bad.csv"]


class TestProcessBranchwaterRow:
    """Test _process_branchwater_row function."""

    def test_process_branchwater_row_basic(self):
        """Test basic row processing."""
        row = pd.Series(
            {
                "acc": "ERR123456",
                "biosample": "SAMN123",
                "organism": "Escherichia coli",
                "geo_loc_name_country_calc": "USA",
            }
        )
        metadata_records = []

        _process_branchwater_row(row, metadata_records)

        assert len(metadata_records) == 1
        record = metadata_records[0]
        assert record["Run_ID"] == "ERR123456"
        assert record["Sample_ID"] == "SAMN123"
        assert record["Sample_Scientific_Name"] == "Escherichia coli"
        assert record["geo_loc_name_country_calc"] == "USA"

    def test_process_branchwater_row_missing_acc(self):
        """Test row processing with missing acc column."""
        row = pd.Series({"other_field": "value"})
        metadata_records = []

        with patch("prepare.data.branchwater.logger") as mock_logger:
            _process_branchwater_row(row, metadata_records)

        assert len(metadata_records) == 0
        mock_logger.warning.assert_called_once()

    def test_process_branchwater_row_partial_fields(self):
        """Test row processing with only some optional fields."""
        row = pd.Series({"acc": "ERR789", "organism": "Salmonella enterica", "unknown_field": "ignored", "cANI": 0.98})
        metadata_records = []

        _process_branchwater_row(row, metadata_records)

        assert len(metadata_records) == 1
        record = metadata_records[0]
        assert record["Run_ID"] == "ERR789"
        assert record["Sample_Scientific_Name"] == "Salmonella enterica"
        assert record["cANI"] == 0.98
        assert "unknown_field" not in record

    def test_process_branchwater_row_with_nan_values(self):
        """Test row processing with NaN values."""
        row = pd.Series({"acc": "ERR999", "organism": "Bacteria", "biosample": pd.NA, "collection_date_sam": None})
        metadata_records = []

        _process_branchwater_row(row, metadata_records)

        assert len(metadata_records) == 1
        record = metadata_records[0]
        assert record["Run_ID"] == "ERR999"
        assert record["Sample_Scientific_Name"] == "Bacteria"
        assert "Sample_ID" not in record
        assert "collection_date_sam" not in record


class TestValidateBranchwaterFile:
    """Test _validate_branchwater_file function."""

    def test_validate_branchwater_file_valid(self, tmp_path):
        """Test validation of valid Branchwater file."""
        csv_file = tmp_path / "valid.csv"
        csv_content = "acc,containment,ksize,scaled\nERR123,0.95,31,1000"
        csv_file.write_text(csv_content)

        result = _validate_branchwater_file(csv_file)

        assert result is True

    def test_validate_branchwater_file_missing_acc(self, tmp_path):
        """Test validation with missing acc column."""
        csv_file = tmp_path / "invalid.csv"
        csv_content = "run_id,containment,ksize,scaled\nERR123,0.95,31,1000"
        csv_file.write_text(csv_content)

        with patch("prepare.data.branchwater.logger") as mock_logger:
            result = _validate_branchwater_file(csv_file)

        assert result is False
        mock_logger.warning.assert_called_once()

    def test_validate_branchwater_file_missing_containment(self, tmp_path):
        """Test validation with missing containment column."""
        csv_file = tmp_path / "invalid.csv"
        csv_content = "acc,similarity,ksize,scaled\nERR123,0.95,31,1000"
        csv_file.write_text(csv_content)

        with patch("prepare.data.branchwater.logger") as mock_logger:
            result = _validate_branchwater_file(csv_file)

        assert result is False
        mock_logger.warning.assert_called_once()


class TestExtractMetadataFromBranchwater:
    """Test extract_metadata_from_branchwater function."""

    def test_extract_metadata_success(self, tmp_path):
        """Test successful metadata extraction."""
        source_dir = tmp_path / "source"
        source_dir.mkdir()
        output_file = tmp_path / "metadata.csv"

        # Create test CSV file
        csv_content = "acc,containment,organism,biosample\nERR123,0.95,E. coli,SAMN123\nERR456,0.87,S. enterica,SAMN456"
        (source_dir / "test.csv").write_text(csv_content)

        result = extract_metadata_from_branchwater(source_dir, output_file)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2
        assert "Run_ID" in result.columns
        assert output_file.exists()

    def test_unreadable_file_is_reported_and_summary_is_warning(self, tmp_path, caplog):
        """A genuinely unreadable file (not a mocked exception) must be named in ``errors``
        and must not stop the other file in the same folder from being extracted; with at
        least one error, the "Processed N files with M errors" summary line must log at
        WARNING, not INFO."""
        import logging

        source_dir = tmp_path / "source"
        source_dir.mkdir()
        output_file = tmp_path / "metadata.csv"

        (source_dir / "good.csv").write_text("acc,containment,organism\nERR123,0.95,E. coli\n")
        (source_dir / "bad.csv").write_bytes(b"\x00\x05\x16\x07\xb0")

        errors: list = []
        # Captured from INFO up, not just WARNING: a regression that logged the summary line
        # at INFO would otherwise be filtered out of caplog entirely (both by caplog's own
        # capture level and by the logger's effective level), so the assertion below would
        # fail on an empty `summary_records` list rather than on the specific level mismatch.
        with caplog.at_level(logging.INFO):
            result = extract_metadata_from_branchwater(source_dir, output_file, errors=errors)

        assert errors == ["bad.csv"]
        assert list(result["Run_ID"]) == ["ERR123"]
        summary_records = [r for r in caplog.records if "Processed" in r.message and "errors" in r.message]
        assert summary_records and summary_records[0].levelno == logging.WARNING

    def test_extract_metadata_api_shaped_csv_warns(self, tmp_path, caplog):
        """A Branchwater-API CSV (acc/containment/cANI only) has no SRA metadata to extract;
        warn and point at download_metadata rather than silently writing a near-empty table."""
        source_dir = tmp_path / "source"
        source_dir.mkdir()
        output_file = tmp_path / "metadata.csv"

        csv_content = "acc,containment,cANI\nERR123,0.95,0.98\nERR456,0.87,0.96\n"
        (source_dir / "test.csv").write_text(csv_content)

        with caplog.at_level("WARNING"):
            result = extract_metadata_from_branchwater(source_dir, output_file)

        assert set(result.columns) == {"Run_ID", "cANI"}
        assert "download_metadata" in caplog.text

    def test_extract_metadata_no_csv_files(self, tmp_path):
        """Test with no CSV files."""
        source_dir = tmp_path / "source"
        source_dir.mkdir()
        output_file = tmp_path / "metadata.csv"

        result = extract_metadata_from_branchwater(source_dir, output_file)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_extract_metadata_nonexistent_folder(self, tmp_path):
        """Test with non-existent source folder."""
        output_file = tmp_path / "metadata.csv"

        with pytest.raises(DataAccessError, match="Branchwater folder does not exist"):
            extract_metadata_from_branchwater("/nonexistent", output_file)

    def test_extract_metadata_invalid_files(self, tmp_path):
        """Test with invalid file formats."""
        source_dir = tmp_path / "source"
        source_dir.mkdir()
        output_file = tmp_path / "metadata.csv"

        # Create invalid CSV file
        invalid_content = "invalid,format\ndata1,data2"
        (source_dir / "invalid.csv").write_text(invalid_content)

        result = extract_metadata_from_branchwater(source_dir, output_file)

        # Should return empty DataFrame but still create output file
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0
        assert output_file.exists()

    def test_extract_metadata_processing_error(self, tmp_path):
        """Test with file processing errors."""
        source_dir = tmp_path / "source"
        source_dir.mkdir()
        output_file = tmp_path / "metadata.csv"

        # Create file with valid header but invalid content
        csv_content = "acc,containment\ninvalid_content_here"
        (source_dir / "error.csv").write_text(csv_content)

        with patch("pandas.read_csv", side_effect=pd.errors.ParserError("Parse error")):
            result = extract_metadata_from_branchwater(source_dir, output_file)

        # Should handle error gracefully
        assert isinstance(result, pd.DataFrame)


class TestFinalizeMetadataExtraction:
    """Test _finalize_metadata_extraction function."""

    def test_finalize_metadata_with_records(self, tmp_path):
        """Test finalization with metadata records."""
        metadata_records = [
            {"Run_ID": "ERR123", "Sample_Scientific_Name": "E. coli"},
            {"Run_ID": "ERR456", "Sample_Scientific_Name": "S. enterica"},
            {"Run_ID": "ERR123", "Sample_Scientific_Name": "E. coli"},  # Duplicate
        ]
        output_file = tmp_path / "metadata.csv"

        result = _finalize_metadata_extraction(metadata_records, output_file, 2, 0)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2  # Duplicates should be removed
        assert output_file.exists()

    def test_finalize_metadata_empty_records(self, tmp_path):
        """Test finalization with no metadata records."""
        metadata_records = []
        output_file = tmp_path / "metadata.csv"

        result = _finalize_metadata_extraction(metadata_records, output_file, 0, 1)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0
        assert output_file.exists()

    def test_finalize_metadata_create_directory(self, tmp_path):
        """Test that output directory is created."""
        metadata_records = [{"Run_ID": "ERR123"}]
        output_file = tmp_path / "nested" / "dir" / "metadata.csv"

        _finalize_metadata_extraction(metadata_records, output_file, 1, 0)

        assert output_file.exists()
        assert output_file.parent.exists()

    def test_finalize_metadata_warns_when_only_run_id_and_cani_populated(self, tmp_path, caplog):
        metadata_records = [
            {"Run_ID": "ERR123", "cANI": 0.98},
            {"Run_ID": "ERR456", "cANI": 0.96},
        ]
        output_file = tmp_path / "metadata.csv"

        with caplog.at_level("WARNING"):
            _finalize_metadata_extraction(metadata_records, output_file, 1, 0)

        assert "download_metadata" in caplog.text

    def test_finalize_metadata_no_warning_with_real_metadata(self, tmp_path, caplog):
        metadata_records = [
            {"Run_ID": "ERR123", "cANI": 0.98, "Sample_Scientific_Name": "E. coli"},
        ]
        output_file = tmp_path / "metadata.csv"

        with caplog.at_level("WARNING"):
            _finalize_metadata_extraction(metadata_records, output_file, 1, 0)

        assert "download_metadata" not in caplog.text


class TestProcessGenomeContainments:
    """Test _process_genome_containments function."""

    def test_process_genome_containments_branchwater(self, tmp_path):
        """Test processing with Branchwater format."""
        csv_file = tmp_path / "genome1.csv"
        containment_data = defaultdict(dict)

        with patch("prepare.data.branchwater.read_csv") as mock_read_csv:
            # Mock header check
            mock_headers = pd.DataFrame(columns=["acc", "containment"])
            mock_read_csv.return_value = mock_headers

            with patch("prepare.data.branchwater.BranchWaterFormatPlugin") as mock_plugin:
                mock_plugin.require_headers.return_value = None
                mock_plugin.parse_file.return_value = [
                    Mock(accession="ERR123", value=0.95),
                    Mock(accession="ERR456", value=0.87),
                ]

                _process_genome_containments(csv_file, "genome1", containment_data)

        assert "ERR123" in containment_data
        assert "ERR456" in containment_data
        assert containment_data["ERR123"]["genome1"] == 0.95
        assert containment_data["ERR456"]["genome1"] == 0.87

    def test_process_genome_containments_unknown_format(self, tmp_path):
        """Test processing with unknown format."""
        csv_file = tmp_path / "genome1.csv"
        containment_data = defaultdict(dict)

        with patch("prepare.data.branchwater.read_csv") as mock_read_csv:
            mock_headers = pd.DataFrame(columns=["unknown", "format"])
            mock_read_csv.return_value = mock_headers

            with patch("prepare.data.branchwater.BranchWaterFormatPlugin") as mock_bw_plugin:
                mock_bw_plugin.require_headers.side_effect = ValidationError("Not branchwater")

                with pytest.raises(DataAccessError, match="Error processing"):
                    _process_genome_containments(csv_file, "genome1", containment_data)

    def test_process_genome_containments_max_value(self, tmp_path):
        """Test that maximum containment value is kept for duplicates."""
        csv_file = tmp_path / "genome1.csv"
        containment_data = defaultdict(dict)

        # Pre-populate with existing value
        containment_data["ERR123"]["genome1"] = 0.5

        with patch("prepare.data.branchwater.read_csv") as mock_read_csv:
            mock_headers = pd.DataFrame(columns=["acc", "containment"])
            mock_read_csv.return_value = mock_headers

            with patch("prepare.data.branchwater.BranchWaterFormatPlugin") as mock_plugin:
                mock_plugin.require_headers.return_value = None
                mock_plugin.parse_file.return_value = [Mock(accession="ERR123", value=0.95)]  # Higher value

                _process_genome_containments(csv_file, "genome1", containment_data)

        # Should keep the higher value
        assert containment_data["ERR123"]["genome1"] == 0.95


class TestGenerateContainmentSummary:
    """Test _generate_containment_summary function."""

    def test_generate_containment_summary_success(self, tmp_path):
        """Test successful summary generation."""
        containment_data = {
            "ERR123": {"genome1": 0.95, "genome2": 0.1},
            "ERR456": {"genome1": 0.87, "genome2": 0.95},
            "ERR789": {"genome1": 0.0, "genome2": 0.5},
        }
        output_file = tmp_path / "parsed.txt"
        summary_file = tmp_path / "summary.txt"

        with patch("prepare.data.branchwater.write_csv") as mock_write:
            result = _generate_containment_summary(containment_data, output_file, summary_file, 0.1)

        assert isinstance(result, ContainmentSummary)
        assert len(result.thresholds) > 0
        assert len(result.counts) > 0
        assert len(result.max_containment) == 3
        assert mock_write.call_count == 2  # parsed and summary files

    def test_generate_containment_summary_empty_data(self, tmp_path):
        """Test summary generation with empty data."""
        containment_data = {}
        output_file = tmp_path / "parsed.txt"
        summary_file = tmp_path / "summary.txt"

        # Empty data will create an empty DataFrame which pandas can handle
        with patch("prepare.data.branchwater.write_csv") as mock_write:
            result = _generate_containment_summary(containment_data, output_file, summary_file, 0.1)

        assert isinstance(result, ContainmentSummary)
        mock_write.assert_called()


class TestParseContainmentData:
    """Test parse_containment_data function."""

    def test_parse_containment_data_success(self, tmp_path):
        """Test successful containment data parsing."""
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()
        output_file = tmp_path / "parsed.txt"
        summary_file = tmp_path / "summary.txt"

        # Create test CSV file
        (matches_dir / "genome1.csv").write_text("acc,containment\nERR123,0.95")

        def mock_process_side_effect(csv_file, genome_id, containment_data, details_rows=None):
            # Simulate adding some data to containment_data
            containment_data["ERR123"]["genome1"] = 0.95

        with patch(
            "prepare.data.branchwater._process_genome_containments", side_effect=mock_process_side_effect
        ) as mock_process:
            with patch("prepare.data.branchwater._generate_containment_summary") as mock_generate:
                mock_summary = ContainmentSummary()
                mock_generate.return_value = mock_summary

                result = parse_containment_data(matches_dir, output_file, summary_file, 0.05)

        assert result == mock_summary
        mock_process.assert_called_once()
        mock_generate.assert_called_once()

    def test_parse_containment_data_processing_errors(self, tmp_path):
        """Every unreadable file is named in one error and nothing is written."""
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()
        output_file = tmp_path / "parsed.txt"
        summary_file = tmp_path / "summary.txt"

        # Create test files
        (matches_dir / "genome1.csv").write_text("data")
        (matches_dir / "genome2.csv").write_text("data")

        with patch("prepare.data.branchwater._process_genome_containments") as mock_process:
            mock_process.side_effect = DataAccessError("Processing error")

            with pytest.raises(DataAccessError, match="2 match file.*genome1.csv.*genome2.csv"):
                parse_containment_data(matches_dir, output_file, summary_file)

        assert not output_file.exists()

    def test_parse_containment_with_an_unreadable_file_writes_nothing(self, tmp_path):
        """One match file with the wrong header: nothing is written and the error names it."""
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()
        output_file = tmp_path / "parsed.txt"
        summary_file = tmp_path / "summary.txt"
        (matches_dir / "good.csv").write_text("acc,containment\nSRR1,0.9\n")
        (matches_dir / "bad.csv").write_text("run,score\nSRR2,0.5\n")

        with pytest.raises(DataAccessError, match="bad.csv"):
            parse_containment_data(matches_dir, output_file, summary_file)

        assert not output_file.exists()
        assert not summary_file.exists()
        assert not (tmp_path / "parsed_details.tsv").exists()

    def test_parse_containment_on_an_empty_matches_folder_raises(self, tmp_path):
        """A folder without match CSVs is an error that names the folder."""
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()
        (matches_dir / "notes.txt").write_text("x")

        with pytest.raises(DataAccessError, match="No match CSV in .*matches"):
            parse_containment_data(matches_dir, tmp_path / "parsed.txt", tmp_path / "summary.txt")

    def test_a_semicolon_separated_csv_names_its_separator(self, tmp_path):
        """A file separated by ';' is reported as such, not as missing columns."""
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()
        (matches_dir / "g.csv").write_text("acc;containment\nSRR1;0.9\n")

        with pytest.raises(DataAccessError, match=r"not comma separated \(found ';'\)"):
            parse_containment_data(matches_dir, tmp_path / "parsed.txt", tmp_path / "summary.txt")

    def test_parse_containment_data_invalid_folder(self):
        """Test with invalid matches folder."""
        with pytest.raises(ValidationError):
            parse_containment_data("/nonexistent", "output.txt", "summary.txt")

    def test_a_binary_match_file_is_named_and_the_good_one_is_not_written(self, tmp_path, caplog):
        """A genuinely unreadable match file (not a mocked exception) is named in the error,
        logged per file, and the good file in the same folder is not written out either."""
        import logging

        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()
        output_file = tmp_path / "parsed.txt"
        summary_file = tmp_path / "summary.txt"

        (matches_dir / "good_genome.csv").write_text(
            "acc,containment,cANI,biosample,bioproject,assay_type,organism,geo_loc_name_country_calc,lat_lon\n"
            "SRR1,0.9,0.98,SAMN1,PRJNA1,WGS,Salmonella enterica,USA,35.7N 100.2W\n"
        )
        (matches_dir / "bad_genome.csv").write_bytes(b"\x00\x05\x16\x07\xb0")

        with caplog.at_level(logging.ERROR):
            with pytest.raises(DataAccessError, match="bad_genome.csv"):
                parse_containment_data(matches_dir, output_file, summary_file)

        assert "bad_genome.csv" in caplog.text
        assert not output_file.exists()


class TestBranchwaterIntegration:
    """Integration tests for branchwater functionality."""

    def test_full_branchwater_workflow(self, tmp_path):
        """Test complete branchwater processing workflow with realistic Salmonella data."""
        # Set up directories
        source_dir = tmp_path / "source"
        target_dir = tmp_path / "target"
        source_dir.mkdir()

        # Use realistic Salmonella test data from salmonella_mini.csv
        real_data_path = Path(__file__).parent.parent / "test_data" / "branchwater" / "salmonella_mini.csv"
        if real_data_path.exists():
            # Copy real data for testing
            import shutil

            shutil.copy(real_data_path, source_dir / "salmonella_genome.csv")
            # Also create a second file for multi-file testing
            shutil.copy(real_data_path, source_dir / "ecoli_genome.csv")
        else:
            # Fallback to minimal synthetic data if real data not available
            csv_content = (
                "acc,containment,cANI,biosample,bioproject,assay_type,"
                "collection_date_sam,geo_loc_name_country_calc,organism,lat_lon\n"
                "ERR123456,0.95,0.98,SAMN123,PRJNA123,WGS,2023-01-01,USA,Salmonella enterica,35.7N 100.2W\n"
                "ERR789012,0.87,0.96,SAMN456,PRJNA456,WGS,2023-02-15,Denmark,gut metagenome,55.67N 12.57E"
            )

            (source_dir / "salmonella_genome.csv").write_text(csv_content)
            (source_dir / "ecoli_genome.csv").write_text(csv_content)

        # Mock the format validation to work with our realistic test data
        with patch("prepare.data.branchwater.read_csv") as mock_read_csv:
            mock_df = pd.DataFrame(
                columns=[
                    "acc",
                    "containment",
                    "cANI",
                    "biosample",
                    "bioproject",
                    "assay_type",
                    "collection_date_sam",
                    "geo_loc_name_country_calc",
                    "organism",
                    "lat_lon",
                ]
            )
            mock_read_csv.return_value = mock_df

            # Process files
            result_files = process_branchwater_files(source_dir, target_dir)

        # Verify file processing
        assert len(result_files) == 2
        assert "ecoli_genome" in result_files
        assert "salmonella_genome" in result_files

        # Test metadata extraction with realistic data
        metadata_file = tmp_path / "metadata.csv"

        with patch("prepare.data.branchwater._validate_branchwater_file", return_value=True):
            with patch("pandas.read_csv") as mock_pd_read:
                # Use realistic values from salmonella_mini.csv for better testing
                mock_pd_read.return_value = pd.DataFrame(
                    {
                        "acc": ["ERR2868097", "SRR22022422", "SRR31320585"],
                        "containment": [0.98, 0.98, 0.97],
                        "cANI": [1.0, 1.0, 1.0],
                        "biosample": ["SAMEA5056563", "SAMN29160454", "SAMN44590404"],
                        "bioproject": ["PRJEB29454", "PRJNA849983", "PRJNA1182286"],
                        "organism": ["metagenome", "food metagenome", "metagenome"],
                        "geo_loc_name_country_calc": ["NP", "United Kingdom", "USA"],
                        "collection_date_sam": ["", "2019-02-07", "2023-09-27"],
                    }
                )

                metadata_df = extract_metadata_from_branchwater(source_dir, metadata_file)

        # Verify metadata extraction
        assert isinstance(metadata_df, pd.DataFrame)
        assert metadata_file.exists()

    def test_error_handling_workflow(self, tmp_path):
        """Test error handling throughout the workflow."""
        source_dir = tmp_path / "source"
        target_dir = tmp_path / "target"
        source_dir.mkdir()

        # Create files that will cause various errors
        (source_dir / "corrupt.csv").write_text("corrupted,data\nwith,errors")
        (source_dir / "empty.csv").write_text("")

        # Should handle errors gracefully
        with patch("prepare.data.branchwater.logger") as mock_logger:
            result = process_branchwater_files(source_dir, target_dir)

        # Should have logged errors but not crashed
        assert isinstance(result, dict)
        mock_logger.error.assert_called()

    def test_realistic_salmonella_data_parsing(self, tmp_path):
        """Test parsing with actual realistic Salmonella data."""
        # Use the realistic salmonella_mini.csv data
        real_data_path = Path(__file__).parent.parent / "test_data" / "branchwater" / "salmonella_mini.csv"
        if not real_data_path.exists():
            pytest.skip("Realistic Salmonella test data not available")

        # Set up test environment
        source_dir = tmp_path / "source"
        source_dir.mkdir()
        output_file = tmp_path / "metadata.csv"

        # Copy real data
        import shutil

        shutil.copy(real_data_path, source_dir / "salmonella.csv")

        # Extract metadata using the real data structure
        result = extract_metadata_from_branchwater(source_dir, output_file)

        # Verify that we can process real data without errors
        assert isinstance(result, pd.DataFrame)
        assert output_file.exists()

        # If data was processed, verify it contains expected fields
        if len(result) > 0:
            assert "Run_ID" in result.columns
            # Verify some realistic entries are present
            run_ids = result["Run_ID"].tolist()
            assert any(run_id.startswith("ERR") or run_id.startswith("SRR") for run_id in run_ids)


class TestProcessGenomeContainmentsDetails:
    """Test that _process_genome_containments can collect a details row per containment."""

    def test_collects_details_rows_when_requested(self, tmp_path):
        csv_file = tmp_path / "GCF_1.csv"
        csv_file.write_text(
            "acc,containment,cANI,biosample,bioproject,assay_type,organism,geo_loc_name_country_calc,lat_lon\n"
            "SRR1,0.9,0.98,SAMN1,PRJNA1,WGS,Salmonella enterica,USA,35.7N 100.2W\n"
        )
        containment_data = defaultdict(dict)
        details_rows: list = []

        _process_genome_containments(csv_file, "GCF_1", containment_data, details_rows=details_rows)

        assert len(details_rows) == 1
        assert details_rows[0] == {
            "accession": "SRR1",
            "genome_id": "GCF_1",
            "containment": 0.9,
            "cANI": "0.98",
            "biosample": "SAMN1",
            "bioproject": "PRJNA1",
            "assay_type": "WGS",
            "organism": "Salmonella enterica",
            "geo_loc_name": "USA",
            "lat_lon": "35.7N 100.2W",
        }

    def test_default_details_rows_is_none_and_backward_compatible(self, tmp_path):
        """Existing callers that do not pass details_rows must keep working unchanged."""
        csv_file = tmp_path / "GCF_1.csv"
        csv_file.write_text("acc,containment\nSRR1,0.9\n")
        containment_data = defaultdict(dict)

        _process_genome_containments(csv_file, "GCF_1", containment_data)

        assert containment_data["SRR1"]["GCF_1"] == 0.9

    def test_missing_fields_become_empty_strings(self, tmp_path):
        csv_file = tmp_path / "GCF_1.csv"
        csv_file.write_text("acc,containment\nSRR1,0.9\n")
        containment_data = defaultdict(dict)
        details_rows: list = []

        _process_genome_containments(csv_file, "GCF_1", containment_data, details_rows=details_rows)

        row = details_rows[0]
        assert row["cANI"] == ""
        assert row["biosample"] == ""
        assert row["geo_loc_name"] == ""


class TestContainmentDetailsTable:
    """Test that parse_containment_data writes a details table alongside the parsed containment table."""

    def test_details_rows_for_two_genome_matches_folder(self, tmp_path):
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()
        (matches_dir / "GCF_1.csv").write_text(
            "acc,containment,cANI,biosample,bioproject,assay_type,organism,geo_loc_name_country_calc,lat_lon\n"
            "SRR1,0.95,0.99,SAMN1,PRJNA1,WGS,Salmonella enterica,USA,35.7N 100.2W\n"
        )
        (matches_dir / "GCF_2.csv").write_text(
            "acc,containment,cANI,biosample,bioproject,assay_type,organism,geo_loc_name_country_calc,lat_lon\n"
            "SRR2,0.80,0.95,SAMN2,PRJNA2,AMPLICON,gut metagenome,Denmark,55.67N 12.57E\n"
        )
        output_file = tmp_path / "parsed.txt"
        summary_file = tmp_path / "summary.txt"

        parse_containment_data(matches_dir, output_file, summary_file, 0.1)

        details_file = tmp_path / "parsed_details.tsv"
        assert details_file.exists()
        details = pd.read_csv(details_file, sep="\t", dtype=str).fillna("")
        assert list(details.columns) == [
            "accession",
            "genome_id",
            "containment",
            "cANI",
            "biosample",
            "bioproject",
            "assay_type",
            "organism",
            "geo_loc_name",
            "lat_lon",
        ]
        rows = {row["accession"]: row for _, row in details.iterrows()}
        assert set(rows) == {"SRR1", "SRR2"}
        assert rows["SRR1"]["genome_id"] == "GCF_1"
        assert rows["SRR1"]["cANI"] == "0.99"
        assert rows["SRR1"]["biosample"] == "SAMN1"
        assert rows["SRR1"]["geo_loc_name"] == "USA"
        assert rows["SRR2"]["genome_id"] == "GCF_2"
        assert rows["SRR2"]["assay_type"] == "AMPLICON"
        assert rows["SRR2"]["geo_loc_name"] == "Denmark"

    def test_details_file_override(self, tmp_path):
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()
        (matches_dir / "GCF_1.csv").write_text("acc,containment,cANI\nSRR1,0.9,0.98\n")
        output_file = tmp_path / "parsed.txt"
        summary_file = tmp_path / "summary.txt"
        custom_details = tmp_path / "custom_details.tsv"

        parse_containment_data(matches_dir, output_file, summary_file, 0.1, details_file=custom_details)

        assert custom_details.exists()
        assert not (tmp_path / "parsed_details.tsv").exists()

    def test_parsed_containment_unchanged_by_details_table(self, tmp_path):
        """Adding the details table must not change parsed_containment.txt's bytes."""
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()
        (matches_dir / "GCF_1.csv").write_text(
            "acc,containment,cANI,biosample\nSRR1,0.9,0.98,SAMN1\nSRR2,0.2,0.9,SAMN2\n"
        )
        output_file = tmp_path / "parsed.txt"
        summary_file = tmp_path / "summary.txt"

        parse_containment_data(matches_dir, output_file, summary_file, 0.1)
        produced = output_file.read_bytes()

        from prepare.data.file_io import write_csv

        containment_data = {"SRR1": {"GCF_1": 0.9}, "SRR2": {"GCF_1": 0.2}}
        df = pd.DataFrame.from_dict(containment_data, orient="index")
        df.fillna(0, inplace=True)
        df["max_containment"] = df.max(axis=1)
        df["max_containment_annotation"] = df.idxmax(axis=1)
        df.sort_values(by="max_containment", ascending=False, inplace=True)
        expected_file = tmp_path / "expected.txt"
        write_csv(df, expected_file, sep="\t")

        assert produced == expected_file.read_bytes()


def test_summary_counts_samples_at_exactly_one(tmp_path):
    """A sample whose containment exactly equals the threshold must be counted (inclusive rule)."""
    from prepare.data.branchwater import parse_containment_data

    matches = tmp_path / "matches"
    matches.mkdir()
    (matches / "GCF_1.csv").write_text("acc,containment,cANI\nSRR1,1.0,1.0\nSRR2,0.5,0.9\n")
    parse_containment_data(matches, tmp_path / "parsed.txt", tmp_path / "summary.txt", 0.5)
    rows = dict(line.split("\t") for line in (tmp_path / "summary.txt").read_text().splitlines()[1:])
    assert rows["1.0"] == "1" and rows["0.5"] == "2"


if __name__ == "__main__":
    pytest.main([__file__])


# --- Narrow exception handling (tech debt): a programming error must not be swallowed ---


def test_unexpected_error_while_parsing_one_containment_file_propagates(tmp_path):
    """Kind (a): a bug is no longer counted as an unreadable file."""
    (tmp_path / "genome1.csv").write_text("acc,containment\nSRR1,0.5\n")
    with patch("prepare.data.branchwater._process_genome_containments", side_effect=TypeError("bug")):
        with pytest.raises(TypeError):
            parse_containment_data(tmp_path, tmp_path / "out.tsv", tmp_path / "summary.tsv")


def test_unreadable_containment_file_is_an_error(tmp_path):
    """Kind (a): a file that cannot be read is named in the error that stops the run."""
    (tmp_path / "genome1.csv").write_text("wrong,header\nx,y\n")
    with pytest.raises(DataAccessError, match="genome1.csv"):
        parse_containment_data(tmp_path, tmp_path / "out.tsv", tmp_path / "summary.tsv")


def test_data_access_error_from_genome_containments_chains_the_cause(tmp_path):
    """Kind (c): the wrapped error keeps its cause."""
    with pytest.raises(DataAccessError) as exc:
        _process_genome_containments(tmp_path / "missing.csv", "g1", defaultdict(dict))
    assert isinstance(exc.value.__cause__, DataAccessError)  # read_csv reports the missing file


def test_data_access_error_from_metadata_extraction_chains_the_cause(tmp_path):
    """Kind (c): the wrapped error keeps its cause."""
    with pytest.raises(DataAccessError) as exc:
        extract_metadata_from_branchwater(tmp_path / "missing", tmp_path / "out.tsv")
    assert isinstance(exc.value.__cause__, DataAccessError)


class TestThresholdCounts:
    """The threshold grid uses exact decimals."""

    def test_maximum_on_a_threshold_is_counted(self):
        thresholds, counts = _threshold_counts(np.array([0.3, 0.6, 0.7]), 0.1)

        assert thresholds == [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
        assert counts == [0, 0, 0, 1, 2, 2, 2, 3, 3, 3, 3]

    def test_labels_are_half_up_decimals(self):
        thresholds, _ = _threshold_counts(np.array([0.5]), 0.005)

        assert thresholds[200 - 191] == 0.96
        assert thresholds[200 - 139] == 0.70
        thresholds, _ = _threshold_counts(np.array([0.5]), 0.1)
        assert thresholds == [round(i / 10, 1) for i in range(10, -1, -1)]

    def test_step_that_does_not_divide_one_uses_the_nearest_whole_number_of_steps(self):
        thresholds, counts = _threshold_counts(np.array([0.5]), 0.007)

        assert len(thresholds) == len(counts) == 144
        assert thresholds[0] == 1.0
        assert thresholds[-1] == 0.0
