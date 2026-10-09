"""
Test processing.containment module functionality.

Tests for containment analysis including test genome download
and containment data analysis functions.
"""

import io
import gzip
from pathlib import Path
from unittest.mock import Mock, patch

import pytest
import pandas as pd

from prepare.processing.containment import (
    download_test_genome,
    count_single_sample,
)
from prepare.core.exceptions import ProcessingError


class TestDownloadTestGenome:
    """Test download_test_genome function."""

    def test_download_test_genome_file_exists(self, tmp_path):
        """Test when test genome file already exists."""
        # Create existing file
        output_path = tmp_path / "GCF_000008985.1.fasta"
        output_path.write_text("existing genome content")

        result = download_test_genome(tmp_path)

        assert result == output_path
        assert output_path.exists()
        assert output_path.read_text() == "existing genome content"

    @patch("prepare.processing.containment.urllib.request.urlopen")
    def test_download_test_genome_success(self, mock_urlopen, tmp_path):
        """Test successful genome download."""
        # Create temporary compressed file content
        genome_content = ">test_genome\nACGTACGTACGT\n"
        compressed = io.BytesIO()
        with gzip.open(compressed, "wt") as f:
            f.write(genome_content)

        # Mock urlopen to return our compressed data
        mock_response = Mock()
        mock_response.read.return_value = compressed.getvalue()
        mock_urlopen.return_value = mock_response

        result = download_test_genome(tmp_path)

        expected_path = tmp_path / "GCF_000008985.1.fasta"
        assert result == expected_path
        assert expected_path.exists()
        assert expected_path.read_text() == genome_content

        # Verify temporary gz file was cleaned up
        temp_gz_path = expected_path.with_suffix(".gz")
        assert not temp_gz_path.exists()

    @patch("prepare.processing.containment.urllib.request.urlopen")
    def test_download_test_genome_download_failure(self, mock_urlopen, tmp_path):
        """Test failure during download."""
        mock_urlopen.side_effect = Exception("Download failed")

        with pytest.raises(ProcessingError, match="Error downloading test genome"):
            download_test_genome(tmp_path)

    @patch("prepare.processing.containment.urllib.request.urlopen")
    def test_download_test_genome_decompression_failure(self, mock_urlopen, tmp_path):
        """Test failure during decompression."""
        # Mock urlopen to return invalid gzip data
        mock_response = Mock()
        mock_response.read.return_value = b"not a valid gzip file"
        mock_urlopen.return_value = mock_response

        with pytest.raises(ProcessingError, match="Error downloading test genome"):
            download_test_genome(tmp_path)

    def test_download_test_genome_invalid_output_folder(self):
        """Test with invalid output folder."""
        # Test with path that cannot be created (e.g., permission denied scenario)
        with patch("prepare.processing.containment.ensure_directory") as mock_ensure:
            mock_ensure.side_effect = PermissionError("Permission denied")

            with pytest.raises(ProcessingError):
                download_test_genome("/invalid/path")

    @patch("prepare.processing.containment.urllib.request.urlretrieve")
    def test_download_test_genome_creates_directory(self, mock_urlretrieve, tmp_path):
        """Test that output directory is created if it doesn't exist."""
        # Use a subdirectory that doesn't exist yet
        output_dir = tmp_path / "subdir" / "genomes"

        # Create temporary compressed file content
        genome_content = ">test_genome\nACGTACGTACGT\n"
        temp_gz_file = tmp_path / "temp.gz"

        with gzip.open(temp_gz_file, "wt") as f:
            f.write(genome_content)

        def mock_download(url, filename):
            temp_gz_path = Path(filename)
            temp_gz_path.parent.mkdir(parents=True, exist_ok=True)
            temp_gz_path.write_bytes(temp_gz_file.read_bytes())

        mock_urlretrieve.side_effect = mock_download

        result = download_test_genome(output_dir)

        expected_path = output_dir / "GCF_000008985.1.fasta"
        assert result == expected_path
        assert expected_path.exists()
        assert output_dir.exists()

    @patch("prepare.processing.containment.urllib.request.urlretrieve")
    @patch("prepare.processing.containment.logger")
    def test_download_test_genome_logging(self, mock_logger, mock_urlretrieve, tmp_path):
        """Test that appropriate logging messages are generated."""
        # Test logging when file exists
        output_path = tmp_path / "GCF_000008985.1.fasta"
        output_path.write_text("existing")

        download_test_genome(tmp_path)
        mock_logger.info.assert_called_with(f"Test genome already exists at {output_path}")

        # Reset mock and test logging during download
        mock_logger.reset_mock()
        output_path.unlink()  # Remove existing file

        # Create temporary compressed file
        genome_content = ">test_genome\nACGT\n"
        temp_gz_file = tmp_path / "temp.gz"

        with gzip.open(temp_gz_file, "wt") as f:
            f.write(genome_content)

        def mock_download(url, filename):
            Path(filename).write_bytes(temp_gz_file.read_bytes())

        mock_urlretrieve.side_effect = mock_download

        download_test_genome(tmp_path)

        # Check that download-related logging occurred
        mock_logger.info.assert_any_call("Downloading test genome")
        mock_logger.info.assert_any_call(f"Downloaded test genome to {output_path}")

    def test_download_test_genome_path_handling(self, tmp_path):
        """Test that function handles both string and Path inputs."""
        # Test with string path
        output_path_str = str(tmp_path)

        with patch("prepare.processing.containment.urllib.request.urlretrieve") as mock_urlretrieve:
            genome_content = ">test\nACGT\n"
            temp_gz_file = tmp_path / "temp.gz"

            with gzip.open(temp_gz_file, "wt") as f:
                f.write(genome_content)

            def mock_download(url, filename):
                Path(filename).write_bytes(temp_gz_file.read_bytes())

            mock_urlretrieve.side_effect = mock_download

            result = download_test_genome(output_path_str)

            expected_path = tmp_path / "GCF_000008985.1.fasta"
            assert result == expected_path
            assert isinstance(result, Path)


class TestCountSingleSample:
    """Test count_single_sample function."""

    def setup_method(self):
        """Set up test data."""
        # Sample summary data
        self.summary_data = pd.DataFrame(
            {
                "GCF_000001.1": [0.8, 0.2, 0.6, 0.05, 0.9],
                "GCF_000002.1": [0.1, 0.7, 0.3, 0.15, 0.2],
                "max_containment": [0.8, 0.7, 0.6, 0.15, 0.9],
            },
            index=["SRR001", "SRR002", "SRR003", "SRR004", "SRR005"],
        )

        # Sample metadata
        self.metadata_data = pd.DataFrame(
            {
                "organism": ["E. coli", "S. aureus", "E. coli", "B. subtilis", "E. coli"],
                "country": ["USA", "UK", "Canada", "Germany", "USA"],
                "year": [2020, 2019, 2021, 2020, 2022],
            },
            index=["SRR001", "SRR002", "SRR003", "SRR004", "SRR005"],
        )

    def test_count_single_sample_success(self, tmp_path):
        """Test successful metadata counting."""
        summary_file = tmp_path / "summary.tsv"
        metadata_file = tmp_path / "metadata.tsv"

        self.summary_data.to_csv(summary_file, sep="\t")
        self.metadata_data.to_csv(metadata_file, sep="\t")

        result = count_single_sample(summary_file, metadata_file, "GCF_000001.1", "organism", threshold=0.5, top_n=10)

        # Expected: SRR001 (0.8), SRR003 (0.6), SRR005 (0.9) pass threshold
        # All have 'E. coli' in metadata
        assert result == {"E. coli": 3}

    def test_count_single_sample_multiple_organisms(self, tmp_path):
        """Test counting with multiple different organisms."""
        summary_file = tmp_path / "summary.tsv"
        metadata_file = tmp_path / "metadata.tsv"

        # Modify metadata to have different organisms for samples above threshold
        metadata_mixed = self.metadata_data.copy()
        metadata_mixed.loc["SRR001", "organism"] = "E. coli"
        metadata_mixed.loc["SRR003", "organism"] = "B. subtilis"
        metadata_mixed.loc["SRR005", "organism"] = "E. coli"

        self.summary_data.to_csv(summary_file, sep="\t")
        metadata_mixed.to_csv(metadata_file, sep="\t")

        result = count_single_sample(summary_file, metadata_file, "GCF_000001.1", "organism", threshold=0.5, top_n=10)

        assert result == {"E. coli": 2, "B. subtilis": 1}

    def test_count_single_sample_no_matches_above_threshold(self, tmp_path):
        """Test behavior when no samples above threshold."""
        summary_file = tmp_path / "summary.tsv"
        metadata_file = tmp_path / "metadata.tsv"

        self.summary_data.to_csv(summary_file, sep="\t")
        self.metadata_data.to_csv(metadata_file, sep="\t")

        result = count_single_sample(summary_file, metadata_file, "GCF_000001.1", "organism", threshold=0.95, top_n=10)

        assert result == {}

    def test_count_single_sample_missing_summary_column(self, tmp_path):
        """Test error when summary column doesn't exist."""
        summary_file = tmp_path / "summary.tsv"
        metadata_file = tmp_path / "metadata.tsv"

        self.summary_data.to_csv(summary_file, sep="\t")
        self.metadata_data.to_csv(metadata_file, sep="\t")

        with pytest.raises(ProcessingError, match="Column MISSING not found in summary file"):
            count_single_sample(summary_file, metadata_file, "MISSING", "organism")

    def test_count_single_sample_missing_metadata_column(self, tmp_path):
        """Test error when metadata column doesn't exist."""
        summary_file = tmp_path / "summary.tsv"
        metadata_file = tmp_path / "metadata.tsv"

        self.summary_data.to_csv(summary_file, sep="\t")
        self.metadata_data.to_csv(metadata_file, sep="\t")

        with pytest.raises(ProcessingError, match="Column missing not found in metadata file"):
            count_single_sample(summary_file, metadata_file, "GCF_000001.1", "missing")

    def test_count_single_sample_no_metadata_matches(self, tmp_path):
        """Test behavior when no metadata matches selected accessions."""
        summary_file = tmp_path / "summary.tsv"
        metadata_file = tmp_path / "metadata.tsv"

        # Create metadata with different index
        different_metadata = pd.DataFrame({"organism": ["X. test"]}, index=["SRR999"])

        self.summary_data.to_csv(summary_file, sep="\t")
        different_metadata.to_csv(metadata_file, sep="\t")

        result = count_single_sample(summary_file, metadata_file, "GCF_000001.1", "organism", threshold=0.5)

        assert result == {}

    def test_count_single_sample_top_n_limit(self, tmp_path):
        """Test top_n parameter limits results."""
        summary_file = tmp_path / "summary.tsv"
        metadata_file = tmp_path / "metadata.tsv"

        # Create data with many different values
        large_metadata = pd.DataFrame(
            {"organism": [f"Species_{i}" for i in range(20)]}, index=[f"SRR{i:03d}" for i in range(1, 21)]
        )

        large_summary = pd.DataFrame(
            {"GCF_000001.1": [0.8] * 20}, index=[f"SRR{i:03d}" for i in range(1, 21)]  # All above threshold
        )

        large_summary.to_csv(summary_file, sep="\t")
        large_metadata.to_csv(metadata_file, sep="\t")

        result = count_single_sample(summary_file, metadata_file, "GCF_000001.1", "organism", threshold=0.5, top_n=5)

        assert len(result) == 5

    def test_count_single_sample_file_read_error(self, tmp_path):
        """Test handling of file reading errors."""
        summary_file = tmp_path / "missing_summary.tsv"
        metadata_file = tmp_path / "metadata.tsv"

        self.metadata_data.to_csv(metadata_file, sep="\t")

        with pytest.raises(ProcessingError, match="Error counting single sample metadata"):
            count_single_sample(summary_file, metadata_file, "GCF_000001.1", "organism")

    def test_count_single_sample_string_paths(self, tmp_path):
        """Test function accepts string paths."""
        summary_file = tmp_path / "summary.tsv"
        metadata_file = tmp_path / "metadata.tsv"

        self.summary_data.to_csv(summary_file, sep="\t")
        self.metadata_data.to_csv(metadata_file, sep="\t")

        result = count_single_sample(str(summary_file), str(metadata_file), "GCF_000001.1", "organism", threshold=0.5)

        assert result == {"E. coli": 3}
