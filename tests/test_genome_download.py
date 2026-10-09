"""Tests for the genome download module."""

import subprocess
import zipfile
from unittest.mock import MagicMock, patch

import pytest

from prepare.core.exceptions import DataAccessError
from prepare.data.genome_download import (
    _validate_genome_accession,
    download_genomes,
    extract_and_organize,
)

# --- _validate_genome_accession ---


class TestValidateGenomeAccession:
    def test_valid_gcf_accession(self):
        assert _validate_genome_accession("GCF_000005845.2") == "GCF_000005845.2"

    def test_valid_gca_accession(self):
        assert _validate_genome_accession("GCA_000005845.2") == "GCA_000005845.2"

    def test_invalid_prefix(self):
        with pytest.raises(DataAccessError, match="Invalid genome accession prefix"):
            _validate_genome_accession("SRR1234567")

    def test_invalid_format_no_version(self):
        with pytest.raises(DataAccessError, match="Invalid genome accession format"):
            _validate_genome_accession("GCF_000005845")

    def test_invalid_format_short_digits(self):
        with pytest.raises(DataAccessError, match="Invalid genome accession format"):
            _validate_genome_accession("GCF_12345.1")

    def test_empty_string(self):
        with pytest.raises(DataAccessError, match="Invalid genome accession prefix"):
            _validate_genome_accession("")


# --- download_genomes ---


class TestDownloadGenomes:
    @patch("prepare.data.genome_download.SecureSubprocess.run_secure")
    def test_single_accession(self, mock_run, tmp_path):
        # Create the zip file that would be produced by the command
        zip_path = tmp_path / "ncbi_dataset.zip"
        zip_path.touch()

        mock_run.return_value = MagicMock(returncode=0)

        result = download_genomes(["GCF_000005845.2"], tmp_path)

        assert result == zip_path
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0]
        assert call_args[0] == "datasets"
        assert "GCF_000005845.2" in call_args[1]
        assert "--include" in call_args[1]
        assert "genome" in call_args[1]

    @patch("prepare.data.genome_download.SecureSubprocess.run_secure")
    def test_multiple_accessions(self, mock_run, tmp_path):
        zip_path = tmp_path / "ncbi_dataset.zip"
        zip_path.touch()

        mock_run.return_value = MagicMock(returncode=0)
        accessions = ["GCF_000005845.2", "GCA_000001405.1"]

        result = download_genomes(accessions, tmp_path)

        assert result == zip_path
        call_args = mock_run.call_args[0][1]
        assert "GCF_000005845.2" in call_args
        assert "GCA_000001405.1" in call_args

    @patch("prepare.data.genome_download.SecureSubprocess.run_secure")
    def test_with_assembly_level(self, mock_run, tmp_path):
        zip_path = tmp_path / "ncbi_dataset.zip"
        zip_path.touch()

        mock_run.return_value = MagicMock(returncode=0)

        download_genomes(["GCF_000005845.2"], tmp_path, assembly_level="complete")

        call_args = mock_run.call_args[0][1]
        assert "--assembly-level" in call_args
        assert "complete" in call_args

    def test_empty_accessions(self, tmp_path):
        with pytest.raises(DataAccessError, match="No accessions provided"):
            download_genomes([], tmp_path)

    def test_invalid_accession(self, tmp_path):
        with pytest.raises(DataAccessError, match="Invalid genome accession"):
            download_genomes(["INVALID_123"], tmp_path)

    @patch("prepare.data.genome_download.SecureSubprocess.run_secure")
    def test_download_failure(self, mock_run, tmp_path):
        mock_run.side_effect = subprocess.CalledProcessError(1, "datasets", stderr="Connection failed")
        with pytest.raises(DataAccessError, match="Failed to download genomes"):
            download_genomes(["GCF_000005845.2"], tmp_path)

    @patch("prepare.data.genome_download.SecureSubprocess.run_secure")
    def test_zip_not_created(self, mock_run, tmp_path):
        # Simulate command success but no zip produced
        mock_run.return_value = MagicMock(returncode=0)

        with pytest.raises(DataAccessError, match="zip file was not created"):
            download_genomes(["GCF_000005845.2"], tmp_path)


# --- extract_and_organize ---


def _create_mock_ncbi_zip(zip_path, accessions_and_files):
    """Helper to create a mock NCBI datasets zip file.

    Args:
        zip_path: Path for the zip file.
        accessions_and_files: Dict mapping accession to FASTA content.
    """
    with zipfile.ZipFile(zip_path, "w") as zf:
        for accession, content in accessions_and_files.items():
            fasta_name = f"ncbi_dataset/data/{accession}/{accession}_genomic.fna"
            zf.writestr(fasta_name, content)


class TestExtractAndOrganize:
    def test_single_genome(self, tmp_path):
        zip_path = tmp_path / "test.zip"
        _create_mock_ncbi_zip(zip_path, {"GCF_000005845.2": ">seq1\nATCG\n"})

        output_dir = tmp_path / "genomes"
        result = extract_and_organize(zip_path, output_dir)

        assert "GCF_000005845.2" in result
        fasta = result["GCF_000005845.2"]
        assert fasta.exists()
        assert fasta.name == "GCF_000005845.2.fna"
        assert ">seq1" in fasta.read_text()

    def test_multiple_genomes(self, tmp_path):
        zip_path = tmp_path / "test.zip"
        _create_mock_ncbi_zip(
            zip_path,
            {
                "GCF_000005845.2": ">ecoli\nATCG\n",
                "GCA_000001405.1": ">human\nGGCC\n",
            },
        )

        output_dir = tmp_path / "genomes"
        result = extract_and_organize(zip_path, output_dir)

        assert len(result) == 2
        assert "GCF_000005845.2" in result
        assert "GCA_000001405.1" in result

    def test_skips_non_accession_dirs(self, tmp_path):
        zip_path = tmp_path / "test.zip"
        with zipfile.ZipFile(zip_path, "w") as zf:
            zf.writestr(
                "ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_genomic.fna",
                ">seq\nATCG\n",
            )
            zf.writestr(
                "ncbi_dataset/data/assembly_data_report.jsonl",
                '{"report": true}\n',
            )

        output_dir = tmp_path / "genomes"
        result = extract_and_organize(zip_path, output_dir)

        assert len(result) == 1
        assert "GCF_000005845.2" in result

    def test_missing_zip(self, tmp_path):
        with pytest.raises(DataAccessError, match="Zip file not found"):
            extract_and_organize(tmp_path / "missing.zip", tmp_path / "out")

    def test_invalid_zip(self, tmp_path):
        bad_zip = tmp_path / "bad.zip"
        bad_zip.write_text("not a zip file")

        with pytest.raises(DataAccessError, match="Not a valid zip file"):
            extract_and_organize(bad_zip, tmp_path / "out")

    def test_corrupt_zip(self, tmp_path):
        bad_zip = tmp_path / "corrupt.zip"
        # Write something that looks like a zip header but is corrupt
        bad_zip.write_bytes(b"PK\x03\x04" + b"\x00" * 100)

        with pytest.raises(DataAccessError):
            extract_and_organize(bad_zip, tmp_path / "out")

    def test_unexpected_structure(self, tmp_path):
        zip_path = tmp_path / "test.zip"
        with zipfile.ZipFile(zip_path, "w") as zf:
            zf.writestr("some_other_dir/file.txt", "data")

        with pytest.raises(DataAccessError, match="Unexpected zip structure"):
            extract_and_organize(zip_path, tmp_path / "out")

    def test_no_fasta_in_accession_dir(self, tmp_path):
        zip_path = tmp_path / "test.zip"
        with zipfile.ZipFile(zip_path, "w") as zf:
            zf.writestr(
                "ncbi_dataset/data/GCF_000005845.2/report.json",
                '{"report": true}',
            )

        output_dir = tmp_path / "genomes"
        result = extract_and_organize(zip_path, output_dir)

        assert len(result) == 0

    def test_cleanup_on_success(self, tmp_path):
        """Verify temp extraction dir is cleaned up after success."""
        zip_path = tmp_path / "test.zip"
        _create_mock_ncbi_zip(zip_path, {"GCF_000005845.2": ">seq\nATCG\n"})

        extract_and_organize(zip_path, tmp_path / "genomes")

        extract_dir = tmp_path / "ncbi_extract_tmp"
        assert not extract_dir.exists()
