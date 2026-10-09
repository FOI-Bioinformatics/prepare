"""
Test CLI enhanced SRA commands functionality.

Tests for the sra_info and sra_validate command classes, focusing on argument parsing,
validation, and proper delegation with mocked dependencies. The statistics that sra_stats
computed before 0.5.0 are sra_profile's, tested in tests/test_cli_sra_profile.py.
"""

import argparse
import json
from unittest.mock import Mock, patch, mock_open

import pytest

from prepare.cli.commands.sra_enhanced import SRAInfoCommand, SRAValidateCommand


class TestSRAInfoCommand:
    """Test SRAInfoCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = SRAInfoCommand()
        assert command.name == "sra_info"
        assert "information" in command.help.lower()
        assert "sra" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = SRAInfoCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test required arguments
        args = parser.parse_args(["--accessions-file", "test.txt", "--email", "test@example.com"])
        assert args.accessions_file == "test.txt"
        assert args.email == "test@example.com"

        # Test optional arguments with defaults
        assert args.output_report == "sra_info_report.csv"
        assert args.bandwidth_mbps == 100.0

    def test_configure_parser_with_optional_args(self):
        """Test parser with optional arguments."""
        command = SRAInfoCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(
            [
                "--accessions-file",
                "test.txt",
                "--email",
                "test@example.com",
                "--api-key",
                "mykey123",
                "--output-report",
                "custom_report.csv",
                "--bandwidth-mbps",
                "250.5",
            ]
        )

        assert args.api_key == "mykey123"
        assert args.output_report == "custom_report.csv"
        assert args.bandwidth_mbps == 250.5

    @patch("prepare.cli.commands.sra_enhanced.save_metadata_report")
    @patch("prepare.cli.commands.sra_enhanced.create_download_preview")
    @patch("prepare.cli.commands.sra_enhanced.SRAMetadataClient")
    def test_execute_success(self, mock_client_class, mock_create_preview, mock_save_report, tmp_path, capsys):
        """Test successful execution."""
        command = SRAInfoCommand()
        args = argparse.Namespace(
            accessions_file="test_accessions.txt",
            email="test@example.com",
            api_key="test_key",
            output_report="test_report.csv",
            bandwidth_mbps=100.0,
        )

        mock_client = Mock()
        mock_client_class.return_value = mock_client

        mock_metadata = {
            "SRR123456": Mock(platform="ILLUMINA", layout="PAIRED", size_mb=1024),
            "SRR789012": Mock(platform="ILLUMINA", layout="SINGLE", size_mb=512),
        }
        mock_tech_counts = {"RNA-Seq": 2}
        mock_total_size = 1.5

        mock_create_preview.return_value = (mock_metadata, mock_tech_counts, mock_total_size)

        # Mock file reading
        mock_file_content = "SRR123456\nSRR789012\n"
        with patch("builtins.open", mock_open(read_data=mock_file_content)):
            result = command.execute(args)

        assert result == 0
        mock_client.close.assert_called_once()
        out = capsys.readouterr().out
        assert "Total accessions: 2" in out
        assert "Metadata fetched: 2" in out
        mock_client_class.assert_called_once_with("test@example.com", "test_key")
        mock_create_preview.assert_called_once_with(["SRR123456", "SRR789012"], mock_client)
        mock_save_report.assert_called_once()

    def test_execute_no_accessions(self, caplog):
        """Test execution with empty accessions file."""
        command = SRAInfoCommand()
        args = argparse.Namespace(
            accessions_file="empty.txt",
            email="test@example.com",
            api_key=None,
            output_report="report.csv",
            bandwidth_mbps=100.0,
        )

        with patch("builtins.open", mock_open(read_data="")):
            result = command.execute(args)

        assert result == 1
        assert "No accessions found in file" in caplog.text

    @patch("prepare.cli.commands.sra_enhanced.create_download_preview")
    @patch("prepare.cli.commands.sra_enhanced.SRAMetadataClient")
    def test_execute_no_metadata(self, mock_client_class, mock_create_preview, caplog):
        """Test execution when no metadata can be fetched."""
        command = SRAInfoCommand()
        args = argparse.Namespace(
            accessions_file="test.txt",
            email="test@example.com",
            api_key=None,
            output_report="report.csv",
            bandwidth_mbps=100.0,
        )

        mock_create_preview.return_value = ({}, {}, 0.0)

        with patch("builtins.open", mock_open(read_data="SRR123456\n")):
            result = command.execute(args)

        assert result == 1
        assert "Could not fetch metadata for any accessions" in caplog.text

    @patch("builtins.open", side_effect=FileNotFoundError())
    def test_execute_file_error(self, mock_open, caplog):
        """Test execution with file reading error: a DataAccessError reaches self.fail (exit 1)."""
        command = SRAInfoCommand()
        args = argparse.Namespace(
            accessions_file="nonexistent.txt",
            email="test@example.com",
            api_key=None,
            output_report="report.csv",
            bandwidth_mbps=100.0,
        )

        result = command.execute(args)

        assert result == 1
        assert "SRA info command failed" in caplog.text

    @patch("prepare.cli.commands.sra_enhanced.create_download_preview")
    @patch("prepare.cli.commands.sra_enhanced.SRAMetadataClient")
    def test_execute_network_error_exits_4(self, mock_client_class, mock_create_preview):
        """A NetworkError from the metadata client reaches ``self.fail`` with exit code 4, not 1.

        Regression test for the broad ``except Exception`` this command used to have, which
        turned a retryable NCBI failure (``NetworkError``, ``prepare.core.exceptions``) into
        the same exit code 1 as any other failure.
        """
        from prepare.core.exceptions import NetworkError

        command = SRAInfoCommand()
        args = argparse.Namespace(
            accessions_file="test.txt",
            email="test@example.com",
            api_key=None,
            output_report="report.csv",
            bandwidth_mbps=100.0,
        )
        mock_create_preview.side_effect = NetworkError("Failed to query NCBI: too many retries")

        with patch("builtins.open", mock_open(read_data="SRR123456\n")):
            result = command.execute(args)

        assert result == 4

    def test_execute_network_error_exits_4_end_to_end(self, monkeypatch):
        """The same exit-4 path, exercised end to end against a real (faked) HTTP layer.

        Unlike ``test_execute_network_error_exits_4`` above, ``SRAMetadataClient`` and
        ``create_download_preview`` are not mocked: ``tests/fake_http.py`` fakes the connection
        pool a persistently-503 NCBI would present, so the real retrying session
        (``prepare.utils.http.retrying_session``), ``SRAMetadataClient._make_request`` and
        ``SRAMetadataClient.get_sra_metadata``'s batch loop are all exercised, confirming a
        ``NetworkError`` actually propagates out of ``get_sra_metadata`` (it is a
        ``DataAccessError`` subclass, which that loop otherwise catches and continues past) all
        the way to ``self.fail``.
        """
        from tests.fake_http import serve

        serve(monkeypatch, lambda path: (503, b"unavailable"))

        command = SRAInfoCommand()
        args = argparse.Namespace(
            accessions_file="test.txt",
            email="test@example.com",
            api_key=None,
            output_report="report.csv",
            bandwidth_mbps=100.0,
        )

        with patch("builtins.open", mock_open(read_data="SRR123456\n")):
            result = command.execute(args)

        assert result == 4


class TestSRAValidateCommand:
    """Test SRAValidateCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = SRAValidateCommand()
        assert command.name == "sra_validate"
        assert "validate" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = SRAValidateCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args([])
        assert args.fastq_folder == "fastq"
        assert args.accession is None and args.accessions_file is None
        assert not args.check_pairs
        assert not args.md5

    def test_configure_parser_with_options(self):
        """Test parser with all options."""
        command = SRAValidateCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(
            [
                "--fastq-folder",
                "custom_fastq",
                "--accession",
                "SRR123",
                "--accession",
                "SRR456",
                "--accessions-file",
                "acc.txt",
                "--check-pairs",
                "--md5",
            ]
        )

        assert args.fastq_folder == "custom_fastq"
        assert args.accession == ["SRR123", "SRR456"]
        assert args.accessions_file == "acc.txt"
        assert args.check_pairs
        assert args.md5

    def test_find_accession_dirs_all(self, tmp_path):
        """Test finding all accession directories."""
        command = SRAValidateCommand()

        # Create test directories
        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        (fastq_folder / "SRR123").mkdir()
        (fastq_folder / "SRR456").mkdir()
        (fastq_folder / "not_accession.txt").touch()  # Should be ignored

        result = command._find_accession_dirs(fastq_folder)

        dir_names = {d.name for d in result}
        assert dir_names == {"SRR123", "SRR456"}

    def test_find_accession_dirs_specific(self, tmp_path):
        """Test finding specific accession directories."""
        command = SRAValidateCommand()

        # Create test directories
        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        (fastq_folder / "SRR123").mkdir()
        (fastq_folder / "SRR456").mkdir()
        (fastq_folder / "SRR789").mkdir()

        result = command._find_accession_dirs(fastq_folder, ["SRR123", "SRR789"])

        dir_names = {d.name for d in result}
        assert dir_names == {"SRR123", "SRR789"}

    def test_find_accession_dirs_skips_hidden_dirs(self, tmp_path):
        """A hidden folder (e.g. ``.Trashes`` or ``._SRR123``) is not an accession."""
        command = SRAValidateCommand()
        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        (fastq_folder / "SRR123").mkdir()
        (fastq_folder / "._SRR123").mkdir()
        (fastq_folder / ".Trashes").mkdir()

        result = command._find_accession_dirs(fastq_folder)

        assert [d.name for d in result] == ["SRR123"]

    def test_validate_directory_ignores_appledouble_files(self, tmp_path):
        """A ``._<name>.fastq.gz`` AppleDouble file is not counted as a FASTQ file."""
        command = SRAValidateCommand()
        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()
        (acc_dir / "test.fastq").write_text("@read1\nACGT\n+\n!!!!\n")
        (acc_dir / "._test.fastq").write_bytes(b"\x00" * 4096)

        result = command._validate_directory(acc_dir)

        assert result["num_files"] == 1

    def test_validate_directory_no_files(self, tmp_path):
        """Test validation of directory with no FASTQ files."""
        command = SRAValidateCommand()

        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()

        result = command._validate_directory(acc_dir)

        assert result["accession"] == "SRR123"
        assert result["status"] == "FAILED"
        assert "No FASTQ files found" in result["issues"]
        assert result["num_files"] == 0

    def test_validate_directory_success(self, tmp_path):
        """Test successful directory validation."""
        command = SRAValidateCommand()

        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()

        # Create test FASTQ file with valid content
        fastq_file = acc_dir / "test.fastq"
        fastq_file.write_text("@read1\nACGT\n+\n!!!!\n")

        result = command._validate_directory(acc_dir)

        assert result["accession"] == "SRR123"
        assert result["status"] == "PASSED"
        assert result["issues"] == "None"
        assert result["issues_list"] == []
        assert result["num_files"] == 1
        assert result["checks"] == ["empty_files", "format", "completeness"]

    def test_validate_directory_empty_files(self, tmp_path):
        """Test validation with empty files."""
        command = SRAValidateCommand()

        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()

        # Create empty FASTQ file
        fastq_file = acc_dir / "test.fastq"
        fastq_file.touch()

        result = command._validate_directory(acc_dir)

        assert result["status"] == "FAILED"
        assert "Empty file: test.fastq" in result["issues"]

    def test_validate_directory_mate_count_mismatch(self, tmp_path):
        """A paired-end dataset whose mates have different read counts is flagged, with the
        read counts named in the message."""
        command = SRAValidateCommand()

        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()

        (acc_dir / "SRR123_1.fastq").write_text("@r\nACGT\n+\nIIII\n" * 2)
        (acc_dir / "SRR123_2.fastq").write_text("@r\nACGT\n+\nIIII\n" * 1)

        result = command._validate_directory(acc_dir, check_pairs=True)

        assert result["status"] == "FAILED"
        assert "mate files differ (2 vs 1)" in result["issues"]
        assert "mate_counts" in result["checks"]

    @patch("prepare.cli.commands.sra_enhanced.cached_stats")
    @patch("prepare.cli.commands.sra_enhanced.count_fastq_reads")
    def test_validate_directory_mate_count_uses_cached_stats(self, mock_count, mock_cached, tmp_path):
        """When a cached stats record is available, mate counts come from its
        ``reads_per_file`` rather than a fresh ``count_fastq_reads`` pass."""
        mock_cached.return_value = {"reads_per_file": {"SRR123_1.fastq": 5, "SRR123_2.fastq": 5}}

        command = SRAValidateCommand()
        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()
        (acc_dir / "SRR123_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (acc_dir / "SRR123_2.fastq").write_text("@r\nACGT\n+\nIIII\n")

        result = command._validate_directory(acc_dir, check_pairs=True)

        assert result["status"] == "PASSED"
        mock_count.assert_not_called()

    def test_validate_directory_format_error_broken_header(self, tmp_path):
        """A file whose first line is not a FASTQ header is caught."""
        command = SRAValidateCommand()

        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()

        fastq_file = acc_dir / "test.fastq"
        fastq_file.write_text("invalid fastq content\nACGT\n+\nIIII\n")

        result = command._validate_directory(acc_dir)

        assert result["status"] == "FAILED"
        assert "FASTQ format error" in result["issues"]
        assert "header does not start with" in result["issues"]

    def test_validate_directory_checks_every_file_not_just_the_first(self, tmp_path):
        """A download can leave one good mate and one broken one; both are checked."""
        command = SRAValidateCommand()

        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()
        (acc_dir / "SRR123_1.fastq").write_text("@read1\nACGT\n+\nIIII\n")
        (acc_dir / "SRR123_2.fastq").write_text("not-a-header\nACGT\n+\nIIII\n")

        result = command._validate_directory(acc_dir)

        assert result["status"] == "FAILED"
        assert "SRR123_2.fastq" in result["issues"]
        assert "SRR123_1.fastq" not in result["issues"]

    def test_validate_directory_reads_a_gz_first_record(self, tmp_path):
        """The first-record check is gzip aware in both directions: a valid gzipped file
        passes and a broken one is caught, without decompressing the whole file."""
        import gzip as gzip_module

        command = SRAValidateCommand()

        good = tmp_path / "SRR1"
        good.mkdir()
        with gzip_module.open(good / "SRR1.fastq.gz", "wt") as handle:
            handle.write("@read1\nACGT\n+\nIIII\n")
        assert command._validate_directory(good)["status"] == "PASSED"

        bad = tmp_path / "SRR2"
        bad.mkdir()
        with gzip_module.open(bad / "SRR2.fastq.gz", "wt") as handle:
            handle.write("@read1\nACGTACGT\n+\nIII\n")
        result = command._validate_directory(bad)
        assert result["status"] == "FAILED"
        assert "sequence/quality length mismatch" in result["issues"]

    def test_validate_directory_md5_without_a_sidecar_is_a_no_op(self, tmp_path):
        """A plain project folder has no recorded md5 to compare against, so --md5 passes
        rather than failing every file."""
        command = SRAValidateCommand()

        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()
        (acc_dir / "SRR123.fastq").write_text("@read1\nACGT\n+\nIIII\n")

        result = command._validate_directory(acc_dir, check_md5=True)

        assert result["status"] == "PASSED"
        assert "md5" in result["checks"]

    def test_validate_directory_reads_the_statistics_record_only_for_check_pairs(self, tmp_path):
        """Without --check-pairs nothing needs the record, so the sidecar is not read."""
        command = SRAValidateCommand()

        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()
        (acc_dir / "SRR123_1.fastq").write_text("@read1\nACGT\n+\nIIII\n")
        (acc_dir / "SRR123_2.fastq").write_text("@read1\nACGT\n+\nIIII\n")

        with patch("prepare.cli.commands.sra_enhanced.cached_stats") as mock_cached:
            command._validate_directory(acc_dir)
            mock_cached.assert_not_called()

            command._validate_directory(acc_dir, check_pairs=True)
            mock_cached.assert_called_once()

    def test_validate_directory_format_error_length_mismatch(self, tmp_path):
        """A first record whose sequence and quality strings differ in length is caught."""
        command = SRAValidateCommand()

        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()

        fastq_file = acc_dir / "test.fastq"
        fastq_file.write_text("@read1\nACGTACGT\n+\nIII\n")

        result = command._validate_directory(acc_dir)

        assert result["status"] == "FAILED"
        assert "sequence/quality length mismatch" in result["issues"]

    def test_validate_directory_format_check_does_not_read_a_large_file_in_full(self, tmp_path):
        """The first-record shape check never falls back to a full-file read: a large file
        with a broken header is rejected without ``count_fastq_reads`` (a full streaming
        pass) ever being called."""
        command = SRAValidateCommand()

        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()

        fastq_file = acc_dir / "test.fastq"
        # A large file (many valid-looking records) but a broken first header.
        fastq_file.write_text("not-a-header\n" + ("@r\nACGT\n+\nIIII\n" * 50000))

        with patch("prepare.cli.commands.sra_enhanced.count_fastq_reads") as mock_count:
            result = command._validate_directory(acc_dir)

        assert result["status"] == "FAILED"
        assert "header does not start with" in result["issues"]
        mock_count.assert_not_called()

    def test_validate_directory_partial_sidecar_reports_spots(self, tmp_path):
        """A store-linked accession whose sidecar records a partial download is failed with
        the reads-on-disk-vs-spots-at-NCBI message."""
        from prepare.store.sidecar import Sidecar, write_sidecar

        store_acc_dir = tmp_path / "store" / "sra" / "SRR123"
        store_acc_dir.mkdir(parents=True)
        (store_acc_dir / "SRR123.fastq").write_text("@r\nACGT\n+\nIIII\n")
        write_sidecar(
            store_acc_dir / "SRR123.json",
            Sidecar(
                accession="SRR123",
                state="partial",
                reads_per_mate=5,
                ncbi={"spots": 100},
            ),
        )

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        acc_dir = fastq_folder / "SRR123"
        acc_dir.symlink_to(store_acc_dir)

        command = SRAValidateCommand()
        result = command._validate_directory(acc_dir)

        assert result["status"] == "FAILED"
        assert "partial: 5 reads on disk vs 100 spots at NCBI" in result["issues"]

    @pytest.mark.parametrize(
        "state, expected",
        [("failed", "store state failed: see store_verify"), ("downloading", "download in progress elsewhere")],
    )
    def test_validate_directory_flags_failed_and_downloading_sidecars(self, tmp_path, state, expected):
        """Only a complete or adopted dataset passes: a failed download, and one another
        project is still downloading, are not finished datasets even when their files parse.

        A "failed" sidecar with no recorded error falls back to pointing at store_verify,
        since there is nothing more specific to tell the reader.
        """
        from prepare.store.sidecar import Sidecar, write_sidecar

        store_acc_dir = tmp_path / "store" / "sra" / "SRR123"
        store_acc_dir.mkdir(parents=True)
        (store_acc_dir / "SRR123.fastq").write_text("@r\nACGT\n+\nIIII\n")
        write_sidecar(store_acc_dir / "SRR123.json", Sidecar(accession="SRR123", state=state))

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        acc_dir = fastq_folder / "SRR123"
        acc_dir.symlink_to(store_acc_dir)

        result = SRAValidateCommand()._validate_directory(acc_dir)

        assert result["status"] == "FAILED"
        assert expected in result["issues"]

    def test_validate_directory_failed_sidecar_reports_its_own_error(self, tmp_path):
        """When the sidecar recorded why the download failed, that reason is surfaced instead
        of the generic 'see store_verify' fallback."""
        from prepare.store.sidecar import Sidecar, write_sidecar

        store_acc_dir = tmp_path / "store" / "sra" / "SRR123"
        store_acc_dir.mkdir(parents=True)
        (store_acc_dir / "SRR123.fastq").write_text("@r\nACGT\n+\nIIII\n")
        write_sidecar(
            store_acc_dir / "SRR123.json",
            Sidecar(accession="SRR123", state="failed", error="connection reset by NCBI"),
        )

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        acc_dir = fastq_folder / "SRR123"
        acc_dir.symlink_to(store_acc_dir)

        result = SRAValidateCommand()._validate_directory(acc_dir)

        assert result["status"] == "FAILED"
        assert "store state failed: connection reset by NCBI" in result["issues"]

    def test_validate_directory_passes_a_complete_sidecar(self, tmp_path):
        from prepare.store.sidecar import Sidecar, write_sidecar

        store_acc_dir = tmp_path / "store" / "sra" / "SRR123"
        store_acc_dir.mkdir(parents=True)
        (store_acc_dir / "SRR123.fastq").write_text("@r\nACGT\n+\nIIII\n")
        write_sidecar(store_acc_dir / "SRR123.json", Sidecar(accession="SRR123", state="complete"))

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        acc_dir = fastq_folder / "SRR123"
        acc_dir.symlink_to(store_acc_dir)

        assert SRAValidateCommand()._validate_directory(acc_dir)["status"] == "PASSED"

    def test_validate_directory_registry_verdict_truncated_reports_spots(self, tmp_path):
        """A plain project directory (no store sidecar) falls back to the registry's own
        download verdict for the same completeness check."""
        from prepare.data.registry import Registry

        acc_dir = tmp_path / "SRR123"
        acc_dir.mkdir()
        (acc_dir / "SRR123.fastq").write_text("@r\nACGT\n+\nIIII\n")

        registry = Registry()
        registry.datasets["SRR123"] = {
            "download": {"complete": {"verdict": "truncated", "reads_r1": 5, "expected_spots": 100}}
        }

        command = SRAValidateCommand()
        result = command._validate_directory(acc_dir, registry)

        assert result["status"] == "FAILED"
        assert "partial: 5 reads on disk vs 100 spots at NCBI" in result["issues"]

    def test_validate_directory_md5_mismatch(self, tmp_path):
        """--md5 fails a file whose content no longer matches the sidecar's recorded md5."""
        from prepare.store.sidecar import Sidecar, write_sidecar

        store_acc_dir = tmp_path / "store" / "sra" / "SRR123"
        store_acc_dir.mkdir(parents=True)
        fastq_path = store_acc_dir / "SRR123.fastq"
        fastq_path.write_text("@r\nACGT\n+\nIIII\n")
        write_sidecar(
            store_acc_dir / "SRR123.json",
            Sidecar(
                accession="SRR123",
                files=[{"name": "SRR123.fastq", "bytes": fastq_path.stat().st_size, "md5": "0" * 32, "reads": 1}],
            ),
        )

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        acc_dir = fastq_folder / "SRR123"
        acc_dir.symlink_to(store_acc_dir)

        command = SRAValidateCommand()
        result = command._validate_directory(acc_dir, check_md5=True)

        assert result["status"] == "FAILED"
        assert "md5 mismatch: SRR123.fastq" in result["issues"]
        assert "md5" in result["checks"]

    def test_validate_directory_md5_match_passes(self, tmp_path):
        """--md5 passes when the file's md5 matches the sidecar's recorded value."""
        from prepare.store.sidecar import Sidecar, md5_file, write_sidecar

        store_acc_dir = tmp_path / "store" / "sra" / "SRR123"
        store_acc_dir.mkdir(parents=True)
        fastq_path = store_acc_dir / "SRR123.fastq"
        fastq_path.write_text("@r\nACGT\n+\nIIII\n")
        write_sidecar(
            store_acc_dir / "SRR123.json",
            Sidecar(
                accession="SRR123",
                files=[
                    {
                        "name": "SRR123.fastq",
                        "bytes": fastq_path.stat().st_size,
                        "md5": md5_file(fastq_path),
                        "reads": 1,
                    }
                ],
            ),
        )

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        acc_dir = fastq_folder / "SRR123"
        acc_dir.symlink_to(store_acc_dir)

        command = SRAValidateCommand()
        result = command._validate_directory(acc_dir, check_md5=True)

        assert result["status"] == "PASSED"

    def test_print_validation_results_success(self):
        """Test printing successful validation results."""
        command = SRAValidateCommand()

        validation_results = [
            {"accession": "SRR123", "status": "PASSED", "issues": "None"},
            {"accession": "SRR456", "status": "PASSED", "issues": "None"},
        ]

        result = command._print_validation_results(validation_results)

        assert result is True

    def test_print_validation_results_with_failures(self, capsys):
        """Test printing validation results with failures."""
        command = SRAValidateCommand()

        validation_results = [
            {"accession": "SRR123", "status": "PASSED", "issues": "None"},
            {"accession": "SRR456", "status": "FAILED", "issues": "Empty files"},
        ]

        result = command._print_validation_results(validation_results)

        assert result is False
        out = capsys.readouterr().out
        assert "Passed: 1\nFailed: 1\n" in out
        assert "  SRR456: Empty files" in out

    def test_execute_success(self, tmp_path):
        """Test successful execution."""
        command = SRAValidateCommand()

        # Create test structure
        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        acc_dir = fastq_folder / "SRR123"
        acc_dir.mkdir()
        (acc_dir / "test.fastq").write_text("content")

        args = argparse.Namespace(
            fastq_folder=str(fastq_folder),
            accession=None,
            check_pairs=False,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        # Mock _validate_directory to return success
        with patch.object(command, "_validate_directory") as mock_validate:
            mock_validate.return_value = {
                "accession": "SRR123",
                "status": "PASSED",
                "issues": "None",
                "num_files": 1,
            }

            result = command.execute(args)

        assert result == 0
        mock_validate.assert_called_once()

    def test_execute_registry_lock_timeout_exits_4(self, tmp_path):
        """A registry lock wait that gave up is retryable: exit 4, not the broad catch's 1."""
        from prepare.core.exceptions import LockTimeoutError

        command = SRAValidateCommand()
        fastq_folder = tmp_path / "fastq"
        (fastq_folder / "SRR123").mkdir(parents=True)
        (fastq_folder / "SRR123" / "test.fastq").write_text("@r\nACGT\n+\n!!!!\n")
        args = argparse.Namespace(
            fastq_folder=str(fastq_folder),
            accession=None,
            check_pairs=False,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )
        with patch("prepare.cli.commands.sra_enhanced.registry_batch", side_effect=LockTimeoutError("registry lock")):
            assert command.execute(args) == 4

    def test_execute_records_analyses_in_registry(self, tmp_path):
        """Each validated accession is recorded with its pass/fail status and file count."""
        command = SRAValidateCommand()

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        acc_dir = fastq_folder / "SRR123"
        acc_dir.mkdir()
        (acc_dir / "test.fastq").write_text("@r\nACGT\n+\n!!!!\n")
        registry_path = tmp_path / "prepare_registry.json"

        args = argparse.Namespace(
            fastq_folder=str(fastq_folder),
            accession=None,
            check_pairs=False,
            registry=str(registry_path),
            data_root=None,
        )

        with patch.object(command, "_validate_directory") as mock_validate:
            mock_validate.return_value = {
                "accession": "SRR123",
                "status": "PASSED",
                "issues": "None",
                "issues_list": [],
                "num_files": 1,
                "checks": ["empty_files", "format", "completeness"],
            }
            result = command.execute(args)

        assert result == 0
        registry = json.loads(registry_path.read_text())
        assert registry["datasets"]["SRR123"]["analyses"]["validate"]["summary"] == {
            "passed": True,
            "files": 1,
            "issues": [],
        }

    def test_execute_records_usage_in_store_catalogue(self, tmp_path):
        """A validated accession is also recorded as 'analysed' usage in the store catalogue."""
        from prepare.data.registry import load_registry as _load, save_registry as _save
        from prepare.store.catalog import Catalog
        from prepare.store.layout import init_store, store_paths

        command = SRAValidateCommand()

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        acc_dir = fastq_folder / "SRR123"
        acc_dir.mkdir()
        (acc_dir / "test.fastq").write_text("@r\nACGT\n+\n!!!!\n")
        registry_path = tmp_path / "prepare_registry.json"
        store_root = tmp_path / "store"
        init_store(store_root)

        registry = _load(registry_path)
        registry.project = {"id": "proj1", "name": "demo", "path": str(tmp_path), "created": "now"}
        _save(registry)

        args = argparse.Namespace(
            fastq_folder=str(fastq_folder),
            accession=None,
            check_pairs=False,
            registry=str(registry_path),
            data_root=str(store_root),
        )

        with patch.object(command, "_validate_directory") as mock_validate:
            mock_validate.return_value = {
                "accession": "SRR123",
                "status": "PASSED",
                "issues": "None",
                "num_files": 1,
            }
            result = command.execute(args)

        assert result == 0
        with Catalog(store_paths(store_root)) as catalog:
            row = catalog.conn.execute(
                "SELECT stage FROM usage WHERE accession = ? AND project_id = ?", ("SRR123", "proj1")
            ).fetchone()
        assert row["stage"] == "analysed"

    def test_catalog_failure_leaves_validate_outcome_unchanged(self, tmp_path):
        """A broken catalogue write never changes validate's registry outcome or exit code."""
        from prepare.data.registry import load_registry as _load, save_registry as _save
        from prepare.store.layout import init_store

        command = SRAValidateCommand()

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        acc_dir = fastq_folder / "SRR123"
        acc_dir.mkdir()
        (acc_dir / "test.fastq").write_text("@r\nACGT\n+\n!!!!\n")
        registry_path = tmp_path / "prepare_registry.json"
        store_root = tmp_path / "store"
        init_store(store_root)

        registry = _load(registry_path)
        registry.project = {"id": "proj1", "name": "demo", "path": str(tmp_path), "created": "now"}
        _save(registry)

        args = argparse.Namespace(
            fastq_folder=str(fastq_folder),
            accession=None,
            check_pairs=False,
            registry=str(registry_path),
            data_root=str(store_root),
        )

        with patch.object(command, "_validate_directory") as mock_validate:
            mock_validate.return_value = {
                "accession": "SRR123",
                "status": "PASSED",
                "issues": "None",
                "num_files": 1,
            }
            with patch("prepare.store.usage.catalog_write", side_effect=RuntimeError("locked")):
                result = command.execute(args)

        assert result == 0
        registry_after = json.loads(registry_path.read_text())
        assert registry_after["datasets"]["SRR123"]["analyses"]["validate"]["summary"] == {
            "passed": True,
            "files": 1,
            "issues": [],
        }

    def test_execute_records_failing_issues_in_registry_summary(self, tmp_path):
        """A failed accession's registry summary carries the raw issue list, not just the
        joined display string, and names the checks that ran."""
        command = SRAValidateCommand()

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        acc_dir = fastq_folder / "SRR123"
        acc_dir.mkdir()
        (acc_dir / "test.fastq").touch()  # empty file: fails the "empty_files" check
        registry_path = tmp_path / "prepare_registry.json"

        args = argparse.Namespace(
            fastq_folder=str(fastq_folder),
            accession=None,
            check_pairs=False,
            md5=False,
            registry=str(registry_path),
            data_root=None,
        )

        result = command.execute(args)

        assert result == 1
        registry = json.loads(registry_path.read_text())
        summary = registry["datasets"]["SRR123"]["analyses"]["validate"]["summary"]
        assert summary["passed"] is False
        assert summary["files"] == 1
        assert summary["issues"] == ["Empty file: test.fastq"]

    def test_execute_folder_not_exists(self, caplog, capsys):
        """Test execution when fastq folder doesn't exist."""
        command = SRAValidateCommand()
        args = argparse.Namespace(fastq_folder="/nonexistent/folder", accession=None, check_pairs=False)

        result = command.execute(args)

        assert result == 1
        assert "FASTQ folder /nonexistent/folder does not exist" in caplog.text
        assert "does not exist" not in capsys.readouterr().out

    def test_execute_no_accession_dirs(self, tmp_path, caplog):
        """Test execution when no accession directories found."""
        command = SRAValidateCommand()

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()

        args = argparse.Namespace(fastq_folder=str(fastq_folder), accession=None, check_pairs=False)

        result = command.execute(args)

        assert result == 1
        assert "No accession directories found" in caplog.text

    @patch("prepare.cli.commands.sra_enhanced.logger")
    def test_execute_exception_handling(self, mock_logger, tmp_path):
        """Test exception handling during execution."""
        command = SRAValidateCommand()

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        acc_dir = fastq_folder / "SRR123"
        acc_dir.mkdir()

        args = argparse.Namespace(
            fastq_folder=str(fastq_folder),
            accession=None,
            check_pairs=False,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        # Mock _validate_directory to raise an exception
        with patch.object(command, "_validate_directory", side_effect=Exception("Test error")):
            result = command.execute(args)

        assert result == 1
        mock_logger.error.assert_called_once()


def test_validate_restricts_to_the_accessions_file_and_flags(tmp_path):
    """--accessions-file and --accession together name the folders to validate (--accessions went in 0.5.0)."""
    fastq_folder = tmp_path / "fastq"
    for accession in ("SRR1", "SRR2", "SRR3"):
        (fastq_folder / accession).mkdir(parents=True)
        (fastq_folder / accession / f"{accession}.fastq").write_text("@r\nACGT\n+\nIIII\n")
    accessions_file = tmp_path / "acc.txt"
    accessions_file.write_text("SRR1\n")
    args = argparse.Namespace(
        fastq_folder=str(fastq_folder),
        accessions_file=str(accessions_file),
        accession=["SRR3"],
        check_pairs=False,
        registry=str(tmp_path / "prepare_registry.json"),
        data_root=None,
    )

    assert SRAValidateCommand().execute(args) == 0

    registry = json.loads((tmp_path / "prepare_registry.json").read_text())
    assert set(registry["datasets"]) == {"SRR1", "SRR3"}
