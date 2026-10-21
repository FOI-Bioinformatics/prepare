"""
Test CLI commands functionality.

This module tests the CLI command classes and their execute methods,
focusing on argument parsing, validation, and proper delegation.
"""

import argparse
import logging
import json
import os
import threading
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from prepare.cli.commands.branchwater import (
    UseBranchwaterCommand,
    ExtractBranchwaterMetadataCommand,
)
from prepare.cli.commands.containment import ParseContainmentCommand, PlotContainmentCommand
from prepare.cli.commands.metadata import (
    DownloadMetadataCommand,
    ParseMetadataCommand,
    CountMetadataCommand,
    PlotMetadataCountsCommand,
    _metadata_fields,
)
from prepare.cli.commands.sra import DownloadSraCommand
from prepare.cli.commands.samples import SingleSampleCommand
from prepare.cli.commands.test_data import DownloadTestGenomeCommand
from prepare.core.constants import DEFAULT_REGISTRY_MAX_SCREENED, FAILED_ACCESSIONS_FILE
from prepare.core.exceptions import PrepareError
from prepare.core.models import ContainmentSummary
from prepare.data import registry_blocks as rb
from prepare.data.registry import load_registry, record_download, record_exclusion, record_metadata, save_registry


class TestUseBranchwaterCommand:
    """Test UseBranchwaterCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = UseBranchwaterCommand()
        assert command.name == "use_branchwater"
        assert "process" in command.help.lower()
        assert "branchwater" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = UseBranchwaterCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Check required arguments are added
        args = parser.parse_args(["--branchwater-folder", "test_folder"])
        assert args.branchwater_folder == "test_folder"
        assert args.matches_folder == "matches"  # default value

        # Test with custom matches folder
        args = parser.parse_args(["--branchwater-folder", "test_folder", "--matches-folder", "custom_matches"])
        assert args.matches_folder == "custom_matches"

    def test_configure_parser_missing_required(self):
        """Test parser with missing required arguments."""
        command = UseBranchwaterCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        with pytest.raises(SystemExit):
            parser.parse_args([])  # Missing required --branchwater-folder

    @patch("prepare.cli.commands.branchwater.process_branchwater_files")
    def test_execute_success(self, mock_command, tmp_path, monkeypatch):
        """Test successful command execution; the matches folder resolves against the project root."""
        monkeypatch.chdir(tmp_path)
        mock_command.return_value = {"file1": Path("test")}
        command = UseBranchwaterCommand()

        args = argparse.Namespace(branchwater_folder="test_folder", matches_folder="matches")

        result = command.execute(args)
        assert result == 0
        mock_command.assert_called_once_with("test_folder", str(tmp_path.resolve() / "matches"), errors=[])

    @patch("prepare.cli.commands.branchwater.process_branchwater_files")
    def test_execute_failure(self, mock_command, tmp_path, monkeypatch):
        """Test command execution failure."""
        monkeypatch.chdir(tmp_path)
        mock_command.side_effect = PrepareError("Test error")
        command = UseBranchwaterCommand()

        args = argparse.Namespace(branchwater_folder="test_folder", matches_folder="matches")

        result = command.execute(args)
        assert result == 1
        mock_command.assert_called_once_with("test_folder", str(tmp_path.resolve() / "matches"), errors=[])

    def test_use_branchwater_returns_1_on_file_error(self, tmp_path):
        source = tmp_path / "bw"
        source.mkdir()
        (source / "bad.csv").write_bytes(b"\x00\x05\x16\x07\xb0")
        args = argparse.Namespace(branchwater_folder=str(source), matches_folder=str(tmp_path / "m"))
        assert UseBranchwaterCommand().execute(args) == 1

    @patch("prepare.cli.commands.branchwater.process_branchwater_files")
    def test_execute_logs_a_next_hint(self, mock_command, caplog):
        """On success, point at the next pipeline step (parse_containment)."""
        mock_command.return_value = {"file1": Path("test")}
        command = UseBranchwaterCommand()
        args = argparse.Namespace(branchwater_folder="test_folder", matches_folder="matches")

        with caplog.at_level("INFO"):
            result = command.execute(args)

        assert result == 0
        assert "Next:" in caplog.text
        assert "parse_containment" in caplog.text
        assert "matches" in caplog.text


class TestExtractBranchwaterMetadataCommand:
    """Test ExtractBranchwaterMetadataCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = ExtractBranchwaterMetadataCommand()
        assert command.name == "extract_branchwater_metadata"
        assert "extract" in command.help.lower()
        assert "metadata" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = ExtractBranchwaterMetadataCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(["--branchwater-folder", "test_folder"])
        assert args.branchwater_folder == "test_folder"
        assert args.metadata_folder == "metadata"

    @patch("prepare.cli.commands.branchwater.extract_metadata_from_branchwater")
    @patch("pathlib.Path.mkdir")
    def test_execute(self, mock_mkdir, mock_command, tmp_path, monkeypatch):
        """Test command execution."""
        import pandas as pd

        monkeypatch.chdir(tmp_path)
        mock_command.return_value = pd.DataFrame()
        command = ExtractBranchwaterMetadataCommand()

        args = argparse.Namespace(branchwater_folder="test_folder", metadata_folder="metadata")

        result = command.execute(args)
        assert result == 0
        mock_mkdir.assert_called_once_with(exist_ok=True)
        expected = tmp_path.resolve() / "metadata" / "branchwater_metadata.txt"
        mock_command.assert_called_once_with("test_folder", expected, errors=[])


class TestParseContainmentCommand:
    """Test ParseContainmentCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = ParseContainmentCommand()
        assert command.name == "parse_containment"
        assert "parse" in command.help.lower()
        assert "containment" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = ParseContainmentCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test with required argument
        args = parser.parse_args(["--matches-folder", "test_matches"])
        assert args.matches_folder == "test_matches"
        assert args.parsed_containment_file == "parsed_containment.txt"
        assert args.summary_containment_file == "top_containments.txt"
        assert args.step_size == 0.1
        assert args.registry_max_screened == DEFAULT_REGISTRY_MAX_SCREENED
        assert args.details_file is None

    def test_configure_parser_with_optional_args(self):
        """Test parser with optional arguments."""
        command = ParseContainmentCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(
            [
                "--matches-folder",
                "test_matches",
                "--parsed-containment-file",
                "custom_parsed.txt",
                "--summary-containment-file",
                "custom_summary.txt",
                "--step-size",
                "0.1",
                "--details-file",
                "custom_details.tsv",
            ]
        )

        assert args.parsed_containment_file == "custom_parsed.txt"
        assert args.summary_containment_file == "custom_summary.txt"
        assert args.step_size == 0.1
        assert args.details_file == "custom_details.tsv"

    @patch("prepare.cli.commands.containment.parse_containment_data")
    def test_execute(self, mock_command, tmp_path):
        """Test command execution."""
        mock_command.return_value = None
        command = ParseContainmentCommand()

        args = argparse.Namespace(
            matches_folder="test_matches",
            parsed_containment_file=str(tmp_path / "parsed.txt"),
            summary_containment_file="summary.txt",
            step_size=0.05,
            details_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            registry_max_screened=DEFAULT_REGISTRY_MAX_SCREENED,
        )

        result = command.execute(args)
        assert result == 0
        mock_command.assert_called_once_with(
            str(tmp_path.resolve() / "test_matches"),
            str(tmp_path / "parsed.txt"),
            "summary.txt",
            0.05,
            details_file=None,
        )

    def test_exit_1_leaves_registry_and_table_unchanged_on_an_unreadable_file(self, tmp_path):
        """A bad match file: exit 1, no table written and the registry's screening blocks kept."""
        matches_folder = tmp_path / "matches"
        matches_folder.mkdir()
        (matches_folder / "GCF_A.csv").write_text("acc,containment\nSRR1,0.9\n")
        registry_path = tmp_path / "prepare_registry.json"
        args = argparse.Namespace(
            matches_folder=str(matches_folder),
            parsed_containment_file=str(tmp_path / "parsed.txt"),
            summary_containment_file=str(tmp_path / "summary.txt"),
            step_size=0.1,
            details_file=None,
            registry=str(registry_path),
            registry_max_screened=DEFAULT_REGISTRY_MAX_SCREENED,
        )
        assert ParseContainmentCommand().execute(args) == 0
        before = registry_path.read_text()
        table_before = (tmp_path / "parsed.txt").read_bytes()

        (matches_folder / "GCF_B.csv").write_text("acc,containment\nSRR3,0.8\n")
        (matches_folder / "bad.csv").write_text("run,score\nSRR2,0.5\n")
        assert ParseContainmentCommand().execute(args) == 1
        assert registry_path.read_text() == before
        assert (tmp_path / "parsed.txt").read_bytes() == table_before

    def test_empty_matches_folder_is_exit_1(self, tmp_path, caplog):
        """No match CSV at all: exit 1 naming the folder, and no registry is created."""
        matches_folder = tmp_path / "matches"
        matches_folder.mkdir()
        args = argparse.Namespace(
            matches_folder=str(matches_folder),
            parsed_containment_file=str(tmp_path / "parsed.txt"),
            summary_containment_file=str(tmp_path / "summary.txt"),
            step_size=0.1,
            details_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            registry_max_screened=DEFAULT_REGISTRY_MAX_SCREENED,
        )
        with caplog.at_level(logging.ERROR):
            assert ParseContainmentCommand().execute(args) == 1
        assert "No match CSV in" in caplog.text and "matches" in caplog.text
        assert not (tmp_path / "prepare_registry.json").exists()

    def test_execute_records_screening_in_registry(self, tmp_path):
        """Every accession in the parsed containment table is recorded as screened from matches."""
        matches_folder = tmp_path / "matches"
        matches_folder.mkdir()
        (matches_folder / "GCF_A.csv").write_text("acc,containment\nSRR1,0.9\nSRR2,0.1\n")

        command = ParseContainmentCommand()
        args = argparse.Namespace(
            matches_folder=str(matches_folder),
            parsed_containment_file=str(tmp_path / "parsed.txt"),
            summary_containment_file=str(tmp_path / "summary.txt"),
            step_size=0.1,
            details_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            registry_max_screened=DEFAULT_REGISTRY_MAX_SCREENED,
        )

        result = command.execute(args)
        assert result == 0

        data = json.loads((tmp_path / "prepare_registry.json").read_text())
        for acc in ("SRR1", "SRR2"):
            screening = data["datasets"][acc]["screening"]
            assert "GCF_A" in screening["genomes"]
            assert screening["genomes"]["GCF_A"]["source"] == "matches"

    @patch("prepare.cli.commands.containment.parse_containment_data", return_value=None)
    def test_a_table_read_from_disk_is_read_outside_the_registry_lock(self, _parse, tmp_path):
        """With no table in the summary, the written file is read before the registry lock is taken."""
        parsed = tmp_path / "parsed.txt"
        parsed.write_text("\tGCF_A\tmax_containment\nSRR1\t0.9\t0.9\n")
        lock = tmp_path / "prepare_registry.json.lock"
        import prepare.data.screening_table as table_module

        real_read = table_module.pd.read_csv
        lock_held_while_reading = []

        def spy_read(*args, **kwargs):
            lock_held_while_reading.append(lock.exists())
            return real_read(*args, **kwargs)

        args = argparse.Namespace(
            matches_folder=str(tmp_path / "matches"),
            parsed_containment_file=str(parsed),
            summary_containment_file=str(tmp_path / "summary.txt"),
            step_size=0.1,
            details_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            registry_max_screened=DEFAULT_REGISTRY_MAX_SCREENED,
        )
        with patch.object(table_module.pd, "read_csv", side_effect=spy_read):
            assert ParseContainmentCommand().execute(args) == 0
        assert lock_held_while_reading == [False]
        data = json.loads((tmp_path / "prepare_registry.json").read_text())
        assert "GCF_A" in data["datasets"]["SRR1"]["screening"]["genomes"]

    def test_execute_writes_details_file_default_path(self, tmp_path):
        """The details table (cANI, sample metadata) is written next to the parsed containment table."""
        matches_folder = tmp_path / "matches"
        matches_folder.mkdir()
        (matches_folder / "GCF_A.csv").write_text("acc,containment,cANI,biosample\nSRR1,0.9,0.98,SAMN1\n")

        command = ParseContainmentCommand()
        parsed_file = tmp_path / "parsed.txt"
        args = argparse.Namespace(
            matches_folder=str(matches_folder),
            parsed_containment_file=str(parsed_file),
            summary_containment_file=str(tmp_path / "summary.txt"),
            step_size=0.1,
            details_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            registry_max_screened=DEFAULT_REGISTRY_MAX_SCREENED,
        )

        result = command.execute(args)
        assert result == 0

        details_file = tmp_path / "parsed_details.tsv"
        assert details_file.exists()
        content = details_file.read_text()
        assert "SRR1" in content and "GCF_A" in content and "0.98" in content

    def test_execute_writes_details_file_explicit_path(self, tmp_path):
        matches_folder = tmp_path / "matches"
        matches_folder.mkdir()
        (matches_folder / "GCF_A.csv").write_text("acc,containment\nSRR1,0.9\n")

        command = ParseContainmentCommand()
        custom_details = tmp_path / "custom_details.tsv"
        args = argparse.Namespace(
            matches_folder=str(matches_folder),
            parsed_containment_file=str(tmp_path / "parsed.txt"),
            summary_containment_file=str(tmp_path / "summary.txt"),
            step_size=0.1,
            details_file=str(custom_details),
            registry=str(tmp_path / "prepare_registry.json"),
            registry_max_screened=DEFAULT_REGISTRY_MAX_SCREENED,
        )

        result = command.execute(args)
        assert result == 0
        assert custom_details.exists()
        assert not (tmp_path / "parsed_details.tsv").exists()

    def test_execute_logs_a_next_hint(self, tmp_path, caplog):
        """On success, point at plot_containment and select_datasets as the next steps."""
        matches_folder = tmp_path / "matches"
        matches_folder.mkdir()
        (matches_folder / "GCF_A.csv").write_text("acc,containment\nSRR1,0.9\nSRR2,0.1\n")

        command = ParseContainmentCommand()
        parsed_file = tmp_path / "parsed.txt"
        args = argparse.Namespace(
            matches_folder=str(matches_folder),
            parsed_containment_file=str(parsed_file),
            summary_containment_file=str(tmp_path / "summary.txt"),
            step_size=0.1,
            details_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            registry_max_screened=DEFAULT_REGISTRY_MAX_SCREENED,
        )

        with caplog.at_level("INFO"):
            result = command.execute(args)

        assert result == 0
        assert "Next:" in caplog.text
        assert "plot_containment" in caplog.text
        assert "select_datasets" in caplog.text
        assert "GCF_A" in caplog.text

    @patch("prepare.cli.commands.containment.parse_containment_data")
    def test_execute_next_hint_names_the_first_genome_column_sorted(self, mock_command, tmp_path, caplog):
        """With several genome columns, the hint's --genome-id must name the alphabetically
        first one, not whichever happened to be inserted first into genome_to_samples."""
        summary = ContainmentSummary()
        summary.genome_to_samples = {"GCF_B": ["SRR1"], "GCF_A": ["SRR2"]}
        mock_command.return_value = summary

        command = ParseContainmentCommand()
        args = argparse.Namespace(
            matches_folder="test_matches",
            parsed_containment_file=str(tmp_path / "parsed.txt"),
            summary_containment_file="summary.txt",
            step_size=0.1,
            details_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            registry_max_screened=DEFAULT_REGISTRY_MAX_SCREENED,
        )

        with caplog.at_level("INFO"):
            result = command.execute(args)

        assert result == 0
        assert "--genome-id GCF_A" in caplog.text
        assert "--genome-id GCF_B" not in caplog.text


class TestDownloadMetadataCommand:
    """Test DownloadMetadataCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = DownloadMetadataCommand()
        assert command.name == "download_metadata"
        assert "download" in command.help.lower()
        assert "metadata" in command.help.lower()

    def test_configure_parser(self, monkeypatch):
        """Test parser configuration."""
        monkeypatch.delenv("NCBI_API_KEY", raising=False)
        command = DownloadMetadataCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test with required email
        args = parser.parse_args(["--email", "test@example.com"])
        assert args.email == "test@example.com"
        assert args.matches_folder == "matches"
        assert args.metadata_folder == "metadata"
        assert args.threshold == 0.0
        assert args.dry_run is False
        assert args.api_key is None
        assert args.batch_size == 200

    @patch("prepare.cli.commands.metadata.download_metadata", return_value={})
    def test_api_key_defaults_from_environment(self, mock_download, monkeypatch):
        """Without --api-key, the NCBI_API_KEY environment variable reaches download_metadata."""
        monkeypatch.setenv("NCBI_API_KEY", "env-key")
        command = DownloadMetadataCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(["--email", "test@example.com", "--dry-run"])
        assert args.api_key is None  # resolved through the runtime settings, not at parse time
        assert command.execute(args) == 0
        assert mock_download.call_args.kwargs["api_key"] == "env-key"

    def test_configure_parser_with_options(self):
        """Test parser with optional arguments."""
        command = DownloadMetadataCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(
            [
                "--email",
                "test@example.com",
                "--matches-folder",
                "custom_matches",
                "--metadata-folder",
                "custom_metadata",
                "--threshold",
                "0.5",
                "--dry-run",
                "--api-key",
                "cli-key",
                "--batch-size",
                "50",
            ]
        )

        assert args.threshold == 0.5
        assert args.dry_run is True
        assert args.api_key == "cli-key"
        assert args.batch_size == 50

    @patch("prepare.cli.commands.metadata.download_metadata")
    def test_execute(self, mock_command, tmp_path):
        """Test command execution."""
        mock_command.return_value = {}
        command = DownloadMetadataCommand()

        args = argparse.Namespace(
            email="test@example.com",
            matches_folder="matches",
            metadata_folder="metadata",
            threshold=0.0,
            dry_run=False,
            accessions_file=None,
            api_key=None,
            batch_size=200,
            registry=str(tmp_path / "prepare_registry.json"),
        )

        result = command.execute(args)
        assert result == 0
        mock_command.assert_called_once_with(
            email="test@example.com",
            matches_folder=str(tmp_path.resolve() / "matches"),
            metadata_folder=str(tmp_path.resolve() / "metadata"),
            threshold=0.0,
            dry_run=False,
            accessions_file=None,
            api_key=None,
            batch_size=200,
        )

    def test_execute_records_metadata_in_registry(self, tmp_path):
        """Every accession download_metadata reports downloaded is recorded with its XML path."""
        metadata_folder = tmp_path / "metadata"
        accessions_file = tmp_path / "accessions.txt"
        accessions_file.write_text("SRR1\nSRR2\n")

        command = DownloadMetadataCommand()
        args = argparse.Namespace(
            email="test@example.com",
            matches_folder=str(tmp_path / "matches"),
            metadata_folder=str(metadata_folder),
            threshold=0.0,
            dry_run=False,
            accessions_file=str(accessions_file),
            api_key=None,
            batch_size=200,
            registry=str(tmp_path / "prepare_registry.json"),
        )

        def fake_batch(batch, metadata_path, email, api_key):
            successes = {}
            for accession in batch:
                xml_path = metadata_path / f"{accession}_metadata.xml"
                xml_path.write_text("<root/>")
                successes[accession] = xml_path
            return successes, {}

        with patch("prepare.data.metadata._download_batch_metadata", side_effect=fake_batch):
            result = command.execute(args)

        assert result == 0
        registry = json.loads((tmp_path / "prepare_registry.json").read_text())
        for acc in ("SRR1", "SRR2"):
            assert registry["datasets"][acc]["metadata"]["xml"].endswith(f"{acc}_metadata.xml")

    def test_execute_dry_run_records_nothing(self, tmp_path):
        """A dry run does not touch the registry."""
        accessions_file = tmp_path / "accessions.txt"
        accessions_file.write_text("SRR1\n")

        command = DownloadMetadataCommand()
        args = argparse.Namespace(
            email="test@example.com",
            matches_folder=str(tmp_path / "matches"),
            metadata_folder=str(tmp_path / "metadata"),
            threshold=0.0,
            dry_run=True,
            accessions_file=str(accessions_file),
            api_key=None,
            batch_size=200,
            registry=str(tmp_path / "prepare_registry.json"),
        )

        result = command.execute(args)
        assert result == 0
        assert not (tmp_path / "prepare_registry.json").exists()

    def test_execute_records_parsed_metadata_fields(self, tmp_path):
        """Downloaded metadata records include parsed fields (run_total_spots, run_md5, etc)."""
        metadata_folder = tmp_path / "metadata"
        metadata_folder.mkdir(parents=True)

        # Create a test XML file with metadata
        xml_content = """<?xml version="1.0"?>
<EXPERIMENT_PACKAGE_SET>
  <EXPERIMENT_PACKAGE>
    <STUDY>
      <IDENTIFIERS><PRIMARY_ID>PRJ123</PRIMARY_ID></IDENTIFIERS>
      <DESCRIPTOR><STUDY_TITLE>Test Study</STUDY_TITLE></DESCRIPTOR>
    </STUDY>
    <SAMPLE>
      <IDENTIFIERS><PRIMARY_ID>SAM123</PRIMARY_ID></IDENTIFIERS>
      <SAMPLE_NAME>
        <SCIENTIFIC_NAME>Escherichia coli</SCIENTIFIC_NAME>
      </SAMPLE_NAME>
    </SAMPLE>
    <EXPERIMENT>
      <IDENTIFIERS><PRIMARY_ID>EXP123</PRIMARY_ID></IDENTIFIERS>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
        <LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT>
      </LIBRARY_DESCRIPTOR>
      <PLATFORM><ILLUMINA/></PLATFORM>
    </EXPERIMENT>
    <RUN_SET>
      <RUN accession="SRR123" total_spots="1234" total_bases="5678" size="999">
        <IDENTIFIERS><PRIMARY_ID>SRR123</PRIMARY_ID></IDENTIFIERS>
        <SRAFiles>
          <SRAFile md5="abc123def456" filename="SRR123.fastq.gz"/>
        </SRAFiles>
      </RUN>
    </RUN_SET>
  </EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>
"""
        xml_path = metadata_folder / "SRR123_metadata.xml"
        xml_path.write_text(xml_content)

        command = DownloadMetadataCommand()
        args = argparse.Namespace(
            email="test@example.com",
            matches_folder=str(tmp_path / "matches"),
            metadata_folder=str(metadata_folder),
            threshold=0.0,
            dry_run=False,
            accessions_file=None,
            api_key=None,
            batch_size=200,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        mock_downloaded = {"SRR123": xml_path}
        with patch("prepare.cli.commands.metadata.download_metadata", return_value=mock_downloaded):
            result = command.execute(args)

        assert result == 0
        registry = json.loads((tmp_path / "prepare_registry.json").read_text())
        assert "SRR123" in registry["datasets"]
        metadata_record = registry["datasets"]["SRR123"]["metadata"]

        # Check that parsed fields are recorded
        assert metadata_record["run_total_spots"] == 1234
        assert isinstance(metadata_record["run_total_spots"], int)
        assert metadata_record["run_md5"] == "abc123def456"
        assert metadata_record["run_size"] == "999"
        assert metadata_record["organism"] == "Escherichia coli"
        assert metadata_record["platform"] == "ILLUMINA"

    def test_execute_records_malformed_xml_with_warning(self, tmp_path, caplog):
        """Malformed XML logs a warning and records empty fields; download continues."""
        metadata_folder = tmp_path / "metadata"
        metadata_folder.mkdir(parents=True)

        # Create a malformed XML file
        xml_path_bad = metadata_folder / "SRR123_metadata.xml"
        xml_path_bad.write_text("<invalid><xml>")

        # Create a valid XML file
        xml_content = """<?xml version="1.0"?>
<EXPERIMENT_PACKAGE_SET>
  <EXPERIMENT_PACKAGE>
    <STUDY>
      <IDENTIFIERS><PRIMARY_ID>PRJ456</PRIMARY_ID></IDENTIFIERS>
      <DESCRIPTOR><STUDY_TITLE>Test Study</STUDY_TITLE></DESCRIPTOR>
    </STUDY>
    <SAMPLE>
      <IDENTIFIERS><PRIMARY_ID>SAM456</PRIMARY_ID></IDENTIFIERS>
      <SAMPLE_NAME>
        <SCIENTIFIC_NAME>Bacillus subtilis</SCIENTIFIC_NAME>
      </SAMPLE_NAME>
    </SAMPLE>
    <EXPERIMENT>
      <IDENTIFIERS><PRIMARY_ID>EXP456</PRIMARY_ID></IDENTIFIERS>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
      </LIBRARY_DESCRIPTOR>
    </EXPERIMENT>
    <RUN_SET>
      <RUN accession="SRR456" total_spots="5000" total_bases="10000" size="2000">
        <IDENTIFIERS><PRIMARY_ID>SRR456</PRIMARY_ID></IDENTIFIERS>
        <SRAFiles>
          <SRAFile md5="xyz789" filename="SRR456.fastq.gz"/>
        </SRAFiles>
      </RUN>
    </RUN_SET>
  </EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>
"""
        xml_path_good = metadata_folder / "SRR456_metadata.xml"
        xml_path_good.write_text(xml_content)

        command = DownloadMetadataCommand()
        args = argparse.Namespace(
            email="test@example.com",
            matches_folder=str(tmp_path / "matches"),
            metadata_folder=str(metadata_folder),
            threshold=0.0,
            dry_run=False,
            accessions_file=None,
            api_key=None,
            batch_size=200,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        mock_downloaded = {"SRR123": xml_path_bad, "SRR456": xml_path_good}
        with patch("prepare.cli.commands.metadata.download_metadata", return_value=mock_downloaded):
            result = command.execute(args)

        assert result == 0
        registry = json.loads((tmp_path / "prepare_registry.json").read_text())

        # Malformed XML: only xml path should be recorded
        assert "SRR123" in registry["datasets"]
        assert registry["datasets"]["SRR123"]["metadata"]["run_total_spots"] is None
        assert registry["datasets"]["SRR123"]["metadata"]["xml"].endswith("SRR123_metadata.xml")

        # Check that warning was logged for malformed XML
        assert "Could not parse metadata for SRR123" in caplog.text and "recorded the file path only" in caplog.text

        # Good XML: all fields should be recorded
        assert "SRR456" in registry["datasets"]
        assert registry["datasets"]["SRR456"]["metadata"]["run_total_spots"] == 5000
        assert registry["datasets"]["SRR456"]["metadata"]["run_md5"] == "xyz789"

    def test_download_and_parse_record_identical_fields(self, tmp_path):
        """DownloadMetadataCommand and ParseMetadataCommand record identical field dicts for same file."""
        metadata_folder = tmp_path / "metadata"
        metadata_folder.mkdir(parents=True)

        # Create the same XML for both workflows
        xml_content = """<?xml version="1.0"?>
<EXPERIMENT_PACKAGE_SET>
  <EXPERIMENT_PACKAGE>
    <STUDY>
      <IDENTIFIERS><PRIMARY_ID>PRJ999</PRIMARY_ID></IDENTIFIERS>
      <DESCRIPTOR><STUDY_TITLE>Consistency Test</STUDY_TITLE></DESCRIPTOR>
    </STUDY>
    <SAMPLE>
      <IDENTIFIERS><PRIMARY_ID>SAM999</PRIMARY_ID></IDENTIFIERS>
      <SAMPLE_NAME>
        <SCIENTIFIC_NAME>Test organism</SCIENTIFIC_NAME>
      </SAMPLE_NAME>
      <SAMPLE_ATTRIBUTES>
        <SAMPLE_ATTRIBUTE>
          <TAG>collection_date</TAG>
          <VALUE>2020-01-01</VALUE>
        </SAMPLE_ATTRIBUTE>
      </SAMPLE_ATTRIBUTES>
    </SAMPLE>
    <EXPERIMENT>
      <IDENTIFIERS><PRIMARY_ID>EXP999</PRIMARY_ID></IDENTIFIERS>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
        <LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT>
      </LIBRARY_DESCRIPTOR>
      <PLATFORM><ILLUMINA/></PLATFORM>
    </EXPERIMENT>
    <RUN_SET>
      <RUN accession="SRR999" total_spots="2000" total_bases="4000" size="500">
        <IDENTIFIERS><PRIMARY_ID>SRR999</PRIMARY_ID></IDENTIFIERS>
        <SRAFiles>
          <SRAFile md5="same123" filename="SRR999.fastq.gz"/>
        </SRAFiles>
      </RUN>
    </RUN_SET>
  </EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>
"""
        xml_path = metadata_folder / "SRR999_metadata.xml"
        xml_path.write_text(xml_content)

        # Test DownloadMetadataCommand
        download_args = argparse.Namespace(
            email="test@example.com",
            matches_folder=str(tmp_path / "matches"),
            metadata_folder=str(metadata_folder),
            threshold=0.0,
            dry_run=False,
            accessions_file=None,
            api_key=None,
            batch_size=200,
            registry=str(tmp_path / "prepare_registry_download.json"),
            data_root=None,
        )

        mock_downloaded = {"SRR999": xml_path}
        download_command = DownloadMetadataCommand()
        with patch("prepare.cli.commands.metadata.download_metadata", return_value=mock_downloaded):
            assert download_command.execute(download_args) == 0

        download_registry = json.loads((tmp_path / "prepare_registry_download.json").read_text())
        download_metadata_rec = download_registry["datasets"]["SRR999"]["metadata"]

        # Test ParseMetadataCommand on the same XML file
        parse_args = argparse.Namespace(
            metadata_folder=str(metadata_folder),
            metadata_table_file=str(tmp_path / "metadata_table.txt"),
            registry=str(tmp_path / "prepare_registry_parse.json"),
        )

        parse_command = ParseMetadataCommand()
        assert parse_command.execute(parse_args) == 0

        parse_registry = json.loads((tmp_path / "prepare_registry_parse.json").read_text())
        parse_metadata_rec = parse_registry["datasets"]["SRR999"]["metadata"]

        # Verify both commands record identical fields except date and xml. field_keys is
        # derived from _metadata_fields' own source-column mapping (fed a row with every
        # source column populated) rather than hand-written, so a field added to that mapping
        # is automatically covered here too.
        all_source_columns = {
            "Run_Size": "x",
            "Run_MD5": "x",
            "Run_Total_Spots": "1",
            "Run_Total_Bases": "1",
            "Experiment_Library_Strategy": "x",
            "Sample_Scientific_Name": "x",
            "collection_date": "x",
            "Experiment_Library_Layout": "x",
            "Platform": "x",
        }
        field_keys = sorted(_metadata_fields(all_source_columns).keys())
        assert field_keys == sorted(
            key for key in download_metadata_rec if key not in ("date", "xml")
        ), "field_keys must cover every key record_metadata writes besides date and xml"
        for key in field_keys:
            assert download_metadata_rec.get(key) == parse_metadata_rec.get(
                key
            ), f"Mismatch for {key}: download={download_metadata_rec.get(key)} vs parse={parse_metadata_rec.get(key)}"


class TestDownloadTestGenomeCommand:
    """Test DownloadTestGenomeCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = DownloadTestGenomeCommand()
        assert command.name == "download_test_genome"
        assert "download" in command.help.lower()
        assert "test" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = DownloadTestGenomeCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test default output folder
        args = parser.parse_args([])
        assert args.output_folder == "genomes"

        # Test custom output folder
        args = parser.parse_args(["--output-folder", "custom_genomes"])
        assert args.output_folder == "custom_genomes"

    @patch("prepare.cli.commands.test_data.download_test_genome")
    def test_execute(self, mock_command, tmp_path, monkeypatch):
        """Test command execution."""
        monkeypatch.chdir(tmp_path)
        mock_command.return_value = None
        command = DownloadTestGenomeCommand()

        args = argparse.Namespace(output_folder="genomes")
        result = command.execute(args)

        assert result == 0
        mock_command.assert_called_once_with(str(tmp_path.resolve() / "genomes"))


class TestDownloadSraCommand:
    """Test DownloadSraCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = DownloadSraCommand()
        assert command.name == "download_sra"
        assert "download" in command.help.lower()
        assert "sra" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = DownloadSraCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test with required accessions file
        args = parser.parse_args(["--accessions-file", "accessions.txt"])
        assert args.accessions_file == "accessions.txt"
        assert args.fastq_folder == "fastq"
        assert args.num_threads == 4
        assert args.max_workers is None
        assert args.dry_run is False
        assert args.force is False
        assert args.verify_downloads is True
        assert args.redownload_truncated is False

    def test_configure_parser_timeout_flag_defaults_to_none_and_parses(self):
        """--timeout defaults to None (run_secure then falls back to the active setting)."""
        command = DownloadSraCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(["--accessions-file", "a.txt"])
        assert args.timeout is None

        args = parser.parse_args(["--accessions-file", "a.txt", "--timeout", "120"])
        assert args.timeout == 120.0

    def test_configure_parser_verify_downloads_flags(self):
        """--no-verify-downloads flips the default; --redownload-truncated is off by default."""
        command = DownloadSraCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(["--accessions-file", "a.txt", "--no-verify-downloads"])
        assert args.verify_downloads is False
        # The no-op --verify-downloads flag is gone; verification is the default.
        with pytest.raises(SystemExit):
            parser.parse_args(["--accessions-file", "a.txt", "--verify-downloads"])

        args = parser.parse_args(["--accessions-file", "a.txt", "--redownload-truncated"])
        assert args.redownload_truncated is True

    def test_configure_parser_rejects_zero_threads_and_workers(self):
        command = DownloadSraCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)
        with pytest.raises(SystemExit):
            parser.parse_args(["--accessions-file", "a.txt", "--num-threads", "0"])
        with pytest.raises(SystemExit):
            parser.parse_args(["--accessions-file", "a.txt", "--max-workers", "0"])

    def test_explicit_max_workers_is_capped(self, caplog):
        from prepare.core.constants import MAX_CONCURRENT_DOWNLOADS

        with caplog.at_level(logging.WARNING, logger="prepare"):
            with patch("prepare.utils.resources.available_cpus", return_value=256):
                workers = DownloadSraCommand()._resolve_max_workers(argparse.Namespace(num_threads=1, max_workers=64))
        assert workers == MAX_CONCURRENT_DOWNLOADS
        assert any("--max-workers 64" in r.getMessage() and "capped" in r.getMessage() for r in caplog.records)

    def test_configure_parser_prefetch_and_compress_flags(self):
        """--sra-cache, --keep-sra, --no-prefetch, --no-compress default sensibly and parse."""
        command = DownloadSraCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(["--accessions-file", "a.txt"])
        assert args.sra_cache is None
        assert args.use_prefetch is True
        assert args.keep_sra is False
        assert args.compress is True

        args = parser.parse_args(
            [
                "--accessions-file",
                "a.txt",
                "--sra-cache",
                "/tmp/cache",
                "--keep-sra",
                "--no-prefetch",
                "--no-compress",
            ]
        )
        assert args.sra_cache == "/tmp/cache"
        assert args.keep_sra is True
        assert args.use_prefetch is False
        assert args.compress is False

    def test_configure_parser_with_options(self):
        """Test parser with optional arguments."""
        command = DownloadSraCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(
            [
                "--accessions-file",
                "accessions.txt",
                "--fastq-folder",
                "custom_fastq",
                "--num-threads",
                "8",
                "--max-workers",
                "2",
                "--max-downloads",
                "10",
                "--dry-run",
                "--force",
            ]
        )

        assert args.fastq_folder == "custom_fastq"
        assert args.num_threads == 8
        assert args.max_workers == 2
        assert args.max_downloads == 10
        assert args.dry_run is True
        assert args.force is True

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_execute(self, mock_command, _which, tmp_path):
        """Test command execution."""
        mock_command.return_value = {
            "total": 1,
            "to_download": 1,
            "already_downloaded": 0,
            "successful": 1,
            "failed": 0,
            "failed_accessions": [],
        }
        command = DownloadSraCommand()

        args = argparse.Namespace(
            accessions_file="accessions.txt",
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        result = command.execute(args)
        assert result == 0
        call_kwargs = dict(mock_command.call_args.kwargs)
        on_result = call_kwargs.pop("on_result")
        assert callable(on_result)
        # Called through execute (not run), the command makes its own stop token for the run.
        assert isinstance(call_kwargs.pop("stop"), threading.Event)
        # The dict download_sra fills with each attempt's timing, read back by the recorder.
        assert call_kwargs.pop("timings") == {}
        assert call_kwargs == {
            "fastq_folder": str(tmp_path / "fastq"),
            "accessions_file": "accessions.txt",
            "max_downloads": None,
            "dry_run": False,
            "num_threads": 4,
            "max_workers": 4,
            "force": False,
            "max_retries": 1,
            "temp_folder": None,
            "blacklist": None,
            "blacklist_accessions": set(),
            "expected_spots": {},
            "redownload_truncated": False,
            "truncated_accessions": set(),
            "sra_cache": None,
            "use_prefetch": True,
            "keep_sra": False,
            "keep_temp": False,
            "compress": True,
            "lock_wait": 0.0,
            "accept_partial": False,
            "run_sizes": {},
            "min_free_gb": 0.0,  # PREPARE_MIN_FREE_GB=0 from conftest
        }

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_execute_forwards_prefetch_and_compress_flags(self, mock_command, _which, tmp_path):
        """--sra-cache, --keep-sra, --no-prefetch and --no-compress reach download_sra unchanged."""
        mock_command.return_value = {
            "total": 1,
            "to_download": 1,
            "already_downloaded": 0,
            "successful": 1,
            "failed": 0,
            "failed_accessions": [],
        }
        command = DownloadSraCommand()

        args = argparse.Namespace(
            accessions_file="accessions.txt",
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
            sra_cache="/tmp/custom-cache",
            use_prefetch=False,
            keep_sra=True,
            compress=False,
        )

        result = command.execute(args)
        assert result == 0
        call_kwargs = mock_command.call_args.kwargs
        assert call_kwargs["sra_cache"] == "/tmp/custom-cache"
        assert call_kwargs["use_prefetch"] is False
        assert call_kwargs["keep_sra"] is True
        assert call_kwargs["compress"] is False

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_execute_computes_default_max_workers_from_cpu_count(self, mock_command, _which, tmp_path):
        """When --max-workers is not given, it is computed from the CPU count and thread count."""
        mock_command.return_value = {
            "total": 1,
            "to_download": 1,
            "already_downloaded": 0,
            "successful": 1,
            "failed": 0,
            "failed_accessions": [],
        }
        command = DownloadSraCommand()

        args = argparse.Namespace(
            accessions_file="accessions.txt",
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=None,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        with patch("prepare.utils.resources.available_cpus", return_value=8):
            result = command.execute(args)

        assert result == 0
        assert mock_command.call_args.kwargs["max_workers"] == 2

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_execute_warns_when_workers_oversubscribe_cpu(self, mock_command, _which, tmp_path, caplog):
        """An explicit --max-workers x --num-threads exceeding the CPU count logs a warning."""
        mock_command.return_value = {
            "total": 1,
            "to_download": 1,
            "already_downloaded": 0,
            "successful": 1,
            "failed": 0,
            "failed_accessions": [],
        }
        command = DownloadSraCommand()

        args = argparse.Namespace(
            accessions_file="accessions.txt",
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        with patch("prepare.utils.resources.available_cpus", return_value=4):
            with caplog.at_level("WARNING"):
                result = command.execute(args)

        assert result == 0
        assert "exceeds" in caplog.text

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_execute_does_not_warn_about_the_derived_default(self, mock_command, _which, tmp_path, caplog):
        """The CPU-derived default is this tool's own choice; only an explicit --max-workers is warned about."""
        mock_command.return_value = {
            "total": 1,
            "to_download": 1,
            "already_downloaded": 0,
            "successful": 1,
            "failed": 0,
            "failed_accessions": [],
        }
        command = DownloadSraCommand()

        args = argparse.Namespace(
            accessions_file="accessions.txt",
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=None,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        # One worker is already the floor of the derived default, so on a 2-CPU machine a
        # single 4-thread download oversubscribes and there is nothing the user could change.
        with patch("prepare.utils.resources.available_cpus", return_value=2):
            with caplog.at_level("WARNING"):
                result = command.execute(args)

        assert result == 0
        assert mock_command.call_args.kwargs["max_workers"] == 1
        assert "exceeds" not in caplog.text

    @patch("prepare.cli.commands.sra.download_sra")
    def test_execute_dry_run(self, mock_command, tmp_path):
        """Dry-run logs the plan (incl. blacklisted and max-downloads branches) and returns 0."""
        mock_command.return_value = {
            "total": 10,
            "to_download": 5,
            "already_downloaded": 3,
            "blacklisted": 2,
            "successful": 0,
            "failed": 0,
        }
        command = DownloadSraCommand()
        args = argparse.Namespace(
            accessions_file="accessions.txt",
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=5,
            num_threads=4,
            max_workers=4,
            dry_run=True,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=["bl.txt"],
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        result = command.execute(args)
        assert result == 0
        assert not (tmp_path / "prepare_registry.json").exists()

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_execute_with_failures_writes_failed_file(self, mock_command, _which, tmp_path, caplog):
        """Failed downloads return 1; the retry hint is logged and the CLI itself writes no file.

        The "Some downloads failed..." warning is the data layer's own (prepare.data.sra's
        download_sra logs it once already); with download_sra mocked out here, as it is in
        this test, the CLI must not log a second copy of that line itself."""
        failed_file = tmp_path / "fastq" / FAILED_ACCESSIONS_FILE

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR999", False, "Download failed: t")
            return {
                "total": 3,
                "already_downloaded": 1,
                "blacklisted": 1,
                "successful": 1,
                "failed": 1,
                "failed_accessions": ["SRR999"],
            }

        mock_command.side_effect = fake_download_sra
        command = DownloadSraCommand()
        args = argparse.Namespace(
            accessions_file="accessions.txt",
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        with caplog.at_level("INFO"):
            result = command.execute(args)
        assert result == 1
        # The data layer (fake download_sra here) is the only writer of this file; the CLI
        # itself never creates it.
        assert not failed_file.exists()
        assert str(failed_file) in caplog.text
        assert "--accessions-file" in caplog.text
        # download_sra is mocked here, so its own "Some downloads failed..." warning never
        # ran; the CLI must not print a second copy of that line itself.
        assert caplog.text.count("Some downloads failed") == 0

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_execute_reports_disk_full_abort(self, mock_command, _which, tmp_path, caplog):
        """A disk-full abort logs the reason distinctly and returns 1, even though outcomes are recorded."""

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", False, "disk-full: No space left on device")
            kwargs["on_result"]("SRR2", False, "disk-full: not attempted")
            return {
                "total": 2,
                "already_downloaded": 0,
                "blacklisted": 0,
                "successful": 0,
                "failed": 2,
                "failed_accessions": ["SRR1", "SRR2"],
                "aborted": "disk-full",
            }

        mock_command.side_effect = fake_download_sra
        command = DownloadSraCommand()
        args = argparse.Namespace(
            accessions_file="accessions.txt",
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=2,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        with caplog.at_level("ERROR"):
            result = command.execute(args)

        assert result == 1
        assert "disk-full" in caplog.text

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_execute_prepare_error(self, mock_command, _which, tmp_path):
        """A PrepareError from the backend is caught and returns 1."""
        mock_command.side_effect = PrepareError("backend boom")
        command = DownloadSraCommand()
        args = argparse.Namespace(
            accessions_file="accessions.txt",
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        result = command.execute(args)
        assert result == 1

    @patch("prepare.utils.tools.shutil.which", return_value=None)
    @patch("prepare.cli.commands.sra.download_sra")
    def test_missing_fasterq_dump_exits_3(self, mock_download, _which, tmp_path, caplog):
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        args = argparse.Namespace(
            fastq_folder=str(tmp_path / "fastq"),
            accessions_file=str(acc),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )
        assert DownloadSraCommand().execute(args) == 3
        mock_download.assert_not_called()
        assert "fasterq-dump not found on PATH" in caplog.text and "sra-tools" in caplog.text

    @patch("prepare.utils.tools.shutil.which", return_value=None)
    @patch("prepare.cli.commands.sra.download_sra")
    def test_dry_run_skips_tool_check(self, mock_download, _which, tmp_path):
        mock_download.return_value = {
            "total": 1,
            "already_downloaded": 0,
            "blacklisted": 0,
            "to_download": 1,
            "successful": 0,
            "failed": 0,
        }
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        args = argparse.Namespace(
            fastq_folder=str(tmp_path / "fastq"),
            accessions_file=str(acc),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=True,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )
        assert DownloadSraCommand().execute(args) == 0

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra", side_effect=KeyboardInterrupt)
    def test_keyboard_interrupt_returns_130(self, _download, _which, tmp_path, caplog):
        args = argparse.Namespace(
            fastq_folder=str(tmp_path / "fastq"),
            accessions_file=str(tmp_path / "acc.txt"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )
        with caplog.at_level(logging.ERROR):
            assert DownloadSraCommand().execute(args) == 130
        assert "Download interrupted by the user" in caplog.text

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_report_file_lists_every_status(self, mock_download, _which, tmp_path):
        def fake_download_sra(**kwargs):
            kwargs["timings"]["SRR1"] = ("2026-10-01T10:00:00+00:00", 12.5)
            return {
                "total": 4,
                "already_downloaded": 1,
                "blacklisted": 1,
                "successful": 1,
                "failed": 1,
                "failed_accessions": ["SRR2"],
                "results": {"SRR1": "Downloaded 2 files", "SRR2": "Download failed: timeout"},
                "already_downloaded_accessions": ["SRR3"],
                "blacklisted_accessions": ["SRR4"],
            }

        mock_download.side_effect = fake_download_sra
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\nSRR2\nSRR3\nSRR4\n")
        report = tmp_path / "reports" / "download_report.csv"
        args = argparse.Namespace(
            fastq_folder=str(tmp_path / "fastq"),
            accessions_file=str(acc),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=str(report),
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )
        DownloadSraCommand().execute(args)
        assert report.read_text().splitlines() == [
            # The first four columns are unchanged; reason and attempts follow. This fake reports no
            # outcome through on_result, so no attempt is counted.
            "accession,status,message,seconds,reason,attempts",
            "SRR1,downloaded,Downloaded 2 files,12.5,,0",
            "SRR2,failed,Download failed: timeout,,network,0",
            "SRR3,already_present,,,,0",
            "SRR4,blacklisted,,,,0",
        ]

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_recorder_writes_download_timing(self, mock_download, _which, tmp_path):
        """A download's start and seconds reach the registry; a failed attempt is timed too."""
        fastq_folder = tmp_path / "fastq"
        (fastq_folder / "SRR1").mkdir(parents=True)
        (fastq_folder / "SRR1" / "SRR1_1.fastq").write_text("@r\nA\n+\nI\n")

        def fake_download_sra(**kwargs):
            kwargs["timings"]["SRR1"] = ("2026-10-01T10:00:00+00:00", 12.5)
            kwargs["on_result"]("SRR1", True, "Downloaded 1 files")
            kwargs["timings"]["SRR2"] = ("2026-10-01T10:00:05+00:00", 3.0)
            kwargs["on_result"]("SRR2", False, "Download failed: t")
            return {
                "total": 2,
                "successful": 1,
                "failed": 1,
                "failed_accessions": ["SRR2"],
                "results": {"SRR1": "Downloaded 1 files", "SRR2": "Download failed: t"},
            }

        mock_download.side_effect = fake_download_sra
        registry_file = tmp_path / "prepare_registry.json"
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_file),
            data_root=None,
        )
        DownloadSraCommand().execute(args)
        datasets = json.loads(registry_file.read_text())["datasets"]
        assert datasets["SRR1"]["download"]["started"] == "2026-10-01T10:00:00+00:00"
        assert datasets["SRR1"]["download"]["seconds"] == 12.5
        assert datasets["SRR2"]["download"]["seconds"] == 3.0

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_failed_redownload_keeps_the_verdict_of_the_folder_on_disk(self, mock_download, _which, tmp_path):
        """A forced redownload that comes back short leaves the complete folder, and its verdict, in place."""
        fastq_folder = tmp_path / "fastq"
        (fastq_folder / "SRR1").mkdir(parents=True)
        (fastq_folder / "SRR1" / "SRR1_1.fastq").write_text("@r\nA\n+\nI\n")
        registry_file = tmp_path / "prepare_registry.json"
        from prepare.data import registry as reg

        registry = reg.load_registry(registry_file)
        reg.record_download(
            registry,
            "SRR1",
            "downloaded",
            fastq_folder,
            complete={"verdict": "complete", "reads_r1": 1000, "expected_spots": 1000, "ratio": 1.0},
        )
        reg.save_registry(registry)

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", False, "Downloaded 1 files, truncated (10 of 1000 spots)")
            kwargs["on_result"]("SRR2", False, "Downloaded 1 files, truncated (10 of 1000 spots)")
            return {"total": 2, "successful": 0, "failed": 2, "failed_accessions": ["SRR1", "SRR2"], "results": {}}

        mock_download.side_effect = fake_download_sra
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=True,
            max_retries=0,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_file),
            data_root=None,
        )
        DownloadSraCommand().execute(args)
        datasets = json.loads(registry_file.read_text())["datasets"]
        # SRR1's complete folder is still on disk: its verdict stands and its files are kept.
        assert datasets["SRR1"]["download"]["state"] == "failed"
        assert datasets["SRR1"]["download"]["complete"]["verdict"] == "complete"
        assert datasets["SRR1"]["download"]["files"]
        # SRR2 has no folder: the failed attempt's truncated verdict is recorded.
        assert datasets["SRR2"]["download"]["complete"]["verdict"] == "truncated"

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_recorder_clears_download_timing_for_a_store_link(self, mock_download, _which, tmp_path):
        """A dataset linked from the store was not downloaded by this run, so an earlier time is removed."""
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        os.symlink(acc_dir, fastq_folder / "SRR1")
        registry_file = tmp_path / "prepare_registry.json"
        seeded = load_registry(registry_file)
        record_download(seeded, "SRR1", "downloaded", fastq_folder)
        seeded.datasets["SRR1"]["download"].update(started="2026-09-01T10:00:00+00:00", seconds=99.0)
        save_registry(seeded)
        message = "linked from store, 1 files"

        def fake_download_sra(**kwargs):
            kwargs["timings"]["SRR1"] = ("2026-10-01T10:00:00+00:00", 0.2)
            kwargs["on_result"]("SRR1", True, message)
            return {"total": 1, "successful": 1, "failed": 0, "failed_accessions": [], "results": {"SRR1": message}}

        mock_download.side_effect = fake_download_sra
        report = tmp_path / "report.csv"
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=str(report),
            registry=str(registry_file),
            data_root=str(store_root),
        )
        assert DownloadSraCommand().execute(args) == 0
        download = json.loads(registry_file.read_text())["datasets"]["SRR1"]["download"]
        assert download["source"] == "store"
        assert "started" not in download and "seconds" not in download
        assert report.read_text().splitlines()[1] == f'SRR1,downloaded,"{message}",,,0'

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_execute_records_outcomes_in_registry(self, mock_download, _which, tmp_path):
        """Every accession's outcome ends up in the registry; the CLI never touches failed_accessions.txt."""
        fastq_folder = tmp_path / "fastq"
        (fastq_folder / "SRR1").mkdir(parents=True)
        (fastq_folder / "SRR1" / "SRR1_1.fastq").write_text("@r\nA\n+\nI\n")
        registry_path = tmp_path / "prepare_registry.json"
        failed_file = fastq_folder / FAILED_ACCESSIONS_FILE

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", True, "Downloaded 2 files")
            kwargs["on_result"]("SRR2", False, "Download failed: t")
            failed_file.parent.mkdir(parents=True, exist_ok=True)
            failed_file.write_text("SRR2\n")
            return {
                "total": 5,
                "already_downloaded": 1,
                "blacklisted": 1,
                "successful": 1,
                "failed": 1,
                "failed_accessions": ["SRR2"],
                "results": {"SRR1": "Downloaded 2 files", "SRR2": "Download failed: t"},
                "already_downloaded_accessions": ["SRR3"],
                "blacklisted_accessions": ["SRR4"],
                "skipped_accessions": ["SRR5"],
            }

        mock_download.side_effect = fake_download_sra

        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\nSRR2\nSRR3\nSRR4\nSRR5\n")
        args = argparse.Namespace(
            accessions_file=str(acc),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_path),
            data_root=None,
        )

        result = DownloadSraCommand().execute(args)
        assert result == 1

        registry = json.loads(registry_path.read_text())
        datasets = registry["datasets"]

        assert datasets["SRR1"]["download"]["state"] == "downloaded"
        assert datasets["SRR1"]["download"]["files"][0]["path"].endswith("SRR1_1.fastq")
        assert "inferred" not in datasets["SRR1"]["download"]

        assert datasets["SRR2"]["download"]["state"] == "failed"
        assert datasets["SRR2"]["download"]["attempts"] == 1

        assert datasets["SRR3"]["download"]["state"] == "downloaded"
        assert "inferred" not in datasets["SRR3"]["download"]
        # Already present on disk, not a real attempt: attempts stays at 0.
        assert datasets["SRR3"]["download"]["attempts"] == 0

        assert datasets["SRR4"]["download"]["state"] == "skipped"
        assert datasets["SRR4"]["download"]["message"] == "blacklisted"

        assert datasets["SRR5"]["download"]["state"] == "skipped"
        assert datasets["SRR5"]["download"]["message"] == "--max-downloads"

        # The failed_accessions.txt written by the fake download_sra (standing in for the
        # data layer) is untouched by the CLI.
        assert failed_file.read_text() == "SRR2\n"

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_verdict_from_message_is_recorded_in_registry(self, mock_download, _which, tmp_path):
        """The verdict encoded in on_result's message ends up as download.complete in the registry."""
        fastq_folder = tmp_path / "fastq"
        registry_path = tmp_path / "prepare_registry.json"

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", True, "Downloaded 1 files, truncated (300000 of 48000000 spots)")
            return {
                "total": 1,
                "already_downloaded": 0,
                "blacklisted": 0,
                "successful": 1,
                "failed": 0,
                "failed_accessions": [],
            }

        mock_download.side_effect = fake_download_sra
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        args = argparse.Namespace(
            accessions_file=str(acc),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_path),
            data_root=None,
        )

        assert DownloadSraCommand().execute(args) == 0
        datasets = json.loads(registry_path.read_text())["datasets"]
        assert datasets["SRR1"]["download"]["complete"]["verdict"] == "truncated"
        assert datasets["SRR1"]["download"]["complete"]["reads_r1"] == 300000
        assert datasets["SRR1"]["download"]["complete"]["expected_spots"] == 48000000

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_expected_spots_read_from_registry_metadata(self, mock_download, _which, tmp_path):
        """--verify-downloads (the default) reads run_total_spots from the registry into expected_spots."""
        registry_path = tmp_path / "prepare_registry.json"
        seeded = load_registry(registry_path)
        record_metadata(seeded, "SRR1", tmp_path / "SRR1.xml", {"run_total_spots": 12345})
        save_registry(seeded)

        mock_download.return_value = {
            "total": 1,
            "already_downloaded": 0,
            "blacklisted": 0,
            "successful": 0,
            "failed": 0,
            "failed_accessions": [],
        }
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        args = argparse.Namespace(
            accessions_file=str(acc),
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_path),
            data_root=None,
            verify_downloads=True,
        )

        DownloadSraCommand().execute(args)
        assert mock_download.call_args.kwargs["expected_spots"] == {"SRR1": 12345}

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_no_verify_downloads_passes_no_expected_spots(self, mock_download, _which, tmp_path):
        registry_path = tmp_path / "prepare_registry.json"
        seeded = load_registry(registry_path)
        record_metadata(seeded, "SRR1", tmp_path / "SRR1.xml", {"run_total_spots": 12345})
        save_registry(seeded)

        mock_download.return_value = {
            "total": 1,
            "already_downloaded": 0,
            "blacklisted": 0,
            "successful": 0,
            "failed": 0,
            "failed_accessions": [],
        }
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        args = argparse.Namespace(
            accessions_file=str(acc),
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_path),
            data_root=None,
            verify_downloads=False,
        )

        DownloadSraCommand().execute(args)
        assert mock_download.call_args.kwargs["expected_spots"] is None

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_redownload_truncated_computes_set_from_registry_verdicts(self, mock_download, _which, tmp_path):
        registry_path = tmp_path / "prepare_registry.json"
        seeded = load_registry(registry_path)
        record_download(seeded, "SRR1", "downloaded", tmp_path / "fastq")
        seeded.datasets["SRR1"]["download"]["complete"] = {"verdict": "truncated"}
        record_download(seeded, "SRR2", "downloaded", tmp_path / "fastq")
        seeded.datasets["SRR2"]["download"]["complete"] = {"verdict": "complete"}
        save_registry(seeded)

        mock_download.return_value = {
            "total": 2,
            "already_downloaded": 0,
            "blacklisted": 0,
            "successful": 0,
            "failed": 0,
            "failed_accessions": [],
        }
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\nSRR2\n")
        args = argparse.Namespace(
            accessions_file=str(acc),
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_path),
            data_root=None,
            redownload_truncated=True,
        )

        DownloadSraCommand().execute(args)
        assert mock_download.call_args.kwargs["truncated_accessions"] == {"SRR1"}
        assert mock_download.call_args.kwargs["redownload_truncated"] is True

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_blacklisted_accession_keeps_its_downloaded_record(self, mock_download, _which, tmp_path):
        """An accession downloaded earlier and blacklisted later keeps its files and sizes."""
        fastq_folder = tmp_path / "fastq"
        (fastq_folder / "SRR1").mkdir(parents=True)
        (fastq_folder / "SRR1" / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        registry_path = tmp_path / "prepare_registry.json"

        seeded = load_registry(registry_path)
        record_download(seeded, "SRR1", "downloaded", fastq_folder)
        record_exclusion(seeded, "SRR1", "16S amplicon")
        save_registry(seeded)
        before = json.loads(registry_path.read_text())["datasets"]["SRR1"]["download"]

        mock_download.return_value = {
            "total": 1,
            "already_downloaded": 0,
            "blacklisted": 1,
            "successful": 0,
            "failed": 0,
            "failed_accessions": [],
            "blacklisted_accessions": ["SRR1"],
        }
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        args = argparse.Namespace(
            accessions_file=str(acc),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_path),
            data_root=None,
        )

        assert DownloadSraCommand().execute(args) == 0
        after = json.loads(registry_path.read_text())["datasets"]["SRR1"]["download"]
        assert after["state"] == "downloaded"
        assert after["files"] == before["files"]
        assert after["bytes_total"] == before["bytes_total"] > 0

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_concurrent_registry_edit_survives_the_run(self, mock_download, _which, tmp_path):
        """An exclusion written by another process mid-run is not reverted by the download's writes."""
        fastq_folder = tmp_path / "fastq"
        registry_path = tmp_path / "prepare_registry.json"

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", True, "Downloaded 1 file")
            # Stand in for a second terminal running `prepare blacklist --add SRR9`.
            other = load_registry(registry_path)
            record_exclusion(other, "SRR9", "16S amplicon")
            save_registry(other)
            kwargs["on_result"]("SRR2", True, "Downloaded 1 file")
            return {
                "total": 2,
                "already_downloaded": 0,
                "blacklisted": 0,
                "successful": 2,
                "failed": 0,
                "failed_accessions": [],
            }

        mock_download.side_effect = fake_download_sra
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\nSRR2\n")
        args = argparse.Namespace(
            accessions_file=str(acc),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_path),
            data_root=None,
        )

        assert DownloadSraCommand().execute(args) == 0
        datasets = json.loads(registry_path.read_text())["datasets"]
        assert datasets["SRR9"]["exclusion"]["excluded"] is True
        assert datasets["SRR1"]["download"]["state"] == "downloaded"
        assert datasets["SRR2"]["download"]["state"] == "downloaded"

    @patch("prepare.cli.commands.sra.download_sra")
    def test_dry_run_resolves_and_logs_data_root(self, mock_download, tmp_path, caplog):
        """--dry-run --data-root resolves the store root (no fasterq-dump check needed) and logs it."""
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        init_store(store_root)
        mock_download.return_value = {
            "total": 1,
            "to_download": 1,
            "already_downloaded": 0,
            "successful": 0,
            "failed": 0,
        }
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        args = argparse.Namespace(
            accessions_file=str(acc),
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=True,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=str(store_root),
        )

        with caplog.at_level("INFO"):
            result = DownloadSraCommand().execute(args)

        assert result == 0
        assert any(str(store_root.resolve()) in message for message in caplog.messages)

    # ------------------------------------------------------------ shared store

    def test_configure_parser_store_flags(self):
        """--link-mode, --accept-partial and --no-resume-partial default sensibly and parse."""
        parser = argparse.ArgumentParser()
        DownloadSraCommand().configure_parser(parser)

        args = parser.parse_args(["--accessions-file", "a.txt"])
        assert args.link_mode == "auto"
        assert args.accept_partial is False
        assert args.resume_partial is True

        args = parser.parse_args(
            ["--accessions-file", "a.txt", "--link-mode", "copy", "--accept-partial", "--no-resume-partial"]
        )
        assert args.link_mode == "copy"
        assert args.accept_partial is True
        assert args.resume_partial is False

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_store_paths_and_link_options_reach_download_sra(self, mock_download, _which, tmp_path):
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        init_store(store_root)
        mock_download.return_value = {
            "total": 1,
            "to_download": 1,
            "already_downloaded": 0,
            "successful": 1,
            "failed": 0,
            "failed_accessions": [],
        }
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        args = argparse.Namespace(
            accessions_file=str(acc),
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=str(store_root),
            link_mode="relative",
            accept_partial=True,
            resume_partial=False,
        )

        assert DownloadSraCommand().execute(args) == 0

        kwargs = mock_download.call_args.kwargs
        assert kwargs["store"].root == store_root.resolve()
        assert kwargs["link_mode"] == "relative"
        assert kwargs["accept_partial"] is True
        assert kwargs["resume_partial"] is False
        assert Path(kwargs["store_metadata"][0]) == tmp_path / "metadata"
        assert Path(kwargs["store_metadata"][1]) == store_root.resolve() / "metadata"

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_no_store_leaves_the_download_call_unchanged(self, mock_download, _which, tmp_path, monkeypatch):
        monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
        monkeypatch.delenv("PREPARE_DATA", raising=False)
        mock_download.return_value = {
            "total": 1,
            "to_download": 1,
            "already_downloaded": 0,
            "successful": 1,
            "failed": 0,
            "failed_accessions": [],
        }
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        args = argparse.Namespace(
            accessions_file=str(acc),
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        assert DownloadSraCommand().execute(args) == 0
        assert "store" not in mock_download.call_args.kwargs

    @pytest.mark.parametrize(
        "message,attempts",
        [("linked from store, 1 files", 0), ("Downloaded 1 files, complete (1 of 1 spots); stored", 1)],
    )
    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_store_backed_results_record_source_and_linked_list(
        self, mock_download, _which, tmp_path, message, attempts
    ):
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        os.symlink(acc_dir, fastq_folder / "SRR1")

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", True, message)
            return {
                "total": 1,
                "to_download": 1,
                "already_downloaded": 0,
                "blacklisted": 0,
                "successful": 1,
                "failed": 0,
                "failed_accessions": [],
                "results": {"SRR1": message},
            }

        mock_download.side_effect = fake_download_sra
        registry_file = tmp_path / "prepare_registry.json"
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_file),
            data_root=str(store_root),
        )

        assert DownloadSraCommand().execute(args) == 0

        written = json.loads(registry_file.read_text())
        download = written["datasets"]["SRR1"]["download"]
        assert download["state"] == "downloaded"
        assert download["source"] == "store"
        assert download["store_name"] == "SRR1"
        assert download["attempts"] == attempts
        assert written["store"]["linked"] == ["SRR1"]

    @pytest.mark.parametrize(
        "previous, sidecar_verdict, sidecar_ratio, expected",
        [
            # complete -> complete: the recorded read count is kept (audit deferred S7-7).
            (
                {"verdict": "complete", "reads_r1": 48000000, "expected_spots": 48000000, "ratio": 1.0},
                "complete",
                1.0,
                {"verdict": "complete", "reads_r1": 48000000, "expected_spots": None, "ratio": 1.0},
            ),
            # truncated project copy relinked to a complete store copy: the truncated count
            # describes other files, so nothing is carried; ratio comes from the new block only.
            (
                {"verdict": "truncated", "reads_r1": 5, "expected_spots": 48000000, "ratio": 0.0001},
                "complete",
                1.0,
                {"verdict": "complete", "reads_r1": None, "expected_spots": None, "ratio": 1.0},
            ),
            # complete -> truncated: nothing carried either.
            (
                {"verdict": "complete", "reads_r1": 48000000, "expected_spots": 48000000, "ratio": 1.0},
                "truncated",
                None,
                {"verdict": "truncated", "reads_r1": None, "expected_spots": None, "ratio": None},
            ),
        ],
        ids=["complete-to-complete", "truncated-to-complete", "complete-to-truncated"],
    )
    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_relink_carries_read_count_only_between_equal_verdicts(
        self, mock_download, _which, tmp_path, previous, sidecar_verdict, sidecar_ratio, expected
    ):
        """A relink to a store copy takes its verdict from the sidecar. A sidecar with no read
        count of its own keeps the registry's previous count only when the previous verdict
        equals the new one; ratio and expected_spots are never carried over, so a relink never
        mixes two downloads' verdict blocks."""
        from prepare.store.layout import init_store, sidecar_path
        from prepare.store.sidecar import Sidecar, write_sidecar

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        write_sidecar(
            sidecar_path(paths, "SRR1"),
            Sidecar(
                accession="SRR1",
                state="complete",
                reads_per_mate=None,
                completeness={"method": "spots", "ratio": sidecar_ratio, "verdict": sidecar_verdict},
            ),
        )

        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        os.symlink(acc_dir, fastq_folder / "SRR1")

        registry_file = tmp_path / "prepare_registry.json"
        seeded = load_registry(registry_file)
        record_download(seeded, "SRR1", "downloaded", fastq_folder)
        seeded.datasets["SRR1"]["download"]["complete"] = dict(previous)
        save_registry(seeded)

        message = "linked from store, 1 files"

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", True, message)
            return {
                "total": 1,
                "to_download": 1,
                "already_downloaded": 0,
                "blacklisted": 0,
                "successful": 1,
                "failed": 0,
                "failed_accessions": [],
                "results": {"SRR1": message},
            }

        mock_download.side_effect = fake_download_sra
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_file),
            data_root=str(store_root),
            redownload_truncated=True,
        )

        assert DownloadSraCommand().execute(args) == 0

        written = json.loads(registry_file.read_text())
        complete = written["datasets"]["SRR1"]["download"]["complete"]
        assert {key: complete.get(key) for key in expected} == expected

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_a_plain_download_is_not_recorded_as_store_backed(self, mock_download, _which, tmp_path):
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        init_store(store_root)
        fastq_folder = tmp_path / "fastq"
        (fastq_folder / "SRR1").mkdir(parents=True)
        (fastq_folder / "SRR1" / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", True, "Downloaded 1 files, unverified")
            return {
                "total": 1,
                "to_download": 1,
                "already_downloaded": 0,
                "blacklisted": 0,
                "successful": 1,
                "failed": 0,
                "failed_accessions": [],
                "results": {"SRR1": "Downloaded 1 files, unverified"},
            }

        mock_download.side_effect = fake_download_sra
        registry_file = tmp_path / "prepare_registry.json"
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_file),
            data_root=str(store_root),
        )

        assert DownloadSraCommand().execute(args) == 0

        written = json.loads(registry_file.read_text())
        assert "source" not in written["datasets"]["SRR1"]["download"]
        assert written["store"].get("linked", []) == []

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_an_already_linked_dataset_is_recorded_as_store_backed(self, mock_download, _which, tmp_path):
        """A link left by an earlier run is recorded with its source, not as a plain download."""
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        os.symlink(acc_dir, fastq_folder / "SRR1")

        mock_download.return_value = {
            "total": 1,
            "to_download": 0,
            "already_downloaded": 1,
            "blacklisted": 0,
            "successful": 0,
            "failed": 0,
            "failed_accessions": [],
            "results": {},
            "already_downloaded_accessions": ["SRR1"],
        }
        registry_file = tmp_path / "prepare_registry.json"
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_file),
            data_root=str(store_root),
        )

        assert DownloadSraCommand().execute(args) == 0

        written = json.loads(registry_file.read_text())
        assert written["datasets"]["SRR1"]["download"]["source"] == "store"
        assert written["store"]["linked"] == ["SRR1"]

    # -------------------------------------------------------- store usage catalogue

    @staticmethod
    def _seed_project(registry_path, project_id="proj1"):
        from prepare.data.registry import load_registry as _load, save_registry as _save

        registry = _load(registry_path)
        registry.project = {"id": project_id, "name": "demo", "path": str(registry_path.parent), "created": "now"}
        _save(registry)

    @staticmethod
    def _usage_rows(store_root):
        from prepare.store.catalog import Catalog
        from prepare.store.layout import store_paths

        with Catalog(store_paths(store_root)) as catalog:
            catalog.migrate()
            return [
                dict(r)
                for r in catalog.conn.execute("SELECT accession, project_id, genome_id, stage FROM usage").fetchall()
            ]

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_fresh_store_download_records_usage_stage_downloaded(self, mock_download, _which, tmp_path):
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        init_store(store_root)
        registry_file = tmp_path / "prepare_registry.json"
        self._seed_project(registry_file)
        message = "Downloaded 1 files, complete (1 of 1 spots); stored"

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", True, message)
            return {
                "total": 1,
                "already_downloaded": 0,
                "blacklisted": 0,
                "successful": 1,
                "failed": 0,
                "failed_accessions": [],
                "results": {"SRR1": message},
            }

        mock_download.side_effect = fake_download_sra
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_file),
            data_root=str(store_root),
        )

        assert DownloadSraCommand().execute(args) == 0

        rows = self._usage_rows(store_root)
        assert rows == [{"accession": "SRR1", "project_id": "proj1", "genome_id": "", "stage": "downloaded"}]

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_linked_from_store_records_usage_stage_linked(self, mock_download, _which, tmp_path):
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        init_store(store_root)
        registry_file = tmp_path / "prepare_registry.json"
        self._seed_project(registry_file)
        message = "linked from store, 1 files"

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", True, message)
            return {
                "total": 1,
                "already_downloaded": 0,
                "blacklisted": 0,
                "successful": 1,
                "failed": 0,
                "failed_accessions": [],
                "results": {"SRR1": message},
            }

        mock_download.side_effect = fake_download_sra
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_file),
            data_root=str(store_root),
        )

        assert DownloadSraCommand().execute(args) == 0

        rows = self._usage_rows(store_root)
        assert rows == [{"accession": "SRR1", "project_id": "proj1", "genome_id": "", "stage": "linked"}]

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_retry_prefixed_store_link_is_recorded_as_a_link(self, mock_download, _which, tmp_path):
        from prepare.data.registry import load_registry
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        init_store(store_root)
        registry_file = tmp_path / "prepare_registry.json"
        self._seed_project(registry_file)
        message = "Retry 1: linked from store, 1 files"

        def fake_download_sra(**kwargs):
            kwargs["timings"]["SRR1"] = ("2026-01-01T00:00:00Z", 3.0)
            kwargs["on_result"]("SRR1", True, message)
            return {
                "total": 1,
                "already_downloaded": 0,
                "blacklisted": 0,
                "successful": 1,
                "failed": 0,
                "failed_accessions": [],
                "results": {"SRR1": message},
            }

        mock_download.side_effect = fake_download_sra
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_file),
            data_root=str(store_root),
        )

        assert DownloadSraCommand().execute(args) == 0

        rows = self._usage_rows(store_root)
        assert rows == [{"accession": "SRR1", "project_id": "proj1", "genome_id": "", "stage": "linked"}]
        block = rb.download_block(load_registry(str(registry_file)), "SRR1")
        assert (block.source, block.attempts, block.message) == ("store", 0, message)
        assert block.seconds is None

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_already_downloaded_store_backed_records_usage_stage_linked(self, mock_download, _which, tmp_path):
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        fastq_folder = tmp_path / "fastq"
        fastq_folder.mkdir()
        os.symlink(acc_dir, fastq_folder / "SRR1")

        registry_file = tmp_path / "prepare_registry.json"
        self._seed_project(registry_file)

        mock_download.return_value = {
            "total": 1,
            "to_download": 0,
            "already_downloaded": 1,
            "blacklisted": 0,
            "successful": 0,
            "failed": 0,
            "failed_accessions": [],
            "results": {},
            "already_downloaded_accessions": ["SRR1"],
        }
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_file),
            data_root=str(store_root),
        )

        assert DownloadSraCommand().execute(args) == 0

        rows = self._usage_rows(store_root)
        assert rows == [{"accession": "SRR1", "project_id": "proj1", "genome_id": "", "stage": "linked"}]

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_catalog_failure_leaves_download_outcome_unchanged(self, mock_download, _which, tmp_path):
        """A broken catalogue write never changes the download's exit code or registry record."""
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        init_store(store_root)
        registry_file = tmp_path / "prepare_registry.json"
        self._seed_project(registry_file)
        message = "Downloaded 1 files, complete (1 of 1 spots); stored"

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", True, message)
            return {
                "total": 1,
                "already_downloaded": 0,
                "blacklisted": 0,
                "successful": 1,
                "failed": 0,
                "failed_accessions": [],
                "results": {"SRR1": message},
            }

        mock_download.side_effect = fake_download_sra
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_file),
            data_root=str(store_root),
        )

        with patch("prepare.store.usage.catalog_write", side_effect=RuntimeError("catalogue is locked")):
            result = DownloadSraCommand().execute(args)

        assert result == 0
        written = json.loads(registry_file.read_text())
        download = written["datasets"]["SRR1"]["download"]
        assert download["state"] == "downloaded"
        assert download["source"] == "store"
        assert download["store_name"] == "SRR1"
        assert written["store"]["linked"] == ["SRR1"]

    # ------------------------------------------------------- download summary logging

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_cli_logs_no_second_download_summary(self, mock_download, _which, tmp_path, caplog):
        """The data layer prints the run summary; the CLI must not print a second block whose
        "Successfully downloaded" line counted store links as downloads (audit S5-6)."""
        mock_download.return_value = {
            "total": 2,
            "already_downloaded": 0,
            "blacklisted": 0,
            "successful": 2,
            "failed": 0,
            "failed_accessions": [],
            "results": {
                "SRR1": "linked from store, 1 files",
                "SRR2": "Downloaded 1 files, complete (1 of 1 spots); stored",
            },
        }
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        with caplog.at_level("INFO"):
            assert DownloadSraCommand().execute(args) == 0

        assert "Download summary:" not in caplog.text
        assert "Successfully downloaded" not in caplog.text

    # ------------------------------------------------------- transient bytes warning

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_warns_when_transient_bytes_exceed_threshold(self, mock_download, _which, tmp_path, caplog, monkeypatch):
        """A kept .sra-cache bigger than the (patched, small) threshold is named in a warning."""
        monkeypatch.setattr("prepare.cli.commands.sra.TRANSIENT_BYTES_WARN_THRESHOLD", 10)
        mock_download.return_value = {
            "total": 1,
            "already_downloaded": 0,
            "blacklisted": 0,
            "successful": 1,
            "failed": 0,
            "failed_accessions": [],
        }
        fastq_folder = tmp_path / "fastq"
        cache = fastq_folder / ".sra-cache"
        cache.mkdir(parents=True)
        (cache / "SRR1.sra").write_bytes(b"x" * 100)

        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        with caplog.at_level("WARNING"):
            assert DownloadSraCommand().execute(args) == 0

        assert str(fastq_folder) in caplog.text
        assert "100 bytes" in caplog.text

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_no_warning_when_transient_bytes_under_threshold(self, mock_download, _which, tmp_path, caplog):
        """The default 1 GB threshold is not tripped by a small leftover cache file."""
        mock_download.return_value = {
            "total": 1,
            "already_downloaded": 0,
            "blacklisted": 0,
            "successful": 1,
            "failed": 0,
            "failed_accessions": [],
        }
        fastq_folder = tmp_path / "fastq"
        cache = fastq_folder / ".sra-cache"
        cache.mkdir(parents=True)
        (cache / "SRR1.sra").write_bytes(b"x" * 100)

        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(fastq_folder),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

        with caplog.at_level("WARNING"):
            assert DownloadSraCommand().execute(args) == 0

        assert "transient" not in caplog.text.lower()


class TestSingleSampleCommand:
    """Test SingleSampleCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = SingleSampleCommand()
        assert command.name == "single_sample"
        assert "single" in command.help.lower()
        assert "sample" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = SingleSampleCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test with required arguments
        args = parser.parse_args(["--summary-column", "test_genome", "--metadata-column", "organism"])

        assert args.summary_column == "test_genome"
        assert args.metadata_column == "organism"
        assert args.summary_file == "parsed_containment.txt"
        assert args.metadata_file is None
        assert args.threshold == 0.1
        assert args.top_n == 100

    def test_configure_parser_with_options(self):
        """Test parser with optional arguments."""
        command = SingleSampleCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(
            [
                "--summary-column",
                "test_genome",
                "--metadata-column",
                "organism",
                "--summary-file",
                "custom_summary.txt",
                "--metadata-file",
                "custom_metadata.txt",
                "--threshold",
                "0.5",
                "--top-n",
                "50",
            ]
        )

        assert args.summary_file == "custom_summary.txt"
        assert args.metadata_file == "custom_metadata.txt"
        assert args.threshold == 0.5
        assert args.top_n == 50

    @patch("prepare.cli.commands.samples.resolve_metadata_table")
    @patch("prepare.cli.commands.samples.count_single_sample")
    def test_execute(self, mock_command, mock_resolve):
        """Test command execution."""
        mock_command.return_value = {"organism1": 5, "organism2": 3}
        mock_resolve.return_value = Path("metadata.txt")
        command = SingleSampleCommand()

        args = argparse.Namespace(
            summary_column="test_genome",
            metadata_column="organism",
            summary_file="summary.txt",
            metadata_file="metadata.txt",
            threshold=0.1,
            top_n=100,
        )

        result = command.execute(args)
        assert result == 0
        mock_resolve.assert_called_once_with("metadata.txt")
        mock_command.assert_called_once_with(
            summary_file="summary.txt",
            metadata_file="metadata.txt",
            summary_column="test_genome",
            metadata_column="organism",
            threshold=0.1,
            top_n=100,
        )


class TestParseMetadataCommand:
    """Test ParseMetadataCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = ParseMetadataCommand()
        assert command.name == "parse_metadata"
        assert "parse" in command.help.lower()
        assert "metadata" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = ParseMetadataCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test default values
        args = parser.parse_args([])
        assert args.metadata_folder == "metadata"
        assert args.metadata_table_file == "metadata_table.txt"

        # Test custom values
        args = parser.parse_args(["--metadata-folder", "custom_metadata", "--metadata-table-file", "custom_table.txt"])
        assert args.metadata_folder == "custom_metadata"
        assert args.metadata_table_file == "custom_table.txt"

    @patch("prepare.cli.commands.metadata.parse_metadata")
    def test_execute(self, mock_command, tmp_path):
        """Test command execution."""
        mock_command.return_value = pd.DataFrame()
        command = ParseMetadataCommand()

        args = argparse.Namespace(
            metadata_folder="metadata",
            metadata_table_file="metadata_table.txt",
            registry=str(tmp_path / "prepare_registry.json"),
        )

        result = command.execute(args)
        assert result == 0
        mock_command.assert_called_once_with(str(tmp_path.resolve() / "metadata"), "metadata_table.txt")

    @patch("prepare.cli.commands.metadata.parse_metadata")
    def test_execute_resolves_the_project_root_once(self, mock_parse, tmp_path):
        """record_metadata gets the root from the command, not one project_root call per row."""
        import prepare.data.registry as registry_mod

        mock_parse.return_value = pd.DataFrame({"Run_ID": [f"SRR{i}" for i in range(50)]})
        args = argparse.Namespace(
            metadata_folder=str(tmp_path / "metadata"),
            metadata_table_file=str(tmp_path / "metadata_table.txt"),
            registry=str(tmp_path / "prepare_registry.json"),
        )
        with patch.object(registry_mod, "project_root", wraps=registry_mod.project_root) as spy:
            assert ParseMetadataCommand().execute(args) == 0
        assert spy.call_count == 0
        registry = load_registry(args.registry)
        assert registry.datasets["SRR49"]["metadata"]["xml"] == "metadata/SRR49_metadata.xml"

    def test_execute_records_metadata_in_registry(self, tmp_path):
        """Parsing a metadata folder records run_size/run_md5 for each accession."""
        metadata_folder = tmp_path / "metadata"
        metadata_folder.mkdir()
        xml_content = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE_SET>
            <EXPERIMENT_PACKAGE>
                <SAMPLE>
                    <IDENTIFIERS>
                        <PRIMARY_ID>SAMN123</PRIMARY_ID>
                    </IDENTIFIERS>
                    <SAMPLE_NAME>
                        <SCIENTIFIC_NAME>Escherichia coli</SCIENTIFIC_NAME>
                    </SAMPLE_NAME>
                </SAMPLE>
                <RUN_SET>
                    <RUN>
                        <IDENTIFIERS>
                            <PRIMARY_ID>SRR123</PRIMARY_ID>
                        </IDENTIFIERS>
                        <size>12345</size>
                        <md5>abcdef0123456789</md5>
                    </RUN>
                </RUN_SET>
            </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""
        (metadata_folder / "SRR123_metadata.xml").write_text(xml_content)

        command = ParseMetadataCommand()
        args = argparse.Namespace(
            metadata_folder=str(metadata_folder),
            metadata_table_file=str(tmp_path / "metadata_table.txt"),
            registry=str(tmp_path / "prepare_registry.json"),
        )

        result = command.execute(args)
        assert result == 0

        registry = json.loads((tmp_path / "prepare_registry.json").read_text())
        metadata = registry["datasets"]["SRR123"]["metadata"]
        assert metadata["run_size"] == "12345"
        assert metadata["run_md5"] == "abcdef0123456789"

    def test_execute_records_ncbi_attribute_spots_in_registry(self, tmp_path):
        """Real NCBI efetch XML holds spots/bases/md5 as attributes; parse_metadata must record them."""
        metadata_folder = tmp_path / "metadata"
        metadata_folder.mkdir()
        xml_content = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE_SET>
            <EXPERIMENT_PACKAGE>
                <EXPERIMENT>
                    <IDENTIFIERS>
                        <PRIMARY_ID>EXP1</PRIMARY_ID>
                    </IDENTIFIERS>
                    <LIBRARY_DESCRIPTOR>
                        <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
                        <LIBRARY_LAYOUT>
                            <PAIRED/>
                        </LIBRARY_LAYOUT>
                    </LIBRARY_DESCRIPTOR>
                    <PLATFORM>
                        <ILLUMINA>
                            <INSTRUMENT_MODEL>Illumina HiSeq 2500</INSTRUMENT_MODEL>
                        </ILLUMINA>
                    </PLATFORM>
                </EXPERIMENT>
                <SAMPLE>
                    <IDENTIFIERS>
                        <PRIMARY_ID>SAMN123</PRIMARY_ID>
                    </IDENTIFIERS>
                </SAMPLE>
                <RUN_SET>
                    <RUN accession="SRR1" total_spots="47964651" total_bases="14389395300" size="4744553813">
                        <IDENTIFIERS>
                            <PRIMARY_ID>SRR1</PRIMARY_ID>
                        </IDENTIFIERS>
                        <SRAFiles>
                            <SRAFile filename="SRR1" md5="abc" semantic_name="run"/>
                        </SRAFiles>
                    </RUN>
                </RUN_SET>
            </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""
        (metadata_folder / "SRR1_metadata.xml").write_text(xml_content)

        command = ParseMetadataCommand()
        args = argparse.Namespace(
            metadata_folder=str(metadata_folder),
            metadata_table_file=str(tmp_path / "metadata_table.txt"),
            registry=str(tmp_path / "prepare_registry.json"),
        )

        result = command.execute(args)
        assert result == 0

        registry = json.loads((tmp_path / "prepare_registry.json").read_text())
        metadata = registry["datasets"]["SRR1"]["metadata"]
        assert metadata["run_total_spots"] == 47964651
        assert metadata["run_total_bases"] == 14389395300
        assert metadata["run_size"] == "4744553813"
        assert metadata["run_md5"] == "abc"
        assert metadata["library_layout"] == "PAIRED"
        assert metadata["platform"] == "ILLUMINA"
        assert metadata["library_strategy"] == "WGS"


class TestCountMetadataCommand:
    """Test CountMetadataCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = CountMetadataCommand()
        assert command.name == "count_metadata"
        assert "count" in command.help.lower()
        assert "metadata" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = CountMetadataCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test with required metadata column
        args = parser.parse_args(["--metadata-column", "organism"])
        assert args.metadata_column == "organism"
        assert args.summary_file == "parsed_containment.txt"
        assert args.metadata_file is None
        assert args.threshold == 0.5
        assert args.output_file == "metadata_counts.txt"
        assert args.stat_file is None

    def test_configure_parser_with_options(self):
        """Test parser with optional arguments."""
        command = CountMetadataCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(
            [
                "--metadata-column",
                "organism",
                "--summary-file",
                "custom_summary.txt",
                "--threshold",
                "0.8",
                "--output-file",
                "custom_counts.txt",
                "--stat-file",
                "statistics.txt",
            ]
        )

        assert args.summary_file == "custom_summary.txt"
        assert args.threshold == 0.8
        assert args.output_file == "custom_counts.txt"
        assert args.stat_file == "statistics.txt"

    @patch("prepare.cli.commands.metadata.resolve_metadata_table")
    @patch("prepare.cli.commands.metadata.count_metadata")
    def test_execute(self, mock_command, mock_resolve):
        """Test command execution."""
        mock_command.return_value = None
        mock_resolve.return_value = Path("metadata.txt")
        command = CountMetadataCommand()

        args = argparse.Namespace(
            metadata_column="organism",
            summary_file="summary.txt",
            metadata_file="metadata.txt",
            threshold=0.5,
            output_file="counts.txt",
            stat_file=None,
        )

        result = command.execute(args)
        assert result == 0
        mock_resolve.assert_called_once_with("metadata.txt")
        mock_command.assert_called_once_with(
            summary_file="summary.txt",
            metadata_file="metadata.txt",
            metadata_column="organism",
            threshold=0.5,
            output_file="counts.txt",
            stat_file=None,
        )


class TestPlotContainmentCommand:
    """Test PlotContainmentCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = PlotContainmentCommand()
        assert command.name == "plot_containment"
        assert "plot" in command.help.lower()
        assert "containment" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = PlotContainmentCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test with required file path
        args = parser.parse_args(["--file-path", "containment.txt"])
        assert args.file_path == "containment.txt"
        assert args.column == "max_containment"
        assert args.plot_type == "rank"
        assert args.title is None
        assert args.colors is None
        assert args.save_format == "png"

    def test_configure_parser_with_options(self):
        """Test parser with optional arguments."""
        command = PlotContainmentCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(
            [
                "--file-path",
                "containment.txt",
                "--column",
                "custom_column",
                "--plot-type",
                "histogram",
                "--title",
                "Custom Title",
                "--colors",
                "red",
                "--save-format",
                "png",
                "--threshold",
                "0.5",
                "--show-title",
            ]
        )

        assert args.column == "custom_column"
        assert args.plot_type == "histogram"
        assert args.title == "Custom Title"
        assert args.colors == "red"
        assert args.save_format == "png"
        assert args.threshold == 0.5
        assert args.show_title is True

    @patch("prepare.cli.commands.containment.viz_plot_containment")
    def test_execute(self, mock_command):
        """Test command execution."""
        import matplotlib.pyplot as plt

        mock_command.return_value = plt.figure()
        command = PlotContainmentCommand()

        args = argparse.Namespace(
            file_path="containment.txt",
            column="max_containment",
            plot_type="rank",
            title=None,
            colors=None,
            save_format=None,
            threshold=None,
            show_title=False,
        )

        result = command.execute(args)
        assert result == 0
        mock_command.assert_called_once_with(
            file_path="containment.txt",
            column="max_containment",
            title=None,
            colors=None,
            show_title=False,
            save_format=None,
            threshold=None,
            plot_type="rank",
        )

    @patch("prepare.cli.commands.containment.viz_plot_containment")
    def test_execute_logs_the_saved_path_and_a_next_hint(self, mock_command, caplog):
        """The saved path the command logs must be the one _save_plot_if_needed actually
        writes to; both are built from the shared plot_output_path helper so there is one
        source of truth for the naming formula, not two copies that can drift apart."""
        import matplotlib.pyplot as plt

        mock_command.return_value = plt.figure()
        command = PlotContainmentCommand()

        args = argparse.Namespace(
            file_path="parsed_containment.txt",
            column="max_containment",
            plot_type="rank",
            title=None,
            colors=None,
            save_format="png",
            threshold=None,
            show_title=False,
        )

        with caplog.at_level("INFO"):
            result = command.execute(args)

        assert result == 0
        assert "parsed_containment_rank_max_containment.png" in caplog.text
        assert "Next:" in caplog.text

    @patch("prepare.cli.commands.containment.viz_plot_containment")
    def test_title_alone_shows_the_title(self, mock_command):
        mock_command.return_value = None
        command = PlotContainmentCommand()
        args = argparse.Namespace(
            file_path="c.txt",
            column="max_containment",
            plot_type="rank",
            title="Mine",
            colors=None,
            save_format="png",
            threshold=None,
            show_title=False,
        )
        assert command.execute(args) == 0
        assert mock_command.call_args.kwargs["show_title"] is True


class TestPlotMetadataCountsCommand:
    """Test PlotMetadataCountsCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = PlotMetadataCountsCommand()
        assert command.name == "plot_metadata_counts"
        assert "plot" in command.help.lower()
        assert "metadata" in command.help.lower()
        assert "counts" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = PlotMetadataCountsCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test with required file path
        args = parser.parse_args(["--file-path", "counts.txt"])
        assert args.file_path == "counts.txt"
        assert args.plot_type == "bar"
        assert args.title is None
        assert args.colors is None
        assert args.save_format is None
        assert args.show_title is False

    def test_configure_parser_with_options(self):
        """Test parser with optional arguments."""
        command = PlotMetadataCountsCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(
            [
                "--file-path",
                "counts.txt",
                "--plot-type",
                "pie",
                "--title",
                "Metadata Counts",
                "--colors",
                "viridis",
                "--save-format",
                "pdf",
                "--show-title",
            ]
        )

        assert args.plot_type == "pie"
        assert args.title == "Metadata Counts"
        assert args.colors == "viridis"
        assert args.save_format == "pdf"
        assert args.show_title is True

    @patch("prepare.cli.commands.metadata.plot_metadata_counts")
    def test_execute(self, mock_command):
        """Test command execution."""
        import matplotlib.pyplot as plt

        mock_command.return_value = plt.figure()
        command = PlotMetadataCountsCommand()

        args = argparse.Namespace(
            file_path="counts.txt", plot_type="bar", title=None, colors=None, save_format=None, show_title=False
        )

        result = command.execute(args)
        assert result == 0
        mock_command.assert_called_once_with(
            file_path="counts.txt", title=None, plot_type="bar", colors=None, show_title=False, save_format=None
        )

    @patch("prepare.cli.commands.metadata.plot_metadata_counts")
    def test_title_alone_shows_the_title(self, mock_command):
        command = PlotMetadataCountsCommand()
        args = argparse.Namespace(
            file_path="counts.txt", plot_type="bar", title="Mine", colors=None, save_format="png", show_title=False
        )
        assert command.execute(args) == 0
        assert mock_command.call_args.kwargs["show_title"] is True

    @patch("prepare.cli.commands.metadata.plot_metadata_counts")
    def test_warns_when_nothing_is_saved(self, mock_command, caplog):
        command = PlotMetadataCountsCommand()
        args = argparse.Namespace(
            file_path="counts.txt", plot_type="bar", title=None, colors=None, save_format=None, show_title=False
        )
        with caplog.at_level("WARNING"):
            assert command.execute(args) == 0
        assert "no --save-format given; nothing was written" in caplog.text


class TestDownloadSraMintsAProjectIdentity:
    """A project that reached the store without store_init still records its usage."""

    @staticmethod
    def _usage_rows(store_root):
        from prepare.store.catalog import Catalog
        from prepare.store.layout import store_paths

        with Catalog(store_paths(store_root)) as catalog:
            catalog.migrate()
            return [dict(r) for r in catalog.conn.execute("SELECT accession, project_id, stage FROM usage").fetchall()]

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_a_project_without_store_init_still_records_usage(self, mock_download, _which, tmp_path, monkeypatch):
        from prepare.data.registry import load_registry
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        init_store(store_root)
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry_file = project_dir / "prepare_registry.json"
        message = "Downloaded 1 files, complete (1 of 1 spots); stored"

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", True, message)
            return {
                "total": 1,
                "already_downloaded": 0,
                "blacklisted": 0,
                "successful": 1,
                "failed": 0,
                "failed_accessions": [],
                "results": {"SRR1": message},
            }

        mock_download.side_effect = fake_download_sra
        args = argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(project_dir / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(registry_file),
            data_root=str(store_root),
        )

        assert DownloadSraCommand().execute(args) == 0

        registry = load_registry(registry_file)
        assert registry.project["id"]
        rows = self._usage_rows(store_root)
        assert [(r["accession"], r["stage"]) for r in rows] == [("SRR1", "downloaded")]
        assert rows[0]["project_id"] == registry.project["id"]


class TestDownloadMetadataSharesWithTheStore:
    """NCBI's spot count is what makes a store dataset verifiable, so the store gets a copy."""

    def _args(self, tmp_path, store_root=None):
        return argparse.Namespace(
            email="test@example.com",
            matches_folder=str(tmp_path / "matches"),
            metadata_folder=str(tmp_path / "metadata"),
            threshold=0.0,
            dry_run=False,
            accessions_file=None,
            api_key=None,
            batch_size=200,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=str(store_root) if store_root else None,
        )

    def test_fetched_xml_is_copied_into_the_store(self, tmp_path):
        from prepare.store.layout import init_store

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        metadata_dir = tmp_path / "metadata"
        metadata_dir.mkdir()
        xml = metadata_dir / "SRR1_metadata.xml"
        xml.write_text('<RunSet><RUN total_spots="10"/></RunSet>')

        with patch("prepare.cli.commands.metadata.download_metadata", return_value={"SRR1": xml}):
            rc = DownloadMetadataCommand().execute(self._args(tmp_path, store_root))

        assert rc == 0
        assert (paths.metadata / "SRR1_metadata.xml").is_file()
        # The project keeps its own copy too.
        assert xml.is_file()

    def test_without_a_store_nothing_is_copied(self, tmp_path):
        metadata_dir = tmp_path / "metadata"
        metadata_dir.mkdir()
        xml = metadata_dir / "SRR1_metadata.xml"
        xml.write_text('<RunSet><RUN total_spots="10"/></RunSet>')

        with patch("prepare.cli.commands.metadata.download_metadata", return_value={"SRR1": xml}):
            assert DownloadMetadataCommand().execute(self._args(tmp_path)) == 0


class TestDownloadSraRegistryWrites:
    """A download run writes the registry in a few transactions, not one per accession."""

    @staticmethod
    def _args(tmp_path):
        return argparse.Namespace(
            accessions_file=str(tmp_path / "acc.txt"),
            fastq_folder=str(tmp_path / "fastq"),
            max_downloads=None,
            num_threads=4,
            max_workers=4,
            dry_run=False,
            force=False,
            max_retries=1,
            temp_folder=None,
            blacklist=None,
            report_file=None,
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
        )

    def test_run_outcomes_write_the_registry_once(self, tmp_path):
        import prepare.data.registry as registry_mod

        args = self._args(tmp_path)
        fastq_dir = tmp_path / "fastq"
        already = [f"SRR{i}" for i in range(300)]
        for accession in already:
            (fastq_dir / accession).mkdir(parents=True)
            (fastq_dir / accession / f"{accession}.fastq.gz").write_bytes(b"x")
        stats = {
            "already_downloaded_accessions": already,
            "blacklisted_accessions": ["SRR900"],
            "skipped_accessions": ["SRR901"],
        }
        with patch.object(registry_mod, "_write_registry", wraps=registry_mod._write_registry) as spy:
            DownloadSraCommand()._record_run_outcomes(args, stats, fastq_dir, store=None)
        assert spy.call_count == 1
        registry = load_registry(args.registry)
        assert rb.download_block(registry, "SRR0").state == "downloaded"
        assert rb.download_block(registry, "SRR299").files[0].bytes == 1
        assert rb.download_block(registry, "SRR900").message == "blacklisted"
        assert rb.download_block(registry, "SRR901").message == "--max-downloads"

    def test_store_sidecars_are_read_outside_the_registry_lock(self, tmp_path):
        """Sidecar verdicts for store-linked accessions are read before the one registry write takes its lock."""
        from prepare.store.layout import init_store

        args = self._args(tmp_path)
        fastq_dir = tmp_path / "fastq"
        store = init_store(tmp_path / "store")
        lock = tmp_path / "prepare_registry.json.lock"
        reads_under_lock = []

        def sidecar(_store, accession):
            reads_under_lock.append(lock.exists())
            return {"verdict": "complete", "ratio": 1.0}

        accessions = [f"SRR{i}" for i in range(5)]
        with (
            patch("prepare.cli.commands.sra.is_store_link", return_value=True),
            patch.object(DownloadSraCommand, "_sidecar_completeness", staticmethod(sidecar)),
            patch("prepare.cli.commands.sra.record_usage_many"),
        ):
            DownloadSraCommand()._record_run_outcomes(
                args, {"already_downloaded_accessions": accessions}, fastq_dir, store=store
            )
        assert reads_under_lock == [False] * 5
        registry = load_registry(args.registry)
        assert rb.download_block(registry, "SRR4").complete.verdict == "complete"
        assert rb.download_block(registry, "SRR4").source == "store"

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_download_results_are_batched(self, mock_download, _which, tmp_path):
        import prepare.data.registry as registry_mod

        def fake_download_sra(**kwargs):
            for index in range(20):
                kwargs["on_result"](f"SRR{index}", False, "Download failed: t")
            return {"total": 20, "successful": 0, "failed": 0, "failed_accessions": []}

        mock_download.side_effect = fake_download_sra
        args = self._args(tmp_path)
        with patch.object(registry_mod, "_write_registry", wraps=registry_mod._write_registry) as spy:
            assert DownloadSraCommand().execute(args) == 0
        assert spy.call_count == 1
        registry = load_registry(args.registry)
        assert all(rb.download_block(registry, f"SRR{i}").state == "failed" for i in range(20))

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_result_recorder_flushes_on_interrupt(self, mock_download, _which, tmp_path):
        """A result recorded before Ctrl-C is in the registry file after the run returns 130."""
        (tmp_path / "fastq" / "SRR1").mkdir(parents=True)
        (tmp_path / "fastq" / "SRR1" / "SRR1.fastq.gz").write_bytes(b"x")

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", True, "Downloaded 1 files")
            raise KeyboardInterrupt

        mock_download.side_effect = fake_download_sra
        args = self._args(tmp_path)
        assert DownloadSraCommand().execute(args) == 130
        assert rb.download_block(load_registry(args.registry), "SRR1").state == "downloaded"
        assert not (tmp_path / "prepare_registry.json.lock").exists()

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_already_exists_after_waiting_is_not_an_attempt(self, mock_download, _which, tmp_path):
        """A run that waited for another run's download and found the files counts no attempt."""
        from prepare.data.registry import record_download, registry_transaction

        args = self._args(tmp_path)
        (tmp_path / "fastq" / "SRR1").mkdir(parents=True)
        with registry_transaction(args.registry) as reg:
            record_download(reg, "SRR1", "downloaded", tmp_path / "fastq", "Downloaded 2 files, unverified")

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", True, "already exists")
            kwargs["on_result"]("SRR2", True, "already exists")
            return {"total": 2, "successful": 2, "failed": 0, "failed_accessions": []}

        mock_download.side_effect = fake_download_sra
        assert DownloadSraCommand().execute(args) == 0
        registry = load_registry(args.registry)
        first = rb.download_block(registry, "SRR1")
        assert (first.state, first.attempts, first.message) == ("downloaded", 1, "Downloaded 2 files, unverified")
        second = rb.download_block(registry, "SRR2")
        assert (second.state, second.attempts) == ("downloaded", 0)

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_retry_prefixed_already_exists_is_not_an_attempt(self, mock_download, _which, tmp_path):
        args = self._args(tmp_path)

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR2", True, "Retry 1: already exists")
            return {"total": 1, "successful": 1, "failed": 0, "failed_accessions": []}

        mock_download.side_effect = fake_download_sra
        assert DownloadSraCommand().execute(args) == 0
        second = rb.download_block(load_registry(args.registry), "SRR2")
        assert (second.state, second.attempts) == ("downloaded", 0)


class TestDownloadSraTerminationAndFinalFlush:
    """SIGTERM/SIGHUP take the interrupt path, and a failed final registry flush is retried once."""

    _args = staticmethod(TestDownloadSraRegistryWrites._args)

    @staticmethod
    def _failing_transactions(monkeypatch, failures):
        """Make the first ``failures`` registry-batch transactions raise DataAccessError."""
        import prepare.data.registry_batch as batch_mod
        from prepare.core.exceptions import DataAccessError

        real = batch_mod.registry_transaction
        calls = []

        def flaky(path):
            calls.append(path)
            if len(calls) <= failures:
                raise DataAccessError(f"Registry is locked by another process: {path}.lock")
            return real(path)

        monkeypatch.setattr(batch_mod, "registry_transaction", flaky)
        monkeypatch.setattr("prepare.cli.commands.sra.FINAL_FLUSH_RETRY_SECONDS", 0)
        return calls

    @staticmethod
    def _three_failed_results(**kwargs):
        for index in range(3):
            kwargs["on_result"](f"SRR{index}", False, "Download failed: t")
        return {"total": 3, "successful": 0, "failed": 0, "failed_accessions": [], "skipped_accessions": ["SRR9"]}

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_handlers_are_installed_for_the_batch_and_restored(self, mock_download, _which, tmp_path):
        import signal

        calls = []
        real_signal = signal.signal

        def spy(signum, handler):
            calls.append((signum, handler))
            return real_signal(signum, handler)

        mock_download.side_effect = self._three_failed_results
        handled = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
        before = {s: signal.getsignal(s) for s in handled}
        with patch.object(signal, "signal", side_effect=spy):
            assert DownloadSraCommand().execute(self._args(tmp_path)) == 0
        assert [c[0] for c in calls] == list(handled) * 2
        installed = calls[0][1]
        assert all(c[1] is installed for c in calls[:3]) and callable(installed)
        assert all(c[1] is before[c[0]] for c in calls[3:])
        assert {s: signal.getsignal(s) for s in before} == before

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_sigterm_during_the_run_flushes_and_returns_130(self, mock_download, _which, tmp_path, caplog):
        import signal

        def fake_download_sra(**kwargs):
            kwargs["on_result"]("SRR1", False, "Download failed: t")
            handler = signal.getsignal(signal.SIGTERM)
            try:
                handler(signal.SIGTERM, None)
            except KeyboardInterrupt:
                # A repeated signal while the interrupt path flushes is logged, not raised again.
                with caplog.at_level(logging.WARNING):
                    handler(signal.SIGTERM, None)
                assert any("SIGTERM again" in r.getMessage() for r in caplog.records)
                raise
            raise AssertionError("the first SIGTERM did not raise KeyboardInterrupt")

        before = signal.getsignal(signal.SIGTERM)
        mock_download.side_effect = fake_download_sra
        args = self._args(tmp_path)
        assert DownloadSraCommand().execute(args) == 130
        assert rb.download_block(load_registry(args.registry), "SRR1").state == "failed"
        assert signal.getsignal(signal.SIGTERM) is before

    def test_handlers_are_not_installed_off_the_main_thread(self):
        import signal
        import threading

        import prepare.cli.commands.sra as sra_mod

        errors = []

        def run():
            try:
                with patch.object(signal, "signal") as spy:
                    with sra_mod._termination_raises_interrupt():
                        pass
                    assert spy.call_count == 0
            except BaseException as e:  # noqa: B902 - reported to the main thread below
                errors.append(e)

        worker = threading.Thread(target=run)
        worker.start()
        worker.join()
        assert errors == []

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_final_flush_failing_once_is_retried(self, mock_download, _which, tmp_path, monkeypatch):
        calls = self._failing_transactions(monkeypatch, failures=1)
        mock_download.side_effect = self._three_failed_results
        args = self._args(tmp_path)
        assert DownloadSraCommand().execute(args) == 0
        registry = load_registry(args.registry)
        assert all(rb.download_block(registry, f"SRR{i}").state == "failed" for i in range(3))
        # The run outcomes (the skipped accession) were recorded after the retried flush.
        assert rb.download_block(registry, "SRR9").message == "--max-downloads"
        assert len(calls) == 3

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_final_flush_failing_twice_logs_the_queued_accessions(
        self, mock_download, _which, tmp_path, monkeypatch, caplog
    ):
        self._failing_transactions(monkeypatch, failures=2)
        mock_download.side_effect = self._three_failed_results
        args = self._args(tmp_path)
        with caplog.at_level(logging.ERROR):
            assert DownloadSraCommand().execute(args) == 1
        errors = [r.getMessage() for r in caplog.records if r.levelno == logging.ERROR]
        assert any("3 queued download outcome(s)" in m and "SRR0, SRR1, SRR2" in m for m in errors)
        assert not Path(args.registry).exists()


class TestDownloadSraMaxDownloadsArgument:
    """--max-downloads takes a positive integer; 0 is refused rather than read as no limit."""

    @staticmethod
    def _parser():
        parser = argparse.ArgumentParser()
        DownloadSraCommand().configure_parser(parser)
        return parser

    @pytest.mark.parametrize("value", ["0", "-3"])
    def test_zero_or_negative_is_rejected(self, value, capsys):
        with pytest.raises(SystemExit):
            self._parser().parse_args(["--accessions-file", "a.txt", "--max-downloads", value])
        assert "--max-downloads must be a positive integer" in capsys.readouterr().err

    def test_positive_value_is_accepted(self):
        assert self._parser().parse_args(["--accessions-file", "a.txt", "--max-downloads", "5"]).max_downloads == 5
        assert self._parser().parse_args(["--accessions-file", "a.txt"]).max_downloads is None


class TestFoldersAnchoredOnTheProject:
    """Run from proj/sub with the registry in proj/, relative folders resolve to proj/ (audit H4 follow-up)."""

    @staticmethod
    def _parse(argv):
        from prepare.cli.main import create_parser

        return create_parser().parse_args(argv)

    @pytest.fixture
    def proj(self, tmp_path, monkeypatch):
        proj = tmp_path / "proj"
        (proj / "sub").mkdir(parents=True)
        save_registry(load_registry(proj / "prepare_registry.json"))
        monkeypatch.chdir(proj / "sub")
        return proj.resolve()

    @patch("prepare.cli.commands.sra.require_tools")
    @patch("prepare.cli.commands.sra.download_sra")
    def test_download_sra_dry_run_uses_the_project_fastq_folder(self, mock_download, _tools, proj, caplog):
        mock_download.return_value = {"total": 1, "to_download": 1, "already_downloaded": 0, "failed": 0}
        (proj / "sub" / "accessions.txt").write_text("SRR1\n")
        args = self._parse(["download_sra", "--accessions-file", "accessions.txt", "--dry-run"])
        with caplog.at_level(logging.INFO):
            assert DownloadSraCommand().execute(args) == 0
        assert mock_download.call_args.kwargs["fastq_folder"] == str(proj / "fastq")
        assert f"Output folder would be: {proj / 'fastq'}" in caplog.text

    @patch("prepare.cli.commands.metadata.download_metadata")
    def test_download_metadata_dry_run_targets_the_project_folders(self, mock_download, proj):
        mock_download.return_value = {}
        args = self._parse(["download_metadata", "--email", "a@b.c", "--dry-run"])
        assert DownloadMetadataCommand().execute(args) == 0
        assert mock_download.call_args.kwargs["metadata_folder"] == str(proj / "metadata")
        assert mock_download.call_args.kwargs["matches_folder"] == str(proj / "matches")

    @patch("prepare.cli.commands.branchwater.process_branchwater_files")
    def test_use_branchwater_writes_the_project_matches_folder(self, mock_process, proj):
        args = self._parse(["use_branchwater", "--branchwater-folder", "bw"])
        assert UseBranchwaterCommand().execute(args) == 0
        assert mock_process.call_args.args == ("bw", str(proj / "matches"))

    def test_extract_branchwater_metadata_writes_the_project_metadata_folder(self, proj):
        (proj / "sub" / "bw").mkdir()
        args = self._parse(["extract_branchwater_metadata", "--branchwater-folder", "bw"])
        assert ExtractBranchwaterMetadataCommand().execute(args) == 0
        assert not (proj / "sub" / "metadata").exists()
        assert (proj / "metadata").is_dir()

    @patch("prepare.cli.commands.containment.parse_containment_data")
    def test_parse_containment_reads_the_project_matches_folder(self, mock_parse, proj):
        table = pd.DataFrame({"GCF_1": [0.9]}, index=pd.Index(["SRR1"]))
        mock_parse.return_value = argparse.Namespace(table=table, genome_to_samples={})
        args = self._parse(["parse_containment"])
        assert ParseContainmentCommand().execute(args) == 0
        assert mock_parse.call_args.args[0] == str(proj / "matches")
        # Recorded relative to the project root, so the registry stays movable with its project.
        csv = load_registry(proj / "prepare_registry.json").datasets["SRR1"]["screening"]["genomes"]["GCF_1"]["csv"]
        assert csv == "matches/GCF_1.csv"

    @patch("prepare.cli.commands.branchwater_search.search_index", return_value=[("SRR1", 0.9)])
    @patch("prepare.cli.commands.branchwater_search.load_signature", return_value=object())
    def test_branchwater_search_writes_the_project_branchwater_folder(self, _sig, _search, proj):
        from prepare.cli.commands.branchwater_search import BranchwaterSearchCommand

        args = self._parse(["branchwater_search", "--signature", "GCF_1.sig", "--no-cache"])
        assert BranchwaterSearchCommand().execute(args) == 0
        assert (proj / "branchwater" / "GCF_1.csv").exists()
        assert not (proj / "sub" / "branchwater").exists()
        csv = load_registry(proj / "prepare_registry.json").datasets["SRR1"]["screening"]["genomes"]["GCF_1"]["csv"]
        assert csv == "branchwater/GCF_1.csv"

    @patch("prepare.cli.commands.test_data.download_test_genome")
    def test_download_test_genome_writes_the_project_genomes_folder(self, mock_download, proj):
        args = self._parse(["download_test_genome"])
        assert DownloadTestGenomeCommand().execute(args) == 0
        mock_download.assert_called_once_with(str(proj / "genomes"))

    def test_sra_validate_reads_the_project_fastq_folder(self, proj, caplog):
        from prepare.cli.commands.sra_enhanced import SRAValidateCommand

        (proj / "fastq" / "SRR1").mkdir(parents=True)
        (proj / "fastq" / "SRR1" / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        args = self._parse(["sra_validate"])
        assert SRAValidateCommand().execute(args) == 0
        assert "does not exist" not in caplog.text

    def test_sra_profile_reads_the_project_fastq_folder(self, proj, caplog):
        from prepare.cli.commands.sra_profile import SRAProfileCommand

        (proj / "fastq" / "SRR1").mkdir(parents=True)
        (proj / "fastq" / "SRR1" / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n" * 4)
        args = self._parse(["sra_profile", "--output-dir", str(proj / "profiles")])
        assert SRAProfileCommand().execute(args) == 0
        assert "does not exist" not in caplog.text
