"""Tests for genome CLI commands (search, download, prepare)."""

import argparse
import csv
import json
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from prepare.cli.commands.genome import (
    GenomeSearchCommand,
    GenomeDownloadCommand,
    GenomePrepareCommand,
)
from prepare.core.exceptions import PrepareError


@pytest.fixture(autouse=True)
def _datasets_present():
    """``download_genomes`` is mocked in these tests, so the ``datasets`` pre-flight check is too;
    tests/test_tools.py covers a missing ``datasets`` with an empty ``PATH``."""
    with patch("prepare.cli.commands.genome.require_tools"):
        yield


class TestGenomeSearchCommand:
    """Tests for GenomeSearchCommand."""

    def test_command_properties(self):
        cmd = GenomeSearchCommand()
        assert cmd.name == "genome_search"
        assert "search" in cmd.help.lower()

    def test_configure_parser(self):
        cmd = GenomeSearchCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        args = parser.parse_args(["--species", "Salmonella enterica"])
        assert args.species == "Salmonella enterica"
        assert args.genus is None
        assert args.representative_only is True
        assert args.output_format == "list"

    def test_configure_parser_genus(self):
        cmd = GenomeSearchCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        args = parser.parse_args(["--genus", "Salmonella"])
        assert args.genus == "Salmonella"
        assert args.species is None

    def test_configure_parser_mutually_exclusive(self):
        cmd = GenomeSearchCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        with pytest.raises(SystemExit):
            parser.parse_args(["--species", "Foo", "--genus", "Bar"])

    def test_configure_parser_missing_required(self):
        cmd = GenomeSearchCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_configure_parser_all_flag(self):
        cmd = GenomeSearchCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        args = parser.parse_args(["--species", "Foo", "--all"])
        assert args.representative_only is False

    @patch("prepare.cli.commands.genome.get_accessions_for_species")
    def test_execute_species_search(self, mock_get):
        mock_get.return_value = ["GCF_000006945.2", "GCF_000007545.1"]
        cmd = GenomeSearchCommand()
        args = argparse.Namespace(
            species="Salmonella enterica",
            genus=None,
            output=None,
            representative_only=True,
            output_format="list",
        )
        result = cmd.execute(args)
        assert result == 0
        mock_get.assert_called_once_with("Salmonella enterica", representative_only=True)

    @patch("prepare.cli.commands.genome.get_accessions_for_genus")
    def test_execute_genus_search(self, mock_get):
        mock_get.return_value = ["GCF_000006945.2"]
        cmd = GenomeSearchCommand()
        args = argparse.Namespace(
            species=None,
            genus="Salmonella",
            output=None,
            representative_only=True,
            output_format="list",
        )
        result = cmd.execute(args)
        assert result == 0
        mock_get.assert_called_once_with("Salmonella", representative_only=True)

    @patch("prepare.cli.commands.genome.get_accessions_for_species")
    def test_execute_output_to_file(self, mock_get):
        mock_get.return_value = ["GCF_000006945.2"]
        cmd = GenomeSearchCommand()
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            outpath = f.name

        try:
            args = argparse.Namespace(
                species="Salmonella enterica",
                genus=None,
                output=outpath,
                representative_only=True,
                output_format="list",
            )
            result = cmd.execute(args)
            assert result == 0
            content = Path(outpath).read_text()
            assert "GCF_000006945.2" in content
        finally:
            Path(outpath).unlink(missing_ok=True)

    @patch("prepare.cli.commands.genome.get_accessions_for_species")
    def test_execute_tsv_format(self, mock_get):
        mock_get.return_value = ["GCF_000006945.2"]
        cmd = GenomeSearchCommand()
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            outpath = f.name

        try:
            args = argparse.Namespace(
                species="Salmonella enterica",
                genus=None,
                output=outpath,
                representative_only=True,
                output_format="tsv",
            )
            result = cmd.execute(args)
            assert result == 0
            content = Path(outpath).read_text()
            assert "GCF_000006945.2\tSalmonella enterica" in content
        finally:
            Path(outpath).unlink(missing_ok=True)

    @patch("prepare.cli.commands.genome.get_accessions_for_species")
    def test_execute_no_results(self, mock_get):
        mock_get.return_value = []
        cmd = GenomeSearchCommand()
        args = argparse.Namespace(
            species="Nonexistent species",
            genus=None,
            output=None,
            representative_only=True,
            output_format="list",
        )
        result = cmd.execute(args)
        assert result == 0

    @patch("prepare.cli.commands.genome.get_accessions_for_species")
    def test_execute_api_failure(self, mock_get):
        mock_get.side_effect = PrepareError("API error")
        cmd = GenomeSearchCommand()
        args = argparse.Namespace(
            species="Salmonella enterica",
            genus=None,
            output=None,
            representative_only=True,
            output_format="list",
        )
        result = cmd.execute(args)
        assert result == 1


class TestGenomeDownloadCommand:
    """Tests for GenomeDownloadCommand."""

    def test_command_properties(self):
        cmd = GenomeDownloadCommand()
        assert cmd.name == "genome_download"
        assert "download" in cmd.help.lower()

    def test_configure_parser(self):
        cmd = GenomeDownloadCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        args = parser.parse_args(["--accessions", "GCF_000006945.2"])
        assert args.accessions == ["GCF_000006945.2"]
        assert args.output_dir == "genomes/"
        assert args.representative_only is True
        assert args.registry is None

    def test_configure_parser_has_registry(self):
        cmd = GenomeDownloadCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        args = parser.parse_args(["--accessions", "GCF_000006945.2", "--registry", "r.json"])
        assert args.registry == "r.json"

    def test_configure_parser_multiple_accessions(self):
        cmd = GenomeDownloadCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        args = parser.parse_args(["--accessions", "GCF_000006945.2", "GCF_000007545.1"])
        assert len(args.accessions) == 2

    def test_execute_missing_all_sources(self):
        cmd = GenomeDownloadCommand()
        args = argparse.Namespace(
            accessions=None,
            accession_file=None,
            species=None,
            genus=None,
            output_dir="genomes/",
            representative_only=True,
            assembly_level=None,
            force=False,
            dry_run=False,
        )
        result = cmd.execute(args)
        assert result == 1

    @patch("prepare.cli.commands.genome.extract_and_organize")
    @patch("prepare.cli.commands.genome.download_genomes")
    def test_execute_with_accessions(self, mock_download, mock_extract):
        mock_download.return_value = Path("genomes/download.zip")
        mock_extract.return_value = {"GCF_000006945.2": Path("genomes/GCF_000006945.2.fna.gz")}
        cmd = GenomeDownloadCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(
                accessions=["GCF_000006945.2"],
                accession_file=None,
                species=None,
                genus=None,
                output_dir=tmpdir,
                representative_only=True,
                assembly_level=None,
                force=False,
                dry_run=False,
                registry=str(Path(tmpdir) / "prepare_registry.json"),
            )
            result = cmd.execute(args)
            assert result == 0
            mock_download.assert_called_once()
            mock_extract.assert_called_once()

    @patch("prepare.cli.commands.genome.extract_and_organize")
    @patch("prepare.cli.commands.genome.download_genomes")
    def test_execute_with_no_extracted_genomes_leaves_registry_untouched(self, mock_download, mock_extract):
        """When extract_and_organize returns nothing (e.g. the zip held no usable genome
        files), no registry transaction should even open: there is nothing to record, so no
        registry file should be created."""
        mock_download.return_value = Path("genomes/download.zip")
        mock_extract.return_value = {}
        cmd = GenomeDownloadCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            registry_file = Path(tmpdir) / "prepare_registry.json"
            args = argparse.Namespace(
                accessions=["GCF_000006945.2"],
                accession_file=None,
                species=None,
                genus=None,
                output_dir=tmpdir,
                representative_only=True,
                assembly_level=None,
                force=False,
                dry_run=False,
                registry=str(registry_file),
            )
            result = cmd.execute(args)
            assert result == 0
            assert not registry_file.exists()

    @patch("prepare.cli.commands.genome.extract_and_organize")
    @patch("prepare.cli.commands.genome.download_genomes")
    def test_execute_with_accession_file(self, mock_download, mock_extract):
        mock_download.return_value = Path("genomes/download.zip")
        mock_extract.return_value = {"GCF_000006945.2": Path("genomes/GCF_000006945.2.fna.gz")}
        cmd = GenomeDownloadCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            acc_file = Path(tmpdir) / "accessions.txt"
            acc_file.write_text("GCF_000006945.2\nGCF_000007545.1\n")
            args = argparse.Namespace(
                accessions=None,
                accession_file=str(acc_file),
                species=None,
                genus=None,
                output_dir=tmpdir,
                representative_only=True,
                assembly_level=None,
                force=False,
                dry_run=False,
                registry=str(Path(tmpdir) / "prepare_registry.json"),
            )
            result = cmd.execute(args)
            assert result == 0

    @patch("prepare.cli.commands.genome.extract_and_organize")
    @patch("prepare.cli.commands.genome.download_genomes")
    @patch("prepare.cli.commands.genome.get_accessions_for_species")
    def test_execute_with_species(self, mock_gtdb, mock_download, mock_extract):
        mock_gtdb.return_value = ["GCF_000006945.2"]
        mock_download.return_value = Path("genomes/download.zip")
        mock_extract.return_value = {"GCF_000006945.2": Path("genomes/GCF_000006945.2.fna.gz")}
        cmd = GenomeDownloadCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(
                accessions=None,
                accession_file=None,
                species="Salmonella enterica",
                genus=None,
                output_dir=tmpdir,
                representative_only=True,
                assembly_level=None,
                force=False,
                dry_run=False,
                registry=str(Path(tmpdir) / "prepare_registry.json"),
            )
            result = cmd.execute(args)
            assert result == 0
            mock_gtdb.assert_called_once_with("Salmonella enterica", representative_only=True)

    def test_execute_accession_file_not_found(self):
        cmd = GenomeDownloadCommand()
        args = argparse.Namespace(
            accessions=None,
            accession_file="/nonexistent/file.txt",
            species=None,
            genus=None,
            output_dir="genomes/",
            representative_only=True,
            assembly_level=None,
            force=False,
            dry_run=False,
        )
        result = cmd.execute(args)
        assert result == 1

    @patch("prepare.cli.commands.genome.download_genomes")
    def test_execute_download_failure(self, mock_download):
        mock_download.side_effect = PrepareError("Download failed")
        cmd = GenomeDownloadCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(
                accessions=["GCF_000006945.2"],
                accession_file=None,
                species=None,
                genus=None,
                output_dir=tmpdir,
                representative_only=True,
                assembly_level=None,
                force=False,
                dry_run=False,
            )
            result = cmd.execute(args)
            assert result == 1

    @patch("prepare.cli.commands.genome.extract_and_organize")
    @patch("prepare.cli.commands.genome.download_genomes")
    def test_execute_skips_present_genome(self, mock_download, mock_extract):
        cmd = GenomeDownloadCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            (Path(tmpdir) / "GCF_000006945.2.fna").write_text(">seq\nACGT\n")
            args = argparse.Namespace(
                accessions=["GCF_000006945.2"],
                accession_file=None,
                species=None,
                genus=None,
                output_dir=tmpdir,
                representative_only=True,
                assembly_level=None,
                force=False,
                dry_run=False,
            )
            result = cmd.execute(args)
            assert result == 0
            mock_download.assert_not_called()
            mock_extract.assert_not_called()

    @patch("prepare.cli.commands.genome.extract_and_organize")
    @patch("prepare.cli.commands.genome.download_genomes")
    def test_execute_force_redownloads_present_genome(self, mock_download, mock_extract):
        mock_download.return_value = Path("genomes/download.zip")
        mock_extract.return_value = {"GCF_000006945.2": Path("genomes/GCF_000006945.2.fna")}
        cmd = GenomeDownloadCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            (Path(tmpdir) / "GCF_000006945.2.fna").write_text(">seq\nACGT\n")
            args = argparse.Namespace(
                accessions=["GCF_000006945.2"],
                accession_file=None,
                species=None,
                genus=None,
                output_dir=tmpdir,
                representative_only=True,
                assembly_level=None,
                force=True,
                dry_run=False,
                registry=str(Path(tmpdir) / "prepare_registry.json"),
            )
            result = cmd.execute(args)
            assert result == 0
            mock_download.assert_called_once()

    @patch("prepare.cli.commands.genome.extract_and_organize")
    @patch("prepare.cli.commands.genome.download_genomes")
    def test_execute_dry_run_downloads_nothing(self, mock_download, mock_extract):
        cmd = GenomeDownloadCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            (Path(tmpdir) / "GCF_000006945.2.fna").write_text(">seq\nACGT\n")
            args = argparse.Namespace(
                accessions=["GCF_000006945.2", "GCF_000007545.1"],
                accession_file=None,
                species=None,
                genus=None,
                output_dir=tmpdir,
                representative_only=True,
                assembly_level=None,
                force=False,
                dry_run=True,
            )
            result = cmd.execute(args)
            assert result == 0
            mock_download.assert_not_called()

    @patch("prepare.cli.commands.genome.extract_and_organize")
    @patch("prepare.cli.commands.genome.download_genomes")
    def test_execute_records_genomes_in_registry(self, mock_download, mock_extract, tmp_path):
        """genome_download must record what it fetched, so later commands (e.g. status) find it."""
        mock_download.return_value = Path("genomes/download.zip")
        mock_extract.return_value = {
            "GCF_000006945.2": tmp_path / "GCF_000006945.2.fna",
            "GCF_000007545.1": tmp_path / "GCF_000007545.1.fna",
        }
        cmd = GenomeDownloadCommand()
        registry_file = tmp_path / "prepare_registry.json"
        args = argparse.Namespace(
            accessions=["GCF_000006945.2", "GCF_000007545.1"],
            accession_file=None,
            species=None,
            genus=None,
            output_dir=str(tmp_path),
            representative_only=True,
            assembly_level=None,
            force=False,
            dry_run=False,
            registry=str(registry_file),
        )

        result = cmd.execute(args)
        assert result == 0

        genomes = json.loads(registry_file.read_text())["genomes"]
        assert set(genomes) == {"GCF_000006945.2", "GCF_000007545.1"}
        assert genomes["GCF_000006945.2"]["fasta"] == "GCF_000006945.2.fna"
        assert genomes["GCF_000006945.2"]["date"]

    @patch("prepare.cli.commands.genome.extract_and_organize")
    @patch("prepare.cli.commands.genome.download_genomes")
    def test_force_redownload_keeps_an_earlier_manifest(self, mock_download, mock_extract, tmp_path, monkeypatch):
        """genome_download writes no manifest of its own; a --force redownload of a genome that
        genome_prepare recorded with a manifest keeps that manifest rather than blanking it."""
        from prepare.data.registry import load_registry, record_genome, save_registry

        monkeypatch.chdir(tmp_path)
        registry_file = tmp_path / "prepare_registry.json"
        fasta = tmp_path / "GCF_000006945.2.fna"
        fasta.write_text(">c\nACGT\n")
        seeded = load_registry(registry_file)
        record_genome(seeded, "GCF_000006945.2", fasta, tmp_path / "genome_manifest.csv")
        save_registry(seeded)

        mock_download.return_value = Path("genomes/download.zip")
        mock_extract.return_value = {"GCF_000006945.2": fasta}
        args = argparse.Namespace(
            accessions=["GCF_000006945.2"],
            accession_file=None,
            species=None,
            genus=None,
            output_dir=str(tmp_path),
            representative_only=True,
            assembly_level=None,
            force=True,
            dry_run=False,
            registry=str(registry_file),
        )
        assert GenomeDownloadCommand().execute(args) == 0

        genome = json.loads(registry_file.read_text())["genomes"]["GCF_000006945.2"]
        assert genome["manifest"] == "genome_manifest.csv"
        assert genome["fasta"] == "GCF_000006945.2.fna"


class TestGenomePrepareCommand:
    """Tests for GenomePrepareCommand."""

    def test_command_properties(self):
        cmd = GenomePrepareCommand()
        assert cmd.name == "genome_prepare"
        assert "manifest" in cmd.help.lower() or "prepare" in cmd.help.lower()

    def test_configure_parser(self):
        cmd = GenomePrepareCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        args = parser.parse_args(["--species", "Salmonella enterica"])
        assert args.species == "Salmonella enterica"
        assert args.manifest_file == "genome_manifest.csv"
        assert args.skip_download is False

    @patch("prepare.cli.commands.genome.extract_and_organize")
    @patch("prepare.cli.commands.genome.download_genomes")
    @patch("prepare.cli.commands.genome.get_accessions_for_species")
    def test_execute_full_pipeline(self, mock_gtdb, mock_download, mock_extract):
        mock_gtdb.return_value = ["GCF_000006945.2"]
        mock_download.return_value = Path("download.zip")
        mock_extract.return_value = {"GCF_000006945.2": Path("GCF_000006945.2.fna.gz")}
        cmd = GenomePrepareCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a fake genome file so manifest picks it up
            genome_file = Path(tmpdir) / "GCF_000006945.2.fna.gz"
            genome_file.touch()
            manifest = Path(tmpdir) / "manifest.csv"
            args = argparse.Namespace(
                species="Salmonella enterica",
                genus=None,
                accession_file=None,
                output_dir=tmpdir,
                representative_only=True,
                manifest_file=str(manifest),
                skip_download=False,
                registry=str(Path(tmpdir) / "prepare_registry.json"),
            )
            result = cmd.execute(args)
            assert result == 0
            assert manifest.exists()
            content = manifest.read_text()
            assert "name,genome_filename,protein_filename" in content
            assert "GCF_000006945.2" in content

    def test_execute_skip_download(self):
        cmd = GenomePrepareCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create fake genome files
            (Path(tmpdir) / "genome1.fna.gz").touch()
            (Path(tmpdir) / "genome2.fna.gz").touch()
            manifest = Path(tmpdir) / "manifest.csv"
            args = argparse.Namespace(
                species=None,
                genus=None,
                accession_file=None,
                output_dir=tmpdir,
                representative_only=True,
                manifest_file=str(manifest),
                skip_download=True,
                registry=str(Path(tmpdir) / "prepare_registry.json"),
            )
            result = cmd.execute(args)
            assert result == 0
            assert manifest.exists()
            with open(manifest) as f:
                reader = csv.reader(f)
                rows = list(reader)
            assert len(rows) == 3  # header + 2 entries

    def test_execute_missing_sources_no_skip(self):
        cmd = GenomePrepareCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(
                species=None,
                genus=None,
                accession_file=None,
                output_dir=tmpdir,
                representative_only=True,
                manifest_file="manifest.csv",
                skip_download=False,
                registry=str(Path(tmpdir) / "prepare_registry.json"),
            )
            result = cmd.execute(args)
            assert result == 1

    @patch("prepare.cli.commands.genome.download_genomes")
    @patch("prepare.cli.commands.genome.get_accessions_for_species")
    def test_execute_api_failure(self, mock_gtdb, mock_download):
        mock_gtdb.side_effect = PrepareError("API error")
        cmd = GenomePrepareCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(
                species="Salmonella enterica",
                genus=None,
                accession_file=None,
                output_dir=tmpdir,
                representative_only=True,
                manifest_file="manifest.csv",
                skip_download=False,
                registry=str(Path(tmpdir) / "prepare_registry.json"),
            )
            result = cmd.execute(args)
            assert result == 1

    def test_execute_skip_download_no_genomes(self):
        cmd = GenomePrepareCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest = Path(tmpdir) / "manifest.csv"
            args = argparse.Namespace(
                species=None,
                genus=None,
                accession_file=None,
                output_dir=tmpdir,
                representative_only=True,
                manifest_file=str(manifest),
                skip_download=True,
                registry=str(Path(tmpdir) / "prepare_registry.json"),
            )
            result = cmd.execute(args)
            assert result == 0

    def test_prepare_with_no_genomes_leaves_registry_untouched(self, tmp_path):
        """With no genome files to manifest, no registry transaction should even open: an
        empty run must not create (or touch) a registry file that never gets anything
        written to it."""
        registry_file = tmp_path / "prepare_registry.json"
        genomes_dir = tmp_path / "genomes"
        genomes_dir.mkdir()
        manifest = tmp_path / "manifest.csv"
        args = argparse.Namespace(
            species=None,
            genus=None,
            accession_file=None,
            output_dir=str(genomes_dir),
            representative_only=True,
            manifest_file=str(manifest),
            skip_download=True,
            registry=str(registry_file),
        )

        result = GenomePrepareCommand().execute(args)

        assert result == 0
        assert not registry_file.exists()

    def test_manifest_lists_plain_fna_files(self):
        from prepare.cli.commands.genome import GenomePrepareCommand

        cmd = GenomePrepareCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            (tmp_path / "GCF_000006945.2.fna").write_text(">c\nACGT\n")
            (tmp_path / "GCF_000006945.2.faa").write_text(">p\nM\n")
            (tmp_path / "other.fasta.gz").write_bytes(b"")
            manifest = tmp_path / "manifest.csv"
            n = cmd._create_manifest(tmp_path, str(manifest), str(tmp_path / "prepare_registry.json"))
            rows = manifest.read_text().splitlines()
            assert n == 2
            assert rows[0] == "name,genome_filename,protein_filename"
            assert rows[1].startswith("GCF_000006945.2,") and rows[1].endswith("GCF_000006945.2.faa")
            assert rows[2].startswith("other,")

    def test_manifest_rows_are_recorded_as_genomes(self, tmp_path):
        """Every manifest row becomes a registry genome, so later commands can find its FASTA."""
        cmd = GenomePrepareCommand()
        (tmp_path / "GCF_000006945.2.fna").write_text(">c\nACGT\n")
        (tmp_path / "wMel.fasta").write_text(">c\nACGT\n")
        manifest = tmp_path / "manifest.csv"
        registry_file = tmp_path / "prepare_registry.json"
        args = argparse.Namespace(
            species=None,
            genus=None,
            accession_file=None,
            output_dir=str(tmp_path),
            representative_only=True,
            manifest_file=str(manifest),
            skip_download=True,
            registry=str(registry_file),
        )

        assert cmd.execute(args) == 0
        genomes = json.loads(registry_file.read_text())["genomes"]
        assert set(genomes) == {"GCF_000006945.2", "wMel"}
        # Both files live under the project root (the registry's own folder), so the registry
        # records them relative to it, which keeps the project movable.
        assert genomes["GCF_000006945.2"]["fasta"] == "GCF_000006945.2.fna"
        assert genomes["GCF_000006945.2"]["manifest"] == "manifest.csv"
        assert genomes["wMel"]["date"]

    def test_configure_parser_has_registry(self):
        cmd = GenomePrepareCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        assert parser.parse_args(["--registry", "r.json"]).registry == "r.json"

    @patch("prepare.cli.commands.genome.extract_and_organize")
    @patch("prepare.cli.commands.genome.download_genomes")
    @patch("prepare.cli.commands.genome.get_accessions_for_genus")
    def test_execute_with_genus(self, mock_gtdb, mock_download, mock_extract):
        mock_gtdb.return_value = ["GCF_000006945.2", "GCF_000007545.1"]
        mock_download.return_value = Path("download.zip")
        mock_extract.return_value = {}
        cmd = GenomePrepareCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest = Path(tmpdir) / "manifest.csv"
            args = argparse.Namespace(
                species=None,
                genus="Salmonella",
                accession_file=None,
                output_dir=tmpdir,
                representative_only=True,
                manifest_file=str(manifest),
                skip_download=False,
                registry=str(Path(tmpdir) / "prepare_registry.json"),
            )
            result = cmd.execute(args)
            assert result == 0
            mock_gtdb.assert_called_once_with("Salmonella", representative_only=True)

    @patch("prepare.cli.commands.genome.extract_and_organize")
    @patch("prepare.cli.commands.genome.download_genomes")
    def test_execute_with_accession_file(self, mock_download, mock_extract):
        mock_download.return_value = Path("download.zip")
        mock_extract.return_value = {}
        cmd = GenomePrepareCommand()
        with tempfile.TemporaryDirectory() as tmpdir:
            acc_file = Path(tmpdir) / "accessions.txt"
            acc_file.write_text("GCF_000006945.2\n# comment\nGCF_000007545.1\n")
            manifest = Path(tmpdir) / "manifest.csv"
            args = argparse.Namespace(
                species=None,
                genus=None,
                accession_file=str(acc_file),
                output_dir=tmpdir,
                representative_only=True,
                manifest_file=str(manifest),
                skip_download=False,
                registry=str(Path(tmpdir) / "prepare_registry.json"),
            )
            result = cmd.execute(args)
            assert result == 0


def test_manifest_ignores_appledouble_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "GCF_1.fna").write_text(">c\nACGT\n")
    (genomes / "._GCF_1.fna").write_bytes(b"\x00\x05\x16\x07")
    count = GenomePrepareCommand()._create_manifest(genomes, "genome_manifest.csv", str(tmp_path / "r.json"))
    assert count == 1
    assert "._GCF_1" not in (tmp_path / "genome_manifest.csv").read_text()
