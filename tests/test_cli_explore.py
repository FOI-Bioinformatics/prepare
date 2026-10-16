"""Tests for prepare/cli/commands/explore.py.

Covers EnrichTaxonomyCommand, ExploreContainmentCommand, and
FindByTaxonomyCommand with mocked backends.

Run: pytest tests/test_cli_explore.py -v
"""

import argparse
import importlib.util
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from prepare.cli.commands.explore import (
    EnrichTaxonomyCommand,
    ExploreContainmentCommand,
    FindByTaxonomyCommand,
)
from prepare.core.models import TaxonomyInfo

# ============================================================================
# Helpers
# ============================================================================


requires_interactive = pytest.mark.skipif(
    not all(importlib.util.find_spec(m) for m in ("plotly", "jinja2")), reason="needs prepare[interactive]"
)


def _write_containment_tsv(path: Path):
    """Write a small containment TSV file for testing."""
    df = pd.DataFrame(
        {"genA": [0.5, 0.1], "genB": [0.8, 0.3]},
        index=["SRR001", "SRR002"],
    )
    df.to_csv(path, sep="\t")


def _write_taxonomy_tsv(path: Path):
    """Write a taxonomy map TSV."""
    with open(path, "w") as f:
        f.write("genome_id\tspecies\tgenus\tfamily\torder\tclass_name\tphylum\torganism\ttax_id\n")
        f.write("genA\tSp1\tGen1\tFam1\t\t\t\t\t\n")
        f.write("genB\tSp2\tGen2\tFam2\t\t\t\t\t\n")


# ============================================================================
# EnrichTaxonomyCommand
# ============================================================================


class TestEnrichTaxonomyCommand:

    def test_properties(self):
        cmd = EnrichTaxonomyCommand()
        assert cmd.name == "enrich_taxonomy"
        assert "taxonomy" in cmd.help.lower()

    def test_configure_parser(self):
        cmd = EnrichTaxonomyCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        args = parser.parse_args(["--parsed-containment", "data.tsv"])
        assert args.parsed_containment == "data.tsv"
        assert args.output == "taxonomy_map.tsv"
        assert args.cache == "taxonomy_cache.tsv"

    @patch("prepare.cli.commands.explore.enrich_genomes_with_taxonomy")
    @patch("prepare.cli.commands.explore.save_taxonomy_cache")
    def test_execute_success(self, mock_save, mock_enrich, tmp_path):
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)
        output = tmp_path / "out.tsv"

        mock_enrich.return_value = {
            "genA": TaxonomyInfo(genome_id="genA", family="Fam1"),
            "genB": TaxonomyInfo(genome_id="genB", family="Fam2"),
        }

        cmd = EnrichTaxonomyCommand()
        args = argparse.Namespace(
            parsed_containment=str(containment),
            output=str(output),
            cache=str(tmp_path / "cache.tsv"),
        )
        result = cmd.execute(args)
        assert result == 0
        mock_enrich.assert_called_once()
        mock_save.assert_called_once()

    def test_execute_missing_file(self, tmp_path):
        cmd = EnrichTaxonomyCommand()
        args = argparse.Namespace(
            parsed_containment=str(tmp_path / "missing.tsv"),
            output="out.tsv",
            cache="cache.tsv",
        )
        result = cmd.execute(args)
        assert result == 1

    @patch("prepare.cli.commands.explore.enrich_genomes_with_taxonomy")
    @patch("prepare.cli.commands.explore.save_taxonomy_cache")
    @patch("prepare.cli.commands.explore.pd.read_csv")
    def test_execute_reads_only_the_header(self, mock_read_csv, mock_save, mock_enrich, tmp_path):
        """Only the containment file's column names are needed, so the data rows are not read."""
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)
        mock_read_csv.return_value = pd.DataFrame(columns=["genA", "genB"])
        mock_enrich.return_value = {}

        cmd = EnrichTaxonomyCommand()
        args = argparse.Namespace(
            parsed_containment=str(containment),
            output=str(tmp_path / "out.tsv"),
            cache=str(tmp_path / "cache.tsv"),
        )
        result = cmd.execute(args)

        assert result == 0
        mock_read_csv.assert_called_once()
        _, kwargs = mock_read_csv.call_args
        assert kwargs.get("nrows") == 0


# ============================================================================
# ExploreContainmentCommand
# ============================================================================


class TestExploreContainmentCommand:

    def test_properties(self):
        cmd = ExploreContainmentCommand()
        assert cmd.name == "explore_containment"
        assert "explorer" in cmd.help.lower()

    def test_configure_parser(self):
        cmd = ExploreContainmentCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        args = parser.parse_args(
            [
                "--parsed-containment",
                "data.tsv",
                "--taxonomy-map",
                "tax.tsv",
                "--min-containment",
                "0.1",
            ]
        )
        assert args.parsed_containment == "data.tsv"
        assert args.taxonomy_map == "tax.tsv"
        assert args.min_containment == 0.1

    @requires_interactive
    @patch("prepare.cli.commands.explore.load_taxonomy_cache")
    @patch("prepare.visualization.explorer.generate_containment_explorer")
    def test_execute_with_taxonomy_map(self, mock_gen, mock_load, tmp_path):
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)
        tax_map = tmp_path / "tax.tsv"
        _write_taxonomy_tsv(tax_map)

        mock_load.return_value = {
            "genA": TaxonomyInfo(genome_id="genA", family="Fam1"),
        }
        mock_gen.return_value = tmp_path / "out.html"

        cmd = ExploreContainmentCommand()
        args = argparse.Namespace(
            parsed_containment=str(containment),
            taxonomy_map=str(tax_map),
            metadata=None,
            output=str(tmp_path / "out.html"),
            min_containment=0.0,
            cache=str(tmp_path / "cache.tsv"),
        )
        result = cmd.execute(args)
        assert result == 0
        mock_gen.assert_called_once()

    @requires_interactive
    @patch("prepare.cli.commands.explore.enrich_genomes_with_taxonomy")
    @patch("prepare.visualization.explorer.generate_containment_explorer")
    def test_execute_without_taxonomy_map(self, mock_gen, mock_enrich, tmp_path):
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)

        mock_enrich.return_value = {}
        mock_gen.return_value = tmp_path / "out.html"

        cmd = ExploreContainmentCommand()
        args = argparse.Namespace(
            parsed_containment=str(containment),
            taxonomy_map=None,
            metadata=None,
            output=str(tmp_path / "out.html"),
            min_containment=0.0,
            cache=str(tmp_path / "cache.tsv"),
        )
        result = cmd.execute(args)
        assert result == 0
        mock_enrich.assert_called_once()

    def test_execute_missing_containment(self, tmp_path):
        cmd = ExploreContainmentCommand()
        args = argparse.Namespace(
            parsed_containment=str(tmp_path / "missing.tsv"),
            taxonomy_map=None,
            metadata=None,
            output="out.html",
            min_containment=0.0,
            cache="cache.tsv",
        )
        result = cmd.execute(args)
        assert result == 1

    @requires_interactive
    @patch("prepare.cli.commands.explore.enrich_genomes_with_taxonomy")
    @patch("prepare.visualization.explorer.generate_containment_explorer")
    def test_execute_logs_a_next_hint_when_not_opened(self, mock_gen, mock_enrich, tmp_path, caplog):
        """Without --open, tell the user to open the HTML file themselves."""
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)

        mock_enrich.return_value = {}
        output = tmp_path / "out.html"
        mock_gen.return_value = output

        cmd = ExploreContainmentCommand()
        args = argparse.Namespace(
            parsed_containment=str(containment),
            taxonomy_map=None,
            metadata=None,
            output=str(output),
            min_containment=0.0,
            cache=str(tmp_path / "cache.tsv"),
            open_browser=False,
        )
        with caplog.at_level("INFO"):
            result = cmd.execute(args)
        assert result == 0
        assert "Next:" in caplog.text
        assert str(output) in caplog.text

    @requires_interactive
    @patch("prepare.utils.browser.open_in_browser")
    @patch("prepare.cli.commands.explore.enrich_genomes_with_taxonomy")
    @patch("prepare.visualization.explorer.generate_containment_explorer")
    def test_execute_no_next_hint_when_opened(self, mock_gen, mock_enrich, mock_open, tmp_path, caplog):
        """--open already opened it, so there is nothing left to hint at."""
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)

        mock_enrich.return_value = {}
        output = tmp_path / "out.html"
        mock_gen.return_value = output
        mock_open.return_value = True

        cmd = ExploreContainmentCommand()
        args = argparse.Namespace(
            parsed_containment=str(containment),
            taxonomy_map=None,
            metadata=None,
            output=str(output),
            min_containment=0.0,
            cache=str(tmp_path / "cache.tsv"),
            open_browser=True,
        )
        with caplog.at_level("INFO"):
            result = cmd.execute(args)
        assert result == 0
        assert "Next:" not in caplog.text

    def test_execute_missing_taxonomy_map(self, tmp_path):
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)

        cmd = ExploreContainmentCommand()
        args = argparse.Namespace(
            parsed_containment=str(containment),
            taxonomy_map=str(tmp_path / "no_such.tsv"),
            metadata=None,
            output="out.html",
            min_containment=0.0,
            cache="cache.tsv",
        )
        result = cmd.execute(args)
        assert result == 1


# ============================================================================
# FindByTaxonomyCommand
# ============================================================================


class TestFindByTaxonomyCommand:

    def test_properties(self):
        cmd = FindByTaxonomyCommand()
        assert cmd.name == "find_by_taxonomy"
        assert "filter" in cmd.help.lower()

    def test_configure_parser(self):
        cmd = FindByTaxonomyCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        args = parser.parse_args(
            [
                "--parsed-containment",
                "data.tsv",
                "--taxonomy-map",
                "tax.tsv",
                "--family",
                "Bacillaceae",
                "--min-containment",
                "0.1",
            ]
        )
        assert args.family == "Bacillaceae"
        assert args.min_containment == 0.1
        assert args.output_format == "table"

    def test_configure_parser_summary_format(self):
        cmd = FindByTaxonomyCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        args = parser.parse_args(
            [
                "--parsed-containment",
                "d.tsv",
                "--taxonomy-map",
                "t.tsv",
                "--genus",
                "Bacillus",
                "--format",
                "summary",
            ]
        )
        assert args.output_format == "summary"
        assert args.genus == "Bacillus"

    @patch("prepare.cli.commands.explore.load_taxonomy_cache")
    @patch("prepare.cli.commands.explore.annotate_containment_with_taxonomy")
    @patch("prepare.cli.commands.explore.filter_by_taxonomy")
    def test_execute_family_filter(self, mock_filter, mock_annotate, mock_load, tmp_path):
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)
        tax_map = tmp_path / "tax.tsv"
        _write_taxonomy_tsv(tax_map)

        mock_load.return_value = {
            "genA": TaxonomyInfo(genome_id="genA", family="Fam1"),
        }
        mock_annotate.return_value = pd.DataFrame(
            {
                "sample": ["S1"],
                "genome": ["genA"],
                "containment": [0.5],
                "family": ["Fam1"],
                "genus": ["G1"],
                "species": ["Sp1"],
            }
        )
        mock_filter.return_value = pd.DataFrame(
            {
                "sample": ["S1"],
                "genome": ["genA"],
                "containment": [0.5],
                "family": ["Fam1"],
                "genus": ["G1"],
                "species": ["Sp1"],
            }
        )

        cmd = FindByTaxonomyCommand()
        output_file = tmp_path / "results.tsv"
        args = argparse.Namespace(
            parsed_containment=str(containment),
            taxonomy_map=str(tax_map),
            family="Fam1",
            genus=None,
            species=None,
            min_containment=0.0,
            output=str(output_file),
            output_format="table",
        )
        result = cmd.execute(args)
        assert result == 0
        mock_filter.assert_called_once_with(
            mock_annotate.return_value,
            family="Fam1",
            genus=None,
            species=None,
            min_containment=0.0,
        )

    def test_execute_no_filter_returns_error(self, tmp_path):
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)
        tax_map = tmp_path / "tax.tsv"
        _write_taxonomy_tsv(tax_map)

        cmd = FindByTaxonomyCommand()
        args = argparse.Namespace(
            parsed_containment=str(containment),
            taxonomy_map=str(tax_map),
            family=None,
            genus=None,
            species=None,
            min_containment=0.0,
            output=None,
            output_format="table",
        )
        result = cmd.execute(args)
        assert result == 2

    @patch("prepare.cli.commands.explore.load_taxonomy_cache")
    @patch("prepare.cli.commands.explore.annotate_containment_with_taxonomy")
    @patch("prepare.cli.commands.explore.filter_by_taxonomy")
    def test_execute_empty_results(self, mock_filter, mock_annotate, mock_load, tmp_path):
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)
        tax_map = tmp_path / "tax.tsv"
        _write_taxonomy_tsv(tax_map)

        mock_load.return_value = {}
        mock_annotate.return_value = pd.DataFrame()
        mock_filter.return_value = pd.DataFrame()

        cmd = FindByTaxonomyCommand()
        args = argparse.Namespace(
            parsed_containment=str(containment),
            taxonomy_map=str(tax_map),
            family="NoMatch",
            genus=None,
            species=None,
            min_containment=0.0,
            output=None,
            output_format="table",
        )
        result = cmd.execute(args)
        assert result == 0  # empty is not an error

    def test_execute_missing_containment(self, tmp_path):
        tax_map = tmp_path / "tax.tsv"
        _write_taxonomy_tsv(tax_map)

        cmd = FindByTaxonomyCommand()
        args = argparse.Namespace(
            parsed_containment=str(tmp_path / "missing.tsv"),
            taxonomy_map=str(tax_map),
            family="Fam1",
            genus=None,
            species=None,
            min_containment=0.0,
            output=None,
            output_format="table",
        )
        result = cmd.execute(args)
        assert result == 1

    def test_execute_missing_taxonomy_map(self, tmp_path):
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)

        cmd = FindByTaxonomyCommand()
        args = argparse.Namespace(
            parsed_containment=str(containment),
            taxonomy_map=str(tmp_path / "missing.tsv"),
            family="Fam1",
            genus=None,
            species=None,
            min_containment=0.0,
            output=None,
            output_format="table",
        )
        result = cmd.execute(args)
        assert result == 1

    def test_table_output_has_no_index_column(self, tmp_path):
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)
        tax_map = tmp_path / "tax.tsv"
        _write_taxonomy_tsv(tax_map)
        annotated = pd.DataFrame(
            {"sample": ["S1"], "genome": ["genA"], "containment": [0.5], "family": ["Fam1"]},
            index=pd.Index(["S1"], name="Run_ID"),
        )
        output_file = tmp_path / "results.tsv"
        args = argparse.Namespace(
            parsed_containment=str(containment),
            taxonomy_map=str(tax_map),
            family="Fam1",
            genus=None,
            species=None,
            min_containment=0.0,
            output=str(output_file),
            output_format="table",
        )
        with (
            patch("prepare.cli.commands.explore.load_taxonomy_cache"),
            patch("prepare.cli.commands.explore.annotate_containment_with_taxonomy", return_value=annotated),
            patch("prepare.cli.commands.explore.filter_by_taxonomy", return_value=annotated),
        ):
            assert FindByTaxonomyCommand().execute(args) == 0
        assert output_file.read_text().splitlines()[0].split("\t")[0] == "sample"

    def test_table_output_to_stdout_has_no_index_column(self, tmp_path, capsys):
        containment = tmp_path / "containment.tsv"
        _write_containment_tsv(containment)
        tax_map = tmp_path / "tax.tsv"
        _write_taxonomy_tsv(tax_map)
        annotated = pd.DataFrame(
            {"sample": ["S1"], "genome": ["genA"], "containment": [0.5], "family": ["Fam1"]},
            index=pd.Index(["S1"], name="Run_ID"),
        )
        args = argparse.Namespace(
            parsed_containment=str(containment),
            taxonomy_map=str(tax_map),
            family="Fam1",
            genus=None,
            species=None,
            min_containment=0.0,
            output=None,
            output_format="table",
        )
        with (
            patch("prepare.cli.commands.explore.load_taxonomy_cache"),
            patch("prepare.cli.commands.explore.annotate_containment_with_taxonomy", return_value=annotated),
            patch("prepare.cli.commands.explore.filter_by_taxonomy", return_value=annotated),
        ):
            assert FindByTaxonomyCommand().execute(args) == 0
        assert capsys.readouterr().out.splitlines()[0].split("\t")[0] == "sample"
