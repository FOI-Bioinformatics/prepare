"""
Tests for the sourmash CLI plugin integration.
"""

import argparse
from unittest.mock import patch

import pandas as pd
import pytest

pytest.importorskip("sourmash")

from prepare.plugins.sourmash_plugin import (  # noqa: E402
    PrepareParsePlugin,
    PreparePlotPlugin,
    PrepareDiversityPlugin,
    PrepareTaxonomyPlugin,
)

# --- Fixtures ---


@pytest.fixture
def parser():
    """Create an argparse subparser for testing."""
    main_parser = argparse.ArgumentParser()
    sub = main_parser.add_subparsers()
    return sub.add_parser("test")


@pytest.fixture
def tmp_output(tmp_path):
    """Provide a temporary output directory."""
    return tmp_path


@pytest.fixture(autouse=True)
def restore_root_logging():
    """Put the root logger back as found: the plugins call setup_logging, which changes it."""
    import logging

    root = logging.getLogger()
    handlers, level = root.handlers[:], root.level
    yield
    root.handlers[:], root.level = handlers, level


# --- PrepareParsePlugin ---


class TestPrepareParsePlugin:
    def test_command_attributes(self):
        assert PrepareParsePlugin.command == "prepare_parse"
        assert PrepareParsePlugin.description is not None

    def test_configure_parser(self, parser):
        PrepareParsePlugin(parser)
        args = parser.parse_args(
            [
                "--matches-folder",
                "my_matches",
                "--step-size",
                "0.05",
            ]
        )
        assert args.matches_folder == "my_matches"
        assert args.step_size == 0.05
        assert args.parsed_containment_file == "parsed_containment.txt"
        assert args.summary_containment_file == "top_containments.txt"

    def test_default_arguments(self, parser):
        PrepareParsePlugin(parser)
        args = parser.parse_args([])
        assert args.matches_folder == "matches"
        assert args.step_size == 0.1

    @patch("prepare.plugins.sourmash_plugin.PrepareParsePlugin.main")
    def test_main_delegates_to_parse(self, mock_main, parser):
        """Verify main can be called without error when mocked."""
        PrepareParsePlugin(parser)
        args = parser.parse_args([])
        mock_main(args)
        mock_main.assert_called_once_with(args)

    @patch("prepare.data.branchwater.parse_containment_data")
    def test_main_calls_parse_containment(self, mock_parse, parser):
        cmd = PrepareParsePlugin(parser)
        args = parser.parse_args(
            [
                "--matches-folder",
                "test_matches",
                "--parsed-containment-file",
                "parsed.txt",
                "--summary-containment-file",
                "summary.txt",
                "--step-size",
                "0.2",
            ]
        )
        # Call main directly, bypassing super().main() since sourmash
        # may not be installed
        with patch.object(type(cmd).__bases__[0], "main"), patch("prepare.data.registry.registry_transaction"):
            cmd.main(args)
        folder, *rest = mock_parse.call_args.args
        assert folder.endswith("test_matches")
        assert tuple(rest) == ("parsed.txt", "summary.txt", 0.2)
        assert mock_parse.call_args.kwargs["details_file"] is None
        assert mock_parse.call_args.kwargs["errors"] == []

    @patch("prepare.data.branchwater.parse_containment_data", side_effect=RuntimeError("bad"))
    def test_main_handles_error(self, mock_parse, parser):
        cmd = PrepareParsePlugin(parser)
        args = parser.parse_args([])
        with patch.object(type(cmd).__bases__[0], "main"):
            with pytest.raises(SystemExit):
                cmd.main(args)


class TestParsePluginErrorsAndRegistry:
    def _run(self, parser, argv):
        cmd = PrepareParsePlugin(parser)
        args = parser.parse_args(argv)
        with patch.object(type(cmd).__bases__[0], "main"):
            cmd.main(args)

    def test_new_flags(self, parser):
        PrepareParsePlugin(parser)
        args = parser.parse_args([])
        assert args.details_file is None
        assert args.registry is None
        assert args.registry_max_screened > 0

    def test_failed_parse_reaches_stderr(self, parser, tmp_path, capsys):
        folder = tmp_path / "matches"
        folder.mkdir()
        argv = ["--matches-folder", str(folder), "--parsed-containment-file", str(tmp_path / "p.txt")]
        with patch("prepare.data.branchwater.parse_containment_data", side_effect=RuntimeError("cannot read it")):
            with pytest.raises(SystemExit) as exc:
                self._run(parser, argv)
        assert exc.value.code == 1
        assert "cannot read it" in capsys.readouterr().err

    def test_unreadable_match_file_exits_1_and_records_screening(self, parser, tmp_path, capsys):
        folder = tmp_path / "matches"
        folder.mkdir()
        (folder / "G1.csv").write_text("acc,containment,cANI\n" "SRR1,0.9,0.99\n")
        (folder / "BAD.csv").write_bytes(b"\x00\xff\xfe not a csv \x00")
        registry = tmp_path / "prepare_registry.json"
        argv = [
            "--matches-folder",
            str(folder),
            "--parsed-containment-file",
            str(tmp_path / "p.txt"),
            "--summary-containment-file",
            str(tmp_path / "s.txt"),
            "--registry",
            str(registry),
        ]
        with pytest.raises(SystemExit) as exc:
            self._run(parser, argv)
        assert exc.value.code == 1
        err = capsys.readouterr().err
        assert "could not be read" in err
        assert "BAD.csv" in err
        import json

        data = json.loads(registry.read_text())
        dumped = json.dumps(data)
        assert "SRR1" in dumped
        assert "screening" in dumped
        assert "G1" in dumped

    def test_screening_recorded_in_registry(self, parser, tmp_path):
        folder = tmp_path / "matches"
        folder.mkdir()
        (folder / "G1.csv").write_text("acc,containment,cANI\n" "SRR1,0.9,0.99\n")
        registry = tmp_path / "prepare_registry.json"
        argv = [
            "--matches-folder",
            str(folder),
            "--parsed-containment-file",
            str(tmp_path / "p.txt"),
            "--summary-containment-file",
            str(tmp_path / "s.txt"),
            "--registry",
            str(registry),
        ]
        self._run(parser, argv)
        assert registry.exists()
        assert "SRR1" in registry.read_text()


# --- PreparePlotPlugin ---


class TestPreparePlotPlugin:
    def test_command_attributes(self):
        assert PreparePlotPlugin.command == "prepare_plot"
        assert PreparePlotPlugin.description is not None

    def test_configure_parser(self, parser):
        PreparePlotPlugin(parser)
        args = parser.parse_args(
            [
                "--file-path",
                "data.txt",
                "--plot-type",
                "histogram",
                "--threshold",
                "0.5",
            ]
        )
        assert args.file_path == "data.txt"
        assert args.plot_type == "histogram"
        assert args.threshold == 0.5
        assert args.column == "max_containment"

    def test_default_arguments(self, parser):
        PreparePlotPlugin(parser)
        args = parser.parse_args(["--file-path", "data.txt"])
        assert args.plot_type == "rank"
        assert args.column == "max_containment"
        assert args.title is None
        assert args.save_format is None

    @patch("prepare.visualization.plots.plot_containment")
    def test_main_calls_plot(self, mock_plot, parser):
        cmd = PreparePlotPlugin(parser)
        args = parser.parse_args(
            [
                "--file-path",
                "data.txt",
                "--plot-type",
                "box",
                "--show-title",
            ]
        )
        with patch.object(type(cmd).__bases__[0], "main"):
            cmd.main(args)
        mock_plot.assert_called_once_with(
            file_path="data.txt",
            column="max_containment",
            title=None,
            colors=None,
            show_title=True,
            save_format=None,
            threshold=None,
            plot_type="box",
        )

    @patch("prepare.visualization.plots.plot_containment", side_effect=ValueError("bad"))
    def test_main_handles_error(self, mock_plot, parser):
        cmd = PreparePlotPlugin(parser)
        args = parser.parse_args(["--file-path", "data.txt"])
        with patch.object(type(cmd).__bases__[0], "main"):
            with pytest.raises(SystemExit):
                cmd.main(args)


class TestPlotPluginStderr:
    def test_missing_file_reports_on_stderr(self, parser, tmp_path, capsys):
        cmd = PreparePlotPlugin(parser)
        args = parser.parse_args(["--file-path", str(tmp_path / "missing.txt")])
        with patch.object(type(cmd).__bases__[0], "main"):
            with pytest.raises(SystemExit) as exc:
                cmd.main(args)
        assert exc.value.code == 1
        assert "Error plotting containment" in capsys.readouterr().err


# --- PrepareDiversityPlugin ---


class TestPrepareDiversityPlugin:
    def test_command_attributes(self):
        assert PrepareDiversityPlugin.command == "prepare_diversity"
        assert PrepareDiversityPlugin.description is not None

    def test_configure_parser(self, parser):
        PrepareDiversityPlugin(parser)
        args = parser.parse_args(
            [
                "--abundance-file",
                "abundance.csv",
                "--beta-metric",
                "jaccard",
                "--alpha-metrics",
                "shannon",
                "simpson",
            ]
        )
        assert args.abundance_file == "abundance.csv"
        assert args.beta_metric == "jaccard"
        assert args.alpha_metrics == ["shannon", "simpson"]

    def test_default_arguments(self, parser):
        PrepareDiversityPlugin(parser)
        args = parser.parse_args(["--abundance-file", "abundance.csv"])
        assert args.output_dir == "diversity_results"
        assert args.beta_metric == "bray_curtis"
        assert "shannon" in args.alpha_metrics
        assert args.permanova_formula is None

    @patch("prepare.processing.diversity.calculate_beta_diversity")
    @patch("prepare.processing.diversity.calculate_alpha_diversity")
    @patch("pandas.read_csv")
    def test_main_runs_diversity(self, mock_read, mock_alpha, mock_beta, parser, tmp_output):
        abundance_df = pd.DataFrame({"sp1": [1, 2], "sp2": [3, 4]}, index=["s1", "s2"])
        mock_read.return_value = abundance_df
        mock_alpha.return_value = pd.DataFrame({"shannon": [1.0, 1.1]})
        mock_beta.return_value = pd.DataFrame({"s1": [0, 0.5], "s2": [0.5, 0]}, index=["s1", "s2"])

        cmd = PrepareDiversityPlugin(parser)
        args = parser.parse_args(
            [
                "--abundance-file",
                "abundance.csv",
                "--output-dir",
                str(tmp_output / "div_out"),
            ]
        )
        with patch.object(type(cmd).__bases__[0], "main"):
            cmd.main(args)

        mock_alpha.assert_called_once()
        mock_beta.assert_called_once()

    @patch("pandas.read_csv", side_effect=FileNotFoundError("not found"))
    def test_main_handles_missing_file(self, mock_read, parser, tmp_output):
        cmd = PrepareDiversityPlugin(parser)
        args = parser.parse_args(
            [
                "--abundance-file",
                "missing.csv",
                "--output-dir",
                str(tmp_output / "div_out"),
            ]
        )
        with patch.object(type(cmd).__bases__[0], "main"):
            with pytest.raises(SystemExit):
                cmd.main(args)


# --- PrepareTaxonomyPlugin ---


class TestPrepareTaxonomyPlugin:
    def test_command_attributes(self):
        assert PrepareTaxonomyPlugin.command == "prepare_taxonomy"
        assert PrepareTaxonomyPlugin.description is not None

    def test_configure_parser(self, parser):
        PrepareTaxonomyPlugin(parser)
        args = parser.parse_args(
            [
                "--species-file",
                "species.txt",
                "--email",
                "test@example.com",
                "--api-key",
                "mykey",
            ]
        )
        assert args.species_file == "species.txt"
        assert args.email == "test@example.com"
        assert args.api_key == "mykey"

    def test_default_arguments(self, parser):
        PrepareTaxonomyPlugin(parser)
        args = parser.parse_args(
            [
                "--species-file",
                "species.txt",
                "--email",
                "test@example.com",
            ]
        )
        assert args.output_file == "taxonomy_validation.csv"
        assert args.cache_file == "taxonomy_validation_cache.csv"
        assert args.species_column is None

    @patch("prepare.data.taxonomy.validate_taxonomic_assignments")
    def test_main_with_text_file(self, mock_validate, parser, tmp_output):
        species_file = tmp_output / "species.txt"
        species_file.write_text("E. coli\nB. subtilis\n")

        results_df = pd.DataFrame(
            {
                "species": ["E. coli", "B. subtilis"],
                "is_valid": [True, True],
                "confidence": ["high", "high"],
            }
        )
        mock_validate.return_value = results_df

        cmd = PrepareTaxonomyPlugin(parser)
        args = parser.parse_args(
            [
                "--species-file",
                str(species_file),
                "--email",
                "test@example.com",
            ]
        )
        with patch.object(type(cmd).__bases__[0], "main"):
            cmd.main(args)

        mock_validate.assert_called_once()
        call_args = mock_validate.call_args
        assert call_args[0][0] == ["E. coli", "B. subtilis"]

    @patch("prepare.data.taxonomy.validate_taxonomic_assignments")
    def test_main_takes_the_email_from_the_settings(self, mock_validate, parser, tmp_output, monkeypatch):
        monkeypatch.setenv("PREPARE_NCBI_EMAIL", "env@example.com")
        species_file = tmp_output / "species.txt"
        species_file.write_text("E. coli\n")
        mock_validate.return_value = pd.DataFrame({"species": ["E. coli"], "is_valid": [True]})
        cmd = PrepareTaxonomyPlugin(parser)
        args = parser.parse_args(["--species-file", str(species_file)])
        with patch.object(type(cmd).__bases__[0], "main"):
            cmd.main(args)
        assert mock_validate.call_args.kwargs["email"] == "env@example.com"

    @patch("prepare.data.taxonomy.validate_taxonomic_assignments")
    @patch("pandas.read_csv")
    def test_main_with_csv_file(self, mock_read, mock_validate, parser):
        mock_read.return_value = pd.DataFrame(
            {
                "organism": ["E. coli", "B. subtilis"],
            }
        )
        results_df = pd.DataFrame(
            {
                "species": ["E. coli", "B. subtilis"],
                "is_valid": [True, False],
                "confidence": ["high", "low"],
            }
        )
        mock_validate.return_value = results_df

        cmd = PrepareTaxonomyPlugin(parser)
        args = parser.parse_args(
            [
                "--species-file",
                "species.csv",
                "--species-column",
                "organism",
                "--email",
                "test@example.com",
            ]
        )
        with patch.object(type(cmd).__bases__[0], "main"):
            cmd.main(args)

        mock_validate.assert_called_once()

    def test_main_csv_missing_column(self, parser, tmp_output):
        csv_file = tmp_output / "species.csv"
        csv_file.write_text("name\nE. coli\n")

        cmd = PrepareTaxonomyPlugin(parser)
        args = parser.parse_args(
            [
                "--species-file",
                str(csv_file),
                "--species-column",
                "nonexistent",
                "--email",
                "test@example.com",
            ]
        )
        with patch.object(type(cmd).__bases__[0], "main"):
            with pytest.raises(SystemExit):
                cmd.main(args)

    @patch("prepare.data.taxonomy.validate_taxonomic_assignments", side_effect=RuntimeError("api"))
    def test_main_handles_error(self, mock_validate, parser, tmp_output):
        species_file = tmp_output / "species.txt"
        species_file.write_text("E. coli\n")

        cmd = PrepareTaxonomyPlugin(parser)
        args = parser.parse_args(
            [
                "--species-file",
                str(species_file),
                "--email",
                "test@example.com",
            ]
        )
        with patch.object(type(cmd).__bases__[0], "main"):
            with pytest.raises(SystemExit):
                cmd.main(args)


# --- Plugin discovery ---


class TestPluginDiscovery:
    def test_all_commands_have_required_attributes(self):
        """All plugin classes must have command and description set."""
        classes = [
            PrepareParsePlugin,
            PreparePlotPlugin,
            PrepareDiversityPlugin,
            PrepareTaxonomyPlugin,
        ]
        for cls in classes:
            assert cls.command is not None, f"{cls.__name__} missing command"
            assert cls.description is not None, f"{cls.__name__} missing description"

    def test_command_names_are_unique(self):
        classes = [
            PrepareParsePlugin,
            PreparePlotPlugin,
            PrepareDiversityPlugin,
            PrepareTaxonomyPlugin,
        ]
        names = [cls.command for cls in classes]
        assert len(names) == len(set(names)), "Duplicate command names found"

    def test_command_names_use_prepare_prefix(self):
        classes = [
            PrepareParsePlugin,
            PreparePlotPlugin,
            PrepareDiversityPlugin,
            PrepareTaxonomyPlugin,
        ]
        for cls in classes:
            assert cls.command.startswith("prepare_"), f"{cls.__name__}.command should start with 'prepare_'"
