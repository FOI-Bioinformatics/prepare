"""
Test CLI advanced analysis commands functionality.

Tests for diversity analysis, interactive plotting, and taxonomy commands.
"""

import argparse
from unittest.mock import patch, mock_open
import pandas as pd
import importlib.util
import pytest

from prepare.cli.commands.advanced_analysis import (
    DiversityAnalysisCommand,
    InteractivePlotCommand,
    TaxonomyValidationCommand,
    TaxonomicSummaryCommand,
)

requires_analysis = pytest.mark.skipif(
    not all(importlib.util.find_spec(m) for m in ("sklearn", "scipy")), reason="needs prepare[analysis]"
)

requires_interactive = pytest.mark.skipif(
    not all(importlib.util.find_spec(m) for m in ("plotly", "jinja2")), reason="needs prepare[interactive]"
)


class TestDiversityAnalysisCommand:
    """Test DiversityAnalysisCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = DiversityAnalysisCommand()
        assert command.name == "diversity_analysis"
        assert "diversity" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = DiversityAnalysisCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test required arguments
        args = parser.parse_args(["--abundance-file", "abundance.csv"])
        assert args.abundance_file == "abundance.csv"
        assert args.output_dir == "diversity_results"
        assert args.alpha_metrics == ["shannon", "simpson", "chao1", "observed_species"]
        assert args.beta_metric == "bray_curtis"

    def test_configure_parser_with_options(self):
        """Test parser with optional arguments."""
        command = DiversityAnalysisCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(
            [
                "--abundance-file",
                "abundance.csv",
                "--metadata-file",
                "metadata.csv",
                "--output-dir",
                "custom_output",
                "--alpha-metrics",
                "shannon",
                "simpson",
                "--beta-metric",
                "jaccard",
                "--permanova-formula",
                "treatment + site",
            ]
        )

        assert args.metadata_file == "metadata.csv"
        assert args.output_dir == "custom_output"
        assert args.alpha_metrics == ["shannon", "simpson"]
        assert args.beta_metric == "jaccard"
        assert args.permanova_formula == "treatment + site"

    @requires_analysis
    @patch("prepare.processing.diversity.calculate_alpha_diversity")
    @patch("prepare.processing.diversity.calculate_beta_diversity")
    @patch("prepare.cli.commands.advanced_analysis.read_matrix")
    def test_execute_success_basic(self, mock_read_matrix, mock_beta_div, mock_alpha_div, tmp_path):
        """Test successful execution without PERMANOVA."""
        # Setup mocks
        mock_abundance_df = pd.DataFrame({"sample1": [1, 2], "sample2": [3, 4]})
        mock_read_matrix.return_value = mock_abundance_df

        mock_alpha_result = pd.DataFrame({"shannon": [1.5, 2.0]})
        mock_alpha_div.return_value = mock_alpha_result

        mock_beta_result = pd.DataFrame([[0, 0.5], [0.5, 0]])
        mock_beta_div.return_value = mock_beta_result

        command = DiversityAnalysisCommand()
        args = argparse.Namespace(
            abundance_file="abundance.csv",
            metadata_file=None,
            output_dir=str(tmp_path / "out"),
            alpha_metrics=["shannon"],
            beta_metric="bray_curtis",
            permanova_formula=None,
        )

        result = command.execute(args)

        assert result == 0
        mock_alpha_div.assert_called_once_with(mock_abundance_df, ["shannon"])
        mock_beta_div.assert_called_once_with(mock_abundance_df, "bray_curtis")
        # Alpha and beta results, and no temporary file left beside them.
        assert sorted(p.name for p in (tmp_path / "out").iterdir()) == [
            "alpha_diversity.csv",
            "beta_diversity_bray_curtis.csv",
        ]

    @requires_analysis
    @patch("prepare.processing.diversity.calculate_alpha_diversity")
    @patch("prepare.processing.diversity.calculate_beta_diversity")
    @patch("prepare.processing.diversity.perform_permanova")
    @patch("pandas.read_csv")
    def test_execute_success_with_permanova(
        self, mock_read_csv, mock_permanova, mock_beta_div, mock_alpha_div, tmp_path
    ):
        """Test successful execution with PERMANOVA."""
        # Setup mocks
        mock_abundance_df = pd.DataFrame({"sample1": [1, 2], "sample2": [3, 4]})
        mock_metadata_df = pd.DataFrame({"treatment": ["A", "B"]})
        mock_read_csv.side_effect = [mock_abundance_df, mock_metadata_df]

        mock_alpha_result = pd.DataFrame({"shannon": [1.5, 2.0]})
        mock_alpha_div.return_value = mock_alpha_result

        mock_beta_result = pd.DataFrame([[0, 0.5], [0.5, 0]])
        mock_beta_div.return_value = mock_beta_result

        mock_permanova_result = {"treatment": {"F_statistic": 2.5, "p_value": 0.03, "significant": True}}
        mock_permanova.return_value = mock_permanova_result

        command = DiversityAnalysisCommand()
        args = argparse.Namespace(
            abundance_file="abundance.csv",
            metadata_file="metadata.csv",
            output_dir=str(tmp_path / "out"),
            alpha_metrics=["shannon"],
            beta_metric="bray_curtis",
            permanova_formula="treatment",
        )

        result = command.execute(args)

        assert result == 0
        mock_permanova.assert_called_once()
        report = (tmp_path / "out" / "permanova_results.txt").read_text()
        assert "Variable: treatment" in report
        assert "P-value: 0.0300" in report

    @requires_analysis
    @patch("pandas.read_csv")
    def test_execute_failure(self, mock_read_csv):
        """Test execution failure."""
        mock_read_csv.side_effect = Exception("File not found")

        command = DiversityAnalysisCommand()
        args = argparse.Namespace(
            abundance_file="nonexistent.csv",
            metadata_file=None,
            output_dir="test_output",
            alpha_metrics=["shannon"],
            beta_metric="bray_curtis",
            permanova_formula=None,
        )

        result = command.execute(args)
        assert result == 1

    @requires_analysis
    def test_execute_permanova_with_text_metadata_column(self, tmp_path):
        """Regression test: a text metadata column must survive the metadata load.

        read_matrix drops non-numeric columns, which previously made
        --permanova-formula fail with a KeyError because the group column
        (e.g. "treatment") was stripped from the metadata before PERMANOVA ran.
        The metadata load must use read_table instead, which keeps all columns.
        Uses real (unmocked) file I/O and diversity/PERMANOVA functions.
        """
        abundance_file = tmp_path / "abundance.txt"
        abundance_file.write_text(
            "\tsp1\tsp2\tsp3\n"
            "SRR1\t10\t5\t2\n"
            "SRR2\t12\t4\t3\n"
            "SRR3\t9\t6\t1\n"
            "SRR4\t2\t8\t12\n"
            "SRR5\t3\t7\t10\n"
            "SRR6\t1\t9\t13\n"
        )

        metadata_file = tmp_path / "metadata.txt"
        metadata_file.write_text(
            "Run_ID\ttreatment\n"
            "SRR1\tcontrol\n"
            "SRR2\tcontrol\n"
            "SRR3\tcontrol\n"
            "SRR4\ttreated\n"
            "SRR5\ttreated\n"
            "SRR6\ttreated\n"
        )

        output_dir = tmp_path / "diversity_results"

        command = DiversityAnalysisCommand()
        args = argparse.Namespace(
            abundance_file=str(abundance_file),
            metadata_file=str(metadata_file),
            output_dir=str(output_dir),
            alpha_metrics=["shannon"],
            beta_metric="bray_curtis",
            permanova_formula="treatment",
        )

        result = command.execute(args)

        assert result == 0
        assert (output_dir / "alpha_diversity.csv").exists()
        assert (output_dir / "beta_diversity_bray_curtis.csv").exists()
        assert (output_dir / "permanova_results.txt").exists()
        permanova_text = (output_dir / "permanova_results.txt").read_text()
        assert "Variable: treatment" in permanova_text


class TestInteractivePlotCommand:
    """Test InteractivePlotCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = InteractivePlotCommand()
        assert command.name == "interactive_plot"
        assert "interactive" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = InteractivePlotCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test required arguments
        args = parser.parse_args(["--data-file", "data.csv", "--plot-type", "pca"])
        assert args.data_file == "data.csv"
        assert args.plot_type == "pca"
        assert args.no_show is False

    def test_configure_parser_with_options(self):
        """Test parser with optional arguments."""
        command = InteractivePlotCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(
            [
                "--data-file",
                "data.csv",
                "--metadata-file",
                "metadata.csv",
                "--plot-type",
                "heatmap",
                "--color-by",
                "treatment",
                "--size-by",
                "weight",
                "--output-file",
                "plot.html",
                "--title",
                "My Plot",
                "--no-show",
            ]
        )

        assert args.metadata_file == "metadata.csv"
        assert args.color_by == "treatment"
        assert args.size_by == "weight"
        assert args.output_file == "plot.html"
        assert args.title == "My Plot"
        assert args.no_show is True

    @requires_interactive
    @patch("prepare.cli.commands.advanced_analysis.read_table")
    @patch("prepare.cli.commands.advanced_analysis.read_matrix")
    @patch("prepare.visualization.interactive.create_interactive_pca")
    def test_execute_pca_plot(self, mock_create_pca, mock_read_matrix, mock_read_table):
        """Test PCA plot creation."""
        mock_data_df = pd.DataFrame({"gene1": [1, 2], "gene2": [3, 4]})
        mock_metadata_df = pd.DataFrame({"treatment": ["A", "B"]})
        mock_read_matrix.return_value = mock_data_df
        mock_read_table.return_value = mock_metadata_df

        command = InteractivePlotCommand()
        args = argparse.Namespace(
            data_file="data.csv",
            metadata_file="metadata.csv",
            plot_type="pca",
            color_by="treatment",
            size_by=None,
            output_file="pca.html",
            title="PCA Test",
            no_show=True,
        )

        result = command.execute(args)

        assert result == 0
        mock_create_pca.assert_called_once_with(
            mock_data_df,
            metadata=mock_metadata_df,
            color_by="treatment",
            size_by=None,
            title="PCA Test",
            output_file="pca.html",
            show_plot=False,
        )

    @requires_interactive
    @patch("prepare.cli.commands.advanced_analysis.read_matrix")
    @patch("prepare.visualization.interactive.create_interactive_heatmap")
    def test_execute_heatmap_plot(self, mock_create_heatmap, mock_read_matrix):
        """Test heatmap plot creation."""
        mock_data_df = pd.DataFrame({"gene1": [1, 2], "gene2": [3, 4]})
        mock_read_matrix.return_value = mock_data_df

        command = InteractivePlotCommand()
        args = argparse.Namespace(
            data_file="data.csv",
            metadata_file=None,
            plot_type="heatmap",
            color_by=None,
            size_by=None,
            output_file=None,
            title=None,
            no_show=False,
        )

        result = command.execute(args)

        assert result == 0
        mock_create_heatmap.assert_called_once_with(
            mock_data_df, sample_metadata=None, title="Interactive Heatmap", output_file=None, show_plot=True
        )

    @requires_interactive
    @patch("prepare.cli.commands.advanced_analysis.read_table")
    @patch("prepare.cli.commands.advanced_analysis.read_matrix")
    @patch("prepare.processing.diversity.calculate_alpha_diversity")
    @patch("prepare.visualization.interactive.create_diversity_comparison_plot")
    def test_execute_diversity_plot(self, mock_create_diversity, mock_calc_alpha, mock_read_matrix, mock_read_table):
        """Test diversity plot creation."""
        mock_data_df = pd.DataFrame({"gene1": [1, 2], "gene2": [3, 4]})
        mock_metadata_df = pd.DataFrame({"treatment": ["A", "B"]})
        mock_alpha_div = pd.DataFrame({"shannon": [1.5, 2.0]})

        mock_read_matrix.return_value = mock_data_df
        mock_read_table.return_value = mock_metadata_df
        mock_calc_alpha.return_value = mock_alpha_div

        command = InteractivePlotCommand()
        args = argparse.Namespace(
            data_file="data.csv",
            metadata_file="metadata.csv",
            plot_type="diversity",
            color_by="treatment",
            size_by=None,
            output_file=None,
            title=None,
            no_show=False,
        )

        result = command.execute(args)

        assert result == 0
        mock_calc_alpha.assert_called_once_with(mock_data_df)
        mock_create_diversity.assert_called_once_with(
            mock_alpha_div,
            mock_metadata_df,
            group_by="treatment",
            title="Diversity Comparison",
            output_file=None,
            show_plot=True,
        )

    @requires_interactive
    @patch("pandas.read_csv")
    def test_execute_diversity_plot_missing_color_by(self, mock_read_csv):
        """Test diversity plot with missing color_by parameter."""
        mock_data_df = pd.DataFrame({"gene1": [1, 2], "gene2": [3, 4]})
        mock_read_csv.return_value = mock_data_df

        command = InteractivePlotCommand()
        args = argparse.Namespace(
            data_file="data.csv",
            metadata_file=None,
            plot_type="diversity",
            color_by=None,
            size_by=None,
            output_file=None,
            title=None,
            no_show=False,
        )

        result = command.execute(args)
        assert result == 1

    @requires_interactive
    @patch("prepare.visualization.interactive.create_interactive_tsne")
    @patch("prepare.cli.commands.advanced_analysis.read_matrix")
    def test_execute_tsne_plot(self, mock_read, mock_tsne, tmp_path):
        import pandas as pd

        mock_read.return_value = pd.DataFrame({"a": [1, 2, 3], "b": [3, 2, 1]}, index=["s1", "s2", "s3"])
        args = argparse.Namespace(
            data_file="m.csv",
            metadata_file=None,
            plot_type="tsne",
            color_by=None,
            size_by=None,
            output_file=str(tmp_path / "t.html"),
            title=None,
            no_show=True,
        )
        assert InteractivePlotCommand().execute(args) == 0
        assert mock_tsne.call_args.kwargs["output_file"] == str(tmp_path / "t.html")
        assert mock_tsne.call_args.kwargs["show_plot"] is False


class TestTaxonomyValidationCommand:
    """Test TaxonomyValidationCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = TaxonomyValidationCommand()
        assert command.name == "validate_taxonomy"
        assert "taxonomy" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = TaxonomyValidationCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test required arguments
        args = parser.parse_args(["--species-file", "species.txt", "--email", "test@example.com"])
        assert args.species_file == "species.txt"
        assert args.email == "test@example.com"
        assert args.output_file == "taxonomy_validation.csv"
        assert args.cache_file == "taxonomy_cache.csv"

    @patch("prepare.cli.commands.advanced_analysis.validate_taxonomic_assignments")
    def test_execute_success_text_file(self, mock_validate, capsys):
        """Test successful execution with text file."""
        mock_results = pd.DataFrame(
            {
                "species": ["Escherichia coli", "Salmonella enterica"],
                "is_valid": [True, True],
                "confidence": ["high", "high"],
            }
        )
        mock_validate.return_value = mock_results

        command = TaxonomyValidationCommand()
        args = argparse.Namespace(
            species_file="species.txt",
            species_column=None,
            email="test@example.com",
            api_key=None,
            output_file="validation.csv",
            cache_file="cache.csv",
        )

        with patch("builtins.open", mock_open(read_data="Escherichia coli\nSalmonella enterica\n")):
            result = command.execute(args)

        assert result == 0
        mock_validate.assert_called_once()
        out = capsys.readouterr().out
        assert "Taxonomy Validation Summary:" in out
        assert "Total species: 2" in out

    @patch("pandas.read_csv")
    @patch("prepare.cli.commands.advanced_analysis.validate_taxonomic_assignments")
    def test_execute_success_csv_file_with_column(self, mock_validate, mock_read_csv, capsys):
        """Test successful execution with CSV file and specified column."""
        mock_df = pd.DataFrame({"species_name": ["Escherichia coli", "Salmonella enterica"], "abundance": [100, 50]})
        mock_read_csv.return_value = mock_df

        mock_results = pd.DataFrame(
            {
                "species": ["Escherichia coli", "Salmonella enterica"],
                "is_valid": [True, False],
                "confidence": ["high", "low"],
            }
        )
        mock_validate.return_value = mock_results

        command = TaxonomyValidationCommand()
        args = argparse.Namespace(
            species_file="species.csv",
            species_column="species_name",
            email="test@example.com",
            api_key="api123",
            output_file="validation.csv",
            cache_file="cache.csv",
        )

        result = command.execute(args)

        assert result == 0
        assert "Valid species: 1 (50.0%)" in capsys.readouterr().out
        mock_validate.assert_called_once_with(
            ["Escherichia coli", "Salmonella enterica"],
            email="test@example.com",
            api_key="api123",
            output_file="validation.csv",
            cache_file="cache.csv",
        )

    @patch("pandas.read_csv")
    def test_execute_csv_missing_column(self, mock_read_csv):
        """Test execution with CSV file missing specified column."""
        mock_df = pd.DataFrame({"other_column": ["data1", "data2"]})
        mock_read_csv.return_value = mock_df

        command = TaxonomyValidationCommand()
        args = argparse.Namespace(
            species_file="species.csv",
            species_column="species_name",
            email="test@example.com",
            api_key=None,
            output_file="validation.csv",
            cache_file="cache.csv",
        )

        result = command.execute(args)
        assert result == 1

    @patch("pandas.read_csv")
    @patch("prepare.cli.commands.advanced_analysis.validate_taxonomic_assignments")
    def test_execute_success_csv_first_column(self, mock_validate, mock_read_csv, capsys):
        """Test successful execution with CSV file using first column."""
        mock_df = pd.DataFrame(
            {"first_col": ["Escherichia coli", "Salmonella enterica"], "other_col": ["data1", "data2"]}
        )
        mock_read_csv.return_value = mock_df

        mock_results = pd.DataFrame(
            {
                "species": ["Escherichia coli", "Salmonella enterica"],
                "is_valid": [True, True],
                "confidence": ["high", "high"],
            }
        )
        mock_validate.return_value = mock_results

        command = TaxonomyValidationCommand()
        args = argparse.Namespace(
            species_file="species.csv",
            species_column=None,  # Will use first column
            email="test@example.com",
            api_key=None,
            output_file="validation.csv",
            cache_file="cache.csv",
        )

        result = command.execute(args)

        assert result == 0
        assert "Total species: 2" in capsys.readouterr().out


class TestTaxonomicSummaryCommand:
    """Test TaxonomicSummaryCommand."""

    def test_command_properties(self):
        """Test command name and help."""
        command = TaxonomicSummaryCommand()
        assert command.name == "taxonomic_summary"
        assert "taxonomic" in command.help.lower()

    def test_configure_parser(self):
        """Test parser configuration."""
        command = TaxonomicSummaryCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        # Test required arguments
        args = parser.parse_args(["--abundance-file", "abundance.csv", "--taxonomy-file", "taxonomy.csv"])
        assert args.abundance_file == "abundance.csv"
        assert args.taxonomy_file == "taxonomy.csv"
        assert args.output_dir == "taxonomic_summaries"
        assert args.levels == ["phylum", "class", "order", "family", "genus"]
        assert args.min_abundance == 0.001

    def test_configure_parser_with_options(self):
        """Test parser with optional arguments."""
        command = TaxonomicSummaryCommand()
        parser = argparse.ArgumentParser()
        command.configure_parser(parser)

        args = parser.parse_args(
            [
                "--abundance-file",
                "abundance.csv",
                "--taxonomy-file",
                "taxonomy.csv",
                "--output-dir",
                "custom_output",
                "--levels",
                "phylum",
                "genus",
                "--min-abundance",
                "0.01",
            ]
        )

        assert args.output_dir == "custom_output"
        assert args.levels == ["phylum", "genus"]
        assert args.min_abundance == 0.01

    @patch("prepare.cli.commands.advanced_analysis.read_records")
    @patch("prepare.cli.commands.advanced_analysis.read_matrix")
    @patch("prepare.cli.commands.advanced_analysis.analyze_taxonomic_composition")
    def test_execute_success(self, mock_analyze, mock_read_matrix, mock_read_records, capsys):
        """Test successful execution."""
        mock_abundance_df = pd.DataFrame({"sample1": [10, 20], "sample2": [15, 25]})
        mock_taxonomy_df = pd.DataFrame({"species": ["Species1", "Species2"], "phylum": ["Phylum1", "Phylum2"]})
        mock_read_matrix.return_value = mock_abundance_df
        mock_read_records.return_value = mock_taxonomy_df

        mock_summaries = {
            "phylum": pd.DataFrame({"Phylum1": [10, 15], "Phylum2": [20, 25]}),
            "genus": pd.DataFrame({"Genus1": [5, 8], "Genus2": [25, 32]}),
        }
        mock_analyze.return_value = mock_summaries

        command = TaxonomicSummaryCommand()
        args = argparse.Namespace(
            abundance_file="abundance.csv",
            taxonomy_file="taxonomy.csv",
            output_dir="summaries",
            levels=["phylum", "genus"],
            min_abundance=0.001,
        )

        result = command.execute(args)

        assert result == 0
        mock_analyze.assert_called_once_with(
            mock_abundance_df, mock_taxonomy_df, levels=["phylum", "genus"], output_dir="summaries"
        )
        out = capsys.readouterr().out
        assert "Taxonomic Summary Results:" in out
        assert "Phylum: 2 taxa across 2 samples" in out

    @patch("pandas.read_csv")
    def test_execute_failure(self, mock_read_csv):
        """Test execution failure."""
        mock_read_csv.side_effect = Exception("File not found")

        command = TaxonomicSummaryCommand()
        args = argparse.Namespace(
            abundance_file="nonexistent.csv",
            taxonomy_file="taxonomy.csv",
            output_dir="summaries",
            levels=["phylum"],
            min_abundance=0.001,
        )

        result = command.execute(args)
        assert result == 1


class TestAdvancedAnalysisIntegration:
    """Integration tests for advanced analysis commands."""

    def test_all_commands_inherit_base_command(self):
        """Test that all command classes inherit from BaseCommand."""
        from prepare.cli.base import BaseCommand

        commands = [
            DiversityAnalysisCommand(),
            InteractivePlotCommand(),
            TaxonomyValidationCommand(),
            TaxonomicSummaryCommand(),
        ]

        for command in commands:
            assert isinstance(command, BaseCommand)
            assert hasattr(command, "name")
            assert hasattr(command, "help")
            assert hasattr(command, "configure_parser")
            assert hasattr(command, "execute")

    def test_command_names_unique(self):
        """Test that all command names are unique."""
        commands = [
            DiversityAnalysisCommand(),
            InteractivePlotCommand(),
            TaxonomyValidationCommand(),
            TaxonomicSummaryCommand(),
        ]

        names = [cmd.name for cmd in commands]
        assert len(names) == len(set(names))

    def test_all_parsers_configured(self):
        """Test that all commands can configure parsers without error."""
        commands = [
            DiversityAnalysisCommand(),
            InteractivePlotCommand(),
            TaxonomyValidationCommand(),
            TaxonomicSummaryCommand(),
        ]

        for command in commands:
            parser = argparse.ArgumentParser()
            # Should not raise exception
            command.configure_parser(parser)

    @requires_analysis
    @patch("prepare.cli.commands.advanced_analysis.logger")
    def test_logging_behavior(self, mock_logger):
        """Test that commands use logging appropriately."""
        command = DiversityAnalysisCommand()

        # Mock the pandas import and other dependencies
        with patch("pandas.read_csv", side_effect=Exception("Test error")):
            args = argparse.Namespace(
                abundance_file="test.csv",
                metadata_file=None,
                output_dir="output",
                alpha_metrics=["shannon"],
                beta_metric="bray_curtis",
                permanova_formula=None,
            )

            result = command.execute(args)

            assert result == 1
            mock_logger.error.assert_called()


class TestTaxonomicSummaryWithEnrichMap:
    def test_parsed_containment_plus_taxonomy_tsv(self, tmp_path):
        from prepare.cli.commands.advanced_analysis import TaxonomicSummaryCommand

        cont = tmp_path / "parsed_containment.txt"
        cont.write_text(
            "\tGCF_1\tGCF_2\tmax_containment\tmax_containment_annotation\nS1\t0.9\t0.1\t0.9\tGCF_1\nS2\t0.0\t0.7\t0.7\tGCF_2\n"  # noqa: E501
        )
        tax = tmp_path / "taxonomy.tsv"
        tax.write_text(
            "genome_id\tspecies\tgenus\tfamily\torder\tclass_name\tphylum\torganism\ttax_id\n"
            "GCF_1\tWolbachia pipientis\tWolbachia\tAnaplasmataceae\tRickettsiales\tAlphaproteobacteria\tPseudomonadota\t\t\n"  # noqa: E501
            "GCF_2\tWolbachia sp947251865\tWolbachia\tAnaplasmataceae\tRickettsiales\tAlphaproteobacteria\tPseudomonadota\t\t\n"  # noqa: E501
        )
        args = argparse.Namespace(
            abundance_file=str(cont),
            taxonomy_file=str(tax),
            output_dir=str(tmp_path / "out"),
            levels=["genus"],
            min_abundance=0.0,
        )
        assert TaxonomicSummaryCommand().execute(args) == 0
        out = (tmp_path / "out" / "taxonomy_summary_genus.csv").read_text().splitlines()
        assert out[0].startswith(",Wolbachia")


if __name__ == "__main__":
    pytest.main([__file__])
