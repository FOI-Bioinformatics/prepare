"""
Advanced analysis commands for PREPARE CLI.

This module provides CLI commands for diversity analysis,
interactive visualization, and taxonomic validation.
"""

import logging
from pathlib import Path

from prepare.cli.base import finite_float, BaseCommand
from prepare.core.exceptions import ConfigurationError, PrepareError, ValidationError
from prepare.core.optional import require
from prepare.core.settings import require_email, setting_for
from prepare.data.file_io import open_atomic, write_csv
from prepare.data.defaults import read_matrix, read_records, read_table
from prepare.data.taxonomy import (
    validate_taxonomic_assignments,
    analyze_taxonomic_composition,
)

logger = logging.getLogger(__name__)


class DiversityAnalysisCommand(BaseCommand):
    """Command for calculating diversity metrics."""

    @property
    def name(self) -> str:
        return "diversity_analysis"

    @property
    def help(self) -> str:
        return "Calculate alpha and beta diversity metrics"

    @property
    def group(self) -> str:
        return "Analysis"

    def configure_parser(self, parser):
        parser.add_argument(
            "--abundance-file",
            required=True,
            help="CSV file with abundance matrix (samples x species)",
        )
        parser.add_argument("--metadata-file", help="CSV file with sample metadata")
        parser.add_argument(
            "--output-dir",
            default="diversity_results",
            help="Output directory for results",
        )
        parser.add_argument(
            "--alpha-metrics",
            nargs="+",
            default=["shannon", "simpson", "chao1", "observed_species"],
            help="Alpha diversity metrics to calculate",
        )
        parser.add_argument(
            "--beta-metric",
            default="bray_curtis",
            choices=["bray_curtis", "jaccard", "euclidean", "manhattan"],
            help="Beta diversity metric",
        )
        parser.add_argument("--permanova-formula", help="PERMANOVA formula (e.g., 'treatment + site')")

    def execute(self, args):
        # Diversity analysis belongs to the analysis extra; check before any output is written.
        require("sklearn", "analysis", "Diversity analysis")
        from prepare.processing.diversity import (
            calculate_alpha_diversity,
            calculate_beta_diversity,
            check_metadata_covers_samples,
            perform_permanova,
        )

        try:
            import pandas as pd

            # Load data
            logger.info("Loading abundance data...")
            abundance_df = read_matrix(args.abundance_file)

            metadata_df = None
            if args.metadata_file:
                logger.info("Loading metadata...")
                metadata_df = read_table(args.metadata_file)
                if args.permanova_formula:
                    # Every sample needs a metadata row; check before any output is written.
                    check_metadata_covers_samples(abundance_df.index, metadata_df)

            # Create output directory
            output_dir = Path(args.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            # Calculate alpha diversity
            logger.info("Calculating alpha diversity...")
            alpha_div = calculate_alpha_diversity(abundance_df, args.alpha_metrics)

            alpha_output = output_dir / "alpha_diversity.csv"
            write_csv(alpha_div, alpha_output)
            logger.info(f"Alpha diversity results saved to {alpha_output}")

            # Calculate beta diversity
            logger.info(f"Calculating beta diversity using {args.beta_metric}...")
            beta_div = calculate_beta_diversity(abundance_df, args.beta_metric)
            assert isinstance(beta_div, pd.DataFrame)

            beta_output = output_dir / f"beta_diversity_{args.beta_metric}.csv"
            write_csv(beta_div, beta_output)
            logger.info(f"Beta diversity results saved to {beta_output}")

            # PERMANOVA if requested
            if args.permanova_formula and metadata_df is not None:
                logger.info("Performing PERMANOVA...")
                permanova_results = perform_permanova(beta_div, metadata_df, args.permanova_formula)

                permanova_output = output_dir / "permanova_results.txt"
                with open_atomic(permanova_output) as f:
                    f.write("PERMANOVA Results\n")
                    f.write("=================\n\n")
                    for variable, results in permanova_results.items():
                        f.write(f"Variable: {variable}\n")
                        f.write(f"F-statistic: {results['F_statistic']:.4f}\n")
                        f.write(f"P-value: {results['p_value']:.4f}\n")
                        f.write(f"Significant: {results['significant']}\n\n")

                logger.info(f"PERMANOVA results saved to {permanova_output}")

            logger.info("Diversity analysis completed successfully!")
            return 0

        except (ConfigurationError, ValidationError):
            raise
        except Exception as e:
            logger.error(f"Diversity analysis failed: {e}")
            return 1


class InteractivePlotCommand(BaseCommand):
    """Command for creating interactive plots."""

    @property
    def name(self) -> str:
        return "interactive_plot"

    @property
    def help(self) -> str:
        return "Create interactive visualizations"

    @property
    def group(self) -> str:
        return "Analysis"

    def configure_parser(self, parser):
        parser.add_argument("--data-file", required=True, help="CSV file with data matrix")
        parser.add_argument("--metadata-file", help="CSV file with sample metadata")
        parser.add_argument(
            "--plot-type",
            required=True,
            choices=["pca", "tsne", "heatmap", "diversity"],
            help="Type of interactive plot to create",
        )
        parser.add_argument("--color-by", help="Metadata column for coloring points")
        parser.add_argument("--size-by", help="Metadata column for sizing points")
        parser.add_argument("--output-file", help="HTML file to save interactive plot")
        parser.add_argument("--title", help="Plot title")
        parser.add_argument("--no-show", action="store_true", help="Don't display plot in browser")

    def execute(self, args):
        require("plotly", "interactive", "An interactive plot")
        if args.plot_type in ("pca", "tsne", "heatmap"):
            require("sklearn" if args.plot_type != "heatmap" else "scipy", "analysis", f"A {args.plot_type} plot")
        from prepare.processing.diversity import calculate_alpha_diversity
        from prepare.visualization.interactive import (
            create_diversity_comparison_plot,
            create_interactive_heatmap,
            create_interactive_pca,
            create_interactive_tsne,
        )

        try:
            # Load data
            logger.info("Loading data...")
            data_df = read_matrix(args.data_file)

            metadata_df = None
            if args.metadata_file:
                logger.info("Loading metadata...")
                metadata_df = read_table(args.metadata_file)

            # Create plot based on type
            show_plot = not args.no_show

            if args.plot_type == "pca":
                logger.info("Creating interactive PCA plot...")
                create_interactive_pca(
                    data_df,
                    metadata=metadata_df,
                    color_by=args.color_by,
                    size_by=args.size_by,
                    title=args.title or "Interactive PCA Plot",
                    output_file=args.output_file,
                    show_plot=show_plot,
                )

            elif args.plot_type == "heatmap":
                logger.info("Creating interactive heatmap...")
                create_interactive_heatmap(
                    data_df,
                    sample_metadata=metadata_df,
                    title=args.title or "Interactive Heatmap",
                    output_file=args.output_file,
                    show_plot=show_plot,
                )

            elif args.plot_type == "diversity":
                if not args.color_by:
                    raise PrepareError("--color-by required for diversity plots")

                logger.info("Creating diversity comparison plot...")
                # Calculate alpha diversity first
                alpha_div = calculate_alpha_diversity(data_df)

                create_diversity_comparison_plot(
                    alpha_div,
                    metadata_df,
                    group_by=args.color_by,
                    title=args.title or "Diversity Comparison",
                    output_file=args.output_file,
                    show_plot=show_plot,
                )

            elif args.plot_type == "tsne":
                logger.info("Creating interactive t-SNE plot...")
                create_interactive_tsne(
                    data_df,
                    metadata_df,
                    color_by=args.color_by,
                    title=args.title or "Interactive t-SNE Plot",
                    output_file=args.output_file,
                    show_plot=show_plot,
                )

            logger.info("Interactive plot created successfully!")
            return 0

        except ConfigurationError:
            raise
        except Exception as e:
            logger.error(f"Interactive plot creation failed: {e}")
            return 1


class TaxonomyValidationCommand(BaseCommand):
    """Command for validating taxonomic assignments."""

    @property
    def name(self) -> str:
        return "validate_taxonomy"

    @property
    def help(self) -> str:
        return "Validate species names against NCBI taxonomy"

    @property
    def group(self) -> str:
        return "Analysis"

    def configure_parser(self, parser):
        parser.add_argument(
            "--species-file",
            required=True,
            help="Text file with species names (one per line) or CSV with species column",  # noqa: E501
        )
        parser.add_argument("--species-column", help="Column name containing species (for CSV files)")
        parser.add_argument(
            "--email",
            default=None,
            help="Email address for NCBI API access (default: PREPARE_NCBI_EMAIL or config [runtime] ncbi_email)",
        )
        parser.add_argument(
            "--api-key", help="NCBI API key for increased rate limits (default: PREPARE_NCBI_API_KEY or NCBI_API_KEY)"
        )
        parser.add_argument(
            "--output-file",
            default="taxonomy_validation.csv",
            help="Output file for validation results",
        )
        parser.add_argument(
            "--cache-file",
            default="taxonomy_validation_cache.csv",
            help="Cache file for the NCBI taxonomy lookups (not the GTDB cache of enrich_taxonomy)",
        )

    def execute(self, args):
        email = require_email(args)
        try:
            import pandas as pd

            # Load species list
            logger.info("Loading species list...")
            if args.species_file.endswith(".csv"):
                df = pd.read_csv(args.species_file)
                if args.species_column:
                    if args.species_column not in df.columns:
                        raise PrepareError(f"Column '{args.species_column}' not found")
                    species_list = df[args.species_column].dropna().tolist()
                else:
                    # Assume first column contains species
                    species_list = df.iloc[:, 0].dropna().tolist()
            else:
                # Text file
                with open(args.species_file, "r") as f:
                    species_list = [line.strip() for line in f if line.strip()]

            logger.info(f"Loaded {len(species_list)} species names")

            # Validate taxonomy
            logger.info("Validating species names against NCBI taxonomy...")
            results_df = validate_taxonomic_assignments(
                species_list,
                email=email,
                api_key=setting_for(args, "ncbi_api_key"),
                output_file=args.output_file,
                cache_file=args.cache_file,
            )

            # Print summary
            valid_count = results_df["is_valid"].sum()
            total_count = len(results_df)

            self.emit("\nTaxonomy Validation Summary:")
            self.emit("============================")
            self.emit(f"Total species: {total_count}")
            self.emit(f"Valid species: {valid_count} " f"({valid_count / total_count * 100:.1f}%)")
            self.emit(
                f"Invalid species: {total_count - valid_count} ({(total_count - valid_count) / total_count * 100:.1f}%)"  # noqa: E501
            )

            # Show confidence distribution
            confidence_counts = results_df["confidence"].value_counts()
            self.emit("\nConfidence distribution:")
            for conf, count in confidence_counts.items():
                self.emit(f"  {conf}: {count}")

            logger.info("Taxonomy validation completed successfully!")
            return 0

        except Exception as e:
            logger.error(f"Taxonomy validation failed: {e}")
            return 1


class TaxonomicSummaryCommand(BaseCommand):
    """Command for creating taxonomic summaries."""

    @property
    def name(self) -> str:
        return "taxonomic_summary"

    @property
    def help(self) -> str:
        return "Create taxonomic summaries at different levels"

    @property
    def group(self) -> str:
        return "Analysis"

    def configure_parser(self, parser):
        parser.add_argument(
            "--abundance-file",
            required=True,
            help="CSV file with abundance matrix (samples x species)",
        )
        parser.add_argument(
            "--taxonomy-file",
            required=True,
            help="Taxonomy table from enrich_taxonomy (TSV, genome_id and ranks) or validate_taxonomy (CSV)",
        )
        parser.add_argument(
            "--output-dir",
            default="taxonomic_summaries",
            help="Output directory for summary files",
        )
        parser.add_argument(
            "--levels",
            nargs="+",
            default=["phylum", "class", "order", "family", "genus"],
            help="Taxonomic levels to summarize",
        )
        parser.add_argument(
            "--min-abundance",
            type=finite_float,
            default=0.001,
            help="Minimum abundance threshold",
        )

    def execute(self, args):
        try:
            # Load data
            logger.info("Loading abundance data...")
            abundance_df = read_matrix(args.abundance_file)

            logger.info("Loading taxonomy data...")
            taxonomy_df = read_records(args.taxonomy_file)

            # Create taxonomic summaries
            logger.info("Creating taxonomic summaries...")
            summaries = analyze_taxonomic_composition(
                abundance_df,
                taxonomy_df,
                levels=args.levels,
                output_dir=args.output_dir,
                min_abundance=args.min_abundance,
            )

            # Print summary statistics
            self.emit("\nTaxonomic Summary Results:")
            self.emit("==========================")
            for level, summary_df in summaries.items():
                n_taxa = summary_df.shape[1]
                n_samples = summary_df.shape[0]
                self.emit(f"{level.title()}: {n_taxa} taxa across {n_samples} samples")

            logger.info("Taxonomic summaries created successfully!")
            return 0

        except Exception as e:
            logger.error(f"Taxonomic summary creation failed: {e}")
            return 1
