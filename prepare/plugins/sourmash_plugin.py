"""
Sourmash CLI plugin integration for PREPARE.

This module provides adapter classes that bridge PREPARE commands to
sourmash's CommandLinePlugin interface, making PREPARE functionality
available as `sourmash scripts prepare_*` subcommands.
"""

import logging
import sys
from typing import TYPE_CHECKING, List

from prepare.cli.base import finite_float
from prepare.core.constants import DEFAULT_REGISTRY_MAX_SCREENED
from prepare.core.optional import require
from prepare.utils.logging import setup_logging

logger = logging.getLogger(__name__)

# sourmash loads this module through its plugin entry points, so sourmash is present whenever
# the plugin is used; importing it elsewhere without sourmash is an error naming the extra.
if TYPE_CHECKING:
    from sourmash.plugins import CommandLinePlugin
else:
    CommandLinePlugin = require("sourmash.plugins", "sourmash", "The PREPARE sourmash plugin").CommandLinePlugin


class PrepareParsePlugin(CommandLinePlugin):
    """Parse containment data from branchwater match files."""

    command = "prepare_parse"
    description = "Parse containment data from branchwater match files"

    def __init__(self, parser):
        super().__init__(parser)
        parser.add_argument(
            "--matches-folder",
            default="matches",
            help="Folder containing containment match files",
        )
        parser.add_argument(
            "--parsed-containment-file",
            default="parsed_containment.txt",
            help="File where the parsed containment will be stored",
        )
        parser.add_argument(
            "--summary-containment-file",
            default="top_containments.txt",
            help="File where the containment summary will be stored",
        )
        parser.add_argument(
            "--step-size",
            default=0.1,
            type=finite_float,
            help="Size of steps for the containment thresholds",
        )
        parser.add_argument(
            "--details-file",
            default=None,
            help="File for the cANI and sample metadata details table "
            "(default: <parsed-containment-file stem>_details.tsv)",
        )
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")
        parser.add_argument(
            "--registry-max-screened",
            type=int,
            default=DEFAULT_REGISTRY_MAX_SCREENED,
            help="Most screening entries to keep in the registry per genome (the match CSVs keep them all)",
        )

    def main(self, args):
        super().main(args)
        setup_logging(level=logging.INFO)
        from prepare.cli.base import anchored_folder
        from prepare.data.branchwater import parse_containment_data
        from prepare.data.registry import record_screening_from_table, registry_transaction
        from prepare.data.screening_table import read_screening_table

        try:
            matches_folder = str(anchored_folder(args.matches_folder, args.registry))
            errors: List[str] = []
            summary = parse_containment_data(
                matches_folder,
                args.parsed_containment_file,
                args.summary_containment_file,
                args.step_size,
                details_file=args.details_file,
                errors=errors,
            )
            # The table just built is used directly; a summary without one falls back to the file.
            table = getattr(summary, "table", None)
            if table is None:
                table = read_screening_table(args.parsed_containment_file)
            with registry_transaction(args.registry) as registry:
                if table is not None:
                    record_screening_from_table(
                        registry, table, matches_folder, max_screened=args.registry_max_screened
                    )
        except Exception as e:
            logger.error(f"Error parsing containment: {e}")
            sys.exit(1)
        if errors:
            logger.error("%d match file(s) could not be read: %s", len(errors), ", ".join(errors))
            sys.exit(1)


class PreparePlotPlugin(CommandLinePlugin):
    """Plot containment data from parsed results."""

    command = "prepare_plot"
    description = "Plot containment data from parsed results"

    def __init__(self, parser):
        super().__init__(parser)
        parser.add_argument(
            "--file-path",
            required=True,
            help="Path to the containment file",
        )
        parser.add_argument(
            "--column",
            default="max_containment",
            help="Column to plot",
        )
        parser.add_argument("--title", default=None, help="Title of the plot")
        parser.add_argument("--colors", default=None, help="Colors to use in the plot")
        parser.add_argument(
            "--show-title",
            action="store_true",
            help="Whether to display the title",
        )
        parser.add_argument(
            "--save-format",
            default=None,
            choices=["png", "jpg", "pdf", "svg"],
            help="Format to save the plot",
        )
        parser.add_argument(
            "--threshold",
            type=finite_float,
            default=None,
            help="Minimum value to include in the plot",
        )
        parser.add_argument(
            "--plot-type",
            default="rank",
            choices=["rank", "histogram", "box", "violin"],
            help="Type of plot to generate",
        )

    def main(self, args):
        super().main(args)
        setup_logging(level=logging.INFO)
        from prepare.visualization.plots import plot_containment

        try:
            plot_containment(
                file_path=args.file_path,
                column=args.column,
                title=args.title,
                colors=args.colors,
                show_title=args.show_title,
                save_format=args.save_format,
                threshold=args.threshold,
                plot_type=args.plot_type,
            )
        except Exception as e:
            logger.error(f"Error plotting containment: {e}")
            sys.exit(1)


class PrepareDiversityPlugin(CommandLinePlugin):
    """Calculate alpha and beta diversity metrics."""

    command = "prepare_diversity"
    description = "Calculate alpha and beta diversity metrics from abundance data"

    def __init__(self, parser):
        super().__init__(parser)
        parser.add_argument(
            "--abundance-file",
            required=True,
            help="CSV file with abundance matrix (samples x species)",
        )
        parser.add_argument(
            "--metadata-file",
            help="CSV file with sample metadata",
        )
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
        parser.add_argument(
            "--permanova-formula",
            help="PERMANOVA formula (e.g., 'treatment + site')",
        )

    def main(self, args):
        super().main(args)
        setup_logging(level=logging.INFO)
        import pandas as pd
        from pathlib import Path

        from prepare.data.file_io import open_atomic, write_csv
        from prepare.processing.diversity import (
            calculate_alpha_diversity,
            calculate_beta_diversity,
            check_metadata_covers_samples,
            perform_permanova,
        )

        try:
            abundance_df = pd.read_csv(args.abundance_file, index_col=0)

            metadata_df = None
            if args.metadata_file:
                metadata_df = pd.read_csv(args.metadata_file, index_col=0)
                if args.permanova_formula:
                    # Every sample needs a metadata row; check before any output is written.
                    check_metadata_covers_samples(abundance_df.index, metadata_df)

            output_dir = Path(args.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            # Alpha diversity
            alpha_div = calculate_alpha_diversity(abundance_df, args.alpha_metrics)
            alpha_output = output_dir / "alpha_diversity.csv"
            write_csv(alpha_div, alpha_output)
            logger.info(f"Alpha diversity results saved to {alpha_output}")

            # Beta diversity
            beta_div = calculate_beta_diversity(abundance_df, args.beta_metric)
            assert isinstance(beta_div, pd.DataFrame)
            beta_output = output_dir / f"beta_diversity_{args.beta_metric}.csv"
            write_csv(beta_div, beta_output)
            logger.info(f"Beta diversity results saved to {beta_output}")

            # PERMANOVA
            if args.permanova_formula and metadata_df is not None:
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

        except Exception as e:
            logger.error(f"Diversity analysis failed: {e}")
            sys.exit(1)


class PrepareTaxonomyPlugin(CommandLinePlugin):
    """Validate species names against NCBI taxonomy."""

    command = "prepare_taxonomy"
    description = "Validate species names against NCBI taxonomy"

    def __init__(self, parser):
        super().__init__(parser)
        parser.add_argument(
            "--species-file",
            required=True,
            help="Text file with species names or CSV with species column",
        )
        parser.add_argument(
            "--species-column",
            help="Column name containing species (for CSV files)",
        )
        parser.add_argument(
            "--email",
            help="Email address for NCBI API access (default: PREPARE_NCBI_EMAIL or [runtime] ncbi_email)",
        )
        parser.add_argument(
            "--api-key",
            help="NCBI API key for increased rate limits (default: PREPARE_NCBI_API_KEY or NCBI_API_KEY)",
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

    def main(self, args):
        super().main(args)
        setup_logging(level=logging.INFO)
        import pandas as pd

        from prepare.core.settings import require_email, setting_for
        from prepare.data.taxonomy import validate_taxonomic_assignments

        try:
            # Load species list
            if args.species_file.endswith(".csv"):
                df = pd.read_csv(args.species_file)
                if args.species_column:
                    if args.species_column not in df.columns:
                        logger.error(f"Column '{args.species_column}' not found")
                        sys.exit(1)
                    species_list = df[args.species_column].dropna().tolist()
                else:
                    species_list = df.iloc[:, 0].dropna().tolist()
            else:
                with open(args.species_file, "r") as f:
                    species_list = [line.strip() for line in f if line.strip()]

            logger.info(f"Loaded {len(species_list)} species names")

            results_df = validate_taxonomic_assignments(
                species_list,
                email=require_email(args),
                api_key=setting_for(args, "ncbi_api_key"),
                output_file=args.output_file,
                cache_file=args.cache_file,
            )

            valid_count = results_df["is_valid"].sum()
            total_count = len(results_df)
            logger.info(f"Valid: {valid_count}/{total_count} ({valid_count / total_count * 100:.1f}%)")

        except Exception as e:
            logger.error(f"Taxonomy validation failed: {e}")
            sys.exit(1)
