"""Taxonomy enrichment and containment exploration CLI commands."""

import argparse
from pathlib import Path

import pandas as pd

from prepare.cli.base import finite_float, BaseCommand
from prepare.core.exceptions import ExitCode, PrepareError
from prepare.core.utils import get_genome_columns
from prepare.data.file_io import write_csv
from prepare.data.genome_taxonomy import (
    annotate_containment_with_taxonomy,
    enrich_genomes_with_taxonomy,
    filter_by_taxonomy,
    load_taxonomy_cache,
    save_taxonomy_cache,
    summarize_by_taxonomy,
)


class EnrichTaxonomyCommand(BaseCommand):
    """Command to enrich genome accessions with GTDB taxonomy."""

    @property
    def name(self) -> str:
        return "enrich_taxonomy"

    @property
    def help(self) -> str:
        return "Enrich genome accessions with taxonomy from GTDB"

    @property
    def group(self) -> str:
        return "Containment"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--parsed-containment",
            required=True,
            help="Path to parsed containment file",
        )
        parser.add_argument(
            "--output",
            default="taxonomy_map.tsv",
            help="Output taxonomy TSV file",
        )
        parser.add_argument(
            "--cache",
            default="taxonomy_cache.tsv",
            help="Cache file path for GTDB lookups",
        )

    def execute(self, args: argparse.Namespace) -> int:
        try:
            containment_path = Path(args.parsed_containment)
            if not containment_path.exists():
                self.logger.error("Containment file not found: %s", containment_path)
                return 1

            # Only the column names are needed here, so the data rows are not read.
            df = pd.read_csv(containment_path, sep="\t", index_col=0, nrows=0)
            genome_cols = get_genome_columns(df)
            self.logger.info("Found %d genome columns in %s", len(genome_cols), containment_path)

            cache_file = Path(args.cache)
            taxonomy = enrich_genomes_with_taxonomy(genome_cols, cache_file=cache_file)

            output_path = Path(args.output)
            save_taxonomy_cache(taxonomy, output_path)
            self.logger.info(
                "Wrote taxonomy map with %d entries to %s",
                len(taxonomy),
                output_path,
            )
            return 0
        except PrepareError as e:
            return self.fail(e, "Error enriching taxonomy")


class ExploreContainmentCommand(BaseCommand):
    """Command to generate an interactive HTML containment explorer."""

    @property
    def name(self) -> str:
        return "explore_containment"

    @property
    def help(self) -> str:
        return "Generate interactive HTML explorer for containment data"

    @property
    def group(self) -> str:
        return "Containment"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--parsed-containment",
            required=True,
            help="Path to parsed containment file",
        )
        parser.add_argument(
            "--taxonomy-map",
            default=None,
            help="Pre-computed taxonomy TSV (enriches on the fly if not given)",
        )
        parser.add_argument(
            "--metadata",
            default=None,
            help="Optional metadata table",
        )
        parser.add_argument(
            "--output",
            default="containment_explorer.html",
            help="Output HTML file",
        )
        parser.add_argument(
            "--min-containment",
            type=finite_float,
            default=0.0,
            help="Minimum containment threshold for inclusion",
        )
        parser.add_argument(
            "--cache",
            default="taxonomy_cache.tsv",
            help="Taxonomy cache file for GTDB lookups",
        )
        parser.add_argument(
            "--open",
            dest="open_browser",
            action="store_true",
            help="Open the generated explorer in a browser",
        )

    def execute(self, args: argparse.Namespace) -> int:
        try:
            from prepare.visualization.explorer import generate_containment_explorer, require_explorer_packages

            # Check before the taxonomy lookup, which may call NCBI.
            require_explorer_packages()

            containment_path = Path(args.parsed_containment)
            if not containment_path.exists():
                self.logger.error("Containment file not found: %s", containment_path)
                return 1

            df = pd.read_csv(containment_path, sep="\t", index_col=0)

            # Load or compute taxonomy
            if args.taxonomy_map:
                taxonomy_path = Path(args.taxonomy_map)
                if not taxonomy_path.exists():
                    self.logger.error("Taxonomy map not found: %s", taxonomy_path)
                    return 1
                taxonomy = load_taxonomy_cache(taxonomy_path)
            else:
                genome_cols = get_genome_columns(df)
                cache_file = Path(args.cache)
                taxonomy = enrich_genomes_with_taxonomy(genome_cols, cache_file=cache_file)

            # Load optional metadata
            metadata = None
            if args.metadata:
                metadata_path = Path(args.metadata)
                if not metadata_path.exists():
                    self.logger.error("Metadata file not found: %s", metadata_path)
                    return 1
                metadata = pd.read_csv(metadata_path, sep="\t", index_col=0)

            output_path = generate_containment_explorer(
                containment_df=df,
                taxonomy=taxonomy,
                metadata=metadata,
                output_file=args.output,
                min_containment=args.min_containment,
            )
            self.logger.info("Generated containment explorer: %s", output_path)

            if getattr(args, "open_browser", False):
                from prepare.utils.browser import open_in_browser

                if not open_in_browser(output_path):
                    self.logger.warning("Could not open a browser automatically for %s", output_path)
            else:
                self.logger.info("Next: open %s in a browser", output_path)

            return 0
        except PrepareError as e:
            return self.fail(e, "Error generating explorer")


class FindByTaxonomyCommand(BaseCommand):
    """Command to filter containment results by taxonomy."""

    @property
    def name(self) -> str:
        return "find_by_taxonomy"

    @property
    def help(self) -> str:
        return "Filter containment results by taxonomic classification"

    @property
    def group(self) -> str:
        return "Containment"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--parsed-containment",
            required=True,
            help="Path to parsed containment file",
        )
        parser.add_argument(
            "--taxonomy-map",
            required=True,
            help="Pre-computed taxonomy TSV",
        )
        parser.add_argument(
            "--family",
            default=None,
            help="Filter by family name",
        )
        parser.add_argument(
            "--genus",
            default=None,
            help="Filter by genus name",
        )
        parser.add_argument(
            "--species",
            default=None,
            help="Filter by species name",
        )
        parser.add_argument(
            "--min-containment",
            type=finite_float,
            default=0.0,
            help="Minimum containment threshold",
        )
        parser.add_argument(
            "--output",
            default=None,
            help="Output TSV file (default: stdout)",
        )
        parser.add_argument(
            "--format",
            dest="output_format",
            choices=["table", "summary"],
            default="table",
            help="Output format: detailed table or taxonomy summary",
        )

    def execute(self, args: argparse.Namespace) -> int:
        try:
            if not any([args.family, args.genus, args.species]):
                self.logger.error("At least one of --family, --genus, or --species is required")
                return int(ExitCode.USAGE)

            containment_path = Path(args.parsed_containment)
            if not containment_path.exists():
                self.logger.error("Containment file not found: %s", containment_path)
                return 1

            taxonomy_path = Path(args.taxonomy_map)
            if not taxonomy_path.exists():
                self.logger.error("Taxonomy map not found: %s", taxonomy_path)
                return 1

            df = pd.read_csv(containment_path, sep="\t", index_col=0)
            taxonomy = load_taxonomy_cache(taxonomy_path)

            annotated = annotate_containment_with_taxonomy(df, taxonomy)
            filtered = filter_by_taxonomy(
                annotated,
                family=args.family,
                genus=args.genus,
                species=args.species,
                min_containment=args.min_containment,
            )

            if filtered.empty:
                self.logger.warning("No matching results found")
                return 0

            if args.output_format == "summary":
                # Determine the most specific taxonomy level requested
                level = "family"
                if args.species:
                    level = "species"
                elif args.genus:
                    level = "genus"
                result = summarize_by_taxonomy(filtered, level=level)
                keep_index = True  # the sample names are the index of the summary table
            else:
                result = filtered
                keep_index = False  # the table has a sample column; its index is a row label

            if args.output:
                write_csv(result, args.output, sep="\t", index=keep_index)
                self.logger.info("Wrote %d results to %s", len(result), args.output)
            else:
                self.emit_raw(result.to_csv(sep="\t", index=keep_index))

            return 0
        except PrepareError as e:
            return self.fail(e, "Error filtering by taxonomy")
