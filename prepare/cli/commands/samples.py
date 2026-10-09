"""
Sample analysis CLI commands.
"""

import argparse

from prepare.cli.base import BaseCommand
from prepare.core.exceptions import PrepareError
from prepare.data.defaults import resolve_metadata_table
from prepare.processing.containment import count_single_sample


class SingleSampleCommand(BaseCommand):
    """Command for analyzing a single sample."""

    @property
    def name(self) -> str:
        return "single_sample"

    @property
    def help(self) -> str:
        return "Analyze a single sample"

    @property
    def group(self) -> str:
        return "Metadata"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--summary-file",
            default="parsed_containment.txt",
            help="Path to the summary file",
        )
        parser.add_argument(
            "--metadata-file",
            default=None,
            help="Parsed metadata table (default: metadata_table.txt, else metadata/branchwater_metadata.txt)",
        )
        parser.add_argument(
            "--summary-column",
            required=True,
            help="Name of the column in the summary file",
        )
        parser.add_argument(
            "--metadata-column",
            required=True,
            help="Name of the column in the metadata file",
        )
        parser.add_argument(
            "--threshold",
            type=float,
            default=0.1,
            help="Threshold for the column in the summary file (inclusive)",
        )
        parser.add_argument("--top-n", type=int, default=100, help="Number of top items to keep")

    def execute(self, args: argparse.Namespace) -> int:
        try:
            count_dict = count_single_sample(
                summary_file=args.summary_file,
                metadata_file=str(resolve_metadata_table(args.metadata_file)),
                summary_column=args.summary_column,
                metadata_column=args.metadata_column,
                threshold=args.threshold,
                top_n=args.top_n,
            )

            if not count_dict:
                self.logger.warning("No data found above threshold")
                return 0

            # Log the top N items
            for key, value in list(count_dict.items())[: args.top_n]:
                self.logger.info(f"{key}: {value}")

            return 0
        except PrepareError as e:
            return self.fail(e, "Error analyzing single sample")
