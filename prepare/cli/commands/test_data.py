"""
Test data CLI commands.
"""

import argparse

from prepare.cli.base import BaseCommand, anchored_folder
from prepare.core.exceptions import PrepareError
from prepare.processing.containment import download_test_genome


class DownloadTestGenomeCommand(BaseCommand):
    """Command for downloading test genome."""

    @property
    def name(self) -> str:
        return "download_test_genome"

    @property
    def help(self) -> str:
        return "Download test genome"

    @property
    def group(self) -> str:
        return "Genomes"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--output-folder",
            default="genomes",
            help="Folder to save the downloaded fasta files",
        )

    def execute(self, args: argparse.Namespace) -> int:
        try:
            args.output_folder = str(anchored_folder(args.output_folder, getattr(args, "registry", None)))
            download_test_genome(args.output_folder)
            return 0
        except PrepareError as e:
            return self.fail(e, "Error downloading test genome")
