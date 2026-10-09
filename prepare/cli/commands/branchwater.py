"""
Branchwater-related CLI commands.
"""

import argparse
from typing import List

from prepare.cli.base import BaseCommand
from pathlib import Path

from prepare.core.exceptions import PrepareError
from prepare.data.branchwater import (
    process_branchwater_files,
    extract_metadata_from_branchwater,
)


class UseBranchwaterCommand(BaseCommand):
    """Command for processing pre-downloaded Branchwater files."""

    @property
    def name(self) -> str:
        return "use_branchwater"

    @property
    def help(self) -> str:
        return "Process pre-downloaded Branchwater files"

    @property
    def group(self) -> str:
        return "Containment"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--branchwater-folder",
            required=True,
            help="Folder containing Branchwater output files",
        )
        parser.add_argument(
            "--matches-folder",
            default="matches",
            help="Folder to save processed matches",
        )

    def execute(self, args: argparse.Namespace) -> int:
        try:
            errors: List[str] = []
            process_branchwater_files(args.branchwater_folder, args.matches_folder, errors=errors)
            if errors:
                self.logger.error("%d match file(s) could not be read: %s", len(errors), ", ".join(errors))
                return 1
            self.logger.info("Next: prepare parse_containment --matches-folder %s", args.matches_folder)
            return 0
        except PrepareError as e:
            return self.fail(e, "Error processing Branchwater files")


class ExtractBranchwaterMetadataCommand(BaseCommand):
    """Command for extracting metadata from Branchwater files."""

    @property
    def name(self) -> str:
        return "extract_branchwater_metadata"

    @property
    def help(self) -> str:
        return "Extract metadata from Branchwater files"

    @property
    def group(self) -> str:
        return "Metadata"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--branchwater-folder",
            required=True,
            help="Folder containing Branchwater output files",
        )
        parser.add_argument(
            "--metadata-folder",
            default="metadata",
            help="Folder to save extracted metadata",
        )

    def execute(self, args: argparse.Namespace) -> int:
        try:
            metadata_folder = Path(args.metadata_folder)
            metadata_folder.mkdir(exist_ok=True)
            output_file = metadata_folder / "branchwater_metadata.txt"

            errors: List[str] = []
            extract_metadata_from_branchwater(args.branchwater_folder, output_file, errors=errors)
            if errors:
                self.logger.error("%d match file(s) could not be read: %s", len(errors), ", ".join(errors))
                return 1
            return 0
        except PrepareError as e:
            return self.fail(e, "Error extracting metadata")
