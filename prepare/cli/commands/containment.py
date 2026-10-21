"""
Containment-related CLI commands.
"""

import argparse
import shlex


from prepare.cli.base import finite_float, BaseCommand, anchored_folder
from prepare.core.exceptions import PrepareError
from prepare.data.branchwater import parse_containment_data
from prepare.core.constants import DEFAULT_REGISTRY_MAX_SCREENED
from prepare.data.registry import record_screening_from_table, registry_transaction
from prepare.data.screening_table import read_screening_table
from prepare.visualization.plots import plot_containment as viz_plot_containment, plot_output_path


class ParseContainmentCommand(BaseCommand):
    """Command for parsing containment data from match files."""

    @property
    def name(self) -> str:
        return "parse_containment"

    @property
    def help(self) -> str:
        return "Parse containment data from match files"

    @property
    def group(self) -> str:
        return "Containment"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--matches-folder",
            default="matches",
            help="Folder containing match files",
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

    def execute(self, args: argparse.Namespace) -> int:
        try:
            args.matches_folder = str(anchored_folder(args.matches_folder, getattr(args, "registry", None)))
            summary = parse_containment_data(
                args.matches_folder,
                args.parsed_containment_file,
                args.summary_containment_file,
                args.step_size,
                details_file=args.details_file,
            )
            # The table just built, so it is not read back from disk. A summary without one (no
            # match data) falls back to the file, read here before the registry lock is taken.
            table = getattr(summary, "table", None)
            if table is None:
                table = read_screening_table(args.parsed_containment_file)
            with registry_transaction(args.registry) as registry:
                if table is not None:
                    record_screening_from_table(
                        registry, table, args.matches_folder, max_screened=args.registry_max_screened
                    )
            hint = f"Next: prepare plot_containment --file-path {shlex.quote(str(args.parsed_containment_file))}"
            genome_id = next(iter(sorted(getattr(summary, "genome_to_samples", None) or {})), None)
            if genome_id:
                hint += (
                    "; then prepare select_datasets "
                    f"--parsed-containment {shlex.quote(str(args.parsed_containment_file))} "
                    f"--genome-id {shlex.quote(str(genome_id))} --threshold 0.5"
                )
            self.logger.info(hint)
            return 0
        except PrepareError as e:
            return self.fail(e, "Error parsing containment")


class PlotContainmentCommand(BaseCommand):
    """Command for plotting containment data."""

    @property
    def name(self) -> str:
        return "plot_containment"

    @property
    def help(self) -> str:
        return "Plot containment data"

    @property
    def group(self) -> str:
        return "Containment"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument("--file-path", required=True, help="Path to the containment file")
        parser.add_argument("--column", default="max_containment", help="Column to plot")
        parser.add_argument("--title", default=None, help="Title of the plot")
        parser.add_argument("--colors", default=None, help="Colors to use in the plot")
        parser.add_argument("--show-title", action="store_true", help="Whether to display the title")
        parser.add_argument(
            "--save-format",
            default="png",
            choices=["png", "jpg", "pdf", "svg"],
            help="Image format of the saved plot (default: png)",
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

    def execute(self, args: argparse.Namespace) -> int:
        try:
            fig = viz_plot_containment(
                file_path=args.file_path,
                column=args.column,
                title=args.title,
                colors=args.colors,
                show_title=args.show_title or args.title is not None,
                save_format=args.save_format,
                threshold=args.threshold,
                plot_type=args.plot_type,
            )
            if fig is not None:
                output_path = plot_output_path(args.file_path, args.plot_type, args.column, args.save_format)
                self.logger.info("Plot saved to %s", output_path)
                self.logger.info("Next: open %s to view the plot", output_path)
            return 0
        except PrepareError as e:
            return self.fail(e, "Error plotting containment")
