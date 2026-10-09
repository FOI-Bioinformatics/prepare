"""
Command Line Interface entry point for PREPARE.

This module provides the main CLI entry point with modular command architecture.
"""

import argparse
import logging
import os
import socket
import sys
import time
from datetime import datetime, timezone
from typing import Dict, List, Optional

from prepare import __version__
from prepare.cli.base import BaseCommand, DefaultsHelpFormatter, add_global_options, command_registry
from prepare.core import settings
from prepare.core.exceptions import ConfigurationError, DataAccessError, PrepareError, exit_code_for
from prepare.data import run_log
from prepare.utils.logging import log_traceback_hint, setup_logging

# Import all command modules to register them
from prepare.cli.commands import (
    UseBranchwaterCommand,
    ExtractBranchwaterMetadataCommand,
    ParseContainmentCommand,
    PlotContainmentCommand,
    DownloadMetadataCommand,
    ParseMetadataCommand,
    CheckMetadataAttributesCommand,
    CountMetadataCommand,
    PlotMetadataCountsCommand,
    SingleSampleCommand,
    DownloadSraCommand,
    StatusCommand,
    ResultsTableCommand,
    ExtractTargetReadsCommand,
    DownloadTestGenomeCommand,
    GenomeSearchCommand,
    GenomeDownloadCommand,
    GenomePrepareCommand,
    EnrichTaxonomyCommand,
    ExploreContainmentCommand,
    FindByTaxonomyCommand,
)
from prepare.cli.commands.select import SelectDatasetsCommand
from prepare.cli.commands.blacklist import BlacklistCommand
from prepare.cli.commands.doctor import STARTUP_ERROR_ATTR, DoctorCommand
from prepare.cli.commands.branchwater_search import BranchwaterSearchCommand
from prepare.cli.commands.store import (
    StoreAdoptCommand,
    StoreGcCommand,
    StoreInitCommand,
    StoreLinkCommand,
    StoreReindexCommand,
    StoreStatusCommand,
    StoreUnlinkCommand,
    StoreUsageCommand,
    StoreVerifyCommand,
)
from prepare.cli.commands.advanced_analysis import (
    DiversityAnalysisCommand,
    InteractivePlotCommand,
    TaxonomyValidationCommand,
    TaxonomicSummaryCommand,
)
from prepare.cli.commands.project_report import ProjectReportCommand
from prepare.cli.commands.renamed import renamed_commands
from prepare.cli.commands.runs import RunsCommand
from prepare.cli.commands.sra_enhanced import SRAInfoCommand, SRAValidateCommand
from prepare.cli.commands.sra_profile import SRAProfileCommand
from prepare.cli.commands.sra_report import SRAReportCommand


def register_all_commands() -> None:
    """Register all available commands with the registry."""
    commands = [
        # Containment commands
        BranchwaterSearchCommand(),
        UseBranchwaterCommand(),
        ParseContainmentCommand(),
        PlotContainmentCommand(),
        ExploreContainmentCommand(),
        EnrichTaxonomyCommand(),
        FindByTaxonomyCommand(),
        # Metadata commands
        ExtractBranchwaterMetadataCommand(),
        DownloadMetadataCommand(),
        ParseMetadataCommand(),
        CheckMetadataAttributesCommand(),
        CountMetadataCommand(),
        SingleSampleCommand(),
        PlotMetadataCountsCommand(),
        # Genome commands
        GenomeSearchCommand(),
        GenomeDownloadCommand(),
        GenomePrepareCommand(),
        DownloadTestGenomeCommand(),
        # Reads commands
        SelectDatasetsCommand(),
        BlacklistCommand(),
        DownloadSraCommand(),
        StatusCommand(),
        ResultsTableCommand(),
        SRAInfoCommand(),
        SRAValidateCommand(),
        SRAProfileCommand(),
        SRAReportCommand(),
        ExtractTargetReadsCommand(),
        # Store commands
        StoreInitCommand(),
        StoreStatusCommand(),
        StoreReindexCommand(),
        StoreAdoptCommand(),
        StoreVerifyCommand(),
        StoreLinkCommand(),
        StoreUnlinkCommand(),
        StoreUsageCommand(),
        StoreGcCommand(),
        # Analysis commands
        DiversityAnalysisCommand(),
        InteractivePlotCommand(),
        TaxonomyValidationCommand(),
        TaxonomicSummaryCommand(),
        # Environment commands
        DoctorCommand(),
        RunsCommand(),
        ProjectReportCommand(),
        # Former names (0.5.0), hidden from the help
        *renamed_commands(),
    ]

    for command in commands:
        command_registry.register(command)


GROUP_ORDER = ["Containment", "Metadata", "Genomes", "Reads", "Store", "Analysis", "Environment", "Other"]


class _HelpFormatter(DefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Show option defaults (when they exist) and keep the epilog's line breaks."""


def _commands_epilog(commands: Dict[str, BaseCommand]) -> str:
    """List commands under their pipeline step for the main --help; hidden commands are left out."""
    by_group: Dict[str, List[BaseCommand]] = {}
    for command in commands.values():
        if command.hidden:
            continue
        by_group.setdefault(command.group, []).append(command)
    lines = ["commands by pipeline step:"]
    for group in GROUP_ORDER + sorted(set(by_group) - set(GROUP_ORDER)):
        if group not in by_group:
            continue
        lines.append(f"  {group}:")
        for command in by_group[group]:
            lines.append(f"    {command.name:<28} {command.help}")
    return "\n".join(lines)


def create_parser() -> argparse.ArgumentParser:
    """
    Create the command line argument parser.

    Returns:
        Configured ArgumentParser instance
    """
    # Register all commands before building the parser so the epilog can list them.
    register_all_commands()

    parser = argparse.ArgumentParser(
        description="PREPARE: A toolkit for analyzing metagenomic datasets based on genome containment.",
        epilog=_commands_epilog(command_registry.get_all_commands()),
        formatter_class=_HelpFormatter,
    )

    parser.add_argument("--version", action="version", version=f"PREPARE v{__version__}")

    add_global_options(parser, suppress_defaults=False)

    command_registry.setup_parsers(parser)

    return parser


def _flag(args: argparse.Namespace, name: str) -> bool:
    """True only when ``args`` really holds ``name`` set to True (a mock namespace invents attributes)."""
    try:
        return vars(args).get(name) is True
    except TypeError:
        return False


def _apply_quiet_and_verbose(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Turn ``-q``/``-v`` into ``args.log_level``; together (in either placement) they are a usage error."""
    quiet, verbose = _flag(args, "log_quiet"), _flag(args, "log_verbose")
    if quiet and verbose:
        parser.error("argument -q/--quiet: not allowed with argument -v/--verbose")
    if quiet:
        args.log_level = "WARNING"
    elif verbose:
        args.log_level = "DEBUG"


# Flags whose value must not reach a log file shared with other users (shared with the run log).
_masked_argv = settings.masked_argv


def _log_run_header(argv: List[str], parsed: Optional[argparse.Namespace] = None) -> None:
    """Log, at DEBUG, what identifies this run in a shared log: version, arguments, host, PID, SLURM IDs."""
    logging.debug("PREPARE v%s: %s", __version__, " ".join(["prepare", *_masked_argv(argv, parsed)]))
    logging.debug("Host %s, process ID %d", socket.gethostname(), os.getpid())
    slurm = [f"{name}={os.environ[name]}" for name in ("SLURM_JOB_ID", "SLURM_ARRAY_TASK_ID") if os.environ.get(name)]
    if slurm:
        logging.debug("SLURM: %s", ", ".join(slurm))


def _configure_logging(runtime: settings.RuntimeSettings) -> Optional[int]:
    """Set up logging from the runtime settings; an exit code when the log file cannot be opened."""
    level = getattr(logging, runtime.log_level)
    try:
        setup_logging(
            level=level,
            log_file=runtime.log_file,
            show_host=runtime.log_host,
            console_traceback=level <= logging.DEBUG,
        )
    except OSError as e:
        setup_logging(level=level, show_host=runtime.log_host)
        error = ConfigurationError(
            f"Cannot open the log file {runtime.log_file} (from {runtime.source('log_file')}): {e}"
        )
        logging.error(f"Error: {error}")
        return exit_code_for(error)
    return None


def _log_failure(message: str, error: BaseException) -> None:
    """Log a failure that left the command: one line on the console, the traceback in the log file."""
    logging.error(message, exc_info=error)
    log_traceback_hint()


def _run_command(parsed_args: argparse.Namespace) -> int:
    """Run the parsed command and map an exception that leaves it to its exit code."""
    try:
        return parsed_args.func(parsed_args)
    except PrepareError as e:
        _log_failure(f"Error: {e}", e)
        return exit_code_for(e)
    except Exception as e:
        _log_failure(f"{type(e).__name__}: {e}", e)
        return exit_code_for(e)
    except KeyboardInterrupt as e:
        # Only a command with graceful_shutdown False gets here; BaseCommand.run handles the rest.
        logging.error("Interrupted")
        return exit_code_for(e)


def _record_run_if_wanted(
    parsed_args: argparse.Namespace, argv: List[str], started: datetime, clock: float, exit_code: int
) -> None:
    """Add this run to the project's run log when its command keeps one; never changes the exit code."""
    command = getattr(getattr(parsed_args, "func", None), "__self__", None)
    try:
        if not isinstance(command, BaseCommand) or not command.records_run(parsed_args):
            return
        project = run_log.project_for(getattr(parsed_args, "registry", None))
        if project is not None:
            seconds = time.monotonic() - clock
            run_log.record_run(project, command.name, argv, parsed_args, started, seconds, exit_code)
    except (OSError, DataAccessError, ValueError, TypeError, AttributeError, KeyError) as e:
        # AttributeError and KeyError: a records_run override reading an argument the parsed
        # namespace lacks must not change the exit code either.
        logging.warning("This run was not recorded in the run log: %s", e)


def main(args: Optional[List[str]] = None) -> int:
    """
    Main entry point for the CLI.

    Args:
        args: Command line arguments (if None, sys.argv[1:] is used)

    Returns:
        Exit code: 0 success, 1 failure, 2 usage error, 3 configuration problem, 4 retryable
        failure (network, a lock wait that gave up), 130 interrupted
        (``prepare.core.exceptions.ExitCode``).
    """
    argv = list(sys.argv[1:] if args is None else args)
    parser = create_parser()
    parsed_args = parser.parse_args(args)
    _apply_quiet_and_verbose(parser, parsed_args)

    # Every runtime setting is resolved once, before logging is set up, since the level is one.
    try:
        runtime = settings.activate(parsed_args)
    except ConfigurationError as e:
        if getattr(parsed_args, "command", None) != "doctor":
            setup_logging(level=logging.INFO)
            logging.error(f"Error: {e}")
            return exit_code_for(e)
        # doctor exists to diagnose exactly this: it runs on the built-in defaults and its
        # config check reports the error (exit code 3 through its own failed check).
        runtime = settings.activate_defaults()
        setattr(parsed_args, STARTUP_ERROR_ATTR, str(e))
    failed = _configure_logging(runtime)
    if failed is not None:
        return failed
    for message in runtime.warnings:
        logging.getLogger("prepare.core.settings").warning(message)
    _log_run_header(argv, parsed_args)
    logging.debug("Runtime settings (value and source):\n  %s", "\n  ".join(runtime.describe()))

    started, clock = datetime.now(timezone.utc), time.monotonic()
    exit_code = 1
    try:
        exit_code = _run_command(parsed_args)
    except SystemExit as e:
        # Recorded with the code the process exits with: an int as given, None as 0, else 1.
        exit_code = 0 if e.code is None else e.code if isinstance(e.code, int) else 1
        raise
    finally:
        _record_run_if_wanted(parsed_args, argv, started, clock, exit_code)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
