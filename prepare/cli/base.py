"""
Base classes for CLI command system.

This module provides the foundation for a modular command architecture.
"""

import argparse
import json
import math
import logging
import re
import sys
from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Sequence, Union

from prepare.core.constants import DEFAULT_LOG_LEVEL, LOG_LEVELS, SRA_ACCESSION_PATTERN
from prepare.core.exceptions import exit_code_for
from prepare.data.accessions import read_accessions_file as _read_accessions
from prepare.data.accessions import require_accession
from prepare.utils.logging import log_traceback_hint
from prepare.utils.security import SecureSubprocess
from prepare.utils.termination import EXIT_INTERRUPTED, graceful_termination

if TYPE_CHECKING:
    from prepare.data.registry import Registry
    from prepare.store.layout import StorePaths


def read_accessions_file(path: Union[str, Path]) -> List[str]:
    """Accessions listed in ``path``, one per line, in file order.

    Surrounding whitespace is stripped; blank lines and lines starting with ``#`` are skipped.
    A file that does not exist or cannot be read raises ``DataAccessError`` naming it.
    """
    return _read_accessions(path)


def refuse_bad_accessions(logger: logging.Logger, accessions: Sequence[str]) -> Optional[int]:
    """Exit code 2 after one error line when any accession is malformed; ``None`` otherwise.

    Commands call this before any filesystem or registry work, so a value such as ``../x`` never becomes a path.
    """
    bad = [acc for acc in accessions if not re.fullmatch(SRA_ACCESSION_PATTERN, acc)]
    if bad:
        logger.error("not an SRA accession: %s", ", ".join(bad))
        return 2
    return None


def accessions_from_args(accessions_file: Optional[str], accessions: Optional[Sequence[str]]) -> List[str]:
    """The accessions named by ``--accessions-file`` and repeated ``--accession`` flags together.

    File entries come first, then the flag values, each accession once in first-seen order.
    Empty when neither was given. A flag value that is not a well-formed SRA accession raises
    ``ValidationError`` (exit code 2), as a bad line in the file does, before it is used as a path.
    """
    named = read_accessions_file(accessions_file) if accessions_file else []
    named += [require_accession(accession) for accession in accessions or []]
    return list(dict.fromkeys(named))


def anchored_folder(folder: Union[str, Path], registry: Optional[str]) -> Path:
    """A project folder (``--fastq-folder`` and similar) resolved against the project root.

    An absolute ``folder`` is returned as is. A relative one is joined onto the folder of the
    registry the command uses (``registry``, or the nearest ``prepare_registry.json`` walking up
    from the working directory), so a command run from a subfolder reads and writes the same
    folders that the registry records. With no registry found above, the root is the working
    directory, as before.
    """
    from prepare.data.registry import project_root, registry_path

    path = Path(folder)
    return path if path.is_absolute() else project_root(registry_path(registry)) / path


def resolve_command_store(args: argparse.Namespace, registry: "Registry") -> Optional["StorePaths"]:
    """The shared data store (if any) for a command's ``--data-root`` and registry.

    ``getattr`` guards ``args.data_root`` so a namespace built without that attribute is not
    broken by it. A store that cannot be reached only costs the usage record, so it is
    logged and skipped rather than failing an analysis the project can run on its own files.
    """
    from prepare.data import registry_blocks as rb
    from prepare.store.resolve import resolve_optional_store

    return resolve_optional_store(getattr(args, "data_root", None), rb.store_block(registry).root)


class DefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Append '(default: ...)' only when the option has a real default.

    A flag (``store_true``, ``store_false``) gets no note: "(default: False)" says nothing a
    reader does not already know, and a flag whose default matters (``--verify-downloads``)
    says so in its own help text.
    """

    def _get_help_string(self, action: argparse.Action) -> str:
        if action.default is None or action.default is argparse.SUPPRESS:
            return action.help or ""
        if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
            return action.help or ""
        return super()._get_help_string(action) or ""


# The level names argparse accepts for --log-level; prepare.core.constants.LOG_LEVELS is the
# one place this list (and its default) is spelled out.
LOG_LEVEL_NAMES = tuple(LOG_LEVELS)
GLOBAL_OPTIONS_TITLE = "logging and progress"


def _non_negative_count(text: str) -> int:
    """An argparse type: a whole number, 0 or more."""
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a whole number, got {text!r}") from None
    if value < 0:
        raise argparse.ArgumentTypeError(f"expected a whole number, 0 or more, got {value}")
    return value


def finite_float(text: str) -> float:
    """An argparse type: a finite number (``nan`` and ``inf`` are refused)."""
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a number, got {text!r}") from None
    if not math.isfinite(value):
        raise argparse.ArgumentTypeError(f"expected a finite number, got {text!r}")
    return value


def add_global_options(parser: argparse.ArgumentParser, suppress_defaults: bool) -> None:
    """Add the logging and progress options to ``parser``.

    They go on the main parser with real defaults (None or False, so the runtime settings
    decide) and on every command's parser with ``argparse.SUPPRESS`` defaults, so a flag given
    after the command name overrides the main parser's value and one left out does not reset
    it. ``-q/--quiet`` and ``-v/--verbose`` are stored as ``log_quiet`` and ``log_verbose``;
    ``main`` turns them into a log level. A command that already has one of these option
    strings keeps its own, and only the clashing flag is left out of its parser, so ``-v``
    always sets the log level (``store_status`` keeps a hidden ``--verbose`` alias for
    ``--details``).
    """

    def default(value: Any) -> Any:
        return argparse.SUPPRESS if suppress_defaults else value

    taken = set(getattr(parser, "_option_string_actions", {}))
    group = parser.add_argument_group(GLOBAL_OPTIONS_TITLE)
    if "--log-level" not in taken:
        group.add_argument(
            "--log-level",
            type=str.upper,
            choices=LOG_LEVEL_NAMES,
            metavar="LEVEL",
            default=default(None),
            help=(
                f"Console logging level (one of {', '.join(LOG_LEVEL_NAMES)}; "
                f"default: PREPARE_LOG_LEVEL, config [runtime] log_level, or {DEFAULT_LOG_LEVEL})"
            ),
        )
    if "--log-file" not in taken:
        group.add_argument(
            "--log-file",
            metavar="PATH",
            default=default(None),
            help=(
                "Append every log line at INFO or above, with host, process ID and tracebacks, to this file "
                "(default: PREPARE_LOG_FILE or config [runtime] log_file; none)"
            ),
        )
    levels = group.add_mutually_exclusive_group()
    if not taken & {"-q", "--quiet"}:
        levels.add_argument(
            "-q",
            "--quiet",
            dest="log_quiet",
            action="store_true",
            default=default(False),
            help="Show only warnings and errors on the console (same as --log-level WARNING)",
        )
    verbose_flags = [flag for flag in ("-v", "--verbose") if flag not in taken]
    if verbose_flags:
        levels.add_argument(
            *verbose_flags,
            dest="log_verbose",
            action="store_true",
            default=default(False),
            help="Show debug lines and tracebacks on the console (same as --log-level DEBUG)",
        )
    if "--progress-every" not in taken:
        group.add_argument(
            "--progress-every",
            type=_non_negative_count,
            metavar="N",
            default=default(None),
            help=(
                "Log a progress summary every N items of a long run; 0 logs one line per item instead "
                "(default: PREPARE_PROGRESS_EVERY, config [runtime] progress_every, or 50)"
            ),
        )


def emit_error_json(message: str, indent: Optional[int] = 2) -> None:
    """Write ``{"error": message}`` to stdout as one JSON document.

    For module-level helpers that have no command instance at hand. The message should also
    be logged by the caller when a human reader needs it; the log goes to stderr.
    """
    print(json.dumps({"error": message}, indent=indent), file=sys.stdout)


class BaseCommand(ABC):
    """Base class for all CLI commands."""

    # False for a command that must keep the default signal behaviour; ``run`` then calls
    # ``execute`` directly.
    graceful_shutdown: bool = True

    def __init__(self):
        """Initialize the command with a logger."""
        self.logger = logging.getLogger(self.__class__.__module__)

    @property
    @abstractmethod
    def name(self) -> str:
        """Return the command name."""
        pass

    @property
    @abstractmethod
    def help(self) -> str:
        """Return the help text for this command."""
        pass

    @property
    def aliases(self) -> List[str]:
        """Return alternate names this command also responds to (default: none)."""
        return []

    @property
    def group(self) -> str:
        """Pipeline step the command belongs to, used to group the main help listing."""
        return "Other"

    @property
    def hidden(self) -> bool:
        """True for a command that parses but is left out of the main help (e.g. a renamed one)."""
        return False

    def records_run(self, args: argparse.Namespace) -> bool:
        """Whether ``main`` adds this run to the project's run log (``prepare.data.run_log``).

        False by default; a command that opts in may still return False for one invocation (a
        dry run, say). The record is written after the command returns and never changes its
        exit code.
        """
        return False

    @abstractmethod
    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Configure the argument parser for this command."""
        pass

    @abstractmethod
    def execute(self, args: argparse.Namespace) -> int:
        """Execute the command with parsed arguments."""
        pass

    def run(self, args: argparse.Namespace) -> int:
        """Run ``execute`` with SIGINT, SIGTERM and SIGHUP handled; the parser calls this, not ``execute``.

        The first signal (or Ctrl-C) raises ``KeyboardInterrupt`` inside ``execute``, so its
        ``finally`` blocks and ``with`` exits still write what the command has done; later
        signals are logged, not raised (see ``prepare.utils.termination``). An interrupt that
        leaves ``execute`` is logged, running tools are stopped and 130 is returned. The
        ``Termination`` is available to ``execute`` as ``args._termination``; its ``stop`` is
        this run's stop token.

        On an interrupt only the tools started under that token are terminated
        (``terminate_children(stop=term.stop)``), so another command running in the same
        process keeps its own. A tool started without a token on the main thread needs no
        such call: ``run_secure`` kills its child when the interrupt reaches the waiting
        thread. The process-wide stopping flag is neither set nor cleared here, so no
        ``clear_stopping`` call is needed at the start of a run.
        """
        if not self.graceful_shutdown:
            return self.execute(args)
        try:
            with graceful_termination() as term:
                args._termination = term
                try:
                    return self.execute(args)
                except KeyboardInterrupt:
                    self.logger.error("Interrupted (%s)", term.cause)
                    SecureSubprocess.terminate_children(stop=term.stop)
                    return EXIT_INTERRUPTED
        except KeyboardInterrupt:
            # A signal while the handlers were being installed or restored, outside execute.
            self.logger.error("Interrupted")
            return EXIT_INTERRUPTED

    def fail(self, error: BaseException, context: str) -> int:
        """Log ``error`` under ``context`` and return its exit code.

        A command's ``except`` path returns this, so the process exits with 3 for a
        configuration problem, 4 for a retryable one (network, a lock wait that gave up) and 1
        for any other failure (see ``prepare.core.exceptions.ExitCode``). The traceback is
        always attached: the console handler shows only the one-line message unless the
        console is at DEBUG, and a log file (``--log-file``) keeps the traceback; without either,
        a line points to ``--log-level DEBUG``.
        """
        self.logger.error("%s: %s", context, error, exc_info=error)
        log_traceback_hint()
        return exit_code_for(error)

    # Output. stdout carries the command's result (tables, JSON); stderr carries logging.
    # These methods and ``emit_error_json`` are the only places in the package that write to
    # stdout (``scripts/check_no_print.sh`` enforces this).

    def emit(self, text: str = "") -> None:
        """Write one newline-terminated line of the command's result to stdout."""
        print(text, file=sys.stdout)

    def emit_raw(self, text: str) -> None:
        """Write ``text`` to stdout as is, without adding a newline (e.g. ``DataFrame.to_csv`` output)."""
        sys.stdout.write(text)

    def emit_json(self, payload: Any) -> None:
        """Write ``payload`` to stdout as exactly one JSON document (indent 2, trailing newline)."""
        print(json.dumps(payload, indent=2), file=sys.stdout)


class CommandRegistry:
    """Registry for managing CLI commands."""

    def __init__(self):
        self._commands: Dict[str, BaseCommand] = {}

    def register(self, command: BaseCommand) -> None:
        """Register a command in the registry."""
        self._commands[command.name] = command

    def get_command(self, name: str) -> Optional[BaseCommand]:
        """Get a command by name."""
        return self._commands.get(name)

    def get_all_commands(self) -> Dict[str, BaseCommand]:
        """Get all registered commands."""
        return self._commands.copy()

    def setup_parsers(self, main_parser: argparse.ArgumentParser) -> None:
        """Add one subparser per command plus a hidden subparser per alias."""
        subparsers = main_parser.add_subparsers(
            title="commands",
            dest="command",
            metavar="COMMAND",
            help="Run 'prepare COMMAND --help' for the options of one command",
        )
        subparsers.required = True

        for command in self._commands.values():
            # Omit the `help` kwarg for the canonical name too: passing it makes argparse
            # print a second, flat listing of every command above the grouped epilog in
            # `prepare.cli.main._commands_epilog`. `description` (set unconditionally in
            # `_add_parser`) still drives the per-command `prepare COMMAND --help` text.
            self._add_parser(subparsers, command, command.name, None)
            for alias in command.aliases:
                # Omit the `help` kwarg entirely rather than passing argparse.SUPPRESS: on
                # some Python versions SUPPRESS is rendered literally as "==SUPPRESS==" for
                # subparser choices instead of being hidden. Not passing `help` at all means
                # argparse never records a choices entry for this alias, so it is left out
                # of the listing while still parsing normally.
                self._add_parser(subparsers, command, alias, None)

    @staticmethod
    def _add_parser(subparsers, command: BaseCommand, name: str, help_text: Optional[str]) -> None:
        kwargs = {} if help_text is None else {"help": help_text}
        subparser = subparsers.add_parser(
            name,
            description=command.help,
            formatter_class=DefaultsHelpFormatter,
            **kwargs,
        )
        command.configure_parser(subparser)
        add_global_options(subparser, suppress_defaults=True)
        subparser.set_defaults(func=command.run)


# Global registry instance
command_registry = CommandRegistry()
