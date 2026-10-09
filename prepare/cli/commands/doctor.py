"""``prepare doctor``: check that this environment can run PREPARE.

The checks themselves live in ``prepare.processing.doctor_report``; this command parses the
options, runs them, and writes either one line per check or, with ``--json``, one JSON document.
It exits 0 when no check failed (warnings allowed) and 3 when one did.
"""

import argparse
from typing import List

from prepare.cli.base import BaseCommand, command_registry
from prepare.core.exceptions import ExitCode
from prepare.processing.doctor_report import FAIL, Check, overall_status, run_checks, summary_counts

# Set by main() on the namespace when the settings could not be activated (a config file or
# PREPARE_* variable that does not parse), so the config check reports it.
STARTUP_ERROR_ATTR = "_settings_error"


class DoctorCommand(BaseCommand):
    """Report tools, configuration, store, disk space, registry, CPUs and memory, and optionally network."""

    @property
    def name(self) -> str:
        """The command name."""
        return "doctor"

    @property
    def help(self) -> str:
        """One line for the command list."""
        return "Check external tools, configuration, store, disk space and resources"

    @property
    def group(self) -> str:
        """Listed under Environment in the main help."""
        return "Environment"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add ``--for``, ``--json``, ``--network``, ``--project`` and ``--data-root``."""
        parser.add_argument(
            "--for",
            dest="for_command",
            metavar="COMMAND",
            help="Fail (rather than warn) when a tool this command needs is missing, e.g. extract_target_reads",
        )
        parser.add_argument("--json", action="store_true", help="Write the report as one JSON document")
        parser.add_argument(
            "--network",
            action="store_true",
            help="Also check that NCBI and the Branchwater server answer (10 s each); off by default",
        )
        parser.add_argument(
            "--project",
            default=".",
            help="Project folder: its registry (or the nearest one above it) and free space are checked",
        )
        parser.add_argument(
            "--data-root",
            default=None,
            help="Shared data store to check (default: PREPARE_DATA, the registry's store, or the config file)",
        )

    def _unknown_command(self, name: str) -> bool:
        """Whether ``name`` is not a command listed in ``prepare --help`` (a hidden former name is not one)."""
        command = command_registry.get_all_commands().get(name)
        return command is None or command.hidden

    def _write_text(self, checks: List[Check]) -> None:
        width = max(len(check.name) for check in checks)
        for check in checks:
            self.emit(f"[{check.status}]{' ' * (5 - len(check.status))} {check.name:<{width}}  {check.detail}")
        counts = summary_counts(checks)
        self.emit(f"Result: {overall_status(checks)} ({counts['ok']} ok, {counts['warn']} warn, {counts['fail']} fail)")

    def execute(self, args: argparse.Namespace) -> int:
        """Run every check and report; 0 without a failed check, 3 with one, 2 for an unknown ``--for``."""
        for_command = getattr(args, "for_command", None)
        if for_command and self._unknown_command(for_command):
            self.logger.error("--for %s: no such command (see 'prepare --help' for the current names)", for_command)
            return int(ExitCode.USAGE)
        checks = run_checks(
            project=getattr(args, "project", "."),
            for_command=for_command,
            network=bool(getattr(args, "network", False)),
            data_root=getattr(args, "data_root", None),
            args=args,
            startup_error=getattr(args, STARTUP_ERROR_ATTR, None),
        )
        status = overall_status(checks)
        if getattr(args, "json", False):
            self.emit_json(
                {
                    "status": status,
                    "for": for_command,
                    "counts": summary_counts(checks),
                    "checks": [check.to_dict() for check in checks],
                }
            )
        else:
            self._write_text(checks)
        return int(ExitCode.CONFIGURATION) if status == FAIL else 0
