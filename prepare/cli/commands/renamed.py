"""Old command names that now point at their replacement.

In 0.5.0 ``sra_stats`` and ``sra_profile_quality`` became ``sra_profile``, and ``sra_dashboard``
and ``sra_compare`` became ``sra_report``. The old names (and their dash aliases) still parse,
with any arguments, so a script that calls one is told the new name instead of receiving
argparse's list of valid commands. They are left out of the main help.
"""

import argparse
from typing import Any, List, Sequence, Tuple

from prepare.cli.base import BaseCommand

# old name -> the command that replaced it
RENAMED_COMMANDS = {
    "sra_stats": "sra_profile",
    "sra_profile_quality": "sra_profile",
    "sra-profile-quality": "sra_profile",
    "sra_dashboard": "sra_report",
    "sra-dashboard": "sra_report",
    "sra_compare": "sra_report",
    "sra-compare": "sra_report",
}


class RenamedCommand(BaseCommand):
    """A former command name: accepts any arguments, logs where the command went and exits 2."""

    def __init__(self, old: str, new: str) -> None:
        """Point ``old`` at ``new``."""
        super().__init__()
        self.old = old
        self.new = new

    @property
    def name(self) -> str:
        """The former command name."""
        return self.old

    @property
    def help(self) -> str:
        """Help text naming the replacement."""
        return f"Renamed: use prepare {self.new}"

    @property
    def group(self) -> str:
        """The group of the replacement command."""
        return "Reads"

    @property
    def hidden(self) -> bool:
        """Always hidden from the main help."""
        return True

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Accept and ignore every argument the old command took, so the pointer is always shown."""
        original = parser.parse_known_args

        def parse_known_args(
            args: Sequence[str] | None = None, namespace: argparse.Namespace | None = None
        ) -> Tuple[argparse.Namespace, List[Any]]:
            parsed, extras = original(args, namespace)
            parsed.ignored_arguments = list(extras)
            return parsed, []

        # argparse hands a subcommand its arguments through parse_known_args and reports any
        # it did not recognise as an error before execute() runs; this parser keeps them.
        setattr(parser, "parse_known_args", parse_known_args)

    def execute(self, args: argparse.Namespace) -> int:
        """Log the new command name and return 2."""
        self.logger.error("%s was renamed: use prepare %s", self.old, self.new)
        return 2


def renamed_commands() -> List[RenamedCommand]:
    """One ``RenamedCommand`` per entry of ``RENAMED_COMMANDS``."""
    return [RenamedCommand(old, new) for old, new in RENAMED_COMMANDS.items()]
