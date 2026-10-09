"""CLI command that records datasets excluded on purpose, with a reason."""

import argparse
from pathlib import Path
from typing import Dict, List

from prepare.cli.base import BaseCommand
from prepare.core.exceptions import PrepareError
from prepare.data import registry_blocks as rb
from prepare.data import run_log
from prepare.data.file_io import write_text_atomic
from prepare.data.registry import Registry, clear_exclusion, load_registry, query, record_exclusion
from prepare.data.registry_batch import registry_update


def read_blacklist_file(path: Path) -> Dict[str, str]:
    """Accession -> reason from a blacklist file (``ACC  # reason`` lines; comments and blanks skipped)."""
    entries: Dict[str, str] = {}
    if not path.exists():
        return entries
    for line in path.read_text().splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        accession, _, reason = text.partition("#")
        entries[accession.strip()] = reason.strip()
    return entries


def write_blacklist_file(path: Path, entries: Dict[str, str]) -> None:
    write_text_atomic(
        path, "".join(f"{acc}  # {reason}\n" if reason else f"{acc}\n" for acc, reason in sorted(entries.items()))
    )


class BlacklistCommand(BaseCommand):
    """Exclude datasets from downloads and extraction, recording why."""

    @property
    def name(self) -> str:
        return "blacklist"

    @property
    def help(self) -> str:
        return "Record datasets to exclude, with a reason; keeps blacklist.txt and the registry in step"

    @property
    def group(self) -> str:
        return "Reads"

    def records_run(self, args: argparse.Namespace) -> bool:
        return not getattr(args, "list", False)

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        action = parser.add_mutually_exclusive_group(required=True)
        action.add_argument("--add", nargs="+", metavar="ACCESSION", help="Accessions to exclude")
        action.add_argument("--remove", nargs="+", metavar="ACCESSION", help="Accessions to allow again")
        action.add_argument("--from-file", help="File of accessions to exclude, one per line")
        action.add_argument("--list", action="store_true", help="Show the excluded accessions and reasons")
        parser.add_argument(
            "--reason", default=None, help="Why the accessions are excluded (required with --add and --from-file)"
        )
        parser.add_argument(
            "--blacklist-file", default="blacklist.txt", help="Plain list kept for download_sra --blacklist"
        )
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")

    def execute(self, args: argparse.Namespace) -> int:
        try:
            if args.list:
                registry = load_registry(args.registry)
                for acc in query(registry, "excluded"):
                    self.emit(f"{acc}\t{(rb.exclusion_block(registry, acc) or rb.ExclusionBlock()).reason}")
                return 0
            accessions: List[str] = []
            if not args.remove:
                if not args.reason:
                    raise PrepareError("--reason is required when adding to the blacklist")
                accessions = args.add or _read_plain_list(Path(args.from_file))
            blacklist_path = Path(args.blacklist_file)

            def update(registry: Registry) -> int:
                # blacklist.txt is read and rewritten inside the registry lock, so two blacklist
                # runs serialise on that lock and neither loses the other's edits to either file.
                entries = read_blacklist_file(blacklist_path)
                for acc in args.remove or []:
                    clear_exclusion(registry, acc)
                    entries.pop(acc, None)
                for acc in accessions:
                    record_exclusion(registry, acc, args.reason)
                    entries[acc] = args.reason
                write_blacklist_file(blacklist_path, entries)
                return len(query(registry, "excluded"))

            excluded = registry_update(args.registry, update)
            run_log.note_run(
                args,
                summary={
                    "action": "remove" if args.remove else "add",
                    "accessions": len(args.remove or accessions),
                    "reason": None if args.remove else args.reason,
                    "excluded": excluded,
                },
            )
            if args.remove:
                self.logger.info("Removed %d accession(s) from the blacklist", len(args.remove))
            else:
                self.logger.info("Excluded %d accession(s): %s", len(accessions), args.reason)
            return 0
        except PrepareError as e:
            return self.fail(e, "Error updating the blacklist")


def _read_plain_list(path: Path) -> List[str]:
    if not path.exists():
        raise PrepareError(f"Accession file not found: {path}")
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
