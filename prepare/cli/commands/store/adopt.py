"""
`store_adopt`: fold an existing project `fastq/` folder into the store.
"""

import argparse
import logging
from pathlib import Path
from typing import Any, List

from prepare.cli.base import BaseCommand
from prepare.cli.commands.sra_verdicts import store_sidecar_verdict, store_spot_count, store_verdict
from prepare.cli.commands.store._shared import _no_store_hint, _gitignore_guard
from prepare.core.exceptions import DataAccessError
from prepare.core.settings import setting_for
from prepare.data import registry_blocks as rb
from prepare.data.registry import load_registry, record_download, registry_transaction, update_linked
from prepare.store.adopt import adopt
from prepare.store.layout import StorePaths, store_paths
from prepare.store.locks import touch_dataset_use
from prepare.store.resolve import resolve_store_root
from prepare.store.usage import ensure_project_identity, record_usage_many

logger = logging.getLogger(__name__)


class StoreAdoptCommand(BaseCommand):
    """Command to fold a project's own downloaded FASTQ folders into the shared store.

    Each accession is staged and moved under its own per-accession lock (see
    ``prepare.store.adopt``), so running this command concurrently against the same store from
    two projects is safe: a shared accession simply serialises rather than racing, and only the
    project's own ``fastq/<ACC>`` folders are ever claimed. Store folders belonging to other
    projects, accessions another run is publishing right now, and accessions the store's
    filesystem has no room to stage are reported and left alone.
    """

    @property
    def name(self) -> str:
        """Return the command name."""
        return "store_adopt"

    @property
    def help(self) -> str:
        """Return the command's help text."""
        return (
            "Move or copy project-owned FASTQ folders into the shared store, then link them back "
            "(each accession is locked, so this is safe to run concurrently from several projects)"
        )

    @property
    def group(self) -> str:
        """Return the pipeline-step group this command is listed under."""
        return "Store"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add the command's arguments."""
        parser.add_argument("--fastq-folder", default="fastq", help="Folder holding per-accession FASTQ downloads")
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")
        parser.add_argument(
            "--registry",
            default=None,
            help="Path to the project registry file (defaults to the nearest prepare_registry.json)",
        )
        mode = parser.add_mutually_exclusive_group()
        mode.add_argument(
            "--move",
            dest="move",
            action="store_true",
            default=True,
            help="Remove each accession's project folder and link it to the store's copy (default)",
        )
        mode.add_argument(
            "--copy",
            dest="move",
            action="store_false",
            help="Leave each accession's project folder as is; the store keeps its own copy, unlinked",
        )
        parser.add_argument(
            "--dry-run", action="store_true", help="Report what would be adopted without changing anything"
        )
        parser.add_argument(
            "--compress",
            dest="compress",
            action="store_true",
            default=True,
            help="Gzip-compress plain FASTQ files while adopting them (default: on)",
        )
        parser.add_argument(
            "--no-compress", dest="compress", action="store_false", help="Leave FASTQ files uncompressed"
        )
        parser.add_argument(
            "--metadata-folder",
            default="metadata",
            help="Folder holding NCBI metadata XML, consulted for each accession's recorded spot count",
        )
        parser.add_argument(
            "--lock-wait",
            dest="lock_wait",
            type=float,
            default=None,
            help=(
                "Seconds to wait for another project's work on the same accession before giving "
                "up on it (default: PREPARE_LOCK_WAIT, config [runtime] lock_wait, or 0, wait for "
                "as long as the other project keeps working)"
            ),
        )

    def execute(self, args: argparse.Namespace) -> int:
        """Run the command; return the exit code."""
        try:
            registry = load_registry(args.registry)
            root = resolve_store_root(args.data_root, rb.store_block(registry).root)
        except DataAccessError as e:
            return self.fail(e, self.name)

        if root is None:
            _no_store_hint()
            return 1

        term = getattr(args, "_termination", None)
        try:
            paths = store_paths(root)
            report = adopt(
                args.fastq_folder,
                paths,
                move=args.move,
                dry_run=args.dry_run,
                compress=args.compress,
                metadata_folders=[Path(args.metadata_folder), paths.metadata],
                lock_wait=setting_for(args, "lock_wait"),
                should_stop=None if term is None else term.stop.is_set,
            )
        except DataAccessError as e:
            return self.fail(e, self.name)

        if args.dry_run:
            self.emit(f"Would adopt {len(report.planned)} dataset(s)")
            if report.planned:
                self.emit("  " + ", ".join(sorted(report.planned)))
            if report.conflicts:
                self.emit(f"Conflicts (left in place): {', '.join(sorted(report.conflicts))}")
            for label, accessions in (
                ("Empty folders, not adopted", report.empty),
                ("Failed (store copy not verified), project copy kept", report.failed),
            ):
                if accessions:
                    self.emit(f"{label}: {', '.join(sorted(accessions))}")
            self._raise_if_stopped(term)
            return 0

        # Only an accession the project now links to needs its download record pointed at the
        # store: one freshly adopted or deduplicated under --move (both replace the project's
        # folder with a link). Under --copy, report.adopted is always empty and a dedup leaves
        # the project's folder exactly as it was, real and unlinked, the same as a fresh --copy
        # adoption (report.copied); those are not linked, but still count as usage (below), so
        # store_gc does not see the store's copy as unused just because this project kept its
        # own copy too.
        newly_linked = sorted(set(report.adopted) | (set(report.deduplicated) if args.move else set()))
        if newly_linked:
            self._record_linked(args, paths, newly_linked)

        copied = sorted(set(report.copied) | (set(report.deduplicated) if not args.move else set()))
        if copied:
            self._record_copied(args, paths, copied)

        self._print_report(report)
        self._raise_if_stopped(term)
        return 0

    @staticmethod
    def _raise_if_stopped(term: Any) -> None:
        """Raise ``KeyboardInterrupt`` when a signal stopped ``adopt`` before every accession.

        Called only after everything ``adopt`` actually finished has been recorded above, so
        ``BaseCommand.run`` converting this to exit 130 never costs the run its bookkeeping.
        """
        if term is not None and term.stop.is_set():
            raise KeyboardInterrupt("store_adopt stopped")

    @staticmethod
    def _record_linked(args: argparse.Namespace, paths: StorePaths, newly_linked: List[str]) -> None:
        """Point the registry's download records at the store and record the linked usage.

        The verdict recorded is ``sra_verdicts.store_verdict``'s, as ``store_link`` and
        ``download_sra`` record it: the new sidecar's verdict judged against the project's spot
        count and merged with the verdict on file, so a recorded ``truncated`` verdict is never
        replaced by ``unverified``. The sidecars and spot counts are read before the lock is taken.
        """
        snapshot = load_registry(args.registry)
        metadata = Path(args.metadata_folder)
        inputs = {
            acc: (store_sidecar_verdict(paths, acc), store_spot_count(snapshot, acc, paths, metadata_folder=metadata))
            for acc in newly_linked
        }
        with registry_transaction(args.registry) as reg:
            ensure_project_identity(reg)
            for acc in newly_linked:
                complete = store_verdict(reg, acc, *inputs[acc])
                record_download(
                    reg,
                    acc,
                    "downloaded",
                    args.fastq_folder,
                    attempt=False,
                    complete=complete,
                    source="store",
                    store_name=acc,
                )
            for acc in newly_linked:
                update_linked(reg, acc, add=True)
            usage_registry = reg
        # Recorded outside the transaction: the catalogue lock is a separate wait, and
        # holding the project's registry lock while queueing for it can time the registry
        # write out.
        record_usage_many(paths, usage_registry, [(acc, "", "linked", "store_adopt") for acc in newly_linked])
        _gitignore_guard(Path.cwd(), logger)

    @staticmethod
    def _record_copied(args: argparse.Namespace, paths: StorePaths, copied: List[str]) -> None:
        """Record usage for accessions left as the project's own, unlinked copy under --copy.

        Neither a fresh --copy adoption (``report.copied``) nor a --copy dedup
        (``report.deduplicated`` when ``args.move`` is False) points the project's download
        record at the store or touches ``registry.store["linked"]``, since the project's
        folder is real, not a link. Without a usage row, store_gc would still see the store's
        copy as unused, even though this project depends on it. These accessions also never
        go through ``link_dataset`` (there is nothing to link under --copy), so their "just
        used" marker is set here instead, the same protection a linked accession gets from
        ``link_dataset`` itself.
        """
        with registry_transaction(args.registry) as reg:
            ensure_project_identity(reg)
            usage_registry = reg
        record_usage_many(paths, usage_registry, [(acc, "", "copied", "store_adopt --copy") for acc in copied])
        for acc in copied:
            touch_dataset_use(paths, acc)

    def _print_report(self, report: Any) -> None:
        self.emit(
            f"Adopted {len(report.adopted)}, copied {len(report.copied)}, "
            f"deduplicated {len(report.deduplicated)}, conflicts {len(report.conflicts)}, "
            f"skipped {len(report.skipped)}, empty {len(report.empty)}, failed {len(report.failed)}"
        )
        for label, accessions in (
            ("Resumed after an interrupted run", report.resumed),
            ("Left to their own project (no sidecar, not ours)", report.foreign),
            ("In progress elsewhere", report.in_progress),
            ("Refused for lack of free space", report.refused),
            ("Empty folders, not adopted", report.empty),
            ("Failed to stage, project copy kept", report.failed),
        ):
            if accessions:
                self.emit(f"{label}: {', '.join(sorted(accessions))}")
        if report.conflicts:
            self.logger.warning("Conflicting accessions left in place: %s", ", ".join(sorted(report.conflicts)))
        if report.failed:
            self.logger.warning(
                "Accessions that failed to stage, project copy kept: %s", ", ".join(sorted(report.failed))
            )
