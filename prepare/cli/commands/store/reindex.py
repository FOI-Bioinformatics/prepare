"""
`store_reindex`: rebuild the SQLite catalogue from the sidecars on disk and replay the journal.
"""

import argparse
import logging
import re
from typing import List, Tuple

from prepare.cli.base import BaseCommand
from prepare.cli.commands.store._shared import _now, _no_store_hint
from prepare.core.constants import SRA_ACCESSION_PATTERN
from prepare.core.exceptions import DataAccessError
from prepare.data import registry_blocks as rb
from prepare.data.file_io import visible_files
from prepare.data.registry import load_registry
from prepare.store import journal
from prepare.store.adopt import rebuild_missing_sidecar
from prepare.store.catalog import REBUILT_WITHOUT_PROJECTS, catalog_write
from prepare.store.layout import StorePaths, sidecar_path, sra_dir, store_paths
from prepare.store.locks import LockHeld, dataset_lock, lock_holder
from prepare.store.resolve import resolve_store_root
from prepare.store.sidecar import Sidecar, read_sidecar

logger = logging.getLogger(__name__)


class StoreReindexCommand(BaseCommand):
    """Command to rebuild the store catalogue from every dataset's sidecar file."""

    @property
    def name(self) -> str:
        """Return the command name."""
        return "store_reindex"

    @property
    def help(self) -> str:
        """Return the command's help text."""
        return (
            "Rebuild the store catalogue from every dataset's sidecar file (a missing sidecar is "
            "rebuilt, unverified, from the files; refuses to run when any sidecar cannot be read)"
        )

    @property
    def group(self) -> str:
        """Return the pipeline-step group this command is listed under."""
        return "Store"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add the command's arguments."""
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")
        parser.add_argument(
            "--registry",
            default=None,
            help="Path to the project registry file (defaults to the nearest prepare_registry.json)",
        )

    @staticmethod
    def _read_all_sidecars(paths: StorePaths) -> Tuple[List[Sidecar], List[str], List[str]]:
        """Every dataset's sidecar, the accessions whose sidecar could not be read, and those with none.

        A reindex rebuilds ``datasets`` from exactly what it reads, so a sidecar missed here
        would look like a dataset that no longer exists. The caller stops rather than acting
        on a partial reading when a sidecar is there but unreadable; a folder with no sidecar
        file at all is rebuilt from its files instead (``_rebuild_missing``), but only when its
        name is an SRA accession (``SRA_ACCESSION_PATTERN``): any other folder in ``sra/`` (notes
        put there by hand, for example) is not a dataset, so it is logged as skipped and gets
        neither a sidecar nor a catalogue row.
        """
        sidecars: List[Sidecar] = []
        unreadable: List[str] = []
        missing: List[str] = []
        for acc_dir in visible_files(paths.sra, dirs=True):
            sc_path = sidecar_path(paths, acc_dir.name)
            if not sc_path.is_file():
                if re.fullmatch(SRA_ACCESSION_PATTERN, acc_dir.name):
                    missing.append(acc_dir.name)
                else:
                    logger.warning("%s in %s is not named like an SRA accession; skipped", acc_dir.name, paths.sra)
                continue
            sidecar = read_sidecar(sc_path)
            if sidecar is None:
                unreadable.append(acc_dir.name)
                continue
            sidecars.append(sidecar)
        return sidecars, unreadable, missing

    def _rebuild_missing(self, paths: StorePaths, missing: List[str]) -> List[Sidecar]:
        """Rebuild the sidecar of each sidecar-less store folder in ``missing``, under its lock.

        The lock is taken without waiting: a held lock means another run is publishing that
        dataset right now, so it is skipped with a warning (its catalogue row, if any, is kept)
        rather than described from a half-written folder.
        """
        rebuilt: List[Sidecar] = []
        for accession in missing:
            try:
                with dataset_lock(paths, accession, blocking=False):
                    sidecar = rebuild_missing_sidecar(accession, sra_dir(paths, accession), paths)
            except LockHeld:
                self.logger.warning(
                    "%s: no sidecar and its lock is held (%s); skipped, rerun store_reindex once that run ends",
                    accession,
                    lock_holder(paths, accession),
                )
                continue
            self.logger.warning(
                "%s: had no sidecar; rebuilt one from the files on disk with state '%s' "
                "(run 'store_verify --rescan --fix-state %s' to check the files)",
                accession,
                sidecar.state,
                accession,
            )
            rebuilt.append(sidecar)
        return rebuilt

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

        try:
            paths = store_paths(root)
            sidecars, unreadable, missing = self._read_all_sidecars(paths)
            if unreadable:
                self.logger.error(
                    "Not reindexing: %d sidecar(s) could not be read, and rebuilding from a "
                    "partial reading would drop those datasets and their usage history: %s",
                    len(unreadable),
                    ", ".join(unreadable),
                )
                return 1
            sidecars.extend(self._rebuild_missing(paths, missing))
            with catalog_write(paths) as catalog:
                # Replay first: a fresh or rebuilt catalogue has no projects yet, so usage rows
                # restored here can insert "unknown" placeholder datasets for accessions the
                # journal references. reindex() then removes any placeholder (and real) row
                # whose accession is not among the sidecars just read and whose folder is gone
                # from disk, so a dataset gc already removed is never resurrected by replay.
                projects, usage = journal.replay(paths, catalog)
                count = catalog.reindex(sidecars)
                has_datasets = catalog.conn.execute("SELECT 1 FROM datasets LIMIT 1").fetchone() is not None
                if projects == 0 and has_datasets:
                    catalog.set_meta(REBUILT_WITHOUT_PROJECTS, _now())
                else:
                    # Either at least one project was restored, or the store now holds no
                    # dataset at all (nothing for gc to mistake as unused either way): a flag
                    # set by an earlier, emptier reindex must not survive as stale.
                    catalog.delete_meta(REBUILT_WITHOUT_PROJECTS)
            if projects == 0:
                self._warn_no_projects_restored(paths, has_datasets)
        except DataAccessError as e:
            return self.fail(e, self.name)

        self.emit(f"Reindexed {count} dataset(s); restored {projects} project(s) and {usage} usage record(s)")
        return 0

    def _warn_no_projects_restored(self, paths: StorePaths, flagged: bool) -> None:
        """Warn that no project was restored. With no datasets there is nothing gc could take for
        unused, so the ``rebuilt_without_projects`` flag is neither set nor cleared in that case."""
        if not flagged:
            self.logger.warning("No project records could be restored from %s", paths.journal)
            return
        self.logger.warning(
            "No project records could be restored from %s; every dataset will look unused until each "
            "project runs store_init or store_link again. store_gc now refuses to run (catalogue flag '%s') "
            "until that is done and it is run once with --accept-rebuilt, or until a later store_reindex "
            "restores at least one project",
            paths.journal,
            REBUILT_WITHOUT_PROJECTS,
        )
