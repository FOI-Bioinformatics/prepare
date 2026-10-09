"""
`store_link` and `store_unlink`: add or remove one accession's project symlink.

`store_link` records the verdict ``sra_verdicts.store_verdict`` gives, the same one ``download_sra``
records for a link: the store sidecar's verdict judged against the project's spot count and merged
with the verdict on file, so a relink never turns a recorded ``truncated`` verdict into
``unverified``. A ready copy whose sidecar is ``unverified`` but short against that spot count is
refused like a ``partial`` one; its sidecar is first rewritten as ``partial`` under the dataset lock,
as the download hand-off does.
"""

import argparse
from pathlib import Path
from typing import List, Optional

from prepare.cli.base import BaseCommand
from prepare.cli.commands.sra_verdicts import store_sidecar_verdict, store_spot_count, store_verdict
from prepare.cli.commands.store._shared import _no_store_hint
from prepare.core.exceptions import DataAccessError
from prepare.data import registry_blocks as rb
from prepare.data.registry import load_registry, record_download, registry_transaction, update_linked
from prepare.data.registry_timing import set_download_timing
from prepare.data.sra.store_handoff import _reverify_sidecar, _short_unverified
from prepare.store.layout import StorePaths, sidecar_path, store_paths
from prepare.store.link import LINK_MODES, link_dataset, unlink_dataset
from prepare.store.locks import dataset_lock
from prepare.store.resolve import resolve_store_root
from prepare.store.sidecar import read_sidecar
from prepare.store.usage import ensure_project_identity, record_usage_many
from prepare.utils.lockfile import LockHeld


class StoreLinkCommand(BaseCommand):
    """Command to link project accessions to the shared store's copies."""

    @property
    def name(self) -> str:
        """Return the command name."""
        return "store_link"

    @property
    def help(self) -> str:
        """Return the command's help text."""
        return "Link project accessions to the shared store's copies (complete datasets only, unless --accept-partial)"

    @property
    def group(self) -> str:
        """Return the pipeline-step group this command is listed under."""
        return "Store"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add the command's arguments."""
        parser.add_argument("accessions", nargs="+", help="Accessions to link from the store")
        parser.add_argument("--fastq-folder", default="fastq", help="Folder holding per-accession FASTQ downloads")
        parser.add_argument(
            "--registry",
            default=None,
            help="Path to the project registry file (defaults to the nearest prepare_registry.json)",
        )
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")
        parser.add_argument(
            "--link-mode",
            choices=list(LINK_MODES),
            default="auto",
            help=(
                "How this project points at the store's copy: a relative or absolute symlink, "
                "a copy of the folder, or auto (relative when the store and the project share a "
                "parent folder)"
            ),
        )
        parser.add_argument(
            "--accept-partial",
            dest="accept_partial",
            action="store_true",
            default=False,
            help="Link a dataset whose store copy is incomplete, failed or undescribed",
        )

    def _refuse_incomplete(
        self, paths: StorePaths, accession: str, accept_partial: bool, expected: Optional[int] = None
    ) -> bool:
        """True when ``accession`` must not be linked as it stands.

        The store's sidecar is the record of whether a dataset is usable. A ``partial`` or
        ``failed`` one, or one with no sidecar at all (an interrupted run, or a folder put
        there by hand), reads as complete once it is linked into ``fastq/``, so linking it
        silently would feed an unfinished dataset into every later analysis. A ready copy
        whose sidecar is ``unverified`` but holds fewer reads than ``expected`` (the project's
        spot count) allows is incomplete too (``_short_unverified``); its sidecar is rewritten
        as ``partial`` under the dataset lock first (``_record_short``).
        """
        sidecar = read_sidecar(sidecar_path(paths, accession))
        if sidecar is None:
            state = "no sidecar"
        elif sidecar.state in ("partial", "failed", "downloading"):
            state = sidecar.state
        elif expected is not None and _short_unverified(sidecar, expected):
            self._record_short(paths, accession, expected)
            state = f"partial (was unverified, holding {sidecar.reads_per_mate} of {expected} spots)"
        else:
            return False

        if accept_partial:
            self.logger.warning("%s: linking a dataset the store records as %s", accession, state)
            return False
        self.logger.error(
            "%s: the store records this dataset as %s; rerun with --accept-partial to link it anyway",
            accession,
            state,
        )
        return True

    def _record_short(self, paths: StorePaths, accession: str, expected: int) -> None:
        """Rewrite a short ``unverified`` sidecar as ``partial``, under the dataset lock.

        The decision is taken again under the lock, since another project may have replaced the
        copy meanwhile. A lock another run holds is not waited for: that run is working on the
        dataset, and the refusal stands either way.
        """
        try:
            with dataset_lock(paths, accession, blocking=False):
                if _short_unverified(read_sidecar(sidecar_path(paths, accession)), expected):
                    _reverify_sidecar(paths, accession, expected)
        except LockHeld:
            self.logger.warning("%s: another run holds the dataset lock; its sidecar is left as it is", accession)

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

        paths = store_paths(root)
        accept_partial = getattr(args, "accept_partial", False)
        linked: List[str] = []
        failed: List[str] = []
        # Read before the registry lock is taken: the spot counts may come from metadata XML files.
        spots = {acc: store_spot_count(registry, acc, paths) for acc in args.accessions}
        for accession in args.accessions:
            if self._refuse_incomplete(paths, accession, accept_partial, spots[accession]):
                failed.append(accession)
                continue
            try:
                link_dataset(args.fastq_folder, accession, paths, mode=args.link_mode)
                linked.append(accession)
            except DataAccessError as e:
                self.logger.error("%s: %s", accession, e)
                failed.append(accession)

        if linked:
            sidecars = {acc: store_sidecar_verdict(paths, acc) for acc in linked}
            with registry_transaction(args.registry) as reg:
                ensure_project_identity(reg)
                for accession in linked:
                    complete = store_verdict(reg, accession, sidecars[accession], spots[accession])
                    record_download(
                        reg,
                        accession,
                        "downloaded",
                        args.fastq_folder,
                        attempt=False,
                        complete=complete,
                        source="store",
                        store_name=accession,
                    )
                    # Linked, not downloaded: an earlier download time of this project is not kept.
                    set_download_timing(reg, accession, None, None)
                for accession in linked:
                    update_linked(reg, accession, add=True)
                usage_registry = reg
            # Outside the transaction: see the note in store_adopt.
            record_usage_many(paths, usage_registry, [(acc, "", "linked", "store_link") for acc in linked])

        self.emit(f"Linked {len(linked)} of {len(args.accessions)} accession(s)")
        return 1 if failed else 0


class StoreUnlinkCommand(BaseCommand):
    """Command to remove a project's store link, without touching the store's copy."""

    @property
    def name(self) -> str:
        """Return the command name."""
        return "store_unlink"

    @property
    def help(self) -> str:
        """Return the command's help text."""
        return "Remove a project's link to the store (never removes a real directory)"

    @property
    def group(self) -> str:
        """Return the pipeline-step group this command is listed under."""
        return "Store"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add the command's arguments."""
        parser.add_argument("accessions", nargs="+", help="Accessions to unlink")
        parser.add_argument("--fastq-folder", default="fastq", help="Folder holding per-accession FASTQ downloads")
        parser.add_argument(
            "--registry",
            default=None,
            help="Path to the project registry file (defaults to the nearest prepare_registry.json)",
        )

    def execute(self, args: argparse.Namespace) -> int:
        """Run the command; return the exit code."""
        removed: List[str] = []
        refused: List[str] = []
        for accession in args.accessions:
            link_path = Path(args.fastq_folder) / accession
            if link_path.exists() and not link_path.is_symlink():
                self.logger.error("%s is a real directory, not a store link; refusing to remove it", accession)
                refused.append(accession)
                continue
            try:
                was_removed = unlink_dataset(args.fastq_folder, accession)
            except DataAccessError as e:
                self.logger.error("%s: %s", self.name, e)
                refused.append(accession)
                continue
            if was_removed:
                removed.append(accession)

        if removed:
            with registry_transaction(args.registry) as reg:
                for accession in removed:
                    record_download(reg, accession, "missing", args.fastq_folder, attempt=False)
                    set_download_timing(reg, accession, None, None)
                for accession in removed:
                    update_linked(reg, accession, add=False)

        self.emit(f"Unlinked {len(removed)} of {len(args.accessions)} accession(s)")
        return 1 if refused else 0
