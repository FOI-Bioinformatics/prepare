"""
The link a project keeps into the shared store.

One dataset is downloaded once, into ``<store>/sra/<ACC>/``. Every project
that uses it gets ``<project>/fastq/<ACC>`` as a symlink to that folder, so
the rest of PREPARE reads the project path exactly as if the reads were
local. This module creates, removes and recognises those links, and reports
links whose target has gone away.
"""

import logging
import os
import shutil
from pathlib import Path
from typing import List, Optional, Union

from prepare.core.exceptions import DataAccessError
from prepare.data.file_io import unique_temp_path
from prepare.store.layout import StorePaths, sra_dir
from prepare.store.locks import touch_dataset_use
from prepare.utils.progress import active_item_level

logger = logging.getLogger(__name__)

# How a project entry points at the store: a symlink whose target is written relative to the
# link ("relative") or as an absolute path ("absolute"), a copy of the folder ("copy"), or
# "auto", which picks relative when the two trees share a parent within one volume and
# absolute otherwise (two separately mounted volumes share no tree that survives a remount).
LINK_MODES = ("auto", "relative", "absolute", "copy")


# Directories that hold one mounted volume each, so two paths sharing only one of these are
# on separate volumes rather than in one tree: macOS mounts under /Volumes/<name>, Linux under
# /mnt/<name> and /media/<user>/<name>.
_MOUNT_PARENTS = (Path("/Volumes"), Path("/mnt"), Path("/media"))


def _is_a_volume_root(common: Path) -> bool:
    """True when ``common`` is only the directory that separate volumes mount into.

    ``/Volumes/A/store`` and ``/Volumes/B/project`` share ``/Volumes``, but they sit on two
    volumes that mount and unmount independently: a relative link between them breaks as soon
    as either moves, so that counts as sharing nothing, the same as sharing only ``/``. Two
    paths inside one volume (``/Volumes/lab/store`` and ``/Volumes/lab/project``) do share a
    tree, and a relative link there is what survives that volume being mounted elsewhere.
    """
    if common in _MOUNT_PARENTS:
        return True
    # /media/<user> holds one directory per volume, the same way /Volumes does.
    return common.parent == Path("/media")


def _shares_a_parent(store_root: Path, project_fastq: Path) -> bool:
    """True when both paths sit under a common directory other than a filesystem or volume root.

    A relative symlink only survives moving the project if the store moves with it, which is
    the case when both live under one shared parent (a lab directory, a scratch mount). Two
    trees whose only common ancestor is ``/``, or a mount directory holding one volume each,
    are unrelated, so an absolute link is safer. ``os.path.commonpath`` raises for paths on
    different drives, which is the same answer.
    """
    try:
        common = Path(os.path.commonpath([str(store_root), str(project_fastq)]))
    except ValueError:
        return False
    if str(common) == common.anchor:
        return False
    return not _is_a_volume_root(common)


def _symlink_target(store_dataset: Path, link_parent: Path, mode: str) -> str:
    """The string to write into a symlink in ``link_parent`` pointing at ``store_dataset``.

    Both paths are already resolved by the caller, so a relative target is computed between
    real locations and a relative ``--fastq-folder`` still gets a link that works from
    anywhere.
    """
    if mode == "relative":
        return os.path.relpath(store_dataset, link_parent)
    return str(store_dataset)


def _is_store_copy(link: Path) -> bool:
    """True when ``link`` is a real folder an earlier copy-mode link made: it holds ``<ACC>.json``.

    A copy of a store dataset carries the store's sidecar with it; a folder the project
    downloaded itself (before the store existed) has none, and holds reads no other copy has.
    """
    return link.is_dir() and not link.is_symlink() and (link / f"{link.name}.json").is_file()


def _refuse_project_data(link: Path, replace_store_copy: bool) -> None:
    """Raise when ``link`` is a real entry that a new link or copy must not replace."""
    if link.is_symlink() or not link.exists():
        return
    if replace_store_copy and _is_store_copy(link):
        return
    raise DataAccessError(
        f"{link} already exists and is not a store link; remove it first if the reads should come from the store"
    )


def _clear_existing(link: Path, replace_store_copy: bool = False) -> Optional[Path]:
    """Make room for a new link at ``link``, refusing to destroy real project data.

    An existing symlink is removed. With ``replace_store_copy``, a real folder that an earlier
    copy-mode link made (it holds ``<ACC>.json``) is renamed aside to a hidden name next to it,
    not removed, and that path is returned so the caller can drop it once the new copy is in
    place (or rename it back if the swap fails). Any other real entry is refused. Returns None
    when nothing was moved aside.
    """
    _refuse_project_data(link, replace_store_copy)
    if link.is_symlink():
        link.unlink()
        return None
    if not link.exists():
        return None
    aside = unique_temp_path(link)
    os.replace(link, aside)
    return aside


def room_shortfall(location: Path, needed: int) -> Optional[int]:
    """The free bytes on ``location``'s filesystem when they are fewer than ``needed``, else None.

    The one free-space rule for staging a copy of a dataset: ``store_adopt`` asks for twice the
    folder's size on the store's filesystem, a copy-mode link for the folder's size on the
    project's. A filesystem that cannot be measured is assumed to have room (logged): refusing on
    an unreadable ``disk_usage`` would be worse than trying.
    """
    try:
        free = shutil.disk_usage(location).free
    except OSError as e:
        logger.warning("Could not check free space on %s: %s", location, e)
        return None
    return free if free < needed else None


def _tree_bytes(folder: Path) -> int:
    """Total bytes of the files under ``folder``; a file that cannot be stat'ed counts as 0."""
    total = 0
    for path in folder.rglob("*"):
        try:
            if path.is_file():
                total += path.stat().st_size
        except OSError:
            continue
    return total


def _copy_into_place(store_dataset: Path, link: Path) -> None:
    """Copy ``store_dataset`` to ``link``, replacing an earlier copy only once the new one is complete.

    The copy is made under a hidden staging name first, so an interrupted copy never leaves a
    partial dataset folder that looks like a download, and an earlier copy at ``link`` is moved
    aside only after the staging copy finished. If the final rename does not happen, for any
    reason including an interrupt, the earlier copy is renamed back, so the project is never left
    with no visible copy at all. The project's filesystem must have room for the staging copy
    (``room_shortfall``); a shortfall, or a copy that fails, raises ``DataAccessError`` naming the
    accession.
    """
    accession = store_dataset.name
    needed = _tree_bytes(store_dataset)
    free = room_shortfall(link.parent, needed)
    if free is not None:
        raise DataAccessError(
            f"Cannot copy {accession} into {link.parent}: {free} bytes free, the copy needs about {needed}"
        )
    staging = unique_temp_path(link)
    try:
        try:
            shutil.copytree(store_dataset, staging)
        except (OSError, shutil.Error) as e:
            raise DataAccessError(f"Cannot copy {accession} into {link.parent}: {e}") from e
        _swap_into_place(staging, link, accession)
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def _swap_into_place(staging: Path, link: Path, accession: str) -> None:
    """Move the finished ``staging`` copy to ``link``, an earlier copy there aside and then away.

    The earlier copy is renamed back whenever the swap did not happen, also on an interrupt
    (``KeyboardInterrupt`` or any other exception), as long as nothing took its place meanwhile.
    """
    aside = _clear_existing(link, replace_store_copy=True)
    swapped = False
    try:
        os.replace(staging, link)
        swapped = True
    except OSError as e:
        raise DataAccessError(f"Cannot copy {accession} into {link.parent}: {e}") from e
    finally:
        if aside is not None:
            if swapped:
                shutil.rmtree(aside, ignore_errors=True)
            elif not os.path.lexists(link):
                os.replace(aside, link)


def link_dataset(
    project_fastq: Union[str, Path],
    accession: str,
    paths: StorePaths,
    mode: str = "auto",
) -> Path:
    """Point ``<project_fastq>/<accession>`` at the store's copy of ``accession``.

    Returns the project-side path. ``mode`` is one of ``auto`` (relative when the store and
    the project share a parent directory, absolute otherwise), ``relative``, ``absolute`` or
    ``copy``, which copies the dataset folder instead of linking it (for a project that must
    keep working when the store is unmounted). An existing symlink at the target is replaced.
    A real directory is replaced only in ``copy`` mode and only when it holds ``<ACC>.json``
    (an earlier copy of the store's dataset), after the new copy is complete; any other real
    directory is refused, since it may hold reads this project downloaded itself. On success,
    records ``accession`` as just used (``prepare.store.locks.touch_dataset_use``), so
    ``store_gc`` keeps it for a grace period even before this project's usage row reaches the
    catalogue.
    """
    if mode not in LINK_MODES:
        raise DataAccessError(f"Unknown link mode '{mode}'; expected one of {', '.join(LINK_MODES)}")

    store_dataset = sra_dir(paths, accession)
    if not store_dataset.is_dir():
        raise DataAccessError(f"Store has no dataset folder for {accession}: {store_dataset}")

    project_path = Path(project_fastq)
    project_path.mkdir(parents=True, exist_ok=True)
    link = project_path / accession
    # Resolved forms: the caller may pass a relative folder, and a relative symlink target
    # can only be worked out between two absolute paths.
    resolved_project = project_path.resolve()
    resolved_dataset = store_dataset.resolve()

    if mode == "copy":
        # Refused before any copying when the entry holds the project's own reads.
        _refuse_project_data(link, replace_store_copy=True)
        _copy_into_place(store_dataset, link)
        touch_dataset_use(paths, accession)
        logger.log(active_item_level(), "Copied %s from the store into %s", accession, link)
        return link

    _clear_existing(link)

    if mode == "auto":
        mode = "relative" if _shares_a_parent(paths.root.resolve(), resolved_project) else "absolute"

    target = _symlink_target(resolved_dataset, resolved_project, mode)
    try:
        os.symlink(target, link, target_is_directory=True)
    except OSError as e:
        raise DataAccessError(f"Cannot link {accession} into {project_path}: {e}") from e
    touch_dataset_use(paths, accession)
    logger.log(active_item_level(), "Linked %s to the store copy at %s", link, target)
    return link


def unlink_dataset(project_fastq: Union[str, Path], accession: str) -> bool:
    """Remove ``<project_fastq>/<accession>`` when it is a symlink; return whether it was removed.

    Only a symlink is ever removed: a real directory holds this project's own reads, and the
    store's copy is never touched either way.
    """
    link = Path(project_fastq) / accession
    if not link.is_symlink():
        return False
    try:
        link.unlink()
    except OSError as e:
        raise DataAccessError(f"Cannot remove the store link {link}: {e}") from e
    return True


def is_store_link(path: Union[str, Path], paths: StorePaths) -> bool:
    """True when ``path`` is a symlink pointing into this store's ``sra`` folder.

    A dangling link still counts: it says where the dataset was meant to come from, which is
    what a caller deciding whether to relink needs to know.
    """
    link = Path(path)
    if not link.is_symlink():
        return False
    target = Path(os.path.realpath(link))
    store_sra = Path(os.path.realpath(paths.sra))
    return target == store_sra or store_sra in target.parents


def dangling_links(project_fastq: Union[str, Path]) -> List[str]:
    """Names in ``project_fastq`` that are symlinks with nothing at the other end, sorted.

    Reported rather than repaired: a missing target usually means the store is unmounted, in
    which case removing the links would lose the record of what the project uses.
    """
    folder = Path(project_fastq)
    if not folder.is_dir():
        return []
    return sorted(entry.name for entry in folder.iterdir() if entry.is_symlink() and not entry.exists())
