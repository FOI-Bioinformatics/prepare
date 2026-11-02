"""Replace an adopted project folder with a store link without ever losing the only copy.

``store_adopt --move`` removes a project's ``fastq/<ACC>`` folder once the store holds the same
data. The helpers here make that removal safe:

- ``link_blocker`` re-checks, under the dataset lock, that the store copy still exists, is in
  a ready state and holds every non-FASTQ file of the project folder (a dedup compares only
  the FASTQ files, so a note or QC output kept beside them would otherwise be lost).
- ``swap_for_link`` renames the project folder aside (``.<ACC>.adopt-old``, a hidden name
  every scanner skips), links, and renames it back when the link fails; the caller removes the
  set-aside copy only once the link exists.
- ``restore_set_aside`` puts back a set-aside folder a crash left behind with no link in its
  place, and removes one whose link was made.
"""

import logging
import shutil
from pathlib import Path
from typing import Callable, List

from prepare.data.sra import STORE_READY_STATES, fastq_files
from prepare.store.layout import StorePaths, sidecar_path, sra_dir
from prepare.store.link import link_dataset
from prepare.store.sidecar import read_sidecar

logger = logging.getLogger(__name__)

ASIDE_SUFFIX = ".adopt-old"


def aside_path(entry: Path) -> Path:
    """Where ``swap_for_link`` keeps ``entry`` while the link is made."""
    return entry.with_name(f".{entry.name}{ASIDE_SUFFIX}")


def store_lacks(entry: Path, store_dir: Path, ignored: Callable[[str], bool]) -> List[str]:
    """Relative paths of the non-FASTQ files in ``entry`` that ``store_dir`` has no copy of.

    FASTQ files are compared by the caller (by content, across compression), so only the other
    files are checked here, by relative path. Names the staging copy skips (``ignored``) are
    temporary files nothing needs.
    """
    reads = set(fastq_files(entry))
    return sorted(
        str(path.relative_to(entry))
        for path in entry.rglob("*")
        if path.is_file()
        and path not in reads
        and not ignored(path.name)
        and not (store_dir / path.relative_to(entry)).exists()
    )


def link_blocker(accession: str, entry: Path, paths: StorePaths, ignored: Callable[[str], bool]) -> str:
    """Why ``entry`` must not be replaced by a link to the store's copy, or "" when it may be.

    The caller holds the dataset lock, so the answer holds until it releases it.
    """
    store_dir = sra_dir(paths, accession)
    if not store_dir.is_dir():
        return "the store copy is gone"
    sidecar = read_sidecar(sidecar_path(paths, accession))
    if sidecar is None or sidecar.state not in STORE_READY_STATES:
        return "the store copy is not a verified copy"
    missing = store_lacks(entry, store_dir, ignored)
    if missing:
        shown = ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else "")
        return f"the store copy lacks {len(missing)} file(s) of the project folder ({shown})"
    return ""


def swap_for_link(accession: str, entry: Path, paths: StorePaths, project_dir: Path) -> Path:
    """Rename ``entry`` aside and link the store copy in its place; return the set-aside path.

    The caller holds the dataset lock. When the link fails (or is interrupted) the folder is
    renamed back before the exception propagates, so the project never ends up with neither.
    """
    aside = aside_path(entry)
    if aside.exists():
        raise OSError(f"{aside} already exists; move it back to {entry.name} or remove it, then rerun")
    entry.rename(aside)
    linked = False
    try:
        link_dataset(project_dir, accession, paths, locked=True)
        linked = True
    finally:
        if not linked:
            if entry.is_symlink():
                entry.unlink()
            if not entry.exists():
                aside.rename(entry)
    return aside


def restore_set_aside(project_dir: Path) -> None:
    """Undo or finish a ``swap_for_link`` that a crash interrupted."""
    if not project_dir.is_dir():
        return
    for aside in project_dir.glob(f".*{ASIDE_SUFFIX}"):
        entry = project_dir / aside.name[1 : -len(ASIDE_SUFFIX)]
        if not aside.is_dir() or aside.is_symlink():
            continue
        if entry.is_symlink() and entry.exists():
            logger.info("%s: removing the copy set aside by an interrupted adopt (linked since)", entry.name)
            shutil.rmtree(aside)
        elif not entry.exists():
            # Nothing there, or a dangling link: the set-aside folder is the project's only copy.
            logger.warning("%s: restoring the project folder an interrupted adopt set aside", entry.name)
            if entry.is_symlink():
                entry.unlink()
            aside.rename(entry)
        else:
            logger.warning("%s: both %s and a set-aside copy exist; left both in place", entry.name, entry)
