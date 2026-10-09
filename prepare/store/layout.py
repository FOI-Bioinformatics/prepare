"""
On-disk layout of the shared data store.

The store is one folder holding a single copy of every downloaded
metagenome, shared across organism projects. This module defines where each
piece lives inside that folder and how to create the folder structure the
first time a root is used.
"""

import json
import logging
import os
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from prepare.core.constants import STORE_LAYOUT, STORE_MARKER
from prepare.data.file_io import unique_temp_path, write_text_atomic

logger = logging.getLogger(__name__)


@dataclass
class StorePaths:
    """Resolved paths for one store root."""

    root: Path
    marker: Path
    catalog: Path
    catalog_lock: Path
    locks: Path
    tmp: Path
    sra: Path
    metadata: Path
    journal: Path


def store_paths(root: Path) -> StorePaths:
    """Compute the layout of a store rooted at ``root`` without touching disk."""
    root = Path(root)
    return StorePaths(
        root=root,
        marker=root / STORE_MARKER,
        catalog=root / "catalog.sqlite",
        catalog_lock=root / "catalog.sqlite.lock",
        locks=root / "locks",
        tmp=root / "tmp",
        sra=root / "sra",
        metadata=root / "metadata",
        journal=root / "journal",
    )


def init_store(root: Path) -> StorePaths:
    """Create the store folders and marker file at ``root`` if not already present.

    Safe to call repeatedly: existing folders are left alone, and an existing
    marker keeps its original id and creation time.
    """
    paths = store_paths(root)

    for directory in (paths.root, paths.locks, paths.tmp, paths.sra, paths.metadata, paths.journal):
        directory.mkdir(parents=True, exist_ok=True)

    if not paths.marker.exists():
        marker = {
            "version": 1,
            "id": str(uuid.uuid4()),
            "created": datetime.now(timezone.utc).isoformat(),
            "layout": STORE_LAYOUT,
        }
        if _publish_marker(paths.marker, json.dumps(marker, indent=2)):
            logger.info("Initialized store at %s", paths.root)

    return paths


def _publish_marker(marker: Path, text: str) -> bool:
    """Create ``marker`` holding ``text`` unless it exists; True when this call created it.

    The complete marker is written to a temporary name and published with ``os.link``, which
    fails when the name exists, so two first-time ``store_init`` runs never overwrite each
    other's marker: the second keeps the first one's id. A filesystem without hard links
    (exFAT, some SMB shares) falls back to checking again and renaming, as before.
    """
    tmp = unique_temp_path(marker)
    try:
        write_text_atomic(tmp, text, fsync=True)
        try:
            os.link(tmp, marker)
        except FileExistsError:
            logger.debug("Store marker %s was created by another process meanwhile", marker)
            return False
        except OSError:
            if marker.exists():
                return False
            os.replace(tmp, marker)
        return True
    finally:
        tmp.unlink(missing_ok=True)


def read_marker(root: Path) -> Optional[dict]:
    """Read the store marker at ``root``, or None when it is missing or unreadable.

    A truncated or corrupt marker answers the same way a missing one does: the caller cannot
    trust this folder to be a store, and its message already says "missing or invalid". A
    warning names the file so the reason is not lost.
    """
    marker_path = store_paths(root).marker
    if not marker_path.exists():
        return None
    try:
        return json.loads(marker_path.read_text())
    except (OSError, ValueError) as e:
        logger.warning("Store marker %s could not be read: %s", marker_path, e)
        return None


def sra_dir(paths: StorePaths, acc: str) -> Path:
    """Directory holding the downloaded files for one SRA accession."""
    return paths.sra / acc


def sidecar_path(paths: StorePaths, acc: str) -> Path:
    """Path to the JSON sidecar recording metadata for one SRA accession."""
    return sra_dir(paths, acc) / f"{acc}.json"


def lock_path(paths: StorePaths, acc: str) -> Path:
    """Path to the lock file guarding concurrent access to one SRA accession."""
    return paths.locks / f"{acc}.lock"
