"""Locks that keep two processes from building the same minimap2 index, or extracting the
same sample, at once.

Both reuse ``prepare.utils.lockfile.held_lock``. The index build lock is blocking: a
second process that asks for the same index while the first is still building it waits,
then reuses what the first one built rather than building it again. The per-sample lock is
not: when it is already held, the caller (``read_extraction.extract_target_reads``) treats
that sample as being worked on elsewhere -- typically an overlapping SLURM shard that
selected the same sample independently -- and reports it skipped rather than waiting or
mapping it twice.
"""

import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Callable, Iterator, Optional

from prepare.core.constants import DATASET_LOCK_STALE_SECONDS, LOCK_HEARTBEAT_SECONDS
from prepare.core.settings import settings_or
from prepare.data.file_io import ensure_directory
from prepare.utils.lockfile import LockHeld, LockPolicy, held_lock

logger = logging.getLogger(__name__)

__all__ = [
    "LockHeld",
    "SAMPLE_LOCKS_DIRNAME",
    "index_build_lock",
    "index_lock_path",
    "sample_extraction_lock",
    "sample_lock_path",
]

#: Subfolder of an extraction output root holding per-sample locks. A leading dot keeps it
#: invisible to every listing built on ``prepare.data.file_io.visible_files``, which is
#: every folder listing in PREPARE, including ``status`` and ``scan_extractions``.
SAMPLE_LOCKS_DIRNAME = ".locks"

#: How often a waiter re-checks a lock it is waiting for; a module-level name (rather than a
#: literal inside the policy below) so a test can shrink it with ``monkeypatch.setattr``.
LOCK_POLL_SECONDS = 1.0


def index_lock_path(index_path: Path) -> Path:
    """Where the build lock for one minimap2 index lives: ``<index>.lock``."""
    return index_path.with_name(index_path.name + ".lock")


def sample_lock_path(output_root: Path, accession: str, genome_id: str) -> Path:
    """Where the per-sample extraction lock for ``accession``/``genome_id`` lives."""
    return Path(output_root) / SAMPLE_LOCKS_DIRNAME / f"{accession}.{genome_id}.lock"


@contextmanager
def index_build_lock(index_path: Path, should_stop: Optional[Callable[[], bool]] = None) -> Iterator[None]:
    """Hold the build lock for one minimap2 index, waiting for another builder to finish.

    A second process asking for the same index while the first is still building it waits
    here; once the first releases the lock, the caller re-checks whether the index is now
    current and reuses it instead of building it again. Blocking, with no wait limit, since
    a build can take minutes on a large genome.
    """
    lock = index_lock_path(index_path)
    ensure_directory(lock.parent)
    stale, heartbeat = settings_or(dataset_lock_stale=DATASET_LOCK_STALE_SECONDS, lock_heartbeat=LOCK_HEARTBEAT_SECONDS)
    policy = LockPolicy(
        what=f"index {index_path.name}",
        stale_seconds=stale,
        heartbeat_seconds=heartbeat,
        poll_seconds=LOCK_POLL_SECONDS,
    )
    with held_lock(lock, policy, should_stop=should_stop):
        yield


@contextmanager
def sample_extraction_lock(output_root: Path, accession: str, genome_id: str) -> Iterator[None]:
    """Hold one sample's extraction lock without waiting.

    Raises ``LockHeld`` at once when another process or thread already holds it, most often
    an overlapping SLURM shard that selected the same sample independently, so the caller
    can report the sample skipped rather than mapping it a second time.
    """
    lock = sample_lock_path(output_root, accession, genome_id)
    ensure_directory(lock.parent)
    stale, heartbeat = settings_or(dataset_lock_stale=DATASET_LOCK_STALE_SECONDS, lock_heartbeat=LOCK_HEARTBEAT_SECONDS)
    policy = LockPolicy(
        what=f"extraction of {accession} against {genome_id}",
        stale_seconds=stale,
        heartbeat_seconds=heartbeat,
    )
    with held_lock(lock, policy, blocking=False):
        yield
