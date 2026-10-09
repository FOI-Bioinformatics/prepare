"""
The per-accession lock that guards one dataset's folder in the shared store.

A store download or adoption owns ``<store>/sra/<ACC>`` for as long as the transfer takes,
which for a large run is minutes to hours. The lock is ``<store>/locks/<ACC>.lock``, taken
through ``prepare.utils.lockfile.held_lock``: the holder's heartbeat refreshes the file's
mtime every ``LOCK_HEARTBEAT_SECONDS``, so a waiter waits for as long as the holder stays
alive, with no overall timeout (giving up on a legitimate multi-hour download would be the
bug, not the wait). A lock whose heartbeat stopped more than ``DATASET_LOCK_STALE_SECONDS``
ago, or whose holder is known to have exited on this host, is taken over with a warning. A
caller that would rather not wait passes ``wait_seconds``; the resulting error names the
accession and the holder. A caller that would rather skip a held lock than wait for it at
all (``store_gc`` and ``store_verify``, which touch the store opportunistically rather than
as a project's main job) passes ``blocking=False``, which raises ``LockHeld`` at once.

``touch_dataset_use``/``last_dataset_use`` are a second, unrelated marker under the same
``<store>/locks/`` folder, ``<ACC>.used``: its mtime records when a project last started
using the dataset (a link, an adopt), independent of whether the accession's own lock is
held right now. ``store_gc`` keeps a dataset touched within ``GC_RECENT_USE_GRACE_SECONDS``
even when nothing currently holds its lock and its catalogue usage row has not caught up yet.

The limits are read from this module at call time, so tests can shrink them.
"""

import logging
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Callable, Iterator, Optional

from prepare.core.constants import DATASET_LOCK_STALE_SECONDS, LOCK_HEARTBEAT_SECONDS
from prepare.core.settings import setting_or, settings_or
from prepare.store.layout import StorePaths, lock_path
from prepare.utils.lockfile import (
    LOCK_WAIT_LOG_SECONDS,
    LockHeld,
    LockPolicy,
    LockWaitStopped,
    describe_holder,
    held_lock,
    holder_is_dead,
    read_holder,
)

logger = logging.getLogger(__name__)

# How often a waiter re-checks a lock it is waiting for.
LOCK_POLL_SECONDS = 1.0

__all__ = [
    "DATASET_LOCK_STALE_SECONDS",
    "LOCK_HEARTBEAT_SECONDS",
    "LOCK_WAIT_LOG_SECONDS",
    "LockHeld",
    "LockWaitStopped",
    "dataset_lock",
    "last_dataset_use",
    "lock_holder",
    "lock_is_held",
    "read_holder",
    "touch_dataset_use",
]


def _dataset_policy(accession: str, wait_seconds: float) -> LockPolicy:
    """The dataset lock's policy: runtime settings when set, else this module's current limits."""
    stale, heartbeat = settings_or(dataset_lock_stale=DATASET_LOCK_STALE_SECONDS, lock_heartbeat=LOCK_HEARTBEAT_SECONDS)
    return LockPolicy(
        what=accession,
        stale_seconds=stale,
        wait_seconds=wait_seconds,
        poll_seconds=LOCK_POLL_SECONDS,
        heartbeat_seconds=heartbeat,
    )


@contextmanager
def dataset_lock(
    paths: StorePaths,
    accession: str,
    wait_seconds: float = 0.0,
    should_stop: Optional[Callable[[], bool]] = None,
    blocking: bool = True,
) -> Iterator[Path]:
    """Hold ``<store>/locks/<ACC>.lock`` for the length of the block.

    ``wait_seconds`` of zero (the default) waits for as long as the current holder stays
    alive, since a download legitimately runs for hours; a positive value gives up after
    that many seconds with ``LockWaitTimeout`` (a ``LockHeld`` and a ``LockTimeoutError``, exit
    code 4) naming the accession and the holder.
    ``should_stop``, when given, ends a wait for a held lock with ``LockWaitStopped`` as
    soon as it returns True; a free lock is taken regardless. ``blocking=False`` raises
    ``LockHeld`` at once instead of waiting (after any stale or dead-holder takeover), for a
    caller such as ``store_gc`` or ``store_verify`` that would rather skip a dataset another
    run is working on than wait for it.
    """
    lock = lock_path(paths, accession)
    lock.parent.mkdir(parents=True, exist_ok=True)
    with held_lock(lock, _dataset_policy(accession, wait_seconds), should_stop=should_stop, blocking=blocking) as held:
        yield held


def lock_is_held(paths: StorePaths, accession: str) -> bool:
    """True when another run is working on ``accession`` right now.

    A lock file whose heartbeat stopped more than ``DATASET_LOCK_STALE_SECONDS`` ago, or
    whose holder has exited on this host, does not count: callers that only inspect (garbage
    collection, adoption) must not treat one as a live download for ever.
    """
    lock = lock_path(paths, accession)
    try:
        age = time.time() - lock.stat().st_mtime
    except OSError:
        return False
    return age <= setting_or("dataset_lock_stale", DATASET_LOCK_STALE_SECONDS) and not holder_is_dead(read_holder(lock))


def lock_holder(paths: StorePaths, accession: str) -> str:
    """One-line description of whoever holds ``accession``'s lock, for a message."""
    return describe_holder(read_holder(lock_path(paths, accession)))


def _use_marker_path(paths: StorePaths, accession: str) -> Path:
    """Path to the file whose mtime records when ``accession`` was last handed to a project."""
    return paths.locks / f"{accession}.used"


def touch_dataset_use(paths: StorePaths, accession: str) -> None:
    """Record that a project just started using ``accession``, for ``store_gc``'s grace period.

    Sets ``<store>/locks/<ACC>.used``'s mtime to now, creating the file if it does not exist
    yet. Called wherever a project is handed a dataset from the store (linked, adopted or
    copied), so a dataset that a project has just started using is never removed by a
    concurrent ``store_gc`` run before that project's own usage row makes it into the
    catalogue. Never raises: this is bookkeeping, and a store on a read-only or otherwise
    uncooperative filesystem must not turn it into a failure of the link or download it
    accompanies.
    """
    marker = _use_marker_path(paths, accession)
    try:
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.touch(exist_ok=True)
    except OSError as e:
        logger.debug("Could not record use of %s at %s: %s", accession, marker, e)


def last_dataset_use(paths: StorePaths, accession: str) -> Optional[float]:
    """When ``touch_dataset_use`` last recorded a use of ``accession``, as a Unix timestamp.

    None when it was never recorded (including a store predating this mechanism) or the
    marker cannot be read; callers then fall back to their other checks alone.
    """
    try:
        return _use_marker_path(paths, accession).stat().st_mtime
    except OSError:
        return None
