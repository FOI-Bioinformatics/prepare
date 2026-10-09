"""Batched writes to the project registry.

Every ``registry_transaction`` reads and rewrites the whole registry file. On a project with
some 15,000 datasets that took about one second per transaction, so a command that records
thousands of accessions one transaction at a time ran for hours. ``registry_batch`` queues
mutations and applies them in a few transactions instead. It lives in its own module, next to
``prepare.data.registry``, which is held at its current size by the module size ceiling.
"""

import logging
import time
from dataclasses import fields
from pathlib import Path
from typing import Any, Callable, List, Optional, Tuple, TypeVar, Union

from prepare.data.registry import Registry, load_registry, registry_path, registry_transaction

logger = logging.getLogger(__name__)

# Indirection so tests can advance the clock a batch uses without touching the lock timing.
_monotonic = time.monotonic

RegistryMutation = Callable[[Registry], Any]
T = TypeVar("T")


class RegistryBatch:
    """Queue registry mutations and apply them in a few transactions instead of one each.

    A batch queues each mutation (a function that changes a loaded ``Registry`` in place) and
    applies the queue, in the order queued, inside one ``registry_transaction``: when
    ``flush_every`` mutations are waiting, when ``flush_seconds`` have passed since the last
    flush, on ``flush()``, and when the ``with`` block exits, including by an exception or
    ``KeyboardInterrupt``. ``flush_seconds`` is checked only when a mutation is queued, so it
    bounds the delay per result, not in wall-clock time: mutations queued before a long quiet
    period stay unwritten until the next one arrives or the batch exits.

    The registry lock is held only while a flush runs, so other writers (and other batches)
    are not blocked between flushes; because each flush loads the file inside the lock, their
    changes are never reverted. A mutation must therefore not assume the registry it receives
    is the one an earlier flush saw. ``flush_every`` or ``flush_seconds`` set to None disables
    that trigger, so a batch with both None writes once, on exit.

    A mutation that raises is dropped together with whatever it had already changed: the
    registry is reloaded from the file, which the transaction has not written yet, and the
    mutations that succeeded before it are replayed, so no accession is left half-updated. The
    replay costs one file load per failing mutation and nothing otherwise.

    Functions passed to ``add_flush_hook`` are called with the written registry after each
    flush has released the lock, so work that takes a different lock (the store catalogue)
    never waits while holding this one. ``registry`` holds the registry the last flush wrote,
    or None before the first flush.
    """

    def __init__(
        self,
        path: Optional[Union[str, Path]] = None,
        flush_every: Optional[int] = 50,
        flush_seconds: Optional[float] = 30.0,
    ) -> None:
        """Bind the batch to a registry file; nothing is read or written until the first flush."""
        self.path = registry_path(path)
        self.flush_every = flush_every
        self.flush_seconds = flush_seconds
        self.registry: Optional[Registry] = None
        self._queue: List[Tuple[RegistryMutation, str]] = []
        self._hooks: List[Callable[[Registry], None]] = []
        self._last_flush = _monotonic()
        # Mutations dropped because they raised, counted once their flush has been written.
        self.dropped = 0

    def __len__(self) -> int:
        """The number of mutations queued and not yet written."""
        return len(self._queue)

    def pending_labels(self) -> List[str]:
        """The labels of the mutations queued and not yet written, in the order queued."""
        return [label or "an unlabelled entry" for _, label in self._queue]

    def add_flush_hook(self, hook: Callable[[Registry], None]) -> None:
        """Call ``hook(registry)`` after every flush that wrote something, once the lock is released."""
        self._hooks.append(hook)

    def apply(self, mutation: RegistryMutation, label: str = "") -> None:
        """Queue ``mutation(registry)``; flush if the size or time limit has been reached.

        ``label`` (typically the accession the mutation records) names the mutation in the
        error logged if it raises.
        """
        self._queue.append((mutation, label))
        if self.flush_every is not None and len(self._queue) >= self.flush_every:
            self.flush()
        elif self.flush_seconds is not None and _monotonic() - self._last_flush >= self.flush_seconds:
            self.flush()

    def flush(self) -> Optional[Registry]:
        """Apply every queued mutation in one transaction and return the written registry.

        Returns None, writing nothing, when the queue is empty. A mutation that raises is
        logged at error, with its label, and dropped; the others are still applied and written,
        so one faulty record neither blocks the rest nor fails every later flush. A failure of
        the transaction itself (a lock timeout, a disk error) raises and keeps the whole queue
        for a later attempt, since the queue is cleared only once the file has been written.
        """
        if not self._queue:
            return None
        pending = list(self._queue)
        with registry_transaction(self.path) as registry:
            dropped = self._apply_all(registry, pending)
        del self._queue[: len(pending)]
        self.dropped += dropped
        self._last_flush = _monotonic()
        self.registry = registry
        for hook in self._hooks:
            hook(registry)
        return registry

    def _apply_all(self, registry: Registry, pending: List[Tuple[RegistryMutation, str]]) -> int:
        """Apply ``pending`` to ``registry`` in order, rolling back and dropping any mutation that raises.

        A rollback reloads the registry from the file (unchanged until the transaction writes it)
        and replays the mutations that had succeeded, so a mutation that raises partway through
        leaves nothing of what it changed. Replaying costs one file load per failure and nothing
        when every mutation succeeds. Returns the number of mutations dropped.
        """
        applied: List[Tuple[RegistryMutation, str]] = []
        for mutation, label in pending:
            try:
                mutation(registry)
            except Exception as e:  # noqa: B902 - one faulty record must not stop the others
                logger.error("Dropped a registry update for %s: %s", label or "an unlabelled entry", e, exc_info=True)
                self._restore(registry, applied)
                continue
            applied.append((mutation, label))
        return len(pending) - len(applied)

    def _restore(self, registry: Registry, applied: List[Tuple[RegistryMutation, str]]) -> None:
        """Reset ``registry`` to the file's state and replay ``applied``, dropping any that now raises."""
        fresh = load_registry(self.path)
        for spec in fields(Registry):
            setattr(registry, spec.name, getattr(fresh, spec.name))
        survivors: List[Tuple[RegistryMutation, str]] = []
        for mutation, label in applied:
            try:
                mutation(registry)
            except Exception as e:  # noqa: B902 - a replay that fails is dropped like a first failure
                logger.error("Dropped a registry update for %s on replay: %s", label or "an unlabelled entry", e)
                self._restore(registry, survivors)
                continue
            survivors.append((mutation, label))
        applied[:] = survivors

    def __enter__(self) -> "RegistryBatch":
        """Return the batch itself."""
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        """Flush what is queued, whether the block finished or raised.

        When the block raised (a ``KeyboardInterrupt`` included), a failing flush is logged
        and the original exception propagates rather than the flush error replacing it.
        """
        if exc_type is None:
            self.flush()
            return
        try:
            self.flush()
        except Exception as e:  # noqa: B902 - the block's own exception must propagate, not this one
            logger.error("Could not write %d queued registry update(s) to %s: %s", len(self._queue), self.path, e)


def registry_batch(
    path: Optional[Union[str, Path]] = None,
    flush_every: Optional[int] = 50,
    flush_seconds: Optional[float] = 30.0,
) -> RegistryBatch:
    """Return a ``RegistryBatch`` for ``path``, used as ``with registry_batch(path) as batch: batch.apply(fn)``.

    Like ``registry_transaction``, a batch must not be flushed while the same process holds
    that registry's lock (inside a ``registry_transaction`` block for the same file), since
    the flush would wait for a lock this process already holds.
    """
    return RegistryBatch(path, flush_every=flush_every, flush_seconds=flush_seconds)


def registry_update(path: Optional[Union[str, Path]], mutation: Callable[[Registry], T]) -> T:
    """Load the registry, apply ``mutation`` to it and write it, all under one lock; return its result.

    For a command that works for a while and then records what it found: do the work on a
    snapshot (``load_registry``, no lock), then pass the recording step here, so the file is
    read again inside the lock and whatever other writers committed meanwhile is kept. The lock
    is held only while ``mutation`` runs, so ``mutation`` must be quick (no read counting, no
    network, no store catalogue). When it raises, nothing is written and the exception
    propagates. Not re-entrant, like ``registry_transaction``.
    """
    with registry_transaction(path) as registry:
        result = mutation(registry)
    return result
