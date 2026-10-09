"""
One lock-file mechanism for every lock PREPARE takes.

A lock is a file created with ``O_EXCL``, which works on local disks, exFAT, NFS and SMB
alike. The file holds a JSON holder record ``{pid, host, started, token, pidns}``: ``token``
is a random value that identifies this one acquisition, and ``pidns`` (Linux only) is the
inode of the holder's pid namespace, so a pid is only interpreted inside the namespace it
belongs to. Each lock is taken under a ``LockPolicy`` that names what it protects and sets
its stale threshold, wait limit, poll interval and heartbeat interval.

While a lock is held, one daemon thread per process refreshes the file's mtime every
``heartbeat_seconds``, so the file's age measures how long ago the holder was last alive,
not how long it has been working. A waiter removes a lock only when that age exceeds
``stale_seconds`` or when the holder is known to be dead (same host, same pid namespace,
and no process with its pid). The removal happens under a second ``O_EXCL`` file,
``<lock>.reclaim``, and only if the lock file is still the one the waiter judged, so two
waiters can never both remove a lock and one of them remove the other's new lock: the lock
is moved aside with ``rename`` and compared again before it is discarded.

Release compares the token under the same guard and never removes a lock that now belongs to
someone else. The heartbeat reads the token and touches the file through one descriptor. A
holder record that is not JSON (a bare pid written by an older version, or a test fixture)
is never judged dead; it follows the age rule alone.

On the main thread, SIGINT, SIGTERM and SIGHUP are blocked (``pthread_sigmask``) from just
before a lock file is created until the acquisition is recorded as held, so a first signal
cannot leave a lock file that nothing releases; a signal arriving meanwhile is delivered as
soon as the lock is recorded. Windows and threads other than the main one get no mask; the
heartbeat thread blocks these signals for good. Remaining windows, all rare and all ending in
a wait rather than in lost data: a signal the kernel delivers to another thread of the
process (a download worker) while the main thread is in that section still reaches Python
there; a ``KeyboardInterrupt`` raised while a release waits for the reclaim guard
(``_wait_for_guard``) leaves the lock file in place until it goes stale; and a process
killed (SIGKILL) while it holds a reclaim guard leaves the guard, so for
``RECLAIM_GUARD_STALE_SECONDS`` every release of that lock waits
``RELEASE_GUARD_WAIT_SECONDS`` and no waiter can reclaim it.
"""

import json
import logging
import os
import secrets
import signal
import socket
import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, List, Optional, Tuple

from prepare.core.exceptions import DataAccessError, LockTimeoutError

logger = logging.getLogger(__name__)

# A waiter says it is waiting once the wait passes FIRST_WAIT_LOG_SECONDS (sub-second
# contention stays quiet), then every LOCK_WAIT_LOG_SECONDS.
FIRST_WAIT_LOG_SECONDS = 1.0
LOCK_WAIT_LOG_SECONDS = 30.0
# A reclaim guard is held for a few file operations; one this old was left by a crash.
RECLAIM_GUARD_STALE_SECONDS = 60.0
# A release compares its token under the reclaim guard and waits this long for a guard in use.
RELEASE_GUARD_WAIT_SECONDS = 5.0
# Errors a heartbeat refresh may raise beyond the OSError it handles itself; any of them marks
# that one lock lost and leaves the heartbeat thread refreshing the others.
_REFRESH_FAILURES = (ArithmeticError, AttributeError, LookupError, OSError, RuntimeError, TypeError, ValueError)
# Signals held back while a lock file is created and recorded as held (see _DeferredSignals).
_DEFERRED_SIGNALS = frozenset(
    getattr(signal, name) for name in ("SIGINT", "SIGTERM", "SIGHUP") if hasattr(signal, name)
)

__all__ = [
    "LockHeld",
    "LockLost",
    "LockPolicy",
    "LockReentry",
    "LockWaitStopped",
    "LockWaitTimeout",
    "describe_holder",
    "held_lock",
    "holder_is_dead",
    "read_holder",
    "verify_held",
]


class LockHeld(DataAccessError):
    """Raised when a lock is taken by another holder and the caller will not wait any longer."""


class LockWaitTimeout(LockHeld, LockTimeoutError):
    """Raised when a blocking wait reaches ``policy.wait_seconds`` with the lock still held.

    Both a ``LockHeld`` (callers that handle a held lock keep doing so) and a
    ``LockTimeoutError``, so a command that stops on it exits with 4 (retryable). A
    non-blocking attempt that finds the lock held raises a plain ``LockHeld`` (exit code 1).
    """


class LockLost(DataAccessError):
    """Raised by ``verify_held`` when a lock this process took is no longer its own."""


class LockReentry(DataAccessError):
    """Raised when a thread asks again for a lock it already holds, which would never be granted."""


class LockWaitStopped(DataAccessError):
    """Raised when a caller's stop predicate ends a wait for a held lock (for example on Ctrl-C)."""


@dataclass(frozen=True)
class LockPolicy:
    """How one kind of lock is waited for, judged stale and kept alive.

    ``what`` names the protected resource in messages. ``wait_seconds`` of zero waits for as
    long as the holder stays alive; a positive value gives up with ``LockHeld`` after that
    long. A lock whose mtime is older than ``stale_seconds`` is taken over.
    """

    what: str
    stale_seconds: float
    wait_seconds: float = 0.0
    poll_seconds: float = 1.0
    heartbeat_seconds: float = 10.0


# ----------------------------------------------------------------- holder record


def _pid_namespace() -> Optional[int]:
    """Inode of this process's pid namespace on Linux, or None where there is none to read."""
    try:
        return os.stat("/proc/self/ns/pid").st_ino
    except OSError:
        return None


def read_holder(lock: Path) -> Dict[str, Any]:
    """The holder record inside ``lock``; ``{"pid": N}`` for a bare pid, ``{}`` when unreadable.

    A lock file caught between creation and its first write, or written in a form this
    version does not know, reads as an empty record: callers then report an unknown holder.
    """
    try:
        data = json.loads(Path(lock).read_text())
    except (OSError, ValueError):
        return {}
    if isinstance(data, dict):
        return data
    if isinstance(data, int) and not isinstance(data, bool):
        return {"pid": data}
    return {}


def describe_holder(holder: Dict[str, Any]) -> str:
    """One-line description of a lock's holder, for a log line or an error message."""
    return (
        f"pid {holder.get('pid', 'unknown')} on host {holder.get('host', 'unknown')} "
        f"since {holder.get('started', 'an unknown time')}"
    )


def holder_is_dead(holder: Dict[str, Any]) -> bool:
    """True only when the holder's process is known to have exited.

    That needs the same host and the same pid namespace as this process and a pid that no
    process has (``ProcessLookupError``). A process of another user (``PermissionError``),
    another host, or a record without a host is treated as alive.
    """
    pid = holder.get("pid")
    if not isinstance(pid, int) or isinstance(pid, bool) or pid <= 0:
        return False
    if holder.get("host") != socket.gethostname() or holder.get("pidns") != _pid_namespace():
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return True
    except OSError:
        return False
    return False


def _create(lock: Path) -> Optional[Dict[str, Any]]:
    """Create ``lock`` exclusively and return the holder record written, or None if it exists."""
    holder: Dict[str, Any] = {
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "started": datetime.now(timezone.utc).isoformat(),
        "token": secrets.token_hex(8),
    }
    pidns = _pid_namespace()
    if pidns is not None:
        holder["pidns"] = pidns
    try:
        handle = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        return None
    except OSError as e:
        raise DataAccessError(f"Cannot create the lock file {lock}: {e}") from e
    try:
        os.write(handle, json.dumps(holder).encode())
    except OSError as e:
        os.close(handle)
        Path(lock).unlink(missing_ok=True)
        raise DataAccessError(f"Cannot write the lock file {lock}: {e}") from e
    os.close(handle)
    return holder


# ----------------------------------------------------------------------- reclaim


@dataclass(frozen=True)
class _Observation:
    """One reading of a lock file: its holder, inode and age in seconds."""

    holder: Dict[str, Any]
    inode: int
    age: float


def _observe(lock: Path) -> Optional[_Observation]:
    """Read ``lock``'s holder, inode and age, or None when the file no longer exists."""
    try:
        status = os.stat(lock)
    except FileNotFoundError:
        return None
    return _Observation(read_holder(lock), status.st_ino, time.time() - status.st_mtime)


def _judged_stale(age: float, policy: LockPolicy) -> bool:
    """Whether a lock whose heartbeat is ``age`` seconds old has outlived ``policy``."""
    return age > policy.stale_seconds


def _take_guard(guard: Path) -> bool:
    """Create the reclaim guard exclusively; remove one left by a crash and try once more."""
    for _ in range(2):
        try:
            os.close(os.open(guard, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644))
            return True
        except FileExistsError:
            try:
                age = time.time() - guard.stat().st_mtime
            except FileNotFoundError:
                continue
            if age <= RECLAIM_GUARD_STALE_SECONDS:
                return False
            logger.warning("Removing the reclaim guard %s (%.0f s old)", guard, age)
            guard.unlink(missing_ok=True)
        except OSError as e:
            raise DataAccessError(f"Cannot create the reclaim guard {guard}: {e}") from e
    return False


def _matches(current: Optional[_Observation], observed: _Observation) -> bool:
    """Whether ``current`` is the same lock file, with the same holder, as ``observed``."""
    return current is not None and current.inode == observed.inode and current.holder == observed.holder


def _restore(aside: Path, lock: Path) -> None:
    """Put a lock moved aside by mistake back under its name, unless a new lock took the name."""
    try:
        os.link(aside, lock)
    except FileExistsError:
        logger.warning(
            "Cannot restore %s: a new lock was created meanwhile; its earlier holder will find it lost", lock
        )
    except OSError as e:
        logger.warning("Cannot restore %s: %s; its holder will find it lost", lock, e)


def _reclaim(lock: Path, observed: _Observation, policy: LockPolicy) -> bool:
    """Remove ``lock`` if it is still the file ``observed`` and is stale or has a dead holder.

    Runs under ``<lock>.reclaim``, so only one waiter decides at a time, and re-reads the
    lock under it: a lock replaced since the observation (another waiter reclaimed it and
    took it) is left alone. The lock is then moved aside with ``rename`` and the moved file
    compared once more, because a guard removed as crashed by two waiters at once can admit
    both; a file that turns out not to be the observed one is put back with ``os.link``.
    Returns True when this call removed the lock.
    """
    guard = lock.with_name(lock.name + ".reclaim")
    if not _take_guard(guard):
        return False
    try:
        current = _observe(lock)
        if current is None or not _matches(current, observed):
            return False
        dead = holder_is_dead(current.holder)
        if not dead and not _judged_stale(current.age, policy):
            return False
        aside = lock.with_name(f"{lock.name}.reclaimed.{secrets.token_hex(4)}")
        try:
            os.rename(lock, aside)
        except FileNotFoundError:
            return False
        except OSError as e:
            raise DataAccessError(f"Cannot move the stale lock {lock} aside: {e}") from e
        try:
            if not _matches(_observe(aside), observed):
                _restore(aside, lock)
                return False
        finally:
            aside.unlink(missing_ok=True)
        if dead:
            logger.warning(
                "Took over the lock on %s: its holder is no longer running (%s)",
                policy.what,
                describe_holder(current.holder),
            )
        else:
            logger.warning(
                "Removed the stale lock on %s (no heartbeat for %.0f s): held by %s",
                policy.what,
                current.age,
                describe_holder(current.holder),
            )
        return True
    finally:
        guard.unlink(missing_ok=True)


def _wait_for_guard(guard: Path) -> bool:
    """Take ``guard`` for a release, waiting up to ``RELEASE_GUARD_WAIT_SECONDS``; False if not taken."""
    deadline = time.monotonic() + RELEASE_GUARD_WAIT_SECONDS
    while not _take_guard(guard):
        if time.monotonic() >= deadline:
            logger.warning("Releasing without the reclaim guard %s, which stays taken", guard)
            return False
        time.sleep(0.005)
    return True


# --------------------------------------------------------------------- heartbeat


@dataclass
class _Held:
    """One acquisition held by this process: its file, resolved path, token and heartbeat state."""

    path: Path
    key: str
    token: str
    interval: float
    refreshed: float = field(default_factory=time.monotonic)
    lost: bool = False
    warned: bool = False


def _read_token(handle: int) -> Optional[str]:
    """The token in the holder record open on ``handle``, or None when it has none."""
    chunks = []
    while True:
        chunk = os.read(handle, 65536)
        if not chunk:
            break
        chunks.append(chunk)
    try:
        data = json.loads(b"".join(chunks))
    except ValueError:
        return None
    token = data.get("token") if isinstance(data, dict) else None
    return str(token) if token is not None else None


def _mark_lost(entry: _Held, reason: str) -> None:
    if not entry.lost:
        logger.warning("Lost the lock %s: %s", entry.path, reason)
    entry.lost = True


def _warn_once(entry: _Held, error: OSError) -> None:
    """Log a refresh failure the first time it happens for ``entry``; later ones go to debug."""
    if entry.warned:
        logger.debug("Cannot refresh the lock %s: %s", entry.path, error)
        return
    entry.warned = True
    logger.warning("Cannot refresh the lock %s: %s (further failures are logged at debug level)", entry.path, error)


def _refresh(entry: _Held) -> None:
    """Touch ``entry``'s lock file, or mark it lost when it is gone or holds another token.

    The token is read, and the mtime set, through one open descriptor, so a lock replaced
    between the two steps is never refreshed on the new holder's behalf.
    """
    entry.refreshed = time.monotonic()
    try:
        handle = os.open(entry.path, os.O_RDONLY)
    except FileNotFoundError:
        _mark_lost(entry, "it was removed")
        return
    except OSError as e:
        _warn_once(entry, e)
        return
    try:
        if _read_token(handle) != entry.token:
            _mark_lost(entry, "it was taken over")
            return
        os.utime(handle if os.utime in os.supports_fd else entry.path, None)
    except OSError as e:
        _warn_once(entry, e)
    finally:
        os.close(handle)


class _Heartbeat:
    """The one daemon thread per process that refreshes every held lock's mtime.

    Entries are keyed by token, not path, so two acquisitions of one path in this process
    (a holder whose lock was taken over by a sibling thread, and that sibling) never share
    or remove each other's entry.
    """

    def __init__(self) -> None:
        self.cond = threading.Condition()
        self.held: Dict[str, _Held] = {}
        self.thread: Optional[threading.Thread] = None

    def _ensure_thread(self) -> None:
        """Start the thread if it is not running; the caller holds ``cond``."""
        if self.thread is None or not self.thread.is_alive():
            self.thread = threading.Thread(target=self._run, name="prepare-lock-heartbeat", daemon=True)
            self.thread.start()

    def add(self, entry: _Held) -> None:
        with self.cond:
            self.held[entry.token] = entry
            self._ensure_thread()
            self.cond.notify_all()

    def remove(self, entry: _Held) -> None:
        with self.cond:
            if self.held.get(entry.token) is entry:
                del self.held[entry.token]

    def ensure_running(self) -> None:
        with self.cond:
            if self.held:
                self._ensure_thread()

    def entries_for(self, key: str) -> List[_Held]:
        with self.cond:
            return [entry for entry in self.held.values() if entry.key == key]

    def held_count(self) -> int:
        with self.cond:
            return len(self.held)

    def _due(self) -> Tuple[List[_Held], Optional[float]]:
        """Entries due for a refresh, and how long to sleep when none is (None: no lock held)."""
        now = time.monotonic()
        live = [entry for entry in self.held.values() if not entry.lost]
        due = [entry for entry in live if now - entry.refreshed >= entry.interval]
        if due or not live:
            return due, None
        return due, min(entry.refreshed + entry.interval for entry in live) - now

    def _run(self) -> None:
        # Python runs signal handlers on the main thread only; blocking them here keeps the kernel
        # from delivering one to this thread while the main thread holds them back (_DeferredSignals).
        if hasattr(signal, "pthread_sigmask"):
            signal.pthread_sigmask(signal.SIG_BLOCK, _DEFERRED_SIGNALS)
        while True:
            with self.cond:
                due, sleep = self._due()
                if not due:
                    self.cond.wait(timeout=sleep)
                    continue
            # File I/O happens outside the condition, so a slow network file system never
            # blocks another thread taking or releasing a lock.
            for entry in due:
                try:
                    _refresh(entry)
                except _REFRESH_FAILURES as e:
                    # One broken entry must not stop the refresh of every other held lock.
                    entry.lost = True
                    logger.warning("Stopped refreshing the lock %s after an unexpected error: %s", entry.path, e)


_HEARTBEAT = _Heartbeat()
_THREAD_STATE = threading.local()


def _thread_entries() -> Dict[str, _Held]:
    """The calling thread's held acquisitions, keyed by the lock's resolved path."""
    held: Optional[Dict[str, _Held]] = getattr(_THREAD_STATE, "held", None)
    if held is None:
        held = {}
        _THREAD_STATE.held = held
    return held


def _after_fork_in_child() -> None:
    """A forked child holds none of its parent's locks and has no heartbeat thread."""
    global _HEARTBEAT, _THREAD_STATE
    _HEARTBEAT = _Heartbeat()
    _THREAD_STATE = threading.local()


if hasattr(os, "register_at_fork"):
    os.register_at_fork(after_in_child=_after_fork_in_child)


# -------------------------------------------------------------- acquire, release


class _DeferredSignals:
    """SIGINT, SIGTERM and SIGHUP blocked on the main thread between creating a lock and recording it.

    ``block`` is called before each attempt to create the lock file and ``restore`` once the
    acquisition is recorded as held (or the attempt failed); a signal that arrives in between
    is delivered by ``restore``. Off the main thread, and where ``pthread_sigmask`` does not
    exist (Windows), both do nothing.
    """

    def __init__(self) -> None:
        self.previous: Optional[Any] = None

    def block(self) -> None:
        if self.previous is not None or not hasattr(signal, "pthread_sigmask"):
            return
        if threading.current_thread() is threading.main_thread():
            self.previous = signal.pthread_sigmask(signal.SIG_BLOCK, _DEFERRED_SIGNALS)

    def restore(self) -> None:
        if self.previous is not None:
            previous, self.previous = self.previous, None
            signal.pthread_sigmask(signal.SIG_SETMASK, previous)


def _held_error(lock: Path, policy: LockPolicy, timed_out: bool = False) -> LockHeld:
    message = f"{policy.what} is locked by {describe_holder(read_holder(lock))}: {lock}"
    return LockWaitTimeout(message) if timed_out else LockHeld(message)


def _check_stop(policy: LockPolicy, should_stop: Optional[Callable[[], bool]], waited: float) -> None:
    if should_stop is not None and should_stop():
        raise LockWaitStopped(f"Stopped waiting for {policy.what} after {waited:.0f} s: a stop was requested")


def _acquire(
    lock: Path,
    key: str,
    policy: LockPolicy,
    should_stop: Optional[Callable[[], bool]],
    blocking: bool,
    deferred: _DeferredSignals,
) -> Dict[str, Any]:
    """Take ``lock`` under ``policy`` and return the holder record written.

    The wait limit is checked before a reclaim, so a lock truly held by another process is
    reported, not taken, when both thresholds pass on the same poll. Signals are deferred
    around each creation attempt and stay deferred when it succeeds: the caller restores them
    once it has recorded the lock as held.
    """
    started = time.monotonic()
    next_log = FIRST_WAIT_LOG_SECONDS
    while True:
        deferred.block()
        holder = _create(lock)
        if holder is not None:
            return holder
        deferred.restore()
        waited = time.monotonic() - started
        if key in _thread_entries():
            # A stop the caller asked for is what it gets; otherwise the wait could never end.
            _check_stop(policy, should_stop, waited)
            raise LockReentry(f"{policy.what} is already locked by this thread: {lock}")
        if blocking and policy.wait_seconds and waited >= policy.wait_seconds:
            raise _held_error(lock, policy, timed_out=True)
        observed = _observe(lock)
        if observed is None:
            continue
        if _judged_stale(observed.age, policy) or holder_is_dead(observed.holder):
            if _reclaim(lock, observed, policy):
                continue
        if not blocking:
            raise _held_error(lock, policy)
        _check_stop(policy, should_stop, waited)
        if waited >= next_log:
            logger.info("waiting for %s: held by %s", policy.what, describe_holder(observed.holder))
            next_log = waited + LOCK_WAIT_LOG_SECONDS
        time.sleep(policy.poll_seconds)


def _release(lock: Path, entry: _Held) -> None:
    """Remove ``lock`` only while it still carries this acquisition's token.

    The token is compared under the reclaim guard, the only place a lock is otherwise
    removed, so a lock reclaimed and taken by someone else after the comparison cannot be
    removed by it.
    """
    guard = lock.with_name(lock.name + ".reclaim")
    guarded = _wait_for_guard(guard)
    try:
        current = read_holder(lock)
        if current.get("token") == entry.token:
            lock.unlink(missing_ok=True)
            return
        if lock.exists():
            logger.warning("Not removing %s: it is now held by %s", lock, describe_holder(current))
        else:
            logger.warning("The lock %s was removed while it was held", lock)
    finally:
        if guarded:
            guard.unlink(missing_ok=True)


@contextmanager
def held_lock(
    lock: Path,
    policy: LockPolicy,
    should_stop: Optional[Callable[[], bool]] = None,
    blocking: bool = True,
) -> Iterator[Path]:
    """Hold ``lock`` for the length of the block, under ``policy``.

    Waits while a live holder keeps the lock, up to ``policy.wait_seconds`` (zero: no
    limit), then raises ``LockWaitTimeout`` (a ``LockHeld``); ``blocking=False`` raises a plain
    ``LockHeld`` at once instead, after any stale or dead-holder takeover. ``should_stop``
    ends a wait with ``LockWaitStopped`` as soon as it returns True; a free lock is taken
    regardless. A thread asking for a lock it already holds (by its resolved path) gets
    ``LockReentry`` at once; other threads contend normally. The parent folder must exist.
    """
    lock = Path(lock)
    key = os.path.realpath(lock)
    owner = os.getpid()
    deferred = _DeferredSignals()
    entry: Optional[_Held] = None
    try:
        holder = _acquire(lock, key, policy, should_stop, blocking, deferred)
        entry = _Held(lock, key, str(holder["token"]), policy.heartbeat_seconds)
        _thread_entries()[key] = entry
        _HEARTBEAT.add(entry)
        deferred.restore()
        yield lock
    finally:
        deferred.restore()
        # A forked child leaves the block too, but the lock is its parent's to release.
        if entry is not None and os.getpid() == owner:
            if _thread_entries().get(key) is entry:
                del _thread_entries()[key]
            _HEARTBEAT.remove(entry)
            _release(lock, entry)


def verify_held(lock: Path) -> None:
    """Raise ``LockLost`` unless this process still holds ``lock`` under its own token.

    The calling thread's own acquisition is checked when it has one, so a holder whose lock
    a sibling thread took over is told so; otherwise any acquisition by this process counts.
    Restarts the heartbeat thread if it has died.
    """
    key = os.path.realpath(lock)
    _HEARTBEAT.ensure_running()
    own = _thread_entries().get(key)
    candidates = [own] if own is not None else _HEARTBEAT.entries_for(key)
    if not candidates:
        raise LockLost(f"{lock} is not held by this process")
    current = read_holder(Path(lock))
    for candidate in candidates:
        if candidate.token != current.get("token"):
            candidate.lost = True
        elif not candidate.lost:
            return
    detail = f"now held by {describe_holder(current)}" if current else "removed"
    raise LockLost(f"Lost the lock {lock}: {detail}")
