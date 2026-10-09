"""Graceful handling of SIGINT, SIGTERM and SIGHUP for the length of a command.

A batch job's time limit, ``kill``, a closed terminal or Ctrl-C would otherwise end the process
at an arbitrary point, or end it again while it is writing what it has done so far.
``graceful_termination`` turns the first of these signals into ``KeyboardInterrupt``, so
``finally`` blocks and ``with`` exits run (a registry batch is flushed, running tools are
stopped). A later signal is logged and counted but not raised, so it cannot cut the final
write short. After ``abandon_after`` signals in total the process stops waiting: running tools
are killed and the process exits with code 130 at once. That is safe for the lock files, whose
holders are recognised as dead by their heartbeat and process id, and for the registry and
store files, which are only ever replaced whole.
"""

import logging
import os
import signal
import sys
import threading
from contextlib import contextmanager
from dataclasses import dataclass, field
from types import FrameType
from typing import Iterator, List, Optional

from prepare.utils.security import SecureSubprocess

logger = logging.getLogger(__name__)

# Signals in total (the first included) after which the process exits without finishing its write.
ABANDON_AFTER_SIGNALS = 3

# Exit code of a command stopped by a signal or Ctrl-C.
EXIT_INTERRUPTED = 130

# The Termination of the context active on the main thread, shared by nested contexts.
_active: Optional["Termination"] = None


@dataclass
class Termination:
    """What a ``graceful_termination`` context has received so far.

    ``stop`` is set by the first signal, ``signum`` records which one it was and ``repeats``
    counts the signals received after it.
    """

    stop: threading.Event = field(default_factory=threading.Event)
    signum: Optional[int] = None
    repeats: int = 0
    _linked: List[threading.Event] = field(default_factory=list, init=False, repr=False, compare=False)

    @property
    def requested(self) -> bool:
        """True once a signal was received or ``stop`` was set."""
        return self.signum is not None or self.stop.is_set()

    @property
    def cause(self) -> str:
        """The first signal's name for a log line, or ``Ctrl-C`` for SIGINT or no recorded signal."""
        if self.signum is None or self.signum == signal.SIGINT:
            return "Ctrl-C"
        return signal_name(self.signum)

    def _set_stop(self) -> None:
        self.stop.set()
        for event in self._linked:
            event.set()


def signal_name(signum: int) -> str:
    """Return the name of signal ``signum`` (``SIGTERM``), or its number for an unknown one."""
    try:
        return signal.Signals(signum).name
    except ValueError:
        return f"signal {signum}"


def handled_signals() -> List[int]:
    """The signals ``graceful_termination`` handles on this platform."""
    names = ("SIGINT", "SIGBREAK") if sys.platform == "win32" else ("SIGINT", "SIGTERM", "SIGHUP")
    return [getattr(signal, name) for name in names if hasattr(signal, name)]


def _make_handler(term: Termination, abandon_after: int):
    """Return the signal handler that implements the first/repeat/abandon behaviour for ``term``."""

    def _handler(signum: int, _frame: Optional[FrameType]) -> None:
        name = signal_name(signum)
        if term.signum is None:
            term.signum = signum
            term._set_stop()
            raise KeyboardInterrupt(f"received {name}")
        term.repeats += 1
        remaining = abandon_after - 1 - term.repeats
        if remaining <= 0:
            logger.error("Received %s %d times; abandoning the write and exiting now", name, term.repeats + 1)
            SecureSubprocess.terminate_children(grace=0)
            os._exit(EXIT_INTERRUPTED)
        else:
            logger.warning("Received %s again; %d more will abandon the write", name, remaining)

    return _handler


@contextmanager
def graceful_termination(
    stop: Optional[threading.Event] = None, abandon_after: int = ABANDON_AFTER_SIGNALS
) -> Iterator[Termination]:
    """Handle SIGINT, SIGTERM and SIGHUP (SIGINT and SIGBREAK on Windows) for the ``with`` block.

    The first signal sets ``stop`` (a new event when None is given), records the signal on the
    yielded ``Termination`` and raises ``KeyboardInterrupt``. Each later one is logged with how
    many more will abandon the write, and not raised. The ``abandon_after``-th signal kills
    running tools and exits the process with code 130 without further cleanup. A signal that
    is ignored when the block starts (SIGHUP under ``nohup``) stays ignored. The previous
    handlers are restored when the block exits, whether it finished or raised.

    Signal handlers can only be installed from the main thread; elsewhere the context installs
    nothing and yields a ``Termination`` that no signal will change. A context opened inside
    another one on the main thread installs nothing either and yields the outer ``Termination``;
    its ``stop``, when given, is set together with the outer one.
    """
    global _active
    if threading.current_thread() is not threading.main_thread():
        yield Termination(stop if stop is not None else threading.Event())
        return
    if _active is not None:
        outer = _active
        if stop is not None and stop is not outer.stop:
            outer._linked.append(stop)
            if outer.requested:
                stop.set()
        try:
            yield outer
        finally:
            if stop is not None and stop in outer._linked:
                outer._linked.remove(stop)
        return

    term = Termination(stop if stop is not None else threading.Event())
    handler = _make_handler(term, abandon_after)
    previous = {}
    try:
        for signum in handled_signals():
            old = signal.getsignal(signum)
            if old is signal.SIG_IGN:
                continue
            # Stored before the new handler goes in: a signal arriving as signal.signal returns
            # runs the handler at once, and the old one must already be known to be restored.
            previous[signum] = old
            signal.signal(signum, handler)
        _active = term
        yield term
    finally:
        _active = None
        for signum, old in previous.items():
            # None means the previous handler was not installed from Python; the default is the best match.
            signal.signal(signum, old if old is not None else signal.SIG_DFL)
