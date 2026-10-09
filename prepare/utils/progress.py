"""
Progress summaries for long runs over many items.

A run over thousands of accessions that logged one INFO line per item would bury its warnings;
``ProgressReporter`` logs one summary line every ``every`` items, or every ``min_interval``
seconds when items are slow, and one line when the run finishes::

    download_sra: 150/2000 done (148 ok, 2 failed), 3.1/min, about 9 h 57 min left
    download_sra: finished 2000/2000 (1990 ok, 10 failed) in 10 h 45 min

The caller logs its own line per item at ``item_level``: DEBUG while summaries are on, INFO when
``every`` is 0, which turns the summaries off and restores one INFO line per item.
"""

import logging
import math
import threading
import time
from typing import Callable, Optional

_module_logger = logging.getLogger(__name__)


def item_level(every: int) -> int:
    """The level for a line about one item: INFO when summaries are off (``every`` 0), else DEBUG."""
    return logging.INFO if every == 0 else logging.DEBUG


def active_item_level() -> int:
    """``item_level`` for the active run's ``progress_every`` setting (``--progress-every``)."""
    from prepare.core import settings

    return item_level(settings.active().progress_every)


class DemoteInfo(logging.Filter):
    """Log the INFO records whose message template is in ``templates`` at ``level`` instead.

    Attached to a module's logger for the length of a run, it moves that module's per-item INFO
    lines to ``item_level`` without editing the module (``data/read_extraction.py`` is held at a
    line ceiling).
    """

    def __init__(self, templates, level: int) -> None:
        """Remember the message templates to match and the level to give them."""
        super().__init__()
        self.templates = frozenset(templates)
        self.level = level

    def filter(self, record: logging.LogRecord) -> bool:
        """Lower a matching INFO record to ``level``; every record is kept for the handlers to judge."""
        if record.levelno == logging.INFO and record.msg in self.templates:
            record.levelno = self.level
            record.levelname = logging.getLevelName(self.level)
        return True


def format_duration(seconds: float) -> str:
    """``seconds`` rounded to the minute, as ``N min`` or ``H h M min`` (``less than 1 min`` below 60 s)."""
    if seconds < 60:
        return "less than 1 min"
    minutes = int(round(seconds / 60))
    if minutes < 60:
        return f"{minutes} min"
    return f"{minutes // 60} h {minutes % 60} min"


class ProgressReporter:
    """Count finished items from any thread and log a summary line now and then.

    Args:
        label: Start of every line, normally the command name
        total: Number of items the run will process
        every: Items between summary lines; 0 turns summaries off (the caller then logs every
            item at INFO, see ``item_level``)
        min_interval: Seconds after which a summary is logged even before the next multiple of
            ``every`` is reached, so a slow run still shows it is alive
        logger: Logger the lines go to (this module's when None)
        clock: Monotonic clock in seconds, replaced by a fake one in tests
    """

    def __init__(
        self,
        label: str,
        total: int,
        every: int,
        min_interval: float = 300.0,
        logger: Optional[logging.Logger] = None,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        """Start counting at the current time of ``clock``."""
        if every < 0:
            raise ValueError(f"every must be 0 or more, got {every}")
        self.label = label
        self.total = total
        self.every = every
        self.min_interval = min_interval
        self.logger = logger or _module_logger
        self.clock = clock
        self.ok = 0
        self.failed = 0
        self._lock = threading.Lock()
        self._started = clock()
        self._last_line_at = self._started
        self._last_line_done = 0
        self._finished = False

    @property
    def done(self) -> int:
        """Items finished so far, successful or not."""
        return self.ok + self.failed

    @property
    def item_level(self) -> int:
        """The level for a caller's line about one item: INFO when summaries are off, else DEBUG."""
        return item_level(self.every)

    def update(self, ok: bool, n: int = 1) -> None:
        """Count ``n`` finished items as successful (``ok``) or failed; safe to call from any thread."""
        self.update_counts(n if ok else 0, 0 if ok else n)

    def update_counts(self, ok: int, failed: int) -> None:
        """Count a batch of ``ok`` successful and ``failed`` failed items as one update (at most one line)."""
        with self._lock:
            self.ok += ok
            self.failed += failed
            if self._due():
                self._last_line_at = self.clock()
                self._last_line_done = self.done
                self._log(self._summary())

    def finish(self, ok: Optional[int] = None, failed: Optional[int] = None) -> None:
        """Log the closing line (once, and only when at least one item was counted).

        ``ok`` and ``failed``, when given, replace the counts first: the final tally after a
        retry pass that turned some failures into successes.
        """
        with self._lock:
            if ok is not None and failed is not None:
                self.ok, self.failed = ok, failed
            if self._finished or self.done == 0:
                return
            self._finished = True
            elapsed = self.clock() - self._started
            self._log(
                f"{self.label}: finished {self.done}/{self.total} ({self.ok} ok, {self.failed} failed) "
                f"in {format_duration(elapsed)}"
            )

    def _due(self) -> bool:
        """True when a summary line is owed: a multiple of ``every`` was passed, or ``min_interval`` went by."""
        if self.every == 0 or self.done >= self.total or self.done == self._last_line_done:
            # The last item is reported by finish(), so it is not logged twice.
            return False
        if self.done // self.every > self._last_line_done // self.every:
            return True
        return self.clock() - self._last_line_at >= self.min_interval

    def _summary(self) -> str:
        """The progress line for the counts as they are now."""
        line = f"{self.label}: {self.done}/{self.total} done ({self.ok} ok, {self.failed} failed)"
        elapsed_minutes = (self.clock() - self._started) / 60
        if elapsed_minutes <= 0:
            return line
        rate = self.done / elapsed_minutes
        line += f", {rate:.1f}/min"
        remaining = self.total - self.done
        if remaining > 0 and rate > 0 and math.isfinite(remaining / rate):
            line += f", about {format_duration(remaining / rate * 60)} left"
        return line

    def _log(self, message: str) -> None:
        """Write one line at INFO."""
        self.logger.info(message)
