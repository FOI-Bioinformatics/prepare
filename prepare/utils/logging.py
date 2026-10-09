"""
Logging utilities for PREPARE.

``setup_logging`` installs a console handler on stderr and, when a log file is given, a file
handler that appends to it. The file receives every line at INFO or above (DEBUG too when the
console is at DEBUG), each stamped with the host name and process ID, and the full traceback
of a failure; the console shows a traceback only when asked to (at DEBUG).
"""

import logging
import socket
import sys
from pathlib import Path
from typing import Optional

CONSOLE_FORMAT = "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
CONSOLE_FORMAT_WITH_HOST = "%(asctime)s - %(hostname)s[%(process)d] - %(levelname)s - %(name)s - %(message)s"
FILE_FORMAT = "%(asctime)s %(hostname)s[%(process)d] %(levelname)s %(name)s: %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

# Looked up once: the host does not change during a run, and a lookup per record would cost a
# system call on every log line.
_HOSTNAME = socket.gethostname()


class HostFilter(logging.Filter):
    """Set ``record.hostname`` so a format string can name the host that wrote the line."""

    def filter(self, record: logging.LogRecord) -> bool:
        """Add the host name to ``record``; never drops a record."""
        record.hostname = _HOSTNAME
        return True


class ConsoleFormatter(logging.Formatter):
    """A formatter that leaves out the traceback and stack of a record unless ``show_traceback``.

    The file handler keeps them. ``logging.Formatter.format`` caches the traceback text on the
    record, so another handler that formatted the record first would otherwise put it on the
    console too; the record's exception fields are hidden for the duration of this call only
    and then restored.
    """

    def __init__(self, fmt: Optional[str] = None, datefmt: Optional[str] = None, show_traceback: bool = False):
        """Format with ``fmt`` and ``datefmt``; ``show_traceback`` keeps tracebacks in the output."""
        super().__init__(fmt, datefmt)
        self.show_traceback = show_traceback

    def format(self, record: logging.LogRecord) -> str:
        """The formatted line, without the traceback unless ``show_traceback``."""
        if self.show_traceback:
            return super().format(record)
        saved = (record.exc_info, record.exc_text, record.stack_info)
        record.exc_info, record.exc_text, record.stack_info = None, None, None
        try:
            return super().format(record)
        finally:
            record.exc_info, record.exc_text, record.stack_info = saved


def _tagged(handler: logging.Handler) -> logging.Handler:
    """Mark ``handler`` as installed by ``setup_logging`` so a later call can remove it."""
    handler._prepare = True  # type: ignore[attr-defined]
    return handler


def _remove_own_handlers(root_logger: logging.Logger) -> None:
    """Remove (and close) only the handlers an earlier ``setup_logging`` call installed."""
    for handler in root_logger.handlers[:]:
        if getattr(handler, "_prepare", False):
            root_logger.removeHandler(handler)
            handler.close()


def log_traceback_hint() -> None:
    """After an error: point to ``--log-level DEBUG`` unless a log file keeps the traceback or DEBUG shows it."""
    root_logger = logging.getLogger()
    if any(isinstance(h, logging.FileHandler) and getattr(h, "_prepare", False) for h in root_logger.handlers):
        return
    if not root_logger.isEnabledFor(logging.DEBUG):
        logging.info("Use --log-level DEBUG for full traceback.")


def _file_handler(log_file: str, level: int) -> logging.Handler:
    """A handler appending to ``log_file`` (its folder created) at ``level`` or INFO, whichever is lower."""
    path = Path(log_file).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(path, mode="a", encoding="utf-8")
    handler.setLevel(min(level, logging.INFO))
    handler.setFormatter(logging.Formatter(FILE_FORMAT, DATE_FORMAT))
    handler.addFilter(HostFilter())
    return _tagged(handler)


def setup_logging(
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    *,
    show_host: bool = False,
    console_traceback: bool = False,
) -> None:
    """
    Set up logging to stderr and, optionally, to a file.

    Calling this twice replaces only the handlers it installed itself, so a host
    application's own handlers (or pytest's caplog handler) are left in place.

    Args:
        level: Console logging level
        log_file: File to append every log line to, at INFO or ``level`` whichever is lower,
            with host name, process ID and full tracebacks; its folder is created. None logs to
            stderr only
        show_host: Put the host name and process ID on every console line as well
        console_traceback: Show tracebacks on the console (the file always has them)

    Raises:
        OSError: If the log file cannot be opened for appending.
    """
    root_logger = logging.getLogger()
    _remove_own_handlers(root_logger)

    # The file handler is built first, so a log file that cannot be opened leaves no half set up
    # console handler behind; the caller then sets up the console alone to report it.
    file_handler = _file_handler(log_file, level) if log_file else None

    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_format = CONSOLE_FORMAT_WITH_HOST if show_host else CONSOLE_FORMAT
    console_handler.setFormatter(ConsoleFormatter(console_format, DATE_FORMAT, show_traceback=console_traceback))
    console_handler.addFilter(HostFilter())
    root_logger.addHandler(_tagged(console_handler))

    if file_handler is not None:
        root_logger.addHandler(file_handler)
        root_logger.setLevel(min(level, logging.INFO))
    else:
        root_logger.setLevel(level)

    # Suppress verbose logging from some libraries
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)

    logging.debug("Logging configured")
