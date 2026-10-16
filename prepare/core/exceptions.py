"""
Custom exception classes for PREPARE, and the exit code each one gives a command.

Exit codes: 0 success, 1 failure, 2 usage error, 3 configuration (the environment lacks
something a command needs), 4 a retryable failure (network, a lock wait that gave up), 130
interrupted. A batch script can resubmit on 4 and stop on 1 or 3.
"""

from enum import IntEnum


class ExitCode(IntEnum):
    """Process exit codes returned by ``prepare``."""

    OK = 0
    FAILURE = 1
    USAGE = 2
    CONFIGURATION = 3
    TRANSIENT = 4
    INTERRUPTED = 130


class PrepareError(Exception):
    """Base exception for all PREPARE errors."""

    exit_code: int = ExitCode.FAILURE


class ValidationError(PrepareError):
    """Exception raised for validation errors in input data; a usage error (exit code 2)."""

    exit_code = ExitCode.USAGE


class FormatError(ValidationError):
    """Exception raised for errors related to file formats."""


class DataAccessError(PrepareError):
    """Exception raised for errors in data access operations."""


class TransientError(DataAccessError):
    """A data access failure that may succeed if the same command is run again later."""

    exit_code = ExitCode.TRANSIENT


class LockTimeoutError(TransientError):
    """Raised when a wait for a lock held by another live process reaches its time limit."""


class NetworkError(TransientError):
    """Raised when a remote service cannot be reached, times out, or answers 429 or 5xx."""


class ProcessingError(PrepareError):
    """Exception raised for errors during data processing."""


class VisualizationError(PrepareError):
    """Exception raised for errors during visualization generation."""


class PluginError(PrepareError):
    """Exception raised for errors related to plugins."""


class SecurityError(PrepareError):
    """Exception raised for security-related errors."""


class ConfigurationError(PrepareError):
    """Exception raised when the environment lacks something a command needs, such as an optional package."""

    exit_code = ExitCode.CONFIGURATION


def exit_code_for(error: BaseException) -> int:
    """The process exit code for ``error``: its class's code, 130 for an interrupt, else 1."""
    if isinstance(error, KeyboardInterrupt):
        return int(ExitCode.INTERRUPTED)
    if isinstance(error, PrepareError):
        return int(error.exit_code)
    return int(ExitCode.FAILURE)
