"""
Security utilities for PREPARE.

This module provides secure subprocess handling and input validation.
"""

import logging
import os
import re
import signal
import subprocess
import tempfile
import threading
import time
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Tuple, Union

from prepare.core.exceptions import SecurityError
from prepare.core.validation import validate_accession
from prepare.core.constants import (
    ALLOWED_BIOINFORMATICS_TOOLS,
    DANGEROUS_ENV_VARS,
    SRA_ACCESSION_PATTERN,
    UNSAFE_SHELL_CHARS,
)
from prepare.core import settings as settings_module

logger = logging.getLogger(__name__)


class _TimeoutNotGiven:
    """Sentinel for ``run_secure``'s ``timeout``: not given means use the active setting."""

    def __repr__(self) -> str:
        return "<not given>"


# A run_secure call that does not pass timeout gets the active run's subprocess_timeout
# setting (--timeout, PREPARE_TIMEOUT, config [runtime] timeout), resolved at call time
# rather than at import time, so a test or a later command sees its own settings.
_TIMEOUT_NOT_GIVEN = _TimeoutNotGiven()

# Per-tool flags that never take a value; any other allowlisted flag for that
# tool consumes the following token as its value.
BOOLEAN_FLAGS: Dict[str, FrozenSet[str]] = {
    "fasterq-dump": frozenset(
        {"--progress", "--split-files", "--split-3", "--skip-technical", "--include-technical", "--force", "--version"}
    ),
    "prefetch": frozenset({"--progress", "--resume", "--version"}),
    "pigz": frozenset({"-f", "-k", "--version"}),
    "minimap2": frozenset({"-a", "--sam-hit-only", "--version"}),
    "samtools": frozenset({"-b", "-c", "--version"}),
    "megahit": frozenset({"--no-mercy", "--version"}),
    "seqkit": frozenset({"-T", "--version"}),
}
# Kept for backward compatibility with any caller importing the old name.
FASTERQ_DUMP_BOOLEAN_FLAGS = BOOLEAN_FLAGS["fasterq-dump"]
# Per-tool flags whose value must be a positive integer (a thread or process count).
POSITIVE_INTEGER_FLAGS = {"fasterq-dump": frozenset({"--threads"}), "pigz": frozenset({"-p"})}
# Kept for callers importing the old name.
FASTERQ_DUMP_INTEGER_FLAGS = POSITIVE_INTEGER_FLAGS["fasterq-dump"]
# Flags (any tool) whose value is a filesystem path and must pass validate_path.
PATH_VALUE_FLAGS = frozenset({"-O", "-o", "-d", "--out-dir", "--temp", "--tmp-dir", "-1", "-2", "-0", "-s"})
# Tools whose positional argument is either an SRA accession or a .sra file path.
SRA_POSITIONAL_TOOLS = frozenset({"fasterq-dump", "prefetch"})
# Return code run_secure reports for a child it did not start because the run is stopping: the
# code of a child killed with SIGKILL, which is what a child started a moment later would get.
_NOT_STARTED_RETURNCODE = -9
# On POSIX each tool runs in its own session, so its process group holds the tool and every
# process it starts (the conda megahit wrapper runs megahit_core as a child that inherits the
# output pipes). A timeout or a termination signals the whole group; otherwise the wrapper
# would die while its child kept the pipes open and ran on.
_PROCESS_GROUPS = os.name == "posix"


def _signal_tool(proc: subprocess.Popen, kill: bool) -> None:
    """Send SIGTERM (SIGKILL with ``kill``) to ``proc``'s process group, or to ``proc`` alone.

    The direct child is signalled instead when the group is gone or not ours to signal, when the
    platform has no process groups (Windows), or when ``proc`` is a test double without a pid.
    """
    if _PROCESS_GROUPS and isinstance(proc.pid, int):
        try:
            os.killpg(proc.pid, signal.SIGKILL if kill else signal.SIGTERM)
            return
        except (ProcessLookupError, PermissionError):
            pass
    if kill:
        proc.kill()
    else:
        proc.terminate()


def _group_running(proc: subprocess.Popen) -> bool:
    """Whether any process of ``proc``'s group is still running (False without process groups)."""
    if not (_PROCESS_GROUPS and isinstance(proc.pid, int)):
        return False
    try:
        os.killpg(proc.pid, 0)
    except (ProcessLookupError, PermissionError):
        return False
    return True


class SecureSubprocess:
    """Secure subprocess wrapper with validation and sanitization."""

    # Get allowed tools from constants
    ALLOWED_EXECUTABLES = set(ALLOWED_BIOINFORMATICS_TOOLS.keys())

    # Get safe parameters from constants
    SAFE_PARAMETERS = {tool: config["safe_params"] for tool, config in ALLOWED_BIOINFORMATICS_TOOLS.items()}

    # Directories registered at runtime from user-supplied output or temp folders.
    _extra_roots: List[Path] = []
    _extra_roots_lock = threading.Lock()

    # Child processes started by run_secure that have not yet finished, each mapped to the stop
    # token of the run that started it (None for a caller that passed none); read by
    # terminate_children so an interrupt can stop tools running in worker threads.
    _children: Dict[subprocess.Popen, Optional[threading.Event]] = {}
    _children_lock = threading.Lock()
    # Process-wide: set by terminate_children called without a token. A child started after
    # that is killed as soon as it is created, whatever its token, so a worker that passed its
    # own stop check just before the interrupt cannot leave a tool running. No run clears it;
    # only clear_stopping does. A run's own token plays the same part for that run alone.
    _stopping = False

    @classmethod
    def clear_stopping(cls) -> None:
        """Allow ``run_secure`` children to run again after ``terminate_children`` without a token."""
        with cls._children_lock:
            cls._stopping = False

    @classmethod
    def _stop_requested(cls, stop: Optional[threading.Event]) -> bool:
        """Whether a child under ``stop`` must not run: the process-wide flag or the token is set."""
        return cls._stopping or (stop is not None and stop.is_set())

    @classmethod
    def _children_snapshot(cls) -> List[Tuple[subprocess.Popen, Optional[threading.Event]]]:
        """The tracked children and their tokens; retried when another thread changes the table meanwhile."""
        for _ in range(100):
            try:
                return list(cls._children.items())
            except RuntimeError:
                continue
        return []

    @classmethod
    def terminate_children(cls, grace: float = 5.0, stop: Optional[threading.Event] = None) -> int:
        """Terminate, then kill, the children started by ``run_secure`` that are still running.

        With ``stop`` None every tracked child is stopped, and any child started after this call
        is killed at once until ``clear_stopping`` is called. With a token only the children
        ``run_secure`` started under that token are stopped; the token is set, so a child that
        run starts afterwards is killed at once (or not started), and other runs are unaffected.
        Each child's process group (the child alone on Windows) is sent SIGTERM; a group with a
        process still running ``grace`` seconds later is sent SIGKILL. Returns the number of
        children stopped.
        """
        # The third-signal handler calls this on the main thread, which may be inside run_secure
        # holding the (non-reentrant) lock: a bounded wait, then an unlocked snapshot, so the
        # handler can never deadlock against the code it interrupted.
        locked = cls._children_lock.acquire(timeout=0.1)
        try:
            if stop is None:
                cls._stopping = True
            else:
                stop.set()
            owners = cls._children_snapshot()
        finally:
            if locked:
                cls._children_lock.release()
        children = [child for child, owner in owners if stop is None or owner is stop]
        for child in children:
            _signal_tool(child, kill=False)
        deadline = time.monotonic() + grace
        for child in children:
            try:
                child.wait(timeout=max(0.0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                _signal_tool(child, kill=True)
                continue
            # The tool itself ended; a process it started may still hold the group.
            while _group_running(child) and time.monotonic() < deadline:
                time.sleep(0.05)
            if _group_running(child):
                _signal_tool(child, kill=True)
        return len(children)

    @classmethod
    def allowed_roots(cls) -> List[Path]:
        """Directories under which validated paths may fall."""
        roots = [Path.cwd(), Path.home(), Path("/tmp"), Path(tempfile.gettempdir())]
        with cls._extra_roots_lock:
            roots.extend(cls._extra_roots)
        return [root.resolve() for root in roots]

    @classmethod
    def add_allowed_root(cls, path: Union[str, Path]) -> Path:
        """Permit paths under a directory the user chose on the command line."""
        resolved = Path(path).resolve()
        with cls._extra_roots_lock:
            if resolved not in cls._extra_roots:
                cls._extra_roots.append(resolved)
        return resolved

    @staticmethod
    def validate_executable(executable: str) -> str:
        """
        Validate and sanitize executable name.

        Args:
            executable: The executable name to validate

        Returns:
            Validated executable name

        Raises:
            SecurityError: If executable is not allowed
        """
        if executable not in SecureSubprocess.ALLOWED_EXECUTABLES:
            raise SecurityError(f"Executable '{executable}' is not allowed")

        # Additional validation - ensure no path traversal
        if "/" in executable or "\\" in executable or ".." in executable:
            raise SecurityError(f"Invalid executable path: {executable}")

        return executable

    @staticmethod
    def validate_parameter(executable: str, param: str) -> str:
        """
        Validate parameter for given executable.

        Args:
            executable: The executable name
            param: The parameter to validate

        Returns:
            Validated parameter

        Raises:
            SecurityError: If parameter is not safe
        """
        if executable not in SecureSubprocess.SAFE_PARAMETERS:
            raise SecurityError(f"No parameter validation defined for {executable}")

        safe_params = SecureSubprocess.SAFE_PARAMETERS[executable]

        # Check if it's a flag parameter
        if param.startswith("-") and param not in safe_params:
            raise SecurityError(f"Parameter '{param}' not allowed for {executable}")

        # Basic sanitization - no shell metacharacters
        if any(char in param for char in UNSAFE_SHELL_CHARS):
            raise SecurityError(f"Parameter contains unsafe characters: {param}")

        return param

    @classmethod
    def validate_path(cls, path: Union[str, Path], allow_creation: bool = True) -> Path:
        """
        Validate and sanitize a file or directory path.

        A path is accepted when it lies under the working directory, the home
        directory, the system temp directory, or a root registered with
        ``add_allowed_root``. Parent-directory segments are rejected when the
        caller does not create the path.

        Raises:
            SecurityError: If the path is unsafe
        """
        path_obj = Path(path).resolve()

        if not any(path_obj == root or root in path_obj.parents for root in cls.allowed_roots()):
            raise SecurityError(f"Path outside allowed directories: {path}")

        if not allow_creation and ".." in Path(path).parts:
            raise SecurityError(f"Unsafe path component in: {path}")

        return path_obj

    @staticmethod
    def validate_accession_for_subprocess(accession: str) -> str:
        """
        Validate SRA accession for subprocess usage.

        Args:
            accession: The SRA accession to validate

        Returns:
            Validated accession

        Raises:
            SecurityError: If accession is invalid or unsafe
        """
        if not validate_accession(accession):
            raise SecurityError(f"Invalid SRA accession format: {accession}")

        # Additional security check - only alphanumeric characters allowed
        if not re.match(SRA_ACCESSION_PATTERN, accession):
            raise SecurityError(f"SRA accession contains invalid characters: {accession}")

        return accession

    @classmethod
    def _build_validated_command(cls, executable: str, args: List[str]) -> List[str]:
        """Validate the executable and each argument, returning the command list to run.

        Flags are checked against the per-tool allowlist; path-valued flags and
        SRA accessions are validated/sanitized. Raises SecurityError on any
        disallowed component.
        """
        cmd = [cls.validate_executable(executable)]

        i = 0
        while i < len(args):
            arg = args[i]

            if arg.startswith("-"):
                cls.validate_parameter(executable, arg)
                cmd.append(arg)

                takes_value = arg not in BOOLEAN_FLAGS.get(executable, frozenset())
                if takes_value and i + 1 < len(args) and not args[i + 1].startswith("-"):
                    i += 1
                    value = args[i]
                    if arg in PATH_VALUE_FLAGS:
                        value = str(cls.validate_path(value))
                    elif arg in POSITIVE_INTEGER_FLAGS.get(executable, frozenset()):
                        if not value.isdigit() or int(value) < 1:
                            raise SecurityError(f"Invalid value for {arg}: {value} (must be a positive integer)")
                    cmd.append(value)
            else:
                # Positional argument: for fasterq-dump and prefetch this is either the SRA
                # accession or a path to an already-downloaded archive. NCBI serves some runs
                # only as .sralite, which fasterq-dump reads like a .sra file.
                if executable in SRA_POSITIONAL_TOOLS:
                    if arg.endswith((".sra", ".sralite")):
                        arg = str(cls.validate_path(arg))
                    else:
                        arg = cls.validate_accession_for_subprocess(arg)
                cmd.append(arg)

            i += 1

        return cmd

    @classmethod
    def run_secure(
        cls,
        executable: str,
        args: List[str],
        cwd: Optional[Union[str, Path]] = None,
        env: Optional[Dict[str, str]] = None,
        timeout: Union[float, int, None, _TimeoutNotGiven] = _TIMEOUT_NOT_GIVEN,
        stop: Optional[threading.Event] = None,
        **kwargs,
    ) -> subprocess.CompletedProcess:
        """
        Run subprocess with security validations.

        Args:
            executable: The executable to run
            args: List of arguments
            cwd: Working directory
            env: Environment variables
            timeout: Seconds before the command is killed. Left unset, the active run's
                ``subprocess_timeout`` setting is used (``--timeout``, ``PREPARE_TIMEOUT``, or
                config ``[runtime] timeout``; see ``prepare.core.settings``). 0 or ``None``,
                given explicitly or resolved from the setting, means no limit.
            **kwargs: Additional subprocess.Popen arguments; ``check`` (default True)
                controls whether a non-zero exit raises, and ``capture_output`` is
                accepted but ignored because output is always captured

        Returns:
            CompletedProcess result

        Raises:
            SecurityError: If any validation fails, or the command times out
            subprocess.CalledProcessError: If the command exits non-zero and ``check`` is True

        The child is started with ``Popen`` in a session of its own (POSIX) and recorded, under
        ``stop``, until it finishes, so ``terminate_children`` can stop it and any process it
        started from another thread; a timeout kills the same process group. When ``stop`` (the calling run's
        token) is already set, or ``terminate_children`` ran without a token, no child is started:
        the call fails as a killed child would, with return code -9 (``CalledProcessError``
        when ``check`` is True).
        """
        cmd = cls._build_validated_command(executable, args)
        resolved_timeout: Optional[float]
        if isinstance(timeout, _TimeoutNotGiven):
            resolved_timeout = settings_module.active().subprocess_timeout
        else:
            resolved_timeout = timeout
        # 0 or None (given or resolved from the setting) means communicate() waits without a limit.
        communicate_timeout = resolved_timeout if resolved_timeout else None

        # Validate working directory if provided
        if cwd:
            cwd = cls.validate_path(cwd)

        # Set safe environment
        safe_env = os.environ.copy() if env is None else env.copy()
        # Remove potentially dangerous environment variables
        for var in DANGEROUS_ENV_VARS:
            safe_env.pop(var, None)

        logger.debug(f"Running secure command: {' '.join(cmd)}")

        check = kwargs.pop("check", True)
        kwargs.pop("capture_output", None)
        if cls._stop_requested(stop):
            logger.debug(f"Not starting {executable}: the run is stopping")
            if check:
                raise subprocess.CalledProcessError(_NOT_STARTED_RETURNCODE, cmd, output="", stderr="")
            return subprocess.CompletedProcess(cmd, _NOT_STARTED_RETURNCODE, "", "")
        try:
            # Use secure defaults
            popen_kwargs = {
                "stdout": subprocess.PIPE,
                "stderr": subprocess.PIPE,
                "text": True,
                "cwd": cwd,
                "env": safe_env,
                **kwargs,
                "start_new_session": _PROCESS_GROUPS,
            }
            proc = subprocess.Popen(cmd, **popen_kwargs)
            with cls._children_lock:
                if cls._stop_requested(stop):
                    # terminate_children ran for this run (or for every run) after the check
                    # above: stop this child too; communicate() below then reaps it and the
                    # non-zero exit is reported as usual.
                    _signal_tool(proc, kill=True)
                else:
                    cls._children[proc] = stop
            try:
                out, err = proc.communicate(timeout=communicate_timeout)
            except subprocess.TimeoutExpired:
                _signal_tool(proc, kill=True)
                proc.communicate()
                raise
            except BaseException:
                # As subprocess.run does: an interrupt in the waiting thread stops the child.
                _signal_tool(proc, kill=True)
                proc.wait()
                raise
            finally:
                with cls._children_lock:
                    cls._children.pop(proc, None)

            if check and proc.returncode != 0:
                raise subprocess.CalledProcessError(proc.returncode, cmd, output=out, stderr=err)
            return subprocess.CompletedProcess(cmd, proc.returncode, out, err)

        except subprocess.TimeoutExpired as e:
            raise SecurityError(
                f"Command timed out after {resolved_timeout:g} s " "(set --timeout or PREPARE_TIMEOUT; 0 disables)"
            ) from e
        except subprocess.CalledProcessError as e:
            # Re-raise as the original exception type for compatibility
            raise e
        except Exception as e:
            raise SecurityError(f"Subprocess execution failed: {e}")
