"""A timeout or a termination stops a tool together with the processes it started.

The conda ``megahit`` is a wrapper that runs ``megahit_core`` as a child holding the inherited
output pipes. The fake tool here does the same with ``sleep 60``: its child records its own
process id and then becomes ``sleep``, so the test can check that this grandchild of PREPARE
is gone, not only the tool itself.
"""

import os
import subprocess
import threading
import time
from pathlib import Path

import pytest

from prepare.core.exceptions import SecurityError
from prepare.utils.security import SecureSubprocess

pytestmark = pytest.mark.skipif(os.name != "posix", reason="process groups are POSIX only")


def _fake_wrapper(tmp_path: Path, ignore_term: bool = False) -> Path:
    """Write a fake ``megahit`` whose child (not an ``exec``) runs ``sleep 60``; return the pid file."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    pid_file = tmp_path / "child.pid"
    trap = "trap '' TERM; " if ignore_term else ""
    script = bin_dir / "megahit"
    script.write_text(f'#!/bin/sh\n/bin/sh -c "{trap}echo \\$\\$ > {pid_file}; exec sleep 60"\necho done\n')
    script.chmod(0o755)
    return pid_file


def _child_pid(pid_file: Path) -> int:
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline:
        text = pid_file.read_text().strip() if pid_file.exists() else ""
        if text:
            return int(text)
        time.sleep(0.05)
    raise AssertionError("the fake tool's child never started")


def _gone(pid: int, within: float = 5.0) -> bool:
    deadline = time.monotonic() + within
    while time.monotonic() < deadline:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return True
        time.sleep(0.05)
    return False


@pytest.fixture
def wrapper_env(tmp_path, monkeypatch):
    """PATH with the fake megahit first; the tracking table and stopping flag are left clean."""
    monkeypatch.setenv("PATH", f"{tmp_path / 'bin'}:/usr/bin:/bin")
    monkeypatch.setattr(SecureSubprocess, "_children", {})
    yield
    SecureSubprocess.terminate_children(grace=0)
    SecureSubprocess.clear_stopping()


def test_timeout_stops_a_tool_whose_child_holds_the_pipes(tmp_path, wrapper_env):
    pid_file = _fake_wrapper(tmp_path)
    started = time.monotonic()
    with pytest.raises(SecurityError, match="timed out after 2 s"):
        SecureSubprocess.run_secure("megahit", ["--version"], timeout=2)
    assert time.monotonic() - started < 15
    assert _gone(_child_pid(pid_file)), "the tool's child is still running after the timeout"


def test_terminate_children_stops_the_tool_and_its_child(tmp_path, wrapper_env):
    pid_file = _fake_wrapper(tmp_path)
    outcome = {}

    def run():
        try:
            SecureSubprocess.run_secure("megahit", ["--version"], timeout=0)
        except subprocess.CalledProcessError as e:
            outcome["error"] = e

    thread = threading.Thread(target=run)
    thread.start()
    child = _child_pid(pid_file)
    assert SecureSubprocess.terminate_children(grace=2.0) == 1
    thread.join(15)
    assert not thread.is_alive(), "run_secure still waits on the pipes after terminate_children"
    assert isinstance(outcome.get("error"), subprocess.CalledProcessError)
    assert _gone(child)


def test_a_child_that_ignores_sigterm_is_killed_after_the_grace(tmp_path, wrapper_env):
    pid_file = _fake_wrapper(tmp_path, ignore_term=True)
    thread = threading.Thread(target=lambda: SecureSubprocess.run_secure("megahit", ["--version"], check=False))
    thread.start()
    child = _child_pid(pid_file)
    started = time.monotonic()
    SecureSubprocess.terminate_children(grace=1.0)
    thread.join(15)
    assert not thread.is_alive()
    assert _gone(child)
    assert time.monotonic() - started < 10
