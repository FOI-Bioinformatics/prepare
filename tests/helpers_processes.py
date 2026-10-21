"""Helpers for tests that run the prepare CLI as separate processes.

The processes find fake ``fasterq-dump``, ``prefetch`` and ``pigz`` scripts first on ``PATH``
(``install_fake_tools``), so no real SRA tool or network is used. Each fake records that it
started and its pid in a barrier folder, and waits while a ``hold-<ACC>`` file exists there, so a
test can keep a download in progress for exactly as long as it needs and release it by removing
that file. ``cli_env`` gives every process its own ``HOME`` and configuration folder under the
test's ``tmp_path``, so no process reads or writes the developer's store or configuration.
"""

import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Union

REPO_ROOT = Path(__file__).resolve().parent.parent

# Read records per FASTQ file a fake fasterq-dump writes.
FAKE_READS = 4

# The body shared by the three fake tools; ``TOOL`` is set per script. A download tool (prefetch,
# fasterq-dump) writes ``<barrier>/<ACC>.<pid>.started``, appends its pid to ``<barrier>/pids``
# and waits while ``<barrier>/hold-<ACC>`` exists; pigz only records its pid. All three exit 143
# on SIGTERM, as a real tool killed by that signal is reported by a shell; while
# ``<barrier>/linger-<ACC>`` exists, a download tool first writes ``<ACC>.<pid>.terminating`` and
# delays that exit until the file is removed, like a tool that takes a while to shut down. While
# ``<barrier>/ignore-term-<ACC>`` exists, a download tool ignores SIGTERM altogether (it writes
# ``<ACC>.<pid>.ignored``), so only SIGKILL ends it.
#
# Failure modes, read from the environment of the fake: ``FAKE_FASTERQ_EXIT="<code>:<stderr>"``
# makes every fasterq-dump call exit with that code and text after writing nothing;
# ``FAKE_FASTERQ_FAIL_ONCE=1`` does so for the first call per accession only (a marker in the
# barrier folder remembers it); ``FAKE_FASTERQ_SPOTS=<n>`` sets the records per mate file;
# ``FAKE_FASTERQ_UNPAIRED=1`` also writes the bare ``<ACC>.fastq`` of unpaired spots that
# ``--split-3`` produces; ``FAKE_PREFETCH_NOOP=1`` makes prefetch exit 0 without an archive.
_FAKE_TOOL_BODY = r"""
import gzip
import os
import shutil
import signal
import sys
import time
from pathlib import Path

BARRIER = Path(BARRIER_DIR)
VALUE_FLAGS = {"-O", "--outdir", "-o", "--outfile", "--temp", "-t", "--threads", "-e", "--max-size", "-p"}


CURRENT = []


def on_term(signum, frame):
    # A test of the SIGKILL that follows the grace period creates ignore-term-<ACC>.
    if CURRENT and (BARRIER / f"ignore-term-{CURRENT[0]}").exists():
        (BARRIER / f"{CURRENT[0]}.{os.getpid()}.ignored").write_text("")
        return
    # A test that needs the run to still be stopping when it sends a second signal creates
    # linger-<ACC>: the tool then marks that it got the signal and exits once the file is gone.
    if CURRENT and (BARRIER / f"linger-{CURRENT[0]}").exists():
        (BARRIER / f"{CURRENT[0]}.{os.getpid()}.terminating").write_text("")
        while (BARRIER / f"linger-{CURRENT[0]}").exists():
            time.sleep(0.05)
    os._exit(143)


def positional(args):
    skip = False
    for token in args:
        if skip:
            skip = False
            continue
        if token in VALUE_FLAGS:
            skip = True
            continue
        if token.startswith("-"):
            continue
        return token
    return None


def value_of(args, *flags):
    for index, token in enumerate(args):
        if token in flags and index + 1 < len(args):
            return args[index + 1]
    return None


def record_pid():
    with open(BARRIER / "pids", "a") as handle:
        handle.write(f"{os.getpid()}\n")


def start_and_hold(accession):
    CURRENT.append(accession)
    record_pid()
    (BARRIER / f"{accession}.{os.getpid()}.started").write_text(TOOL + "\n")
    while (BARRIER / f"hold-{accession}").exists():
        time.sleep(0.05)


def fastq_record(accession, mate, index):
    return f"@{accession}.{index} read {index}/{mate}\nACGTACGTAC\n+\nIIIIIIIIII\n"


def fail_once(accession):
    # The first fasterq-dump call per accession under FAKE_FASTERQ_FAIL_ONCE.
    if not os.environ.get("FAKE_FASTERQ_FAIL_ONCE"):
        return False
    marker = BARRIER / f"failed-once-{accession}"
    if marker.exists():
        return False
    marker.write_text("")
    return True


def fasterq_dump(args):
    source = positional(args)
    accession = Path(source).name.split(".")[0]
    outdir = Path(value_of(args, "-O", "--outdir") or ".")
    start_and_hold(accession)
    outdir.mkdir(parents=True, exist_ok=True)
    exit_spec = os.environ.get("FAKE_FASTERQ_EXIT")
    if exit_spec or fail_once(accession):
        code, _, text = (exit_spec or "1:err: connection reset by peer").partition(":")
        sys.stderr.write(text + "\n")
        sys.exit(int(code))
    reads = int(os.environ.get("FAKE_FASTERQ_SPOTS", READS))
    for mate in (1, 2):
        text = "".join(fastq_record(accession, mate, index) for index in range(1, reads + 1))
        (outdir / f"{accession}_{mate}.fastq").write_text(text)
    if os.environ.get("FAKE_FASTERQ_UNPAIRED"):
        (outdir / f"{accession}.fastq").write_text(fastq_record(accession, 0, reads + 1))


def prefetch(args):
    accession = positional(args)
    outdir = Path(value_of(args, "-O", "--output-directory") or ".")
    start_and_hold(accession)
    if os.environ.get("FAKE_PREFETCH_NOOP"):
        return
    folder = outdir / accession
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{accession}.sra").write_bytes(b"fake sra archive\n")


def pigz(args):
    record_pid()
    path = Path(positional(args))
    target = path.with_name(path.name + ".gz")
    with open(path, "rb") as source, gzip.open(target, "wb") as sink:
        shutil.copyfileobj(source, sink)
    if "-k" not in args:
        path.unlink()


def main():
    signal.signal(signal.SIGTERM, on_term)
    args = sys.argv[1:]
    if "--version" in args:
        print(f"{TOOL} : 3.0.0-fake")
        return 0
    {"fasterq-dump": fasterq_dump, "prefetch": prefetch, "pigz": pigz}[TOOL](args)
    return 0


sys.exit(main())
"""

FAKE_TOOLS = ("fasterq-dump", "prefetch", "pigz")


def install_fake_tools(bin_dir: Union[str, Path], barrier_dir: Union[str, Path]) -> None:
    """Write executable fake ``fasterq-dump``, ``prefetch`` and ``pigz`` scripts into ``bin_dir``.

    The scripts run under this interpreter (``sys.executable`` in the shebang) and accept the
    argument shapes ``prepare.data.sra.accession`` passes: fasterq-dump with ``-O <folder>``
    and the accession or a ``.sra`` path, writing ``<ACC>_1.fastq`` and ``<ACC>_2.fastq`` with
    ``FAKE_READS`` reads each; prefetch with ``-O <cache>``, writing ``<cache>/<ACC>/<ACC>.sra``;
    pigz with ``-p N -f <file>``, replacing the file with ``<file>.gz``.
    """
    bin_path = Path(bin_dir)
    barrier = Path(barrier_dir)
    bin_path.mkdir(parents=True, exist_ok=True)
    barrier.mkdir(parents=True, exist_ok=True)
    for tool in FAKE_TOOLS:
        header = f"#!{sys.executable}\nTOOL = {tool!r}\nBARRIER_DIR = {str(barrier)!r}\nREADS = {FAKE_READS}\n"
        script = bin_path / tool
        script.write_text(header + _FAKE_TOOL_BODY)
        script.chmod(0o755)


def cli_env(
    tmp_path: Union[str, Path], fake_bin: Union[str, Path], store: Optional[Union[str, Path]] = None, **extra: str
) -> Dict[str, str]:
    """The environment for a CLI child process: this checkout, an isolated home, only fake tools.

    ``PYTHONPATH`` is the repository root, so the child imports this checkout's package. ``HOME``
    and ``XDG_CONFIG_HOME`` point under ``tmp_path``, so no configuration file of the developer's
    is read. ``PREPARE_DATA`` is ``store`` when given and removed otherwise. ``PATH`` holds
    ``fake_bin`` and a folder with nothing but ``python``/``python3`` links to this interpreter:
    not the interpreter's own folder, which in a conda environment also holds the real SRA tools,
    samtools and minimap2. ``extra`` entries are added last.
    """
    base = Path(tmp_path)
    home = base / "home"
    (home / ".config").mkdir(parents=True, exist_ok=True)
    python_bin = base / "python-bin"
    python_bin.mkdir(exist_ok=True)
    for name in ("python", "python3"):
        link = python_bin / name
        if not link.exists():
            link.symlink_to(sys.executable)
    env = dict(os.environ)
    for name in ("PREPARE_DATA", "PYTHONSTARTUP", "PYTHONHOME", "VIRTUAL_ENV"):
        env.pop(name, None)
    env.update(
        {
            "PYTHONPATH": str(REPO_ROOT),
            "HOME": str(home),
            "XDG_CONFIG_HOME": str(home / ".config"),
            "PATH": os.pathsep.join([str(fake_bin), str(python_bin)]),
            "PYTHONUNBUFFERED": "1",
        }
    )
    if store is not None:
        env["PREPARE_DATA"] = str(store)
    env.update(extra)
    return env


def _log_paths(cwd: Path, log_dir: Optional[Path]) -> tuple:
    """Unused stdout and stderr file paths for one child, in ``log_dir`` (default: beside ``cwd``)."""
    folder = Path(log_dir) if log_dir is not None else cwd.parent / "process-logs"
    folder.mkdir(parents=True, exist_ok=True)
    index = 0
    while (folder / f"{cwd.name}-{index}.stderr").exists():
        index += 1
    return folder / f"{cwd.name}-{index}.stdout", folder / f"{cwd.name}-{index}.stderr"


def spawn_cli(
    args: Sequence[str], cwd: Union[str, Path], env: Dict[str, str], log_dir: Optional[Union[str, Path]] = None
) -> subprocess.Popen:
    """Start ``python -m prepare.cli.main <args>`` in ``cwd`` in a session of its own.

    ``start_new_session`` keeps a signal sent to the child (or its process group) away from the
    test runner. Output goes to files rather than pipes, so a child that writes a lot never
    blocks and a test can read its stderr while it runs: the paths are on the returned object
    as ``stdout_path`` and ``stderr_path``.
    """
    work = Path(cwd)
    stdout_path, stderr_path = _log_paths(work, Path(log_dir) if log_dir is not None else None)
    stdout_path.touch()
    stderr_path.touch()
    with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
        proc = subprocess.Popen(
            [sys.executable, "-m", "prepare.cli.main", *args],
            cwd=str(work),
            env=env,
            stdin=subprocess.DEVNULL,
            stdout=out,
            stderr=err,
            start_new_session=True,
        )
    proc.stdout_path = stdout_path  # type: ignore[attr-defined]
    proc.stderr_path = stderr_path  # type: ignore[attr-defined]
    return proc


def stderr_of(proc: subprocess.Popen) -> str:
    """What a child started by ``spawn_cli`` has written to stderr so far."""
    try:
        return Path(proc.stderr_path).read_text(errors="replace")  # type: ignore[attr-defined]
    except OSError:
        return ""


def stdout_of(proc: subprocess.Popen) -> str:
    """What a child started by ``spawn_cli`` has written to stdout so far."""
    try:
        return Path(proc.stdout_path).read_text(errors="replace")  # type: ignore[attr-defined]
    except OSError:
        return ""


def run_cli(
    args: Sequence[str], cwd: Union[str, Path], env: Dict[str, str], timeout: float = 20.0
) -> subprocess.CompletedProcess:
    """Run one CLI command to completion and return its exit code, stdout and stderr."""
    proc = spawn_cli(args, cwd, env)
    try:
        returncode = proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        stop_process(proc)
        raise AssertionError(f"prepare {' '.join(args)} did not finish in {timeout} s:\n{stderr_of(proc)}")
    return subprocess.CompletedProcess(list(args), returncode, stdout_of(proc), stderr_of(proc))


def wait_for(
    predicate: Callable[[], bool],
    timeout: float = 10.0,
    interval: float = 0.05,
    what: str = "the condition",
    detail: Optional[Callable[[], str]] = None,
) -> None:
    """Poll ``predicate`` until it is true; raise ``AssertionError`` naming ``what`` after ``timeout``.

    ``detail``, when given, is called on a timeout and its text (for example a child's stderr)
    is added to the message.
    """
    deadline = time.monotonic() + timeout
    while True:
        if predicate():
            return
        if time.monotonic() >= deadline:
            extra = f"\n{detail()}" if detail is not None else ""
            raise AssertionError(f"timed out after {timeout} s waiting for {what}{extra}")
        time.sleep(interval)


def alive(pid: int) -> bool:
    """Whether a process with ``pid`` exists (a zombie not yet reaped counts as existing)."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def stop_process(proc: subprocess.Popen, grace: float = 3.0) -> None:
    """End a child and its process group: SIGTERM, then SIGKILL after ``grace`` seconds; reap it."""
    if proc.poll() is None:
        _signal_group(proc.pid, signal.SIGTERM)
        try:
            proc.wait(timeout=grace)
        except subprocess.TimeoutExpired:
            pass
    # Whatever is left in the group is killed outright (a tool run_secure started is in a group of
    # its own; the harness kills any such fake tool by pid).
    _signal_group(proc.pid, signal.SIGKILL)
    try:
        proc.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        pass


def _signal_group(pgid: int, signum: int) -> None:
    """Send ``signum`` to process group ``pgid``, ignoring a group that no longer exists."""
    try:
        os.killpg(pgid, signum)
    except (ProcessLookupError, PermissionError):
        pass


def started_files(barrier_dir: Union[str, Path], accession: str) -> List[Path]:
    """The ``<ACC>.<pid>.started`` files the fake download tools wrote for ``accession``."""
    return sorted(Path(barrier_dir).glob(f"{accession}.*.started"))


def fake_pids(barrier_dir: Union[str, Path]) -> List[int]:
    """Every pid a fake tool recorded in ``<barrier>/pids``."""
    path = Path(barrier_dir) / "pids"
    if not path.exists():
        return []
    return [int(line) for line in path.read_text().split() if line.strip()]


def pid_of_started(path: Path) -> int:
    """The fake tool's pid encoded in a ``<ACC>.<pid>.started`` file name."""
    return int(path.name.split(".")[1])


def is_fake_tool(pid: int, bin_dir: Union[str, Path]) -> bool:
    """Whether ``pid`` is still one of the fake tools in ``bin_dir``, judged by its command line.

    A recorded pid can be reused by an unrelated process once the fake has exited; checking the
    command line (which names the script under ``bin_dir``) keeps cleanup from killing it.
    """
    ps = "/bin/ps" if Path("/bin/ps").exists() else "ps"
    try:
        result = subprocess.run(
            [ps, "-o", "command=", "-p", str(pid)], capture_output=True, text=True, timeout=5, check=False
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return str(Path(bin_dir)) in result.stdout
