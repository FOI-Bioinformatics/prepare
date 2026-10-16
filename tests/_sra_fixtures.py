"""Helpers and fixtures shared by the test_data_sra_*.py files."""

import gzip
import shutil
import time
import pytest
from pathlib import Path
from unittest.mock import Mock
from prepare.data.sra import accession as accession_mod


def _staging_download(*results):
    """A ``download_accession`` stand-in that writes one FASTQ file into its output folder.

    Without a store, ``download_sra`` publishes the folder ``download_accession`` built, so a
    stand-in that reports success must also leave that folder behind. ``results`` are the
    ``(success, message)`` pairs returned, one per call in call order.
    """
    remaining = list(results)

    def _download(accession, output_folder, *args, **kwargs):
        acc_dir = Path(output_folder) / accession
        acc_dir.mkdir(parents=True, exist_ok=True)
        (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        return remaining.pop(0) if len(remaining) > 1 else remaining[0]

    return _download


def _slow_download_tools(state, dump_seconds=0.2, pigz_seconds=0.05):
    """A ``run_secure`` stand-in for fasterq-dump (run directly, no prefetch) and pigz.

    fasterq-dump writes a two-file pair into its ``-O`` folder after ``dump_seconds``; pigz
    compresses one file in place after ``pigz_seconds``, which is the window in which a
    folder compressed where it is published would show a half-compressed pair.
    """
    import threading as threading_module

    lock = threading_module.Lock()

    def run(executable, args, **kwargs):
        with lock:
            state.setdefault("calls", []).append((executable, list(args)))
        if executable == "fasterq-dump":
            time.sleep(dump_seconds)
            out_dir = Path(args[args.index("-O") + 1])
            out_dir.mkdir(parents=True, exist_ok=True)
            for mate in ("1", "2"):
                (out_dir / f"SRR1_{mate}.fastq").write_text("@r\nACGT\n+\nIIII\n")
        elif executable == "pigz":
            time.sleep(pigz_seconds)
            plain = Path(args[-1])
            with open(plain, "rb") as src, gzip.open(str(plain) + ".gz", "wb") as dest:
                shutil.copyfileobj(src, dest)
            plain.unlink()
        return Mock(returncode=0, stdout="", stderr="")

    return run


def _no_prefetch_with_pigz(tool):
    """``shutil.which`` stand-in: no prefetch (fasterq-dump runs directly), pigz present."""
    return "/usr/bin/pigz" if tool == "pigz" else None


@pytest.fixture
def fast_project_lock_poll(monkeypatch):
    """Poll a held project lock every 0.05 s so contention tests stay short."""
    monkeypatch.setattr(accession_mod, "PROJECT_LOCK_POLL_SECONDS", 0.05)


def _with_prefetch(tool):
    """``shutil.which`` stand-in: prefetch and fasterq-dump are on PATH."""
    return f"/usr/bin/{tool}"


def _scratch_recording_tools(state):
    """A ``run_secure`` stand-in for fasterq-dump that records its ``--temp`` folder.

    It notes whether that folder existed and whether an earlier run's file was still in it,
    then writes a FASTQ pair into its ``-O`` folder, or fails when ``state["fail"]`` is set.
    """
    import subprocess as subprocess_module

    def run(executable, args, **kwargs):
        state.setdefault("calls", []).append((executable, list(args)))
        if executable == "fasterq-dump":
            temp = Path(args[args.index("--temp") + 1])
            state["temp"] = str(temp)
            state["temp_existed"] = temp.is_dir()
            state["left_over_present"] = (temp / "left-over").exists()
            (temp / "scratch-file").parent.mkdir(parents=True, exist_ok=True)
            (temp / "scratch-file").write_text("scratch")
            if state.get("fail"):
                raise subprocess_module.CalledProcessError(1, ["fasterq-dump"], stderr="connection reset")
            out_dir = Path(args[args.index("-O") + 1])
            out_dir.mkdir(parents=True, exist_ok=True)
            for mate in ("1", "2"):
                (out_dir / f"SRR1_{mate}.fastq").write_text("@r\nACGT\n+\nIIII\n")
        return Mock(returncode=0, stdout="", stderr="")

    return run
