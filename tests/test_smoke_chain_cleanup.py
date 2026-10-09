"""scripts/smoke_chain.sh removes the scratch project directory it created, but only on success
outside CI; a failed run, a CI run, or a directory the caller named is left in place.

A stub ``python`` on PATH stands in for every prepare step, so these tests need no network or
bioinformatics tools.
"""

import os
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "smoke_chain.sh"

STUB = """#!/usr/bin/env bash
# Stand-in for "python -m prepare.cli.main <command> ...".
if [ -n "${SMOKE_STUB_FAIL:-}" ]; then exit 1; fi
if [ "$3" = "results_table" ]; then printf 'accession\\tgenome_id\\nSRR2517620\\tGCF\\n' >results.tsv; fi
exit 0
"""


@pytest.fixture
def stub_env(tmp_path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    stub = bin_dir / "python"
    stub.write_text(STUB)
    stub.chmod(0o755)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    env = {k: v for k, v in os.environ.items() if k not in ("CI", "SMOKE_STUB_FAIL")}
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    env["TMPDIR"] = str(scratch)
    return env, scratch


def _run(env, *args):
    return subprocess.run(["bash", str(SCRIPT), *args], env=env, capture_output=True, text=True, timeout=60)


def _smoke_dirs(scratch):
    return [p for p in scratch.rglob("prepare-smoke-*") if p.is_dir()]


def test_success_outside_ci_removes_the_scratch_directory(stub_env):
    env, scratch = stub_env
    result = _run(env)
    assert result.returncode == 0, result.stdout + result.stderr
    assert _smoke_dirs(scratch) == []


def test_failure_keeps_the_scratch_directory(stub_env):
    env, scratch = stub_env
    env["SMOKE_STUB_FAIL"] = "1"
    result = _run(env)
    assert result.returncode != 0
    assert len(_smoke_dirs(scratch)) == 1


def test_success_in_ci_keeps_the_scratch_directory(stub_env):
    env, scratch = stub_env
    env["CI"] = "true"
    result = _run(env)
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(_smoke_dirs(scratch)) == 1


def test_a_directory_the_caller_names_is_kept(stub_env, tmp_path):
    env, _ = stub_env
    project = tmp_path / "project"
    result = _run(env, str(project))
    assert result.returncode == 0, result.stdout + result.stderr
    assert (project / "results.tsv").is_file()
