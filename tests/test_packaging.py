"""Packaging consistency: the typed marker ships and the version is declared consistently."""

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "check_version.py"


def _run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)


def test_py_typed_marker_exists():
    assert (REPO_ROOT / "prepare" / "py.typed").is_file()


def test_versions_agree():
    result = _run()
    assert result.returncode == 0, result.stdout


def test_matching_tag_passes_and_other_tag_fails():
    import prepare

    assert _run("--tag", f"v{prepare.__version__}").returncode == 0
    bad = _run("--tag", "v0.0.0")
    assert bad.returncode == 1
    assert len(bad.stdout.strip().splitlines()) == 1
