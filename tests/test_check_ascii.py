"""scripts/check_ascii.sh: non-ASCII bytes fail the gate unless the line carries an ``# ascii-ok``
marker; the exemption is per line, not per file.
"""

import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_ascii.sh"
# Built with chr() so this file itself stays ASCII-only.
E_ACUTE = chr(0xE9)


def _repo(tmp_path, files):
    for name, text in files.items():
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "add", "-A"], cwd=tmp_path, check=True)
    return tmp_path


def _export(tmp_path, files):
    """A source tree that is not a git checkout, as `git archive` or an unpacked sdist leaves it."""
    for name, text in files.items():
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    return tmp_path


def _gate(root):
    return subprocess.run(["bash", str(SCRIPT), str(root)], capture_output=True, text=True, timeout=60)


def test_an_export_without_git_is_still_checked(tmp_path):
    root = _export(tmp_path, {"prepare/x.py": f'LABEL = "caf{E_ACUTE}"\n', "Makefile": "all:\n\ttrue\n"})
    result = _gate(root)
    assert result.returncode == 1
    assert "prepare/x.py" in result.stdout


def test_a_clean_export_without_git_passes(tmp_path):
    root = _export(tmp_path, {"tests/test_x.py": "X = 1\n", "scripts/run.sh": "echo ok\n", "setup.cfg": "[flake8]\n"})
    result = _gate(root)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "No non-ASCII" in result.stdout


def test_a_root_with_no_source_files_fails_rather_than_passing_silently(tmp_path):
    result = _gate(_export(tmp_path, {"README.md": "text\n"}))
    assert result.returncode == 1
    assert "nothing was checked" in result.stdout


def test_unmarked_non_ascii_line_in_a_test_file_fails(tmp_path):
    root = _repo(tmp_path, {"tests/test_x.py": f'LABEL = "caf{E_ACUTE}"\n'})
    result = _gate(root)
    assert result.returncode == 1
    assert "tests/test_x.py" in result.stdout


def test_a_line_marked_ascii_ok_passes(tmp_path):
    root = _repo(tmp_path, {"tests/test_x.py": f'LABEL = "caf{E_ACUTE}"  # ascii-ok: fixture\n'})
    result = _gate(root)
    assert result.returncode == 0, result.stdout + result.stderr


def test_the_marker_covers_only_its_own_line(tmp_path):
    text = f'A = "{E_ACUTE}"  # ascii-ok: fixture\nB = "{E_ACUTE}"\n'
    root = _repo(tmp_path, {"tests/test_security_comprehensive.py": text})
    result = _gate(root)
    assert result.returncode == 1
    assert "2:B = " in result.stdout
    assert "1:A = " not in result.stdout


def test_pyproject_is_checked_line_by_line_too(tmp_path):
    text = f'authors = [\n    {{name = "{E_ACUTE}"}}  # ascii-ok: author name\n]\ndescription = "{E_ACUTE}"\n'
    root = _repo(tmp_path, {"pyproject.toml": text})
    result = _gate(root)
    assert result.returncode == 1
    assert "4:description" in result.stdout
    assert "2:" not in result.stdout
