"""Tests for scripts/check_atomic_writes.sh, the gate against direct file writes under prepare/."""

import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
GATE = REPO_ROOT / "scripts" / "check_atomic_writes.sh"


def _run_gate(root: Path) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", str(GATE), str(root)], capture_output=True, text=True)


def _fixture_tree(tmp_path: Path, body: str) -> Path:
    pkg = tmp_path / "prepare" / "sub"
    pkg.mkdir(parents=True)
    (pkg / "mod.py").write_text(body)
    scripts = tmp_path / "scripts"
    scripts.mkdir()
    (scripts / "atomic_writes_allowlist.txt").write_text("prepare/allowed.py  test fixture\n")
    (tmp_path / "prepare" / "allowed.py").write_text('open(p, "w")\n')
    return tmp_path


@pytest.mark.parametrize(
    "line",
    [
        "x.write_text(data)\n",
        "x.write_bytes(data)\n",
        'df.to_csv("out.tsv", sep="\\t")\n',
        "df.to_csv(path)\n",
        "df.to_csv(path_or_buf=path)\n",
        'with open(path, "w") as f:\n',
        "with open(path, 'wb') as f:\n",
        'with gzip.open(path, mode="wt") as f:\n',
        "shutil.copy2(src, dest)\n",
        "shutil.copyfile(src, dest)\n",
        "fig.write_html(str(path))\n",
        "plt.savefig(path, dpi=300)\n",
        "with path.open('w') as f:\n",
        'with path.open(mode="wb") as f:\n',
        "with open(p, 'x') as f:\n",
        "with open(p, mode) as f:\n",
        "with open(p, mode=write_mode, encoding='utf-8') as f:\n",
    ],
)
def test_gate_fails_on_direct_write(tmp_path, line):
    result = _run_gate(_fixture_tree(tmp_path, line))
    assert result.returncode == 1, result.stdout + result.stderr
    assert "prepare/sub/mod.py" in result.stdout


@pytest.mark.parametrize(
    "line",
    [
        'text = df.to_csv(sep="\\t")\n',
        "text = df.to_csv(index=False)\n",
        "text = df.to_csv()\n",
        'with open(path, "rb") as f:\n',
        'with open(path, "a") as f:\n',
        "write_text_atomic(path, text)\n",
        'with open_atomic(path, "w") as f:\n',
        "# a comment about x.write_text(data)\n",
        "with path.open('rb') as f:\n",
        "with open(p, encoding='utf-8') as f:\n",
        "fd = os.open(lock, os.O_CREAT | os.O_EXCL)\n",
        "webbrowser.open(uri)\n",
    ],
)
def test_gate_passes_on_allowed_forms(tmp_path, line):
    result = _run_gate(_fixture_tree(tmp_path, line))
    assert result.returncode == 0, result.stdout + result.stderr


def test_gate_passes_on_the_real_tree():
    result = _run_gate(REPO_ROOT)
    assert result.returncode == 0, result.stdout + result.stderr


def test_gate_rejects_a_stale_allowlist_entry(tmp_path):
    root = _fixture_tree(tmp_path, "pass\n")
    (root / "scripts" / "atomic_writes_allowlist.txt").write_text(
        "prepare/allowed.py  fixture\nprepare/gone.py  removed module\n"
    )
    result = _run_gate(root)
    assert result.returncode == 1
    assert "prepare/gone.py" in result.stdout


def test_gate_rejects_an_entry_without_reason(tmp_path):
    root = _fixture_tree(tmp_path, "pass\n")
    (root / "scripts" / "atomic_writes_allowlist.txt").write_text("prepare/allowed.py\n")
    result = _run_gate(root)
    assert result.returncode == 1
    assert "no reason" in result.stdout
