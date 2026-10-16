"""Tests for the output helpers on `BaseCommand` and the one-stdout-channel rule.

stdout carries a command's result (tables, JSON); stderr carries logging. The JSON-mode tests
run each ``--json`` command through ``prepare.cli.main.main``, so ``setup_logging`` binds its
real stderr handler, and then parse the whole of stdout as one JSON document: a stray print or
a log line on stdout would fail them.

Every store test runs under tmp_path and monkeypatches HOME/XDG_CONFIG_HOME/PREPARE_DATA so
nothing here reads or writes the real user config.
"""

import argparse
import json
import logging
import subprocess
from pathlib import Path
from typing import List

import pytest

from prepare.cli.base import BaseCommand, emit_error_json
from prepare.cli.commands.store import StoreInitCommand
from prepare.cli.main import main
from prepare.core.constants import STORE_ENV
from prepare.store.catalog import catalog_write
from prepare.store.layout import init_store, store_paths
from prepare.store.sidecar import Sidecar

REPO_ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(autouse=True)
def isolated_env(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    monkeypatch.delenv(STORE_ENV, raising=False)
    yield


class _EmitCommand(BaseCommand):
    """A minimal command that writes whatever it is given through the helpers."""

    def __init__(self, action):
        super().__init__()
        self._action = action

    @property
    def name(self) -> str:
        return "x"

    @property
    def help(self) -> str:
        return "x"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        pass

    def execute(self, args: argparse.Namespace) -> int:
        self._action(self)
        return 0


class TestEmitHelpers:
    def test_emit_json_writes_exactly_one_document(self, capsys):
        _EmitCommand(lambda cmd: cmd.emit_json({"a": 1, "b": [1, 2]})).execute(argparse.Namespace())

        captured = capsys.readouterr()
        assert json.loads(captured.out) == {"a": 1, "b": [1, 2]}
        assert captured.out.endswith("}\n")
        assert captured.err == ""

    def test_emit_json_is_indented_by_two(self, capsys):
        _EmitCommand(lambda cmd: cmd.emit_json({"a": 1})).execute(argparse.Namespace())

        assert capsys.readouterr().out == '{\n  "a": 1\n}\n'

    def test_emit_writes_one_newline_terminated_line(self, capsys):
        _EmitCommand(lambda cmd: cmd.emit("hello")).execute(argparse.Namespace())

        assert capsys.readouterr().out == "hello\n"

    def test_emit_raw_writes_the_text_as_is(self, capsys):
        _EmitCommand(lambda cmd: cmd.emit_raw("a\tb\n1\t2")).execute(argparse.Namespace())

        assert capsys.readouterr().out == "a\tb\n1\t2"

    def test_emit_without_text_writes_an_empty_line(self, capsys):
        _EmitCommand(lambda cmd: cmd.emit()).execute(argparse.Namespace())

        assert capsys.readouterr().out == "\n"

    def test_emit_error_json_writes_one_error_document(self, capsys):
        emit_error_json("something went wrong")

        captured = capsys.readouterr()
        assert json.loads(captured.out) == {"error": "something went wrong"}
        assert captured.err == ""


# ---------------------------------------------------------------- JSON commands


@pytest.fixture
def restore_root_logging():
    """``main`` replaces the root logger's handlers with a stderr handler bound to the stream
    capsys installed; put the previous handlers and level back so later tests do not log to a
    closed stream."""
    root = logging.getLogger()
    handlers, level = root.handlers[:], root.level
    yield
    for handler in root.handlers[:]:
        root.removeHandler(handler)
    for handler in handlers:
        root.addHandler(handler)
    root.setLevel(level)


def _corrupt_project_registry(tmp_path: Path) -> Path:
    """A store with one registered project whose registry no longer parses.

    ``stale_projects`` (used by ``store_status`` and ``store_gc``) logs a warning for such a
    project instead of failing.
    """
    root = tmp_path / "store"
    project_dir = tmp_path / "project"
    project_dir.mkdir()
    registry = project_dir / "prepare_registry.json"
    args = argparse.Namespace(data_root=str(root), project_name=None, set_default=False, registry=str(registry))
    assert StoreInitCommand().execute(args) == 0
    registry.write_text("{ not json")
    return root


def _store_status_argv(tmp_path: Path) -> List[str]:
    root = _corrupt_project_registry(tmp_path)
    return ["store_status", "--data-root", str(root), "--details", "--json"]


def _store_gc_argv(tmp_path: Path) -> List[str]:
    root = _corrupt_project_registry(tmp_path)
    return ["store_gc", "--data-root", str(root), "--dry-run", "--json"]


def _store_usage_argv(tmp_path: Path) -> List[str]:
    root = tmp_path / "store"
    paths = init_store(root)
    with catalog_write(paths) as cat:
        cat.upsert_project("proja", "Wolbachia", str(tmp_path), str(tmp_path / "prepare_registry.json"))
        cat.upsert_dataset(Sidecar(accession="SRR1", state="complete"))
        cat.record_usage("SRR1", "proja", "wMel", "downloaded")
    return ["store_usage", "--data-root", str(store_paths(root).root), "--accession", "SRR404", "--json"]


def _status_argv(tmp_path: Path) -> List[str]:
    # The store this project names is not there: status warns and still reports.
    for folder in ("fastq", "metadata", "genomes"):
        (tmp_path / folder).mkdir()
    return [
        "status",
        "--fastq-folder",
        str(tmp_path / "fastq"),
        "--metadata-folder",
        str(tmp_path / "metadata"),
        "--genomes-folder",
        str(tmp_path / "genomes"),
        "--targeted-folder",
        str(tmp_path / "targeted"),
        "--matches-folder",
        str(tmp_path / "matches"),
        "--registry",
        str(tmp_path / "prepare_registry.json"),
        "--data-root",
        str(tmp_path / "unmounted"),
        "--json",
    ]


# (argv builder, a WARNING the state provokes or None). store_usage has no warning path of its
# own, so for it the test only shows that its INFO line lands on stderr.
JSON_COMMANDS = [
    pytest.param(_store_status_argv, "Could not check registry", id="store_status"),
    pytest.param(_store_gc_argv, "Could not check registry", id="store_gc"),
    pytest.param(_store_usage_argv, None, id="store_usage"),
    pytest.param(_status_argv, "store unavailable", id="status"),
]


@pytest.mark.parametrize("argv_builder,warning", JSON_COMMANDS)
def test_json_mode_stdout_is_a_single_document_and_logs_go_to_stderr(
    argv_builder, warning, tmp_path, monkeypatch, capsys, restore_root_logging
):
    monkeypatch.chdir(tmp_path)
    argv = argv_builder(tmp_path)
    capsys.readouterr()

    rc = main(argv)

    captured = capsys.readouterr()
    assert rc == 0, captured.err
    json.loads(captured.out)  # the whole stream parses as one document
    assert "Resolved store root" in captured.err
    assert "Resolved store root" not in captured.out
    if warning:
        assert "WARNING" in captured.err
        assert warning in captured.err


def test_no_library_module_prints():
    """The gate's own pattern, limited to the library packages."""
    hits = subprocess.run(
        [
            "grep",
            "-rlnE",
            "--include=*.py",
            r"(^|[^.a-zA-Z0-9_])print\(|sys\.stdout\.write",
            "prepare/data",
            "prepare/store",
            "prepare/processing",
            "prepare/plugins",
            "prepare/sra",
            "prepare/visualization",
        ],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert hits.stdout.strip() == ""


def test_print_gate_script_passes_on_the_tree():
    result = subprocess.run(["bash", "scripts/check_no_print.sh"], capture_output=True, text=True, cwd=REPO_ROOT)
    assert result.returncode == 0, result.stdout + result.stderr


def test_print_gate_script_fails_on_a_print_outside_base(tmp_path):
    (tmp_path / "prepare" / "cli").mkdir(parents=True)
    (tmp_path / "prepare" / "cli" / "base.py").write_text('print("allowed")\n')
    (tmp_path / "prepare" / "library.py").write_text('def f():\n    print("not allowed")\n')
    (tmp_path / "prepare" / "writer.py").write_text('import sys\nsys.stdout.write("not allowed")\n')
    (tmp_path / "prepare" / "other.py").write_text("import pprint\npprint.pprint(1)\nself._print_table(1)\n")

    result = subprocess.run(
        ["bash", str(REPO_ROOT / "scripts" / "check_no_print.sh"), str(tmp_path)], capture_output=True, text=True
    )

    assert result.returncode == 1
    assert "prepare/library.py:2:" in result.stdout
    assert "prepare/writer.py:2:" in result.stdout
    assert "prepare/cli/base.py:1:" not in result.stdout
    assert "prepare/other.py:" not in result.stdout


class TestAnchoredFolder:
    """Relative project folders live beside the registry the command uses (audit H4 follow-up)."""

    def test_no_registry_above_resolves_against_the_working_directory(self, tmp_path, monkeypatch):
        from prepare.cli.base import anchored_folder

        monkeypatch.chdir(tmp_path)
        assert anchored_folder("fastq", None) == tmp_path.resolve() / "fastq"

    def test_registry_above_the_working_directory_anchors_on_its_folder(self, tmp_path, monkeypatch):
        from prepare.cli.base import anchored_folder

        proj = tmp_path / "proj"
        (proj / "sub").mkdir(parents=True)
        (proj / "prepare_registry.json").write_text("{}")
        monkeypatch.chdir(proj / "sub")
        assert anchored_folder("fastq", None) == proj.resolve() / "fastq"

    def test_explicit_registry_anchors_on_its_folder(self, tmp_path, monkeypatch):
        from prepare.cli.base import anchored_folder

        monkeypatch.chdir(tmp_path)
        explicit = tmp_path / "other" / "prepare_registry.json"
        assert anchored_folder(Path("metadata"), str(explicit)) == (tmp_path / "other").resolve() / "metadata"

    def test_absolute_folder_is_returned_unchanged(self, tmp_path, monkeypatch):
        from prepare.cli.base import anchored_folder

        proj = tmp_path / "proj"
        (proj / "sub").mkdir(parents=True)
        (proj / "prepare_registry.json").write_text("{}")
        monkeypatch.chdir(proj / "sub")
        elsewhere = tmp_path / "elsewhere" / "fastq"
        assert anchored_folder(str(elsewhere), None) == elsewhere
