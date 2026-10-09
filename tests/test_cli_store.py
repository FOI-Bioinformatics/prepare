"""Tests for the `store_init`, `store_status` and `store_reindex` CLI commands.

Every test runs under tmp_path and monkeypatches HOME/XDG_CONFIG_HOME/PREPARE_DATA so
nothing here reads or writes the real user config, matching the isolation pattern used in
tests/test_store_resolve.py.
"""

import argparse
import gzip
import json
import logging
import subprocess
import threading
from unittest.mock import patch

import pytest

from prepare.cli.commands.store import (
    StoreAdoptCommand,
    StoreGcCommand,
    StoreInitCommand,
    StoreLinkCommand,
    StoreReindexCommand,
    StoreStatusCommand,
    StoreUnlinkCommand,
    StoreUsageCommand,
    StoreVerifyCommand,
)
from prepare.core.constants import STORE_ENV
from prepare.data.registry import load_registry, save_registry
from prepare.store.catalog import Catalog, catalog_write
from prepare.store.layout import init_store, lock_path, read_marker, sidecar_path, sra_dir, store_paths
from prepare.store.locks import dataset_lock
from prepare.store.sidecar import Sidecar, read_sidecar, write_sidecar
from prepare.utils.security import SecureSubprocess


@pytest.fixture(autouse=True)
def isolated_env(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    monkeypatch.delenv(STORE_ENV, raising=False)
    yield


@pytest.fixture(autouse=True)
def reset_allowed_roots():
    SecureSubprocess._extra_roots = []
    yield
    SecureSubprocess._extra_roots = []


def _init_args(root, project_dir, **overrides):
    base = dict(
        data_root=str(root),
        project_name=None,
        set_default=False,
        registry=str(project_dir / "prepare_registry.json"),
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _status_args(**overrides):
    base = dict(data_root=None, registry=None, json=False, verbose=False)
    base.update(overrides)
    return argparse.Namespace(**base)


def _reindex_args(**overrides):
    base = dict(data_root=None, registry=None)
    base.update(overrides)
    return argparse.Namespace(**base)


def _adopt_args(**overrides):
    base = dict(
        fastq_folder="fastq",
        data_root=None,
        registry=None,
        move=True,
        dry_run=False,
        compress=True,
        metadata_folder="metadata",
        lock_wait=0.0,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _verify_args(accessions=None, **overrides):
    base = dict(
        accessions=list(accessions or []),
        data_root=None,
        registry=None,
        md5=False,
        spots=False,
        fix_state=False,
        rescan=False,
        json=False,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _link_args(accessions, **overrides):
    base = dict(
        accessions=list(accessions),
        fastq_folder="fastq",
        registry=None,
        data_root=None,
        link_mode="auto",
        accept_partial=False,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _unlink_args(accessions, **overrides):
    base = dict(accessions=list(accessions), fastq_folder="fastq", registry=None)
    base.update(overrides)
    return argparse.Namespace(**base)


def _gc_args(**overrides):
    base = dict(
        data_root=None,
        registry=None,
        dry_run=False,
        yes=False,
        older_than=None,
        keep_partial=False,
        include_stale=False,
        accept_rebuilt=False,
        json=False,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _usage_args(**overrides):
    base = dict(
        data_root=None,
        registry=None,
        accession=None,
        project=None,
        organism=None,
        unused=False,
        bytes_by_organism=False,
        json=False,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _sidecar(accession="SRR1", state="complete"):
    return Sidecar(
        accession=accession,
        state=state,
        layout="SINGLE",
        downloaded="2026-09-06T00:00:00+00:00",
        tool="fasterq-dump",
        tool_version="3.0.0",
        compression="gzip",
        files=[{"name": f"{accession}.fastq.gz", "bytes": 100, "md5": "a" * 32, "reads": 5}],
        reads_per_mate=5,
        bases_total=500,
        ncbi={"spots": 5, "bases": 500, "size": 1000, "layout": "SINGLE", "files": []},
        completeness={"method": "spots", "ratio": 1.0, "verdict": "complete"},
    )


class TestStoreInitCommand:
    def test_command_properties(self):
        cmd = StoreInitCommand()
        assert cmd.name == "store_init"
        assert cmd.group == "Store"
        assert "store" in cmd.help.lower()

    def test_configure_parser_requires_data_root(self):
        cmd = StoreInitCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        with pytest.raises(SystemExit):
            parser.parse_args([])
        args = parser.parse_args(["--data-root", "/x"])
        assert args.data_root == "/x"
        assert args.project_name is None
        assert args.set_default is False

    def test_execute_creates_store_and_records_project(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreInitCommand().execute(_init_args(root, project_dir))
        assert rc == 0

        marker = read_marker(root)
        assert marker is not None

        paths = store_paths(root)
        with Catalog(paths) as cat:
            tables = {
                row["name"]
                for row in cat.conn.execute("SELECT name FROM sqlite_master WHERE type IN ('table', 'view')").fetchall()
            }
        assert {"datasets", "projects", "usage"} <= tables

        registry = load_registry(project_dir / "prepare_registry.json")
        assert registry.store["root"] == str(root.resolve())
        assert registry.store["mode"] == "symlink"
        assert registry.store["linked"] == []
        assert registry.project["name"]
        assert registry.project["id"]
        assert registry.project["path"] == str(project_dir.resolve())

    def test_execute_uses_project_name_flag(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreInitCommand().execute(_init_args(root, project_dir, project_name="Wolbachia"))
        assert rc == 0

        registry = load_registry(project_dir / "prepare_registry.json")
        assert registry.project["name"] == "Wolbachia"

    def test_execute_defaults_project_name_to_cwd_name(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        project_dir = tmp_path / "my-project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreInitCommand().execute(_init_args(root, project_dir))
        assert rc == 0

        registry = load_registry(project_dir / "prepare_registry.json")
        assert registry.project["name"] == "my-project"

    def test_execute_keeps_recorded_exports(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry = load_registry(project_dir / "prepare_registry.json")
        registry.project = {
            "exports": {"results_table": {"output": "results.tsv", "summary": {"rows": 1}}},
            "note": "kept",
        }
        save_registry(registry)

        rc = StoreInitCommand().execute(_init_args(root, project_dir))
        assert rc == 0

        registry = load_registry(project_dir / "prepare_registry.json")
        assert registry.project["exports"]["results_table"]["output"] == "results.tsv"
        assert registry.project["note"] == "kept"
        assert registry.project["id"]

    def test_execute_keeps_project_id_on_second_run(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        StoreInitCommand().execute(_init_args(root, project_dir))
        first_id = load_registry(project_dir / "prepare_registry.json").project["id"]

        StoreInitCommand().execute(_init_args(root, project_dir))
        second_id = load_registry(project_dir / "prepare_registry.json").project["id"]

        assert first_id == second_id

    def test_store_status_counts_project_after_init(self, tmp_path, monkeypatch, capsys):
        root = tmp_path / "store"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"

        assert StoreInitCommand().execute(_init_args(root, project_dir, registry=str(registry_path))) == 0

        rc = StoreStatusCommand().execute(_status_args(data_root=str(root), json=True, registry=str(registry_path)))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["projects"] == 1

    def test_execute_set_default_writes_config(self, tmp_path, monkeypatch):
        from prepare.store.resolve import read_config

        root = tmp_path / "store"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreInitCommand().execute(_init_args(root, project_dir, set_default=True))
        assert rc == 0
        assert read_config()["store"]["data_root"] == root.resolve().as_posix()

    def test_execute_without_set_default_does_not_write_config(self, tmp_path, monkeypatch):
        from prepare.store.resolve import read_config

        root = tmp_path / "store"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        StoreInitCommand().execute(_init_args(root, project_dir, set_default=False))
        assert read_config() == {}

    def test_gitignore_guard_appends_fastq_when_git_present(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        (project_dir / ".git").mkdir()
        monkeypatch.chdir(project_dir)

        with patch("prepare.cli.commands.store._shared.subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(args=[], returncode=0, stdout="", stderr="")
            rc = StoreInitCommand().execute(_init_args(root, project_dir))
        assert rc == 0

        gitignore = (project_dir / ".gitignore").read_text()
        assert "fastq/" in gitignore.splitlines()
        mock_run.assert_called_once()
        assert mock_run.call_args.args[0] == ["git", "ls-files", "fastq"]

    def test_gitignore_guard_skips_when_already_present(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        (project_dir / ".git").mkdir()
        (project_dir / ".gitignore").write_text("fastq/\nother\n")
        monkeypatch.chdir(project_dir)

        with patch("prepare.cli.commands.store._shared.subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(args=[], returncode=0, stdout="", stderr="")
            StoreInitCommand().execute(_init_args(root, project_dir))

        gitignore = (project_dir / ".gitignore").read_text()
        assert gitignore.count("fastq/") == 1

    def test_gitignore_guard_warns_when_fastq_tracked(self, tmp_path, monkeypatch, caplog):
        root = tmp_path / "store"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        (project_dir / ".git").mkdir()
        monkeypatch.chdir(project_dir)

        with patch("prepare.cli.commands.store._shared.subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                args=[], returncode=0, stdout="fastq/SRR1/SRR1.fastq\n", stderr=""
            )
            with caplog.at_level("WARNING"):
                StoreInitCommand().execute(_init_args(root, project_dir))

        assert any("git" in message.lower() and "fastq" in message.lower() for message in caplog.messages)

    def test_gitignore_guard_never_runs_when_no_git_dir(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        with patch("prepare.cli.commands.store._shared.subprocess.run") as mock_run:
            rc = StoreInitCommand().execute(_init_args(root, project_dir))
        assert rc == 0
        mock_run.assert_not_called()
        assert not (project_dir / ".gitignore").exists()


class TestStoreStatusCommand:
    def test_command_properties(self):
        cmd = StoreStatusCommand()
        assert cmd.name == "store_status"
        assert cmd.group == "Store"

    def test_no_store_configured_logs_hint_and_returns_1(self, tmp_path, monkeypatch, capsys, caplog):
        monkeypatch.chdir(tmp_path)
        rc = StoreStatusCommand().execute(_status_args(registry=str(tmp_path / "prepare_registry.json")))
        assert rc == 1
        assert "store_init" in caplog.text
        # Without --json the hint is an error on stderr; stdout carries no result.
        assert capsys.readouterr().out == ""

    def test_no_store_configured_with_json_prints_json_error(self, tmp_path, monkeypatch, capsys):
        monkeypatch.chdir(tmp_path)
        rc = StoreStatusCommand().execute(_status_args(json=True, registry=str(tmp_path / "prepare_registry.json")))
        payload = json.loads(capsys.readouterr().out)
        assert rc == 1
        assert payload["error"].startswith("No store configured")

    def test_reports_counts_bytes_and_projects(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        project_dir = tmp_path / "proj1"
        project_dir.mkdir()
        registry_file = project_dir / "prepare_registry.json"
        registry_file.write_text(json.dumps({"project": {"id": "proj1"}}))
        with catalog_write(paths) as cat:
            cat.upsert_project("proj1", "Wolbachia", str(project_dir), str(registry_file))
            cat.upsert_dataset(_sidecar("SRR1", state="complete"))
            cat.upsert_dataset(_sidecar("SRR2", state="partial"))
            cat.record_usage("SRR1", "proj1", "wMel", "downloaded")

        rc = StoreStatusCommand().execute(_status_args(data_root=str(root), json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["root"] == str(root.resolve())
        assert report["datasets"] == {"complete": 1, "partial": 1}
        assert report["bytes_total"] == 200
        assert report["projects"] == 1
        assert report["stale_projects"] == []
        assert "datasets_list" not in report

    def test_verbose_lists_datasets(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        with catalog_write(paths) as cat:
            cat.upsert_project("proj1", "Wolbachia", str(tmp_path / "proj1"), "reg1")
            cat.upsert_dataset(_sidecar("SRR1"))
            cat.record_usage("SRR1", "proj1", "wMel", "downloaded")

        rc = StoreStatusCommand().execute(_status_args(data_root=str(root), json=True, verbose=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["datasets_list"] == [{"accession": "SRR1", "state": "complete", "bytes": 100, "projects": 1}]

    def test_stale_project_detected_when_registry_missing(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        missing_path = tmp_path / "gone"
        missing_registry = tmp_path / "gone" / "prepare_registry.json"
        with catalog_write(paths) as cat:
            cat.upsert_project("proj1", "Wolbachia", str(missing_path), str(missing_registry))

        rc = StoreStatusCommand().execute(_status_args(data_root=str(root), json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert len(report["stale_projects"]) == 1
        row = report["stale_projects"][0]
        assert (row["project_id"], row["name"], row["registry"]) == ("proj1", "Wolbachia", str(missing_registry))
        assert row["reason"] == "registry missing"

    def test_stale_project_detected_when_registry_id_differs(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        project_dir = tmp_path / "proj1"
        project_dir.mkdir()
        registry_file = project_dir / "prepare_registry.json"
        registry_file.write_text(json.dumps({"project": {"id": "some-other-id"}}))
        with catalog_write(paths) as cat:
            cat.upsert_project("proj1", "Wolbachia", str(project_dir), str(registry_file))

        rc = StoreStatusCommand().execute(_status_args(data_root=str(root), json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert len(report["stale_projects"]) == 1
        row = report["stale_projects"][0]
        assert (row["project_id"], row["name"], row["registry"]) == ("proj1", "Wolbachia", str(registry_file))
        # A registry is there; it simply belongs to a different project now. Saying "registry
        # missing" for that case sends the reader looking for a file that exists.
        assert row["reason"] == "project id differs"

    def test_text_output_includes_store_header(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        with catalog_write(paths):
            pass

        rc = StoreStatusCommand().execute(_status_args(data_root=str(root)))
        out = capsys.readouterr().out

        assert rc == 0
        assert "Store" in out


class TestStoreReindexCommand:
    def test_command_properties(self):
        cmd = StoreReindexCommand()
        assert cmd.name == "store_reindex"
        assert cmd.group == "Store"

    def test_no_store_configured_returns_1(self, tmp_path, monkeypatch, capsys):
        monkeypatch.chdir(tmp_path)
        rc = StoreReindexCommand().execute(_reindex_args(registry=str(tmp_path / "prepare_registry.json")))
        assert rc == 1

    def test_rebuilds_from_sidecars_and_prints_count(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        for acc in ("SRR1", "SRR2"):
            sc_path = sidecar_path(paths, acc)
            write_sidecar(sc_path, _sidecar(acc))

        rc = StoreReindexCommand().execute(_reindex_args(data_root=str(root)))
        out = capsys.readouterr().out

        assert rc == 0
        assert "2" in out

        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1") is not None
            assert cat.get_dataset("SRR2") is not None

    def test_reindex_drops_removed_accessions(self, tmp_path):
        root = tmp_path / "store"
        paths = init_store(root)
        write_sidecar(sidecar_path(paths, "SRR1"), _sidecar("SRR1"))
        StoreReindexCommand().execute(_reindex_args(data_root=str(root)))

        with catalog_write(paths) as cat:
            cat.upsert_dataset(_sidecar("SRR-stale"))

        StoreReindexCommand().execute(_reindex_args(data_root=str(root)))

        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR-stale") is None
            assert cat.get_dataset("SRR1") is not None


def _write_fastq_gz(path, text="@r\nACGT\n+\nIIII\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def _sidecar_matching_disk(acc_dir, accession, state="complete", reads=5):
    """A sidecar whose file records match the FASTQ files actually on disk for this accession
    (`<accession>.fastq.gz` by default, but also a mate-suffixed name such as
    `<accession>_1.fastq.gz`), so a plain `store_verify` (no --md5) reports it as healthy."""
    file_paths = sorted(p for p in acc_dir.glob("*.fastq.gz") if not p.name.startswith("."))
    if not file_paths:
        file_paths = [acc_dir / f"{accession}.fastq.gz"]
    sidecar = _sidecar(accession, state=state)
    sidecar.files = [
        {"name": file_path.name, "bytes": file_path.stat().st_size, "md5": "ignored", "reads": reads}
        for file_path in file_paths
    ]
    return sidecar


class TestStoreAdoptCommand:
    def test_command_properties(self):
        cmd = StoreAdoptCommand()
        assert cmd.name == "store_adopt"
        assert cmd.group == "Store"

    def test_no_store_configured_returns_1(self, tmp_path, monkeypatch, capsys, caplog):
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreAdoptCommand().execute(_adopt_args(registry=str(project_dir / "prepare_registry.json")))

        assert rc == 1
        assert "store_init" in caplog.text
        assert capsys.readouterr().out == ""

    def test_adopts_moves_links_and_records_registry(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        init_store(root)
        project_dir = tmp_path / "project"
        (project_dir / "fastq" / "SRR1").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1" / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        monkeypatch.chdir(project_dir)

        registry_path = project_dir / "prepare_registry.json"
        rc = StoreAdoptCommand().execute(_adopt_args(data_root=str(root), registry=str(registry_path)))

        assert rc == 0
        assert (project_dir / "fastq" / "SRR1").is_symlink()

        paths = store_paths(root)
        assert sra_dir(paths, "SRR1").is_dir()
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1") is not None

        registry = load_registry(registry_path)
        assert registry.datasets["SRR1"]["download"]["state"] == "downloaded"
        assert registry.datasets["SRR1"]["download"]["source"] == "store"
        assert registry.store["linked"] == ["SRR1"]

    def test_adopt_copies_sidecar_completeness_and_records_usage(self, tmp_path, monkeypatch):
        """A newly adopted accession's registry record carries its sidecar's completeness
        verdict, and the store catalogue gets a 'linked' usage row for this project."""
        root = tmp_path / "store"
        init_store(root)
        project_dir = tmp_path / "project"
        (project_dir / "fastq" / "SRR1").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1" / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"

        assert StoreInitCommand().execute(_init_args(root, project_dir, registry=str(registry_path))) == 0

        rc = StoreAdoptCommand().execute(_adopt_args(data_root=str(root), registry=str(registry_path)))
        assert rc == 0

        registry = load_registry(registry_path)
        complete = registry.datasets["SRR1"]["download"]["complete"]
        # No metadata XML is present, so the sidecar's own verify against NCBI's spot count
        # is unverified; the registry record must carry exactly that verdict, not silently
        # invent one.
        assert complete["verdict"] == "unverified"
        assert complete["expected_spots"] is None
        assert complete["reads_r1"] == 1

        paths = store_paths(root)
        with Catalog(paths) as cat:
            row = cat.conn.execute(
                "SELECT stage FROM usage WHERE accession = ? AND project_id = ?",
                ("SRR1", registry.project["id"]),
            ).fetchone()
        assert row["stage"] == "linked"

    def test_catalog_failure_leaves_adopt_outcome_unchanged(self, tmp_path, monkeypatch):
        """A broken catalogue write never changes store_adopt's registry outcome or exit code."""
        root = tmp_path / "store"
        init_store(root)
        project_dir = tmp_path / "project"
        (project_dir / "fastq" / "SRR1").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1" / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"

        assert StoreInitCommand().execute(_init_args(root, project_dir, registry=str(registry_path))) == 0

        with patch("prepare.store.usage.catalog_write", side_effect=RuntimeError("locked")):
            rc = StoreAdoptCommand().execute(_adopt_args(data_root=str(root), registry=str(registry_path)))

        assert rc == 0
        registry = load_registry(registry_path)
        assert registry.datasets["SRR1"]["download"]["state"] == "downloaded"
        assert registry.datasets["SRR1"]["download"]["source"] == "store"
        assert registry.store["linked"] == ["SRR1"]

    def test_dry_run_changes_nothing(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        init_store(root)
        project_dir = tmp_path / "project"
        (project_dir / "fastq" / "SRR1").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1" / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        monkeypatch.chdir(project_dir)

        registry_path = project_dir / "prepare_registry.json"
        rc = StoreAdoptCommand().execute(_adopt_args(data_root=str(root), registry=str(registry_path), dry_run=True))

        assert rc == 0
        assert (project_dir / "fastq" / "SRR1" / "SRR1.fastq").is_file()
        assert not (project_dir / "fastq" / "SRR1").is_symlink()
        assert not registry_path.exists()

    def test_dry_run_reports_empty_and_failed_accessions(self, tmp_path, monkeypatch, capsys):
        """A dry run lists the folders it would leave alone as empty or failed, the same as a
        real run, rather than only the planned adoptions and conflicts."""
        root = tmp_path / "store"
        paths = init_store(root)
        project_dir = tmp_path / "project"
        (project_dir / "fastq" / "SRR1").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1" / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (project_dir / "fastq" / "SRR2").mkdir(parents=True)
        (project_dir / "fastq" / "SRR3").mkdir(parents=True)
        (project_dir / "fastq" / "SRR3" / "SRR3.fastq").write_text("@r\nACGT\n+\nIIII\n")
        store_copy = sra_dir(paths, "SRR3")
        store_copy.mkdir(parents=True)
        (store_copy / "SRR3.fastq").write_text("@r\nACGT\n+\nIIII\n")
        write_sidecar(sidecar_path(paths, "SRR3"), _sidecar("SRR3", state="failed"))
        monkeypatch.chdir(project_dir)

        registry_path = project_dir / "prepare_registry.json"
        rc = StoreAdoptCommand().execute(_adopt_args(data_root=str(root), registry=str(registry_path), dry_run=True))
        out = capsys.readouterr().out

        assert rc == 0
        assert "Would adopt 1 dataset(s)" in out
        assert "Empty folders, not adopted: SRR2" in out
        assert "SRR3" in [line for line in out.splitlines() if line.startswith("Failed")][0]

    def test_copy_mode_leaves_project_folder_untouched_and_unlinked(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        init_store(root)
        project_dir = tmp_path / "project"
        entry = project_dir / "fastq" / "SRR1"
        entry.mkdir(parents=True)
        (entry / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        monkeypatch.chdir(project_dir)

        registry_path = project_dir / "prepare_registry.json"
        rc = StoreAdoptCommand().execute(_adopt_args(data_root=str(root), registry=str(registry_path), move=False))

        assert rc == 0
        # --copy: the project's own folder is left exactly as it was, not linked.
        assert entry.is_dir() and not entry.is_symlink()
        assert (entry / "SRR1.fastq").is_file()

        paths = store_paths(root)
        assert sra_dir(paths, "SRR1").is_dir()

        # The project's own download record is untouched (--copy never points it at the
        # store), but the registry does now get a project identity, so the "copied" usage
        # recorded below has a project to belong to and store_gc does not see this dataset
        # as unused just because it was never linked.
        registry = load_registry(registry_path)
        assert registry.project.get("id")
        assert "SRR1" not in registry.datasets
        assert not registry.store.get("linked")
        with Catalog(paths) as cat:
            row = cat.conn.execute(
                "SELECT stage FROM usage WHERE accession = ? AND project_id = ?",
                ("SRR1", registry.project["id"]),
            ).fetchone()
        assert row is not None and row["stage"] == "copied"

    def test_copy_mode_dedup_does_not_record_link(self, tmp_path, monkeypatch):
        """A --copy run that finds the store already holds an identical copy (a dedup) must
        not be recorded as linked: the project's folder is left exactly as it was, a real
        directory, the same as a fresh --copy adoption, so no store-sourced download record
        or store["linked"] entry should appear."""
        root = tmp_path / "store"
        init_store(root)
        project_dir = tmp_path / "project"
        entry = project_dir / "fastq" / "SRR1"
        entry.mkdir(parents=True)
        (entry / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"

        assert StoreInitCommand().execute(_init_args(root, project_dir, registry=str(registry_path))) == 0

        first_rc = StoreAdoptCommand().execute(
            _adopt_args(data_root=str(root), registry=str(registry_path), move=False)
        )
        assert first_rc == 0

        second_rc = StoreAdoptCommand().execute(
            _adopt_args(data_root=str(root), registry=str(registry_path), move=False)
        )
        assert second_rc == 0

        assert entry.is_dir() and not entry.is_symlink()
        assert (entry / "SRR1.fastq").is_file()

        registry = load_registry(registry_path)
        assert not registry.store.get("linked")
        assert "SRR1" not in registry.datasets

    def test_copy_adoption_records_copied_usage_so_gc_keeps_it(self, tmp_path, monkeypatch, capsys):
        """A --copy adoption must record a 'copied' usage row for the project, so store_gc does
        not treat the store's copy as unused just because the project kept its own folder too
        (see test_copy_mode_dedup_does_not_record_link for the same fixture shape)."""
        root = tmp_path / "store"
        init_store(root)
        project_dir = tmp_path / "project"
        entry = project_dir / "fastq" / "SRR1"
        entry.mkdir(parents=True)
        (entry / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"

        assert StoreInitCommand().execute(_init_args(root, project_dir, registry=str(registry_path))) == 0

        rc = StoreAdoptCommand().execute(_adopt_args(data_root=str(root), registry=str(registry_path), move=False))
        assert rc == 0

        registry = load_registry(registry_path)
        paths = store_paths(root)
        with Catalog(paths) as cat:
            row = cat.conn.execute(
                "SELECT stage FROM usage WHERE accession = ? AND project_id = ?",
                ("SRR1", registry.project["id"]),
            ).fetchone()
        assert row is not None and row["stage"] == "copied"

        capsys.readouterr()
        gc_rc = StoreGcCommand().execute(_gc_args(data_root=str(root), dry_run=True, json=True))
        report = json.loads(capsys.readouterr().out)
        assert gc_rc == 0
        # The recorded usage keeps this dataset out of gc's removal candidates.
        assert report["datasets"] == []

    def test_signal_between_accessions_stops_before_the_next_one(self, tmp_path, monkeypatch):
        """A stop noticed after SRR1 is adopted ends the run before SRR2 is ever staged."""
        import importlib

        # prepare.store re-exports the ``adopt`` function under the same name, which
        # shadows the submodule on ``prepare.store.adopt``; import_module bypasses that.
        adopt_mod = importlib.import_module("prepare.store.adopt")

        root = tmp_path / "store"
        init_store(root)
        project_dir = tmp_path / "project"
        for acc in ("SRR1", "SRR2"):
            entry = project_dir / "fastq" / acc
            entry.mkdir(parents=True)
            (entry / f"{acc}.fastq").write_text("@r\nACGT\n+\nIIII\n")
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"

        cmd = StoreAdoptCommand()
        args = _adopt_args(data_root=str(root), registry=str(registry_path))

        original = adopt_mod._adopt_one
        calls = []

        def spy(accession, *a, **kw):
            result = original(accession, *a, **kw)
            calls.append(accession)
            if len(calls) == 1:
                args._termination.stop.set()
            return result

        with patch.object(adopt_mod, "_adopt_one", side_effect=spy):
            rc = cmd.run(args)

        assert rc == 130
        assert calls == ["SRR1"]

        paths = store_paths(root)
        assert sra_dir(paths, "SRR1").is_dir()
        assert not sra_dir(paths, "SRR2").exists()
        assert (project_dir / "fastq" / "SRR1").is_symlink()
        assert (project_dir / "fastq" / "SRR2").is_dir()
        assert not (project_dir / "fastq" / "SRR2").is_symlink()


class TestStoreVerifyCommand:
    def test_command_properties(self):
        cmd = StoreVerifyCommand()
        assert cmd.name == "store_verify"
        assert cmd.group == "Store"

    def test_no_store_configured_returns_1(self, tmp_path, monkeypatch, capsys):
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreVerifyCommand().execute(_verify_args(registry=str(project_dir / "prepare_registry.json")))
        assert rc == 1

    def test_verify_healthy_dataset_reports_ok(self, tmp_path, monkeypatch, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz")
        sidecar = _sidecar_matching_disk(acc_dir, "SRR1")
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreVerifyCommand().execute(
            _verify_args(data_root=str(root), registry=str(project_dir / "prepare_registry.json"))
        )
        out = capsys.readouterr().out

        assert rc == 0
        assert "SRR1" in out
        assert "ok" in out

    def test_verify_flags_truncated_dataset_and_fix_state_updates_sidecar(self, tmp_path, monkeypatch, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        # Only 1 read on disk, but the sidecar's recorded NCBI spot count says there should be 5.
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz", text="@r\nACGT\n+\nIIII\n")
        sidecar = _sidecar("SRR1", state="complete")
        sidecar.files = [
            {
                "name": "SRR1.fastq.gz",
                "bytes": acc_dir.joinpath("SRR1.fastq.gz").stat().st_size,
                "md5": "ignored",
                "reads": 1,
            }
        ]
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreVerifyCommand().execute(
            _verify_args(
                data_root=str(root), registry=str(project_dir / "prepare_registry.json"), spots=True, fix_state=True
            )
        )
        out = capsys.readouterr().out

        assert rc == 1
        assert "truncated" in out

        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "partial"
        assert fixed.completeness["verdict"] == "truncated"

        with Catalog(paths) as cat:
            row = cat.get_dataset("SRR1")
            assert row["state"] == "partial"

    def test_verify_flags_bad_md5_as_failed_regardless_of_spots(self, tmp_path, monkeypatch, capsys):
        """A bytes/md5 mismatch always wins over the spots check: a corrupted file can still
        happen to contain the right number of reads, so --fix-state must mark it failed, not
        complete."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz", text="@r\nACGT\n+\nIIII\n")
        sidecar = _sidecar("SRR1", state="complete")
        sidecar.files = [
            {
                "name": "SRR1.fastq.gz",
                "bytes": acc_dir.joinpath("SRR1.fastq.gz").stat().st_size,
                # Wrong md5 on purpose; read count/spots still match (ncbi.spots defaults to 5).
                "md5": "0" * 32,
                "reads": 5,
            }
        ]
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreVerifyCommand().execute(
            _verify_args(
                data_root=str(root),
                registry=str(project_dir / "prepare_registry.json"),
                md5=True,
                spots=True,
                fix_state=True,
            )
        )
        out = capsys.readouterr().out

        assert rc == 1
        assert "corrupt" in out

        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "failed"
        assert fixed.error and "md5 mismatch" in fixed.error

        with Catalog(paths) as cat:
            row = cat.get_dataset("SRR1")
            assert row["state"] == "failed"

    def test_verify_specific_accessions_only(self, tmp_path, monkeypatch, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        for acc in ("SRR1", "SRR2"):
            acc_dir = sra_dir(paths, acc)
            _write_fastq_gz(acc_dir / f"{acc}.fastq.gz")
            sidecar = _sidecar_matching_disk(acc_dir, acc)
            write_sidecar(sidecar_path(paths, acc), sidecar)
            with catalog_write(paths) as cat:
                cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreVerifyCommand().execute(
            _verify_args(
                accessions=["SRR2"], data_root=str(root), registry=str(project_dir / "prepare_registry.json")
            )
        )
        out = capsys.readouterr().out

        assert rc == 0
        assert "SRR2" in out
        assert "SRR1" not in out

    def test_missing_dataset_fails(self, tmp_path, monkeypatch, capsys):
        root = tmp_path / "store"
        init_store(root)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreVerifyCommand().execute(
            _verify_args(
                accessions=["SRR-ghost"],
                data_root=str(root),
                registry=str(project_dir / "prepare_registry.json"),
            )
        )
        out = capsys.readouterr().out

        assert rc == 1
        assert "missing" in out

    def test_corrupt_gzip_reported_as_corrupt_without_crashing(self, tmp_path, monkeypatch, capsys):
        """A truncated/corrupt gzip file makes read-count verification raise inside the data
        layer; --spots must catch that and report it, not crash the whole command."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        acc_dir.mkdir(parents=True)
        # Looks like a gzip file (magic bytes) but is not valid gzip data.
        (acc_dir / "SRR1.fastq.gz").write_bytes(b"\x1f\x8b\x00not-really-gzip")
        sidecar = _sidecar_matching_disk(acc_dir, "SRR1")
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreVerifyCommand().execute(
            _verify_args(data_root=str(root), registry=str(project_dir / "prepare_registry.json"), spots=True)
        )
        out = capsys.readouterr().out

        assert rc == 1
        assert "corrupt" in out

    def test_spots_from_store_metadata_xml(self, tmp_path, monkeypatch):
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        _write_fastq_gz(acc_dir / "SRR1_1.fastq.gz", "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n")
        sc = _sidecar_matching_disk(acc_dir, "SRR1", state="failed", reads=2)
        sc.ncbi = {}
        sc.error = "._SRR1_1.fastq.gz: missing"
        write_sidecar(sidecar_path(paths, "SRR1"), sc)
        (paths.metadata / "SRR1_metadata.xml").write_text(
            "<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><RUN_SET>"
            '<RUN accession="SRR1" total_spots="2" total_bases="8" size="100">'
            '<SRAFiles><SRAFile filename="SRR1_1.fastq.gz" md5="00"/></SRAFiles></RUN></RUN_SET>'
            "</EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>"
        )
        rc = StoreVerifyCommand().execute(_verify_args(data_root=str(paths.root), spots=True, fix_state=True))
        assert rc == 0
        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "complete"
        assert fixed.completeness["verdict"] == "complete"
        assert fixed.ncbi.get("spots") == 2
        assert fixed.error is None

    def test_fix_state_promotes_readable_dataset_without_spots(self, tmp_path, monkeypatch):
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        _write_fastq_gz(acc_dir / "SRR1_1.fastq.gz")
        sc = _sidecar_matching_disk(acc_dir, "SRR1", state="failed", reads=1)
        sc.ncbi = {}
        sc.error = "old error"
        write_sidecar(sidecar_path(paths, "SRR1"), sc)
        rc = StoreVerifyCommand().execute(_verify_args(data_root=str(paths.root), spots=True, fix_state=True))
        assert rc == 0
        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "complete" and fixed.error is None
        assert fixed.completeness["verdict"] == "unverified"

    def test_rescan_rebuilds_file_list_from_disk(self, tmp_path, monkeypatch):
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        _write_fastq_gz(acc_dir / "SRR1_1.fastq.gz")
        sc = _sidecar_matching_disk(acc_dir, "SRR1", state="failed", reads=1)
        sc.files.insert(0, {"name": "._SRR1_1.fastq.gz", "bytes": 4096, "md5": "x", "reads": None})
        # The recorded spot count matches the one read on disk: --fix-state compares against it.
        sc.ncbi = {"spots": 1}
        write_sidecar(sidecar_path(paths, "SRR1"), sc)
        rc = StoreVerifyCommand().execute(_verify_args(data_root=str(paths.root), rescan=True, fix_state=True))
        assert rc == 0
        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert [f["name"] for f in fixed.files] == ["SRR1_1.fastq.gz"]
        assert fixed.state == "complete"

    def test_fix_state_without_spots_requires_a_read_through(self, tmp_path, monkeypatch, capsys):
        """A file that matches its recorded size can still be a corrupt gzip stream that was
        never opened; --fix-state without --spots must read it before promoting to complete,
        not promote on size alone."""
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        # Looks like a gzip file (magic bytes) but is not valid gzip data.
        (acc_dir / "SRR1.fastq.gz").write_bytes(b"\x1f\x8b\x00not-really-gzip")
        sc = _sidecar_matching_disk(acc_dir, "SRR1", state="failed", reads=1)
        sc.error = "old error"
        write_sidecar(sidecar_path(paths, "SRR1"), sc)

        rc = StoreVerifyCommand().execute(_verify_args(data_root=str(paths.root), fix_state=True))
        out = capsys.readouterr().out

        assert rc == 1
        assert "corrupt" in out
        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "failed"
        assert fixed.error and "SRR1.fastq.gz" in fixed.error

    def test_fix_state_never_promotes_a_dataset_with_no_files_recorded(self, tmp_path, monkeypatch, capsys):
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        sc = _sidecar("SRR1", state="failed")
        sc.files = []
        sc.error = "old error"
        write_sidecar(sidecar_path(paths, "SRR1"), sc)

        rc = StoreVerifyCommand().execute(_verify_args(data_root=str(paths.root), fix_state=True))
        out = capsys.readouterr().out

        # A dataset --fix-state leaves (or marks) failed is never reported as ok.
        assert rc == 1
        assert "corrupt" in out
        assert "failed" in out
        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "failed"
        assert fixed.error == "no files recorded"

    def test_fix_state_without_spots_keeps_a_known_truncated_dataset_partial(self, tmp_path, monkeypatch, capsys):
        """The sidecar already records NCBI's spot count, so --fix-state alone must compare
        against it rather than promote a truncated dataset to complete as unverified."""
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        acc_dir = paths.sra / "SRR1"
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz", "@r\nACGT\n+\nIIII\n" * 10)
        sc = _sidecar_matching_disk(acc_dir, "SRR1", state="partial", reads=10)
        sc.ncbi = {"spots": 100}
        sc.completeness = {"method": "spots", "ratio": 0.1, "verdict": "truncated"}
        write_sidecar(sidecar_path(paths, "SRR1"), sc)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sc)

        rc = StoreVerifyCommand().execute(_verify_args(data_root=str(paths.root), fix_state=True))
        out = capsys.readouterr().out

        assert rc == 1
        assert "truncated" in out
        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "partial"
        assert fixed.completeness == {"method": "spots", "ratio": 0.1, "verdict": "truncated"}
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1")["state"] == "partial"

    def test_fix_state_without_spots_promotes_a_complete_dataset_by_its_known_spot_count(self, tmp_path, monkeypatch):
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        acc_dir = paths.sra / "SRR1"
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz", "@r\nACGT\n+\nIIII\n" * 10)
        sc = _sidecar_matching_disk(acc_dir, "SRR1", state="failed", reads=10)
        sc.ncbi = {"spots": 10}
        sc.error = "old error"
        write_sidecar(sidecar_path(paths, "SRR1"), sc)

        rc = StoreVerifyCommand().execute(_verify_args(data_root=str(paths.root), fix_state=True))

        assert rc == 0
        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "complete" and fixed.error is None
        assert fixed.completeness == {"method": "spots", "ratio": 1.0, "verdict": "complete"}

    def test_fix_state_with_spots_still_reads_the_mate_2_file(self, tmp_path, monkeypatch, capsys):
        """--spots counts mate 1 only, so --fix-state must read mate 2 before promoting."""
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        acc_dir = paths.sra / "SRR1"
        _write_fastq_gz(acc_dir / "SRR1_1.fastq.gz", "@r\nACGT\n+\nIIII\n" * 10)
        (acc_dir / "SRR1_2.fastq.gz").write_bytes(b"\x1f\x8b\x00not-really-gzip")
        sc = _sidecar_matching_disk(acc_dir, "SRR1", state="failed", reads=10)
        sc.ncbi = {"spots": 10}
        write_sidecar(sidecar_path(paths, "SRR1"), sc)

        rc = StoreVerifyCommand().execute(_verify_args(data_root=str(paths.root), spots=True, fix_state=True))
        out = capsys.readouterr().out

        assert rc == 1
        assert "corrupt" in out
        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "failed"
        assert fixed.error and "SRR1_2.fastq.gz" in fixed.error

    def test_fix_state_with_a_known_spot_count_still_reads_the_mate_2_file(self, tmp_path, monkeypatch, capsys):
        """The spot comparison counts mate 1 only; a corrupt mate 2 file must still be found."""
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        acc_dir = paths.sra / "SRR1"
        _write_fastq_gz(acc_dir / "SRR1_1.fastq.gz", "@r\nACGT\n+\nIIII\n" * 10)
        (acc_dir / "SRR1_2.fastq.gz").write_bytes(b"\x1f\x8b\x00not-really-gzip")
        sc = _sidecar_matching_disk(acc_dir, "SRR1", state="failed", reads=10)
        sc.ncbi = {"spots": 10}
        write_sidecar(sidecar_path(paths, "SRR1"), sc)

        rc = StoreVerifyCommand().execute(_verify_args(data_root=str(paths.root), fix_state=True))
        out = capsys.readouterr().out

        assert rc == 1
        assert "corrupt" in out
        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "failed"
        assert fixed.error and "SRR1_2.fastq.gz" in fixed.error

    def test_rescan_no_fastq_left_reports_missing_and_leaves_sidecar_untouched(self, tmp_path, monkeypatch, capsys):
        """The accession's folder (and its sidecar) still exist, but every FASTQ file was
        removed from it by hand; --rescan must not silently rebuild an empty file list and let
        --fix-state promote that to "complete", so it reports the dataset missing instead and
        leaves the sidecar exactly as it was."""
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        _write_fastq_gz(acc_dir / "SRR1_1.fastq.gz")
        sc = _sidecar_matching_disk(acc_dir, "SRR1", state="complete", reads=1)
        write_sidecar(sidecar_path(paths, "SRR1"), sc)
        original_files = [dict(f) for f in sc.files]
        (acc_dir / "SRR1_1.fastq.gz").unlink()

        rc = StoreVerifyCommand().execute(_verify_args(data_root=str(paths.root), rescan=True, fix_state=True))
        out = capsys.readouterr().out

        assert rc == 1
        assert "missing" in out
        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "complete"
        assert fixed.files == original_files

    def test_fix_state_skipped_when_the_dataset_is_locked(self, tmp_path, monkeypatch, capsys, caplog):
        """--fix-state's write-back takes the accession's lock without waiting; another run
        already holding it (a real dataset_lock here, not a fabricated lock file) means the
        sidecar is left exactly as it was, and the accession is logged as in use."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz", text="@r\nACGT\n+\nIIII\n")
        sidecar = _sidecar("SRR1", state="complete")
        sidecar.files = [
            {
                "name": "SRR1.fastq.gz",
                "bytes": acc_dir.joinpath("SRR1.fastq.gz").stat().st_size,
                "md5": "ignored",
                "reads": 1,
            }
        ]
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        locked = threading.Event()
        release = threading.Event()

        def _hold():
            with dataset_lock(paths, "SRR1"):
                locked.set()
                assert release.wait(timeout=5), "the main thread never released the holder"

        holder = threading.Thread(target=_hold)
        holder.start()
        try:
            assert locked.wait(timeout=5), "the holder thread never took the lock"
            with caplog.at_level(logging.WARNING):
                rc = StoreVerifyCommand().execute(
                    _verify_args(
                        data_root=str(root),
                        registry=str(project_dir / "prepare_registry.json"),
                        spots=True,
                        fix_state=True,
                    )
                )
        finally:
            release.set()
            holder.join(timeout=5)
        out = capsys.readouterr().out

        assert rc == 1  # the truncated verdict is still reported; only the write-back is skipped
        assert "truncated" in out
        assert any("in use, not updated" in record.message for record in caplog.records)
        unchanged = read_sidecar(sidecar_path(paths, "SRR1"))
        assert unchanged.state == "complete"
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1")["state"] == "complete"

    def test_fix_state_skipped_when_the_sidecar_changed_since_the_check(self, tmp_path, caplog):
        """A sidecar replaced by another run between _verify_one's read and the locked
        write-back (patched in here to stand in for that race) must not be clobbered by a
        write-back computed from the stale copy this run read."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz", text="@r\nACGT\n+\nIIII\n")
        sidecar = _sidecar("SRR1", state="complete")
        sidecar.files = [
            {
                "name": "SRR1.fastq.gz",
                "bytes": acc_dir.joinpath("SRR1.fastq.gz").stat().st_size,
                "md5": "ignored",
                "reads": 1,
            }
        ]
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        cmd = StoreVerifyCommand()
        result = cmd._verify_one("SRR1", paths, check_md5=False, check_spots=True, metadata_dirs=[], rescan=False)

        # Another run replaces the sidecar between the check above and the write-back below.
        replaced = _sidecar("SRR1", state="failed")
        replaced.error = "replaced by another run"
        write_sidecar(sidecar_path(paths, "SRR1"), replaced)

        with caplog.at_level(logging.WARNING):
            cmd._fix_state_locked(result, paths)

        assert any("changed since the check" in record.message for record in caplog.records)
        current = read_sidecar(sidecar_path(paths, "SRR1"))
        assert current.state == "failed"
        assert current.error == "replaced by another run"

    def test_fix_state_not_skipped_when_nothing_changed(self, tmp_path):
        """The baseline: with no lock held and no change to the sidecar, _fix_state_locked
        still applies the write-back exactly as the unlocked _fix_state would."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz", text="@r\nACGT\n+\nIIII\n")
        sidecar = _sidecar("SRR1", state="complete")
        sidecar.files = [
            {
                "name": "SRR1.fastq.gz",
                "bytes": acc_dir.joinpath("SRR1.fastq.gz").stat().st_size,
                "md5": "ignored",
                "reads": 1,
            }
        ]
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        cmd = StoreVerifyCommand()
        result = cmd._verify_one("SRR1", paths, check_md5=False, check_spots=True, metadata_dirs=[], rescan=False)
        cmd._fix_state_locked(result, paths)

        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "partial"
        assert fixed.completeness["verdict"] == "truncated"

    def test_no_store_configured_with_json_prints_json_error(self, tmp_path, monkeypatch, capsys):
        monkeypatch.chdir(tmp_path)
        rc = StoreVerifyCommand().execute(_verify_args(registry=str(tmp_path / "prepare_registry.json"), json=True))
        payload = json.loads(capsys.readouterr().out)
        assert rc == 1
        assert payload["error"].startswith("No store configured")

    def test_json_reports_corrupt_dataset_with_counts_and_no_sidecar_object(self, tmp_path, monkeypatch, capsys):
        """A size mismatch reports one "corrupt" dataset; --json prints exactly one document
        whose dataset row drops the Sidecar object and ncbi_found lookup _json_row is meant
        to drop, and whose counts tally the verdict."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz", text="@r\nACGT\n+\nIIII\n")
        sidecar = _sidecar("SRR1", state="complete")
        sidecar.files = [
            {
                "name": "SRR1.fastq.gz",
                # Wrong recorded size on purpose: the file on disk does not match.
                "bytes": acc_dir.joinpath("SRR1.fastq.gz").stat().st_size + 1,
                "md5": "ignored",
                "reads": 1,
            }
        ]
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreVerifyCommand().execute(
            _verify_args(data_root=str(root), registry=str(project_dir / "prepare_registry.json"), json=True)
        )
        doc = json.loads(capsys.readouterr().out)

        assert rc == 1
        assert doc["root"] == str(root.resolve())
        assert doc["checks"] == {"md5": False, "spots": False, "rescan": False, "fix_state": False}
        assert doc["counts"] == {"corrupt": 1}
        assert doc["fixed"] == []
        assert len(doc["datasets"]) == 1
        row = doc["datasets"][0]
        assert row["accession"] == "SRR1"
        assert row["verdict"] == "corrupt"
        assert "sidecar" not in row
        assert "ncbi_found" not in row
        assert "size mismatch" in row["mismatch_detail"]

    def test_fix_state_json_records_before_and_after_states(self, tmp_path, monkeypatch, capsys):
        """--fix-state --json: the document's "fixed" list and each dataset's nested "fix"
        both carry the state the sidecar had before the write-back and the state it was
        promoted/demoted to, matching what was actually written to disk."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        # Only 1 read on disk, but the sidecar's recorded NCBI spot count says there should be 5.
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz", text="@r\nACGT\n+\nIIII\n")
        sidecar = _sidecar("SRR1", state="complete")
        sidecar.files = [
            {
                "name": "SRR1.fastq.gz",
                "bytes": acc_dir.joinpath("SRR1.fastq.gz").stat().st_size,
                "md5": "ignored",
                "reads": 1,
            }
        ]
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreVerifyCommand().execute(
            _verify_args(
                data_root=str(root),
                registry=str(project_dir / "prepare_registry.json"),
                spots=True,
                fix_state=True,
                json=True,
            )
        )
        doc = json.loads(capsys.readouterr().out)

        assert rc == 1
        assert doc["checks"]["fix_state"] is True
        assert doc["counts"] == {"truncated": 1}
        assert doc["fixed"] == [
            {"accession": "SRR1", "state_before": "complete", "state_after": "partial", "error": None}
        ]
        row = doc["datasets"][0]
        assert row["state"] == "partial"
        assert row["state_before"] == "complete"
        assert row["fix"] == {
            "action": "updated",
            "state_before": "complete",
            "state_after": "partial",
            "error": None,
        }

        fixed = read_sidecar(sidecar_path(paths, "SRR1"))
        assert fixed.state == "partial"

    def test_fix_state_json_reports_skipped_in_use_when_locked(self, tmp_path, monkeypatch, capsys, caplog):
        """--fix-state --json: an accession another run holds the lock for is reported in the
        document as "skipped-in-use", not folded into "fixed" (nothing was actually written)."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz", text="@r\nACGT\n+\nIIII\n")
        sidecar = _sidecar("SRR1", state="complete")
        sidecar.files = [
            {
                "name": "SRR1.fastq.gz",
                "bytes": acc_dir.joinpath("SRR1.fastq.gz").stat().st_size,
                "md5": "ignored",
                "reads": 1,
            }
        ]
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        locked = threading.Event()
        release = threading.Event()

        def _hold():
            with dataset_lock(paths, "SRR1"):
                locked.set()
                assert release.wait(timeout=5), "the main thread never released the holder"

        holder = threading.Thread(target=_hold)
        holder.start()
        try:
            assert locked.wait(timeout=5), "the holder thread never took the lock"
            with caplog.at_level(logging.WARNING):
                rc = StoreVerifyCommand().execute(
                    _verify_args(
                        data_root=str(root),
                        registry=str(project_dir / "prepare_registry.json"),
                        spots=True,
                        fix_state=True,
                        json=True,
                    )
                )
        finally:
            release.set()
            holder.join(timeout=5)
        doc = json.loads(capsys.readouterr().out)

        assert rc == 1
        assert doc["fixed"] == []
        row = doc["datasets"][0]
        assert row["fix"]["action"] == "skipped-in-use"
        assert row["fix"]["state_before"] == "complete"
        assert row["fix"]["state_after"] == "complete"

        unchanged = read_sidecar(sidecar_path(paths, "SRR1"))
        assert unchanged.state == "complete"


class TestStoreLinkCommand:
    def test_command_properties(self):
        cmd = StoreLinkCommand()
        assert cmd.name == "store_link"
        assert cmd.group == "Store"

    def test_no_store_configured_returns_1(self, tmp_path, monkeypatch, capsys):
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreLinkCommand().execute(_link_args(["SRR1"], registry=str(project_dir / "prepare_registry.json")))
        assert rc == 1

    def test_links_accession_and_updates_registry(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz")
        sidecar = _sidecar("SRR1")
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"

        rc = StoreLinkCommand().execute(
            _link_args(
                ["SRR1"], data_root=str(root), registry=str(registry_path), fastq_folder=str(project_dir / "fastq")
            )
        )

        assert rc == 0
        assert (project_dir / "fastq" / "SRR1").is_symlink()

        registry = load_registry(registry_path)
        assert registry.datasets["SRR1"]["download"]["state"] == "downloaded"
        assert registry.datasets["SRR1"]["download"]["source"] == "store"
        assert registry.store["linked"] == ["SRR1"]

    def test_link_and_unlink_clear_an_earlier_download_time(self, tmp_path, monkeypatch):
        """A linked or removed dataset was not downloaded by this project: no download time is kept."""
        from prepare.data.registry import record_download, registry_transaction
        from prepare.data.registry_timing import set_download_timing

        root = tmp_path / "store"
        paths = init_store(root)
        _write_fastq_gz(sra_dir(paths, "SRR1") / "SRR1.fastq.gz")
        sidecar = _sidecar("SRR1")
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"
        fastq_dir = project_dir / "fastq"

        def timed_failure():
            with registry_transaction(registry_path) as reg:
                record_download(reg, "SRR1", "failed", fastq_dir, "network: reset")
                set_download_timing(reg, "SRR1", "2026-10-01T10:00:00+00:00", 42.0)

        timed_failure()
        link_args = _link_args(["SRR1"], data_root=str(root), registry=str(registry_path), fastq_folder=str(fastq_dir))
        assert StoreLinkCommand().execute(link_args) == 0
        assert "seconds" not in load_registry(registry_path).datasets["SRR1"]["download"]

        with registry_transaction(registry_path) as reg:
            set_download_timing(reg, "SRR1", "2026-10-01T10:00:00+00:00", 42.0)
        unlink_args = _unlink_args(["SRR1"], registry=str(registry_path), fastq_folder=str(fastq_dir))
        assert StoreUnlinkCommand().execute(unlink_args) == 0
        download = load_registry(registry_path).datasets["SRR1"]["download"]
        assert download["state"] == "missing" and "seconds" not in download and "started" not in download

    def test_link_copies_sidecar_completeness_and_records_usage(self, tmp_path, monkeypatch):
        """A linked accession's registry record carries the store sidecar's completeness
        verdict, and the store catalogue gets a 'linked' usage row for this project."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz")
        sidecar = _sidecar("SRR1")
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"

        assert StoreInitCommand().execute(_init_args(root, project_dir, registry=str(registry_path))) == 0

        rc = StoreLinkCommand().execute(
            _link_args(
                ["SRR1"], data_root=str(root), registry=str(registry_path), fastq_folder=str(project_dir / "fastq")
            )
        )
        assert rc == 0

        registry = load_registry(registry_path)
        complete = registry.datasets["SRR1"]["download"]["complete"]
        # Recomputed from the sidecar's read and spot counts (sra_verdicts.store_verdict), the same
        # verdict download_sra records for a link.
        assert complete == {"method": "spots", "verdict": "complete", "ratio": 1.0, "expected_spots": 5, "reads_r1": 5}

        with Catalog(paths) as cat:
            row = cat.conn.execute(
                "SELECT stage FROM usage WHERE accession = ? AND project_id = ?",
                ("SRR1", registry.project["id"]),
            ).fetchone()
        assert row["stage"] == "linked"

    def test_catalog_failure_leaves_link_outcome_unchanged(self, tmp_path, monkeypatch):
        """A broken catalogue write never changes store_link's registry outcome or exit code."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz")
        sidecar = _sidecar("SRR1")
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"

        assert StoreInitCommand().execute(_init_args(root, project_dir, registry=str(registry_path))) == 0

        with patch("prepare.store.usage.catalog_write", side_effect=RuntimeError("locked")):
            rc = StoreLinkCommand().execute(
                _link_args(
                    ["SRR1"],
                    data_root=str(root),
                    registry=str(registry_path),
                    fastq_folder=str(project_dir / "fastq"),
                )
            )

        assert rc == 0
        registry = load_registry(registry_path)
        assert registry.datasets["SRR1"]["download"]["state"] == "downloaded"
        assert registry.datasets["SRR1"]["download"]["source"] == "store"
        assert registry.store["linked"] == ["SRR1"]

    def test_link_failure_reports_nonzero(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        init_store(root)  # no SRR1 dataset in the store

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"

        rc = StoreLinkCommand().execute(
            _link_args(
                ["SRR1"], data_root=str(root), registry=str(registry_path), fastq_folder=str(project_dir / "fastq")
            )
        )

        assert rc == 1
        assert not registry_path.exists()


class TestStoreUnlinkCommand:
    def test_command_properties(self):
        cmd = StoreUnlinkCommand()
        assert cmd.name == "store_unlink"
        assert cmd.group == "Store"

    def test_unlinks_symlink_and_updates_registry(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        _write_fastq_gz(acc_dir / "SRR1.fastq.gz")
        sidecar = _sidecar("SRR1")
        write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(sidecar)

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"
        fastq_dir = project_dir / "fastq"

        StoreLinkCommand().execute(
            _link_args(["SRR1"], data_root=str(root), registry=str(registry_path), fastq_folder=str(fastq_dir))
        )

        rc = StoreUnlinkCommand().execute(
            _unlink_args(["SRR1"], registry=str(registry_path), fastq_folder=str(fastq_dir))
        )

        assert rc == 0
        assert not (fastq_dir / "SRR1").exists()

        registry = load_registry(registry_path)
        assert registry.datasets["SRR1"]["download"]["state"] == "missing"
        assert registry.store["linked"] == []

    def test_refuses_a_real_directory(self, tmp_path, monkeypatch):
        project_dir = tmp_path / "project"
        (project_dir / "fastq" / "SRR1").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1" / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"
        fastq_dir = project_dir / "fastq"

        rc = StoreUnlinkCommand().execute(
            _unlink_args(["SRR1"], registry=str(registry_path), fastq_folder=str(fastq_dir))
        )

        assert rc == 1
        assert (fastq_dir / "SRR1").is_dir() and not (fastq_dir / "SRR1").is_symlink()
        assert (fastq_dir / "SRR1" / "SRR1.fastq").is_file()


class TestStoreUsageCommand:
    def test_command_properties(self):
        cmd = StoreUsageCommand()
        assert cmd.name == "store_usage"
        assert cmd.group == "Store"

    def test_configure_parser_requires_one_selector(self):
        cmd = StoreUsageCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_configure_parser_rejects_two_selectors(self):
        cmd = StoreUsageCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        with pytest.raises(SystemExit):
            parser.parse_args(["--accession", "SRR1", "--unused"])

    def test_no_store_configured_returns_1(self, tmp_path, monkeypatch, capsys):
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreUsageCommand().execute(
            _usage_args(unused=True, registry=str(project_dir / "prepare_registry.json"))
        )
        assert rc == 1

    def test_no_store_configured_with_json_prints_json_error(self, tmp_path, monkeypatch, capsys):
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreUsageCommand().execute(
            _usage_args(unused=True, json=True, registry=str(project_dir / "prepare_registry.json"))
        )
        payload = json.loads(capsys.readouterr().out)
        assert rc == 1
        assert payload["error"].startswith("No store configured")

    def _seed(self, tmp_path):
        root = tmp_path / "store"
        paths = init_store(root)
        proj_a_dir = tmp_path / "proja"
        proj_a_dir.mkdir()
        proj_b_dir = tmp_path / "projb"
        proj_b_dir.mkdir()
        with catalog_write(paths) as cat:
            cat.upsert_project("proja", "Wolbachia", str(proj_a_dir), str(proj_a_dir / "prepare_registry.json"))
            cat.upsert_project("projb", "Rickettsia", str(proj_b_dir), str(proj_b_dir / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1"))
            cat.upsert_dataset(_sidecar("SRR2"))
            cat.record_usage("SRR1", "proja", "wMel", "downloaded")
            cat.record_usage("SRR1", "projb", "wMel", "analysed")
            cat.record_usage("SRR2", "proja", "wRi", "downloaded")
        return root

    def test_selector_by_accession(self, tmp_path, capsys):
        root = self._seed(tmp_path)
        rc = StoreUsageCommand().execute(_usage_args(data_root=str(root), accession="SRR1", json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["selector"] == "accession"
        project_ids = {row["project_id"] for row in report["rows"]}
        assert project_ids == {"proja", "projb"}
        for row in report["rows"]:
            assert set(row) == {"project_name", "project_id", "genome_id", "stage", "first_used", "last_used"}

    def test_selector_by_project_name(self, tmp_path, capsys):
        root = self._seed(tmp_path)
        rc = StoreUsageCommand().execute(_usage_args(data_root=str(root), project="Wolbachia", json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["selector"] == "project"
        accessions = {row["accession"] for row in report["rows"]}
        assert accessions == {"SRR1", "SRR2"}

    def test_selector_by_project_id(self, tmp_path, capsys):
        root = self._seed(tmp_path)
        rc = StoreUsageCommand().execute(_usage_args(data_root=str(root), project="projb", json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        accessions = {row["accession"] for row in report["rows"]}
        assert accessions == {"SRR1"}

    def test_selector_by_ambiguous_project_name_exits_1(self, tmp_path, caplog):
        root = self._seed(tmp_path)
        paths = store_paths(root)
        with catalog_write(paths) as cat:
            proj_c_dir = tmp_path / "projc"
            proj_c_dir.mkdir()
            cat.upsert_project("projc", "Wolbachia", str(proj_c_dir), str(proj_c_dir / "prepare_registry.json"))

        with caplog.at_level("ERROR"):
            rc = StoreUsageCommand().execute(_usage_args(data_root=str(root), project="Wolbachia"))

        assert rc == 1
        assert any("proja" in message and "projc" in message for message in caplog.messages)

    def test_selector_by_organism(self, tmp_path, capsys):
        root = self._seed(tmp_path)
        rc = StoreUsageCommand().execute(_usage_args(data_root=str(root), organism="wMel", json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["selector"] == "organism"
        accessions = {row["accession"] for row in report["rows"]}
        assert accessions == {"SRR1"}

    def test_selector_unused(self, tmp_path, capsys):
        root = self._seed(tmp_path)
        paths = store_paths(root)
        with catalog_write(paths) as cat:
            cat.upsert_dataset(_sidecar("SRR3"))

        rc = StoreUsageCommand().execute(_usage_args(data_root=str(root), unused=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["selector"] == "unused"
        assert report["rows"] == [{"accession": "SRR3", "state": "complete", "bytes": 100}]

    def test_selector_bytes_by_organism(self, tmp_path, capsys):
        root = self._seed(tmp_path)
        rc = StoreUsageCommand().execute(_usage_args(data_root=str(root), bytes_by_organism=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["selector"] == "bytes-by-organism"
        by_genome = {row["genome_id"]: row for row in report["rows"]}
        assert by_genome["wMel"]["datasets"] == 1
        assert by_genome["wMel"]["bytes"] == 100
        assert by_genome["wRi"]["datasets"] == 1
        assert by_genome["wRi"]["bytes"] == 100

    def test_text_output_lists_projects(self, tmp_path, capsys):
        root = self._seed(tmp_path)
        rc = StoreUsageCommand().execute(_usage_args(data_root=str(root), accession="SRR1"))
        out = capsys.readouterr().out

        assert rc == 0
        assert "proja" in out and "projb" in out


class TestStoreInitRebinding:
    """store_init is run again: keep what the project already links, and say what changed."""

    def test_rebinding_to_another_root_preserves_linked_and_warns(self, tmp_path, monkeypatch, caplog):
        import logging

        first = tmp_path / "store-a"
        second = tmp_path / "store-b"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"

        assert StoreInitCommand().execute(_init_args(first, project_dir)) == 0
        registry = load_registry(registry_path)
        registry.store["linked"] = ["SRR1", "SRR2"]
        from prepare.data.registry import save_registry

        save_registry(registry, registry_path)

        with caplog.at_level(logging.WARNING):
            assert StoreInitCommand().execute(_init_args(second, project_dir)) == 0

        registry = load_registry(registry_path)
        assert registry.store["root"] == str(second.resolve())
        # The list of what this project links is not something store_init knows how to rebuild.
        assert registry.store["linked"] == ["SRR1", "SRR2"]
        assert any("was bound to" in record.message for record in caplog.records)

    def test_rerun_keeps_keys_it_does_not_know(self, tmp_path, monkeypatch):
        from prepare.data.registry import save_registry

        root = tmp_path / "store"
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"
        assert StoreInitCommand().execute(_init_args(root, project_dir)) == 0
        registry = load_registry(registry_path)
        registry.store["comment"] = "written by a newer version"
        registry.project["notes"] = "also kept"
        save_registry(registry, registry_path)

        assert StoreInitCommand().execute(_init_args(root, project_dir)) == 0

        registry = load_registry(registry_path)
        assert registry.store["comment"] == "written by a newer version"
        assert registry.project["notes"] == "also kept"
        assert set(registry.store) == {"root", "mode", "linked", "comment"}

    def test_refuses_a_non_empty_folder_that_is_not_a_store(self, tmp_path, monkeypatch):
        home_like = tmp_path / "documents"
        home_like.mkdir()
        (home_like / "thesis.txt").write_text("chapter one")
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreInitCommand().execute(_init_args(home_like, project_dir))

        assert rc == 1
        assert not (home_like / "sra").exists()
        assert not (home_like / "prepare_store.json").exists()

    def test_accepts_an_empty_folder(self, tmp_path, monkeypatch):
        empty = tmp_path / "new-store"
        empty.mkdir()
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        assert StoreInitCommand().execute(_init_args(empty, project_dir)) == 0
        assert read_marker(empty) is not None


class TestStoreLinkGuardsIncompleteDatasets:
    def _store_with(self, tmp_path, state):
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        acc_dir.mkdir(parents=True, exist_ok=True)
        with gzip.open(acc_dir / "SRR1.fastq.gz", "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n")
        if state is not None:
            write_sidecar(sidecar_path(paths, "SRR1"), Sidecar(accession="SRR1", state=state))
        return root, paths

    @pytest.mark.parametrize("state", ["partial", "failed", None])
    def test_an_incomplete_dataset_needs_accept_partial(self, tmp_path, monkeypatch, capsys, state):
        root, paths = self._store_with(tmp_path, state)
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreLinkCommand().execute(_link_args(["SRR1"], data_root=str(root)))

        assert rc == 1
        assert not (project_dir / "fastq" / "SRR1").exists()

        rc = StoreLinkCommand().execute(_link_args(["SRR1"], data_root=str(root), accept_partial=True))

        assert rc == 0
        assert (project_dir / "fastq" / "SRR1").is_symlink()

    def test_a_complete_dataset_links_without_the_flag(self, tmp_path, monkeypatch):
        root, paths = self._store_with(tmp_path, "complete")
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        assert StoreLinkCommand().execute(_link_args(["SRR1"], data_root=str(root))) == 0
        assert (project_dir / "fastq" / "SRR1").is_symlink()


class TestLinkersAlwaysHaveAProjectIdentity:
    """A project that never ran store_init still records who it is, so gc can see its links."""

    def test_store_link_mints_a_project_identity(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        acc_dir.mkdir(parents=True, exist_ok=True)
        with gzip.open(acc_dir / "SRR1.fastq.gz", "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n")
        write_sidecar(sidecar_path(paths, "SRR1"), Sidecar(accession="SRR1", state="complete"))

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"

        # No store_init: the store was found through --data-root alone.
        assert StoreLinkCommand().execute(_link_args(["SRR1"], data_root=str(root))) == 0

        registry = load_registry(registry_path)
        assert registry.project["id"]
        assert registry.project["name"] == project_dir.name
        with Catalog(paths) as catalog:
            rows = catalog.conn.execute("SELECT project_id FROM projects").fetchall()
        assert [row["project_id"] for row in rows] == [registry.project["id"]]

    def test_store_adopt_mints_a_project_identity(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        paths = init_store(root)
        project_dir = tmp_path / "project"
        (project_dir / "fastq" / "SRR1").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1" / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        monkeypatch.chdir(project_dir)

        rc = StoreAdoptCommand().execute(_adopt_args(data_root=str(root)))

        assert rc == 0
        registry = load_registry(project_dir / "prepare_registry.json")
        assert registry.project["id"]
        with Catalog(paths) as catalog:
            assert catalog.conn.execute("SELECT COUNT(*) AS n FROM usage").fetchone()["n"] == 1

    def test_store_link_keeps_recorded_exports_and_other_project_keys(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        acc_dir.mkdir(parents=True, exist_ok=True)
        with gzip.open(acc_dir / "SRR1.fastq.gz", "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n")
        write_sidecar(sidecar_path(paths, "SRR1"), Sidecar(accession="SRR1", state="complete"))

        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)
        registry_path = project_dir / "prepare_registry.json"
        registry = load_registry(registry_path)
        registry.project = {
            "exports": {"results_table": {"output": "results.tsv", "summary": {"rows": 1}}},
            "note": "kept",
        }
        save_registry(registry)

        assert StoreLinkCommand().execute(_link_args(["SRR1"], data_root=str(root))) == 0

        registry = load_registry(registry_path)
        assert registry.project["id"]
        assert registry.project["exports"]["results_table"]["output"] == "results.tsv"
        assert registry.project["note"] == "kept"


class TestEnsureProjectIdentity:
    def test_minting_an_identity_keeps_unknown_project_keys(self, tmp_path, monkeypatch):
        from prepare.store.usage import ensure_project_identity

        monkeypatch.chdir(tmp_path)
        registry = load_registry(tmp_path / "prepare_registry.json")
        registry.project = {"exports": {"results_table": {"output": "results.tsv"}}, "note": "kept"}

        project = ensure_project_identity(registry)

        assert project["id"]
        assert registry.project["exports"] == {"results_table": {"output": "results.tsv"}}
        assert registry.project["note"] == "kept"
        assert registry.project["name"] == tmp_path.name


class TestStoreCommandsWithoutAMarker:
    """A store root that is not a store: the error names the rule that produced the root."""

    def test_store_status_reports_the_missing_marker(self, tmp_path, monkeypatch, caplog):
        import logging

        not_a_store = tmp_path / "not-a-store"
        not_a_store.mkdir()
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        with caplog.at_level(logging.ERROR):
            rc = StoreStatusCommand().execute(_status_args(data_root=str(not_a_store)))

        assert rc == 1
        assert any("no store marker" in record.message for record in caplog.records)


class TestStoreAdoptGitignoreGuard:
    def test_adopt_adds_fastq_to_gitignore(self, tmp_path, monkeypatch):
        root = tmp_path / "store"
        init_store(root)
        project_dir = tmp_path / "project"
        (project_dir / ".git").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1" / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        monkeypatch.chdir(project_dir)

        with patch("subprocess.run", return_value=subprocess.CompletedProcess([], 0, "", "")):
            assert StoreAdoptCommand().execute(_adopt_args(data_root=str(root))) == 0

        assert "fastq/" in (project_dir / ".gitignore").read_text()


class TestStoreReindexNeverLosesHistory:
    """A reindex rebuilds from what it reads, so it must read everything or stop."""

    def _store_with_dataset(self, tmp_path):
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        acc_dir.mkdir(parents=True, exist_ok=True)
        with gzip.open(acc_dir / "SRR1.fastq.gz", "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n")
        write_sidecar(sidecar_path(paths, "SRR1"), Sidecar(accession="SRR1", state="complete"))
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "P1", str(tmp_path / "p1"), str(tmp_path / "p1" / "prepare_registry.json"))
            cat.upsert_dataset(read_sidecar(sidecar_path(paths, "SRR1")))
            cat.record_usage("SRR1", "p1", "wMel", "downloaded")
        return root, paths

    def test_an_unreadable_sidecar_aborts_the_reindex(self, tmp_path, caplog):
        import logging

        root, paths = self._store_with_dataset(tmp_path)
        sidecar_path(paths, "SRR1").write_text("{ this is not json")

        with caplog.at_level(logging.ERROR):
            rc = StoreReindexCommand().execute(_reindex_args(data_root=str(root)))

        assert rc == 1
        assert any("SRR1" in record.message for record in caplog.records)
        # The usage history, which lives only in the catalogue, is untouched.
        with Catalog(paths) as catalog:
            assert catalog.conn.execute("SELECT COUNT(*) AS n FROM usage").fetchone()["n"] == 1

    def test_a_sidecar_holding_a_json_list_aborts_the_reindex_as_unreadable(self, tmp_path, caplog):
        import logging

        root, paths = self._store_with_dataset(tmp_path)
        sidecar_path(paths, "SRR1").write_text('[{"accession": "SRR1"}]')

        with caplog.at_level(logging.ERROR):
            rc = StoreReindexCommand().execute(_reindex_args(data_root=str(root)))

        assert rc == 1
        assert any("SRR1" in record.message for record in caplog.records)
        with Catalog(paths) as catalog:
            assert catalog.conn.execute("SELECT COUNT(*) AS n FROM usage").fetchone()["n"] == 1

    def test_reindex_after_catalog_loss_restores_projects_and_usage(self, tmp_path, monkeypatch):
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        write_sidecar(sidecar_path(paths, "SRR1"), _sidecar("SRR1"))
        with catalog_write(paths) as c:
            c.upsert_dataset(_sidecar("SRR1"))
            c.upsert_project("pid1", "proj", str(tmp_path), "r.json")
            c.record_usage("SRR1", "pid1", "", "linked", "")
        (paths.root / "catalog.sqlite").unlink()
        assert StoreReindexCommand().execute(_reindex_args(data_root=str(paths.root))) == 0
        with catalog_write(paths) as c:
            assert c.conn.execute("SELECT COUNT(*) FROM projects").fetchone()[0] == 1
            assert c.conn.execute("SELECT COUNT(*) FROM usage WHERE accession='SRR1'").fetchone()[0] == 1

    def test_reindex_survives_a_corrupt_project_journal_line(self, tmp_path, caplog):
        """A damaged project line must not sink the whole reindex: the dataset rebuild (already
        done in the same catalog_write transaction) must still be kept, and the command must
        still return 0, with the usage row that referenced the lost project simply skipped."""
        import logging

        root, paths = self._store_with_dataset(tmp_path)
        (paths.journal / "projects.jsonl").write_text("{ this is not json\n")
        (paths.root / "catalog.sqlite").unlink()

        with caplog.at_level(logging.WARNING):
            rc = StoreReindexCommand().execute(_reindex_args(data_root=str(root)))

        assert rc == 0
        with Catalog(paths) as catalog:
            assert catalog.conn.execute("SELECT COUNT(*) FROM datasets").fetchone()[0] == 1
            assert catalog.conn.execute("SELECT COUNT(*) FROM projects").fetchone()[0] == 0
            assert catalog.conn.execute("SELECT COUNT(*) FROM usage").fetchone()[0] == 0
        assert any("Skipping unreadable journal line" in r.message for r in caplog.records)
        assert any("Skipped 1 usage record" in r.message for r in caplog.records)

    def test_reindex_recovers_a_pre_journal_project_once_backfilled(self, tmp_path):
        """A project and usage row written before the journal existed (journaling disabled)
        are copied into the journal the next time the store is opened for writing, while
        catalog.sqlite still holds them; only then does losing catalog.sqlite stay recoverable."""
        root = tmp_path / "store"
        paths = init_store(root)
        # SRR1's folder and sidecar must actually be on disk, or the later reindex this test
        # exercises (which now replays before rebuilding, see test_reindex_does_not_resurrect_a_
        # removed_dataset) would correctly treat this usage row's accession as gone rather than
        # resurrect a placeholder for it.
        write_sidecar(sidecar_path(paths, "SRR1"), _sidecar("SRR1"))
        with catalog_write(paths) as c:
            c.journal_enabled = False
            c.upsert_project("pid1", "proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            c.upsert_dataset(_sidecar("SRR1"))
            c.record_usage("SRR1", "pid1", "", "linked", "")
        assert not (paths.journal / "projects.jsonl").exists()

        with catalog_write(paths):
            pass  # any ordinary write session backfills the journal from catalog.sqlite
        assert (paths.journal / "projects.jsonl").is_file()

        (paths.root / "catalog.sqlite").unlink()
        assert StoreReindexCommand().execute(_reindex_args(data_root=str(root))) == 0
        with catalog_write(paths) as c:
            assert c.conn.execute("SELECT COUNT(*) FROM projects").fetchone()[0] == 1
            assert c.conn.execute("SELECT COUNT(*) FROM usage WHERE accession='SRR1'").fetchone()[0] == 1

    def test_reindex_does_not_resurrect_a_removed_dataset(self, tmp_path, monkeypatch):
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        write_sidecar(sidecar_path(paths, "SRR1"), _sidecar("SRR1"))
        with catalog_write(paths) as c:
            c.upsert_dataset(_sidecar("SRR1"))
            c.upsert_project("pid1", "proj", str(tmp_path), "r.json")
            c.record_usage("SRR1", "pid1", "", "linked", "")
            c.record_usage("SRR9", "pid1", "", "linked", "")  # SRR9 was removed by gc later
        assert StoreReindexCommand().execute(_reindex_args(data_root=str(paths.root))) == 0
        with catalog_write(paths) as c:
            rows = {r["accession"] for r in c.conn.execute("SELECT accession FROM datasets").fetchall()}
        assert rows == {"SRR1"}

    def test_replay_restores_hostname_and_first_used(self, tmp_path, monkeypatch):
        paths = init_store(tmp_path / "store")
        monkeypatch.chdir(tmp_path)
        write_sidecar(sidecar_path(paths, "SRR1"), _sidecar("SRR1"))
        with catalog_write(paths) as c:
            c.upsert_dataset(_sidecar("SRR1"))
            c.upsert_project("pid1", "proj", str(tmp_path), "r.json", hostname="other-host")
            c.record_usage("SRR1", "pid1", "", "linked", "", at="2026-01-01T00:00:00+00:00")
        (paths.root / "catalog.sqlite").unlink()
        assert StoreReindexCommand().execute(_reindex_args(data_root=str(paths.root))) == 0
        with catalog_write(paths) as c:
            row = c.conn.execute("SELECT hostname FROM projects WHERE project_id='pid1'").fetchone()
            used = c.conn.execute("SELECT first_used FROM usage WHERE accession='SRR1'").fetchone()
        assert row["hostname"] == "other-host"
        assert used["first_used"].startswith("2026-01-01")


class TestCorruptStoreMarker:
    def test_a_corrupt_marker_reads_as_missing_with_a_warning(self, tmp_path, caplog):
        import logging

        from prepare.store.layout import store_paths as _store_paths

        root = tmp_path / "store"
        init_store(root)
        _store_paths(root).marker.write_text("{ truncated")

        with caplog.at_level(logging.WARNING):
            assert read_marker(root) is None

        assert any("marker" in record.message for record in caplog.records)


def test_reindex_and_verify_ignore_hidden_entries(tmp_path, monkeypatch):
    paths = init_store(tmp_path / "store")
    monkeypatch.chdir(tmp_path)
    acc_dir = sra_dir(paths, "SRR1")
    _write_fastq_gz(acc_dir / "SRR1.fastq.gz")
    write_sidecar(sidecar_path(paths, "SRR1"), _sidecar_matching_disk(acc_dir, "SRR1"))
    (paths.sra / "._SRR1").write_bytes(b"\x00\x05")
    (paths.sra / ".sra-cache").mkdir()
    (paths.sra / ".hidden_dir").mkdir()
    assert StoreReindexCommand().execute(_reindex_args(data_root=str(paths.root))) == 0
    rc = StoreVerifyCommand().execute(_verify_args(data_root=str(paths.root)))
    assert rc == 0


def test_every_store_command_failure_starts_with_the_command_name():
    """Store command errors read ``store_<name>: ...``; a literal context would break that."""
    import re

    from pathlib import Path

    import prepare.cli.commands.store as store_pkg

    folder = Path(store_pkg.__file__).parent
    contexts = []
    for module in sorted(folder.glob("*.py")):
        contexts += re.findall(r"self\.fail\(e, ([^)]*)\)", module.read_text())
    assert contexts and set(contexts) == {"self.name"}


class TestStoreReindexMissingSidecar:
    """A store folder with no sidecar at all is described from disk rather than blocking the reindex."""

    def _sidecar_less_store(self, tmp_path):
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = sra_dir(paths, "SRR1")
        acc_dir.mkdir(parents=True, exist_ok=True)
        with gzip.open(acc_dir / "SRR1.fastq.gz", "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n")
        return root, paths

    def test_a_missing_sidecar_is_rebuilt_and_catalogued(self, tmp_path, caplog):
        root, paths = self._sidecar_less_store(tmp_path)

        with caplog.at_level(logging.WARNING):
            rc = StoreReindexCommand().execute(_reindex_args(data_root=str(root)))

        assert rc == 0
        sidecar = read_sidecar(sidecar_path(paths, "SRR1"))
        assert sidecar is not None and sidecar.state == "failed"
        assert any("SRR1" in record.message and record.levelno == logging.WARNING for record in caplog.records)
        with Catalog(paths) as catalog:
            row = catalog.get_dataset("SRR1")
        assert row is not None and row["state"] == "failed"

    def test_a_missing_sidecar_whose_lock_is_held_is_skipped_with_a_warning(self, tmp_path, caplog):
        root, paths = self._sidecar_less_store(tmp_path)
        # A fresh lock file naming a run on another host: alive as far as this host can tell.
        holder = {"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00", "token": "feedbeef"}
        paths.locks.mkdir(parents=True, exist_ok=True)
        lock_path(paths, "SRR1").write_text(json.dumps(holder))

        with caplog.at_level(logging.WARNING):
            rc = StoreReindexCommand().execute(_reindex_args(data_root=str(root)))

        assert rc == 0
        assert not sidecar_path(paths, "SRR1").exists()
        assert any(
            "SRR1" in record.message and record.levelno == logging.WARNING and "lock" in record.message
            for record in caplog.records
        )

    def test_an_unreadable_sidecar_still_refuses_and_rebuilds_nothing(self, tmp_path, caplog):
        root, paths = self._sidecar_less_store(tmp_path)
        bad_dir = sra_dir(paths, "SRR2")
        bad_dir.mkdir(parents=True)
        sidecar_path(paths, "SRR2").write_text("{ not json")

        with caplog.at_level(logging.ERROR):
            rc = StoreReindexCommand().execute(_reindex_args(data_root=str(root)))

        assert rc == 1
        assert any("SRR2" in record.message for record in caplog.records)
        assert not sidecar_path(paths, "SRR1").exists()

    def test_a_folder_not_named_like_an_accession_is_skipped(self, tmp_path, caplog):
        root, paths = self._sidecar_less_store(tmp_path)
        stray = paths.sra / "notes"
        stray.mkdir()
        (stray / "readme.txt").write_text("kept by hand\n")

        with caplog.at_level(logging.WARNING):
            rc = StoreReindexCommand().execute(_reindex_args(data_root=str(root)))

        assert rc == 0
        assert not (stray / "notes.json").exists()
        assert sidecar_path(paths, "SRR1").is_file()
        assert any("notes" in record.message and "skipped" in record.message for record in caplog.records)
        with Catalog(paths) as catalog:
            assert catalog.get_dataset("notes") is None
            assert catalog.get_dataset("SRR1") is not None

    def test_a_rebuilt_sidecar_is_refused_by_store_link_until_verified(self, tmp_path, monkeypatch):
        """The repair path the rebuilt sidecar names: reindex, a refused link, verify --rescan --fix-state."""
        root, paths = self._sidecar_less_store(tmp_path)
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        assert StoreReindexCommand().execute(_reindex_args(data_root=str(root))) == 0
        assert read_sidecar(sidecar_path(paths, "SRR1")).state == "failed"

        assert StoreLinkCommand().execute(_link_args(["SRR1"], data_root=str(root))) == 1
        assert not (project_dir / "fastq" / "SRR1").exists()

        rc = StoreVerifyCommand().execute(_verify_args(["SRR1"], data_root=str(root), rescan=True, fix_state=True))
        assert rc == 0
        assert read_sidecar(sidecar_path(paths, "SRR1")).state == "complete"

        assert StoreLinkCommand().execute(_link_args(["SRR1"], data_root=str(root))) == 0
        assert (project_dir / "fastq" / "SRR1").is_symlink()
