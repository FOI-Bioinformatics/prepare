"""Tests for the `store_gc` CLI command.

Every test runs under tmp_path and monkeypatches HOME/XDG_CONFIG_HOME/PREPARE_DATA so
nothing here reads or writes the real user config, matching the isolation pattern used in
tests/test_cli_store.py and tests/test_store_resolve.py.
"""

import argparse
import json
import logging
import os
import threading

import pytest

from prepare.cli.commands.store import (
    StoreAdoptCommand,
    StoreGcCommand,
    StoreInitCommand,
    StoreLinkCommand,
    StoreReindexCommand,
)
from prepare.core.constants import STORE_ENV
from prepare.core.exceptions import DataAccessError
from prepare.data.registry import load_registry
from prepare.store.catalog import REBUILT_WITHOUT_PROJECTS, Catalog, catalog_write
from prepare.store.layout import init_store, sidecar_path, sra_dir, store_paths
from prepare.store.locks import dataset_lock, touch_dataset_use
from prepare.store.sidecar import Sidecar, write_sidecar


@pytest.fixture(autouse=True)
def isolated_env(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    monkeypatch.delenv(STORE_ENV, raising=False)
    yield


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


def _sidecar(accession="SRR1", state="complete", downloaded="2026-09-06T00:00:00+00:00"):
    return Sidecar(
        accession=accession,
        state=state,
        layout="SINGLE",
        downloaded=downloaded,
        tool="fasterq-dump",
        tool_version="3.0.0",
        compression="gzip",
        files=[{"name": f"{accession}.fastq.gz", "bytes": 100, "md5": "a" * 32, "reads": 5}],
        reads_per_mate=5,
        bases_total=500,
        ncbi={"spots": 5, "bases": 500, "size": 1000, "layout": "SINGLE", "files": []},
        completeness={"method": "spots", "ratio": 1.0, "verdict": "complete"},
    )


def _write_dataset_dir(paths, accession):
    acc_dir = sra_dir(paths, accession)
    acc_dir.mkdir(parents=True, exist_ok=True)
    (acc_dir / f"{accession}.fastq.gz").write_bytes(b"x" * 100)
    return acc_dir


def _init_args(root, project_dir, **overrides):
    base = dict(
        data_root=str(root),
        project_name=None,
        set_default=False,
        registry=str(project_dir / "prepare_registry.json"),
    )
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


class TestStoreGcCommand:
    def test_command_properties(self):
        cmd = StoreGcCommand()
        assert cmd.name == "store_gc"
        assert cmd.group == "Store"

    def test_no_store_configured_returns_1(self, tmp_path, monkeypatch, capsys):
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreGcCommand().execute(_gc_args(registry=str(project_dir / "prepare_registry.json")))
        assert rc == 1

    def test_no_store_configured_with_json_prints_json_error(self, tmp_path, monkeypatch, capsys):
        """No store, --json requested: the hint must be valid JSON on stdout, same shape as
        the refusal in test_json_refusal_prints_an_error_object, not the plain-text hint."""
        project_dir = tmp_path / "project"
        project_dir.mkdir()
        monkeypatch.chdir(project_dir)

        rc = StoreGcCommand().execute(_gc_args(json=True, registry=str(project_dir / "prepare_registry.json")))
        payload = json.loads(capsys.readouterr().out)

        assert rc == 1
        assert payload["error"].startswith("No store configured")

    def test_json_refusal_prints_an_error_object(self, tmp_path, capsys):
        """A refusal before --json can even build a report must still be visible on stdout as
        JSON, not only logged, so a script driving store_gc --json can parse it."""
        root = tmp_path / "store"
        paths = init_store(root)
        _write_dataset_dir(paths, "SRR1")
        with catalog_write(paths) as cat:
            cat.upsert_dataset(_sidecar("SRR1"))  # no project recorded at all

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), json=True))
        out = capsys.readouterr().out

        assert rc == 1
        payload = json.loads(out)
        assert "records no project at all" in payload["error"]

    def test_dry_run_refusal_under_no_project_returns_1_and_touches_nothing(self, tmp_path, capsys):
        """--dry-run must not bypass the no-project refusal: the refusal is checked before
        candidates are even built, so it applies whether or not --yes would have followed."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = _write_dataset_dir(paths, "SRR1")
        with catalog_write(paths) as cat:
            cat.upsert_dataset(_sidecar("SRR1"))  # no project recorded at all

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), dry_run=True))
        out = capsys.readouterr().out

        assert rc == 1
        assert out == ""  # no --json: the refusal is only logged, nothing on stdout
        assert acc_dir.is_dir()
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1") is not None

    def test_dry_run_and_json_refusal_prints_an_error_object(self, tmp_path, capsys):
        """--dry-run combined with --json under the same no-project refusal must still print
        the JSON error object (not an empty report) and remove nothing."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = _write_dataset_dir(paths, "SRR1")
        with catalog_write(paths) as cat:
            cat.upsert_dataset(_sidecar("SRR1"))  # no project recorded at all

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), dry_run=True, json=True))
        out = capsys.readouterr().out

        assert rc == 1
        payload = json.loads(out)
        assert "records no project at all" in payload["error"]
        assert acc_dir.is_dir()
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1") is not None

    def test_unused_dataset_listed_with_bytes(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        _write_dataset_dir(paths, "SRR1")
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1"))

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        accessions = {row["accession"]: row for row in report["datasets"]}
        assert accessions["SRR1"]["bytes"] == 100
        assert accessions["SRR1"]["reason"] == "unused"

    def test_dataset_used_only_by_stale_project_needs_include_stale(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        _write_dataset_dir(paths, "SRR1")
        stale_dir = tmp_path / "stale_proj"
        stale_dir.mkdir()
        with catalog_write(paths) as cat:
            cat.upsert_project("stale1", "OldProject", str(stale_dir), str(stale_dir / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1"))
            cat.record_usage("SRR1", "stale1", "wMel", "downloaded")

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), json=True))
        report = json.loads(capsys.readouterr().out)

        # A project is "stale" whenever its registry is not a file on this machine, which is
        # true of every project on another workstation of a shared store: by default such a
        # dataset is kept and only reported.
        assert rc == 0
        assert report["datasets"] == []
        kept = {row["accession"]: row for row in report["kept_stale"]}
        assert "OldProject" in kept["SRR1"]["reason"]

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), include_stale=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        accessions = {row["accession"]: row for row in report["datasets"]}
        assert "OldProject" in accessions["SRR1"]["reason"]
        assert "stale" in accessions["SRR1"]["reason"]

    def test_live_linked_dataset_never_listed(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        _write_dataset_dir(paths, "SRR1")
        live_dir = tmp_path / "live_proj"
        live_dir.mkdir()
        registry_file = live_dir / "prepare_registry.json"
        registry_file.write_text(json.dumps({"project": {"id": "live1"}}))
        with catalog_write(paths) as cat:
            cat.upsert_project("live1", "LiveProject", str(live_dir), str(registry_file))
            cat.upsert_dataset(_sidecar("SRR1"))
            cat.record_usage("SRR1", "live1", "wMel", "downloaded")

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["datasets"] == []

    def test_keep_partial_spares_partial_datasets(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        _write_dataset_dir(paths, "SRR1")
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1", state="partial"))

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), keep_partial=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["datasets"] == []

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), keep_partial=False, json=True))
        report = json.loads(capsys.readouterr().out)
        assert len(report["datasets"]) == 1

    def test_older_than_filters_recent_datasets(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        _write_dataset_dir(paths, "SRR1")
        from datetime import datetime, timezone

        recent = datetime.now(timezone.utc).isoformat()
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1", downloaded=recent))

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), older_than=30, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["datasets"] == []

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), json=True))
        report = json.loads(capsys.readouterr().out)
        assert len(report["datasets"]) == 1

    def test_leftovers_listed_with_sizes(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        with catalog_write(paths):
            pass
        (paths.tmp / "SRR9_temp").mkdir(parents=True)
        (paths.tmp / "SRR9_temp" / "partial.fastq").write_bytes(b"x" * 50)
        (paths.tmp / "SRR8_adopt").mkdir(parents=True)
        (paths.tmp / "SRR8_adopt" / "SRR8.fastq").write_bytes(b"y" * 30)
        (paths.tmp / ".sra-cache" / "SRR7.sra").parent.mkdir(parents=True, exist_ok=True)
        (paths.tmp / ".sra-cache" / "SRR7.sra").write_bytes(b"z" * 20)
        (paths.sra / ".sra-cache").mkdir(parents=True, exist_ok=True)
        (paths.sra / ".sra-cache" / "leftover.sra").write_bytes(b"w" * 10)

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        total_bytes = sum(entry["bytes"] for entry in report["leftovers"])
        assert total_bytes == 50 + 30 + 20 + 10
        assert len(report["leftovers"]) == 4
        assert all(entry["reason"] == "leftover" for entry in report["leftovers"])

    def test_dry_run_deletes_nothing(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = _write_dataset_dir(paths, "SRR1")
        (paths.tmp / "SRR9_temp").mkdir(parents=True)
        (paths.tmp / "SRR9_temp" / "f").write_bytes(b"x" * 10)
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1"))

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root)))
        assert rc == 0

        assert acc_dir.is_dir()
        assert (paths.tmp / "SRR9_temp").is_dir()
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1") is not None

    def test_yes_removes_folders_catalogue_rows_and_leftovers(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = _write_dataset_dir(paths, "SRR1")
        (paths.tmp / "SRR9_temp").mkdir(parents=True)
        (paths.tmp / "SRR9_temp" / "f").write_bytes(b"x" * 10)
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1"))

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), dry_run=False, yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert not acc_dir.exists()
        assert not (paths.tmp / "SRR9_temp").exists()
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1") is None
        assert "SRR1" in report["removed_datasets"]
        assert any("SRR9_temp" in entry for entry in report["removed_leftovers"])

    def test_dry_run_and_yes_together_refuses_and_removes_nothing(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = _write_dataset_dir(paths, "SRR1")
        with catalog_write(paths) as cat:
            cat.upsert_dataset(_sidecar("SRR1"))

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), dry_run=True, yes=True))

        assert rc == 1
        assert acc_dir.is_dir()
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1") is not None

    def test_live_project_symlink_without_usage_row_keeps_dataset(self, tmp_path, capsys):
        """A dataset a live project still symlinks to must never be removed, even when it has
        no usage row at all (usage tracking predating that project, or a failed hook)."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = _write_dataset_dir(paths, "SRR1")

        project_dir = tmp_path / "live_proj"
        (project_dir / "fastq").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1").symlink_to(acc_dir, target_is_directory=True)
        registry_file = project_dir / "prepare_registry.json"
        registry_file.write_text(json.dumps({"project": {"id": "live1"}}))

        with catalog_write(paths) as cat:
            cat.upsert_project("live1", "LiveProject", str(project_dir), str(registry_file))
            cat.upsert_dataset(_sidecar("SRR1"))
            # Deliberately no usage row recorded for SRR1.

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["datasets"] == []
        still_linked = {row["accession"]: row for row in report["still_linked"]}
        assert "LiveProject" in still_linked["SRR1"]["reason"]

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert acc_dir.is_dir()
        assert "SRR1" not in report["removed_datasets"]
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1") is not None

    def test_stale_projects_named_in_report(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        stale_dir = tmp_path / "stale_proj"
        stale_dir.mkdir()
        missing_registry = stale_dir / "prepare_registry.json"
        with catalog_write(paths) as cat:
            cat.upsert_project("stale1", "OldProject", str(stale_dir), str(missing_registry))

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        row = report["stale_projects"][0]
        assert len(report["stale_projects"]) == 1
        assert row["project_id"] == "stale1"
        assert row["name"] == "OldProject"
        assert row["registry"] == str(missing_registry)
        assert row["reason"] == "registry missing"

    def test_gc_refuses_when_catalog_has_no_projects(self, tmp_path, caplog):
        import logging

        paths = init_store(tmp_path / "store")
        _write_dataset_dir(paths, "SRR1")
        write_sidecar(sidecar_path(paths, "SRR1"), _sidecar("SRR1"))
        with catalog_write(paths) as c:
            c.upsert_dataset(_sidecar("SRR1"))
        with caplog.at_level(logging.ERROR):
            rc = StoreGcCommand().execute(_gc_args(data_root=str(paths.root), yes=True))
        assert rc == 1
        assert (paths.sra / "SRR1" / "SRR1.fastq.gz").exists()
        assert any("no project" in record.message.lower() for record in caplog.records)

    def test_a_second_project_linking_the_same_accession_also_keeps_it_from_gc(self, tmp_path, monkeypatch, capsys):
        """A --copy adoption records a 'copied' usage row for the adopting project
        (test_copy_adoption_records_copied_usage_so_gc_keeps_it in tests/test_cli_store.py);
        a second, independent project that later links the same accession through
        store_link must record its own usage row too. Both are live (neither stale), so
        _classify_dataset's 'keep' branch applies and the accession is dropped from the
        report entirely, the same non-appearance the single-project case already pins via
        report["datasets"] == [] -- not surfaced under any bucket, since 'in_use' there is
        reserved for a download lock currently held, not for multi-project usage."""
        root = tmp_path / "store"
        init_store(root)

        project1_dir = tmp_path / "project1"
        entry = project1_dir / "fastq" / "SRR1"
        entry.mkdir(parents=True)
        (entry / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        registry1_path = project1_dir / "prepare_registry.json"
        monkeypatch.chdir(project1_dir)
        assert StoreInitCommand().execute(_init_args(root, project1_dir, registry=str(registry1_path))) == 0
        assert (
            StoreAdoptCommand().execute(_adopt_args(data_root=str(root), registry=str(registry1_path), move=False)) == 0
        )

        project2_dir = tmp_path / "project2"
        project2_dir.mkdir()
        registry2_path = project2_dir / "prepare_registry.json"
        monkeypatch.chdir(project2_dir)
        assert StoreInitCommand().execute(_init_args(root, project2_dir, registry=str(registry2_path))) == 0
        assert StoreLinkCommand().execute(_link_args(["SRR1"], data_root=str(root), registry=str(registry2_path))) == 0

        registry1 = load_registry(registry1_path)
        registry2 = load_registry(registry2_path)
        paths = store_paths(root)
        with Catalog(paths) as cat:
            rows = {
                row["project_id"]: row["stage"]
                for row in cat.conn.execute(
                    "SELECT project_id, stage FROM usage WHERE accession = ?", ("SRR1",)
                ).fetchall()
            }
        assert rows == {registry1.project["id"]: "copied", registry2.project["id"]: "linked"}

        capsys.readouterr()
        gc_rc = StoreGcCommand().execute(_gc_args(data_root=str(root), dry_run=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert gc_rc == 0
        assert report["datasets"] == []
        assert report["still_linked"] == []
        assert report["in_use"] == []
        assert report["kept_stale"] == []


class TestStoreGcRespectsLocksAndPlaceholders:
    """Never remove what another run is working on, or a row that stands for no files."""

    @staticmethod
    def _hold(paths, accession):
        import json as json_module

        from prepare.store.layout import lock_path

        lock = lock_path(paths, accession)
        lock.parent.mkdir(parents=True, exist_ok=True)
        lock.write_text(json_module.dumps({"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00"}))
        return lock

    def test_a_locked_dataset_is_kept_and_reported_as_in_use(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        _write_dataset_dir(paths, "SRR1")
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1"))
        self._hold(paths, "SRR1")

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["datasets"] == []
        assert [row["accession"] for row in report["in_use"]] == ["SRR1"]
        assert sra_dir(paths, "SRR1").is_dir()

    def test_leftovers_of_a_locked_accession_are_kept(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        with catalog_write(paths):
            pass
        building = paths.tmp / "SRR1_temp"
        building.mkdir(parents=True)
        (building / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        cached = paths.tmp / ".sra-cache" / "SRR1"
        cached.mkdir(parents=True)
        (cached / "SRR1.sra").write_bytes(b"x" * 10)
        self._hold(paths, "SRR1")

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["leftovers"] == []
        assert building.is_dir()
        assert cached.is_dir()

    def test_fasterq_dump_scratch_of_a_locked_accession_is_kept(self, tmp_path, capsys):
        """``<store>/tmp/<ACC>_fqtmp`` is fasterq-dump's live scratch folder while the
        accession's lock is held; gc must map it back to the accession and leave it alone."""
        root = tmp_path / "store"
        paths = init_store(root)
        with catalog_write(paths):
            pass
        scratch = paths.tmp / "SRR1_fqtmp"
        scratch.mkdir(parents=True)
        (scratch / "fasterq.tmp.part").write_bytes(b"x" * 10)
        self._hold(paths, "SRR1")

        assert StoreGcCommand._accession_of_leftover("SRR1_fqtmp") == "SRR1"
        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["leftovers"] == []
        assert report["removed_leftovers"] == []
        assert scratch.is_dir()

    def test_fasterq_dump_scratch_without_a_lock_is_a_leftover(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        with catalog_write(paths):
            pass
        scratch = paths.tmp / "SRR1_fqtmp"
        scratch.mkdir(parents=True)

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert any("SRR1_fqtmp" in entry for entry in report["removed_leftovers"])
        assert not scratch.exists()

    def test_bare_staged_download_folder_is_listed_and_removed_with_yes(self, tmp_path, capsys):
        """A killed download can leave a bare `<store>/tmp/<ACC>` staging folder behind
        (``store_handoff._store_fetch`` stages into ``store.tmp / accession`` before the
        dataset is published); once its lock is gone this is a plain leftover."""
        root = tmp_path / "store"
        paths = init_store(root)
        with catalog_write(paths):
            pass
        staged = paths.tmp / "SRR1"
        staged.mkdir(parents=True)
        (staged / "SRR1.sra").write_bytes(b"x" * 10)

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert [entry["reason"] for entry in report["leftovers"] if "SRR1" in entry["path"]] == [
            "leftover (staged download)"
        ]
        assert any("SRR1" in entry for entry in report["removed_leftovers"])
        assert not staged.exists()

    def test_bare_staged_download_folder_is_kept_while_its_lock_is_live(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        with catalog_write(paths):
            pass
        staged = paths.tmp / "SRR1"
        staged.mkdir(parents=True)
        (staged / "SRR1.sra").write_bytes(b"x" * 10)
        self._hold(paths, "SRR1")

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["leftovers"] == []
        assert report["removed_leftovers"] == []
        assert staged.is_dir()

    def test_other_tmp_folder_names_are_untouched(self, tmp_path, capsys):
        """A folder under `tmp/` that is neither a known transient suffix nor a bare
        accession name (e.g. not matching SRA_ACCESSION_PATTERN) is left alone entirely,
        not even listed."""
        root = tmp_path / "store"
        paths = init_store(root)
        with catalog_write(paths):
            pass
        other = paths.tmp / "notanaccession"
        other.mkdir(parents=True)
        (other / "file.txt").write_text("keep me")

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["leftovers"] == []
        assert report["removed_leftovers"] == []
        assert other.is_dir()

    def test_a_placeholder_row_is_never_a_candidate(self, tmp_path, capsys):
        """A usage row for an accession that is not catalogued yet inserts a state="unknown"
        placeholder; it stands for no files, so gc must not offer to remove it."""
        root = tmp_path / "store"
        paths = init_store(root)
        proj = tmp_path / "proj"
        proj.mkdir()
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(proj), str(proj / "prepare_registry.json"))
            cat.record_usage("SRR-PLACEHOLDER", "p1", "wMel", "linked")

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert [row["accession"] for row in report["datasets"]] == []

    def test_stale_projects_report_their_host_and_reason(self, tmp_path, capsys):
        root = tmp_path / "store"
        paths = init_store(root)
        gone = tmp_path / "gone"
        gone.mkdir()
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Gone", str(gone), str(gone / "prepare_registry.json"))

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        row = report["stale_projects"][0]
        assert row["reason"] == "registry missing"
        assert row["hostname"]


class TestStoreGcAfterARebuildWithoutProjects:
    """A reindex that restores no project leaves every dataset looking unused. The refusal must
    outlast the first project to register again, until the user says every project is back."""

    @staticmethod
    def _store_with_unjournaled_dataset(tmp_path):
        paths = init_store(tmp_path / "store")
        _write_dataset_dir(paths, "SRR1")
        write_sidecar(sidecar_path(paths, "SRR1"), _sidecar("SRR1"))
        return paths

    @staticmethod
    def _flag(paths):
        with Catalog(paths) as catalog:
            return catalog.get_meta(REBUILT_WITHOUT_PROJECTS)

    @staticmethod
    def _register_project(paths, tmp_path, project_id="p1"):
        proj = tmp_path / project_id
        proj.mkdir(exist_ok=True)
        with catalog_write(paths) as cat:
            cat.upsert_project(project_id, project_id, str(proj), str(proj / "prepare_registry.json"))

    def test_reindex_without_a_journal_sets_the_flag(self, tmp_path, caplog):
        import logging

        paths = self._store_with_unjournaled_dataset(tmp_path)

        with caplog.at_level(logging.WARNING):
            rc = StoreReindexCommand().execute(argparse.Namespace(data_root=str(paths.root), registry=None))

        assert rc == 0
        assert self._flag(paths)
        assert any("--accept-rebuilt" in record.message for record in caplog.records)

    def test_gc_refuses_while_the_flag_is_set_even_after_one_project_registers(self, tmp_path, caplog):
        import logging

        paths = self._store_with_unjournaled_dataset(tmp_path)
        StoreReindexCommand().execute(argparse.Namespace(data_root=str(paths.root), registry=None))
        self._register_project(paths, tmp_path)

        with caplog.at_level(logging.ERROR):
            rc = StoreGcCommand().execute(_gc_args(data_root=str(paths.root), yes=True))

        assert rc == 1
        assert sra_dir(paths, "SRR1").is_dir()
        message = " ".join(record.message for record in caplog.records)
        assert "--accept-rebuilt" in message
        assert "store_init" in message

    def test_gc_json_refuses_when_rebuilt_without_projects(self, tmp_path, capsys):
        """The refusal set by an unjournaled reindex must also reach a --json caller: a
        script parsing stdout needs the same error a human sees in the log, not empty
        output (same fixture as test_gc_refuses_while_the_flag_is_set_even_after_one_project
        _registers, with --json in place of --yes)."""
        paths = self._store_with_unjournaled_dataset(tmp_path)
        StoreReindexCommand().execute(argparse.Namespace(data_root=str(paths.root), registry=None))
        self._register_project(paths, tmp_path)
        capsys.readouterr()

        rc = StoreGcCommand().execute(_gc_args(data_root=str(paths.root), json=True))

        assert rc == 1
        assert sra_dir(paths, "SRR1").is_dir()
        report = json.loads(capsys.readouterr().out)
        assert "--accept-rebuilt" in report["error"]

    def test_accept_rebuilt_clears_the_flag_and_gc_proceeds(self, tmp_path, capsys):
        paths = self._store_with_unjournaled_dataset(tmp_path)
        StoreReindexCommand().execute(argparse.Namespace(data_root=str(paths.root), registry=None))
        self._register_project(paths, tmp_path)
        capsys.readouterr()

        rc = StoreGcCommand().execute(_gc_args(data_root=str(paths.root), accept_rebuilt=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert [row["accession"] for row in report["datasets"]] == ["SRR1"]
        assert self._flag(paths) is None
        # Nothing was removed: without --yes this is still a dry run.
        assert sra_dir(paths, "SRR1").is_dir()

    def test_accept_rebuilt_is_refused_until_a_project_registers(self, tmp_path, caplog):
        """--accept-rebuilt with no project recorded must not clear the flag: otherwise the first
        project to register would let a plain gc treat every other project's data as unused."""
        import logging

        paths = self._store_with_unjournaled_dataset(tmp_path)
        StoreReindexCommand().execute(argparse.Namespace(data_root=str(paths.root), registry=None))

        with caplog.at_level(logging.ERROR):
            rc = StoreGcCommand().execute(_gc_args(data_root=str(paths.root), accept_rebuilt=True, yes=True))

        assert rc == 1
        assert self._flag(paths)
        assert sra_dir(paths, "SRR1").is_dir()
        assert any("store_init" in record.message for record in caplog.records)

        self._register_project(paths, tmp_path)
        assert StoreGcCommand().execute(_gc_args(data_root=str(paths.root), yes=True)) == 1
        assert sra_dir(paths, "SRR1").is_dir()

        rc = StoreGcCommand().execute(_gc_args(data_root=str(paths.root), accept_rebuilt=True))

        assert rc == 0
        assert self._flag(paths) is None

    def test_reindex_that_restores_a_project_clears_the_flag(self, tmp_path):
        paths = self._store_with_unjournaled_dataset(tmp_path)
        StoreReindexCommand().execute(argparse.Namespace(data_root=str(paths.root), registry=None))
        self._register_project(paths, tmp_path)
        assert self._flag(paths)

        rc = StoreReindexCommand().execute(argparse.Namespace(data_root=str(paths.root), registry=None))

        assert rc == 0
        assert self._flag(paths) is None

    def test_empty_reindex_clears_a_stale_rebuilt_flag(self, tmp_path):
        """Once every dataset is removed from disk, a reindex restores no project and finds no
        datasets either; the flag must not survive as stale forever in that case."""
        import shutil

        paths = self._store_with_unjournaled_dataset(tmp_path)
        StoreReindexCommand().execute(argparse.Namespace(data_root=str(paths.root), registry=None))
        assert self._flag(paths)

        shutil.rmtree(sra_dir(paths, "SRR1"))

        rc = StoreReindexCommand().execute(argparse.Namespace(data_root=str(paths.root), registry=None))

        assert rc == 0
        assert self._flag(paths) is None


class TestStoreGcRemovalUnderTheLock:
    """--yes only ever removes a dataset under its own per-accession lock, re-checked right
    before removal: something that starts using an accession after the read-only
    classification pass, but before removal reaches it, must still keep the dataset."""

    def test_dataset_locked_between_classification_and_removal_is_kept_and_reported(
        self, tmp_path, capsys, monkeypatch
    ):
        """A real dataset_lock, taken by another thread only after classification has
        already finished (so the read-only pass sees it free), must still stop removal: the
        lock is re-taken, non-blocking, right before a candidate's folder is moved."""
        root = tmp_path / "store"
        paths = init_store(root)
        _write_dataset_dir(paths, "SRR1")
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1"))

        classified = threading.Event()
        locked = threading.Event()
        release = threading.Event()
        original = StoreGcCommand._dataset_candidates.__func__

        def _classify_then_wait(cls, *args, **kwargs):
            result = original(cls, *args, **kwargs)
            classified.set()
            assert locked.wait(timeout=5), "the holder thread never took the lock"
            return result

        monkeypatch.setattr(StoreGcCommand, "_dataset_candidates", classmethod(_classify_then_wait))

        def _hold():
            assert classified.wait(timeout=5), "classification never ran"
            with dataset_lock(paths, "SRR1"):
                locked.set()
                assert release.wait(timeout=5), "the main thread never released the holder"

        holder = threading.Thread(target=_hold)
        holder.start()
        try:
            rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        finally:
            release.set()
            holder.join(timeout=5)
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["removed_datasets"] == []
        assert [row["accession"] for row in report["in_use"]] == ["SRR1"]
        assert sra_dir(paths, "SRR1").is_dir()
        assert (sra_dir(paths, "SRR1") / "SRR1.fastq.gz").is_file()
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1") is not None

    def test_a_recently_touched_dataset_is_kept(self, tmp_path, capsys):
        """touch_dataset_use (set by a link, an adopt, or a store_link) keeps a dataset out
        of removal for GC_RECENT_USE_GRACE_SECONDS, even with no usage row recorded yet."""
        root = tmp_path / "store"
        paths = init_store(root)
        _write_dataset_dir(paths, "SRR1")
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1"))
        touch_dataset_use(paths, "SRR1")

        # The read-only classification pass does not know about touch_dataset_use, so SRR1
        # is still offered as an ordinary candidate here.
        dry = StoreGcCommand().execute(_gc_args(data_root=str(root), json=True))
        dry_report = json.loads(capsys.readouterr().out)
        assert dry == 0
        assert [row["accession"] for row in dry_report["datasets"]] == ["SRR1"]

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["removed_datasets"] == []
        in_use = {row["accession"]: row["reason"] for row in report["in_use"]}
        assert "recently used" in in_use["SRR1"]
        assert sra_dir(paths, "SRR1").is_dir()
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1") is not None

    def test_an_old_touch_no_longer_keeps_a_dataset(self, tmp_path, capsys):
        """Once the grace period has elapsed, a touch from long ago no longer protects the
        dataset from an otherwise ordinary removal."""
        import os

        root = tmp_path / "store"
        paths = init_store(root)
        _write_dataset_dir(paths, "SRR1")
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1"))
        touch_dataset_use(paths, "SRR1")
        marker = paths.locks / "SRR1.used"
        old = marker.stat().st_mtime - 90000  # just over a day ago
        os.utime(marker, (old, old))

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["removed_datasets"] == ["SRR1"]
        assert not sra_dir(paths, "SRR1").exists()

    def test_catalog_write_failure_after_the_rename_restores_the_folder(self, tmp_path, monkeypatch, caplog):
        """A catalogue write failure that happens after a candidate's folder has already been
        moved aside for removal must not leave that folder gone with no record of it: it is
        put back under its original name and the failure is logged."""
        root = tmp_path / "store"
        paths = init_store(root)
        acc_dir = _write_dataset_dir(paths, "SRR1")
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1"))

        def _boom(paths):
            raise DataAccessError("catalogue exploded")

        monkeypatch.setattr("prepare.cli.commands.store.gc.catalog_write", _boom)

        with caplog.at_level(logging.ERROR):
            rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True))

        assert rc == 1
        assert acc_dir.is_dir()
        assert (acc_dir / "SRR1.fastq.gz").is_file()
        assert not (paths.tmp / "SRR1_gc").exists()
        assert any("catalogue exploded" in record.message for record in caplog.records)
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1") is not None

    def test_a_dataset_republished_after_its_rename_keeps_its_catalogue_row(self, tmp_path, capsys, monkeypatch):
        """Between a candidate's rename (lock released) and the batched catalogue delete, a
        download publishes the accession again: the new folder and its row must stay."""
        root = tmp_path / "store"
        paths = init_store(root)
        _write_dataset_dir(paths, "SRR1")
        _write_dataset_dir(paths, "SRR2")
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1"))
            cat.upsert_dataset(_sidecar("SRR2"))
        original = StoreGcCommand._remove_one_dataset

        def _remove_then_republish(self, paths_arg, candidate, report, aside_by_accession, *rest):
            original(self, paths_arg, candidate, report, aside_by_accession, *rest)
            if candidate["accession"] == "SRR1":
                # A download_sra of SRR1 takes the now free lock and publishes a new copy.
                with dataset_lock(paths, "SRR1"):
                    _write_dataset_dir(paths, "SRR1")
                    with catalog_write(paths) as cat:
                        cat.upsert_dataset(_sidecar("SRR1", downloaded="2026-09-30T00:00:00+00:00"))

        monkeypatch.setattr(StoreGcCommand, "_remove_one_dataset", _remove_then_republish)
        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["removed_datasets"] == ["SRR2"]
        assert {row["accession"]: row["reason"] for row in report["in_use"]} == {
            "SRR1": "in use: published again during removal"
        }
        assert (sra_dir(paths, "SRR1") / "SRR1.fastq.gz").is_file()
        assert not (paths.tmp / "SRR1_gc").exists()
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR1")["downloaded"] == "2026-09-30T00:00:00+00:00"
            assert cat.get_dataset("SRR2") is None

    def test_a_catalogue_row_whose_folder_is_gone_is_removed(self, tmp_path, capsys, caplog):
        """A catalogue row with no dataset folder (removed by hand, or by an interrupted run)
        is dropped by --yes and reported, rather than staying a candidate on every run."""
        root = tmp_path / "store"
        paths = init_store(root)
        _write_dataset_dir(paths, "SRR1")
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "Proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "prepare_registry.json"))
            cat.upsert_dataset(_sidecar("SRR1"))
            cat.upsert_dataset(_sidecar("SRR9"))
        assert not sra_dir(paths, "SRR9").exists()

        with caplog.at_level(logging.INFO):
            rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert sorted(report["removed_datasets"]) == ["SRR1", "SRR9"]
        assert any("SRR9" in r.message and "folder already gone" in r.message for r in caplog.records)
        with Catalog(paths) as cat:
            assert cat.get_dataset("SRR9") is None
            assert cat.get_dataset("SRR1") is None
        assert not (paths.tmp / "SRR9_gc").exists()

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)
        assert rc == 0
        assert report["datasets"] == []

    def test_leftover_aside_from_an_interrupted_gc_is_swept_up_next_run(self, tmp_path, capsys):
        """A `<ACC>_gc` folder left behind by an interrupted removal (the rename succeeded,
        the catalogue delete or final cleanup did not) is picked up as an ordinary leftover
        by the next run, via the added `_gc` suffix in _accession_of_leftover."""
        root = tmp_path / "store"
        paths = init_store(root)
        with catalog_write(paths):
            pass
        aside = paths.tmp / "SRR1_gc"
        aside.mkdir(parents=True)
        (aside / "SRR1.fastq.gz").write_bytes(b"x" * 10)

        assert StoreGcCommand._accession_of_leftover("SRR1_gc") == "SRR1"

        rc = StoreGcCommand().execute(_gc_args(data_root=str(root), yes=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert any("SRR1_gc" in entry for entry in report["removed_leftovers"])
        assert not aside.exists()


def test_gc_include_stale_keeps_a_dataset_linked_by_a_relative_registry_project(tmp_path, monkeypatch, capsys):
    """Audit H2: a project bound with a relative ``--registry`` must not look stale from
    another working directory, and its live symlink must keep the dataset out of gc."""
    root = tmp_path / "store"
    init_store(root)
    paths = store_paths(root)
    acc_dir = _write_dataset_dir(paths, "SRR2")
    write_sidecar(sidecar_path(paths, "SRR2"), _sidecar("SRR2"))
    with catalog_write(paths) as cat:
        cat.upsert_dataset(_sidecar("SRR2"))

    proj = tmp_path / "proj"
    proj.mkdir()
    monkeypatch.chdir(proj)
    assert StoreInitCommand().execute(_init_args(root, proj, registry="prepare_registry.json")) == 0
    assert StoreLinkCommand().execute(_link_args(["SRR2"], data_root=str(root), registry="prepare_registry.json")) == 0
    used = root / "locks" / "SRR2.used"
    assert used.is_file()
    os.utime(used, (1, 1))  # past the 24 h grace

    capsys.readouterr()
    monkeypatch.chdir(tmp_path)
    rc = StoreGcCommand().execute(_gc_args(data_root=str(root), include_stale=True, yes=True, json=True))
    assert rc == 0
    assert acc_dir.is_dir()
    report = json.loads(capsys.readouterr().out)
    assert "SRR2" not in report["removed_datasets"]


class TestStoreGcLoadsRegistriesOnceAndPrunesTheJournal:
    """Classification reads every project's registry once per run, not once per dataset, and a
    removed accession leaves no usage line in the journal for a later replay to bring back."""

    @staticmethod
    def _count_stale_calls(monkeypatch):
        from prepare.cli.commands.store import gc as gc_module
        from prepare.store import usage as usage_module

        real = usage_module.stale_projects
        calls = []

        def counting(catalog):
            calls.append(1)
            return real(catalog)

        monkeypatch.setattr(usage_module, "stale_projects", counting)
        monkeypatch.setattr(gc_module, "stale_projects", counting)
        return calls

    @staticmethod
    def _store_with_datasets(tmp_path, count):
        paths = init_store(tmp_path / "store")
        stale_dir = tmp_path / "stale_proj"
        stale_dir.mkdir()
        with catalog_write(paths) as cat:
            cat.upsert_project("stale1", "OldProject", str(stale_dir), str(stale_dir / "prepare_registry.json"))
            for number in range(1, count + 1):
                accession = f"SRR{number}"
                _write_dataset_dir(paths, accession)
                cat.upsert_dataset(_sidecar(accession))
                if number % 2:
                    cat.record_usage(accession, "stale1", "wMel", "downloaded")
        return paths

    def test_stale_projects_is_called_once_for_a_50_dataset_classification(self, tmp_path, monkeypatch, capsys):
        paths = self._store_with_datasets(tmp_path, 50)
        calls = self._count_stale_calls(monkeypatch)

        rc = StoreGcCommand().execute(_gc_args(data_root=str(paths.root), json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert len(report["datasets"]) == 25
        assert len(report["kept_stale"]) == 25
        assert len(calls) == 1

    def test_removal_rechecks_staleness_once_per_batch(self, tmp_path, monkeypatch, capsys):
        paths = self._store_with_datasets(tmp_path, 10)
        calls = self._count_stale_calls(monkeypatch)

        rc = StoreGcCommand().execute(_gc_args(data_root=str(paths.root), yes=True, include_stale=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert len(report["removed_datasets"]) == 10
        # One for the classification and report, one for the removal batch's re-check.
        assert len(calls) == 2

    def test_gc_prunes_the_usage_journal_for_removed_accessions(self, tmp_path, capsys, caplog):
        paths = init_store(tmp_path / "store")
        stale_dir = tmp_path / "stale_proj"
        stale_dir.mkdir()
        live_dir = tmp_path / "live_proj"
        live_dir.mkdir()
        registry_file = live_dir / "prepare_registry.json"
        registry_file.write_text(json.dumps({"project": {"id": "live1"}}))
        for accession in ("SRR1", "SRR2"):
            _write_dataset_dir(paths, accession)
        with catalog_write(paths) as cat:
            cat.upsert_project("stale1", "OldProject", str(stale_dir), str(stale_dir / "prepare_registry.json"))
            cat.upsert_project("live1", "LiveProject", str(live_dir), str(registry_file))
            cat.upsert_dataset(_sidecar("SRR1"))
            cat.upsert_dataset(_sidecar("SRR2"))
            cat.record_usage("SRR1", "stale1", "wMel", "downloaded")
            cat.record_usage("SRR1", "stale1", "wMel", "linked")
            cat.record_usage("SRR2", "live1", "wMel", "downloaded")
        usage_file = paths.journal / "usage.jsonl"
        with open(usage_file, "a", encoding="utf-8") as handle:
            handle.write("{not json\n")

        with caplog.at_level(logging.WARNING):
            rc = StoreGcCommand().execute(_gc_args(data_root=str(paths.root), yes=True, include_stale=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["removed_datasets"] == ["SRR1"]
        lines = usage_file.read_text(encoding="utf-8").splitlines()
        assert "{not json" in lines
        records = [json.loads(line) for line in lines if line != "{not json"]
        assert [record["accession"] for record in records] == ["SRR2"]
        assert sum("could not be parsed" in record.message for record in caplog.records) == 1

    def test_reindex_compacts_the_projects_journal_to_the_latest_line_per_project(self, tmp_path):
        paths = init_store(tmp_path / "store")
        with catalog_write(paths) as cat:
            cat.upsert_project("p1", "first", str(tmp_path / "a"), "r1.json")
            cat.upsert_project("p2", "other", str(tmp_path / "b"), "r2.json")
            cat.upsert_project("p1", "second", str(tmp_path / "a"), "r1.json")
            cat.upsert_project("p1", "third", str(tmp_path / "a"), "r1.json")

        rc = StoreReindexCommand().execute(argparse.Namespace(data_root=str(paths.root), registry=None))

        assert rc == 0
        records = [
            json.loads(line) for line in (paths.journal / "projects.jsonl").read_text(encoding="utf-8").splitlines()
        ]
        assert [(record["project_id"], record["name"]) for record in records] == [("p2", "other"), ("p1", "third")]
        with Catalog(paths) as cat:
            names = dict(cat.conn.execute("SELECT project_id, name FROM projects").fetchall())
        assert names == {"p1": "third", "p2": "other"}
