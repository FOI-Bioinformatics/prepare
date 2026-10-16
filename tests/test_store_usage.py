"""
Tests for prepare.store.usage: the safe wrapper that records dataset usage in the
shared store's catalogue without ever letting a catalogue failure reach the caller, plus the
``stale_projects`` and ``linked_by`` helpers ``store_status``/``store_gc`` build on.
"""

import json
import os
from pathlib import Path
from unittest.mock import patch

import pytest

from prepare.data.registry import Registry
from prepare.store.catalog import Catalog, catalog_write
from prepare.store.layout import init_store
from prepare.store.usage import _upsert_project_row, linked_by, record_usage_many, record_usage_safe, stale_projects


@pytest.fixture
def paths(tmp_path):
    return init_store(tmp_path / "store")


def _registry(tmp_path, project_id="proj1", name="demo", path=None):
    registry = Registry(path=tmp_path / "prepare_registry.json")
    if project_id is not None:
        registry.project = {"id": project_id, "name": name, "path": path or str(tmp_path), "created": "now"}
    return registry


class TestRecordUsageSafe:
    def test_no_store_returns_false(self, tmp_path):
        registry = _registry(tmp_path)
        assert record_usage_safe(None, registry, "SRR1", "GCF_1", "downloaded") is False

    def test_no_project_id_returns_false(self, paths, tmp_path):
        registry = _registry(tmp_path, project_id=None)
        assert registry.project == {}
        assert record_usage_safe(paths, registry, "SRR1", "GCF_1", "downloaded") is False

    def test_no_project_id_logs_nothing_bad_and_writes_nothing(self, paths, tmp_path, caplog):
        registry = _registry(tmp_path, project_id=None)
        with caplog.at_level("WARNING"):
            result = record_usage_safe(paths, registry, "SRR1", "GCF_1", "downloaded")
        assert result is False
        with Catalog(paths, create=True) as catalog:
            catalog.migrate()
            rows = catalog.conn.execute("SELECT COUNT(*) AS n FROM usage").fetchone()
        assert rows["n"] == 0

    def test_records_row_with_stage_and_genome(self, paths, tmp_path):
        registry = _registry(tmp_path)
        assert record_usage_safe(paths, registry, "SRR1", "GCF_1", "downloaded", detail="stored") is True

        with Catalog(paths, create=True) as catalog:
            catalog.migrate()
            row = catalog.conn.execute(
                "SELECT * FROM usage WHERE accession = ? AND project_id = ?", ("SRR1", "proj1")
            ).fetchone()
        assert row["genome_id"] == "GCF_1"
        assert row["stage"] == "downloaded"
        assert row["detail"] == "stored"

    def test_default_genome_id_is_empty_string(self, paths, tmp_path):
        registry = _registry(tmp_path)
        assert record_usage_safe(paths, registry, "SRR1", "", "linked") is True

        with Catalog(paths, create=True) as catalog:
            catalog.migrate()
            row = catalog.conn.execute("SELECT genome_id FROM usage WHERE accession = ?", ("SRR1",)).fetchone()
        assert row["genome_id"] == ""

    def test_upserts_project_row(self, paths, tmp_path):
        registry = _registry(tmp_path, name="my-project", path=str(tmp_path / "proj"))
        record_usage_safe(paths, registry, "SRR1", "GCF_1", "extracted")

        with Catalog(paths, create=True) as catalog:
            catalog.migrate()
            row = catalog.conn.execute("SELECT * FROM projects WHERE project_id = ?", ("proj1",)).fetchone()
        assert row["name"] == "my-project"
        assert row["path"] == str(tmp_path / "proj")

    def test_moved_project_refreshes_the_projects_row(self, paths, tmp_path):
        """A project whose path changed between two calls has its catalogue row refreshed."""
        registry = _registry(tmp_path, path=str(tmp_path / "old"))
        record_usage_safe(paths, registry, "SRR1", "GCF_1", "downloaded")

        registry.project["path"] = str(tmp_path / "new")
        record_usage_safe(paths, registry, "SRR2", "GCF_1", "downloaded")

        with Catalog(paths, create=True) as catalog:
            catalog.migrate()
            row = catalog.conn.execute("SELECT path FROM projects WHERE project_id = ?", ("proj1",)).fetchone()
        assert row["path"] == str(tmp_path / "new")

    def test_catalog_write_raising_is_caught_and_logged(self, paths, tmp_path, caplog):
        registry = _registry(tmp_path)
        with patch("prepare.store.usage.catalog_write", side_effect=RuntimeError("locked")):
            with caplog.at_level("WARNING"):
                result = record_usage_safe(paths, registry, "SRR1", "GCF_1", "downloaded")
        assert result is False
        assert "SRR1" in caplog.text
        assert "downloaded" in caplog.text

    def test_second_call_updates_last_used_not_first_used(self, paths, tmp_path):
        registry = _registry(tmp_path)
        record_usage_safe(paths, registry, "SRR1", "GCF_1", "downloaded", detail="first")
        with Catalog(paths, create=True) as catalog:
            catalog.migrate()
            first = catalog.conn.execute(
                "SELECT first_used, last_used FROM usage WHERE accession = ?", ("SRR1",)
            ).fetchone()

        record_usage_safe(paths, registry, "SRR1", "GCF_1", "downloaded", detail="second")
        with Catalog(paths, create=True) as catalog:
            catalog.migrate()
            second = catalog.conn.execute(
                "SELECT first_used, last_used, detail FROM usage WHERE accession = ?", ("SRR1",)
            ).fetchone()

        assert second["first_used"] == first["first_used"]
        assert second["detail"] == "second"


class TestRecordUsageMany:
    def test_no_store_returns_false(self, tmp_path):
        registry = _registry(tmp_path)
        assert record_usage_many(None, registry, [("SRR1", "", "linked", "")]) is False

    def test_no_project_id_returns_false(self, paths, tmp_path):
        registry = _registry(tmp_path, project_id=None)
        assert record_usage_many(paths, registry, [("SRR1", "", "linked", "")]) is False

    def test_empty_rows_is_a_noop_success(self, paths, tmp_path):
        registry = _registry(tmp_path)
        assert record_usage_many(paths, registry, []) is True

    def test_records_every_row_in_one_write(self, paths, tmp_path):
        registry = _registry(tmp_path)
        rows = [
            ("SRR1", "GCF_1", "linked", "already downloaded"),
            ("SRR2", "GCF_1", "linked", "already downloaded"),
            ("SRR3", "", "linked", "already downloaded"),
        ]
        assert record_usage_many(paths, registry, rows) is True

        with Catalog(paths, create=True) as catalog:
            catalog.migrate()
            recorded = catalog.conn.execute(
                "SELECT accession, genome_id, stage FROM usage ORDER BY accession"
            ).fetchall()
        assert [(r["accession"], r["genome_id"], r["stage"]) for r in recorded] == [
            ("SRR1", "GCF_1", "linked"),
            ("SRR2", "GCF_1", "linked"),
            ("SRR3", "", "linked"),
        ]

    def test_one_catalog_write_call_for_the_whole_batch(self, paths, tmp_path):
        registry = _registry(tmp_path)
        rows = [("SRR1", "", "linked", ""), ("SRR2", "", "linked", "")]

        with patch(
            "prepare.store.usage.catalog_write",
            wraps=__import__("prepare.store.catalog", fromlist=["catalog_write"]).catalog_write,
        ) as wrapped:
            record_usage_many(paths, registry, rows)

        assert wrapped.call_count == 1

    def test_catalog_write_raising_is_caught_and_logged(self, paths, tmp_path, caplog):
        registry = _registry(tmp_path)
        with patch("prepare.store.usage.catalog_write", side_effect=RuntimeError("locked")):
            with caplog.at_level("WARNING"):
                result = record_usage_many(paths, registry, [("SRR1", "", "linked", "")])
        assert result is False


class TestStaleProjects:
    def test_unreadable_registry_directory_marks_project_stale_without_crashing(self, paths, tmp_path, caplog):
        """A registry file this process cannot even stat (e.g. a directory permissions
        problem on a shared multi-user store) must be treated as stale, not raise."""
        if os.name != "posix":
            pytest.skip("directory permission bits are not enforced the same way outside POSIX")
        if hasattr(os, "geteuid") and os.geteuid() == 0:
            pytest.skip("running as root bypasses directory permission checks")

        blocked_dir = tmp_path / "blocked"
        blocked_dir.mkdir()
        registry_file = blocked_dir / "prepare_registry.json"
        registry_file.write_text(json.dumps({"project": {"id": "proj1"}}))

        with catalog_write(paths) as cat:
            cat.upsert_project("proj1", "Blocked", str(blocked_dir), str(registry_file))

        original_mode = blocked_dir.stat().st_mode
        blocked_dir.chmod(0o000)
        try:
            with Catalog(paths) as cat:
                cat.migrate()
                with caplog.at_level("WARNING"):
                    stale = stale_projects(cat)
        finally:
            blocked_dir.chmod(original_mode)

        assert [row["project_id"] for row in stale] == ["proj1"]
        assert any(str(registry_file) in message for message in caplog.messages)


def test_relative_registry_is_stored_absolute(tmp_path, paths, monkeypatch):
    monkeypatch.chdir(tmp_path)
    registry = _registry(tmp_path, "p1", "proj", str(tmp_path))
    registry.path = Path("prepare_registry.json")
    with catalog_write(paths) as cat:
        _upsert_project_row(cat, registry, "p1")
        row = cat.conn.execute("SELECT registry FROM projects").fetchone()
    assert row["registry"] == str((tmp_path / "prepare_registry.json").resolve())


class TestLinkedBy:
    def test_no_projects_returns_empty(self, paths):
        with Catalog(paths, create=True) as cat:
            cat.migrate()
            assert linked_by(paths, cat, "SRR1") == []

    def test_live_project_with_matching_symlink_is_returned(self, paths, tmp_path):
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR1.fastq.gz").write_bytes(b"x")

        project_dir = tmp_path / "proj"
        (project_dir / "fastq").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1").symlink_to(acc_dir, target_is_directory=True)
        registry_file = project_dir / "prepare_registry.json"
        registry_file.write_text(json.dumps({"project": {"id": "proj1"}}))

        with catalog_write(paths) as cat:
            cat.upsert_project("proj1", "Live", str(project_dir), str(registry_file))

        with Catalog(paths, create=True) as cat:
            cat.migrate()
            assert linked_by(paths, cat, "SRR1") == ["Live"]

    def test_a_stale_projects_live_symlink_still_counts(self, paths, tmp_path):
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)

        project_dir = tmp_path / "proj"
        (project_dir / "fastq").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1").symlink_to(acc_dir, target_is_directory=True)
        # No registry file is ever written at this path, so the project is stale.
        with catalog_write(paths) as cat:
            cat.upsert_project("proj1", "proj", str(project_dir), str(project_dir / "prepare_registry.json"))

        with Catalog(paths, create=True) as cat:
            cat.migrate()
            assert linked_by(paths, cat, "SRR1") == ["proj (stale)"]

    def test_symlink_to_a_different_accession_is_ignored(self, paths, tmp_path):
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        other_dir = paths.sra / "SRR2"
        other_dir.mkdir(parents=True)

        project_dir = tmp_path / "proj"
        (project_dir / "fastq").mkdir(parents=True)
        (project_dir / "fastq" / "SRR1").symlink_to(other_dir, target_is_directory=True)
        registry_file = project_dir / "prepare_registry.json"
        registry_file.write_text(json.dumps({"project": {"id": "proj1"}}))

        with catalog_write(paths) as cat:
            cat.upsert_project("proj1", "Live", str(project_dir), str(registry_file))

        with Catalog(paths, create=True) as cat:
            cat.migrate()
            assert linked_by(paths, cat, "SRR1") == []
