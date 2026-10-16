"""Tests for the per-project dataset registry (prepare.data.registry)."""

import gzip
import json
import os
from pathlib import Path
from unittest.mock import patch

import pytest

from prepare.core.exceptions import DataAccessError
from prepare.data import registry as reg
from prepare.data import registry_batch as batch_mod
from prepare.data import registry_blocks as rb
from prepare.data.read_extraction import summarise_contigs
from prepare.processing.status_report import to_dataframes


def _fastq(path: Path, reads: int = 2, gz: bool = False) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(reads))
    if gz:
        with gzip.open(path, "wt") as handle:
            handle.write(body)
    else:
        path.write_text(body)
    return path


def _project(tmp_path: Path) -> reg.ProjectPaths:
    return reg.ProjectPaths(
        fastq=tmp_path / "fastq",
        metadata=tmp_path / "metadata",
        genomes=tmp_path / "genomes",
        targeted=tmp_path / "targeted",
        matches=tmp_path / "matches",
    )


class TestLoadSave:
    def test_registry_path_resolves_an_explicit_relative_path(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        assert reg.registry_path("prepare_registry.json") == (tmp_path / "prepare_registry.json").resolve()

    def test_absent_registry_is_empty(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        assert r.datasets == {} and r.version == reg.SCHEMA_VERSION and r.path == tmp_path / "prepare_registry.json"

    def test_round_trip_and_atomic_write(self, tmp_path):
        target = tmp_path / "prepare_registry.json"
        r = reg.load_registry(target)
        reg.record_exclusion(r, "SRR1", "amplicon", source="user")
        reg.save_registry(r)
        assert target.exists()
        assert not list(tmp_path.glob("*.tmp.*")) and not list(tmp_path.glob("*.lock"))
        again = reg.load_registry(target)
        assert again.datasets["SRR1"]["exclusion"]["reason"] == "amplicon"
        assert json.loads(target.read_text())["version"] == reg.SCHEMA_VERSION

    @pytest.mark.parametrize("text", ["[]", '"x"', "3"])
    def test_non_object_json_raises_data_access_error(self, tmp_path, text):
        target = tmp_path / "prepare_registry.json"
        target.write_text(text)
        with pytest.raises(DataAccessError, match="not a JSON object"):
            reg.load_registry(target)

    def test_non_numeric_version_raises_data_access_error(self, tmp_path):
        target = tmp_path / "prepare_registry.json"
        target.write_text('{"version": "x"}')
        with pytest.raises(DataAccessError, match="version is not a number"):
            reg.load_registry(target)

    def test_newer_schema_is_refused(self, tmp_path):
        target = tmp_path / "prepare_registry.json"
        target.write_text('{"version": %d}' % (reg.SCHEMA_VERSION + 1))
        with pytest.raises(DataAccessError, match="newer"):
            reg.load_registry(target)

    def test_unknown_top_level_keys_survive_a_round_trip(self, tmp_path):
        target = tmp_path / "prepare_registry.json"
        target.write_text('{"version": 2, "notes": {"a": 1}, "datasets": {}}')
        loaded = reg.load_registry(target)
        assert loaded.extra == {"notes": {"a": 1}}
        reg.save_registry(loaded)
        written = json.loads(target.read_text())
        assert written["notes"] == {"a": 1} and written["version"] == reg.SCHEMA_VERSION

    def test_invalid_json_raises(self, tmp_path):
        target = tmp_path / "prepare_registry.json"
        target.write_text("{not json")
        with pytest.raises(DataAccessError, match="not valid JSON"):
            reg.load_registry(target)

    def test_load_registry_wraps_os_error(self, tmp_path):
        """A read failure (permissions, I/O error) surfaces as DataAccessError, not a raw traceback."""
        target = tmp_path / "prepare_registry.json"
        target.write_text("{}")
        with patch.object(Path, "read_text", side_effect=OSError("permission denied")):
            with pytest.raises(DataAccessError, match="Cannot read registry"):
                reg.load_registry(target)

    def test_save_registry_wraps_os_error(self, tmp_path):
        """A write failure (disk full, permissions) surfaces as DataAccessError, not a raw traceback."""
        target = tmp_path / "prepare_registry.json"
        r = reg.load_registry(target)
        with patch("os.replace", side_effect=OSError("disk full")):
            with pytest.raises(DataAccessError, match="Cannot write registry"):
                reg.save_registry(r)
        # The lock and any temp file must not be left behind after the failure.
        assert not list(tmp_path.glob("*.lock"))
        assert not list(tmp_path.glob("*.tmp.*"))

    def test_registry_path_walks_up(self, tmp_path, monkeypatch):
        (tmp_path / "prepare_registry.json").write_text("{}")
        sub = tmp_path / "targeted" / "SRR1"
        sub.mkdir(parents=True)
        assert reg.registry_path(start=sub).resolve() == (tmp_path / "prepare_registry.json").resolve()
        assert reg.registry_path(explicit=str(tmp_path / "x.json")) == tmp_path / "x.json"

    def test_stale_lock_is_removed_and_fresh_lock_times_out(self, tmp_path, monkeypatch):
        target = tmp_path / "prepare_registry.json"
        lock = tmp_path / "prepare_registry.json.lock"
        monkeypatch.setattr(reg, "LOCK_STALE_SECONDS", 5.0)  # well apart from the wait window
        monkeypatch.setattr(reg, "LOCK_WAIT_SECONDS", 0.3)
        lock.write_text("1")
        os.utime(lock, (0, 0))  # ancient -> stale -> removed
        reg.save_registry(reg.load_registry(target))
        assert target.exists() and not lock.exists()
        lock.write_text("1")  # fresh lock held by "another process"
        with pytest.raises(DataAccessError, match="locked"):
            reg.save_registry(reg.load_registry(target))
        lock.unlink()

    def test_wait_deadline_is_its_own_constant(self, tmp_path, monkeypatch):
        """A long staleness window still gives up after the (separate) wait deadline."""
        target = tmp_path / "prepare_registry.json"
        lock = tmp_path / "prepare_registry.json.lock"
        monkeypatch.setattr(reg, "LOCK_STALE_SECONDS", 600.0)
        monkeypatch.setattr(reg, "LOCK_WAIT_SECONDS", 0.2)
        lock.write_text("1")
        with pytest.raises(DataAccessError, match="locked"):
            reg.save_registry(reg.load_registry(target))
        assert lock.exists()  # not reclaimed as stale
        lock.unlink()


class TestRegistrySchemaV2:
    def test_fresh_registry_has_empty_project_and_store(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        assert r.project == {}
        assert r.store == {}
        assert r.version == reg.SCHEMA_VERSION == 2

    def test_v1_file_loads_with_empty_project_and_store(self, tmp_path):
        target = tmp_path / "prepare_registry.json"
        v1_payload = {
            "version": 1,
            "created": "2026-01-01T00:00:00",
            "updated": "2026-01-01T00:00:00",
            "genomes": {},
            "datasets": {},
        }
        target.write_text(json.dumps(v1_payload))

        r = reg.load_registry(target)

        assert r.version == 1
        assert r.project == {}
        assert r.store == {}

    def test_v1_file_is_upgraded_to_v2_on_save(self, tmp_path):
        target = tmp_path / "prepare_registry.json"
        v1_payload = {
            "version": 1,
            "created": "2026-01-01T00:00:00",
            "updated": "2026-01-01T00:00:00",
            "genomes": {},
            "datasets": {},
        }
        target.write_text(json.dumps(v1_payload))

        r = reg.load_registry(target)
        reg.save_registry(r)

        written = json.loads(target.read_text())
        assert written["version"] == reg.SCHEMA_VERSION
        assert written["project"] == {}
        assert written["store"] == {}

    def test_project_and_store_round_trip(self, tmp_path):
        target = tmp_path / "prepare_registry.json"
        r = reg.load_registry(target)
        r.project = {"id": "abc-123", "name": "Wolbachia", "path": "/projects/wolbachia", "created": "2026-01-01"}
        r.store = {"root": "/data/store", "mode": "symlink", "linked": []}
        reg.save_registry(r)

        again = reg.load_registry(target)

        assert again.project == r.project
        assert again.store == r.store


class TestRegistryTransaction:
    def test_transaction_reloads_before_writing(self, tmp_path):
        """A change made between transactions survives, because each one loads from disk."""
        target = tmp_path / "prepare_registry.json"
        reg.save_registry(reg.load_registry(target))
        other = reg.load_registry(target)
        reg.record_exclusion(other, "SRR9", "amplicon")
        reg.save_registry(other)

        with reg.registry_transaction(target) as r:
            reg.record_download(r, "SRR1", "failed", tmp_path / "fastq", message="timeout")

        data = json.loads(target.read_text())
        assert data["datasets"]["SRR9"]["exclusion"]["excluded"] is True
        assert data["datasets"]["SRR1"]["download"]["state"] == "failed"
        assert not list(tmp_path.glob("*.lock")) and not list(tmp_path.glob("*.tmp.*"))

    def test_a_lock_taken_over_during_the_transaction_writes_nothing(self, tmp_path):
        from prepare.utils.lockfile import LockLost

        target = tmp_path / "prepare_registry.json"
        reg.save_registry(reg.load_registry(target))
        before = target.read_text()
        lock = target.with_name(target.name + ".lock")
        other = {"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00", "token": "feedbeef"}
        with pytest.raises(LockLost):
            with reg.registry_transaction(target) as r:
                reg.record_exclusion(r, "SRR1", "amplicon")
                # Another process reclaimed the lock while this holder was stalled.
                lock.write_text(json.dumps(other))
        assert target.read_text() == before
        assert json.loads(lock.read_text())["token"] == "feedbeef"

    def test_a_batch_keeps_its_queue_when_the_lock_was_lost(self, tmp_path):
        from prepare.utils.lockfile import LockLost

        target = tmp_path / "prepare_registry.json"
        lock = target.with_name(target.name + ".lock")
        other = {"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00", "token": "feedbeef"}
        batch = batch_mod.RegistryBatch(target, flush_every=None, flush_seconds=None)

        def _mutation(registry):
            reg.record_exclusion(registry, "SRR1", "amplicon")
            lock.write_text(json.dumps(other))

        batch.apply(_mutation, "SRR1")
        with pytest.raises(LockLost):
            batch.flush()
        assert batch.pending_labels() == ["SRR1"]
        assert not target.exists()

    def test_transaction_writes_nothing_on_error(self, tmp_path):
        target = tmp_path / "prepare_registry.json"
        with pytest.raises(ValueError, match="boom"):
            with reg.registry_transaction(target) as r:
                reg.record_exclusion(r, "SRR1", "amplicon")
                raise ValueError("boom")
        assert not target.exists()
        assert not list(tmp_path.glob("*.lock"))


class TestNanToNone:
    def test_nan_to_none(self):
        assert reg.nan_to_none(float("nan")) is None
        assert reg.nan_to_none("SRR1") == "SRR1"


class TestProjectRelativePaths:
    """A registry's paths are recorded relative to its own file, so a project can be moved."""

    def test_project_root_is_the_registry_files_parent(self, tmp_path):
        r = reg.load_registry(tmp_path / "sub" / "prepare_registry.json")
        assert reg.project_root(r) == (tmp_path / "sub").resolve()

    def test_project_root_falls_back_to_cwd_when_unbound(self, monkeypatch, tmp_path):
        monkeypatch.chdir(tmp_path)
        r = reg.Registry()
        assert r.path is None
        assert reg.project_root(r) == tmp_path.resolve()

    def test_resolve_project_path_round_trip(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        inside = tmp_path / "genomes" / "GCF_1.fna"
        stored = reg._project_relative(inside, reg.project_root(r))
        assert not Path(stored).is_absolute()
        assert reg.resolve_project_path(r, stored) == inside.resolve()

    def test_resolve_project_path_leaves_an_absolute_value_unchanged(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        outside = tmp_path.parent / "elsewhere" / "GCF_2.fna"
        assert reg.resolve_project_path(r, str(outside)) == outside

    def test_record_genome_stores_a_relative_path_under_the_project_root(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        fasta = tmp_path / "genomes" / "GCF_1.fna"
        manifest = tmp_path / "manifest.csv"
        reg.record_genome(r, "GCF_1", fasta, manifest)
        recorded = r.genomes["GCF_1"]["fasta"]
        assert recorded == "genomes/GCF_1.fna"
        assert reg.resolve_project_path(r, recorded) == fasta.resolve()
        assert r.genomes["GCF_1"]["manifest"] == "manifest.csv"

    def test_record_genome_with_empty_manifest_stores_empty_string(self, tmp_path):
        """A genome fetched without a manifest (e.g. genome_download) must not record a
        meaningless path resolved from wherever the process happened to run."""
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        fasta = tmp_path / "GCF_1.fna"
        reg.record_genome(r, "GCF_1", fasta, "")
        assert r.genomes["GCF_1"]["manifest"] == ""

    def test_record_genome_keeps_an_absolute_path_outside_the_project_root(self, tmp_path):
        r = reg.load_registry(tmp_path / "project" / "prepare_registry.json")
        outside_fasta = tmp_path / "shared_refs" / "GCF_1.fna"
        reg.record_genome(r, "GCF_1", outside_fasta, outside_fasta)
        recorded = r.genomes["GCF_1"]["fasta"]
        assert Path(recorded).is_absolute()
        assert Path(recorded) == outside_fasta.resolve()

    def test_recorded_download_files_are_project_relative(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        _fastq(tmp_path / "fastq" / "SRR1" / "SRR1_1.fastq")
        reg.record_download(r, "SRR1", "downloaded", tmp_path / "fastq")
        path = r.datasets["SRR1"]["download"]["files"][0]["path"]
        assert path == "fastq/SRR1/SRR1_1.fastq"

    def test_recorded_metadata_analysis_extraction_assembly_are_project_relative(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_metadata(r, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {})
        reg.record_analysis(r, "SRR1", "sra_stats", tmp_path / "reports" / "sra_statistics.csv", {})
        reg.record_extraction(
            r,
            "SRR1",
            "GCF_1",
            [tmp_path / "targeted" / "SRR1" / "GCF_1_1.fastq.gz"],
            10,
            False,
            {"genome_fasta": tmp_path / "genomes" / "GCF_1.fna", "preset": "sr", "threshold": 0.1},
        )
        reg.record_assembly(
            r, "SRR1", "GCF_1", tmp_path / "targeted" / "SRR1" / "GCF_1_assembly", {"contigs": 1}, "v1.2.9", {}
        )
        assert r.datasets["SRR1"]["metadata"]["xml"] == "metadata/SRR1_metadata.xml"
        assert r.datasets["SRR1"]["analyses"]["sra_stats"]["output"] == "reports/sra_statistics.csv"
        extraction = r.datasets["SRR1"]["extractions"]["GCF_1"]
        assert extraction["genome_fasta"] == "genomes/GCF_1.fna"
        assert extraction["files"] == ["targeted/SRR1/GCF_1_1.fastq.gz"]
        assert extraction["assembly"]["dir"] == "targeted/SRR1/GCF_1_assembly"

    def test_absolute_entries_from_an_older_registry_still_resolve(self, tmp_path):
        """A pre-existing registry with absolute paths (written before this change) is
        unaffected: resolve_project_path returns an absolute value unchanged."""
        target = tmp_path / "prepare_registry.json"
        absolute_fasta = tmp_path / "genomes" / "GCF_1.fna"
        target.write_text(
            json.dumps(
                {
                    "version": 1,
                    "created": "2026-01-01T00:00:00+00:00",
                    "updated": "2026-01-01T00:00:00+00:00",
                    "genomes": {"GCF_1": {"fasta": str(absolute_fasta), "manifest": ""}},
                    "datasets": {},
                }
            )
        )
        r = reg.load_registry(target)
        assert reg.resolve_project_path(r, r.genomes["GCF_1"]["fasta"]) == absolute_fasta


class TestRecords:
    def test_screening_selection_exclusion(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_screening(r, "SRR1", "GCF_1", 0.94611, 0.9974, "branchwater", 0.1, tmp_path / "bw.csv")
        reg.record_screening(r, "SRR1", "GCF_2", 0.2, None, "branchwater", 0.1, None)
        reg.record_selection(r, ["SRR1"], {"column": "GCF_1", "threshold": 0.5}, tmp_path / "accessions.txt")
        reg.record_selection(r, ["SRR2"], {"column": "GCF_1", "threshold": 0.9}, tmp_path / "accessions.txt")
        assert r.datasets["SRR1"]["screening"]["genomes"]["GCF_1"]["containment"] == 0.9461
        assert r.datasets["SRR1"]["screening"]["genomes"]["GCF_2"]["cani"] is None
        assert r.datasets["SRR1"]["selection"]["selected"] is False
        assert r.datasets["SRR2"]["selection"]["selected"] is True
        assert r.datasets["SRR2"]["selection"]["criteria"]["threshold"] == 0.9
        reg.record_exclusion(r, "SRR2", "16S amplicon")
        assert reg.query(r, "excluded") == ["SRR2"]
        reg.clear_exclusion(r, "SRR2")
        assert reg.query(r, "excluded") == []

    def test_download_records_files_and_attempts(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        _fastq(tmp_path / "fastq" / "SRR1" / "SRR1_1.fastq")
        _fastq(tmp_path / "fastq" / "SRR1" / "SRR1_2.fastq")
        reg.record_download(r, "SRR1", "failed", tmp_path / "fastq", message="timeout")
        reg.record_download(r, "SRR1", "downloaded", tmp_path / "fastq")
        dl = r.datasets["SRR1"]["download"]
        assert dl["state"] == "downloaded" and dl["attempts"] == 2
        assert [Path(f["path"]).name for f in dl["files"]] == ["SRR1_1.fastq", "SRR1_2.fastq"]
        assert dl["bytes_total"] == sum(f["bytes"] for f in dl["files"]) > 0
        reg.record_download(r, "SRR3", "skipped", tmp_path / "fastq", message="--max-downloads")
        assert r.datasets["SRR3"]["download"]["attempts"] == 0

        # A zero-byte file (e.g. left by an interrupted download) must never appear in
        # download["files"]; fastq_files() is the single source of truth here too.
        (tmp_path / "fastq" / "SRR2").mkdir(parents=True)
        (tmp_path / "fastq" / "SRR2" / "SRR2_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (tmp_path / "fastq" / "SRR2" / "SRR2_2.fastq").write_text("")
        reg.record_download(r, "SRR2", "downloaded", tmp_path / "fastq")
        names = [Path(f["path"]).name for f in r.datasets["SRR2"]["download"]["files"]]
        assert names == ["SRR2_1.fastq"]

        # A "downloaded" state recorded without an actual attempt (e.g. already present on
        # disk) must not bump attempts.
        _fastq(tmp_path / "fastq" / "SRR4" / "SRR4_1.fastq")
        reg.record_download(r, "SRR4", "downloaded", tmp_path / "fastq", attempt=False)
        assert r.datasets["SRR4"]["download"]["state"] == "downloaded"
        assert r.datasets["SRR4"]["download"]["attempts"] == 0

    def test_record_download_stores_and_preserves_the_completeness_verdict(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        verdict = {"verdict": "truncated", "reads_r1": 300000, "expected_spots": 48000000, "ratio": 0.0063}
        reg.record_download(r, "SRR1", "downloaded", tmp_path / "fastq", complete=verdict)
        assert r.datasets["SRR1"]["download"]["complete"] == verdict

        # A later call with no verdict (e.g. a skip recorded with no real attempt) must not
        # erase the verdict already on file.
        reg.record_download(r, "SRR1", "downloaded", tmp_path / "fastq", attempt=False)
        assert r.datasets["SRR1"]["download"]["complete"] == verdict

    def test_record_download_drops_the_cached_mate_read_counts(self, tmp_path):
        """The mate counts extract_target_reads caches describe the files this outcome
        replaces. A pair counted from a truncated download surviving the re-download that
        fixed it would force single-end mapping of a now-complete pair."""
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        download = reg.upsert_dataset(r, "SRR1").setdefault("download", {"attempts": 0})
        download["mate_reads"] = [10, 9]
        download["mate_reads_signature"] = [["SRR1_1.fastq.gz", 1, 1.0], ["SRR1_2.fastq.gz", 1, 1.0]]

        reg.record_download(r, "SRR1", "downloaded", tmp_path / "fastq")

        assert "mate_reads" not in r.datasets["SRR1"]["download"]
        assert "mate_reads_signature" not in r.datasets["SRR1"]["download"]

    def test_set_download_verdict_only_touches_complete(self, tmp_path):
        """Unlike record_download, set_download_verdict must not reset date, files, bytes_total,
        message, source or store_name; it only sets the completeness verdict."""
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_download(
            r, "SRR1", "downloaded", tmp_path / "fastq", message="original", source="store", store_name="SRR1"
        )
        original = dict(r.datasets["SRR1"]["download"])

        verdict = {"method": "spots", "ratio": 0.5, "verdict": "truncated"}
        reg.set_download_verdict(r, "SRR1", verdict)

        download = r.datasets["SRR1"]["download"]
        assert download["complete"] == verdict
        assert download["date"] == original["date"]
        assert download["files"] == original["files"]
        assert download["bytes_total"] == original["bytes_total"]
        assert download["message"] == original["message"]
        assert download["source"] == original["source"]
        assert download["store_name"] == original["store_name"]
        assert download["attempts"] == original["attempts"]

    def test_screening_source_and_threshold_are_per_genome(self, tmp_path):
        """A later screening from another source must not overwrite the first genome's provenance."""
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_screening(r, "SRR1", "GCF_1", 0.9, 0.99, "branchwater", 0.1, tmp_path / "bw.csv")
        reg.record_screening(r, "SRR1", "GCF_2", 0.3, None, "matches", 0.0, tmp_path / "m.csv")
        screening = r.datasets["SRR1"]["screening"]
        assert screening["genomes"]["GCF_1"]["source"] == "branchwater"
        assert screening["genomes"]["GCF_1"]["query_threshold"] == 0.1
        assert screening["genomes"]["GCF_2"]["source"] == "matches"
        assert "source" not in screening and "query_threshold" not in screening
        assert screening["date"]

    def test_record_extraction_stores_mapped_total(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_extraction(r, "SRR1", "GCF_1", [], 80, False, {}, mapped_total=100)
        assert r.datasets["SRR1"]["extractions"]["GCF_1"]["mapped_total"] == 100

    def test_record_extraction_mapped_total_defaults_to_none(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_extraction(r, "SRR1", "GCF_1", [], 80, False, {})
        assert r.datasets["SRR1"]["extractions"]["GCF_1"]["mapped_total"] is None

    def test_record_extraction_stores_reference_coverage(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        tsv = tmp_path / "targeted" / "SRR1" / "GCF_1_coverage.tsv"
        coverage = {"breadth": 0.875, "mean_depth": 8.0, "covered_bases": 350, "reference_bp": 400, "coverage_tsv": tsv}
        reg.record_extraction(r, "SRR1", "GCF_1", [], 80, False, {}, coverage=coverage)
        entry = r.datasets["SRR1"]["extractions"]["GCF_1"]
        assert entry["breadth"] == 0.875 and entry["mean_depth"] == 8.0
        assert entry["coverage_tsv"] == "targeted/SRR1/GCF_1_coverage.tsv"

    def test_record_extraction_coverage_defaults_to_none(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_extraction(r, "SRR1", "GCF_1", [], 80, False, {})
        entry = r.datasets["SRR1"]["extractions"]["GCF_1"]
        assert entry["breadth"] is None and entry["mean_depth"] is None and entry["coverage_tsv"] is None

    def test_record_assembly_passes_the_whole_stats_dict_through(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_extraction(r, "SRR1", "GCF_1", [], 100, False, {})
        reg.record_assembly(
            r,
            "SRR1",
            "GCF_1",
            tmp_path / "asm",
            {
                "contigs": 5,
                "total_bp": 5000,
                "n50": 800,
                "largest": 1200,
                "n90": 400,
                "gc": 0.512,
                "contigs_ge_1kb": 2,
                "genome_fraction_estimate": 0.75,
                "reads_mapped": 90,
                "mapping_rate": 0.9,
                "mean_depth_estimate": 12.3,
            },
            "v1.2.9",
            {"preset": "meta-sensitive"},
        )
        assembly = r.datasets["SRR1"]["extractions"]["GCF_1"]["assembly"]
        assert assembly["contigs"] == 5 and isinstance(assembly["contigs"], int)
        assert assembly["total_bp"] == 5000
        assert assembly["n50"] == 800
        assert assembly["largest"] == 1200
        assert assembly["n90"] == 400
        assert assembly["gc"] == 0.512
        assert assembly["contigs_ge_1kb"] == 2
        assert assembly["genome_fraction_estimate"] == 0.75
        assert assembly["reads_mapped"] == 90
        assert assembly["mapping_rate"] == 0.9
        assert assembly["mean_depth_estimate"] == 12.3
        assert assembly["params"]["preset"] == "meta-sensitive"

    def test_record_assembly_defaults_the_four_required_keys(self, tmp_path):
        """A stats dict missing some of the four required keys still gets them, defaulted to 0."""
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_extraction(r, "SRR1", "GCF_1", [], 0, False, {})
        reg.record_assembly(r, "SRR1", "GCF_1", tmp_path / "asm", {"contigs": 1}, "v1.2.9", {})
        assembly = r.datasets["SRR1"]["extractions"]["GCF_1"]["assembly"]
        assert assembly["contigs"] == 1
        assert assembly["total_bp"] == 0 and assembly["n50"] == 0 and assembly["largest"] == 0

    def test_clear_assembly_removes_the_block_but_keeps_the_extraction(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_extraction(r, "SRR1", "GCF_1", [tmp_path / "r.fastq.gz"], 10, False, {})
        reg.record_assembly(r, "SRR1", "GCF_1", tmp_path / "asm", {"contigs": 1}, "v1.2.9", {})
        reg.clear_assembly(r, "SRR1", "GCF_1")
        entry = r.datasets["SRR1"]["extractions"]["GCF_1"]
        assert entry["assembly"] is None
        assert entry["mapped_reads"] == 10  # the extraction record itself is untouched

    def test_clear_assembly_is_a_no_op_without_an_existing_record(self, tmp_path):
        """Nothing to clear (no extraction record at all yet) must not raise or create one."""
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.clear_assembly(r, "SRR1", "GCF_1")
        assert "SRR1" not in r.datasets

    def test_real_records_clear_the_inferred_flag(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_screening(r, "SRR1", "GCF_1", 0.5, None, "matches", 0.0, None)
        r.datasets["SRR1"]["screening"]["inferred"] = True
        reg.record_screening(r, "SRR1", "GCF_1", 0.9, 0.99, "branchwater", 0.1, None)
        assert "inferred" not in r.datasets["SRR1"]["screening"]

        reg.record_extraction(r, "SRR1", "GCF_1", [], 3, False, {})
        r.datasets["SRR1"]["extractions"]["GCF_1"]["inferred"] = True
        reg.record_assembly(r, "SRR1", "GCF_1", Path("d"), {"contigs": 1}, "v1.2.9", {})
        assert "inferred" not in r.datasets["SRR1"]["extractions"]["GCF_1"]

    def test_screening_cap_keeps_the_best(self, tmp_path, caplog):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        table = tmp_path / "parsed_containment.txt"
        table.write_text("\tGCF_A\nSRR1\t0.9\nSRR2\t0.2\nSRR3\t0.5\n")
        with caplog.at_level("WARNING"):
            reg.record_screening_from_table(r, table, tmp_path / "matches", max_screened=2)
        assert sorted(reg.query(r, "screened")) == ["SRR1", "SRR3"]
        assert "SRR2" not in r.datasets  # nothing else was recorded for it
        assert "2" in caplog.text and "screen" in caplog.text.lower()

    def test_screening_cap_keeps_accessions_with_other_records(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_exclusion(r, "SRR2", "16S amplicon")
        table = tmp_path / "parsed_containment.txt"
        table.write_text("\tGCF_A\nSRR1\t0.9\nSRR2\t0.2\n")
        reg.record_screening_from_table(r, table, tmp_path / "matches", max_screened=1)
        assert reg.query(r, "screened") == ["SRR1"]
        assert r.datasets["SRR2"]["exclusion"]["excluded"] is True and "screening" not in r.datasets["SRR2"]

    def test_screening_from_table_skips_non_numeric_cells(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        table = tmp_path / "parsed_containment.txt"
        table.write_text("\tGCF_A\tGCF_B\nSRR1\t0.9\tnot-a-number\n")
        count = reg.record_screening_from_table(r, table, tmp_path / "matches")
        assert count == 1
        assert list(r.datasets["SRR1"]["screening"]["genomes"]) == ["GCF_A"]

    def test_registry_path_logs_a_registry_above_the_working_directory(self, tmp_path, caplog):
        (tmp_path / "prepare_registry.json").write_text("{}")
        sub = tmp_path / "targeted"
        sub.mkdir()
        with caplog.at_level("INFO"):
            found = reg.registry_path(start=sub)
        assert found.resolve() == (tmp_path / "prepare_registry.json").resolve()
        assert str(found) in caplog.text

    def test_record_screening_from_table(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        table = tmp_path / "parsed_containment.txt"
        table.write_text(
            "\tGCF_A\tGCF_B\tmax_containment\tmax_containment_annotation\n"
            "SRR1\t0.9\t0.0\t0.9\tGCF_A\n"
            "SRR2\t0.0\t0.3\t0.3\tGCF_B\n"
        )
        matches_folder = tmp_path / "matches"
        count = reg.record_screening_from_table(r, table, matches_folder)
        assert count == 2
        assert r.datasets["SRR1"]["screening"]["genomes"]["GCF_A"]["containment"] == 0.9
        assert r.datasets["SRR1"]["screening"]["genomes"]["GCF_A"]["csv"] == "matches/GCF_A.csv"
        assert r.datasets["SRR1"]["screening"]["genomes"]["GCF_A"]["source"] == "matches"
        assert "GCF_B" not in r.datasets["SRR1"]["screening"]["genomes"]
        assert r.datasets["SRR2"]["screening"]["genomes"]["GCF_B"]["containment"] == 0.3
        assert "GCF_A" not in r.datasets["SRR2"]["screening"]["genomes"]

    def test_record_screening_from_table_missing_returns_zero(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        count = reg.record_screening_from_table(r, tmp_path / "missing.txt", tmp_path / "matches")
        assert count == 0
        assert r.datasets == {}

    def test_metadata_analysis_extraction_assembly(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_metadata(
            r, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {"run_size": "1234", "organism": "x"}
        )
        reg.record_analysis(r, "SRR1", "sra_stats", tmp_path / "sra_statistics.csv", {"total_reads": 10})
        reg.record_extraction(
            r,
            "SRR1",
            "GCF_1",
            [tmp_path / "t" / "GCF_1_1.fastq.gz"],
            239464,
            False,
            {"genome_fasta": "genomes/GCF_1.fna", "preset": "sr", "threshold": 0.1},
        )
        reg.record_assembly(
            r,
            "SRR1",
            "GCF_1",
            tmp_path / "t" / "GCF_1_assembly",
            {"contigs": 188, "total_bp": 1209849, "n50": 15400, "largest": 57491},
            "v1.2.9",
            {"threads": 1},
        )
        rec = rb.extraction_block(r, "SRR1", "GCF_1")
        assert rec.mapped_reads == 239464 and rec.assembly.n50 == 15400 and rec.assembly.version == "v1.2.9"
        assert r.datasets["SRR1"]["metadata"]["run_size"] == "1234"
        assert r.datasets["SRR1"]["analyses"]["sra_stats"]["summary"]["total_reads"] == 10
        # a new extraction record drops the assembly block: it was built from the replaced reads
        reg.record_extraction(r, "SRR1", "GCF_1", [], 0, True, {"genome_fasta": "g", "preset": "sr", "threshold": 0.1})
        assert rb.extraction_block(r, "SRR1", "GCF_1").assembly is None

    def test_record_metadata_keeps_spots_bases_layout_platform_as_int(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_metadata(
            r,
            "SRR1",
            tmp_path / "metadata" / "SRR1_metadata.xml",
            {
                "run_size": "4744553813",
                "run_md5": "abc",
                "run_total_spots": "47964651",
                "run_total_bases": "14389395300",
                "library_layout": "PAIRED",
                "platform": "ILLUMINA",
                "library_strategy": "WGS",
            },
        )
        metadata = r.datasets["SRR1"]["metadata"]
        assert metadata["run_total_spots"] == 47964651
        assert isinstance(metadata["run_total_spots"], int)
        assert metadata["run_total_bases"] == 14389395300
        assert isinstance(metadata["run_total_bases"], int)
        assert metadata["library_layout"] == "PAIRED"
        assert metadata["platform"] == "ILLUMINA"
        assert metadata["library_strategy"] == "WGS"

    def test_record_metadata_with_a_given_root_records_the_same_path(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        xml = tmp_path / "metadata" / "SRR1_metadata.xml"
        reg.record_metadata(r, "SRR1", xml, {})
        reg.record_metadata(r, "SRR2", xml, {}, root=reg.project_root(r))
        assert r.datasets["SRR1"]["metadata"]["xml"] == r.datasets["SRR2"]["metadata"]["xml"]
        assert r.datasets["SRR2"]["metadata"]["xml"] == "metadata/SRR1_metadata.xml"

    def test_record_metadata_spots_none_when_not_numeric(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_metadata(r, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {})
        metadata = r.datasets["SRR1"]["metadata"]
        assert metadata["run_total_spots"] is None
        assert metadata["run_total_bases"] is None
        assert metadata["library_layout"] is None
        assert metadata["platform"] is None
        assert metadata["library_strategy"] is None

    def test_query_and_stage_counts(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        for acc, cont in (("SRR1", 0.9), ("SRR2", 0.5), ("SRR3", 0.05)):
            reg.record_screening(r, acc, "GCF_1", cont, None, "branchwater", 0.01, None)
        reg.record_selection(r, ["SRR1", "SRR2"], {"threshold": 0.1}, Path("accessions.txt"))
        reg.record_download(r, "SRR1", "downloaded", tmp_path / "fastq")
        reg.record_extraction(r, "SRR1", "GCF_1", [Path("a")], 100, False, {})
        reg.record_extraction(r, "SRR2", "GCF_1", [], 0, False, {})
        reg.record_assembly(
            r, "SRR1", "GCF_1", Path("d"), {"contigs": 3, "total_bp": 9, "n50": 3, "largest": 4}, "v", {}
        )
        assert reg.query(r, "screened") == ["SRR1", "SRR2", "SRR3"]
        assert reg.query(r, "selected") == ["SRR1", "SRR2"]
        assert reg.query(r, "downloaded") == ["SRR1"]
        assert reg.query(r, "extracted", "GCF_1") == ["SRR1"]
        assert reg.query(r, "assembled", "GCF_1") == ["SRR1"]
        counts = reg.stage_counts(r)
        assert counts["stages"]["selected"] == 2 and counts["genomes"]["GCF_1"]["zero_mapped"] == ["SRR2"]
        with pytest.raises(DataAccessError, match="Unknown stage"):
            reg.query(r, "bogus")


class TestScanners:
    def test_summarise_contigs_uses_len_headers(self, tmp_path):
        fa = tmp_path / "final.contigs.fa"
        fa.write_text(
            ">k141_1 flag=1 multi=2.0 len=10\nACGTACGTAC\n>k141_2 flag=1 multi=2.0 len=4\nACGT\n>k141_3 len=6\nACGTAC\n"
        )
        assert summarise_contigs(fa) == {
            "contigs": 3,
            "total_bp": 20,
            "n50": 10,
            "n90": 4,
            "largest": 10,
            "gc": 0.5,
            "contigs_ge_1kb": 0,
        }
        fa.write_text(">a\nACGTACGT\n>b\nAC\n")
        assert summarise_contigs(fa) == {
            "contigs": 2,
            "total_bp": 10,
            "n50": 8,
            "n90": 2,
            "largest": 8,
            "gc": 0.5,
            "contigs_ge_1kb": 0,
        }
        assert summarise_contigs(tmp_path / "missing.fa") == {
            "contigs": 0,
            "total_bp": 0,
            "n50": 0,
            "n90": 0,
            "largest": 0,
            "gc": 0.0,
            "contigs_ge_1kb": 0,
        }

    def test_scan_downloads_ignores_transient_temp_folder(self, tmp_path):
        """A <acc>_temp folder holding a partial FASTQ must not count as a downloaded accession."""
        fastq_dir = tmp_path / "fastq"
        _fastq(fastq_dir / "SRR1" / "SRR1_1.fastq")
        _fastq(fastq_dir / "SRR2_temp" / "SRR2_temp_1.fastq")

        found = reg.scan_downloads(fastq_dir)

        assert list(found.keys()) == ["SRR1"]

    def test_scan_downloads_ignores_empty_and_partial_files(self, tmp_path):
        """A zero-byte file and a .gz.tmp.<pid> leftover are not downloaded reads."""
        fastq_dir = tmp_path / "fastq"
        good = _fastq(fastq_dir / "SRR1" / "SRR1_1.fastq")
        (fastq_dir / "SRR1" / "SRR1_2.fastq").write_text("")
        (fastq_dir / "SRR1" / "SRR1_1.fastq.gz.tmp.4242").write_text("partial")

        found = reg.scan_downloads(fastq_dir)

        assert found == {"SRR1": (1, good.stat().st_size)}

    def test_count_fastq_reads_plain_and_gz(self, tmp_path):
        assert reg.count_fastq_reads(_fastq(tmp_path / "a.fastq", reads=3)) == 3
        assert reg.count_fastq_reads(_fastq(tmp_path / "b.fastq.gz", reads=5, gz=True)) == 5

    def test_count_fastq_reads_across_chunk_boundaries(self, tmp_path):
        """The chunked binary reader must still count correctly for a file over 1 MiB."""
        record = "@r{0}\n" + "A" * 32 + "\n+\n" + "I" * 32 + "\n"
        path = tmp_path / "big.fastq"
        with open(path, "w") as handle:
            for i in range(20000):
                handle.write(record.format(i))
        assert path.stat().st_size > 1024 * 1024
        assert reg.count_fastq_reads(path) == 20000

    def test_split_extract_filename(self):
        genomes = ["GCF_000008025.1", "wMel_ref_1"]
        assert reg.split_extract_filename("GCF_000008025.1_1.fastq.gz", genomes) == ("GCF_000008025.1", "_1")
        assert reg.split_extract_filename("GCF_000008025.1.fastq.gz", genomes) == ("GCF_000008025.1", "")
        assert reg.split_extract_filename("wMel_ref_1_s.fastq.gz", genomes) == ("wMel_ref_1", "_s")
        assert reg.split_extract_filename("GCF_9_1.fastq.gz", genomes) == ("GCF_9", "_1")  # suffix fallback
        assert reg.split_extract_filename("notes.txt", genomes) is None

    def test_empty_assembly_dirs(self, tmp_path):
        paths = _project(tmp_path)
        empty_asm = paths.targeted / "SRR1" / "GCF_1_assembly"
        empty_asm.mkdir(parents=True)  # interrupted, no contigs file
        full_asm = paths.targeted / "SRR2" / "GCF_1_assembly"
        full_asm.mkdir(parents=True)
        (full_asm / "final.contigs.fa").write_text(">k len=5\nACGTA\n")
        (paths.targeted / "SRR3").mkdir(parents=True)  # no assembly folder at all

        pairs = reg.empty_assembly_dirs(paths.targeted, ["GCF_1"])
        assert pairs == [("SRR1", "GCF_1")]

    def test_bootstrap_handles_non_numeric_cani(self, tmp_path):
        paths = _project(tmp_path)
        paths.matches.mkdir(parents=True)
        (paths.matches / "GCF_1.csv").write_text("acc,containment,cANI\nSRR1,0.5,NA\nSRR2,0.6,\n")
        r = reg.bootstrap_from_disk(paths)
        assert r.datasets["SRR1"]["screening"]["genomes"]["GCF_1"]["containment"] == 0.5
        assert r.datasets["SRR1"]["screening"]["genomes"]["GCF_1"]["cani"] is None
        assert r.datasets["SRR2"]["screening"]["genomes"]["GCF_1"]["containment"] == 0.6
        assert r.datasets["SRR2"]["screening"]["genomes"]["GCF_1"]["cani"] is None

    def test_bootstrap_and_reconcile(self, tmp_path):
        paths = _project(tmp_path)
        _fastq(paths.fastq / "SRR1" / "SRR1_1.fastq")
        _fastq(paths.fastq / "SRR2" / "SRR2_1.fastq")
        (paths.metadata).mkdir()
        (paths.metadata / "SRR1_metadata.xml").write_text("<x/>")
        paths.genomes.mkdir()
        (paths.genomes / "GCF_1.fna").write_text(">c\nACGT\n")
        _fastq(paths.targeted / "SRR1" / "GCF_1_1.fastq.gz", gz=True)
        _fastq(paths.targeted / "SRR1" / "GCF_1_2.fastq.gz", gz=True)
        asm = paths.targeted / "SRR1" / "GCF_1_assembly"
        asm.mkdir(parents=True)
        (asm / "final.contigs.fa").write_text(">k len=5\nACGTA\n")
        (paths.targeted / "SRR2" / "GCF_1_assembly").mkdir(parents=True)  # interrupted, empty
        paths.matches.mkdir()
        (paths.matches / "GCF_1.csv").write_text("acc,containment,cANI\nSRR1,0.9,0.99\nSRR2,0.3,0.9\nSRR9,0.2,0.8\n")
        acc_file = tmp_path / "accessions.txt"
        acc_file.write_text("SRR1\nSRR2\n")

        r = reg.bootstrap_from_disk(paths, accessions_file=acc_file)
        assert r.datasets["SRR9"]["screening"]["genomes"]["GCF_1"]["containment"] == 0.2
        assert (
            r.datasets["SRR1"]["selection"]["selected"] is True and r.datasets["SRR1"]["selection"]["inferred"] is True
        )
        assert (
            r.datasets["SRR1"]["download"]["state"] == "downloaded"
            and r.datasets["SRR1"]["download"]["inferred"] is True
        )
        assert r.datasets["SRR1"]["metadata"]["inferred"] is True and "metadata" not in r.datasets["SRR2"]
        ext = rb.extraction_block(r, "SRR1", "GCF_1")
        # both mates are counted, so the inferred count matches what a real extraction records
        assert ext.mapped_reads == 4 and ext.assembly.contigs == 1 and ext.inferred is True
        assert "GCF_1" in r.genomes

        (paths.fastq / "SRR2" / "SRR2_1.fastq").unlink()
        _fastq(paths.fastq / "SRR7" / "SRR7_1.fastq")
        report = reg.reconcile(r, paths)
        assert report.recorded_missing == ["SRR2"]
        assert report.untracked_fastq == ["SRR7"]
        assert report.empty_assembly_dirs == [("SRR2", "GCF_1")]
        assert r.datasets["SRR2"]["download"]["state"] == "missing"

    def test_bootstrap_binds_the_registry_path_before_recording_so_paths_survive_a_different_cwd(
        self, tmp_path, monkeypatch
    ):
        """status --init with --registry pointing elsewhere must not record paths relative to
        the working directory: bootstrap_from_disk binds the registry to its intended file
        (via target_path) before any writer runs, so every recorded path is already relative to
        the right project root."""
        paths = _project(tmp_path)
        _fastq(paths.fastq / "SRR1" / "SRR1_1.fastq")
        elsewhere = tmp_path / "elsewhere"
        elsewhere.mkdir()
        monkeypatch.chdir(elsewhere)

        registry_file = tmp_path / "prepare_registry.json"
        r = reg.bootstrap_from_disk(paths, target_path=registry_file)
        reg.save_registry(r)

        loaded = reg.load_registry(registry_file)
        recorded = loaded.datasets["SRR1"]["download"]["files"][0]["path"]
        resolved = reg.resolve_project_path(loaded, recorded)
        assert resolved == (paths.fastq / "SRR1" / "SRR1_1.fastq").resolve()
        assert resolved.exists()

    def test_bootstrap_fills_coverage_from_the_coverage_table_on_disk(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        paths = _project(tmp_path)
        paths.genomes.mkdir()
        (paths.genomes / "GCF_1.fna").write_text(">c\nACGT\n")
        _fastq(paths.targeted / "SRR1" / "GCF_1_1.fastq.gz", gz=True)
        (paths.targeted / "SRR1" / "GCF_1_coverage.tsv").write_text(
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
            "chr\t1\t300\t10\t240\t80.0\t4.0\t30\t60\n"
            "plasmid\t1\t100\t2\t10\t10.0\t1.0\t30\t60\n"
        )
        _fastq(paths.targeted / "SRR2" / "GCF_1_1.fastq.gz", gz=True)

        r = reg.bootstrap_from_disk(paths)

        ext = rb.extraction_block(r, "SRR1", "GCF_1")
        assert ext.breadth == 0.625 and ext.mean_depth == 3.25
        assert ext.coverage_tsv == "targeted/SRR1/GCF_1_coverage.tsv"
        assert ext.inferred is True
        without = rb.extraction_block(r, "SRR2", "GCF_1")
        assert without.breadth is None and without.mean_depth is None and without.coverage_tsv is None

    def test_bootstrap_with_an_unreadable_coverage_table_records_no_coverage(self, tmp_path, monkeypatch, caplog):
        monkeypatch.chdir(tmp_path)
        paths = _project(tmp_path)
        paths.genomes.mkdir()
        (paths.genomes / "GCF_1.fna").write_text(">c\nACGT\n")
        _fastq(paths.targeted / "SRR1" / "GCF_1_1.fastq.gz", gz=True)
        (paths.targeted / "SRR1" / "GCF_1_coverage.tsv").write_text("not\ta coverage table\n")

        with caplog.at_level("WARNING"):
            r = reg.bootstrap_from_disk(paths)

        ext = rb.extraction_block(r, "SRR1", "GCF_1")
        assert ext.breadth is None and ext.mean_depth is None and ext.coverage_tsv is None
        assert ext.mapped_reads == 2
        assert "GCF_1_coverage.tsv" in caplog.text

    def test_bootstrap_from_disk_ignores_transient_temp_folder(self, tmp_path):
        """A <acc>_temp folder must never be recorded as a downloaded accession."""
        paths = _project(tmp_path)
        _fastq(paths.fastq / "SRR1" / "SRR1_1.fastq")
        _fastq(paths.fastq / "SRR2_temp" / "SRR2_temp_1.fastq")

        r = reg.bootstrap_from_disk(paths)

        assert r.datasets["SRR1"]["download"]["state"] == "downloaded"
        assert "SRR2_temp" not in r.datasets

    def test_to_dataframes(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_selection(r, ["SRR1"], {"threshold": 0.1}, Path("a.txt"))
        reg.record_extraction(r, "SRR1", "GCF_1", [], 7, False, {"preset": "sr"})
        datasets, extractions = to_dataframes(r)
        assert list(datasets.index) == ["SRR1"] and bool(datasets.loc["SRR1", "selected"]) is True
        assert extractions.loc[0, "genome_id"] == "GCF_1" and int(extractions.loc[0, "mapped_reads"]) == 7
        assert "breadth" in extractions.columns and "mean_depth" in extractions.columns

    def test_to_dataframes_carries_breadth_and_mean_depth(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        coverage = {"breadth": 0.5, "mean_depth": 3.25, "coverage_tsv": tmp_path / "c.tsv"}
        reg.record_extraction(r, "SRR1", "GCF_1", [], 7, False, {}, coverage=coverage)
        _, extractions = to_dataframes(r)
        assert extractions.loc[0, "breadth"] == 0.5 and extractions.loc[0, "mean_depth"] == 3.25

    def test_to_dataframes_reads_the_profile_or_an_old_quality_analysis(self, tmp_path):
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_analysis(r, "SRR1", "profile", tmp_path / "p.json", {"total_reads": 10, "gc_percent": 40.0})
        reg.record_analysis(r, "SRR2", "quality", tmp_path / "q.json", {"grade": "good", "gc_content": 0.5})
        datasets, _ = to_dataframes(r)
        assert datasets.loc["SRR1", "gc_percent"] == 40.0 and datasets.loc["SRR1", "total_reads"] == 10
        assert datasets.loc["SRR2", "gc_percent"] == 50.0 and datasets.loc["SRR2", "quality_grade"] == "good"


class TestStoreLinksInTheRegistry:
    """Store-backed downloads: recorded source, and dangling links reported as drift."""

    def _store_dataset(self, tmp_path, accession="SRR1"):
        from prepare.store.layout import init_store

        paths = init_store(tmp_path / "store")
        _fastq(paths.sra / accession / f"{accession}_1.fastq")
        return paths

    def test_record_download_keeps_source_and_store_name(self, tmp_path):
        paths = self._store_dataset(tmp_path)
        project = _project(tmp_path)
        from prepare.store.link import link_dataset

        link_dataset(project.fastq, "SRR1", paths)

        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_download(
            r,
            "SRR1",
            "downloaded",
            project.fastq,
            "linked from store, 1 files",
            attempt=False,
            source="store",
            store_name="SRR1",
        )

        download = r.datasets["SRR1"]["download"]
        assert download["source"] == "store"
        assert download["store_name"] == "SRR1"
        assert download["attempts"] == 0
        # The link resolves into the store, so the recorded path is the store's file.
        assert reg.resolve_project_path(r, download["files"][0]["path"]) == paths.sra / "SRR1" / "SRR1_1.fastq"

    def test_record_download_without_a_source_clears_a_stale_one(self, tmp_path):
        project = _project(tmp_path)
        _fastq(project.fastq / "SRR1" / "SRR1_1.fastq")
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_download(r, "SRR1", "downloaded", project.fastq, source="store", store_name="SRR1")

        reg.record_download(r, "SRR1", "downloaded", project.fastq)

        assert "source" not in r.datasets["SRR1"]["download"]
        assert "store_name" not in r.datasets["SRR1"]["download"]

    def test_reconcile_reports_dangling_links(self, tmp_path):
        import shutil

        from prepare.store.link import link_dataset

        paths = self._store_dataset(tmp_path)
        project = _project(tmp_path)
        link_dataset(project.fastq, "SRR1", paths)
        r = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_download(r, "SRR1", "downloaded", project.fastq, source="store", store_name="SRR1")
        shutil.rmtree(paths.sra / "SRR1")

        report = reg.reconcile(r, project)

        assert report.dangling_links == ["SRR1"]
        assert report.recorded_missing == ["SRR1"]

    def test_reconcile_reports_no_dangling_links_for_a_healthy_link(self, tmp_path):
        from prepare.store.link import link_dataset

        paths = self._store_dataset(tmp_path)
        project = _project(tmp_path)
        link_dataset(project.fastq, "SRR1", paths)
        r = reg.load_registry(tmp_path / "prepare_registry.json")

        report = reg.reconcile(r, project)

        assert report.dangling_links == []
        assert report.untracked_fastq == ["SRR1"]


class TestStoreBackedVerdictBackfill:
    """A record linked from the store gets its verdict from the store's own sidecar."""

    def test_a_store_backed_record_is_filled_in_from_the_sidecar(self, tmp_path):
        from prepare.data.registry import ProjectPaths, Registry
        from prepare.data.registry_reconcile import _fill_missing_download_verdicts
        from prepare.store.layout import init_store, sidecar_path
        from prepare.store.sidecar import Sidecar, write_sidecar

        paths_store = init_store(tmp_path / "store")
        write_sidecar(
            sidecar_path(paths_store, "SRR1"),
            Sidecar(
                accession="SRR1",
                state="complete",
                reads_per_mate=100,
                ncbi={"spots": 100},
                completeness={"method": "spots", "ratio": 1.0, "verdict": "complete"},
            ),
        )

        registry = Registry(path=tmp_path / "prepare_registry.json")
        registry.store = {"root": str(tmp_path / "store")}
        registry.datasets = {"SRR1": {"download": {"state": "downloaded", "source": "store"}}}

        _fill_missing_download_verdicts(registry, ProjectPaths(fastq=tmp_path / "fastq"))

        complete = registry.datasets["SRR1"]["download"]["complete"]
        assert complete["verdict"] == "complete"
        assert complete["expected_spots"] == 100
        assert complete["reads_r1"] == 100

    def test_a_store_backed_record_without_a_sidecar_is_left_alone(self, tmp_path):
        from prepare.data.registry import ProjectPaths, Registry
        from prepare.data.registry_reconcile import _fill_missing_download_verdicts
        from prepare.store.layout import init_store

        init_store(tmp_path / "store")
        registry = Registry(path=tmp_path / "prepare_registry.json")
        registry.store = {"root": str(tmp_path / "store")}
        registry.datasets = {"SRR1": {"download": {"state": "downloaded", "source": "store"}}}

        _fill_missing_download_verdicts(registry, ProjectPaths(fastq=tmp_path / "fastq"))

        assert "complete" not in registry.datasets["SRR1"]["download"]


class TestScannersIgnoreHiddenEntries:
    """``._*`` AppleDouble files and hidden folders are never extractions or assemblies."""

    def test_scan_extractions_ignores_appledouble_files_and_hidden_dirs(self, tmp_path):
        targeted = tmp_path / "targeted"
        acc_dir = targeted / "SRR1"
        acc_dir.mkdir(parents=True)
        (acc_dir / "GCF_1_1.fastq.gz").write_bytes(b"x")
        (acc_dir / "._GCF_1_1.fastq.gz").write_bytes(b"x" * 4096)
        hidden_dir = targeted / ".Trashes"
        hidden_dir.mkdir()
        (hidden_dir / "GCF_1_1.fastq.gz").write_bytes(b"x")

        found = reg.scan_extractions(targeted, [])

        assert list(found) == ["SRR1"]
        assert list(found["SRR1"]) == ["GCF_1"]
        assert [p.name for p in found["SRR1"]["GCF_1"]] == ["GCF_1_1.fastq.gz"]

    def test_scan_assemblies_ignores_hidden_entries(self, tmp_path):
        targeted = tmp_path / "targeted"
        (targeted / "SRR1" / "GCF_1_assembly").mkdir(parents=True)
        (targeted / "SRR1" / "._GCF_1_assembly").mkdir()
        (targeted / ".hidden" / "GCF_2_assembly").mkdir(parents=True)

        found = reg.scan_assemblies(targeted, [])

        assert found == {"SRR1": {"GCF_1": targeted / "SRR1" / "GCF_1_assembly"}}

    def test_genome_ids_on_disk_ignores_hidden_assembly_entries(self, tmp_path):
        paths = _project(tmp_path)
        (paths.targeted / "SRR1" / "GCF_1_assembly").mkdir(parents=True)
        (paths.targeted / "SRR1" / "._GCF_9_assembly").write_bytes(b"x" * 4096)

        ids = reg._genome_ids_on_disk(paths, None)

        assert "GCF_1" in ids
        assert not any(genome_id.startswith(".") for genome_id in ids)


class TestRecordExport:
    def test_record_export_under_project_exports(self, tmp_path):
        path = tmp_path / "prepare_registry.json"
        r = reg.load_registry(path)
        r.project = {"id": "p1", "name": "demo"}
        reg.record_export(r, "results_table", tmp_path / "out" / "results.tsv", {"rows": 3})
        reg.save_registry(r)
        loaded = reg.load_registry(path)
        export = loaded.project["exports"]["results_table"]
        assert export["output"] == "out/results.tsv"
        assert export["summary"] == {"rows": 3}
        assert export["date"]
        assert loaded.project["id"] == "p1"

    def test_to_int_or_none_is_public(self):
        assert reg.to_int_or_none("12") == 12
        assert reg.to_int_or_none("1.2G") is None


class TestUpdateLinked:
    """``update_linked`` is the one place the registry's ``store["linked"]`` list is changed."""

    def test_add_keeps_the_list_sorted_and_free_of_duplicates(self):
        registry = reg.Registry()
        reg.update_linked(registry, "SRR2", add=True)
        reg.update_linked(registry, "SRR1", add=True)
        reg.update_linked(registry, "SRR2", add=True)
        assert registry.store["linked"] == ["SRR1", "SRR2"]

    def test_remove_drops_the_accession_and_is_idempotent(self):
        registry = reg.Registry()
        registry.store["linked"] = ["SRR1", "SRR2"]
        reg.update_linked(registry, "SRR1", add=False)
        reg.update_linked(registry, "SRR1", add=False)
        assert registry.store["linked"] == ["SRR2"]

    def test_remove_from_an_absent_list_leaves_an_empty_list(self):
        registry = reg.Registry()
        reg.update_linked(registry, "SRR1", add=False)
        assert registry.store["linked"] == []


class TestStageCountsSpeed:
    """Stage checks read single fields rather than building typed blocks; catch a regression.

    The timing assertion for this moved to ``tests/test_performance_regressions.py`` (Task 9); this
    class keeps the non-timing regression guard, that ``stage_counts`` never converts a dataset into
    a typed ``RegistryBlock``.
    """

    def test_stage_counts_does_not_build_typed_blocks(self):
        r = reg.Registry()
        for i in range(100):
            accession = f"SRR{i:07d}"
            reg.record_screening(r, accession, "G1", 0.5, None, "matches", 0.0, None)
            r.datasets[accession]["download"] = {
                "attempts": 1,
                "state": "downloaded",
                "date": "d",
                "files": [{"path": "p", "bytes": 1, "mtime": "m"}, {"path": "q", "bytes": 1, "mtime": "m"}],
                "bytes_total": 2,
                "message": "",
                "complete": {"verdict": "complete", "ratio": 1.0},
            }
            if i % 10 == 0:
                r.datasets[accession]["extractions"] = {"G1": {"mapped_reads": i % 20, "files": [], "assembly": None}}
        conversions = []
        original = rb.RegistryBlock.from_dict.__func__

        def counting(cls, data):
            conversions.append(cls.__name__)
            return original(cls, data)

        with patch.object(rb.RegistryBlock, "from_dict", classmethod(counting)):
            counts = reg.stage_counts(r)
        assert counts["stages"]["downloaded"] == 100 and counts["genomes"]["G1"]["extracted"] == 5
        assert len(counts["genomes"]["G1"]["zero_mapped"]) == 5
        assert conversions == []


class TestCompactRegistryFile:
    def test_registry_is_written_compact_and_reloads(self, tmp_path):
        """The registry file carries no indentation, ends in a newline, and loads back unchanged."""
        path = tmp_path / reg.REGISTRY_FILENAME
        registry = reg.Registry(path=path)
        reg.record_screening(registry, "SRR1", "G1", 0.5, 0.9, "branchwater", 0.1, None)
        reg.record_exclusion(registry, "SRR2", "low quality")
        reg.save_registry(registry, path)
        text = path.read_text()
        assert "\n  " not in text
        assert text.endswith("}\n")
        assert text.count("\n") == 1
        assert reg.load_registry(path).datasets == registry.datasets


class TestRegistryBatch:
    def test_queued_mutations_are_written_in_one_transaction_on_exit(self, tmp_path):
        path = tmp_path / reg.REGISTRY_FILENAME
        with patch.object(reg, "_write_registry", wraps=reg._write_registry) as spy:
            with batch_mod.registry_batch(path) as batch:
                for index in range(10):
                    batch.apply(lambda r, i=index: reg.record_exclusion(r, f"SRR{i}", "test"))
                assert not path.exists()
        assert spy.call_count == 1
        assert sorted(reg.load_registry(path).datasets) == sorted(f"SRR{i}" for i in range(10))

    def test_mutations_apply_in_the_order_queued(self, tmp_path):
        path = tmp_path / reg.REGISTRY_FILENAME
        with batch_mod.registry_batch(path) as batch:
            batch.apply(lambda r: reg.record_exclusion(r, "SRR1", "first"))
            batch.apply(lambda r: reg.record_exclusion(r, "SRR1", "second"))
        assert rb.exclusion_block(reg.load_registry(path), "SRR1").reason == "second"

    def test_flushes_every_n_mutations(self, tmp_path):
        path = tmp_path / reg.REGISTRY_FILENAME
        with patch.object(reg, "_write_registry", wraps=reg._write_registry) as spy:
            with batch_mod.registry_batch(path, flush_every=3) as batch:
                for index in range(7):
                    batch.apply(lambda r, i=index: reg.record_exclusion(r, f"SRR{i}", "test"))
                assert spy.call_count == 2
                assert len(reg.load_registry(path).datasets) == 6
        assert spy.call_count == 3

    def test_flushes_when_the_interval_has_elapsed(self, tmp_path, monkeypatch):
        path = tmp_path / reg.REGISTRY_FILENAME
        clock = {"now": 100.0}
        with batch_mod.registry_batch(path, flush_every=1000, flush_seconds=30.0) as batch:
            monkeypatch.setattr(batch_mod, "_monotonic", lambda: clock["now"])
            batch.apply(lambda r: reg.record_exclusion(r, "SRR1", "test"))
            assert not path.exists()
            clock["now"] += 1e9
            batch.apply(lambda r: reg.record_exclusion(r, "SRR2", "test"))
            assert sorted(reg.load_registry(path).datasets) == ["SRR1", "SRR2"]

    def test_flushes_on_keyboard_interrupt_and_lets_it_propagate(self, tmp_path):
        path = tmp_path / reg.REGISTRY_FILENAME
        with pytest.raises(KeyboardInterrupt):
            with batch_mod.registry_batch(path) as batch:
                batch.apply(lambda r: reg.record_exclusion(r, "SRR1", "test"))
                raise KeyboardInterrupt
        assert "SRR1" in reg.load_registry(path).datasets
        assert not path.with_name(path.name + ".lock").exists()

    def test_a_failed_flush_on_interrupt_does_not_replace_the_interrupt(self, tmp_path, caplog):
        path = tmp_path / reg.REGISTRY_FILENAME
        with patch.object(reg, "_write_registry", side_effect=DataAccessError("disk gone")):
            with pytest.raises(KeyboardInterrupt):
                with batch_mod.registry_batch(path) as batch:
                    batch.apply(lambda r: reg.record_exclusion(r, "SRR1", "test"))
                    raise KeyboardInterrupt
        assert "disk gone" in caplog.text

    def test_a_failing_mutation_is_logged_and_dropped(self, tmp_path, caplog):
        """The mutations around a failing one are written, the queue empties, and one error names the label."""
        path = tmp_path / reg.REGISTRY_FILENAME

        def broken(registry):
            raise KeyError("no such block")

        with patch.object(reg, "_write_registry", wraps=reg._write_registry) as spy:
            with batch_mod.registry_batch(path) as batch:
                batch.apply(lambda r: reg.record_exclusion(r, "SRR1", "test"), label="SRR1")
                batch.apply(broken, label="SRR2")
                batch.apply(lambda r: reg.record_exclusion(r, "SRR3", "test"), label="SRR3")
                batch.flush()
                assert len(batch) == 0
                # A later mutation does not retry the dropped one.
                batch.apply(lambda r: reg.record_exclusion(r, "SRR4", "test"), label="SRR4")
        assert spy.call_count == 2
        assert batch.dropped == 1
        assert sorted(reg.load_registry(path).datasets) == ["SRR1", "SRR3", "SRR4"]
        errors = [r for r in caplog.records if r.levelname == "ERROR"]
        assert len(errors) == 1
        assert "SRR2" in errors[0].getMessage() and "no such block" in errors[0].getMessage()

    def test_a_mutation_failing_partway_leaves_no_trace(self, tmp_path, caplog):
        """What a failing mutation changed before raising is rolled back; the other mutations still land."""
        path = tmp_path / reg.REGISTRY_FILENAME
        reg.save_registry(reg.Registry(project={"name": "p"}), path)

        def half_applied(registry):
            reg.record_exclusion(registry, "SRR2", "first field")
            registry.project["touched"] = True
            raise KeyError("second field")

        with batch_mod.registry_batch(path) as batch:
            batch.apply(lambda r: reg.record_exclusion(r, "SRR1", "test"), label="SRR1")
            batch.apply(half_applied, label="SRR2")
            batch.apply(lambda r: reg.record_exclusion(r, "SRR3", "test"), label="SRR3")
        written = reg.load_registry(path)
        assert sorted(written.datasets) == ["SRR1", "SRR3"]
        assert written.project == {"name": "p"}
        errors = [r for r in caplog.records if r.levelname == "ERROR"]
        assert len(errors) == 1 and "SRR2" in errors[0].getMessage()

    def test_rollback_reapplies_the_earlier_mutations_only_once(self, tmp_path):
        """Mutations before the failing one are replayed from the file state, not applied twice."""
        path = tmp_path / reg.REGISTRY_FILENAME

        def count(registry):
            registry.project["count"] = registry.project.get("count", 0) + 1

        def broken(registry):
            registry.project["count"] = 99
            raise ValueError("boom")

        with batch_mod.registry_batch(path) as batch:
            batch.apply(count, label="a")
            batch.apply(count, label="b")
            batch.apply(broken, label="c")
            batch.apply(count, label="d")
        assert reg.load_registry(path).project["count"] == 3

    def test_a_failed_transaction_keeps_the_queue(self, tmp_path):
        path = tmp_path / reg.REGISTRY_FILENAME
        batch = batch_mod.registry_batch(path)
        batch.apply(lambda r: reg.record_exclusion(r, "SRR1", "test"))
        with patch.object(reg, "_write_registry", side_effect=DataAccessError("disk gone")):
            with pytest.raises(DataAccessError):
                batch.flush()
        assert len(batch) == 1
        batch.flush()
        assert "SRR1" in reg.load_registry(path).datasets

    def test_an_empty_batch_writes_nothing(self, tmp_path):
        path = tmp_path / reg.REGISTRY_FILENAME
        with batch_mod.registry_batch(path):
            pass
        assert not path.exists()

    def test_the_lock_is_held_only_during_a_flush(self, tmp_path):
        path = tmp_path / reg.REGISTRY_FILENAME
        with batch_mod.registry_batch(path) as batch:
            batch.apply(lambda r: reg.record_exclusion(r, "SRR1", "test"))
            # Another writer can still use the registry while mutations are queued.
            with reg.registry_transaction(path) as other:
                reg.record_exclusion(other, "SRR2", "other writer")
        assert sorted(reg.load_registry(path).datasets) == ["SRR1", "SRR2"]

    def test_flush_hooks_receive_the_written_registry(self, tmp_path):
        path = tmp_path / reg.REGISTRY_FILENAME
        seen = []
        with batch_mod.registry_batch(path) as batch:
            batch.add_flush_hook(lambda r: seen.append(sorted(r.datasets)))
            batch.apply(lambda r: reg.record_exclusion(r, "SRR1", "test"))
        assert seen == [["SRR1"]]
        assert batch.registry is not None and "SRR1" in batch.registry.datasets


class TestProjectFolders:
    """Project folders are anchored on the registry's folder (audit H4)."""

    def test_anchored_joins_relative_folders_onto_the_root(self, tmp_path):
        outside = tmp_path / "shared" / "genomes"
        paths = reg.ProjectPaths.anchored(tmp_path / "proj", fastq="reads", genomes=str(outside))
        assert paths.fastq == tmp_path / "proj" / "reads"
        assert paths.genomes == outside
        assert paths.metadata == tmp_path / "proj" / "metadata"
        assert paths.targeted == tmp_path / "proj" / "targeted"
        assert paths.matches == tmp_path / "proj" / "matches"

    def test_project_root_of_a_registry_file_is_its_folder(self, tmp_path, monkeypatch):
        (tmp_path / "elsewhere").mkdir()
        monkeypatch.chdir(tmp_path / "elsewhere")
        assert reg.project_root(tmp_path / "proj" / "prepare_registry.json") == (tmp_path / "proj").resolve()
        assert reg.project_root("prepare_registry.json") == (tmp_path / "elsewhere").resolve()
