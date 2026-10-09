"""The commands that keep a run log, each run through ``main`` in a project with the run log on.

Every test runs on ``tmp_path`` with the external tools and downloads replaced by fakes, and reads
back what ``main`` recorded through ``run_log.read_runs`` and ``run_log.read_detail``. A detail
holds rows only: mappings of plain values keyed by accession or ``accession/genome``, under one
section per kind of result (the form ``runs --diff`` and ``runs --accession`` compare).
"""

import argparse
import gzip
import json
from typing import Any, Dict
from unittest.mock import patch

import pytest

from prepare.cli.main import main
from prepare.core import settings
from prepare.data import run_log
from prepare.data.read_extraction import ExtractionResult
from prepare.data.registry import REGISTRY_FILENAME, load_registry, save_registry
from prepare.processing.doctor_report import OK, Check
from prepare.store.catalog import catalog_write
from prepare.store.layout import init_store, sidecar_path, sra_dir
from prepare.store.resolve import STORE_ENV
from prepare.store.sidecar import Sidecar, write_sidecar

from tests.helpers_extraction import tools_present  # noqa: F401 - autouse fixture

FASTQ = "@r1\nACGTACGTGC\n+\nIIIIIIIIII\n@r2\nGGCCAATTGC\n+\nIIIIIIIIII\n"


@pytest.fixture
def run_log_on(tmp_path, monkeypatch):
    """The run log on, and HOME and the store variable pointed away from the developer's own."""
    monkeypatch.setenv(settings.SETTINGS["run_log"].env, "true")
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    monkeypatch.delenv(STORE_ENV, raising=False)
    settings.reset_for_tests()
    yield
    settings.reset_for_tests()


@pytest.fixture
def project(tmp_path, monkeypatch, run_log_on):
    """A project folder with an empty registry, as the working directory."""
    folder = tmp_path / "project"
    folder.mkdir()
    save_registry(load_registry(folder / REGISTRY_FILENAME))
    monkeypatch.chdir(folder)
    return folder


def _only_run(project, command):
    records = run_log.read_runs(project, command)
    assert len(records) == 1, [r.command for r in run_log.read_runs(project)]
    return records[0]


def _no_runs(project):
    return not (run_log.runs_dir(project) / run_log.RUNS_FILE).exists()


def _assert_rows(section: Dict[str, Any]) -> None:
    """Every entry is a row: a non-empty mapping of plain values."""
    assert section
    for row in section.values():
        assert isinstance(row, dict) and row
        assert all(not isinstance(value, (dict, list)) for value in row.values()), row


def _write_fastq(path, text=FASTQ):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(text)


# ---------------------------------------------------------------- download_sra

NETWORK = "Network error: connection reset by peer"


def _fake_download(**kwargs):
    kwargs["on_result"]("SRR2", False, NETWORK)
    return {
        "total": 1,
        "successful": 0,
        "failed": 1,
        "failed_accessions": ["SRR2"],
        "results": {"SRR2": NETWORK},
        "aborted": None,
    }


def _download(project, *extra):
    (project / "acc.txt").write_text("SRR2\n")
    with (
        patch("prepare.cli.commands.sra.require_tools"),
        patch("prepare.cli.commands.sra.download_sra", side_effect=_fake_download),
    ):
        return main(["download_sra", "--accessions-file", "acc.txt", "--max-retries", "0", *extra])


def test_download_sra_records_its_run_document_totals(project):
    assert _download(project) == 4
    record = _only_run(project, "download_sra")
    assert record.exit_code == 4
    assert record.summary["accessions"] == 1
    assert record.summary["failed"] == 1
    assert record.summary["downloaded"] == 0
    assert record.summary["failures_by_reason"] == {"network": 1}
    assert record.summary["aborted"] is None
    detail = run_log.read_detail(project, record)
    assert set(detail) == {"failed"}
    _assert_rows(detail["failed"])
    assert detail["failed"]["SRR2"]["reason"] == "network"
    assert detail["failed"]["SRR2"]["attempts"] == 1


def _fake_mixed_download(**kwargs):
    kwargs["timings"]["SRR1"] = ("2026-10-01T12:00:00+00:00", 2.5)
    kwargs["on_result"]("SRR1", True, "Downloaded SRR1")
    kwargs["on_result"]("SRR2", False, NETWORK)
    return {
        "total": 2,
        "successful": 1,
        "failed": 1,
        "failed_accessions": ["SRR2"],
        "results": {"SRR1": "Downloaded SRR1", "SRR2": NETWORK},
        "aborted": None,
    }


def test_download_sra_detail_has_a_row_per_downloaded_accession(project, capsys):
    (project / "acc.txt").write_text("SRR1\nSRR2\n")
    with (
        patch("prepare.cli.commands.sra.require_tools"),
        patch("prepare.cli.commands.sra.download_sra", side_effect=_fake_mixed_download),
    ):
        main(["download_sra", "--accessions-file", "acc.txt", "--max-retries", "0"])
    record = _only_run(project, "download_sra")
    detail = run_log.read_detail(project, record)
    assert set(detail) == {"failed", "downloaded"}
    _assert_rows(detail["downloaded"])
    assert detail["downloaded"] == {"SRR1": {"status": "downloaded", "seconds": 2.5, "attempts": 1}}
    assert set(detail["failed"]) == {"SRR2"}

    # runs --accession now finds the run that fetched SRR1, with unprefixed field names.
    capsys.readouterr()
    assert main(["runs", "--accession", "SRR1", "--json"]) == 0
    history = json.loads(capsys.readouterr().out)["runs"]
    assert [entry["command"] for entry in history] == ["download_sra"]
    assert history[0]["values"] == {"SRR1": {"status": "downloaded", "seconds": 2.5, "attempts": 1}}


def test_a_dry_run_records_nothing(project):
    def fake(**kwargs):
        return {"total": 1, "to_download": 1, "already_downloaded": 0, "blacklisted": 0}

    (project / "acc.txt").write_text("SRR2\n")
    with patch("prepare.cli.commands.sra.download_sra", side_effect=fake):
        assert main(["download_sra", "--accessions-file", "acc.txt", "--dry-run"]) == 0
    assert _no_runs(project)


# ---------------------------------------------------------------- sra_profile and sra_report


@pytest.fixture
def fastq_project(project):
    _write_fastq(project / "fastq" / "SRR1" / "SRR1.fastq.gz")
    return project


def test_sra_profile_records_profile_rows(fastq_project):
    args = ["sra_profile", "--accession", "SRR1", "--output-report", "stats.csv", "--output-dir", "profiles"]
    assert main([*args, "--summary-only"]) == 0
    record = _only_run(fastq_project, "sra_profile")
    assert record.summary["accessions"] == 1
    assert record.summary["profiled"] == 1
    assert record.summary["failed"] == 0
    assert record.summary["total_reads"] == 2
    detail = run_log.read_detail(fastq_project, record)
    assert set(detail) == {"analyses"} and set(detail["analyses"]) == {"profile"}
    rows = detail["analyses"]["profile"]
    _assert_rows(rows)
    assert rows["SRR1"]["total_reads"] == 2
    assert rows["SRR1"]["gc_percent"] == pytest.approx(60.0)


def test_sra_report_records_report_rows(fastq_project):
    (fastq_project / "acc.txt").write_text("SRR1\n")
    argv = ["sra_report", "--accessions-file", "acc.txt", "--no-report", "--no-open", "--output-dir", "report"]
    assert main(argv) == 0
    record = _only_run(fastq_project, "sra_report")
    assert record.summary["accessions"] == 1
    assert record.summary["reported"] == 1
    assert record.summary["failed"] == 0
    assert record.summary["html"] is False
    detail = run_log.read_detail(fastq_project, record)
    assert set(detail) == {"analyses"} and set(detail["analyses"]) == {"report"}
    rows = detail["analyses"]["report"]
    _assert_rows(rows)
    assert rows["SRR1"]["total_reads"] == 2
    assert "quality_grade" in rows["SRR1"]


# ---------------------------------------------------------------- results_table


def test_results_table_records_the_export_summary(project):
    assert main(["results_table", "--output", "results.tsv"]) == 0
    record = _only_run(project, "results_table")
    assert record.summary["rows"] == 0
    assert record.summary["output"] == "results.tsv"
    assert {"accessions", "genomes", "genome_id", "min_containment"} <= set(record.summary)
    assert record.detail is None


def test_results_table_with_no_record_records_nothing(project):
    assert main(["results_table", "--output", "results.tsv", "--no-record"]) == 0
    assert _no_runs(project)


# ---------------------------------------------------------------- extract_target_reads


def _extract(project, *extra):
    (project / "parsed_containment.txt").write_text("accession\tG1\nSRR1\t0.0\n")
    (project / "g1.fasta").write_text(">c\nACGT\n")
    reads = project / "targeted" / "G1" / "SRR1_1.fastq"
    reads.parent.mkdir(parents=True)
    reads.write_text(FASTQ)

    def fake(**kwargs):
        result = ExtractionResult(
            files=[reads],
            mapped_records=7,
            mapped_total=8,
            coverage={"breadth": 0.5, "mean_depth": 2.25, "covered_bases": 2, "reference_bp": 4},
        )
        if not kwargs["dry_run"]:
            kwargs["on_result"]("SRR1", result)
        return {"SRR1": result}

    argv = ["extract_target_reads", "--parsed-containment", "parsed_containment.txt", "--genome-id", "G1"]
    argv += ["--genome-fasta", "g1.fasta", "--threshold", "0.1"]
    with patch("prepare.cli.commands.read_extraction.extract_target_reads", side_effect=fake):
        return main([*argv, *extra])


def test_extract_target_reads_records_one_row_per_accession_and_genome(project):
    assert _extract(project) == 0
    record = _only_run(project, "extract_target_reads")
    assert record.summary == {
        "genome_id": "G1",
        "samples": 1,
        "extracted": 1,
        "with_reads": 1,
        "skipped": 0,
    }
    detail = run_log.read_detail(project, record)
    assert detail == {"extractions": {"SRR1/G1": {"mapped_reads": 7, "breadth": 0.5, "mean_depth": 2.25}}}


def test_extract_target_reads_notes_a_skipped_sample_with_its_recorded_values(project, capsys):
    assert _extract(project) == 0
    seen = {}

    def skipping(**kwargs):
        seen.update(kwargs["already_done"])
        result = ExtractionResult(files=[project / "targeted" / "G1" / "SRR1_1.fastq"], mapped_records=7, skipped=True)
        kwargs["on_result"]("SRR1", result)
        return {"SRR1": result}

    argv = ["extract_target_reads", "--parsed-containment", "parsed_containment.txt", "--genome-id", "G1"]
    argv += ["--genome-fasta", "g1.fasta", "--threshold", "0.1"]
    with patch("prepare.cli.commands.read_extraction.extract_target_reads", side_effect=skipping):
        assert main(argv) == 0
    assert "SRR1" in seen
    latest = run_log.read_runs(project, "extract_target_reads")[-1]
    assert latest.summary["skipped"] == 1
    detail = run_log.read_detail(project, latest)
    assert detail == {
        "extractions": {"SRR1/G1": {"mapped_reads": 7, "breadth": 0.5, "mean_depth": 2.25, "skipped": True}}
    }
    # runs --diff no longer reads the skipped sample as removed.
    capsys.readouterr()
    assert main(["runs", "--command", "extract_target_reads", "--diff", "previous", "latest", "--json"]) == 0
    diff = json.loads(capsys.readouterr().out)["detail"]
    assert diff["removed"] == [] and diff["added"] == []
    assert diff["changed"] == {"SRR1/G1": {"skipped": [None, True]}}


def test_extract_target_reads_dry_run_records_nothing(project):
    assert _extract(project, "--dry-run") == 0
    assert _no_runs(project)


# ---------------------------------------------------------------- store_verify


def test_store_verify_records_counts_and_fixed(project, tmp_path):
    root = tmp_path / "store"
    paths = init_store(root)
    fastq = sra_dir(paths, "SRR1") / "SRR1.fastq.gz"
    _write_fastq(fastq)
    sidecar = Sidecar(
        accession="SRR1",
        state="partial",
        layout="SINGLE",
        downloaded="2026-09-06T00:00:00+00:00",
        tool="fasterq-dump",
        tool_version="3.0.0",
        compression="gzip",
        files=[{"name": "SRR1.fastq.gz", "bytes": fastq.stat().st_size, "md5": "ignored", "reads": 2}],
        reads_per_mate=2,
        bases_total=20,
        ncbi={"spots": 2, "bases": 20, "size": 100, "layout": "SINGLE", "files": []},
        completeness={"method": "spots", "ratio": 1.0, "verdict": "complete"},
    )
    write_sidecar(sidecar_path(paths, "SRR1"), sidecar)
    with catalog_write(paths) as cat:
        cat.upsert_dataset(sidecar)

    assert main(["store_verify", "--data-root", str(root), "--spots", "--fix-state"]) == 0
    record = _only_run(project, "store_verify")
    assert record.summary["datasets"] == 1
    assert record.summary["counts"] == {"ok": 1}
    assert record.summary["fixed"] == 1
    detail = run_log.read_detail(project, record)
    assert set(detail) == {"verify"}
    _assert_rows(detail["verify"])
    row = detail["verify"]["SRR1"]
    assert row["verdict"] == "ok"
    assert row["state_before"] == "partial"
    assert row["state_after"] == "complete"
    assert row["fix"] == "updated"


# ---------------------------------------------------------------- project_report


def test_project_report_records_its_summary(project):
    argv = ["project_report", "--html", "never", "--no-environment", "--output-dir", "report"]
    assert main(argv) == 0
    record = _only_run(project, "project_report")
    assert record.summary["files"] == ["project_report.md", "project_report.json"]
    assert record.summary["html"] is False
    assert {"datasets", "extraction_rows", "failed_downloads", "environment"} <= set(record.summary)


def test_project_report_with_no_record_records_nothing(project):
    argv = ["project_report", "--html", "never", "--no-environment", "--output-dir", "report", "--no-record"]
    assert main(argv) == 0
    assert _no_runs(project)


# ---------------------------------------------------------------- select_datasets


def test_select_datasets_records_the_selection(project):
    (project / "parsed_containment.txt").write_text("accession\tG1\nSRR1\t0.9\nSRR2\t0.5\nSRR3\t0.01\n")
    assert main(["select_datasets", "--genome-id", "G1", "--threshold", "0.1"]) == 0
    record = _only_run(project, "select_datasets")
    assert record.summary["selected"] == 2
    assert record.summary["column"] == "G1"
    assert record.summary["threshold"] == 0.1
    detail = run_log.read_detail(project, record)
    assert set(detail) == {"selection"}
    _assert_rows(detail["selection"])
    assert detail["selection"]["SRR1"] == {"rank": 1, "column": "G1", "value": 0.9}
    assert detail["selection"]["SRR2"]["rank"] == 2


def test_select_datasets_with_no_record_records_nothing(project):
    (project / "parsed_containment.txt").write_text("accession\tG1\nSRR1\t0.9\n")
    argv = ["select_datasets", "--genome-id", "G1", "--threshold", "0.1", "--output", "scratch.txt", "--no-record"]
    assert main(argv) == 0
    assert _no_runs(project)


# ---------------------------------------------------------------- blacklist


def test_blacklist_records_what_it_changed(project):
    assert main(["blacklist", "--add", "SRR9", "SRR8", "--reason", "host reads"]) == 0
    assert main(["blacklist", "--remove", "SRR8"]) == 0
    added, removed = run_log.read_runs(project, "blacklist")
    assert added.summary == {"action": "add", "accessions": 2, "reason": "host reads", "excluded": 2}
    assert removed.summary == {"action": "remove", "accessions": 1, "reason": None, "excluded": 1}


def test_blacklist_list_records_nothing(project):
    assert main(["blacklist", "--list"]) == 0
    assert _no_runs(project)


# ---------------------------------------------------------------- status


def test_status_init_records_the_new_registry(tmp_path, monkeypatch, run_log_on):
    folder = tmp_path / "fresh"
    folder.mkdir()
    monkeypatch.chdir(folder)
    _write_fastq(folder / "fastq" / "SRR1" / "SRR1.fastq.gz")
    assert main(["status", "--init"]) == 0
    record = _only_run(folder, "status")
    assert record.summary["action"] == "init"
    assert record.summary["datasets"] == 1
    assert "drift" not in record.summary


def test_status_reconcile_records_drift_counts(project):
    _write_fastq(project / "fastq" / "SRR1" / "SRR1.fastq.gz")
    assert main(["status", "--reconcile"]) == 0
    record = _only_run(project, "status")
    assert record.summary["action"] == "reconcile"
    assert record.summary["drift"]["untracked_fastq"] == 1
    assert all(isinstance(count, int) for count in record.summary["drift"].values())


# ---------------------------------------------------------------- read-only commands


def test_read_only_invocations_record_nothing(project, capsys):
    assert main(["status"]) == 0
    assert main(["status", "--json"]) == 0
    with patch("prepare.cli.commands.doctor.run_checks", return_value=[Check("python", OK, "3.12")]):
        assert main(["doctor"]) == 0
    assert _no_runs(project)
    # runs exits 1 without a run log, and still writes none.
    assert main(["runs"]) == 1
    assert _no_runs(project)


@pytest.mark.parametrize("error", [AttributeError("no such argument"), KeyError("missing")])
@pytest.mark.parametrize(
    "argv, code",
    [(["blacklist", "--add", "SRR9", "--reason", "host reads"], 0), (["blacklist", "--add", "SRR9"], 1)],
    ids=["exits-0", "fails"],
)
def test_a_records_run_override_that_raises_keeps_the_exit_code(project, caplog, error, argv, code):
    def broken(self, args):
        raise error

    with patch("prepare.cli.commands.blacklist.BlacklistCommand.records_run", broken):
        assert main(argv) == code
    assert _no_runs(project)
    assert "not recorded in the run log" in caplog.text


# ---------------------------------------------------------------- note_rows


def test_note_rows_keeps_earlier_rows_and_other_sections():
    args = argparse.Namespace()
    run_log.note_run(args, detail={"other": {"SRR0": {"x": 1}}})
    run_log.note_rows(args, {"SRR1/G1": {"mapped_reads": 1}}, "extractions")
    run_log.note_rows(args, {"SRR2/G1": {"mapped_reads": 2}}, "extractions")
    run_log.note_rows(args, {"SRR1": {"gc_percent": 41.2}}, "analyses", "profile")
    run_log.note_rows(args, {}, "ignored")
    assert args._run_detail == {
        "other": {"SRR0": {"x": 1}},
        "extractions": {"SRR1/G1": {"mapped_reads": 1}, "SRR2/G1": {"mapped_reads": 2}},
        "analyses": {"profile": {"SRR1": {"gc_percent": 41.2}}},
    }
