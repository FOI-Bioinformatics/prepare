"""Tests for the per-run download summary: report CSV columns and ``download_run.json``.

The command tests run ``DownloadSraCommand.execute`` with ``download_sra`` replaced by a fake
that reports outcomes through ``on_result``, so no download tool, store or network is used.
"""

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch

import pytest

from prepare.cli.commands.sra import DownloadSraCommand
from prepare.data.sra.retry import DISK_FULL_NOT_ATTEMPTED
from prepare.data.sra.run_report import (
    DOWNLOAD_RUN_FILE,
    REASONS,
    REPORT_HEADER,
    RunOutcomes,
    failure_reason,
    report_rows,
    run_document,
    stats_from_outcomes,
    write_report_csv,
    write_run_document,
)
from prepare.data.sra.space import INSUFFICIENT_SPACE_PREFIX

NETWORK = "network: Download failed: connection timed out"


@pytest.mark.parametrize(
    "message, reason",
    [
        (NETWORK, "network"),
        ("not-found: failed to resolve accession SRR1: no data ( 404 )", "not-found"),
        ("disk-full: fasterq-dump: storage exhausted", "disk-full"),
        (DISK_FULL_NOT_ATTEMPTED, "disk-full"),
        (f"{INSUFFICIENT_SPACE_PREFIX} for SRR1: 12 GB needed on /data, 3 GB free", "insufficient-space"),
        ("locked: download of SRR1 is locked by pid 12 on host since 404: /x/.locks/SRR1.lock", "locked"),
        ("lock lost: download of SRR1 was taken over", "locked"),
        ("store error: dataset SRR1 is locked by pid 7 on node1: /store/locks/SRR1.lock", "locked"),
        ("interrupted", "interrupted"),
        ("Retry 1: interrupted", "interrupted"),
        ("Download failed: fasterq-dump exited with 3", "unknown"),
        ("", "unknown"),
        # A store copy kept but not linked is a settled outcome with no cause of its own.
        ("incomplete: 10 of 20 spots in store copy (connection reset)", "unknown"),
        ("partial in store; not found complete, refetch limit reached", "unknown"),
    ],
)
def test_failure_reason(message, reason):
    assert failure_reason(message) == reason
    assert reason in REASONS


def test_failure_reason_reads_through_a_retry_prefix():
    assert failure_reason(f"Retry 2: {INSUFFICIENT_SPACE_PREFIX} for SRR1") == "insufficient-space"
    assert failure_reason("Retry 1 error: [Errno 28] No space left on device") == "disk-full"


class TestRunOutcomes:
    def test_a_retried_network_failure_counts_two_attempts(self):
        outcomes = RunOutcomes()
        outcomes.observe("SRR1", False, NETWORK)
        outcomes.observe("SRR1", False, f"Retry 1: {NETWORK}")
        assert outcomes.attempts == {"SRR1": 2}
        assert outcomes.last["SRR1"] == (False, f"Retry 1: {NETWORK}")

    @pytest.mark.parametrize(
        "success, message",
        [
            (False, DISK_FULL_NOT_ATTEMPTED),
            (False, "interrupted"),
            (False, "Retry 1: interrupted"),
            (False, f"{INSUFFICIENT_SPACE_PREFIX} for SRR1"),
            (True, "already exists"),
            (True, "linked from store, 2 files, complete"),
            (True, "Retry 1: linked from store, 2 files"),
            (False, "locked: download of SRR1 is locked by pid 12: /x.lock"),
            (False, "partial in store; refetch limit reached"),
            (
                False,
                "incomplete: 2 refetches of SRR1 gained no reads (1 of 100 spots); NCBI's count may not be "
                "reachable; use --accept-partial, or --force to fetch again",
            ),
        ],
    )
    def test_outcomes_that_started_no_download_count_zero(self, success, message):
        outcomes = RunOutcomes()
        outcomes.observe("SRR1", success, message)
        assert outcomes.attempts == {"SRR1": 0}

    def test_the_store_refusal_of_an_exhausted_copy_counts_zero(self):
        # The message as the store precheck writes it, before any lock or download is taken.
        from types import SimpleNamespace

        from prepare.data.sra.store_handoff import _exhausted_message

        sidecar = SimpleNamespace(refetch={"unchanged": 2}, reads_per_mate=1, ncbi={"spots": 100})
        message = _exhausted_message("SRR1", sidecar, 100)
        outcomes = RunOutcomes()
        outcomes.observe("SRR1", False, message)
        outcomes.observe("SRR2", False, f"Retry 1: {message}")
        assert outcomes.attempts == {"SRR1": 0, "SRR2": 0}
        assert failure_reason(message) == "unknown"

    @pytest.mark.parametrize(
        "success, message",
        [
            (True, "Downloaded 2 files, complete"),
            (True, "Downloaded 2 files, complete; stored"),
            (False, "lock lost: taken over"),
            (False, "incomplete: 10 of 20 spots"),
            (False, "not-found: no data ( 404 )"),
        ],
    )
    def test_outcomes_that_started_a_download_count_one(self, success, message):
        outcomes = RunOutcomes()
        outcomes.observe("SRR1", success, message)
        assert outcomes.attempts == {"SRR1": 1}

    def test_wrap_observes_then_calls_the_callback(self):
        seen = []
        outcomes = RunOutcomes()
        wrapped = outcomes.wrap(lambda acc, ok, msg: seen.append((acc, ok, msg, dict(outcomes.attempts))))
        wrapped("SRR1", True, "Downloaded 1 files")
        assert seen == [("SRR1", True, "Downloaded 1 files", {"SRR1": 1})]
        RunOutcomes().wrap(None)("SRR2", False, NETWORK)

    def test_stats_from_outcomes(self):
        outcomes = RunOutcomes()
        outcomes.observe("SRR1", True, "Downloaded 1 files")
        outcomes.observe("SRR2", False, NETWORK)
        stats = stats_from_outcomes(outcomes)
        assert stats["successful"] == 1 and stats["failed"] == 1
        assert stats["failed_accessions"] == ["SRR2"]
        assert stats["results"] == {"SRR1": "Downloaded 1 files", "SRR2": NETWORK}
        assert stats["aborted"] is None


def test_report_rows_keep_the_first_four_columns_and_add_reason_and_attempts():
    stats = {
        "failed_accessions": ["SRR2", "SRR5"],
        "results": {"SRR1": "Downloaded 2 files", "SRR2": f"Retry 1: {NETWORK}", "SRR5": DISK_FULL_NOT_ATTEMPTED},
        "already_downloaded_accessions": ["SRR3"],
        "blacklisted_accessions": ["SRR4"],
        "skipped_accessions": ["SRR6"],
    }
    timings = {"SRR1": ("2026-10-01T10:00:00+00:00", 12.5)}
    rows = report_rows(stats, timings, {"SRR1": 1, "SRR2": 2, "SRR5": 0})
    assert rows == [
        ("SRR1", "downloaded", "Downloaded 2 files", "12.5", "", "1"),
        ("SRR2", "failed", f"Retry 1: {NETWORK}", "", "network", "2"),
        ("SRR3", "already_present", "", "", "", "0"),
        ("SRR4", "blacklisted", "", "", "", "0"),
        ("SRR5", "failed", DISK_FULL_NOT_ATTEMPTED, "", "disk-full", "0"),
        ("SRR6", "skipped", "--max-downloads", "", "", "0"),
    ]
    assert REPORT_HEADER[:4] == ("accession", "status", "message", "seconds")


def test_write_report_csv(tmp_path):
    path = tmp_path / "reports" / "report.csv"
    write_report_csv(path, [("SRR1", "failed", "a, b", "", "unknown", "1")])
    assert list(csv.reader(path.open())) == [list(REPORT_HEADER), ["SRR1", "failed", "a, b", "", "unknown", "1"]]


def test_run_document_and_its_file(tmp_path):
    outcomes = RunOutcomes()
    outcomes.observe("SRR1", False, NETWORK)
    outcomes.observe("SRR1", False, f"Retry 1: {NETWORK}")
    outcomes.observe("SRR2", True, "Downloaded 1 files")
    stats = dict(stats_from_outcomes(outcomes), already_downloaded_accessions=["SRR3"])
    started = datetime(2026, 10, 1, 10, 0, 0, tzinfo=timezone.utc)
    finished = datetime(2026, 10, 1, 10, 0, 30, tzinfo=timezone.utc)
    document = run_document(
        stats=stats,
        outcomes=outcomes,
        started=started,
        finished=finished,
        exit_code=4,
        aborted=None,
        settings={"workers": 2},
        paths={"fastq": str(tmp_path)},
    )
    assert document["started"] == "2026-10-01T10:00:00+00:00"
    assert document["seconds"] == 30.0
    assert document["exit_code"] == 4 and document["aborted"] is None
    assert document["host"]
    assert document["totals"] == {
        "accessions": 3,
        "downloaded": 1,
        "failed": 1,
        "already_present": 1,
        "blacklisted": 0,
        "skipped": 0,
        "attempts": 3,
    }
    assert document["failures_by_reason"] == {"network": 1}
    assert document["failed"] == [
        {"accession": "SRR1", "reason": "network", "attempts": 2, "message": f"Retry 1: {NETWORK}"}
    ]
    path = write_run_document(tmp_path / "fastq", document)
    assert path == tmp_path / "fastq" / DOWNLOAD_RUN_FILE
    assert json.loads(path.read_text()) == document


# ---------------------------------------------------------------- through the command


def _args(tmp_path, **overrides):
    values = dict(
        accessions_file=str(tmp_path / "acc.txt"),
        fastq_folder=str(tmp_path / "fastq"),
        max_downloads=None,
        num_threads=4,
        max_workers=2,
        dry_run=False,
        force=False,
        max_retries=1,
        temp_folder=None,
        blacklist=None,
        report_file=str(tmp_path / "report.csv"),
        registry=str(tmp_path / "prepare_registry.json"),
        data_root=None,
    )
    values.update(overrides)
    return argparse.Namespace(**values)


def _execute(fake, args):
    with (
        patch("prepare.cli.commands.sra.require_tools"),
        patch("prepare.cli.commands.sra.download_sra", side_effect=fake),
    ):
        return DownloadSraCommand().execute(args)


def _document(tmp_path):
    return json.loads((tmp_path / "fastq" / DOWNLOAD_RUN_FILE).read_text())


def _report(tmp_path):
    return list(csv.reader((tmp_path / "report.csv").open()))


def _preflight_failure(args):
    from prepare.core.exceptions import ConfigurationError

    fake = patch("prepare.cli.commands.sra.download_sra", side_effect=AssertionError("no download expected"))
    with patch("prepare.cli.commands.sra.require_tools", side_effect=ConfigurationError("fasterq-dump not found")):
        with fake:
            return DownloadSraCommand().execute(args)


def test_a_preflight_failure_creates_no_fastq_folder_and_writes_no_document(tmp_path):
    assert _preflight_failure(_args(tmp_path)) == 3
    assert not (tmp_path / "fastq").exists()
    # The report CSV the user asked for is still written, with no rows.
    assert _report(tmp_path) == [list(REPORT_HEADER)]


def test_a_preflight_failure_in_an_existing_fastq_folder_still_writes_the_document(tmp_path):
    (tmp_path / "fastq").mkdir()
    assert _preflight_failure(_args(tmp_path)) == 3
    document = _document(tmp_path)
    assert document["exit_code"] == 3
    assert document["totals"]["accessions"] == 0


def test_a_retried_network_failure_is_reported_with_two_attempts(tmp_path):
    def fake(**kwargs):
        kwargs["on_result"]("SRR1", False, NETWORK)
        kwargs["on_result"]("SRR1", False, f"Retry 1: {NETWORK}")
        return {
            "total": 1,
            "successful": 0,
            "failed": 1,
            "failed_accessions": ["SRR1"],
            "results": {"SRR1": f"Retry 1: {NETWORK}"},
            "aborted": None,
        }

    assert _execute(fake, _args(tmp_path)) == 4
    assert _report(tmp_path)[1] == ["SRR1", "failed", f"Retry 1: {NETWORK}", "", "network", "2"]
    document = _document(tmp_path)
    assert document["exit_code"] == 4
    assert document["failures_by_reason"] == {"network": 1}
    assert document["settings"]["retries"] == 1 and document["settings"]["workers"] == 2
    # download_sra writes the retry file when it returns with failures, so the document names it.
    assert document["paths"]["failed_accessions"] == str(tmp_path / "fastq" / "failed_accessions.txt")


def test_a_disk_full_abort_reports_not_attempted_rows_with_zero_attempts(tmp_path):
    def fake(**kwargs):
        kwargs["on_result"]("SRR1", False, "disk-full: No space left on device")
        kwargs["on_result"]("SRR2", False, DISK_FULL_NOT_ATTEMPTED)
        return {
            "total": 2,
            "successful": 0,
            "failed": 2,
            "failed_accessions": ["SRR1", "SRR2"],
            "results": {"SRR1": "disk-full: No space left on device", "SRR2": DISK_FULL_NOT_ATTEMPTED},
            "aborted": "disk-full",
        }

    assert _execute(fake, _args(tmp_path)) == 1
    rows = _report(tmp_path)
    assert [row[4:] for row in rows[1:]] == [["disk-full", "1"], ["disk-full", "0"]]
    document = _document(tmp_path)
    assert document["aborted"] == "disk-full"
    assert document["failures_by_reason"] == {"disk-full": 2}


def test_an_interrupt_mid_run_writes_both_files_and_exits_130(tmp_path):
    def fake(**kwargs):
        kwargs["timings"]["SRR1"] = ("2026-10-01T10:00:00+00:00", 4.0)
        kwargs["on_result"]("SRR1", True, "Downloaded 1 files")
        kwargs["on_result"]("SRR2", False, "interrupted")
        raise KeyboardInterrupt

    assert _execute(fake, _args(tmp_path)) == 130
    assert _report(tmp_path) == [
        list(REPORT_HEADER),
        ["SRR1", "downloaded", "Downloaded 1 files", "4.0", "", "1"],
        ["SRR2", "failed", "interrupted", "", "interrupted", "0"],
    ]
    document = _document(tmp_path)
    assert document["aborted"] == "interrupted"
    # An interrupted run writes no retry file; one left by an earlier run is not named.
    assert document["paths"]["failed_accessions"] is None
    assert document["exit_code"] == 130
    assert document["failed"] == [
        {"accession": "SRR2", "reason": "interrupted", "attempts": 0, "message": "interrupted"}
    ]


def test_a_failed_final_flush_still_writes_both_files(tmp_path, monkeypatch):
    import prepare.data.registry_batch as batch_mod
    from prepare.core.exceptions import DataAccessError

    def locked(path):
        raise DataAccessError(f"Registry is locked by another process: {path}.lock")

    monkeypatch.setattr(batch_mod, "registry_transaction", locked)
    monkeypatch.setattr("prepare.cli.commands.sra.FINAL_FLUSH_RETRY_SECONDS", 0)

    def fake(**kwargs):
        kwargs["on_result"]("SRR1", False, "Download failed: exit 3")
        return {
            "total": 2,
            "successful": 0,
            "failed": 1,
            "failed_accessions": ["SRR1"],
            "results": {"SRR1": "Download failed: exit 3"},
            "skipped_accessions": ["SRR9"],
            "aborted": None,
        }

    assert _execute(fake, _args(tmp_path)) == 1
    assert [row[:2] for row in _report(tmp_path)[1:]] == [["SRR1", "failed"], ["SRR9", "skipped"]]
    document = _document(tmp_path)
    assert document["exit_code"] == 1
    assert document["totals"]["skipped"] == 1


def test_a_dry_run_writes_neither_file(tmp_path):
    def fake(**kwargs):
        return {"total": 1, "to_download": 1, "already_downloaded": 0, "successful": 0, "failed": 0, "aborted": None}

    assert _execute(fake, _args(tmp_path, dry_run=True)) == 0
    assert not (tmp_path / "report.csv").exists()
    assert not (tmp_path / "fastq" / DOWNLOAD_RUN_FILE).exists()


def test_the_run_document_is_written_without_a_report_file(tmp_path):
    def fake(**kwargs):
        kwargs["on_result"]("SRR1", True, "Downloaded 1 files")
        return {"total": 1, "successful": 1, "failed": 0, "failed_accessions": [], "results": {"SRR1": "x"}}

    assert _execute(fake, _args(tmp_path, report_file=None)) == 0
    document = _document(tmp_path)
    assert document["exit_code"] == 0 and document["totals"]["downloaded"] == 1
    assert Path(document["paths"]["fastq"]) == tmp_path / "fastq"
    assert not (tmp_path / "report.csv").exists()


def test_a_summary_that_cannot_be_written_does_not_change_the_exit_code(tmp_path, caplog):
    def fake(**kwargs):
        return {"total": 0, "successful": 0, "failed": 0, "failed_accessions": [], "results": {}}

    with patch("prepare.cli.commands.sra.write_run_document", side_effect=OSError("read-only")):
        assert _execute(fake, _args(tmp_path)) == 0
    assert "Could not write the download run summary" in caplog.text
