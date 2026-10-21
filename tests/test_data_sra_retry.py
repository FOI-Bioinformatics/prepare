"""Tests for prepare.data.sra: processing download results, retries and failure handling."""

import logging
import pytest
from unittest.mock import Mock, call, patch
from prepare.data.sra.retry import (
    _handle_download_failure,
    _process_download_results,
    _retry_failed_downloads,
)
from prepare.data.sra.download import download_sra
from prepare.data.sra import accession as accession_mod
from prepare.data.sra import retry as retry_mod


class TestProcessDownloadResults:
    """Test _process_download_results function."""

    def test_process_download_results_all_success(self):
        """Test processing all successful results."""
        # Create mock futures
        future1 = Mock()
        future1.result.return_value = (True, "Success message 1")
        future2 = Mock()
        future2.result.return_value = (True, "Success message 2")

        futures_results = [("SRR123", (True, "Success message 1")), ("SRR456", (True, "Success message 2"))]
        accessions_to_download = ["SRR123", "SRR456"]
        download_results = {}
        failed_accessions = []

        _process_download_results(futures_results, accessions_to_download, download_results, failed_accessions)

        assert download_results == {"SRR123": "Success message 1", "SRR456": "Success message 2"}
        assert failed_accessions == []

    def test_process_download_results_some_failures(self):
        """Test processing mixed success/failure results."""
        future1 = Mock()
        future1.result.return_value = (True, "Success message")
        future2 = Mock()
        future2.result.return_value = (False, "Error message")

        futures_results = [("SRR123", (True, "Success message")), ("SRR456", (False, "Error message"))]
        accessions_to_download = ["SRR123", "SRR456"]
        download_results = {}
        failed_accessions = []

        _process_download_results(futures_results, accessions_to_download, download_results, failed_accessions)

        assert download_results == {"SRR123": "Success message", "SRR456": "Error message"}
        assert failed_accessions == ["SRR456"]

    def test_collect_result_records_the_raised_error_text(self):
        """A worker that raised is tallied as (False, "Error: <text>"), never as an unpack error."""
        from concurrent.futures import Future

        from prepare.data.sra.retry import _collect_result, _progress_reporter

        future: Future = Future()
        future.set_exception(OSError("Future failed"))
        futures_results = []
        seen = []

        with patch("prepare.data.sra.retry.logger"):
            outcome = _collect_result(
                future, "SRR456", futures_results, lambda *a: seen.append(a), _progress_reporter(1)
            )

        assert outcome == (False, "Future failed")
        assert futures_results == [("SRR456", (False, "Error: Future failed"))]
        assert seen == [("SRR456", False, "Future failed")]

        download_results = {}
        failed_accessions = []
        _process_download_results(futures_results, ["SRR456"], download_results, failed_accessions)
        assert download_results == {"SRR456": "Error: Future failed"}
        assert failed_accessions == ["SRR456"]

    def test_raised_network_error_is_classified_network(self, tmp_path, monkeypatch):
        """The exception text a worker raised is the recorded result, so its class survives."""
        from prepare.core.exceptions import NetworkError
        from prepare.data.sra.accession import classify_download_error
        from prepare.data.sra.retry import _download_with_retries

        def raising_downloader(accession, *args, **kwargs):
            raise NetworkError("Connection timed out while resolving SRR1")

        seen = []
        _successful, _failed, failed_accessions, download_results, _abort = _download_with_retries(
            ["SRR1"],
            tmp_path,
            1,
            1,
            False,
            None,
            0,
            on_result=lambda acc, ok, msg: seen.append((acc, ok, msg)),
            downloader=raising_downloader,
        )

        assert failed_accessions == ["SRR1"]
        assert download_results["SRR1"] == "Error: Connection timed out while resolving SRR1"
        assert classify_download_error(download_results["SRR1"]) == "network"
        assert seen == [("SRR1", False, "Connection timed out while resolving SRR1")]

    def test_raised_network_error_without_keyword_is_classified_network(self, tmp_path):
        """A NetworkError whose text names no network symptom still keeps its class."""
        from prepare.core.exceptions import NetworkError
        from prepare.data.sra.accession import classify_download_error
        from prepare.data.sra.retry import _download_with_retries

        def raising_downloader(accession, *args, **kwargs):
            raise NetworkError("boom")

        _successful, _failed, failed_accessions, download_results, _abort = _download_with_retries(
            ["SRR1"], tmp_path, 1, 1, False, None, 0, downloader=raising_downloader
        )

        assert failed_accessions == ["SRR1"]
        assert download_results["SRR1"].startswith("Error: network: ")
        assert classify_download_error(download_results["SRR1"]) == "network"

    def test_final_count_logged_once(self, caplog):
        """The final tally is logged exactly once, after the loop, not every fifth success."""
        results = [("SRR1", (True, "ok")), ("SRR2", (False, "boom")), ("SRR3", (True, "linked from store, 2 files"))]

        with caplog.at_level(logging.INFO, logger="prepare.data.sra"):
            ok, failed = _process_download_results(results, ["SRR1", "SRR2", "SRR3"], {}, [])

        assert (ok, failed) == (2, 1)
        assert "Downloaded 2 of 3 (1 failed)" in caplog.text
        assert caplog.text.count("Downloaded 2 of 3") == 1


class TestRetryFailedDownloads:
    """Test _retry_failed_downloads function."""

    def test_retry_failed_downloads_success(self, tmp_path):
        """Test successful retry of failed downloads."""
        failed_accessions = ["SRR123", "SRR456"]
        download_results = {}

        with patch("prepare.data.sra.accession.download_accession") as mock_download:
            mock_download.side_effect = [(True, "Retry success 1"), (True, "Retry success 2")]

            retried_successful, updated_failed, abort_reason = _retry_failed_downloads(
                failed_accessions,
                max_retries=2,
                fastq_path=tmp_path,
                num_threads=4,
                temp_folder=None,
                download_results=download_results,
            )

        assert download_results == {"SRR123": "Retry 1: Retry success 1", "SRR456": "Retry 1: Retry success 2"}
        assert len(updated_failed) == 0
        assert retried_successful == 2
        assert abort_reason is None

    def test_retry_passes_the_callers_force(self, tmp_path):
        """The retry pass hands the caller's force to the downloader; it defaults to False."""
        for kwargs, expected in (({"force": True}, True), ({}, False)):
            with patch("prepare.data.sra.accession.download_accession") as mock_download:
                mock_download.return_value = (True, "Retry success")

                _retry_failed_downloads(
                    ["SRR123"],
                    max_retries=1,
                    fastq_path=tmp_path,
                    num_threads=4,
                    temp_folder=None,
                    download_results={},
                    **kwargs,
                )

            mock_download.assert_called_once()
            assert mock_download.call_args.kwargs["force"] is expected

    def test_download_sra_with_force_keeps_it_on_the_retry(self, tmp_path):
        """A forced run that fails once then succeeds retries with force and never reports 'already exists'."""
        accessions_file = tmp_path / "accessions.txt"
        accessions_file.write_text("SRR123\n")
        calls = []

        def flaky(accession, *args, **kwargs):
            force = kwargs["force"] if "force" in kwargs else args[2]
            calls.append(force)
            if len(calls) == 1:
                return False, "Download failed: timeout"
            return (True, "already exists") if not force else (True, "downloaded")

        with patch("prepare.data.sra.accession.download_accession", side_effect=flaky):
            with patch("prepare.data.sra.retry._handle_download_failure"):
                with patch("prepare.data.sra.retry.time.sleep"):
                    download_sra(tmp_path / "out", accessions_file, force=True, max_retries=1, max_workers=1)

        assert calls == [True, True]

    def test_retry_failed_downloads_partial_success(self, tmp_path):
        """Test partial success in retry."""
        failed_accessions = ["SRR123", "SRR456"]
        download_results = {}

        with patch("prepare.data.sra.accession.download_accession") as mock_download:
            with patch("prepare.data.sra.retry.time.sleep") as mock_sleep:
                mock_download.side_effect = [(True, "Retry success"), (False, "Retry failed"), OSError("still failing")]

                retried_successful, updated_failed, abort_reason = _retry_failed_downloads(
                    failed_accessions,
                    max_retries=2,
                    fastq_path=tmp_path,
                    num_threads=4,
                    temp_folder=None,
                    download_results=download_results,
                )

        assert download_results == {"SRR123": "Retry 1: Retry success", "SRR456": "Retry 2 error: still failing"}
        assert updated_failed == ["SRR456"]
        assert retried_successful == 1
        assert abort_reason is None
        mock_sleep.assert_called_once_with(1)

    def test_split_not_found_reads_the_message_under_the_retry_prefix(self):
        results = {"SRR1": "Retry 1: incomplete: 3 of 9 spots", "SRR2": "Retry 2: Download failed: timeout"}
        retry_batch, not_found = retry_mod._split_not_found(["SRR1", "SRR2"], results)
        assert not_found == ["SRR1"]
        assert retry_batch == ["SRR2"]

    def test_result_body_strips_one_retry_prefix(self):
        assert accession_mod.result_body("Retry 3: already exists") == "already exists"
        assert accession_mod.result_body("Retry 2 error: boom") == "boom"
        assert accession_mod.result_body("plain") == "plain"
        assert accession_mod.result_body("") == ""

    def test_retry_failed_downloads_skips_not_found(self, tmp_path):
        """An accession already classified as not-found must not be retried again."""
        failed_accessions = ["SRR123", "SRR404"]
        download_results = {"SRR123": "network: Connection timed out", "SRR404": "not-found: 404 Not Found"}

        with patch("prepare.data.sra.accession.download_accession") as mock_download:
            mock_download.return_value = (True, "Retry success")

            retried_successful, updated_failed, abort_reason = _retry_failed_downloads(
                failed_accessions,
                max_retries=1,
                fastq_path=tmp_path,
                num_threads=4,
                temp_folder=None,
                download_results=download_results,
            )

        mock_download.assert_called_once_with(
            "SRR123",
            tmp_path,
            4,
            force=False,
            temp_folder=None,
            expected_spots=None,
            redownload_truncated=False,
            sra_cache=None,
            use_prefetch=True,
            keep_sra=False,
            compress=True,
            stop=None,
        )
        assert updated_failed == ["SRR404"]
        assert retried_successful == 1
        assert abort_reason is None

    def test_retry_failed_downloads_sleeps_exponentially_between_rounds(self, tmp_path):
        """Between retry rounds, sleep 2**attempt seconds so a flaky network gets a backoff."""
        failed_accessions = ["SRR123"]
        download_results = {}

        with patch("prepare.data.sra.accession.download_accession") as mock_download:
            with patch("prepare.data.sra.retry.time.sleep") as mock_sleep:
                mock_download.return_value = (False, "network: Connection timed out")

                _retry_failed_downloads(
                    failed_accessions,
                    max_retries=3,
                    fastq_path=tmp_path,
                    num_threads=4,
                    temp_folder=None,
                    download_results=download_results,
                )

        assert mock_sleep.call_args_list == [call(1), call(2)]

    def test_retry_pass_drops_interrupted_results(self, tmp_path):
        """A worker stopped by the token is not recorded as a failed attempt; the rest are not tried."""
        import threading

        stop = threading.Event()
        attempted = []
        seen = []

        def downloader(accession, *args, **kwargs):
            attempted.append(accession)
            stop.set()
            return False, "interrupted"

        download_results = {"SRR1": "network: timed out", "SRR2": "network: timed out"}
        retried, failed, abort = _retry_failed_downloads(
            ["SRR1", "SRR2"],
            1,
            tmp_path,
            1,
            None,
            download_results,
            on_result=lambda acc, ok, msg: seen.append((acc, ok, msg)),
            downloader=downloader,
            stop=stop,
        )

        assert attempted == ["SRR1"]
        assert seen == []
        assert download_results == {"SRR1": "network: timed out", "SRR2": "network: timed out"}
        assert sorted(failed) == ["SRR1", "SRR2"]
        assert retried == 0 and abort is None

    def test_retry_backoff_ends_when_the_token_is_set(self, tmp_path):
        """The wait between rounds ends as soon as the run is stopped."""
        import threading
        import time

        stop = threading.Event()

        def downloader(accession, *args, **kwargs):
            return False, "network: Connection timed out"

        threading.Timer(0.1, stop.set).start()
        started = time.monotonic()
        _retry_failed_downloads(["SRR1"], 3, tmp_path, 1, None, {}, downloader=downloader, stop=stop)
        assert time.monotonic() - started < 1.0

    def test_retry_failed_downloads_disk_full_aborts_without_raising(self, tmp_path):
        """A disk-full failure aborts the run but must not raise, so the caller keeps its state."""
        failed_accessions = ["SRR1", "SRR2", "SRR3"]
        download_results = {}
        notified = []

        def _on_result(accession, success, message):
            notified.append((accession, success, message))

        with patch("prepare.data.sra.accession.download_accession") as mock_download:
            with patch("prepare.data.sra.retry.time.sleep") as mock_sleep:
                mock_download.return_value = (False, "disk-full: No space left on device")

                retried_successful, updated_failed, abort_reason = _retry_failed_downloads(
                    failed_accessions,
                    max_retries=2,
                    fastq_path=tmp_path,
                    num_threads=4,
                    temp_folder=None,
                    download_results=download_results,
                    on_result=_on_result,
                )

        # Only the triggering accession is actually attempted; the rest of the round is
        # marked failed without ever calling download_accession.
        mock_download.assert_called_once()
        assert retried_successful == 0
        assert abort_reason == "disk-full"
        assert updated_failed == ["SRR1", "SRR2", "SRR3"]
        assert download_results["SRR1"] == "Retry 1: disk-full: No space left on device"
        assert download_results["SRR2"] == "disk-full: not attempted"
        assert download_results["SRR3"] == "disk-full: not attempted"
        # Every accession in the round is still notified, including the ones never attempted.
        assert notified == [
            ("SRR1", False, "Retry 1: disk-full: No space left on device"),
            ("SRR2", False, "disk-full: not attempted"),
            ("SRR3", False, "disk-full: not attempted"),
        ]
        # No further retry round runs, so no backoff sleep either.
        mock_sleep.assert_not_called()


class TestHandleDownloadFailure:
    """Test _handle_download_failure function."""

    def test_handle_download_failure_no_failures(self, tmp_path):
        """Test when no failures occurred."""
        fastq_path = tmp_path / "downloads"
        failed_accessions = []

        with patch("prepare.data.sra.retry.logger") as mock_logger:
            _handle_download_failure(fastq_path, failed_accessions)

        # No logging should occur when there are no failures
        mock_logger.info.assert_not_called()

    def test_clean_run_removes_a_stale_failed_file(self, tmp_path):
        fastq_path = tmp_path / "downloads"
        fastq_path.mkdir()
        stale = fastq_path / "failed_accessions.txt"
        stale.write_text("SRR999\n")

        _handle_download_failure(fastq_path, [])

        assert not stale.exists()

    def test_retry_hint_is_left_to_the_command(self, tmp_path):
        """The data layer writes the file; the command logs the one hint naming the user's folder."""
        fastq_path = tmp_path / "downloads"
        fastq_path.mkdir()

        with patch("prepare.data.sra.retry.logger") as mock_logger:
            _handle_download_failure(fastq_path, ["SRR123"])

        messages = [str(c.args[0]) for c in mock_logger.info.call_args_list]
        assert not any("To retry only failed" in m for m in messages)

    def test_handle_download_failure_with_failures(self, tmp_path):
        """Test when failures occurred."""
        fastq_path = tmp_path / "downloads"
        fastq_path.mkdir(parents=True, exist_ok=True)
        failed_accessions = ["SRR123", "SRR456"]

        with patch("prepare.data.sra.retry.logger") as mock_logger:
            _handle_download_failure(fastq_path, failed_accessions)

        mock_logger.info.assert_called()
        # Check that failed accessions file was created
        failed_file = fastq_path / "failed_accessions.txt"
        assert failed_file.exists()
        content = failed_file.read_text()
        assert "SRR123" in content
        assert "SRR456" in content


def test_unexpected_error_in_retry_downloader_propagates(tmp_path):
    """Kind (a): a downloader bug is no longer counted as a failed accession."""

    def buggy_downloader(*args, **kwargs):
        raise TypeError("bug")

    with pytest.raises(TypeError):
        _retry_failed_downloads(["SRR1"], 1, tmp_path, 1, None, {}, downloader=buggy_downloader)


def test_retry_downloader_os_error_still_marks_the_accession_failed(tmp_path):
    """Kind (a): an expected error still marks the item failed and continues."""

    def failing_downloader(*args, **kwargs):
        raise OSError("disk full")

    download_results = {}
    result = _retry_failed_downloads(["SRR1"], 1, tmp_path, 1, None, download_results, downloader=failing_downloader)
    assert "SRR1" in result[1]
    assert "disk full" in download_results["SRR1"]
