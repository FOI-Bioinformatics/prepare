"""Tests for ProgressReporter and the progress summaries of download_sra and download_metadata.

A long run logs one summary line every ``every`` items (or every ``min_interval`` seconds,
whichever comes first) and one at the end; the line for each item moves to DEBUG, unless
``every`` is 0, which keeps one INFO line per item and no summaries.
"""

import logging
import threading
from unittest.mock import MagicMock, patch

import pytest

from prepare.data import metadata as metadata_mod
from prepare.data.sra import download as download_mod
from prepare.utils.progress import ProgressReporter, format_duration

LOGGER = "prepare.test_progress"


class FakeClock:
    """A monotonic clock the test moves by hand."""

    def __init__(self) -> None:
        self.now = 1000.0

    def __call__(self) -> float:
        return self.now

    def advance(self, seconds: float) -> None:
        self.now += seconds


def _lines(caplog, level=logging.INFO):
    return [r.getMessage() for r in caplog.records if r.name == LOGGER and r.levelno == level]


def _reporter(total, every, clock, **kwargs):
    return ProgressReporter("download_sra", total, every, logger=logging.getLogger(LOGGER), clock=clock, **kwargs)


# --- the reporter -----------------------------------------------------------------


def test_lines_at_n_and_2n_and_finish(caplog):
    clock = FakeClock()
    reporter = _reporter(25, 10, clock, min_interval=1e9)
    with caplog.at_level(logging.INFO, logger=LOGGER):
        for index in range(25):
            clock.advance(60.0)
            reporter.update(ok=index != 3)
        reporter.finish()
    lines = _lines(caplog)
    assert lines == [
        "download_sra: 10/25 done (9 ok, 1 failed), 1.0/min, about 15 min left",
        "download_sra: 20/25 done (19 ok, 1 failed), 1.0/min, about 5 min left",
        "download_sra: finished 25/25 (24 ok, 1 failed) in 25 min",
    ]


def test_the_line_format_in_hours():
    clock = FakeClock()
    reporter = _reporter(2000, 150, clock, min_interval=1e9)
    clock.advance(150 / 3.1 * 60)
    with patch.object(reporter, "_log") as log:
        reporter.update(ok=True, n=148)
        reporter.update(ok=False, n=2)
    log.assert_called_once_with("download_sra: 150/2000 done (148 ok, 2 failed), 3.1/min, about 9 h 57 min left")


def test_a_batch_update_that_crosses_a_mark_logs_once(caplog):
    clock = FakeClock()
    reporter = _reporter(1000, 50, clock)
    with caplog.at_level(logging.INFO, logger=LOGGER):
        clock.advance(10)
        reporter.update(ok=True, n=200)
        reporter.update(ok=True, n=10)
    assert len(_lines(caplog)) == 1
    assert _lines(caplog)[0].startswith("download_sra: 200/1000 done (200 ok, 0 failed)")


def test_min_interval_triggers_a_line_between_marks(caplog):
    clock = FakeClock()
    reporter = _reporter(1000, 100, clock, min_interval=300.0)
    with caplog.at_level(logging.INFO, logger=LOGGER):
        reporter.update(ok=True)
        clock.advance(299)
        reporter.update(ok=True)
        assert _lines(caplog) == []
        clock.advance(2)
        reporter.update(ok=True)
        assert len(_lines(caplog)) == 1 and _lines(caplog)[0].startswith("download_sra: 3/1000 done")
        clock.advance(10)
        reporter.update(ok=True)
        assert len(_lines(caplog)) == 1  # the interval counts from the last line


def test_every_zero_logs_no_summaries_and_items_at_info(caplog):
    clock = FakeClock()
    reporter = _reporter(500, 0, clock)
    assert reporter.item_level == logging.INFO
    assert _reporter(500, 10, clock).item_level == logging.DEBUG
    with caplog.at_level(logging.INFO, logger=LOGGER):
        for _ in range(500):
            clock.advance(1000)
            reporter.update(ok=True)
    assert _lines(caplog) == []
    with caplog.at_level(logging.INFO, logger=LOGGER):
        reporter.finish()
    assert _lines(caplog) == ["download_sra: finished 500/500 (500 ok, 0 failed) in 138 h 53 min"]


def test_finish_logs_once_and_not_at_all_with_nothing_done(caplog):
    clock = FakeClock()
    with caplog.at_level(logging.INFO, logger=LOGGER):
        _reporter(10, 5, clock).finish()
        reporter = _reporter(10, 5, clock)
        reporter.update(ok=True)
        reporter.finish()
        reporter.finish()
    assert len(_lines(caplog)) == 1


def test_no_rate_before_any_time_has_passed(caplog):
    clock = FakeClock()
    reporter = _reporter(10, 1, clock)
    with caplog.at_level(logging.INFO, logger=LOGGER):
        reporter.update(ok=True)
    assert _lines(caplog) == ["download_sra: 1/10 done (1 ok, 0 failed)"]


def test_a_batch_of_successes_and_failures_is_one_update(caplog):
    clock = FakeClock()
    reporter = _reporter(1000, 100, clock)
    with caplog.at_level(logging.INFO, logger=LOGGER):
        clock.advance(60)
        reporter.update_counts(180, 20)
    assert [line.split(" (")[0] for line in _lines(caplog)] == ["download_sra: 200/1000 done"]
    assert (reporter.ok, reporter.failed) == (180, 20)


def test_a_negative_every_is_refused():
    with pytest.raises(ValueError):
        ProgressReporter("x", 10, -1)


@pytest.mark.parametrize(
    "seconds, text",
    [(0, "less than 1 min"), (59, "less than 1 min"), (60, "1 min"), (3569, "59 min"), (3599, "1 h 0 min")],
)
def test_format_duration(seconds, text):
    assert format_duration(seconds) == text


def test_eight_threads_times_a_thousand_updates(caplog):
    clock = FakeClock()
    reporter = _reporter(8000, 1000, clock)

    def work(ok: bool) -> None:
        for _ in range(1000):
            reporter.update(ok=ok)

    threads = [threading.Thread(target=work, args=(index % 2 == 0,)) for index in range(8)]
    with caplog.at_level(logging.INFO, logger=LOGGER):
        clock.advance(60)
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        reporter.finish()
    assert (reporter.done, reporter.ok, reporter.failed) == (8000, 4000, 4000)
    lines = _lines(caplog)
    assert len(lines) == 8  # 1000, 2000, ... 7000, then the finish line
    assert lines[-1].startswith("download_sra: finished 8000/8000 (4000 ok, 4000 failed)")
    assert [line.split(" ")[1] for line in lines[:-1]] == [f"{k * 1000}/8000" for k in range(1, 8)]


# --- download_sra -----------------------------------------------------------------


def _fake_project_download(accession, *args, **kwargs):
    ok = not accession.endswith("7")
    return ok, "Downloaded 2 files" if ok else "fasterq-dump failed: exit status 3"


@pytest.mark.parametrize("every", [25, 0])
def test_download_sra_logs_summaries_not_one_line_per_accession(tmp_path, caplog, monkeypatch, every):
    monkeypatch.setenv("PREPARE_PROGRESS_EVERY", str(every))
    accessions = [f"SRR{100000 + index}" for index in range(200)]
    accessions_file = tmp_path / "accessions.txt"
    accessions_file.write_text("".join(f"{acc}\n" for acc in accessions))
    with (
        patch("prepare.data.sra.accession._project_download", side_effect=_fake_project_download),
        caplog.at_level(logging.DEBUG),
    ):
        stats = download_mod.download_sra(tmp_path / "fastq", accessions_file, max_workers=4, max_retries=0)
    assert stats["successful"] == 180 and stats["failed"] == 20

    info = [r for r in caplog.records if r.levelno == logging.INFO]
    per_item = [r for r in caplog.records if r.getMessage().endswith(": Downloaded 2 files")]
    assert len(per_item) == 180
    summaries = [r for r in info if r.getMessage().startswith("download_sra: ")]
    if every:
        fixed = 15
        assert len(info) <= len(accessions) // every + fixed
        assert all(r.levelno == logging.DEBUG for r in per_item)
        assert len(summaries) == len(accessions) // every
        assert summaries[-1].getMessage().startswith("download_sra: finished 200/200 (180 ok, 20 failed)")
    else:
        assert all(r.levelno == logging.INFO for r in per_item)
        assert [r.getMessage().split(" (")[0] for r in summaries] == ["download_sra: finished 200/200"]
    # A failure stays a warning whatever the progress interval.
    assert sum(1 for r in caplog.records if r.levelno == logging.WARNING and "SRR100007" in r.getMessage()) == 1


# --- download_metadata --------------------------------------------------------------


def _fake_batch(batch, metadata_path, email, api_key):
    successes = {acc: metadata_path / f"{acc}_metadata.xml" for acc in batch if not acc.endswith("9")}
    failures = {acc: "not in the NCBI response" for acc in batch if acc.endswith("9")}
    return successes, failures


@pytest.mark.parametrize("every", [100, 0])
def test_download_metadata_logs_summaries_per_batch(tmp_path, caplog, monkeypatch, every):
    monkeypatch.setenv("PREPARE_PROGRESS_EVERY", str(every))
    accessions = [f"SRR{200000 + index}" for index in range(1000)]
    with (
        patch("prepare.data.metadata._download_batch_metadata", side_effect=_fake_batch),
        caplog.at_level(logging.DEBUG, logger="prepare.data.metadata"),
    ):
        result = metadata_mod._download_accessions_metadata(accessions, tmp_path, "a@b.c", 1000, batch_size=200)
    assert len(result) == 900
    summaries = [r.getMessage() for r in caplog.records if r.getMessage().startswith("download_metadata: ")]
    assert all(r.levelno == logging.INFO for r in caplog.records if r.getMessage() in summaries)
    if every:
        assert [line.split(" ")[1] for line in summaries] == [
            "200/1000",
            "400/1000",
            "600/1000",
            "800/1000",
            "finished",
        ]
    else:
        assert [line.split(" (")[0] for line in summaries] == ["download_metadata: finished 1000/1000"]
    assert summaries[-1].startswith("download_metadata: finished 1000/1000 (900 ok, 100 failed)")
    # Each failed accession is still named at ERROR.
    assert sum(1 for r in caplog.records if r.levelno == logging.ERROR) == 100


@pytest.mark.parametrize("every, level", [(50, logging.DEBUG), (0, logging.INFO)])
def test_download_metadata_request_lines_follow_the_progress_setting(tmp_path, caplog, monkeypatch, every, level):
    monkeypatch.setenv("PREPARE_PROGRESS_EVERY", str(every))

    class Handle:
        def read(self):
            return b"<EXPERIMENT_PACKAGE_SET/>"

        def close(self):
            pass

    with (
        patch("prepare.data.metadata.Entrez.efetch", return_value=Handle()),
        patch("prepare.data.metadata._pace_requests"),
        caplog.at_level(logging.DEBUG, logger="prepare.data.metadata"),
    ):
        metadata_mod._download_batch_metadata(["SRR1", "SRR2"], tmp_path, "a@b.c", None)
        metadata_mod._download_single_metadata("SRR3", tmp_path, "a@b.c", None)
    requests = [r for r in caplog.records if r.getMessage().startswith("Downloading metadata for")]
    assert len(requests) == 2 and all(r.levelno == level for r in requests)


# --- fix wave: settings read before work starts, closing line after the retry pass ------------


def test_a_bad_progress_setting_stops_downloads_before_any_starts(tmp_path, monkeypatch):
    from prepare.core.exceptions import ConfigurationError
    from prepare.data.sra import retry as retry_mod

    monkeypatch.setenv("PREPARE_PROGRESS_EVERY", "abc")
    worker = MagicMock(return_value=(True, "ok"))
    with pytest.raises(ConfigurationError, match="PREPARE_PROGRESS_EVERY"):
        retry_mod._download_with_retries(["SRR1", "SRR2"], tmp_path, 1, 2, False, None, 1, downloader=worker)
    worker.assert_not_called()


def test_a_bad_progress_setting_stops_metadata_before_any_request(tmp_path, monkeypatch):
    from prepare.core.exceptions import ConfigurationError

    monkeypatch.setenv("PREPARE_PROGRESS_EVERY", "abc")
    with patch("prepare.data.metadata._download_batch_metadata") as batch:
        with pytest.raises(ConfigurationError):
            metadata_mod._download_accessions_metadata(["SRR1"], tmp_path, "a@b.c", 1)
    batch.assert_not_called()


def test_download_closing_line_follows_the_retry_pass(tmp_path, caplog):
    from prepare.data.sra import retry as retry_mod

    attempts = {}

    def flaky(accession, *args, **kwargs):
        attempts[accession] = attempts.get(accession, 0) + 1
        if accession == "SRR2" and attempts[accession] == 1:
            return False, "network: connection reset"
        return True, "Downloaded 2 files"

    with caplog.at_level(logging.INFO):
        successful, failed_count, *_ = retry_mod._download_with_retries(
            ["SRR1", "SRR2"], tmp_path, 1, 2, False, None, 1, downloader=flaky
        )
    assert (successful, failed_count) == (2, 0)
    messages = [r.getMessage() for r in caplog.records]
    closing = [i for i, m in enumerate(messages) if m.startswith("download_sra: finished")]
    assert len(closing) == 1
    assert messages[closing[0]].startswith("download_sra: finished 2/2 (2 ok, 0 failed)")
    retry_line = next(i for i, m in enumerate(messages) if m.startswith("Retry attempt"))
    assert retry_line < closing[0]
