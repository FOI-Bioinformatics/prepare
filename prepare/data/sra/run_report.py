"""The per-run summary of a ``download_sra`` run: the report CSV and ``<fastq>/download_run.json``.

``RunOutcomes`` wraps the run's result callback and keeps, per accession, the last outcome and
the number of attempts that started a download. From those (or from ``download_sra``'s
statistics when the run returned them) ``report_rows`` builds the rows of the ``--report-file``
CSV and ``run_document`` the JSON document written to the FASTQ folder after every run that is
not a dry run, an interrupted one included.

Each failed accession gets one of ``REASONS``. ``failure_reason`` checks the messages the
download layer writes itself (an interruption, a free-space refusal, a lock) before falling back
to ``classify_download_error``, which is left unchanged. A store copy kept but not linked
(``SETTLED_PREFIXES``) is a settled outcome rather than an error with a cause of its own, so it
is ``"unknown"`` whatever words its message contains.
"""

import argparse
import csv
import json
import socket
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

from prepare.core.constants import FAILED_ACCESSIONS_FILE
from prepare.core.settings import setting_for
from prepare.data.file_io import open_atomic, write_text_atomic
from prepare.data.sra import accession as accession_mod
from prepare.data.sra.download import default_max_workers
from prepare.data.sra.retry import DISK_FULL_NOT_ATTEMPTED
from prepare.data.sra.space import INSUFFICIENT_SPACE_PREFIX
from prepare.data.sra.store_handoff import (
    SETTLED_PREFIXES,
    STORE_EXHAUSTED_SUFFIX,
    STORE_INCOMPLETE_PREFIX,
    STORE_LINKED_PREFIX,
    STORE_PARTIAL_PREFIX,
    STORE_PROJECT_FOLDER_PREFIX,
)

REASONS = ("network", "not-found", "disk-full", "insufficient-space", "locked", "interrupted", "unknown")
REPORT_HEADER = ("accession", "status", "message", "seconds", "reason", "attempts")
DOWNLOAD_RUN_FILE = "download_run.json"

# Version of the download_run.json layout, raised when a key changes meaning or is removed.
DOCUMENT_VERSION = 1

ResultCallback = Callable[[str, bool, str], None]
# Accession -> (start as ISO 8601 UTC, seconds) of its last download attempt, filled by download_sra.
Timings = Dict[str, Tuple[str, float]]

_INTERRUPTED = "interrupted"
_LOCK_LOST = "lock lost: "


def failure_reason(message: str) -> str:
    """The ``REASONS`` value for a failed accession's last result message.

    An interruption, a free-space refusal and a lock (one not taken, or lost while the download
    ran) are recognised first, then a settled store outcome (``"unknown"``); anything else is
    ``classify_download_error``'s class of the message.
    """
    body = accession_mod.result_body(message)
    if body == _INTERRUPTED:
        return "interrupted"
    if body.startswith(INSUFFICIENT_SPACE_PREFIX):
        return "insufficient-space"
    if accession_mod._LOCK_MESSAGE_RE.search(body):
        return "locked"
    if body.startswith(SETTLED_PREFIXES):
        return "unknown"
    return accession_mod.classify_download_error(body)


def started_download(success: bool, message: str) -> bool:
    """Whether a result came from an attempt that started a download.

    Not one: a disk-full "not attempted" mark, an interruption, a free-space refusal, a wait for
    a lock that gave up, the store's refusal of a partial copy, of an incomplete copy whose
    refetches are used up or of a project folder that is not a link, files found already in
    place, and a link to a dataset the store already held. A lock lost during the download
    did start one, and so did every other ``incomplete:`` outcome (it follows a fetch).
    """
    body = accession_mod.result_body(message)
    if success:
        return body != accession_mod.ALREADY_EXISTS and not body.startswith(STORE_LINKED_PREFIX)
    if body == DISK_FULL_NOT_ATTEMPTED or body.startswith((STORE_PARTIAL_PREFIX, STORE_PROJECT_FOLDER_PREFIX)):
        return False
    if body.startswith(STORE_INCOMPLETE_PREFIX) and body.endswith(STORE_EXHAUSTED_SUFFIX):
        return False
    reason = failure_reason(body)
    if reason == "locked":
        return body.startswith(_LOCK_LOST)
    return reason not in ("interrupted", "insufficient-space")


class RunOutcomes:
    """The outcomes one download run reported: last result and attempt count per accession.

    ``timings`` is the dict handed to ``download_sra`` for each attempt's start and seconds;
    ``stats`` is ``download_sra``'s statistics once it returned (None after an interrupt);
    ``reached_download`` is set when the run got as far as calling ``download_sra`` (a run that
    stopped in its preflight checks leaves it False).
    The callback is only called on the run's main thread, so no lock is needed.
    """

    def __init__(self) -> None:
        """Start with no outcome observed."""
        self.attempts: Dict[str, int] = {}
        self.last: Dict[str, Tuple[bool, str]] = {}
        self.timings: Timings = {}
        self.stats: Optional[Dict[str, Any]] = None
        self.reached_download = False

    def observe(self, accession: str, success: bool, message: str) -> None:
        """Record one result for ``accession``, counting it as an attempt when it started a download."""
        self.last[accession] = (bool(success), message)
        self.attempts[accession] = self.attempts.get(accession, 0) + int(started_download(success, message))

    def wrap(self, on_result: Optional[ResultCallback]) -> ResultCallback:
        """``on_result`` (which may be None) preceded by ``observe``."""

        def _observed(accession: str, success: bool, message: str) -> None:
            self.observe(accession, success, message)
            if on_result is not None:
                on_result(accession, success, message)

        return _observed


def stats_from_outcomes(outcomes: RunOutcomes) -> Dict[str, Any]:
    """Statistics in ``download_sra``'s shape built from the observed results alone.

    Used when the run returned none (an interrupt): only accessions that reported a result are
    included, and none is listed as already present, blacklisted or skipped.
    """
    failed = [acc for acc, (success, _) in outcomes.last.items() if not success]
    return {
        "total": len(outcomes.last),
        "successful": len(outcomes.last) - len(failed),
        "failed": len(failed),
        "failed_accessions": failed,
        "results": {acc: message for acc, (_, message) in outcomes.last.items()},
        "already_downloaded_accessions": [],
        "blacklisted_accessions": [],
        "skipped_accessions": [],
        "aborted": None,
    }


def attempt_timing(accession: str, success: bool, message: str, timings: Timings) -> Optional[Tuple[str, float]]:
    """``(started, seconds)`` of this run's download of ``accession``, or None when it ran none.

    A dataset linked from the store was not downloaded, so it has no time even when the worker
    that linked it was timed.
    """
    if success and accession_mod.result_body(message).startswith(STORE_LINKED_PREFIX):
        return None
    return timings.get(accession)


def report_rows(stats: Dict[str, Any], timings: Optional[Timings], attempts: Dict[str, int]) -> List[Tuple[str, ...]]:
    """One sorted row per accession: ``REPORT_HEADER``'s columns, all strings.

    ``seconds`` is filled for a download this run timed; ``reason`` for a failed accession only.
    """
    failed = set(stats.get("failed_accessions") or [])
    rows = []
    for accession, message in (stats.get("results") or {}).items():
        success = accession not in failed
        timing = attempt_timing(accession, success, message, timings or {})
        seconds = str(timing[1]) if timing is not None else ""
        reason = "" if success else failure_reason(message)
        status = "downloaded" if success else "failed"
        rows.append((accession, status, message, seconds, reason, str(attempts.get(accession, 0))))
    for key, status, message in (
        ("already_downloaded_accessions", "already_present", ""),
        ("blacklisted_accessions", "blacklisted", ""),
        ("skipped_accessions", "skipped", "--max-downloads"),
    ):
        rows.extend((acc, status, message, "", "", str(attempts.get(acc, 0))) for acc in stats.get(key) or [])
    return sorted(rows)


def write_report_csv(path: Union[str, Path], rows: Sequence[Sequence[str]]) -> None:
    """Write ``REPORT_HEADER`` and ``rows`` to ``path`` atomically, creating its folder."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with open_atomic(target, newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(REPORT_HEADER)
        writer.writerows(rows)


def run_settings(args: argparse.Namespace) -> Dict[str, Any]:
    """The settings of a ``download_sra`` run that a reader of its document needs to compare runs."""
    workers = args.max_workers if args.max_workers is not None else default_max_workers(args.num_threads)
    return {
        "workers": workers,
        "threads": args.num_threads,
        "retries": args.max_retries,
        "max_downloads": args.max_downloads,
        "force": bool(args.force),
        "verify": bool(getattr(args, "verify_downloads", True)),
        "redownload": bool(getattr(args, "redownload_truncated", False)),
        "prefetch": bool(getattr(args, "use_prefetch", True)),
        "compress": bool(getattr(args, "compress", True)),
        "min_free_gb": setting_for(args, "min_free_gb"),
        "link_mode": getattr(args, "link_mode", "auto"),
    }


def run_paths(args: argparse.Namespace, returned: Optional[Dict[str, Any]] = None) -> Dict[str, Optional[str]]:
    """The files a ``download_sra`` run read or wrote.

    ``returned`` is the statistics ``download_sra`` returned (None after an interrupt).
    ``failed_accessions`` names the retry file only when this run wrote it, which ``download_sra``
    does when it returns with failures; a file left by an earlier run is not named.
    """
    fastq = Path(args.fastq_folder)
    stats = returned or {}
    wrote_failed = int(stats.get("failed") or 0) > 0 and bool(stats.get("failed_accessions"))
    return {
        "fastq": str(fastq),
        "accessions_file": args.accessions_file,
        "report": args.report_file,
        "failed_accessions": str(fastq / FAILED_ACCESSIONS_FILE) if wrote_failed else None,
        "registry": args.registry,
        "data_root": getattr(args, "data_root", None),
    }


def run_document(
    *,
    stats: Dict[str, Any],
    outcomes: RunOutcomes,
    started: datetime,
    finished: datetime,
    exit_code: int,
    aborted: Optional[str],
    settings: Dict[str, Any],
    paths: Dict[str, Optional[str]],
) -> Dict[str, Any]:
    """The ``download_run.json`` document of one run.

    ``totals.accessions`` counts the accessions with a report row (after an interrupt, only the
    ones that reported a result); ``totals.attempts`` sums the attempts that started a download.
    ``failed`` holds one entry per failed accession, sorted, and ``failures_by_reason`` their
    count per ``REASONS`` value.
    """
    rows = report_rows(stats, outcomes.timings, outcomes.attempts)
    failed = [
        {"accession": acc, "reason": reason, "attempts": int(attempts), "message": message}
        for acc, status, message, _, reason, attempts in rows
        if status == "failed"
    ]
    statuses = Counter(row[1] for row in rows)
    return {
        "version": DOCUMENT_VERSION,
        "command": "download_sra",
        "started": started.isoformat(timespec="seconds"),
        "finished": finished.isoformat(timespec="seconds"),
        "seconds": round((finished - started).total_seconds(), 3),
        "exit_code": int(exit_code),
        "aborted": aborted,
        "host": socket.gethostname(),
        "totals": {
            "accessions": len(rows),
            "downloaded": statuses["downloaded"],
            "failed": statuses["failed"],
            "already_present": statuses["already_present"],
            "blacklisted": statuses["blacklisted"],
            "skipped": statuses["skipped"],
            "attempts": sum(outcomes.attempts.values()),
        },
        "failures_by_reason": dict(sorted(Counter(entry["reason"] for entry in failed).items())),
        "failed": failed,
        "settings": settings,
        "paths": paths,
    }


def run_log_entries(
    document: Dict[str, Any], rows: Sequence[Sequence[str]] = ()
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """The run-log summary and detail of one ``run_document``.

    The summary holds the totals, the failure count per reason and the abort cause. The detail
    holds rows only: one per failed accession (reason, attempts, message) under ``"failed"``, and
    one per accession downloaded in this run (status, seconds, attempts) under ``"downloaded"``,
    taken from ``rows`` (``report_rows``); a section without a row is left out.
    """
    summary = {
        **document["totals"],
        "failures_by_reason": document["failures_by_reason"],
        "aborted": document["aborted"],
    }
    detail: Dict[str, Any] = {}
    failed = {entry["accession"]: {k: v for k, v in entry.items() if k != "accession"} for entry in document["failed"]}
    if failed:
        detail["failed"] = failed
    downloaded = {
        acc: {"status": status, "seconds": float(seconds) if seconds else None, "attempts": int(attempts)}
        for acc, status, _, seconds, _, attempts in rows
        if status == "downloaded"
    }
    if downloaded:
        detail["downloaded"] = downloaded
    return summary, detail


def write_run_document(fastq_dir: Union[str, Path], document: Dict[str, Any]) -> Path:
    """Write ``document`` to ``<fastq_dir>/download_run.json`` atomically and return that path."""
    target = Path(fastq_dir) / DOWNLOAD_RUN_FILE
    target.parent.mkdir(parents=True, exist_ok=True)
    write_text_atomic(target, json.dumps(document, indent=2, sort_keys=False) + "\n")
    return target
