"""
The per-project run log: one record per run of a command that keeps one.

Layout, under ``<project>/.prepare/runs/`` (the project is the folder holding the registry file):

- ``runs.jsonl``: one JSON object per line, one line per run, in the order the runs finished. Lines
  are appended under ``runs.jsonl.lock``; a line cut short by an interrupted append is skipped on
  reading, with one warning.
- ``<run_id>.json``: the detail of one run (larger tables a command passes through ``note_run``),
  written atomically. Only the last ``DETAILS_KEPT_PER_COMMAND`` detail files of each command are
  kept; a pruned run keeps its line in ``runs.jsonl`` with ``"detail": null`` (and any key a later
  schema added to the line). A detail file written for an append that then failed names no line
  and so is never pruned; it is left in the folder.

``main`` records a run after the command returns, when the command's ``records_run`` is true, the
``run_log`` setting is on and a registry file exists; a failure to record is logged as a warning and
never changes the command's exit code. The registry itself holds no list of runs.
"""

import argparse
import json
import logging
import os
import secrets
import socket
from dataclasses import asdict, dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

from prepare import __version__
from prepare.core import settings
from prepare.core.exceptions import ValidationError
from prepare.data.file_io import write_text_atomic
from prepare.data.registry import REGISTRY_FILENAME
from prepare.utils.lockfile import LockPolicy, held_lock

logger = logging.getLogger(__name__)

RUN_LOG_SCHEMA = 1
RUNS_FILE = "runs.jsonl"
DETAILS_KEPT_PER_COMMAND = 10
RUN_LOG_LOCK_POLICY = LockPolicy("Run log", stale_seconds=60, wait_seconds=10)
MASK = "***"

PathLike = Union[str, Path]


@dataclass
class RunRecord:
    """One run of one command, as stored on one line of ``runs.jsonl``.

    ``started`` and ``finished`` are UTC times (``2026-10-01T12:05:01Z``); ``seconds`` is the measured
    run time; ``argv`` and ``args`` have secret values replaced by ``***``; ``detail`` is the name of
    the run's detail file inside the runs folder, or None when there is none (or it was pruned).
    """

    run_id: str
    command: str
    started: str
    finished: str
    seconds: float
    exit_code: int
    argv: List[str]
    args: Dict[str, Any]
    version: str
    host: str
    pid: int
    summary: Dict[str, Any] = field(default_factory=dict)
    detail: Optional[str] = None
    schema: int = RUN_LOG_SCHEMA

    def to_dict(self) -> Dict[str, Any]:
        """The record as a JSON-ready dict."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RunRecord":
        """A record from ``to_dict`` output; unknown keys are ignored, a wrong shape raises ValueError."""
        if not isinstance(data, dict) or not isinstance(data.get("run_id"), str) or not data.get("run_id"):
            raise ValueError("not a run record")
        if not isinstance(data.get("command"), str):
            raise ValueError("run record without a command")
        detail = data.get("detail")
        return cls(
            run_id=data["run_id"],
            command=data["command"],
            started=str(data.get("started", "")),
            finished=str(data.get("finished", "")),
            seconds=float(data.get("seconds", 0.0)),
            exit_code=int(data.get("exit_code", 1)),
            argv=[str(item) for item in data.get("argv") or []],
            args=dict(data.get("args") or {}),
            version=str(data.get("version", "")),
            host=str(data.get("host", "")),
            pid=int(data.get("pid", 0)),
            summary=dict(data.get("summary") or {}),
            detail=detail if isinstance(detail, str) and detail else None,
            schema=int(data.get("schema", RUN_LOG_SCHEMA)),
        )


def runs_dir(project: PathLike) -> Path:
    """The run log folder of ``project``: ``<project>/.prepare/runs``."""
    return Path(project) / ".prepare" / "runs"


def project_for(registry: Optional[PathLike] = None, start: PathLike = ".") -> Optional[Path]:
    """The project folder a run belongs to: the parent of its registry file, or None without one.

    ``registry`` is the command's ``--registry`` value; without it the registry is looked for in
    ``start`` and its parents, as ``registry_path`` does (without logging where it was found).
    """
    if registry:
        path = Path(registry)
        return path.resolve().parent if path.is_file() else None
    current = Path(start).resolve()
    for folder in (current, *current.parents):
        if (folder / REGISTRY_FILENAME).is_file():
            return folder
    return None


def _utc_text(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def new_run_id(command: str, started: datetime) -> str:
    """A run identifier: UTC start time, command and four random hex digits (``20261001T120501Z-sra_profile-3fa2``)."""
    stamp = started.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{stamp}-{command}-{secrets.token_hex(2)}"


def _json_safe(value: Any) -> Any:
    """``value`` with paths as strings, sequences and sets as lists and anything unknown as its text."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted((_json_safe(item) for item in value), key=str)
    return str(value)


def json_safe_args(args: argparse.Namespace) -> Dict[str, Any]:
    """The parsed arguments as JSON-ready values, without ``_*`` attributes or ``func``.

    The value of a secret setting's flag (the NCBI API key) is replaced by ``***``.
    """
    try:
        values = dict(vars(args))
    except TypeError:
        return {}
    secret_dests = {spec.cli_dest for spec in settings.SETTINGS.values() if spec.secret and spec.cli_dest}
    safe: Dict[str, Any] = {}
    for key in sorted(values):
        if key.startswith("_") or key == "func":
            continue
        value = values[key]
        safe[key] = MASK if key in secret_dests and value else _json_safe(value)
    return safe


def note_run(
    args: argparse.Namespace, summary: Optional[Dict[str, Any]] = None, detail: Optional[Dict[str, Any]] = None
) -> None:
    """Add ``summary`` and ``detail`` entries to what the run log records for this run.

    Entries are merged into ``args._run_summary`` and ``args._run_detail``, a later key replacing an
    earlier one. The summary goes on the run's line in ``runs.jsonl``; the detail, when there is any,
    goes to the run's own detail file.
    """
    for attribute, entries in (("_run_summary", summary), ("_run_detail", detail)):
        if entries:
            merged = dict(getattr(args, attribute, None) or {})
            merged.update(entries)
            setattr(args, attribute, merged)


def note_rows(args: argparse.Namespace, rows: Dict[str, Dict[str, Any]], *section: str) -> None:
    """Add per-item ``rows`` under the nested ``section`` of this run's detail, keeping earlier rows.

    ``note_rows(args, {"SRR1/G1": {...}}, "extractions")`` gives ``{"extractions": {"SRR1/G1": {...}}}``.
    A row is a mapping of plain values keyed by accession (or ``accession/genome``), the form ``runs``
    compares; a command that notes one row at a time keeps every row noted before it.
    """
    if not rows or not section:
        return
    detail = dict(getattr(args, "_run_detail", None) or {})
    node = detail
    for name in section:
        child = node.get(name)
        node[name] = dict(child) if isinstance(child, dict) else {}
        node = node[name]
    node.update(rows)
    args._run_detail = detail


def _parse_lines(text: str) -> Tuple[List[Tuple[str, Optional[RunRecord]]], int]:
    """Each non-empty line of ``text`` with its record (None for an unreadable one), and that count."""
    parsed: List[Tuple[str, Optional[RunRecord]]] = []
    unreadable = 0
    for line in text.splitlines():
        if not line.strip():
            continue
        try:
            record: Optional[RunRecord] = RunRecord.from_dict(json.loads(line))
        except (ValueError, TypeError, KeyError):
            record = None
            unreadable += 1
        parsed.append((line, record))
    return parsed, unreadable


def _read_text(path: Path) -> str:
    """The log's text, "" when it does not exist yet.

    A byte that is not UTF-8 is replaced (U+FFFD) rather than raised, so one damaged line never
    stops a reader or the next append: outside a JSON string it makes the line fail parsing, and
    the line is skipped like one cut short; inside a string the line is read with the replacement
    character, which a later rewrite of the log by pruning keeps in place of the byte.
    """
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except FileNotFoundError:
        return ""


def _append_line(log: Path, line: str, after_partial: bool) -> None:
    """Append one line to ``log`` and flush it to disk; start a new line after a partial one."""
    with open(log, "a", encoding="utf-8") as handle:
        if after_partial:
            handle.write("\n")
        handle.write(line + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def _prune_details(folder: Path, log: Path, command: str) -> None:
    """Keep the last ``DETAILS_KEPT_PER_COMMAND`` detail files of ``command``; null the older lines' detail."""
    parsed, _ = _parse_lines(_read_text(log))
    with_detail = [record for _, record in parsed if record and record.command == command and record.detail]
    pruned = with_detail[: max(0, len(with_detail) - DETAILS_KEPT_PER_COMMAND)]
    if not pruned:
        return
    pruned_ids: Set[str] = {record.run_id for record in pruned}
    lines = []
    for line, record in parsed:
        if record is not None and record.run_id in pruned_ids:
            # The line's own JSON is edited, not re-serialised through RunRecord, so keys a later
            # schema adds to a line are carried through.
            data = json.loads(line)
            data["detail"] = None
            line = json.dumps(data, sort_keys=True)
        lines.append(line)
    # The log is rewritten before the files go, so no line ever names a removed file for long.
    write_text_atomic(log, "\n".join(lines) + "\n", fsync=True)
    for record in pruned:
        (folder / f"{record.run_id}.json").unlink(missing_ok=True)


def record_run(
    project: PathLike,
    command: str,
    argv: Sequence[str],
    args: argparse.Namespace,
    started: datetime,
    seconds: float,
    exit_code: int,
) -> Optional[RunRecord]:
    """Append one run of ``command`` to ``project``'s run log and return its record.

    Returns None, writing nothing, when the ``run_log`` setting is off. ``argv`` and ``args`` are
    stored with secret values masked; ``args._run_summary`` and ``args._run_detail`` (see
    ``note_run``) become the record's summary and detail file. Raises ``OSError`` when the folder
    cannot be written and ``LockWaitTimeout`` (a ``DataAccessError``) when another run holds the
    log's lock for longer than ``RUN_LOG_LOCK_POLICY.wait_seconds``.
    """
    if not settings.active().run_log:
        return None
    folder = runs_dir(project)
    folder.mkdir(parents=True, exist_ok=True)
    log = folder / RUNS_FILE
    summary = _json_safe(dict(getattr(args, "_run_summary", None) or {}))
    detail = _json_safe(dict(getattr(args, "_run_detail", None) or {}))
    with held_lock(folder / f"{RUNS_FILE}.lock", RUN_LOG_LOCK_POLICY):
        text = _read_text(log)
        taken = {record.run_id for _, record in _parse_lines(text)[0] if record}
        run_id = new_run_id(command, started)
        while run_id in taken or (folder / f"{run_id}.json").exists():
            run_id = new_run_id(command, started)
        record = RunRecord(
            run_id=run_id,
            command=command,
            started=_utc_text(started),
            finished=_utc_text(started + timedelta(seconds=seconds)),
            seconds=round(float(seconds), 3),
            exit_code=int(exit_code),
            argv=settings.masked_argv(list(argv), args),
            args=json_safe_args(args),
            version=__version__,
            host=socket.gethostname(),
            pid=os.getpid(),
            summary=summary,
            detail=f"{run_id}.json" if detail else None,
        )
        if detail:
            payload = {"schema": RUN_LOG_SCHEMA, "run_id": run_id, "command": command, "detail": detail}
            write_text_atomic(folder / f"{run_id}.json", json.dumps(payload, indent=2, sort_keys=True) + "\n")
        _append_line(log, json.dumps(record.to_dict(), sort_keys=True), bool(text) and not text.endswith("\n"))
        _prune_details(folder, log, command)
    return record


def read_runs(project: PathLike, command: Optional[str] = None) -> List[RunRecord]:
    """The runs recorded for ``project`` (only ``command``'s when given), oldest first.

    Unreadable lines (a line cut short by an interrupted append) are skipped with one warning.
    Empty when the project has no run log.
    """
    log = runs_dir(project) / RUNS_FILE
    parsed, unreadable = _parse_lines(_read_text(log))
    if unreadable:
        logger.warning("Skipped %d unreadable line(s) in %s", unreadable, log)
    return [record for _, record in parsed if record is not None and (command is None or record.command == command)]


def read_detail(project: PathLike, record: RunRecord) -> Optional[Dict[str, Any]]:
    """The detail ``record``'s command noted for that run, or None when it has none or it was pruned."""
    if not record.detail:
        return None
    path = runs_dir(project) / Path(record.detail).name
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None
    except (OSError, ValueError) as e:
        logger.warning("Cannot read the detail of run %s from %s: %s", record.run_id, path, e)
        return None
    detail = payload.get("detail") if isinstance(payload, dict) else None
    return detail if isinstance(detail, dict) else None


def resolve_run(records: Sequence[RunRecord], selector: str) -> RunRecord:
    """The record ``selector`` names among ``records`` (oldest first).

    ``latest`` is the last record and ``previous`` the one before it; otherwise the selector is a
    run identifier or a prefix that matches exactly one. Raises ``ValidationError`` when nothing,
    or more than one run, matches.
    """
    if not records:
        raise ValidationError(f"Cannot select run '{selector}': no runs are recorded")
    if selector == "latest":
        return records[-1]
    if selector == "previous":
        if len(records) < 2:
            raise ValidationError("No previous run: only one run is recorded")
        return records[-2]
    for record in records:
        if record.run_id == selector:
            return record
    matches = [record for record in records if record.run_id.startswith(selector)]
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise ValidationError(f"No run matches '{selector}'")
    shown = ", ".join(record.run_id for record in matches[:5])
    raise ValidationError(f"'{selector}' matches {len(matches)} runs ({shown}); give more of the run identifier")
