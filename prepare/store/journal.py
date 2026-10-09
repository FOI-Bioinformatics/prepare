"""Append-only journal of the store's project and usage records.

The SQLite catalogue is the working copy of "which project used which dataset", but it is the
only copy: sidecars know nothing about projects. ``store_reindex`` rebuilds the catalogue from
sidecars and would otherwise come back with no projects and no usage, after which ``store_gc``
sees every dataset as unused. Every ``upsert_project`` and ``record_usage`` therefore also
appends one JSON line here, and ``replay`` feeds those lines back into a rebuilt catalogue.
Appends happen inside ``catalog_write``, which holds the store-wide write lock.
"""

import json
import logging
import socket
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterator, Optional, Set, Tuple

from prepare.core.exceptions import DataAccessError
from prepare.store.layout import StorePaths

logger = logging.getLogger(__name__)

PROJECTS_FILE = "projects.jsonl"
USAGE_FILE = "usage.jsonl"


def _append(paths: StorePaths, name: str, record: Dict[str, Any], at: Optional[str] = None) -> None:
    """Append one record, stamped with ``at`` (or now, when the caller has no timestamp of its
    own to preserve)."""
    paths.journal.mkdir(parents=True, exist_ok=True)
    record = dict(record, at=at or datetime.now(timezone.utc).isoformat())
    with open(paths.journal / name, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def append_project(
    paths: StorePaths, project_id: str, name: str, path: str, registry: str, hostname: Optional[str] = None
) -> None:
    """Append one project record (mirroring ``Catalog.upsert_project``'s arguments).

    Records the host that wrote this row, so replaying the line later (possibly on a
    different machine, rebuilding a lost catalogue) can restore the host that actually wrote
    the project rather than the host doing the rebuild. ``hostname`` should be the same value
    the caller just wrote to ``catalog.sqlite`` (``Catalog.upsert_project`` passes its own
    resolved host through); it defaults to this machine's own host only when the caller has
    none to pass (e.g. a direct call outside ``Catalog``).
    """
    _append(
        paths,
        PROJECTS_FILE,
        {
            "project_id": project_id,
            "name": name,
            "path": path,
            "registry": registry,
            "hostname": hostname if hostname is not None else socket.gethostname(),
        },
    )


def append_usage(
    paths: StorePaths,
    accession: str,
    project_id: str,
    genome_id: str,
    stage: str,
    detail: str,
    at: Optional[str] = None,
    last_used: Optional[str] = None,
) -> None:
    """Append one usage record (mirroring ``Catalog.record_usage``'s arguments).

    ``at`` should be the same timestamp the caller just wrote to ``first_used`` in
    ``catalog.sqlite`` (``Catalog.record_usage`` passes it through), so replaying this line
    later restores the date the usage actually happened rather than the date it was replayed.
    ``last_used`` is written as its own key only when given (an ordinary call has none, and a
    single ``at`` already covers both bounds); a backfilled or explicitly out-of-order line
    that carries a later ``last_used`` than ``at`` passes it so replay restores both bounds.
    """
    record = {
        "accession": accession,
        "project_id": project_id,
        "genome_id": genome_id,
        "stage": stage,
        "detail": detail,
    }
    if last_used is not None:
        record["last_used"] = last_used
    _append(paths, USAGE_FILE, record, at=at)


def _lines(path: Path) -> Iterator[Dict[str, Any]]:
    """Every parseable JSON object on its own line in ``path``; skips and warns on the rest.

    A line that is not valid JSON, or that parses to something other than a JSON object (a
    list, a number, ``null``), is not a usable record: it is skipped with a warning naming the
    file and line number rather than raised, so one damaged line never blocks every record
    after it.
    """
    if not path.is_file():
        return
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                logger.warning("Skipping unreadable journal line %s:%d", path, number)
                continue
            if not isinstance(record, dict):
                logger.warning("Skipping journal line %s:%d: not a JSON object", path, number)
                continue
            yield record


def _has_records(path: Path) -> bool:
    """Whether ``path`` holds at least one parseable JSON-object line."""
    return any(True for _ in _lines(path))


def replay(paths: StorePaths, catalog: Any) -> Tuple[int, int]:
    """Feed every journaled project and usage record into ``catalog``; returns the two counts.

    The catalogue calls back into ``append_*`` when it writes, so replay runs with journaling
    suspended (``catalog.journal_enabled = False``) to avoid duplicating the file.

    The project count is the number of *distinct* ``project_id`` values upserted, not the
    number of journal lines replayed: ``upsert_project`` appends one line per call, including
    a call that only refreshes an existing project's ``last_seen``/``hostname``, so a project
    touched more than once would otherwise be over-counted.

    A usage line whose ``project_id`` was never restored (its project's line is missing,
    corrupt, or absent from the journal entirely) is skipped rather than raised: without this,
    one such line would make ``record_usage`` raise ``DataAccessError("Unknown project_id")``,
    aborting the whole replay, including every dataset row ``reindex`` already rebuilt in the
    same transaction. Skipped lines are counted and reported in one warning, not one per line.
    A journal line missing ``project_id`` (for projects) or ``accession``/``project_id`` (for
    usage) is skipped the same way, since it cannot be replayed either.
    """
    project_ids: Set[str] = set()
    usage = 0
    skipped_usage = 0
    previous = getattr(catalog, "journal_enabled", True)
    catalog.journal_enabled = False
    try:
        for record in _lines(paths.journal / PROJECTS_FILE):
            if "project_id" not in record:
                logger.warning("Skipping project journal line without project_id: %r", record)
                continue
            catalog.upsert_project(
                record["project_id"],
                record.get("name", ""),
                record.get("path", ""),
                record.get("registry", ""),
                hostname=record.get("hostname"),
            )
            project_ids.add(record["project_id"])
        for record in _lines(paths.journal / USAGE_FILE):
            if "accession" not in record or "project_id" not in record:
                logger.warning("Skipping usage journal line without accession/project_id: %r", record)
                continue
            try:
                catalog.record_usage(
                    record["accession"],
                    record["project_id"],
                    record.get("genome_id", ""),
                    record.get("stage", ""),
                    record.get("detail", ""),
                    at=record.get("at"),
                    last_used=record.get("last_used"),
                )
                usage += 1
            except DataAccessError:
                skipped_usage += 1
        if skipped_usage:
            logger.warning(
                "Skipped %d usage record(s) referencing a project that could not be restored; "
                "those datasets will look unused until that project runs store_init or store_link again",
                skipped_usage,
            )
    finally:
        catalog.journal_enabled = previous
    return len(project_ids), usage


def backfill_from_catalog(paths: StorePaths, catalog: Any) -> Tuple[int, int]:
    """Copy an existing catalogue's ``projects``/``usage`` rows into the journal, once.

    A store that started using PREPARE before this journal existed has its project and usage
    history only in ``catalog.sqlite``; the journal starts empty for it. Without this, losing
    that database would still lose every such project (``store_reindex`` has nothing to
    replay), and ``store_gc`` would then see every one of its datasets as unused. Called from
    ``catalog_write`` right after ``migrate()``, under the store-wide write lock, so it runs at
    most once per store: once ``projects.jsonl`` holds at least one project line (whether from
    a real write or from this backfill), this is a no-op. Writes projects before usage, since
    usage rows reference a project by id. Uses ``append_project``/``append_usage`` directly
    (not ``catalog.upsert_project``/``record_usage``), so nothing already in ``catalog.sqlite``
    is written back to it and nothing is appended twice. Each line carries the catalogue row's
    own ``hostname`` and ``first_used`` (else ``last_used``) time as ``at``, plus the row's own
    ``last_used`` (which can be later than ``first_used`` for a row touched more than once
    before the journal existed), rather than the host and time of the backfill, so a later
    replay restores both bounds the catalogue held; a row with no recorded host falls back to
    this machine, as ``append_project`` does.
    """
    projects_path = paths.journal / PROJECTS_FILE
    if _has_records(projects_path):
        return 0, 0

    try:
        project_rows = catalog.conn.execute(
            "SELECT project_id, name, path, registry, hostname FROM projects"
        ).fetchall()
    except sqlite3.Error as e:
        raise DataAccessError(str(e)) from e
    if not project_rows:
        return 0, 0

    for row in project_rows:
        append_project(
            paths,
            row["project_id"],
            row["name"] or "",
            row["path"] or "",
            row["registry"] or "",
            hostname=row["hostname"],
        )

    try:
        usage_rows = catalog.conn.execute(
            "SELECT accession, project_id, genome_id, stage, detail, first_used, last_used FROM usage"
        ).fetchall()
    except sqlite3.Error as e:
        raise DataAccessError(str(e)) from e
    for row in usage_rows:
        append_usage(
            paths,
            row["accession"],
            row["project_id"],
            row["genome_id"] or "",
            row["stage"] or "",
            row["detail"] or "",
            at=row["first_used"] or row["last_used"],
            last_used=row["last_used"],
        )

    return len(project_rows), len(usage_rows)
