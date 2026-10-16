"""
SQLite catalogue over the shared data store: ``<root>/catalog.sqlite``.

The catalogue is a rebuildable index over every sidecar in the store, plus two
cross-project tables (``projects`` and ``usage``) that record which project used
which dataset, for which target genome, at which pipeline stage. It answers
questions no single sidecar or project registry can: which projects used one
accession, which downloaded datasets no project references any more, how many
bytes are held per target genome, and which datasets were analysed for one
genome across every project. Because it is rebuildable from the sidecars
(``reindex``), losing or corrupting ``catalog.sqlite`` never loses data; it
only loses the cross-project ``usage`` history, which lives only here.

Writes are serialised across processes with the lock-file mechanism every
PREPARE lock uses (``prepare.utils.lockfile``), on a
``catalog.sqlite.lock`` file next to the database, so two projects on one
machine or a shared network volume never interleave writes. Reads may open the
database directly without taking the lock, but only one that already exists:
creating and migrating are ``catalog_write``'s job, so a read never runs DDL
outside the lock. The database journals in DELETE mode, never WAL, because a
shared store is expected to live on a NAS and SQLite documents WAL as unsafe
over NFS and SMB while accepting the pragma there regardless.
"""

import functools
import logging
import socket
import sqlite3
from contextlib import contextmanager
from datetime import datetime, timezone
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Tuple, TypeVar

from prepare.core.constants import (
    CATALOG_BUSY_TIMEOUT_SECONDS,
    CATALOG_LOCK_STALE_SECONDS,
    CATALOG_LOCK_WAIT_SECONDS,
    SHORT_LOCK_HEARTBEAT_SECONDS,
    SHORT_LOCK_POLL_SECONDS,
)
from prepare.core.exceptions import DataAccessError
from prepare.core.settings import setting_or
from prepare.store.layout import StorePaths
from prepare.store.sidecar import Sidecar
from prepare.utils.lockfile import LockPolicy, held_lock, verify_held

logger = logging.getLogger(__name__)

# Bump when the catalogue's schema changes shape.
SCHEMA_VERSION = 1

# ``store_meta`` key set by ``store_reindex`` when it rebuilds a catalogue that holds datasets
# but restores no project: every dataset then looks unused, so ``store_gc`` refuses to run
# until the user confirms with ``--accept-rebuilt`` (or a later reindex restores a project).
REBUILT_WITHOUT_PROJECTS = "rebuilt_without_projects"

_SCHEMA_STATEMENTS = (
    """
    CREATE TABLE IF NOT EXISTS store_meta (
        key TEXT PRIMARY KEY,
        value TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS datasets (
        accession TEXT PRIMARY KEY,
        state TEXT,
        layout TEXT,
        compression TEXT,
        reads_per_mate INTEGER,
        bases_total INTEGER,
        bytes_total INTEGER,
        ncbi_spots INTEGER,
        ncbi_bases INTEGER,
        completeness TEXT,
        downloaded TEXT,
        tool_version TEXT,
        updated TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS files (
        accession TEXT NOT NULL REFERENCES datasets(accession) ON DELETE CASCADE,
        name TEXT NOT NULL,
        bytes INTEGER,
        md5 TEXT,
        reads INTEGER,
        PRIMARY KEY (accession, name)
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS projects (
        project_id TEXT PRIMARY KEY,
        name TEXT,
        path TEXT,
        registry TEXT,
        created TEXT,
        last_seen TEXT
    )
    """,
    """
    CREATE TABLE IF NOT EXISTS usage (
        accession TEXT NOT NULL REFERENCES datasets(accession) ON DELETE CASCADE,
        project_id TEXT NOT NULL REFERENCES projects(project_id) ON DELETE CASCADE,
        genome_id TEXT NOT NULL DEFAULT '',
        stage TEXT NOT NULL,
        first_used TEXT,
        last_used TEXT,
        detail TEXT,
        PRIMARY KEY (accession, project_id, genome_id, stage)
    )
    """,
    "CREATE INDEX IF NOT EXISTS idx_usage_project ON usage(project_id)",
    "CREATE INDEX IF NOT EXISTS idx_usage_genome ON usage(genome_id)",
    """
    CREATE VIEW IF NOT EXISTS unused_datasets AS
    SELECT d.accession AS accession
    FROM datasets d
    LEFT JOIN usage u ON u.accession = d.accession
    WHERE u.accession IS NULL
    """,
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


_F = TypeVar("_F", bound=Callable[..., Any])


def _wrap_sqlite_errors(func: _F) -> _F:
    """Translate any ``sqlite3.Error`` the wrapped call raises into ``DataAccessError``.

    Every public method of ``Catalog`` (and ``catalog_write``) is wrapped with this so a
    caller never has to catch ``sqlite3.Error`` directly; a corrupt database, a locked file,
    or a disk I/O error all surface the same way as every other PREPARE data-access failure.
    """

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        try:
            return func(*args, **kwargs)
        except sqlite3.Error as e:
            raise DataAccessError(f"Catalog operation '{func.__name__}' failed: {e}") from e

    return wrapper  # type: ignore[return-value]


class Catalog:
    """One open connection to ``catalog.sqlite``.

    Use as a context manager. Opening never migrates the schema and never takes
    the write lock; both are the caller's job (``catalog_write`` does both for a
    write session). A plain ``with Catalog(paths) as catalog:`` is for reads
    against a catalogue that already exists, and refuses to open one that does
    not: ``create=True`` (what ``catalog_write`` passes) is the only way to
    bring a new database into being.
    """

    def __init__(self, paths: StorePaths, create: bool = False):
        self.paths = paths
        self.create = create
        self._conn: Optional[sqlite3.Connection] = None
        # Suspended by journal.replay() while it feeds journaled records back in, so a
        # rebuilt catalogue does not duplicate the lines it is replaying.
        self.journal_enabled: bool = True

    @property
    def conn(self) -> sqlite3.Connection:
        """The open connection; raises if used outside a ``with Catalog(...) as catalog:`` block."""
        if self._conn is None:
            raise RuntimeError("Catalog is not open; use it as a context manager")
        return self._conn

    @_wrap_sqlite_errors
    def __enter__(self) -> "Catalog":
        if not self.create and not self.paths.catalog.is_file():
            raise DataAccessError(f"No store catalogue at {self.paths.catalog}; build one with: prepare store_reindex")
        self.paths.root.mkdir(parents=True, exist_ok=True)
        # SQLite's own busy handler, on top of catalog_write's O_EXCL lock file: the implicit
        # default (5 s) is shorter than a store_reindex replay can hold the database, and
        # shorter than the window this connection may need to wait out another connection's
        # transaction.
        conn = sqlite3.connect(str(self.paths.catalog), timeout=CATALOG_BUSY_TIMEOUT_SECONDS)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON")
        # DELETE, never WAL: SQLite documents WAL as unsafe over NFS and SMB (its shared-memory
        # index needs coherent mmap across hosts) and the pragma succeeds there anyway, so a
        # shared store on a NAS would silently run in an unsupported mode. Every write here is
        # already serialised by catalog.sqlite.lock, so WAL would buy nothing.
        conn.execute("PRAGMA journal_mode = DELETE")
        self._conn = conn
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    # ------------------------------------------------------------------ schema

    def _add_missing_columns(self) -> None:
        """Add columns later schema versions introduced to a database created by an earlier one.

        ``projects.hostname`` records which machine wrote a project's row. A project's registry
        living on another workstation's disk looks missing from here, so staleness cannot be
        read the same way for a row written elsewhere; recording the host at least makes that
        visible in ``store_status`` and ``store_gc``.
        """
        existing = {row["name"] for row in self.conn.execute("PRAGMA table_info(projects)").fetchall()}
        if "hostname" not in existing:
            self.conn.execute("ALTER TABLE projects ADD COLUMN hostname TEXT")

    def _refuse_newer_schema(self) -> None:
        """Raise ``DataAccessError`` when the catalogue records a schema version above ``SCHEMA_VERSION``.

        A catalogue written by a later release may hold tables or columns this release does not
        know, and migrating or writing it here could damage them. An unreadable version value is
        refused the same way, since it cannot be shown to be one this release understands.
        """
        has_meta = self.conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'store_meta'"
        ).fetchone()
        if has_meta is None:
            return
        row = self.conn.execute("SELECT value FROM store_meta WHERE key = 'schema_version'").fetchone()
        if row is None:
            return
        try:
            version = int(row["value"])
        except (TypeError, ValueError):
            version = None
        if version is None or version > SCHEMA_VERSION:
            raise DataAccessError(
                f"The store catalogue at {self.paths.catalog} has schema version {row['value']}, but this "
                f"prepare supports schema version {SCHEMA_VERSION} at most; upgrade prepare to use this store"
            )

    @_wrap_sqlite_errors
    def migrate(self) -> None:
        """Create the catalogue schema if it does not already exist. Safe to call repeatedly.

        Refuses, before any migration step, a catalogue whose recorded schema version is newer
        than this release's (``_refuse_newer_schema``).
        """
        self._refuse_newer_schema()
        for statement in _SCHEMA_STATEMENTS:
            self.conn.execute(statement)
        self._add_missing_columns()
        self.conn.execute(
            "INSERT OR IGNORE INTO store_meta (key, value) VALUES ('schema_version', ?)",
            (str(SCHEMA_VERSION),),
        )

    # -------------------------------------------------------------------- meta

    @_wrap_sqlite_errors
    def get_meta(self, key: str) -> Optional[str]:
        """The ``store_meta`` value for ``key``, or None when it is not set."""
        row = self.conn.execute("SELECT value FROM store_meta WHERE key = ?", (key,)).fetchone()
        return None if row is None else row["value"]

    @_wrap_sqlite_errors
    def set_meta(self, key: str, value: str) -> None:
        """Set the ``store_meta`` value for ``key``, replacing any earlier value."""
        self.conn.execute("INSERT OR REPLACE INTO store_meta (key, value) VALUES (?, ?)", (key, value))

    @_wrap_sqlite_errors
    def delete_meta(self, key: str) -> None:
        """Remove ``key`` from ``store_meta``; a key that is not set is left as it is."""
        self.conn.execute("DELETE FROM store_meta WHERE key = ?", (key,))

    # ---------------------------------------------------------------- datasets

    @_wrap_sqlite_errors
    def upsert_dataset(self, sidecar: Sidecar) -> None:
        """Insert or update the ``datasets``/``files`` rows for one sidecar.

        ``bytes_total`` is the sum of the sidecar's file sizes; ``completeness``
        stores the verdict string only (the method and ratio are not indexed).
        The file rows for this accession are replaced wholesale, since a sidecar
        always lists every file currently on disk for the accession.
        """
        files = sidecar.files or []
        bytes_total = sum(int(f.get("bytes") or 0) for f in files)
        now = _now()

        self.conn.execute(
            """
            INSERT INTO datasets (
                accession, state, layout, compression, reads_per_mate, bases_total,
                bytes_total, ncbi_spots, ncbi_bases, completeness, downloaded,
                tool_version, updated
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(accession) DO UPDATE SET
                state=excluded.state,
                layout=excluded.layout,
                compression=excluded.compression,
                reads_per_mate=excluded.reads_per_mate,
                bases_total=excluded.bases_total,
                bytes_total=excluded.bytes_total,
                ncbi_spots=excluded.ncbi_spots,
                ncbi_bases=excluded.ncbi_bases,
                completeness=excluded.completeness,
                downloaded=excluded.downloaded,
                tool_version=excluded.tool_version,
                updated=excluded.updated
            """,
            (
                sidecar.accession,
                sidecar.state,
                sidecar.layout,
                sidecar.compression,
                sidecar.reads_per_mate,
                sidecar.bases_total,
                bytes_total,
                sidecar.ncbi.get("spots"),
                sidecar.ncbi.get("bases"),
                sidecar.completeness.get("verdict"),
                sidecar.downloaded,
                sidecar.tool_version,
                now,
            ),
        )
        self.conn.execute("DELETE FROM files WHERE accession = ?", (sidecar.accession,))
        self.conn.executemany(
            "INSERT INTO files (accession, name, bytes, md5, reads) VALUES (?, ?, ?, ?, ?)",
            [(sidecar.accession, f.get("name"), f.get("bytes"), f.get("md5"), f.get("reads")) for f in files],
        )

    @_wrap_sqlite_errors
    def get_dataset(self, accession: str) -> Optional[Dict[str, Any]]:
        """The ``datasets`` row for ``accession`` as a plain dict (with its ``files``), or None."""
        row = self.conn.execute("SELECT * FROM datasets WHERE accession = ?", (accession,)).fetchone()
        if row is None:
            return None
        result = dict(row)
        file_rows = self.conn.execute(
            "SELECT name, bytes, md5, reads FROM files WHERE accession = ? ORDER BY name",
            (accession,),
        ).fetchall()
        result["files"] = [dict(f) for f in file_rows]
        return result

    @_wrap_sqlite_errors
    def delete_dataset(self, accession: str) -> None:
        """Remove ``accession``'s row from ``datasets``.

        The foreign key cascade (``ON DELETE CASCADE``) removes its ``files`` and ``usage``
        rows along with it. Used by ``store_gc --yes`` once the accession's folder has been
        removed from disk; a dataset the catalogue no longer backs with real data must never
        keep a row here. A no-op when ``accession`` is not present.
        """
        self.conn.execute("DELETE FROM datasets WHERE accession = ?", (accession,))

    # ---------------------------------------------------------------- projects

    @_wrap_sqlite_errors
    def upsert_project(
        self, project_id: str, name: str, path: str, registry: str, hostname: Optional[str] = None
    ) -> None:
        """Insert or update one project's row, keeping its original ``created`` timestamp.

        The row records the host that wrote it, so a shared store used from two machines can
        say where a project it cannot see from here was last written. ``hostname`` defaults to
        this machine (``socket.gethostname()``); ``journal.replay`` passes the host a journal
        line actually recorded, so rebuilding a lost catalogue on a different machine does not
        overwrite every project's real host with the rebuilding machine's.
        """
        now = _now()
        host = hostname if hostname is not None else socket.gethostname()
        self.conn.execute(
            """
            INSERT INTO projects (project_id, name, path, registry, created, last_seen, hostname)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(project_id) DO UPDATE SET
                name=excluded.name,
                path=excluded.path,
                registry=excluded.registry,
                last_seen=excluded.last_seen,
                hostname=excluded.hostname
            """,
            (project_id, name, path, registry, now, now, host),
        )
        if self.journal_enabled:
            from prepare.store import journal

            journal.append_project(self.paths, project_id, name, path, registry, hostname=host)

    # ------------------------------------------------------------------- usage

    @_wrap_sqlite_errors
    def record_usage(
        self,
        accession: str,
        project_id: str,
        genome_id: str,
        stage: str,
        detail: str = "",
        at: Optional[str] = None,
        last_used: Optional[str] = None,
    ) -> None:
        """Record that ``project_id`` used ``accession`` (for ``genome_id``, at ``stage``).

        ``first_used`` is set once and kept at the earliest timestamp seen for the same
        (accession, project_id, genome_id, stage) key; ``last_used`` and ``detail`` advance to
        the latest. Calls can arrive out of chronological order (``journal.replay`` feeds lines
        back in file order, not time order), so both bounds are computed with SQL ``min``/``max``
        across the existing row and the incoming call rather than blindly overwritten. If
        ``accession`` is not yet in ``datasets`` (usage can be recorded before a dataset is
        catalogued, e.g. a link recorded ahead of the next reindex), a minimal placeholder row
        with ``state="unknown"`` is inserted first so the foreign key from ``usage`` to
        ``datasets`` is satisfied; a later ``upsert_dataset`` or ``reindex`` fills it in properly.

        ``at`` sets this call's own timestamp (``first_used`` candidate) instead of now;
        ``journal.replay`` passes the journal line's own ``at`` so a rebuilt catalogue keeps
        the date usage actually happened rather than the date it was replayed. ``last_used``
        additionally sets this call's ``last_used`` candidate when the caller already knows a
        later timestamp than ``at`` (a backfilled or replayed line that itself spans a
        first/last pair); it defaults to ``at`` (or now) when not given, which is what an
        ordinary, non-replayed call wants.

        Raises ``DataAccessError`` if ``project_id`` is not a known project (from
        ``upsert_project``); usage is never recorded against a fabricated project.
        """
        known = self.conn.execute("SELECT 1 FROM projects WHERE project_id = ?", (project_id,)).fetchone()
        if known is None:
            raise DataAccessError(f"Unknown project_id: {project_id}")

        # A live call (neither ``at`` nor ``last_used`` given) is this host's own clock
        # reporting a real event happening now; it always wins on detail, even when the
        # row already stored a later last_used from a host whose clock ran ahead. Only a
        # replay/backfill call (one that passes ``at``) defers to a stored, later timestamp.
        is_live_call = at is None and last_used is None
        detail_clause = (
            "detail=excluded.detail"
            if is_live_call
            else "detail=CASE WHEN last_used IS NULL OR excluded.last_used >= last_used "
            "THEN excluded.detail ELSE detail END"
        )

        now = at or _now()
        latest = last_used or now
        genome_id = genome_id or ""

        self.conn.execute(
            "INSERT OR IGNORE INTO datasets (accession, state, updated) VALUES (?, ?, ?)",
            (accession, "unknown", now),
        )
        self.conn.execute(
            f"""
            INSERT INTO usage (accession, project_id, genome_id, stage, first_used, last_used, detail)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(accession, project_id, genome_id, stage) DO UPDATE SET
                first_used=min(coalesce(first_used, excluded.first_used), excluded.first_used),
                last_used=max(coalesce(last_used, excluded.last_used), excluded.last_used),
                {detail_clause}
            """,
            (accession, project_id, genome_id, stage, now, latest, detail),
        )
        if self.journal_enabled:
            from prepare.store import journal

            journal.append_usage(
                self.paths,
                accession,
                project_id,
                genome_id,
                stage,
                detail,
                at=now,
                last_used=latest if last_used is not None else None,
            )

    # ----------------------------------------------------------------- queries

    @_wrap_sqlite_errors
    def projects_for(self, accession: str) -> List[Dict[str, Any]]:
        """Every project that has recorded usage of ``accession``, sorted by project id."""
        rows = self.conn.execute(
            """
            SELECT DISTINCT p.*
            FROM projects p
            JOIN usage u ON u.project_id = p.project_id
            WHERE u.accession = ?
            ORDER BY p.project_id
            """,
            (accession,),
        ).fetchall()
        return [dict(row) for row in rows]

    @_wrap_sqlite_errors
    def datasets_for_project(self, project_id: str) -> List[Dict[str, Any]]:
        """Every dataset ``project_id`` has recorded usage of, sorted by accession."""
        rows = self.conn.execute(
            """
            SELECT DISTINCT d.*
            FROM datasets d
            JOIN usage u ON u.accession = d.accession
            WHERE u.project_id = ?
            ORDER BY d.accession
            """,
            (project_id,),
        ).fetchall()
        return [dict(row) for row in rows]

    @_wrap_sqlite_errors
    def unused(self) -> List[str]:
        """Accessions in ``datasets`` with no usage row at all, sorted."""
        rows = self.conn.execute("SELECT accession FROM unused_datasets ORDER BY accession").fetchall()
        return [row["accession"] for row in rows]

    @_wrap_sqlite_errors
    def bytes_by_genome(self) -> Dict[str, int]:
        """Total ``bytes_total`` per genome, summed over distinct (accession, genome_id) usage pairs.

        A dataset used by several projects for the same genome is counted once,
        since it occupies the bytes only once on disk.
        """
        rows = self.conn.execute("""
            SELECT pair.genome_id AS genome_id, SUM(d.bytes_total) AS total
            FROM (SELECT DISTINCT genome_id, accession FROM usage) pair
            JOIN datasets d ON d.accession = pair.accession
            GROUP BY pair.genome_id
            ORDER BY pair.genome_id
            """).fetchall()
        return {row["genome_id"]: (row["total"] or 0) for row in rows}

    @_wrap_sqlite_errors
    def datasets_for_genome(self, genome_id: str) -> List[Tuple[str, str]]:
        """Distinct (accession, project_id) pairs with recorded usage for ``genome_id``."""
        rows = self.conn.execute(
            """
            SELECT DISTINCT accession, project_id
            FROM usage
            WHERE genome_id = ?
            ORDER BY accession, project_id
            """,
            (genome_id,),
        ).fetchall()
        return [(row["accession"], row["project_id"]) for row in rows]

    # ----------------------------------------------------------------- rebuild

    @_wrap_sqlite_errors
    def reindex(self, sidecars: Iterable[Sidecar]) -> int:
        """Rebuild ``datasets``/``files`` from ``sidecars``; keep ``projects`` and ``usage``.

        Every accession not present in ``sidecars`` is removed from ``datasets``;
        the foreign key cascade then removes its ``files`` and ``usage`` rows too.
        Every accession that is present is upserted, so its ``usage`` history
        (referencing the same, still-present, accession) is left untouched.
        Returns how many sidecars were indexed.
        """
        sidecar_list = list(sidecars)
        keep = {sidecar.accession for sidecar in sidecar_list}

        existing = {row["accession"] for row in self.conn.execute("SELECT accession FROM datasets").fetchall()}
        for accession in existing - keep:
            if (self.paths.sra / accession).is_dir():
                # The folder is there; only its sidecar is missing or unreadable. Dropping the
                # row would take the cross-project usage history with it (the foreign key
                # cascades), and that history exists nowhere else.
                logger.warning("%s: keeping its catalogue row, the files are still in the store", accession)
                continue
            self.conn.execute("DELETE FROM datasets WHERE accession = ?", (accession,))

        for sidecar in sidecar_list:
            self.upsert_dataset(sidecar)

        return len(sidecar_list)


@contextmanager
def catalog_write(paths: StorePaths) -> Iterator[Catalog]:
    """Open the catalogue for a write session, serialised across processes.

    Acquires ``paths.catalog_lock`` (waiting up to the ``catalog_lock_wait`` runtime setting, by
    default ``CATALOG_LOCK_WAIT_SECONDS``, then raising ``LockHeld`` naming the holder), opens the
    catalogue, migrates its schema, backfills the journal from any ``projects``/``usage`` rows
    that predate it (``journal.backfill_from_catalog``, a no-op once the journal already has
    project lines), yields it for the caller to write through, commits on a clean exit, and
    always releases the lock. If the block raises, the connection is closed without committing
    (uncommitted changes are discarded) and the lock is still released. The commit happens only
    while the lock is still this process's (``verify_held``): a holder stalled past the stale
    window gets ``LockLost``, a ``DataAccessError``, and its changes are discarded.

    Not re-entrant: nesting a second ``catalog_write`` against the same store root inside
    this block's body, in the same thread, raises ``LockReentry`` at once.
    """
    paths.root.mkdir(parents=True, exist_ok=True)
    policy = LockPolicy(
        what="Store catalogue",
        stale_seconds=CATALOG_LOCK_STALE_SECONDS,
        wait_seconds=setting_or("catalog_lock_wait", CATALOG_LOCK_WAIT_SECONDS),
        poll_seconds=SHORT_LOCK_POLL_SECONDS,
        heartbeat_seconds=SHORT_LOCK_HEARTBEAT_SECONDS,
    )
    with held_lock(paths.catalog_lock, policy) as lock:
        with Catalog(paths, create=True) as catalog:
            catalog.migrate()
            from prepare.store import journal

            journal.backfill_from_catalog(paths, catalog)
            yield catalog
            verify_held(lock)
            try:
                catalog.conn.commit()
            except sqlite3.Error as e:
                raise DataAccessError(f"Cannot commit catalog write: {e}") from e
