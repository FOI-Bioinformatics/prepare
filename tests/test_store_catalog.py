"""
Tests for prepare.store.catalog: the rebuildable SQLite index over sidecars,
plus the cross-project ``projects`` and ``usage`` tables.
"""

import sqlite3
import threading
import time

import pytest

from prepare.core.exceptions import DataAccessError
from prepare.store.catalog import Catalog, catalog_write
from prepare.store.layout import init_store
from prepare.store.sidecar import Sidecar


def _sidecar(accession="SRR1", state="complete", bytes_r1=1000, bytes_r2=1000, spots=10, bases=1400):
    return Sidecar(
        accession=accession,
        state=state,
        layout="PAIRED",
        downloaded="2026-09-06T00:00:00+00:00",
        tool="fasterq-dump",
        tool_version="3.0.0",
        compression="gzip",
        files=[
            {"name": f"{accession}_1.fastq.gz", "bytes": bytes_r1, "md5": "a" * 32, "reads": 10},
            {"name": f"{accession}_2.fastq.gz", "bytes": bytes_r2, "md5": "b" * 32, "reads": 10},
        ],
        reads_per_mate=10,
        bases_total=bases,
        ncbi={"spots": spots, "bases": bases, "size": 12345, "layout": "PAIRED", "files": []},
        completeness={"method": "spots", "ratio": 1.0, "verdict": "complete"},
        stats={},
        stats_computed=None,
    )


@pytest.fixture
def paths(tmp_path):
    return init_store(tmp_path / "store")


def test_migrate_creates_schema_and_is_idempotent(paths):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.migrate()  # idempotent, no error

        tables = {
            row["name"]
            for row in catalog._conn.execute(
                "SELECT name FROM sqlite_master WHERE type IN ('table', 'view')"
            ).fetchall()
        }
    assert {"store_meta", "datasets", "files", "projects", "usage", "unused_datasets"} <= tables


def test_upsert_dataset_then_get(paths):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_dataset(_sidecar())

        row = catalog.get_dataset("SRR1")

    assert row is not None
    assert row["accession"] == "SRR1"
    assert row["state"] == "complete"
    assert row["layout"] == "PAIRED"
    assert row["compression"] == "gzip"
    assert row["reads_per_mate"] == 10
    assert row["bases_total"] == 1400
    assert row["bytes_total"] == 2000
    assert row["ncbi_spots"] == 10
    assert row["ncbi_bases"] == 1400
    assert row["completeness"] == "complete"
    assert row["tool_version"] == "3.0.0"


def test_get_dataset_returns_none_when_absent(paths):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        assert catalog.get_dataset("SRR-missing") is None


def test_upsert_dataset_replaces_files_and_totals(paths):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_dataset(_sidecar(bytes_r1=1000, bytes_r2=1000))
        catalog.upsert_dataset(_sidecar(state="partial", bytes_r1=500, bytes_r2=500))

        row = catalog.get_dataset("SRR1")
        file_count = catalog._conn.execute("SELECT COUNT(*) AS n FROM files WHERE accession = ?", ("SRR1",)).fetchone()[
            "n"
        ]

    assert row["state"] == "partial"
    assert row["bytes_total"] == 1000
    assert file_count == 2


def test_record_usage_keeps_first_used_updates_last_used(paths):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_project("proj1", "Wolbachia", "/projects/wolbachia", "prepare_registry.json")
        catalog.upsert_dataset(_sidecar())

        catalog.record_usage("SRR1", "proj1", "wMel", "downloaded", detail="first")
        first = dict(
            catalog._conn.execute(
                "SELECT * FROM usage WHERE accession=? AND project_id=? AND genome_id=? AND stage=?",
                ("SRR1", "proj1", "wMel", "downloaded"),
            ).fetchone()
        )

        time.sleep(0.01)
        catalog.record_usage("SRR1", "proj1", "wMel", "downloaded", detail="second")
        second = dict(
            catalog._conn.execute(
                "SELECT * FROM usage WHERE accession=? AND project_id=? AND genome_id=? AND stage=?",
                ("SRR1", "proj1", "wMel", "downloaded"),
            ).fetchone()
        )

    assert first["first_used"] == second["first_used"]
    assert second["last_used"] >= first["last_used"]
    assert second["detail"] == "second"


def test_record_usage_with_explicit_at_keeps_earliest_first_used_latest_last_used(paths):
    """journal.replay feeds record_usage calls back in journal order, which need not be
    chronological (a journal can be replayed out of order, or rebuilt from a backfill that
    carries an explicit first/last pair). first_used must land on the earliest ``at`` seen
    and last_used on the latest, regardless of which call happens first. A replay/backfill
    call (one that passes ``at``) whose own timestamp is older than what is already stored
    must also leave the stored detail alone, since it is not the most recent real event."""
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_project("proj1", "Wolbachia", "/projects/wolbachia", "prepare_registry.json")

        catalog.record_usage("SRR1", "proj1", "wMel", "downloaded", detail="newer", at="2026-09-02T00:00:00+00:00")
        catalog.record_usage("SRR1", "proj1", "wMel", "downloaded", detail="older", at="2026-09-01T00:00:00+00:00")

        row = dict(
            catalog._conn.execute(
                "SELECT * FROM usage WHERE accession=? AND project_id=? AND genome_id=? AND stage=?",
                ("SRR1", "proj1", "wMel", "downloaded"),
            ).fetchone()
        )

    assert row["first_used"] == "2026-09-01T00:00:00+00:00"
    assert row["last_used"] == "2026-09-02T00:00:00+00:00"
    assert row["detail"] == "newer"


def test_record_usage_live_call_refreshes_detail_despite_lagging_clock(paths):
    """A live call (neither ``at`` nor ``last_used`` passed) reports this host's own clock;
    if an earlier call's ``at`` came from a host whose clock ran ahead, the stored last_used
    can already sit in the future relative to a live call made afterwards. The live call is
    still the most recent real event and must win on detail regardless, while last_used
    keeps the later (future) timestamp already on the row. Only a replay/backfill call
    (which passes ``at``) should defer to a stored, later timestamp instead."""
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_project("proj1", "Wolbachia", "/projects/wolbachia", "prepare_registry.json")
        catalog.upsert_dataset(_sidecar())

        catalog.record_usage("SRR1", "proj1", "wMel", "downloaded", detail="a", at="2099-01-01T00:00:00+00:00")
        catalog.record_usage("SRR1", "proj1", "wMel", "downloaded", detail="b")

        row = dict(
            catalog._conn.execute(
                "SELECT * FROM usage WHERE accession=? AND project_id=? AND genome_id=? AND stage=?",
                ("SRR1", "proj1", "wMel", "downloaded"),
            ).fetchone()
        )

    assert row["detail"] == "b"
    assert row["last_used"] == "2099-01-01T00:00:00+00:00"


def test_record_usage_updates_a_hand_edited_row_with_null_dates(paths):
    """first_used/last_used are nullable; SQLite's scalar min()/max() return NULL when either
    side is NULL, so a row with a NULL date (e.g. hand-edited, or migrated from an older
    schema) would never advance past NULL without coalescing against the incoming value.
    record_usage must still set both dates and refresh detail for such a row."""
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_project("proj1", "Wolbachia", "/projects/wolbachia", "prepare_registry.json")
        catalog.upsert_dataset(_sidecar())
        catalog._conn.execute(
            "INSERT INTO usage (accession, project_id, genome_id, stage, first_used, last_used, detail) "
            "VALUES (?, ?, ?, ?, NULL, NULL, ?)",
            ("SRR1", "proj1", "wMel", "downloaded", "stale"),
        )

        catalog.record_usage("SRR1", "proj1", "wMel", "downloaded", detail="fresh", at="2026-09-05T00:00:00+00:00")

        row = dict(
            catalog._conn.execute(
                "SELECT * FROM usage WHERE accession=? AND project_id=? AND genome_id=? AND stage=?",
                ("SRR1", "proj1", "wMel", "downloaded"),
            ).fetchone()
        )

    assert row["first_used"] == "2026-09-05T00:00:00+00:00"
    assert row["last_used"] == "2026-09-05T00:00:00+00:00"
    assert row["detail"] == "fresh"


def test_record_usage_inserts_placeholder_dataset_when_missing(paths):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_project("proj1", "Wolbachia", "/projects/wolbachia", "prepare_registry.json")

        catalog.record_usage("SRR-unseen", "proj1", "", "linked")

        row = catalog.get_dataset("SRR-unseen")

    assert row is not None
    assert row["state"] == "unknown"


def test_upsert_project_keeps_created_updates_last_seen(paths):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_project("proj1", "Wolbachia", "/projects/wolbachia", "prepare_registry.json")
        created_first = catalog._conn.execute("SELECT created FROM projects WHERE project_id=?", ("proj1",)).fetchone()[
            "created"
        ]

        time.sleep(0.01)
        catalog.upsert_project("proj1", "Wolbachia renamed", "/projects/wolbachia2", "prepare_registry.json")
        row = dict(catalog._conn.execute("SELECT * FROM projects WHERE project_id=?", ("proj1",)).fetchone())

    assert row["created"] == created_first
    assert row["name"] == "Wolbachia renamed"
    assert row["path"] == "/projects/wolbachia2"


def test_projects_for_and_datasets_for_project(paths):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_project("proj1", "Wolbachia", "/p1", "reg1")
        catalog.upsert_project("proj2", "Anopheles", "/p2", "reg2")
        catalog.upsert_dataset(_sidecar("SRR1"))
        catalog.upsert_dataset(_sidecar("SRR2"))

        catalog.record_usage("SRR1", "proj1", "wMel", "downloaded")
        catalog.record_usage("SRR1", "proj2", "wMel", "downloaded")
        catalog.record_usage("SRR2", "proj1", "wMel", "downloaded")

        projects_for_srr1 = catalog.projects_for("SRR1")
        datasets_for_proj1 = catalog.datasets_for_project("proj1")

    assert [p["project_id"] for p in projects_for_srr1] == ["proj1", "proj2"]
    assert [d["accession"] for d in datasets_for_proj1] == ["SRR1", "SRR2"]


def test_unused_lists_datasets_without_usage(paths):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_project("proj1", "Wolbachia", "/p1", "reg1")
        catalog.upsert_dataset(_sidecar("SRR1"))
        catalog.upsert_dataset(_sidecar("SRR2"))
        catalog.record_usage("SRR1", "proj1", "wMel", "downloaded")

        unused = catalog.unused()

    assert unused == ["SRR2"]


def test_bytes_by_genome_sums_distinct_accession_genome_pairs(paths):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_project("proj1", "P1", "/p1", "reg1")
        catalog.upsert_project("proj2", "P2", "/p2", "reg2")
        catalog.upsert_dataset(_sidecar("SRR1", bytes_r1=1000, bytes_r2=1000))  # bytes_total 2000
        catalog.upsert_dataset(_sidecar("SRR2", bytes_r1=500, bytes_r2=500))  # bytes_total 1000

        # Same (accession, genome_id) pair used in two projects: must count once.
        catalog.record_usage("SRR1", "proj1", "wMel", "downloaded")
        catalog.record_usage("SRR1", "proj2", "wMel", "downloaded")
        # Different genome for SRR2.
        catalog.record_usage("SRR2", "proj1", "wAlbB", "downloaded")

        totals = catalog.bytes_by_genome()

    assert totals == {"wMel": 2000, "wAlbB": 1000}


def test_datasets_for_genome_returns_accession_project_tuples(paths):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_project("proj1", "P1", "/p1", "reg1")
        catalog.upsert_project("proj2", "P2", "/p2", "reg2")
        catalog.upsert_dataset(_sidecar("SRR1"))
        catalog.upsert_dataset(_sidecar("SRR2"))

        catalog.record_usage("SRR1", "proj1", "wMel", "analysed")
        catalog.record_usage("SRR2", "proj2", "wMel", "analysed")
        catalog.record_usage("SRR1", "proj1", "wAlbB", "analysed")

        result = catalog.datasets_for_genome("wMel")

    assert result == [("SRR1", "proj1"), ("SRR2", "proj2")]


def test_reindex_rebuilds_datasets_and_keeps_usage(paths):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_project("proj1", "P1", "/p1", "reg1")
        catalog.upsert_dataset(_sidecar("SRR1", state="partial"))
        catalog.upsert_dataset(_sidecar("SRR2", state="complete"))
        catalog.record_usage("SRR1", "proj1", "wMel", "downloaded")
        catalog.record_usage("SRR2", "proj1", "wMel", "downloaded")

        count = catalog.reindex([_sidecar("SRR1", state="complete")])

        srr1 = catalog.get_dataset("SRR1")
        srr2 = catalog.get_dataset("SRR2")
        remaining_usage = catalog._conn.execute("SELECT accession FROM usage ORDER BY accession").fetchall()
        projects_still_there = catalog._conn.execute("SELECT project_id FROM projects").fetchall()

    assert count == 1
    assert srr1["state"] == "complete"
    assert srr2 is None
    assert [r["accession"] for r in remaining_usage] == ["SRR1"]
    assert [p["project_id"] for p in projects_still_there] == ["proj1"]


def test_journal_mode_is_delete_never_wal(paths):
    """SQLite documents WAL as unsafe over NFS and SMB, and the pragma succeeds there anyway,
    so a store on a NAS would run in an unsupported mode without saying so. Writes are already
    serialised by catalog.sqlite.lock, so DELETE costs nothing here."""
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        mode = catalog.conn.execute("PRAGMA journal_mode").fetchone()[0]

    assert str(mode).lower() == "delete"


def test_catalog_write_commits_and_releases_lock(paths):
    with catalog_write(paths) as catalog:
        catalog.upsert_dataset(_sidecar())

    assert not paths.catalog_lock.exists()

    with Catalog(paths, create=True) as catalog:
        row = catalog.get_dataset("SRR1")
    assert row is not None


def test_catalog_write_rolls_back_and_releases_lock_on_exception(paths):
    with pytest.raises(RuntimeError):
        with catalog_write(paths) as catalog:
            catalog.upsert_dataset(_sidecar())
            raise RuntimeError("boom")

    assert not paths.catalog_lock.exists()

    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        row = catalog.get_dataset("SRR1")
    assert row is None


def test_catalog_write_serializes_two_writers(paths):
    order = []
    barrier_entered = threading.Event()

    def hold_lock():
        with catalog_write(paths) as catalog:
            barrier_entered.set()
            time.sleep(0.3)
            catalog.upsert_dataset(_sidecar("SRR1"))
            order.append("thread")

    thread = threading.Thread(target=hold_lock)
    thread.start()
    barrier_entered.wait(timeout=2)

    with catalog_write(paths) as catalog:
        order.append("main")
        catalog.upsert_dataset(_sidecar("SRR2"))

    thread.join(timeout=5)

    assert order == ["thread", "main"]


def test_catalog_write_raises_when_lock_never_released(paths, monkeypatch):
    monkeypatch.setattr("prepare.store.catalog.CATALOG_LOCK_WAIT_SECONDS", 0.2)

    # Simulate another process holding the lock: a fresh, non-stale lock file.
    paths.catalog_lock.parent.mkdir(parents=True, exist_ok=True)
    paths.catalog_lock.write_text("999999")

    with pytest.raises(DataAccessError):
        with catalog_write(paths):
            pass


def test_record_usage_unknown_project_raises(paths):
    """record_usage never fabricates a project row for an id it does not recognise."""
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        with pytest.raises(DataAccessError, match="Unknown project_id"):
            catalog.record_usage("SRR1", "no-such-project", "wMel", "downloaded")

        # No placeholder dataset or usage row was left behind by the failed call.
        assert catalog.get_dataset("SRR1") is None
        assert catalog._conn.execute("SELECT COUNT(*) AS n FROM usage").fetchone()["n"] == 0


class _FlagConnection(sqlite3.Connection):
    """A real sqlite3.Connection subclass whose statement methods can be told to fail.

    sqlite3.Connection is an immutable C type: its bound methods cannot be
    monkeypatched on an instance. Subclassing it (as ``test_wal_fallback_logs_and_continues``
    above already does) and flipping a plain instance attribute is the way to make a real,
    already-open connection start raising ``sqlite3.Error`` on demand, after any setup calls
    that must succeed have already run.
    """

    _boom = False

    def execute(self, sql, *args, **kwargs):
        if self._boom:
            raise sqlite3.OperationalError("simulated disk I/O error")
        return super().execute(sql, *args, **kwargs)

    def executemany(self, sql, *args, **kwargs):
        if self._boom:
            raise sqlite3.OperationalError("simulated disk I/O error")
        return super().executemany(sql, *args, **kwargs)

    def commit(self):
        if self._boom:
            raise sqlite3.OperationalError("simulated disk I/O error")
        return super().commit()


@pytest.fixture
def flag_connect(monkeypatch):
    """Route sqlite3.connect through _FlagConnection for the duration of one test."""
    real_connect = sqlite3.connect

    def fake_connect(database, *args, **kwargs):
        kwargs["factory"] = _FlagConnection
        return real_connect(database, *args, **kwargs)

    monkeypatch.setattr(sqlite3, "connect", fake_connect)


@pytest.mark.parametrize(
    "method_name, call_args",
    [
        ("migrate", ()),
        ("upsert_dataset", (_sidecar(),)),
        ("get_dataset", ("SRR1",)),
        ("upsert_project", ("proj1", "P1", "/p1", "reg1")),
        ("projects_for", ("SRR1",)),
        ("datasets_for_project", ("proj1",)),
        ("unused", ()),
        ("bytes_by_genome", ()),
        ("datasets_for_genome", ("wMel",)),
        ("reindex", ([],)),
    ],
)
def test_public_methods_wrap_sqlite_errors(paths, flag_connect, method_name, call_args):
    """A raw sqlite3.Error out of the connection surfaces as DataAccessError, never raw."""
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.conn._boom = True

        method = getattr(catalog, method_name)
        with pytest.raises(DataAccessError):
            method(*call_args)


def test_record_usage_wraps_sqlite_error(paths, flag_connect):
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        catalog.upsert_project("proj1", "P1", "/p1", "reg1")
        catalog.conn._boom = True

        with pytest.raises(DataAccessError):
            catalog.record_usage("SRR1", "proj1", "wMel", "downloaded")


def test_catalog_open_wraps_sqlite_error(paths, monkeypatch):
    def boom(*_args, **_kwargs):
        raise sqlite3.OperationalError("unable to open database file")

    monkeypatch.setattr(sqlite3, "connect", boom)
    with pytest.raises(DataAccessError):
        with Catalog(paths):
            pass


def test_catalog_write_wraps_commit_sqlite_error(paths, flag_connect):
    with pytest.raises(DataAccessError):
        with catalog_write(paths) as catalog:
            catalog.upsert_dataset(_sidecar())
            catalog.conn._boom = True


def test_reindex_keeps_a_row_whose_folder_is_still_in_the_store(paths):
    """A sidecar missed by one reading must not take the dataset's usage history with it.

    The cascade from ``datasets`` removes ``files`` and ``usage`` too, and usage is the one
    record that exists nowhere but the catalogue.
    """
    (paths.sra / "SRR2").mkdir(parents=True)
    (paths.sra / "SRR2" / "SRR2.fastq.gz").write_bytes(b"x")

    with catalog_write(paths) as catalog:
        catalog.upsert_project("p1", "P1", "/p1", "reg1")
        catalog.upsert_dataset(_sidecar("SRR1"))
        catalog.upsert_dataset(_sidecar("SRR2"))
        catalog.record_usage("SRR2", "p1", "wMel", "downloaded")

    with catalog_write(paths) as catalog:
        catalog.reindex([_sidecar("SRR1")])
        accessions = [
            row["accession"] for row in catalog.conn.execute("SELECT accession FROM datasets ORDER BY accession")
        ]
        usage = catalog.conn.execute("SELECT COUNT(*) AS n FROM usage").fetchone()["n"]

    assert accessions == ["SRR1", "SRR2"]
    assert usage == 1


def test_reindex_drops_a_row_whose_folder_is_gone(paths):
    with catalog_write(paths) as catalog:
        catalog.upsert_dataset(_sidecar("SRR1"))
        catalog.upsert_dataset(_sidecar("SRR2"))

    with catalog_write(paths) as catalog:
        catalog.reindex([_sidecar("SRR1")])
        accessions = [
            row["accession"] for row in catalog.conn.execute("SELECT accession FROM datasets ORDER BY accession")
        ]

    assert accessions == ["SRR1"]


def test_catalog_read_waits_out_a_brief_exclusive_lock(paths):
    """Every catalogue connection's own SQLite busy handler, not just catalog_write's lock
    file, must let a reader wait out a short-lived exclusive lock from another connection
    instead of failing at once with 'database is locked'."""
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()

    entered = threading.Event()

    def hold_exclusive():
        conn = sqlite3.connect(str(paths.catalog), timeout=5.0)
        try:
            conn.execute("BEGIN EXCLUSIVE")
            entered.set()
            time.sleep(1.0)
            conn.rollback()
        finally:
            conn.close()

    thread = threading.Thread(target=hold_exclusive)
    thread.start()
    assert entered.wait(timeout=2)

    with Catalog(paths) as catalog:
        row = catalog.get_dataset("SRR1")

    thread.join(timeout=5)
    assert row is None


def test_catalog_busy_timeout_constant_is_wired_to_connect(paths, monkeypatch):
    """A short CATALOG_BUSY_TIMEOUT_SECONDS makes a new connection give up quickly against a
    held exclusive lock, which only happens if the constant actually reaches
    sqlite3.connect's ``timeout`` rather than a hard-coded default."""
    monkeypatch.setattr("prepare.store.catalog.CATALOG_BUSY_TIMEOUT_SECONDS", 0.05)

    with Catalog(paths, create=True) as catalog:
        catalog.migrate()

    entered = threading.Event()
    release = threading.Event()

    def hold_exclusive():
        conn = sqlite3.connect(str(paths.catalog), timeout=5.0)
        try:
            conn.execute("BEGIN EXCLUSIVE")
            entered.set()
            release.wait(timeout=5)
            conn.rollback()
        finally:
            conn.close()

    thread = threading.Thread(target=hold_exclusive)
    thread.start()
    assert entered.wait(timeout=2)

    try:
        with pytest.raises(DataAccessError):
            with Catalog(paths) as catalog:
                catalog.get_dataset("SRR1")
    finally:
        release.set()
        thread.join(timeout=5)


def test_a_catalogue_lock_taken_over_during_a_write_commits_nothing(paths):
    import json

    from prepare.utils.lockfile import LockLost

    other = {"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00", "token": "feedbeef"}
    with pytest.raises(LockLost):
        with catalog_write(paths) as catalog:
            catalog.upsert_dataset(_sidecar("SRR1"))
            # Another process reclaimed the lock while this holder was stalled.
            paths.catalog_lock.write_text(json.dumps(other))
    with Catalog(paths, create=True) as catalog:
        catalog.migrate()
        assert catalog.get_dataset("SRR1") is None
    assert json.loads(paths.catalog_lock.read_text())["token"] == "feedbeef"
