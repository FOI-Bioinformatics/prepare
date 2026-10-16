import json
import sqlite3

import pytest

from prepare.core.exceptions import DataAccessError
from prepare.store.catalog import Catalog, catalog_write
from prepare.store.layout import init_store
from prepare.store import journal


def test_upsert_project_and_usage_are_journaled(tmp_path):
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), str(tmp_path / "proj" / "r.json"))
        c.record_usage("SRR1", "pid1", "GCF_1", "linked", "test")
    projects = [json.loads(line) for line in (paths.journal / "projects.jsonl").read_text().splitlines()]
    usage = [json.loads(line) for line in (paths.journal / "usage.jsonl").read_text().splitlines()]
    assert projects[0]["project_id"] == "pid1" and projects[0]["name"] == "proj"
    assert usage[0]["accession"] == "SRR1" and usage[0]["stage"] == "linked"


def test_replay_restores_projects_and_usage_into_a_fresh_catalog(tmp_path):
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
        c.record_usage("SRR1", "pid1", "", "downloaded", "")
    (paths.root / "catalog.sqlite").unlink()
    with catalog_write(paths) as c:
        restored = journal.replay(paths, c)
        assert restored == (1, 1)
        assert c.conn.execute("SELECT COUNT(*) FROM projects").fetchone()[0] == 1
        assert c.conn.execute("SELECT COUNT(*) FROM usage").fetchone()[0] == 1


def test_replaying_the_journal_twice_is_idempotent(tmp_path):
    """A second replay (e.g. a repeated ``store_reindex``) must not duplicate rows: the
    catalogue's own upsert semantics (ON CONFLICT DO UPDATE for projects, DO UPDATE for the
    unique (accession, project_id, genome_id, stage) usage key) make replaying the same
    journal lines again a no-op on row counts, whichever catalogue it is replayed into."""
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
        c.record_usage("SRR1", "pid1", "GCF_1", "downloaded", "first pass")
    (paths.root / "catalog.sqlite").unlink()

    with catalog_write(paths) as c:
        first = journal.replay(paths, c)
        second = journal.replay(paths, c)

        assert first == (1, 1)
        assert second == (1, 1)
        assert c.conn.execute("SELECT COUNT(*) FROM projects").fetchone()[0] == 1
        assert c.conn.execute("SELECT COUNT(*) FROM usage").fetchone()[0] == 1
        project_row = c.conn.execute("SELECT name FROM projects WHERE project_id = ?", ("pid1",)).fetchone()
        assert project_row["name"] == "proj"


def test_replay_without_journal_returns_zero(tmp_path):
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        assert journal.replay(paths, c) == (0, 0)


def test_replay_skips_usage_whose_project_was_not_restored(tmp_path, caplog):
    """A usage line for a project whose own line is missing must not abort the whole replay."""
    import logging

    paths = init_store(tmp_path / "store")
    journal.append_usage(paths, "SRR1", "pid-missing", "", "linked", "")

    with caplog.at_level(logging.WARNING):
        with catalog_write(paths) as c:
            restored = journal.replay(paths, c)
            assert restored == (0, 0)
            assert c.conn.execute("SELECT COUNT(*) FROM usage").fetchone()[0] == 0
    assert any("Skipped 1 usage record" in r.message for r in caplog.records)


def test_replay_skips_malformed_journal_lines_without_raising(tmp_path):
    """A line that parses but is not a usable record (not an object, or missing a required
    key) is skipped rather than raised through KeyError/TypeError."""
    paths = init_store(tmp_path / "store")
    paths.journal.mkdir(parents=True, exist_ok=True)
    (paths.journal / "projects.jsonl").write_text('[1, 2, 3]\n{"name": "no id here"}\n')
    (paths.journal / "usage.jsonl").write_text('"just a string"\n{"accession": "SRR1"}\n')

    with catalog_write(paths) as c:
        assert journal.replay(paths, c) == (0, 0)


def test_replay_out_of_order_keeps_latest_last_used_and_earliest_first_used(tmp_path):
    """Journal lines are replayed in file order, which is not necessarily chronological (e.g.
    a line appended late by a slow writer). The restored row must still end up with first_used
    at the earliest ``at`` and last_used at the latest ``at``, not just the last line replayed."""
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
    journal.append_usage(paths, "SRR1", "pid1", "", "linked", "", at="2026-09-02T00:00:00+00:00")
    journal.append_usage(paths, "SRR1", "pid1", "", "linked", "", at="2026-09-01T00:00:00+00:00")

    with catalog_write(paths) as c:
        journal.replay(paths, c)
        row = c.conn.execute("SELECT first_used, last_used FROM usage WHERE accession='SRR1'").fetchone()

    assert row["first_used"] == "2026-09-01T00:00:00+00:00"
    assert row["last_used"] == "2026-09-02T00:00:00+00:00"


def test_replay_restores_an_explicit_last_used_from_a_single_line(tmp_path):
    """A single journal line can itself carry a ``last_used`` later than its own ``at``
    (e.g. a backfilled row touched more than once before the journal existed). Replay must
    pass that ``last_used`` through to ``record_usage`` rather than collapsing both bounds to
    ``at``."""
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
    journal.append_usage(
        paths,
        "SRR1",
        "pid1",
        "",
        "linked",
        "",
        at="2026-09-01T00:00:00+00:00",
        last_used="2026-09-03T00:00:00+00:00",
    )

    with catalog_write(paths) as c:
        journal.replay(paths, c)
        row = c.conn.execute("SELECT first_used, last_used FROM usage WHERE accession='SRR1'").fetchone()

    assert row["first_used"] == "2026-09-01T00:00:00+00:00"
    assert row["last_used"] == "2026-09-03T00:00:00+00:00"


def test_replay_counts_distinct_projects_not_journal_lines(tmp_path):
    """``upsert_project`` appends one journal line per call, even when it is the same project
    updated again (e.g. a ``store_link`` run touching the same project's ``last_seen``
    repeatedly). Replay's project count must reflect distinct projects restored, not the
    number of lines replayed, or ``store_reindex`` misreports how many projects came back."""
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
    (paths.root / "catalog.sqlite").unlink()

    with catalog_write(paths) as c:
        projects, _usage = journal.replay(paths, c)

    assert projects == 1


def test_catalog_write_backfills_a_pre_journal_store_once(tmp_path):
    """A pre-journal store's ``projects``/``usage`` rows are copied into the journal the first
    time ``catalog_write`` opens it afterwards (``catalog_write`` calls ``backfill_from_catalog``
    itself, right after ``migrate()``), and never duplicated on a later open."""
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        c.journal_enabled = False  # simulate a store written before the journal existed
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
        c.record_usage("SRR1", "pid1", "", "linked", "")
    assert not (paths.journal / "projects.jsonl").exists()

    with catalog_write(paths):
        pass  # backfill runs on entry, before this block's body

    projects = [json.loads(line) for line in (paths.journal / "projects.jsonl").read_text().splitlines()]
    usage = [json.loads(line) for line in (paths.journal / "usage.jsonl").read_text().splitlines()]
    assert [p["project_id"] for p in projects] == ["pid1"]
    assert projects[0]["name"] == "proj" and projects[0]["registry"] == "r.json"
    assert [line["accession"] for line in usage] == ["SRR1"]

    with catalog_write(paths):
        pass  # a second open must not duplicate the lines already backfilled
    assert len((paths.journal / "projects.jsonl").read_text().splitlines()) == 1
    assert len((paths.journal / "usage.jsonl").read_text().splitlines()) == 1


def test_backfill_keeps_the_catalogue_rows_hostname_and_time(tmp_path):
    """Backfilled lines carry the host and time the catalogue row recorded, not the host and
    time of the backfill, so a later replay restores the true values."""
    paths = init_store(tmp_path / "store")
    when = "2020-01-02T03:04:05+00:00"
    with catalog_write(paths) as c:
        c.journal_enabled = False  # simulate a store written before the journal existed
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json", hostname="other-host")
        c.record_usage("SRR1", "pid1", "", "linked", "", at=when)

    with catalog_write(paths):
        pass  # backfill runs on entry

    project = json.loads((paths.journal / "projects.jsonl").read_text().splitlines()[0])
    usage = json.loads((paths.journal / "usage.jsonl").read_text().splitlines()[0])
    assert project["hostname"] == "other-host"
    assert usage["at"] == when


def test_backfill_line_carries_first_and_last_used(tmp_path):
    """A pre-journal usage row can have been touched more than once before the journal
    existed, so first_used and last_used differ; the backfilled line must carry both, not
    collapse to a single timestamp (which would make a later replay lose the earlier date)."""
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        c.journal_enabled = False  # simulate a store written before the journal existed
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
        c.record_usage("SRR1", "pid1", "", "linked", "", at="2026-09-01T00:00:00+00:00")
        c.record_usage("SRR1", "pid1", "", "linked", "", at="2026-09-03T00:00:00+00:00")
    assert not (paths.journal / "projects.jsonl").exists()

    with catalog_write(paths):
        pass  # backfill runs on entry

    usage = json.loads((paths.journal / "usage.jsonl").read_text().splitlines()[0])
    assert usage["at"] == "2026-09-01T00:00:00+00:00"
    assert usage["last_used"] == "2026-09-03T00:00:00+00:00"


def test_backfill_from_catalog_is_a_noop_on_an_empty_catalog(tmp_path):
    """No project rows yet (a freshly initialised store) means nothing to backfill."""
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        assert journal.backfill_from_catalog(paths, c) == (0, 0)
    assert not (paths.journal / "projects.jsonl").exists()


def test_backfill_sqlite_error_becomes_data_access_error(tmp_path, monkeypatch):
    """backfill_from_catalog's own sqlite calls must be wrapped like every other catalogue
    method: a raw sqlite3.Error must never reach the caller of catalog_write."""
    paths = init_store(tmp_path / "store")
    with Catalog(paths, create=True) as c:
        c.migrate()
        c.journal_enabled = False  # a pre-journal project: nothing appended yet
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
        c.conn.commit()
    assert not (paths.journal / "projects.jsonl").exists()

    class _RaisingConn:
        def execute(self, *args, **kwargs):
            raise sqlite3.OperationalError("boom")

    # sqlite3.Connection is a C extension type: its methods cannot be monkeypatched directly
    # (setattr on an instance or the class both raise), so the ``conn`` property itself is
    # patched to hand back a stand-in whose ``execute`` raises, same effect as if the real
    # connection had failed mid-query.
    monkeypatch.setattr(Catalog, "conn", property(lambda self: _RaisingConn()))
    with Catalog(paths) as c:
        with pytest.raises(DataAccessError):
            journal.backfill_from_catalog(paths, c)


def test_prune_usage_drops_only_the_named_accessions_and_keeps_corrupt_lines(tmp_path, caplog):
    import logging

    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
        c.record_usage("SRR1", "pid1", "", "downloaded", "")
        c.record_usage("SRR2", "pid1", "", "downloaded", "")
        c.record_usage("SRR3", "pid1", "", "downloaded", "")
    usage_file = paths.journal / "usage.jsonl"
    with open(usage_file, "a", encoding="utf-8") as handle:
        handle.write("{broken\n[1, 2]\n")

    with caplog.at_level(logging.WARNING):
        dropped = journal.prune_usage(paths, {"SRR1", "SRR3"})

    assert dropped == 2
    lines = usage_file.read_text(encoding="utf-8").splitlines()
    assert lines[-2:] == ["{broken", "[1, 2]"]
    assert [json.loads(line)["accession"] for line in lines[:-2]] == ["SRR2"]
    assert sum("could not be parsed" in record.message for record in caplog.records) == 1


def test_prune_usage_leaves_the_file_untouched_when_nothing_matches(tmp_path):
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
        c.record_usage("SRR2", "pid1", "", "downloaded", "")
    usage_file = paths.journal / "usage.jsonl"
    before = usage_file.stat().st_mtime_ns, usage_file.read_text(encoding="utf-8")

    assert journal.prune_usage(paths, {"SRR1"}) == 0
    assert (usage_file.stat().st_mtime_ns, usage_file.read_text(encoding="utf-8")) == before


def test_prune_usage_without_a_journal_is_a_noop(tmp_path):
    paths = init_store(tmp_path / "store")
    assert journal.prune_usage(paths, {"SRR1"}) == 0
    assert not (paths.journal / "usage.jsonl").exists()


def test_compact_projects_keeps_the_latest_line_per_project_and_corrupt_lines(tmp_path):
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        c.upsert_project("p1", "first", "/a", "r1.json")
        c.upsert_project("p2", "other", "/b", "r2.json")
        c.upsert_project("p1", "second", "/a", "r1.json")
    projects_file = paths.journal / "projects.jsonl"
    with open(projects_file, "a", encoding="utf-8") as handle:
        handle.write('{oops\n{"name": "no id"}\n')

    assert journal.compact_projects(paths) == 1

    lines = projects_file.read_text(encoding="utf-8").splitlines()
    assert json.loads(lines[0])["project_id"] == "p2"
    assert json.loads(lines[1])["name"] == "second"
    assert lines[2:] == ["{oops", '{"name": "no id"}']


def test_prune_usage_writes_kept_lines_back_byte_identical(tmp_path):
    """Kept lines are written exactly as read: no whitespace stripping, and no split on the
    other separators ``str.splitlines`` honours (U+2028, U+0085) inside a foreign line."""
    paths = init_store(tmp_path / "store")
    paths.journal.mkdir(parents=True, exist_ok=True)
    keep_record = '{"accession": "SRR2", "project_id": "p1", "stage": "downloaded"}  '
    corrupt = "{broken \u2028 still the same line \x85 and more   "
    drop_record = '{"accession": "SRR1", "project_id": "p1", "stage": "downloaded"}'
    usage_file = paths.journal / "usage.jsonl"
    with open(usage_file, "w", encoding="utf-8", newline="") as handle:
        handle.write(keep_record + "\n" + corrupt + "\n" + drop_record + "\n")

    assert journal.prune_usage(paths, {"SRR1"}) == 1

    expected = (keep_record + "\n" + corrupt + "\n").encode("utf-8")
    assert usage_file.read_bytes() == expected


def test_journal_rewrites_are_fsynced(tmp_path, monkeypatch):
    paths = init_store(tmp_path / "store")
    with catalog_write(paths) as c:
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
        c.upsert_project("pid1", "proj", str(tmp_path / "proj"), "r.json")
        c.record_usage("SRR1", "pid1", "", "downloaded", "")
    calls = []
    real = journal.write_text_atomic

    def recording(path, text, *args, **kwargs):
        calls.append(kwargs.get("fsync"))
        return real(path, text, *args, **kwargs)

    monkeypatch.setattr(journal, "write_text_atomic", recording)

    assert journal.prune_usage(paths, {"SRR1"}) == 1
    assert journal.compact_projects(paths) == 1
    assert calls == [True, True]
