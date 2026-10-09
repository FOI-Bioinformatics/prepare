"""Tests for comparing recorded runs (``prepare/processing/run_diff.py``)."""

import argparse
from datetime import datetime, timedelta, timezone

import pytest

from prepare.core import settings
from prepare.data import run_log
from prepare.data.registry import REGISTRY_FILENAME
from prepare.processing.run_diff import (
    accession_history,
    detail_kept,
    detail_rows,
    diff_details,
    diff_summaries,
)

STARTED = datetime(2026, 10, 1, 12, 5, 1, tzinfo=timezone.utc)


@pytest.fixture
def project(tmp_path, monkeypatch):
    """A project folder holding a registry file, with the run log on."""
    monkeypatch.setenv(settings.SETTINGS["run_log"].env, "true")
    settings.reset_for_tests()
    folder = tmp_path / "project"
    folder.mkdir()
    (folder / REGISTRY_FILENAME).write_text("{}")
    return folder


def _profile_detail(rows):
    return {"analyses": {"profile": rows}}


def _record(project, command, minutes, argv=None, summary=None, detail=None):
    args = argparse.Namespace(registry=None)
    run_log.note_run(args, summary=summary, detail=detail)
    return run_log.record_run(project, command, argv or [command], args, STARTED + timedelta(minutes=minutes), 1.0, 0)


# --- diff_summaries ------------------------------------------------------------------------------


def test_diff_summaries_reports_every_key_with_numeric_delta():
    rows = diff_summaries({"profiled": 2, "failed": 0, "grade": "A"}, {"profiled": 3, "failed": 0, "grade": "B"})
    by_key = {row["key"]: row for row in rows}
    assert [row["key"] for row in rows] == ["failed", "grade", "profiled"]
    assert by_key["profiled"] == {"key": "profiled", "before": 2, "after": 3, "delta": 1, "changed": True}
    assert by_key["failed"]["changed"] is False and by_key["failed"]["delta"] == 0
    assert by_key["grade"]["delta"] is None and by_key["grade"]["changed"] is True


def test_diff_summaries_flattens_nested_values_and_handles_missing_keys():
    rows = diff_summaries({"counts": {"ok": 1}}, {"counts": {"ok": 1, "failed": 2}, "new": True})
    by_key = {row["key"]: row for row in rows}
    assert by_key["counts.ok"]["changed"] is False
    assert by_key["counts.failed"] == {
        "key": "counts.failed",
        "before": None,
        "after": 2,
        "delta": None,
        "changed": True,
    }
    # A boolean is not treated as a number.
    assert by_key["new"]["delta"] is None


def test_diff_summaries_rounds_float_delta():
    rows = diff_summaries({"gc_percent": 41.2}, {"gc_percent": 43.0})
    assert rows[0]["delta"] == pytest.approx(1.8)
    assert rows[0]["delta"] == round(43.0 - 41.2, 6)


# --- detail_rows / diff_details ------------------------------------------------------------------


def test_detail_rows_finds_rows_under_nested_sections():
    rows = detail_rows(_profile_detail({"SRR1": {"gc_percent": 41.2}, "SRR2": {"gc_percent": 50.0}}))
    assert rows == {"SRR1": {"gc_percent": 41.2}, "SRR2": {"gc_percent": 50.0}}


def test_detail_rows_prefixes_fields_of_a_key_found_in_several_sections():
    detail = {"analyses": {"profile": {"SRR1": {"gc_percent": 41.2}}, "report": {"SRR1": {"grade": "A"}}}}
    assert detail_rows(detail) == {"SRR1": {"analyses.profile.gc_percent": 41.2, "analyses.report.grade": "A"}}


def test_detail_rows_ignores_a_flat_mapping_beside_the_rows():
    # The Task 18 review's reproduction: a settings mapping beside the rows is not a row and
    # does not rename the row's fields, so a diff against the same detail without it is empty.
    with_settings = {"analyses": {"profile": {"SRR1": {"gc_percent": 41.2}}}, "settings": {"threads": 2}}
    assert detail_rows(with_settings) == {"SRR1": {"gc_percent": 41.2}}
    diff = diff_details(with_settings, _profile_detail({"SRR1": {"gc_percent": 41.2}}))
    assert diff == {"added": [], "removed": [], "changed": {}, "unchanged": 1}


def test_detail_rows_keeps_sections_with_the_same_last_name_apart():
    detail = {"x": {"profile": {"SRR1": {"a": 1}}}, "y": {"profile": {"SRR1": {"a": 2}}}}
    assert detail_rows(detail) == {"SRR1": {"x.profile.a": 1, "y.profile.a": 2}}


def test_detail_rows_does_not_prefix_keys_found_in_one_section_only():
    detail = {"failed": {"SRR1": {"reason": "network"}}, "downloaded": {"SRR2": {"seconds": 3.0}}}
    assert detail_rows(detail) == {"SRR1": {"reason": "network"}, "SRR2": {"seconds": 3.0}}


def test_detail_rows_of_flat_rows_and_of_none():
    assert detail_rows({"SRR1/genome": {"mapped_reads": 5}}) == {"SRR1/genome": {"mapped_reads": 5}}
    assert detail_rows(None) == {}
    assert detail_rows({"count": 3}) == {}


def test_diff_details_added_removed_changed():
    a = _profile_detail({"SRR1": {"gc_percent": 41.2, "total_reads": 100}, "SRR2": {"gc_percent": 50.0}})
    b = _profile_detail({"SRR1": {"gc_percent": 43.0, "total_reads": 100}, "SRR3": {"gc_percent": 39.0}})
    diff = diff_details(a, b)
    assert diff["added"] == ["SRR3"]
    assert diff["removed"] == ["SRR2"]
    assert diff["changed"] == {"SRR1": {"gc_percent": [41.2, 43.0]}}
    assert diff["unchanged"] == 0


def test_diff_details_identical_details():
    detail = _profile_detail({"SRR1": {"gc_percent": 41.2}})
    assert diff_details(detail, detail) == {"added": [], "removed": [], "changed": {}, "unchanged": 1}


def test_diff_details_reports_field_present_on_one_side_only():
    diff = diff_details({"SRR1": {"a": 1}}, {"SRR1": {"a": 1, "b": 2}})
    assert diff["changed"] == {"SRR1": {"b": [None, 2]}}


# --- detail_kept / accession_history -------------------------------------------------------------


def test_detail_kept_false_after_pruning(project, monkeypatch):
    monkeypatch.setattr(run_log, "DETAILS_KEPT_PER_COMMAND", 1)
    first = _record(project, "sra_profile", 0, detail=_profile_detail({"SRR1": {"gc_percent": 41.2}}))
    second = _record(project, "sra_profile", 1, detail=_profile_detail({"SRR1": {"gc_percent": 43.0}}))
    records = run_log.read_runs(project)
    assert [detail_kept(project, r) for r in records] == [False, True]
    assert first.run_id == records[0].run_id and second.run_id == records[1].run_id


def test_detail_kept_false_when_file_is_gone(project):
    record = _record(project, "sra_profile", 0, detail=_profile_detail({"SRR1": {"gc_percent": 41.2}}))
    (run_log.runs_dir(project) / record.detail).unlink()
    assert detail_kept(project, record) is False


def test_accession_history_lists_values_per_run(project):
    _record(project, "sra_profile", 0, detail=_profile_detail({"SRR1": {"gc_percent": 41.2}, "SRR2": {"x": 1}}))
    _record(project, "sra_profile", 1, detail=_profile_detail({"SRR2": {"x": 2}}))
    _record(project, "extract_target_reads", 2, detail={"SRR1/genomeA": {"mapped_reads": 7}})
    _record(project, "sra_profile", 3, detail=_profile_detail({"SRR1": {"gc_percent": 43.0}}))
    records = run_log.read_runs(project)

    history = accession_history(project, records, "SRR1")
    assert [entry["command"] for entry in history] == ["sra_profile", "extract_target_reads", "sra_profile"]
    assert history[0]["values"] == {"SRR1": {"gc_percent": 41.2}}
    assert history[1]["values"] == {"SRR1/genomeA": {"mapped_reads": 7}}
    assert history[2]["values"] == {"SRR1": {"gc_percent": 43.0}}
    assert all(entry["detail_kept"] for entry in history)
    assert history[0]["run_id"] == records[0].run_id


def test_accession_history_includes_pruned_run_named_in_argv(project, monkeypatch):
    monkeypatch.setattr(run_log, "DETAILS_KEPT_PER_COMMAND", 1)
    _record(project, "sra_profile", 0, argv=["sra_profile", "SRR1"], detail=_profile_detail({"SRR1": {"a": 1}}))
    _record(project, "sra_profile", 1, argv=["sra_profile", "SRR9"], detail=_profile_detail({"SRR9": {"a": 1}}))
    history = accession_history(project, run_log.read_runs(project), "SRR1")
    assert len(history) == 1
    assert history[0]["detail_kept"] is False
    assert history[0]["values"] is None


def test_accession_history_includes_pruned_run_naming_the_accession_as_option_value(project, monkeypatch):
    monkeypatch.setattr(run_log, "DETAILS_KEPT_PER_COMMAND", 1)
    _record(project, "sra_profile", 0, argv=["sra_profile", "--accession=SRR1"], detail=_profile_detail({"SRR1": {}}))
    _record(project, "sra_profile", 1, argv=["sra_profile", "--accession=SRR11"], detail=_profile_detail({"S": {}}))
    _record(project, "sra_profile", 2, argv=["sra_profile", "SRR9"], detail=_profile_detail({"SRR9": {"a": 1}}))
    history = accession_history(project, run_log.read_runs(project), "SRR1")
    assert [entry["started"] for entry in history] == [run_log.read_runs(project)[0].started]
    assert history[0]["values"] is None


def test_accession_history_does_not_match_accession_prefixes(project):
    _record(project, "sra_profile", 0, detail=_profile_detail({"SRR11": {"a": 1}}))
    assert accession_history(project, run_log.read_runs(project), "SRR1") == []
