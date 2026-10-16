"""Tests for the ``runs`` command (``prepare/cli/commands/runs.py``)."""

import argparse
import json
import logging
from datetime import datetime, timedelta, timezone

import pytest

from prepare.cli.main import main
from prepare.core import settings
from prepare.data import run_log
from prepare.data.registry import REGISTRY_FILENAME

STARTED = datetime(2026, 10, 1, 12, 5, 1, tzinfo=timezone.utc)


@pytest.fixture
def project(tmp_path, monkeypatch):
    """A project folder holding a registry file, with the run log on while fixtures are written."""
    monkeypatch.setenv(settings.SETTINGS["run_log"].env, "true")
    settings.reset_for_tests()
    folder = tmp_path / "project"
    folder.mkdir()
    (folder / REGISTRY_FILENAME).write_text("{}")
    return folder


def _registry(project):
    return str(project / REGISTRY_FILENAME)


def _profile_run(project, minutes, gc_srr1, argv=None):
    """Record one fake ``sra_profile`` run over SRR1 and SRR2."""
    args = argparse.Namespace(registry=None, threads=2)
    run_log.note_run(
        args,
        summary={"profiled": 2, "failed": 0},
        detail={
            "analyses": {
                "profile": {
                    "SRR1": {"gc_percent": gc_srr1, "total_reads": 1000},
                    "SRR2": {"gc_percent": 50.0, "total_reads": 2000},
                }
            }
        },
    )
    return run_log.record_run(
        project,
        "sra_profile",
        argv or ["sra_profile", "--threads", "2"],
        args,
        STARTED + timedelta(minutes=minutes),
        3.25,
        0,
    )


@pytest.fixture
def two_profiles(project):
    """Two ``sra_profile`` runs that differ only in SRR1's gc_percent."""
    first = _profile_run(project, 0, 41.2)
    second = _profile_run(project, 5, 43.0)
    return first, second


def _json(capsys):
    out = capsys.readouterr().out
    return json.loads(out)


def test_runs_is_listed_in_environment_group(capsys):
    with pytest.raises(SystemExit):
        main(["--help"])
    out = capsys.readouterr().out
    environment = out.split("Environment:")[1]
    assert "runs" in environment


def test_list_shows_runs_newest_first(project, two_profiles, capsys):
    first, second = two_profiles
    assert main(["runs", "--registry", _registry(project)]) == 0
    out = capsys.readouterr().out
    assert out.index(second.run_id) < out.index(first.run_id)
    assert "sra_profile" in out
    assert "profiled=2" in out


def test_list_found_from_working_directory(project, two_profiles, capsys, monkeypatch):
    nested = project / "sub"
    nested.mkdir()
    monkeypatch.chdir(nested)
    assert main(["runs"]) == 0
    assert two_profiles[1].run_id in capsys.readouterr().out


def test_list_limit_and_command_filter(project, two_profiles, capsys):
    args = argparse.Namespace(registry=None)
    run_log.record_run(project, "blacklist", ["blacklist", "--list"], args, STARTED + timedelta(hours=1), 0.1, 0)
    assert main(["runs", "--registry", _registry(project), "--limit", "1", "--json"]) == 0
    document = _json(capsys)
    assert [run["command"] for run in document["runs"]] == ["blacklist"]
    assert document["total"] == 3

    assert main(["runs", "--registry", _registry(project), "--command", "sra_profile", "--json"]) == 0
    document = _json(capsys)
    assert [run["run_id"] for run in document["runs"]] == [two_profiles[1].run_id, two_profiles[0].run_id]
    assert document["total"] == 2


def test_list_limit_zero_lists_every_run(project, two_profiles, capsys):
    assert main(["runs", "--registry", _registry(project), "--limit", "0", "--json"]) == 0
    assert len(_json(capsys)["runs"]) == 2


def test_list_default_limit_is_twenty_newest_runs(project, capsys):
    args = argparse.Namespace(registry=None)
    records = [
        run_log.record_run(project, "blacklist", ["blacklist", "--list"], args, STARTED + timedelta(minutes=i), 0.1, 0)
        for i in range(25)
    ]
    assert main(["runs", "--registry", _registry(project), "--json"]) == 0
    document = _json(capsys)
    assert document["total"] == 25
    assert [run["run_id"] for run in document["runs"]] == [record.run_id for record in reversed(records[5:])]


def test_list_json_is_one_document(project, two_profiles, capsys):
    assert main(["runs", "--registry", _registry(project), "--json"]) == 0
    document = _json(capsys)
    assert document["project"] == str(project.resolve())
    assert document["runs"][0]["detail_kept"] is True
    assert document["runs"][0]["summary"] == {"profiled": 2, "failed": 0}


def test_list_command_without_runs(project, two_profiles, capsys):
    assert main(["runs", "--registry", _registry(project), "--command", "download_sra"]) == 0
    assert "No runs of download_sra" in capsys.readouterr().out


def test_diff_latest_previous_reports_changed_gc_percent(project, two_profiles, capsys):
    assert main(["runs", "--registry", _registry(project), "--diff", "previous", "latest"]) == 0
    out = capsys.readouterr().out
    assert two_profiles[0].run_id in out and two_profiles[1].run_id in out
    srr1 = [line for line in out.splitlines() if "SRR1" in line and "gc_percent" in line]
    assert len(srr1) == 1
    assert "41.2" in srr1[0] and "43.0" in srr1[0]
    assert "SRR2" not in out


def test_diff_json_document(project, two_profiles, capsys):
    assert main(["runs", "--registry", _registry(project), "--diff", "latest", "previous", "--json"]) == 0
    document = _json(capsys)
    assert document["run_a"]["run_id"] == two_profiles[1].run_id
    assert document["run_b"]["run_id"] == two_profiles[0].run_id
    assert document["detail"]["changed"] == {"SRR1": {"gc_percent": [43.0, 41.2]}}
    assert document["detail"]["added"] == [] and document["detail"]["removed"] == []
    assert document["detail_kept"] == {"run_a": True, "run_b": True}
    assert all(row["changed"] is False for row in document["summary"])


def test_diff_by_run_id_prefix(project, two_profiles, capsys):
    first, second = two_profiles
    prefix = first.run_id[: len("20261001T120501Z")]
    assert main(["runs", "--registry", _registry(project), "--diff", prefix, second.run_id, "--json"]) == 0
    assert _json(capsys)["run_a"]["run_id"] == first.run_id


def test_diff_with_pruned_detail_reports_not_kept(project, capsys, monkeypatch):
    monkeypatch.setattr(run_log, "DETAILS_KEPT_PER_COMMAND", 1)
    first = _profile_run(project, 0, 41.2)
    _profile_run(project, 5, 43.0)
    assert run_log.read_runs(project)[0].detail is None

    assert main(["runs", "--registry", _registry(project), "--diff", "previous", "latest"]) == 0
    out = capsys.readouterr().out
    assert f"not kept for run {first.run_id}" in out

    assert main(["runs", "--registry", _registry(project), "--diff", "previous", "latest", "--json"]) == 0
    document = _json(capsys)
    assert document["detail"] is None
    assert document["detail_kept"] == {"run_a": False, "run_b": True}


def test_show_run_with_detail(project, two_profiles, capsys):
    assert main(["runs", "--registry", _registry(project), "--show", "latest"]) == 0
    out = capsys.readouterr().out
    assert two_profiles[1].run_id in out
    assert "sra_profile --threads 2" in out
    assert "Detail: kept" in out
    assert "SRR1" in out and "43.0" in out


def test_show_kept_detail_without_rows_prints_it_as_json(project, capsys):
    args = argparse.Namespace(registry=None)
    run_log.note_run(args, summary={"removed": 0}, detail={"count": 3})
    record = run_log.record_run(project, "blacklist", ["blacklist"], args, STARTED, 0.1, 0)
    assert main(["runs", "--registry", _registry(project), "--show", record.run_id]) == 0
    out = capsys.readouterr().out
    assert "Detail: kept" in out
    assert "  {" in out and '"count": 3' in out


def test_show_pruned_detail_reported_as_not_kept(project, capsys, monkeypatch):
    monkeypatch.setattr(run_log, "DETAILS_KEPT_PER_COMMAND", 1)
    first = _profile_run(project, 0, 41.2)
    _profile_run(project, 5, 43.0)

    assert main(["runs", "--registry", _registry(project), "--show", first.run_id]) == 0
    assert "Detail: not kept" in capsys.readouterr().out

    assert main(["runs", "--registry", _registry(project), "--show", first.run_id, "--json"]) == 0
    document = _json(capsys)
    assert document["run"]["run_id"] == first.run_id
    assert document["detail_kept"] is False
    assert document["detail"] is None


def test_accession_traces_values_across_runs(project, two_profiles, capsys):
    assert main(["runs", "--registry", _registry(project), "--accession", "SRR1"]) == 0
    out = capsys.readouterr().out
    assert out.index(two_profiles[0].run_id) < out.index(two_profiles[1].run_id)
    assert "41.2" in out and "43.0" in out

    assert main(["runs", "--registry", _registry(project), "--accession", "SRR1", "--json"]) == 0
    document = _json(capsys)
    assert document["accession"] == "SRR1"
    assert [entry["values"]["SRR1"]["gc_percent"] for entry in document["runs"]] == [41.2, 43.0]
    assert document["runs_without_detail"] == 0


def test_accession_not_found(project, two_profiles, capsys):
    assert main(["runs", "--registry", _registry(project), "--accession", "SRR999"]) == 0
    assert "No recorded run" in capsys.readouterr().out


def test_accession_counts_runs_whose_detail_was_pruned(project, capsys, monkeypatch):
    monkeypatch.setattr(run_log, "DETAILS_KEPT_PER_COMMAND", 1)
    _profile_run(project, 0, 41.2)
    _profile_run(project, 5, 43.0)
    assert main(["runs", "--registry", _registry(project), "--accession", "SRR1", "--json"]) == 0
    document = _json(capsys)
    assert len(document["runs"]) == 1
    assert document["runs_without_detail"] == 1


def test_no_registry_exits_1(tmp_path, capsys, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    with caplog.at_level(logging.ERROR):
        assert main(["runs"]) == 1
    assert "no run log" in caplog.text.lower()


def test_no_run_log_exits_1_with_json_error(project, capsys):
    assert main(["runs", "--registry", _registry(project), "--json"]) == 1
    assert "error" in _json(capsys)


def test_bad_selector_is_validation_error(project, two_profiles, capsys, caplog):
    with caplog.at_level(logging.ERROR):
        assert main(["runs", "--registry", _registry(project), "--show", "nosuchrun"]) == 2
    assert "No run matches 'nosuchrun'" in caplog.text
    assert main(["runs", "--registry", _registry(project), "--show", "nosuchrun", "--json"]) == 2
    assert "No run matches" in _json(capsys)["error"]


def test_ambiguous_prefix_is_validation_error(project, two_profiles, caplog):
    with caplog.at_level(logging.ERROR):
        assert main(["runs", "--registry", _registry(project), "--show", "2026"]) == 2
    assert "matches 2 runs" in caplog.text


def test_previous_with_one_run_is_validation_error(project, caplog):
    _profile_run(project, 0, 41.2)
    with caplog.at_level(logging.ERROR):
        assert main(["runs", "--registry", _registry(project), "--diff", "previous", "latest"]) == 2
    assert "No previous run" in caplog.text


def test_show_diff_accession_are_mutually_exclusive(project, two_profiles):
    with pytest.raises(SystemExit) as excinfo:
        main(["runs", "--registry", _registry(project), "--show", "latest", "--accession", "SRR1"])
    assert excinfo.value.code == 2


def test_runs_records_nothing(project, two_profiles, capsys):
    before = run_log.read_runs(project)
    assert main(["runs", "--registry", _registry(project)]) == 0
    assert run_log.read_runs(project) == before
