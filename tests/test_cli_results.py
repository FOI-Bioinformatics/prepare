"""Tests for the results_table CLI command."""

import argparse
import json
import logging

from prepare.cli.commands.results import ResultsTableCommand
from prepare.data.registry import load_registry, record_extraction, record_screening, save_registry
from prepare.processing.results import RESULTS_COLUMNS


def _args(tmp_path, **kwargs):
    base = dict(
        output=str(tmp_path / "results.tsv"),
        genome_id=None,
        parsed_containment=str(tmp_path / "parsed_containment.txt"),
        min_containment=0.0,
        registry=str(tmp_path / "prepare_registry.json"),
        no_record=False,
    )
    base.update(kwargs)
    return argparse.Namespace(**base)


def _seed(tmp_path):
    registry = load_registry(tmp_path / "prepare_registry.json")
    record_screening(registry, "SRR1", "GCF_A", 0.9, None, "matches", 0.0, None)
    record_extraction(registry, "SRR1", "GCF_A", [], 10, False, {})
    save_registry(registry)
    return tmp_path / "prepare_registry.json"


def test_writes_header_and_rows(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _seed(tmp_path)
    (tmp_path / "parsed_containment.txt").write_text("\tGCF_A\tmax_containment\nSRR1\t0.95\t0.95\nSRR2\t0.4\t0.4\n")
    rc = ResultsTableCommand().execute(_args(tmp_path))
    assert rc == 0
    lines = (tmp_path / "results.tsv").read_text().splitlines()
    assert lines[0].split("\t") == RESULTS_COLUMNS
    assert [line.split("\t")[:3] for line in lines[1:]] == [["SRR1", "GCF_A", "0.95"], ["SRR2", "GCF_A", "0.4"]]


def test_missing_parsed_table_falls_back_to_registry(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    _seed(tmp_path)
    with caplog.at_level(logging.INFO):
        rc = ResultsTableCommand().execute(_args(tmp_path))
    assert rc == 0
    assert any(r.levelno == logging.INFO and "parsed_containment.txt" in r.message for r in caplog.records)
    lines = (tmp_path / "results.tsv").read_text().splitlines()
    assert len(lines) == 2 and lines[1].startswith("SRR1\tGCF_A\t0.9\t")


def test_export_is_recorded(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    registry_path = _seed(tmp_path)
    rc = ResultsTableCommand().execute(_args(tmp_path, output="results.tsv", genome_id="GCF_A"))
    assert rc == 0
    export = json.loads(registry_path.read_text())["project"]["exports"]["results_table"]
    assert export["output"] == "results.tsv"
    assert export["date"]
    assert export["summary"]["rows"] == 1
    assert export["summary"]["accessions"] == 1 and export["summary"]["genomes"] == 1
    assert export["summary"]["genome_id"] == "GCF_A"


def test_no_record_leaves_registry_untouched(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    registry_path = _seed(tmp_path)
    before = registry_path.read_text()
    rc = ResultsTableCommand().execute(_args(tmp_path, no_record=True))
    assert rc == 0
    assert (tmp_path / "results.tsv").exists()
    assert registry_path.read_text() == before


def test_zero_rows_warns_and_writes_header_only(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    with caplog.at_level(logging.INFO):
        rc = ResultsTableCommand().execute(_args(tmp_path))
    assert rc == 0
    assert (tmp_path / "results.tsv").read_text().splitlines() == ["\t".join(RESULTS_COLUMNS)]
    assert any(r.levelno == logging.WARNING and "0 row(s)" in r.message for r in caplog.records)


def test_no_registry_is_not_created_and_the_table_is_still_written(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "parsed_containment.txt").write_text("\tGCF_A\tmax_containment\nSRR1\t0.95\t0.95\n")
    registry_path = tmp_path / "prepare_registry.json"
    assert not registry_path.exists()
    with caplog.at_level(logging.INFO):
        rc = ResultsTableCommand().execute(_args(tmp_path))
    assert rc == 0
    assert not registry_path.exists()
    assert not (tmp_path / "prepare_registry.json.lock").exists()
    assert (tmp_path / "results.tsv").read_text().splitlines()[1].startswith("SRR1\tGCF_A\t0.95\t")
    assert f"no project registry at {registry_path}; export not recorded" in caplog.text


def test_min_containment_limits_rows_and_is_recorded(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    registry_path = _seed(tmp_path)
    (tmp_path / "parsed_containment.txt").write_text("\tGCF_A\tmax_containment\nSRR1\t0.95\t0.95\nSRR2\t0.4\t0.4\n")
    rc = ResultsTableCommand().execute(_args(tmp_path, min_containment=0.5))
    assert rc == 0
    lines = (tmp_path / "results.tsv").read_text().splitlines()
    assert [line.split("\t")[:2] for line in lines[1:]] == [["SRR1", "GCF_A"]]
    summary = json.loads(registry_path.read_text())["project"]["exports"]["results_table"]["summary"]
    assert summary["min_containment"] == 0.5 and summary["rows"] == 1


def test_unparseable_parsed_table_returns_1(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    registry_path = _seed(tmp_path)
    before = registry_path.read_text()
    (tmp_path / "parsed_containment.txt").write_text("\tGCF_A\nSRR1\t0.9\nSRR2\t0.1\t0.2\t0.3\t0.4\n")
    with caplog.at_level(logging.ERROR):
        rc = ResultsTableCommand().execute(_args(tmp_path))
    assert rc == 1
    assert "Error writing the results table" in caplog.text
    assert not (tmp_path / "results.tsv").exists()
    assert registry_path.read_text() == before


def test_command_is_registered():
    from prepare.cli.main import create_parser, register_all_commands

    register_all_commands()
    parser = create_parser()
    action = next(a for a in parser._subparsers._group_actions if getattr(a, "choices", None))
    assert "results_table" in action.choices
    args = action.choices["results_table"].parse_args([])
    assert args.output == "results.tsv" and args.min_containment == 0.0 and args.no_record is False
    assert args.parsed_containment == "parsed_containment.txt"
