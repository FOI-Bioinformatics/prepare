"""Tests for the ``project_report`` command (``prepare/cli/commands/project_report.py``)."""

import builtins
import gzip
import io
import json
import logging
import os
import sys

import pytest

from helpers_project_report import build_project
from prepare.cli.commands.project_report import ProjectReportCommand
from prepare.cli.main import create_parser, main
from prepare.data.registry import load_registry
from prepare.processing import project_report as pr
from prepare.processing.doctor_report import Check, run_checks

FILES = ("project_report.md", "project_report.json", "project_report.html")


@pytest.fixture
def project(tmp_path, monkeypatch):
    """A project folder with a registry, as the working directory; the environment checks are faked."""
    folder = tmp_path / "project"
    folder.mkdir()
    build_project(folder)
    monkeypatch.chdir(folder)
    monkeypatch.setattr(pr, "run_checks", lambda *a, **k: [Check("python", "ok", "Python 3.12")])
    return folder


def _run(*argv):
    args = create_parser().parse_args(["project_report", *argv])
    return ProjectReportCommand().execute(args)


def _block(monkeypatch, *modules):
    for name in modules:
        for loaded in [m for m in sys.modules if m == name or m.startswith(name + ".")]:
            monkeypatch.delitem(sys.modules, loaded)
        monkeypatch.setitem(sys.modules, name, None)


def _exports(folder):
    return load_registry(folder / "prepare_registry.json").project.get("exports", {})


def test_project_report_is_listed_in_environment_group(capsys):
    with pytest.raises(SystemExit):
        main(["--help"])
    environment = capsys.readouterr().out.split("Environment:")[1]
    assert "project_report" in environment


def test_without_a_registry_exits_1_pointing_to_status_init(tmp_path, monkeypatch, caplog):
    monkeypatch.chdir(tmp_path)
    caplog.set_level(logging.INFO)
    assert _run("--html", "never") == 1
    assert "status --init" in caplog.text
    assert not (tmp_path / "project_report").exists()
    assert not (tmp_path / "prepare_registry.json").exists()


def test_writes_markdown_and_json_and_records_the_export(project):
    assert _run("--html", "never") == 0
    out = project / "project_report"
    assert (out / "project_report.md").is_file()
    assert not (out / "project_report.html").exists()
    report = json.loads((out / "project_report.json").read_text())
    assert [key for key in report if key in pr.SECTIONS] == list(pr.SECTIONS)
    assert report["environment"]["included"] is True
    export = _exports(project)["project_report"]
    assert export["output"] == "project_report/project_report.md"
    assert export["summary"]["files"] == ["project_report.md", "project_report.json"]
    assert export["summary"]["html"] is False


def test_no_record_leaves_the_registry_alone(project):
    before = (project / "prepare_registry.json").read_text()
    assert _run("--html", "never", "--no-record") == 0
    assert (project / "prepare_registry.json").read_text() == before


def test_output_dir_and_registry_flags(project, tmp_path, monkeypatch):
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    registry = str(project / "prepare_registry.json")
    assert _run("--html", "never", "--registry", registry, "--output-dir", "out") == 0
    assert (elsewhere / "out" / "project_report.md").is_file()


def test_max_rows_cuts_the_tables(project):
    assert _run("--html", "never", "--max-rows", "1") == 0
    out = project / "project_report"
    report = json.loads((out / "project_report.json").read_text())
    assert len(report["extractions"]["rows"]) == 1
    assert report["extractions"]["rows_total"] == 4
    assert "prepare results_table" in (out / "project_report.md").read_text()


def test_max_rows_rejects_a_negative_number(project, capsys):
    with pytest.raises(SystemExit):
        create_parser().parse_args(["project_report", "--max-rows", "-1"])


def test_no_environment_skips_run_checks(project, monkeypatch):
    def refuse(*_args, **_kwargs):
        raise AssertionError("run_checks must not be called")

    monkeypatch.setattr(pr, "run_checks", refuse)
    assert _run("--html", "never", "--no-environment") == 0
    report = json.loads((project / "project_report" / "project_report.json").read_text())
    assert report["environment"] == {"included": False}


def test_json_prints_one_document(project, capsys):
    assert _run("--html", "never", "--json") == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["output_dir"] == "project_report"
    assert payload["files"]["markdown"] == "project_report/project_report.md"
    assert payload["files"]["html"] is None
    assert payload["recorded"] is True
    assert payload["sections"] == list(pr.SECTIONS)


def test_json_without_a_registry_prints_an_error_document(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    assert _run("--json") == 1
    assert "status --init" in json.loads(capsys.readouterr().out)["error"]


# --- the interactive extra ----------------------------------------------------------------------


def test_auto_without_jinja2_writes_markdown_and_json_only(project, monkeypatch, caplog):
    _block(monkeypatch, "jinja2")
    caplog.set_level(logging.INFO)
    assert _run() == 0
    out = project / "project_report"
    assert sorted(p.name for p in out.iterdir()) == ["project_report.json", "project_report.md"]
    skipped = [r for r in caplog.records if "prepare[interactive]" in r.getMessage()]
    assert len(skipped) == 1 and skipped[0].levelno == logging.INFO
    assert _exports(project)["project_report"]["summary"]["html"] is False


def test_always_without_jinja2_exits_3_and_writes_nothing(project, monkeypatch, caplog):
    _block(monkeypatch, "jinja2")
    before = (project / "prepare_registry.json").read_text()
    assert _run("--html", "always") == 3
    assert "prepare[interactive]" in caplog.text
    assert not (project / "project_report").exists()
    assert (project / "prepare_registry.json").read_text() == before


def test_always_without_plotly_exits_3(project, monkeypatch):
    _block(monkeypatch, "plotly")
    assert _run("--html", "always") == 3
    assert not (project / "project_report").exists()


def test_auto_with_the_extra_writes_html_with_every_section(project):
    pytest.importorskip("plotly")
    pytest.importorskip("jinja2")
    assert _run() == 0
    out = project / "project_report"
    assert sorted(p.name for p in out.iterdir()) == sorted(FILES)
    html = (out / "project_report.html").read_text()
    for section in pr.SECTIONS:
        assert f'id="{section}"' in html
    assert _exports(project)["project_report"]["summary"]["files"] == list(FILES)


# --- reads no FASTQ -----------------------------------------------------------------------------


def test_no_fastq_file_is_opened(project, monkeypatch):
    opened = []
    real_open, real_io_open, real_gzip_open = builtins.open, io.open, gzip.open

    def watch(real):
        def wrapper(file, *args, **kwargs):
            opened.append(str(file))
            return real(file, *args, **kwargs)

        return wrapper

    # The real environment checks run here too: they must not open a FASTQ file either.
    monkeypatch.setattr(pr, "run_checks", run_checks)
    monkeypatch.setattr(builtins, "open", watch(real_open))
    monkeypatch.setattr(io, "open", watch(real_io_open))
    monkeypatch.setattr(gzip, "open", watch(real_gzip_open))
    assert _run() == 0
    assert opened, "the watch saw no file at all"
    assert not [path for path in opened if ".fastq" in path or path.endswith((".fq", ".fq.gz"))]
    assert any(path.endswith("project_report.json") or "prepare_registry.json" in path for path in opened)


def test_no_fastq_folder_is_listed(project, monkeypatch):
    # Path.iterdir reaches os.listdir and glob reaches os.scandir: neither may list the FASTQ folder.
    listed = []
    real_scandir, real_listdir = os.scandir, os.listdir

    def watch(real):
        def wrapper(path=".", *args, **kwargs):
            if isinstance(path, (str, bytes, os.PathLike)):
                listed.append(os.path.abspath(os.fsdecode(path)))
            return real(path, *args, **kwargs)

        return wrapper

    monkeypatch.setattr(pr, "run_checks", run_checks)
    monkeypatch.setattr(os, "scandir", watch(real_scandir))
    monkeypatch.setattr(os, "listdir", watch(real_listdir))
    fastq = os.path.abspath(project / "fastq")
    assert os.path.isdir(fastq)
    os.listdir(project)  # the watch is in place
    assert listed == [os.path.abspath(project)]
    assert _run() == 0
    assert not [path for path in listed if path == fastq or path.startswith(fastq + os.sep)]


def test_a_byte_that_is_not_utf8_in_the_run_log_still_gives_a_report(project, caplog):
    log = project / ".prepare" / "runs" / "runs.jsonl"
    assert log.is_file()
    with open(log, "ab") as handle:
        handle.write(b"\xff damaged line\n")
    with caplog.at_level(logging.WARNING):
        assert _run("--html", "never") == 0
    report = json.loads((project / "project_report" / "project_report.json").read_text())
    assert report["runs"]["runs"]
    assert "unreadable line" in caplog.text


class TestStageChartIsABarChart:
    """G13: the stages are not nested, so they are drawn as bars from zero, not as a funnel."""

    def test_stage_chart_is_a_bar_chart(self):
        go = pytest.importorskip("plotly.graph_objects")
        from prepare.visualization.project_report import _stage_figure_object

        figure = _stage_figure_object({"screened": {"accessions": 10}, "downloaded": {"accessions": 4}})
        assert isinstance(figure.data[0], go.Bar)
        assert list(figure.data[0].x) == ["screened", "downloaded"] and list(figure.data[0].y) == [10, 4]
        assert figure.layout.title.text == "Datasets per stage"
        assert figure.layout.yaxis.rangemode == "tozero"

    def test_genome_chart_axis_starts_at_zero(self):
        pytest.importorskip("plotly.graph_objects")
        from prepare.visualization.project_report import _genome_figure_object

        figure = _genome_figure_object({"G1": {"extracted": 0, "assembled": 0}})
        assert figure.layout.yaxis.rangemode == "tozero"
