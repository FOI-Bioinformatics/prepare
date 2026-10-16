"""Tests for the project report builder and its Markdown and HTML renderers."""

import pytest

from helpers_project_report import FAILED_MESSAGE, build_project
from prepare.data.registry import load_registry
from prepare.processing import project_report as pr
from prepare.processing.doctor_report import Check
from prepare.processing.project_funnel import funnel
from prepare.processing.project_report_markdown import escape_cell, render_markdown


@pytest.fixture
def registry(tmp_path):
    return load_registry(build_project(tmp_path))


def _fake_checks(*_args, **_kwargs):
    return [
        Check("python", "ok", "Python 3.12"),
        Check("fasterq-dump", "warn", "not found on PATH"),
        Check("registry", "ok", "readable"),
    ]


@pytest.fixture
def no_checks(monkeypatch):
    def refuse(*_args, **_kwargs):
        raise AssertionError("run_checks must not be called")

    monkeypatch.setattr(pr, "run_checks", refuse)


def _report(registry, **kwargs):
    kwargs.setdefault("include_environment", False)
    return pr.build_project_report(registry, **kwargs)


# --- builder ------------------------------------------------------------------------------------


def test_report_has_every_section_in_order(registry, no_checks):
    report = _report(registry)
    assert list(pr.SECTIONS) == [
        "project",
        "funnel",
        "genomes",
        "extractions",
        "downloads",
        "failures",
        "timing",
        "environment",
        "outputs",
        "runs",
    ]
    assert [key for key in report if key in pr.SECTIONS] == list(pr.SECTIONS)
    assert report["report_version"] == pr.REPORT_VERSION
    assert report["generated"].endswith("Z")


def test_project_section_names_registry_and_versions(registry, tmp_path, no_checks):
    project = _report(registry)["project"]
    assert project["registry"] == str(tmp_path / "prepare_registry.json")
    assert project["root"] == str(tmp_path.resolve())
    assert project["registry_version"] == registry.version
    assert project["updated"] == registry.updated
    assert project["datasets"] == 6
    assert project["prepare_version"]


def test_funnel_section_is_the_status_funnel(registry, no_checks):
    assert _report(registry)["funnel"] == funnel(registry)


def test_genomes_section_counts_and_medians(registry, no_checks):
    genomes = _report(registry)["genomes"]
    assert genomes["GCF_A"] == {
        "extracted": 2,
        "assembled": 1,
        "zero_mapped": 1,
        "median_breadth": 0.8,
        "median_depth": 8.0,
    }
    assert genomes["GCF_B"] == {
        "extracted": 1,
        "assembled": 0,
        "zero_mapped": 0,
        "median_breadth": 0.5,
        "median_depth": 2.0,
    }


def test_extraction_rows_sorted_by_genome_then_mapped_reads(registry, no_checks):
    section = _report(registry)["extractions"]
    order = [(row["genome_id"], row["accession"], row["mapped_reads"]) for row in section["rows"]]
    assert order == [("GCF_A", "SRR1", 500), ("GCF_A", "SRR3", 70), ("GCF_A", "SRR2", 0), ("GCF_B", "SRR2", 40)]
    assert section["columns"] == list(pr.EXTRACTION_COLUMNS)
    assert set(section["rows"][0]) == set(pr.EXTRACTION_COLUMNS)
    assert section["rows_total"] == 4
    assert section["truncated"] is False
    assert section["note"] is None


def test_extraction_rows_cut_to_max_rows_with_a_pointer_to_results_table(registry, no_checks):
    section = _report(registry, max_rows=2)["extractions"]
    assert len(section["rows"]) == 2
    assert section["rows_total"] == 4
    assert section["truncated"] is True
    assert "2 of 4" in section["note"] and "results_table" in section["note"]


def test_max_rows_zero_keeps_every_row(registry, no_checks):
    section = _report(registry, max_rows=0)["extractions"]
    assert len(section["rows"]) == section["rows_total"] == 4


def test_downloads_section_counts_verdicts_and_lists(registry, no_checks):
    downloads = _report(registry)["downloads"]
    assert downloads["states"] == {"downloaded": 3, "failed": 1}
    assert downloads["verdicts"] == {"complete": 1, "truncated": 1, "unverified": 1, "none": 0}
    assert downloads["truncated"] == ["SRR2"]
    assert downloads["unverified"] == ["SRR3"]
    assert (downloads["truncated_total"], downloads["unverified_total"]) == (1, 1)
    assert downloads["truncated_note"] is None and downloads["unverified_note"] is None


def test_downloads_lists_and_counts_cover_the_same_downloaded_datasets(registry, no_checks):
    # The Task 25 review's reproduction: SRR3 given a truncated verdict, then SRR2 (truncated)
    # recorded as missing. The verdict table counts SRR3 only, so the list must name SRR3 too.
    from prepare.data.registry import record_download, set_download_verdict

    set_download_verdict(registry, "SRR3", {"method": "spots", "verdict": "truncated", "ratio": 0.4})
    record_download(registry, "SRR2", "missing", registry.path.parent / "fastq")
    downloads = _report(registry, max_rows=1)["downloads"]
    assert downloads["verdicts"]["truncated"] == 1
    assert downloads["truncated"] == ["SRR3"]
    assert downloads["truncated_total"] == 1
    assert downloads["truncated_note"] is None


def test_downloads_list_cut_by_max_rows_carries_its_total_and_a_note(registry, no_checks):
    from prepare.data.registry import set_download_verdict

    set_download_verdict(registry, "SRR1", {"method": "spots", "verdict": "truncated", "ratio": 0.4})
    downloads = _report(registry, max_rows=1)["downloads"]
    assert downloads["truncated"] == ["SRR1"]
    assert downloads["truncated_total"] == 2
    assert downloads["truncated_note"] == "1 of 2 truncated datasets shown"
    markdown = render_markdown(_report(registry, max_rows=1))
    assert "Truncated (first 1 of 2): SRR1" in markdown


def test_failures_section_gives_reason_attempts_and_date(registry, no_checks):
    failures = _report(registry)["failures"]
    assert failures["rows_total"] == 1
    assert failures["by_reason"] == {"network": 1}
    row = failures["rows"][0]
    assert row["accession"] == "SRR4"
    assert row["reason"] == "network"
    assert row["attempts_total"] == 2
    assert "attempts" not in row
    assert row["message"] == FAILED_MESSAGE
    assert row["date"]


def test_timing_section_adds_the_distribution_per_kind(registry, no_checks):
    timing = _report(registry)["timing"]
    assert timing["summary"]["downloads_timed"] == 4
    assert timing["distribution"]["download"] == {
        "count": 4,
        "min": 10.0,
        "p25": 17.5,
        "median": 25.0,
        "p75": 32.5,
        "p90": 37.0,
        "max": 40.0,
    }
    assert timing["distribution"]["extraction"]["median"] == 5.0
    assert timing["distribution"]["assembly"]["p90"] == 60.0


def test_distribution_of_nothing_is_all_none():
    assert pr.distribution([]) == {
        "count": 0,
        "min": None,
        "p25": None,
        "median": None,
        "p75": None,
        "p90": None,
        "max": None,
    }


def test_distribution_of_one_value_is_that_value_everywhere():
    assert pr.distribution([7.0]) == {
        "count": 1,
        "min": 7.0,
        "p25": 7.0,
        "median": 7.0,
        "p75": 7.0,
        "p90": 7.0,
        "max": 7.0,
    }


def test_distribution_of_two_values_interpolates_linearly():
    assert pr.distribution([7.0, 3.0]) == {
        "count": 2,
        "min": 3.0,
        "p25": 4.0,
        "median": 5.0,
        "p75": 6.0,
        "p90": 6.6,
        "max": 7.0,
    }


def test_environment_section_from_run_checks(registry, monkeypatch):
    calls = []

    def fake(*args, **kwargs):
        calls.append(kwargs)
        return _fake_checks()

    monkeypatch.setattr(pr, "run_checks", fake)
    environment = pr.build_project_report(registry)["environment"]
    assert calls and calls[0]["network"] is False
    assert environment["included"] is True
    assert environment["status"] == "warn"
    assert environment["counts"] == {"ok": 2, "warn": 1, "fail": 0}
    assert environment["problems"] == [{"name": "fasterq-dump", "status": "warn", "detail": "not found on PATH"}]


def test_environment_section_skipped_without_run_checks(registry, no_checks):
    assert _report(registry)["environment"] == {"included": False}


def test_outputs_section_lists_exports_and_analysis_outputs(registry, no_checks):
    outputs = _report(registry)["outputs"]
    assert outputs["exports"]["results_table"]["output"] == "results.tsv"
    assert outputs["exports"]["results_table"]["summary"] == {"rows": 4}
    profile = outputs["analyses"]["profile"]
    assert profile["datasets"] == 2
    assert profile["outputs"] == ["profiles/SRR1.json", "profiles/SRR2.json"]
    assert profile["outputs_total"] == 2
    assert profile["latest"]


def test_runs_section_is_the_last_runs_newest_first(registry, no_checks):
    runs = _report(registry, runs_limit=2)["runs"]
    assert runs["total"] == 3
    assert [run["command"] for run in runs["runs"]] == ["results_table", "extract_target_reads"]
    assert set(runs["runs"][0]) == {"run_id", "command", "started", "seconds", "exit_code", "summary"}


def test_runs_section_empty_without_a_run_log(tmp_path, no_checks):
    registry = load_registry(build_project(tmp_path, with_runs=False))
    assert _report(registry)["runs"] == {"total": 0, "runs": []}


def test_report_is_json_serialisable(registry, no_checks):
    import json

    json.dumps(_report(registry))


# --- Markdown -----------------------------------------------------------------------------------


def test_escape_cell_escapes_pipes_and_newlines():
    assert escape_cell("a|b\nc\r\nd") == "a\\|b<br>c<br>d"
    assert escape_cell(None) == ""
    assert escape_cell(0.123456789) == "0.1235"
    assert escape_cell(True) == "yes"


def test_markdown_has_a_heading_per_section(registry, no_checks):
    text = render_markdown(_report(registry))
    for title in pr.SECTION_TITLES.values():
        assert f"\n## {title}\n" in text
    assert text.startswith("# PREPARE project report")
    assert text.isascii()


def test_markdown_table_rows_keep_their_columns_with_pipes_in_cells(registry, no_checks):
    report = _report(registry)
    report["failures"]["rows"][0]["message"] = "bad | line\nsecond"
    text = render_markdown(report)
    line = next(line for line in text.splitlines() if line.startswith("| SRR4 "))
    assert "bad \\| line<br>second" in line
    header = next(line for line in text.splitlines() if line.startswith("| accession | reason"))
    assert line.replace("\\|", "").count("|") == header.count("|")


def test_markdown_truncation_note(registry, no_checks):
    text = render_markdown(_report(registry, max_rows=1))
    assert "1 of 4" in text and "prepare results_table" in text


def test_markdown_failures_name_attempts_over_all_runs(registry, no_checks):
    text = render_markdown(_report(registry))
    assert "| accession | reason | attempts_total | date | message |" in text
    assert "attempts_total counts the attempts over all runs" in text


def test_markdown_environment_not_included(registry, no_checks):
    assert "not included (--no-environment)" in render_markdown(_report(registry))


# --- HTML ---------------------------------------------------------------------------------------


def test_html_has_every_section_id_and_escapes_text(registry, no_checks):
    pytest.importorskip("plotly")
    pytest.importorskip("jinja2")
    from prepare.visualization.project_report import render_html

    report = _report(registry)
    report["failures"]["rows"][0]["message"] = "<script>alert(1)</script>"
    html = render_html(report)
    for section in pr.SECTIONS:
        assert f'id="{section}"' in html
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html
    assert "attempts_total" in html and "over all runs" in html


def test_html_stat_cards_include_the_analysed_stage(registry, no_checks):
    pytest.importorskip("plotly")
    pytest.importorskip("jinja2")
    from prepare.visualization.project_report import render_html

    assert '<p class="k">analysed</p>' in render_html(_report(registry))
