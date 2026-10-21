"""Tests for ``sra_report``, which replaces ``sra_dashboard`` and ``sra_compare``.

Each accession is profiled once (or read back from a saved profile), and that one set of
profiles feeds the quality section and, with ``--groups-file``, the comparison of groups.
"""

import argparse
import json
import logging
import re
from unittest.mock import patch

import numpy as np
import pytest

from prepare.cli.commands.sra_report import SRAReportCommand, load_groups
from prepare.core.exceptions import DataAccessError, ValidationError
from prepare.sra.analytics import ComparativeAnalysis, QualityProfile, SRADatasetAnalyzer
from prepare.sra.profiles import write_profile_json
from prepare.utils.logging import ConsoleFormatter


def _write_fastq(path, reads):
    path.write_text("".join(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(reads)))
    return path


def _args(tmp_path, **overrides):
    parser = argparse.ArgumentParser()
    SRAReportCommand().configure_parser(parser)
    args = parser.parse_args([])
    args.fastq_folder = str(tmp_path / "fastq")
    args.output_dir = str(tmp_path / "reports")
    args.registry = str(tmp_path / "prepare_registry.json")
    args.no_open = True
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def _datasets(tmp_path, accessions=("SRR1", "SRR2", "SRR3", "SRR4")):
    """One single-end dataset per accession; read count, length and GC differ between them so the
    group tests see real variation."""
    bases = ["GGCCATATCG", "GGGGCCATTA", "ATATATGCAA", "GCGCGCATGG", "ACGTACGATT", "TTTTGGCCAG"]
    for i, (accession, seq) in enumerate(zip(accessions, bases)):
        acc_dir = tmp_path / "fastq" / accession
        acc_dir.mkdir(parents=True)
        reads = [(seq * (i + 1))[: 8 + 3 * i + j] for j in range(3 + i)]
        _write_fastq(acc_dir / f"{accession}.fastq", reads)


def _groups_file(tmp_path, groups=None):
    path = tmp_path / "groups.json"
    path.write_text(json.dumps(groups or {"A": ["SRR1", "SRR2"], "B": ["SRR3", "SRR4"]}))
    return str(path)


def make_profile(accession, gc_percent=45.0, grade="good"):
    return QualityProfile(
        accession=accession,
        total_reads=1000,
        reads_sampled=1000,
        total_bases=150000,
        avg_read_length=150.0,
        read_length_distribution={},
        gc_percent=gc_percent,
        gc_histogram={},
        quality_distribution={"excellent_q30+": 0.9},
        n_content=0.0,
        contamination_indicators={"adapter_contamination": 0.0},
        complexity_score=0.85,
        duplication_rate=None,
        technology_confidence=0.8,
        quality_grade=grade,
        recommendations=[],
    )


def _saved_profiles(tmp_path, accessions=("SRR1", "SRR2")):
    folder = tmp_path / "profiles"
    folder.mkdir()
    for i, accession in enumerate(accessions):
        write_profile_json(make_profile(accession, gc_percent=40.0 + i), folder)
    return str(folder)


# ---------------------------------------------------------------------------
# Profiled once
# ---------------------------------------------------------------------------


def test_groups_file_profiles_each_accession_once(tmp_path, monkeypatch):
    _datasets(tmp_path)
    groups = _groups_file(tmp_path)

    calls = []
    original = SRADatasetAnalyzer.profile_dataset_quality

    def spy(self, accession, *args, **kwargs):
        calls.append(accession)
        return original(self, accession, *args, **kwargs)

    monkeypatch.setattr(SRADatasetAnalyzer, "profile_dataset_quality", spy)

    assert SRAReportCommand().execute(_args(tmp_path, groups_file=groups)) == 0

    assert sorted(calls) == ["SRR1", "SRR2", "SRR3", "SRR4"]
    assert (tmp_path / "reports" / "sra_report.html").is_file()
    registry = json.loads((tmp_path / "prepare_registry.json").read_text())
    assert registry["datasets"]["SRR1"]["analyses"]["report"]["output"] == "reports/sra_report.html"


def test_one_report_holds_the_quality_and_the_comparison_sections(tmp_path):
    _datasets(tmp_path)
    assert SRAReportCommand().execute(_args(tmp_path, groups_file=_groups_file(tmp_path), title="Wolbachia")) == 0

    html = (tmp_path / "reports" / "sra_report.html").read_text()
    assert "Wolbachia" in html
    assert "Average GC content" in html  # quality section
    assert "Statistical tests" in html  # comparison section
    assert sorted(p.name for p in (tmp_path / "reports").iterdir()) == ["sra_report.html", "sra_report.json"]


def test_a_group_with_no_profiled_accession_is_dropped_with_a_warning(tmp_path, caplog):
    """A group whose accessions have no FASTQ files does not crash the report; it is left out and said so."""
    _datasets(tmp_path, ("SRR1", "SRR2", "SRR3"))
    groups = _groups_file(tmp_path, {"A": ["SRR1", "SRR2"], "B": ["SRR3"], "C": ["SRR9"]})

    with caplog.at_level(logging.WARNING):
        assert SRAReportCommand().execute(_args(tmp_path, groups_file=groups)) == 0

    assert (tmp_path / "reports" / "sra_report.html").is_file()
    assert "C" in caplog.text and "SRR9" in caplog.text
    html = (tmp_path / "reports" / "sra_report.html").read_text()
    assert "Statistical tests" in html and "p-value nan" not in html and "p-value <b>nan" not in html


def test_a_single_group_left_gives_no_comparison(tmp_path, caplog):
    _datasets(tmp_path, ("SRR1", "SRR2"))
    groups = _groups_file(tmp_path, {"A": ["SRR1", "SRR2"], "B": ["SRR9"]})

    with caplog.at_level(logging.WARNING):
        assert SRAReportCommand().execute(_args(tmp_path, groups_file=groups)) == 0

    html = (tmp_path / "reports" / "sra_report.html").read_text()
    assert "Statistical tests" not in html
    assert "Fewer than two groups" in caplog.text


def test_total_reads_counts_mates_and_leaves_sampled_profiles_out(tmp_path):
    """The card sums exact totals only, says mates are counted, and names the sampled datasets."""
    from prepare.sra.reporting import SRAReportGenerator

    generator = SRAReportGenerator(tmp_path / "reports")
    full = make_profile("SRR1")
    sampled = make_profile("SRR2")
    sampled.sampled = True
    sampled.total_reads = 10
    summary = generator._calculate_quality_summary({"SRR1": full, "SRR2": sampled})

    assert summary["total_reads"] == 1000
    assert summary["sampled_datasets"] == 1
    html = generator.generate_report({"SRR1": full, "SRR2": sampled})
    text = html.read_text()
    assert "mates counted" in text
    assert re.search(r"1\s+dataset\(s\) sampled only, not counted", text)


def test_without_groups_there_is_no_comparison(tmp_path):
    _datasets(tmp_path, ("SRR1", "SRR2"))
    accessions = tmp_path / "acc.txt"
    accessions.write_text("SRR1\nSRR2\n")

    assert SRAReportCommand().execute(_args(tmp_path, accessions_file=str(accessions))) == 0

    html = (tmp_path / "reports" / "sra_report.html").read_text()
    assert "Statistical tests" not in html
    assert json.loads((tmp_path / "reports" / "sra_report.json").read_text())["comparison"] is None


def test_saved_profiles_are_reused_and_not_profiled_again(tmp_path):
    profiles = _saved_profiles(tmp_path, ("SRR1", "SRR2", "SRR3", "SRR4"))
    with patch.object(SRADatasetAnalyzer, "profile_dataset_quality") as profile_call:
        args = _args(tmp_path, quality_profiles=profiles, groups_file=_groups_file(tmp_path))
        assert SRAReportCommand().execute(args) == 0
    profile_call.assert_not_called()
    saved = json.loads((tmp_path / "reports" / "sra_report.json").read_text())
    assert saved["gc_percent"] == {"SRR1": 40.0, "SRR2": 41.0, "SRR3": 42.0, "SRR4": 43.0}
    assert set(saved["comparison"]["groups"]) == {"A", "B"}


def test_saved_profiles_alone_name_the_accessions(tmp_path, capsys):
    profiles = _saved_profiles(tmp_path)
    assert SRAReportCommand().execute(_args(tmp_path, quality_profiles=profiles)) == 0
    assert json.loads((tmp_path / "reports" / "sra_report.json").read_text())["accessions"] == ["SRR1", "SRR2"]
    assert "Reusing 2 saved quality profile(s)" in capsys.readouterr().out


def test_an_accession_without_fastq_is_left_out_with_a_warning(tmp_path, caplog):
    _datasets(tmp_path, ("SRR1", "SRR2"))
    accessions = tmp_path / "acc.txt"
    accessions.write_text("SRR1\nSRR2\nSRR404\n")
    with caplog.at_level(logging.WARNING):
        assert SRAReportCommand().execute(_args(tmp_path, accessions_file=str(accessions))) == 0
    assert "SRR404" in caplog.text
    saved = json.loads((tmp_path / "reports" / "sra_report.json").read_text())
    assert saved["failed_accessions"] == ["SRR404"]


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def test_the_comparison_and_tests_are_printed(tmp_path, capsys):
    _datasets(tmp_path)
    assert SRAReportCommand().execute(_args(tmp_path, groups_file=_groups_file(tmp_path))) == 0
    out = capsys.readouterr().out
    assert "Comparison of groups:" in out
    assert "Mean total reads (mates counted):" in out
    assert "gc_percent: t-test p=" in out
    assert out.isascii()


def test_print_comparison_marks_significant_tests(capsys):
    comparison = ComparativeAnalysis(
        dataset_groups={"A": ["SRR1"], "B": ["SRR2"]},
        summary_statistics={"A": {"gc_percent": {"mean": 45.0}, "total_reads": {"mean": 2000000.0}}},
        statistical_tests={"gc_percent": {"test": "t-test", "p_value": 0.03, "significant": True}},
        outlier_datasets=[],
        clustering_results=None,
        batch_effects={},
        recommendations=[],
        visualization_data={},
    )
    SRAReportCommand()._print_comparison(comparison)
    out = capsys.readouterr().out
    assert "Mean GC content: 45.0%" in out
    assert "Mean total reads (mates counted): 2,000,000" in out
    assert "gc_percent: t-test p=0.0300 *" in out


def test_print_comparison_shows_the_adjusted_p_and_follows_it(capsys):
    comparison = ComparativeAnalysis(
        dataset_groups={"A": ["SRR1"], "B": ["SRR2"]},
        summary_statistics={},
        statistical_tests={"gc_percent": {"test": "t-test", "p_value": 0.03, "p_adjusted": 0.12, "significant": False}},
        outlier_datasets=[],
        clustering_results=None,
        batch_effects={},
        recommendations=[],
        visualization_data={},
    )
    SRAReportCommand()._print_comparison(comparison)
    out = capsys.readouterr().out
    assert "gc_percent: t-test p=0.0300 adjusted=0.1200" in out and "*" not in out


def test_report_shows_raw_and_adjusted_p(tmp_path):
    """The page carries the raw and the Holm-adjusted p-value of each test that ran."""
    from prepare.sra.reporting import _render_page

    comparison = ComparativeAnalysis(
        dataset_groups={"A": ["SRR1", "SRR2"], "B": ["SRR3", "SRR4"]},
        summary_statistics={},
        statistical_tests={
            "gc_percent": {"test": "t-test", "p_value": 0.01, "p_adjusted": 0.04, "significant": True},
            "total_reads": {"test": "skipped (identical values)", "p_value": 1.0, "significant": False, "note": "x"},
        },
        outlier_datasets=[],
        clustering_results=None,
        batch_effects={},
        recommendations=[],
        visualization_data={},
    )
    html = _render_page(
        "t", "t", "now", quality=None, comparative={"comparison": comparison, "plots": {}, "group_counts": {}}
    )
    assert "p 0.0100, adjusted 0.0400 (Holm)" in html
    assert "p 1.0000" in html and "adjusted" not in html.split("p 1.0000")[1].split("</p>")[0]


def test_report_json_is_numpy_safe(tmp_path):
    comparison = ComparativeAnalysis(
        dataset_groups={"a": ["A1"], "b": ["B1"]},
        summary_statistics={},
        statistical_tests={
            "gc_percent": {
                "test": "t-test",
                "statistic": np.float64(2.0),
                "p_value": np.float64(0.03),
                "significant": np.bool_(True),
            }
        },
        outlier_datasets=[],
        clustering_results=None,
        batch_effects={},
        recommendations=[],
        visualization_data={},
    )
    from prepare.sra.analytics import AnomalyReport

    anomalies = AnomalyReport([], {}, {}, {}, {})
    path = SRAReportCommand()._write_json(tmp_path, {"A1": make_profile("A1")}, [], anomalies, comparison)

    saved = json.loads(path.read_text())
    assert saved["comparison"]["statistical_tests"]["gc_percent"]["significant"] is True
    assert saved["comparison"]["significant_differences"] == ["gc_percent"]


def test_no_report_writes_only_the_json_and_needs_no_plotting_packages(tmp_path, monkeypatch):
    import sys

    _datasets(tmp_path, ("SRR1",))
    monkeypatch.setitem(sys.modules, "plotly", None)
    args = _args(tmp_path, no_report=True, no_open=False, groups_file=_groups_file(tmp_path, {"A": ["SRR1"]}))
    with patch("prepare.cli.commands.sra_report.open_in_browser") as opener:
        assert SRAReportCommand().execute(args) == 0
    opener.assert_not_called()
    assert sorted(p.name for p in (tmp_path / "reports").iterdir()) == ["sra_report.json"]
    registry = json.loads((tmp_path / "prepare_registry.json").read_text())
    assert registry["datasets"]["SRR1"]["analyses"]["report"]["output"] == "reports/sra_report.json"


@pytest.mark.parametrize("package", ["plotly", "jinja2"])
def test_a_missing_interactive_extra_fails_before_any_work_or_output(tmp_path, monkeypatch, caplog, package):
    import sys

    _datasets(tmp_path, ("SRR1", "SRR2"))
    monkeypatch.setitem(sys.modules, package, None)
    with patch.object(SRADatasetAnalyzer, "profile_dataset_quality") as profile_call:
        with caplog.at_level(logging.ERROR):
            assert SRAReportCommand().execute(_args(tmp_path, groups_file=_groups_file(tmp_path))) == 3
    profile_call.assert_not_called()
    assert "prepare[interactive]" in caplog.text
    assert not (tmp_path / "reports").exists()
    assert not (tmp_path / "prepare_registry.json").exists()


@pytest.mark.parametrize("no_report", [False, True])
def test_a_missing_analysis_extra_fails_before_any_work_with_groups_file(tmp_path, monkeypatch, caplog, no_report):
    """--groups-file needs scipy (for the statistical tests); checked up front, not after profiling."""
    import sys

    _datasets(tmp_path, ("SRR1", "SRR2"))
    monkeypatch.setitem(sys.modules, "scipy", None)
    with patch.object(SRADatasetAnalyzer, "profile_dataset_quality") as profile_call:
        with caplog.at_level(logging.ERROR):
            args = _args(tmp_path, groups_file=_groups_file(tmp_path), no_report=no_report)
            assert SRAReportCommand().execute(args) == 3
    profile_call.assert_not_called()
    assert "prepare[analysis]" in caplog.text
    assert not (tmp_path / "reports").exists()
    assert not (tmp_path / "prepare_registry.json").exists()


def test_without_groups_file_a_missing_scipy_does_not_matter(tmp_path, monkeypatch):
    """scipy is only needed for the group comparison; a plain report must not require it."""
    import sys

    _datasets(tmp_path, ("SRR1",))
    monkeypatch.setitem(sys.modules, "scipy", None)
    profiles = _saved_profiles(tmp_path, ("SRR1",))
    assert SRAReportCommand().execute(_args(tmp_path, quality_profiles=profiles)) == 0


def test_the_report_is_opened_unless_no_open(tmp_path):
    _datasets(tmp_path, ("SRR1",))
    accessions = tmp_path / "acc.txt"
    accessions.write_text("SRR1\n")
    with patch("prepare.cli.commands.sra_report.open_in_browser", return_value=True) as opener:
        assert SRAReportCommand().execute(_args(tmp_path, accessions_file=str(accessions), no_open=False)) == 0
    opener.assert_called_once_with(tmp_path / "reports" / "sra_report.html")


def test_registry_summary_names_group_grade_and_gc_percent(tmp_path):
    _datasets(tmp_path)
    assert SRAReportCommand().execute(_args(tmp_path, groups_file=_groups_file(tmp_path))) == 0
    summary = json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]["SRR3"]["analyses"]["report"][
        "summary"
    ]
    assert summary["group"] == "B"
    assert set(summary) == {"quality_grade", "gc_percent", "group", "anomalous", "total_reads", "total_bases"}


# ---------------------------------------------------------------------------
# Accession sources and errors
# ---------------------------------------------------------------------------


def test_parser_defaults():
    parser = argparse.ArgumentParser()
    SRAReportCommand().configure_parser(parser)
    args = parser.parse_args([])
    assert (args.fastq_folder, args.output_dir, args.title) == ("fastq", "sra_reports", "SRA Report")
    assert args.groups_file is None and args.quality_profiles is None and not args.no_report
    for removed in ("--dashboard-type", "--statistical-tests", "--generate-report", "--fastq-dir"):
        with pytest.raises(SystemExit):
            parser.parse_args([removed, "x"])


def test_no_accession_source_is_an_error(tmp_path, caplog):
    with caplog.at_level(logging.ERROR):
        assert SRAReportCommand().execute(_args(tmp_path)) == 2
    assert any("Give --accessions-file, --groups-file or a --quality-profiles" in r.message for r in caplog.records)
    assert not any("Traceback" in ConsoleFormatter().format(r) for r in caplog.records)


def test_no_fastq_for_any_accession_is_an_error(tmp_path, caplog):
    with caplog.at_level(logging.ERROR):
        assert SRAReportCommand().execute(_args(tmp_path, groups_file=_groups_file(tmp_path))) == 1
    assert "No dataset could be profiled" in caplog.text


@pytest.mark.parametrize("kind", ["missing", "file"])
def test_a_quality_profiles_path_that_is_not_a_folder_is_logged(tmp_path, caplog, kind):
    path = tmp_path / "profiles"
    if kind == "file":
        path.write_text("not a folder")
    with caplog.at_level(logging.WARNING):
        assert SRAReportCommand().execute(_args(tmp_path, quality_profiles=str(path))) == 2
    expected = "is not a directory" if kind == "file" else "not found"
    assert str(path) in caplog.text and expected in caplog.text


def test_load_groups_success(tmp_path):
    assert load_groups(_groups_file(tmp_path)) == {"A": ["SRR1", "SRR2"], "B": ["SRR3", "SRR4"]}


@pytest.mark.parametrize(
    "content,message",
    [("{invalid json", "Invalid JSON"), ('["SRR1"]', "must map group names")],
)
def test_load_groups_rejects_bad_files(tmp_path, content, message):
    path = tmp_path / "groups.json"
    path.write_text(content)
    with pytest.raises(ValidationError, match=message):
        load_groups(str(path))


def test_a_missing_groups_file_is_a_data_access_error(tmp_path):
    with pytest.raises(DataAccessError, match="Cannot read groups file") as caught:
        load_groups(str(tmp_path / "groups.json"))
    assert not isinstance(caught.value, ValidationError)


def test_a_missing_groups_file_exits_1(tmp_path, caplog):
    with caplog.at_level(logging.ERROR):
        assert SRAReportCommand().execute(_args(tmp_path, groups_file=str(tmp_path / "groups.json"))) == 1
    assert "Cannot read groups file" in caplog.text


def test_a_bad_groups_file_is_an_error_on_stderr_not_stdout(tmp_path, caplog, capsys):
    path = tmp_path / "groups.json"
    path.write_text("{invalid json")
    with caplog.at_level(logging.ERROR):
        assert SRAReportCommand().execute(_args(tmp_path, groups_file=str(path))) == 2
    assert "Invalid JSON" in caplog.text
    assert "Invalid JSON" not in capsys.readouterr().out


def test_an_expected_error_is_one_line_without_traceback(caplog, tmp_path):
    args = _args(tmp_path, accessions_file=str(tmp_path / "missing.txt"))
    with caplog.at_level("ERROR"):
        assert SRAReportCommand().execute(args) == 1
    assert any("Report failed" in r.message and "missing.txt" in r.message for r in caplog.records)
    assert not any("Traceback" in ConsoleFormatter().format(r) for r in caplog.records)


def test_cli_execute_logs_traceback_for_unexpected_error(caplog, monkeypatch, tmp_path):
    """An unexpected error still returns 1, and its traceback is logged."""
    cmd = SRAReportCommand()

    def buggy(*args, **kwargs):
        raise RuntimeError("bug")

    monkeypatch.setattr(cmd, "_saved_profiles", buggy)
    with caplog.at_level("ERROR"):
        assert cmd.execute(_args(tmp_path)) == 1
    assert any(r.exc_info for r in caplog.records)
