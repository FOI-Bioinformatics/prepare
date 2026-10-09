"""Tests for ``sra_profile``, which replaces ``sra_stats`` and ``sra_profile_quality``.

One statistics path: read and base totals and GC come from the dataset's shared statistics
record (``prepare.store.stats.compute_dataset_stats``), per-read quality from a sample of
every mate file. GC is reported in percent (0-100) in the CSV, the per-accession JSON and the
registry summary.
"""

import argparse
import json
import logging
from unittest.mock import patch

import pandas as pd
import pytest

from prepare.cli.commands.sra_profile import SRAProfileCommand
from prepare.sra.analytics import QualityProfile, SRADatasetAnalyzer
from prepare.sra.profiles import load_quality_profiles
from prepare.store.stats import compute_dataset_stats
from prepare.utils.logging import ConsoleFormatter


def _write_fastq(path, reads):
    """Write ``reads`` (sequences) as a FASTQ file with constant quality."""
    path.write_text("".join(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(reads)))
    return path


def _args(tmp_path, **overrides):
    parser = argparse.ArgumentParser()
    SRAProfileCommand().configure_parser(parser)
    args = parser.parse_args([])
    args.fastq_folder = str(tmp_path / "fastq")
    args.output_report = str(tmp_path / "sra_statistics.csv")
    args.output_dir = str(tmp_path / "profiles")
    args.registry = str(tmp_path / "prepare_registry.json")
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def _paired_dataset(tmp_path, accession="SRR1"):
    acc_dir = tmp_path / "fastq" / accession
    acc_dir.mkdir(parents=True)
    r1 = _write_fastq(acc_dir / f"{accession}_1.fastq", ["GGGGCCCCAT", "GCGCATATAT", "ATATATATAT"])
    r2 = _write_fastq(acc_dir / f"{accession}_2.fastq", ["CCCCCCCCCC", "GGGGGGGGGG", "ATATATATAT"])
    return [r1, r2]


def test_gc_percent_and_totals_come_from_the_shared_statistics_record(tmp_path):
    files = _paired_dataset(tmp_path)
    record = compute_dataset_stats(files, use_seqkit=False)

    assert SRAProfileCommand().execute(_args(tmp_path)) == 0

    row = pd.read_csv(tmp_path / "sra_statistics.csv").set_index("accession").loc["SRR1"]
    assert row["gc_percent"] == record["gc_content"] * 100
    assert row["total_reads"] == record["reads_total"] == 6
    assert "gc_content" not in row.index

    profile = json.loads((tmp_path / "profiles" / "SRR1_quality_profile.json").read_text())
    assert profile["gc_percent"] == record["gc_content"] * 100
    assert profile["total_reads"] == record["reads_total"]
    assert "gc_content" not in profile
    assert "sequence_complexity" not in profile and "complexity_score" in profile

    registry = json.loads((tmp_path / "prepare_registry.json").read_text())
    summary = registry["datasets"]["SRR1"]["analyses"]["profile"]["summary"]
    assert summary["gc_percent"] == record["gc_content"] * 100
    assert summary["total_reads"] == 6


def test_per_read_quality_samples_every_mate_not_mate_1_only(tmp_path):
    _paired_dataset(tmp_path)

    assert SRAProfileCommand().execute(_args(tmp_path)) == 0

    profile = json.loads((tmp_path / "profiles" / "SRR1_quality_profile.json").read_text())
    assert profile["reads_sampled"] == 6


def test_accession_restricts_the_run(tmp_path):
    """The old sra_stats parsed --accessions but profiled every folder anyway."""
    _paired_dataset(tmp_path, "SRR1")
    _paired_dataset(tmp_path, "SRR2")

    assert SRAProfileCommand().execute(_args(tmp_path, accession=["SRR1"])) == 0

    table = pd.read_csv(tmp_path / "sra_statistics.csv")
    assert list(table["accession"]) == ["SRR1"]
    registry = json.loads((tmp_path / "prepare_registry.json").read_text())
    assert set(registry["datasets"]) == {"SRR1"}


def test_signal_between_accessions_stops_before_the_next_one(tmp_path):
    """A stop noticed after SRR1 is profiled ends the run before SRR2 is ever profiled."""
    _paired_dataset(tmp_path, "SRR1")
    _paired_dataset(tmp_path, "SRR2")
    cmd = SRAProfileCommand()
    args = _args(tmp_path, accession=["SRR1", "SRR2"])

    original = cmd._profile_one
    seen = []

    def spy(analyzer, args_, accession):
        result = original(analyzer, args_, accession)
        seen.append(accession)
        if len(seen) == 1:
            args_._termination.stop.set()
        return result

    with patch.object(cmd, "_profile_one", side_effect=spy):
        rc = cmd.run(args)

    assert rc == 130
    assert seen == ["SRR1"]
    registry = json.loads((tmp_path / "prepare_registry.json").read_text())
    assert set(registry["datasets"]) == {"SRR1"}
    assert not (tmp_path / "profiles" / "SRR2_quality_profile.json").exists()


def test_cli_execute_logs_traceback_for_unexpected_error(caplog, monkeypatch, tmp_path):
    """An unexpected error still returns 1, and its traceback is logged."""
    cmd = SRAProfileCommand()

    def buggy(*args, **kwargs):
        raise RuntimeError("bug")

    monkeypatch.setattr(cmd, "_resolve_accessions", buggy)
    (tmp_path / "fastq").mkdir()
    with caplog.at_level("ERROR"):
        assert cmd.execute(_args(tmp_path)) == 1
    assert any(r.exc_info for r in caplog.records)


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def test_parser_defaults_and_repeatable_accession():
    parser = argparse.ArgumentParser()
    SRAProfileCommand().configure_parser(parser)
    args = parser.parse_args([])
    assert (args.fastq_folder, args.output_report, args.output_dir) == (
        "fastq",
        "sra_statistics.csv",
        "sra_quality_profiles",
    )
    assert (args.sample_size, args.sampler) == (10000, "uniform")
    assert args.accession is None and args.accessions_file is None
    assert parser.parse_args(["--accession", "SRR1", "--accession", "SRR2"]).accession == ["SRR1", "SRR2"]


def test_sample_size_rejects_non_positive_values():
    parser = argparse.ArgumentParser()
    SRAProfileCommand().configure_parser(parser)
    assert parser.parse_args(["--sample-size", "500"]).sample_size == 500
    with pytest.raises(SystemExit):
        parser.parse_args(["--sample-size", "0"])


def test_removed_flags_are_gone():
    # The old --accessions still parses, as argparse's abbreviation of --accessions-file.
    parser = argparse.ArgumentParser()
    SRAProfileCommand().configure_parser(parser)
    for flag in ("--include-contamination", "--fastq-dir"):
        with pytest.raises(SystemExit):
            parser.parse_args([flag, "x"])


# ---------------------------------------------------------------------------
# Accession selection and folder layout
# ---------------------------------------------------------------------------


def test_accessions_file_and_accession_flags_combine(tmp_path):
    for accession in ("SRR1", "SRR2", "SRR3"):
        _paired_dataset(tmp_path, accession)
    accessions_file = tmp_path / "acc.txt"
    accessions_file.write_text("# selected\nSRR2\n\n")

    args = _args(tmp_path, accessions_file=str(accessions_file), accession=["SRR3", "SRR2"])
    assert SRAProfileCommand().execute(args) == 0

    assert list(pd.read_csv(tmp_path / "sra_statistics.csv")["accession"]) == ["SRR2", "SRR3"]


def test_every_accession_folder_is_profiled_by_default_and_hidden_ones_are_skipped(tmp_path):
    _paired_dataset(tmp_path, "SRR1")
    hidden = tmp_path / "fastq" / "._SRR1"
    hidden.mkdir()
    _write_fastq(hidden / "SRR1.fastq", ["ACGT"])
    transient = tmp_path / "fastq" / "SRR9_temp"
    transient.mkdir()
    _write_fastq(transient / "SRR9.fastq", ["ACGT"])

    assert SRAProfileCommand().execute(_args(tmp_path)) == 0

    assert list(pd.read_csv(tmp_path / "sra_statistics.csv")["accession"]) == ["SRR1"]


def test_paired_and_single_layouts_are_detected(tmp_path):
    _paired_dataset(tmp_path, "SRR1")
    single = tmp_path / "fastq" / "SRR2"
    single.mkdir()
    _write_fastq(single / "SRR2.fastq", ["ACGTACGT"])

    assert SRAProfileCommand().execute(_args(tmp_path)) == 0

    layouts = pd.read_csv(tmp_path / "sra_statistics.csv").set_index("accession")["layout"]
    assert layouts.to_dict() == {"SRR1": "PAIRED", "SRR2": "SINGLE"}


def test_flat_files_and_prefix_collisions(tmp_path):
    """SRR1 and SRR10 must not collide: each accession profiles its own files only."""
    fastq = tmp_path / "fastq"
    (fastq / "SRR10").mkdir(parents=True)
    _write_fastq(fastq / "SRR10" / "SRR10_1.fastq", ["ACGT"] * 3)
    _write_fastq(fastq / "SRR10" / "SRR10_2.fastq", ["ACGT"] * 3)
    _write_fastq(fastq / "SRR1_1.fastq", ["ACGT"])

    assert SRAProfileCommand().execute(_args(tmp_path, accession=["SRR10", "SRR1"])) == 0

    totals = pd.read_csv(tmp_path / "sra_statistics.csv").set_index("accession")["total_reads"]
    assert totals.to_dict() == {"SRR10": 6, "SRR1": 1}


def test_missing_folder_is_an_error_on_stderr(tmp_path, caplog, capsys):
    args = _args(tmp_path, fastq_folder=str(tmp_path / "nowhere"))
    assert SRAProfileCommand().execute(args) == 1
    assert "does not exist" in caplog.text
    assert "does not exist" not in capsys.readouterr().out


def test_a_folder_without_accession_folders_is_an_error(tmp_path, caplog):
    (tmp_path / "fastq").mkdir()
    (tmp_path / "fastq" / "not_a_dir.txt").write_text("x")
    with caplog.at_level("ERROR"):
        assert SRAProfileCommand().execute(_args(tmp_path)) == 1
    assert "No accession folders" in caplog.text
    assert not any("Traceback" in ConsoleFormatter().format(r) for r in caplog.records)


def test_an_accession_without_fastq_files_fails_and_is_summarised(tmp_path, caplog):
    _paired_dataset(tmp_path, "SRR1")
    empty = tmp_path / "fastq" / "SRR404"
    empty.mkdir()
    (empty / "notes.txt").write_text("x")

    assert SRAProfileCommand().execute(_args(tmp_path)) == 1

    assert "No FASTQ files found for SRR404" in caplog.text
    summary = json.loads((tmp_path / "profiles" / "quality_summary.json").read_text())
    assert summary["total_analyzed"] == 1
    assert summary["failed_accessions"] == ["SRR404"]
    assert list(pd.read_csv(tmp_path / "sra_statistics.csv")["accession"]) == ["SRR1"]


def test_nothing_profiled_writes_no_table_and_a_null_summary(tmp_path):
    (tmp_path / "fastq" / "SRR404").mkdir(parents=True)

    assert SRAProfileCommand().execute(_args(tmp_path)) == 1

    assert not (tmp_path / "sra_statistics.csv").exists()
    summary = json.loads((tmp_path / "profiles" / "quality_summary.json").read_text())
    assert summary["summary_stats"] is None and summary["failed_accessions"] == ["SRR404"]


# ---------------------------------------------------------------------------
# The one statistics path
# ---------------------------------------------------------------------------


def _store_linked(tmp_path, accession="SRR1", reads=("ATCG", "GCTA"), record=None):
    """A dataset in a store folder with a sidecar, linked from the project's fastq folder."""
    from prepare.store.sidecar import Sidecar, write_sidecar

    store_acc_dir = tmp_path / "store" / "sra" / accession
    store_acc_dir.mkdir(parents=True)
    _write_fastq(store_acc_dir / f"{accession}.fastq", list(reads))
    sidecar_path = store_acc_dir / f"{accession}.json"
    write_sidecar(sidecar_path, Sidecar(accession=accession, stats=record))
    (tmp_path / "fastq").mkdir(exist_ok=True)
    (tmp_path / "fastq" / accession).symlink_to(store_acc_dir)
    return store_acc_dir, sidecar_path


def test_the_totals_are_the_cached_record_not_the_sample(tmp_path, capsys):
    """A store link with a matching cached record reuses it: its totals are reported."""
    from prepare.store.sidecar import Sidecar, write_sidecar

    store_acc_dir, sidecar_path = _store_linked(tmp_path, reads=["ACGT"] * 4)
    record = compute_dataset_stats([store_acc_dir / "SRR1.fastq"], use_seqkit=False)
    record.update(reads_total=1724338, bases_total=258650700)
    write_sidecar(sidecar_path, Sidecar(accession="SRR1", stats=record))

    assert SRAProfileCommand().execute(_args(tmp_path, sample_size=2, sampler="head")) == 0

    assert "Total reads (mates counted): 1,724,338 (sampled 2)" in capsys.readouterr().out
    profile = json.loads((tmp_path / "profiles" / "SRR1_quality_profile.json").read_text())
    assert (profile["total_reads"], profile["total_bases"]) == (1724338, 258650700)
    assert (profile["reads_sampled"], profile["sampled"]) == (2, False)
    registry = json.loads((tmp_path / "prepare_registry.json").read_text())
    summary = registry["datasets"]["SRR1"]["analyses"]["profile"]["summary"]
    assert (summary["total_reads"], summary["reads_sampled"]) == (1724338, 2)


def test_a_computed_record_is_written_back_to_the_sidecar(tmp_path):
    from prepare.store.sidecar import read_sidecar

    _, sidecar_path = _store_linked(tmp_path)

    assert SRAProfileCommand().execute(_args(tmp_path)) == 0

    sidecar = read_sidecar(sidecar_path)
    assert sidecar.stats["reads_total"] == 2 and sidecar.stats_computed is not None


def test_the_cache_survives_a_zero_byte_extra_file(tmp_path):
    """The record's signature covers the same file list cached_stats checks, so a zero-byte
    mate left by an interrupted download does not force a recompute on every run."""
    store_acc_dir, _ = _store_linked(tmp_path)
    (store_acc_dir / "SRR1_2.fastq").write_text("")

    assert SRAProfileCommand().execute(_args(tmp_path)) == 0
    with patch("prepare.sra.dataset_stats.compute_dataset_stats") as recompute:
        assert SRAProfileCommand().execute(_args(tmp_path)) == 0
    recompute.assert_not_called()


def test_a_cache_that_cannot_be_written_is_reported_once(tmp_path, caplog):
    _store_linked(tmp_path)
    with caplog.at_level(logging.WARNING):
        with patch("prepare.sra.dataset_stats.store_stats", side_effect=OSError("read-only store")):
            assert SRAProfileCommand().execute(_args(tmp_path)) == 0
    warnings = [r for r in caplog.records if r.levelno == logging.WARNING and "Could not cache" in r.message]
    assert len(warnings) == 1
    assert "SRR1" in caplog.text and "read-only store" in caplog.text


def test_a_plain_project_folder_gets_no_sidecar(tmp_path):
    _paired_dataset(tmp_path)
    assert SRAProfileCommand().execute(_args(tmp_path)) == 0
    assert not (tmp_path / "fastq" / "SRR1" / "SRR1.json").exists()


def test_totals_stay_exact_when_the_per_read_metrics_are_sampled(tmp_path):
    acc_dir = tmp_path / "fastq" / "SRR1"
    acc_dir.mkdir(parents=True)
    for name in ("SRR1_1.fastq", "SRR1_2.fastq"):
        _write_fastq(acc_dir / name, ["ACGT"] * 5)

    assert SRAProfileCommand().execute(_args(tmp_path, sample_size=2)) == 0

    row = pd.read_csv(tmp_path / "sra_statistics.csv").loc[0]
    assert (row["total_reads"], row["total_bases"], row["reads_sampled"]) == (10, 40, 2)
    assert not bool(row["sampled"])


_ORIGINAL_PROFILE = SRADatasetAnalyzer.profile_dataset_quality


def _real_profile(self, *args, **kwargs):
    return _ORIGINAL_PROFILE(self, *args, **kwargs)


def test_sample_size_and_sampler_reach_the_analyzer(tmp_path):
    _paired_dataset(tmp_path)
    with patch.object(
        SRADatasetAnalyzer, "profile_dataset_quality", autospec=True, side_effect=_real_profile
    ) as profile_call:
        assert SRAProfileCommand().execute(_args(tmp_path, sample_size=500, sampler="head")) == 0
    kwargs = profile_call.call_args.kwargs
    assert (kwargs["sample_size"], kwargs["sampler"]) == (500, "head")
    assert [p.name for p in kwargs["fastq_path"]] == ["SRR1_1.fastq", "SRR1_2.fastq"]


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def _profile(accession="SRR1", n_content=0.0, duplication_rate=None, adapter=0.0):
    return QualityProfile(
        accession=accession,
        total_reads=1000,
        reads_sampled=1000,
        total_bases=150000,
        avg_read_length=150.0,
        read_length_distribution={},
        gc_percent=45.0,
        gc_histogram={},
        quality_distribution={"excellent_q30+": 0.9},
        n_content=n_content,
        contamination_indicators={"adapter_contamination": adapter},
        complexity_score=0.85,
        duplication_rate=duplication_rate,
        technology_confidence=0.8,
        quality_grade="good",
        recommendations=[],
    )


def test_printed_profile_flags_problems_in_plain_ascii(capsys):
    SRAProfileCommand()._print_profile(_profile(n_content=0.02, duplication_rate=0.25, adapter=0.08))
    out = capsys.readouterr().out
    assert "Total reads (mates counted): 1,000 (sampled 1,000)" in out
    assert "GC content: 45.0%" in out
    assert "WARNING: high N content: 2.0%" in out
    assert "WARNING: high duplicate rate: 25.0%" in out
    assert "WARNING: adapter contamination: 8.0%" in out
    assert out.isascii()


def test_printed_profile_handles_an_unknown_duplicate_rate(capsys):
    SRAProfileCommand()._print_profile(_profile(duplication_rate=None))
    assert "duplicate rate" not in capsys.readouterr().out


def test_summary_lines_go_to_stdout_and_label_mates(tmp_path, capsys):
    _paired_dataset(tmp_path)
    assert SRAProfileCommand().execute(_args(tmp_path, summary_only=True)) == 0
    out = capsys.readouterr().out
    assert "Total reads (mates counted): 6" in out
    assert "Average GC content:" in out
    assert "Quality profile: SRR1" not in out  # --summary-only
    assert out.isascii()


def test_detailed_reports_prints_each_json_path(tmp_path, capsys):
    _paired_dataset(tmp_path)
    assert SRAProfileCommand().execute(_args(tmp_path, detailed_reports=True)) == 0
    assert "Profile saved:" in capsys.readouterr().out


def test_profile_json_round_trips_through_load_quality_profiles(tmp_path):
    _paired_dataset(tmp_path)
    assert SRAProfileCommand().execute(_args(tmp_path)) == 0

    reloaded = load_quality_profiles(tmp_path / "profiles")
    assert set(reloaded) == {"SRR1"}
    assert reloaded["SRR1"].total_reads == 6 and reloaded["SRR1"].sampled is False
    raw = json.loads((tmp_path / "profiles" / "SRR1_quality_profile.json").read_text())
    for key in ("read_length_distribution", "gc_histogram", "technology_confidence", "mean_quality"):
        assert key in raw


# ---------------------------------------------------------------------------
# Registry and store
# ---------------------------------------------------------------------------


def test_registry_records_the_profile_relative_to_the_project(tmp_path):
    _paired_dataset(tmp_path)
    assert SRAProfileCommand().execute(_args(tmp_path)) == 0

    analysis = json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]["SRR1"]["analyses"]
    assert set(analysis) == {"profile"}
    assert analysis["profile"]["output"] == "profiles/SRR1_quality_profile.json"
    assert set(analysis["profile"]["summary"]) == {
        "total_reads",
        "total_bases",
        "reads_sampled",
        "sampled",
        "gc_percent",
        "avg_read_length",
        "quality_grade",
    }


def _project_with_store(tmp_path):
    from prepare.data.registry import load_registry, save_registry
    from prepare.store.layout import init_store

    store_root = tmp_path / "store_root"
    init_store(store_root)
    registry = load_registry(tmp_path / "prepare_registry.json")
    registry.project = {"id": "proj1", "name": "demo", "path": str(tmp_path), "created": "now"}
    save_registry(registry)
    return store_root


def test_usage_is_recorded_in_the_store_catalogue(tmp_path):
    from prepare.store.catalog import Catalog
    from prepare.store.layout import store_paths

    _paired_dataset(tmp_path)
    store_root = _project_with_store(tmp_path)

    assert SRAProfileCommand().execute(_args(tmp_path, data_root=str(store_root))) == 0

    with Catalog(store_paths(store_root)) as catalog:
        catalog.migrate()
        row = catalog.conn.execute(
            "SELECT stage FROM usage WHERE accession = ? AND project_id = ?", ("SRR1", "proj1")
        ).fetchone()
    assert row["stage"] == "analysed"


def test_a_catalogue_failure_leaves_the_outcome_unchanged(tmp_path):
    _paired_dataset(tmp_path)
    store_root = _project_with_store(tmp_path)
    with patch("prepare.store.usage.catalog_write", side_effect=RuntimeError("locked")):
        assert SRAProfileCommand().execute(_args(tmp_path, data_root=str(store_root))) == 0
    registry = json.loads((tmp_path / "prepare_registry.json").read_text())
    assert registry["datasets"]["SRR1"]["analyses"]["profile"]["summary"]["total_reads"] == 6


def test_an_unreachable_store_costs_the_usage_record_only(tmp_path, caplog):
    _paired_dataset(tmp_path)
    with caplog.at_level(logging.WARNING):
        assert SRAProfileCommand().execute(_args(tmp_path, data_root=str(tmp_path / "unmounted"))) == 0
    assert (tmp_path / "sra_statistics.csv").is_file()
    assert any("store unavailable" in r.message for r in caplog.records)


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


def test_an_expected_error_is_one_line_without_traceback(caplog, tmp_path):
    (tmp_path / "fastq").mkdir()
    args = _args(tmp_path, accessions_file=str(tmp_path / "missing.txt"))
    with caplog.at_level("ERROR"):
        assert SRAProfileCommand().execute(args) == 1
    assert any("Profiling failed" in r.message and "missing.txt" in r.message for r in caplog.records)
    assert not any("Traceback" in ConsoleFormatter().format(r) for r in caplog.records)


def test_a_bug_while_profiling_one_accession_is_not_counted_as_a_failed_accession(tmp_path, caplog):
    _paired_dataset(tmp_path)
    with patch("prepare.cli.commands.sra_profile.profile_accession", side_effect=TypeError("bug")):
        with caplog.at_level("ERROR"):
            assert SRAProfileCommand().execute(_args(tmp_path)) == 1
    assert any(r.exc_info for r in caplog.records)
    assert "Failed to profile" not in caplog.text


def _write_corrupt_gzip_fastq(path):
    """A gzip FASTQ whose deflate stream has one flipped byte, so reading it raises zlib.error."""
    import gzip as _gzip

    records = "".join(f"@r{i}\n{'ACGT' * 25}\n+\n{'I' * 100}\n" for i in range(200))
    data = bytearray(_gzip.compress(records.encode()))
    data[20] ^= 0xFF
    path.write_bytes(bytes(data))
    return path


def test_a_corrupt_gzip_fails_its_accession_and_the_run_continues(tmp_path):
    fastq = tmp_path / "fastq"
    (fastq / "SRR1").mkdir(parents=True)
    _write_corrupt_gzip_fastq(fastq / "SRR1" / "SRR1.fastq.gz")
    _paired_dataset(tmp_path, "SRR2")

    with patch("prepare.store.stats.shutil.which", return_value=None):
        assert SRAProfileCommand().execute(_args(tmp_path)) == 1

    summary = json.loads((tmp_path / "profiles" / "quality_summary.json").read_text())
    assert summary["failed_accessions"] == ["SRR1"] and summary["total_analyzed"] == 1


def test_load_dataset_stats_is_none_for_unreadable_files_and_a_bug_propagates(tmp_path):
    from prepare.sra.dataset_stats import load_dataset_stats

    acc_dir = tmp_path / "SRR1"
    acc_dir.mkdir()
    corrupt = _write_corrupt_gzip_fastq(acc_dir / "SRR1.fastq.gz")
    with patch("prepare.store.stats.shutil.which", return_value=None):
        assert load_dataset_stats([corrupt]) is None
    good = _write_fastq(acc_dir / "SRR1_1.fastq", ["ACGT"])
    with patch("prepare.sra.dataset_stats.compute_dataset_stats", side_effect=TypeError("bug")):
        with pytest.raises(TypeError):
            load_dataset_stats([good])


def test_load_dataset_stats_cache_write_failure_is_logged_and_a_bug_propagates(tmp_path, caplog):
    from prepare.sra import dataset_stats

    acc_dir = tmp_path / "SRR1"
    acc_dir.mkdir()
    fastq = _write_fastq(acc_dir / "SRR1.fastq", ["ACGT"])
    with (
        patch.object(dataset_stats, "_resolved_sidecar_path", return_value=tmp_path / "SRR1.json"),
        patch.object(dataset_stats, "cached_stats", return_value=None),
    ):
        with patch.object(dataset_stats, "store_stats", side_effect=PermissionError("read-only")):
            with caplog.at_level("WARNING"):
                assert dataset_stats.load_dataset_stats([fastq]) is not None
        assert "Could not cache statistics" in caplog.text
        with patch.object(dataset_stats, "store_stats", side_effect=TypeError("bug")):
            with pytest.raises(TypeError):
                dataset_stats.load_dataset_stats([fastq])
