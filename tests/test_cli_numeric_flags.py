"""Numeric flags reject values that mean nothing (zero or negative counts, out-of-range sizes) as usage errors."""

import pytest

from prepare.cli.main import create_parser

EXTRACT = ["extract_target_reads", "--parsed-containment", "p.txt", "--genome-id", "G", "--genome-fasta", "g.fna"]
SAMPLE = ["single_sample", "--summary-column", "s", "--metadata-column", "m"]


@pytest.mark.parametrize(
    "argv",
    [
        ["store_gc", "--older-than", "-1"],
        EXTRACT + ["--threads", "0"],
        EXTRACT + ["--assembly-threads", "-2"],
        EXTRACT + ["--min-contig-len", "0"],
        SAMPLE + ["--top-n", "-5"],
        ["branchwater_search", "--signature", "g.sig", "--registry-max-screened", "-1"],
        ["parse_containment", "--registry-max-screened", "0"],
        ["download_metadata", "--batch-size", "0"],
        ["download_metadata", "--batch-size", "600"],
    ],
)
def test_bad_numeric_values_are_usage_errors(argv, capsys):
    with pytest.raises(SystemExit) as info:
        create_parser().parse_args(argv)
    assert info.value.code == 2
    assert "argument" in capsys.readouterr().err


def test_older_than_zero_is_accepted():
    assert create_parser().parse_args(["store_gc", "--older-than", "0"]).older_than == 0


def test_batch_size_bounds_are_accepted():
    parser = create_parser()
    assert parser.parse_args(["download_metadata", "--batch-size", "1"]).batch_size == 1
    assert parser.parse_args(["download_metadata", "--batch-size", "500"]).batch_size == 500


EXTRACT_THRESHOLD = ["extract_target_reads", "--parsed-containment", "p.txt", "--genome-id", "G", "--genome-fasta", "g"]


@pytest.mark.parametrize("value", ["1.5", "-0.1"])
@pytest.mark.parametrize(
    "argv",
    [
        ["branchwater_search", "--signature", "g.sig", "--threshold"],
        ["select_datasets", "--threshold"],
        EXTRACT_THRESHOLD + ["--threshold"],
        ["download_metadata", "--threshold"],
        ["count_metadata", "--metadata-column", "m", "--threshold"],
        ["results_table", "--min-containment"],
        ["explore_containment", "--parsed-containment", "p.txt", "--min-containment"],
    ],
)
def test_a_containment_threshold_outside_0_to_1_is_a_usage_error(argv, value, capsys):
    with pytest.raises(SystemExit) as info:
        create_parser().parse_args(argv + [value])
    assert info.value.code == 2


def test_max_cache_age_days_must_not_be_negative():
    with pytest.raises(SystemExit):
        create_parser().parse_args(["branchwater_search", "--signature", "g.sig", "--max-cache-age-days", "-1"])


@pytest.mark.parametrize(
    "argv",
    [
        ["explore_containment", "--parsed-containment", "missing.txt"],
        ["store_gc", "--dry-run", "--yes"],
        ["genome_download"],
    ],
)
def test_usage_errors_go_through_fail(argv, tmp_path, monkeypatch, caplog):
    """A usage or input error is one ERROR line with the exception attached (BaseCommand.fail): exit 2 for a
    usage error (ValidationError), 1 for a missing file (DataAccessError)."""
    import logging

    from prepare.cli.main import create_parser

    monkeypatch.chdir(tmp_path)
    args = create_parser().parse_args(argv)
    with caplog.at_level(logging.ERROR):
        rc = args.func(args) if hasattr(args, "func") else None
    errors = [r for r in caplog.records if r.levelno >= logging.ERROR]
    assert rc == (1 if argv[0] == "explore_containment" else 2)
    assert len(errors) == 1 and errors[0].exc_info is not None
