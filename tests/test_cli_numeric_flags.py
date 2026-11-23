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
