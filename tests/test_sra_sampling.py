"""Tests for prepare.data.sra.sampling: the index-based FASTQ record sampler shared by
quality profiling and the store statistics record, and quality scores kept as histograms."""

import random
import statistics
from unittest.mock import patch

import numpy as np
import pytest

from prepare.data.sra import sample_records
from prepare.data.sra.sampling import distribution_from_histogram, quality_histogram
from tests.perf_fixtures import write_fastq


def test_sample_records_matches_brute_force_selection(tmp_path):
    path = write_fastq(tmp_path / "r.fastq.gz", n=1000)
    got = sample_records([path], 50, seed=3)
    idx = sorted(random.Random(3).sample(range(1000), 50))
    assert [s.decode() for s, _ in got] == [f"A{i}" for i in idx]
    assert all(q == b"I" * len(s) for s, q in got)


def test_sample_records_with_a_given_total_matches_the_counted_one(tmp_path):
    path = write_fastq(tmp_path / "r.fastq", n=1000)
    assert sample_records([path], 50, total_records=1000, seed=7) == sample_records([path], 50, seed=7)


def test_sample_records_returns_all_when_fewer_than_sample_size(tmp_path):
    path = write_fastq(tmp_path / "r.fastq.gz", n=30)
    got = sample_records([path], 100)
    assert [s.decode() for s, _ in got] == [f"A{i}" for i in range(30)]


def test_sample_records_spans_both_mates(tmp_path):
    r1 = write_fastq(tmp_path / "S_1.fastq.gz", n=100, seq_for=lambda i: f"A{i}")
    r2 = write_fastq(tmp_path / "S_2.fastq.gz", n=100, seq_for=lambda i: f"C{i}")
    got = sample_records([r1, r2], 60, seed=11)
    idx = sorted(random.Random(11).sample(range(200), 60))
    expected = [f"A{i}" if i < 100 else f"C{i - 100}" for i in idx]
    assert [s.decode() for s, _ in got] == expected
    assert any(i >= 100 for i in idx)


def test_sample_records_tolerates_overstated_total(tmp_path):
    path = write_fastq(tmp_path / "r.fastq.gz", n=1000)
    got = sample_records([path], 50, total_records=2000, seed=5)
    idx = sorted(i for i in random.Random(5).sample(range(2000), 50) if i < 1000)
    assert [s.decode() for s, _ in got] == [f"A{i}" for i in idx]


def test_sample_records_handles_zero_length_reads_crlf_and_no_final_newline(tmp_path):
    path = tmp_path / "r.fastq"
    path.write_bytes(b"@r0\r\nACGT\r\n+\r\nIIII\r\n@r1\n\n+\n\n@r2\nGG\n+\n#I")
    assert sample_records([path], 10) == [(b"ACGT", b"IIII"), (b"", b""), (b"GG", b"#I")]


@pytest.mark.parametrize("chunk_size", [1, 7, 64, 1000, 1 << 20])
@pytest.mark.parametrize("suffix", [".fastq", ".fastq.gz"])
def test_sample_records_record_across_a_chunk_boundary(tmp_path, monkeypatch, chunk_size, suffix):
    from prepare.data.sra import sampling

    monkeypatch.setattr(sampling, "CHUNK_SIZE", chunk_size)
    r1 = write_fastq(tmp_path / f"S_1{suffix}", n=200)
    r2 = write_fastq(tmp_path / f"S_2{suffix}", n=50, seq_for=lambda i: f"C{i}")
    got = sample_records([r1, r2], 60, seed=2)
    idx = sorted(random.Random(2).sample(range(250), 60))
    assert [s.decode() for s, _ in got] == [f"A{i}" if i < 200 else f"C{i - 200}" for i in idx]


def test_sample_records_rejects_a_file_that_is_not_fastq(tmp_path):
    path = tmp_path / "r.fastq"
    path.write_text(">r0\nACGT\n>r1\nACGT\n")
    with pytest.raises(ValueError, match="FASTQ"):
        sample_records([path], 1, total_records=1)


def test_sample_records_rejects_a_truncated_selected_record(tmp_path):
    path = tmp_path / "r.fastq"
    path.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGT\n")
    with pytest.raises(ValueError, match="Truncated FASTQ record"):
        sample_records([path], 2, total_records=2)


def test_sample_records_rejects_a_malformed_selected_record(tmp_path):
    path = tmp_path / "r.fastq"
    path.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGT\nIIII\nACGT\n")
    with pytest.raises(ValueError, match="Malformed FASTQ record 2 of"):
        sample_records([path], 2, total_records=2)


def test_a_malformed_record_in_mate_2_is_numbered_within_its_file(tmp_path):
    r1 = write_fastq(tmp_path / "S_1.fastq", n=5)
    r2 = tmp_path / "S_2.fastq"
    r2.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n@r3\nACGT\nIIII\nACGT\n")
    with pytest.raises(ValueError, match=r"Malformed FASTQ record 3 of .*S_2\.fastq"):
        sample_records([r1, r2], 8, total_records=8)


def test_sample_records_empty_inputs(tmp_path):
    path = write_fastq(tmp_path / "r.fastq", n=0)
    assert sample_records([path], 10) == []
    assert sample_records([write_fastq(tmp_path / "s.fastq", n=5)], 0) == []


def _list_based(scores):
    """The list-based figures the quality profile reported before histograms."""
    total = len(scores)
    return {
        "mean": statistics.mean(scores),
        "median": statistics.median(scores),
        "q25": float(np.percentile(scores, 25)),
        "q75": float(np.percentile(scores, 75)),
        "excellent_q30+": sum(1 for q in scores if q >= 30) / total,
        "good_q20-29": sum(1 for q in scores if 20 <= q < 30) / total,
        "fair_q10-19": sum(1 for q in scores if 10 <= q < 20) / total,
        "poor_q0-9": sum(1 for q in scores if q < 10) / total,
    }


@pytest.mark.parametrize("n_reads", [1, 2, 7, 500])
def test_quality_histogram_matches_list_based_distribution(n_reads):
    rng = random.Random(n_reads)
    quals = [bytes(33 + rng.randint(0, 41) for _ in range(rng.randint(1, 150))) for _ in range(n_reads)]
    scores = [b - 33 for q in quals for b in q]

    hist = quality_histogram(quals)
    assert hist.shape == (94,) and int(hist.sum()) == len(scores)

    got = distribution_from_histogram(hist)
    for key, value in _list_based(scores).items():
        assert got[key] == pytest.approx(value, abs=1e-9), key


def test_quality_histogram_clamps_out_of_range_characters():
    hist = quality_histogram([b" \x7f~!"])
    assert hist.shape == (94,)
    assert hist[0] == 2 and hist[93] == 2


def test_distribution_from_an_empty_histogram():
    got = distribution_from_histogram(np.zeros(94, dtype=np.int64))
    assert got["mean"] == 0.0 and got["excellent_q30+"] == 0.0


def test_compute_dataset_stats_samples_from_its_own_count(tmp_path):
    """The statistics record passes the count it already has to the sampler, so the sample
    file is not counted a second time."""
    from prepare.store.stats import compute_dataset_stats

    path = write_fastq(tmp_path / "SRR1.fastq.gz", n=300, qual_for=lambda i: "5" * len(f"A{i}"))
    with patch("prepare.data.sra.sampling.count_fastq_reads", side_effect=AssertionError("counted twice")):
        stats = compute_dataset_stats([path], sample_size=100, use_seqkit=False)
    assert stats["reads_total"] == 300 and stats["sampled"] is True
    assert stats["quality_summary"] == {"mean": 20.0, "median": 20.0, "q25": 20.0, "q75": 20.0}


def test_profile_passes_the_cached_record_count_to_the_sampler(tmp_path):
    from prepare.sra.analytics import SRADatasetAnalyzer

    r1 = write_fastq(tmp_path / "SRR1_1.fastq", n=40)
    r2 = write_fastq(tmp_path / "SRR1_2.fastq", n=40)
    analyzer = SRADatasetAnalyzer()
    record = {"reads_per_file": {r1.name: 40, r2.name: 40}, "reads_total": 80}
    with patch("prepare.sra.quality.sample_records", wraps=sample_records) as sampler:
        profile = analyzer.profile_dataset_quality("SRR1", fastq_path=[r1, r2], dataset_stats=record, sample_size=10)
    assert sampler.call_args.kwargs["total_records"] == 80
    assert profile.reads_sampled == 10

    partial = {"reads_per_file": {r1.name: 40}, "reads_total": 40}
    with patch("prepare.sra.quality.sample_records", wraps=sample_records) as sampler:
        analyzer.profile_dataset_quality("SRR1", fastq_path=[r1, r2], dataset_stats=partial, sample_size=10)
    assert sampler.call_args.kwargs["total_records"] is None


def test_compressed_input_is_read_through_a_large_buffer(tmp_path, monkeypatch):
    """A .gz file's compressed bytes are read with the COMPRESSED_READ_BUFFER buffer size."""
    import builtins
    import gzip

    from prepare.data.sra import sampling

    path = tmp_path / "SRR1_1.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@r1\nACGT\n+\nIIII\n@r2\nGGGG\n+\nIIII\n")
    buffers = []
    real_open = builtins.open

    def spy_open(file, mode="r", buffering=-1, *args, **kwargs):
        if str(file) == str(path):
            buffers.append(buffering)
        return real_open(file, mode, buffering, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", spy_open)
    got = sampling.sample_records([path], 2, total_records=2)
    assert got == [(b"ACGT", b"IIII"), (b"GGGG", b"IIII")]
    assert buffers == [sampling.COMPRESSED_READ_BUFFER]


def test_a_much_overstated_total_can_return_no_records(tmp_path):
    """As the docstring says: indices past the end are dropped, possibly all of them."""
    from prepare.data.sra import sampling

    path = tmp_path / "reads.fastq"
    path.write_text("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(7)))
    assert sampling.sample_records([path], 3, total_records=1000, seed=0) == []
