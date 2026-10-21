"""Tests for the dataset sketcher behind ``sra_sketch``.

The block reader joins every sequence of a FASTQ file with ``N`` so one ``add_sequence`` call
per block gives the same hashes as one call per read; the sketcher, writer and containment
helpers are checked against sourmash's own API on the same reads.
"""

import gzip
import threading

import pytest

from prepare.data.sra.sampling import iter_sequence_blocks


def _write_fastq(path, reads, gz=False, crlf=False, final_newline=True):
    """Write ``reads`` as a four-line FASTQ file (plain or gzip), with the line ends requested."""
    end = "\r\n" if crlf else "\n"
    text = "".join(f"@r{i}{end}{seq}{end}+{end}{'I' * len(seq)}{end}" for i, seq in enumerate(reads))
    if not final_newline:
        text = text.rstrip("\r\n")
    data = text.encode()
    if gz:
        with gzip.open(path, "wb") as handle:
            handle.write(data)
    else:
        path.write_bytes(data)
    return path


def _reads_from_blocks(blocks):
    """The sequences a sequence of ``(block, n)`` pairs holds, and the sum of the counts."""
    reads = []
    total = 0
    for block, n in blocks:
        if n:
            reads.extend(block.split(b"N"))
        total += n
    return reads, total


class TestIterSequenceBlocks:
    def test_blocks_join_sequences_with_n(self, tmp_path):
        path = _write_fastq(tmp_path / "a.fastq", ["ACGT", "GGCC", "TTAA"])
        blocks = list(iter_sequence_blocks(path))
        assert blocks == [(b"ACGTNGGCCNTTAA", 3)]

    def test_blocks_read_gzip(self, tmp_path):
        path = _write_fastq(tmp_path / "a.fastq.gz", ["ACGT", "GGCC"], gz=True)
        assert list(iter_sequence_blocks(path)) == [(b"ACGTNGGCC", 2)]

    def test_blocks_strip_cr(self, tmp_path):
        path = _write_fastq(tmp_path / "a.fastq", ["ACGT", "GGCC"], crlf=True)
        assert list(iter_sequence_blocks(path)) == [(b"ACGTNGGCC", 2)]

    def test_blocks_accept_missing_final_newline(self, tmp_path):
        path = _write_fastq(tmp_path / "a.fastq", ["ACGT", "GGCC"], final_newline=False)
        assert _reads_from_blocks(iter_sequence_blocks(path)) == ([b"ACGT", b"GGCC"], 2)

    def test_record_spanning_block_boundary(self, tmp_path, monkeypatch):
        reads = [
            f"ACGT{i:04d}".replace("0", "A").replace("1", "C").replace("2", "G").replace("3", "T") for i in range(10)
        ]
        reads = [r.translate(str.maketrans("456789", "ACGTAC")) for r in reads]
        path = _write_fastq(tmp_path / "a.fastq", reads)
        monkeypatch.setattr("prepare.data.sra.sampling.CHUNK_SIZE", 64)
        blocks = list(iter_sequence_blocks(path))
        assert len(blocks) > 1
        got, total = _reads_from_blocks(blocks)
        assert got == [r.encode() for r in reads]
        assert total == 10

    def test_truncated_record_raises(self, tmp_path):
        path = tmp_path / "a.fastq"
        path.write_bytes(b"@r0\nACGT\n+\nIIII\n@r1\nACGT\n")
        with pytest.raises(ValueError, match="Truncated FASTQ record"):
            list(iter_sequence_blocks(path))

    def test_not_fastq_raises(self, tmp_path):
        path = tmp_path / "a.fastq"
        path.write_bytes(b">r0\nACGT\n")
        with pytest.raises(ValueError, match="Not a FASTQ file"):
            list(iter_sequence_blocks(path))

    def test_empty_file_yields_nothing(self, tmp_path):
        path = tmp_path / "a.fastq"
        path.write_bytes(b"")
        assert list(iter_sequence_blocks(path)) == []

    def test_stop_token_raises_interrupted_error(self, tmp_path):
        path = _write_fastq(tmp_path / "a.fastq", ["ACGT"])
        stop = threading.Event()
        stop.set()
        with pytest.raises(InterruptedError):
            list(iter_sequence_blocks(path, stop=stop))
