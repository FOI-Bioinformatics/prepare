"""
Index-based sampling of FASTQ records, and per-base quality scores kept as histograms.

``sample_records`` is the one sampler behind the quality profile
(``prepare.sra.quality.SequenceQualityAnalyzer``) and the store statistics record
(``prepare.store.stats.compute_dataset_stats``). It draws the record indices to keep up
front, from the dataset's record count, and then streams the files once in large binary
blocks, locating records by counting newlines and slicing out only the selected ones. No
per-record Python work is done for the records that are skipped, which is what made the
earlier reservoir sampler slow on a large run.

The sampler assumes the four-line FASTQ layout that ``fasterq-dump`` and ``count_fastq_reads``
assume: header, sequence, ``+`` separator and quality, each on exactly one line.
"""

import contextlib
import gzip
import io
import random
import threading
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union

import numpy as np

from prepare.data.sra.fastq import count_fastq_reads

# Bytes decompressed per read call. 1 MiB was the fastest of 256 KiB to 4 MiB on synthetic
# files. Larger reads from gzip.GzipFile are slower because CPython assembles them from smaller
# decompressed pieces, and on a real 427 MB mate file 4 MiB blocks spent most of the time copying.
CHUNK_SIZE = 1024 * 1024

# Read buffer for the compressed input of a ``.gz`` file. gzip.GzipFile otherwise reads its input
# in pieces of about 128 KiB, which on a real run's external volume meant thousands of small reads.
COMPRESSED_READ_BUFFER = 8 * 1024 * 1024

# Phred scores 0-93 (``!`` to ``~`` in Phred+33), the range of the quality histograms.
QUALITY_BINS = 94
PHRED_OFFSET = 33

Record = Tuple[bytes, bytes]


@contextlib.contextmanager
def _open_binary(path: Path) -> Iterator[io.BufferedIOBase]:
    """Open ``path`` for binary reading, decompressing a ``.gz`` file read through a large buffer."""
    if not str(path).endswith(".gz"):
        with open(path, "rb") as handle:
            yield handle
        return
    with open(path, "rb", buffering=COMPRESSED_READ_BUFFER) as raw, gzip.GzipFile(fileobj=raw, mode="rb") as handle:
        yield handle


def _strip_cr(line: bytes) -> bytes:
    """``line`` without a trailing carriage return (a file written with CRLF line ends)."""
    return line[:-1] if line.endswith(b"\r") else line


def _take_records(
    buf: bytes,
    newline_at: np.ndarray,
    wanted: List[int],
    pos: int,
    base: int,
    out: List[Record],
    path: Path,
    file_start: int,
) -> int:
    """Append the selected complete records in ``buf`` to ``out``; return the next ``wanted`` position.

    ``base`` is the dataset-wide index of the first record in ``buf``, which starts on a record
    boundary, and ``newline_at`` holds the offsets of every newline in ``buf``. ``file_start`` is
    the dataset-wide index of the first record of ``path``, used to report a malformed record by
    its 1-based number within the file.
    """
    last = base + len(newline_at) // 4
    while pos < len(wanted) and wanted[pos] < last:
        line = 4 * (wanted[pos] - base)
        header_start = int(newline_at[line - 1]) + 1 if line else 0
        seq_start, seq_end = int(newline_at[line]) + 1, int(newline_at[line + 1])
        qual_start, qual_end = int(newline_at[line + 2]) + 1, int(newline_at[line + 3])
        if buf[header_start : header_start + 1] != b"@" or buf[seq_end + 1 : seq_end + 2] != b"+":
            header = buf[header_start : seq_start - 1].strip()
            number = wanted[pos] - file_start + 1
            raise ValueError(f"Malformed FASTQ record {number} of {path}: {header!r}")
        out.append((_strip_cr(buf[seq_start:seq_end]), _strip_cr(buf[qual_start:qual_end])))
        pos += 1
    return pos


def _sample_file(path: Path, wanted: List[int], pos: int, base: int, out: List[Record]) -> Tuple[int, int]:
    """Collect the records of ``path`` whose dataset-wide index is in ``wanted[pos:]``.

    Returns the next position in ``wanted`` and the dataset-wide index just past the last
    record read, so the next file continues the numbering. Stops reading as soon as every
    wanted index has been taken. A file that ends before an index it was expected to hold
    simply ends; a trailing record cut short (fewer than four lines) raises ``ValueError``
    only when it is one of the selected records.
    """
    carry = b""
    checked_first = False
    file_start = base
    with _open_binary(path) as handle:
        while pos < len(wanted):
            chunk = handle.read(CHUNK_SIZE)
            at_eof = not chunk
            buf = carry + chunk if carry else chunk
            if at_eof and buf and not buf.endswith(b"\n"):
                buf += b"\n"
            if not checked_first and buf:
                if not buf.startswith(b"@"):
                    raise ValueError(f"Not a FASTQ file (the first record does not start with '@'): {path}")
                checked_first = True
            # One vectorised pass gives both the record count and every line boundary; it is
            # cheaper than bytes.count followed by a second pass over blocks holding a sample.
            newline_at = np.flatnonzero(np.frombuffer(buf, dtype=np.uint8) == 10)
            complete = len(newline_at) // 4
            cut = int(newline_at[4 * complete - 1]) + 1 if complete else 0
            if complete and wanted[pos] < base + complete:
                pos = _take_records(buf, newline_at, wanted, pos, base, out, path, file_start)
            base += complete
            carry = buf[cut:]
            if at_eof:
                if carry.strip() and pos < len(wanted) and wanted[pos] == base:
                    header = carry.split(b"\n", 1)[0].strip()
                    raise ValueError(f"Truncated FASTQ record after header: {header!r} in {path}")
                break
    return pos, base


def sample_records(
    paths: Sequence[Union[str, Path]],
    sample_size: int,
    total_records: Optional[int] = None,
    seed: int = 0,
) -> List[Record]:
    """Draw up to ``sample_size`` records uniformly from the concatenation of ``paths``.

    The records of all files are numbered in order (mate 1, then mate 2), and
    ``k = min(sample_size, total_records)`` distinct indices are drawn with
    ``random.Random(seed).sample(range(total_records), k)``, so every record has the same
    chance of being included and a given seed and total always select the same records.
    Returns ``(sequence, quality)`` byte strings in file order, without line ends.

    Args:
        paths: The dataset's FASTQ files, plain or gzip-compressed
        sample_size: Number of records to return at most
        total_records: Records across all ``paths``, e.g. from the cached statistics
            record; counted here with ``count_fastq_reads`` (one extra pass) when None. A
            total larger than the files hold is tolerated: indices past the end are
            dropped, so fewer records are returned, and with a much overstated total
            possibly none at all. A total smaller than the true count
            makes only the first ``total_records`` records eligible, so the later ones are
            never sampled.
        seed: Seed of the index draw

    Raises:
        ValueError: when a file does not start with ``@``, or a selected record is
            malformed (no ``+`` separator line) or cut short at the end of the file
    """
    file_paths = [Path(p) for p in paths]
    if sample_size <= 0 or not file_paths:
        return []
    if total_records is None:
        total_records = sum(count_fastq_reads(p) for p in file_paths)
    if total_records <= 0:
        return []
    k = min(sample_size, total_records)
    wanted = sorted(random.Random(seed).sample(range(total_records), k))

    out: List[Record] = []
    pos = base = 0
    for path in file_paths:
        if pos >= len(wanted):
            break
        pos, base = _sample_file(path, wanted, pos, base, out)
    return out


def iter_sequence_blocks(path: Union[str, Path], stop: Optional[threading.Event] = None) -> Iterator[Tuple[bytes, int]]:
    """Yield ``(sequences joined by b"N", record count)`` for each block of the FASTQ file ``path``.

    The file is read in ``CHUNK_SIZE`` blocks (gzip-compressed or plain, through
    ``_open_binary``); the sequence line of every complete record in a block is sliced out and
    the slices are joined with a single ``N``, so a sketcher can hash a whole block with one
    ``add_sequence(block, force=True)`` call: every k-mer that spans a junction contains the
    ``N`` and is skipped, which leaves exactly the k-mers of the individual reads. The partial
    record at the end of a block is carried into the next one. A trailing carriage return is
    removed from each sequence; a file without a final newline is accepted.

    Raises:
        ValueError: when the file does not start with ``@`` or ends in a cut-short record
        InterruptedError: when ``stop`` is set before a block is read
    """
    file_path = Path(path)
    carry = b""
    checked_first = False
    with _open_binary(file_path) as handle:
        while True:
            if stop is not None and stop.is_set():
                raise InterruptedError(f"Reading {file_path} stopped: the run was interrupted")
            chunk = handle.read(CHUNK_SIZE)
            at_eof = not chunk
            buf = carry + chunk if carry else chunk
            if at_eof and buf and not buf.endswith(b"\n"):
                buf += b"\n"
            if not checked_first and buf:
                if not buf.startswith(b"@"):
                    raise ValueError(f"Not a FASTQ file (the first record does not start with '@'): {file_path}")
                checked_first = True
            newline_at = np.flatnonzero(np.frombuffer(buf, dtype=np.uint8) == 10)
            complete = len(newline_at) // 4
            if complete:
                yield _joined_sequences(buf, newline_at, complete), complete
            cut = int(newline_at[4 * complete - 1]) + 1 if complete else 0
            carry = buf[cut:]
            if at_eof:
                if carry.strip():
                    header = carry.split(b"\n", 1)[0].strip()
                    raise ValueError(f"Truncated FASTQ record after header: {header!r} in {file_path}")
                return


def _joined_sequences(buf: bytes, newline_at: np.ndarray, complete: int) -> bytes:
    """The sequence lines of the first ``complete`` records in ``buf``, joined with ``N``."""
    starts = newline_at[0 : 4 * complete : 4] + 1
    ends = newline_at[1 : 4 * complete : 4]
    if ends.size and buf[int(ends[0]) - 1 : int(ends[0])] == b"\r":
        ends = ends - 1
    return b"N".join(buf[int(start) : int(end)] for start, end in zip(starts, ends))


def _fold_raw_counts(raw: np.ndarray) -> np.ndarray:
    """Fold a 256-bin byte-value count into Phred bins 0-93, clamping values outside ``!``-``~``."""
    hist = raw[PHRED_OFFSET : PHRED_OFFSET + QUALITY_BINS].astype(np.int64)
    hist[0] += int(raw[:PHRED_OFFSET].sum())
    hist[-1] += int(raw[PHRED_OFFSET + QUALITY_BINS :].sum())
    return hist


def quality_histogram(quals: Iterable[bytes]) -> np.ndarray:
    """Count of each Phred+33 quality score (0-93) over all quality strings in ``quals``.

    Returns an int64 array of length 94; a character below ``!`` counts as 0 and one above
    ``~`` as 93.
    """
    data = b"".join(quals)
    raw = np.bincount(np.frombuffer(data, dtype=np.uint8), minlength=256)
    return _fold_raw_counts(raw)


def histogram_from_scores(scores: Iterable[int]) -> np.ndarray:
    """The ``quality_histogram`` of already-decoded Phred scores (clamped to 0-93)."""
    values = np.fromiter(scores, dtype=np.int64)
    return np.bincount(np.clip(values, 0, QUALITY_BINS - 1), minlength=QUALITY_BINS).astype(np.int64)


def _score_at_rank(cumulative: np.ndarray, rank: int) -> int:
    """The score of the ``rank``-th (0-based) value in sorted order, from cumulative bin counts."""
    return int(np.searchsorted(cumulative, rank, side="right"))


def _histogram_percentile(cumulative: np.ndarray, total: int, pct: float) -> float:
    """Linear-interpolated percentile (numpy's default method) of the values a histogram counts."""
    k = (total - 1) * (pct / 100)
    lo = int(k)
    hi = min(lo + 1, total - 1)
    v_lo = _score_at_rank(cumulative, lo)
    if hi == lo:
        return float(v_lo)
    v_hi = _score_at_rank(cumulative, hi)
    return float(v_lo + (v_hi - v_lo) * (k - lo))


def distribution_from_histogram(hist: np.ndarray) -> Dict[str, float]:
    """Summary figures of the per-base quality scores a ``quality_histogram`` counts.

    Returns ``mean``, ``median``, ``q25`` and ``q75`` (linear-interpolated percentiles, as
    ``numpy.percentile``) and the fractions of bases at Q30 or more (``excellent_q30+``),
    Q20-29 (``good_q20-29``), Q10-19 (``fair_q10-19``) and below Q10 (``poor_q0-9``). All
    figures are 0.0 for an empty histogram.
    """
    counts = np.asarray(hist, dtype=np.int64)
    total = int(counts.sum())
    keys = ("mean", "median", "q25", "q75", "excellent_q30+", "good_q20-29", "fair_q10-19", "poor_q0-9")
    if total == 0:
        return {key: 0.0 for key in keys}
    cumulative = np.cumsum(counts)
    return {
        "mean": int(np.dot(np.arange(len(counts), dtype=np.int64), counts)) / total,
        "median": _histogram_percentile(cumulative, total, 50),
        "q25": _histogram_percentile(cumulative, total, 25),
        "q75": _histogram_percentile(cumulative, total, 75),
        "excellent_q30+": int(counts[30:].sum()) / total,
        "good_q20-29": int(counts[20:30].sum()) / total,
        "fair_q10-19": int(counts[10:20].sum()) / total,
        "poor_q0-9": int(counts[:10].sum()) / total,
    }
