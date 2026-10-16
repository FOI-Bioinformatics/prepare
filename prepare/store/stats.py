"""
One cached FASTQ statistics record per dataset.

``compute_dataset_stats`` computes the record that ``sra_profile`` and ``sra_report`` take a
dataset's read and base totals and GC content from (through
``prepare.sra.dataset_stats.load_dataset_stats``); ``cached_stats`` and ``store_stats`` are
the sidecar-backed cache around it (see ``prepare.store.sidecar.Sidecar.stats``), so a
dataset already profiled once is not re-read for the next command.

Read counts and (when available) bases/min/max/avg length come from an exact pass: a
streaming newline count (``prepare.data.sra.count_fastq_reads``), or ``seqkit stats -T``
when the tool is installed, which is both exact and faster. Everything that needs per-read
content (GC, quality, complexity, length distribution, N50) comes from a uniform sample over
the first file (``prepare.data.sra.sample_records``, which draws record indices from the
exact count and slices only those records out of the stream), so a large dataset is never
parsed record by record in Python.
"""

import logging
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from prepare.core.constants import DEFAULT_NUM_THREADS, SRA_ACCESSION_PATTERN
from prepare.core.exceptions import DataAccessError
from prepare.data.sra import count_fastq_reads, fastq_files, sample_records
from prepare.data.sra.sampling import distribution_from_histogram, quality_histogram
from prepare.store.layout import StorePaths, sidecar_path as dataset_sidecar_path, store_paths
from prepare.store.locks import dataset_lock
from prepare.store.sidecar import read_sidecar, write_sidecar
from prepare.utils import resources
from prepare.utils.security import SecureSubprocess

logger = logging.getLogger(__name__)

# Default sample size for the quality/GC/complexity part of the stats record.
DEFAULT_SAMPLE_SIZE = 10000

# How long ``store_stats`` waits for the dataset lock before giving up on the cache. A
# download or adoption holding the lock runs for minutes to hours, and the statistics record
# is only a cache, so waiting that out would stall an analysis for no benefit.
STATS_LOCK_WAIT_SECONDS = 30.0

# Read-length histogram buckets (inclusive bounds), matching the convention already used by
# prepare.sra.analytics's per-read length distribution; a length past the last bound falls
# into "1000+".
_LENGTH_BUCKETS = [(0, 50), (51, 100), (101, 150), (151, 250), (251, 500), (501, 1000)]

_REQUIRED_SEQKIT_COLUMNS = ("num_seqs", "sum_len", "min_len", "avg_len", "max_len")


def _file_signature(files: List[Path]) -> Dict[str, List[Union[int, float]]]:
    """``{file name: [size bytes, mtime]}`` snapshot used to detect a changed dataset."""
    signature: Dict[str, List[Union[int, float]]] = {}
    for f in files:
        st = f.stat()
        signature[f.name] = [st.st_size, st.st_mtime]
    return signature


def _sample_reads(path: Path, sample_size: int, reads_in_file: Optional[int]) -> Tuple[List[str], List[bytes], int]:
    """Sample up to ``sample_size`` (sequence, quality) pairs uniformly from ``path``.

    Returns ``(sequences, qualities, reads_in_file)``. ``reads_in_file`` is the count the
    caller already has (counted here when None); the record indices are drawn from it with a
    fixed seed, so repeated runs on the same file give the same sample.
    """
    if reads_in_file is None:
        reads_in_file = count_fastq_reads(path)
    records = sample_records([path], sample_size, total_records=reads_in_file, seed=0)
    return [seq.decode("utf-8") for seq, _qual in records], [qual for _seq, qual in records], reads_in_file


def _gc_fraction(seq: str) -> float:
    """Fraction of ``seq`` that is G or C; 0.0 for an empty sequence."""
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _n50(lengths: List[int]) -> int:
    """N50 of ``lengths``: the length at which the cumulative sum first reaches half the total."""
    if not lengths:
        return 0
    sorted_lengths = sorted(lengths, reverse=True)
    target = sum(sorted_lengths) / 2
    cumulative = 0
    for length in sorted_lengths:
        cumulative += length
        if cumulative >= target:
            return length
    return 0


def _length_histogram(lengths: List[int]) -> Dict[str, int]:
    """Bucket read lengths into named ranges, e.g. ``{"101-150": 42, ...}``."""
    histogram = {f"{lo}-{hi}": 0 for lo, hi in _LENGTH_BUCKETS}
    histogram["1000+"] = 0
    for length in lengths:
        for lo, hi in _LENGTH_BUCKETS:
            if lo <= length <= hi:
                histogram[f"{lo}-{hi}"] += 1
                break
        else:
            histogram["1000+"] += 1
    return histogram


def _quality_summary(qualities: List[bytes]) -> Dict[str, float]:
    """Mean/median/q25/q75 of the per-base Phred+33 quality scores of the sampled reads.

    Computed from a histogram of the scores (``quality_histogram``), not a per-base list.
    """
    figures = distribution_from_histogram(quality_histogram(qualities))
    return {key: figures[key] for key in ("mean", "median", "q25", "q75")}


def _parse_seqkit_table(stdout: str, files: List[Path]) -> Dict[str, Dict[str, Any]]:
    """Parse ``seqkit stats -T`` tab-delimited output into ``{file name: {column: value}}``.

    Columns are located by name in the header row rather than by fixed position, since a
    seqkit version may add columns; a comma in a numeric field (some builds format large
    counts that way even in ``-T`` mode) is stripped before conversion. Rows are matched to
    the input files by the ``file`` column rather than by output order, which ``-j`` makes no
    promises about. Raises ``ValueError`` when the expected columns or row count are not
    found, or when a row names a file that was not asked for, which the caller treats as a
    seqkit failure and falls back to the streaming count.
    """
    lines = [line for line in stdout.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError("seqkit stats produced no data rows")

    header = lines[0].split("\t")
    index = {name: i for i, name in enumerate(header)}
    if not all(name in index for name in (*_REQUIRED_SEQKIT_COLUMNS, "file")):
        raise ValueError(f"seqkit stats output missing expected columns: {header}")

    data_lines = lines[1:]
    if len(data_lines) != len(files):
        raise ValueError(f"seqkit stats row count ({len(data_lines)}) did not match input file count ({len(files)})")

    by_name = {file_path.name: file_path for file_path in files}
    rows: Dict[str, Dict[str, Any]] = {}
    for line in data_lines:
        cols = line.split("\t")
        name = Path(cols[index["file"]]).name
        if name not in by_name:
            raise ValueError(f"seqkit stats reported a file that was not requested: {name}")
        rows[name] = {
            "num_seqs": int(cols[index["num_seqs"]].replace(",", "")),
            "sum_len": int(cols[index["sum_len"]].replace(",", "")),
            "min_len": int(cols[index["min_len"]].replace(",", "")),
            "avg_len": float(cols[index["avg_len"]].replace(",", "")),
            "max_len": int(cols[index["max_len"]].replace(",", "")),
        }
    return rows


def _run_seqkit_stats(files: List[Path]) -> Dict[str, Dict[str, Any]]:
    """Run ``seqkit stats -T -j <threads> <files...>`` and parse its output table."""
    threads = min(DEFAULT_NUM_THREADS, resources.available_cpus())
    args = ["stats", "-T", "-j", str(threads), *(str(f) for f in files)]
    result = SecureSubprocess.run_secure("seqkit", args)
    return _parse_seqkit_table(result.stdout or "", files)


def _try_seqkit_stats(file_paths: List[Path], use_seqkit: bool) -> Optional[Dict[str, Dict[str, Any]]]:
    """``_run_seqkit_stats`` when requested and installed, else None; a failure logs and falls back."""
    if not (use_seqkit and shutil.which("seqkit")):
        return None
    try:
        return _run_seqkit_stats(file_paths)
    except Exception as e:
        logger.warning("seqkit stats failed for %s, falling back to streaming counts: %s", file_paths, e)
        return None


def _counts_from_seqkit(seqkit_rows: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Exact reads/bases/length figures straight from a parsed seqkit stats table."""
    reads_per_file = {name: row["num_seqs"] for name, row in seqkit_rows.items()}
    reads_total = sum(reads_per_file.values())
    bases_total = sum(row["sum_len"] for row in seqkit_rows.values())
    return {
        "reads_per_file": reads_per_file,
        "reads_total": reads_total,
        "bases_total": bases_total,
        "min_read_length": min(row["min_len"] for row in seqkit_rows.values()),
        "max_read_length": max(row["max_len"] for row in seqkit_rows.values()),
        "avg_read_length": (bases_total / reads_total) if reads_total else 0.0,
    }


def _counts_from_streaming(file_paths: List[Path]) -> Dict[str, Any]:
    """Exact read counts from a streaming newline count; length/base figures filled in later
    from the sample, since a full pass just to measure lengths would defeat the
    point of not using seqkit."""
    reads_per_file = {f.name: count_fastq_reads(f) for f in file_paths}
    return {
        "reads_per_file": reads_per_file,
        "reads_total": sum(reads_per_file.values()),
        "bases_total": None,
        "min_read_length": None,
        "max_read_length": None,
        "avg_read_length": None,
    }


def _approximate_lengths_from_sample(counts: Dict[str, Any], sample_lengths: List[int]) -> Dict[str, Any]:
    """Fill in length/base figures from the sample when there is no exact count."""
    avg_read_length = _mean(sample_lengths)
    counts["avg_read_length"] = avg_read_length
    counts["min_read_length"] = min(sample_lengths) if sample_lengths else 0
    counts["max_read_length"] = max(sample_lengths) if sample_lengths else 0
    counts["bases_total"] = int(round(avg_read_length * counts["reads_total"]))
    return counts


def compute_dataset_stats(
    files: List[Union[str, Path]], sample_size: int = DEFAULT_SAMPLE_SIZE, use_seqkit: bool = True
) -> Dict[str, Any]:
    """Compute the shared statistics record for one dataset's FASTQ files.

    ``reads_per_file``/``reads_total`` (and, when available, ``bases_total`` and
    ``min``/``avg``/``max_read_length``) are exact: from ``seqkit stats -T`` when
    ``use_seqkit`` and the tool is installed (a seqkit failure logs a warning and falls back
    rather than raising), otherwise from a streaming newline count, with the length figures
    then approximated from the sample below. Everything that needs per-read content
    (``gc_content``, ``n_content``, ``length_histogram``, ``n50``, ``quality_summary``,
    ``duplication_rate``) comes from a uniform sample of ``sample_size`` reads over the first
    file, drawn by record index from the count above; ``sampled`` is True when that file holds
    more reads than the sample.

    ``signature`` records each file's size and mtime so ``cached_stats`` can tell whether the
    dataset on disk still matches a cached copy of this result.
    """
    if not files:
        raise ValueError("compute_dataset_stats requires at least one FASTQ file")
    file_paths = [Path(f) for f in files]
    signature = _file_signature(file_paths)

    seqkit_rows = _try_seqkit_stats(file_paths, use_seqkit)
    counts = _counts_from_seqkit(seqkit_rows) if seqkit_rows is not None else _counts_from_streaming(file_paths)

    sequences, qualities, seen_in_sample_file = _sample_reads(
        file_paths[0], sample_size, counts["reads_per_file"].get(file_paths[0].name)
    )
    sampled = seen_in_sample_file > sample_size
    sample_lengths = [len(s) for s in sequences]

    if seqkit_rows is None:
        counts = _approximate_lengths_from_sample(counts, sample_lengths)

    return {
        **counts,
        "n50": _n50(sample_lengths),
        "gc_content": _mean([_gc_fraction(s) for s in sequences]),
        "length_histogram": _length_histogram(sample_lengths),
        "quality_summary": _quality_summary(qualities),
        "duplication_rate": (1.0 - len(set(sequences)) / len(sequences)) if sequences else 0.0,
        "n_content": _mean([(s.count("N") / len(s)) if s else 0.0 for s in sequences]),
        "sample_size": sample_size,
        "sampled": sampled,
        "computed": datetime.now(timezone.utc).isoformat(),
        "signature": signature,
    }


def cached_stats(acc_dir: Union[str, Path], sidecar_path: Optional[Union[str, Path]]) -> Optional[Dict[str, Any]]:
    """The sidecar's cached stats when its recorded signature matches the files in ``acc_dir`` now.

    Returns None when there is no sidecar path, no readable sidecar, no stats recorded yet,
    or any file's size/mtime has changed since the stats were computed (a file added,
    removed, resized, or rewritten in place).
    """
    if sidecar_path is None:
        return None
    sidecar = read_sidecar(sidecar_path)
    if sidecar is None or not sidecar.stats:
        return None
    signature = sidecar.stats.get("signature")
    if not signature:
        return None
    current = _file_signature(fastq_files(acc_dir))
    if current != signature:
        return None
    return sidecar.stats


def _store_location(path: Path) -> Optional[Tuple[StorePaths, str]]:
    """``(StorePaths, accession)`` when ``path`` is a store sidecar, else None.

    A store sidecar lives at ``<root>/sra/<ACC>/<ACC>.json``. The candidate root is taken
    from the path's own position and the layout helpers are asked where that store would
    keep this accession's sidecar; only an exact match is treated as a store, so a sidecar
    written beside a plain project folder is left to the unlocked write below.
    """
    acc_dir = path.parent
    accession = acc_dir.name
    if not re.fullmatch(SRA_ACCESSION_PATTERN, accession):
        return None
    paths = store_paths(acc_dir.parent.parent)
    if dataset_sidecar_path(paths, accession) != path:
        return None
    return paths, accession


def _write_stats(sidecar_path: Union[str, Path], stats: Dict[str, Any]) -> None:
    """Re-read the sidecar and write back ``stats`` and a fresh ``stats_computed``.

    The sidecar is read here rather than reused from the caller because computing the
    statistics takes seconds to minutes on a real dataset, during which another project may
    have published a new state or file list for the same accession. Only the two statistics
    fields are assigned, so everything else keeps whatever is on disk now.
    """
    sidecar = read_sidecar(sidecar_path)
    if sidecar is None:
        return
    sidecar.stats = stats
    sidecar.stats_computed = datetime.now(timezone.utc).isoformat()
    write_sidecar(sidecar_path, sidecar)


def store_stats(sidecar_path: Optional[Union[str, Path]], stats: Dict[str, Any]) -> None:
    """Write ``stats`` into the sidecar at ``sidecar_path``, with a fresh ``stats_computed``.

    A no-op when there is no sidecar to update (``sidecar_path`` is None, or names a sidecar
    that does not exist) -- a project without a shared store then keeps the stats only in its
    own analysis record rather than in a sidecar that was never created.

    For a sidecar inside a store, the read-modify-write runs under the accession's dataset
    lock, the same lock every other sidecar writer takes (download, adoption,
    ``store_verify --fix-state``), so a statistics write cannot overwrite a state change
    another project published meanwhile. The wait is bounded: a dataset locked by a running
    download only loses the cache, which is recomputed on the next command, so the write is
    skipped with a warning rather than blocking the analysis that asked for it.
    """
    if sidecar_path is None:
        return
    located = _store_location(Path(sidecar_path))
    if located is None:
        _write_stats(sidecar_path, stats)
        return
    paths, accession = located
    try:
        with dataset_lock(paths, accession, wait_seconds=STATS_LOCK_WAIT_SECONDS):
            _write_stats(sidecar_path, stats)
    except DataAccessError as e:
        logger.warning("Not caching statistics for %s: %s", accession, e)
