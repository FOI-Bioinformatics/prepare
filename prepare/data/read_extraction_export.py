"""Mapped-read export and coverage helpers for targeted read extraction.

This module holds the samtools and minimap2 argument builders, the export of a filtered BAM back
to FASTQ, the per-sample reference coverage table and the coverage estimate for an assembly.
``prepare.data.read_extraction`` imports and re-exports these names, so existing imports keep working.
"""

import gzip
import logging
import re
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import pandas as pd

from prepare.core.exceptions import ConfigurationError, DataAccessError, SecurityError
from prepare.data.assembly import is_singleton_file, summarise_contigs
from prepare.utils.security import SecureSubprocess

# The same logger as prepare.data.read_extraction: the code here ran there before the split, and the
# command layer filters and the tests capture that logger by name.
logger = logging.getLogger("prepare.data.read_extraction")

# samtools -F flags dropped by the alignment filter: unmapped (0x4), secondary (0x100)
# and supplementary (0x800) alignments.
FILTER_FLAGS = "0x904"

# minimap2 prints this when the two mate files differ in length; it then maps single-end.
UNEQUAL_MATES_MARKER = "different number of records"

MINIMAP2_MIN_VERSION = "2.17"  # the first minimap2 release with --sam-hit-only (minimap2 NEWS.md)


def _count_records(path: Path, *filter_args: str) -> int:
    """Number of records in a SAM/BAM file, via ``samtools view -c``.

    ``filter_args`` are extra ``samtools view`` flags (e.g. ``"-F", "4"``) applied before
    the count; with none, every record in the file is counted.
    """
    result = SecureSubprocess.run_secure("samtools", _samtools_count_args(filter_args, path))
    text = (result.stdout or "").strip()
    return int(text) if text.isdigit() else 0


def _fastq_is_empty(path: Path) -> bool:
    """True when the file is missing or has no records (gzip or plain)."""
    if not path.exists():
        return True
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return handle.readline() == ""


def _filter_and_merge_bam(
    sam_paths: List[Path], filter_args: List[str], threads: int, out_dir: Path, genome_id: str, bam_path: Path
) -> None:
    """Filter each SAM into a BAM (``-F 0x904 [-q min_mapq]``); merge with ``samtools cat``
    when there is more than one (the deliberate per-mate single-end fallback)."""
    if len(sam_paths) == 1:
        SecureSubprocess.run_secure("samtools", _samtools_view_args(filter_args, threads, bam_path, sam_paths[0]))
        return
    part_bams = [out_dir / f"{genome_id}.mate{i}.bam" for i in range(1, len(sam_paths) + 1)]
    for sam_path, part_bam in zip(sam_paths, part_bams):
        SecureSubprocess.run_secure("samtools", _samtools_view_args(filter_args, threads, part_bam, sam_path))
    SecureSubprocess.run_secure("samtools", _samtools_cat_args(bam_path, part_bams))
    for part_bam in part_bams:
        part_bam.unlink(missing_ok=True)


def _export_mapped_fastq(reads: List[Path], bam_path: Path, out_dir: Path, genome_id: str, threads: int) -> List[Path]:
    """Export a filtered BAM back to FASTQ.

    A paired mapping gives the mate files ``_1``/``_2`` and, when some reads mapped without
    their mate, a singleton file ``_s`` that is part of the extraction. Reads without a mate
    flag (a single-end mapping, or the unequal-mates fallback) go to the ``_0`` file.
    """
    if len(reads) < 2:
        out0 = out_dir / f"{genome_id}.fastq.gz"
        SecureSubprocess.run_secure("samtools", _samtools_fastq_single_args(threads, out0, bam_path))
        return [] if _fastq_is_empty(out0) else [out0]

    out1 = out_dir / f"{genome_id}_1.fastq.gz"
    out2 = out_dir / f"{genome_id}_2.fastq.gz"
    singles = out_dir / f"{genome_id}_s.fastq.gz"
    orphans = out_dir / f"{genome_id}_0.fastq.gz"
    SecureSubprocess.run_secure(
        "samtools", _samtools_fastq_paired_args(threads, out1, out2, singles, orphans, bam_path)
    )
    for path in (out1, out2, singles, orphans):
        if _fastq_is_empty(path):
            path.unlink(missing_ok=True)
    if out1.exists() and out2.exists():
        return [out1, out2, singles] if singles.exists() else [out1, out2]
    if orphans.exists():
        return [orphans]
    return [singles] if singles.exists() else []


def coverage_table_path(out_dir: Path, genome_id: str) -> Path:
    """Where ``reference_coverage`` writes one sample's ``samtools coverage`` table."""
    return out_dir / f"{genome_id}_coverage.tsv"


def summarise_coverage_table(path: Union[str, Path]) -> Dict[str, Any]:
    """Aggregate a ``samtools coverage`` table over every reference sequence.

    The table counts the kept alignments of the filtered BAM; ``samtools coverage`` also skips
    duplicate and QC-fail reads by default (its ``--ff`` default). Breadth is the fraction of
    reference bases covered by at least one read (samtools coverage has no minimum-depth
    option, so >= 1x is the only threshold available); mean depth is each sequence's
    ``meandepth`` weighted by its length. Both are rounded to four decimals and are None when
    the table holds no reference bases.

    Returns:
        ``{"breadth", "mean_depth", "covered_bases", "reference_bp"}``.
    """
    table = pd.read_csv(path, sep="\t")
    lengths = table["endpos"] - table["startpos"] + 1
    reference_bp = int(lengths.sum())
    covered_bases = int(table["covbases"].sum())
    if reference_bp == 0:
        breadth: Optional[float] = None
        mean_depth: Optional[float] = None
    else:
        breadth = round(covered_bases / reference_bp, 4)
        mean_depth = round(float((table["meandepth"] * lengths).sum()) / reference_bp, 4)
    return {
        "breadth": breadth,
        "mean_depth": mean_depth,
        "covered_bases": covered_bases,
        "reference_bp": reference_bp,
    }


def reference_coverage(
    accession: str,
    genome_id: str,
    bam_path: Path,
    out_dir: Path,
    sam_root: Path,
    threads: int,
) -> Optional[Dict[str, Any]]:
    """Breadth and mean depth of the target genome covered by one sample's kept alignments.

    The filtered BAM is coordinate-sorted to ``sam_root/<genome>.mapped.sorted.bam`` (beside
    the SAM files, so ``--temp-folder`` applies) and ``samtools coverage`` writes the
    per-sequence table to ``out_dir/<genome>_coverage.tsv``. The sorted BAM is always removed.
    Only the kept alignments count (unmapped, secondary and supplementary records, and any
    below ``--min-mapq``, were filtered out earlier), and ``samtools coverage`` additionally
    skips duplicate and QC-fail reads by default. Coverage is supplementary to the extracted
    reads, so a tool or parsing failure is logged as a warning, any partial table is removed,
    and None is returned rather than failing the sample.

    Returns:
        The ``summarise_coverage_table`` summary plus ``coverage_tsv`` (the table's path), or
        None when the coverage could not be computed.
    """
    sorted_bam = sam_root / f"{genome_id}.mapped.sorted.bam"
    tsv_path = coverage_table_path(out_dir, genome_id)
    try:
        SecureSubprocess.run_secure("samtools", _samtools_sort_args(threads, sorted_bam, bam_path))
        SecureSubprocess.run_secure("samtools", _samtools_coverage_args(tsv_path, sorted_bam))
        summary = summarise_coverage_table(tsv_path)
    except (subprocess.CalledProcessError, SecurityError, DataAccessError, OSError, ValueError, KeyError) as exc:
        logger.warning("%s: reference coverage skipped (%s)", accession, exc)
        tsv_path.unlink(missing_ok=True)
        return None
    finally:
        sorted_bam.unlink(missing_ok=True)
    summary["coverage_tsv"] = tsv_path
    return summary


def _minimap2_index_args(preset: str, index_path: Path, genome_path: Path) -> List[str]:
    """Build the minimap2 index-build argument list used by ``build_index``."""
    return ["-x", preset, "-d", str(index_path), str(genome_path)]


def _minimap2_map_args(
    preset: str, threads: int, sam_path: Path, reference: Union[str, Path], reads: Sequence[Union[str, Path]]
) -> List[str]:
    """minimap2 mapping args for ``_run_minimap2``/``assembly_coverage``; ``--sam-hit-only`` skips unmapped reads."""
    args = ["-a", "-x", preset, "-t", str(threads), "-o", str(sam_path), str(reference), *(str(r) for r in reads)]
    return args + ["--sam-hit-only"]


def _samtools_count_args(filter_args: Sequence[str], path: Path) -> List[str]:
    """Build the ``samtools view -c`` argument list used by ``_count_records``."""
    return ["view", "-c", *filter_args, str(path)]


def _samtools_view_args(
    filter_args: Sequence[str], threads: int, out_path: Path, in_path: Union[str, Path]
) -> List[str]:
    """Build the ``samtools view -b`` filter argument list shared by ``_filter_and_merge_bam``
    (the single-SAM case and the per-mate loop) and ``assembly_coverage``."""
    return ["view", "-b", *filter_args, "-@", str(threads), "-o", str(out_path), str(in_path)]


def _samtools_cat_args(out_path: Path, part_paths: Sequence[Path]) -> List[str]:
    """Build the ``samtools cat`` argument list ``_filter_and_merge_bam`` uses to merge the
    per-mate BAMs of the unequal-mates single-end fallback."""
    return ["cat", "-o", str(out_path), *(str(p) for p in part_paths)]


def _samtools_sort_args(threads: int, out_path: Path, in_path: Path) -> List[str]:
    """Build the ``samtools sort`` argument list ``reference_coverage`` uses to coordinate-sort
    the filtered BAM before ``samtools coverage``."""
    return ["sort", "-@", str(threads), "-o", str(out_path), str(in_path)]


def _samtools_coverage_args(out_path: Path, bam_path: Path) -> List[str]:
    """Build the ``samtools coverage`` argument list ``reference_coverage`` uses to write the
    per-sequence coverage table of the sorted BAM."""
    return ["coverage", "-o", str(out_path), str(bam_path)]


def _samtools_fastq_single_args(threads: int, out_path: Path, bam_path: Path) -> List[str]:
    """Build the ``samtools fastq`` argument list ``_export_mapped_fastq`` uses for single-
    end (or already-collapsed) output."""
    return ["fastq", "-@", str(threads), "-0", str(out_path), str(bam_path)]


def _samtools_fastq_paired_args(
    threads: int, out1: Path, out2: Path, singles: Path, orphans: Path, bam_path: Path
) -> List[str]:
    """Build the ``samtools fastq`` argument list ``_export_mapped_fastq`` uses for paired
    output."""
    args = ["fastq", "-@", str(threads), "-1", str(out1), "-2", str(out2)]
    return args + ["-s", str(singles), "-0", str(orphans), str(bam_path)]


def _average_read_length(path: Union[str, Path], sample_size: int = 1000) -> float:
    """Average sequence length over the first ``sample_size`` records of a FASTQ file.

    Gzip aware; returns 0.0 for a file with no records.
    """
    opener = gzip.open if Path(path).suffix == ".gz" else open
    total = 0
    count = 0
    with opener(path, "rt") as handle:
        for i, line in enumerate(handle):
            if i % 4 == 1:
                total += len(line.strip())
                count += 1
                if count >= sample_size:
                    break
    return total / count if count else 0.0


def assembly_coverage(
    contigs: Union[str, Path],
    reads: Sequence[Union[str, Path]],
    preset: str,
    threads: int,
    work_dir: Union[str, Path],
    mapped_reads: Optional[int],
) -> Dict[str, Any]:
    """Estimate assembly coverage by mapping the extracted reads back onto its contigs.

    Aligns ``reads`` against ``contigs`` with minimap2, filters unmapped/secondary/supplementary
    alignments the same way extraction does, and counts what remains. ``mapped_reads`` is the
    sample's mapped-read count from extraction, used for what fraction of those reads the assembly
    recruits; pass ``None`` (or 0) when that count is not known, and ``mapping_rate`` comes back
    ``None``. The SAM and BAM written under ``work_dir`` are removed before returning either way.

    Returns:
        ``{"reads_mapped": int, "mapping_rate": Optional[float], "mean_depth_estimate":
        Optional[float]}``, the last None when the assembly has no contigs to spread reads over.
    """
    work_root = Path(work_dir)
    reads_mapped = sum(_count_back_mapped(contigs, group, preset, threads, work_root) for group in _query_groups(reads))

    mapping_rate = reads_mapped / mapped_reads if mapped_reads else None
    avg_read_len = _average_read_length(reads[0]) if reads else 0.0
    total_bp = summarise_contigs(contigs)["total_bp"]
    # An assembly with no contigs has no depth to report; 0.0 would read as a measured
    # depth of zero, which is a different statement. ``mapping_rate`` says None for the
    # same reason when the mapped-read count is unknown.
    mean_depth_estimate = (reads_mapped * avg_read_len / total_bp) if total_bp else None

    return {
        "reads_mapped": reads_mapped,
        "mapping_rate": mapping_rate,
        "mean_depth_estimate": mean_depth_estimate,
    }


def _query_groups(reads: Sequence[Union[str, Path]]) -> List[List[Path]]:
    """The read files in the groups minimap2 maps together: the mate files, then the singleton file.

    minimap2 pairs its query files by position, so a singleton file next to the mate files is
    mapped in a run of its own rather than paired with a mate file.
    """
    paths = [Path(r) for r in reads]
    if len(paths) <= 2:
        return [paths]
    mates = [r for r in paths if not is_singleton_file(r)]
    singles = [r for r in paths if is_singleton_file(r)]
    return [group for group in (mates, singles) if group]


def _count_back_mapped(
    contigs: Union[str, Path], queries: Sequence[Union[str, Path]], preset: str, threads: int, work_root: Path
) -> int:
    """Map ``queries`` onto ``contigs`` and count the kept records; the SAM and BAM are removed either way."""
    sam_path = work_root / "coverage.sam"
    bam_path = work_root / "coverage.bam"
    try:
        SecureSubprocess.run_secure("minimap2", _minimap2_map_args(preset, threads, sam_path, contigs, queries))
        SecureSubprocess.run_secure("samtools", _samtools_view_args(["-F", FILTER_FLAGS], threads, bam_path, sam_path))
        return _count_records(bam_path)
    finally:
        sam_path.unlink(missing_ok=True)
        bam_path.unlink(missing_ok=True)


def _run_minimap2(
    accession: str,
    preset: str,
    threads: int,
    sam_path: Path,
    reference: Path,
    genome_fasta: Path,
    reads: List[Path],
) -> subprocess.CompletedProcess:
    """Align ``reads`` against the prebuilt index, retrying once against the FASTA directly
    if the index fails (e.g. it was built by an incompatible minimap2 version).

    Only the failures a tool run can produce are retried: a non-zero exit
    (``CalledProcessError``), a rejected or timed-out command (``SecurityError``), and a
    filesystem error. A ``KeyboardInterrupt`` or a programming error is not a reason to run
    minimap2 a second time, so it propagates.
    """
    args = _minimap2_map_args(preset, threads, sam_path, reference, reads)
    try:
        return SecureSubprocess.run_secure("minimap2", args)
    except (subprocess.CalledProcessError, SecurityError, DataAccessError, OSError) as exc:
        _refuse_old_minimap2(exc)
        logger.warning(
            "%s: minimap2 failed against the prebuilt index %s (%s); retrying against the FASTA directly",
            accession,
            reference,
            exc,
        )
        fallback_args = _minimap2_map_args(preset, threads, sam_path, genome_fasta, reads)
        return SecureSubprocess.run_secure("minimap2", fallback_args)


def _refuse_old_minimap2(exc: Exception) -> None:
    """Raise ConfigurationError when minimap2's stderr rejects ``--sam-hit-only`` as an unknown option."""
    stderr = str(getattr(exc, "stderr", None) or "")
    if "sam-hit-only" in stderr and re.search(r"unknown|unrecogni[sz]ed|invalid", stderr, re.IGNORECASE):
        message = f"minimap2 does not support --sam-hit-only; install minimap2 {MINIMAP2_MIN_VERSION} or later"
        raise ConfigurationError(f"{message} ({stderr.strip()})") from exc


def _align_and_count_mapped(
    accession: str,
    preset: str,
    threads: int,
    reference: Path,
    genome_fasta: Path,
    sam_paths: List[Path],
    mate_groups: List[List[Path]],
) -> Tuple[int, bool]:
    """Align each read group to its own SAM and return (total mapped records, saw the
    stderr marker). Mapped records are counted with ``-F 4`` (drops unmapped), before any
    other filtering."""
    mapped_total = 0
    stderr_marker_seen = False
    for sam_path, group in zip(sam_paths, mate_groups):
        aligned = _run_minimap2(accession, preset, threads, sam_path, reference, genome_fasta, group)
        if UNEQUAL_MATES_MARKER in (aligned.stderr or ""):
            stderr_marker_seen = True
        mapped_total += _count_records(sam_path, "-F", "4")
    return mapped_total, stderr_marker_seen
