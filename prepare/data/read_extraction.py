"""Targeted read extraction prior to assembly.

Given a parsed containment table (samples x target genomes), this module maps
each sample's reads against a chosen target genome and keeps only the reads that
align. The resulting reduced FASTQ set supports a small, targeted assembly rather
than a whole-metagenome assembly.

The mapping uses minimap2 (one aligner for both short and long reads) and samtools
to filter and export the mapped reads. External tools run through
``SecureSubprocess.run_secure`` and must be present on the system.
"""

import gzip
import json
import logging
import platform  # noqa: F401 - tests patch platform.system through this module
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple, Union

import pandas as pd

from prepare.core.exceptions import ConfigurationError, DataAccessError, ProcessingError, SecurityError
from prepare.data.assembly import (  # noqa: F401 - re-exported; the assembly code moved there
    _megahit_args,
    assemble_extracted_reads,
    fasta_length,
    megahit_version,
    resolve_assembly_threads,
    summarise_contigs,
)
from prepare.data.extraction_locks import LockHeld, index_build_lock, sample_extraction_lock
from prepare.data.extraction_staging import publish_staged, require_nonempty, staged_sample_outputs
from prepare.data.file_io import atomic_path, ensure_directory, write_text_atomic
from prepare.data.sra import MATE1_SUFFIXES, fastq_files, fastq_stem, orphan_fastq, primary_fastq
from prepare.utils.security import SecureSubprocess

logger = logging.getLogger(__name__)

# minimap2 presets by read type; keys are accepted on the CLI via --preset
MINIMAP2_PRESETS = {
    "sr": "sr",  # short paired/single Illumina reads
    "map-ont": "map-ont",  # Oxford Nanopore
    "map-pb": "map-pb",  # PacBio CLR
    "map-hifi": "map-hifi",  # PacBio HiFi
}

# minimap2 prints this when the two mate files differ in length; it then maps single-end.
UNEQUAL_MATES_MARKER = "different number of records"

# samtools -F flags dropped by the alignment filter: unmapped (0x4), secondary (0x100)
# and supplementary (0x800) alignments.
FILTER_FLAGS = "0x904"

MINIMAP2_MIN_VERSION = "2.17"  # the first minimap2 release with --sam-hit-only (minimap2 NEWS.md)


@dataclass
class ExtractionResult:
    """Outcome of mapping one sample against the target genome."""

    files: List[Path]
    mapped_records: int
    unequal_mates: bool = False
    skipped: bool = False
    # Mapped records before the secondary/supplementary/MAPQ filter (0 for a skipped result).
    mapped_total: int = 0
    # Coverage of the reference genome by the kept alignments (``breadth``, ``mean_depth``,
    # ``covered_bases``, ``reference_bp``, ``coverage_tsv``); None when nothing mapped or
    # ``samtools coverage`` could not be run.
    coverage: Optional[Dict[str, Any]] = None


def _notify_result(
    on_result: Optional[Callable[[str, "ExtractionResult"], None]], accession: str, result: "ExtractionResult"
) -> None:
    """Call ``on_result`` in isolation so a failing callback never stops the run.

    The callback records the result (a registry write, which can raise, for example on a
    lock timeout). One sample's bookkeeping failing must not cost the remaining samples
    their extraction.
    """
    if on_result is None:
        return
    try:
        on_result(accession, result)
    except Exception as e:
        logger.warning("Recording the extraction result for %s failed: %s", accession, e)


def _count_records(path: Path, *filter_args: str) -> int:
    """Number of records in a SAM/BAM file, via ``samtools view -c``.

    ``filter_args`` are extra ``samtools view`` flags (e.g. ``"-F", "4"``) applied before
    the count; with none, every record in the file is counted.
    """
    result = SecureSubprocess.run_secure("samtools", _samtools_count_args(filter_args, path))
    text = (result.stdout or "").strip()
    return int(text) if text.isdigit() else 0


def resolve_index_path(genome_fasta: Union[str, Path], preset: str, index_dir: Union[str, Path]) -> Path:
    """The minimap2 index path ``build_index`` reads or writes for this genome and preset."""
    return Path(index_dir) / f"{Path(genome_fasta).stem}.{preset}.mmi"


def _index_source(genome_path: Path) -> Dict[str, Any]:
    """The identity of the FASTA an index was built from: resolved path, size and mtime."""
    stat = genome_path.stat()
    return {"fasta": str(genome_path.resolve()), "bytes": stat.st_size, "mtime": stat.st_mtime}


def _index_is_current(index_path: Path, source: Dict[str, Any]) -> bool:
    """True when ``index_path`` is non-empty and its record says it was built from ``source``.

    The record is compared in full because the index is named after the FASTA's stem only:
    a different genome with the same file name (a second assembly called ``wMel.fna``, a copy
    restored from a tarball) would otherwise be mapped against the wrong index; ``cp -p``,
    ``mv``, ``rsync -a`` and tar all preserve an older mtime, so an mtime alone does not tell.
    """
    if not index_path.exists() or index_path.stat().st_size == 0:
        return False
    record_path = index_path.with_suffix(index_path.suffix + ".json")
    try:
        recorded = json.loads(record_path.read_text())
    except (OSError, ValueError):
        return False
    return isinstance(recorded, dict) and all(recorded.get(key) == value for key, value in source.items())


def build_index(genome_fasta: Union[str, Path], preset: str, index_dir: Union[str, Path]) -> Path:
    """Build (or reuse) a minimap2 index for the target genome, shared across every sample.

    The index lives at ``resolve_index_path(genome_fasta, preset, index_dir)``, with a
    ``<index>.json`` record beside it naming the FASTA it was built from (resolved path,
    size and mtime). It is reused only when all three still match, and rebuilt otherwise.

    The build writes to a unique temporary name and is moved into place, so a second
    run against the same output folder either sees the previous index or none at all, never
    a half-written one; it also runs under a lock, so two processes never build it at once.
    """
    genome_path = Path(genome_fasta)
    index_root = ensure_directory(index_dir)
    index_path = resolve_index_path(genome_path, preset, index_root)
    source = _index_source(genome_path)
    if _index_is_current(index_path, source):
        return index_path

    with index_build_lock(index_path):
        if _index_is_current(index_path, source):
            return index_path
        with atomic_path(index_path) as staged:
            SecureSubprocess.run_secure("minimap2", _minimap2_index_args(preset, staged, genome_path))
            require_nonempty(staged, "The minimap2 index")
        write_text_atomic(index_path.with_suffix(index_path.suffix + ".json"), json.dumps(source, indent=2))
        return index_path


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


def _fastq_is_empty(path: Path) -> bool:
    """True when the file is missing or has no records (gzip or plain)."""
    if not path.exists():
        return True
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return handle.readline() == ""


def select_samples_for_genome(containment: pd.DataFrame, genome_id: str, threshold: float) -> List[str]:
    """Return the sample accessions whose containment for a genome meets a threshold.

    Args:
        containment: Parsed containment table indexed by sample accession.
        genome_id: Target genome column to select on.
        threshold: Minimum containment value (inclusive).

    Raises:
        ProcessingError: If the genome column is absent from the table.
    """
    if genome_id not in containment.columns:
        raise ProcessingError(
            f"Genome '{genome_id}' not found in the containment table. "
            f"Available columns include: {', '.join(map(str, containment.columns[:5]))}"
        )
    values = pd.to_numeric(containment[genome_id], errors="coerce")
    return [str(acc) for acc in containment.index[values >= threshold]]


def selected_samples(parsed_containment: Union[str, Path], genome_id: str, threshold: float) -> List[str]:
    """Return the sample accessions that meet containment >= threshold for a genome.

    Reads the parsed containment table and delegates to ``select_samples_for_genome``.
    Exposed so a caller (the CLI) can learn which samples were selected without
    duplicating that logic, or re-running the full extraction.

    Raises:
        DataAccessError: If the table is missing.
        ProcessingError: If the genome column is absent from the table.
    """
    table_path = Path(parsed_containment)
    if not table_path.exists():
        raise DataAccessError(f"Parsed containment table not found: {table_path}")
    containment = pd.read_csv(table_path, sep="\t", index_col=0)
    return select_samples_for_genome(containment, genome_id, threshold)


#: Why a registry extraction record was rejected: the call's own parameters (genome, preset,
#: threshold, min_mapq) changed since the record was written, versus the parameters still
#: matching but one of the recorded output files being gone from disk.
_PARAMETERS_DIFFER = "recorded parameters differ from this call"
_OUTPUT_FILE_MISSING = "recorded output file missing"


def _record_mismatch_reason(
    record: Dict[str, Any], genome_path: Path, preset: str, threshold: float, min_mapq: int = 0
) -> Optional[str]:
    """Why a registry extraction record does not match the current call, or None when it does.

    A parameter that is absent or None is a wildcard. Records rebuilt from disk by
    ``status --init`` carry no parameters at all, and rejecting them would remap every
    sample of an already extracted project. Distinguishes a parameter mismatch from a
    recorded output file that has since been removed from disk, so a caller can log which
    one triggered the redo instead of always blaming the parameters.
    """
    recorded_fasta = record.get("genome_fasta")
    if recorded_fasta is not None and Path(recorded_fasta).resolve() != genome_path.resolve():
        return _PARAMETERS_DIFFER
    recorded_preset = record.get("preset")
    if recorded_preset is not None and recorded_preset != preset:
        return _PARAMETERS_DIFFER
    recorded_threshold = record.get("threshold")
    if recorded_threshold is not None and float(recorded_threshold) != float(threshold):
        return _PARAMETERS_DIFFER
    recorded_mapq = record.get("min_mapq")
    if recorded_mapq is not None and int(recorded_mapq) != int(min_mapq):
        return _PARAMETERS_DIFFER
    if record.get("mapped_reads") == 0:
        return None
    if not all(Path(p).exists() for p in record.get("files", [])):
        return _OUTPUT_FILE_MISSING
    return None


def _record_matches(
    record: Dict[str, Any], genome_path: Path, preset: str, threshold: float, min_mapq: int = 0
) -> bool:
    """True when a registry extraction record matches the current call and its files are still usable.

    A parameter that is absent or None is a wildcard. Records rebuilt from disk by
    ``status --init`` carry no parameters at all, and rejecting them would remap every
    sample of an already extracted project.
    """
    return _record_mismatch_reason(record, genome_path, preset, threshold, min_mapq=min_mapq) is None


def _skipped_result(record: Dict[str, Any]) -> ExtractionResult:
    """The result of a sample whose recorded extraction still stands."""
    return ExtractionResult(
        files=[Path(p) for p in record.get("files", [])],
        mapped_records=int(record.get("mapped_reads") or 0),
        unequal_mates=bool(record.get("unequal_mates", False)),
        skipped=True,
    )


def _sample_reads(fastq_folder: Path, accession: str) -> List[Path]:
    """Return the FASTQ files for one accession that minimap2 should map, R1 before R2.

    For paired data that is the mate pair alone: the bare ``<acc>.fastq`` file that
    ``--split-3`` writes for unpaired spots would make minimap2 read three files as an
    interleaved set and report mismatched mate counts, so it is left out. Single-end data
    gives the bare file. Zero-byte files and hidden temporary leftovers of an interrupted
    download are excluded by ``fastq_files``.
    """
    acc_dir = fastq_folder / accession
    files = fastq_files(acc_dir)
    if not files:
        return []

    primary = primary_fastq(acc_dir)
    if primary is None:
        return files

    stem = fastq_stem(primary)
    for marker in MATE1_SUFFIXES:
        if stem.endswith(marker):
            mate_stem = stem[: -len(marker)] + marker[:-1] + "2"
            mate_two = next((p for p in files if fastq_stem(p) == mate_stem), None)
            if mate_two is None:
                return files
            orphan = orphan_fastq(acc_dir)
            if orphan is not None:
                logger.debug("%s: mapping the mate files only; %s holds unpaired reads", accession, orphan.name)
            return [primary, mate_two]

    return [primary]


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
    """Export a filtered BAM back to FASTQ. Reads without a mate flag (single-end mapping,
    or the unequal-mates fallback) go to the ``-0`` file."""
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
        return [out1, out2]
    if orphans.exists():
        return [orphans]
    return [singles] if singles.exists() else []


def _map_and_extract(
    accession: str,
    reads: List[Path],
    reference: Path,
    genome_fasta: Path,
    out_dir: Path,
    genome_id: str,
    preset: str,
    threads: int,
    min_mapq: int = 0,
    sam_dir: Optional[Path] = None,
    keep_sam: bool = False,
    force_single_end: bool = False,
) -> ExtractionResult:
    """Map one sample's reads to the target genome and write the mapped reads.

    ``reference`` is the prebuilt minimap2 index shared across every sample of this call;
    ``genome_fasta`` is the FASTA it was built from, used as a one-time fallback if aligning
    against the index fails. ``force_single_end`` (set when the caller already knows the two
    mate files disagree in read count) maps each mate file in its own minimap2 run and merges
    the two alignments as single-end, rather than relying on minimap2's stderr warning.

    Returns the FASTQ files written (two for paired input, one otherwise, none when nothing
    mapped) with the mapped-record counts and the reference coverage (its table is written
    beside the FASTQ, which is staged: ``extraction_staging``). SAM files (under ``sam_dir`` or
    ``out_dir``) go once the BAM exists unless ``keep_sam``; the BAM goes after the export."""
    ensure_directory(out_dir)
    sam_root = ensure_directory(sam_dir) if sam_dir is not None else out_dir
    bam_path = out_dir / f"{genome_id}.mapped.bam"
    filter_args = ["-F", FILTER_FLAGS]
    if min_mapq:
        filter_args += ["-q", str(min_mapq)]

    unequal = force_single_end and len(reads) >= 2
    if unequal:
        mate_groups: List[List[Path]] = [[read] for read in reads[:2]]
        sam_paths = [sam_root / f"{genome_id}.mate1.sam", sam_root / f"{genome_id}.mate2.sam"]
        logger.warning(
            "%s: the mate files are recorded with different read counts, so each is mapped "
            "independently as single-end reads; re-download with fasterq-dump (which keeps mates "
            "in step) for paired extraction",
            accession,
        )
    else:
        mate_groups = [reads]
        sam_paths = [sam_root / f"{genome_id}.sam"]

    mapped_total, stderr_marker_seen = _align_and_count_mapped(
        accession, preset, threads, reference, genome_fasta, sam_paths, mate_groups
    )
    if stderr_marker_seen and not unequal:
        unequal = True
        logger.warning(
            "%s: the mate files have different read counts, so minimap2 mapped them as single-end reads; "
            "re-download with fasterq-dump (which keeps mates in step) for paired extraction",
            accession,
        )

    if mapped_total == 0:
        logger.warning("No reads from %s mapped to %s; nothing written", accession, genome_id)
        if not keep_sam:
            for sam_path in sam_paths:
                sam_path.unlink(missing_ok=True)
        # A forced rerun that now keeps nothing must not leave an earlier run's table behind.
        coverage_table_path(out_dir, genome_id).unlink(missing_ok=True)
        return ExtractionResult([], 0, unequal, mapped_total=0)

    _filter_and_merge_bam(sam_paths, filter_args, threads, out_dir, genome_id, bam_path)
    if not keep_sam:
        for sam_path in sam_paths:
            sam_path.unlink(missing_ok=True)

    mapped = _count_records(bam_path)
    if mapped == 0:
        logger.warning("No reads from %s mapped to %s; nothing written", accession, genome_id)
        bam_path.unlink(missing_ok=True)
        coverage_table_path(out_dir, genome_id).unlink(missing_ok=True)
        return ExtractionResult([], 0, unequal, mapped_total=mapped_total)
    mapq_clause = f" and MAPQ below {min_mapq}" if min_mapq > 0 else ""
    logger.info(
        "%s: kept %d of %d mapped records (secondary/supplementary%s removed: %d)",
        accession,
        mapped,
        mapped_total,
        mapq_clause,
        mapped_total - mapped,
    )

    with staged_sample_outputs(out_dir, genome_id) as staging:
        staged = _export_mapped_fastq(reads, bam_path, staging, genome_id, threads)
        coverage = reference_coverage(accession, genome_id, bam_path, out_dir, sam_root, threads)
        written = publish_staged(staging, staged, out_dir, genome_id)
    bam_path.unlink(missing_ok=True)

    return ExtractionResult(written, mapped, unequal, mapped_total=mapped_total, coverage=coverage)


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


def _resolve_done_state(
    accession: str,
    record: Optional[Dict[str, Any]],
    genome_path: Path,
    preset: str,
    threshold: float,
    force: bool,
    min_mapq: int = 0,
) -> bool:
    """True when a recorded extraction still stands and can be skipped.

    Logs the reason for a redo (parameters differ vs. a recorded output file missing) when
    a record exists, is not forced, and does not match. Split out of ``_extract_one_sample``
    to keep that function's branching simple.
    """
    if not record or force:
        return False
    mismatch_reason = _record_mismatch_reason(record, genome_path, preset, threshold, min_mapq=min_mapq)
    if mismatch_reason is not None:
        logger.info("%s: redoing extraction, %s", accession, mismatch_reason)
        return False
    return True


def _skip_if_truncated(accession: str, verdict: Optional[Dict[str, Any]], allow_truncated: bool) -> bool:
    """True (after logging) when this sample's truncated download or unusable store copy should be skipped."""
    if not verdict or allow_truncated:
        return False
    spots = f"download truncated ({verdict.get('reads_r1', 0)} of {verdict.get('expected_spots', 0)} spots)"
    logger.warning("skipped %s: %s; use --allow-truncated", accession, verdict.get("reason", spots))
    return True


def _extract_one_sample(
    accession: str,
    record: Optional[Dict[str, Any]],
    genome_path: Path,
    preset: str,
    threshold: float,
    force: bool,
    dry_run: bool,
    fastq_root: Path,
    genome_id: str,
    output_root: Path,
    threads: int,
    min_mapq: int,
    sam_dir: Optional[Path],
    keep_sam: bool,
    mate_counts: Dict[str, Tuple[int, int]],
    reference_holder: Dict[str, Path],
) -> Optional[ExtractionResult]:
    """Handle one selected sample: an already-done skip, missing reads, dry-run reporting,
    or a real mapping run. Returns None when the sample contributes nothing to ``results``
    (no FASTQ files found for it). The shared minimap2 index is built on first use here and
    cached in ``reference_holder["reference"]`` for every later sample of this call. The
    mapping run itself holds this sample's own lock; one already held elsewhere (an
    overlapping SLURM shard extracting the same sample) is reported skipped, not run twice.
    """
    done = _resolve_done_state(accession, record, genome_path, preset, threshold, force, min_mapq=min_mapq)
    if done and not dry_run:
        result = _skipped_result(record or {})
        logger.info(
            "%s already extracted against %s (%d mapped reads); use --force to redo",
            accession,
            genome_id,
            result.mapped_records,
        )
        return result

    reads = _sample_reads(fastq_root, accession)
    if not reads:
        logger.warning("No FASTQ files found for %s under %s; skipping", accession, fastq_root)
        return None

    if dry_run:
        result = _skipped_result(record or {}) if done else ExtractionResult([], 0)
        if done:
            logger.info(
                "would skip %s (already extracted, %d mapped reads); use --force to redo",
                accession,
                result.mapped_records,
            )
        return result

    if "reference" not in reference_holder:
        index_dir = output_root / ".index"
        SecureSubprocess.add_allowed_root(index_dir)
        reference_holder["reference"] = build_index(genome_path, preset, index_dir)

    counts = mate_counts.get(accession)
    force_single_end = bool(counts is not None and len(reads) >= 2 and counts[0] != counts[1])
    try:
        with sample_extraction_lock(output_root, accession, genome_id):
            outcome = _map_and_extract(
                accession,
                reads,
                reference_holder["reference"],
                genome_path,
                output_root / accession,
                genome_id,
                preset,
                threads,
                min_mapq=min_mapq,
                sam_dir=sam_dir,
                keep_sam=keep_sam,
                force_single_end=force_single_end,
            )
    except LockHeld:
        logger.info("%s: extraction against %s is in progress elsewhere; skipped", accession, genome_id)
        return ExtractionResult(files=[], mapped_records=0, skipped=True)
    if outcome.files:
        logger.info(
            "Extracted %d mapped records for %s -> %s",
            outcome.mapped_records,
            accession,
            ", ".join(str(p) for p in outcome.files),
        )
    return outcome


def extract_target_reads(
    parsed_containment: Union[str, Path],
    genome_id: str,
    genome_fasta: Union[str, Path],
    fastq_folder: Union[str, Path] = "fastq",
    output_folder: Union[str, Path] = "targeted",
    threshold: float = 0.1,
    preset: str = "sr",
    threads: int = 4,
    dry_run: bool = False,
    force: bool = False,
    already_done: Optional[Dict[str, Dict[str, Any]]] = None,
    on_result: Optional[Callable[[str, ExtractionResult], None]] = None,
    min_mapq: int = 0,
    temp_folder: Optional[Union[str, Path]] = None,
    allow_truncated: bool = False,
    mate_counts: Optional[Dict[str, Tuple[int, int]]] = None,
    truncated_downloads: Optional[Dict[str, Dict[str, Any]]] = None,
    keep_sam: bool = False,
    available: Optional[Set[str]] = None,
) -> Dict[str, ExtractionResult]:
    """Extract reads mapping to a target genome for every qualifying sample.

    A sample already locked by another process or thread (an overlapping SLURM shard that
    selected it independently) is reported skipped rather than mapped twice; dry runs, which
    touch no file, take no lock.

    Args:
        parsed_containment: Parsed containment table (samples x genomes).
        genome_id: Target genome to extract against (a column in the table).
        genome_fasta: FASTA file for the target genome.
        fastq_folder: Root folder with per-accession FASTQ subdirectories.
        output_folder: Root folder for the extracted reads (per-accession subdirs).
        threshold: Minimum containment for a sample to be included.
        preset: minimap2 preset (sr, map-ont, map-pb, map-hifi).
        threads: Threads for minimap2 and samtools.
        dry_run: If True, report the qualifying samples without running any tool.
        force: If True, redo extraction even for a sample already recorded in ``already_done``.
        already_done: Accession -> the registry's extraction record for this genome; skipped
            (unless ``force``) when its genome FASTA, preset and threshold still match this
            call and its files are still on disk (or it mapped 0).
        on_result: Called on the main thread with (accession, result) after each real-run
            sample, skipped ones included, so the caller can checkpoint as it lands. A raising
            callback is logged, not fatal. Dry runs never call it.
        min_mapq: Minimum mapping quality (samtools ``-q``) kept, beyond dropping
            secondary/supplementary alignments. 0 keeps every mapped record, since a divergent
            strain of the target genome can map with a genuinely low MAPQ.
        temp_folder: Where the intermediate SAM alignment(s) are written; defaults to the
            sample's own output folder.
        allow_truncated: If True, extract even a sample listed in ``truncated_downloads``;
            otherwise that sample is skipped.
        mate_counts: Accession -> (mate 1 reads, mate 2 reads). A mismatch maps the two mate
            files independently as single-end, rather than relying on minimap2's own warning.
        truncated_downloads: Accession -> download verdict, or a ``reason`` for an unusable store
            copy, for every accession to skip unless ``allow_truncated``.
        keep_sam: If True, keep the intermediate SAM alignment(s) for debugging.
        available: Accessions with FASTQ files on disk. A selected sample not in it is left
            out and counted rather than logged individually, so a dry run over many thousands
            of screened-but-not-downloaded samples prints one summary line. ``None`` checks
            each selected sample's FASTQ files individually instead.

    Returns:
        Mapping of accession to an ``ExtractionResult`` (empty files in dry-run).

    Raises:
        DataAccessError: If inputs are missing.
        ProcessingError: If the genome column or preset is invalid.
    """
    table_path = Path(parsed_containment)
    if not table_path.exists():
        raise DataAccessError(f"Parsed containment table not found: {table_path}")
    if preset not in MINIMAP2_PRESETS:
        raise ProcessingError(f"Unknown minimap2 preset '{preset}'. Choose one of: {', '.join(MINIMAP2_PRESETS)}")

    fastq_root = Path(fastq_folder)
    genome_path = Path(genome_fasta)
    if not dry_run and not genome_path.exists():
        raise DataAccessError(f"Target genome FASTA not found: {genome_path}")

    containment = pd.read_csv(table_path, sep="\t", index_col=0)
    samples = select_samples_for_genome(containment, genome_id, threshold)
    logger.info("%d sample(s) meet containment >= %.3f for %s", len(samples), threshold, genome_id)
    if available is not None:
        # A sample already recorded (``already_done``) is kept even when it is not in
        # ``available``: its FASTQ input may have been removed since a successful
        # extraction (e.g. by store_gc/store_unlink), and it must still reach the
        # "already extracted" fast path below rather than being counted as missing.
        known = set(available) | set((already_done or {}).keys())
        missing = [acc for acc in samples if acc not in known]
        samples = [acc for acc in samples if acc in known]
        if missing:
            logger.info(
                "%d of the %d selected sample(s) have no FASTQ under %s; skipped",
                len(missing),
                len(missing) + len(samples),
                fastq_root,
            )

    results: Dict[str, ExtractionResult] = {}
    output_root = Path(output_folder)
    SecureSubprocess.add_allowed_root(output_root)
    already_done = already_done or {}
    mate_counts = mate_counts or {}
    truncated_downloads = truncated_downloads or {}

    # The index is built lazily, the first time a sample actually needs mapping, so a run
    # where every sample is already done (or dry) never touches minimap2 at all. Held in a
    # dict (rather than a local reassigned across the loop) so the helper below can cache
    # it too.
    reference_holder: Dict[str, Path] = {}
    sam_dir: Optional[Path] = None
    if temp_folder is not None:
        sam_dir = Path(temp_folder)
        SecureSubprocess.add_allowed_root(sam_dir)

    for accession in samples:
        if _skip_if_truncated(accession, truncated_downloads.get(accession), allow_truncated):
            continue
        result = _extract_one_sample(
            accession,
            already_done.get(accession),
            genome_path,
            preset,
            threshold,
            force,
            dry_run,
            fastq_root,
            genome_id,
            output_root,
            threads,
            min_mapq,
            sam_dir,
            keep_sam,
            mate_counts,
            reference_holder,
        )
        if result is None:
            continue
        results[accession] = result
        if not dry_run:
            _notify_result(on_result, accession, result)

    return results


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
    sam_path = work_root / "coverage.sam"
    bam_path = work_root / "coverage.bam"
    try:
        SecureSubprocess.run_secure("minimap2", _minimap2_map_args(preset, threads, sam_path, contigs, reads))
        SecureSubprocess.run_secure("samtools", _samtools_view_args(["-F", FILTER_FLAGS], threads, bam_path, sam_path))
        reads_mapped = _count_records(bam_path)
    finally:
        sam_path.unlink(missing_ok=True)
        bam_path.unlink(missing_ok=True)

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
