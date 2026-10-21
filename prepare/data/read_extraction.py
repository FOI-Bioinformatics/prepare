"""Targeted read extraction prior to assembly.

Given a parsed containment table (samples x target genomes), this module maps
each sample's reads against a chosen target genome and keeps only the reads that
align. The resulting reduced FASTQ set supports a small, targeted assembly rather
than a whole-metagenome assembly.

The mapping uses minimap2 (one aligner for both short and long reads) and samtools
to filter and export the mapped reads. External tools run through
``SecureSubprocess.run_secure`` and must be present on the system.
"""

import json
import logging
import platform  # noqa: F401 - tests patch platform.system through this module
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Set, Tuple, Union

import pandas as pd

from prepare.core.exceptions import DataAccessError, ProcessingError
from prepare.data.assembly import (  # noqa: F401 - re-exported; the assembly code moved there
    _megahit_args,
    assemble_extracted_reads,
    fasta_length,
    megahit_version,
    resolve_assembly_threads,
    summarise_contigs,
)
from prepare.data.extraction_locks import LockHeld, index_build_lock, sample_extraction_lock
from prepare.data.extraction_staging import publish_staged, remove_published, require_nonempty, staged_sample_outputs
from prepare.data.file_io import atomic_path, ensure_directory, write_text_atomic
from prepare.data.read_extraction_export import (  # noqa: F401 - re-exported; the export and coverage code moved there
    FILTER_FLAGS,
    MINIMAP2_MIN_VERSION,
    UNEQUAL_MATES_MARKER,
    _align_and_count_mapped,
    _average_read_length,
    _count_records,
    _export_mapped_fastq,
    _fastq_is_empty,
    _filter_and_merge_bam,
    _run_minimap2,
    _minimap2_index_args,
    _minimap2_map_args,
    _refuse_old_minimap2,
    _samtools_cat_args,
    _samtools_count_args,
    _samtools_coverage_args,
    _samtools_fastq_paired_args,
    _samtools_fastq_single_args,
    _samtools_sort_args,
    _samtools_view_args,
    assembly_coverage,
    coverage_table_path,
    reference_coverage,
    summarise_coverage_table,
)
from prepare.data.sra import MATE1_SUFFIXES, fastq_files, fastq_stem, orphan_fastq, primary_fastq
from prepare.utils.security import SecureSubprocess

if TYPE_CHECKING:
    from prepare.data.registry import Registry

logger = logging.getLogger(__name__)

# minimap2 presets by read type; keys are accepted on the CLI via --preset
MINIMAP2_PRESETS = {
    "sr": "sr",  # short paired/single Illumina reads
    "map-ont": "map-ont",  # Oxford Nanopore
    "map-pb": "map-pb",  # PacBio CLR
    "map-hifi": "map-hifi",  # PacBio HiFi
}


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


def selected_samples(
    parsed_containment: Union[str, Path],
    genome_id: str,
    threshold: float,
    registry: Optional["Registry"] = None,
    ignore_selection: bool = False,
) -> List[str]:
    """Return the sample accessions to extract for a genome.

    The table's samples at containment >= ``threshold``, minus the accessions the registry
    excludes, and, when the registry holds a selection and ``ignore_selection`` is off, only
    those the selection names. Without a registry the table alone decides.

    Raises:
        DataAccessError: If the table is missing.
        ProcessingError: If the genome column is absent from the table.
    """
    table_path = Path(parsed_containment)
    if not table_path.exists():
        raise DataAccessError(f"Parsed containment table not found: {table_path}")
    containment = pd.read_csv(table_path, sep="\t", index_col=0)
    samples = select_samples_for_genome(containment, genome_id, threshold)
    if registry is None:
        return samples
    from prepare.data.registry import query  # registry_scan imports this module

    excluded = set(query(registry, "excluded"))
    selection = set() if ignore_selection else set(query(registry, "selected"))
    kept = [acc for acc in samples if acc not in excluded and (not selection or acc in selection)]
    logger.info(
        "%d sample(s) in the table at containment >= %.3f for %s; %d excluded; %s; %d to extract",
        len(samples),
        threshold,
        genome_id,
        sum(acc in excluded for acc in samples),
        f"{len(selection)} in the recorded selection" if selection else "no recorded selection applied",
        len(kept),
    )
    return kept


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
        coverage_table_path(out_dir, genome_id).unlink(missing_ok=True)
        remove_published(out_dir, genome_id)  # a forced rerun that keeps nothing leaves no earlier output
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
        remove_published(out_dir, genome_id)
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


def _samples_to_extract(
    table_path: Path,
    genome_id: str,
    threshold: float,
    samples: Optional[List[str]],
    available: Optional[Set[str]],
    already_done: Dict[str, Dict[str, Any]],
    fastq_root: Path,
) -> List[str]:
    """The samples ``extract_target_reads`` loops over: the given list or the table's, minus those without FASTQ."""
    if samples is None:
        containment = pd.read_csv(table_path, sep="\t", index_col=0)
        samples = select_samples_for_genome(containment, genome_id, threshold)
        logger.info("%d sample(s) meet containment >= %.3f for %s", len(samples), threshold, genome_id)
    if available is None:
        return samples
    # A sample already recorded (``already_done``) is kept even when it is not in
    # ``available``: its FASTQ input may have been removed since a successful
    # extraction (e.g. by store_gc/store_unlink), and it must still reach the
    # "already extracted" fast path rather than being counted as missing.
    known = set(available) | set(already_done.keys())
    missing = [acc for acc in samples if acc not in known]
    if missing:
        logger.info(
            "%d of the %d selected sample(s) have no FASTQ under %s; skipped", len(missing), len(samples), fastq_root
        )
    return [acc for acc in samples if acc in known]


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
    samples: Optional[List[str]] = None,
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
        samples: The accessions to extract, as ``selected_samples`` returns them (the table
            bounded by the registry's exclusions and selection). ``None`` takes every sample
            of the table at ``threshold``.

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

    already_done = already_done or {}
    samples = _samples_to_extract(table_path, genome_id, threshold, samples, available, already_done, fastq_root)

    results: Dict[str, ExtractionResult] = {}
    output_root = Path(output_folder)
    SecureSubprocess.add_allowed_root(output_root)
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
