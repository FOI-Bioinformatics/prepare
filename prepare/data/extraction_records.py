"""Whether a recorded extraction still stands for the current call.

An extraction record names the genome FASTA (and its size and modification time), the minimap2
preset and the mapping quality floor it was made with, the name, size and modification time of
the read files it mapped, and the files it wrote. A rerun reuses it only while all of those still
match; ``read_extraction`` decides what to redo from this.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional

#: Why a registry extraction record was rejected: the call's own parameters (genome, preset,
#: min_mapq) changed since the record was written, versus the parameters still matching but
#: one of the recorded output files being gone from disk.
_PARAMETERS_DIFFER = "recorded parameters differ from this call"
_OUTPUT_FILE_MISSING = "recorded output file missing"
#: The sample's reads differ from the ones the record was made from (e.g. ``download_sra --force``).
_READS_REPLACED = "the reads were replaced since"


def reads_signature(reads: List[Path]) -> List[List[Any]]:
    """``[[name, size bytes, mtime], ...]`` of a sample's FASTQ files, recorded with its extraction."""
    signature: List[List[Any]] = []
    for path in reads:
        st = Path(path).stat()
        signature.append([Path(path).name, st.st_size, st.st_mtime])
    return signature


def _reads_replaced(record: Dict[str, Any], reads: Optional[List[Path]]) -> bool:
    """Whether ``reads`` differ from the read files the record names; False when either is unknown."""
    recorded = record.get("reads_signature")
    if not recorded or not reads:
        return False
    try:
        current = reads_signature(reads)
    except OSError:
        return False
    return [list(entry) for entry in recorded] != current


def genome_signature(genome_path: Path) -> List[float]:
    """``[size, mtime]`` of the genome FASTA, recorded with an extraction to notice a replaced file."""
    st = Path(genome_path).stat()
    return [st.st_size, st.st_mtime]


def _other_genome(record: Dict[str, Any], genome_path: Path) -> bool:
    """Whether the record names another FASTA, or the same path holding another assembly."""
    fasta, signature = record.get("genome_fasta"), record.get("genome_signature")
    if fasta is not None and Path(fasta).resolve() != genome_path.resolve():
        return True
    # A missing FASTA is reported by the caller with its own error; it cannot be compared here.
    return signature is not None and genome_path.is_file() and list(signature) != genome_signature(genome_path)


def _record_mismatch_reason(
    record: Dict[str, Any],
    genome_path: Path,
    preset: str,
    threshold: float,
    min_mapq: int = 0,
    reads: Optional[List[Path]] = None,
) -> Optional[str]:
    """Why a registry extraction record does not match the current call, or None when it does.

    The identity of an extraction is the genome FASTA, the minimap2 preset and the mapping
    quality floor; the containment threshold only decides which samples are tried, so it is
    recorded for information and never compared. A parameter that is absent or None is a
    wildcard. Records rebuilt from disk by ``status --init`` carry no parameters at all, and
    rejecting them would remap every sample of an already extracted project. Distinguishes a
    parameter mismatch from a recorded output file that has since been removed from disk, so
    a caller can log which one triggered the redo instead of always blaming the parameters.

    ``reads`` are the sample's read files as they are now; when given, a record that names other
    read files (another size or modification time: the sample was downloaded again) is stale. A
    record written before read signatures were recorded matches any reads.
    """
    if _other_genome(record, genome_path):
        return _PARAMETERS_DIFFER
    recorded_preset = record.get("preset")
    if recorded_preset is not None and recorded_preset != preset:
        return _PARAMETERS_DIFFER
    recorded_mapq = record.get("min_mapq")
    if recorded_mapq is not None and int(recorded_mapq) != int(min_mapq):
        return _PARAMETERS_DIFFER
    if _reads_replaced(record, reads):
        return _READS_REPLACED
    if record.get("mapped_reads") == 0:
        return None
    if not all(Path(p).exists() for p in record.get("files", [])):
        return _OUTPUT_FILE_MISSING
    return None


def _record_matches(
    record: Dict[str, Any],
    genome_path: Path,
    preset: str,
    threshold: float,
    min_mapq: int = 0,
    reads: Optional[List[Path]] = None,
) -> bool:
    """True when a registry extraction record matches the current call and its files are still usable.

    The key is the genome FASTA, preset and ``min_mapq``; ``threshold`` is accepted for the
    call signature and ignored (see ``_record_mismatch_reason``).
    """
    return _record_mismatch_reason(record, genome_path, preset, threshold, min_mapq=min_mapq, reads=reads) is None
