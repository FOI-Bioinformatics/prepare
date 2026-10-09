"""The FASTQ files of one accession and its shared statistics record.

``sra_profile`` and ``sra_report`` both profile a dataset from the same two inputs: every mate
file of the accession, and the statistics record ``prepare.store.stats.compute_dataset_stats``
computes over them (exact read and base totals, GC). The record is cached in the store sidecar
when the accession is a link into the shared store, so a dataset is read in full once.
"""

import logging
import zlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

from prepare.core.exceptions import DataAccessError, PrepareError
from prepare.data.file_io import visible_files
from prepare.data.sra import fastq_files
from prepare.data.sra_metadata import _resolved_sidecar_path
from prepare.store.stats import DEFAULT_SAMPLE_SIZE, cached_stats, compute_dataset_stats, store_stats

logger = logging.getLogger(__name__)

# What reading a dataset's FASTQ files raises for a file that is missing, unreadable,
# truncated or malformed (a corrupt gzip stream raises zlib.error, which is not an OSError);
# PrepareError covers the store and analyzer layers.
DATASET_READ_ERRORS = (OSError, EOFError, ValueError, UnicodeDecodeError, zlib.error, PrepareError)


def accession_fastq_files(root: Union[str, Path], accession: str) -> List[Path]:
    """Every non-empty FASTQ file of ``accession`` under ``root``, mates in name order.

    Downloads live in ``<root>/<accession>/`` (fasterq-dump layout, possibly a link into the
    shared store); flat ``<root>/<accession>.fastq.gz`` and ``<root>/<accession>_1.fastq``
    files are accepted too. The flat patterns require an exact accession match or an
    underscore right after it, so ``SRR1`` never picks up ``SRR10``'s files.
    """
    base = Path(root)
    acc_dir = base / accession
    if acc_dir.is_dir():
        files = fastq_files(acc_dir)
        if files:
            return files
    flat = visible_files(base, f"{accession}.fastq*", f"{accession}.fq*", f"{accession}_*.fastq*", f"{accession}_*.fq*")
    return [f for f in flat if f.stat().st_size > 0]


def load_dataset_stats(
    files: Sequence[Union[str, Path]], sample_size: int = DEFAULT_SAMPLE_SIZE
) -> Optional[Dict[str, Any]]:
    """The shared statistics record for one dataset's ``files``, or None when it cannot be read.

    Reuses the store sidecar's cached record when its signature still matches the files on
    disk, otherwise computes it here. The record is written back to the sidecar only when
    there is one (a plain project folder keeps it for this run only); a write that fails is a
    lost cache, not a failed profile, so it is reported as a warning and the record is used.
    """
    paths = [Path(f) for f in files]
    if not paths:
        return None
    acc_dir = paths[0].parent
    sidecar_path = _resolved_sidecar_path(acc_dir)
    try:
        cached = cached_stats(acc_dir, sidecar_path)
        if cached is not None:
            return cached
        stats = compute_dataset_stats(list(paths), sample_size=sample_size)
    except DATASET_READ_ERRORS as e:
        logger.debug("Could not compute dataset stats for %s: %s", acc_dir.name, e)
        return None
    if sidecar_path is not None:
        try:
            store_stats(sidecar_path, stats)
        except (OSError, DataAccessError) as e:
            logger.warning("Could not cache statistics for %s: %s", acc_dir.name, e)
    return stats
