"""
Sourmash sketches of downloaded datasets, for ``sra_sketch``.

``sketch_dataset`` reads every record of an accession's FASTQ files through
``prepare.data.sra.sampling.iter_sequence_blocks`` and feeds each block to one FracMinHash
sketch per k-mer size in ``KSIZES``, all at ``SCALED`` with abundance tracked. The k = 21
sketch is built exactly as ``prepare.data.branchwater_search.sketch_fasta`` builds the
Branchwater query, so the two can be compared. ``write_signature`` stores the sketches as one
gzip-compressed JSON signature file, the format ``sourmash sketch -o x.sig.gz`` writes, with
the k = 21 sketch first so ``load_signature`` finds it. sourmash is imported at the point of use
(the ``sourmash`` extra).
"""

import logging
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

from prepare.core.exceptions import DataAccessError
from prepare.core.optional import require
from prepare.data.branchwater_search import KSIZE, SCALED
from prepare.data.file_io import open_atomic
from prepare.data.sra.sampling import iter_sequence_blocks
from prepare.store.stats import _file_signature
from prepare.utils.progress import active_item_level

logger = logging.getLogger(__name__)

# The k-mer sizes of every dataset sketch, ascending; Branchwater's k must be among them.
KSIZES = (21, 31, 51)
assert KSIZE in KSIZES

SIGNATURE_SUFFIX = ".sig.gz"
_PURPOSE = "Sketching reads"
# Blocks (of CHUNK_SIZE bytes) between progress lines while one file is sketched.
PROGRESS_BLOCKS = 200


@dataclass
class DatasetSketch:
    """One accession's sketches and the figures gathered while building them."""

    accession: str
    files: List[Path]
    reads: int
    bases: int
    hashes: Dict[int, int]
    signature: Dict[str, List[Union[int, float]]]
    minhashes: Dict[int, Any] = field(default_factory=dict, repr=False)


def signature_path(accession: str, output_dir: Union[str, Path]) -> Path:
    """``<output_dir>/<accession>.sig.gz``."""
    return Path(output_dir) / f"{accession}{SIGNATURE_SUFFIX}"


def sketch_dataset(
    accession: str, files: Sequence[Union[str, Path]], stop: Optional[threading.Event] = None
) -> DatasetSketch:
    """Sketch every read of ``files`` (an accession's mate and orphan files) at each k in ``KSIZES``.

    Each block of joined sequences goes into one ``add_sequence`` call per sketch; sourmash
    releases the GIL in that call, so several accessions can be sketched on threads at once.
    ``stop`` ends the read with ``InterruptedError`` between blocks. Read errors propagate as
    ``prepare.sra.dataset_stats.DATASET_READ_ERRORS``.
    """
    sourmash = require("sourmash", "sourmash", _PURPOSE)
    paths = [Path(f) for f in files]
    minhashes = {k: sourmash.MinHash(n=0, ksize=k, scaled=SCALED, track_abundance=True) for k in KSIZES}
    reads = bases = 0
    for path in paths:
        file_reads = 0
        for index, (block, n) in enumerate(iter_sequence_blocks(path, stop), 1):
            file_reads += n
            bases += len(block) - (n - 1)
            for minhash in minhashes.values():
                minhash.add_sequence(block, force=True)
            if index % PROGRESS_BLOCKS == 0:
                logger.log(active_item_level(), "Sketching %s: %d reads so far", path.name, file_reads)
        reads += file_reads
        logger.log(active_item_level(), "Sketched %s: %d reads", path.name, file_reads)
    return DatasetSketch(
        accession=accession,
        files=paths,
        reads=reads,
        bases=bases,
        hashes={k: len(mh) for k, mh in minhashes.items()},
        signature=_file_signature(paths),
        minhashes=minhashes,
    )


def write_signature(sketch: DatasetSketch, path: Union[str, Path]) -> Path:
    """Write ``sketch`` as a gzip-compressed sourmash signature file at ``path``, atomically.

    One signature per k, in ascending k, each named by the accession. The file appears under
    its final name only once it is complete, so a failure or an interrupt leaves nothing behind.
    """
    sourmash = require("sourmash", "sourmash", _PURPOSE)
    signature_module = require("sourmash.signature", "sourmash", _PURPOSE)
    signatures = [sourmash.SourmashSignature(sketch.minhashes[k], name=sketch.accession) for k in KSIZES]
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with open_atomic(target, "wb") as handle:
        signature_module.save_signatures_to_json(signatures, handle, compression=1)
    return target


def sketch_summary(sketch: DatasetSketch, source: str = "computed") -> Dict[str, Any]:
    """The registry summary of ``sketch``: parameters, totals, hashes per k and the file snapshot."""
    return {
        "ksizes": list(KSIZES),
        "scaled": SCALED,
        "track_abundance": True,
        "reads": sketch.reads,
        "bases": sketch.bases,
        "hashes": {str(k): n for k, n in sketch.hashes.items()},
        "files": sketch.signature,
        "source": source,
    }


def sketch_is_current(summary: Optional[Dict[str, Any]], files: Sequence[Path], path: Union[str, Path]) -> bool:
    """Whether a recorded sketch (``summary``, its file at ``path``) still describes ``files``."""
    if not summary or not Path(path).is_file():
        return False
    return (
        summary.get("files") == _file_signature([Path(f) for f in files])
        and summary.get("ksizes") == list(KSIZES)
        and summary.get("scaled") == SCALED
    )


def load_dataset_minhash(path: Union[str, Path], ksize: int = KSIZE) -> Any:
    """The flattened (abundance-free) sketch at ``ksize`` from a signature file, checked for ``SCALED``."""
    sourmash = require("sourmash", "sourmash", _PURPOSE)
    try:
        loaded = next(iter(sourmash.load_file_as_signatures(str(path), ksize=ksize)))
    except (OSError, ValueError, StopIteration) as e:
        raise DataAccessError(f"No k={ksize} sketch in {path}: {e}") from e
    minhash = loaded.minhash
    if minhash.scaled != SCALED:
        raise DataAccessError(f"{path} has scaled={minhash.scaled}; the dataset sketches use scaled={SCALED}")
    return minhash.flatten()
