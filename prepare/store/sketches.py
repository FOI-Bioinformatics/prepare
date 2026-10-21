"""
One cached sourmash signature per store dataset.

``sra_sketch`` builds a signature from a dataset's reads; when the dataset is a copy in the
shared store, ``store_sketch`` keeps the signature beside it (``<store>/sra/<ACC>/<ACC>.sig.gz``)
and records its description in the sidecar's ``sketch`` block, so another project linking the
same dataset reuses it through ``cached_sketch`` instead of reading the reads again. The
pattern is that of the statistics record in ``prepare.store.stats``: a file snapshot decides
whether the cache is current, the dataset lock guards the sidecar write for a bounded wait, and
a held lock only costs the cache. The hashes themselves stay in the signature file, never in
the sidecar, which every other command reads.
"""

import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Sequence, Tuple, Union

from prepare.core.exceptions import DataAccessError
from prepare.data.file_io import copy_file
from prepare.store.locks import dataset_lock
from prepare.store.sidecar import read_sidecar, write_sidecar
from prepare.store.stats import STATS_LOCK_WAIT_SECONDS, _file_signature, _store_location

logger = logging.getLogger(__name__)

SIGNATURE_SUFFIX = ".sig.gz"


def store_signature_path(sidecar_path: Union[str, Path]) -> Path:
    """``<ACC>.sig.gz`` beside the sidecar at ``sidecar_path``."""
    path = Path(sidecar_path)
    return path.parent / f"{path.parent.name}{SIGNATURE_SUFFIX}"


def cached_sketch(
    acc_dir: Union[str, Path], sidecar_path: Optional[Union[str, Path]], files: Sequence[Path]
) -> Optional[Tuple[Path, Dict[str, Any]]]:
    """The cached signature path and its sidecar block when they still describe ``files``.

    None when there is no sidecar path, no readable sidecar, no ``sketch`` block, the block's
    file snapshot differs from the files now, the block was built with other k sizes or
    another scaled value than ``sra_sketch`` uses now, or the signature file is missing.
    """
    # Imported here: the sketch parameters live with the sketching code in prepare.sra.
    from prepare.sra.sketch import KSIZES, SCALED

    if sidecar_path is None:
        return None
    sidecar = read_sidecar(sidecar_path)
    block = sidecar.sketch if sidecar is not None else None
    if not block or not block.get("file"):
        return None
    path = Path(sidecar_path).parent / str(block["file"])
    if not _block_describes(block, files, KSIZES, SCALED) or not path.is_file():
        return None
    return path, block


def _block_describes(block: Dict[str, Any], files: Sequence[Path], ksizes: Sequence[int], scaled: int) -> bool:
    """Whether a cached ``sketch`` block was built from ``files`` as they are now, with these parameters."""
    return (
        block.get("files") == _file_signature([Path(f) for f in files])
        and list(block.get("ksizes") or []) == list(ksizes)
        and block.get("scaled") == scaled
    )


def store_sketch(
    sidecar_path: Optional[Union[str, Path]],
    signature_file: Union[str, Path],
    summary: Dict[str, Any],
    should_stop: Optional[Callable[[], bool]] = None,
) -> None:
    """Copy ``signature_file`` beside the sidecar and record ``summary`` as its ``sketch`` block.

    A no-op without a sidecar (``sidecar_path`` None or not a file). Inside a store the copy
    and the read-modify-write of the sidecar run under the dataset lock with the statistics
    wait; a lock held by a running download leaves the cache unwritten, with a warning, since
    the next ``sra_sketch`` run computes the sketch again. ``should_stop`` (the run's stop
    token) ends that wait early, with the same warning.
    """
    if sidecar_path is None or not Path(sidecar_path).is_file():
        return
    located = _store_location(Path(sidecar_path))
    if located is None:
        _write_sketch(Path(sidecar_path), Path(signature_file), summary)
        return
    paths, accession = located
    try:
        with dataset_lock(paths, accession, wait_seconds=STATS_LOCK_WAIT_SECONDS, should_stop=should_stop):
            _write_sketch(Path(sidecar_path), Path(signature_file), summary)
    except DataAccessError as e:
        logger.warning("Not caching sketch for %s: %s", accession, e)


def _write_sketch(sidecar_path: Path, signature_file: Path, summary: Dict[str, Any]) -> None:
    """Copy the signature into place, then re-read the sidecar and set only its ``sketch`` block."""
    target = store_signature_path(sidecar_path)
    copy_file(signature_file, target)
    sidecar = read_sidecar(sidecar_path)
    if sidecar is None:
        return
    sidecar.sketch = {**summary, "file": target.name, "computed": datetime.now(timezone.utc).isoformat()}
    write_sidecar(sidecar_path, sidecar)
