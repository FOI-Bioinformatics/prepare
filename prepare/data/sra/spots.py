"""Expected spot counts of an SRA accession and the completeness verdict a read count gives against them.

A download's completeness verdict compares the reads on disk with NCBI's recorded spot count
(``Run_Total_Spots``). That count can be known in several places: the project registry's metadata
block, the store sidecar's ``ncbi`` block, a metadata XML file on disk, or an earlier download
verdict. ``expected_spots`` is the one lookup that consults them in a fixed order, so every caller
gets the same count for the same accession, and ``merged_verdict`` is the one rule for combining a
new verdict with the recorded one: a known verdict (``complete`` or ``truncated``) is never replaced
by ``unverified``, since the absence of a count is not evidence about the reads either way.
"""

import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Iterable, Mapping, Optional, Union

from prepare.core.exceptions import DataAccessError
from prepare.data import registry_blocks as rb
from prepare.data.sra.fastq import COMPLETE_RATIO_THRESHOLD

if TYPE_CHECKING:
    from prepare.data.registry import Registry
    from prepare.store.layout import StorePaths

logger = logging.getLogger(__name__)

_UNVERIFIED = "unverified"


def positive_int(value: Any) -> Optional[int]:
    """``value`` as an int when it is a positive whole number (or its string form), else None.

    A spot count of zero is treated as unknown, as ``verify_download`` has always done.
    """
    if isinstance(value, bool):
        return None
    try:
        number = int(value)
    except (TypeError, ValueError):
        return None
    return number if number > 0 else None


def verdict_for_count(reads_r1: Optional[int], expected_spots: Optional[int]) -> Dict[str, Any]:
    """The completeness verdict of ``reads_r1`` reads against ``expected_spots`` NCBI spots.

    ``"complete"`` when the ratio of reads to spots is at least ``COMPLETE_RATIO_THRESHOLD``,
    ``"truncated"`` below it, and ``"unverified"`` when either count is unknown (None) or the spot
    count is zero. The ratio is rounded to four decimals. Returns a dict with the keys
    ``method``, ``ratio``, ``verdict``, ``expected_spots`` and ``reads_r1``.
    """
    ratio: Optional[float]
    if expected_spots and reads_r1 is not None:
        ratio = round(reads_r1 / expected_spots, 4)
        verdict = "complete" if ratio >= COMPLETE_RATIO_THRESHOLD else "truncated"
        method = "spots"
    else:
        ratio = None
        verdict = _UNVERIFIED
        method = _UNVERIFIED
    return {
        "method": method,
        "ratio": ratio,
        "verdict": verdict,
        "expected_spots": expected_spots,
        "reads_r1": reads_r1,
    }


def spots_from_xml(accession: str, folders: Iterable[Union[str, Path]]) -> Optional[int]:
    """``Run_Total_Spots`` from the first ``<accession>_metadata.xml`` in ``folders`` that records one.

    A folder without the file is skipped. A file that cannot be read or parsed (``OSError``,
    ``ValueError``, ``xml.etree.ElementTree.ParseError``) is logged and skipped. Returns None when
    no folder yields a positive spot count.
    """
    # Imported here: prepare.data.metadata pulls in pandas, Biopython and lxml, which the
    # download path does not otherwise need at import time.
    from prepare.data.metadata import parse_metadata_xml

    for folder in folders:
        xml_path = Path(folder) / f"{accession}_metadata.xml"
        if not xml_path.is_file():
            continue
        try:
            parsed = parse_metadata_xml(xml_path)
        except (OSError, ValueError, ET.ParseError) as e:
            logger.warning("Could not read the spot count of %s from %s: %s", accession, xml_path, e)
            continue
        spots = positive_int((parsed or {}).get("Run_Total_Spots"))
        if spots is not None:
            return spots
    return None


def _sidecar_spots(accession: str, store: Union["StorePaths", str, Path]) -> Optional[int]:
    """``ncbi.spots`` from ``accession``'s store sidecar, or None when there is no readable sidecar."""
    # Imported here: prepare.store imports prepare.data.sra, so a module-level import would
    # form a cycle.
    from prepare.store.layout import StorePaths, sidecar_path, store_paths
    from prepare.store.sidecar import read_sidecar

    paths = store if isinstance(store, StorePaths) else store_paths(Path(store))
    path = sidecar_path(paths, accession)
    if not path.is_file():
        return None
    try:
        sidecar = read_sidecar(path)
    except DataAccessError as e:
        logger.warning("Ignoring the sidecar of %s: %s", accession, e)
        return None
    if sidecar is None or not isinstance(sidecar.ncbi, dict):
        return None
    return positive_int(sidecar.ncbi.get("spots"))


def expected_spots(
    registry: Optional["Registry"],
    accession: str,
    store: Optional[Union["StorePaths", str, Path]] = None,
    xml_folders: Iterable[Union[str, Path]] = (),
) -> Optional[int]:
    """NCBI's recorded spot count for ``accession``, or None when no source records one.

    The sources are consulted in this order, and the first positive count wins:

    1. the registry's metadata block (``run_total_spots``);
    2. ``ncbi.spots`` in the store sidecar, when ``store`` (a ``StorePaths`` or a store root) is given;
    3. ``Run_Total_Spots`` in ``<accession>_metadata.xml`` in each of ``xml_folders``, in order;
    4. the ``expected_spots`` of the download verdict the registry already records.

    The first three read NCBI's count directly; the recorded verdict holds a copy made by an earlier
    run, so it comes last and a fresh sidecar or XML corrects a wrong copy.
    """
    if registry is not None:
        metadata = rb.metadata_block(registry, accession)
        if metadata is not None:
            spots = positive_int(metadata.run_total_spots)
            if spots is not None:
                return spots
    if store is not None:
        spots = _sidecar_spots(accession, store)
        if spots is not None:
            return spots
    spots = spots_from_xml(accession, xml_folders)
    if spots is not None or registry is None:
        return spots
    previous = rb.download_verdict(registry, accession)
    return positive_int(previous.expected_spots) if previous is not None else None


def _as_dict(verdict: Any) -> Optional[Dict[str, Any]]:
    """A plain-dict copy of ``verdict`` (a mapping or a registry ``Verdict``), or None."""
    if verdict is None:
        return None
    if isinstance(verdict, rb.RegistryBlock):
        return verdict.to_dict()
    if isinstance(verdict, Mapping):
        return dict(verdict)
    raise TypeError(f"Not a verdict: {verdict!r}")


def merged_verdict(
    previous: Optional[Union[Mapping[str, Any], "rb.Verdict"]],
    new: Optional[Union[Mapping[str, Any], "rb.Verdict"]],
    reads_r1: Optional[int],
    expected: Optional[int],
) -> Optional[Dict[str, Any]]:
    """The completeness verdict to record, given the ``previous`` one and a ``new`` one.

    When both ``reads_r1`` and a spot count ``expected`` are known, the verdict is recomputed from
    them with ``verdict_for_count``, which may turn a recorded verdict into the other known one
    (``truncated`` into ``complete`` or the reverse). Otherwise ``new`` is taken, except that a known
    ``previous`` verdict (``complete`` or ``truncated``) is never replaced by an ``unverified`` (or
    missing) one: without a count there is no evidence about the reads either way, e.g. a relink to
    a store copy whose sidecar records no read count. A missing ``new`` keeps ``previous``. Returns
    None when neither is known. The result is always a new dict; the arguments are not modified.
    """
    if reads_r1 is not None and positive_int(expected) is not None:
        return verdict_for_count(reads_r1, positive_int(expected))
    before = _as_dict(previous)
    after = _as_dict(new)
    if after is None:
        return before
    known_before = before is not None and before.get("verdict") in ("complete", "truncated")
    if known_before and after.get("verdict") in (None, _UNVERIFIED):
        return before
    return after


def read_count(value: Any) -> Optional[int]:
    """``value`` as a read count: a whole number of zero or more, else None (unknown).

    A bool is not a count, and neither is a float or a numeric string with a fraction.
    """
    if isinstance(value, bool) or isinstance(value, float):
        return None
    try:
        number = int(value)
    except (TypeError, ValueError):
        return None
    return number if number >= 0 else None
