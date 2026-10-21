"""
Containment of a genome in local dataset sketches, written in Branchwater's CSV layout.

``containment_rows`` compares the genome's k = 21, scaled 1000 sketch (the one
``prepare.data.branchwater_search.sketch_fasta`` or ``load_signature`` returns) with each
dataset signature ``sra_sketch`` wrote, as the fraction of the genome's hashes present in the
dataset. ``merge_containment_csv`` writes the rows into ``<matches>/<genome_id>.csv`` with the
columns of a Branchwater result, keeping the rows of other accessions, so ``parse_containment``
and ``select_datasets`` read the file unchanged. The file's stem stays the genome id because
``parse_containment_data`` keys genomes on it and keeps the maximum of duplicate rows.
"""

import csv
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple, Union

from prepare.core.exceptions import DataAccessError
from prepare.core.optional import require
from prepare.data.branchwater_search import BRANCHWATER_COLUMNS, KSIZE
from prepare.data.file_io import open_atomic
from prepare.sra.sketch import load_dataset_minhash
from prepare.utils.progress import active_item_level

logger = logging.getLogger(__name__)


def genome_minhash(signature: Mapping[str, Any]) -> Any:
    """The MinHash of a sourmash signature given as the JSON object ``sketch_fasta`` returns."""
    signature_module = require("sourmash.signature", "sourmash", "Comparing sketches")
    return signature_module.load_one_signature_from_json(json.dumps([signature])).minhash


def containment_rows(genome_mh: Any, signatures: Mapping[str, Union[str, Path]]) -> List[Tuple[str, float]]:
    """``(accession, containment)`` for each dataset signature, in the order given.

    A signature that cannot be read is logged and left out, so one bad file does not stop the
    comparison of the others.
    """
    rows: List[Tuple[str, float]] = []
    for accession, path in signatures.items():
        try:
            dataset_mh = load_dataset_minhash(path, ksize=genome_mh.ksize)
        except DataAccessError as e:
            logger.warning("No containment for %s: %s", accession, e)
            continue
        value = float(genome_mh.contained_by(dataset_mh))
        logger.log(active_item_level(), "%s: containment %.4f", accession, value)
        rows.append((accession, value))
    return rows


def merge_containment_csv(path: Union[str, Path], rows: Sequence[Tuple[str, float]], ksize: int = KSIZE) -> Path:
    """Write ``rows`` into the matches CSV at ``path``, replacing the rows of the same accessions.

    Rows of other accessions, and every column they carry, are kept; the header is the existing
    one extended with any missing Branchwater column, or ``BRANCHWATER_COLUMNS`` for a new file.
    New rows carry ``containment`` to four decimals and ``cANI`` as containment to the power
    ``1 / ksize``, as ``write_branchwater_csv`` does; their metadata columns are empty.
    """
    target = Path(path)
    fieldnames, existing = _read_existing(target)
    replaced = {accession for accession, _ in rows}
    kept = [row for row in existing if row.get("acc") not in replaced]
    for accession, containment in rows:
        cani = round(containment ** (1 / ksize), 4) if containment > 0 else 0.0
        kept.append({"acc": accession, "containment": f"{containment:.4f}", "cANI": str(cani)})
    kept.sort(key=lambda row: -_as_float(row.get("containment")))
    target.parent.mkdir(parents=True, exist_ok=True)
    with open_atomic(target, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, restval="")
        writer.writeheader()
        writer.writerows(kept)
    return target


def _read_existing(target: Path) -> Tuple[List[str], List[Dict[str, Any]]]:
    """The header (extended with the Branchwater columns) and rows of ``target``, or the defaults."""
    if not target.is_file():
        return list(BRANCHWATER_COLUMNS), []
    with open(target, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    if "acc" not in fieldnames or "containment" not in fieldnames:
        raise DataAccessError(f"{target} is not a matches file (needs 'acc' and 'containment' columns)")
    for column in BRANCHWATER_COLUMNS:
        if column not in fieldnames:
            fieldnames.append(column)
    return fieldnames, rows


def _as_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0
