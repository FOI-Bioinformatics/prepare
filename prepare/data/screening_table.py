"""Screening entries recorded from a parsed containment table in one pass.

``parse_containment`` writes a wide table (one row per accession, one column per genome) and
records every positive cell as a screening entry in the project registry. Recording the cells one
at a time rebuilt an accession's screening block once per genome and stamped each cell with its
own time; capping afterwards rebuilt every block once per genome again. Here the positive cells
are found with numpy, the per-genome cap is applied before anything is written, and each
accession's block is converted and written once, with one timestamp for the whole run. The result
equals recording the cells one by one and then calling ``cap_screening`` for each genome.

It lives next to ``prepare.data.registry``, which is held at its current size by the module
size ceiling; ``registry.record_screening_from_table`` delegates here.
"""

import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd

from prepare.data import registry as reg
from prepare.data import registry_blocks as rb

logger = logging.getLogger(__name__)

_NOT_GENOMES = ("max_containment", "max_containment_annotation")

Cell = Tuple[int, int]


def _positive_cells(table: pd.DataFrame) -> Tuple[List[str], List[str], np.ndarray, np.ndarray, List[float]]:
    """Accessions, genomes, and the row, column and rounded value of every cell above 0.

    Cells that do not hold a number are skipped. Cells come in row order, then column order.
    """
    genomes = [str(column) for column in table.columns if column not in _NOT_GENOMES]
    frame = table[[column for column in table.columns if column not in _NOT_GENOMES]]
    values = frame.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    rows, cols = np.nonzero(values > 0)
    rounded = [round(value, 4) for value in values[rows, cols].tolist()]
    return [str(accession) for accession in table.index], genomes, rows, cols, rounded


def _earlier_entries(registry: reg.Registry, genomes: Set[str]) -> Dict[str, Dict[str, float]]:
    """Genome -> {accession: containment} for the screening entries already in the registry."""
    earlier: Dict[str, Dict[str, float]] = defaultdict(dict)
    for accession, record in registry.datasets.items():
        recorded = (record.get("screening") or {}).get("genomes") or {}
        for genome in genomes.intersection(recorded):
            earlier[genome][accession] = float((recorded[genome] or {}).get("containment") or 0.0)
    return earlier


def _apply_cap(
    registry: reg.Registry,
    accessions: List[str],
    genomes: List[str],
    cells: Dict[int, List[Tuple[int, float]]],
    max_screened: int,
) -> Tuple[Set[Cell], List[Tuple[str, str]]]:
    """Which new cells to drop, and which earlier (accession, genome) entries to drop, per genome.

    Ranks new and earlier entries together as ``cap_screening`` does: by containment, highest
    first, ties in registry order, with accessions new to the registry after the existing ones
    in table order. Logs the same warning ``cap_screening`` does for each genome it trims.
    """
    earlier = _earlier_entries(registry, set(genomes))
    position = {accession: index for index, accession in enumerate(registry.datasets)}
    dropped_cells: Set[Cell] = set()
    dropped_earlier: List[Tuple[str, str]] = []
    for col, genome in enumerate(genomes):
        new = cells.get(col, [])
        new_accessions = {accessions[row] for row, _ in new}
        kept_earlier = {acc: value for acc, value in earlier.get(genome, {}).items() if acc not in new_accessions}
        total = len(new_accessions) + len(kept_earlier)
        if total <= max_screened:
            continue
        ranked: List[Tuple[float, int, Any]] = [
            (-value, position.get(accessions[row], len(position) + row), (row, col)) for row, value in new
        ]
        ranked += [(-value, position[acc], acc) for acc, value in kept_earlier.items()]
        ranked.sort(key=lambda item: (item[0], item[1]))
        for _, _, entry in ranked[max_screened:]:
            if isinstance(entry, tuple):
                dropped_cells.add(entry)
            else:
                dropped_earlier.append((entry, genome))
        logger.warning(
            "Kept the %d highest containments for %s in the registry and dropped %d more "
            "(the limit is --registry-max-screened); the match CSV keeps them all",
            max_screened,
            genome,
            total - max_screened,
        )
    return dropped_cells, dropped_earlier


def _store_screening(registry: reg.Registry, accession: str, screening: rb.ScreeningBlock) -> None:
    """Write ``screening`` back; with no genome left, remove the block, and the record if that empties it."""
    if screening.genomes:
        rb.set_screening_block(registry, accession, screening)
        return
    record = registry.datasets.get(accession)
    if record is None:
        return
    record.pop("screening", None)
    if not record:
        del registry.datasets[accession]


def _drop_screening(registry: reg.Registry, accession: str, genome: str) -> None:
    """Remove one genome's screening entry, through the typed block."""
    screening = rb.screening_block(registry, accession) or rb.ScreeningBlock()
    screening.genomes.pop(genome, None)
    _store_screening(registry, accession, screening)


def _write_row(
    registry: reg.Registry,
    accession: str,
    entries: List[Tuple[str, float, bool]],
    date: str,
    csv_paths: Dict[str, str],
) -> None:
    """Write one accession's screening block once: kept entries set, dropped ones removed."""
    existing = accession in registry.datasets
    if not existing and not any(keep for _, _, keep in entries):
        return
    screening = rb.screening_block(registry, accession) or rb.ScreeningBlock()
    screening.date = date
    screening.discard("inferred")
    for genome, value, keep in entries:
        if keep:
            screening.genomes[genome] = rb.ScreeningEntry(
                containment=value, cani=None, csv=csv_paths[genome], source="matches"
            )
        else:
            screening.genomes.pop(genome, None)
    _store_screening(registry, accession, screening)


def read_screening_table(path: Union[str, Path]) -> Optional[pd.DataFrame]:
    """The parsed containment table written at ``path``, or None when there is no such file.

    Read before a registry transaction is opened, so the lock is never held while a large
    table is parsed; a missing file is logged at debug level.
    """
    if not Path(path).exists():
        logger.debug("Parsed containment table %s does not exist; nothing to record", path)
        return None
    return pd.read_csv(Path(path), sep="\t", index_col=0)


def record_screening_table(
    registry: reg.Registry,
    table: Union[pd.DataFrame, str, Path],
    matches_folder: Union[str, Path],
    max_screened: int,
) -> int:
    """Record every positive cell of ``table`` as a screening entry; return how many were recorded.

    ``table`` is the parsed containment DataFrame, or the path it was written to, read here (a
    path that does not exist is logged at debug level and records nothing). Keeps at most
    ``max_screened`` accessions per genome, counting entries already in the registry. The count
    returned includes the cells the cap then dropped, as recording them one by one did. Every
    block written gets the same timestamp.
    """
    if not isinstance(table, pd.DataFrame):
        loaded = read_screening_table(table)
        if loaded is None:
            return 0
        table = loaded
    accessions, genomes, rows, cols, rounded = _positive_cells(table)
    if not len(rows):
        return 0
    cells: Dict[int, List[Tuple[int, float]]] = defaultdict(list)
    for row, col, value in zip(rows.tolist(), cols.tolist(), rounded):
        cells[col].append((row, value))
    dropped_cells, dropped_earlier = _apply_cap(registry, accessions, genomes, cells, max_screened)

    date = reg._now()
    csv_paths = {genome: str(Path(matches_folder) / f"{genome}.csv") for genome in genomes}
    row_list, col_list = rows.tolist(), cols.tolist()
    bounds = np.searchsorted(rows, np.arange(len(accessions) + 1)).tolist()
    for row in dict.fromkeys(row_list):
        entries = [
            (genomes[col_list[i]], rounded[i], (row, col_list[i]) not in dropped_cells)
            for i in range(bounds[row], bounds[row + 1])
        ]
        _write_row(registry, accessions[row], entries, date, csv_paths)
        for genome, _, _ in entries:
            registry.genomes.setdefault(genome, {})
    for accession, genome in dropped_earlier:
        _drop_screening(registry, accession, genome)
    return len(row_list)
