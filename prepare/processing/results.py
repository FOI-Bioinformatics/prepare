"""The consolidated results table: one row per screened (accession, genome) pair.

Joins what the project registry records for each pair (screening, selection, exclusion,
download, run metadata, the dataset quality summary of ``sra_profile``/``sra_report``, read
extraction, reference coverage and assembly) with the containment values of the parsed
containment table, which are unrounded and not limited by ``cap_screening``. Three columns
before the end are the seconds the download, the extraction and the assembly took, when the
registry records them; ``download_verdict`` is the accession's recorded download completeness
verdict (``complete``, ``truncated`` or ``unverified``), empty when none was recorded.

The last eight columns were added after ``download_verdict`` and are empty wherever the
registry has nothing to report: ``quality_source`` names which analysis supplied the
``total_reads``/``gc_percent``/``quality_grade`` columns (``profile``, ``report``, ``legacy``,
or empty when neither ran); ``assembly_largest``, ``assembly_n90``, ``assembly_gc_percent``
(the assembly's GC fraction as a percent, two decimals) and ``assembly_contigs_ge_1kb`` are
extra contig statistics beyond ``contigs``/``total_bp``/``n50``, empty without an assembly;
``assembly_mean_depth_estimate`` is empty unless the assembly's coverage mapping was run;
``assembly_dir`` and ``coverage_tsv`` are the assembly folder and the coverage table's path,
each relative to the project root (the registry file's folder) like every other path column.
"""

from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from prepare.core.utils import _KNOWN_METADATA_COLUMNS
from prepare.data import registry_blocks as rb
from prepare.data.registry import Registry, to_int_or_none

RESULTS_COLUMNS = [
    "accession",
    "genome_id",
    "containment",
    "selected",
    "excluded",
    "exclusion_reason",
    "download_state",
    "run_total_spots",
    "run_size",
    "total_reads",
    "gc_percent",
    "quality_grade",
    "mapped_reads",
    "mapping_rate_to_reference",
    "breadth",
    "mean_depth",
    "contigs",
    "total_bp",
    "n50",
    "genome_fraction_estimate",
    "assembly_mapping_rate",
    "download_seconds",
    "extraction_seconds",
    "assembly_seconds",
    "download_verdict",
    "quality_source",
    "assembly_largest",
    "assembly_n90",
    "assembly_gc_percent",
    "assembly_contigs_ge_1kb",
    "assembly_mean_depth_estimate",
    "assembly_dir",
    "coverage_tsv",
]

Pair = Tuple[str, str]

# Empty blocks for datasets that lack one, built once and only read: most pairs have no
# extraction, and a fresh block per pair costs more than the rest of the row together.
_NO_EXCLUSION = rb.ExclusionBlock()
_NO_METADATA = rb.MetadataBlock()
_NO_SELECTION = rb.SelectionBlock()
_NO_EXTRACTION = rb.ExtractionBlock()


def _positive_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if pd.isna(number) or number <= 0:
        return None
    return number


def _assembly_gc_percent(gc: Any) -> Optional[float]:
    """``gc`` (the assembly's 0-1 GC fraction from ``extra``) as a percent, two decimals, or None."""
    try:
        return round(float(gc) * 100, 2)
    except (TypeError, ValueError):
        return None


def screened_pairs(registry: Registry, parsed_table: Optional[pd.DataFrame] = None) -> Dict[Pair, Optional[float]]:
    """Every (accession, genome) pair with its containment, or None when it was never screened.

    Pairs come from three sources: the registry's ``screening.genomes`` entries; every genome
    column of ``parsed_table`` (all columns except ``max_containment`` and
    ``max_containment_annotation``) with a value above 0; and extraction records with no
    screening entry (containment None). Where the registry and the table both hold a value,
    the table's is used, since the registry rounds it to four decimals.
    """
    pairs: Dict[Pair, Optional[float]] = {}
    for accession in registry.datasets:
        # One field per entry, so it is read straight from the registry dict (see ``rb.raw``).
        for genome_id, entry in (rb.raw(registry, accession, "screening", "genomes") or {}).items():
            value = entry.get("containment") if isinstance(entry, dict) else None
            pairs[(accession, genome_id)] = float(value) if value is not None else None
    if parsed_table is not None:
        # Genome column name -> position of its first column: a repeated name is read once, by
        # position, since selecting it by name would return a frame rather than one column.
        genome_columns: Dict[Any, int] = {}
        for index, name in enumerate(parsed_table.columns):
            if name not in _KNOWN_METADATA_COLUMNS:
                genome_columns.setdefault(name, index)
        # Plain column lists, walked row by row as before, instead of one pandas Series per row.
        columns = [(str(name), parsed_table.iloc[:, index].tolist()) for name, index in genome_columns.items()]
        for position, accession in enumerate(parsed_table.index):
            for genome_id, values in columns:
                value = _positive_float(values[position])
                if value is not None:
                    pairs[(str(accession), genome_id)] = value
    for accession in registry.datasets:
        for genome_id, entry in (registry.datasets[accession].get("extractions") or {}).items():
            if isinstance(entry, dict):
                pairs.setdefault((accession, genome_id), None)
    return pairs


def _mapping_rate(mapped_reads: Optional[int], spots: Optional[int]) -> Optional[float]:
    """Mapped reads over the run's spot count, or None when either is missing or zero.

    ``mapped_reads`` counts BAM records that passed the extraction filters, while a spot is one
    read or one read pair, so paired-end data can give a value above 1. The ratio is reported
    as is, without halving, since whether both mates mapped is not known here. It is rounded
    to four decimals, like the other ratios the registry records.
    """
    if not mapped_reads or not spots or mapped_reads <= 0 or spots <= 0:
        return None
    return round(mapped_reads / spots, 4)


def _dataset_fields(registry: Registry, accession: str) -> Tuple[Dict[str, Any], Optional[str]]:
    """The columns of a results row that depend on the accession only, shared by all its genomes.

    Returned as ``(fields, quality_source)`` rather than one dict with ``quality_source`` folded
    in, so that ``_row`` can place it after ``download_verdict`` and a raw row's key order stays
    the same as ``RESULTS_COLUMNS``. This is for the readability of the raw rows only:
    ``results_dataframe`` orders the output columns by ``RESULTS_COLUMNS`` in any case.
    """
    exclusion = rb.exclusion_block(registry, accession) or _NO_EXCLUSION
    excluded = bool(exclusion.excluded)
    metadata = rb.metadata_block(registry, accession) or _NO_METADATA
    profile, quality_source = rb.quality_summary(registry, accession)
    fields = {
        "selected": bool((rb.selection_block(registry, accession) or _NO_SELECTION).selected),
        "excluded": excluded,
        "exclusion_reason": (exclusion.reason or None) if excluded else None,
        # Read raw, like download_verdict: a hand-edited verdict that is not a mapping must not
        # stop the table, and the state is the only field of the download block used here.
        "download_state": rb.raw(registry, accession, "download", "state") or None,
        "run_total_spots": to_int_or_none(metadata.run_total_spots),
        "run_size": to_int_or_none(metadata.run_size),
        "total_reads": to_int_or_none(profile["total_reads"]),
        "gc_percent": profile["gc_percent"],
        "quality_grade": profile["quality_grade"],
    }
    return fields, quality_source


def _row(
    registry: Registry,
    accession: str,
    genome_id: str,
    containment: Optional[float],
    dataset: Dict[str, Any],
    quality_source: Optional[str],
) -> Dict[str, Any]:
    extraction = rb.extraction_block(registry, accession, genome_id) or _NO_EXTRACTION
    # A missing assembly reads as None in every column, unlike a recorded one whose stats are 0.
    assembly = extraction.assembly.to_dict() if extraction.assembly is not None else {}
    mapped_reads = to_int_or_none(extraction.mapped_reads)
    return {
        "accession": accession,
        "genome_id": genome_id,
        "containment": containment,
        **dataset,
        "mapped_reads": mapped_reads,
        "mapping_rate_to_reference": _mapping_rate(mapped_reads, dataset["run_total_spots"]),
        "breadth": extraction.breadth,
        "mean_depth": extraction.mean_depth,
        "contigs": assembly.get("contigs"),
        "total_bp": assembly.get("total_bp"),
        "n50": assembly.get("n50"),
        "genome_fraction_estimate": assembly.get("genome_fraction_estimate"),
        "assembly_mapping_rate": assembly.get("mapping_rate"),
        # One field, read straight from the registry dict (see rb.raw) rather than cached per accession.
        "download_seconds": rb.raw(registry, accession, "download", "seconds"),
        "extraction_seconds": extraction.seconds,
        "assembly_seconds": assembly.get("seconds"),
        # Same rb.raw pattern as download_seconds; appended at the end, not grouped with the
        # other download fields, so existing column positions are kept.
        "download_verdict": _download_verdict(registry, accession),
        # Per accession, but placed here rather than in `dataset` so the raw row's key order
        # matches RESULTS_COLUMNS (results_dataframe orders the output columns explicitly anyway).
        "quality_source": quality_source,
        "assembly_largest": assembly.get("largest"),
        "assembly_n90": assembly.get("n90"),
        "assembly_gc_percent": _assembly_gc_percent(assembly.get("gc")),
        "assembly_contigs_ge_1kb": assembly.get("contigs_ge_1kb"),
        "assembly_mean_depth_estimate": assembly.get("mean_depth_estimate"),
        "assembly_dir": assembly.get("dir"),
        "coverage_tsv": extraction.coverage_tsv,
    }


def _download_verdict(registry: Registry, accession: str) -> Optional[str]:
    """The recorded download verdict of ``accession``, or None when none is recorded.

    A ``complete`` value that is not a mapping (a hand-edited registry) reads as no verdict, as in
    ``prepare.processing.status_report.download_verdicts``.
    """
    complete = rb.raw(registry, accession, "download", "complete")
    return complete.get("verdict") if isinstance(complete, dict) else None


def results_rows(
    registry: Registry,
    parsed_table: Optional[pd.DataFrame] = None,
    genome_id: Optional[str] = None,
    min_containment: float = 0.0,
) -> List[Dict[str, Any]]:
    """One dict per (accession, genome) pair, keyed by ``RESULTS_COLUMNS`` in that order.

    ``genome_id`` keeps only that genome's pairs. ``min_containment`` keeps pairs whose
    containment is at least that value; above 0 it also drops pairs with no containment
    (extracted but never screened), since those cannot be shown to meet it. Rows are sorted
    by decreasing containment, with unknown containment last, then by accession and genome.
    """
    rows = []
    # Built once per accession and reused for each of its genomes.
    datasets: Dict[str, Tuple[Dict[str, Any], Optional[str]]] = {}
    for (accession, genome), containment in screened_pairs(registry, parsed_table).items():
        if genome_id is not None and genome != genome_id:
            continue
        if min_containment > 0 and (containment is None or containment < min_containment):
            continue
        if accession not in datasets:
            datasets[accession] = _dataset_fields(registry, accession)
        dataset, quality_source = datasets[accession]
        rows.append(_row(registry, accession, genome, containment, dataset, quality_source))
    rows.sort(
        key=lambda row: (
            row["containment"] is None,
            -(row["containment"] or 0.0),
            row["accession"],
            row["genome_id"],
        )
    )
    return rows


def results_dataframe(rows: List[Dict[str, Any]]) -> pd.DataFrame:
    """The rows as a DataFrame with exactly ``RESULTS_COLUMNS`` (header only when ``rows`` is empty).

    Built with object dtype so integer columns with missing values stay integers in the
    written table instead of turning into floats.
    """
    return pd.DataFrame(rows, columns=RESULTS_COLUMNS, dtype=object)
