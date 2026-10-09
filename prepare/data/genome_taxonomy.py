"""Genome taxonomy enrichment via the GTDB API.

This module maps genome accessions (GCF_*/GCA_*) to their taxonomic
classification by querying the GTDB API, and provides utilities for
annotating containment data with taxonomy information.
"""

import csv
import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd
import requests

from prepare.core.exceptions import DataAccessError
from prepare.core.models import TaxonomyInfo
from prepare.data.file_io import open_atomic
from prepare.data.gtdb import GTDB_API_BASE, REQUEST_TIMEOUT, get_session

logger = logging.getLogger(__name__)

# Mapping from GTDB rank prefixes to TaxonomyInfo field names
_GTDB_RANK_MAP = {
    "d": None,  # domain -- not stored
    "p": "phylum",
    "c": "class_name",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}

_CACHE_COLUMNS = [
    "genome_id",
    "species",
    "genus",
    "family",
    "order",
    "class_name",
    "phylum",
    "organism",
    "tax_id",
]


def parse_gtdb_taxonomy_string(gtdb_string: str, genome_id: str) -> TaxonomyInfo:
    """Parse a GTDB taxonomy string into a TaxonomyInfo dataclass.

    Expects the standard GTDB format, e.g.
    ``d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;
    g__Bacillus;s__Bacillus subtilis``
    """
    info = TaxonomyInfo(genome_id=genome_id)
    if not gtdb_string:
        return info

    for token in gtdb_string.split(";"):
        token = token.strip()
        if "__" not in token:
            continue
        prefix, value = token.split("__", 1)
        field = _GTDB_RANK_MAP.get(prefix)
        if field and value:
            setattr(info, field, value)

    return info


def _taxonomy_from_record(record: dict, accession: str) -> Optional[TaxonomyInfo]:
    """Build a TaxonomyInfo from a GTDB record dict, or None if it has no taxonomy."""
    gtdb_taxonomy = record.get("gtdb_taxonomy") or record.get("gtdbTaxonomy") or record.get("taxonomy") or ""
    if not gtdb_taxonomy:
        return None

    info = parse_gtdb_taxonomy_string(gtdb_taxonomy, accession)
    info.organism = record.get("organism_name") or record.get("organismName")
    info.tax_id = str(record.get("ncbi_taxid") or record.get("taxId") or "") or None
    return info


def _extract_search_rows(data) -> list:
    """Normalize a GTDB search response into a list of record dicts.

    The endpoint may return a bare list, a dict under a ``rows`` key, or a
    dict under a ``results`` key.
    """
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        return data.get("rows") or data.get("results") or []
    return []


def _lookup_genome_taxonomy_gtdb(accession: str) -> Optional[TaxonomyInfo]:
    """Query the GTDB API for a single genome's taxonomy.

    Tries the genome detail endpoint first, then falls back to the search
    endpoint. Returns None when no taxonomy can be resolved.
    """
    quoted = requests.utils.quote(accession)
    genome_url = f"{GTDB_API_BASE}/genome/{quoted}"
    search_url = f"{GTDB_API_BASE}/search/gtdb"
    logger.debug("Querying GTDB genome endpoint: %s", genome_url)

    try:
        session = get_session()
        response = session.get(genome_url, timeout=REQUEST_TIMEOUT)
        if response.status_code == 200:
            info = _taxonomy_from_record(response.json(), accession)
            if info:
                return info

        # Fallback: search endpoint
        params: Dict[str, Union[str, int]] = {"search": accession, "page": 1, "itemsPerPage": 1}
        response = session.get(search_url, params=params, timeout=REQUEST_TIMEOUT)
        if response.status_code == 200:
            rows = _extract_search_rows(response.json())
            if rows:
                return _taxonomy_from_record(rows[0], accession)

    except requests.exceptions.RequestException as e:
        raise DataAccessError(f"GTDB API error looking up genome '{accession}': {e}")

    return None


def load_taxonomy_cache(cache_file: Path) -> Dict[str, TaxonomyInfo]:
    """Load a genome-to-taxonomy cache from a TSV file."""
    cache: Dict[str, TaxonomyInfo] = {}
    if not cache_file.exists():
        return cache

    with open(cache_file, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid = row.get("genome_id", "")
            if not gid:
                continue
            cache[gid] = TaxonomyInfo(
                genome_id=gid,
                species=row.get("species") or None,
                genus=row.get("genus") or None,
                family=row.get("family") or None,
                order=row.get("order") or None,
                class_name=row.get("class_name") or None,
                phylum=row.get("phylum") or None,
                organism=row.get("organism") or None,
                tax_id=row.get("tax_id") or None,
            )
    logger.info("Loaded %d entries from taxonomy cache %s", len(cache), cache_file)
    return cache


def save_taxonomy_cache(taxonomy: Dict[str, TaxonomyInfo], cache_file: Path) -> None:
    """Save a genome-to-taxonomy mapping to a TSV file."""
    with open_atomic(cache_file, newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_CACHE_COLUMNS, delimiter="\t")
        writer.writeheader()
        for info in taxonomy.values():
            writer.writerow(
                {
                    "genome_id": info.genome_id,
                    "species": info.species or "",
                    "genus": info.genus or "",
                    "family": info.family or "",
                    "order": info.order or "",
                    "class_name": info.class_name or "",
                    "phylum": info.phylum or "",
                    "organism": info.organism or "",
                    "tax_id": info.tax_id or "",
                }
            )
    logger.info("Saved %d entries to taxonomy cache %s", len(taxonomy), cache_file)


def _append_taxonomy_cache_row(cache_file: Path, info: TaxonomyInfo) -> None:
    """Append one genome's taxonomy to the TSV cache, writing the header only when the file is new.

    Only a successfully resolved lookup reaches here; a miss or a failed API call is not written,
    so a later run retries it instead of caching the gap forever.
    """
    cache_path = Path(cache_file)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    is_new = not cache_path.exists() or cache_path.stat().st_size == 0
    with open(cache_path, "a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_CACHE_COLUMNS, delimiter="\t")
        if is_new:
            writer.writeheader()
        writer.writerow(
            {
                "genome_id": info.genome_id,
                "species": info.species or "",
                "genus": info.genus or "",
                "family": info.family or "",
                "order": info.order or "",
                "class_name": info.class_name or "",
                "phylum": info.phylum or "",
                "organism": info.organism or "",
                "tax_id": info.tax_id or "",
            }
        )


def enrich_genomes_with_taxonomy(
    genome_ids: List[str],
    cache_file: Optional[Path] = None,
) -> Dict[str, TaxonomyInfo]:
    """Map genome accessions to taxonomy via the GTDB API.

    Loads cached results when available. Each newly resolved genome is appended to the cache
    file right after its lookup succeeds, so an interrupted run keeps every row it already
    fetched; a genome with no taxonomy found, or whose lookup failed, is not written to the
    cache and is retried on the next run (though this call still returns an empty TaxonomyInfo
    for it, so the caller always gets one entry per requested genome).
    """
    taxonomy: Dict[str, TaxonomyInfo] = {}

    if cache_file:
        taxonomy = load_taxonomy_cache(cache_file)

    missing = [gid for gid in genome_ids if gid not in taxonomy]
    if missing:
        logger.info(
            "Enriching %d genome(s) with taxonomy (%d cached)",
            len(missing),
            len(taxonomy),
        )
        for accession in missing:
            try:
                info = _lookup_genome_taxonomy_gtdb(accession)
                if info:
                    taxonomy[accession] = info
                    if cache_file:
                        _append_taxonomy_cache_row(cache_file, info)
                else:
                    logger.warning("No taxonomy found for genome '%s'", accession)
                    taxonomy[accession] = TaxonomyInfo(genome_id=accession)
            except DataAccessError:
                logger.warning("Failed to retrieve taxonomy for '%s'", accession)
                taxonomy[accession] = TaxonomyInfo(genome_id=accession)

    return taxonomy


_ANNOTATED_COLUMNS = ["sample", "genome", "containment", "species", "genus", "family"]


def annotate_containment_with_taxonomy(
    containment_df: pd.DataFrame,
    taxonomy: Dict[str, TaxonomyInfo],
) -> pd.DataFrame:
    """Add taxonomy columns to a parsed containment DataFrame.

    Converts the wide-format containment table (samples as rows, genomes as columns) into a
    long-format DataFrame with columns ``sample, genome, containment, species, genus, family``.

    Melts the genome columns instead of iterating row by row, and keeps only the sample/genome
    pairs with a positive containment value; a zero (or missing) containment cell is dropped
    rather than kept as a zero-valued row. This is a behaviour change from the previous
    row-by-row implementation, which kept every cell including zeros. Rows keep that
    implementation's order: sample by sample, and within a sample in genome column order.
    """
    from prepare.core.utils import get_genome_columns

    genome_cols = get_genome_columns(containment_df)

    melted = containment_df[genome_cols].melt(ignore_index=False, var_name="genome", value_name="containment")
    # melt is genome-major; a stable sort on the source row position restores sample-major order
    # while keeping genome column order within each sample.
    melted["_row"] = list(range(len(containment_df))) * len(genome_cols)
    melted = melted[melted["containment"] > 0].sort_values("_row", kind="stable")
    if melted.empty:
        return pd.DataFrame(columns=_ANNOTATED_COLUMNS)

    melted["sample"] = [name if isinstance(name, str) else str(name) for name in melted.index]
    melted["containment"] = melted["containment"].astype(float)

    infos = melted["genome"].map(lambda genome: taxonomy.get(genome, TaxonomyInfo(genome_id=genome)))
    melted["species"] = infos.map(lambda info: info.species)
    melted["genus"] = infos.map(lambda info: info.genus)
    melted["family"] = infos.map(lambda info: info.family)

    return melted[_ANNOTATED_COLUMNS].reset_index(drop=True)


def filter_by_taxonomy(
    annotated_df: pd.DataFrame,
    family: Optional[str] = None,
    genus: Optional[str] = None,
    species: Optional[str] = None,
    min_containment: float = 0.0,
) -> pd.DataFrame:
    """Filter an annotated containment DataFrame by taxonomy criteria."""
    df = annotated_df.copy()

    if family:
        df = df[df["family"].str.lower() == family.lower()]
    if genus:
        df = df[df["genus"].str.lower() == genus.lower()]
    if species:
        df = df[df["species"].str.lower() == species.lower()]
    if min_containment > 0:
        df = df[df["containment"] >= min_containment]

    return df


def summarize_by_taxonomy(
    annotated_df: pd.DataFrame,
    level: str = "family",
) -> pd.DataFrame:
    """Aggregate containment by taxonomic level.

    For each sample, computes the maximum containment per taxonomic group
    and returns a pivot table with samples as rows and taxonomy groups as
    columns.
    """
    valid_levels = {"family", "genus", "species"}
    if level not in valid_levels:
        raise ValueError(f"level must be one of {valid_levels}, got '{level}'")

    df = annotated_df.dropna(subset=[level])
    if df.empty:
        return pd.DataFrame()

    grouped = df.groupby(["sample", level])["containment"].max().reset_index()
    pivot = grouped.pivot_table(index="sample", columns=level, values="containment", fill_value=0.0)
    return pivot
