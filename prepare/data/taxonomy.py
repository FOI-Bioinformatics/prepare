"""
Taxonomy integration and validation for PREPARE.

This module provides functions for integrating with NCBI taxonomy,
validating taxonomic assignments, and creating taxonomic summaries.
"""

import csv
import logging
import requests
import time
import xml.etree.ElementTree as ET
import pandas as pd
from typing import Dict, List, Optional, Union
from pathlib import Path
import re

from prepare.core.exceptions import ProcessingError
from prepare.data.file_io import open_atomic, write_csv
from prepare.utils.http import retrying_session

logger = logging.getLogger(__name__)

# Columns of the enrich_taxonomy map and the lineage rank each one holds.
MAP_RANK_COLUMNS = {
    "species": "species",
    "genus": "genus",
    "family": "family",
    "order": "order",
    "class_name": "class",
    "phylum": "phylum",
}

RETRY_TOTAL = 3
RETRY_BACKOFF_FACTOR = 0.5
RETRY_STATUS_FORCELIST = [429, 500, 502, 503, 504]

# Columns of the incremental species-validation cache, in write order.
_VALIDATION_CACHE_COLUMNS = [
    "original_name",
    "cleaned_name",
    "validated_name",
    "tax_id",
    "rank",
    "lineage",
    "is_valid",
    "confidence",
]


def _build_retrying_session() -> requests.Session:
    """Build a Session that retries a 429/5xx NCBI response with backoff.

    Each NCBITaxonomyClient owns one of these (it already keeps per-instance rate-limit
    state), reused across every esearch/efetch call that client makes. When every retry fails the
    last response is returned, so raise_for_status() reports the HTTP error rather than a RetryError.

    A thin wrapper (re-exported under its historical name) over
    ``prepare.utils.http.retrying_session``, which every other HTTP-polling client shares.
    """
    return retrying_session(
        total=RETRY_TOTAL,
        backoff_factor=RETRY_BACKOFF_FACTOR,
        status_forcelist=RETRY_STATUS_FORCELIST,
        allowed_methods=["GET"],
    )


class NCBITaxonomyClient:
    """Client for interacting with NCBI Taxonomy database."""

    def __init__(self, email: str, api_key: Optional[str] = None):
        """
        Initialize NCBI client.

        Args:
            email: Email address for NCBI API access
            api_key: Optional API key for increased rate limits
        """
        self.email = email
        self.api_key = api_key
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.last_request_time: float = 0
        self.request_delay = 1.0 / 3 if api_key else 1.0 / 10  # Rate limiting
        self.session = _build_retrying_session()

    def close(self) -> None:
        """Close the client's HTTP session and its pooled connections."""
        self.session.close()

    def __enter__(self) -> "NCBITaxonomyClient":
        """Return the client itself; the session is closed when the block exits."""
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        """Close the session."""
        self.close()

    def _make_request(self, url: str, params: Dict[str, str]) -> str:
        """Make rate-limited request to NCBI API."""
        # Rate limiting
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        if time_since_last < self.request_delay:
            time.sleep(self.request_delay - time_since_last)

        # Add required parameters
        params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key

        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            self.last_request_time = time.time()
            return response.text
        except requests.RequestException as e:
            logger.error(f"NCBI API request failed: {e}")
            raise ProcessingError(f"Failed to query NCBI: {e}")

    def search_taxonomy(self, query: str) -> List[Dict[str, str]]:
        """
        Search for taxonomic entries by name.

        Args:
            query: Scientific name or partial name to search

        Returns:
            List of matching taxonomy entries
        """
        search_url = f"{self.base_url}esearch.fcgi"
        params = {"db": "taxonomy", "term": query, "retmax": "10", "retmode": "xml"}

        response = self._make_request(search_url, params)

        # Parse XML response. xml.etree (expat), not prepare.utils.xml's lxml parser:
        # ElementTree never resolves an external entity or fetches a DTD, so XXE does not apply,
        # and expat 2.4.1 and later rejects entity amplification such as billion laughs. This is
        # NCBI's own small API response.
        root = ET.fromstring(response)
        tax_ids = [id_elem.text for id_elem in root.findall(".//Id") if id_elem.text is not None]

        if not tax_ids:
            return []

        # Get details for found IDs
        return self.get_taxonomy_details(tax_ids[:5])  # Limit to first 5 results

    @staticmethod
    def _element_text(parent, tag: str) -> str:
        """Return the text of a child element, or '' when missing."""
        el = parent.find(tag)
        return (el.text or "") if el is not None else ""

    @staticmethod
    def _parse_lineage(taxon) -> str:
        """Build a 'rank:name;...' lineage string from a <Taxon>'s LineageEx."""
        lineage = []
        lineage_elem = taxon.find(".//LineageEx")
        if lineage_elem is not None:
            for taxon_elem in lineage_elem.findall("Taxon"):
                name = taxon_elem.find("ScientificName")
                rank_elem = taxon_elem.find("Rank")
                if name is not None and rank_elem is not None:
                    lineage.append(f"{rank_elem.text}:{name.text}")
        return ";".join(lineage)

    @classmethod
    def _parse_taxon_element(cls, taxon) -> Dict[str, str]:
        """Parse a single <Taxon> element into a taxonomy detail dict."""
        return {
            "tax_id": cls._element_text(taxon, "TaxId"),
            "scientific_name": cls._element_text(taxon, "ScientificName"),
            "rank": cls._element_text(taxon, "Rank"),
            "lineage": cls._parse_lineage(taxon),
        }

    def get_taxonomy_details(self, tax_ids: List[str]) -> List[Dict[str, str]]:
        """
        Get detailed taxonomy information for given taxonomy IDs.

        Args:
            tax_ids: List of NCBI taxonomy IDs

        Returns:
            List of taxonomy detail dictionaries
        """
        if not tax_ids:
            return []

        fetch_url = f"{self.base_url}efetch.fcgi"
        params = {"db": "taxonomy", "id": ",".join(tax_ids), "retmode": "xml"}

        response = self._make_request(fetch_url, params)

        # xml.etree (expat): see the comment in search_taxonomy above.
        root = ET.fromstring(response)
        return [self._parse_taxon_element(taxon) for taxon in root.findall("./Taxon")]

    def validate_species_name(self, species_name: str) -> Dict[str, Union[str, bool]]:
        """
        Validate a species name against NCBI taxonomy.

        Args:
            species_name: Species name to validate

        Returns:
            Dictionary with validation results
        """
        try:
            # Clean the species name
            cleaned_name = self._clean_species_name(species_name)

            # Search for exact match first
            results = self.search_taxonomy(f'"{cleaned_name}"[Scientific Name]')

            if not results:
                # Try broader search
                results = self.search_taxonomy(cleaned_name)

            if results:
                best_match = results[0]  # Take first/best result
                return {
                    "original_name": species_name,
                    "cleaned_name": cleaned_name,
                    "validated_name": best_match["scientific_name"],
                    "tax_id": best_match["tax_id"],
                    "rank": best_match["rank"],
                    "lineage": best_match["lineage"],
                    "is_valid": True,
                    "confidence": (
                        "high" if cleaned_name.lower() == best_match["scientific_name"].lower() else "medium"
                    ),
                }
            else:
                return {
                    "original_name": species_name,
                    "cleaned_name": cleaned_name,
                    "validated_name": "",
                    "tax_id": "",
                    "rank": "",
                    "lineage": "",
                    "is_valid": False,
                    "confidence": "low",
                }

        except Exception as e:
            logger.warning(f"Failed to validate species name '{species_name}': {e}")
            return {
                "original_name": species_name,
                "cleaned_name": species_name,
                "validated_name": "",
                "tax_id": "",
                "rank": "",
                "lineage": "",
                "is_valid": False,
                "confidence": "error",
            }

    def _clean_species_name(self, name: str) -> str:
        """Clean and standardize species name for taxonomy lookup."""
        # Remove common prefixes/suffixes
        cleaned = re.sub(r"^(candidate|unclassified|uncultured)\s+", "", name, flags=re.IGNORECASE)
        cleaned = re.sub(r"\s+(sp\.|strain|isolate|clone)\s+.*$", "", cleaned, flags=re.IGNORECASE)

        # Remove brackets and their contents
        cleaned = re.sub(r"\[.*?\]", "", cleaned)
        cleaned = re.sub(r"\(.*?\)", "", cleaned)

        # Standardize whitespace
        cleaned = " ".join(cleaned.split())

        return cleaned.strip()


def _append_validation_cache_row(
    cache_file: Union[str, Path],
    result: Dict[str, Union[str, bool]],
    reset: bool = False,
) -> None:
    """Append one species validation to the CSV cache, writing the header when the file is new.

    ``reset`` truncates the file first, for the one write that follows a cache file that failed
    to load: a corrupt cache is discarded rather than appended to.
    """
    cache_path = Path(cache_file)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    row = {column: result.get(column, "") for column in _VALIDATION_CACHE_COLUMNS}
    if reset:
        # A fresh file replaces the corrupt one in a single rename.
        with open_atomic(cache_path, newline="") as f:
            writer = csv.DictWriter(f, fieldnames=_VALIDATION_CACHE_COLUMNS)
            writer.writeheader()
            writer.writerow(row)
        return
    is_new = not cache_path.exists() or cache_path.stat().st_size == 0
    with open(cache_path, "a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_VALIDATION_CACHE_COLUMNS)
        if is_new:
            writer.writeheader()
        writer.writerow(row)


def validate_taxonomic_assignments(
    species_list: List[str],
    email: str,
    api_key: Optional[str] = None,
    output_file: Optional[Union[str, Path]] = None,
    cache_file: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Validate a list of species names against NCBI taxonomy.

    Args:
        species_list: List of species names to validate
        email: Email for NCBI API access
        api_key: Optional API key
        output_file: File to save validation results
        cache_file: File to cache results for reuse

    Returns:
        DataFrame with validation results
    """
    try:
        # Load cache if available
        cached_results = {}
        cache_needs_reset = False
        if cache_file and Path(cache_file).exists():
            try:
                cache_df = pd.read_csv(cache_file, encoding="utf-8")
                cached_results = dict(zip(cache_df["original_name"], cache_df.to_dict("records")))
                logger.info(f"Loaded {len(cached_results)} cached taxonomy validations")
            except Exception as e:
                logger.warning(f"Failed to load cache file: {e}")
                # The file is unreadable; the first successful lookup below replaces it
                # instead of appending after its corrupt content.
                cache_needs_reset = True

        # One client, and so one HTTP session, for the whole list; closed when done.
        client = NCBITaxonomyClient(email, api_key)
        try:
            # Validate species
            results = []
            unique_species = list(set(species_list))

            logger.info(f"Validating {len(unique_species)} unique species names")

            for i, species in enumerate(unique_species, 1):
                if species in cached_results:
                    results.append(cached_results[species])
                    logger.debug(f"Using cached result for {species}")
                else:
                    logger.info(f"Validating {i}/{len(unique_species)}: {species}")
                    result = client.validate_species_name(species)
                    results.append(result)

                    # A confidence of "error" means the lookup itself failed (e.g. a network
                    # error); leave it out of the cache so a later run retries it, instead of
                    # caching the failure forever.
                    if cache_file and result.get("confidence") != "error":
                        _append_validation_cache_row(cache_file, result, reset=cache_needs_reset)
                        cache_needs_reset = False

                    # Rate limiting for API calls
                    if i % 10 == 0:
                        logger.info(f"Progress: {i}/{len(unique_species)} species validated")
        finally:
            client.close()

        # Create results DataFrame
        if results:
            results_df = pd.DataFrame(results)
        else:
            # Create empty DataFrame with expected columns
            results_df = pd.DataFrame(
                columns=[
                    "original_name",
                    "validated_name",
                    "is_valid",
                    "tax_id",
                    "rank",
                    "lineage",
                    "confidence",
                ]
            )

        # Save results
        if output_file:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            write_csv(results_df, output_path, index=False)
            logger.info(f"Taxonomy validation results saved to {output_path}")

        # The cache (if requested) was already updated incrementally above, one row per
        # successful lookup, so an interrupted run keeps whatever it fetched.

        # Log summary
        valid_count = results_df["is_valid"].sum() if not results_df.empty else 0
        logger.info(f"Taxonomy validation complete: {valid_count}/{len(results_df)} species validated")  # noqa: E501

        return results_df

    except Exception as e:
        raise ProcessingError(f"Failed to validate taxonomic assignments: {e}")


def _build_taxonomy_lineage_map(taxonomy_data: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    """Map a species name or genome id to {rank: name}, from either taxonomy table PREPARE writes."""
    columns = set(taxonomy_data.columns)
    if "genome_id" in columns:
        return _lineages_from_taxonomy_map(taxonomy_data)
    if {"original_name", "is_valid", "lineage"}.issubset(columns):
        return _lineages_from_validation_table(taxonomy_data)
    raise ProcessingError(
        "Taxonomy file must come from validate_taxonomy (columns original_name, is_valid, lineage) "
        "or enrich_taxonomy (column genome_id with rank columns)"
    )


def _lineages_from_validation_table(taxonomy_data: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    taxonomy_dict: Dict[str, Dict[str, str]] = {}
    for _, row in taxonomy_data.iterrows():
        if not (row["is_valid"] and row["lineage"]):
            continue
        lineage_dict = {}
        for part in str(row["lineage"]).split(";"):
            if ":" in part:
                rank, name = part.split(":", 1)
                lineage_dict[rank.lower()] = name
        taxonomy_dict[row["original_name"]] = lineage_dict
    return taxonomy_dict


def _lineages_from_taxonomy_map(taxonomy_data: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    lineages: Dict[str, Dict[str, str]] = {}
    for _, row in taxonomy_data.iterrows():
        ranks = {}
        for column, rank in MAP_RANK_COLUMNS.items():
            value = row.get(column)
            if isinstance(value, str) and value.strip():
                ranks[rank] = value.strip()
        lineages[str(row["genome_id"])] = ranks
    return lineages


def _summarize_sample_taxonomy(sample_abundances, taxonomy_dict, level, min_abundance) -> Dict[str, float]:
    """Aggregate one sample's abundances by taxon at `level`, bucketing unknowns."""
    sample_summary: Dict[str, float] = {}
    for species, abundance in sample_abundances.items():
        if abundance < min_abundance:
            continue
        lineage = taxonomy_dict.get(species, {})
        taxon_name = lineage.get(level, f"Unclassified_{level}")
        sample_summary[taxon_name] = sample_summary.get(taxon_name, 0) + abundance
    return sample_summary


def create_taxonomic_summary(
    abundance_data: pd.DataFrame,
    taxonomy_data: pd.DataFrame,
    level: str = "genus",
    min_abundance: float = 0.001,
) -> pd.DataFrame:
    """
    Create taxonomic summary at specified level.

    Args:
        abundance_data: Sample x Species abundance matrix
        taxonomy_data: Taxonomy information with lineage
        level: Taxonomic level for summary (genus, family, order, etc.)
        min_abundance: Minimum abundance threshold

    Returns:
        Sample x Taxon summary matrix
    """
    try:
        taxonomy_dict = _build_taxonomy_lineage_map(taxonomy_data)

        summary_data = [
            _summarize_sample_taxonomy(abundance_data.loc[sample_name], taxonomy_dict, level, min_abundance)
            for sample_name in abundance_data.index
        ]

        summary_df = pd.DataFrame(summary_data, index=abundance_data.index).fillna(0)

        logger.info(f"Created taxonomic summary at {level} level: {summary_df.shape}")
        return summary_df

    except ProcessingError:
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to create taxonomic summary: {e}")


def analyze_taxonomic_composition(
    abundance_data: pd.DataFrame,
    taxonomy_data: pd.DataFrame,
    levels: List[str] = ["phylum", "class", "order", "family", "genus"],
    output_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Analyze taxonomic composition at multiple levels.

    Args:
        abundance_data: Sample x Species abundance matrix
        taxonomy_data: Taxonomy information
        levels: List of taxonomic levels to analyze
        output_dir: Directory to save results

    Returns:
        Dictionary of summary matrices for each level
    """
    try:
        results = {}

        for level in levels:
            logger.info(f"Creating taxonomic summary at {level} level")
            summary_df = create_taxonomic_summary(abundance_data, taxonomy_data, level)
            results[level] = summary_df

            if output_dir:
                output_path = Path(output_dir) / f"taxonomy_summary_{level}.csv"
                output_path.parent.mkdir(parents=True, exist_ok=True)
                write_csv(summary_df, output_path)
                logger.info(f"Saved {level} summary to {output_path}")

        return results

    except ProcessingError:
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to analyze taxonomic composition: {e}")
