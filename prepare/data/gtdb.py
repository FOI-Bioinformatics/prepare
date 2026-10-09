"""GTDB API client for genome taxonomy database queries."""

import logging
import threading
from typing import Dict, List, Optional

import requests

from prepare.core.exceptions import DataAccessError
from prepare.utils.http import retrying_session

logger = logging.getLogger(__name__)

GTDB_API_BASE = "https://gtdb-api.ecogenomic.org"
REQUEST_TIMEOUT = 30
RETRY_TOTAL = 3
RETRY_BACKOFF_FACTOR = 0.5
RETRY_STATUS_FORCELIST = [429, 500, 502, 503, 504]

_session: Optional[requests.Session] = None
_SESSION_LOCK = threading.Lock()


def get_session() -> requests.Session:
    """Return the process-wide GTDB HTTP session, building it on first use.

    Every genome, species and taxon lookup shares this one session, so TCP/TLS connections are
    pooled across a whole run instead of rebuilt per request. A 429 or 5xx GTDB response is
    retried up to RETRY_TOTAL times with exponential backoff; if every attempt fails, the last
    response is returned (not raised), so status checks and raise_for_status() see it as before.

    Double-checked locking: the cheap unlocked check on the fast path skips the lock once the
    session exists; two threads racing on the first call both wait on _SESSION_LOCK and only
    the first one to acquire it builds the session, so every caller ends up sharing one.
    """
    global _session
    if _session is None:
        with _SESSION_LOCK:
            if _session is None:
                _session = retrying_session(
                    total=RETRY_TOTAL,
                    backoff_factor=RETRY_BACKOFF_FACTOR,
                    status_forcelist=RETRY_STATUS_FORCELIST,
                    allowed_methods=["GET"],
                )
    return _session


def _api_error_message(kind: str, name: str, e: requests.exceptions.RequestException) -> str:
    """Build the DataAccessError message for a failed GTDB request.

    A 400 response usually means the name is not a current GTDB name: GTDB periodically
    reclassifies genera (for example many former Lactobacillus species now sit in
    Lacticaseibacillus, Limosilactobacillus or another genus), so the message suggests
    searching the current genus instead. Other HTTP and connection errors keep the plain
    message, since a status other than 400 is not evidence of a stale name.
    """
    message = f"GTDB API error searching {kind} '{name}': {e}"
    response = getattr(e, "response", None)
    if response is not None and getattr(response, "status_code", None) == 400:
        message += (
            f". '{name}' is not a current GTDB name (GTDB periodically reclassifies genera, "
            "for example former Lactobacillus species now placed in Lacticaseibacillus, "
            "Limosilactobacillus or another genus); try searching the current genus name instead."
        )
    return message


def search_species(species_name: str) -> List[Dict]:
    """Search GTDB for a species, return list of genome records with accessions."""
    url = f"{GTDB_API_BASE}/species/search/{requests.utils.quote(species_name)}"
    logger.debug("Searching GTDB species: %s", species_name)

    try:
        response = get_session().get(url, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as e:
        raise DataAccessError(_api_error_message("species", species_name, e))

    if not data:
        logger.warning("No results found for species '%s'", species_name)
        return []

    # The API may return a list or a single object with a genomes array
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        return data.get("genomes", [data])

    return []


def search_taxon(taxon_name: str, limit: int = 100) -> List[Dict]:
    """Search GTDB by any taxonomic name (genus, family, etc.)."""
    url = f"{GTDB_API_BASE}/taxon/search/{requests.utils.quote(taxon_name)}"
    params = {"limit": limit}
    logger.debug("Searching GTDB taxon: %s (limit=%d)", taxon_name, limit)

    try:
        response = get_session().get(url, params=params, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as e:
        raise DataAccessError(_api_error_message("taxon", taxon_name, e))

    if not data:
        logger.warning("No results found for taxon '%s'", taxon_name)
        return []

    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        if "matches" in data:
            return [{"name": str(name)} for name in data["matches"]]
        return data.get("results", [data])

    return []


def _record_accession(record: Dict) -> Optional[str]:
    """Pull the assembly accession from a GTDB record under any of its key spellings."""
    return record.get("accession") or record.get("gid") or record.get("ncbi_accession")


def _is_representative(record: Dict) -> bool:
    """Whether a GTDB record is flagged as a representative genome (any key spelling)."""
    return bool(record.get("isRep") or record.get("is_representative") or record.get("gtdb_species_rep"))


def _keep_accession(record: Dict, representative_only: bool) -> Optional[str]:
    """Return the record's accession if it should be kept under the rep filter, else None."""
    accession = _record_accession(record)
    if accession and (not representative_only or _is_representative(record)):
        return accession
    return None


def get_accessions_for_species(species_name: str, representative_only: bool = True) -> List[str]:
    """Get assembly accessions (GCF_*/GCA_*) for a species.

    If representative_only is True, return only the representative genome.
    """
    results = search_species(species_name)
    if not results:
        return []

    accessions = [acc for record in results if (acc := _keep_accession(record, representative_only))]

    # If representative_only but none flagged as representative, return the first accession
    if representative_only and not accessions:
        first = _record_accession(results[0])
        if first:
            accessions.append(first)

    return accessions


def _append_unique(accessions: List[str], accession: Optional[str]) -> None:
    """Append an accession to the list if present and not already there."""
    if accession and accession not in accessions:
        accessions.append(accession)


def _resolve_genus_species(
    species_name: str, genus_name: str, representative_only: bool, accessions: List[str]
) -> bool:
    """Fan out to the species endpoint for one taxon match.

    Appends any resolved accessions to ``accessions`` and returns whether the
    lookup succeeded (False after logging a warning on failure).
    """
    logger.debug("Resolving species %s for genus %s", species_name, genus_name)
    try:
        for species_accession in get_accessions_for_species(species_name, representative_only):
            _append_unique(accessions, species_accession)
        return True
    except DataAccessError as e:
        logger.warning("Skipping species %s: %s", species_name, e)
        return False


def get_accessions_for_genus(genus_name: str, representative_only: bool = True) -> List[str]:
    """Get representative accessions for all species in a genus.

    The live taxon endpoint returns taxon names only; each ``s__`` species name is
    resolved through the species endpoint. Older record-shaped responses (with an
    accession per record) are still handled.
    """
    taxon_results = search_taxon(genus_name)
    if not taxon_results:
        return []

    accessions: List[str] = []
    seen_species = set()
    attempted = 0
    failed = 0
    skipped_other_genus = 0
    genus_prefix = f"{genus_name} "
    for record in taxon_results:
        name = str(record.get("species") or record.get("name", ""))
        if name in seen_species:
            continue
        seen_species.add(name)

        if name.startswith("s__") and _record_accession(record) is None:
            species_name = name[3:]
            if not species_name.startswith(genus_prefix):
                # /taxon/search is a name-contains search, so a genus query also
                # matches species of other genera (e.g. "Bacillus" also matches
                # "Paenibacillus larvae"); keep only this genus's own species.
                skipped_other_genus += 1
                continue
            attempted += 1
            if not _resolve_genus_species(species_name, genus_name, representative_only, accessions):
                failed += 1
            continue

        _append_unique(accessions, _keep_accession(record, representative_only))

    if skipped_other_genus:
        logger.info(
            "Skipped %d taxon match(es) outside genus '%s'",
            skipped_other_genus,
            genus_name,
        )

    if attempted and failed == attempted:
        raise DataAccessError(f"All {attempted} species lookups for genus '{genus_name}' failed")

    if representative_only and not accessions:
        for record in taxon_results:
            _append_unique(accessions, _record_accession(record))

    return accessions
