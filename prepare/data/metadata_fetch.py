"""Entrez download of SRA run metadata: single and batched efetch calls, pacing and credentials.

``prepare.data.metadata`` imports and re-exports these names, so existing imports keep working.
The orchestration (which accessions to fetch, the batch loop and the parsing of the XML) stays there.
"""

import contextlib
import copy
import http.client
import logging
import threading
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple
from urllib.error import HTTPError, URLError

from Bio import Entrez

from prepare.core import settings
from prepare.data.file_io import write_text_atomic
from prepare.data.ncbi_pacing import RequestPacer
from prepare.utils.progress import item_level

# The same logger as prepare.data.metadata: the code here ran there before the split.
logger = logging.getLogger("prepare.data.metadata")

# Maximum number of retries for failed downloads
MAX_RETRIES = 3

# Maximum accessions per single Entrez.efetch call; NCBI's own guidance caps URL-based
# id lists well below this, so this is a defensive ceiling rather than a tuned optimum.
MAX_BATCH_SIZE = 500

# One RequestPacer per rate (with/without API key), shared across batches; built under _PACERS_LOCK.
_PACERS: Dict[bool, RequestPacer] = {}
_PACERS_LOCK = threading.Lock()
# Serializes setting Bio.Entrez's email/api_key and the efetch call that uses them; see _entrez_credentials.
_ENTREZ_LOCK = threading.Lock()


def _pace_requests(api_key: Optional[str]) -> None:
    """Wait out NCBI's minimum gap since the previous request at this key's rate."""
    with _PACERS_LOCK:
        pacer = _PACERS.get(key := bool(api_key)) or _PACERS.setdefault(key, RequestPacer(api_key))
    pacer.wait()


@contextlib.contextmanager
def _entrez_credentials(email: str, api_key: Optional[str]) -> Iterator[None]:
    """Hold Bio.Entrez's email/api_key module attributes at the given values for one call.

    Biopython keeps ``Entrez.email``/``Entrez.api_key`` as plain module attributes, so two
    threads fetching with different credentials (different projects, or a library host running
    more than one download at once) race on them: one thread's efetch call can pick up the
    other's email or key. This holds _ENTREZ_LOCK for the duration, so only one thread has the
    attributes set to its own values while it calls efetch, and restores whatever was set
    before on exit (via ``finally``, so an exception from the call does not leave the wrong
    credentials in place for the next caller).
    """
    with _ENTREZ_LOCK:
        previous_email = Entrez.email
        previous_api_key = Entrez.api_key
        Entrez.email = email  # type: ignore[assignment]
        Entrez.api_key = api_key  # type: ignore[assignment]
        try:
            yield
        finally:
            Entrez.email = previous_email
            Entrez.api_key = previous_api_key


def _write_metadata_file(metadata_path: Path, accession: str, content: str) -> Path:
    """Write one accession's metadata XML atomically (``write_text_atomic``)."""
    return write_text_atomic(metadata_path / f"{accession}_metadata.xml", content)


def _split_efetch_packages(xml_text: str, wanted: Set[str]) -> Dict[str, str]:
    """Split a multi-accession efetch response into one XML document per wanted run accession.

    NCBI's batched efetch response is one ``EXPERIMENT_PACKAGE_SET`` holding one
    ``EXPERIMENT_PACKAGE`` per experiment; a package's ``RUN_SET`` can carry more than one
    ``RUN`` when several requested accessions share an experiment. For each wanted run, this
    returns a deep copy of its package whose ``RUN_SET`` keeps only that run, wrapped back in
    an ``EXPERIMENT_PACKAGE_SET`` with an XML declaration, so the result parses exactly like
    today's single-accession response (``parse_metadata_xml`` / ``_extract_metadata_fields``
    read ``.//RUN`` and its attributes unchanged).

    Returns ``{accession: xml_string}`` for the wanted accessions actually found; accessions
    absent from the response are simply missing from the result.
    """
    # xml.etree (expat), not prepare.utils.xml's lxml parser: ElementTree never resolves an
    # external entity or fetches a DTD, so XXE does not apply, and expat 2.4.1 and later rejects
    # entity amplification such as billion laughs. This is NCBI's own small efetch response.
    root = ET.fromstring(xml_text)
    found: Dict[str, str] = {}

    for package in root.findall("EXPERIMENT_PACKAGE"):
        for run_element in package.findall(".//RUN"):
            accession = run_element.get("accession")
            if not accession or accession not in wanted or accession in found:
                continue

            package_copy = copy.deepcopy(package)
            run_set = package_copy.find("RUN_SET")
            if run_set is not None:
                for run in list(run_set.findall("RUN")):
                    if run.get("accession") != accession:
                        run_set.remove(run)

            package_set = ET.Element("EXPERIMENT_PACKAGE_SET")
            package_set.append(package_copy)
            xml_bytes = ET.tostring(package_set, encoding="UTF-8", xml_declaration=True)
            found[accession] = xml_bytes.decode("utf-8")

    return found


def _download_single_metadata(
    accession: str, metadata_path: Path, entrez_email: str, api_key: Optional[str] = None
) -> Tuple[bool, Any]:
    """Download metadata for one accession; the file keeps only the requested run.

    Used directly and as the fallback when a batch response is rejected. A body that is not XML,
    or that holds no package for the accession, is reported as a failure and nothing is written.

    Returns:
        Tuple of (success, path or error message)
    """
    last_error_message = f"Failed after {MAX_RETRIES} attempts"

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            _pace_requests(api_key)
            logger.log(item_level(settings.active().progress_every), f"Downloading metadata for {accession}")
            with _entrez_credentials(entrez_email, api_key):
                handle = Entrez.efetch(db="sra", id=accession, retmode="xml")
                try:
                    metadata = handle.read().decode()
                finally:
                    handle.close()
            try:
                package = _split_efetch_packages(metadata, {accession}).get(accession)
            except ET.ParseError:
                return False, f"NCBI returned a non-XML body for {accession}"
            if package is None:
                return False, f"NCBI returned no package for {accession}"
            return True, _write_metadata_file(metadata_path, accession, package)

        except HTTPError as e:
            if e.code in (400, 404):
                return False, f"HTTP {e.code}: not found at NCBI"
            transient = "network: " if e.code == 429 or e.code >= 500 else ""
            last_error_message = f"{transient}HTTP {e.code} after {MAX_RETRIES} attempts"
            logger.warning(f"Error downloading {accession}, retrying ({attempt}/{MAX_RETRIES}): {e}")
            time.sleep(2**attempt)

        # A truncated NCBI reply raises http.client.IncompleteRead, an HTTPException, not an OSError.
        except (URLError, OSError, http.client.HTTPException) as e:
            last_error_message = f"network: {e}"
            logger.warning(f"Error downloading {accession}, retrying ({attempt}/{MAX_RETRIES}): {e}")
            time.sleep(2**attempt)

    return False, last_error_message


def _download_accessions_individually(
    batch: List[str], metadata_path: Path, email: str, api_key: Optional[str]
) -> Tuple[Dict[str, Path], Dict[str, str]]:
    """Fetch each accession in ``batch`` with its own ``efetch`` call.

    Used when a batched request comes back rejected outright (NCBI does not recognize one of
    the accessions in it), so each accession is re-tried on its own to isolate the bad one.
    """
    successes: Dict[str, Path] = {}
    failures: Dict[str, str] = {}
    for accession in batch:
        success, result = _download_single_metadata(accession, metadata_path, email, api_key)
        if success:
            successes[accession] = result
        else:
            failures[accession] = result
    return successes, failures


def _download_batch_metadata(
    batch: List[str], metadata_path: Path, email: str, api_key: Optional[str]
) -> Tuple[Dict[str, Path], Dict[str, str]]:
    """Fetch one batch of accessions with a single ``efetch`` call and split the response.

    Args:
        batch: Accessions to fetch together, joined into one comma-separated ``id`` parameter.
        metadata_path: Folder to write each accession's split-out XML file into.
        email: Email for NCBI API, forwarded to the single-accession fallback.
        api_key: Optional NCBI API key, forwarded to the pacing helper.

    Returns:
        Tuple of (``{accession: path}`` for successes, ``{accession: reason}`` for failures).
    """
    wanted = set(batch)
    id_string = ",".join(batch)
    last_failures = {accession: f"Failed after {MAX_RETRIES} attempts" for accession in batch}

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            _pace_requests(api_key)
            logger.log(
                item_level(settings.active().progress_every), f"Downloading metadata for {len(batch)} accession(s)"
            )
            with _entrez_credentials(email, api_key):
                handle = Entrez.efetch(db="sra", id=id_string, retmode="xml")
                try:
                    xml_text = handle.read().decode()
                finally:
                    handle.close()

            packages = _split_efetch_packages(xml_text, wanted)
            successes = {
                accession: _write_metadata_file(metadata_path, accession, packages[accession])
                for accession in batch
                if accession in packages
            }
            failures = {accession: "not in the NCBI response" for accession in batch if accession not in packages}
            return successes, failures

        except HTTPError as e:
            if e.code in (400, 404):
                if len(batch) > 1:
                    return _download_accessions_individually(batch, metadata_path, email, api_key)
                return {}, {batch[0]: f"HTTP {e.code}: not found at NCBI"}
            if e.code == 429 or 500 <= e.code < 600:
                last_failures = {acc: f"network: HTTP {e.code} after {MAX_RETRIES} attempts" for acc in batch}
                logger.warning(f"Error fetching batch, retrying ({attempt}/{MAX_RETRIES}): {e}")
                time.sleep(2**attempt)
                continue
            return {}, {accession: f"HTTP {e.code}: {e.reason}" for accession in batch}

        # A truncated NCBI reply raises http.client.IncompleteRead, an HTTPException, not an OSError.
        except (URLError, OSError, http.client.HTTPException) as e:
            last_failures = {accession: f"network: {e}" for accession in batch}
            logger.warning(f"Error fetching batch, retrying ({attempt}/{MAX_RETRIES}): {e}")
            time.sleep(2**attempt)

    return {}, last_failures
