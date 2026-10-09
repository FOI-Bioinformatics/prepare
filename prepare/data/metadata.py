"""
Metadata handling for PREPARE.

This module provides functions for downloading and processing metadata from NCBI.
"""

import contextlib
import copy
import http.client
import logging
import os
import threading
import xml.etree.ElementTree as ET
from pathlib import Path
import time
import pandas as pd
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple, Union

from Bio import Entrez
from lxml import etree
from urllib.error import HTTPError, URLError

from prepare.core import settings
from prepare.core.exceptions import DataAccessError, PrepareError, NetworkError
from prepare.core.validation import validate_folder
from prepare.data.file_io import ensure_directory, list_files, write_csv, write_text_atomic
from prepare.utils.progress import ProgressReporter, item_level
from prepare.utils.xml import parse_xml_file

logger = logging.getLogger(__name__)

# Maximum number of retries for failed downloads
MAX_RETRIES = 3

# Maximum accessions per single Entrez.efetch call; NCBI's own guidance caps URL-based
# id lists well below this, so this is a defensive ceiling rather than a tuned optimum.
MAX_BATCH_SIZE = 500

# Minimum seconds between successive NCBI requests, without and with an API key
# (NCBI allows roughly 3 requests/second without a key and 10/second with one).
_RATE_LIMIT_DELAY_NO_KEY = 0.34
_RATE_LIMIT_DELAY_WITH_KEY = 0.1

# Monotonic timestamp of the previous NCBI request, module-level so pacing holds across
# batches within one process.
_last_request_time = 0.0
# Held across the sleep in _pace_requests, so two threads racing on a request each wait
# out the full gap rather than both sleeping the same (too-short) remainder concurrently.
_PACE_LOCK = threading.Lock()

# Serializes the window between setting Bio.Entrez's email/api_key module attributes and
# making the efetch call that depends on them; see _entrez_credentials below.
_ENTREZ_LOCK = threading.Lock()


def _pace_requests(api_key: Optional[str]) -> None:
    """Sleep only long enough to respect NCBI's rate limit since the previous request.

    This replaces a flat per-accession sleep: batching already cuts the number of requests,
    so the only sleep left on the success path is the minimum gap NCBI expects between calls.
    """
    global _last_request_time
    delay = _RATE_LIMIT_DELAY_WITH_KEY if api_key else _RATE_LIMIT_DELAY_NO_KEY
    with _PACE_LOCK:
        elapsed = time.monotonic() - _last_request_time
        if elapsed < delay:
            time.sleep(delay - elapsed)
        _last_request_time = time.monotonic()


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


def _get_unique_accessions(matches_folder, threshold):
    """
    Extract unique accessions from match files.

    Args:
        matches_folder: Folder containing match files
        threshold: Minimum containment threshold

    Returns:
        Set of unique accessions
    """
    unique_accessions: set = set()
    matches_path = Path(matches_folder)

    # Find all CSV files in matches folder
    csv_files = list_files(matches_path, "*.csv")

    if not csv_files:
        logger.warning(f"No CSV files found in {matches_path}")
        return unique_accessions

    for csv_file in csv_files:
        # Read CSV file
        try:
            df = pd.read_csv(csv_file)

            # Determine accession column name based on format
            if "acc" in df.columns:
                accession_col = "acc"
                containment_col = "containment"
            elif "SRA accession" in df.columns:
                accession_col = "SRA accession"
                containment_col = "containment"
            else:
                logger.warning(f"Unknown file format: {csv_file}")
                continue

            # Filter by threshold
            if threshold > 0:
                df = df[df[containment_col] >= threshold]

            # Add accessions to set
            unique_accessions.update(df[accession_col].tolist())

        # pandas' ParserError and EmptyDataError are ValueErrors; KeyError is a missing
        # containment column; TypeError is a containment column holding text.
        except (OSError, ValueError, KeyError, TypeError) as e:
            logger.warning(f"Error reading {csv_file}: {e}")

    return unique_accessions


def _read_accessions_file(path: Union[str, Path]) -> List[str]:
    """Read non-empty, non-comment accession lines from a file, one accession per line."""
    return [ln.strip() for ln in Path(path).read_text().splitlines() if ln.strip() and not ln.strip().startswith("#")]


def _download_single_metadata(
    accession: str, metadata_path: Path, entrez_email: str, api_key: Optional[str] = None
) -> Tuple[bool, Any]:
    """
    Download metadata for a single accession.

    Used both directly (one accession requested on its own) and as the fallback path when a
    batch response is rejected outright and each of its accessions must be re-fetched alone.

    Args:
        accession: SRA accession
        metadata_path: Path to save metadata
        entrez_email: Email for NCBI API
        api_key: Optional NCBI API key, for the pacing delay and higher rate limits

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
            return True, _write_metadata_file(metadata_path, accession, metadata)

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


def download_metadata(
    email: str,
    matches_folder: Union[str, Path],
    metadata_folder: Union[str, Path],
    threshold: float = 0.0,
    dry_run: bool = False,
    accessions_file: Optional[Union[str, Path]] = None,
    api_key: Optional[str] = None,
    batch_size: int = 200,
) -> Dict[str, Path]:
    """
    Download metadata for SRA accessions found in match files, or from an explicit list.

    Accessions are fetched from NCBI in batches (one ``Entrez.efetch`` call per batch) rather
    than one request per accession, which is both faster and gentler on NCBI's rate limits.

    Args:
        email: Email address for NCBI API
        matches_folder: Folder containing match files
        metadata_folder: Folder to save metadata files
        threshold: Minimum containment threshold
        dry_run: If True, only count accessions without downloading
        accessions_file: When given, the wanted accessions come from this file's
            non-empty, non-comment lines and the matches folder is not read.
        api_key: Optional NCBI API key; raises the rate limit and is used for every request.
        batch_size: Accessions per ``Entrez.efetch`` call, from 1 to 500.

    Returns:
        Dictionary mapping accessions to metadata file paths

    Raises:
        ValueError: If batch_size is outside 1 to 500.
        DataAccessError: If the download fails
    """
    if not 1 <= batch_size <= MAX_BATCH_SIZE:
        raise ValueError(f"batch_size must be between 1 and {MAX_BATCH_SIZE}, got {batch_size}")

    try:
        metadata_path = ensure_directory(metadata_folder)

        if accessions_file:
            unique_accessions: set = set(_read_accessions_file(accessions_file))
        else:
            matches_path = validate_folder(matches_folder)
            unique_accessions = _get_unique_accessions(matches_path, threshold)

        total_accessions = len(unique_accessions)
        logger.info(f"Found {total_accessions} unique accessions")

        # Check which accessions need downloading, from one listing of the folder
        existing = set(os.listdir(metadata_path))
        accessions_to_download = [
            accession for accession in unique_accessions if f"{accession}_metadata.xml" not in existing
        ]

        to_download_count = len(accessions_to_download)
        logger.info(f"Need to download {to_download_count} accessions")

        if dry_run:
            logger.info("Dry run, not downloading metadata")
            return {}

        # Download metadata for each accession, in batches
        return _download_accessions_metadata(
            accessions_to_download, metadata_path, email, to_download_count, api_key=api_key, batch_size=batch_size
        )

    except NetworkError:
        raise
    except (OSError, ValueError, PrepareError) as e:
        raise DataAccessError(f"Error downloading metadata: {e}") from e


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


def _download_accessions_metadata(
    accessions_to_download: List[str],
    metadata_path: Path,
    email: str,
    total_count: int,
    api_key: Optional[str] = None,
    batch_size: int = 200,
) -> Dict[str, Path]:
    """
    Download metadata for multiple accessions, fetched in batches.

    Args:
        accessions_to_download: List of accessions to download
        metadata_path: Path to save metadata
        email: Email for NCBI API
        total_count: Total number of accessions to download
        api_key: Optional NCBI API key
        batch_size: Accessions per ``Entrez.efetch`` call

    Returns:
        Dictionary mapping accessions to metadata file paths, for successes only.
    """
    # Each efetch call sets Entrez.email/api_key itself under _ENTREZ_LOCK (_entrez_credentials).
    result_files: Dict[str, Path] = {}
    failures: Dict[str, str] = {}
    batches = [accessions_to_download[i : i + batch_size] for i in range(0, len(accessions_to_download), batch_size)]
    progress = ProgressReporter(
        "download_metadata", len(accessions_to_download), settings.active().progress_every, logger=logger
    )

    for batch in batches:
        batch_successes, batch_failures = _download_batch_metadata(batch, metadata_path, email, api_key)
        result_files.update(batch_successes)
        failures.update(batch_failures)
        progress.update_counts(len(batch_successes), len(batch_failures))
    progress.finish()

    logger.info(f"Fetched {len(batches)} batches: {len(result_files)} metadata files, {len(failures)} failures")
    for accession, reason in failures.items():
        logger.error(f"Failed to download {accession}: {reason}")
    if failures and not result_files and all(reason.startswith("network: ") for reason in failures.values()):
        raise NetworkError(f"NCBI could not be reached for any of {len(failures)} accession(s); try again later")
    return result_files


def _run_attr(run_element, attr_name):
    """Read an attribute off a ``<RUN>`` element, or ``None`` when absent or unset."""
    if run_element is None:
        return None
    return run_element.get(attr_name)


def _first_srafile(srafile_elements):
    """Return the ``<SRAFile>`` element whose ``semantic_name`` is ``run``, else the first one.

    NCBI's efetch XML lists every file associated with a run (the sequencing data plus any
    reference or index files); the actual run data is the entry marked ``semantic_name="run"``.
    """
    for element in srafile_elements:
        if element.get("semantic_name") == "run":
            return element
    return srafile_elements[0] if srafile_elements else None


def _first_child_tag(element):
    """Return the tag name of the first child of ``element``, or ``None``."""
    if element is None:
        return None
    children = list(element)
    return children[0].tag if children else None


# The fixed columns of the metadata table, in table order; sample attribute columns follow them.
METADATA_FIXED_COLUMNS = tuple(
    "Run_ID Run_Total_Spots Run_Total_Bases Run_Size Run_Download_Path Run_MD5 Run_Filename Run_Spot_Length"
    " Run_Reads Run_FTP Run_Aspera Run_Galaxy Project_ID Project_Title Project_Abstract Sample_ID Sample_External_ID"
    " Sample_Name Sample_Scientific_Name Sample_Title Experiment_ID Experiment_Title Experiment_Design"
    " Experiment_Library_Name Experiment_Library_Strategy Experiment_Library_Source Experiment_Library_Selection"
    " Experiment_Library_Layout Platform SRA_Normalized_URL".split()
)


def _run_spot_length(tree):
    """Sum the ``average`` read length reported for each read in ``Statistics/Read``.

    Returns ``None`` when no such statistics are present, so the caller can fall back to
    the older ``<RUN/spot_length>`` child-element read.
    """
    read_elements = tree.findall(".//RUN/Statistics/Read")
    averages = []
    for read_element in read_elements:
        average = read_element.get("average")
        if average is None:
            continue
        try:
            averages.append(float(average))
        except ValueError:
            continue
    if not averages:
        return None
    total = sum(averages)
    return str(int(total)) if total.is_integer() else str(total)


def _extract_metadata_fields(tree, xml_file):
    """
    Extract metadata fields from an XML tree.

    Args:
        tree: XML tree to extract data from
        xml_file: Path to XML file (for error reporting)

    Returns:
        Dictionary with extracted metadata
    """
    try:
        # Extract project information
        project_id = tree.findtext(".//STUDY/IDENTIFIERS/PRIMARY_ID")
        project_title = tree.findtext(".//STUDY/DESCRIPTOR/STUDY_TITLE")
        project_abstract = tree.findtext(".//STUDY/DESCRIPTOR/STUDY_ABSTRACT")

        # Extract sample information
        sample_id = tree.findtext(".//SAMPLE/IDENTIFIERS/PRIMARY_ID")
        sample_external_id = tree.findtext(".//SAMPLE/IDENTIFIERS/EXTERNAL_ID")
        sample_name = tree.findtext(".//SAMPLE/SAMPLE_NAME/TAXON_ID")
        sample_scientific_name = tree.findtext(".//SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME")
        sample_title = tree.findtext(".//SAMPLE/TITLE")

        # Extract run information.
        # NCBI's efetch XML carries spots, bases, size and md5 as attributes on <RUN> and
        # <SRAFile>, not as child elements; the old child-element reads are kept as a
        # fallback for XML that predates this (or comes from a different source).
        run_id = tree.findtext(".//RUN/IDENTIFIERS/PRIMARY_ID")
        run_element = tree.find(".//RUN")
        run_total_spots = _run_attr(run_element, "total_spots") or tree.findtext(".//RUN/Total_spots")
        run_total_bases = _run_attr(run_element, "total_bases") or tree.findtext(".//RUN/Total_bases")
        run_size = _run_attr(run_element, "size") or tree.findtext(".//RUN/size")
        run_download_path = tree.findtext(".//RUN/download_path")
        srafile_elements = tree.findall(".//RUN/SRAFiles/SRAFile")
        srafile_element = _first_srafile(srafile_elements)
        run_md5 = (srafile_element.get("md5") if srafile_element is not None else None) or tree.findtext(".//RUN/md5")
        run_filename = (srafile_element.get("filename") if srafile_element is not None else None) or tree.findtext(
            ".//RUN/filename"
        )
        run_spot_length = _run_spot_length(tree) or tree.findtext(".//RUN/spot_length")
        run_reads = tree.findtext(".//RUN/reads")
        run_ftp = tree.findtext(".//RUN/ftp")
        run_aspera = tree.findtext(".//RUN/aspera")
        run_galaxy = tree.findtext(".//RUN/galaxy")

        # Extract experiment information
        experiment_id = tree.findtext(".//EXPERIMENT/IDENTIFIERS/PRIMARY_ID")
        experiment_title = tree.findtext(".//EXPERIMENT/TITLE")
        experiment_design = tree.findtext(".//EXPERIMENT/DESIGN/DESIGN_DESCRIPTION")
        experiment_library_name = tree.findtext(".//EXPERIMENT//LIBRARY_DESCRIPTOR/LIBRARY_NAME")
        experiment_library_strategy = tree.findtext(".//EXPERIMENT//LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY")
        experiment_library_source = tree.findtext(".//EXPERIMENT//LIBRARY_DESCRIPTOR/LIBRARY_SOURCE")
        experiment_library_selection = tree.findtext(".//EXPERIMENT//LIBRARY_DESCRIPTOR/LIBRARY_SELECTION")
        experiment_library_layout = _first_child_tag(tree.find(".//LIBRARY_LAYOUT"))
        platform = _first_child_tag(tree.find(".//PLATFORM"))

        # Extract SRA URL
        sra_normalized_url = None
        if len(srafile_elements) > 1:
            sra_normalized_url = srafile_elements[1].get("url")

        # Create metadata record
        metadata_dict = {
            "Run_ID": run_id,
            "Run_Total_Spots": run_total_spots,
            "Run_Total_Bases": run_total_bases,
            "Run_Size": run_size,
            "Run_Download_Path": run_download_path,
            "Run_MD5": run_md5,
            "Run_Filename": run_filename,
            "Run_Spot_Length": run_spot_length,
            "Run_Reads": run_reads,
            "Run_FTP": run_ftp,
            "Run_Aspera": run_aspera,
            "Run_Galaxy": run_galaxy,
            "Project_ID": project_id,
            "Project_Title": project_title,
            "Project_Abstract": project_abstract,
            "Sample_ID": sample_id,
            "Sample_External_ID": sample_external_id,
            "Sample_Name": sample_name,
            "Sample_Scientific_Name": sample_scientific_name,
            "Sample_Title": sample_title,
            "Experiment_ID": experiment_id,
            "Experiment_Title": experiment_title,
            "Experiment_Design": experiment_design,
            "Experiment_Library_Name": experiment_library_name,
            "Experiment_Library_Strategy": experiment_library_strategy,
            "Experiment_Library_Source": experiment_library_source,
            "Experiment_Library_Selection": experiment_library_selection,
            "Experiment_Library_Layout": experiment_library_layout,
            "Platform": platform,
            "SRA_Normalized_URL": sra_normalized_url,
        }

        return metadata_dict

    except ValueError as e:
        logger.error(f"Error extracting fields from {xml_file}: {e}")
        return {}


def _extract_sample_attributes(tree, unique_attributes=None):
    """
    Extract sample attributes (SAMPLE_ATTRIBUTE TAG/VALUE pairs) from an XML tree.

    Each attribute's tag and value are the text of its first TAG and first VALUE child (``""``
    for an empty element, as ``findtext`` gives). A tag repeated within one file keeps its last
    value; an attribute without VALUE maps to ``None``, and one without a TAG text is skipped.

    Args:
        tree: XML tree or element (lxml or stdlib ElementTree) to extract data from
        unique_attributes: Attribute names to keep (any container supporting ``in``), or ``None``
            to keep every attribute in the file

    Returns:
        Dictionary with extracted attributes
    """
    sample_attributes = {}
    for attribute in tree.findall(".//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
        # One pass over the children instead of two findtext calls: several times faster in lxml.
        tag = value = None
        for child in attribute:
            if child.tag == "TAG":
                if tag is None:
                    tag = child.text or ""
            elif child.tag == "VALUE" and value is None:
                value = child.text or ""
        if tag and (unique_attributes is None or tag in unique_attributes):
            sample_attributes[tag] = value

    return sample_attributes


def _parse_metadata_file(xml_file: Union[str, Path]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Parse one metadata XML file once, returning its fixed fields and its own sample attributes.

    Raises OSError when the file cannot be read and ``lxml.etree.XMLSyntaxError`` on malformed XML.
    """
    root = parse_xml_file(xml_file)
    return _extract_metadata_fields(root, xml_file), _extract_sample_attributes(root)


def parse_metadata_xml(path: Union[str, Path]) -> Dict[str, Any]:
    """Parse one NCBI efetch metadata XML file into the same fields ``parse_metadata`` extracts per record.

    The single-file counterpart to ``parse_metadata``, used to read one already-downloaded metadata
    file. Raises OSError when the file is missing, ET.ParseError on XML syntax errors, or ValueError
    on extraction errors, so the caller can handle them with appropriate logging. Also folds in this
    file's own SAMPLE_ATTRIBUTE tags (e.g. ``collection_date``), the same way ``parse_metadata``
    does per file, so a single-file parse and a folder-wide parse of the same file agree.
    """
    xml_path = Path(path)
    if not xml_path.is_file():
        raise OSError(f"Metadata XML file not found: {xml_path}")
    try:
        metadata_dict, sample_attributes = _parse_metadata_file(xml_path)
    except etree.XMLSyntaxError as e:
        raise ET.ParseError(f"{xml_path}: {e}") from e
    metadata_dict.update(sample_attributes)
    return metadata_dict


def _metadata_table(field_rows: List[Dict[str, Any]], attribute_rows: List[Dict[str, Any]]) -> pd.DataFrame:
    """Build the metadata table from per-file fixed fields and per-file (sparse) sample attributes.

    Columns are the fixed fields in ``METADATA_FIXED_COLUMNS`` order, then every attribute tag seen
    in any file, sorted. A file without a given attribute gets ``None`` in that column. An attribute
    tag equal to a fixed column name replaces that column's values for every row (with ``None`` for
    files lacking the attribute), as the table has always done.
    """
    row_count = len(field_rows)
    # Attribute columns start as all None and only the cells a file actually has are filled in,
    # so the work grows with the number of attributes present rather than rows times tags.
    attribute_columns: Dict[str, List[Any]] = {}
    for index, attributes in enumerate(attribute_rows):
        for tag, value in attributes.items():
            column = attribute_columns.get(tag)
            if column is None:
                column = attribute_columns[tag] = [None] * row_count
            column[index] = value

    columns: Dict[str, List[Any]] = {}
    for name in METADATA_FIXED_COLUMNS:
        replaced = attribute_columns.get(name)
        columns[name] = replaced if replaced is not None else [row.get(name) for row in field_rows]
    for tag in sorted(attribute_columns.keys() - set(METADATA_FIXED_COLUMNS)):
        columns[tag] = attribute_columns[tag]
    logger.info(f"Found {len(attribute_columns)} unique sample attributes")
    return pd.DataFrame(columns)


def parse_metadata(metadata_folder: Union[str, Path], output_file: Union[str, Path]) -> pd.DataFrame:
    """
    Parse metadata XML files and create a consolidated table.

    Each file is parsed once; its fixed fields and sample attributes are read from the same tree.

    Args:
        metadata_folder: Folder containing metadata XML files
        output_file: Path to save the consolidated metadata table

    Returns:
        DataFrame containing the parsed metadata

    Raises:
        DataAccessError: If parsing fails
    """
    metadata_path = validate_folder(metadata_folder)

    field_rows: List[Dict[str, Any]] = []
    attribute_rows: List[Dict[str, Any]] = []
    error_count = 0

    try:
        xml_files = list_files(metadata_path, "*.xml")

        if not xml_files:
            logger.warning(f"No XML files found in {metadata_path}")
            return pd.DataFrame()

        logger.info(f"Processing {len(xml_files)} metadata files")

        for xml_file in xml_files:
            try:
                fields, sample_attributes = _parse_metadata_file(xml_file)
            except (OSError, etree.XMLSyntaxError, ValueError) as e:
                error_count += 1
                logger.error(f"Error parsing {xml_file}: {e}")
                continue

            field_rows.append(fields)
            attribute_rows.append(sample_attributes)

            # Log progress periodically
            if len(field_rows) % 100 == 0:
                logger.info(f"Processed {len(field_rows)} metadata files")

        if not field_rows:
            logger.warning("No metadata records created")
            return pd.DataFrame()

        metadata_df = _metadata_table(field_rows, attribute_rows)

        # Save to file
        # One chunk: pandas' default chunk size shrinks with the column count, and with one column
        # per sample attribute its per-chunk, per-column conversion dominated the write time.
        write_csv(metadata_df, output_file, sep="\t", index=False, chunksize=len(metadata_df))
        logger.info(f"Saved metadata table with {len(metadata_df)} records to {output_file}")

        return metadata_df

    except (OSError, ValueError, PrepareError) as e:
        raise DataAccessError(f"Error parsing metadata: {e}") from e


def check_metadata_attributes(file_path: Union[str, Path], output_file: Union[str, Path]) -> Dict[str, int]:
    """
    Count occurrences of metadata attributes and save to file.

    Args:
        file_path: Path to the metadata table file
        output_file: Path to save the attribute counts

    Returns:
        Dictionary mapping attribute names to counts

    Raises:
        DataAccessError: If the operation fails
    """
    try:
        # Read metadata table
        df = pd.read_csv(file_path, sep="\t")

        # Filter columns to exclude standard categories
        excluded_prefixes = ("Run_", "Project_", "Sample_", "Experiment_")
        filtered_columns = [col for col in df.columns if not col.startswith(excluded_prefixes)]

        # Count non-null values for each column
        counts = {col: df[col].count() for col in filtered_columns}

        # Sort by count
        sorted_counts = {k: v for k, v in sorted(counts.items(), key=lambda item: item[1], reverse=True)}

        # Save to file
        write_text_atomic(output_file, "".join(f"{key}\t{value}\n" for key, value in sorted_counts.items()))

        logger.info(f"Saved attribute counts to {output_file}")
        return sorted_counts

    # pandas' ParserError and EmptyDataError are ValueErrors.
    except (OSError, ValueError) as e:
        raise DataAccessError(f"Error checking metadata attributes: {e}") from e
