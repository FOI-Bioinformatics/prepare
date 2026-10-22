"""
Metadata handling for PREPARE.

This module provides functions for downloading and processing metadata from NCBI.
"""

import logging
import os
import xml.etree.ElementTree as ET
from pathlib import Path
import time  # noqa: F401 - tests patch time.sleep through this module
import pandas as pd
from typing import Any, Dict, List, Optional, Tuple, Union

from Bio import Entrez  # noqa: F401 - tests patch Entrez.efetch and Entrez.email through this module
from lxml import etree

from prepare.core import settings
from prepare.core.exceptions import DataAccessError, PrepareError, NetworkError
from prepare.core.validation import validate_folder
from prepare.data.accessions import read_accessions_file
from prepare.data.metadata_fetch import (  # noqa: F401 - re-exported; the Entrez download code moved there
    _ENTREZ_LOCK,
    _PACERS,
    _PACERS_LOCK,
    MAX_BATCH_SIZE,
    MAX_RETRIES,
    _download_accessions_individually,
    _download_batch_metadata,
    _download_single_metadata,
    _entrez_credentials,
    _pace_requests,
    _split_efetch_packages,
    _write_metadata_file,
)
from prepare.data.file_io import ensure_directory, list_files, write_csv, write_text_atomic
from prepare.utils.progress import ProgressReporter
from prepare.utils.xml import parse_xml_file

logger = logging.getLogger(__name__)


def _get_unique_accessions(matches_folder, threshold):
    """
    Extract unique accessions from match files.

    Args:
        matches_folder: Folder containing match files
        threshold: Minimum containment threshold

    Returns:
        Set of unique accessions

    Raises:
        DataAccessError: If the folder holds no match CSV or a match file cannot be read
    """
    unique_accessions: set = set()
    matches_path = Path(matches_folder)

    # Find all CSV files in matches folder
    csv_files = list_files(matches_path, "*.csv")

    if not csv_files:
        raise DataAccessError(f"No match CSV in {matches_path}")

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
                raise ValueError("no 'acc' or 'SRA accession' column")

            # Filter by threshold
            if threshold > 0:
                df = df[df[containment_col] >= threshold]

            # Add accessions to set
            unique_accessions.update(df[accession_col].tolist())

        # pandas' ParserError and EmptyDataError are ValueErrors; KeyError is a missing
        # containment column; TypeError is a containment column holding text.
        except (OSError, ValueError, KeyError, TypeError) as e:
            raise DataAccessError(f"Cannot read match file {csv_file}: {e}") from e

    return unique_accessions


def _read_accessions_file(path: Union[str, Path]) -> List[str]:
    """Read non-empty, non-comment accession lines from a file, one accession per line."""
    return read_accessions_file(path)


def download_metadata(
    email: str,
    matches_folder: Union[str, Path],
    metadata_folder: Union[str, Path],
    threshold: float = 0.0,
    dry_run: bool = False,
    accessions_file: Optional[Union[str, Path]] = None,
    api_key: Optional[str] = None,
    batch_size: int = 200,
    failures: Optional[Dict[str, str]] = None,
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
        failures: Optional mapping the caller inspects afterwards; every accession that
            could not be fetched is added with the reason. Left untouched when None.

    Returns:
        Dictionary mapping accessions to metadata file paths

    Raises:
        ValueError: If batch_size is outside 1 to 500.
        DataAccessError: If the download fails
        NetworkError: If NCBI could not be reached for any accession
    """
    if not 1 <= batch_size <= MAX_BATCH_SIZE:
        raise ValueError(f"batch_size must be between 1 and {MAX_BATCH_SIZE}, got {batch_size}")

    try:
        # A dry run writes nothing, so it leaves a missing metadata folder missing.
        metadata_path = Path(metadata_folder) if dry_run else ensure_directory(metadata_folder)

        if accessions_file:
            unique_accessions: set = set(_read_accessions_file(accessions_file))
        else:
            matches_path = validate_folder(matches_folder)
            unique_accessions = _get_unique_accessions(matches_path, threshold)

        total_accessions = len(unique_accessions)
        logger.info(f"Found {total_accessions} unique accessions")

        # Check which accessions need downloading, from one listing of the folder
        existing = set(os.listdir(metadata_path)) if metadata_path.is_dir() else set()
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
            accessions_to_download,
            metadata_path,
            email,
            to_download_count,
            api_key=api_key,
            batch_size=batch_size,
            failures=failures,
        )

    except NetworkError:
        raise
    except (OSError, ValueError, PrepareError) as e:
        raise DataAccessError(f"Error downloading metadata: {e}") from e


def _download_accessions_metadata(
    accessions_to_download: List[str],
    metadata_path: Path,
    email: str,
    total_count: int,
    api_key: Optional[str] = None,
    batch_size: int = 200,
    failures: Optional[Dict[str, str]] = None,
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
        failures: Optional mapping that receives every failed accession and its reason

    Returns:
        Dictionary mapping accessions to metadata file paths, for successes only.
    """
    # Each efetch call sets Entrez.email/api_key itself under _ENTREZ_LOCK (_entrez_credentials).
    result_files: Dict[str, Path] = {}
    if failures is None:
        failures = {}
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

        # SRA URL: the "SRA Normalized" file wherever listed; with several files and none so named, the last
        named = [e for e in srafile_elements if e.get("semantic_name") == "SRA Normalized"]
        normalized = named[0] if named else (srafile_elements[-1] if len(srafile_elements) > 1 else None)
        sra_normalized_url = normalized.get("url") if normalized is not None else None

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


def parse_metadata(
    metadata_folder: Union[str, Path], output_file: Union[str, Path], errors: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Parse metadata XML files and create a consolidated table.

    Each file is parsed once; its fixed fields and sample attributes are read from the same tree.
    A file that cannot be parsed is left out of the table and named in ``errors``.

    Args:
        metadata_folder: Folder containing metadata XML files
        output_file: Path to save the consolidated metadata table
        errors: Optional list the caller inspects afterwards; the name of each file that
            could not be parsed is appended to it. Left untouched when None.

    Returns:
        DataFrame containing the parsed metadata

    Raises:
        DataAccessError: If the folder holds no metadata XML or the table cannot be written
    """
    metadata_path = validate_folder(metadata_folder)

    field_rows: List[Dict[str, Any]] = []
    attribute_rows: List[Dict[str, Any]] = []

    try:
        xml_files = list_files(metadata_path, "*.xml")

        if not xml_files:
            raise DataAccessError(f"No metadata XML in {metadata_path}")

        logger.info(f"Processing {len(xml_files)} metadata files")

        for xml_file in xml_files:
            try:
                fields, sample_attributes = _parse_metadata_file(xml_file)
            except (OSError, etree.XMLSyntaxError, ValueError) as e:
                logger.error(f"Error parsing {xml_file}: {e}")
                if errors is not None:
                    errors.append(xml_file.name)
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

    except DataAccessError:
        raise
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
