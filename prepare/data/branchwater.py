"""
Branchwater data handling module for PREPARE.

This module provides functions for processing Branchwater containment files.
"""

import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, DefaultDict, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from prepare.core.exceptions import DataAccessError, PrepareError
from prepare.core.models import ContainmentSummary
from prepare.core.validation import validate_folder
from prepare.data.file_io import (
    ensure_directory,
    list_files,
    copy_file,
    read_csv,
    visible_files,
    write_csv,
)
from prepare.plugins.base import format_registry
from prepare.plugins.formats.branchwater import BranchWaterFormatPlugin

logger = logging.getLogger(__name__)

# Register format plugins
format_registry.register(BranchWaterFormatPlugin)

# Columns of the details table written alongside the parsed containment table. It carries the
# cANI and sample metadata a Branchwater match brings, which the wide containment matrix (one
# column per genome) has no room for.
DETAILS_COLUMNS = [
    "accession",
    "genome_id",
    "containment",
    "cANI",
    "biosample",
    "bioproject",
    "assay_type",
    "organism",
    "geo_loc_name",
    "lat_lon",
]

# Maps a details column to the key it is read from in a Containment's additional_data
# (the Branchwater CSV column name, as BranchWaterFormatPlugin.parse_file exposes it).
_DETAILS_FIELD_SOURCE = {
    "cANI": "cANI",
    "biosample": "biosample",
    "bioproject": "bioproject",
    "assay_type": "assay_type",
    "organism": "organism",
    "geo_loc_name": "geo_loc_name_country_calc",
    "lat_lon": "lat_lon",
}


def _containment_details_row(accession, genome_id, value, containment) -> Dict[str, Any]:
    """Build one details row from a parsed Containment; a field absent from the CSV is an empty string."""
    additional_data = getattr(containment, "additional_data", None)
    if not isinstance(additional_data, dict):
        additional_data = {}
    row: Dict[str, Any] = {"accession": accession, "genome_id": genome_id, "containment": value}
    for column, source_field in _DETAILS_FIELD_SOURCE.items():
        field_value = additional_data.get(source_field, "")
        row[column] = "" if field_value is None else field_value
    return row


def _resolve_details_file(output_file: Union[str, Path], details_file: Optional[Union[str, Path]]) -> Path:
    """Return the details table path: the given one, or ``<parsed file stem>_details.tsv`` next to it."""
    if details_file is not None:
        return Path(details_file)
    output_path = Path(output_file)
    return output_path.parent / f"{output_path.stem}_details.tsv"


def _write_details_table(details_rows: List[Dict[str, Any]], details_file: Union[str, Path]) -> None:
    """Write the collected details rows as a TSV, in the fixed DETAILS_COLUMNS order."""
    df = pd.DataFrame(details_rows, columns=DETAILS_COLUMNS)
    write_csv(df, details_file, sep="\t", index=False)
    logger.info(f"Containment details saved to {details_file}")


def process_branchwater_files(
    source_folder: Union[str, Path],
    target_folder: Union[str, Path],
    errors: Optional[List[str]] = None,
) -> Dict[str, Path]:
    """
    Process Branchwater files from source folder and copy to target folder.

    Args:
        source_folder: Folder containing Branchwater CSV files
        target_folder: Folder to save processed CSV files
        errors: Optional list the caller inspects afterwards; the name of each file
            that could not be read is appended to it. Left untouched when None.

    Returns:
        Dictionary mapping genome IDs to processed file paths

    Raises:
        DataAccessError: If the process fails
    """
    source_path = validate_folder(source_folder)
    target_path = ensure_directory(target_folder)

    result_files: Dict[str, Path] = {}
    processed_count = 0
    error_count = 0

    # Get all CSV files in the source folder
    csv_files = list_files(source_path, "*.csv")
    if not csv_files:
        logger.warning(f"No CSV files found in {source_path}")
        return result_files

    for csv_file in csv_files:
        output_file = target_path / csv_file.name
        genome_id = csv_file.stem

        # Skip if the output file already exists
        if output_file.exists():
            logger.info(f"Skipping {csv_file.name} as output file already exists")
            result_files[genome_id] = output_file
            processed_count += 1
            continue

        try:
            # Validate branchwater format
            BranchWaterFormatPlugin.validate_header(read_csv(csv_file, nrows=0).columns.tolist())

            # Copy the file
            result_file = copy_file(csv_file, output_file)
            logger.info(f"Processed branchwater file: {csv_file.name}")
            result_files[genome_id] = result_file
            processed_count += 1

        except (OSError, PrepareError) as e:
            error_count += 1
            logger.error(f"Error processing {csv_file}: {e}")
            if errors is not None:
                errors.append(csv_file.name)

    level = logger.warning if error_count else logger.info
    level(f"Processed {processed_count} files with {error_count} errors")
    return result_files


def _process_branchwater_row(row, metadata_records):
    """
    Process a single row from a Branchwater CSV and extract metadata.

    Args:
        row: Row data from Branchwater CSV
        metadata_records: List to append metadata records to
    """
    if "acc" not in row:
        logger.warning("Row missing 'acc' column. Skipping.")
        return

    # Basic metadata record
    metadata_record = {"Run_ID": row["acc"]}

    # Add any other available metadata fields
    field_mapping = {
        "biosample": "Sample_ID",
        "bioproject": "Project_ID",
        "organism": "Sample_Scientific_Name",
        "geo_loc_name_country_calc": "geo_loc_name_country_calc",
        "collection_date_sam": "collection_date_sam",
        "assay_type": "assay_type",
        "lat_lon": "lat_lon",
        "cANI": "cANI",
    }

    for source_field, target_field in field_mapping.items():
        if source_field in row and not pd.isna(row[source_field]):
            metadata_record[target_field] = row[source_field]

    metadata_records.append(metadata_record)


def _validate_branchwater_file(csv_file):
    """
    Validate if a file is in Branchwater format.

    Args:
        csv_file: Path to CSV file

    Returns:
        bool: True if file is valid, False otherwise
    """
    # Read the CSV file directly to check its format
    with open(csv_file, "r") as f:
        header_line = f.readline().strip()
        columns = [col.strip() for col in header_line.split(",")]

    # Check if this is a valid Branchwater format file
    if "acc" not in columns or "containment" not in columns:
        logger.warning(f"File {csv_file} does not appear to be in Branchwater format. Skipping.")
        return False
    return True


def extract_metadata_from_branchwater(
    branchwater_folder: Union[str, Path],
    output_file: Union[str, Path],
    errors: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Extract metadata from Branchwater CSV files and save to a file.

    Args:
        branchwater_folder: Path to the folder containing Branchwater CSV files
        output_file: Path to save metadata CSV file
        errors: Optional list the caller inspects afterwards; the name of each file
            that could not be read is appended to it. Left untouched when None.

    Returns:
        DataFrame containing extracted metadata

    Raises:
        DataAccessError: If the extraction fails
    """
    try:
        source_path = Path(branchwater_folder)
        metadata_records: List[Dict[str, Any]] = []
        processed_count = 0
        error_count = 0

        if not source_path.exists():
            raise DataAccessError(f"Branchwater folder does not exist: {source_path}")

        # Get all CSV files in the source folder
        csv_files = visible_files(source_path, "*.csv")
        if not csv_files:
            logger.warning(f"No CSV files found in {source_path}")
            return pd.DataFrame()

        for csv_file in csv_files:
            try:
                logger.info(f"Processing file: {csv_file.name}")

                # Validate file format
                if not _validate_branchwater_file(csv_file):
                    continue

                # Read the CSV with pandas
                df = pd.read_csv(csv_file)

                # Extract metadata from each row
                for _, row in df.iterrows():
                    _process_branchwater_row(row, metadata_records)

                processed_count += 1
                logger.info(f"Extracted metadata from {csv_file}")

            # pandas' ParserError and EmptyDataError, and UnicodeDecodeError, are ValueErrors.
            except (OSError, ValueError) as e:
                error_count += 1
                logger.error(f"Error extracting metadata from {csv_file}: {e}")
                if errors is not None:
                    errors.append(csv_file.name)

        return _finalize_metadata_extraction(metadata_records, output_file, processed_count, error_count)

    except (OSError, ValueError, PrepareError) as e:
        raise DataAccessError(f"Error extracting metadata from Branchwater files: {e}") from e


def _finalize_metadata_extraction(metadata_records, output_file, processed_count, error_count):
    """
    Create and save the final metadata DataFrame.

    Args:
        metadata_records: List of metadata records
        output_file: Path to save the output
        processed_count: Number of processed files
        error_count: Number of errors encountered

    Returns:
        pd.DataFrame: The metadata DataFrame
    """
    level = logger.warning if error_count else logger.info
    level(f"Processed {processed_count} files with {error_count} errors")

    if not metadata_records:
        logger.warning("No metadata records extracted")
        # Create an empty DataFrame but save it anyway
        metadata_df = pd.DataFrame()
    else:
        # Create DataFrame from records
        metadata_df = pd.DataFrame(metadata_records)
        # Remove duplicates
        if "Run_ID" in metadata_df.columns:
            metadata_df = metadata_df.drop_duplicates(subset=["Run_ID"])
        # A Branchwater-API-shaped CSV (acc/containment/cANI columns only) carries no SRA
        # metadata, so every record here has nothing beyond Run_ID and cANI; that looks like
        # a successful extraction unless we say so explicitly.
        extra_columns = set(metadata_df.columns) - {"Run_ID", "cANI"}
        if not extra_columns:
            logger.warning(
                "No metadata beyond Run_ID and cANI was found; the source CSVs look like "
                "Branchwater API output (acc/containment/cANI only), not an SRA metadata download. "
                "Run download_metadata to fetch biosample, bioproject and other SRA fields."
            )

    # Save to output file
    # Ensure directory exists
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    write_csv(metadata_df, output_file, sep="\t", index=False)

    if not metadata_records:
        logger.info(f"Saved empty metadata file to {output_file}")
    else:
        logger.info(f"Extracted metadata saved to {output_file} with {len(metadata_df)} records")

    return metadata_df


def _process_genome_containments(csv_file, genome_id, containment_data, details_rows=None):
    """
    Process containment data from a single genome file.

    Args:
        csv_file: Path to CSV file
        genome_id: ID of the genome
        containment_data: Dictionary to store containment data
        details_rows: Optional list to append a details row (cANI and sample metadata) to,
            one per containment. Left untouched when None (the default), so existing
            callers keep their prior behaviour.

    Raises:
        DataAccessError: If processing fails
    """
    try:
        # Determine file format and use appropriate plugin
        df_headers = read_csv(csv_file, nrows=0)
        headers = df_headers.columns.tolist()

        BranchWaterFormatPlugin.validate_header(headers)

        # Parse file to get containments
        containments = BranchWaterFormatPlugin.parse_file(csv_file, genome_id)

        # Add to containment data dictionary
        for containment in containments:
            if genome_id not in containment_data[containment.accession]:
                containment_data[containment.accession][genome_id] = containment.value
            else:
                # Keep maximum containment value if multiple entries
                containment_data[containment.accession][genome_id] = max(
                    containment_data[containment.accession][genome_id],
                    containment.value,
                )

            if details_rows is not None:
                details_rows.append(
                    _containment_details_row(containment.accession, genome_id, containment.value, containment)
                )

    except (OSError, ValueError, PrepareError) as e:
        raise DataAccessError(f"Error processing {csv_file}: {e}") from e


def _threshold_counts(max_values: np.ndarray, step_size: float) -> Tuple[List[float], List[int]]:
    """Thresholds from 1.0 down to 0.0 in ``step_size`` steps, and how many rows reach each.

    A row counts at a threshold when its maximum containment is at or above it. One sorted copy
    of the maxima answers every threshold with a binary search.
    """
    raw = np.arange(int(1 / step_size), -1, -1) * step_size
    ordered = np.sort(max_values)
    reached = len(ordered) - np.searchsorted(ordered, raw, side="left")
    return [round(float(threshold), 2) for threshold in raw], [int(count) for count in reached]


def _positive_mappings(df: pd.DataFrame) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """Genome -> accessions and accession -> genomes for every cell above 0, both in table order."""
    genomes = [column for column in df.columns if column not in ("max_containment", "max_containment_annotation")]
    accessions = df.index.to_numpy(dtype=object)
    positive = df[genomes].to_numpy(dtype=float) > 0
    genome_to_samples = {genome: accessions[positive[:, j]].tolist() for j, genome in enumerate(genomes)}
    rows, cols = np.nonzero(positive)
    per_row = np.split(np.asarray(genomes, dtype=object)[cols], np.searchsorted(rows, np.arange(1, len(df))))
    sample_to_genomes = {accession: names.tolist() for accession, names in zip(accessions.tolist(), per_row)}
    return genome_to_samples, sample_to_genomes


def _generate_containment_summary(containment_data, output_file, summary_file, step_size):
    """
    Generate summary data from containment data.

    Args:
        containment_data: Dictionary of containment data
        output_file: Path to save parsed containment data
        summary_file: Path to save containment summary
        step_size: Step size for threshold calculation

    Returns:
        ContainmentSummary object

    Raises:
        DataAccessError: If generation fails
    """
    if step_size <= 0 or step_size > 1.0:
        raise DataAccessError("step_size must be between 0 and 1.0")
    if step_size < 0.001:
        raise DataAccessError("step_size too small (minimum 0.001)")

    try:
        # Create DataFrame from containment data
        df = pd.DataFrame.from_dict(containment_data, orient="index")

        # Fill NA values with 0
        df.fillna(0, inplace=True)

        # Add max_containment column
        df["max_containment"] = df.max(axis=1)

        # Add max_containment_annotation column
        df["max_containment_annotation"] = df.idxmax(axis=1)

        # Sort by max_containment. A stable sort keeps tied accessions in their input order, so the
        # table (and everything recorded from it) is the same on every platform; the default
        # quicksort orders ties differently depending on the CPU's sort kernel.
        df.sort_values(by="max_containment", ascending=False, inplace=True, kind="stable")

        # Save parsed containment data
        write_csv(df, output_file, sep="\t")
        logger.info(f"Parsed containment data saved to {output_file}")

        thresholds, counts = _threshold_counts(df["max_containment"].to_numpy(dtype=float), step_size)

        # Save summary data
        summary_df = pd.DataFrame({"Threshold": thresholds, "Count": counts})
        write_csv(summary_df, summary_file, sep="\t", index=False)
        logger.info(f"Containment summary saved to {summary_file}")

        genome_to_samples, sample_to_genomes = _positive_mappings(df)
        return ContainmentSummary(
            thresholds=thresholds,
            counts=counts,
            max_containment=dict(zip(df.index, df["max_containment"].tolist())),
            genome_to_samples=genome_to_samples,
            sample_to_genomes=sample_to_genomes,
            table=df,
        )

    except (OSError, ValueError, PrepareError) as e:
        raise DataAccessError(f"Error generating containment summary: {e}") from e


def parse_containment_data(
    matches_folder: Union[str, Path],
    output_file: Union[str, Path],
    summary_file: Union[str, Path],
    step_size: float = 0.1,
    details_file: Optional[Union[str, Path]] = None,
    errors: Optional[List[str]] = None,
) -> ContainmentSummary:
    """
    Parse containment data from match files and generate summary.

    Alongside the parsed containment table, writes a details table (default
    ``<output_file stem>_details.tsv``) with one row per (accession, genome_id)
    carrying the cANI and sample metadata a Branchwater match brings; the wide
    containment matrix itself is unchanged.

    Args:
        matches_folder: Folder containing match files
        output_file: Path to save parsed containment data
        summary_file: Path to save containment summary
        step_size: Step size for threshold calculation
        details_file: Path to save the details table (default: derived from output_file)
        errors: Optional list the caller inspects afterwards; the name of each file
            that could not be read is appended to it. Left untouched when None.

    Returns:
        ContainmentSummary object

    Raises:
        DataAccessError: If the parsing fails
    """
    matches_path = validate_folder(matches_folder)

    # Dictionary to store containment data
    containment_data: DefaultDict[str, Dict[str, float]] = defaultdict(dict)
    details_rows: List[Dict[str, Any]] = []

    # Get all CSV files in the matches folder
    csv_files = list_files(matches_path, "*.csv")
    if not csv_files:
        logger.warning(f"No CSV files found in {matches_path}")
        return ContainmentSummary()

    processed_count = 0
    error_count = 0

    for csv_file in csv_files:
        genome_id = csv_file.stem

        try:
            _process_genome_containments(csv_file, genome_id, containment_data, details_rows=details_rows)
            processed_count += 1
        except PrepareError as e:
            error_count += 1
            logger.error(f"Error parsing containment from {csv_file}: {e}")
            if errors is not None:
                errors.append(csv_file.name)

    level = logger.warning if error_count else logger.info
    level(f"Processed {processed_count} files with {error_count} errors")

    _write_details_table(details_rows, _resolve_details_file(output_file, details_file))

    if not containment_data:
        logger.warning("No valid containment data found")
        return ContainmentSummary()

    return _generate_containment_summary(containment_data, output_file, summary_file, step_size)
