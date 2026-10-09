"""
Metadata counting for PREPARE.

This module provides functions for counting and analyzing metadata.
"""

import logging
import os
import pandas as pd
from pathlib import Path
from typing import Optional, Union

from prepare.core.exceptions import ProcessingError
from prepare.core.utils import get_genome_columns as _get_genome_columns
from prepare.data.file_io import write_csv

logger = logging.getLogger(__name__)


def _validate_metadata_column(metadata_df, metadata_column):
    """
    Validate that metadata column exists.

    Args:
        metadata_df: Metadata DataFrame
        metadata_column: Column name to validate

    Raises:
        ProcessingError: If column doesn't exist
    """
    if metadata_column not in metadata_df.columns:
        available_columns = ", ".join(metadata_df.columns)
        raise ProcessingError(
            f"Column '{metadata_column}' not found in metadata file. " f"Available columns: {available_columns}"
        )


def _process_genome_accessions(genome_column, summary_df, threshold, metadata_df, metadata_column, df_list):
    """
    Process accessions for a single genome.

    Args:
        genome_column: Genome column name
        summary_df: Summary DataFrame
        threshold: Containment threshold
        metadata_df: Metadata DataFrame
        metadata_column: Metadata column to count
        df_list: List to append count DataFrames to

    Returns:
        Number of samples processed
    """
    try:
        # Find accessions with containment at or above threshold
        selected_accessions = summary_df[summary_df[genome_column] >= threshold].index

        # Skip if no matching accessions
        if len(selected_accessions) == 0:
            logger.warning(f"No accessions found for {genome_column} at or above threshold {threshold}")
            return 0

        # Filter metadata to selected accessions
        selected_metadata = metadata_df[metadata_df.index.isin(selected_accessions)]

        # Skip if no matching metadata
        if selected_metadata.empty:
            logger.warning(f"No metadata found for {genome_column} accessions")
            return 0

        # Count values in metadata column
        count_series = selected_metadata[metadata_column].value_counts()

        # Create DataFrame with counts
        count_df = pd.DataFrame({genome_column: count_series})
        df_list.append(count_df)

        return len(selected_accessions)

    except Exception as e:
        logger.error(f"Error processing genome {genome_column}: {e}")
        return 0


def count_metadata(
    summary_file: Union[str, Path],
    metadata_file: Union[str, Path],
    metadata_column: str,
    threshold: float,
    output_file: Union[str, Path],
    stat_file: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """
    Count occurrences of metadata values by genome.

    Args:
        summary_file: Path to the summary file
        metadata_file: Path to the metadata file
        metadata_column: Column in metadata file to count
        threshold: Minimum containment threshold
        output_file: Path to save the result table
        stat_file: Path to save the statistics file

    Returns:
        DataFrame with count results

    Raises:
        ProcessingError: If the operation fails
    """
    try:
        # Generate statistics file path if not provided
        if stat_file is None:
            base_name, ext = os.path.splitext(output_file)
            stat_file = f"{base_name}_stats{ext}"

        # Load data
        summary_df = pd.read_csv(summary_file, sep="\t", index_col=0)
        metadata_df = pd.read_csv(metadata_file, sep="\t", index_col=0)

        # Validate metadata column exists
        _validate_metadata_column(metadata_df, metadata_column)

        # Get genome columns
        genome_columns = _get_genome_columns(summary_df)
        logger.info(f"Processing {len(genome_columns)} genome columns")

        # Create result dataframes for each genome
        df_list: list = []
        unique_accessions = set()
        processed_count = 0

        for genome_column in genome_columns:
            selected_count = _process_genome_accessions(
                genome_column,
                summary_df,
                threshold,
                metadata_df,
                metadata_column,
                df_list,
            )
            unique_accessions.update(summary_df[summary_df[genome_column] >= threshold].index)
            if selected_count > 0:
                processed_count += 1

        # Skip if no data frames created
        if not df_list:
            logger.warning("No count data generated")
            return pd.DataFrame()

        # Combine all dataframes
        result_df = pd.concat(df_list, axis=1)

        # Fill NA values with 0 and convert to int
        result_df = result_df.fillna(0).astype(int)

        # Save to file
        write_csv(result_df, output_file, sep="\t")
        logger.info(f"Saved count table to {output_file}")

        # Calculate and save column sums
        column_sums = result_df.sum().sort_values(ascending=False)
        write_csv(column_sums, stat_file, sep="\t", header=False)
        logger.info(f"Saved column statistics to {stat_file}")

        # Log summary information
        total_counts = result_df.values.sum()
        logger.info(f"Total counts in table: {total_counts}")
        logger.info(f"Unique accessions after filtering: {len(unique_accessions)}")

        return result_df

    except Exception as e:
        if isinstance(e, ProcessingError):
            raise
        raise ProcessingError(f"Error counting metadata: {e}")
