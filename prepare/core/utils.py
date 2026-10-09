"""
Core utility functions for PREPARE.
"""

import logging

from prepare.core.exceptions import ProcessingError

logger = logging.getLogger(__name__)

# Known metadata columns that are not genome identifiers
_KNOWN_METADATA_COLUMNS = {"max_containment", "max_containment_annotation"}


def get_genome_columns(df):
    """
    Return genome columns from a DataFrame by excluding known metadata columns.

    Args:
        df: DataFrame with containment data

    Returns:
        List of column names that represent genomes

    Raises:
        ProcessingError: If no genome columns found
    """
    genome_columns = [col for col in df.columns if col not in _KNOWN_METADATA_COLUMNS]

    if not genome_columns:
        raise ProcessingError("No genome columns found in summary file")

    return genome_columns
