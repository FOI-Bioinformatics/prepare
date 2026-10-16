"""
Validation utilities for PREPARE.

This module provides functions for validating input data and configurations.
"""

import logging
from pathlib import Path
from typing import Optional, Union

from prepare.core.exceptions import DataAccessError, ValidationError

logger = logging.getLogger(__name__)


def validate_accession(accession: str) -> bool:
    """
    Validate if a string looks like a valid SRA accession.

    Args:
        accession: The accession string to validate

    Returns:
        True if the accession is valid, False otherwise
    """
    # Basic validation - SRA accessions follow specific patterns
    from prepare.core.constants import SRA_VALID_PREFIXES

    valid_prefixes = SRA_VALID_PREFIXES

    if not accession:
        return False

    if not accession.startswith(valid_prefixes):
        return False

    # Check if the rest is numeric
    suffix = accession[3:]
    if not suffix.isdigit():
        return False

    return True


def validate_containment_value(value: Union[str, float]) -> Optional[float]:
    """
    Validate and convert a containment value.

    Args:
        value: The containment value to validate

    Returns:
        The validated containment value as a float, or None if invalid
    """
    try:
        float_value = float(value)

        # Containment should be between 0 and 1
        if 0 <= float_value <= 1:
            return float_value
        else:
            logger.warning(f"Containment value {float_value} outside expected range [0, 1]")
            return None
    except (ValueError, TypeError):
        logger.warning(f"Invalid containment value: {value}")
        return None


def validate_folder(folder_path: Union[str, Path], create: bool = False) -> Path:
    """
    Validate that a folder exists and is accessible.

    Args:
        folder_path: Path to the folder to validate
        create: If True, create the folder if it doesn't exist

    Returns:
        The validated folder path as a Path object

    Raises:
        ValidationError: If the path is not a directory, or does not exist and create is False
        DataAccessError: If the folder cannot be created (an operating-system failure, not invalid input)
    """
    path = Path(folder_path)

    if path.exists():
        if not path.is_dir():
            raise ValidationError(f"{path} exists but is not a directory")
    elif create:
        try:
            path.mkdir(parents=True)
            logger.info(f"Created directory: {path}")
        except Exception as e:
            raise DataAccessError(f"Failed to create directory {path}: {e}") from e
    else:
        raise ValidationError(f"Directory does not exist: {path}")

    return path
