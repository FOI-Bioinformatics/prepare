"""Genome download functionality using NCBI datasets CLI.

This module provides functions for downloading genome assemblies from NCBI
using the datasets command-line tool, extracting FASTA files, and creating
manifests suitable for downstream analysis tools such as sourmash.
"""

import logging
import os
import re
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from prepare.core.constants import GENOME_ACCESSION_PATTERN, GENOME_ACCESSION_PREFIXES
from prepare.core.exceptions import DataAccessError
from prepare.data.file_io import copy_file
from prepare.utils.security import SecureSubprocess

logger = logging.getLogger(__name__)


def _validate_genome_accession(accession: str) -> str:
    """Validate a genome assembly accession (GCF_/GCA_ format).

    Args:
        accession: The accession string to validate.

    Returns:
        The validated accession.

    Raises:
        DataAccessError: If the accession format is invalid.
    """
    if not accession.startswith(GENOME_ACCESSION_PREFIXES):
        raise DataAccessError(f"Invalid genome accession prefix: {accession}. " f"Expected GCF_ or GCA_ prefix.")
    if not re.match(GENOME_ACCESSION_PATTERN, accession):
        raise DataAccessError(
            f"Invalid genome accession format: {accession}. " f"Expected format: GCF_000000000.0 or GCA_000000000.0"
        )
    return accession


def genome_fasta_path(accession: str, output_dir: Path) -> Path:
    """Return the organized FASTA path an accession is extracted to."""
    return Path(output_dir) / f"{accession}.fna"


def partition_present_genomes(accessions: List[str], output_dir: Path) -> Tuple[List[str], List[str]]:
    """Split accessions into (already present on disk, missing) by their .fna file."""
    present: List[str] = []
    missing: List[str] = []
    for acc in accessions:
        (present if genome_fasta_path(acc, output_dir).exists() else missing).append(acc)
    return present, missing


def download_genomes(
    accessions: List[str],
    output_dir: Path,
    include: str = "genome",
    assembly_level: Optional[str] = None,
) -> Path:
    """Download genome FASTA files using NCBI datasets CLI.

    Args:
        accessions: List of genome accessions (GCF_/GCA_ format).
        output_dir: Directory for the downloaded zip file.
        include: Data type to include (default: genome).
        assembly_level: Optional assembly level filter (e.g., complete, chromosome).

    Returns:
        Path to the downloaded zip file.

    Raises:
        DataAccessError: If validation or download fails.
    """
    if not accessions:
        raise DataAccessError("No accessions provided for download.")

    for acc in accessions:
        _validate_genome_accession(acc)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    zip_path = output_dir / "ncbi_dataset.zip"

    args = ["download", "genome", "accession"]
    args.extend(accessions)
    args.extend(["--include", include])
    args.extend(["--filename", str(zip_path)])

    if assembly_level:
        args.extend(["--assembly-level", assembly_level])

    logger.info(f"Downloading {len(accessions)} genome(s) from NCBI datasets")

    try:
        SecureSubprocess.run_secure("datasets", args)
    except Exception as e:
        raise DataAccessError(f"Failed to download genomes: {e}")

    if not zip_path.exists():
        raise DataAccessError("Download completed but zip file was not created.")

    logger.info(f"Downloaded genomes to {zip_path}")
    return zip_path


class MissingGenomesError(DataAccessError):
    """Requested genomes were absent from the download; ``genome_paths`` holds those that were extracted."""

    def __init__(self, message: str, genome_paths: Dict[str, Path]):
        super().__init__(message)
        self.genome_paths = genome_paths


def _require_expected(genome_paths: Dict[str, Path], expected: Optional[Sequence[str]]) -> None:
    """Raise DataAccessError naming requested accessions that the zip did not contain."""
    if expected is None:
        return
    missing = [acc for acc in expected if acc not in genome_paths]
    if missing:
        raise MissingGenomesError(
            f"{len(missing)} requested genome(s) missing from the NCBI download: {', '.join(missing)}", genome_paths
        )


def extract_and_organize(
    zip_path: Path,
    output_dir: Path,
    expected: Optional[Sequence[str]] = None,
) -> Dict[str, Path]:
    """Extract downloaded zip and organize FASTA files by accession.

    The NCBI datasets zip typically contains files at:
    ncbi_dataset/data/{accession}/{accession}_{assembly}_genomic.fna

    Each FASTA is copied to output_dir/{accession}.fna for convenience.

    Args:
        zip_path: Path to the downloaded NCBI datasets zip file.
        output_dir: Directory to place organized FASTA files.
        expected: Accessions that were requested. When given, accessions absent from the
            zip are reported in a DataAccessError after the others are extracted.

    Returns:
        Dict mapping accession to the organized FASTA file path.

    Raises:
        DataAccessError: If the zip is invalid, extraction fails, or an expected
            accession is missing from the zip.
    """
    zip_path = Path(zip_path)
    output_dir = Path(output_dir)

    if not zip_path.exists():
        raise DataAccessError(f"Zip file not found: {zip_path}")

    if not zipfile.is_zipfile(zip_path):
        raise DataAccessError(f"Not a valid zip file: {zip_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    genome_paths: Dict[str, Path] = {}

    extract_dir = zip_path.parent / f"ncbi_extract_{os.getpid()}"

    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(extract_dir)

        # Find FASTA files in the extracted data
        data_dir = extract_dir / "ncbi_dataset" / "data"
        if not data_dir.exists():
            raise DataAccessError("Unexpected zip structure: ncbi_dataset/data/ not found.")

        for accession_dir in sorted(data_dir.iterdir()):
            if not accession_dir.is_dir():
                continue
            accession = accession_dir.name
            # Skip non-accession directories (e.g., assembly_data_report)
            if not accession.startswith(GENOME_ACCESSION_PREFIXES):
                continue

            # Find FASTA files (.fna or .fasta)
            fasta_files = sorted(accession_dir.glob("*.fna")) + sorted(accession_dir.glob("*.fasta"))
            if not fasta_files:
                logger.warning(f"No FASTA files found for {accession}")
                continue

            # Use the first (typically only) FASTA file
            source_fasta = fasta_files[0]
            dest_fasta = output_dir / f"{accession}.fna"
            copy_file(source_fasta, dest_fasta)
            genome_paths[accession] = dest_fasta
            logger.debug(f"Extracted {accession} -> {dest_fasta}")

    except zipfile.BadZipFile as e:
        raise DataAccessError(f"Corrupt zip file {zip_path}: {e}")
    except DataAccessError:
        raise
    except Exception as e:
        raise DataAccessError(f"Error extracting genomes: {e}")
    finally:
        # Clean up temporary extraction directory
        if extract_dir.exists():
            shutil.rmtree(extract_dir, ignore_errors=True)

    logger.info(f"Extracted and organized {len(genome_paths)} genome(s)")
    _require_expected(genome_paths, expected)
    return genome_paths
