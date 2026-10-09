"""Mapping from parsed NCBI metadata XML fields to the registry's metadata block fields.

Moved out of ``cli.commands.metadata`` so a registry can be filled in from an already-downloaded
metadata XML folder (``fill_metadata_from_xml``) without going through the CLI layer. The CLI
commands import ``FIELD_COLUMNS`` and ``metadata_fields`` from here so both still agree on one
column-to-field mapping.
"""

import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Mapping, Tuple, Union

from prepare.data import registry_blocks as rb
from prepare.data.metadata import parse_metadata_xml
from prepare.data.registry import Registry, nan_to_none, record_metadata

logger = logging.getLogger(__name__)


# Registry field name and the metadata table column it is read from, as used by metadata_fields.
FIELD_COLUMNS: Tuple[Tuple[str, str], ...] = (
    ("run_size", "Run_Size"),
    ("run_md5", "Run_MD5"),
    ("run_total_spots", "Run_Total_Spots"),
    ("run_total_bases", "Run_Total_Bases"),
    ("assay_type", "Experiment_Library_Strategy"),
    ("organism", "Sample_Scientific_Name"),
    ("collection_date", "collection_date"),
    ("library_layout", "Experiment_Library_Layout"),
    ("platform", "Platform"),
    ("library_strategy", "Experiment_Library_Strategy"),
)


def metadata_fields(row: Mapping[str, Any]) -> Dict[str, Any]:
    """Map one parsed metadata row to registry field names.

    ``row`` is a mapping with keys like ``Run_Total_Spots`` and ``Run_MD5``: the dict
    ``parse_metadata_xml`` returns, or a pandas row from the metadata table ``parse_metadata``
    writes. Returns a dict with registry field names such as ``run_total_spots`` and
    ``run_md5``, as read by ``record_metadata``. Shared by ``DownloadMetadataCommand``,
    ``ParseMetadataCommand`` and ``fill_metadata_from_xml`` so all three agree on the same
    mapping.
    """
    fields: Dict[str, Any] = {}
    for field_name, column in FIELD_COLUMNS:
        value = row.get(column)
        if value is not None:
            fields[field_name] = nan_to_none(value)
    return fields


def metadata_fields_from_xml(xml_path: Union[str, Path]) -> Dict[str, Any]:
    """Parse one metadata XML file into registry field names.

    Wraps ``parse_metadata_xml``. A missing file (``OSError``), malformed XML
    (``xml.etree.ElementTree.ParseError``) or an extraction error (``ValueError``) is logged as a
    warning and treated as no fields (``{}``) rather than raised, so a caller filling many
    accessions can skip one bad file and continue with the rest.
    """
    try:
        row = parse_metadata_xml(xml_path)
    except (OSError, ValueError, ET.ParseError) as error:
        logger.warning("Could not read metadata XML %s: %s", xml_path, error)
        return {}
    return metadata_fields(row)


def _has_positive_spots(fields: Mapping[str, Any]) -> bool:
    """True when ``fields`` carries a ``run_total_spots`` that is a positive whole number."""
    value: Any = fields.get("run_total_spots")
    try:
        return int(value) > 0
    except (TypeError, ValueError):
        return False


def fill_metadata_from_xml(registry: Registry, metadata_folder: Union[str, Path]) -> List[str]:
    """Fill in metadata for every dataset recorded without a spot count, from its XML file.

    For each accession in ``registry`` whose metadata block exists but has ``run_total_spots``
    still ``None``, and whose ``<accession>_metadata.xml`` exists in ``metadata_folder``, parses
    that file and records its fields with ``record_metadata``, keeping the block's ``inferred``
    mark if it was already set. Returns the accessions filled, in registry order. An accession
    with no metadata block, a block already holding a spot count, a missing XML file, an XML
    file that fails to parse, or one that records no positive spot count (e.g. a run NCBI has not
    loaded) is left untouched, so filling the same folder again writes nothing new.
    """
    folder = Path(metadata_folder)
    filled: List[str] = []
    for accession in registry.datasets:
        block = rb.metadata_block(registry, accession)
        if block is None or block.run_total_spots is not None:
            continue
        xml_path = folder / f"{accession}_metadata.xml"
        if not xml_path.is_file():
            continue
        fields = metadata_fields_from_xml(xml_path)
        if not _has_positive_spots(fields):
            continue
        was_inferred = block.inferred
        record_metadata(registry, accession, xml_path, fields)
        if was_inferred:
            rb.mark_inferred(registry, accession, "metadata")
        filled.append(accession)
    return filled
