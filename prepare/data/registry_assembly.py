"""What an assembly was built from, and whether a recorded assembly still matches its reads.

An assembly block may carry an optional ``inputs`` key describing the extracted reads and
settings it was built from. It is left out while unset, so a registry written before the key
existed is rewritten unchanged. ``record_extraction`` drops a recorded assembly, because a new
extraction replaces the reads the assembly was built from.

For an assembly recorded without ``inputs`` (one made by an earlier version), two checks give a
conservative view of whether it is stale: ``assembly_predates_extraction`` compares the two
recorded dates, and ``legacy_assembly_current`` also compares the recorded megahit preset and
minimum contig length with the requested ones (a value missing on either side matches anything).
"""

from datetime import datetime
from typing import Any, Dict, Optional, Tuple

from prepare.data import registry_blocks as rb
from prepare.data.registry import Registry


def _blocks(
    registry: Registry, accession: str, genome_id: str
) -> Tuple[Optional[rb.ExtractionBlock], Optional[rb.AssemblyBlock]]:
    """The extraction block of ``accession`` against ``genome_id`` and its assembly, either may be None."""
    extraction = rb.extraction_block(registry, accession, genome_id)
    if extraction is None:
        return None, None
    return extraction, extraction.assembly


def set_assembly_inputs(registry: Registry, accession: str, genome_id: str, inputs: Optional[Dict[str, Any]]) -> None:
    """Record what the assembly of ``accession``'s reads for ``genome_id`` was built from.

    None removes the key. A no-op when that extraction, or its assembly, is not recorded.
    """
    extraction, assembly = _blocks(registry, accession, genome_id)
    if extraction is None or assembly is None:
        return
    if inputs is None:
        assembly.discard("inputs")
    else:
        assembly.inputs = dict(inputs)
    rb.set_extraction_block(registry, accession, genome_id, extraction)


def _earlier(first: str, second: str) -> bool:
    """True when timestamp ``first`` is earlier than ``second``.

    Both are parsed as ISO 8601, so two times written with different UTC offsets (before and
    after a daylight saving change) compare correctly; text that does not parse, or a mix of
    times with and without an offset, is compared as plain strings.
    """
    try:
        return datetime.fromisoformat(first) < datetime.fromisoformat(second)
    except (ValueError, TypeError):
        return str(first) < str(second)


def assembly_predates_extraction(registry: Registry, accession: str, genome_id: str) -> bool:
    """True when the recorded assembly is dated earlier than the recorded extraction.

    False when either block, or either date, is missing or empty.
    """
    extraction, assembly = _blocks(registry, accession, genome_id)
    if extraction is None or assembly is None:
        return False
    if not extraction.date or not assembly.date:
        return False
    return _earlier(assembly.date, extraction.date)


def _matches(recorded: Any, requested: Any) -> bool:
    """True when the values are equal, or when either is missing (a missing value matches anything)."""
    return recorded is None or requested is None or recorded == requested


def legacy_assembly_current(
    registry: Registry, accession: str, genome_id: str, preset: Optional[str], min_contig_len: Optional[int]
) -> bool:
    """True when a recorded assembly without ``inputs`` can be taken as current for these settings.

    That is: an assembly is recorded, it is not dated earlier than the extraction (see
    ``assembly_predates_extraction``), and its recorded ``params.preset`` and
    ``params.min_contig_len`` equal ``preset`` and ``min_contig_len``. Only a value missing on
    either side (None) matches anything: a recorded preset, ``"default"`` included, that differs
    from the requested one is not current, so the assembly is built again with the requested one.
    """
    _, assembly = _blocks(registry, accession, genome_id)
    if assembly is None or assembly_predates_extraction(registry, accession, genome_id):
        return False
    params = assembly.params if isinstance(assembly.params, dict) else {}
    return _matches(params.get("preset"), preset) and _matches(params.get("min_contig_len"), min_contig_len)
