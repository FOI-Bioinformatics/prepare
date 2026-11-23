"""Which recorded selection applies to a genome.

``select_datasets`` records one project-wide ``selected`` flag per accession with the criteria of
the run that set it. Extraction and the ``status --next`` suggestions both need to know which of
the selected accessions were selected for the genome at hand, so the rule lives here, once.
"""

import logging
from typing import TYPE_CHECKING, Any, Dict, Optional, Set

from prepare.core.constants import DEFAULT_CONTAINMENT_THRESHOLD
from prepare.data import registry_blocks as rb

if TYPE_CHECKING:
    from prepare.data.registry import Registry

logger = logging.getLogger(__name__)


def _covers(criteria: Dict[str, Any], genome_id: str) -> bool:
    """Whether a selection made with ``criteria`` was made for ``genome_id``.

    True when the criteria name no column, name this genome's column or ``max_containment``, or
    name this genome in a ``--genome-ids`` list.
    """
    column = criteria.get("column")
    return column in (None, genome_id, "max_containment") or genome_id in (criteria.get("genome_ids") or [])


def selection_for_genome(registry: "Registry", genome_id: str) -> Set[str]:
    """The selected accessions whose recorded criteria cover ``genome_id``; an empty set when none do.

    The selection is project-wide (one flag per accession), so one made for another genome's
    column must not narrow this genome's extraction. Each selected accession is judged by the
    criteria of its own selection block, so accessions selected by different runs (for example
    downloads kept by ``--skip-downloaded``) give the same answer on every call.
    """
    from prepare.data.registry import query  # registry_scan imports read_extraction, which imports this module

    selected = set(query(registry, "selected"))
    covered: Set[str] = set()
    for acc in selected:
        block = rb.selection_block(registry, acc)
        criteria = dict(block.criteria) if block is not None and block.criteria else {}
        if _covers(criteria, genome_id):
            covered.add(acc)
    if len(covered) < len(selected):
        logger.info(
            "%d of %d selected accession(s) were selected for a genome other than %s; left out",
            len(selected) - len(covered),
            len(selected),
            genome_id,
        )
    return covered


def selection_threshold(registry: "Registry", genome_id: str) -> Optional[float]:
    """The lowest containment threshold recorded by the selections that cover ``genome_id``, or None.

    ``extract_target_reads`` uses it when no ``--threshold`` is given, so a sample selected below
    the default threshold is still extracted. None when nothing is selected for the genome or no
    covering selection recorded a threshold.
    """
    from prepare.data.registry import query  # registry_scan imports read_extraction, which imports this module

    thresholds = []
    for acc in query(registry, "selected"):
        block = rb.selection_block(registry, acc)
        value = _threshold_for(dict(block.criteria) if block is not None and block.criteria else {}, genome_id)
        if value is not None:
            thresholds.append(value)
    return min(thresholds) if thresholds else None


def _threshold_for(criteria: Dict[str, Any], genome_id: str) -> Optional[float]:
    """The threshold one selection's ``criteria`` set for ``genome_id``, or None when they set none.

    A threshold on max_containment or on any of --genome-ids bounds another genome's value, not
    this one's, so it may lower the default but never raise it.
    """
    value = criteria.get("threshold")
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not _covers(criteria, genome_id):
        return None
    bounds_this = criteria.get("column") == genome_id or (
        genome_id in (criteria.get("genome_ids") or []) and criteria.get("require") == "all"
    )
    return float(value) if bounds_this else min(float(value), DEFAULT_CONTAINMENT_THRESHOLD)
