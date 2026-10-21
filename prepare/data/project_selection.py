"""Which recorded selection applies to a genome.

``select_datasets`` records one project-wide ``selected`` flag per accession with the criteria of
the run that set it. Extraction and the ``status --next`` suggestions both need to know whether
that selection was made for the genome at hand, so the rule lives here, once.
"""

import logging
from typing import TYPE_CHECKING, Any, Dict, Set

from prepare.data import registry_blocks as rb

if TYPE_CHECKING:
    from prepare.data.registry import Registry

logger = logging.getLogger(__name__)


def selection_for_genome(registry: "Registry", genome_id: str) -> Set[str]:
    """The recorded selection when its criteria cover ``genome_id``, else an empty set.

    The selection is project-wide (one flag per accession), so one made for another genome's
    column must not narrow this genome's extraction; it applies when it was made on this
    genome's column, on ``max_containment``, or on a ``--genome-ids`` list that names it.
    """
    from prepare.data.registry import query  # registry_scan imports read_extraction, which imports this module

    selected = set(query(registry, "selected"))
    if not selected:
        return set()
    criteria: Dict[str, Any] = {}
    for acc in selected:
        block = rb.selection_block(registry, acc)
        if block is not None and block.criteria:
            criteria = dict(block.criteria)
            break
    column = criteria.get("column")
    covered = column in (None, genome_id, "max_containment") or genome_id in (criteria.get("genome_ids") or [])
    if covered:
        return selected
    logger.info("The recorded selection was made for %s, not %s; the table decides", column, genome_id)
    return set()
