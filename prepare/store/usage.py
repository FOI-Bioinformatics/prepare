"""
Safe usage recording into the shared store's catalogue.

Every pipeline stage that touches a dataset (a fresh download, a link from the store, an
analysis, a targeted extraction, an assembly) can record that fact in the catalogue's
``usage`` table, so the catalogue answers "which projects used this accession, for which
genome, at which stage" across every project sharing the store. That history is a
convenience, never a dependency: a project must keep working exactly as before when there
is no store, when the project has never been bound to one (no ``store_init``), or when the
catalogue write itself fails (a lock timeout, a corrupt database, a disk error). The two
functions here are the only way any CLI command touches the catalogue for this purpose, and
neither ever raises.
"""

import logging
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from prepare.core.exceptions import DataAccessError
from prepare.data.registry import Registry, load_registry
from prepare.data.registry_blocks import project_block, set_project_block
from prepare.store.catalog import Catalog, catalog_write
from prepare.store.layout import StorePaths

logger = logging.getLogger(__name__)


def _upsert_project_row(catalog, registry: Registry, project_id: str) -> None:
    project = project_block(registry)
    catalog.upsert_project(project_id, project.name or "", project.path or "", str(registry.path))


def ensure_project_identity(registry: Registry) -> Dict[str, Any]:
    """Give ``registry`` a project identity if it has none, and return it.

    A project reaches the store through four discovery rules, and three of them (the
    ``PREPARE_DATA`` environment variable, ``--data-root`` and the user config's default)
    never run ``store_init``, so the registry carries no ``project.id``. Without one the
    project records no usage, and ``store_gc`` then sees a dataset it links as used by nobody.
    Minting the identity on first use closes that gap: every project that links from the store
    has a row in the catalogue, whichever rule found the store.

    The identity is written into ``registry.project`` only; the caller's ``registry_transaction``
    persists it with whatever else that transaction records, and the catalogue row follows from
    the next ``record_usage_safe``/``record_usage_many`` call.
    """
    project = project_block(registry)
    if project.id:
        return project.to_dict()

    cwd = Path.cwd()
    # Keys this function does not own (e.g. "exports" from results_table) are kept.
    project.id = str(uuid.uuid4())
    project.name = project.name or cwd.name
    project.path = project.path or str(cwd.resolve())
    project.created = project.created or datetime.now(timezone.utc).isoformat()
    set_project_block(registry, project)
    logger.info("Recorded this project's identity for the shared store: %s (%s)", project.id, project.name)
    return project.to_dict()


def record_usage_safe(
    paths: Optional[StorePaths],
    registry: Registry,
    accession: str,
    genome_id: str,
    stage: str,
    detail: str = "",
) -> bool:
    """Record one dataset's use in the store catalogue; never raises.

    Returns False, doing nothing else, when there is no store (``paths`` is None) or this
    project has never been bound to one (``registry.project`` carries no ``id``, i.e.
    ``store_init`` was never run for it). Otherwise upserts the project's row (so a project
    that has moved, or one seen for the first time, is always current) and records the
    usage row, both inside one ``catalog_write`` session. Every exception the write can
    raise (``DataAccessError`` from a corrupt or locked catalogue, ``OSError`` from the
    filesystem, or anything else) is caught, logged at warning naming the accession and
    stage, and turned into False: the pipeline stage that triggered this call must never
    fail because the catalogue could not be written.
    """
    project_id = project_block(registry).id
    if paths is None or not project_id:
        return False

    try:
        with catalog_write(paths) as catalog:
            _upsert_project_row(catalog, registry, project_id)
            catalog.record_usage(accession, project_id, genome_id, stage, detail)
        return True
    except Exception as e:
        logger.warning("Could not record usage for %s (stage=%s): %s", accession, stage, e)
        return False


def record_usage_many(
    paths: Optional[StorePaths],
    registry: Registry,
    rows: Iterable[Tuple[str, str, str, str]],
) -> bool:
    """Like ``record_usage_safe``, for several ``(accession, genome_id, stage, detail)`` rows.

    All rows are written inside one ``catalog_write`` session, so a batch (e.g. every
    already-downloaded accession a download run found linked from the store) takes one
    lock instead of one per accession. Same safety as ``record_usage_safe``: no store or no
    bound project returns False without writing anything; any exception during the write
    is caught, logged at warning, and turned into False.
    """
    project_id = project_block(registry).id
    if paths is None or not project_id:
        return False

    rows = list(rows)
    if not rows:
        return True

    try:
        with catalog_write(paths) as catalog:
            _upsert_project_row(catalog, registry, project_id)
            for accession, genome_id, stage, detail in rows:
                catalog.record_usage(accession, project_id, genome_id, stage, detail)
        return True
    except Exception as e:
        logger.warning("Could not record usage batch (%d row(s)): %s", len(rows), e)
        return False


def stale_projects(catalog: Catalog) -> List[Dict[str, Any]]:
    """Projects in the catalogue whose registry no longer keeps them alive.

    A project row is stale when the pipeline run that recorded it will never write to this
    store again: its registry file (the ``registry`` column, written by ``store_init``) has
    been removed, or the registry now found at that path belongs to a different project (its
    ``project.id`` no longer matches this row's ``project_id``, e.g. the project folder was
    reused for a fresh ``store_init``). A registry that exists but fails to parse, or whose
    path cannot even be checked (e.g. a permissions error on a shared multi-user store), is
    treated the same as missing: it cannot vouch for this project either, and the check must
    never crash the caller over one unreadable project.

    Shared by ``store_status`` (which only reports stale projects) and ``store_gc`` (which
    also uses this to decide that a dataset's only usage rows no longer keep it alive).
    Returns each stale row as a dict with ``project_id``, ``name``, ``path``, ``registry``,
    ``hostname`` and a ``reason``, sorted by ``project_id``. The reason separates the two
    cases a reader must tell apart: ``registry missing`` (nothing at that path on this
    machine, which is also what every project on another workstation of a shared store looks
    like from here) and ``project id differs`` (a registry is there, but it now belongs to a
    different project).
    """
    # SELECT *, not a named column list: a catalogue created before ``hostname`` existed has
    # no such column, and a read path must not run DDL to add one.
    rows = [dict(row) for row in catalog.conn.execute("SELECT * FROM projects ORDER BY project_id").fetchall()]
    stale: List[Dict[str, Any]] = []
    for row in rows:
        registry_path = row.get("registry")
        try:
            if not registry_path or not Path(registry_path).is_file():
                stale.append({**row, "reason": "registry missing"})
                continue
            registry = load_registry(registry_path)
        except (OSError, DataAccessError, ValueError) as e:
            logger.warning("Could not check registry %s while checking staleness: %s", registry_path, e)
            stale.append({**row, "reason": "registry unreadable"})
            continue
        if project_block(registry).id != row["project_id"]:
            stale.append({**row, "reason": "project id differs"})
    return stale


def linked_by(paths: StorePaths, catalog: Catalog, accession: str) -> List[str]:
    """Names of every non-stale project that still symlinks ``accession`` to this store.

    Usage rows are the normal record of "a project uses this dataset", but a project can end
    up linking an accession without one: usage tracking predates that project's ``store_link``
    call, or a hook that would have recorded it failed silently. Before ``store_gc`` removes a
    dataset with no (live) usage rows, it must also check the filesystem directly: for every
    project ``stale_projects`` does not consider gone, does its default ``fastq/<accession>``
    folder (relative to the project's recorded ``path``) exist as a symlink resolving to this
    store's copy of the dataset? A project whose path or link cannot be checked (missing,
    unreadable, broken symlink) is silently treated as not linking it, never as an error.
    Returns the matching project names, sorted, or an empty list when none link it.
    """
    stale_ids = {row["project_id"] for row in stale_projects(catalog)}
    target = (paths.sra / accession).resolve()

    rows = catalog.conn.execute("SELECT project_id, name, path FROM projects ORDER BY project_id").fetchall()
    linked: List[str] = []
    for row in rows:
        if row["project_id"] in stale_ids:
            continue
        project_path = row["path"]
        if not project_path:
            continue
        link_path = Path(project_path) / "fastq" / accession
        try:
            if link_path.is_symlink() and link_path.resolve() == target:
                linked.append(row["name"] or row["project_id"])
        except OSError:
            continue
    return sorted(linked)
