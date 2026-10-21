"""Downloads routed through the shared data store: the store keeps the one copy, the project gets a link.

The ``prepare.store`` imports stay inside the functions because ``prepare.store`` itself
imports this package.

The precheck's refusal of an incomplete copy (``_store_precheck`` without ``--resume-partial`` or
with exhausted refetches) removes the project's store link without holding the dataset lock, on a
sidecar read that may already be stale. This is tolerated: it only matters when two runs of the
same project handle the same accession at once and the other one has just linked a better copy,
and the next run links that copy again.
"""

import functools
import logging
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from prepare.core.exceptions import DataAccessError
from prepare.data.sra import accession as accession_mod
from prepare.data.sra import cleanup as cleanup_mod
from prepare.data.sra import fastq as fastq_mod
from prepare.data.sra import spots as spots_mod
from prepare.utils.lockfile import LockLost, verify_held

logger = logging.getLogger(__name__)

# The publish step moved to ``cleanup.publish_folder`` so plain projects share it; the old
# name stays for callers and tests that use it.
_publish_store_dataset = cleanup_mod.publish_folder


# Sidecar states that mean the store's copy is usable as it stands: verified against NCBI's
# spot count, or downloaded without a spot count to check it against.
STORE_READY_STATES = ("complete", "unverified")

# Prefix ``_link_result`` puts on the message for a dataset the store already held (linked
# rather than downloaded); the single place this string is produced, so every caller that
# needs to tell a link apart from a fresh download (the CLI's result recorder, this module's
# own run summary) matches against this constant rather than a copy of the literal.
STORE_LINKED_PREFIX = "linked from store"

# Starts of the failure messages for a store copy that is kept but not linked: a fetch or kept
# refetch that left the copy incomplete, and the precheck's refusal of an incomplete copy under
# --no-resume-partial. These are settled outcomes, not transient errors: the retry pass does not
# retry them (a second fetch in the same run would gain nothing and use up a refetch), and
# download_sra records the store sidecar's verdict for them.
STORE_INCOMPLETE_PREFIX = "incomplete:"
STORE_PARTIAL_PREFIX = "partial in store;"
# The project's own reads sit at fastq/<ACC> (a real folder, not a link): a store link can
# never replace them, so no fetch is started and the outcome is settled.
STORE_PROJECT_FOLDER_PREFIX = "project folder exists:"
SETTLED_PREFIXES = (STORE_INCOMPLETE_PREFIX, STORE_PARTIAL_PREFIX, STORE_PROJECT_FOLDER_PREFIX)
# The end of the precheck's refusal of an incomplete copy whose refetches are used up: an
# ``incomplete:`` message that started no download (the other ``incomplete:`` messages follow one).
STORE_EXHAUSTED_SUFFIX = "use --accept-partial, or --force to fetch again"

# Refetches of an incomplete store copy that may gain no reads before the copy is no longer
# fetched again on its own: it is then linked with --accept-partial, or refused, until a run
# passes --force. NCBI's spot count is sometimes not reachable by fasterq-dump at all.
STORE_PARTIAL_REFETCH_LIMIT = 2


def _xml_folders(folders, store) -> List[Union[str, Path]]:
    """The folders that may hold ``<ACC>_metadata.xml``: ``folders`` (one, several or None), then the store's."""
    if folders is None:
        candidates: List[Union[str, Path]] = []
    elif isinstance(folders, (str, Path)):
        candidates = [folders]
    else:
        candidates = list(folders)
    if store is not None:
        candidates.append(store.metadata)
    return candidates


def _metadata_xml(folders, accession: str) -> Optional[Path]:
    """The first ``<accession>_metadata.xml`` found in ``folders``, or None."""
    for folder in _xml_folders(folders, None):
        xml = Path(folder) / f"{accession}_metadata.xml"
        if xml.is_file():
            return xml
    return None


def _resolve_expected_spots(accession: str, store, store_metadata, given, verify: bool = True) -> Optional[int]:
    """The spot count a store download of ``accession`` is judged against.

    ``given`` is the caller's count (the project registry's ``run_total_spots``); without one,
    ``prepare.data.sra.spots.expected_spots`` looks in the store sidecar and then in the
    metadata XML folders, the same lookup order every other caller uses. With ``verify``
    off (``--no-verify-downloads``) no count is used at all, so the copy is ``unverified``.
    """
    if not verify:
        return None
    spots = spots_mod.positive_int(given)
    if spots is not None:
        return spots
    return spots_mod.expected_spots(None, accession, store=store, xml_folders=_xml_folders(store_metadata, store))


def _read_store_sidecar(store, accession: str):
    """The sidecar of the store's copy of ``accession``, or None when there is none to read."""
    from prepare.store.layout import sidecar_path
    from prepare.store.sidecar import read_sidecar

    path = sidecar_path(store, accession)
    return read_sidecar(path) if path.is_file() else None


def _short_unverified(sidecar, expected_spots: Optional[int]) -> bool:
    """True when a ready but unverified sidecar records fewer reads than ``expected_spots`` allows.

    Uses the sidecar's ``reads_per_mate`` only, so no FASTQ file is read. A copy downloaded
    before its spot count was known is recorded as ``unverified``; once a count is known, the
    same verdict computation as a fresh download (``verdict_for_count``) decides whether it is short.
    """
    if sidecar is None or sidecar.state not in STORE_READY_STATES:
        return False
    if (sidecar.completeness or {}).get("verdict") not in (None, "unverified"):
        return False
    return spots_mod.verdict_for_count(sidecar.reads_per_mate, expected_spots)["verdict"] == "truncated"


def _store_state(store, accession: str, expected_spots: Optional[int] = None) -> str:
    """What the store holds for ``accession``: ``ready``, ``incomplete`` or ``absent``.

    ``incomplete`` covers a sidecar recording a partial, failed or in-progress download, and
    also files sitting there with no sidecar at all, which is what an interrupted run leaves
    behind and cannot be trusted without re-fetching. With ``expected_spots``, a ready copy
    whose sidecar is ``unverified`` but records fewer reads than that count allows is
    ``incomplete`` too.
    """
    from prepare.store.layout import sidecar_path, sra_dir

    if sidecar_path(store, accession).is_file():
        sidecar = _read_store_sidecar(store, accession)
        state = sidecar.state if sidecar else None
        if state not in STORE_READY_STATES or _short_unverified(sidecar, expected_spots):
            return "incomplete"
        return "ready"
    return "incomplete" if fastq_mod.fastq_files(sra_dir(store, accession)) else "absent"


def _store_settles(store, accession: str, expected_spots: Optional[int] = None) -> bool:
    """True when a run without ``--force`` links or refuses ``accession`` without downloading it.

    That is a ready copy, or an incomplete one whose refetches stopped gaining reads
    (``_refetch_exhausted``). The free-space guard exempts these accessions.
    """
    state = _store_state(store, accession, expected_spots)
    return state == "ready" or (state == "incomplete" and _refetch_exhausted(_read_store_sidecar(store, accession)))


def _reverify_sidecar(store, accession: str, expected_spots: int) -> None:
    """Record a short ``unverified`` store copy as ``partial`` against ``expected_spots``.

    The caller holds the dataset lock. The verdict comes from the sidecar's ``reads_per_mate``
    (no FASTQ file is read); ``ncbi.spots`` is filled in when the sidecar had none. The sidecar
    is rewritten atomically by ``write_sidecar`` and the catalogue row updated (never raises).
    """
    from prepare.store.layout import sidecar_path
    from prepare.store.sidecar import write_sidecar

    sidecar = _read_store_sidecar(store, accession)
    if sidecar is None:
        return
    verdict = spots_mod.verdict_for_count(sidecar.reads_per_mate, expected_spots)
    if not sidecar.ncbi.get("spots"):
        sidecar.ncbi = {**sidecar.ncbi, "spots": expected_spots}
    sidecar.completeness = {"method": verdict["method"], "ratio": verdict["ratio"], "verdict": verdict["verdict"]}
    sidecar.state = "partial" if verdict["verdict"] == "truncated" else "complete"
    write_sidecar(sidecar_path(store, accession), sidecar)
    _catalogue_published(store, sidecar)
    logger.warning(
        "%s: the store copy holds %s of %s spots; recorded as %s",
        accession,
        sidecar.reads_per_mate,
        expected_spots,
        sidecar.state,
    )


def _spot_text(count: Optional[int]) -> str:
    """A spot or read count for a message, ``unknown`` when it is not recorded."""
    return "unknown" if count is None else str(count)


def _incomplete_message(accession: str, sidecar) -> str:
    """The result message for a store copy that is kept for resuming but not linked."""
    reads = sidecar.reads_per_mate if sidecar is not None else None
    spots = spots_mod.positive_int(sidecar.ncbi.get("spots")) if sidecar is not None else None
    return (
        f"{STORE_INCOMPLETE_PREFIX} {accession} store copy holds {_spot_text(reads)} of {_spot_text(spots)} spots; "
        "kept for --resume-partial; rerun with --accept-partial to use it"
    )


def _read_count(value: Any) -> Optional[int]:
    """``value`` as a read count: a whole number of zero or more, else None (unknown)."""
    if isinstance(value, bool):
        return None
    try:
        number = int(value)
    except (TypeError, ValueError):
        return None
    return number if number >= 0 else None


def _refetch_record(previous, new) -> Optional[Dict[str, Any]]:
    """The ``refetch`` record of the store copy after ``new`` (a refetch) was built over ``previous``.

    ``unchanged`` counts consecutive refetches whose ``reads_per_mate`` did not grow beyond the
    larger of the previous copy's and its record's count; a gain resets it to 0. A count of 0 is
    a known count, so refetches that keep returning no reads are counted as unchanged and stop.
    A ``complete`` result, or a first download (no ``previous``), has no record (None).
    ``reads_per_mate`` is the largest count seen so far, ``expected_spots`` the count the
    refetch was judged against.
    """
    if previous is None or new.state == "complete":
        return None
    prior = previous.refetch or {}
    counts = (_read_count(previous.reads_per_mate), _read_count(prior.get("reads_per_mate")))
    known = [n for n in counts if n is not None]
    best = max(known) if known else None
    gained = new.reads_per_mate is not None and (best is None or new.reads_per_mate > best)
    unchanged = 0 if gained else int(prior.get("unchanged") or 0) + 1
    expected = spots_mod.positive_int(new.ncbi.get("spots")) or spots_mod.positive_int(previous.ncbi.get("spots"))
    return {
        "unchanged": unchanged,
        "reads_per_mate": new.reads_per_mate if gained else best,
        "expected_spots": expected,
        "last": datetime.now(timezone.utc).isoformat(),
    }


def _refetch_exhausted(sidecar) -> bool:
    """True when ``STORE_PARTIAL_REFETCH_LIMIT`` refetches in a row gained no reads for this copy."""
    if sidecar is None or not sidecar.refetch:
        return False
    return int(sidecar.refetch.get("unchanged") or 0) >= STORE_PARTIAL_REFETCH_LIMIT


def _exhausted_message(accession: str, sidecar, expected_spots: Optional[int]) -> str:
    """The result message for an incomplete store copy that is no longer refetched without --force."""
    spots = expected_spots or spots_mod.positive_int(sidecar.ncbi.get("spots"))
    return (
        f"{STORE_INCOMPLETE_PREFIX} {sidecar.refetch.get('unchanged')} refetches of {accession} gained no reads "
        f"({_spot_text(sidecar.reads_per_mate)} of {_spot_text(spots)} spots); NCBI's count may not be "
        f"reachable; {STORE_EXHAUSTED_SUFFIX}"
    )


def _drop_store_link(project_fastq: Path, accession: str, store) -> None:
    """Remove the project's symlink to ``accession`` when it points into the store.

    Called when the store's copy is not complete enough to link, so an earlier link never
    presents partial data as present. A real folder (the project's own reads) is never removed.
    """
    from prepare.store.link import is_store_link, unlink_dataset

    if is_store_link(Path(project_fastq) / accession, store) and unlink_dataset(project_fastq, accession):
        logger.info("Removed the project link to the incomplete store copy of %s", accession)


def _link_result(accession: str, project_fastq: Path, store, link_mode: str, note: str) -> Tuple[bool, str]:
    """Link the store's copy into the project and report how many files it holds.

    The caller holds ``accession``'s dataset lock.
    """
    from prepare.store.link import link_dataset

    link = link_dataset(project_fastq, accession, store, mode=link_mode, locked=True)
    return True, f"{STORE_LINKED_PREFIX}{note}, {len(fastq_mod.fastq_files(link))} files"


def _try_link(
    accession: str, project_fastq: Path, store, link_mode: str, note: str, locked: bool
) -> Optional[Tuple[bool, str]]:
    """``_link_result``, taking the dataset lock without waiting unless ``locked`` says it is held.

    Returns None when another run holds the lock, so that ``_store_download`` waits for it and
    decides again under it instead of reporting a ready copy as an error.
    """
    from prepare.store.locks import LockHeld, dataset_lock

    if locked:
        return _link_result(accession, project_fastq, store, link_mode, note)
    try:
        with dataset_lock(store, accession, blocking=False):
            return _link_result(accession, project_fastq, store, link_mode, note)
    except LockHeld:
        logger.debug("%s: dataset lock held by another run; deciding again under the lock", accession)
        return None


def _store_precheck(
    accession: str,
    project_fastq: Path,
    store,
    link_mode: str,
    accept_partial: bool,
    resume_partial: bool,
    expected_spots: Optional[int] = None,
    locked: bool = False,
) -> Optional[Tuple[bool, str]]:
    """Decide what the store alone can settle for ``accession``, without any network call.

    Returns a result when the dataset is already usable (linked), or when it is incomplete
    and the caller asked not to resume it; returns None when it must be downloaded. A ready
    copy that is ``unverified`` but short against ``expected_spots`` is incomplete: with
    ``locked`` (the caller holds the dataset lock) its sidecar is rewritten as ``partial``
    first; without the lock this returns None so the decision is taken again under it. A link
    made without ``locked`` takes the lock without waiting and also returns None when another
    run holds it (``_try_link``).
    """
    if expected_spots is not None and _short_unverified(_read_store_sidecar(store, accession), expected_spots):
        if not locked:
            return None
        _reverify_sidecar(store, accession, expected_spots)
    state = _store_state(store, accession, expected_spots)
    if state == "ready":
        return _try_link(accession, project_fastq, store, link_mode, "", locked)
    if state == "incomplete" and not resume_partial:
        if accept_partial:
            return _try_link(accession, project_fastq, store, link_mode, " (partial)", locked)
        _drop_store_link(project_fastq, accession, store)
        return False, f"{STORE_PARTIAL_PREFIX} rerun with --resume-partial to finish {accession}"
    if state == "incomplete":
        return _exhausted_result(accession, project_fastq, store, link_mode, accept_partial, expected_spots, locked)
    return None


def _exhausted_result(
    accession: str,
    project_fastq: Path,
    store,
    link_mode: str,
    accept_partial: bool,
    expected_spots: Optional[int],
    locked: bool = False,
) -> Optional[Tuple[bool, str]]:
    """The result for an incomplete copy whose refetches stopped gaining reads, or None to fetch it again.

    Once ``STORE_PARTIAL_REFETCH_LIMIT`` refetches in a row gained nothing, the copy is linked
    as partial with ``accept_partial`` and otherwise refused (any store link removed); only
    ``--force``, which skips this precheck, fetches it again.
    """
    sidecar = _read_store_sidecar(store, accession)
    if not _refetch_exhausted(sidecar):
        return None
    if accept_partial:
        return _try_link(accession, project_fastq, store, link_mode, " (partial)", locked)
    _drop_store_link(project_fastq, accession, store)
    return False, _exhausted_message(accession, sidecar, expected_spots)


def _publish_decision(previous, new) -> str:
    """Whether a freshly built store copy (``new`` sidecar) replaces the published one (``previous``).

    Returns ``"publish"`` or ``"keep_previous"``: a refetch that comes back worse never replaces
    a better copy. A ``failed`` result never replaces a copy that is not itself ``failed``; a
    ``complete`` result always publishes; otherwise a copy with fewer reads per mate (judged
    against the same spot count) than the published one is kept out. Without a previous
    sidecar, or over a ``failed`` one, the new copy is published.
    """
    if previous is None:
        return "publish"
    if new.state == "failed":
        return "publish" if previous.state == "failed" else "keep_previous"
    if previous.state == "failed" or new.state == "complete":
        return "publish"
    if previous.reads_per_mate is not None and new.reads_per_mate is not None:
        if new.reads_per_mate < previous.reads_per_mate:
            return "keep_previous"
    return "publish"


def _record_kept_refetch(store, accession: str, previous, new, expected_spots: Optional[int]):
    """Update the published sidecar after a refetch (``new``) that did not replace it; returns it re-read.

    The caller holds the dataset lock. A copy that is ``unverified`` but short against
    ``expected_spots`` is first recorded as ``partial`` (``_reverify_sidecar``); an incomplete
    copy then gets its ``refetch`` counter advanced, written atomically by ``write_sidecar``.
    """
    from prepare.store.layout import sidecar_path
    from prepare.store.sidecar import write_sidecar

    if expected_spots is not None and _short_unverified(previous, expected_spots):
        _reverify_sidecar(store, accession, expected_spots)
        previous = _read_store_sidecar(store, accession) or previous
    if previous.state in STORE_READY_STATES:
        return previous
    previous.refetch = _refetch_record(previous, new)
    write_sidecar(sidecar_path(store, accession), previous)
    return previous


def _keep_previous(
    accession: str,
    project_fastq: Path,
    store,
    link_mode: str,
    accept_partial: bool,
    previous,
    expected_spots: Optional[int] = None,
) -> Tuple[bool, str]:
    """Hand over the published store copy that a worse refetch did not replace.

    A copy that is ``unverified`` but short against ``expected_spots`` is not ready, so it is
    linked only with ``accept_partial``, like any other incomplete copy.
    """
    if previous.state in STORE_READY_STATES and not _short_unverified(previous, expected_spots):
        return _link_result(accession, project_fastq, store, link_mode, " (kept the earlier copy)")
    if accept_partial:
        return _link_result(accession, project_fastq, store, link_mode, " (partial)")
    _drop_store_link(project_fastq, accession, store)
    return False, _incomplete_message(accession, previous)


def _catalogue_published(store, sidecar) -> bool:
    """Record a just-published dataset's sidecar in the store catalogue; never raises.

    ``_store_fetch`` calls this only after the dataset is already complete in ``sra/<ACC>``,
    so the catalogue write is bookkeeping, not part of what makes the download succeed: a
    locked or otherwise unwritable catalogue (``DataAccessError``, e.g. another writer
    holding the catalogue lock, or a store_reindex replay in progress) must not turn a
    correctly published dataset into a reported download failure. On that failure this logs
    a warning naming ``prepare store_reindex`` (which rebuilds the catalogue from the
    sidecars already on disk, including this one) and returns False; a real success returns
    True. Same shape as ``record_usage_safe``/``record_usage_many`` in
    ``prepare.store.usage``, which protect usage recording the same way.
    """
    from prepare.store.catalog import catalog_write

    try:
        with catalog_write(store) as catalog:
            catalog.upsert_dataset(sidecar)
        return True
    except DataAccessError as e:
        logger.warning(
            "Could not record %s in the store catalogue at %s: %s; run prepare store_reindex "
            "--data-root %s to repair it",
            sidecar.accession,
            store.root,
            e,
            store.root,
        )
        return False


def _sidecar_ncbi(accession: str, store, store_metadata, expected_spots, verify: bool) -> Dict[str, Any]:
    """The ``ncbi`` block of a fresh store sidecar: the metadata XML's figures, with ``spots`` resolved.

    ``expected_spots`` (the resolved count of ``_resolve_expected_spots``) fills in ``spots``
    when the XML records none. With ``verify`` off no count is recorded at all: one on file
    would let ``store_verify`` judge a copy the user asked not to judge.
    """
    from prepare.store.sidecar import ncbi_from_metadata_xml

    xml = _metadata_xml(store_metadata, accession) or _metadata_xml(store.metadata, accession)
    ncbi: Dict[str, Any] = ncbi_from_metadata_xml(xml) if xml is not None else {}
    expected = spots_mod.positive_int(expected_spots)
    if not ncbi.get("spots") and expected is not None:
        ncbi = {**ncbi, "spots": expected}
    if not verify:
        return {key: value for key, value in ncbi.items() if key != "spots"}
    return ncbi


def _store_fetch(
    accession: str,
    project_fastq: Path,
    store,
    link_mode: str,
    store_metadata,
    accept_partial: bool = False,
    verify: bool = True,
    **download_kwargs,
) -> Tuple[bool, str]:
    """Download ``accession`` into the store, describe it, and link the project to it.

    Everything happens in the store's ``tmp`` folder: fasterq-dump builds
    ``tmp/<ACC>_temp``, the files are verified and compressed into ``tmp/<ACC>``, and the
    sidecar describing them is written there too. Only then is the finished folder published
    into ``sra/<ACC>`` with one rename, so that folder never holds an unverified or
    sidecar-less dataset for another project to find, and an existing copy is replaced only
    once its replacement is complete.

    The sidecar's ``ncbi.spots`` comes from the metadata XML, or from ``expected_spots`` (the
    resolved spot count of ``_resolve_expected_spots``, in ``download_kwargs``: the registry's,
    the previous sidecar's or the XML's) when the XML records none. A refetch that comes back
    worse than the published copy is discarded (``_publish_decision``), and the published
    sidecar's ``refetch`` counter is advanced instead; a published refetch carries the counter
    in its own sidecar (``_refetch_record``). A copy that is not
    ready (``partial`` or ``failed``) is published so ``--resume-partial`` can finish it, but it
    is linked only with ``accept_partial``; otherwise any project link to it is removed and the
    result is ``(False, "incomplete: ...")``.
    """
    from prepare.store.layout import lock_path, sra_dir
    from prepare.store.link import link_dataset
    from prepare.store.sidecar import build_sidecar, write_sidecar

    target = sra_dir(store, accession)
    replacing = target.exists()
    staged = store.tmp / accession
    # A forced refetch stays forced even when nothing is in the store yet: it must not reuse
    # a cached archive from an earlier attempt.
    download_kwargs["force"] = bool(download_kwargs.get("force")) or replacing
    download_kwargs.setdefault("sra_cache", None)
    if download_kwargs["sra_cache"] is None:
        download_kwargs["sra_cache"] = store.tmp / ".sra-cache"

    # A caller that gave no temp_folder gets one under the store's own tmp, not the system
    # temp directory: fasterq-dump's scratch space can run to several gigabytes per accession,
    # and download_accession only cleans up a temp_folder it created itself (its finally block
    # skips a folder the caller supplied), so this scratch folder is ours to remove afterwards.
    own_scratch: Optional[Path] = None
    if download_kwargs.get("temp_folder") is None:
        own_scratch = store.tmp / f"{accession}_fqtmp"
        download_kwargs["temp_folder"] = own_scratch

    # Whatever an earlier interrupted attempt left staged is not a resume point: the download
    # would otherwise be skipped as "already exists" and that partial copy published.
    cleanup_mod._safe_rmtree(staged)

    try:
        # The sidecar states decide what an incomplete copy means in the store, so the
        # download itself reports a truncated run as done rather than failed here.
        success, message = accession_mod.download_accession(
            accession, store.tmp, staging_folder=store.tmp, accept_partial=True, **download_kwargs
        )
    finally:
        if own_scratch is not None:
            cleanup_mod._safe_rmtree(own_scratch)

    if not success:
        cleanup_mod._safe_rmtree(staged)
        return False, message

    files = fastq_mod.fastq_files(staged)
    compression = "gzip" if any(path.name.endswith(".gz") for path in files) else "none"
    ncbi = _sidecar_ncbi(accession, store, store_metadata, download_kwargs.get("expected_spots"), verify)
    expected = ncbi.get("spots")

    previous = _read_store_sidecar(store, accession) if replacing else None
    sidecar = build_sidecar(accession, staged, ncbi, accession_mod.fasterq_dump_version(), compression)
    sidecar.refetch = _refetch_record(previous, sidecar)
    write_sidecar(staged / f"{accession}.json", sidecar)
    # A holder stalled past the stale window may have lost the lock to another project, which
    # then owns the staged folder too: publish nothing, and leave that folder to its new owner.
    verify_held(lock_path(store, accession))

    if previous is not None and _publish_decision(previous, sidecar) == "keep_previous":
        logger.warning(
            "%s: the new download (%s, %s reads per mate) is not more complete than the store copy "
            "(%s, %s reads per mate); keeping the store copy",
            accession,
            sidecar.state,
            _spot_text(sidecar.reads_per_mate),
            previous.state,
            _spot_text(previous.reads_per_mate),
        )
        cleanup_mod._safe_rmtree(staged)
        previous = _record_kept_refetch(store, accession, previous, sidecar, expected)
        return _keep_previous(accession, project_fastq, store, link_mode, accept_partial, previous, expected)

    cleanup_mod.publish_folder(staged, target, store.tmp)

    catalogued = _catalogue_published(store, sidecar)

    if sidecar.state not in STORE_READY_STATES and not accept_partial:
        _drop_store_link(project_fastq, accession, store)
        return False, _incomplete_message(accession, sidecar)

    link_dataset(project_fastq, accession, store, mode=link_mode, locked=True)
    suffix = "; stored" if catalogued else "; catalogue pending; stored"
    return True, f"{message}{suffix}"


def _store_download(
    accession: str,
    project_fastq: Union[str, Path],
    store,
    link_mode: str = "auto",
    accept_partial: bool = False,
    resume_partial: bool = True,
    store_metadata=None,
    force: bool = False,
    lock_wait: float = 0.0,
    stop: Optional[threading.Event] = None,
    verify: bool = True,
    **download_kwargs,
) -> Tuple[bool, str]:
    """Get ``accession`` for this project through the shared store.

    The store holds one copy of every dataset; the project only ever gets a link to it. A
    copy that is already there is linked without touching the network. Otherwise the store's
    per-accession lock is taken, the sidecar is checked again (another project may have
    finished the download while this one waited), and the download runs into the store.

    The lock is held for the whole download (``prepare.store.locks.dataset_lock``, which
    heartbeats while held), so a second project waits rather than writing the same folder.
    ``lock_wait`` of zero waits for as long as the other project keeps working; a positive
    value gives up after that many seconds, naming the accession and the holder. ``stop`` is
    the run's stop token: when it (or the process-wide ``STOP``) is set, a wait for the lock
    ends with ``(False, "interrupted")``; it is passed on to ``download_accession``. A lock
    lost while the download ran (``verify_held`` just before the publish) ends with
    ``(False, "lock lost: ...")`` and nothing published.
    """
    from prepare.store.locks import LockWaitStopped, dataset_lock

    project_path = Path(project_fastq)
    own = project_path / accession
    if own.exists() and not own.is_symlink():
        # link_dataset would refuse this after the fetch; refusing here starts no fetch, and
        # the outcome is settled (not retried) because a rerun alone cannot change it.
        return False, (
            f"{STORE_PROJECT_FOLDER_PREFIX} {own} holds the project's own reads and is not a store link; "
            f"run store_adopt to move them into the store, or remove the folder"
        )
    for directory in (store.sra, store.tmp, store.locks):
        directory.mkdir(parents=True, exist_ok=True)

    try:
        # One spot count decides completeness everywhere below: the precheck, the download's
        # own verdict message, and the sidecar of a fresh copy.
        expected = _resolve_expected_spots(
            accession, store, store_metadata, download_kwargs.get("expected_spots"), verify=verify
        )
        download_kwargs["expected_spots"] = expected
        precheck = functools.partial(
            _store_precheck,
            accession,
            project_path,
            store,
            link_mode,
            accept_partial,
            resume_partial,
            expected_spots=expected,
        )
        if not force:
            settled = precheck()
            if settled is not None:
                return settled

        should_stop = functools.partial(accession_mod.stop_requested, stop)
        with dataset_lock(store, accession, wait_seconds=lock_wait, should_stop=should_stop):
            if not force:
                settled = precheck(locked=True)
                if settled is not None:
                    return settled
            return _store_fetch(
                accession,
                project_path,
                store,
                link_mode,
                store_metadata,
                accept_partial=accept_partial,
                verify=verify,
                force=force,
                stop=stop,
                **download_kwargs,
            )
    except LockWaitStopped:
        logger.info(f"Stopped waiting for the store lock on {accession}: the run was interrupted")
        return False, "interrupted"
    except LockLost as e:
        logger.error(f"Store download of {accession} not published: {e}")
        return False, f"lock lost: {e}"
    except DataAccessError as e:
        logger.error(f"Store download failed for {accession}: {e}")
        return False, f"store error: {e}"


def _store_downloader(
    store,
    link_mode: str,
    accept_partial: bool,
    resume_partial: bool,
    store_metadata,
    lock_wait: float = 0.0,
    verify: bool = True,
):
    """A ``download_accession``-shaped callable that routes every download through the store.

    The download loops call it with the project's FASTQ folder as ``output_folder``, which is
    where the link is created; the files themselves land in the store. Returns None when no
    store is configured, which leaves the loops downloading into the project folder as before.
    """
    if store is None:
        return None

    def _download(accession, output_folder, num_threads=4, force=False, temp_folder=None, **kwargs):
        return _store_download(
            accession,
            output_folder,
            store,
            link_mode=link_mode,
            accept_partial=accept_partial,
            resume_partial=resume_partial,
            store_metadata=store_metadata,
            force=force,
            lock_wait=lock_wait,
            verify=verify,
            num_threads=num_threads,
            temp_folder=temp_folder,
            **kwargs,
        )

    return _download
