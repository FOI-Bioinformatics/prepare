"""Recording of download outcomes in the project registry, for the ``sra-download`` command.

The functions here are module-level so the command class stays small. The command keeps thin
wrappers with the former method names and passes them in, so the recorder looks up its helpers
through the command, as it did before the split.
"""

import functools
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from prepare.cli.commands.sra_verdicts import recorded_verdict, store_sidecar_verdict
from prepare.data import registry_blocks as rb
from prepare.data.registry import Registry, record_download, update_linked
from prepare.data.registry_batch import RegistryBatch
from prepare.data.registry_timing import set_download_timing
from prepare.data.sra import ALREADY_EXISTS, SETTLED_PREFIXES, STORE_LINKED_PREFIX, parse_verdict_message, result_body
from prepare.data.sra.run_report import Timings, attempt_timing
from prepare.store.layout import StorePaths
from prepare.store.usage import ensure_project_identity, record_usage_many

# Marker the data layer puts in a result message for a dataset this run downloaded and
# saved into the store (as opposed to STORE_LINKED_PREFIX, imported above, for one the
# store already held).
STORE_SAVED_SUFFIX = "; stored"


def sidecar_completeness(store: Optional[StorePaths], accession: str) -> Optional[dict]:
    """The completeness verdict recorded in the store's sidecar for ``accession``, or None.

    None when there is no store, or no sidecar (not yet catalogued, or unreadable);
    ``read_sidecar`` already logs a warning for the latter case.
    """
    return store_sidecar_verdict(store, accession)


def record_present(
    reg: Registry,
    acc: str,
    fastq_dir: Path,
    from_store: bool,
    complete: Optional[dict],
    spots: Optional[int] = None,
) -> None:
    """Record an accession found already downloaded, unless its record already says so.

    ``complete`` is the store sidecar's completeness verdict for a store-linked accession, or
    ``present_verdicts``'s for one in a folder of its own, merged with the one on file by
    ``recorded_verdict`` against ``spots``; the caller reads either before the lock is taken.
    """
    if from_store:
        ensure_project_identity(reg)
    if (rb.download_block(reg, acc) or rb.DownloadBlock()).state == "downloaded":
        return
    complete = recorded_verdict(reg, acc, complete, spots, linked=from_store)
    record_download(
        reg,
        acc,
        "downloaded",
        fastq_dir,
        attempt=False,
        complete=complete,
        source="store" if from_store else None,
        store_name=acc if from_store else None,
    )
    if from_store:
        update_linked(reg, acc, add=True)


def result_recorder(
    fastq_dir: Path,
    store: Optional[StorePaths],
    batch: RegistryBatch,
    timings: Optional[Timings] = None,
    expected_spots: Optional[Dict[str, int]] = None,
    record_found: Callable[..., None] = record_present,
    sidecar_of: Callable[[Optional[StorePaths], str], Optional[dict]] = sidecar_completeness,
) -> Callable[[str, bool, str], None]:
    """The callback the download loop uses to record each accession's outcome.

    Each outcome is queued on ``batch`` (which ``_run`` opens around the download call and
    flushes every 50 results, every 30 seconds, and on exit including an interrupt) rather
    than written in a transaction of its own.

    A dataset linked from the shared store is recorded as downloaded without counting an
    attempt against it, since no download ran; one this run downloaded into the store
    counts as an attempt like any other. Both are added to the project's list of linked
    datasets, and recorded as store catalogue usage: ``"linked"`` for a dataset the store
    already held, ``"downloaded"`` for one this run saved into it. The usage rows are
    written after each registry flush, once its lock is released: the catalogue has its own
    lock, and waiting for it while holding the project's registry lock can time a
    concurrent writer's registry write out.

    ``timings`` is the dict ``download_sra`` fills before each ``on_result`` call; the
    attempt's start and seconds are recorded with the outcome, and a dataset linked from
    the store, which no download produced, has any earlier time removed. The verdict recorded
    is ``recorded_verdict`` of the message's (or, for a store dataset whose message carries
    none, the sidecar's) and the one on file, with ``expected_spots`` as the spot count.
    """
    timings = timings if timings is not None else {}
    spot_counts = expected_spots or {}
    usage_rows: List[Tuple[str, str, str, str]] = []

    def _record_usage(reg: Registry) -> None:
        if usage_rows:
            record_usage_many(store, reg, list(usage_rows))
            usage_rows.clear()

    batch.add_flush_hook(_record_usage)

    def _record_result(accession: str, success: bool, message: str) -> None:
        body = result_body(message)
        if success and body == ALREADY_EXISTS:
            # Found in place (typically after waiting for another run that downloaded it):
            # not an attempt, and an existing record keeps the other run's message.
            found = functools.partial(record_found, acc=accession, fastq_dir=fastq_dir, from_store=False, complete=None)
            batch.apply(found, label=accession)
            return
        linked = bool(success) and body.startswith(STORE_LINKED_PREFIX)
        from_store = linked or (bool(success) and body.endswith(STORE_SAVED_SUFFIX))
        if from_store:
            usage_rows.append((accession, "", "linked" if linked else "downloaded", message))
        complete = parse_verdict_message(message) if success else None
        # A store link or a store copy kept unlinked (SETTLED_PREFIXES) records its sidecar's verdict
        # (sra_verdicts.store_verdict); the sidecar is read here, outside the flush's lock.
        use_sidecar = complete is None and (from_store or (not success and body.startswith(SETTLED_PREFIXES)))
        sidecar = sidecar_of(store, accession) if use_sidecar else None
        started, seconds = attempt_timing(accession, success, message, timings) or (None, None)

        def _mutation(reg: Registry) -> None:
            if from_store:
                ensure_project_identity(reg)
            verdict = recorded_verdict(reg, accession, sidecar or complete, spot_counts.get(accession), bool(sidecar))
            record_download(
                reg,
                accession,
                "downloaded" if success else "failed",
                fastq_dir,
                message,
                attempt=not linked,
                complete=verdict,
                source="store" if from_store else None,
                store_name=accession if from_store else None,
            )
            set_download_timing(reg, accession, started, seconds)
            if from_store:
                update_linked(reg, accession, add=True)

        batch.apply(_mutation, label=accession)

    return _record_result
