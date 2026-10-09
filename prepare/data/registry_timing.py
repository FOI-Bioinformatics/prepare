"""How long downloads, extractions and assemblies took, recorded in the project registry.

Each of the download, extraction and assembly blocks carries two optional keys: ``started``, when
the step began (ISO 8601 UTC, to the second), and ``seconds``, the wall-clock time it took. Both
are left out of a block that has no timing, so a registry written before these keys existed is
rewritten unchanged. The setters here are called by the commands right after the step's own
``record_*`` call; passing None for both values removes the keys, so a step that did not run (a
dataset linked from the shared store, for example) does not keep the time of an earlier one.

``timing_summary`` gives the counts, totals and medians ``status`` reports.
"""

import statistics
from datetime import datetime, timezone
from time import monotonic
from typing import Any, Dict, List, Optional, Tuple

from prepare.data import registry_blocks as rb
from prepare.data.registry import Registry

# Download states whose recorded time describes an attempt that ran; status counts only these.
TIMED_DOWNLOAD_STATES = ("downloaded", "failed")


def _apply(block: rb.RegistryBlock, started: Optional[str], seconds: Optional[float]) -> None:
    """Set ``started`` and ``seconds`` on ``block``, or remove both when both are None."""
    if started is None and seconds is None:
        block.discard("started")
        block.discard("seconds")
        return
    setattr(block, "started", started)
    setattr(block, "seconds", seconds)


def set_download_timing(registry: Registry, accession: str, started: Optional[str], seconds: Optional[float]) -> None:
    """Record when ``accession``'s download started and how long it took; None for both clears them.

    A no-op when no download is recorded for ``accession``.
    """
    download = rb.download_block(registry, accession)
    if download is None:
        return
    _apply(download, started, seconds)
    rb.set_download_block(registry, accession, download)


def set_extraction_timing(
    registry: Registry, accession: str, genome_id: str, started: Optional[str], seconds: Optional[float]
) -> None:
    """Record the timing of ``accession``'s extraction against ``genome_id``; None for both clears it.

    A no-op when that extraction is not recorded.
    """
    extraction = rb.extraction_block(registry, accession, genome_id)
    if extraction is None:
        return
    _apply(extraction, started, seconds)
    rb.set_extraction_block(registry, accession, genome_id, extraction)


def set_assembly_timing(
    registry: Registry, accession: str, genome_id: str, started: Optional[str], seconds: Optional[float]
) -> None:
    """Record the timing of the assembly of ``accession``'s reads for ``genome_id``; None for both clears it.

    A no-op when that extraction, or its assembly, is not recorded.
    """
    extraction = rb.extraction_block(registry, accession, genome_id)
    if extraction is None or extraction.assembly is None:
        return
    _apply(extraction.assembly, started, seconds)
    rb.set_extraction_block(registry, accession, genome_id, extraction)


class Stopwatch:
    """Start time and elapsed seconds of consecutive steps, in the form the registry records."""

    def __init__(self) -> None:
        """Start timing the first step now."""
        self.restart()

    def restart(self) -> None:
        """Start timing the next step now."""
        self._started = datetime.now(timezone.utc)
        self._clock = monotonic()

    def lap(self) -> Tuple[str, float]:
        """``(started, seconds)`` of the step now ending, then start timing the next one."""
        result = (self._started.isoformat(timespec="seconds"), round(monotonic() - self._clock, 3))
        self.restart()
        return result


def _seconds(block: Any) -> Optional[float]:
    """The ``seconds`` a raw registry block records, or None when it has no numeric value."""
    if not isinstance(block, dict):
        return None
    value = block.get("seconds")
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _summarise(prefix: str, plural: str, values: List[float]) -> Dict[str, Any]:
    """The count, total and median of ``values`` under the keys ``status`` reports."""
    return {
        f"{plural}_timed": len(values),
        f"{prefix}_seconds_total": round(float(sum(values)), 3),
        f"{prefix}_seconds_median": round(statistics.median(values), 3) if values else None,
    }


def timing_summary(registry: Registry) -> Dict[str, Any]:
    """How many downloads, extractions and assemblies have a recorded time, with the total and median.

    Keys are ``downloads_timed``, ``download_seconds_total`` and ``download_seconds_median``, and
    the same for ``extractions`` and ``assemblies``. A median is None when nothing is timed. Only a
    download block whose state is ``downloaded`` or ``failed`` counts: a ``seconds`` left on a
    ``missing`` or ``skipped`` block by a registry written before the writers cleared it is not.
    """
    downloads: List[float] = []
    extractions: List[float] = []
    assemblies: List[float] = []
    for record in registry.datasets.values():
        download = record.get("download")
        timed = isinstance(download, dict) and download.get("state") in TIMED_DOWNLOAD_STATES
        value = _seconds(download) if timed else None
        if value is not None:
            downloads.append(value)
        for extraction in (record.get("extractions") or {}).values():
            value = _seconds(extraction)
            if value is not None:
                extractions.append(value)
            value = _seconds(extraction.get("assembly") if isinstance(extraction, dict) else None)
            if value is not None:
                assemblies.append(value)
    return {
        **_summarise("download", "downloads", downloads),
        **_summarise("extraction", "extractions", extractions),
        **_summarise("assembly", "assemblies", assemblies),
    }
