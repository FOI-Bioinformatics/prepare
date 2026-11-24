"""Cross-stage funnel counts for `status`: how far a project's datasets got through the pipeline.

`funnel` reduces a registry to one dict per pipeline stage (screened, selected, downloaded,
analysed, extracted, assembled), each holding at least an ``accessions`` count that matches
`prepare.data.registry.stage_members`; the download, extraction and assembly stages also carry
bytes, seconds and/or pairs read from the registry's raw blocks in the same pass.
`prepare.cli.commands.status.render_text` renders the result as one text line, and the 0.9.0
`project_report` command reads it to build the project-level report.
"""

from typing import Any, Callable, Dict, List, Optional

from prepare.data.registry import Registry, stage_members
from prepare.data.registry_timing import TIMED_DOWNLOAD_STATES


def _seconds(value: Any) -> Optional[float]:
    """A raw registry field as a float ``seconds`` value, or None when it is missing or not numeric."""
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return None
    return float(value)


def _total_seconds(values: List[float]) -> Optional[float]:
    """The rounded sum of ``values``, or None when there are none (nothing in the stage was timed)."""
    return round(sum(values), 3) if values else None


# A check whether one raw extraction record's work for a stage is still on disk (see status_report).
OnDisk = Optional[Callable[[Dict[str, Any], str], bool]]


def _tally_download(download: Dict[str, Any], tally: Dict[str, Any], counted: bool = True) -> None:
    """Add one dataset's raw download block to ``tally`` (bytes, failed count and the two time lists).

    A finished download is added only when ``counted`` (it is among the downloaded members); a
    failed one always is.
    """
    state = download.get("state")
    if state == "downloaded" and not counted:
        return
    if state == "downloaded":
        tally["bytes"] += int(download.get("bytes_total") or 0)
    elif state == "failed":
        tally["failed"] += 1
    if state in TIMED_DOWNLOAD_STATES:
        seconds = _seconds(download.get("seconds"))
        if seconds is not None:
            tally["seconds" if state == "downloaded" else "failed_seconds"].append(seconds)


def _tally_extractions(
    record: Dict[str, Any], work: Dict[str, Any], extracted: bool, assembled: bool, on_disk: OnDisk = None
) -> None:
    """Add one dataset's extractions with mapped reads (when ``extracted``) and assemblies with
    contigs (when ``assembled``) to ``work``: pair counts, assembled bases and the time lists.

    ``on_disk(extraction, stage)``, when given, leaves out each pair whose own work is gone from disk,
    so an accession kept in a stage by one genome does not count another genome's deleted pair.
    """
    for extraction in (record.get("extractions") or {}).values():
        if not isinstance(extraction, dict):
            continue
        if on_disk is not None:
            extracted_here = extracted and on_disk(extraction, "extracted")
            assembled_here = assembled and on_disk(extraction, "assembled")
        else:
            extracted_here, assembled_here = extracted, assembled
        if extracted_here and (extraction.get("mapped_reads") or 0) > 0:
            work["extracted"] += 1
            _append_seconds(work["extracted_seconds"], extraction.get("seconds"))
        assembly = extraction.get("assembly")
        if assembled_here and isinstance(assembly, dict) and (assembly.get("contigs") or 0) > 0:
            work["assembled"] += 1
            work["total_bp"] += int(assembly.get("total_bp") or 0)
            _append_seconds(work["assembled_seconds"], assembly.get("seconds"))


def _append_seconds(times: List[float], value: Any) -> None:
    """Append ``value`` to ``times`` when it is a usable number of seconds."""
    seconds = _seconds(value)
    if seconds is not None:
        times.append(seconds)


def funnel(
    registry: Registry, members: Optional[Dict[str, List[str]]] = None, on_disk: OnDisk = None
) -> Dict[str, Any]:
    """How many of ``registry``'s datasets made it through each stage of the pipeline.

    ``members`` is ``stage_members(registry)`` when the caller already built it (``status``
    reuses the one it built for the report's ``stages`` block); every ``accessions`` count below
    is ``len(members[stage])``, so the funnel's counts never drift from ``status``'s own stage
    counts. One further pass over ``registry.datasets`` reads the rest, straight from each
    dataset's raw ``download``/``extractions`` dict rather than building a typed block per
    dataset:

    - ``downloaded``: ``bytes`` sums ``bytes_total`` over the datasets in the "downloaded"
      stage (a download that is not in that state never has files on record, so its
      ``bytes_total`` is always 0); ``seconds`` sums the recorded download time of the same
      datasets, so all three figures describe one set; ``failed`` counts datasets whose download
      state is "failed" and ``failed_seconds`` sums their recorded time (the other state a
      download's ``seconds`` is recorded for, see
      ``prepare.data.registry_timing.TIMED_DOWNLOAD_STATES``).
    - ``extracted``/``assembled``: ``pairs`` counts the (accession, genome) extractions and
      assemblies that meet the same mapped-reads/contigs threshold
      ``prepare.data.registry._extraction_stage`` uses to decide stage membership (one
      accession can hold more than one such pair, so ``pairs`` can exceed ``accessions``);
      ``seconds`` and (for assemblies) ``total_bp`` are summed over that same set of pairs, not
      every recorded extraction/assembly attempt, so an extraction with zero mapped reads (or an
      assembly with zero contigs) contributes to neither.

    A stage's ``seconds`` is None, not 0, when nothing in it was ever timed.
    """
    if members is None:
        members = stage_members(registry)

    download_tally: Dict[str, Any] = {"bytes": 0, "failed": 0, "seconds": [], "failed_seconds": []}
    work: Dict[str, Any] = {
        "extracted": 0,
        "extracted_seconds": [],
        "assembled": 0,
        "total_bp": 0,
        "assembled_seconds": [],
    }

    # The figures are summed over the same accessions the counts take from ``members``, so
    # work that status found gone from disk adds to neither.
    downloaded, extracted, assembled = (set(members[stage]) for stage in ("downloaded", "extracted", "assembled"))
    for acc, record in registry.datasets.items():
        download = record.get("download")
        if isinstance(download, dict):
            _tally_download(download, download_tally, counted=acc in downloaded)
        if acc in extracted or acc in assembled:
            _tally_extractions(record, work, acc in extracted, acc in assembled, on_disk)

    return {
        "screened": {"accessions": len(members["screened"])},
        "selected": {"accessions": len(members["selected"]), "excluded": len(members["excluded"])},
        "downloaded": {
            "accessions": len(members["downloaded"]),
            "bytes": download_tally["bytes"],
            "seconds": _total_seconds(download_tally["seconds"]),
            "failed": download_tally["failed"],
            "failed_seconds": _total_seconds(download_tally["failed_seconds"]),
        },
        "analysed": {"accessions": len(members["analysed"])},
        "extracted": {
            "accessions": len(members["extracted"]),
            "pairs": work["extracted"],
            "seconds": _total_seconds(work["extracted_seconds"]),
        },
        "assembled": {
            "accessions": len(members["assembled"]),
            "pairs": work["assembled"],
            "total_bp": work["total_bp"],
            "seconds": _total_seconds(work["assembled_seconds"]),
        },
    }
