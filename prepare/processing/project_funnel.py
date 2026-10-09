"""Cross-stage funnel counts for `status`: how far a project's datasets got through the pipeline.

`funnel` reduces a registry to one dict per pipeline stage (screened, selected, downloaded,
analysed, extracted, assembled), each holding at least an ``accessions`` count that matches
`prepare.data.registry.stage_members`; the download, extraction and assembly stages also carry
bytes, seconds and/or pairs read from the registry's raw blocks in the same pass.
`prepare.cli.commands.status.render_text` renders the result as one text line, and the 0.9.0
`project_report` command reads it to build the project-level report.
"""

from typing import Any, Dict, List, Optional

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


def _tally_download(download: Dict[str, Any], tally: Dict[str, Any]) -> None:
    """Add one dataset's raw download block to ``tally`` (bytes, failed count and the two time lists)."""
    state = download.get("state")
    if state == "downloaded":
        tally["bytes"] += int(download.get("bytes_total") or 0)
    elif state == "failed":
        tally["failed"] += 1
    if state in TIMED_DOWNLOAD_STATES:
        seconds = _seconds(download.get("seconds"))
        if seconds is not None:
            tally["seconds" if state == "downloaded" else "failed_seconds"].append(seconds)


def funnel(registry: Registry, members: Optional[Dict[str, List[str]]] = None) -> Dict[str, Any]:
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
    extracted_pairs = 0
    extraction_seconds: List[float] = []
    assembled_pairs = 0
    assembled_total_bp = 0
    assembly_seconds: List[float] = []

    for record in registry.datasets.values():
        download = record.get("download")
        if isinstance(download, dict):
            _tally_download(download, download_tally)
        for extraction in (record.get("extractions") or {}).values():
            if not isinstance(extraction, dict):
                continue
            if (extraction.get("mapped_reads") or 0) > 0:
                extracted_pairs += 1
                seconds = _seconds(extraction.get("seconds"))
                if seconds is not None:
                    extraction_seconds.append(seconds)
            assembly = extraction.get("assembly")
            if isinstance(assembly, dict) and (assembly.get("contigs") or 0) > 0:
                assembled_pairs += 1
                assembled_total_bp += int(assembly.get("total_bp") or 0)
                seconds = _seconds(assembly.get("seconds"))
                if seconds is not None:
                    assembly_seconds.append(seconds)

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
            "pairs": extracted_pairs,
            "seconds": _total_seconds(extraction_seconds),
        },
        "assembled": {
            "accessions": len(members["assembled"]),
            "pairs": assembled_pairs,
            "total_bp": assembled_total_bp,
            "seconds": _total_seconds(assembly_seconds),
        },
    }
