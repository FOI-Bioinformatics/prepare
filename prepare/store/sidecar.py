"""
Per-dataset sidecar: ``<root>/sra/<ACC>/<ACC>.json``.

The sidecar makes the store self-describing: what files exist for one
downloaded accession, their size, md5 and read count, the run layout
(``PAIRED``/``SINGLE``), NCBI's recorded facts about the run, whether the
download is complete against those facts, and whatever per-run statistics a
later step computes. Nothing here talks to NCBI or the network; it only
reads what is already on disk plus a caller-supplied ``ncbi`` dict.
"""

import hashlib
import json
import logging
import xml.etree.ElementTree as ET
import zlib
from dataclasses import asdict, dataclass, field, fields as dataclass_fields
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from prepare.core.exceptions import DataAccessError
from prepare.data.file_io import write_text_atomic
from prepare.data.metadata import parse_metadata_xml
from prepare.data.sra import (
    MATE1_SUFFIXES,
    MATE_SUFFIXES,
    fastq_digest,
    fastq_files,
    fastq_stem,
    orphan_fastq,
    primary_fastq,
    verify_download,
)

logger = logging.getLogger(__name__)

# Bump when the sidecar's fields change shape; a sidecar written by an older schema is still
# read back by Sidecar.from_dict, which fills in defaults for whatever is missing.
SIDECAR_SCHEMA = 1

# Suffixes marking the second mate of a pair, i.e. MATE_SUFFIXES minus MATE1_SUFFIXES.
_MATE2_SUFFIXES = tuple(suffix for suffix in MATE_SUFFIXES if suffix not in MATE1_SUFFIXES)


@dataclass
class Sidecar:
    """One dataset's sidecar record.

    Every field has a default so ``from_dict`` tolerates a sidecar written by an older or
    partial schema, or one missing keys because a later processing step has not run yet.
    """

    accession: str = ""
    state: str = "unknown"
    layout: str = "SINGLE"
    downloaded: Optional[str] = None
    tool: str = "fasterq-dump"
    tool_version: str = ""
    compression: str = "none"
    files: List[Dict[str, Any]] = field(default_factory=list)
    reads_per_mate: Optional[int] = None
    bases_total: Optional[int] = None
    ncbi: Dict[str, Any] = field(default_factory=dict)
    completeness: Dict[str, Any] = field(default_factory=dict)
    stats: Dict[str, Any] = field(default_factory=dict)
    stats_computed: Optional[str] = None
    schema: int = SIDECAR_SCHEMA
    error: Optional[str] = None
    # Optional record of a re-download of this dataset; left out of the JSON while None, so a
    # sidecar written before the field existed is rewritten byte-identical (schema stays 1).
    refetch: Optional[Dict[str, Any]] = None
    # Optional description of the sourmash signature cached beside the dataset
    # (``prepare.store.sketches``); left out of the JSON while None, like ``refetch``.
    sketch: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Plain-dict form suitable for JSON serialization; ``refetch`` and ``sketch`` are left out while None."""
        data = asdict(self)
        for optional in ("refetch", "sketch"):
            if data.get(optional) is None:
                data.pop(optional, None)
        return data

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Sidecar":
        """Build a ``Sidecar`` from a plain dict, ignoring unknown keys and defaulting missing ones."""
        known_names = {f.name for f in dataclass_fields(cls)}
        filtered = {key: value for key, value in (data or {}).items() if key in known_names}
        return cls(**filtered)


def md5_file(path: Union[str, Path]) -> str:
    """MD5 hex digest of the file at ``path``, read in 1 MiB chunks so a large file is never
    loaded into memory at once.

    The one home for this: adoption's dedup check, ``store_verify`` and sidecar building all
    compare against the same digest, so they must compute it the same way."""
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def detect_layout(files: List[Path]) -> str:
    """``PAIRED`` when a mate-1 (``_1``/``_R1``) and a mate-2 (``_2``/``_R2``) file both exist,
    else ``SINGLE``."""
    stems = [fastq_stem(p) for p in files]
    has_mate1 = any(stem.endswith(MATE1_SUFFIXES) for stem in stems)
    has_mate2 = any(stem.endswith(_MATE2_SUFFIXES) for stem in stems)
    return "PAIRED" if has_mate1 and has_mate2 else "SINGLE"


def build_sidecar(
    accession: str,
    acc_dir: Union[str, Path],
    ncbi: Dict[str, Any],
    tool_version: str,
    compression: str,
    tool: str = "fasterq-dump",
    downloaded: Optional[str] = None,
    verify: bool = True,
) -> Sidecar:
    """Build a sidecar describing the files already downloaded for ``accession`` in ``acc_dir``.

    Reads every FASTQ file present once, through ``prepare.data.sra.fastq_digest``, which
    yields its size, md5 (of the bytes as stored) and read count together; the layout comes
    from the file names, and completeness reuses ``prepare.data.sra.verify_download`` against
    ``ncbi.get("spots")`` with those read counts, so no file is read a second time. A truncated
    gzip file raises ``EOFError`` (a corrupt or non-gzip one ``gzip.BadGzipFile``, an ``OSError``);
    either is caught per file and turns the whole result into ``state="failed"`` with the
    error recorded, since a corrupt file cannot be verified against NCBI's spot count. Such a
    file still records the md5 of its stored bytes, which takes a second read of that file only.

    ``tool`` and ``downloaded`` default to a fresh fasterq-dump download happening now;
    adoption passes ``tool="adopted"`` and the files' own age instead, since it did not
    download them. With ``verify`` off (``--no-verify-downloads``) the files are not judged
    against ``ncbi["spots"]``: the count stays on file as NCBI's figure, the verdict is
    ``unverified`` and the state ``complete``.
    """
    acc_path = Path(acc_dir)
    files = fastq_files(acc_path)
    layout = detect_layout(files)

    reads_by_path: Dict[Path, Optional[int]] = {}
    file_records = []
    error: Optional[str] = None
    for file_path in files:
        try:
            digest = fastq_digest(file_path)
            reads_by_path[file_path] = digest.records
            file_records.append(
                {"name": file_path.name, "bytes": digest.size, "md5": digest.md5, "reads": digest.records}
            )
        # A corrupt gzip stream raises zlib.error, which is not an OSError.
        except (EOFError, OSError, zlib.error) as exc:
            reads_by_path[file_path] = None
            error = f"{file_path.name}: {exc}"
            file_records.append(
                {"name": file_path.name, "bytes": file_path.stat().st_size, "md5": md5_file(file_path), "reads": None}
            )

    reads_per_mate: Optional[int]
    if error is not None:
        completeness = {"method": "unverified", "ratio": None, "verdict": "unverified"}
        reads_per_mate = None
        state = "failed"
    else:
        primary = primary_fastq(acc_path)
        orphan = orphan_fastq(acc_path)
        expected = ncbi.get("spots") if verify else None
        verdict = verify_download(
            accession,
            acc_path,
            expected,
            reads_r1=reads_by_path.get(primary, 0) if primary is not None else 0,
            reads_orphan=reads_by_path.get(orphan) if orphan is not None else None,
        )
        reads_per_mate = verdict["reads_r1"]
        method = "spots" if expected else "unverified"
        completeness = {"method": method, "ratio": verdict["ratio"], "verdict": verdict["verdict"]}
        state = "partial" if verdict["verdict"] == "truncated" else "complete"

    return Sidecar(
        accession=accession,
        state=state,
        layout=layout,
        downloaded=downloaded or datetime.now(timezone.utc).isoformat(),
        tool=tool,
        tool_version=tool_version,
        compression=compression,
        files=file_records,
        reads_per_mate=reads_per_mate,
        bases_total=None,
        ncbi=dict(ncbi),
        completeness=completeness,
        stats={},
        stats_computed=None,
        schema=SIDECAR_SCHEMA,
        error=error,
    )


def _to_int(value: Any) -> Optional[int]:
    """Best-effort int conversion; ``None`` for anything not numeric (including ``None`` itself)."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def ncbi_from_metadata_xml(xml_path: Union[str, Path]) -> Dict[str, Any]:
    """Map one NCBI metadata XML file's fields into the flat ``ncbi`` dict a sidecar records.

    ``spots``/``bases``/``size`` come from ``Run_Total_Spots``/``Run_Total_Bases``/``Run_Size``
    (int when present and numeric, else None); ``layout`` from ``Experiment_Library_Layout``;
    ``files`` holds one ``{"name", "md5"}`` entry built from ``Run_Filename``/``Run_MD5`` when
    either is present, else an empty list. Returns ``{}`` when the XML file is missing or
    unparsable, logging a warning here since ``parse_metadata_xml`` itself only raises.
    """
    try:
        parsed = parse_metadata_xml(xml_path)
    except (OSError, ValueError, ET.ParseError) as e:
        logger.warning(f"Could not read NCBI metadata from {xml_path}: {e}")
        return {}
    if not parsed:
        return {}

    files: List[Dict[str, Any]] = []
    name = parsed.get("Run_Filename")
    md5 = parsed.get("Run_MD5")
    if name is not None or md5 is not None:
        files.append({"name": name, "md5": md5})

    return {
        "spots": _to_int(parsed.get("Run_Total_Spots")),
        "bases": _to_int(parsed.get("Run_Total_Bases")),
        "size": _to_int(parsed.get("Run_Size")),
        "layout": parsed.get("Experiment_Library_Layout"),
        "files": files,
    }


def write_sidecar(path: Union[str, Path], sidecar: Sidecar) -> Path:
    """Write ``sidecar`` to ``path`` atomically and flushed to disk (``write_text_atomic``), sorted keys, indent 2."""
    target = Path(path)
    try:
        write_text_atomic(target, json.dumps(sidecar.to_dict(), indent=2, sort_keys=True) + "\n", fsync=True)
    except OSError as e:
        raise DataAccessError(f"Cannot write sidecar {target}: {e}") from e
    return target


def read_sidecar(path: Union[str, Path]) -> Optional[Sidecar]:
    """Read the sidecar at ``path``, returning None with a logged warning when missing or unparseable.

    Raises ``DataAccessError`` when the file parses but holds something other than a JSON object,
    so the caller reports a store error rather than failing with an ``AttributeError``.
    """
    sidecar_path = Path(path)
    if not sidecar_path.is_file():
        logger.warning(f"Sidecar file not found: {sidecar_path}")
        return None
    try:
        data = json.loads(sidecar_path.read_text())
    except (OSError, json.JSONDecodeError) as e:
        logger.warning(f"Could not read sidecar {sidecar_path}: {e}")
        return None
    if not isinstance(data, dict):
        raise DataAccessError(f"{sidecar_path}: sidecar is not a JSON object")
    return Sidecar.from_dict(data)


def sidecar_completeness(path: Union[str, Path]) -> Optional[Dict[str, Any]]:
    """The completeness verdict a store sidecar records, shaped for a registry download record.

    Returns None when there is no sidecar to read (``read_sidecar`` logs the reason). Shared by
    every command that points a project's download record at the store's copy: ``download_sra``,
    ``store_link`` and ``store_adopt`` must all record the same verdict for the same dataset.
    """
    sidecar = read_sidecar(path)
    if sidecar is None:
        return None
    return {
        "verdict": sidecar.completeness.get("verdict"),
        "ratio": sidecar.completeness.get("ratio"),
        "expected_spots": sidecar.ncbi.get("spots"),
        "reads_r1": sidecar.reads_per_mate,
    }
