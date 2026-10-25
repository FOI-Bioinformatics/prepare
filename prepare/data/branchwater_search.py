"""Query the Branchwater index of SRA metagenomes with a genome sketch.

The public Branchwater service (https://branchwater.sourmash.bio) indexes about
1.1 million SRA metagenomes as FracMinHash sketches (k=21, scaled=1000). Its
search API takes one sourmash signature and a minimum containment and returns
the matching SRA accessions with their containment. This module builds the
sketch, runs the search and writes the result in the CSV layout the rest of
PREPARE reads (see ``use_branchwater``). The metadata columns are left empty;
``download_metadata`` fills them from NCBI.
"""

import csv
import gzip
import hashlib
import io
import json
import logging
import zlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from prepare.core.exceptions import DataAccessError
from prepare.core.validation import genome_name
from prepare.core.optional import require
from prepare.data.file_io import open_atomic, write_text_atomic

logger = logging.getLogger(__name__)

DEFAULT_SERVER = "https://api.branchwater.sourmash.bio"
KSIZE = 21
SCALED = 1000
RETRY_TOTAL = 4
RETRY_BACKOFF_FACTOR = 2
RETRY_STATUS_FORCELIST = [429, 500, 502, 503, 504]
BRANCHWATER_COLUMNS = [
    "acc",
    "containment",
    "cANI",
    "biosample",
    "bioproject",
    "assay_type",
    "collection_date_sam",
    "geo_loc_name_country_calc",
    "organism",
    "lat_lon",
]


def _open_text(path: Path):
    """Open a plain or gzipped (by magic number) text file for reading."""
    with open(path, "rb") as probe:
        packed = probe.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if packed else open(path, "r")


def sketch_fasta(fasta_path: Union[str, Path]) -> Dict[str, Any]:
    """Build the k=21, scaled=1000 sourmash signature Branchwater expects for one FASTA file.

    Returns the signature as the JSON object sourmash writes (one element of a
    ``.sig`` file), ready to be posted to the search API.

    Raises:
        ConfigurationError: If sourmash is not installed.
        DataAccessError: If the file is missing, cannot be read (a truncated or corrupt
            gzip among the reasons) or holds no sequences.
    """
    sourmash = require("sourmash", "sourmash", "Sketching a genome")
    signature_module = require("sourmash.signature", "sourmash", "Sketching a genome")
    from Bio import SeqIO

    path = Path(fasta_path)
    if not path.exists():
        raise DataAccessError(f"Genome FASTA not found: {path}")

    minhash = sourmash.MinHash(n=0, ksize=KSIZE, scaled=SCALED)
    n_records = 0
    try:
        with _open_text(path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                minhash.add_sequence(str(record.seq).upper(), force=True)
                n_records += 1
    except (OSError, EOFError, zlib.error) as e:  # gzip.BadGzipFile is an OSError
        raise DataAccessError(f"Cannot read {path}: {e}") from e
    if n_records == 0:
        raise DataAccessError(f"No sequences found in {path}")

    signature = sourmash.SourmashSignature(minhash, name=genome_name(path), filename=path.name)
    buffer = io.BytesIO()
    signature_module.save_signatures_to_json([signature], buffer)
    payload = json.loads(buffer.getvalue().decode("utf-8"))
    logger.info(
        "Sketched %s: %d sequence(s), %d hashes (k=%d, scaled=%d)",
        path.name,
        n_records,
        len(minhash.hashes),
        KSIZE,
        SCALED,
    )
    return payload[0]


def load_signature(sig_path: Union[str, Path]) -> Dict[str, Any]:
    """Read one sourmash JSON signature (plain or gzipped) and select its k=21, scaled=1000 sketch.

    A signature holding several sketches (for example one file with k=21 and k=31
    sketches) is common when a .sig was built with multiple k-mer sizes; only the
    k=21, scaled=1000 sketch is what Branchwater expects, so the returned object's
    ``signatures`` list is narrowed to that one sketch.
    """
    path = Path(sig_path)
    if not path.exists():
        raise DataAccessError(f"Signature file not found: {path}")
    try:
        with open(path, "rb") as handle:
            raw = handle.read()
        if raw[:2] == b"\x1f\x8b":
            # A .sig.gz, as sourmash's own -o x.sig.gz and sra_sketch write it.
            raw = gzip.decompress(raw)
        data = json.loads(raw.decode("utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError, OSError, EOFError, zlib.error) as e:
        raise DataAccessError(f"Signature file is not valid JSON: {path} ({e})") from e

    signature = data[0] if isinstance(data, list) and data else data
    if not isinstance(signature, dict) or not signature.get("signatures"):
        raise DataAccessError(f"Not a sourmash signature file: {path}")
    sketch = _select_sketch(signature, path)
    return {**signature, "signatures": [sketch]}


def _select_sketch(signature: Dict[str, Any], path: Path) -> Dict[str, Any]:
    """Return the k=21, scaled=1000 sketch from a signature's (possibly multi-sketch) list."""
    for sketch in signature["signatures"]:
        max_hash = sketch.get("max_hash") or 0
        scaled = round(2**64 / max_hash) if max_hash else None
        if sketch.get("ksize") == KSIZE and scaled == SCALED:
            return sketch

    sketch = signature["signatures"][0]
    ksize = sketch.get("ksize")
    max_hash = sketch.get("max_hash") or 0
    scaled = round(2**64 / max_hash) if max_hash else None
    raise DataAccessError(f"Branchwater needs k={KSIZE}, scaled={SCALED}; {path.name} has k={ksize}, scaled={scaled}")


def _session() -> requests.Session:
    """Build a requests.Session that retries transient Branchwater failures with backoff.

    A 429 or 5xx response is retried up to ``RETRY_TOTAL`` times with exponential
    backoff; any other status (a 4xx client error, for example) is returned as-is
    and left for the caller to raise on.
    """
    session = requests.Session()
    retry = Retry(
        total=RETRY_TOTAL,
        backoff_factor=RETRY_BACKOFF_FACTOR,
        status_forcelist=RETRY_STATUS_FORCELIST,
        allowed_methods=["POST"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def _sketch_mins(signature: Dict[str, Any]) -> List[int]:
    """Return the sorted mins of signature's k=21, scaled=1000 sketch, used as the cache key input."""
    sketch = _select_sketch(signature, Path("signature"))
    return sorted(sketch.get("mins", []))


def _cache_key(signature: Dict[str, Any], threshold: float, server: str) -> str:
    """Derive a stable cache key from the sketch content and the search parameters."""
    payload = {
        "mins": _sketch_mins(signature),
        "ksize": KSIZE,
        "scaled": SCALED,
        "threshold": threshold,
        "server": server,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def _cache_read(cache_dir: Path, key: str, max_cache_age_days: Optional[int]) -> Optional[str]:
    """Return the cached CSV text for ``key``, or None if there is no usable cache entry."""
    csv_path = cache_dir / f"{key}.csv"
    meta_path = cache_dir / f"{key}.json"
    if not csv_path.exists() or not meta_path.exists():
        return None
    try:
        meta = json.loads(meta_path.read_text())
        fetched = datetime.fromisoformat(meta["fetched"])
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return None
    if max_cache_age_days is not None:
        age_days = (datetime.now(timezone.utc) - fetched).total_seconds() / 86400.0
        if age_days > max_cache_age_days:
            return None
    try:
        text = csv_path.read_text()
    except OSError:
        return None
    logger.info("using cached Branchwater result from %s", fetched.date().isoformat())
    return text


def _cache_write(cache_dir: Path, key: str, text: str, server: str, threshold: float, rows: int) -> None:
    """Write a cache entry; a failure (for example a read-only cache directory) is logged, not raised."""
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        write_text_atomic(cache_dir / f"{key}.csv", text)
        meta = {
            "fetched": datetime.now(timezone.utc).isoformat(),
            "server": server,
            "threshold": threshold,
            "rows": rows,
        }
        write_text_atomic(cache_dir / f"{key}.json", json.dumps(meta))
    except OSError as e:
        logger.warning("Could not write Branchwater cache entry %s: %s", key, e)


def search_index(
    signature: Dict[str, Any],
    threshold: float,
    server: str = DEFAULT_SERVER,
    timeout: int = 600,
    cache_dir: Optional[Union[str, Path]] = None,
    refresh: bool = False,
    max_cache_age_days: Optional[int] = None,
) -> List[Tuple[str, float]]:
    """Post a signature to the Branchwater search API and return (accession, containment) pairs, best first.

    When ``cache_dir`` is given, a prior result for the same sketch, threshold and
    server is reused instead of querying the server again, unless ``refresh`` is set
    or the cached entry is older than ``max_cache_age_days``.
    """
    cache_path = Path(cache_dir) if cache_dir is not None else None
    cache_key = _cache_key(signature, threshold, server) if cache_path is not None else None

    text: Optional[str] = None
    if cache_path is not None and not refresh:
        assert cache_key is not None  # set together with cache_path above
        text = _cache_read(cache_path, cache_key, max_cache_age_days)

    fetched_now = False
    if text is None:
        url = f"{server.rstrip('/')}/search"
        session = _session()
        try:
            response = session.post(
                url, json={"threshold": threshold, "signature": signature}, timeout=timeout, stream=True
            )
        except requests.exceptions.RequestException as e:
            raise DataAccessError(f"Branchwater search failed: {e}") from e
        if response.status_code != 200:
            raise DataAccessError(f"Branchwater search returned HTTP {response.status_code}: {response.text[:200]}")
        # A response with no charset in its Content-Type (or a non-text content type) leaves
        # encoding unset, and iter_lines then yields bytes, which csv.reader rejects with an
        # error that is not a PrepareError. The server sends UTF-8 CSV either way.
        response.encoding = response.encoding or "utf-8"
        # decode_unicode=True yields str given the encoding set above, but iter_lines is typed
        # to allow bytes when encoding cannot be determined; decode defensively so "\n".join
        # below and _parse_search_rows always see str.
        lines: List[str] = [
            line if isinstance(line, str) else line.decode(response.encoding or "utf-8")
            for line in response.iter_lines(decode_unicode=True)
        ]
        matches = _parse_search_rows(lines)
        text = "\n".join(lines)
        fetched_now = True
    else:
        matches = _parse_search_rows(text.splitlines())

    kept = [match for match in matches if match[1] >= threshold]
    n_dropped = len(matches) - len(kept)
    if n_dropped:
        logger.info("Dropped %d match(es) below containment %.2f reported by the server", n_dropped, threshold)

    if fetched_now and cache_path is not None:
        assert cache_key is not None  # set together with cache_path above
        _cache_fresh_result(cache_path, cache_key, text, server, threshold, len(kept), len(matches))

    return kept


def _cache_fresh_result(
    cache_path: Path, cache_key: str, text: str, server: str, threshold: float, n_kept: int, n_rows: int
) -> None:
    """Cache a fetched reply, except a reply with no rows at all.

    An empty reply is often an outage or a transient server state; caching it would
    hide later results for the lifetime of the cache entry.
    """
    if n_rows == 0:
        logger.info("Branchwater returned no rows: empty result not cached")
        return
    _cache_write(cache_path, cache_key, text, server, threshold, n_kept)


def _parse_search_rows(lines) -> List[Tuple[str, float]]:
    """Parse Branchwater search result lines (a header row plus data rows) into (accession, containment) pairs."""
    reader = csv.reader(lines)
    header = next(reader, [])
    if "SRA accession" not in header or "containment" not in header:
        raise DataAccessError(f"Unexpected Branchwater response header: {header}")
    acc_idx = header.index("SRA accession")
    cont_idx = header.index("containment")
    matches: List[Tuple[str, float]] = []
    for row in reader:
        if not row:
            continue
        try:
            matches.append((row[acc_idx].strip(), float(row[cont_idx])))
        except (IndexError, ValueError, AttributeError):
            logger.warning("Skipping malformed Branchwater row: %s", row)
    matches.sort(key=lambda match: match[1], reverse=True)
    return matches


def _parse_search_csv(text: str) -> List[Tuple[str, float]]:
    """Parse a full Branchwater CSV response body (kept for a cached result read back as one string)."""
    return _parse_search_rows(text.splitlines())


def write_branchwater_csv(matches: List[Tuple[str, float]], output_path: Union[str, Path], ksize: int = KSIZE) -> Path:
    """Write matches in Branchwater CSV layout; cANI is derived from containment, metadata columns stay empty."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    empty_metadata = [""] * (len(BRANCHWATER_COLUMNS) - 3)
    with open_atomic(path, newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(BRANCHWATER_COLUMNS)
        for accession, containment in matches:
            cani = round(containment ** (1 / ksize), 4) if containment > 0 else 0.0
            writer.writerow([accession, f"{containment:.4f}", f"{cani:.4f}"] + empty_metadata)
    logger.info("Wrote %d match(es) to %s", len(matches), path)
    return path
