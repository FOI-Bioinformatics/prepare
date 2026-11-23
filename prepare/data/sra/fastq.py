"""FASTQ files of a downloaded SRA accession: discovery, read counting, compression and verification.

Every function here reads or rewrites files already on disk; none of them starts a download.
"""

import gzip
import hashlib
import json
import logging
import re
import subprocess
import threading
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from prepare.core.constants import FASTQ_GLOBS
from prepare.core.exceptions import SecurityError
from prepare.data.file_io import open_atomic, visible_files
from prepare.utils.security import SecureSubprocess
from prepare.utils.tools import usable_optional_tool

logger = logging.getLogger(__name__)


# What SecureSubprocess.run_secure raises for a tool that fails, times out, is refused or
# cannot be started.
_TOOL_ERRORS = (subprocess.CalledProcessError, subprocess.TimeoutExpired, SecurityError, OSError)


# Ratio of downloaded reads to NCBI's recorded run_total_spots at or above which a download
# counts as complete rather than truncated.
COMPLETE_RATIO_THRESHOLD = 0.99

# FASTQ file-name suffixes marking a mate of a pair: fasterq-dump's --split-3 output uses
# _1/_2, while files from other sources often use _R1/_R2.
MATE1_SUFFIXES = ("_1", "_R1")
MATE2_SUFFIXES = ("_2", "_R2")
MATE_SUFFIXES = ("_1", "_2", "_R1", "_R2")


def fastq_files(acc_dir: Union[str, Path]) -> List[Path]:
    """Non-empty, visible FASTQ files directly in ``acc_dir``, sorted by name (see ``visible_files``).

    Matches ``FASTQ_GLOBS`` (plain and gzipped ``.fastq``/``.fq``). A directory that does not
    exist (or a dangling symlink) yields an empty list; a symlinked directory is followed
    since ``Path.is_dir``/``Path.glob`` already resolve it transparently. A zero-byte file
    (e.g. left behind by an interrupted download) is never returned, and neither is a hidden
    name such as the ``._<name>`` AppleDouble files macOS writes next to every file on a
    volume without native extended attributes.
    """
    found = []
    for candidate in visible_files(acc_dir, *FASTQ_GLOBS):
        try:
            if candidate.stat().st_size > 0:
                found.append(candidate)
        except OSError:
            continue
    return found


def fastq_stem(path: Union[str, Path]) -> str:
    """The file name of ``path`` without its FASTQ extension (``.fastq``/``.fq``, plain or gzipped)."""
    name = Path(path).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def primary_fastq(acc_dir: Union[str, Path]) -> Optional[Path]:
    """The FASTQ file in ``acc_dir`` holding one record per sequenced spot, or None if there is none.

    A mate-1 file (``_1``/``_R1``) is preferred, since for paired data it holds exactly one
    record per spot; otherwise the bare ``<acc>.fastq`` file written for single-end data.
    Sorting by name is not enough here: ``<acc>.fastq`` sorts before ``<acc>_1.fastq``, and
    for paired data that bare file holds only the unpaired leftovers of ``--split-3``.
    """
    files = fastq_files(acc_dir)
    if not files:
        return None
    mate_one = next((p for p in files if fastq_stem(p).endswith(MATE1_SUFFIXES)), None)
    if mate_one is not None:
        return mate_one
    bare = next((p for p in files if not fastq_stem(p).endswith(MATE_SUFFIXES)), None)
    return bare if bare is not None else files[0]


def mate2_fastq(acc_dir: Union[str, Path]) -> Optional[Path]:
    """The mate-2 file (``_2``/``_R2``) in ``acc_dir``, or None when there is none."""
    return next((p for p in fastq_files(acc_dir) if fastq_stem(p).endswith(MATE2_SUFFIXES)), None)


def orphan_fastq(acc_dir: Union[str, Path]) -> Optional[Path]:
    """The bare ``<acc>.fastq`` file of unpaired spots ``--split-3`` writes beside a mate pair.

    Returns None when the folder holds no mate-1 file, since then the bare file is the
    single-end data itself rather than a set of leftovers (``primary_fastq`` returns it).
    """
    files = fastq_files(acc_dir)
    if not any(fastq_stem(p).endswith(MATE1_SUFFIXES) for p in files):
        return None
    return next((p for p in files if not fastq_stem(p).endswith(MATE_SUFFIXES)), None)


def accession_has_fastq(acc_dir: Union[str, Path]) -> bool:
    """Return True if the per-accession directory holds at least one usable FASTQ file.

    This is the single source of truth for "this accession is already
    downloaded" used across the download and status paths. A sidecar
    ``<acc_dir>/<acc_dir.name>.json`` recording an in-progress or failed state
    (``partial``, ``failed`` or ``downloading``) overrides an otherwise
    present-looking directory, so a half-written or restarted download is not
    mistaken for a finished one.
    """
    acc_path = Path(acc_dir)
    if not fastq_files(acc_path):
        return False

    sidecar = acc_path / f"{acc_path.name}.json"
    if sidecar.exists():
        try:
            data = json.loads(sidecar.read_text())
        except (OSError, ValueError):
            data = None
        # A sidecar that is unreadable, not JSON, or not a JSON object carries no state.
        state = data.get("state") if isinstance(data, dict) else None
        if state in ("partial", "failed", "downloading"):
            return False

    return True


# Block size for streaming reads of FASTQ files (``count_fastq_reads`` and ``fastq_digest``).
_BLOCK_SIZE = 1024 * 1024


def count_fastq_reads(path: Union[str, Path]) -> int:
    """Count FASTQ records in ``path`` via a chunked binary newline count (4 lines/record).

    Reads in 1 MiB blocks so a large FASTQ file is never loaded into memory, and counts
    raw ``b"\\n"`` bytes rather than decoding text, which is both faster and immune to
    encoding errors in a corrupted file. Works for both gzip-compressed and plain files.

    A truncated gzip raises ``EOFError``; a corrupt deflate stream raises ``gzip.BadGzipFile``
    (an ``OSError``), as ``_GzipInflater`` does, rather than the bare ``zlib.error`` that
    ``gzip`` lets through and that no caller expects.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    total_newlines = 0
    last_byte = b""
    with opener(path, "rb") as handle:
        while True:
            try:
                block = handle.read(_BLOCK_SIZE)
            except zlib.error as exc:
                raise gzip.BadGzipFile(f"{path}: {exc}") from exc
            if not block:
                break
            total_newlines += block.count(b"\n")
            last_byte = block[-1:]
    # A file whose last line has no trailing newline still ends a record; count it too.
    if last_byte and last_byte != b"\n":
        total_newlines += 1
    return total_newlines // 4


# Buffer size of the file handle under ``fastq_digest``'s block reads; large sequential reads
# suit an external or network volume.
_READ_BUFFER = 8 * 1024 * 1024
# Upper bound on decompressed output per decompress call, so one highly compressible block
# never expands into an unbounded buffer.
_MAX_INFLATE = 4 * 1024 * 1024


@dataclass(frozen=True)
class FastqDigest:
    """What one pass over a stored FASTQ file yields.

    ``records`` follows ``count_fastq_reads``; ``md5`` is the digest of the file as stored (the
    compressed bytes for a ``.gz`` file), matching ``prepare.store.sidecar.md5_file``; ``size``
    is the number of bytes read. ``content_md5`` is the digest of the decompressed content, and
    is only computed when asked for.
    """

    records: int
    md5: str
    size: int
    content_md5: Optional[str] = None


class _NewlineCounter:
    """Counts newlines in a byte stream fed in pieces and remembers its last byte."""

    def __init__(self, content_digest: Optional[Any]) -> None:
        self.newlines = 0
        self.last_byte = b""
        self._content_digest = content_digest

    def feed(self, data: bytes) -> None:
        """Account for one piece of (decompressed) content."""
        if not data:
            return
        self.newlines += data.count(b"\n")
        self.last_byte = data[-1:]
        if self._content_digest is not None:
            self._content_digest.update(data)

    def content_md5(self) -> Optional[str]:
        """Hex md5 of the content fed so far, or None when it was not asked for."""
        return self._content_digest.hexdigest() if self._content_digest is not None else None

    def records(self) -> int:
        """Records under the rule ``count_fastq_reads`` uses: an unterminated last line still counts."""
        newlines = self.newlines
        if self.last_byte and self.last_byte != b"\n":
            newlines += 1
        return newlines // 4


class _GzipInflater:
    """Inflates a gzip stream fed in pieces, including files of several concatenated members.

    Errors are raised as ``gzip.open`` raises them: ``gzip.BadGzipFile`` (an ``OSError``) for a
    member that does not start with the gzip magic bytes, with gzip's own wording, or for a
    corrupt deflate stream, and ``EOFError`` from ``finish`` for a stream cut short.
    """

    def __init__(self, sink: _NewlineCounter) -> None:
        self._sink = sink
        self._inflater = zlib.decompressobj(16 + zlib.MAX_WBITS)
        self._in_member = False

    def feed(self, data: bytes) -> None:
        """Inflate ``data``, starting a new member whenever the previous one has ended."""
        try:
            self._feed(data)
        except zlib.error as exc:
            raise gzip.BadGzipFile(str(exc)) from exc

    def _feed(self, data: bytes) -> None:
        """``feed`` without the translation of ``zlib.error``."""
        while data:
            if self._inflater.eof:
                # Zero bytes after a member are padding, as gzip.open also accepts.
                data = data.lstrip(b"\x00")
                if not data:
                    return
                self._inflater = zlib.decompressobj(16 + zlib.MAX_WBITS)
            if not self._in_member and len(data) >= 2 and data[:2] != b"\x1f\x8b":
                raise gzip.BadGzipFile(f"Not a gzipped file ({data[:2]!r})")
            self._in_member = True
            self._sink.feed(self._inflater.decompress(data, _MAX_INFLATE))
            while self._inflater.unconsumed_tail and not self._inflater.eof:
                self._sink.feed(self._inflater.decompress(self._inflater.unconsumed_tail, _MAX_INFLATE))
            if self._inflater.eof:
                self._in_member = False
                data = self._inflater.unused_data
            else:
                data = b""

    def finish(self) -> None:
        """Raise ``EOFError`` when the stream stopped inside a member, as ``gzip.open`` does."""
        self._sink.feed(self._inflater.flush())
        if self._in_member and not self._inflater.eof:
            raise EOFError("Compressed file ended before the end-of-stream marker was reached")


def fastq_digest(path: Union[str, Path], content_md5: bool = False) -> FastqDigest:
    """Record count, md5 and size of the FASTQ file at ``path`` from a single read of its bytes.

    The stored bytes are read once in 1 MiB blocks; each block updates the md5 and, for a
    ``.gz`` file, is inflated on the fly (several concatenated gzip members are handled) so
    newlines are counted in the decompressed content. A plain file counts newlines in the
    blocks directly. The result equals ``count_fastq_reads(path)``, ``md5_file(path)`` and the
    file size, at the cost of one pass instead of two. A truncated gzip stream raises
    ``EOFError``, and a corrupt one or a file that is not gzip ``gzip.BadGzipFile``, as
    ``count_fastq_reads`` would for the latter.

    With ``content_md5`` the md5 of the decompressed content is also computed, so a plain and
    a gzipped copy of the same reads can be compared without a second pass.
    """
    raw_digest = hashlib.md5()
    counter = _NewlineCounter(hashlib.md5() if content_md5 else None)
    inflater = _GzipInflater(counter) if str(path).endswith(".gz") else None
    size = 0
    with open(path, "rb", buffering=_READ_BUFFER) as handle:
        while True:
            block = handle.read(_BLOCK_SIZE)
            if not block:
                break
            size += len(block)
            raw_digest.update(block)
            if inflater is not None:
                inflater.feed(block)
            else:
                counter.feed(block)
    if inflater is not None:
        inflater.finish()
    return FastqDigest(
        records=counter.records(), md5=raw_digest.hexdigest(), size=size, content_md5=counter.content_md5()
    )


def iter_fastq_records(path: Union[str, Path]):
    """Yield ``(sequence, quality)`` string pairs for each record in ``path``, streaming.

    Works for both gzip-compressed and plain files. A raw four-line reader for a caller that
    needs every record without Biopython's slower per-record parser; a caller that needs only
    a sample of records uses ``prepare.data.sra.sampling.sample_records`` instead, which
    does no per-record work for the records it skips.

    Raises ``ValueError`` when a header line is not followed by a complete
    sequence/plus/quality triplet, since a truncated trailing record cannot be trusted.
    Truncation is decided on the raw lines: a line that is missing entirely (``readline``
    returns the empty string at end of file) ends the record early, whereas a present but
    empty line does not. A record whose sequence and quality lines are both empty is a
    zero-length read, which is legal FASTQ (a read trimmed to nothing, or the empty mate
    ``fasterq-dump --split-files`` writes for a half-empty spot) and is yielded as such.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq_line = handle.readline()
            plus_line = handle.readline()
            qual_line = handle.readline()
            if not seq_line or not plus_line or not qual_line:
                raise ValueError(f"Truncated FASTQ record after header: {header.strip()!r}")
            if not plus_line.startswith("+"):
                raise ValueError(f"Malformed FASTQ record after header: {header.strip()!r}")
            yield seq_line.rstrip("\r\n"), qual_line.rstrip("\r\n")


def verify_download(
    accession: str,
    acc_dir: Union[str, Path],
    expected_spots: Optional[int],
    reads_r1: Optional[int] = None,
    reads_orphan: Optional[int] = None,
    reads_r2: Optional[int] = None,
) -> Dict[str, Any]:
    """Compare what actually downloaded for ``accession`` against NCBI's recorded spot count.

    ``reads_r1`` counts one record per sequenced spot: the mate-1 file (or the single-end
    file) plus the bare ``<acc>.fastq`` file of unpaired spots that ``--split-3`` writes
    beside a mate pair, since NCBI's ``total_spots`` covers those too. The verdict is
    ``"complete"`` when the ratio of downloaded reads to ``expected_spots`` is at least
    ``COMPLETE_RATIO_THRESHOLD``, ``"truncated"`` below that, and ``"unverified"`` when
    ``expected_spots`` is unknown (e.g. NCBI metadata was never fetched for this accession).

    For paired data the mate-2 file must hold as many records as the mate-1 file (the orphan
    file aside). When it does not, ``mates_match`` is False, ``reads_r1`` counts only the spots
    with both mates (the smaller mate count, plus the orphan file) and is judged against
    ``expected_spots`` as usual, and without a spot count the verdict is ``"truncated"``. This
    reads the mate-2 file once unless ``reads_r2`` is given. A mate-1 file with no mate-2 file
    beside it is not judged on mates, since some runs yield only ``_1``.

    A caller that has already counted the records of the primary file (``reads_r1``), the
    orphan file (``reads_orphan``) or the mate-2 file (``reads_r2``) passes them in, and that
    file is then not read again.
    """
    files = fastq_files(acc_dir)
    primary = primary_fastq(acc_dir)
    orphan = orphan_fastq(acc_dir)
    if reads_r1 is None:
        reads_r1 = _count_named(primary) if primary is not None else 0
    mate1_reads = reads_r1
    if orphan is not None:
        reads_r1 += reads_orphan if reads_orphan is not None else _count_named(orphan)
    bytes_total = sum(p.stat().st_size for p in files)

    reads_r2 = _mate2_reads(acc_dir, primary, orphan, reads_r2)
    return {**_judged(reads_r1, mate1_reads, reads_r2, expected_spots), "bytes_total": bytes_total}


def _judged(reads_r1: int, mate1_reads: int, reads_r2: Optional[int], expected_spots: Optional[int]) -> Dict[str, Any]:
    """The verdict fields of ``verify_download`` for these counts (see there)."""
    mates_match = reads_r2 is None or reads_r2 == mate1_reads
    if reads_r2 is not None and not mates_match:
        # Only the spots with both mates are complete: every later recompute from reads_r1 (a
        # registry merge, a reverify) then reaches the same verdict.
        reads_r1 -= mate1_reads - min(mate1_reads, reads_r2)

    # Imported here: the spots module imports COMPLETE_RATIO_THRESHOLD from this one.
    from prepare.data.sra.spots import verdict_for_count

    verdict = verdict_for_count(reads_r1, expected_spots)
    if not mates_match and verdict["verdict"] == "unverified":
        verdict["verdict"] = "truncated"
    return {
        "reads_r1": reads_r1,
        "reads_mate1": mate1_reads,
        "reads_r2": reads_r2,
        "mates_match": mates_match,
        "expected_spots": expected_spots,
        "ratio": verdict["ratio"],
        "verdict": verdict["verdict"],
    }


def _count_named(path: Path) -> int:
    """``count_fastq_reads(path)``, with the file name added to the message of an unreadable file."""
    try:
        return count_fastq_reads(path)
    except (EOFError, OSError) as e:
        raise type(e)(f"{path.name}: {e}") from e


def _mate2_reads(
    acc_dir: Union[str, Path], primary: Optional[Path], orphan: Optional[Path], given: Optional[int]
) -> Optional[int]:
    """The read count of the mate-2 file beside a mate-1 file (``given`` when known), or None when there is none."""
    paired = orphan is not None or (primary is not None and fastq_stem(primary).endswith(MATE1_SUFFIXES))
    mate2 = mate2_fastq(acc_dir) if paired else None
    if mate2 is None:
        return None
    return given if given is not None else _count_named(mate2)


_VERDICT_MESSAGE_RE = re.compile(r", (complete|truncated) \((\d+) of (\d+) spots\)")
_MATES_MESSAGE_RE = re.compile(r", truncated \(mate files differ: (\d+) and (\d+) reads\)")


def mates_message(reads_r1: int, reads_r2: int) -> str:
    """The verdict text for a download whose mate files differ, in the form ``parse_verdict_message`` reads."""
    return f"truncated (mate files differ: {reads_r1} and {reads_r2} reads)"


def parse_verdict_message(message: str) -> Optional[Dict[str, Any]]:
    """Recover the verdict dict encoded in a download result message by ``_handle_download_output``.

    The ``on_result`` callback contract only carries a plain message string across the
    download/registry boundary, so the CLI recovers the verdict from it rather than the data
    layer reaching into the registry directly. Returns ``None`` for a message that carries no
    verdict (a failure message, or "already exists").

    Both patterns are anchored on the ", <verdict>" separator ``_handle_download_output``
    writes, so a failure message that merely contains one of these words is not read as a
    verdict.
    """
    mates = _MATES_MESSAGE_RE.search(message)
    if mates:
        reads_r1, reads_r2 = int(mates.group(1)), int(mates.group(2))
        # The complete spots are those with both mates: the smaller of the two counts.
        return {
            "verdict": "truncated",
            "reads_r1": min(reads_r1, reads_r2),
            "reads_r2": reads_r2,
            "expected_spots": None,
            "ratio": None,
        }
    match = _VERDICT_MESSAGE_RE.search(message)
    if match:
        verdict, reads_r1, expected_spots = match.group(1), int(match.group(2)), int(match.group(3))
        ratio = round(reads_r1 / expected_spots, 4) if expected_spots else 0.0
        return {"verdict": verdict, "reads_r1": reads_r1, "expected_spots": expected_spots, "ratio": ratio}
    if ", unverified" in message:
        return {"verdict": "unverified"}
    return None


def compress_fastq(path: Path, threads: int, stop: Optional[threading.Event] = None) -> Path:
    """Gzip-compress ``path`` in place, returning the path to the compressed file.

    Uses ``pigz`` (parallel gzip) when it is on PATH, since it is substantially faster
    than single-threaded gzip on the multi-core machines this tool typically runs on;
    otherwise falls back to Python's ``gzip`` module, streaming the file through in
    1 MiB blocks so a large FASTQ file is never fully loaded into memory. Either way
    the uncompressed source is removed and only the ``.gz`` file remains.

    On failure, no partial ``.gz`` is left behind and the uncompressed source survives
    untouched: the Python fallback writes to a unique hidden temp file and only
    renames it onto the final ``.gz`` name (and only then unlinks the source)
    once the gzip stream has closed cleanly; a failing ``pigz`` invocation has any
    ``.gz`` it managed to write before dying removed. Either way the original
    exception propagates to the caller.

    ``stop`` is the calling download run's stop token; ``pigz`` is recorded under it so
    ``SecureSubprocess.terminate_children(stop=stop)`` stops it with the rest of that run.
    The Python fallback checks it before every block and raises ``InterruptedError`` (an
    ``OSError``) once it is set, removing its temporary file and keeping the source.
    """
    target = path.with_suffix(path.suffix + ".gz")

    # An optional pigz older than its floor (or not runnable) is skipped for Python's gzip.
    if usable_optional_tool("pigz"):
        try:
            SecureSubprocess.run_secure("pigz", ["-p", str(threads), "-f", str(path)], stop=stop)
        except _TOOL_ERRORS:
            if target.exists():
                try:
                    target.unlink()
                except OSError as cleanup_error:
                    logger.warning(f"Could not remove partial {target} after a failed pigz run: {cleanup_error}")
            raise
        return target

    block_size = 1024 * 1024
    # open_atomic writes a dot-prefixed temporary name (unique_temp_path), which fastq_files and
    # every other folder listing skip, and removes it on any exit that did not publish it. The
    # gzip header records the real file name rather than the temporary one.
    with open(path, "rb") as source, open_atomic(target, "wb") as raw:
        with gzip.GzipFile(filename=target.name, mode="wb", fileobj=raw, compresslevel=6) as dest:
            while True:
                if stop is not None and stop.is_set():
                    raise InterruptedError(f"Compression of {path} stopped: the run was interrupted")
                block = source.read(block_size)
                if not block:
                    break
                dest.write(block)
    path.unlink()
    return target
