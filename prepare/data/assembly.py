"""Assembly of extracted reads with megahit, and summaries of the resulting contigs.

The reads ``prepare.data.read_extraction`` keeps for one target genome are assembled here into
a small, targeted assembly. megahit runs through ``SecureSubprocess.run_secure`` and must be
present on the system. This module does not import ``read_extraction``; that module re-exports
the names defined here so existing imports keep working.
"""

import logging
import platform
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple, Union

from prepare.core.exceptions import ProcessingError
from prepare.data.assembly_identity import (
    CONTIGS_NAME,
    MARKER_NAME,
    AssemblyInputs,
    assembly_state,
    publish_assembly,
    write_marker,
)
from prepare.data.file_io import unique_temp_path
from prepare.utils import resources, tools
from prepare.utils.security import SecureSubprocess

logger = logging.getLogger(__name__)


def assembly_memory(value: str) -> Optional[Union[int, float]]:
    """The megahit ``--memory`` value for an ``--assembly-memory`` setting, or None to leave it out.

    ``auto`` gives 90% of the memory limit detected for this process (cgroup or SLURM) in
    bytes, and None when there is none (megahit then uses its own default); a fraction is
    passed through, and a size becomes bytes (see ``resources.parse_memory``).

    Raises:
        ValueError: If ``value`` is not a memory value.
    """
    limit = resources.memory_limit_bytes() if str(value).strip().lower() == "auto" else None
    return resources.parse_memory(value, limit)


def resolve_assembly_threads(requested: Optional[int], fallback: int) -> int:
    """Pick the megahit thread count, defaulting to single-thread on macOS.

    megahit 1.2.9's multithreaded k-mer sorting step segfaults on recent macOS
    (the prebuilt binary predates the current runtime), so when the user has not
    asked for a specific value we cap the assembly at one thread on Darwin.
    minimap2/samtools are unaffected and keep using ``fallback`` threads.

    Args:
        requested: An explicit thread count from the user, or None for the default.
        fallback: The thread count to use when no override is given (non-macOS).

    Returns:
        The number of threads to run megahit with.
    """
    if requested is not None:
        return requested
    if platform.system() == "Darwin":
        return 1
    return fallback


def is_singleton_file(path: Union[str, Path]) -> bool:
    """Whether ``path`` is the ``<genome>_s.fastq.gz`` file of reads that mapped without their mate."""
    return Path(path).name.endswith("_s.fastq.gz")


def _megahit_args(
    reads: List[Path],
    out_dir: Path,
    threads: int,
    min_contig_len: Optional[int],
    preset: Optional[str],
    k_flags: Optional[Dict[str, int]],
    tmp_dir: Optional[Path],
    memory: Optional[Union[int, float]] = None,
) -> List[str]:
    """Build the megahit argument list: input reads, thread count, k-mer/preset choice, an
    optional ``--memory`` (bytes, or a fraction of the node's memory) and an optional
    ``--tmp-dir`` (created and allow-listed here when given).

    Raises:
        ProcessingError: If the number of reads is unsupported, ``k_flags`` is given
            together with a preset, or ``tmp_dir`` is ``out_dir`` or a folder inside it
            (megahit refuses to run when its ``-o`` directory already exists).
    """
    args: List[str] = []
    if len(reads) == 3 and is_singleton_file(reads[2]):
        args += ["-1", str(reads[0]), "-2", str(reads[1]), "-r", str(reads[2])]
    elif len(reads) == 2:
        args += ["-1", str(reads[0]), "-2", str(reads[1])]
    elif len(reads) == 1:
        args += ["-r", str(reads[0])]
    else:
        raise ProcessingError(f"Expected 1 or 2 FASTQ files (plus a singleton file) to assemble, got {len(reads)}")

    args += ["--num-cpu-threads", str(threads), "-o", str(out_dir)]
    if memory is not None:
        args += ["--memory", str(memory)]
    if min_contig_len is not None:
        args += ["--min-contig-len", str(min_contig_len)]

    preset_active = preset not in (None, "default")
    if k_flags and preset_active:
        raise ProcessingError("megahit presets and explicit k values cannot be combined")
    if k_flags:
        for key, value in k_flags.items():
            args += [f"--{key}", str(value)]
    elif preset_active:
        args += ["--presets", str(preset)]

    if tmp_dir is not None:
        tmp_dir = Path(tmp_dir)
        _refuse_tmp_dir_inside(tmp_dir, out_dir)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        SecureSubprocess.add_allowed_root(tmp_dir)
        args += ["--tmp-dir", str(tmp_dir)]

    return args


def _refuse_tmp_dir_inside(tmp_dir: Path, out_dir: Path) -> None:
    """Raise ``ProcessingError`` when ``tmp_dir`` is ``out_dir`` or a folder inside it."""
    if tmp_dir.resolve().is_relative_to(out_dir.resolve()):
        raise ProcessingError(
            f"tmp_dir ({tmp_dir}) must not be output_dir or a folder inside it ({out_dir}); "
            "megahit refuses to run when its -o directory already exists"
        )


def _reuse_existing(
    out_dir: Path, expected: Optional[AssemblyInputs], accept_unmarked: bool, version: str, params: Any
) -> bool:
    """True when the assembly in ``out_dir`` can be kept; logs why it is redone otherwise.

    An unmarked folder that is accepted gets a marker for ``expected`` written in place.
    """
    state, reason = assembly_state(out_dir, expected, accept_unmarked)
    if state == "current":
        if expected is not None and not (out_dir / MARKER_NAME).exists():
            write_marker(out_dir, expected, version, params)
        logger.info("%s already assembled (%s); use --force to redo", out_dir, reason)
        return True
    if state == "incomplete":
        logger.warning("%s holds no %s (interrupted run?); assembling again", out_dir, CONTIGS_NAME)
    elif state == "stale":
        logger.info("%s is out of date (%s); assembling again", out_dir, reason)
    return False


def assemble_extracted_reads(
    reads: List[Path],
    output_dir: Union[str, Path],
    threads: int = 4,
    min_contig_len: Optional[int] = None,
    force: bool = False,
    preset: Optional[str] = "meta-sensitive",
    keep_intermediate: bool = False,
    k_flags: Optional[Dict[str, int]] = None,
    tmp_dir: Optional[Path] = None,
    memory: Optional[Union[int, float]] = None,
    *,
    expected: Optional[AssemblyInputs] = None,
    accept_unmarked: bool = True,
    version: str = "",
    params: Optional[Mapping[str, Any]] = None,
) -> Tuple[Path, bool]:
    """Assemble a set of extracted FASTQ files with megahit.

    Single-end input (one file) uses megahit ``-r``; paired input (two files) uses ``-1``/``-2``. The
    reads are expected to be a small, target-filtered set. Unless ``force`` is set, an existing
    assembly is kept when ``assembly_state`` reports it current: with ``expected`` None any folder
    holding contigs is kept, otherwise its marker must match ``expected`` (an unmarked folder is kept,
    and given a marker, when ``accept_unmarked`` is True). A folder without contigs (an interrupted
    run) is logged as a warning and assembled again.

    megahit writes into a hidden staging folder beside ``output_dir``; the marker (when ``expected``
    is given) is written into it, and the finished folder replaces ``output_dir`` by rename
    (``publish_assembly``). If megahit fails, the previous ``output_dir`` is left as it was and the
    staging folder is removed.

    Args:
        reads: One or two FASTQ files to assemble.
        output_dir: megahit output directory.
        threads: CPU threads.
        min_contig_len: Optional minimum contig length.
        force: If True, redo the assembly even if it already ran.
        preset: megahit ``--presets`` value (e.g. ``meta-sensitive``, ``meta-large``); ``None`` or
            ``"default"`` omits the flag and uses megahit's own defaults.
        keep_intermediate: If True, keep megahit's ``intermediate_contigs/`` folder instead of
            removing it once the assembly succeeds (useful for debugging a k-mer step, at the cost
            of extra disk space).
        k_flags: Explicit ``--k-min``/``--k-max``/``--k-step`` values (keys without the leading
            dashes, e.g. ``{"k-min": 21}``), used instead of a preset; combining the two is rejected.
        tmp_dir: Where megahit writes its scratch files (``--tmp-dir``); defaults to megahit's own
            choice under the staging folder when not given. megahit needs FIFOs for its scratch files, so
            a default landing on a filesystem without them (e.g. ExFAT) fails; a POSIX filesystem
            works around it.
        memory: megahit ``--memory``: bytes, or a fraction of the node's memory; None leaves the
            flag out (see ``assembly_memory``).
        expected: The inputs the assembly should be built from (``assembly_identity.assembly_inputs``);
            None keeps any folder that holds contigs and writes no marker.
        accept_unmarked: Whether a folder with contigs but no marker counts as current.
        version: The megahit version line, recorded in the marker (not compared).
        params: Further settings recorded in the marker (not compared).

    Returns:
        The megahit output directory and whether megahit actually ran (False when the assembly was
        already there), so the caller can leave an existing record alone.

    Raises:
        ProcessingError: If the number of reads is unsupported, ``k_flags`` is given together with a
            preset, ``tmp_dir`` is ``output_dir`` or a folder inside it, megahit itself fails, or
            megahit finishes without writing ``final.contigs.fa``.
    """
    out_dir = Path(output_dir)
    if not force and _reuse_existing(out_dir, expected, accept_unmarked, version, params):
        return out_dir, False
    if tmp_dir is not None:
        _refuse_tmp_dir_inside(Path(tmp_dir), out_dir)

    out_dir.parent.mkdir(parents=True, exist_ok=True)
    SecureSubprocess.add_allowed_root(out_dir.parent)
    staging = unique_temp_path(out_dir)
    args = _megahit_args(reads, staging, threads, min_contig_len, preset, k_flags, tmp_dir, memory)

    try:
        try:
            SecureSubprocess.run_secure("megahit", args)
        except subprocess.CalledProcessError as exc:
            tail = "\n".join((exc.stderr or "").strip().splitlines()[-5:])
            raise ProcessingError(f"megahit failed (exit {exc.returncode}) for {out_dir}:\n{tail}") from exc
        if not (staging / CONTIGS_NAME).is_file():
            raise ProcessingError(f"megahit finished without writing {CONTIGS_NAME} for {out_dir}")
        if not keep_intermediate:
            shutil.rmtree(staging / "intermediate_contigs", ignore_errors=True)
        if expected is not None:
            write_marker(staging, expected, version, params)
        publish_assembly(staging, out_dir)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    logger.info("Assembly written to %s", out_dir)
    return out_dir, True


def summarise_contigs(contigs: Union[str, Path]) -> Dict[str, Any]:
    """Contig count, total length, N50/N90, largest contig, GC fraction and the number of
    contigs at least 1 kb long, for a FASTA file.

    megahit headers carry ``len=<bp>``; when present that value is used for contig length,
    so the scan reads only header lines for length. GC content is always computed from the
    actual sequence lines, regardless of whether a header length is present.
    """
    path = Path(contigs)
    empty: Dict[str, Any] = {
        "contigs": 0,
        "total_bp": 0,
        "n50": 0,
        "n90": 0,
        "largest": 0,
        "gc": 0.0,
        "contigs_ge_1kb": 0,
    }
    if not path.exists():
        return empty
    lengths: List[int] = []
    current = 0
    have_current = False
    header_len = False
    gc_count = 0
    bases_seen = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if have_current:
                    lengths.append(current)
                have_current = True
                match = re.search(r"\blen=(\d+)", line)
                current = int(match.group(1)) if match else 0
                header_len = match is not None
            elif have_current:
                seq = line.strip()
                if not header_len:
                    current += len(seq)
                bases_seen += len(seq)
                gc_count += sum(1 for base in seq.upper() if base in "GC")
    if have_current:
        lengths.append(current)
    lengths.sort(reverse=True)
    total = sum(lengths)

    def _n_stat(numerator: int, denominator: int) -> int:
        """The shortest contig in the sorted-by-length prefix whose cumulative length is at
        least ``numerator/denominator`` of the total (N50 is numerator=1, denominator=2;
        N90 is numerator=9, denominator=10). Integer arithmetic only, to avoid floating-point
        error on a large total."""
        running = 0
        for length in lengths:
            running += length
            if running * denominator >= total * numerator:
                return length
        return 0

    n50 = _n_stat(1, 2)
    n90 = _n_stat(9, 10)
    gc = round(gc_count / bases_seen, 4) if bases_seen else 0.0
    contigs_ge_1kb = sum(1 for length in lengths if length >= 1000)
    return {
        "contigs": len(lengths),
        "total_bp": total,
        "n50": n50,
        "n90": n90,
        "largest": lengths[0] if lengths else 0,
        "gc": gc,
        "contigs_ge_1kb": contigs_ge_1kb,
    }


def fasta_length(path: Union[str, Path]) -> int:
    """Total base count of a FASTA file (sequence lines only, headers excluded)."""
    total = 0
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def megahit_version() -> str:
    """The installed megahit's version line (``MEGAHIT v1.2.9``), or an empty string if it cannot be read."""
    return tools.probe_tool("megahit").version_text
