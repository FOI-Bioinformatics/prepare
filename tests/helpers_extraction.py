"""Shared test double for minimap2/samtools/megahit/prefetch/fasterq-dump/pigz.

Used by read-extraction tests and by the download tests for the prefetch/split-3/
compression sequence in ``download_accession``.
"""

import functools
import gzip
import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from prepare.utils import tools

UNEQUAL_WARNING = "[W::mm_bseq_read_frag2] query files have different number of records; extra records skipped."

# Flags (across the tools this fake stands in for) whose next token is a value, not
# another flag or the positional argument.
_VALUE_FLAGS = frozenset({"-O", "--temp", "--threads", "--max-size", "-p"})

# ``samtools coverage`` header, as samtools 1.10+ writes it.
COVERAGE_HEADER = "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"

# Two contigs: 100 bp with 50 covered at mean depth 2, and 300 bp fully covered at mean
# depth 10 (breadth 350/400 = 0.875, length-weighted mean depth 3200/400 = 8.0).
DEFAULT_COVERAGE_ROWS = [("contig1", 100, 50, 2.0), ("contig2", 300, 300, 10.0)]


@pytest.fixture(autouse=True)
def tools_present():
    """Let ``extract_target_reads``' pre-flight tool check pass without minimap2/samtools/megahit.

    Autouse in every test module that imports it. These tests mock ``SecureSubprocess.run_secure``
    and may run where the tools are genuinely absent (the pip-only CI job). A test that exercises a
    missing tool applies its own, more specific ``shutil.which`` patch, which takes precedence for
    its duration. Version floors are checked in tests/test_tools.py with fake tools on PATH; here
    the check looks tools up without running them, so the mocked ``run_secure`` sees only the
    command's own calls.
    """
    without_versions = functools.partial(tools.require_tools, check_versions=False)
    with (
        patch("prepare.utils.tools.shutil.which", return_value="/usr/bin/tool"),
        patch("prepare.cli.commands.read_extraction.require_tools", without_versions),
    ):
        yield


def sam_record(qname: str, flag: int, rname: str = "chr1", pos: int = 100) -> str:
    """One minimal SAM alignment line for ``qname``/``flag``; the rest are placeholders."""
    mapq = "0" if flag & 0x100 else "60"
    return f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t4M\t*\t0\t0\tACGT\tIIII"


def sam_text(records) -> str:
    """SAM body text for ``(qname, flag)`` tuples, one ``sam_record`` line each."""
    return "".join(sam_record(qname, flag) + "\n" for qname, flag in records)


def _read_sam_records(path):
    """Parse the ``(qname, flag)`` pairs back out of a SAM file ``sam_text`` wrote."""
    records = []
    for line in Path(path).read_text().splitlines():
        if not line:
            continue
        fields = line.split("\t")
        records.append((fields[0], int(fields[1])))
    return records


def _filter_mask(args) -> int:
    """The bitmask a ``-F`` argument passes to ``samtools view``/``view -c``, or 0."""
    if "-F" in args:
        value = args[args.index("-F") + 1]
        return int(value, 16) if value.lower().startswith("0x") else int(value)
    return 0


def _kept_count(records, mask: int) -> int:
    """Count of ``(qname, flag)`` records samtools would keep after ``-F <mask>``."""
    return sum(1 for _, flag in records if (flag & mask) == 0)


def coverage_table(rows):
    """``samtools coverage`` output for ``rows`` of ``(rname, length, covbases, meandepth)``."""
    lines = [COVERAGE_HEADER]
    for rname, length, covbases, meandepth in rows:
        percent = 100.0 * covbases / length if length else 0.0
        lines.append(f"{rname}\t1\t{length}\t10\t{covbases}\t{percent:.4f}\t{meandepth}\t36\t60")
    return "\n".join(lines) + "\n"


def _positional(args):
    """Return the first token in ``args`` that is neither a flag nor a flag's value.

    Works regardless of where the positional argument (an SRA accession, or a path to
    a ``.sra``/FASTQ file) sits in the list, since real callers order it differently
    for the prefetch and the direct fasterq-dump command shapes.
    """
    skip_next = False
    for token in args:
        if skip_next:
            skip_next = False
            continue
        if token in _VALUE_FLAGS:
            skip_next = True
            continue
        if token.startswith("-"):
            continue
        return token
    return None


def _fake_tools(state):
    """Stand in for minimap2/samtools/megahit/prefetch/fasterq-dump/pigz.

    Records every call and creates the files each real tool would write.

    state keys: mapped (int, default 10, the post-filter/kept count reported for a BAM
    count), mapped_total (int, the pre-filter mapped count reported for a SAM count;
    defaults to ``mapped``), nonempty (flags whose output file gets a read, default
    ("-1", "-2")), unequal (bool, emit minimap2's mate-count warning), single (bool,
    fasterq-dump writes a single ``<acc>.fastq`` instead of a ``<acc>_1.fastq``/
    ``<acc>_2.fastq`` pair), reads (int, records per FASTQ file fasterq-dump writes,
    default 4), coverage_rows (``(rname, length, covbases, meandepth)`` tuples written by
    ``samtools coverage``; defaults to ``DEFAULT_COVERAGE_ROWS``), coverage_fail (bool,
    ``samtools coverage`` exits non-zero).

    sam_records (list of ``(qname, flag)``, default None): when given, minimap2 writes
    real SAM lines for these records instead of an empty file, honouring ``--sam-hit-only``
    (records with the unmapped flag, bit 0x4, are dropped before writing, as the real tool
    does), and every later ``samtools view -c``/``view -b -F ...`` call counts or filters
    them for real rather than reading ``mapped``/``mapped_total``/``coverage_mapped`` back
    off the state dict. Set it to pin a mapped/kept count to real filter semantics instead
    of a hard-coded number.

    minimap2_fail_stderr (str, default None): when given, every minimap2 alignment exits
    non-zero (``CalledProcessError``) with this text on stderr, as an old or broken minimap2 would.
    """

    def run(executable, args, **kwargs):
        state.setdefault("calls", []).append((executable, list(args)))
        result = MagicMock(returncode=0, stdout="", stderr="")
        if executable == "minimap2":
            if "-d" in args:
                # Building an index: create the .mmi the real tool would write (empty when the
                # state asks for a tool that exits 0 without writing).
                index_path = Path(args[args.index("-d") + 1])
                index_path.parent.mkdir(parents=True, exist_ok=True)
                index_path.write_bytes(b"" if state.get("empty_index") else b"MMI\x02")
            elif "-o" in args and state.get("minimap2_fail_stderr") is not None:
                raise subprocess.CalledProcessError(1, ["minimap2", *args], stderr=state["minimap2_fail_stderr"])
            elif "-o" in args:
                # Aligning: create the SAM output the real tool would write.
                sam_path = Path(args[args.index("-o") + 1])
                sam_path.parent.mkdir(parents=True, exist_ok=True)
                records = state.get("sam_records")
                if records is not None:
                    if "--sam-hit-only" in args:
                        records = [r for r in records if not (r[1] & 0x4)]
                    sam_path.write_text(sam_text(records))
                else:
                    sam_path.write_text("")
                if state.get("unequal"):
                    result.stderr = UNEQUAL_WARNING
        if executable == "samtools" and args[:2] == ["view", "-c"]:
            target = args[-1]
            records = state.get("sam_records")
            if records is not None and str(target).endswith(".sam"):
                actual = _read_sam_records(target) if Path(target).exists() else []
                result.stdout = f"{_kept_count(actual, _filter_mask(args))}\n"
            elif (
                records is not None
                and str(target).endswith(".bam")
                and Path(target).exists()
                and Path(target).read_bytes().strip().isdigit()
            ):
                result.stdout = f"{int(Path(target).read_bytes())}\n"
            elif str(target).endswith(".sam"):
                result.stdout = f"{state.get('mapped_total', state.get('mapped', 10))}\n"
            elif Path(target).name == "coverage.bam":
                result.stdout = f"{state.get('coverage_mapped', 5)}\n"
            else:
                result.stdout = f"{state.get('mapped', 10)}\n"
        if executable == "samtools" and args[0] == "view" and "-b" in args:
            state.setdefault("view_filter_calls", []).append(list(args))
            out_path = Path(args[args.index("-o") + 1])
            out_path.parent.mkdir(parents=True, exist_ok=True)
            records = state.get("sam_records")
            if records is not None:
                in_path = Path(args[-1])
                actual = _read_sam_records(in_path) if in_path.exists() else []
                out_path.write_bytes(str(_kept_count(actual, _filter_mask(args))).encode())
            else:
                out_path.write_bytes(b"")
        if executable == "samtools" and args[0] == "cat":
            out_path = Path(args[args.index("-o") + 1])
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_bytes(b"")
        if executable == "samtools" and args[0] == "sort":
            out_path = Path(args[args.index("-o") + 1])
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_bytes(b"")
        if executable == "samtools" and args[0] == "coverage":
            out_path = Path(args[args.index("-o") + 1])
            if state.get("coverage_fail"):
                # A real failure can leave a partial table behind.
                out_path.write_text(COVERAGE_HEADER + "\n")
                raise subprocess.CalledProcessError(1, ["samtools", *args], stderr="coverage failed")
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(coverage_table(state.get("coverage_rows", DEFAULT_COVERAGE_ROWS)))
        if executable == "samtools" and args[0] == "fastq":
            for flag in ("-1", "-2", "-0", "-s"):
                if flag in args and flag in state.get("nonempty", ("-1", "-2")):
                    path = Path(args[args.index(flag) + 1])
                    path.parent.mkdir(parents=True, exist_ok=True)
                    with gzip.open(path, "wt") as handle:
                        handle.write("@r1\nACGT\n+\nIIII\n")
        if executable == "megahit":
            if "--version" in args:
                result.stdout = ""
            elif "-o" in args:
                out_dir = Path(args[args.index("-o") + 1])
                out_dir.mkdir(parents=True, exist_ok=True)
                (out_dir / "final.contigs.fa").write_text(">c1 len=100\nACGT\n>c2 len=50\nACGT\n")
                (out_dir / "intermediate_contigs").mkdir(parents=True, exist_ok=True)
        if executable == "prefetch":
            accession = _positional(args)
            cache_dir = Path(args[args.index("-O") + 1])
            acc_dir = cache_dir / accession
            acc_dir.mkdir(parents=True, exist_ok=True)
            # NCBI serves some runs only as the smaller .sralite format.
            suffix = ".sralite" if state.get("sralite") else ".sra"
            (acc_dir / f"{accession}{suffix}").write_bytes(b"")
        if executable == "fasterq-dump":
            positional = _positional(args)
            accession = Path(positional).stem if positional.endswith((".sra", ".sralite")) else positional
            out_dir = Path(args[args.index("-O") + 1])
            out_dir.mkdir(parents=True, exist_ok=True)
            record = "@r\nACGT\n+\nIIII\n" * state.get("reads", 4)
            if state.get("single"):
                (out_dir / f"{accession}.fastq").write_text(record)
            else:
                (out_dir / f"{accession}_1.fastq").write_text(record)
                (out_dir / f"{accession}_2.fastq").write_text(record)
        if executable == "pigz":
            path = Path(_positional(args))
            target = path.with_suffix(path.suffix + ".gz")
            with open(path, "rb") as src, gzip.open(target, "wb") as dst:
                dst.write(src.read())
            path.unlink()
        return result

    return run
