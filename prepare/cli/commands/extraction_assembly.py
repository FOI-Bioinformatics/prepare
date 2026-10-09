"""The assembly step of ``extract_target_reads --assemble``.

Each sample with extracted reads is assembled while its extraction lock is held, so a sample
being extracted or assembled by another process is skipped rather than assembled twice. An
existing assembly is kept only when its identity marker (``data/assembly_identity.py``) matches
the reads and settings of this run; an assembly built from other reads or settings, or a folder
an interrupted megahit left without contigs, is assembled again. A sample whose assembly fails is
reported and the remaining samples are still assembled.
"""

import argparse
import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from prepare.cli.base import BaseCommand
from prepare.core.exceptions import ConfigurationError, ProcessingError
from prepare.core.settings import setting_for
from prepare.data import registry_blocks as rb
from prepare.data.assembly import assembly_memory
from prepare.data.assembly_identity import (
    MARKER_NAME,
    AssemblyInputs,
    assembly_inputs,
    assembly_state,
    read_marker,
    sweep_staging,
)
from prepare.data.extraction_locks import LockHeld, sample_extraction_lock
from prepare.data.read_extraction import (
    ExtractionResult,
    assemble_extracted_reads,
    assembly_coverage,
    fasta_length,
    megahit_version,
    resolve_assembly_threads,
    summarise_contigs,
)
from prepare.data.registry import Registry, clear_assembly, load_registry, record_assembly
from prepare.data.registry_assembly import (
    assembly_predates_extraction,
    legacy_assembly_current,
    set_assembly_inputs,
)
from prepare.data.registry_batch import registry_update
from prepare.data.registry_timing import Stopwatch, set_assembly_timing
from prepare.store.layout import StorePaths
from prepare.store.usage import record_usage_safe

# Recorded in the marker written into an unmarked folder that is accepted as it is: the megahit
# that built it is not known, so the version installed now is not recorded as its provenance.
_ACCEPTED_PARAMS = {"accepted_without_marker": True}


@dataclass
class AssemblyOutcome:
    """What the assembly step did with each sample.

    ``assembled``: megahit ran and the result was recorded. ``reused``: an existing assembly
    matched this run's inputs and was kept. ``failed``: sample -> reason, for a sample whose
    assembly could not be made. ``busy``: samples whose lock another process held.
    """

    assembled: List[str] = field(default_factory=list)
    reused: List[str] = field(default_factory=list)
    failed: Dict[str, str] = field(default_factory=dict)
    busy: List[str] = field(default_factory=list)


@dataclass
class _RunSettings:
    """The megahit settings shared by every sample of one run."""

    threads: int
    memory: Optional[Any]
    version: str
    genome_length: int

    def params(self, args: argparse.Namespace) -> Dict[str, Any]:
        """The parameters recorded in the registry and in the marker of a new assembly."""
        return {"threads": self.threads, "min_contig_len": args.min_contig_len, "preset": args.assembly_preset}


def resolve_assembly_memory(args: argparse.Namespace, log: logging.Logger) -> Optional[Any]:
    """The megahit ``--memory`` value for ``--assembly-memory`` (or its setting), or None to omit it.

    Raises:
        ConfigurationError: If the value given is not a memory value.
    """
    value = setting_for(args, "assembly_memory")
    try:
        memory = assembly_memory(value)
    except ValueError as e:
        raise ConfigurationError(f"--assembly-memory: {e}") from e
    if memory is None:
        log.info("No memory limit detected for this job; megahit uses its own default")
    else:
        log.info("megahit --memory %s (--assembly-memory %s)", memory, value)
    return memory


def _run_settings(args: argparse.Namespace, log: logging.Logger) -> _RunSettings:
    """Threads, memory, megahit version and reference length for this run."""
    threads = resolve_assembly_threads(args.assembly_threads, args.threads)
    if args.assembly_threads is None and threads < args.threads:
        log.info(
            "Running megahit single-threaded on macOS (its parallel sort is unstable here); "
            "override with --assembly-threads"
        )
    memory = resolve_assembly_memory(args, log)
    return _RunSettings(threads, memory, megahit_version(), fasta_length(args.genome_fasta))


def _megahit_scratch(args: argparse.Namespace) -> Tuple[Path, bool]:
    """megahit's ``--tmp-dir`` for one sample, and whether it was created here (and is removed after).

    megahit needs FIFOs for its scratch files, which some filesystems (e.g. ExFAT) do not provide;
    ``--temp-folder`` points it elsewhere when given. Otherwise a folder unique to this run is made
    under the output root, beside every per-accession folder and never inside an assembly folder,
    so two runs sharing one output folder never remove each other's scratch.
    """
    temp_folder = setting_for(args, "temp_folder")
    if temp_folder:
        return Path(temp_folder), False
    Path(args.output_folder).mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(dir=args.output_folder, prefix=".megahit-tmp-")), True


def _assembly_stats(
    args: argparse.Namespace, run: _RunSettings, out_dir: Path, reads: List[Path], mapped_records: int
) -> Dict[str, Any]:
    """Contig statistics of ``out_dir``, with the read coverage of the contigs unless ``--no-coverage``."""
    contigs_path = out_dir / "final.contigs.fa"
    stats: Dict[str, Any] = dict(summarise_contigs(contigs_path))
    stats["genome_fraction_estimate"] = (stats["total_bp"] / run.genome_length) if run.genome_length else None
    if not args.no_coverage:
        stats.update(
            assembly_coverage(contigs_path, reads, args.preset, run.threads, out_dir, mapped_reads=mapped_records)
        )
    return stats


class _SampleAssembly:
    """The assembly of one sample, run while its extraction lock is held."""

    def __init__(
        self, args: argparse.Namespace, run: _RunSettings, accession: str, reads: List[Path], mapped_records: int
    ) -> None:
        """Bind one sample's reads to the run's settings; nothing is read or run yet."""
        self.args = args
        self.run = run
        self.accession = accession
        self.reads = reads
        self.mapped_records = mapped_records
        self.out_dir = Path(args.output_folder) / accession / f"{args.genome_id}_assembly"

    def _identity(self) -> Tuple[AssemblyInputs, bool, bool]:
        """This run's inputs, whether an unmarked folder is accepted, and whether a record exists.

        The registry is read once per sample, under the sample lock, for the extraction date.
        """
        registry = load_registry(self.args.registry)
        accession, genome_id = self.accession, self.args.genome_id
        block = rb.extraction_block(registry, accession, genome_id)
        expected = assembly_inputs(
            self.reads, block.date if block else None, self.args.assembly_preset, self.args.min_contig_len
        )
        accept_unmarked = legacy_assembly_current(
            registry, accession, genome_id, self.args.assembly_preset, self.args.min_contig_len
        ) and not assembly_predates_extraction(registry, accession, genome_id)
        return expected, accept_unmarked, block is not None and block.assembly is not None

    def run_megahit(self, expected: AssemblyInputs, accept_unmarked: bool, accepted_as_is: bool) -> Tuple[bool, Any]:
        """Call megahit unless the folder is current; returns (ran, (started, seconds) or None)."""
        version, params = self.run.version, self.run.params(self.args)
        if accepted_as_is:
            version, params = "", _ACCEPTED_PARAMS
        tmp_dir, own_tmp_dir = _megahit_scratch(self.args)
        watch = Stopwatch()
        try:
            _, ran = assemble_extracted_reads(
                self.reads,
                self.out_dir,
                threads=self.run.threads,
                min_contig_len=self.args.min_contig_len,
                force=self.args.force,
                preset=self.args.assembly_preset,
                keep_intermediate=self.args.keep_intermediate,
                tmp_dir=tmp_dir,
                memory=self.run.memory,
                expected=expected,
                accept_unmarked=accept_unmarked,
                version=version,
                params=params,
            )
        finally:
            if own_tmp_dir:
                shutil.rmtree(tmp_dir, ignore_errors=True)
        # megahit's own run time (and its scratch removal), not the statistics or coverage below.
        return ran, (watch.lap() if ran else None)

    def _record(
        self, expected: AssemblyInputs, version: str, params: Dict[str, Any], timing: Optional[Tuple[str, float]]
    ) -> Tuple[Registry, int]:
        """Record the assembly in ``out_dir`` with its inputs and timing in one registry update.

        Returns the registry as written and the number of contigs.
        """
        stats = _assembly_stats(self.args, self.run, self.out_dir, self.reads, self.mapped_records)
        accession, genome_id = self.accession, self.args.genome_id
        started, seconds = timing or (None, None)

        def write(reg: Registry) -> Registry:
            record_assembly(reg, accession, genome_id, self.out_dir, stats, version, params)
            set_assembly_inputs(reg, accession, genome_id, expected.to_dict())
            set_assembly_timing(reg, accession, genome_id, started, seconds)
            return reg

        return registry_update(self.args.registry, write), int(stats.get("contigs", 0))

    def assemble(self, store: Optional[StorePaths]) -> str:
        """Assemble (or keep) this sample's assembly and record it; returns "assembled" or "reused"."""
        sweep_staging(self.out_dir)
        expected, accept_unmarked, has_record = self._identity()
        state, _ = assembly_state(self.out_dir, expected, accept_unmarked)
        accepted_as_is = state == "current" and not (self.out_dir / MARKER_NAME).exists()
        if has_record and (self.args.force or state != "current"):
            # The record describes an assembly about to be replaced. It is cleared before megahit
            # runs: a failed run keeps the previous folder on disk, but no record of it.
            registry_update(self.args.registry, lambda reg: clear_assembly(reg, self.accession, self.args.genome_id))
            has_record = False
        ran, timing = self.run_megahit(expected, accept_unmarked, accepted_as_is and not self.args.force)
        accession, genome_id = self.accession, self.args.genome_id
        if ran:
            reg, contigs = self._record(expected, self.run.version, self.run.params(self.args), timing)
        elif not has_record:
            # Current by its marker, but the registry lost the record: describe the folder from its marker.
            marker = read_marker(self.out_dir) or {}
            version, params = str(marker.get("megahit_version") or ""), dict(marker.get("params") or {})
            reg, contigs = self._record(expected, version, params, None)
        else:
            if accepted_as_is:
                # The legacy record now also states the inputs the marker was just written with.
                registry_update(
                    self.args.registry, lambda r: set_assembly_inputs(r, accession, genome_id, expected.to_dict())
                )
            return "reused"
        # After the registry lock is released: the store catalogue has its own lock.
        record_usage_safe(store, reg, accession, genome_id, "assembled", detail=f"{contigs} contigs")
        return "assembled" if ran else "reused"


def _assemble_one(
    command: BaseCommand,
    args: argparse.Namespace,
    run: _RunSettings,
    sample: Tuple[str, List[Path], int],
    store: Optional[StorePaths],
    outcome: AssemblyOutcome,
) -> None:
    """Assemble one sample under its lock and file it under ``outcome``."""
    accession, reads, mapped_records = sample
    try:
        with sample_extraction_lock(Path(args.output_folder), accession, args.genome_id):
            result = _SampleAssembly(args, run, accession, reads, mapped_records).assemble(store)
    except LockHeld:
        command.logger.info("%s: assembly against %s is in progress elsewhere; skipped", accession, args.genome_id)
        outcome.busy.append(accession)
        return
    except (ProcessingError, OSError, subprocess.SubprocessError) as exc:
        # megahit, publishing the folder, or the coverage mapping (minimap2/samtools) failed.
        command.logger.error("%s: assembly against %s failed: %s", accession, args.genome_id, exc)
        outcome.failed[accession] = str(exc)
        return
    (outcome.assembled if result == "assembled" else outcome.reused).append(accession)


def assemble_samples(
    command: BaseCommand,
    args: argparse.Namespace,
    with_reads: Dict[str, List[Path]],
    results: Dict[str, ExtractionResult],
    store: Optional[StorePaths] = None,
) -> AssemblyOutcome:
    """Assemble every sample that has mapped reads, recording each assembly as it lands.

    ``with_reads`` holds the read files this run's extraction step returned, so an assembly is
    compared with the reads on disk now. A signal (``args._termination.stop``) is checked before
    each sample and stops the loop there; the caller raises afterwards, so the run exits 130
    without losing what was already assembled. A sample whose assembly fails (megahit, or an
    I/O error while publishing it) is listed in ``failed`` and the loop goes on.
    """
    run = _run_settings(args, command.logger)
    outcome = AssemblyOutcome()
    term = getattr(args, "_termination", None)
    for accession, reads in with_reads.items():
        if term is not None and term.stop.is_set():
            command.logger.warning("Stopping before %s: interrupted", accession)
            break
        result = results.get(accession)
        mapped = result.mapped_records if result is not None else 0
        _assemble_one(command, args, run, (accession, reads, mapped), store, outcome)
    command.logger.info(
        "Assembly: %d sample(s) assembled, reused %d, failed %d, busy elsewhere %d",
        len(outcome.assembled),
        len(outcome.reused),
        len(outcome.failed),
        len(outcome.busy),
    )
    return outcome
