"""The ``sra_sketch`` command: one sourmash signature per downloaded dataset, with local containment.

Each accession's FASTQ files are read once and sketched at k = 21, 31 and 51 (scaled 1000,
abundance tracked) into ``<output-dir>/<ACC>.sig.gz`` (``prepare.sra.sketch``). A sketch whose
recorded file snapshot still matches is skipped; a signature cached beside a store copy is
copied instead of recomputed (``prepare.store.sketches``). With ``--genome-fasta`` or
``--signature`` the containment of the genome in every dataset is merged into
``<matches-folder>/<genome_id>.csv`` in Branchwater's layout (``prepare.sra.sketch_containment``),
so ``parse_containment`` and ``select_datasets`` read it unchanged. A ``"sketch"`` analysis is
recorded per accession in the project registry.
"""

import argparse
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from prepare.cli.base import BaseCommand, accessions_from_args, anchored_folder, positive_int, resolve_command_store
from prepare.core import settings
from prepare.core.exceptions import DataAccessError, PrepareError, ValidationError
from prepare.core.optional import require
from prepare.data import run_log
from prepare.data.branchwater_search import load_signature, sketch_fasta
from prepare.data.file_io import copy_file, visible_files
from prepare.data.registry import load_registry, record_analysis
from prepare.data.registry_batch import registry_batch
from prepare.data.registry_blocks import _analysis_summary
from prepare.data.sra import is_transient_folder
from prepare.data.sra_metadata import _resolved_sidecar_path
from prepare.sra.dataset_stats import DATASET_READ_ERRORS, accession_fastq_files
from prepare.sra.sketch import (
    KSIZES,
    SCALED,
    DatasetSketch,
    signature_path,
    sketch_dataset,
    sketch_is_current,
    sketch_summary,
    write_signature,
)
from prepare.sra.sketch_containment import containment_rows, genome_minhash, merge_containment_csv
from prepare.store.sketches import cached_sketch, store_sketch
from prepare.store.usage import record_usage_many
from prepare.utils import resources
from prepare.utils.progress import ProgressReporter

ANALYSIS_NAME = "sketch"
_PURPOSE = "Sketching reads"


@dataclass
class Planned:
    """One accession and what the run will do for it."""

    accession: str
    files: List[Path] = field(default_factory=list)
    acc_dir: Optional[Path] = None
    action: str = "sketch"  # sketch, reuse, skip or failed
    cached: Optional[Path] = None
    cached_block: Optional[Dict[str, Any]] = None


@dataclass
class Done:
    """One accession's signature on disk and the summary to record."""

    accession: str
    path: Path
    summary: Dict[str, Any]
    recorded: bool = True


class SRASketchCommand(BaseCommand):
    """Build sourmash sketches of downloaded datasets and compare them with a genome."""

    @property
    def name(self) -> str:
        """Command name."""
        return "sra_sketch"

    @property
    def help(self) -> str:
        """One-line help."""
        return "Sourmash sketches (k=21,31,51, scaled 1000) of downloaded datasets, with local containment of a genome"

    @property
    def group(self) -> str:
        """Pipeline step."""
        return "Reads"

    def records_run(self, args: argparse.Namespace) -> bool:
        """Every run but a dry run is added to the project's run log."""
        return not getattr(args, "dry_run", False)

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add the command's options."""
        parser.add_argument("--fastq-folder", default="fastq", help="Folder holding one folder per accession")
        parser.add_argument("--accessions-file", default=None, help="Sketch the accessions listed here, one per line")
        parser.add_argument(
            "--accession",
            action="append",
            default=None,
            help="Sketch this accession (repeatable). With neither this nor --accessions-file, every "
            "accession folder in --fastq-folder is sketched",
        )
        parser.add_argument("--output-dir", default="sketches", help="Folder for the <ACC>.sig.gz signature files")
        genome = parser.add_mutually_exclusive_group()
        genome.add_argument("--genome-fasta", default=None, help="Genome to compare with each dataset (FASTA)")
        genome.add_argument(
            "--signature", default=None, help="Existing plain sourmash signature of the genome (k=21, scaled=1000)"
        )
        parser.add_argument(
            "--genome-id", default=None, help="Name of the matches CSV (default: the genome file's stem)"
        )
        parser.add_argument("--matches-folder", default="matches", help="Folder the containment CSV is written to")
        parser.add_argument(
            "--max-workers",
            type=positive_int,
            default=None,
            help="Accessions sketched at once (default: the CPUs available to this job, at most max_workers_cap; "
            "each worker holds three sketches of its dataset in memory)",
        )
        parser.add_argument("--force", action="store_true", help="Sketch again even when the sketch is up to date")
        parser.add_argument("--dry-run", action="store_true", help="List what would be sketched and write nothing")
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")

    # ------------------------------------------------------------------------------------ planning

    def _resolve_accessions(self, args: argparse.Namespace) -> List[str]:
        """The accessions named on the command line, else every accession folder in --fastq-folder."""
        named = accessions_from_args(args.accessions_file, args.accession)
        if named:
            return named
        folder = Path(args.fastq_folder)
        found = [d.name for d in visible_files(folder, dirs=True) if not is_transient_folder(d.name)]
        if not found:
            raise ValidationError(f"No accession folders in {folder}; give --accession or --accessions-file")
        return found

    def _plan(self, args: argparse.Namespace, accessions: List[str]) -> List[Planned]:
        """Decide per accession whether to sketch, reuse the store's signature, skip, or fail."""
        registry = load_registry(args.registry)
        output_dir = Path(args.output_dir)
        planned: List[Planned] = []
        for accession in accessions:
            item = Planned(accession)
            item.files = accession_fastq_files(args.fastq_folder, accession)
            if not item.files:
                self.logger.warning("No FASTQ files found for %s under %s", accession, args.fastq_folder)
                item.action = "failed"
                planned.append(item)
                continue
            item.acc_dir = Path(args.fastq_folder) / accession
            target = signature_path(accession, output_dir)
            summary = _analysis_summary(registry, accession, ANALYSIS_NAME)
            if not args.force and sketch_is_current(summary, item.files, target):
                item.action = "skip"
            elif not args.force:
                cached = cached_sketch(item.acc_dir, _resolved_sidecar_path(item.acc_dir), item.files)
                if cached is not None:
                    item.action = "reuse"
                    item.cached, item.cached_block = cached
            planned.append(item)
        return planned

    def _dry_run(self, planned: List[Planned]) -> int:
        for item in planned:
            self.emit(f"{item.accession}: {item.action} ({len(item.files)} file(s))")
        counts = {action: sum(1 for p in planned if p.action == action) for action in ("sketch", "reuse", "skip")}
        self.emit(f"Would sketch {counts['sketch']}, reuse {counts['reuse']} from the store, skip {counts['skip']}")
        return 0

    # ------------------------------------------------------------------------------------- sketching

    @staticmethod
    def _max_workers(args: argparse.Namespace) -> int:
        if args.max_workers:
            return int(args.max_workers)
        return max(1, min(resources.available_cpus(), settings.active().max_workers_cap))

    def _sketch_one(self, item: Planned, output_dir: Path, stop) -> Done:
        """Sketch or copy one accession's signature into the project and cache a fresh one in the store."""
        target = signature_path(item.accession, output_dir)
        if item.action == "reuse" and item.cached is not None and item.cached_block is not None:
            copy_file(item.cached, target)
            block = {k: v for k, v in item.cached_block.items() if k not in ("file", "computed")}
            return Done(item.accession, target, {**block, "source": "store"})
        sketch: DatasetSketch = sketch_dataset(item.accession, item.files, stop=stop)
        write_signature(sketch, target)
        summary = sketch_summary(sketch)
        if item.acc_dir is not None:
            store_sketch(_resolved_sidecar_path(item.acc_dir), target, summary)
        return Done(item.accession, target, summary)

    def _sketch_all(self, args: argparse.Namespace, planned: List[Planned]) -> Tuple[List[Done], List[str], bool]:
        """Run the sketches on a thread pool; returns (done, failed accessions, interrupted)."""
        output_dir = Path(args.output_dir)
        term = getattr(args, "_termination", None)
        stop = term.stop if term is not None else None
        work = [p for p in planned if p.action in ("sketch", "reuse")]
        done: List[Done] = []
        failed: List[str] = [p.accession for p in planned if p.action == "failed"]
        interrupted = False
        progress = ProgressReporter("sra_sketch", len(work), settings.active().progress_every, logger=self.logger)
        executor = ThreadPoolExecutor(max_workers=self._max_workers(args))
        futures: Dict[Future, Planned] = {}
        try:
            for item in work:
                futures[executor.submit(self._sketch_one, item, output_dir, stop)] = item
            for future in as_completed(futures):
                item = futures[future]
                outcome = self._collect(future, item, progress)
                if outcome is None:
                    failed.append(item.accession)
                else:
                    done.append(outcome)
                if stop is not None and stop.is_set():
                    interrupted = True
                    break
        except KeyboardInterrupt:
            interrupted = True
            if stop is not None:
                stop.set()
        finally:
            executor.shutdown(wait=True, cancel_futures=True)
        if interrupted:
            # Workers that finished after the stop are kept; a cancelled one left no file.
            for future, item in futures.items():
                if future.done() and not future.cancelled() and item.accession not in {d.accession for d in done}:
                    outcome = self._collect(future, item, progress)
                    if outcome is not None:
                        done.append(outcome)
        return done, failed, interrupted

    def _collect(self, future: Future, item: Planned, progress: ProgressReporter) -> Optional[Done]:
        """The finished future's result, or None (logged) when the accession could not be sketched."""
        try:
            outcome = future.result()
        except (*DATASET_READ_ERRORS, DataAccessError) as e:
            if isinstance(e, InterruptedError):
                return None
            self.logger.warning("Failed to sketch %s: %s", item.accession, e)
            progress.update(ok=False)
            return None
        progress.update(ok=True)
        self.logger.log(
            progress.item_level,
            "%s: %s reads, %s hashes at k=%d",
            item.accession,
            outcome.summary.get("reads"),
            outcome.summary.get("hashes", {}).get("21"),
            KSIZES[0],
        )
        return outcome

    # ----------------------------------------------------------------------------------- containment

    @staticmethod
    def _genome(args: argparse.Namespace) -> Optional[Tuple[str, Dict[str, Any]]]:
        """``(genome_id, signature)`` from --genome-fasta or --signature, or None."""
        source = args.genome_fasta or args.signature
        if not source:
            return None
        signature = sketch_fasta(args.genome_fasta) if args.genome_fasta else load_signature(args.signature)
        return args.genome_id or Path(source).stem, signature

    def _containment(
        self, args: argparse.Namespace, genome: Tuple[str, Dict[str, Any]], planned: List[Planned]
    ) -> Dict[str, Tuple[float, float]]:
        """Write the matches CSV for every accession with a signature on disk; returns the values by accession."""
        genome_id, signature = genome
        output_dir = Path(args.output_dir)
        paths = {
            p.accession: signature_path(p.accession, output_dir)
            for p in planned
            if p.action != "failed" and signature_path(p.accession, output_dir).is_file()
        }
        rows = containment_rows(genome_minhash(signature), paths)
        target = merge_containment_csv(Path(args.matches_folder) / f"{genome_id}.csv", rows)
        self.emit(f"Containment of {genome_id} written to {target} ({len(rows)} dataset(s))")
        for accession, value in sorted(rows, key=lambda r: -r[1])[:10]:
            self.emit(f"  {accession}: {value:.4f}")
        ksize = KSIZES[0]
        return {acc: (value, round(value ** (1 / ksize), 4) if value > 0 else 0.0) for acc, value in rows}

    # ------------------------------------------------------------------------------------- recording

    def _record(self, args: argparse.Namespace, done: List[Done]) -> None:
        """Record a "sketch" analysis and an "analysed" store usage for every sketched or reused accession."""
        if not done:
            return
        store = resolve_command_store(args, load_registry(args.registry))
        with registry_batch(args.registry, flush_every=None, flush_seconds=None) as batch:
            for item in done:
                batch.apply(
                    partial(
                        record_analysis,
                        accession=item.accession,
                        analysis=ANALYSIS_NAME,
                        output=item.path,
                        summary=item.summary,
                    ),
                    item.accession,
                )
        if batch.registry is not None:
            record_usage_many(store, batch.registry, [(d.accession, "", "analysed", ANALYSIS_NAME) for d in done])

    def _note_run(
        self,
        args: argparse.Namespace,
        planned: List[Planned],
        done: List[Done],
        failed: List[str],
        genome: Optional[Tuple[str, Dict[str, Any]]],
        containment: Dict[str, Tuple[float, float]],
    ) -> None:
        sources = {d.accession: d.summary.get("source") for d in done}
        summary = {
            "accessions": len(planned),
            "sketched": sum(1 for s in sources.values() if s == "computed"),
            "reused": sum(1 for s in sources.values() if s == "store"),
            "skipped": sum(1 for p in planned if p.action == "skip"),
            "failed": len(failed),
            "reads": sum(int(d.summary.get("reads") or 0) for d in done),
            "bases": sum(int(d.summary.get("bases") or 0) for d in done),
            "ksizes": ",".join(str(k) for k in KSIZES),
            "scaled": SCALED,
            "max_workers": self._max_workers(args),
            "genome_id": genome[0] if genome else None,
            "containment_rows": len(containment),
        }
        run_log.note_run(args, summary=summary)
        rows: Dict[str, Dict[str, Any]] = {}
        for item in done:
            hashes = item.summary.get("hashes", {})
            row: Dict[str, Any] = {
                "reads": item.summary.get("reads"),
                "bases": item.summary.get("bases"),
                "hashes_k21": hashes.get("21"),
                "hashes_k31": hashes.get("31"),
                "hashes_k51": hashes.get("51"),
                "source": item.summary.get("source"),
            }
            if item.accession in containment:
                row["containment"], row["cani"] = containment[item.accession]
            rows[item.accession] = row
        if rows:
            run_log.note_rows(args, rows, "analyses", ANALYSIS_NAME)

    # ------------------------------------------------------------------------------------------ run

    def _run(self, args: argparse.Namespace) -> int:
        folder = Path(args.fastq_folder)
        if not folder.is_dir():
            raise ValidationError(f"FASTQ folder {folder} does not exist")
        accessions = self._resolve_accessions(args)
        planned = self._plan(args, accessions)
        if args.dry_run:
            return self._dry_run(planned)
        require("sourmash", "sourmash", _PURPOSE)
        genome = self._genome(args)

        done, failed, interrupted = self._sketch_all(args, planned)
        containment: Dict[str, Tuple[float, float]] = {}
        if genome is not None and not interrupted:
            containment = self._containment(args, genome, planned)
        self._record(args, done)
        self._note_run(args, planned, done, failed, genome, containment)
        self.emit(
            f"Sketched {sum(1 for d in done if d.summary.get('source') == 'computed')}, "
            f"reused {sum(1 for d in done if d.summary.get('source') == 'store')}, "
            f"skipped {sum(1 for p in planned if p.action == 'skip')}, failed {len(failed)} of {len(planned)}"
        )
        if interrupted:
            raise KeyboardInterrupt("sra_sketch stopped")
        if failed:
            self.logger.warning("%d of %d accession(s) could not be sketched", len(failed), len(planned))
            return 1
        return 0

    def execute(self, args: argparse.Namespace) -> int:
        """Sketch the datasets; 1 when any accession could not be sketched."""
        try:
            registry = getattr(args, "registry", None)
            args.fastq_folder = str(anchored_folder(args.fastq_folder, registry))
            args.output_dir = str(anchored_folder(args.output_dir, registry))
            args.matches_folder = str(anchored_folder(args.matches_folder, registry))
            return self._run(args)
        except (PrepareError, OSError) as e:
            return self.fail(e, "Sketching failed")
        except Exception as e:  # noqa: B902 - top-level catch: keep the traceback, return 1
            self.logger.exception("Sketching failed: %s", e)
            return 1
