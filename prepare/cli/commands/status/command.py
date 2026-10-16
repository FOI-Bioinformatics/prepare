"""
The `status` command: local inventory and registry stages.

Reports what PREPARE has already downloaded locally (SRA reads, NCBI metadata,
genome assemblies) so a user can see what is available without re-downloading,
and reports where every accession sits in the registry (screened, selected,
excluded, downloaded, analysed, extracted, assembled). When no registry file
exists yet, the report is reconstructed in memory from what is on disk.

The report is built by `prepare.processing.status_report`, the suggested next steps by
`suggest` and the text output by `render_text`.

`status --init` is the one place that parses metadata XML files under the registry lock
(`fill_metadata_from_xml` inside the `registry_update` that writes the bootstrapped registry):
the registry file did not exist before that write, so no other writer is waiting on it, and one
write leaves no half-filled registry behind.
"""

import argparse
from functools import partial
from pathlib import Path
from typing import Any, Dict, Tuple

from prepare.cli.base import BaseCommand
from prepare.cli.commands.status.render_text import print_report
from prepare.cli.commands.status.suggest import next_steps
from prepare.core.exceptions import DataAccessError, PrepareError
from prepare.data import run_log
from prepare.data.file_io import write_csv
from prepare.data.metadata_fields import fill_metadata_from_xml
from prepare.data.registry import (
    ProjectPaths,
    Registry,
    STAGES,
    bootstrap_from_disk,
    load_registry,
    project_root,
    registry_path,
)
from prepare.data.registry_batch import registry_update
from prepare.data.registry_reconcile import ReconcilePlan, StoreReconcileReport, apply_reconcile, scan_reconcile
from prepare.processing.status_report import build_report, to_dataframes


def _adopt_bootstrap(registry: Registry, built: Registry, metadata_folder: Path) -> Registry:
    """Fill the empty ``registry`` loaded under the lock with what ``status --init`` rebuilt from disk.

    Raises ``DataAccessError`` when the registry file exists by now: another process created it
    after this one checked, and overwriting it would lose what that process recorded. Once the
    bootstrapped fields are copied in, fills in every metadata block bootstrap recorded without a
    spot count from its ``<accession>_metadata.xml`` file, in this same registry write, so a kill
    partway through cannot leave the registry persisted with metadata half-filled.
    """
    if registry.path is not None and registry.path.exists():
        raise DataAccessError(
            f"Registry {registry.path} was created by another process while this one rebuilt it from disk; "
            "run status --reconcile to update it instead"
        )
    for name in ("created", "genomes", "datasets", "project", "store"):
        setattr(registry, name, getattr(built, name))
    fill_metadata_from_xml(registry, metadata_folder)
    return registry


def _apply_plan(registry: Registry, plan: ReconcilePlan) -> Tuple[StoreReconcileReport, Registry]:
    """Apply a reconcile plan to the registry loaded under the lock; return the report and that registry."""
    return apply_reconcile(registry, plan), registry


class StatusCommand(BaseCommand):
    """Command to report locally available data and the registry's per-accession stages."""

    @property
    def name(self) -> str:
        """Return the command name."""
        return "status"

    @property
    def help(self) -> str:
        """Return the command's help text."""
        return "Report which SRA reads, metadata, and genomes are already available locally"

    @property
    def group(self) -> str:
        """Return the pipeline-step group this command is listed under."""
        return "Reads"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        """Add the command's arguments."""
        parser.add_argument("--fastq-folder", default="fastq", help="Folder holding per-accession FASTQ downloads")
        parser.add_argument("--metadata-folder", default="metadata", help="Folder holding NCBI metadata XML")
        parser.add_argument("--genomes-folder", default="genomes", help="Folder holding genome FASTA files")
        parser.add_argument(
            "--targeted-folder", default="targeted", help="Root folder of extracted reads and assemblies"
        )
        parser.add_argument("--matches-folder", default="matches", help="Folder of Branchwater match CSVs")
        parser.add_argument(
            "--accessions-file",
            help="Optional file of SRA accessions (one per line) to reconcile against local FASTQ/metadata",
        )
        parser.add_argument(
            "--parsed-containment",
            help="Optional parsed containment table; its sample accessions are the wanted list",
        )
        parser.add_argument(
            "--registry",
            default=None,
            help="Registry file (default: prepare_registry.json found upwards from here)",
        )
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")
        parser.add_argument("--stage", choices=list(STAGES), default=None, help="List the accessions in one stage")
        parser.add_argument(
            "--genome",
            action="append",
            default=None,
            help="Restrict extraction and assembly stages to a genome id (repeatable)",
        )
        parser.add_argument("--init", action="store_true", help="Create the registry from what is on disk")
        parser.add_argument(
            "--reconcile",
            action="store_true",
            help="Compare the registry with the disk and record missing downloads",
        )
        parser.add_argument(
            "--export-tsv", default=None, help="Write <PREFIX>_datasets.tsv and <PREFIX>_extractions.tsv"
        )
        parser.add_argument("--next", action="store_true", help="Suggest the commands that advance the most accessions")
        parser.add_argument("--list-missing", action="store_true", help="Also print the accessions that are missing")
        parser.add_argument("--json", action="store_true", help="Emit the report as JSON")

    def records_run(self, args: argparse.Namespace) -> bool:
        """Only ``--init`` and ``--reconcile``, which write the registry, are added to the run log."""
        return bool(getattr(args, "init", False) or getattr(args, "reconcile", False))

    @staticmethod
    def _note_run(args: argparse.Namespace, registry: Registry, report: Dict[str, Any]) -> None:
        """What ``--init`` or ``--reconcile`` did: the dataset count and, for a reconcile, each finding's count."""
        if not (args.init or args.reconcile):
            return
        summary: Dict[str, Any] = {"action": "reconcile" if args.reconcile else "init"}
        summary["datasets"] = len(registry.datasets)
        if args.reconcile:
            summary["drift"] = {key: len(found) for key, found in (report.get("drift") or {}).items()}
        run_log.note_run(args, summary=summary)

    def _export_tsv(self, registry: Registry, prefix: str) -> None:
        datasets_df, extractions_df = to_dataframes(registry)
        datasets_path = f"{prefix}_datasets.tsv"
        extractions_path = f"{prefix}_extractions.tsv"
        write_csv(datasets_df, datasets_path, sep="\t")
        write_csv(extractions_df, extractions_path, sep="\t", index=False)
        self.logger.info("Wrote %s and %s", datasets_path, extractions_path)

    def _emit(self, args: argparse.Namespace, report: Dict[str, Any], registry: Registry) -> None:
        if args.json:
            self.emit_json(report)
        else:
            print_report(args, report, registry, self.emit)

    def execute(self, args: argparse.Namespace) -> int:
        """Build and emit the status report; return 0 on success or 1 on a handled error.

        Loads the registry (or bootstraps and, with ``--init``, writes one from what is on disk),
        optionally reconciles it against the filesystem with ``--reconcile``, and assembles a
        report of registry stages, downloads, genomes, store usage, and drift, emitted as JSON or
        as text depending on ``args.json``. Returns 1 rather than raising when ``--init`` is asked
        for an already-existing registry, when ``--reconcile`` is asked with no registry yet, or
        when report assembly raises a :class:`~prepare.core.exceptions.PrepareError` (logged
        and swallowed here so the CLI exits cleanly instead of printing a traceback).
        """
        try:
            registry_file = registry_path(args.registry)
            # Relative folders belong to the project that owns the registry, not the working directory.
            paths = ProjectPaths.anchored(
                project_root(registry_file),
                fastq=args.fastq_folder,
                metadata=args.metadata_folder,
                genomes=args.genomes_folder,
                targeted=args.targeted_folder,
                matches=args.matches_folder,
            )
            existed = registry_file.exists()
            if existed and args.init:
                self.logger.error(
                    "Registry already exists at %s; run status --reconcile to update it from disk, "
                    "or remove the file to rebuild it",
                    registry_file,
                )
                return 1
            if args.reconcile and not existed:
                self.logger.error(
                    "No registry at %s to reconcile; create one first with: prepare status --init",
                    registry_file,
                )
                return 1
            if existed:
                registry = load_registry(registry_file)
            else:
                registry = bootstrap_from_disk(paths, args.accessions_file, args.parsed_containment, registry_file)
                if args.init:
                    registry = registry_update(
                        registry_file, partial(_adopt_bootstrap, built=registry, metadata_folder=paths.metadata)
                    )
                    self.logger.info("Registry written to %s", registry_file)

            drift = None
            if args.reconcile:
                # Reads are counted on a snapshot, without the lock; only the apply runs under it,
                # on the registry as it is then.
                plan = scan_reconcile(registry, paths)
                drift, registry = registry_update(registry_file, partial(_apply_plan, plan=plan))

            report = build_report(registry, args, paths, registry_file, existed, drift)
            self._note_run(args, registry, report)
            if args.next:
                report["next"] = next_steps(registry, paths)
            if args.export_tsv:
                self._export_tsv(registry, args.export_tsv)

            if not existed and not args.init:
                self.logger.info(
                    "No registry yet; the stages above were reconstructed from disk. Run: prepare status --init"
                )

            self._emit(args, report, registry)
            return 0
        except PrepareError as e:
            return self.fail(e, "Error building status report")
