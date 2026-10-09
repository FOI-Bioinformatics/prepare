"""CLI command that turns containment results into an accession list for downloads."""

import argparse
from functools import partial
from pathlib import Path
from typing import Any, Dict, List

from prepare.cli.base import BaseCommand
from prepare.core.constants import DEFAULT_CONTAINMENT_THRESHOLD, DEFAULT_TOP_N
from prepare.core.exceptions import PrepareError
from prepare.data import registry_blocks as rb
from prepare.data import run_log
from prepare.data.file_io import write_text_atomic
from prepare.data.defaults import resolve_metadata_table
from prepare.data.registry import Registry, load_registry, project_root, query, record_selection
from prepare.data.registry_batch import registry_update
from prepare.processing.selection import RunFilters, parse_size, select_accessions_ranked


def _positive_int(value: str) -> int:
    """argparse type for --top-n: rejects zero and negative values with a clear message."""
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError(f"--top-n must be a positive integer, got {value!r}")
    return parsed


def _size_bytes(value: str) -> int:
    """argparse type for --max-run-size: a byte count with an optional K/M/G/T suffix, above zero."""
    try:
        parsed = parse_size(value)
    except PrepareError as e:
        raise argparse.ArgumentTypeError(str(e)) from e
    if parsed < 1:
        raise argparse.ArgumentTypeError(f"--max-run-size must be above zero, got {value!r}")
    return parsed


def _non_negative_int(value: str) -> int:
    """argparse type for --min-spots and --max-spots: an integer of zero or more."""
    try:
        parsed = int(value)
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"expected a whole number, got {value!r}") from e
    if parsed < 0:
        raise argparse.ArgumentTypeError(f"expected zero or a positive whole number, got {value!r}")
    return parsed


def _non_blank(value: str) -> str:
    """argparse type for --platform: rejects an empty or whitespace-only value, which would match no run."""
    if not value.strip():
        raise argparse.ArgumentTypeError("--platform needs a platform name such as ILLUMINA, got a blank value")
    return value


def _recorded_selection_files(registry: Registry) -> set:
    """Resolved paths of every file a recorded selection names.

    A relative recorded path is resolved against both the project root and the current
    directory, since the recording run may have been started from either.
    """
    files = set()
    for acc in registry.datasets:
        recorded = (rb.selection_block(registry, acc) or rb.SelectionBlock()).output
        if not recorded:
            continue
        path = Path(recorded)
        if path.is_absolute():
            files.add(path.resolve())
        else:
            files.add((project_root(registry) / path).resolve())
            files.add(path.resolve())
    return files


def _note_run(
    args: argparse.Namespace, ranked: List[Dict[str, Any]], downloaded_count: int, excluded_count: int
) -> None:
    """Counts and criteria for the run-log summary; rank, column and value per accession for its detail."""
    summary: Dict[str, Any] = {"selected": len(ranked), "already_downloaded": downloaded_count}
    summary["excluded"] = excluded_count
    summary["column"] = ranked[0]["column"] if ranked else (args.genome_id or "max_containment")
    run_log.note_run(args, summary={**summary, "threshold": args.threshold, "top_n": args.top_n})
    rows = {r["accession"]: {"rank": r["rank"], "column": r["column"], "value": r["value"]} for r in ranked}
    run_log.note_rows(args, rows, "selection")


class SelectDatasetsCommand(BaseCommand):
    """Write the accessions that meet a containment threshold (and optional metadata filter)."""

    @property
    def name(self) -> str:
        return "select_datasets"

    @property
    def help(self) -> str:
        return "Write an accessions file from parsed containment, for download_sra"

    @property
    def group(self) -> str:
        return "Reads"

    def records_run(self, args: argparse.Namespace) -> bool:
        return not getattr(args, "no_record", False)

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--parsed-containment", default="parsed_containment.txt", help="Table from parse_containment"
        )
        genome_group = parser.add_mutually_exclusive_group()
        genome_group.add_argument(
            "--genome-id", default=None, help="Genome column to rank on (default: max_containment)"
        )
        genome_group.add_argument(
            "--genome-ids", nargs="+", default=None, help="Multiple genome columns to rank on together"
        )
        parser.add_argument(
            "--require",
            choices=["any", "all"],
            default="any",
            help="With --genome-ids, require any or all columns to meet the threshold",
        )
        parser.add_argument(
            "--threshold", type=float, default=DEFAULT_CONTAINMENT_THRESHOLD, help="Minimum containment, inclusive"
        )
        parser.add_argument(
            "--top-n",
            type=_positive_int,
            default=None,
            help=f"Keep only the top N accessions after filtering (suggested default: {DEFAULT_TOP_N})",
        )
        parser.add_argument(
            "--metadata-file",
            default=None,
            help="Metadata table for filtering (default: metadata_table.txt, else metadata/branchwater_metadata.txt)",
        )
        parser.add_argument("--metadata-column", default=None, help="Metadata column to filter on")
        parser.add_argument("--metadata-value", default=None, help="Required value in that column")
        parser.add_argument(
            "--max-run-size",
            type=_size_bytes,
            default=None,
            metavar="BYTES",
            help="Drop runs whose Run_Size exceeds BYTES (suffixes K, M, G, T are powers of 10, e.g. 500M); "
            "needs the NCBI metadata table, and a run with no Run_Size value is dropped",
        )
        parser.add_argument(
            "--min-spots",
            type=_non_negative_int,
            default=None,
            metavar="N",
            help="Drop runs with fewer than N spots (Run_Total_Spots); a run with no value is dropped",
        )
        parser.add_argument(
            "--max-spots",
            type=_non_negative_int,
            default=None,
            metavar="N",
            help="Drop runs with more than N spots (Run_Total_Spots); a run with no value is dropped",
        )
        parser.add_argument(
            "--platform",
            type=_non_blank,
            default=None,
            metavar="NAME",
            help="Keep only runs from this sequencing platform, e.g. ILLUMINA (case-insensitive); "
            "a run with no Platform value is dropped",
        )
        parser.add_argument("--output", default="accessions.txt", help="Output file, one accession per line")
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")
        parser.add_argument(
            "--skip-excluded",
            action=argparse.BooleanOptionalAction,
            default=True,
            help="Drop accessions already marked excluded in the registry",
        )
        parser.add_argument(
            "--skip-downloaded",
            action="store_true",
            help="Drop accessions already marked downloaded in the registry",
        )
        parser.add_argument(
            "--no-record",
            dest="no_record",
            action="store_true",
            help="Write the output file and log the counts, but do not record the selection in "
            "the registry or the run in the run log (for an exploratory run that should not "
            "redefine the target list); refuses an --output that a recorded selection names, so "
            "give it a scratch file",
        )

    def execute(self, args: argparse.Namespace) -> int:
        try:
            run_filters = RunFilters(
                max_run_size=args.max_run_size,
                min_spots=args.min_spots,
                max_spots=args.max_spots,
                platform=args.platform,
            )
            # Checked before the metadata table is looked for, so a project without one still sees
            # that the spot bounds cannot both hold rather than a missing-table error.
            run_filters.validate()
            metadata_file = None
            if args.metadata_column or run_filters.active():
                metadata_file = resolve_metadata_table(args.metadata_file)

            # A snapshot, read without the lock: ranking can take a while on a large table, and the
            # selection is recorded afterwards against the registry as it is then (see below).
            registry = load_registry(args.registry)
            if args.no_record and Path(args.output).resolve() in _recorded_selection_files(registry):
                self.logger.error(
                    "%s is the recorded selection's file; give --no-record a different --output "
                    "or run without --no-record",
                    args.output,
                )
                return 1
            exclude = set(query(registry, "excluded")) if args.skip_excluded else set()
            excluded_count = len(exclude)
            downloaded = set(query(registry, "downloaded")) if args.skip_downloaded else set()
            downloaded_count = len(downloaded)
            exclude |= downloaded

            ranked = select_accessions_ranked(
                parsed_containment=args.parsed_containment,
                genome_id=args.genome_id,
                threshold=args.threshold,
                metadata_file=metadata_file,
                metadata_column=args.metadata_column,
                metadata_value=args.metadata_value,
                top_n=args.top_n,
                exclude=exclude,
                genome_ids=args.genome_ids,
                require=args.require,
                run_filters=run_filters,
            )
            accessions = [accession for accession, _, _ in ranked]

            output = Path(args.output)
            output.parent.mkdir(parents=True, exist_ok=True)
            write_text_atomic(output, "".join(f"{acc}\n" for acc in accessions))
            self.logger.info("Wrote %d accession(s) to %s", len(accessions), output)
            if not accessions:
                self.logger.warning("No accessions met the criteria; %s is empty", output)
            self.logger.info(
                "%d already downloaded, %d excluded, %d selected",
                downloaded_count,
                excluded_count,
                len(accessions),
            )

            if args.no_record:
                self.logger.info("Not recording this selection in the registry (--no-record)")
                return 0

            ranked_records = [
                {"accession": accession, "rank": i + 1, "column": column, "value": value}
                for i, (accession, column, value) in enumerate(ranked)
            ]
            _note_run(args, ranked_records, downloaded_count, excluded_count)
            record = partial(
                record_selection,
                accessions=accessions,
                criteria={
                    # The column the ranking actually used, which for --genome-ids is the
                    # combined label select_accessions_ranked builds (e.g. "A+B"), not the
                    # single --genome-id this falls back to when nothing was selected.
                    "column": ranked[0][1] if ranked else (args.genome_id or "max_containment"),
                    "threshold": args.threshold,
                    # The metadata table actually used, whether it was --metadata-file or the
                    # default resolve_metadata_table autodetected (see the top of execute()),
                    # so a reselect built from these criteria carries the real path instead of
                    # relying on cwd-relative autodetection happening again.
                    "metadata_file": (str(metadata_file.resolve()) if metadata_file is not None else None),
                    "metadata_column": args.metadata_column,
                    "metadata_value": args.metadata_value,
                    "table": str(args.parsed_containment),
                    "top_n": args.top_n,
                    "require": args.require,
                    "genome_ids": args.genome_ids,
                    "skip_excluded": args.skip_excluded,
                    "skip_downloaded": args.skip_downloaded,
                    "max_run_size": args.max_run_size,
                    "min_spots": args.min_spots,
                    "max_spots": args.max_spots,
                    "platform": args.platform,
                },
                output=output,
                ranked=ranked_records,
            )
            # Loaded again inside the lock, so a download or exclusion another process recorded
            # while this one was ranking is kept; record_selection unselects earlier selections
            # against that registry, not against the snapshot.
            registry_update(registry.path, record)
            return 0
        except PrepareError as e:
            return self.fail(e, "Error selecting datasets")
