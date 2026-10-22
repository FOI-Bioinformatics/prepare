"""
Metadata-related CLI commands.
"""

import argparse
import logging
from functools import partial
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd

from prepare.cli.base import finite_float, BaseCommand, anchored_folder
from prepare.core.exceptions import DataAccessError, PrepareError
from prepare.core.settings import require_email, setting_for
from prepare.data import registry_blocks as rb
from prepare.data.defaults import resolve_metadata_table
from prepare.data.file_io import atomic_path, copy_file
from prepare.data.metadata import (
    check_metadata_attributes,
    download_metadata,
    parse_metadata,
    parse_metadata_xml,
)
from prepare.data.metadata_fields import FIELD_COLUMNS, metadata_fields
from prepare.data.registry import Registry, record_metadata, registry_path
from prepare.data.registry_batch import registry_batch, registry_update
from prepare.processing.counts import count_metadata
from prepare.plugins.base import visualizer_registry
from prepare.store.resolve import resolve_optional_store
from prepare.visualization.plots import plot_metadata_counts

logger = logging.getLogger(__name__)


# The field-to-column mapping and the row-to-fields function now live in data.metadata_fields,
# shared with fill_metadata_from_xml; kept as aliases here for the rest of this module and tests.
_FIELD_COLUMNS = FIELD_COLUMNS
_metadata_fields = metadata_fields

# The table columns ParseMetadataCommand reads per row: the accession and the _FIELD_COLUMNS sources.
_ROW_COLUMNS = tuple(dict.fromkeys(("Run_ID",) + tuple(column for _, column in _FIELD_COLUMNS)))


class DownloadMetadataCommand(BaseCommand):
    """Command for downloading metadata for SRA accessions."""

    @property
    def name(self) -> str:
        return "download_metadata"

    @property
    def help(self) -> str:
        return "Download metadata for SRA accessions"

    @property
    def group(self) -> str:
        return "Metadata"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--email",
            default=None,
            help="Email address for NCBI API access (default: PREPARE_NCBI_EMAIL or config [runtime] ncbi_email)",
        )
        parser.add_argument("--matches-folder", default="matches", help="Folder containing match files")
        parser.add_argument(
            "--metadata-folder",
            default="metadata",
            help="Folder to save downloaded metadata",
        )
        parser.add_argument(
            "--threshold",
            type=finite_float,
            default=0.0,
            help="Threshold for containment values (inclusive)",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Calculate number of accessions without downloading",
        )
        parser.add_argument(
            "--accessions-file",
            default=None,
            help="File of accessions to fetch metadata for; replaces the matches folder scan",
        )
        parser.add_argument(
            "--api-key",
            default=None,
            help="NCBI API key for a higher rate limit (default: PREPARE_NCBI_API_KEY or NCBI_API_KEY)",
        )
        parser.add_argument(
            "--batch-size",
            type=int,
            default=200,
            help="Accessions per NCBI request, from 1 to 500 (default: 200)",
        )
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")

    def _share_with_store(self, args: argparse.Namespace, registry, downloaded: dict) -> None:
        """Copy each fetched XML into the store's ``metadata/`` folder, when a store resolves.

        NCBI's spot count for a run is what makes a store dataset's sidecar verifiable, and the
        store branch of ``download_sra`` already reads that folder. Copying is best effort: a
        store that cannot be written to costs the sharing, never the metadata this project just
        fetched.
        """
        paths = resolve_optional_store(getattr(args, "data_root", None), rb.store_block(registry).root)
        if paths is None:
            return
        try:
            paths.metadata.mkdir(parents=True, exist_ok=True)
            for xml_path in downloaded.values():
                # One rename, so another project or store_verify never reads a half-copied XML.
                copy_file(xml_path, paths.metadata / Path(xml_path).name)
        except (OSError, DataAccessError) as e:
            self.logger.warning("Could not copy metadata into the store at %s: %s", paths.metadata, e)
            return
        self.logger.info("Copied %d metadata file(s) into the store at %s", len(downloaded), paths.metadata)

    def _parse_downloaded(self, downloaded: Dict[str, Path], term) -> Dict[str, Tuple[Any, Dict[str, Any]]]:
        """Parse each fetched XML into its registry fields, without the registry lock.

        Checked before each accession: a signal stops the loop there, leaving every later
        accession's XML on disk (already written atomically by ``download_metadata``) but not
        yet recorded; what was parsed so far is still recorded by the caller.
        """
        parsed: Dict[str, Tuple[Any, Dict[str, Any]]] = {}
        for accession, xml_path in downloaded.items():
            if term is not None and term.stop.is_set():
                self.logger.warning("Stopping before %s: interrupted", accession)
                break
            try:
                fields = _metadata_fields(parse_metadata_xml(xml_path))
            except (PrepareError, ValueError, OSError, ET.ParseError) as e:
                self.logger.warning(f"Could not parse metadata for {accession}: {e}; recorded the file path only")
                fields = {}
            parsed[accession] = (xml_path, fields)
        return parsed

    def execute(self, args: argparse.Namespace) -> int:
        email = require_email(args)
        try:
            args.matches_folder = str(anchored_folder(args.matches_folder, getattr(args, "registry", None)))
            args.metadata_folder = str(anchored_folder(args.metadata_folder, getattr(args, "registry", None)))
            failures: Dict[str, str] = {}
            downloaded = download_metadata(
                email=email,
                matches_folder=args.matches_folder,
                metadata_folder=args.metadata_folder,
                threshold=args.threshold,
                dry_run=args.dry_run,
                accessions_file=args.accessions_file,
                api_key=setting_for(args, "ncbi_api_key"),
                batch_size=args.batch_size,
                failures=failures,
            )
            term = getattr(args, "_termination", None)
            if not args.dry_run and downloaded:
                # Every XML is parsed before the registry lock is taken; the records are then
                # written in one transaction on the registry as it is at that point.
                registry = _record_all_metadata(args.registry, self._parse_downloaded(downloaded, term))
                self._share_with_store(args, registry, downloaded)
            if term is not None and term.stop.is_set():
                raise KeyboardInterrupt("download_metadata stopped")
            if failures:
                self.logger.error("No metadata for %d accession(s): %s", len(failures), ", ".join(sorted(failures)))
                return 1
            return 0
        except PrepareError as e:
            return self.fail(e, "Error downloading metadata")


def _record_all_metadata(registry_arg: Optional[str], parsed: Dict[str, Tuple[Any, Dict[str, Any]]]) -> Registry:
    """Record every ``accession -> (xml_path, fields)`` in one registry transaction; return the registry written.

    Each record is queued on a batch that writes once, on exit, so the registry is read inside
    the lock and a download or selection another process recorded meanwhile is kept. A record
    that fails is logged and dropped; the summary line counts the dropped records, so a run
    that still exits 0 says how many are missing.
    """
    target = registry_path(registry_arg)
    # Paths are recorded relative to the registry's folder, as project_root would give.
    root = target.parent.resolve()
    with registry_batch(target, flush_every=None, flush_seconds=None) as batch:
        for accession, (xml_path, fields) in parsed.items():
            batch.apply(
                partial(record_metadata, accession=accession, xml_path=xml_path, fields=fields, root=root), accession
            )
    if batch.dropped:
        logger.warning(
            "Recorded metadata for %d accession(s) in the registry; %d dropped (see the errors above)",
            len(parsed) - batch.dropped,
            batch.dropped,
        )
    elif parsed:
        logger.info("Recorded metadata for %d accession(s) in the registry", len(parsed))
    if batch.registry is None:
        # Nothing to record: still write the registry, as a run with an empty table always has.
        return registry_update(target, lambda registry: registry)
    return batch.registry


class ParseMetadataCommand(BaseCommand):
    """Command for parsing downloaded metadata files."""

    @property
    def name(self) -> str:
        return "parse_metadata"

    @property
    def help(self) -> str:
        return "Parse downloaded metadata files"

    @property
    def group(self) -> str:
        return "Metadata"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--metadata-folder",
            default="metadata",
            help="Folder containing metadata files",
        )
        parser.add_argument(
            "--metadata-table-file",
            default="metadata_table.txt",
            help="File where the parsed metadata will be stored",
        )
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")

    @staticmethod
    def _row_record(metadata_folder: Path, row: Mapping[str, Any]) -> Optional[Tuple[str, Path, Dict[str, Any]]]:
        accession = row.get("Run_ID")
        if accession is None or pd.isna(accession):
            return None
        return str(accession), metadata_folder / f"{accession}_metadata.xml", _metadata_fields(row)

    def _table_records(self, df: pd.DataFrame, metadata_folder: Path, term) -> Dict[str, Tuple[Any, Dict[str, Any]]]:
        """The registry record of every row of the parsed table, stopping early on a signal."""
        # Only the columns the registry records, as plain dicts: a wide table (one column per
        # sample attribute) makes a pandas Series per row costly.
        columns = [column for column in _ROW_COLUMNS if column in df.columns]
        parsed: Dict[str, Tuple[Any, Dict[str, Any]]] = {}
        for row in df[columns].to_dict("records"):
            if term is not None and term.stop.is_set():
                self.logger.warning("Stopping metadata parsing: interrupted")
                break
            record = self._row_record(metadata_folder, row)
            if record is not None:
                parsed[record[0]] = (record[1], record[2])
        return parsed

    def execute(self, args: argparse.Namespace) -> int:
        try:
            args.metadata_folder = str(anchored_folder(args.metadata_folder, getattr(args, "registry", None)))
            errors: List[str] = []
            df = parse_metadata(args.metadata_folder, args.metadata_table_file, errors=errors)
            term = getattr(args, "_termination", None)
            # Parsed without the registry lock; recorded in one transaction.
            _record_all_metadata(args.registry, self._table_records(df, Path(args.metadata_folder), term))
            if term is not None and term.stop.is_set():
                raise KeyboardInterrupt("parse_metadata stopped")
            if errors:
                self.logger.error("%d metadata file(s) could not be parsed: %s", len(errors), ", ".join(errors))
                return 1
            return 0
        except PrepareError as e:
            return self.fail(e, "Error parsing metadata")


class CheckMetadataAttributesCommand(BaseCommand):
    """Command for counting how often each sample-attribute column is populated."""

    @property
    def name(self) -> str:
        return "check_metadata_attributes"

    @property
    def help(self) -> str:
        return "Count how often each metadata attribute is populated"

    @property
    def group(self) -> str:
        return "Metadata"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--file-path",
            default=None,
            help="Parsed metadata table (default: metadata_table.txt, else metadata/branchwater_metadata.txt)",
        )
        parser.add_argument(
            "--output-file",
            default="metadata_attribute_counts.txt",
            help="Path to save the attribute counts",
        )

    def execute(self, args: argparse.Namespace) -> int:
        try:
            check_metadata_attributes(str(resolve_metadata_table(args.file_path)), args.output_file)
            return 0
        except PrepareError as e:
            return self.fail(e, "Error checking metadata attributes")


class CountMetadataCommand(BaseCommand):
    """Command for counting metadata values by genome."""

    @property
    def name(self) -> str:
        return "count_metadata"

    @property
    def help(self) -> str:
        return "Count metadata values by genome"

    @property
    def group(self) -> str:
        return "Metadata"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--summary-file",
            default="parsed_containment.txt",
            help="Summary file",
        )
        parser.add_argument(
            "--metadata-file",
            default=None,
            help="Parsed metadata table (default: metadata_table.txt, else metadata/branchwater_metadata.txt)",
        )
        parser.add_argument(
            "--metadata-column",
            required=True,
            help="Name of the column in the metadata file",
        )
        parser.add_argument(
            "--threshold",
            type=finite_float,
            default=0.5,
            help="Threshold for containment values (inclusive)",
        )
        parser.add_argument(
            "--output-file",
            default="metadata_counts.txt",
            help="Path to the output file",
        )
        parser.add_argument("--stat-file", default=None, help="Path to the statistics file")

    def execute(self, args: argparse.Namespace) -> int:
        try:
            count_metadata(
                summary_file=args.summary_file,
                metadata_file=str(resolve_metadata_table(args.metadata_file)),
                metadata_column=args.metadata_column,
                threshold=args.threshold,
                output_file=args.output_file,
                stat_file=args.stat_file,
            )
            return 0
        except PrepareError as e:
            return self.fail(e, "Error counting metadata")


class PlotMetadataCountsCommand(BaseCommand):
    """Command for plotting metadata counts."""

    @property
    def name(self) -> str:
        return "plot_metadata_counts"

    @property
    def help(self) -> str:
        return "Plot metadata counts"

    @property
    def group(self) -> str:
        return "Metadata"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument("--file-path", required=True, help="Counts table from count_metadata (or its _stats file)")
        parser.add_argument("--title", default=None, help="Title for the plot")
        parser.add_argument(
            "--plot-type",
            default="bar",
            choices=["bar", "pie", "radar"],
            help="Type of plot to generate",
        )
        parser.add_argument("--colors", default=None, help="Colors or colormap name")
        parser.add_argument("--show-title", action="store_true", help="Whether to display the title")
        parser.add_argument(
            "--save-format",
            default=None,
            choices=["png", "jpg", "pdf", "svg"],
            help="Format to save the figure",
        )

    def execute(self, args: argparse.Namespace) -> int:
        try:
            plot_metadata_counts(
                file_path=args.file_path,
                title=args.title,
                plot_type=args.plot_type,
                colors=args.colors,
                show_title=args.show_title or args.title is not None,
                save_format=args.save_format,
            )
            if args.save_format is None:
                self.logger.warning("no --save-format given; nothing was written")
            return 0
        except PrepareError as e:
            return self.fail(e, "Error plotting metadata counts")


class PlotMapCommand(BaseCommand):
    """Command for drawing sample locations from a metadata table on a map."""

    @property
    def name(self) -> str:
        return "plot_map"

    @property
    def help(self) -> str:
        return "Draw sample locations from a metadata table on a map"

    @property
    def group(self) -> str:
        return "Metadata"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument("--metadata-file", required=True, help="Metadata table (tab-separated)")
        parser.add_argument("--output", required=True, help="Image file to write (png, jpg, pdf or svg)")
        parser.add_argument(
            "--lat-lon-column", default="lat_lon", help="Column holding coordinates, e.g. '38.98 N 77.10 W'"
        )

    def execute(self, args: argparse.Namespace) -> int:
        fig = None
        try:
            data = pd.read_csv(resolve_metadata_table(args.metadata_file), sep="\t")
            plugin = visualizer_registry.get("map")
            suffix = Path(args.output).suffix.lstrip(".").lower()
            with atomic_path(args.output) as tmp:
                fig = plugin.create_plot(
                    data,
                    lat_lon_column=args.lat_lon_column,
                    output_file=tmp,
                    output_format=suffix if suffix in ("png", "jpg", "pdf", "svg") else "png",
                )
            self.logger.info("Saved map to %s", args.output)
            return 0
        except PrepareError as e:
            return self.fail(e, "Error plotting map")
        finally:
            if fig is not None:
                plt.close(fig)
