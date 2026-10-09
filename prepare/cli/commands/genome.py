"""Genome search, download, and preparation CLI commands."""

import argparse
import csv
from pathlib import Path
from typing import Optional

from prepare.cli.base import BaseCommand
from prepare.core.constants import GENOME_FASTA_GLOBS
from prepare.core.exceptions import PrepareError
from prepare.data.file_io import open_atomic, visible_files, write_text_atomic
from prepare.data.genome_download import download_genomes, extract_and_organize, partition_present_genomes
from prepare.data.gtdb import (
    get_accessions_for_genus,
    get_accessions_for_species,
)
from prepare.data.registry import record_genome, registry_transaction, resolve_project_path
from prepare.utils.tools import require_tools


def _genome_name(filename: str) -> str:
    """Strip the FASTA suffix (and .gz) from a genome file name."""
    for suffix in (".fna.gz", ".fasta.gz", ".fa.gz", ".fna", ".fasta", ".fa"):
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return filename


class GenomeSearchCommand(BaseCommand):
    """Command to search GTDB for genome accessions."""

    @property
    def name(self) -> str:
        return "genome_search"

    @property
    def help(self) -> str:
        return "Search GTDB for genome accessions by species or genus"

    @property
    def group(self) -> str:
        return "Genomes"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument(
            "--species",
            help="Species name to search (e.g. 'Salmonella enterica')",
        )
        group.add_argument(
            "--genus",
            help="Genus name to search (e.g. 'Salmonella')",
        )
        parser.add_argument(
            "--output",
            default=None,
            help="Output file path (default: print to stdout)",
        )
        parser.add_argument(
            "--all",
            dest="representative_only",
            action="store_false",
            default=True,
            help="Include all genomes, not just representatives",
        )
        parser.add_argument(
            "--format",
            dest="output_format",
            choices=["list", "tsv"],
            default="list",
            help="Output format",
        )

    def execute(self, args: argparse.Namespace) -> int:
        try:
            if args.species:
                accessions = get_accessions_for_species(args.species, representative_only=args.representative_only)
                label = args.species
            else:
                accessions = get_accessions_for_genus(args.genus, representative_only=args.representative_only)
                label = args.genus

            if not accessions:
                self.logger.warning("No accessions found for '%s'", label)
                return 0

            self.logger.info("Found %d accessions for '%s'", len(accessions), label)

            if args.output_format == "tsv":
                lines = [f"{acc}\t{label}" for acc in accessions]
            else:
                lines = accessions

            output_text = "\n".join(lines) + "\n"

            if args.output:
                write_text_atomic(args.output, output_text)
                self.logger.info("Wrote accessions to %s", args.output)
            else:
                self.emit_raw(output_text)

            return 0
        except PrepareError as e:
            return self.fail(e, "Error searching genomes")


class GenomeDownloadCommand(BaseCommand):
    """Command to download genome assemblies from NCBI."""

    @property
    def name(self) -> str:
        return "genome_download"

    @property
    def help(self) -> str:
        return "Download genome assemblies from NCBI by accession, species, or genus"

    @property
    def group(self) -> str:
        return "Genomes"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--accessions",
            nargs="+",
            help="One or more assembly accessions (e.g. GCF_000006945.2)",
        )
        parser.add_argument(
            "--accession-file",
            help="Text file with one accession per line",
        )
        parser.add_argument(
            "--species",
            help="Species name to search GTDB and download genomes",
        )
        parser.add_argument(
            "--genus",
            help="Genus name to search GTDB and download genomes",
        )
        parser.add_argument(
            "--output-dir",
            default="genomes/",
            help="Directory to save downloaded genomes",
        )
        parser.add_argument(
            "--all",
            dest="representative_only",
            action="store_false",
            default=True,
            help="Include all genomes, not just representatives",
        )
        parser.add_argument(
            "--assembly-level",
            choices=["complete", "chromosome", "scaffold", "contig"],
            default=None,
            help="Filter by assembly level",
        )
        parser.add_argument(
            "--force",
            action="store_true",
            help="Re-download genomes even if their FASTA already exists",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Report which genomes are already present vs would be downloaded, without downloading",
        )
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")

    def _collect_accessions(self, args: argparse.Namespace) -> list:
        """Collect accessions from all input sources."""
        accessions = []

        if args.accessions:
            accessions.extend(args.accessions)

        if args.accession_file:
            path = Path(args.accession_file)
            if not path.exists():
                raise PrepareError(f"Accession file not found: {path}")
            for line in path.read_text().strip().splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    accessions.append(line)

        if args.species:
            accessions.extend(get_accessions_for_species(args.species, representative_only=args.representative_only))

        if args.genus:
            accessions.extend(get_accessions_for_genus(args.genus, representative_only=args.representative_only))

        return accessions

    def execute(self, args: argparse.Namespace) -> int:
        try:
            if not any([args.accessions, args.accession_file, args.species, args.genus]):
                self.logger.error(
                    "At least one of --accessions, --accession-file, --species, " "or --genus is required"
                )
                return 1

            accessions = self._collect_accessions(args)
            if not accessions:
                self.logger.warning("No accessions to download")
                return 0

            output_dir = Path(args.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            present, missing = partition_present_genomes(accessions, output_dir)
            to_download = accessions if args.force else missing

            if present and not args.force:
                self.logger.info("%d genome(s) already present, skipping", len(present))
            if args.dry_run:
                self.logger.info("Dry run: %d already present, %d would be downloaded", len(present), len(to_download))
                return 0
            if not to_download:
                self.logger.info("All %d genome(s) already present", len(accessions))
                return 0

            require_tools(["datasets"])
            self.logger.info("Downloading %d genome(s) to %s", len(to_download), output_dir)
            zip_path = download_genomes(
                to_download,
                output_dir,
                assembly_level=args.assembly_level,
            )

            genome_paths = extract_and_organize(zip_path, output_dir)
            self.logger.info("Downloaded and extracted %d genome(s)", len(genome_paths))

            if genome_paths:
                with registry_transaction(args.registry) as reg:
                    for genome_id, path in genome_paths.items():
                        # genome_download writes no manifest; keep one an earlier genome_prepare
                        # recorded for this genome rather than blanking it on a --force redownload.
                        earlier = (reg.genomes.get(genome_id) or {}).get("manifest")
                        manifest = resolve_project_path(reg, earlier) if earlier else ""
                        record_genome(reg, genome_id, path, manifest)

            return 0
        except PrepareError as e:
            return self.fail(e, "Error downloading genomes")


class GenomePrepareCommand(BaseCommand):
    """Command to search, download, and prepare a genome manifest."""

    @property
    def name(self) -> str:
        return "genome_prepare"

    @property
    def help(self) -> str:
        return "Search GTDB, download genomes, and create a manifest CSV"

    @property
    def group(self) -> str:
        return "Genomes"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        group = parser.add_mutually_exclusive_group()
        group.add_argument(
            "--species",
            help="Species name to search GTDB",
        )
        group.add_argument(
            "--genus",
            help="Genus name to search GTDB",
        )
        parser.add_argument(
            "--accession-file",
            help="Text file with one accession per line",
        )
        parser.add_argument(
            "--output-dir",
            default="genomes/",
            help="Directory for downloaded genomes",
        )
        parser.add_argument(
            "--all",
            dest="representative_only",
            action="store_false",
            default=True,
            help="Include all genomes, not just representatives",
        )
        parser.add_argument(
            "--manifest-file",
            default="genome_manifest.csv",
            help="Output manifest file path",
        )
        parser.add_argument(
            "--skip-download",
            action="store_true",
            help="Only create manifest from existing files in output-dir",
        )
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")

    def _collect_accessions(self, args: argparse.Namespace) -> list:
        """Collect accessions from GTDB or accession file."""
        accessions = []

        if args.accession_file:
            path = Path(args.accession_file)
            if not path.exists():
                raise PrepareError(f"Accession file not found: {path}")
            for line in path.read_text().strip().splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    accessions.append(line)

        if args.species:
            accessions.extend(get_accessions_for_species(args.species, representative_only=args.representative_only))

        if args.genus:
            accessions.extend(get_accessions_for_genus(args.genus, representative_only=args.representative_only))

        return accessions

    def _create_manifest(self, output_dir: Path, manifest_file: str, registry: Optional[str] = None) -> int:
        """Create a manifest CSV from the genome FASTA files in output_dir and record its rows."""
        genome_files = visible_files(output_dir, *GENOME_FASTA_GLOBS)
        if not genome_files:
            self.logger.warning("No genome files found in %s", output_dir)
            return 0

        manifest_path = Path(manifest_file)
        rows = []
        with open_atomic(manifest_path, newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["name", "genome_filename", "protein_filename"])
            for gf in genome_files:
                name = _genome_name(gf.name)
                protein = next(
                    (p for p in (gf.with_name(name + ".faa"), gf.with_name(name + ".faa.gz")) if p.exists()),
                    None,
                )
                writer.writerow([name, str(gf), str(protein) if protein else ""])
                rows.append((name, gf))

        if rows:
            with registry_transaction(registry) as reg:
                for name, gf in rows:
                    record_genome(reg, name, gf, manifest_path)

        self.logger.info("Created manifest with %d entries: %s", len(genome_files), manifest_path)
        return len(genome_files)

    def execute(self, args: argparse.Namespace) -> int:
        try:
            output_dir = Path(args.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            if not args.skip_download:
                if not any([args.species, args.genus, args.accession_file]):
                    self.logger.error(
                        "At least one of --species, --genus, or --accession-file "
                        "is required when not using --skip-download"
                    )
                    return 1

                accessions = self._collect_accessions(args)
                if not accessions:
                    self.logger.warning("No accessions found")
                    return 0

                require_tools(["datasets"])
                self.logger.info("Downloading %d genome(s) to %s", len(accessions), output_dir)

                zip_path = download_genomes(accessions, output_dir)
                extract_and_organize(zip_path, output_dir)

            self._create_manifest(output_dir, args.manifest_file, args.registry)
            return 0
        except PrepareError as e:
            return self.fail(e, "Error preparing genomes")
