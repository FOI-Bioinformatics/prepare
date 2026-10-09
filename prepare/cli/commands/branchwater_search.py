"""CLI command that searches the Branchwater index with a genome and writes a Branchwater CSV."""

import argparse
from pathlib import Path

from prepare.cli.base import BaseCommand
from prepare.core.constants import DEFAULT_REGISTRY_MAX_SCREENED
from prepare.core.exceptions import PrepareError
from prepare.data.branchwater_search import (
    DEFAULT_SERVER,
    KSIZE,
    load_signature,
    search_index,
    sketch_fasta,
    write_branchwater_csv,
)
from prepare.data.registry import cap_screening, record_screening, registry_transaction


class BranchwaterSearchCommand(BaseCommand):
    """Search Branchwater with a genome FASTA or a sourmash signature."""

    @property
    def name(self) -> str:
        return "branchwater_search"

    @property
    def help(self) -> str:
        return "Search the Branchwater index of SRA metagenomes with a genome and write a Branchwater CSV"

    @property
    def group(self) -> str:
        return "Containment"

    def configure_parser(self, parser: argparse.ArgumentParser) -> None:
        source = parser.add_mutually_exclusive_group(required=True)
        source.add_argument("--genome-fasta", help="Genome FASTA to sketch (needs the sourmash extra)")
        source.add_argument("--signature", help="Existing sourmash signature file (k=21, scaled=1000)")
        parser.add_argument("--threshold", type=float, default=0.1, help="Minimum containment to keep (inclusive)")
        parser.add_argument("--branchwater-folder", default="branchwater", help="Folder for the Branchwater CSV")
        parser.add_argument(
            "--output", default=None, help="Output CSV (default: <branchwater-folder>/<input stem>.csv)"
        )
        parser.add_argument("--server", default=DEFAULT_SERVER, help="Branchwater search API base URL")
        parser.add_argument(
            "--no-cache", action="store_true", help="Do not read or write the Branchwater response cache"
        )
        parser.add_argument("--refresh", action="store_true", help="Bypass a cached result and query the server again")
        parser.add_argument(
            "--max-cache-age-days",
            type=int,
            default=None,
            help="Refetch if the cached result is older than this many days (default: no age limit)",
        )
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")
        parser.add_argument(
            "--registry-max-screened",
            type=int,
            default=DEFAULT_REGISTRY_MAX_SCREENED,
            help="Most screening entries to keep in the registry per genome (the CSV keeps them all)",
        )

    def execute(self, args: argparse.Namespace) -> int:
        try:
            if args.genome_fasta:
                signature = sketch_fasta(args.genome_fasta)
                source = Path(args.genome_fasta)
            else:
                signature = load_signature(args.signature)
                source = Path(args.signature)

            output = Path(args.output) if args.output else Path(args.branchwater_folder) / f"{source.stem}.csv"
            cache_dir = None if args.no_cache else output.parent / ".branchwater-cache"
            matches = search_index(
                signature,
                args.threshold,
                server=args.server,
                cache_dir=cache_dir,
                refresh=args.refresh,
                max_cache_age_days=args.max_cache_age_days,
            )
            write_branchwater_csv(matches, output)

            if matches:
                # Break ties on containment by accession, so the reported "best" hit is
                # deterministic across runs instead of depending on the server's row order.
                best_accession, best_containment = min(matches, key=lambda match: (-match[1], match[0]))
                self.logger.info(
                    "%d metagenome(s) contain %s at >= %.2f; best containment %.4f (%s)",
                    len(matches),
                    source.stem,
                    args.threshold,
                    best_containment,
                    best_accession,
                )
            else:
                self.logger.warning(
                    "No metagenome reached containment %.2f. The public index has returned nothing for known "
                    "positives before; try a control genome (for example Salmonella LT2, GCF_000006945.2) to "
                    "check that the index is answering.",
                    args.threshold,
                )
            self.logger.info("Next: prepare use_branchwater --branchwater-folder %s", output.parent)

            with registry_transaction(args.registry) as registry:
                for accession, containment in matches:
                    cani = containment ** (1 / KSIZE) if containment > 0 else 0.0
                    record_screening(
                        registry, accession, source.stem, containment, cani, "branchwater", args.threshold, output
                    )
                cap_screening(registry, source.stem, args.registry_max_screened)
            return 0
        except PrepareError as e:
            return self.fail(e, "Error searching Branchwater")
