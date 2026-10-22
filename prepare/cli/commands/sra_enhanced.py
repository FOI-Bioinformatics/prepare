"""
Enhanced SRA CLI commands for PREPARE.

This module provides the sra_info and sra_validate commands for previewing NCBI metadata
and validating downloaded datasets. Dataset statistics and quality profiles are
``sra_profile`` (``prepare.cli.commands.sra_profile``).
"""

import gzip
import logging
from functools import partial
from pathlib import Path
from typing import Any, Dict, Optional

from prepare.cli.base import (
    BaseCommand,
    accessions_from_args,
    finite_float,
    anchored_folder,
    read_accessions_file,
    resolve_command_store,
)
from prepare.core.exceptions import DataAccessError, PrepareError
from prepare.core.settings import require_email, setting_for
from prepare.core.validation import validate_accession
from prepare.data import registry_blocks as rb
from prepare.data.file_io import visible_files
from prepare.data.registry import Registry, load_registry, record_analysis
from prepare.data.registry_batch import registry_batch
from prepare.data.sra.cleanup import is_transient_folder
from prepare.data.sra import (
    MATE1_SUFFIXES,
    MATE_SUFFIXES,
    count_fastq_reads,
    fastq_files,
    fastq_stem,
)
from prepare.data.sra_metadata import (
    SRAMetadataClient,
    _resolved_sidecar_path,
    create_download_preview,
    estimate_download_time,
    save_metadata_report,
)
from prepare.store.sidecar import md5_file, read_sidecar
from prepare.store.stats import cached_stats
from prepare.store.usage import record_usage_many

logger = logging.getLogger(__name__)

# Suffixes marking the second mate of a pair, i.e. MATE_SUFFIXES minus MATE1_SUFFIXES.
_MATE2_SUFFIXES = tuple(suffix for suffix in MATE_SUFFIXES if suffix not in MATE1_SUFFIXES)


class SRAInfoCommand(BaseCommand):
    """Command for getting SRA dataset information before downloading."""

    @property
    def name(self) -> str:
        return "sra_info"

    @property
    def help(self) -> str:
        return "Get detailed information about SRA datasets before downloading"

    @property
    def group(self) -> str:
        return "Reads"

    def configure_parser(self, parser):
        parser.add_argument(
            "--accessions-file",
            required=True,
            help="File containing SRA accessions, one per line",
        )
        parser.add_argument(
            "--email",
            default=None,
            help="Email address for NCBI API access (default: PREPARE_NCBI_EMAIL or config [runtime] ncbi_email)",
        )
        parser.add_argument(
            "--api-key",
            help="NCBI API key for increased rate limits (default: PREPARE_NCBI_API_KEY or NCBI_API_KEY)",
        )
        parser.add_argument(
            "--output-report",
            default="sra_info_report.csv",
            help="Output file for detailed report",
        )
        parser.add_argument(
            "--bandwidth-mbps",
            type=finite_float,
            default=100.0,
            help="Estimated bandwidth in Mbps for download time estimation",
        )

    def _print_analysis_summary(self, accessions, metadata, tech_counts, total_size_gb, bandwidth_mbps):
        """Print the SRA dataset analysis summary (counts, distributions, size, ETA)."""
        self.emit("\nSRA Dataset Analysis:")
        self.emit("===================")
        self.emit(f"Total accessions: {len(accessions)}")
        self.emit(f"Metadata fetched: {len(metadata)}")
        self.emit(f"Total estimated size: {total_size_gb:.2f} GB")

        if tech_counts:
            self.emit("\nTechnology distribution:")
            for tech, count in tech_counts.items():
                self.emit(f"  {tech}: {count} datasets")

        platforms: dict = {}
        layouts: dict = {}
        for info in metadata.values():
            platforms[info.platform] = platforms.get(info.platform, 0) + 1
            layouts[info.layout] = layouts.get(info.layout, 0) + 1

        if platforms:
            self.emit("\nPlatform distribution:")
            for platform, count in platforms.items():
                self.emit(f"  {platform}: {count}")
        if layouts:
            self.emit("\nLayout distribution:")
            for layout, count in layouts.items():
                self.emit(f"  {layout}: {count}")

        sizes = [info.size_mb / 1024 for info in metadata.values()]  # Convert to GB
        if sizes:
            self.emit("\nSize statistics:")
            self.emit(f"  Average size per dataset: {sum(sizes)/len(sizes):.2f} GB")
            self.emit(f"  Largest dataset: {max(sizes):.2f} GB")
            self.emit(f"  Smallest dataset: {min(sizes):.2f} GB")

        estimated_hours = estimate_download_time(total_size_gb, bandwidth_mbps, 4)
        if estimated_hours < 1:
            self.emit(f"  Estimated download time: {estimated_hours*60:.0f} minutes")
        else:
            self.emit(f"  Estimated download time: {estimated_hours:.1f} hours")

    def execute(self, args):
        email = require_email(args)
        try:
            accessions = read_accessions_file(args.accessions_file)

            if not accessions:
                self.logger.error("No accessions found in file")
                return 1

            self.emit(f"Analyzing {len(accessions)} SRA accessions...")

            client = SRAMetadataClient(email, setting_for(args, "ncbi_api_key"))
            try:
                metadata, tech_counts, total_size_gb = create_download_preview(accessions, client)
            finally:
                client.close()

            if not metadata:
                self.logger.error("Could not fetch metadata for any accessions")
                return 1

            self._print_analysis_summary(accessions, metadata, tech_counts, total_size_gb, args.bandwidth_mbps)

            save_metadata_report(metadata, args.output_report)
            self.emit(f"\nDetailed report saved to: {args.output_report}")

            return 0

        except (PrepareError, OSError) as e:
            return self.fail(e, "SRA info command failed")


class SRAValidateCommand(BaseCommand):
    """Command for validating downloaded SRA datasets."""

    @property
    def name(self) -> str:
        return "sra_validate"

    @property
    def help(self) -> str:
        return "Validate integrity of downloaded SRA datasets"

    @property
    def group(self) -> str:
        return "Reads"

    def configure_parser(self, parser):
        parser.add_argument(
            "--fastq-folder",
            default="fastq",
            help="Folder containing downloaded FASTQ files",
        )
        parser.add_argument("--accessions-file", default=None, help="Validate the accessions listed here")
        parser.add_argument(
            "--accession",
            action="append",
            default=None,
            help="Validate this accession (repeatable). With neither this nor --accessions-file, every "
            "accession folder in --fastq-folder is validated",
        )
        parser.add_argument(
            "--check-pairs",
            action="store_true",
            help="Check that paired-end files have matching read counts",
        )
        parser.add_argument(
            "--md5",
            action="store_true",
            help="Re-hash each FASTQ file and compare it with the md5 the store recorded for "
            "that file when it was stored, which detects a file changed or corrupted since "
            "(no-op for a dataset without a store sidecar)",
        )
        parser.add_argument("--registry", default=None, help="Registry file (default: found upwards from here)")
        parser.add_argument("--data-root", default=None, help="Shared data store root (overrides discovery)")

    def _find_accession_dirs(self, fastq_folder, specific_accessions=None):
        """Find accession directories to validate.

        Only folders named like an accession count; a download-in-progress folder (``<acc>_temp``,
        ``.prepare-tmp``) or any other folder under the FASTQ root is neither validated nor recorded.
        """
        accession_dirs = [
            d
            for d in visible_files(fastq_folder, dirs=True)
            if validate_accession(d.name) and not is_transient_folder(d.name)
        ]
        if specific_accessions:
            accession_dirs = [d for d in accession_dirs if d.name in specific_accessions]
        return accession_dirs

    def _missing_results(self, wanted, accession_dirs):
        """A FAILED result for every named accession that has no folder, each logged as an error."""
        found = {d.name for d in accession_dirs}
        results = []
        for accession in wanted or ():
            if accession in found:
                continue
            self.logger.error("%s: no folder under the FASTQ folder; nothing to validate", accession)
            results.append(
                {"accession": accession, "status": "FAILED", "issues": "no folder", "issues_list": ["no folder"]}
            )
        return results

    @staticmethod
    def _empty_file_issues(fastq_files) -> list:
        """Issues for any zero-byte FASTQ files."""
        return [f"Empty file: {f.name}" for f in fastq_files if f.stat().st_size == 0]

    @staticmethod
    def _first_record_issue(path: Path) -> Optional[str]:
        """The problem with ``path``'s first FASTQ record, or None when it looks right.

        Opens the file once (gzip aware) and reads four lines: the '@' header, the sequence,
        the '+' separator and a quality string of the same length. A large corrupted file is
        never parsed beyond that.
        """
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as handle:
            header = handle.readline()
            if not header:
                return f"No valid FASTQ records in {path.name}"
            if not header.startswith("@"):
                return f"FASTQ format error in {path.name}: header does not start with '@'"
            seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
        if not seq or not plus or not qual:
            return f"FASTQ format error in {path.name}: the first record is incomplete"
        if not plus.startswith("+"):
            return f"FASTQ format error in {path.name}: the third line does not start with '+'"
        if len(seq.rstrip("\r\n")) != len(qual.rstrip("\r\n")):
            return f"FASTQ format error in {path.name}: sequence/quality length mismatch"
        return None

    @staticmethod
    def _fastq_format_issues(acc_dir: Path) -> list:
        """Issue for any FASTQ file in ``acc_dir`` whose first record is malformed.

        Every file is checked (gzip aware, via ``prepare.data.sra.fastq_files``), not only
        the first, since a download can leave one good mate and one broken one.
        """
        issues = []
        for f in fastq_files(acc_dir):
            try:
                issue = SRAValidateCommand._first_record_issue(f)
            except (ValueError, OSError) as e:
                issues.append(f"FASTQ format error in {f.name}: {e}")
                continue
            if issue is not None:
                issues.append(issue)
        return issues

    @staticmethod
    def _mate_count_issues(acc_dir: Path, cached: Optional[Dict[str, Any]]) -> list:
        """Issue when a paired-end dataset's two mate files have different read counts.

        Read counts come from a cached stats record's ``reads_per_file`` (already computed,
        no file I/O) when available, else a fresh ``count_fastq_reads`` per mate file. A
        dataset with no complete mate-1/mate-2 pair (single-end, or an incomplete pair) is
        not flagged here.
        """
        files = fastq_files(acc_dir)
        mate1 = next((f for f in files if fastq_stem(f).endswith(MATE1_SUFFIXES)), None)
        mate2 = next((f for f in files if fastq_stem(f).endswith(_MATE2_SUFFIXES)), None)
        if mate1 is None or mate2 is None:
            return []
        reads_per_file = (cached or {}).get("reads_per_file") or {}
        n1 = reads_per_file.get(mate1.name)
        if n1 is None:
            n1 = count_fastq_reads(mate1)
        n2 = reads_per_file.get(mate2.name)
        if n2 is None:
            n2 = count_fastq_reads(mate2)
        if n1 != n2:
            return [f"mate files differ ({n1} vs {n2})"]
        return []

    @staticmethod
    def _completeness_issues(acc_dir: Path, verdict: Optional[rb.Verdict]) -> list:
        """Issue when this accession's download did not complete against NCBI's spot count.

        Prefers the store sidecar (freshest, when ``acc_dir`` is a store link) over the
        registry's own download verdict (``"truncated"`` in its own vocabulary). Silent when
        neither source has ever verified this accession against NCBI.

        Three sidecar states are reported: ``"partial"`` (fewer reads on disk than NCBI's
        spot count), ``"failed"`` (the download did not finish) and ``"downloading"`` (a
        download is in progress, so the files on disk are not the finished dataset). Only
        ``"complete"`` and ``"adopted"`` datasets pass.
        """
        sidecar_path = _resolved_sidecar_path(acc_dir)
        if sidecar_path is not None:
            sidecar = read_sidecar(sidecar_path)
            if sidecar is None:
                return []
            if sidecar.state == "partial":
                reads = sidecar.reads_per_mate
                spots = sidecar.ncbi.get("spots")
                return [f"partial: {reads} reads on disk vs {spots} spots at NCBI"]
            if sidecar.state == "failed":
                detail = sidecar.error or "see store_verify"
                return [f"store state failed: {detail}"]
            if sidecar.state == "downloading":
                return ["download in progress elsewhere"]
            return []

        if verdict is not None and verdict.verdict == "truncated":
            reads = verdict.reads_r1
            spots = verdict.expected_spots
            return [f"partial: {reads} reads on disk vs {spots} spots at NCBI"]
        return []

    @staticmethod
    def _md5_issues(acc_dir: Path) -> list:
        """Issue for any FASTQ file whose md5 no longer matches the one the store recorded.

        The comparison is against the sidecar's own ``files[].md5``, computed over the stored
        FASTQ file when it was downloaded or adopted, so it detects a file changed or
        corrupted since. It is not NCBI's md5, which covers the ``.sra`` archive rather than
        the FASTQ files extracted from it and so can never match one. A no-op when there is
        no sidecar to compare against."""
        sidecar_path = _resolved_sidecar_path(acc_dir)
        if sidecar_path is None:
            return []
        sidecar = read_sidecar(sidecar_path)
        if sidecar is None:
            return []
        recorded = {entry.get("name"): entry.get("md5") for entry in sidecar.files}
        issues = []
        for f in fastq_files(acc_dir):
            expected = recorded.get(f.name)
            if expected is None:
                continue
            if md5_file(f) != expected:
                issues.append(f"md5 mismatch: {f.name}")
        return issues

    def _validate_directory(
        self,
        acc_dir,
        registry: Optional[Registry] = None,
        check_pairs: bool = False,
        check_md5: bool = False,
    ):
        """Validate a single accession directory."""
        self.emit(f"Validating {acc_dir.name}...")

        raw_files = visible_files(acc_dir, "*.fastq*")
        if not raw_files:
            return {
                "accession": acc_dir.name,
                "status": "FAILED",
                "issues": "No FASTQ files found",
                "issues_list": ["No FASTQ files found"],
                "num_files": 0,
                "checks": [],
            }

        checks = ["empty_files", "format"]
        issues = self._empty_file_issues(raw_files)
        issues += self._fastq_format_issues(acc_dir)

        if check_pairs:
            # Only the mate-count check reads the statistics record, so a run without
            # --check-pairs does not stat the files or read the sidecar for nothing.
            checks.append("mate_counts")
            issues += self._mate_count_issues(acc_dir, cached_stats(acc_dir, _resolved_sidecar_path(acc_dir)))

        checks.append("completeness")
        verdict = rb.download_verdict(registry, acc_dir.name) if registry is not None else None
        issues += self._completeness_issues(acc_dir, verdict)

        if check_md5:
            checks.append("md5")
            issues += self._md5_issues(acc_dir)

        return {
            "accession": acc_dir.name,
            "status": "PASSED" if not issues else "FAILED",
            "issues": "; ".join(issues) if issues else "None",
            "issues_list": issues,
            "num_files": len(raw_files),
            "checks": checks,
        }

    def _print_validation_results(self, validation_results):
        """Print validation results summary."""
        self.emit("\nValidation Results:")
        self.emit("=================")

        passed = [r for r in validation_results if r["status"] == "PASSED"]
        failed = [r for r in validation_results if r["status"] == "FAILED"]

        self.emit(f"Total validated: {len(validation_results)}")
        self.emit(f"Passed: {len(passed)}")
        self.emit(f"Failed: {len(failed)}")

        if failed:
            self.emit("\nFailed validations:")
            for result in failed:
                self.emit(f"  {result['accession']}: {result['issues']}")

        return len(failed) == 0

    def execute(self, args):
        try:
            args.fastq_folder = str(anchored_folder(args.fastq_folder, getattr(args, "registry", None)))
            fastq_folder = Path(args.fastq_folder)
            if not fastq_folder.exists():
                self.logger.error(f"FASTQ folder {fastq_folder} does not exist")
                return 1

            self.emit("Validating downloaded SRA datasets...")

            wanted = accessions_from_args(getattr(args, "accessions_file", None), getattr(args, "accession", None))
            accession_dirs = self._find_accession_dirs(fastq_folder, wanted)
            missing = self._missing_results(wanted, accession_dirs)
            if not accession_dirs and not missing:
                raise DataAccessError(f"No accession folder under {fastq_folder}")

            # A snapshot, read without the lock: validation reads every file (and md5s it with
            # --md5), which can take minutes. The results are recorded afterwards in one
            # transaction on the registry as it is then, and the catalogue after that lock is released.
            registry = load_registry(args.registry)
            store = resolve_command_store(args, registry)
            check_md5 = getattr(args, "md5", False)
            validation_results = [
                self._validate_directory(acc_dir, registry, args.check_pairs, check_md5) for acc_dir in accession_dirs
            ]
            with registry_batch(args.registry, flush_every=None, flush_seconds=None) as batch:
                for result in validation_results:  # a named accession with no folder is not recorded
                    summary = {
                        "passed": result["status"] == "PASSED",
                        "files": result.get("num_files", 0),
                        "issues": result.get("issues_list", []),
                    }
                    batch.apply(
                        partial(
                            record_analysis,
                            accession=result["accession"],
                            analysis="validate",
                            output="",
                            summary=summary,
                        ),
                        result["accession"],
                    )
            if batch.registry is not None:
                rows = [(r["accession"], "", "analysed", "validate") for r in validation_results]
                record_usage_many(store, batch.registry, rows)

            success = self._print_validation_results(validation_results + missing)
            return 0 if success else 1

        except (PrepareError, OSError) as e:
            # A registry lock wait that gave up (4) or a configuration problem (3) keeps its code.
            return self.fail(e, "SRA validation failed")
