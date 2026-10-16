"""
SRA metadata and statistics handling for PREPARE.

This module provides comprehensive functionality for fetching SRA metadata,
detecting sequencing technologies, and writing the dataset statistics table of sra_profile.
"""

import json
import logging
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

import pandas as pd
import requests

from prepare.core.exceptions import DataAccessError, NetworkError
from prepare.data.file_io import write_csv
from prepare.data.ncbi_pacing import RequestPacer
from prepare.utils.http import retrying_session

logger = logging.getLogger(__name__)

_RUN_ACCESSION_RE = re.compile(r"^[SED]RR\d+$", re.IGNORECASE)


@dataclass
class SRADatasetInfo:
    """Information about an SRA dataset."""

    accession: str
    title: str
    organism: str
    platform: str
    instrument: str
    strategy: str
    layout: str
    spots: int
    bases: int
    avg_length: float
    size_mb: float
    release_date: str
    bioproject: str
    biosample: str
    library_selection: str
    library_source: str


def _is_transient(error: "requests.RequestException") -> bool:
    """True for a failure worth retrying later: no connection, a timeout, HTTP 429 or a 5xx answer."""
    if isinstance(error, (requests.ConnectionError, requests.Timeout)):
        return True
    status = getattr(getattr(error, "response", None), "status_code", None)
    return isinstance(status, int) and (status == 429 or status >= 500)


class SRAMetadataClient:
    """Client for fetching SRA metadata from NCBI."""

    def __init__(self, email: str, api_key: Optional[str] = None):
        """
        Initialize SRA metadata client.

        Args:
            email: Email address for NCBI API access
            api_key: Optional API key for increased rate limits
        """
        self.email = email
        self.api_key = api_key
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.pacer = RequestPacer(api_key)
        # One retrying session per client, reused across every esearch/efetch call it makes: a
        # connection failure, a 429 or a 5xx from NCBI is retried with backoff before
        # _make_request ever sees it; a NetworkError below is raised only once retries (and the
        # transport-level connection attempts they cover) are exhausted.
        self.session = retrying_session()

    def close(self) -> None:
        """Close the client's HTTP session and its pooled connections."""
        self.session.close()

    def __enter__(self) -> "SRAMetadataClient":
        """Return the client itself; the session is closed when the block exits."""
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        """Close the session."""
        self.close()

    def _make_request(self, url: str, params: Dict[str, str]) -> str:
        """Make rate-limited request to NCBI API."""
        self.pacer.wait()

        # Add required parameters
        params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key

        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            return response.text
        except requests.RequestException as e:
            logger.error(f"NCBI API request failed: {e}")
            if _is_transient(e):
                raise NetworkError(f"Failed to query NCBI: {e}") from e
            raise DataAccessError(f"Failed to query NCBI: {e}") from e

    def get_sra_metadata(self, accessions: List[str]) -> Dict[str, SRADatasetInfo]:
        """
        Fetch metadata for SRA accessions.

        Args:
            accessions: List of SRA accessions

        Returns:
            Dictionary mapping accessions to SRADatasetInfo objects

        Raises:
            NetworkError: a batch could not be reached (connection failure, timeout, or a
                persistent 429/5xx after retries). This is deliberately not caught here and
                continued past: it is a subclass of DataAccessError (below), but the caller
                needs to tell "NCBI is unreachable right now" (worth retrying the whole run
                later, see the ``ExitCode`` docs) apart from "this one batch's data was bad"
                (worth skipping and continuing with the rest).
        """
        if not accessions:
            return {}

        logger.info(f"Fetching metadata for {len(accessions)} SRA accessions")
        results = {}
        batch_size = 200  # Process in batches to avoid URL length limits

        for i in range(0, len(accessions), batch_size):
            batch = accessions[i : i + batch_size]
            batch_num = i // batch_size + 1
            total_batches = (len(accessions) + batch_size - 1) // batch_size
            logger.info(f"Processing batch {batch_num}/{total_batches}")

            try:
                batch_results = self._fetch_batch_metadata(batch)
                results.update(batch_results)
            except NetworkError:
                # A retryable failure, not a bad batch: let it propagate so the caller (and
                # ultimately the CLI's exit code) sees it, instead of silently returning
                # whatever partial results were collected so far as if nothing went wrong.
                raise
            except (DataAccessError, json.JSONDecodeError, KeyError) as e:
                logger.error(f"Failed to fetch metadata for batch: {e}")
                # Continue with other batches
                continue

        logger.info(f"Successfully fetched metadata for {len(results)} accessions")
        return results

    def _fetch_batch_metadata(self, accessions: List[str]) -> Dict[str, SRADatasetInfo]:
        """Fetch metadata for a batch of accessions."""
        # Search for accessions
        search_url = f"{self.base_url}esearch.fcgi"
        search_query = " OR ".join(f"{acc}[ACCN]" for acc in accessions)
        params = {
            "db": "sra",
            "term": search_query,
            "retmax": str(len(accessions) * 2),
            "retmode": "json",
        }

        search_response = self._make_request(search_url, params)
        search_data = json.loads(search_response)

        if "esearchresult" not in search_data or not search_data["esearchresult"]["idlist"]:
            logger.warning(f"No results found for batch: {accessions[:3]}...")
            return {}

        # Fetch detailed information
        ids = search_data["esearchresult"]["idlist"]
        fetch_url = f"{self.base_url}efetch.fcgi"
        params = {
            "db": "sra",
            "id": ",".join(ids),
            "retmode": "xml",
        }

        fetch_response = self._make_request(fetch_url, params)
        # Restrict the result to packages this batch actually asked for: eSearch can match an
        # accession at any level (run, experiment, study, sample), and the matching efetch
        # package for one requested run can bundle other runs (e.g. other lanes/replicates of
        # the same experiment) alongside it, which are kept too, not dropped.
        results = self._parse_sra_xml(fetch_response, requested=set(accessions))
        present = {acc.upper() for acc in results}
        missing = [acc for acc in accessions if _RUN_ACCESSION_RE.match(acc) and acc.upper() not in present]
        if missing:
            logger.warning(
                "No metadata returned for %d requested run accession(s): %s", len(missing), ", ".join(missing)
            )
        return results

    def _parse_sra_xml(self, xml_content: str, requested: Optional[Set[str]] = None) -> Dict[str, SRADatasetInfo]:
        """Parse SRA XML response to extract metadata, one entry per RUN accession.

        ``requested``, when given, is the set of accessions the caller actually asked for.
        Filtering happens per EXPERIMENT_PACKAGE, not per RUN: a package is kept (every RUN
        in it returned) when any requested accession matches, case-insensitively, that
        package's EXPERIMENT, STUDY or SAMPLE accession, any of its RUN accessions, or an
        identifier listed under one of those elements (a BioProject or BioSample accession); a
        package matching none of those is dropped. When a non-empty reply would be filtered
        down to nothing, every package is kept and a WARNING is logged, so the caller sees the
        runs that were returned rather than a false report that nothing could be fetched. A request list may reasonably
        hold an accession from any of those levels (nothing about how it is built rules out
        an experiment, study or sample accession alongside RUN accessions), and one package
        can bundle several RUNs (e.g. other lanes/replicates of the same experiment) that a
        RUN-only match would otherwise have dropped even though the package was asked for.
        Called directly with no ``requested`` set (a script, a REPL, or a test working with
        raw XML), every RUN in every package is returned, matching the historical behaviour.
        """
        try:
            # xml.etree (expat), not prepare.utils.xml's lxml parser: ElementTree never resolves
            # an external entity or fetches a DTD, so XXE does not apply, and expat 2.4.1 and later
            # rejects entity amplification such as billion laughs. This is NCBI's own API response.
            root = ET.fromstring(xml_content)
        except ET.ParseError as e:
            logger.error(f"Failed to parse SRA XML: {e}")
            return {}

        results = {}
        requested_upper = {r.upper() for r in requested} if requested is not None else None

        packages = root.findall(".//EXPERIMENT_PACKAGE")

        matched = []
        uninspected = 0
        for package in packages:
            try:
                if requested_upper is None or self._package_matches_requested(package, requested_upper):
                    matched.append(package)
            except ValueError as e:
                uninspected += 1
                logger.warning(f"Failed to inspect dataset package: {e}")
        if requested_upper is not None and packages and not matched:
            if uninspected:
                # At least one package could not be checked at all: falling back to "keep
                # everything" here would resurrect a package that was never confirmed to
                # match, so only the (empty) set of cleanly-inspected matches is kept.
                logger.warning("%d package(s) could not be inspected for the requested accessions", uninspected)
            else:
                logger.warning("requested accessions matched no package in the reply; listing every run returned")
                matched = packages

        # _extract_dataset_info logs and skips a package whose values cannot be read.
        for package in matched:
            for info in self._extract_dataset_info(package):
                results[info.accession] = info

        return results

    @staticmethod
    def _package_matches_requested(package, requested_upper: Set[str]) -> bool:
        """True when any accession this EXPERIMENT_PACKAGE carries is in ``requested_upper``
        (already uppercased).

        The accessions compared are the ``accession`` attribute of its EXPERIMENT, STUDY and
        SAMPLE and of each RUN, plus the text of every IDENTIFIERS/PRIMARY_ID, EXTERNAL_ID and
        SECONDARY_ID under those elements, which is where a BioProject (PRJNA...) or BioSample
        (SAMN...) accession appears. Comparison is case-insensitive on this side too, since an
        accession's own casing in the XML is not guaranteed to match how a caller wrote it.
        """
        elements = [package.find(".//EXPERIMENT"), package.find(".//STUDY"), package.find(".//SAMPLE")]
        elements.extend(package.findall(".//RUN_SET/RUN"))
        candidates = []
        for element in elements:
            if element is None:
                continue
            candidates.append(element.get("accession", ""))
            for tag in ("PRIMARY_ID", "EXTERNAL_ID", "SECONDARY_ID"):
                candidates.extend((node.text or "").strip() for node in element.findall(f"IDENTIFIERS/{tag}"))
        return any(candidate and candidate.upper() in requested_upper for candidate in candidates)

    def _extract_platform(self, experiment) -> Tuple[str, str]:
        """Return (platform, instrument) from the first PLATFORM child of an experiment."""
        platform_elem = experiment.find(".//PLATFORM")
        if platform_elem is None:
            return "", ""
        for child in platform_elem:
            instrument_elem = child.find(".//INSTRUMENT_MODEL")
            instrument = (instrument_elem.text or "") if instrument_elem is not None else ""
            return child.tag, instrument
        return "", ""

    @staticmethod
    def _run_numbers(run) -> Tuple[int, int, float, str]:
        """spots, bases, size in MB and published date of one ``<RUN>`` element.

        Spots/bases are normally the RUN's own ``total_spots``/``total_bases`` attributes.
        When those are missing (0), the RUN's nested ``<Statistics nspots="..."
        nbases="...">`` child, where present, is used instead, rather than reporting a
        dataset that has real reads as having none.
        """

        def _int_attr(element, name: str) -> int:
            value = element.get(name) if element is not None else None
            return int(value) if value and value.isdigit() else 0

        spots = _int_attr(run, "total_spots")
        bases = _int_attr(run, "total_bases")
        if spots == 0 or bases == 0:
            statistics_elem = run.find("./Statistics")
            if statistics_elem is not None:
                if spots == 0:
                    spots = _int_attr(statistics_elem, "nspots")
                if bases == 0:
                    bases = _int_attr(statistics_elem, "nbases")

        size_bytes = _int_attr(run, "size")
        return spots, bases, size_bytes / (1024 * 1024), run.get("published", "") or ""

    def _extract_biosample(self, package) -> str:
        """Return the BioSample accession from SAMPLE_ATTRIBUTE tags, or ''."""
        for attr in package.findall(".//SAMPLE_ATTRIBUTE"):
            if self._get_text(attr, ".//TAG", "").lower() == "biosample":
                return self._get_text(attr, ".//VALUE", "")
        return ""

    def _extract_dataset_info(self, package) -> List[SRADatasetInfo]:
        """Extract dataset information from an XML package, one entry per RUN.

        A package with no RUN_SET/RUN elements still yields a single entry keyed by
        the experiment accession, with zeroed run-level numbers.
        """
        try:
            experiment = package.find(".//EXPERIMENT")
            if experiment is None:
                return []

            accession = experiment.get("accession", "")
            title = self._get_text(experiment, ".//TITLE", "")
            platform, instrument = self._extract_platform(experiment)

            # Get library info
            library_descriptor = experiment.find(".//LIBRARY_DESCRIPTOR")
            strategy = self._get_text(library_descriptor, ".//LIBRARY_STRATEGY", "")
            selection = self._get_text(library_descriptor, ".//LIBRARY_SELECTION", "")
            source = self._get_text(library_descriptor, ".//LIBRARY_SOURCE", "")
            layout_elem = library_descriptor.find(".//LIBRARY_LAYOUT") if library_descriptor is not None else None
            layout = "PAIRED" if layout_elem is not None and layout_elem.find(".//PAIRED") is not None else "SINGLE"

            # Get sample / study / submission info
            sample = package.find(".//SAMPLE")
            organism = self._get_text(sample, ".//SCIENTIFIC_NAME", "") if sample is not None else ""
            study = package.find(".//STUDY")
            bioproject = (
                self._get_text(study, ".//EXTERNAL_ID[@namespace='BioProject']", "") if study is not None else ""
            )
            submission = package.find(".//SUBMISSION")
            submission_received = submission.get("received", "") if submission is not None else ""
            biosample = self._extract_biosample(package)

            runs = package.findall(".//RUN_SET/RUN")
            infos = []
            for run in runs or [None]:
                if run is None:
                    run_accession, spots, bases, size_mb, published = accession, 0, 0, 0.0, ""
                else:
                    run_accession = run.get("accession", "") or accession
                    spots, bases, size_mb, published = self._run_numbers(run)
                infos.append(
                    SRADatasetInfo(
                        accession=run_accession,
                        title=title,
                        organism=organism,
                        platform=platform,
                        instrument=instrument,
                        strategy=strategy,
                        layout=layout,
                        spots=spots,
                        bases=bases,
                        avg_length=bases / spots if spots else 0.0,
                        size_mb=size_mb,
                        release_date=published or submission_received,
                        bioproject=bioproject,
                        biosample=biosample,
                        library_selection=selection,
                        library_source=source,
                    )
                )
            return infos

        except ValueError as e:
            logger.warning(f"Failed to extract dataset info: {e}")
            return []

    def _get_text(self, element, xpath: str, default: str = "") -> str:
        """Safely extract text from XML element."""
        if element is None:
            return default

        try:
            if xpath.startswith("./@"):
                # Attribute
                attr_name = xpath[3:]
                return element.get(attr_name, default)
            elif "/@" in xpath:
                # Nested attribute
                elem_path, attr_name = xpath.rsplit("/@", 1)
                elem = element.find(elem_path)
                return elem.get(attr_name, default) if elem is not None else default
            else:
                # Element text
                elem = element.find(xpath)
                return elem.text or default if elem is not None else default
        except SyntaxError:
            # ElementPath rejects a path it cannot evaluate with SyntaxError.
            return default


def detect_sequencing_technology(dataset_info: SRADatasetInfo) -> str:
    """
    Detect sequencing technology from SRA metadata.

    Args:
        dataset_info: SRA dataset information

    Returns:
        Technology type: 'illumina', 'nanopore', 'pacbio', or 'unknown'
    """
    platform = dataset_info.platform.lower()
    instrument = dataset_info.instrument.lower()

    # Illumina detection
    if platform == "illumina" or "illumina" in instrument:
        return "illumina"

    # Nanopore detection
    if platform == "oxford_nanopore" or "nanopore" in instrument or "minion" in instrument or "gridion" in instrument:
        return "nanopore"

    # PacBio detection
    if platform == "pacbio_smrt" or "pacbio" in instrument or "sequel" in instrument or "rs ii" in instrument:
        return "pacbio"

    # Check strategy for additional hints
    strategy = dataset_info.strategy.lower()
    if "nanopore" in strategy:
        return "nanopore"
    if "pacbio" in strategy:
        return "pacbio"

    logger.warning(f"Unknown sequencing technology for {dataset_info.accession}: " f"{platform}/{instrument}")
    return "unknown"


def create_download_preview(
    accessions: List[str], metadata_client: SRAMetadataClient
) -> Tuple[Dict[str, SRADatasetInfo], Dict[str, int], float]:
    """
    Create a preview of what would be downloaded.

    Args:
        accessions: List of SRA accessions
        metadata_client: SRA metadata client

    Returns:
        Tuple of (metadata_dict, technology_counts, total_size_gb)
    """
    logger.info("Creating download preview...")

    # Fetch metadata
    metadata = metadata_client.get_sra_metadata(accessions)

    # Count technologies
    tech_counts: dict = {}
    total_size_mb = 0.0

    for acc, info in metadata.items():
        tech = detect_sequencing_technology(info)
        tech_counts[tech] = tech_counts.get(tech, 0) + 1
        total_size_mb += info.size_mb

    total_size_gb = total_size_mb / 1024

    logger.info(f"Preview: {len(metadata)} datasets, {total_size_gb:.2f} GB total")
    return metadata, tech_counts, total_size_gb


def save_metadata_report(metadata: Dict[str, SRADatasetInfo], output_file: Union[str, Path]) -> None:
    """
    Save metadata report to CSV file.

    Args:
        metadata: Dictionary of SRA metadata
        output_file: Output CSV file path
    """
    if not metadata:
        logger.warning("No metadata to save")
        return

    # Convert to DataFrame
    records = []
    for acc, info in metadata.items():
        records.append(
            {
                "accession": info.accession,
                "title": info.title,
                "organism": info.organism,
                "platform": info.platform,
                "instrument": info.instrument,
                "strategy": info.strategy,
                "layout": info.layout,
                "spots": info.spots,
                "bases": info.bases,
                "avg_length": info.avg_length,
                "size_mb": info.size_mb,
                "release_date": info.release_date,
                "bioproject": info.bioproject,
                "biosample": info.biosample,
                "library_selection": info.library_selection,
                "library_source": info.library_source,
                "technology": detect_sequencing_technology(info),
            }
        )

    df = pd.DataFrame(records)
    write_csv(df, output_file, index=False)
    logger.info(f"Metadata report saved to {output_file}")


def _resolved_sidecar_path(acc_dir: Path) -> Optional[Path]:
    """Sidecar next to ``acc_dir``'s resolved target, when ``acc_dir`` is a store link.

    Returns None for a plain directory (nothing to resolve to) or when the resolved
    directory holds no sidecar of its own.
    """
    if not acc_dir.is_symlink():
        return None
    resolved = acc_dir.resolve()
    candidate = resolved / f"{resolved.name}.json"
    return candidate if candidate.is_file() else None


def format_statistics_summary(df: pd.DataFrame) -> List[str]:
    """The aggregate statistics and layout distribution for a statistics table, as text lines.

    The caller decides where the lines go (``sra_profile`` writes them to stdout). Read totals
    count mates (both ends of a pair), not NCBI spots.
    """
    sampled = " (lower bound: some totals are sample counts)" if bool(df["sampled"].any()) else ""
    lines = [
        "\nDataset Statistics Summary:",
        "==========================",
        f"Total datasets: {len(df)}",
        f"Total reads (mates counted): {df['total_reads'].sum():,}{sampled}",
        f"Total bases: {df['total_bases'].sum():,}",
        f"Average read length: {df['avg_read_length'].mean():.1f}",
        f"Average GC content: {df['gc_percent'].mean():.1f}%",
        "\nLayout distribution:",
    ]
    lines.extend(f"  {layout}: {count}" for layout, count in df["layout"].value_counts().items())
    return lines


def generate_statistics_report(rows: Sequence[Dict[str, Any]], output_file: Union[str, Path]) -> List[str]:
    """Write the ``sra_profile`` statistics table, one row per profiled accession.

    ``rows`` come from ``prepare.sra.profiles.statistics_row``; GC is in percent under
    ``gc_percent``. Returns the summary lines from ``format_statistics_summary`` for the
    caller to show, or an empty list (and no file) when there are no rows.
    """
    if not rows:
        logger.warning("No statistics to write")
        return []
    df = pd.DataFrame(list(rows))
    write_csv(df, output_file, index=False)
    logger.info(f"Statistics report saved to {output_file} ({len(df)} datasets)")
    return format_statistics_summary(df)


def estimate_download_time(total_size_gb: float, bandwidth_mbps: float = 100, num_parallel: int = 4) -> float:
    """
    Estimate download time based on size and bandwidth.

    Args:
        total_size_gb: Total size in GB
        bandwidth_mbps: Bandwidth in Mbps
        num_parallel: Number of parallel downloads

    Returns:
        Estimated time in hours
    """
    # Convert GB to Mb
    total_size_mb = total_size_gb * 1024 * 8

    # Account for parallel downloads (with some efficiency loss)
    effective_bandwidth = bandwidth_mbps * num_parallel * 0.8

    # Calculate time in seconds, convert to hours
    time_seconds = total_size_mb / effective_bandwidth
    time_hours = time_seconds / 3600

    return time_hours
