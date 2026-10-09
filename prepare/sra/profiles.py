"""Per-accession quality profiles on disk, and the statistics table row built from one.

``sra_profile`` writes one ``<accession>_quality_profile.json`` per accession into its
``--output-dir``; ``sra_report --quality-profiles`` reads the same files back with
``load_quality_profiles`` instead of profiling those accessions again.
"""

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

from prepare.data.file_io import visible_files, write_text_atomic
from prepare.data.sra import MATE_SUFFIXES, MATE1_SUFFIXES, fastq_stem
from prepare.sra.analytics import QualityProfile, SRADatasetAnalyzer, _gc_histogram
from prepare.sra.dataset_stats import load_dataset_stats
from prepare.store.stats import DEFAULT_SAMPLE_SIZE

logger = logging.getLogger(__name__)

PROFILE_JSON_SUFFIX = "_quality_profile.json"

# Suffixes marking the second mate of a pair, i.e. MATE_SUFFIXES minus MATE1_SUFFIXES.
_MATE2_SUFFIXES = tuple(suffix for suffix in MATE_SUFFIXES if suffix not in MATE1_SUFFIXES)


def profile_json_path(accession: str, output_dir: Union[str, Path]) -> Path:
    """Where ``accession``'s quality profile JSON lives under ``output_dir``."""
    return Path(output_dir) / f"{accession}{PROFILE_JSON_SUFFIX}"


def profile_to_dict(profile: QualityProfile) -> Dict[str, Any]:
    """The JSON payload of one quality profile. GC is in percent under ``gc_percent``."""
    return {
        "accession": profile.accession,
        "total_reads": profile.total_reads,
        "reads_sampled": profile.reads_sampled,
        "sampled": profile.sampled,
        "total_bases": profile.total_bases,
        "avg_read_length": profile.avg_read_length,
        "read_length_distribution": profile.read_length_distribution,
        "gc_percent": profile.gc_percent,
        "gc_histogram": profile.gc_histogram,
        "mean_quality": profile.mean_quality,
        "quality_grade": profile.quality_grade,
        "quality_distribution": profile.quality_distribution,
        "complexity_score": profile.complexity_score,
        "n_content": profile.n_content,
        "duplication_rate": profile.duplication_rate,
        "technology_confidence": profile.technology_confidence,
        "contamination_indicators": profile.contamination_indicators,
        "recommendations": profile.recommendations,
    }


def write_profile_json(profile: QualityProfile, output_dir: Union[str, Path]) -> Path:
    """Write ``profile`` as JSON into ``output_dir`` and return the file's path."""
    path = profile_json_path(profile.accession, output_dir)
    write_text_atomic(path, json.dumps(profile_to_dict(profile), indent=2))
    return path


def _gc_percent(data: Dict[str, Any]) -> float:
    """GC in percent from a profile JSON; one written before 0.5.0 holds a 0-1 ``gc_content``."""
    if "gc_percent" in data:
        return float(data["gc_percent"] or 0.0)
    return float(data.get("gc_content") or 0.0) * 100


def _profile_from_dict(accession: str, data: Dict[str, Any]) -> QualityProfile:
    """A QualityProfile from one profile JSON's payload, old or new."""
    return QualityProfile(
        accession=accession,
        total_reads=data.get("total_reads", 0),
        total_bases=data.get("total_bases", 0),
        avg_read_length=data.get("avg_read_length", 0.0),
        read_length_distribution=data.get("read_length_distribution", {}),
        gc_percent=_gc_percent(data),
        gc_histogram=data.get("gc_histogram") or _gc_histogram(data.get("gc_distribution", [])),
        quality_distribution=data.get("quality_distribution", {}),
        n_content=data.get("n_content", 0.0),
        contamination_indicators=data.get("contamination_indicators", {}),
        complexity_score=data.get("complexity_score", data.get("sequence_complexity", 0.0)),
        duplication_rate=data.get("duplication_rate"),
        technology_confidence=data.get("technology_confidence", 0.0),
        quality_grade=data.get("quality_grade", ""),
        recommendations=data.get("recommendations", []),
        reads_sampled=data.get("reads_sampled", data.get("total_reads", 0)),
        sampled=data.get("sampled", "reads_sampled" not in data),
        mean_quality=data.get("mean_quality", 0.0),
    )


def load_quality_profiles(profiles_dir: Union[str, Path]) -> Dict[str, QualityProfile]:
    """Load previously saved per-accession quality profile JSONs.

    Returns an accession -> QualityProfile mapping for every readable
    ``*_quality_profile.json`` file found directly under ``profiles_dir``. A missing
    directory yields an empty mapping; a file that is not valid JSON is skipped with a
    warning rather than raising.

    Files written before 0.5.0 are read too: their ``gc_content`` fraction becomes
    ``gc_percent``; a per-read ``gc_distribution`` list is bucketed into the 5-percent-wide
    ``gc_histogram``; a ``sequence_complexity`` key stands in for ``complexity_score``; and a
    file without ``reads_sampled`` (whose ``total_reads`` was the sample size) loads with
    ``reads_sampled`` equal to ``total_reads`` and ``sampled`` True.
    """
    profiles: Dict[str, QualityProfile] = {}
    for path in visible_files(Path(profiles_dir), f"*{PROFILE_JSON_SUFFIX}"):
        try:
            data = json.loads(path.read_text())
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as e:
            logger.warning("Could not read quality profile %s: %s", path, e)
            continue
        accession = data.get("accession") or path.name[: -len(PROFILE_JSON_SUFFIX)]
        profiles[accession] = _profile_from_dict(accession, data)
    return profiles


def _layout(files: Sequence[Path]) -> str:
    """``PAIRED`` when a second-mate file is among ``files``, else ``SINGLE``."""
    return "PAIRED" if any(fastq_stem(f).endswith(_MATE2_SUFFIXES) for f in files) else "SINGLE"


def statistics_row(
    profile: QualityProfile, files: Sequence[Path], dataset_stats: Optional[Dict[str, Any]]
) -> Dict[str, Any]:
    """One row of the ``sra_profile`` statistics table for ``profile``.

    Totals, mean length and GC are the profile's (which come from the shared statistics
    record); minimum and maximum length and N50 are the record's own figures, None when it
    could not be computed.
    """
    stats = dataset_stats or {}
    return {
        "accession": profile.accession,
        "num_files": len(files),
        "layout": _layout(files),
        "total_reads": profile.total_reads,
        "total_bases": profile.total_bases,
        "avg_read_length": profile.avg_read_length,
        "min_read_length": stats.get("min_read_length"),
        "max_read_length": stats.get("max_read_length"),
        "n50": stats.get("n50"),
        "gc_percent": profile.gc_percent,
        "avg_quality": profile.mean_quality,
        "complexity_score": profile.complexity_score,
        "duplication_rate": profile.duplication_rate,
        "quality_grade": profile.quality_grade,
        "reads_sampled": profile.reads_sampled,
        "sampled": profile.sampled,
    }


@dataclass
class ProfiledDataset:
    """One accession's quality profile with the files and statistics record it was built from."""

    profile: QualityProfile
    files: List[Path]
    stats: Optional[Dict[str, Any]]


def profile_accession(
    analyzer: SRADatasetAnalyzer, accession: str, sample_size: int = DEFAULT_SAMPLE_SIZE, sampler: str = "uniform"
) -> Optional[ProfiledDataset]:
    """Profile ``accession`` once from its FASTQ files under the analyzer's folder; None without files.

    The statistics record is loaded first (cached in the store sidecar, or computed) and
    handed to ``profile_dataset_quality``, so the files are not read for it a second time.
    A dataset whose record cannot be computed is still profiled, with its totals taken
    from the sample and flagged as such. Read errors from the files propagate as
    ``DATASET_READ_ERRORS``.
    """
    files = analyzer.find_fastq_files(accession)
    if not files:
        return None
    stats = load_dataset_stats(files, sample_size=sample_size)
    profile = analyzer.profile_dataset_quality(
        accession,
        fastq_path=files,
        sample_size=sample_size,
        sampler=sampler,
        dataset_stats=stats if stats is not None else {},
    )
    return ProfiledDataset(profile=profile, files=files, stats=stats)
