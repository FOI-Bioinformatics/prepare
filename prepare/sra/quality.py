"""Per-read quality figures from a sample of a dataset's FASTQ files.

``SequenceQualityAnalyzer.analyze_fastq_quality`` samples reads from every mate file of a
dataset (uniformly over the whole of each file by default, with
``prepare.data.sra.sample_records``) and returns their length, GC (in percent), base
quality, N content, complexity, duplication and adapter figures. Base qualities are kept as a
histogram of Phred scores rather than a per-base list. The dataset totals come from the
shared statistics record instead; see
``prepare.sra.analytics.SRADatasetAnalyzer.profile_dataset_quality``.
"""

import itertools
import logging
import statistics
import zlib
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from prepare.core.exceptions import DataAccessError
from prepare.data.sra import sample_records
from prepare.data.sra.sampling import distribution_from_histogram, histogram_from_scores, quality_histogram

logger = logging.getLogger(__name__)

# What reading a FASTQ file raises for a file that is unreadable, truncated or malformed (a
# corrupt gzip stream raises zlib.error, which is not an OSError); anything else is a bug.
_FASTQ_READ_ERRORS = (OSError, EOFError, ValueError, UnicodeDecodeError, zlib.error)

# One sample: per-read lengths, GC fractions, the Phred score histogram of all bases, per-read
# N fractions, and the sequences.
Sample = Tuple[List[int], List[float], np.ndarray, List[float], List[str]]

_QUALITY_CLASSES = ("excellent_q30+", "good_q20-29", "fair_q10-19", "poor_q0-9")


def _gc_histogram(gc_contents: List[float]) -> Dict[str, int]:
    """Bucket per-read GC fractions into 5-percent-wide bins, e.g. ``{"40-45": 12, ...}``."""
    histogram: Dict[str, int] = defaultdict(int)
    for gc in gc_contents:
        start = min(95, max(0, int((gc * 100) // 5) * 5))
        histogram[f"{start}-{start + 5}"] += 1
    return dict(histogram)


class SequenceQualityAnalyzer:
    """Analyzes sequence quality metrics from FASTQ files."""

    def __init__(self):
        self.quality_encodings = {"sanger": 33, "illumina_1.3": 64, "illumina_1.5": 64, "solexa": 64}

    def analyze_fastq_quality(
        self,
        fastq_path: Union[str, Path, Sequence[Union[str, Path]]],
        sample_size: int = 10000,
        sampler: str = "uniform",
        total_records: Optional[int] = None,
    ) -> Dict[str, Any]:
        """
        Analyze quality metrics from the FASTQ file(s) of one dataset.

        Args:
            fastq_path: One FASTQ file, or every mate file of a dataset; the sample is then
                drawn from all of them, so mate 2 is represented as well as mate 1
            sample_size: Number of reads to sample for analysis, across all files
            sampler: ``"uniform"`` (default) samples reads uniformly across the whole of
                every file (``prepare.data.sra.sample_records``), so reads past the first
                ``sample_size`` are represented too, not just the head. ``"head"`` takes the
                first reads of each file only (an equal share of ``sample_size`` per file),
                for a caller that wants the cheapest possible read of files it already knows
                to be homogeneous.
            total_records: Reads across all the files, e.g. from the dataset's cached
                statistics record; lets the uniform sampler skip its own counting pass.

        Returns:
            Dictionary of quality metrics; GC figures are in percent (0-100)
        """
        paths = [Path(fastq_path)] if isinstance(fastq_path, (str, Path)) else [Path(p) for p in fastq_path]
        for path in paths:
            if not path.exists():
                raise DataAccessError(f"FASTQ file not found: {path}")

        try:
            sample = self._sample(paths, sample_size, sampler, total_records)
        except _FASTQ_READ_ERRORS as e:
            logger.error(f"Error analyzing FASTQ file {paths[0] if paths else ''}: {e}")
            raise DataAccessError(f"Failed to analyze FASTQ file: {e}")

        if not sample[0]:
            raise DataAccessError("No valid reads found in FASTQ file")
        return self._summarise(*sample)

    def _sample(self, paths: List[Path], sample_size: int, sampler: str, total_records: Optional[int] = None) -> Sample:
        """Per-read lengths, GC fractions, base quality histogram, N fractions and sequences of the sample."""
        if sampler != "head":
            return self._sample_uniform(paths, sample_size, total_records)
        share = -(-sample_size // max(len(paths), 1))
        lengths: List[int] = []
        gc_fractions: List[float] = []
        quality_hist = histogram_from_scores([])
        n_fractions: List[float] = []
        sequences: List[str] = []
        for path in paths:
            file_lengths, file_gc, file_hist, file_n, file_sequences = self._sample_head(path, share)
            lengths.extend(file_lengths)
            gc_fractions.extend(file_gc)
            quality_hist = quality_hist + file_hist
            n_fractions.extend(file_n)
            sequences.extend(file_sequences)
        return lengths, gc_fractions, quality_hist, n_fractions, sequences

    def _summarise(
        self,
        read_lengths: List[int],
        gc_contents: List[float],
        quality_hist: np.ndarray,
        n_contents: List[float],
        sequences: List[str],
    ) -> Dict[str, Any]:
        """The quality metrics of one sample; GC figures in percent."""
        gc_percents = [gc * 100 for gc in gc_contents]
        quality = distribution_from_histogram(quality_hist)
        return {
            "total_reads_sampled": len(read_lengths),
            "read_length_stats": {
                "mean": statistics.mean(read_lengths),
                "median": statistics.median(read_lengths),
                "std": statistics.stdev(read_lengths) if len(read_lengths) > 1 else 0,
                "min": min(read_lengths),
                "max": max(read_lengths),
                "distribution": self._get_length_distribution(read_lengths),
            },
            "gc_content_stats": {
                "mean": statistics.mean(gc_percents),
                "median": statistics.median(gc_percents),
                "std": statistics.stdev(gc_percents) if len(gc_percents) > 1 else 0,
                # A histogram, not the raw per-read list: the list used to make every
                # profile JSON grow with the dataset instead of staying a fixed size.
                "histogram": _gc_histogram(gc_contents),
            },
            "quality_stats": {
                "mean": quality["mean"],
                "median": quality["median"],
                "q25": quality["q25"],
                "q75": quality["q75"],
                "distribution": {key: quality[key] for key in _QUALITY_CLASSES} if quality_hist.sum() else {},
            },
            "n_content_stats": {
                "mean": statistics.mean(n_contents),
                "max": max(n_contents),
                "high_n_reads": sum(1 for n in n_contents if n > 0.05),
            },
            "complexity_metrics": self._analyze_sequence_complexity(sequences),
            "contamination_indicators": self._detect_contamination_indicators(sequences),
            "duplication_rate": self._calculate_duplication_rate(sequences),
        }

    def _sample_head(self, fastq_path: Path, sample_size: int) -> Sample:
        """First ``sample_size`` reads of ``fastq_path`` (the original, head-only sampling).

        Quality strings are counted into a histogram (``quality_histogram``) rather than
        decoded into one Phred score per base.
        """
        import gzip

        read_lengths: List[int] = []
        gc_contents: List[float] = []
        qualities: List[bytes] = []
        n_contents: List[float] = []
        sequences: List[str] = []

        opener = gzip.open if fastq_path.suffix.endswith(".gz") else open
        with opener(fastq_path, "rt") as handle:
            for _title, sequence, quality in itertools.islice(FastqGeneralIterator(handle), sample_size):
                read_lengths.append(len(sequence))
                gc_contents.append(self._calculate_gc_content(sequence))
                qualities.append(quality.encode("ascii"))
                n_contents.append(sequence.count("N") / len(sequence) if sequence else 0.0)
                sequences.append(sequence)

        return read_lengths, gc_contents, quality_histogram(qualities), n_contents, sequences

    def _sample_uniform(
        self,
        fastq_path: Union[Path, Sequence[Path]],
        sample_size: int,
        total_records: Optional[int] = None,
    ) -> Sample:
        """Sample ``sample_size`` reads uniformly over the whole of every file.

        Unlike ``_sample_head``, a read from anywhere in the files has an equal chance of
        being included, so a dataset whose later reads (or second mates) differ from its
        first ones is still represented in the quality metrics. The record indices are drawn
        up front with a fixed seed (``prepare.data.sra.sample_records``), and the files are
        streamed once in binary blocks; ``total_records``, when known, spares the counting
        pass that the draw otherwise needs.
        """
        paths = [fastq_path] if isinstance(fastq_path, (str, Path)) else list(fastq_path)
        records = sample_records(paths, sample_size, total_records=total_records, seed=0)

        read_lengths: List[int] = []
        gc_contents: List[float] = []
        n_contents: List[float] = []
        sequences: List[str] = []
        for seq_bytes, _qual in records:
            seq = seq_bytes.decode("utf-8")
            read_lengths.append(len(seq))
            gc_contents.append(self._calculate_gc_content(seq))
            n_contents.append((seq.count("N") / len(seq)) if seq else 0.0)
            sequences.append(seq)

        return read_lengths, gc_contents, quality_histogram(q for _s, q in records), n_contents, sequences

    def _calculate_duplication_rate(self, sequences: List[str]) -> float:
        """Fraction of sampled reads that are exact duplicates of another read.

        Computed over the sampled sequences as ``1 - unique / total``; 0.0 when
        every read is distinct and rising toward 1.0 as duplicates dominate.
        """
        if not sequences:
            return 0.0
        return 1.0 - len(set(sequences)) / len(sequences)

    def _calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content of sequence."""
        gc_count = sequence.count("G") + sequence.count("C")
        return gc_count / len(sequence) if len(sequence) > 0 else 0.0

    def _get_length_distribution(self, lengths: List[int]) -> Dict[str, int]:
        """Get read length distribution in ranges."""
        distribution = {"0-50": 0, "51-100": 0, "101-150": 0, "151-250": 0, "251-500": 0, "501-1000": 0, "1000+": 0}

        for length in lengths:
            if length <= 50:
                distribution["0-50"] += 1
            elif length <= 100:
                distribution["51-100"] += 1
            elif length <= 150:
                distribution["101-150"] += 1
            elif length <= 250:
                distribution["151-250"] += 1
            elif length <= 500:
                distribution["251-500"] += 1
            elif length <= 1000:
                distribution["501-1000"] += 1
            else:
                distribution["1000+"] += 1

        return distribution

    def _get_quality_distribution(self, quality_scores: List[int]) -> Dict[str, float]:
        """Fractions of ``quality_scores`` at Q30+, Q20-29, Q10-19 and below Q10; empty for no scores."""
        if not quality_scores:
            return {}
        quality = distribution_from_histogram(histogram_from_scores(quality_scores))
        return {key: quality[key] for key in _QUALITY_CLASSES}

    def _analyze_sequence_complexity(self, sequences: List[str]) -> Dict[str, float]:
        """Analyze sequence complexity indicators."""
        if not sequences:
            return {}

        # Calculate various complexity metrics
        kmer_diversities = []
        homopolymer_rates = []

        for seq in sequences[:1000]:  # Sample first 1000 sequences
            # K-mer diversity (using 3-mers)
            kmers = [seq[i : i + 3] for i in range(len(seq) - 2)]
            unique_kmers = len(set(kmers))
            total_kmers = len(kmers)
            kmer_diversities.append(unique_kmers / total_kmers if total_kmers > 0 else 0)

            # Homopolymer run detection
            homopolymer_count = 0
            current_run = 1
            for i in range(1, len(seq)):
                if seq[i] == seq[i - 1]:
                    current_run += 1
                else:
                    if current_run >= 4:  # Runs of 4+ same base
                        homopolymer_count += 1
                    current_run = 1
            if current_run >= 4:  # a run that reaches the end of the read
                homopolymer_count += 1
            homopolymer_rates.append(homopolymer_count / len(seq) if len(seq) > 0 else 0)

        return {
            "kmer_diversity_mean": statistics.mean(kmer_diversities) if kmer_diversities else 0,
            "homopolymer_rate_mean": statistics.mean(homopolymer_rates) if homopolymer_rates else 0,
            "complexity_score": statistics.mean(kmer_diversities) if kmer_diversities else 0,
        }

    def _detect_contamination_indicators(self, sequences: List[str]) -> Dict[str, float]:
        """Detect potential contamination indicators."""
        # Simple contamination indicators
        indicators = {"adapter_contamination": 0.0, "vector_contamination": 0.0, "overrepresented_sequences": 0.0}

        if not sequences:
            return indicators

        # Common adapter sequences
        adapters = [
            "AGATCGGAAGAGC",  # Illumina universal adapter
            "CTGTCTCTTATACACATCT",  # Illumina adapter
            "AATGATACGGCGACCACCGAGATCTACAC",  # Illumina P5
        ]

        adapter_hits = 0
        for seq in sequences:
            for adapter in adapters:
                if adapter in seq:
                    adapter_hits += 1
                    break

        indicators["adapter_contamination"] = adapter_hits / len(sequences)

        # Check for overrepresented sequences
        sequence_counts = Counter(sequences)
        most_common = sequence_counts.most_common(10)
        if most_common:
            max_count = most_common[0][1]
            indicators["overrepresented_sequences"] = max_count / len(sequences)

        return indicators
