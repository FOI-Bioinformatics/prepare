"""
Comprehensive tests for SRA dataset analytics and quality profiling.

Tests cover:
- Sequence quality analysis from FASTQ files
- Dataset quality profiling and grading
- Comparative analysis across dataset groups
- Anomaly detection in datasets
- Processing parameter recommendations
- Statistical testing and visualization data preparation
"""

from pathlib import Path
from unittest.mock import Mock, patch
import importlib.util
import pytest
import pandas as pd

# Biopython is a hard project dependency (pyproject.toml), always installed alongside
# numpy/scipy -- no need to fake it out at import time. A prior version of this file wrapped
# this import in `with patch.dict("sys.modules", {"Bio": Mock(), ...}):`, which is actively
# harmful: unittest.mock's dict-patching restores sys.modules to an exact snapshot on exit,
# silently unregistering every module imported for the first time inside the block
# (including numpy/scipy, transitively). A later fresh import of anything under prepare.sra
# in the same process then hits numpy's C extensions with "ImportError: cannot load module
# more than once per process". Tests below that read FASTQ records write small files to tmp_path.
from prepare.sra.analytics import (
    SRADatasetAnalyzer,
    SequenceQualityAnalyzer,
    QualityProfile,
    ComparativeAnalysis,
    AnomalyReport,
    ProcessingRecommendations,
)
from prepare.sra.profiles import load_quality_profiles

from prepare.core.exceptions import DataAccessError

requires_analysis = pytest.mark.skipif(
    not all(importlib.util.find_spec(m) for m in ("sklearn", "scipy")), reason="needs prepare[analysis]"
)


class TestSequenceQualityAnalyzer:
    """Test sequence quality analysis functionality."""

    def setup_method(self):
        """Set up test fixtures."""
        self.analyzer = SequenceQualityAnalyzer()

    def test_analyzer_initialization(self):
        """Test SequenceQualityAnalyzer initialization."""
        assert hasattr(self.analyzer, "quality_encodings")
        assert "sanger" in self.analyzer.quality_encodings
        assert self.analyzer.quality_encodings["sanger"] == 33

    def test_calculate_gc_content(self):
        """Test GC content calculation."""
        # Test normal sequence
        gc_content = self.analyzer._calculate_gc_content("ATGCGCAT")
        assert gc_content == 0.5  # 4 G/C out of 8 bases

        # Test all GC
        gc_content_all = self.analyzer._calculate_gc_content("GGGGCCCC")
        assert gc_content_all == 1.0

        # Test no GC
        gc_content_none = self.analyzer._calculate_gc_content("ATATATAT")
        assert gc_content_none == 0.0

        # Test empty sequence
        gc_content_empty = self.analyzer._calculate_gc_content("")
        assert gc_content_empty == 0.0

    def test_get_length_distribution(self):
        """Test read length distribution calculation."""
        lengths = [25, 75, 125, 175, 300, 750, 1500]
        distribution = self.analyzer._get_length_distribution(lengths)

        assert distribution["0-50"] == 1  # 25
        assert distribution["51-100"] == 1  # 75
        assert distribution["101-150"] == 1  # 125
        assert distribution["151-250"] == 1  # 175
        assert distribution["251-500"] == 1  # 300
        assert distribution["501-1000"] == 1  # 750
        assert distribution["1000+"] == 1  # 1500

    def test_get_quality_distribution(self):
        """Test quality score distribution calculation."""
        quality_scores = [5, 15, 25, 35, 8, 22, 32, 12]
        distribution = self.analyzer._get_quality_distribution(quality_scores)

        assert distribution["excellent_q30+"] == 2 / 8  # 35, 32
        assert distribution["good_q20-29"] == 2 / 8  # 25, 22
        assert distribution["fair_q10-19"] == 2 / 8  # 15, 12
        assert distribution["poor_q0-9"] == 2 / 8  # 5, 8

    def test_analyze_sequence_complexity(self):
        """Test sequence complexity analysis."""
        # High complexity sequence
        high_complexity = ["ATGCGTACGTAGC", "CGTACGTAGCTAG", "TACGTAGCTAGCT"]
        complexity_high = self.analyzer._analyze_sequence_complexity(high_complexity)

        assert "kmer_diversity_mean" in complexity_high
        assert "homopolymer_rate_mean" in complexity_high
        assert "complexity_score" in complexity_high

        # Low complexity sequence (homopolymer runs)
        low_complexity = ["AAAAATTTTTGGGG", "CCCCCCAAAAAAA", "TTTTTTTTGGGGGG"]
        complexity_low = self.analyzer._analyze_sequence_complexity(low_complexity)

        assert complexity_low["homopolymer_rate_mean"] > 0

    def test_trailing_homopolymer_run_is_counted(self):
        """A run of four or more bases that ends the read counts like one inside it."""
        for seq in ("ACGTAAAA", "AAAACGTA"):
            result = self.analyzer._analyze_sequence_complexity([seq])
            assert result["homopolymer_rate_mean"] == pytest.approx(1 / 8)

    def test_head_sampler_tolerates_an_empty_read(self, tmp_path):
        """A zero-length record has an N fraction of 0.0 rather than raising."""
        fastq = tmp_path / "reads.fastq"
        fastq.write_text("@empty\n\n+\n\n@r1\nACGT\n+\nIIII\n")
        _lengths, _gc, _hist, n_contents, _seqs = self.analyzer._sample_head(fastq, 10)
        assert n_contents == [0.0, 0.0]

    def test_detect_contamination_indicators(self):
        """Test contamination detection."""
        # Sequences with adapter contamination
        contaminated_sequences = [
            "ATGCAGATCGGAAGAGCGTACG",  # Contains Illumina adapter
            "CGTACGTAGCTAG",
            "TACGTAGCTAGCT",
        ]

        indicators = self.analyzer._detect_contamination_indicators(contaminated_sequences)

        assert "adapter_contamination" in indicators
        assert "vector_contamination" in indicators
        assert "overrepresented_sequences" in indicators
        assert indicators["adapter_contamination"] > 0

    def test_analyze_fastq_quality_success(self, tmp_path):
        """Head sampling of a two-record file: every section present, quality from the histogram."""
        fastq = tmp_path / "test.fastq"
        # Phred+33: '?' is 30, 'D' is 35, 'I' is 40, ':' is 25, '5' is 20.
        fastq.write_text("@seq1\nATGCGTACGT\n+\n?DIDN:?DID\n@seq2\nCGTACGTAGC\n+\n:?D?:5:?D?\n")

        result = self.analyzer.analyze_fastq_quality(fastq, sample_size=100, sampler="head")

        assert "total_reads_sampled" in result
        assert "read_length_stats" in result
        assert "gc_content_stats" in result
        assert "quality_stats" in result
        assert "complexity_metrics" in result
        assert result["total_reads_sampled"] == 2
        assert result["quality_stats"]["median"] == 30.0

    def test_head_sampler_takes_only_the_first_reads(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text("".join(f"@r{i}\n{'G' if i < 3 else 'A'}AAA\n+\nIIII\n" for i in range(10)))

        result = self.analyzer.analyze_fastq_quality(fastq, sample_size=3, sampler="head")

        assert result["total_reads_sampled"] == 3
        assert result["gc_content_stats"]["mean"] == 25.0
        assert result["quality_stats"]["mean"] == 40.0

    def test_analyze_fastq_quality_file_not_found(self):
        """Test FASTQ analysis with missing file."""
        with pytest.raises(DataAccessError, match="FASTQ file not found"):
            self.analyzer.analyze_fastq_quality("/nonexistent/file.fastq")

    @patch("builtins.open", side_effect=IOError("Permission denied"))
    def test_analyze_fastq_quality_io_error(self, mock_file):
        """Test FASTQ analysis with IO error."""
        with patch("pathlib.Path.exists", return_value=True):
            with pytest.raises(DataAccessError, match="Failed to analyze FASTQ file"):
                self.analyzer.analyze_fastq_quality("test.fastq")

    def test_analyze_fastq_quality_lets_a_bug_propagate(self, tmp_path):
        fastq = tmp_path / "reads.fastq"
        fastq.write_text("@r\nACGT\n+\nIIII\n")
        with patch.object(self.analyzer, "_sample_uniform", side_effect=TypeError("bug")):
            with pytest.raises(TypeError):
                self.analyzer.analyze_fastq_quality(fastq)

    def test_analyze_fastq_quality_samples_every_mate(self, tmp_path):
        r1 = tmp_path / "SRR1_1.fastq"
        r2 = tmp_path / "SRR1_2.fastq"
        r1.write_text("@r\nAAAA\n+\nIIII\n")
        r2.write_text("@r\nGGGG\n+\nIIII\n")
        for sampler in ("uniform", "head"):
            result = self.analyzer.analyze_fastq_quality([r1, r2], sample_size=10, sampler=sampler)
            assert result["total_reads_sampled"] == 2
            assert result["gc_content_stats"]["mean"] == 50.0  # percent

    def test_calculate_duplication_rate(self):
        """Duplication rate is 1 - unique/total over the sampled reads."""
        assert self.analyzer._calculate_duplication_rate([]) == 0.0
        assert self.analyzer._calculate_duplication_rate(["A", "C", "G", "T"]) == 0.0
        assert self.analyzer._calculate_duplication_rate(["A", "A", "C", "G"]) == 0.25
        assert self.analyzer._calculate_duplication_rate(["A", "A", "A", "A"]) == 0.75

    def test_analyze_fastq_quality_reports_duplication(self, tmp_path):
        """analyze_fastq_quality includes a real duplication_rate for duplicate reads."""
        fastq = tmp_path / "test.fastq"
        # The second read is identical to the first.
        fastq.write_text("@r1\nATGCATGCAT\n+\n??????????\n@r2\nATGCATGCAT\n+\n??????????\n")

        result = self.analyzer.analyze_fastq_quality(fastq, sample_size=100, sampler="head")

        assert result["duplication_rate"] == 0.5  # one of two reads is a duplicate

    def test_analyze_fastq_quality_uniform_default_holds_histogram_not_per_read_list(self, tmp_path):
        """The default (uniform) sampler's output holds a GC histogram, not a raw per-read list."""
        fastq = tmp_path / "reads.fastq"
        lines = []
        for i in range(20):
            seq = "GGGGGGGGGG" if i % 2 == 0 else "AAAAAAAAAA"
            lines.extend([f"@r{i}", seq, "+", "I" * len(seq)])
        fastq.write_text("\n".join(lines) + "\n")

        result = self.analyzer.analyze_fastq_quality(fastq, sample_size=100)

        assert "histogram" in result["gc_content_stats"]
        assert "distribution" not in result["gc_content_stats"]
        assert isinstance(result["gc_content_stats"]["histogram"], dict)
        assert sum(result["gc_content_stats"]["histogram"].values()) == 20

    def test_analyze_fastq_quality_uniform_samples_tail_of_large_file(self, tmp_path):
        """The uniform sampler must not be biased toward the head of a large file."""
        fastq = tmp_path / "reads.fastq"
        lines = []
        for _ in range(4000):
            lines.extend(["@head", "AAAAAAAA", "+", "IIIIIIII"])
        for _ in range(1000):
            lines.extend(["@tail", "CCCCCCCC", "+", "IIIIIIII"])
        fastq.write_text("\n".join(lines) + "\n")

        result = self.analyzer.analyze_fastq_quality(fastq, sample_size=200, sampler="uniform")

        assert result["total_reads_sampled"] == 200
        # A head-biased sample would show a GC mean of 0.0; the tail is 20% of the file.
        # GC is in percent.
        assert result["gc_content_stats"]["mean"] > 2


class TestLoadQualityProfiles:
    """load_quality_profiles must read both the old (gc_distribution list, gc_content fraction)
    and the new (gc_histogram dict, gc_percent) on-disk JSON shape."""

    def _base_profile_json(self, accession):
        return {
            "accession": accession,
            "total_reads": 100,
            "total_bases": 15000,
            "avg_read_length": 150.0,
            "read_length_distribution": {},
            "gc_content": 0.5,
            "quality_distribution": {},
            "n_content": 0.0,
            "contamination_indicators": {},
            "complexity_score": 0.8,
            "duplication_rate": 0.1,
            "technology_confidence": 0.9,
            "quality_grade": "good",
            "recommendations": [],
        }

    def test_reads_a_profile_written_before_the_sampling_fields_existed(self, tmp_path):
        """An old JSON has neither reads_sampled nor sampled; its total_reads was the sample
        size, so it loads as a sampled profile of that size rather than failing."""
        import json

        data = self._base_profile_json("SRR_OLD")
        (tmp_path / "SRR_OLD_quality_profile.json").write_text(json.dumps(data))

        profile = load_quality_profiles(tmp_path)["SRR_OLD"]

        assert profile.total_reads == 100
        assert profile.reads_sampled == 100
        assert profile.sampled is True

    def test_reads_the_sampling_fields_when_present(self, tmp_path):
        import json

        data = self._base_profile_json("SRR_NEW")
        data["reads_sampled"] = 10
        data["sampled"] = False
        (tmp_path / "SRR_NEW_quality_profile.json").write_text(json.dumps(data))

        profile = load_quality_profiles(tmp_path)["SRR_NEW"]

        assert profile.total_reads == 100
        assert profile.reads_sampled == 10
        assert profile.sampled is False

    def test_reads_old_format_gc_distribution_key(self, tmp_path):
        import json

        data = self._base_profile_json("SRR_OLD")
        data["gc_distribution"] = [0.4, 0.5, 0.6]  # legacy per-read list
        (tmp_path / "SRR_OLD_quality_profile.json").write_text(json.dumps(data))

        profiles = load_quality_profiles(tmp_path)

        assert "SRR_OLD" in profiles
        # Old data is not discarded: the legacy per-read list is bucketed into the same
        # 5-percent-wide histogram the new field holds.
        assert profiles["SRR_OLD"].gc_histogram == {"40-45": 1, "50-55": 1, "60-65": 1}

    def test_old_format_with_no_gc_distribution_defaults_to_empty_histogram(self, tmp_path):
        import json

        data = self._base_profile_json("SRR_OLD_EMPTY")
        # Neither gc_histogram nor gc_distribution present at all.
        (tmp_path / "SRR_OLD_EMPTY_quality_profile.json").write_text(json.dumps(data))

        profiles = load_quality_profiles(tmp_path)

        assert profiles["SRR_OLD_EMPTY"].gc_histogram == {}

    def test_reads_new_format_gc_histogram_key(self, tmp_path):
        import json

        data = self._base_profile_json("SRR_NEW")
        data["gc_histogram"] = {"45-50": 2, "50-55": 1}
        (tmp_path / "SRR_NEW_quality_profile.json").write_text(json.dumps(data))

        profiles = load_quality_profiles(tmp_path)

        assert "SRR_NEW" in profiles
        assert profiles["SRR_NEW"].gc_histogram == {"45-50": 2, "50-55": 1}

    def test_prefers_complexity_score_key_when_both_present(self, tmp_path):
        """complexity_score is the only key written since 0.5.0; when an older profile JSON
        carries both it and the sequence_complexity alias, complexity_score takes priority."""
        import json

        data = self._base_profile_json("SRR_BOTH")
        data["complexity_score"] = 0.8
        data["sequence_complexity"] = 0.5
        (tmp_path / "SRR_BOTH_quality_profile.json").write_text(json.dumps(data))

        profile = load_quality_profiles(tmp_path)["SRR_BOTH"]

        assert profile.complexity_score == 0.8

    def test_reads_sequence_complexity_key_when_complexity_score_absent(self, tmp_path):
        """A profile JSON that only carries the newer sequence_complexity key (no
        complexity_score at all) still loads a real score instead of the 0.0 default."""
        import json

        data = self._base_profile_json("SRR_ALIAS_ONLY")
        del data["complexity_score"]
        data["sequence_complexity"] = 0.63
        (tmp_path / "SRR_ALIAS_ONLY_quality_profile.json").write_text(json.dumps(data))

        profile = load_quality_profiles(tmp_path)["SRR_ALIAS_ONLY"]

        assert profile.complexity_score == 0.63

    def test_reads_an_old_gc_content_fraction_as_percent(self, tmp_path):
        import json

        data = self._base_profile_json("SRR_OLD_GC")  # "gc_content": 0.5, as written before 0.5.0
        (tmp_path / "SRR_OLD_GC_quality_profile.json").write_text(json.dumps(data))

        assert load_quality_profiles(tmp_path)["SRR_OLD_GC"].gc_percent == 50.0

    def test_reads_gc_percent_directly(self, tmp_path):
        import json

        data = self._base_profile_json("SRR_NEW_GC")
        del data["gc_content"]
        data["gc_percent"] = 41.5
        (tmp_path / "SRR_NEW_GC_quality_profile.json").write_text(json.dumps(data))

        assert load_quality_profiles(tmp_path)["SRR_NEW_GC"].gc_percent == 41.5


class TestSRADatasetAnalyzer:
    """Test main SRA dataset analyzer functionality."""

    def setup_method(self):
        """Set up test fixtures."""
        self.analyzer = SRADatasetAnalyzer()

    def test_analyzer_initialization(self):
        """Test SRADatasetAnalyzer initialization."""
        assert hasattr(self.analyzer, "quality_analyzer")
        assert isinstance(self.analyzer.quality_analyzer, SequenceQualityAnalyzer)

    def test_find_fastq_in_per_accession_folder(self, tmp_path):
        acc_dir = tmp_path / "fastq" / "SRR123456"
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR123456_2.fastq").write_text("@r\nA\n+\nI\n")
        r1 = acc_dir / "SRR123456_1.fastq"
        r1.write_text("@r\nA\n+\nI\n")
        analyzer = SRADatasetAnalyzer(fastq_dir=tmp_path / "fastq")
        assert analyzer.find_fastq("SRR123456") == r1

    def test_find_fastq_flat_layout(self, tmp_path):
        flat = tmp_path / "fastq" / "SRR123456.fastq.gz"
        flat.parent.mkdir(parents=True)
        flat.write_text("@r\nA\n+\nI\n")
        analyzer = SRADatasetAnalyzer(fastq_dir=tmp_path / "fastq")
        assert analyzer.find_fastq("SRR123456") == flat

    def test_find_fastq_defaults_to_fastq_folder_in_cwd(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        r1 = tmp_path / "fastq" / "SRR1" / "SRR1_1.fastq"
        r1.parent.mkdir(parents=True)
        r1.write_text("@r\nA\n+\nI\n")
        assert SRADatasetAnalyzer().find_fastq("SRR1") == Path("fastq/SRR1/SRR1_1.fastq")

    def test_find_fastq_missing_returns_none(self, tmp_path):
        assert SRADatasetAnalyzer(fastq_dir=tmp_path).find_fastq("NONE") is None

    def test_find_fastq_ignores_prefix_collisions(self, tmp_path):
        other = tmp_path / "SRR10_1.fastq"
        other.write_text("@r\nA\n+\nI\n")
        (tmp_path / "SRR10").mkdir()
        (tmp_path / "SRR10" / "SRR10_1.fastq").write_text("@r\nA\n+\nI\n")
        analyzer = SRADatasetAnalyzer(fastq_dir=tmp_path)
        assert analyzer.find_fastq("SRR1") is None
        assert analyzer.find_fastq("SRR10") == tmp_path / "SRR10" / "SRR10_1.fastq"

    def test_find_fastq_files_returns_every_mate_and_skips_zero_byte_files(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        for name in ("SRR1_1.fastq", "SRR1_2.fastq"):
            (acc_dir / name).write_text("@r\nA\n+\nI\n")
        (acc_dir / "SRR1_3.fastq").write_text("")  # left by an interrupted download
        analyzer = SRADatasetAnalyzer(fastq_dir=tmp_path)
        assert analyzer.find_fastq_files("SRR1") == [acc_dir / "SRR1_1.fastq", acc_dir / "SRR1_2.fastq"]

    def test_calculate_quality_grade(self):
        """Test quality grade calculation."""
        # Excellent quality
        excellent_grade = self.analyzer._calculate_quality_grade(
            quality_stats={"mean": 35},
            contamination={"adapter_contamination": 0.005},
            complexity={"complexity_score": 0.8},
        )
        assert excellent_grade == "excellent"

        # Poor quality
        poor_grade = self.analyzer._calculate_quality_grade(
            quality_stats={"mean": 15},
            contamination={"adapter_contamination": 0.1},
            complexity={"complexity_score": 0.2},
        )
        assert poor_grade == "poor"

        # Fair quality (score = 0.2 + 0.1 + 0.15 = 0.45)
        fair_grade = self.analyzer._calculate_quality_grade(
            quality_stats={"mean": 25},
            contamination={"adapter_contamination": 0.02},
            complexity={"complexity_score": 0.6},
        )
        assert fair_grade == "fair"

    def test_generate_quality_recommendations(self):
        """Test quality recommendation generation."""
        # Low quality scenario
        recommendations = self.analyzer._generate_quality_recommendations(
            quality_stats={"mean": 15},
            contamination={"adapter_contamination": 0.1},
            complexity={"complexity_score": 0.2},
            metadata=None,
        )

        assert len(recommendations) > 0
        assert any("quality trimming" in rec.lower() for rec in recommendations)
        assert any("adapter" in rec.lower() for rec in recommendations)
        assert any("duplicates" in rec.lower() for rec in recommendations)

        # High quality scenario
        good_recommendations = self.analyzer._generate_quality_recommendations(
            quality_stats={"mean": 35},
            contamination={"adapter_contamination": 0.01},
            complexity={"complexity_score": 0.8},
            metadata=None,
        )

        assert any("high quality" in rec.lower() for rec in good_recommendations)

    def test_estimate_technology_confidence(self):
        """Test technology detection confidence estimation."""
        # Complete metadata
        complete_metadata = Mock()
        complete_metadata.platform = "ILLUMINA"
        complete_metadata.instrument = "HiSeq 2500"

        confidence = self.analyzer._estimate_technology_confidence(complete_metadata)
        assert confidence == 1.0

        # Partial metadata
        partial_metadata = Mock()
        partial_metadata.platform = "ILLUMINA"
        partial_metadata.instrument = None

        confidence_partial = self.analyzer._estimate_technology_confidence(partial_metadata)
        assert 0.5 < confidence_partial < 1.0

        # No metadata
        confidence_none = self.analyzer._estimate_technology_confidence(None)
        assert confidence_none == 0.5

    @patch.object(SequenceQualityAnalyzer, "analyze_fastq_quality")
    def test_profile_dataset_quality(self, mock_analyze):
        """Test dataset quality profiling."""
        # Mock FASTQ analysis results
        mock_analyze.return_value = {
            "total_reads_sampled": 10000,
            "read_length_stats": {"mean": 150, "distribution": {"101-150": 10000}},
            "gc_content_stats": {"mean": 45.0},
            "quality_stats": {"mean": 30, "distribution": {"excellent_q30+": 0.8}},
            "n_content_stats": {"mean": 0.01},
            "complexity_metrics": {"complexity_score": 0.7},
            "contamination_indicators": {"adapter_contamination": 0.02},
            "duplication_rate": 0.3,
        }

        with patch.object(self.analyzer, "find_fastq_files", return_value=[Path("test.fastq")]):
            with patch("pathlib.Path.exists", return_value=True):
                profile = self.analyzer.profile_dataset_quality("SRR123456")

        assert isinstance(profile, QualityProfile)
        assert profile.accession == "SRR123456"
        assert profile.total_reads == 10000
        assert profile.avg_read_length == 150
        assert profile.gc_percent == 45.0
        assert profile.complexity_score == 0.7
        assert profile.duplication_rate == 0.3
        assert profile.quality_grade in ["excellent", "good", "fair", "poor"]
        assert len(profile.recommendations) > 0
        # Without a shared statistics record the totals can only be the sample's own.
        assert profile.reads_sampled == 10000
        assert profile.sampled is True

    @patch.object(SequenceQualityAnalyzer, "analyze_fastq_quality")
    def test_profile_dataset_quality_prefers_the_shared_statistics_record(self, mock_analyze):
        """Given the dataset's statistics record, the profile reports the dataset totals and
        keeps the sample size separately."""
        mock_analyze.return_value = {
            "total_reads_sampled": 10000,
            "read_length_stats": {"mean": 150, "distribution": {}},
            "gc_content_stats": {"mean": 45.0},
            "quality_stats": {"mean": 30, "distribution": {}},
            "n_content_stats": {"mean": 0.01},
            "complexity_metrics": {"complexity_score": 0.7},
            "contamination_indicators": {},
            "duplication_rate": 0.3,
        }
        record = {"reads_total": 1724338, "bases_total": 258650700}

        with patch.object(self.analyzer, "find_fastq_files", return_value=[Path("test.fastq")]):
            with patch("pathlib.Path.exists", return_value=True):
                profile = self.analyzer.profile_dataset_quality("SRR123456", dataset_stats=record)

        assert profile.total_reads == 1724338
        assert profile.total_bases == 258650700
        assert profile.reads_sampled == 10000
        assert profile.sampled is False

    @patch.object(SequenceQualityAnalyzer, "analyze_fastq_quality")
    def test_profile_dataset_quality_scales_bases_when_the_record_has_none(self, mock_analyze):
        """A record from the streaming path carries no exact base total; it is scaled from
        the sampled mean read length instead."""
        mock_analyze.return_value = {
            "total_reads_sampled": 100,
            "read_length_stats": {"mean": 150, "distribution": {}},
            "gc_content_stats": {"mean": 45.0},
            "quality_stats": {"mean": 30, "distribution": {}},
            "n_content_stats": {"mean": 0.0},
            "complexity_metrics": {"complexity_score": 0.7},
            "contamination_indicators": {},
            "duplication_rate": 0.0,
        }

        with patch.object(self.analyzer, "find_fastq_files", return_value=[Path("test.fastq")]):
            with patch("pathlib.Path.exists", return_value=True):
                profile = self.analyzer.profile_dataset_quality(
                    "SRR123456", dataset_stats={"reads_total": 1000, "bases_total": None}
                )

        assert profile.total_reads == 1000
        assert profile.total_bases == 150000
        assert profile.sampled is False

    @requires_analysis
    def test_compare_datasets(self):
        """Test dataset comparison functionality."""
        groups = {"group1": ["SRR123456", "SRR123457"], "group2": ["SRR789012", "SRR789013"]}

        # Mock profile generation
        mock_profiles = {
            "SRR123456": QualityProfile(
                accession="SRR123456",
                total_reads=10000,
                total_bases=1500000,
                avg_read_length=150,
                read_length_distribution={},
                gc_percent=45.0,
                gc_histogram={},
                quality_distribution={},
                n_content=0.01,
                contamination_indicators={"adapter_contamination": 0.02},
                complexity_score=0.7,
                duplication_rate=None,
                technology_confidence=0.9,
                quality_grade="good",
                recommendations=[],
            ),
            "SRR123457": QualityProfile(
                accession="SRR123457",
                total_reads=12000,
                total_bases=1800000,
                avg_read_length=150,
                read_length_distribution={},
                gc_percent=48.0,
                gc_histogram={},
                quality_distribution={},
                n_content=0.015,
                contamination_indicators={"adapter_contamination": 0.025},
                complexity_score=0.75,
                duplication_rate=None,
                technology_confidence=0.85,
                quality_grade="good",
                recommendations=[],
            ),
            "SRR789012": QualityProfile(
                accession="SRR789012",
                total_reads=8000,
                total_bases=1200000,
                avg_read_length=150,
                read_length_distribution={},
                gc_percent=42.0,
                gc_histogram={},
                quality_distribution={},
                n_content=0.008,
                contamination_indicators={"adapter_contamination": 0.03},
                complexity_score=0.65,
                duplication_rate=None,
                technology_confidence=0.88,
                quality_grade="fair",
                recommendations=[],
            ),
            "SRR789013": QualityProfile(
                accession="SRR789013",
                total_reads=9000,
                total_bases=1350000,
                avg_read_length=150,
                read_length_distribution={},
                gc_percent=40.0,
                gc_histogram={},
                quality_distribution={},
                n_content=0.012,
                contamination_indicators={"adapter_contamination": 0.035},
                complexity_score=0.68,
                duplication_rate=None,
                technology_confidence=0.82,
                quality_grade="fair",
                recommendations=[],
            ),
        }

        with patch.object(self.analyzer, "profile_dataset_quality") as mock_profile:

            def profile_side_effect(accession):
                return mock_profiles.get(accession)

            mock_profile.side_effect = profile_side_effect

            comparison = self.analyzer.compare_datasets(groups)

        assert isinstance(comparison, ComparativeAnalysis)
        assert comparison.dataset_groups == groups
        assert "group1" in comparison.summary_statistics
        assert len(comparison.recommendations) > 0

    def test_detect_dataset_anomalies(self):
        """Test anomaly detection functionality."""
        accessions = ["SRR123456", "SRR789012", "SRR555555"]

        # Mock profiles with different quality characteristics
        mock_profiles = {
            "SRR123456": QualityProfile(  # Normal dataset
                accession="SRR123456",
                total_reads=10000,
                total_bases=1500000,
                avg_read_length=150,
                read_length_distribution={},
                gc_percent=45.0,
                gc_histogram={},
                quality_distribution={},
                n_content=0.01,
                contamination_indicators={"adapter_contamination": 0.02},
                complexity_score=0.7,
                duplication_rate=None,
                technology_confidence=0.9,
                quality_grade="good",
                recommendations=[],
            ),
            "SRR789012": QualityProfile(  # Anomalous dataset - high contamination
                accession="SRR789012",
                total_reads=8000,
                total_bases=1200000,
                avg_read_length=150,
                read_length_distribution={},
                gc_percent=15.0,  # Unusual GC
                gc_histogram={},
                quality_distribution={},
                n_content=0.08,  # High N
                contamination_indicators={"adapter_contamination": 0.15},
                complexity_score=0.05,  # Low complexity
                duplication_rate=None,
                technology_confidence=0.6,
                quality_grade="poor",
                recommendations=[],
            ),
            "SRR555555": QualityProfile(  # Normal dataset
                accession="SRR555555",
                total_reads=11000,
                total_bases=1650000,
                avg_read_length=150,
                read_length_distribution={},
                gc_percent=42.0,
                gc_histogram={},
                quality_distribution={},
                n_content=0.012,
                contamination_indicators={"adapter_contamination": 0.018},
                complexity_score=0.72,
                duplication_rate=None,
                technology_confidence=0.88,
                quality_grade="good",
                recommendations=[],
            ),
        }

        with patch.object(self.analyzer, "profile_dataset_quality") as mock_profile:

            def profile_side_effect(accession):
                return mock_profiles.get(accession)

            mock_profile.side_effect = profile_side_effect

            anomaly_report = self.analyzer.detect_dataset_anomalies(accessions)

        assert isinstance(anomaly_report, AnomalyReport)
        assert "SRR789012" in anomaly_report.anomalous_datasets  # Should be flagged
        assert "SRR123456" not in anomaly_report.anomalous_datasets  # Should be normal
        assert len(anomaly_report.anomaly_types) > 0
        assert "SRR789012" in anomaly_report.severity_scores
        assert anomaly_report.severity_scores["SRR789012"] > 0.2

    @requires_analysis
    def test_compare_datasets_with_profiles_does_not_call_profile_dataset_quality(self):
        """Supplying profiles reuses them; profile_dataset_quality must not be called."""
        # Two accessions per group: a group of one has no variance and is not tested.
        groups = {"group1": ["SRR1", "SRR2"], "group2": ["SRR3", "SRR4"]}
        profiles = {
            accession: QualityProfile(
                accession=accession,
                total_reads=10000 - 500 * i,
                total_bases=1500000 - 70000 * i,
                avg_read_length=150 - 3 * i,
                read_length_distribution={},
                gc_percent=45.0 - 2.5 * i,
                gc_histogram={},
                quality_distribution={},
                n_content=0.01,
                contamination_indicators={"adapter_contamination": 0.02},
                complexity_score=0.7 - 0.05 * i,
                duplication_rate=None,
                technology_confidence=0.9,
                quality_grade="good",
                recommendations=[],
            )
            for i, accession in enumerate(["SRR1", "SRR2", "SRR3", "SRR4"])
        }

        with patch.object(self.analyzer, "profile_dataset_quality") as mock_profile:
            comparison = self.analyzer.compare_datasets(groups, profiles=profiles)

        mock_profile.assert_not_called()
        assert isinstance(comparison, ComparativeAnalysis)
        assert "group1" in comparison.summary_statistics

    def test_a_group_of_one_is_not_tested(self, caplog):
        """A t-test with a single dataset in a group gives NaN and scipy warnings; it is skipped."""
        import logging
        import math
        import warnings

        df = pd.DataFrame(
            {
                "accession": ["SRR1", "SRR2", "SRR3"],
                "group": ["A", "A", "B"],
                "avg_read_length": [150.0, 152.0, 149.0],
                "gc_percent": [45.0, 47.0, 52.0],
                "total_reads": [10000, 12000, 8000],
                "complexity_score": [0.7, 0.75, 0.6],
            }
        )
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with caplog.at_level(logging.INFO):
                tests = self.analyzer._perform_statistical_tests(df, {"A": ["SRR1", "SRR2"], "B": ["SRR3"]})
        assert tests["gc_percent"]["test"] == "skipped (a group has fewer than 2 values)"
        assert math.isnan(tests["gc_percent"]["p_value"]) and tests["gc_percent"]["significant"] is False
        assert "a group has fewer than 2 values" in caplog.text

    def test_compare_datasets_profiles_a_group_accession_missing_from_supplied_profiles(self):
        """A partial ``profiles`` dict must not silently drop the accessions it lacks:
        the gap is profiled fresh, and the comparison holds every accession."""
        groups = {"group1": ["SRR1", "SRR2"]}
        profile_srr1 = QualityProfile(
            accession="SRR1",
            total_reads=10000,
            total_bases=1500000,
            avg_read_length=150,
            read_length_distribution={},
            gc_percent=45.0,
            gc_histogram={},
            quality_distribution={},
            n_content=0.01,
            contamination_indicators={"adapter_contamination": 0.02},
            complexity_score=0.7,
            duplication_rate=None,
            technology_confidence=0.9,
            quality_grade="good",
            recommendations=[],
        )
        profile_srr2 = QualityProfile(
            accession="SRR2",
            total_reads=20000,
            total_bases=3000000,
            avg_read_length=150,
            read_length_distribution={},
            gc_percent=50.0,
            gc_histogram={},
            quality_distribution={},
            n_content=0.01,
            contamination_indicators={"adapter_contamination": 0.01},
            complexity_score=0.8,
            duplication_rate=None,
            technology_confidence=0.9,
            quality_grade="excellent",
            recommendations=[],
        )
        # Only SRR1 is supplied; SRR2 is a gap that must still get profiled.
        partial_profiles = {"SRR1": profile_srr1}

        with patch.object(self.analyzer, "profile_dataset_quality", return_value=profile_srr2) as mock_profile:
            comparison = self.analyzer.compare_datasets(groups, profiles=partial_profiles)

        mock_profile.assert_called_once_with("SRR2")
        assert "group1" in comparison.summary_statistics
        # The mean total_reads over the group must reflect BOTH profiles (15000), not just
        # the one that was supplied (10000) -- proof SRR2 was folded into the comparison.
        assert comparison.summary_statistics["group1"]["total_reads"]["mean"] == pytest.approx(15000.0)

    def test_compare_datasets_profiling_failure_for_missing_accession_is_logged_and_excluded(self):
        """An accession that cannot be profiled either way is excluded, not silently
        dropped without a trace: it is logged and noted in recommendations."""
        groups = {"group1": ["SRR1", "SRR_MISSING"]}
        profile_srr1 = QualityProfile(
            accession="SRR1",
            total_reads=10000,
            total_bases=1500000,
            avg_read_length=150,
            read_length_distribution={},
            gc_percent=45.0,
            gc_histogram={},
            quality_distribution={},
            n_content=0.01,
            contamination_indicators={"adapter_contamination": 0.02},
            complexity_score=0.7,
            duplication_rate=None,
            technology_confidence=0.9,
            quality_grade="good",
            recommendations=[],
        )
        partial_profiles = {"SRR1": profile_srr1}

        with patch.object(
            self.analyzer, "profile_dataset_quality", side_effect=DataAccessError("no FASTQ found")
        ) as mock_profile:
            comparison = self.analyzer.compare_datasets(groups, profiles=partial_profiles)

        mock_profile.assert_called_once_with("SRR_MISSING")
        assert any("SRR_MISSING" in note for note in comparison.recommendations)
        # The group's stats still reflect the one profile that was available.
        assert comparison.summary_statistics["group1"]["total_reads"]["mean"] == pytest.approx(10000.0)

    def test_detect_dataset_anomalies_with_profiles_does_not_call_profile_dataset_quality(self):
        """Supplying profiles reuses them; profile_dataset_quality must not be called."""
        accessions = ["SRR1", "SRR2"]
        profiles = {
            "SRR1": QualityProfile(
                accession="SRR1",
                total_reads=10000,
                total_bases=1500000,
                avg_read_length=150,
                read_length_distribution={},
                gc_percent=45.0,
                gc_histogram={},
                quality_distribution={},
                n_content=0.01,
                contamination_indicators={"adapter_contamination": 0.02},
                complexity_score=0.7,
                duplication_rate=None,
                technology_confidence=0.9,
                quality_grade="good",
                recommendations=[],
            ),
            "SRR2": QualityProfile(  # anomalous: high contamination, low complexity
                accession="SRR2",
                total_reads=8000,
                total_bases=1200000,
                avg_read_length=150,
                read_length_distribution={},
                gc_percent=15.0,
                gc_histogram={},
                quality_distribution={},
                n_content=0.08,
                contamination_indicators={"adapter_contamination": 0.15},
                complexity_score=0.05,
                duplication_rate=None,
                technology_confidence=0.6,
                quality_grade="poor",
                recommendations=[],
            ),
        }

        with patch.object(self.analyzer, "profile_dataset_quality") as mock_profile:
            anomaly_report = self.analyzer.detect_dataset_anomalies(accessions, profiles=profiles)

        mock_profile.assert_not_called()
        assert isinstance(anomaly_report, AnomalyReport)
        assert "SRR2" in anomaly_report.anomalous_datasets

    def test_recommend_processing_params(self):
        """Test processing parameter recommendations."""
        # Create test profile
        test_profile = QualityProfile(
            accession="SRR123456",
            total_reads=10000,
            total_bases=1500000,
            avg_read_length=150,
            read_length_distribution={},
            gc_percent=45.0,
            gc_histogram={},
            quality_distribution={},
            n_content=0.01,
            contamination_indicators={"adapter_contamination": 0.08},
            complexity_score=0.7,
            duplication_rate=None,
            technology_confidence=0.9,
            quality_grade="fair",
            recommendations=[],
        )

        recommendations = self.analyzer.recommend_processing_params("SRR123456", test_profile)

        assert isinstance(recommendations, ProcessingRecommendations)
        assert recommendations.accession == "SRR123456"
        assert recommendations.recommended_pipeline in ["short_read", "paired_end", "long_read"]
        assert "enabled" in recommendations.quality_trimming
        assert "enabled" in recommendations.adapter_removal
        assert recommendations.adapter_removal["enabled"] is True  # High contamination
        assert "memory_gb" in recommendations.computational_requirements
        assert recommendations.estimated_processing_time is not None


class TestStatisticalAnalysis:
    """Test statistical analysis functionality."""

    def setup_method(self):
        """Set up test fixtures."""
        self.analyzer = SRADatasetAnalyzer()

    @requires_analysis
    def test_perform_statistical_tests(self):
        """Test statistical test performance."""
        # Create test comparison data
        comparison_df = pd.DataFrame(
            {
                "accession": ["SRR1", "SRR2", "SRR3", "SRR4"],
                "group": ["A", "A", "B", "B"],
                "avg_read_length": [150, 155, 250, 245],
                "gc_percent": [45.0, 47.0, 52.0, 50.0],
                "total_reads": [10000, 12000, 8000, 9000],
                "complexity_score": [0.7, 0.75, 0.6, 0.65],
            }
        )

        groups = {"A": ["SRR1", "SRR2"], "B": ["SRR3", "SRR4"]}

        tests = self.analyzer._perform_statistical_tests(comparison_df, groups)

        assert "avg_read_length" in tests
        assert "gc_percent" in tests
        assert tests["avg_read_length"]["test"] == "t-test"  # Two groups
        assert "statistic" in tests["avg_read_length"]
        assert "p_value" in tests["avg_read_length"]
        assert "significant" in tests["avg_read_length"]

    def test_detect_outliers(self):
        """Test outlier detection."""
        # Create test data with obvious outlier
        comparison_df = pd.DataFrame(
            {
                "accession": ["SRR1", "SRR2", "SRR3", "SRR4", "SRR5"],
                "avg_read_length": [150, 155, 148, 152, 500],  # SRR5 is outlier
                "gc_percent": [45.0, 47.0, 44.0, 46.0, 48.0],
                "total_reads": [10000, 12000, 11000, 10500, 9500],
            }
        )

        outliers = self.analyzer._detect_outliers(comparison_df)

        assert "SRR5" in outliers  # Should detect read length outlier

    def test_holm_adjusted_p_values_decide_significance(self):
        """G5: four tests at raw 0.01, 0.02, 0.03, 0.2 adjust to 0.04, 0.06, 0.06, 0.2; one is significant."""
        comparison_df = pd.DataFrame(
            {
                "accession": ["SRR1", "SRR2", "SRR3", "SRR4"],
                "group": ["A", "A", "B", "B"],
                "avg_read_length": [150, 155, 250, 245],
                "gc_percent": [45.0, 47.0, 52.0, 50.0],
                "total_reads": [10000, 12000, 8000, 9000],
                "complexity_score": [0.7, 0.75, 0.6, 0.65],
            }
        )
        raw = iter([0.01, 0.02, 0.03, 0.2])
        with patch("prepare.sra.analytics.compare_group_means", side_effect=lambda data: ("t-test", 1.0, next(raw))):
            tests = self.analyzer._perform_statistical_tests(
                comparison_df, {"A": ["SRR1", "SRR2"], "B": ["SRR3", "SRR4"]}
            )

        columns = ["avg_read_length", "gc_percent", "total_reads", "complexity_score"]
        assert [tests[c]["p_value"] for c in columns] == [0.01, 0.02, 0.03, 0.2]
        assert [round(tests[c]["p_adjusted"], 6) for c in columns] == [0.04, 0.06, 0.06, 0.2]
        assert [tests[c]["significant"] for c in columns] == [True, False, False, False]

    def test_holm_adjusted_helper(self):
        from prepare.sra.analytics import holm_adjusted

        assert [round(p, 6) for p in holm_adjusted([0.2, 0.01, 0.03, 0.02])] == [0.2, 0.04, 0.06, 0.06]
        assert holm_adjusted([0.6, 0.7]) == [1.0, 1.0]
        assert holm_adjusted([]) == []

    def test_recommendation_names_the_correction(self):
        tests = {"avg_read_length": {"significant": True, "p_value": 0.001, "p_adjusted": 0.004}}
        recommendations = self.analyzer._generate_comparative_recommendations({}, tests, [])
        assert any("after Holm correction" in rec for rec in recommendations)

    def test_generate_comparative_recommendations(self):
        """Test comparative analysis recommendations."""
        summary_stats = {"group_A": {"avg_read_length": {"mean": 150}}, "group_B": {"avg_read_length": {"mean": 250}}}

        statistical_tests = {"avg_read_length": {"significant": True, "p_value": 0.001}}

        outliers = ["SRR123456"]

        recommendations = self.analyzer._generate_comparative_recommendations(
            summary_stats, statistical_tests, outliers
        )

        assert len(recommendations) > 0
        assert any("significant" in rec.lower() for rec in recommendations)
        assert any("outlier" in rec.lower() for rec in recommendations)

    def test_prepare_visualization_data(self):
        """Test visualization data preparation."""
        comparison_df = pd.DataFrame(
            {
                "accession": ["SRR1", "SRR2", "SRR3"],
                "group": ["A", "A", "B"],
                "avg_read_length": [150, 155, 250],
                "gc_percent": [45.0, 47.0, 52.0],
            }
        )

        viz_data = self.analyzer._prepare_visualization_data(comparison_df)

        assert "boxplot_data" in viz_data
        assert "summary_table" in viz_data
        assert "correlation_matrix" in viz_data
        assert len(viz_data["boxplot_data"]) == 3


class TestDataStructures:
    """Test data structure creation and validation."""

    def test_quality_profile_creation(self):
        """Test QualityProfile creation."""
        profile = QualityProfile(
            accession="SRR123456",
            total_reads=10000,
            total_bases=1500000,
            avg_read_length=150.0,
            read_length_distribution={"101-150": 8000, "151-250": 2000},
            gc_percent=45.0,
            gc_histogram={},
            quality_distribution={"excellent_q30+": 0.8},
            n_content=0.01,
            contamination_indicators={"adapter_contamination": 0.02},
            complexity_score=0.7,
            duplication_rate=0.15,
            technology_confidence=0.9,
            quality_grade="good",
            recommendations=["High quality dataset"],
        )

        assert profile.accession == "SRR123456"
        assert profile.total_reads == 10000
        assert profile.quality_grade == "good"
        assert len(profile.recommendations) == 1

    def test_comparative_analysis_creation(self):
        """Test ComparativeAnalysis creation."""
        analysis = ComparativeAnalysis(
            dataset_groups={"A": ["SRR1", "SRR2"], "B": ["SRR3", "SRR4"]},
            summary_statistics={"A": {"mean_gc": 0.45}, "B": {"mean_gc": 0.52}},
            statistical_tests={"gc_content": {"p_value": 0.001, "significant": True}},
            outlier_datasets=["SRR5"],
            clustering_results=None,
            batch_effects={"batch1": 0.1},
            recommendations=["Significant differences detected"],
            visualization_data={},
        )

        assert len(analysis.dataset_groups) == 2
        assert "A" in analysis.dataset_groups
        assert analysis.statistical_tests["gc_content"]["significant"] is True
        assert "SRR5" in analysis.outlier_datasets

    def test_anomaly_report_creation(self):
        """Test AnomalyReport creation."""
        report = AnomalyReport(
            anomalous_datasets=["SRR123456", "SRR789012"],
            anomaly_types={"high_contamination": ["SRR123456"], "low_complexity": ["SRR789012"]},
            severity_scores={"SRR123456": 0.8, "SRR789012": 0.6},
            explanations={"SRR123456": "High adapter contamination", "SRR789012": "Low complexity"},
            recommended_actions={"SRR123456": ["Trim adapters"], "SRR789012": ["Check PCR artifacts"]},
        )

        assert len(report.anomalous_datasets) == 2
        assert "high_contamination" in report.anomaly_types
        assert report.severity_scores["SRR123456"] > report.severity_scores["SRR789012"]
        assert len(report.recommended_actions["SRR123456"]) > 0

    def test_processing_recommendations_creation(self):
        """Test ProcessingRecommendations creation."""
        recommendations = ProcessingRecommendations(
            accession="SRR123456",
            recommended_pipeline="paired_end",
            quality_trimming={"enabled": True, "quality_threshold": 20},
            adapter_removal={"enabled": True, "stringency": "high"},
            contamination_filtering={"enabled": True, "check_vector": True},
            assembly_parameters={"kmer_size": 31, "coverage_cutoff": 5},
            expected_coverage=50.0,
            computational_requirements={"memory_gb": 16, "cpu_cores": 4},
            estimated_processing_time="2.5 hours",
        )

        assert recommendations.accession == "SRR123456"
        assert recommendations.recommended_pipeline == "paired_end"
        assert recommendations.quality_trimming["enabled"] is True
        assert recommendations.computational_requirements["memory_gb"] == 16


def test_statistical_tests_are_json_serialisable():
    """json_safe converts numpy scalars, sets and non-finite floats to JSON-safe values."""
    import json

    import numpy as np

    from prepare.sra.analytics import json_safe

    payload = {"significant": np.bool_(True), "p": np.float64(0.01), "nan": float("nan"), "s": {1, 2}}
    text = json.dumps(json_safe(payload))
    assert json.loads(text) == {"significant": True, "p": 0.01, "nan": None, "s": [1, 2]}


def test_json_safe_converts_a_pandas_series_to_a_dict():
    """A pandas Series (e.g. a DataFrame column pulled out on its own) has no json.dumps
    support; json_safe must convert it via to_dict() rather than leaving it as-is."""
    import json

    import pandas as pd

    from prepare.sra.analytics import json_safe

    series = pd.Series([1, 2, 3], index=["a", "b", "c"])
    safe = json_safe(series)
    assert safe == {"a": 1, "b": 2, "c": 3}
    json.dumps(safe)  # must not raise


def test_json_safe_converts_a_pandas_dataframe_to_a_dict():
    """A pandas DataFrame (e.g. comparison.summary_statistics built with pd.DataFrame) is
    likewise converted via to_dict() rather than passed through unchanged."""
    import json

    import pandas as pd

    from prepare.sra.analytics import json_safe

    df = pd.DataFrame({"gc_content": [0.4, 0.5]})
    safe = json_safe(df)
    # json_safe's dict branch stringifies every key it recurses through, including the row
    # index DataFrame.to_dict() nests under each column.
    assert safe == {"gc_content": {"0": 0.4, "1": 0.5}}
    json.dumps(safe)  # must not raise


def _profile_with_gc(accession: str, gc_percent: float) -> QualityProfile:
    """Build a real QualityProfile with a chosen GC percent, otherwise matching the
    ``make_profile`` helper in tests/test_cli_sra_report.py."""
    return QualityProfile(
        accession=accession,
        total_reads=1000,
        reads_sampled=1000,
        total_bases=150000,
        avg_read_length=150.0,
        read_length_distribution={},
        gc_percent=gc_percent,
        gc_histogram={},
        quality_distribution={"excellent_q30+": 0.9},
        n_content=0.0,
        contamination_indicators={"adapter_contamination": 0.0},
        complexity_score=0.85,
        duplication_rate=None,
        technology_confidence=0.8,
        quality_grade="good",
        recommendations=[],
    )


@requires_analysis
def test_compare_datasets_result_dumps(tmp_path):
    """compare_datasets' real statistical_tests output (numpy bool included) must be
    JSON-serialisable once passed through json_safe, matching what sra_report does before
    writing sra_report.json."""
    import json

    from prepare.sra.analytics import SRADatasetAnalyzer, json_safe

    analyzer = SRADatasetAnalyzer(tmp_path)
    groups = {"a": ["A1", "A2", "A3"], "b": ["B1", "B2"]}
    profiles = {}
    for i, accession in enumerate(groups["a"]):
        profiles[accession] = _profile_with_gc(accession, 40.0 + i)
    for i, accession in enumerate(groups["b"]):
        profiles[accession] = _profile_with_gc(accession, 55.0 + i)

    result = analyzer.compare_datasets(groups, profiles=profiles)

    assert result.statistical_tests["gc_percent"]["test"] == "t-test"
    json.dumps(json_safe(result.statistical_tests))
