"""Tests for the staged FASTQ export of one sample (``prepare.data.extraction_staging``)."""

import gzip
from pathlib import Path
from unittest.mock import patch

import pytest

from prepare.data import registry as reg
from prepare.data.extraction_staging import publish_staged, remove_published, staged_sample_outputs
from prepare.data.file_io import visible_files
from prepare.data.read_extraction import _map_and_extract
from prepare.data.registry_reconcile import apply_reconcile, scan_reconcile
from helpers_extraction import _fake_tools

READ = "@r1\nACGT\n+\nIIII\n"


def _extract(tmp_path: Path, out_dir: Path, state, genome_id="GCF_1", paired=True):
    reads = [tmp_path / "SRR1_1.fastq.gz", tmp_path / "SRR1_2.fastq.gz"] if paired else [tmp_path / "SRR1.fastq.gz"]
    with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
        return _map_and_extract(
            accession="SRR1",
            reads=reads,
            reference=tmp_path / "ref.mmi",
            genome_fasta=tmp_path / "genome.fna",
            out_dir=out_dir,
            genome_id=genome_id,
            preset="sr",
            threads=1,
        )


def _interrupting_tools(state):
    """The fake tools, except that ``samtools fastq`` writes part of mate 1 and is interrupted."""
    fake = _fake_tools(state)

    def run(executable, args, **kwargs):
        if executable == "samtools" and args[0] == "fastq":
            out1 = Path(args[args.index("-1") + 1])
            with gzip.open(out1, "wt") as handle:
                handle.write("@r1\nAC")
            raise KeyboardInterrupt
        return fake(executable, args, **kwargs)

    return run


class TestCompleteSample:
    def test_files_are_published_under_their_final_names_with_the_exported_reads(self, tmp_path):
        out_dir = tmp_path / "targeted" / "SRR1"

        result = _extract(tmp_path, out_dir, {})

        assert result.files == [out_dir / "GCF_1_1.fastq.gz", out_dir / "GCF_1_2.fastq.gz"]
        for path in result.files:
            with gzip.open(path, "rt") as handle:
                assert handle.read() == READ
        # Nothing staged is left behind; only the record of what was published is hidden there.
        assert [p.name for p in out_dir.iterdir() if p.name.startswith(".")] == [".GCF_1.extracted.json"]

    def test_a_rerun_that_writes_fewer_files_removes_the_older_ones(self, tmp_path):
        out_dir = tmp_path / "targeted" / "SRR1"
        first = _extract(tmp_path, out_dir, {"nonempty": ("-1", "-2", "-s")})
        assert first.files == [out_dir / "GCF_1_1.fastq.gz", out_dir / "GCF_1_2.fastq.gz"]
        # The singletons are kept beside the pairs, as the in-place export kept them.
        assert (out_dir / "GCF_1_s.fastq.gz").is_file()

        result = _extract(tmp_path, out_dir, {"nonempty": ("-0",)})

        assert result.files == [out_dir / "GCF_1_0.fastq.gz"]
        names = sorted(p.name for p in visible_files(out_dir, "*.fastq.gz"))
        assert names == ["GCF_1_0.fastq.gz"]


class TestRerunCleanupStaysWithinTheGenome:
    """Genome ids can look like another genome's mate names (``G1`` paired writes ``G1_1.fastq.gz``)."""

    def test_a_single_end_rerun_of_g1_leaves_the_files_of_genome_g1_1(self, tmp_path):
        out_dir = tmp_path / "targeted" / "SRR1"
        other = _extract(tmp_path, out_dir, {"nonempty": ("-0",)}, genome_id="G1_1", paired=False)
        assert other.files == [out_dir / "G1_1.fastq.gz"]
        _extract(tmp_path, out_dir, {"nonempty": ("-0",)}, genome_id="G1", paired=False)

        result = _extract(tmp_path, out_dir, {"nonempty": ("-0",)}, genome_id="G1", paired=False)

        assert result.files == [out_dir / "G1.fastq.gz"]
        with gzip.open(out_dir / "G1_1.fastq.gz", "rt") as handle:
            assert handle.read() == READ
        names = sorted(p.name for p in visible_files(out_dir, "*.fastq.gz"))
        assert names == ["G1.fastq.gz", "G1_1.fastq.gz"]

    def test_a_paired_rerun_of_g1_replaces_its_own_single_end_file(self, tmp_path):
        out_dir = tmp_path / "targeted" / "SRR1"
        _extract(tmp_path, out_dir, {"nonempty": ("-0",)}, genome_id="G1", paired=False)
        assert (out_dir / "G1.fastq.gz").is_file()

        result = _extract(tmp_path, out_dir, {}, genome_id="G1")

        assert result.files == [out_dir / "G1_1.fastq.gz", out_dir / "G1_2.fastq.gz"]
        names = sorted(p.name for p in visible_files(out_dir, "*.fastq.gz"))
        assert names == ["G1_1.fastq.gz", "G1_2.fastq.gz"]

    def test_a_folder_without_records_loses_nothing_but_overwritten_names(self, tmp_path):
        out_dir = tmp_path / "targeted" / "SRR1"
        out_dir.mkdir(parents=True)
        (out_dir / "G1_s.fastq.gz").write_bytes(b"written before 0.6.0")

        _extract(tmp_path, out_dir, {}, genome_id="G1")

        assert (out_dir / "G1_s.fastq.gz").read_bytes() == b"written before 0.6.0"


class TestInterruptedSample:
    def test_an_interrupted_export_leaves_no_visible_fastq(self, tmp_path):
        out_dir = tmp_path / "targeted" / "SRR1"
        with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=_interrupting_tools({})):
            with pytest.raises(KeyboardInterrupt):
                _map_and_extract(
                    accession="SRR1",
                    reads=[tmp_path / "SRR1_1.fastq.gz", tmp_path / "SRR1_2.fastq.gz"],
                    reference=tmp_path / "ref.mmi",
                    genome_fasta=tmp_path / "genome.fna",
                    out_dir=out_dir,
                    genome_id="GCF_1",
                    preset="sr",
                    threads=1,
                )

        assert visible_files(out_dir, "*.fastq.gz") == []
        assert reg.scan_extractions(tmp_path / "targeted", ["GCF_1"]) == {}

    def test_a_staging_folder_left_by_a_kill_is_ignored_by_scan_and_reconcile(self, tmp_path):
        paths = reg.ProjectPaths(
            fastq=tmp_path / "fastq",
            metadata=tmp_path / "metadata",
            genomes=tmp_path / "genomes",
            targeted=tmp_path / "targeted",
            matches=tmp_path / "matches",
        )
        paths.genomes.mkdir()
        (paths.genomes / "GCF_1.fna").write_text(">c\nACGT\n")
        out_dir = paths.targeted / "SRR1"
        # What a SIGKILL during the export leaves: the staging folder and its partial files.
        leftover = out_dir / ".GCF_1.node7.4242.0badf00d.tmp"
        leftover.mkdir(parents=True)
        for name in ("GCF_1_1.fastq.gz", "GCF_1_2.fastq.gz"):
            with gzip.open(leftover / name, "wt") as handle:
                handle.write("@r1\nAC")

        registry = reg.bootstrap_from_disk(paths, target_path=tmp_path / reg.REGISTRY_FILENAME)
        assert reg.scan_extractions(paths.targeted, ["GCF_1"]) == {}
        report = apply_reconcile(registry, scan_reconcile(registry, paths))
        assert report.untracked_extractions == []
        assert "extractions" not in registry.datasets.get("SRR1", {})


class TestPublishStaged:
    def test_publishes_each_file_with_one_rename(self, tmp_path):
        out_dir = tmp_path / "out"
        with staged_sample_outputs(out_dir, "G") as staging:
            staged = [staging / "G_1.fastq.gz", staging / "G_2.fastq.gz"]
            for path in staged:
                path.write_bytes(b"x")
            published = publish_staged(staging, staged, out_dir, "G")
        assert published == [out_dir / "G_1.fastq.gz", out_dir / "G_2.fastq.gz"]
        assert sorted(p.name for p in visible_files(out_dir)) == ["G_1.fastq.gz", "G_2.fastq.gz"]
        assert sorted(p.name for p in out_dir.iterdir()) == [".G.extracted.json", "G_1.fastq.gz", "G_2.fastq.gz"]


class TestRemovePublished:
    def test_removes_the_recorded_files_and_the_record(self, tmp_path):
        out_dir = tmp_path / "out"
        with staged_sample_outputs(out_dir, "G") as staging:
            staged = [staging / "G_1.fastq.gz", staging / "G_2.fastq.gz"]
            for path in staged:
                path.write_bytes(b"x")
            publish_staged(staging, staged, out_dir, "G")
        (out_dir / "other.txt").write_text("keep")

        removed = remove_published(out_dir, "G")

        assert sorted(p.name for p in removed) == ["G_1.fastq.gz", "G_2.fastq.gz"]
        assert sorted(p.name for p in out_dir.iterdir()) == ["other.txt"]

    def test_is_a_no_op_without_a_record(self, tmp_path):
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        (out_dir / "G_1.fastq.gz").write_bytes(b"older")
        assert remove_published(out_dir, "G") == []
        assert (out_dir / "G_1.fastq.gz").exists()

    def test_keeps_a_file_another_genomes_record_lists(self, tmp_path):
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        (out_dir / "G1_1.fastq.gz").write_bytes(b"x")
        (out_dir / ".G1.extracted.json").write_text('["G1_1.fastq.gz"]\n')
        (out_dir / ".G1_1.extracted.json").write_text('["G1_1.fastq.gz"]\n')
        assert remove_published(out_dir, "G1") == []
        assert (out_dir / "G1_1.fastq.gz").exists()
        assert not (out_dir / ".G1.extracted.json").exists()


class TestZeroMappingForcedRerun:
    def _first_run(self, tmp_path, out_dir):
        _extract(tmp_path, out_dir, {})
        assert (out_dir / "GCF_1_1.fastq.gz").is_file()
        assert (out_dir / ".GCF_1.extracted.json").is_file()

    def test_nothing_mapped_in_the_alignment_removes_earlier_reads(self, tmp_path):
        out_dir = tmp_path / "targeted" / "SRR1"
        self._first_run(tmp_path, out_dir)
        result = _extract(tmp_path, out_dir, {"mapped_total": 0, "mapped": 0})
        assert result.files == []
        assert not (out_dir / "GCF_1_1.fastq.gz").exists()
        assert not (out_dir / "GCF_1_2.fastq.gz").exists()
        assert not (out_dir / ".GCF_1.extracted.json").exists()

    def test_nothing_kept_after_filtering_removes_earlier_reads(self, tmp_path):
        out_dir = tmp_path / "targeted" / "SRR1"
        self._first_run(tmp_path, out_dir)
        result = _extract(tmp_path, out_dir, {"mapped_total": 5, "mapped": 0})
        assert result.files == []
        assert not (out_dir / "GCF_1_1.fastq.gz").exists()
        assert not (out_dir / ".GCF_1.extracted.json").exists()
