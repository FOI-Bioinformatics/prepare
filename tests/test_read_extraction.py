"""Tests for targeted read extraction (prepare.data.read_extraction)."""

import gzip
import json
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from prepare.core.exceptions import DataAccessError, ProcessingError, SecurityError
from prepare.data.read_extraction import (
    ExtractionResult,
    _map_and_extract,
    _minimap2_map_args,
    _record_matches,
    _run_minimap2,
    _sample_reads,
    assemble_extracted_reads,
    assembly_coverage,
    build_index,
    extract_target_reads,
    fasta_length,
    megahit_version,
    resolve_assembly_threads,
    resolve_index_path,
    reference_coverage,
    select_samples_for_genome,
    summarise_contigs,
    summarise_coverage_table,
)
from prepare.data.registry import load_registry, record_extraction, resolve_project_path, save_registry
from helpers_extraction import _fake_tools, coverage_table


def _make_tree(tmp, paired=True):
    """Build a fixture: containment table, per-accession FASTQ, and a genome FASTA."""
    root = Path(tmp)
    table = root / "parsed_containment.txt"
    table.write_text(
        "\tGCF_1\tGCF_2\tmax_containment\n"
        "SRR1\t0.90\t0.0\t0.90\n"
        "SRR2\t0.40\t0.85\t0.85\n"
        "SRR3\t0.05\t0.0\t0.05\n"
    )
    for acc in ("SRR1", "SRR2", "SRR3"):
        d = root / "fastq" / acc
        d.mkdir(parents=True)
        (d / f"{acc}_1.fastq.gz").write_text("x")
        if paired:
            (d / f"{acc}_2.fastq.gz").write_text("x")
    genome = root / "GCF_1.fna"
    genome.write_text(">s\nACGT\n")
    return root, table, genome


class TestSelectSamples:
    def test_threshold_selection(self):
        df = pd.DataFrame({"GCF_1": [0.9, 0.4, 0.05]}, index=["SRR1", "SRR2", "SRR3"])
        assert select_samples_for_genome(df, "GCF_1", 0.5) == ["SRR1"]
        assert select_samples_for_genome(df, "GCF_1", 0.3) == ["SRR1", "SRR2"]

    def test_missing_genome_column(self):
        df = pd.DataFrame({"GCF_1": [0.9]}, index=["SRR1"])
        with pytest.raises(ProcessingError):
            select_samples_for_genome(df, "GCF_missing", 0.1)

    def test_selected_samples_applies_exclusions_and_the_selection(self, tmp_path, caplog):
        from prepare.data.read_extraction import selected_samples
        from prepare.data.registry import record_exclusion, record_selection

        table = tmp_path / "parsed_containment.txt"
        table.write_text("\tGCF_1\nSRR1\t0.9\nSRR2\t0.6\nSRR3\t0.5\nSRR4\t0.3\n")
        registry = load_registry(tmp_path / "registry.json")
        record_exclusion(registry, "SRR2", "blacklisted")
        record_selection(registry, ["SRR2", "SRR3"], {"genome_id": "GCF_1"}, "accessions.txt")

        assert selected_samples(table, "GCF_1", 0.1) == ["SRR1", "SRR2", "SRR3", "SRR4"]
        with caplog.at_level(logging.INFO):
            assert selected_samples(table, "GCF_1", 0.1, registry=registry) == ["SRR3"]
        assert "4 sample(s) in the table" in caplog.text and "1 excluded" in caplog.text
        assert selected_samples(table, "GCF_1", 0.1, registry=registry, ignore_selection=True) == [
            "SRR1",
            "SRR3",
            "SRR4",
        ]


class TestSampleReads:
    """_sample_reads: which files of an accession folder are handed to minimap2."""

    def _acc_dir(self, tmp_path, *names):
        acc_dir = tmp_path / "fastq" / "SRR1"
        acc_dir.mkdir(parents=True, exist_ok=True)
        for name in names:
            (acc_dir / name).write_text("@r\nACGT\n+\nIIII\n")
        return tmp_path / "fastq"

    def test_paired_folder_returns_mates_in_order(self, tmp_path):
        folder = self._acc_dir(tmp_path, "SRR1_2.fastq.gz", "SRR1_1.fastq.gz")

        assert [p.name for p in _sample_reads(folder, "SRR1")] == ["SRR1_1.fastq.gz", "SRR1_2.fastq.gz"]

    def test_orphan_file_is_not_passed_to_the_aligner(self, tmp_path):
        """--split-3 writes unpaired spots to a bare <acc>.fastq; minimap2 gets the mates only."""
        folder = self._acc_dir(tmp_path, "SRR1.fastq", "SRR1_1.fastq", "SRR1_2.fastq")

        assert [p.name for p in _sample_reads(folder, "SRR1")] == ["SRR1_1.fastq", "SRR1_2.fastq"]

    def test_single_end_folder_returns_the_bare_file(self, tmp_path):
        folder = self._acc_dir(tmp_path, "SRR1.fastq")

        assert [p.name for p in _sample_reads(folder, "SRR1")] == ["SRR1.fastq"]

    def test_partial_download_leftovers_are_ignored(self, tmp_path):
        """Zero-byte files and .gz.tmp.<pid> leftovers are not reads."""
        folder = self._acc_dir(tmp_path, "SRR1_1.fastq")
        (folder / "SRR1" / "SRR1_2.fastq").write_text("")
        (folder / "SRR1" / "SRR1_1.fastq.gz.tmp.4242").write_text("junk")

        assert [p.name for p in _sample_reads(folder, "SRR1")] == ["SRR1_1.fastq"]

    def test_missing_folder_returns_empty(self, tmp_path):
        assert _sample_reads(tmp_path / "fastq", "SRR1") == []


class TestExtractTargetReads:
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_paired_extraction_command_construction(self, mock_run):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            results = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
            )
            assert list(results) == ["SRR1"]
            assert [p.name for p in results["SRR1"].files] == ["GCF_1_1.fastq.gz", "GCF_1_2.fastq.gz"]
            # The orphan file (-0) is removed when empty.
            assert not (root / "targeted" / "SRR1" / "GCF_1_0.fastq.gz").exists()

        tools = [c[0] for c in state["calls"]]
        # First minimap2 call builds the shared index (-d); the second aligns the sample.
        # The last two samtools calls sort the filtered BAM and write its reference coverage.
        assert tools == ["minimap2", "minimap2"] + ["samtools"] * 6
        assert "-d" in state["calls"][0][1]
        assert state["calls"][1][1][:2] == ["-a", "-x"]
        assert state["calls"][2][1][:2] == ["view", "-c"]
        assert state["calls"][3][1][:4] == ["view", "-b", "-F", "0x904"]
        assert state["calls"][4][1][:2] == ["view", "-c"]
        fastq_args = state["calls"][5][1]
        assert fastq_args[0] == "fastq"
        assert all(flag in fastq_args for flag in ("-1", "-2", "-s", "-0", "-@"))
        assert state["calls"][6][1][0] == "sort" and state["calls"][7][1][0] == "coverage"

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_single_end_uses_flag_0(self, mock_run):
        state = {"nonempty": ("-0",)}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=False)
            results = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
            )
            assert [p.name for p in results["SRR1"].files] == ["GCF_1.fastq.gz"]
        fastq_args = state["calls"][5][1]
        assert "-0" in fastq_args
        assert "-@" in fastq_args

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_zero_mapped_records_writes_nothing(self, mock_run, caplog):
        state = {"mapped": 0}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            with caplog.at_level("WARNING"):
                results = extract_target_reads(
                    parsed_containment=table,
                    genome_id="GCF_1",
                    genome_fasta=genome,
                    fastq_folder=root / "fastq",
                    output_folder=root / "targeted",
                    threshold=0.5,
                )
            assert results == {"SRR1": ExtractionResult([], 0)}
            assert list((root / "targeted" / "SRR1").glob("*.fastq.gz")) == []
        assert "No reads from SRR1 mapped to GCF_1" in caplog.text
        # Nothing mapped (per the SAM count), so filtering/counting the BAM and running
        # samtools fastq are all skipped.
        assert [c[0] for c in state["calls"]] == ["minimap2", "minimap2", "samtools"]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_unequal_mates_fall_back_to_orphan_file(self, mock_run, caplog):
        state = {"unequal": True, "nonempty": ("-0",)}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            with caplog.at_level("WARNING"):
                results = extract_target_reads(
                    parsed_containment=table,
                    genome_id="GCF_1",
                    genome_fasta=genome,
                    fastq_folder=root / "fastq",
                    output_folder=root / "targeted",
                    threshold=0.5,
                )
            assert [p.name for p in results["SRR1"].files] == ["GCF_1_0.fastq.gz"]
        assert "different read counts" in caplog.text

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_paired_export_keeps_the_singleton_file(self, mock_run):
        """P10: reads whose mate did not map (samtools fastq -s) are part of the extraction."""
        state = {"nonempty": ("-1", "-2", "-s")}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            results = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
            )
            assert [p.name for p in results["SRR1"].files] == [
                "GCF_1_1.fastq.gz",
                "GCF_1_2.fastq.gz",
                "GCF_1_s.fastq.gz",
            ]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_mapped_count_is_logged(self, mock_run, caplog):
        state = {"mapped": 1234}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            with caplog.at_level("INFO"):
                extract_target_reads(
                    parsed_containment=table,
                    genome_id="GCF_1",
                    genome_fasta=genome,
                    fastq_folder=root / "fastq",
                    output_folder=root / "targeted",
                    threshold=0.5,
                )
        assert "1234 mapped records for SRR1" in caplog.text

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_dry_run_runs_no_tools(self, mock_run):
        seen = []
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp)
            results = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.3,
                dry_run=True,
                on_result=lambda acc, result: seen.append(acc),
            )
        assert set(results) == {"SRR1", "SRR2"}
        mock_run.assert_not_called()
        assert seen == []  # a dry run records nothing

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_on_result_runs_per_sample_and_failures_are_isolated(self, mock_run, caplog):
        """A callback that raises for one sample must not stop the next sample's extraction."""
        mock_run.side_effect = _fake_tools({})
        seen = []

        def on_result(accession, result):
            seen.append((accession, result.mapped_records))
            if accession == "SRR1":
                raise DataAccessError("registry file could not be written")

        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            with caplog.at_level("WARNING"):
                results = extract_target_reads(
                    parsed_containment=table,
                    genome_id="GCF_1",
                    genome_fasta=genome,
                    fastq_folder=root / "fastq",
                    output_folder=root / "targeted",
                    threshold=0.3,
                    on_result=on_result,
                )
            assert [p.name for p in results["SRR2"].files] == ["GCF_1_1.fastq.gz", "GCF_1_2.fastq.gz"]
        assert [acc for acc, _ in seen] == ["SRR1", "SRR2"]
        assert set(results) == {"SRR1", "SRR2"}
        assert "Recording the extraction result for SRR1 failed" in caplog.text

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_a_lock_timeout_while_recording_is_a_transient_error(self, mock_run):
        """P12: the run finishes the remaining samples, then the lock timeout surfaces (exit 4)."""
        from prepare.core.exceptions import LockTimeoutError, TransientError

        mock_run.side_effect = _fake_tools({})
        seen = []

        def on_result(accession, result):
            seen.append(accession)
            if accession == "SRR1":
                raise LockTimeoutError("registry lock held by another process")

        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            with pytest.raises(TransientError, match="registry lock"):
                extract_target_reads(
                    parsed_containment=table,
                    genome_id="GCF_1",
                    genome_fasta=genome,
                    fastq_folder=root / "fastq",
                    output_folder=root / "targeted",
                    threshold=0.3,
                    on_result=on_result,
                )
            assert (root / "targeted" / "SRR2" / "GCF_1_1.fastq.gz").exists()  # the run went on
        assert seen == ["SRR1", "SRR2"]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_an_unexpected_callback_error_propagates(self, mock_run):
        """A bug in the callback is not hidden behind a warning."""
        mock_run.side_effect = _fake_tools({})

        def on_result(accession, result):
            raise RuntimeError("bug")

        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            with pytest.raises(RuntimeError):
                extract_target_reads(
                    parsed_containment=table,
                    genome_id="GCF_1",
                    genome_fasta=genome,
                    fastq_folder=root / "fastq",
                    output_folder=root / "targeted",
                    threshold=0.3,
                    on_result=on_result,
                )

    def test_missing_table_raises(self):
        with pytest.raises(DataAccessError):
            extract_target_reads("/no/such/table.txt", "GCF_1", "/no/genome.fna")

    def test_bad_preset_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp)
            with pytest.raises(ProcessingError):
                extract_target_reads(table, "GCF_1", genome, preset="not-a-preset")

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_missing_genome_fasta_raises(self, mock_run):
        with tempfile.TemporaryDirectory() as tmp:
            root, table, _ = _make_tree(tmp)
            with pytest.raises(DataAccessError):
                extract_target_reads(table, "GCF_1", root / "absent.fna", fastq_folder=root / "fastq")


class TestBuildIndex:
    """build_index/resolve_index_path: one shared minimap2 index per genome/preset."""

    def test_index_built_once_and_reused(self, tmp_path):
        genome = tmp_path / "g.fna"
        genome.write_text(">s\nACGT\n")
        index_dir = tmp_path / ".index"
        state = {}
        with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
            first = build_index(genome, "sr", index_dir)
            assert len(state["calls"]) == 1
            executable, args = state["calls"][0]
            assert executable == "minimap2"
            # Built under a temporary name and moved into place, so a concurrent run never
            # reads a half-written index.
            assert args[:3] == ["-x", "sr", "-d"]
            assert Path(args[3]).parent == first.parent and Path(args[3]).name.startswith(f".{first.name}.")
            assert args[4] == str(genome)
            assert first.is_file()
            again = build_index(genome, "sr", index_dir)
        assert again == first
        assert len(state["calls"]) == 1  # the FASTA has not changed, so the index is reused
        assert not list(index_dir.glob("*.tmp"))

    def test_an_empty_index_is_never_published(self, tmp_path):
        """A minimap2 that exits 0 without writing (or after a full disk truncated it) fails the build."""
        genome = tmp_path / "g.fna"
        genome.write_text(">s\nACGT\n")
        index_dir = tmp_path / ".index"
        state = {"empty_index": True}
        with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
            with pytest.raises(ProcessingError, match="empty"):
                build_index(genome, "sr", index_dir)
        assert [p.name for p in index_dir.iterdir() if not p.name.endswith(".lock")] == []

    def test_an_empty_index_left_by_an_older_version_is_rebuilt(self, tmp_path):
        genome = tmp_path / "g.fna"
        genome.write_text(">s\nACGT\n")
        index_dir = tmp_path / ".index"
        state = {}
        with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
            index = build_index(genome, "sr", index_dir)
            index.write_bytes(b"")
            again = build_index(genome, "sr", index_dir)
        assert again == index
        assert len(state["calls"]) == 2
        assert index.stat().st_size > 0

    def test_index_rebuilt_for_a_different_genome_with_the_same_name(self, tmp_path):
        """Two genome files can share a stem, so the index is keyed on the FASTA's identity.

        The second FASTA deliberately carries an older mtime, which is what a copy restored
        from a tarball or moved with ``cp -p`` looks like; an mtime comparison alone would
        map it against the first genome's index.
        """
        old = tmp_path / "genomes_v1" / "wMel.fna"
        new = tmp_path / "genomes_v2" / "wMel.fna"
        old.parent.mkdir(parents=True)
        new.parent.mkdir(parents=True)
        old.write_text(">old\nAAAACCCC\n")
        new.write_text(">new\nCCCC\n")
        past = old.stat().st_mtime - 3600
        os.utime(new, (past, past))

        index_dir = tmp_path / ".index"
        state = {}
        with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
            first = build_index(old, "sr", index_dir)
            second = build_index(new, "sr", index_dir)

        assert first == second  # same stem, so the same index path
        assert len([c for c in state["calls"] if "-d" in c[1]]) == 2
        record = json.loads(index_dir.joinpath("wMel.sr.mmi.json").read_text())
        assert record["fasta"] == str(new.resolve())

    def test_a_failed_build_leaves_no_index(self, tmp_path):
        genome = tmp_path / "g.fna"
        genome.write_text(">s\nACGT\n")
        index_dir = tmp_path / ".index"

        def failing_build(executable, args, **kwargs):
            Path(args[args.index("-d") + 1]).write_bytes(b"half")  # a partial index
            raise SecurityError("minimap2 refused the command")

        with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=failing_build):
            with pytest.raises(SecurityError):
                build_index(genome, "sr", index_dir)

        assert not (index_dir / "g.sr.mmi").exists()
        assert not list(index_dir.glob("*.tmp"))

    def test_index_rebuilt_when_fasta_touched_newer(self, tmp_path):
        genome = tmp_path / "g.fna"
        genome.write_text(">s\nACGT\n")
        index_dir = tmp_path / ".index"
        state = {}
        with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
            index_path = build_index(genome, "sr", index_dir)
            assert len(state["calls"]) == 1

            newer = index_path.stat().st_mtime + 10
            os.utime(genome, (newer, newer))

            build_index(genome, "sr", index_dir)
        assert len(state["calls"]) == 2

    def test_resolve_index_path_names_the_file_by_genome_stem_and_preset(self, tmp_path):
        path = resolve_index_path(tmp_path / "genomes" / "GCF_1.fna", "map-ont", tmp_path / ".index")
        assert path == tmp_path / ".index" / "GCF_1.map-ont.mmi"


class TestIndexReuseAcrossSamples:
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_index_built_once_for_three_samples_and_reused(self, mock_run):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            results = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.05,
            )
            assert set(results) == {"SRR1", "SRR2", "SRR3"}
        index_calls = [c for c in state["calls"] if c[0] == "minimap2" and "-d" in c[1]]
        assert len(index_calls) == 1


class TestMapqAndSamLifetime:
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_min_mapq_passed_to_the_filter_and_logged(self, mock_run, caplog):
        state = {"mapped_total": 100, "mapped": 80}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            with caplog.at_level("INFO"):
                extract_target_reads(
                    parsed_containment=table,
                    genome_id="GCF_1",
                    genome_fasta=genome,
                    fastq_folder=root / "fastq",
                    output_folder=root / "targeted",
                    threshold=0.5,
                    min_mapq=20,
                )
        filter_call = next(c for c in state["calls"] if c[0] == "samtools" and c[1][:2] == ["view", "-b"])
        args = filter_call[1]
        assert args[args.index("-F") + 1] == "0x904"
        assert args[args.index("-q") + 1] == "20"
        assert "kept 80 of 100 mapped records (secondary/supplementary and MAPQ below 20 removed: 20)" in caplog.text

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_no_min_mapq_omits_the_q_flag(self, mock_run):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
            )
        filter_call = next(c for c in state["calls"] if c[0] == "samtools" and c[1][:2] == ["view", "-b"])
        assert "-q" not in filter_call[1]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_sam_is_removed_after_the_run(self, mock_run):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
            )
            assert list((root / "targeted" / "SRR1").glob("*.sam")) == []

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_sam_is_kept_with_the_debug_flag(self, mock_run):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
                keep_sam=True,
            )
            assert len(list((root / "targeted" / "SRR1").glob("*.sam"))) == 1

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_sam_written_under_temp_folder_when_given(self, mock_run):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            temp_folder = root / "scratch"
            extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
                temp_folder=temp_folder,
                keep_sam=True,
            )
            assert list(temp_folder.glob("*.sam"))
            assert list((root / "targeted" / "SRR1").glob("*.sam")) == []


class TestSamHitOnly:
    """``--sam-hit-only`` regression: minimap2 must emit only aligned records, and every
    count derived from the SAM/BAM afterwards must stay exactly what it was before the
    flag was added. Two unmapped reads and three mapped (one of them secondary) pin
    mapped_total (the ``-F 4`` SAM count, computed on the pre-change code) at 3 and the
    kept count after the ``0x904`` filter at 2; since ``--sam-hit-only`` only removes the
    unmapped records minimap2 itself would otherwise have written, both filters see the
    same records either way and these numbers must not move.
    """

    RECORDS = [
        ("read1", 4),  # unmapped
        ("read2", 4),  # unmapped
        ("read3", 0),  # mapped, primary
        ("read4", 0),  # mapped, primary
        ("read5", 256),  # mapped, secondary
    ]

    def test_builder_emits_the_flag(self):
        args = _minimap2_map_args("sr", 4, Path("out.sam"), Path("ref.mmi"), [Path("r.fastq.gz")])
        assert "--sam-hit-only" in args

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_mapped_total_and_kept_count_are_pinned(self, mock_run, tmp_path):
        state = {"sam_records": self.RECORDS}
        mock_run.side_effect = _fake_tools(state)

        result = _map_and_extract(
            accession="SRR1",
            reads=[tmp_path / "SRR1_1.fastq.gz", tmp_path / "SRR1_2.fastq.gz"],
            reference=tmp_path / "ref.mmi",
            genome_fasta=tmp_path / "genome.fna",
            out_dir=tmp_path / "out",
            genome_id="GCF_1",
            preset="sr",
            threads=1,
        )

        assert result.mapped_total == 3
        assert result.mapped_records == 2

        align_call = next(c for c in state["calls"] if c[0] == "minimap2")
        assert "--sam-hit-only" in align_call[1]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_assembly_coverage_mapping_rate_is_pinned(self, mock_run, tmp_path):
        state = {"sam_records": self.RECORDS}
        mock_run.side_effect = _fake_tools(state)
        contigs = tmp_path / "final.contigs.fa"
        contigs.write_text(">c1 len=10\nACGT\n")
        reads_file = tmp_path / "r.fastq.gz"
        with gzip.open(reads_file, "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n")

        result = assembly_coverage(contigs, [reads_file], "sr", 1, tmp_path, mapped_reads=2)

        assert result["reads_mapped"] == 2
        assert result["mapping_rate"] == 1.0


class TestTruncatedDownloadSkip:
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_truncated_sample_is_skipped_by_default(self, mock_run, caplog):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            truncated = {"SRR1": {"verdict": "truncated", "reads_r1": 5, "expected_spots": 20}}
            with caplog.at_level("WARNING"):
                results = extract_target_reads(
                    parsed_containment=table,
                    genome_id="GCF_1",
                    genome_fasta=genome,
                    fastq_folder=root / "fastq",
                    output_folder=root / "targeted",
                    threshold=0.5,
                    truncated_downloads=truncated,
                )
        assert "SRR1" not in results
        assert not mock_run.called
        assert "skipped SRR1: download truncated (5 of 20 spots); use --allow-truncated" in caplog.text

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_allow_truncated_extracts_anyway(self, mock_run):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            truncated = {"SRR1": {"verdict": "truncated", "reads_r1": 5, "expected_spots": 20}}
            results = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
                truncated_downloads=truncated,
                allow_truncated=True,
            )
        assert "SRR1" in results
        assert mock_run.called


class TestMateCountMismatch:
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_mismatched_mate_counts_map_single_end_without_the_stderr_marker(self, mock_run, caplog):
        """mate_counts alone (no minimap2 stderr marker) is enough to trigger the fallback."""
        state = {"nonempty": ("-0",)}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            with caplog.at_level("WARNING"):
                results = extract_target_reads(
                    parsed_containment=table,
                    genome_id="GCF_1",
                    genome_fasta=genome,
                    fastq_folder=root / "fastq",
                    output_folder=root / "targeted",
                    threshold=0.5,
                    mate_counts={"SRR1": (100, 80)},
                )
        assert results["SRR1"].unequal_mates is True
        assert [p.name for p in results["SRR1"].files] == ["GCF_1_0.fastq.gz"]
        assert "mapped independently as single-end" in caplog.text

        minimap2_calls = [c for c in state["calls"] if c[0] == "minimap2"]
        # one shared index build, plus one alignment per mate file
        assert len(minimap2_calls) == 3
        alignment_calls = [c for c in minimap2_calls if "-a" in c[1]]
        assert all(len([r for r in c[1] if str(r).endswith(".fastq.gz")]) == 1 for c in alignment_calls)

        cat_calls = [c for c in state["calls"] if c[0] == "samtools" and c[1][:1] == ["cat"]]
        assert len(cat_calls) == 1

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_matching_mate_counts_map_as_a_pair(self, mock_run):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            results = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
                mate_counts={"SRR1": (100, 100)},
            )
        assert results["SRR1"].unequal_mates is False
        alignment_calls = [c for c in state["calls"] if c[0] == "minimap2" and "-a" in c[1]]
        assert len(alignment_calls) == 1


class TestRecordMatches:
    def test_a_different_threshold_does_not_remap(self, tmp_path):
        """P3: the threshold decides which samples are tried, not what one sample's extraction is."""
        record = {
            "genome_fasta": str(tmp_path / "g.fna"),
            "preset": "sr",
            "threshold": 0.5,
            "min_mapq": 0,
            "mapped_reads": 0,
        }
        (tmp_path / "g.fna").write_text(">a\nA\n")
        assert _record_matches(record, tmp_path / "g.fna", "sr", 0.1, min_mapq=0)

    def test_a_different_mapq_still_remaps(self, tmp_path):
        record = {"genome_fasta": str(tmp_path / "g.fna"), "preset": "sr", "min_mapq": 20, "mapped_reads": 0}
        (tmp_path / "g.fna").write_text(">a\nA\n")
        assert not _record_matches(record, tmp_path / "g.fna", "sr", 0.1, min_mapq=0)
        assert not _record_matches(record, tmp_path / "g.fna", "map-ont", 0.1, min_mapq=20)

    def test_record_matches_compares_min_mapq(self, tmp_path):
        """A record's min_mapq is compared like every other extraction parameter: a
        mismatch means the sample must be redone rather than skipped."""
        record = {
            "genome_fasta": str(tmp_path / "g.fna"),
            "preset": "sr",
            "threshold": 0.1,
            "min_mapq": 20,
            "mapped_reads": 0,
        }
        (tmp_path / "g.fna").write_text(">a\nA\n")
        assert _record_matches(record, tmp_path / "g.fna", "sr", 0.1, min_mapq=20)
        assert not _record_matches(record, tmp_path / "g.fna", "sr", 0.1, min_mapq=0)


class TestExtractionIdempotency:
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_recorded_extraction_is_skipped_unless_forced(self, mock_run):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            first = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
            )
            record = {
                "SRR1": {
                    "genome_fasta": str(genome),
                    "preset": "sr",
                    "threshold": 0.5,
                    "mapped_reads": first["SRR1"].mapped_records,
                    "unequal_mates": False,
                    "files": [str(p) for p in first["SRR1"].files],
                }
            }
            calls_before = len(state["calls"])
            notified = []
            again = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
                already_done=record,
                on_result=lambda acc, result: notified.append((acc, result.skipped)),
            )
            assert again["SRR1"].skipped is True and len(state["calls"]) == calls_before
            assert notified == [("SRR1", True)]  # skips are reported too
            forced = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
                already_done=record,
                force=True,
            )
            assert forced["SRR1"].skipped is False and len(state["calls"]) > calls_before

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_assembly_skips_existing_contigs_and_redoes_empty_dir(self, mock_run, caplog):
        mock_run.side_effect = _fake_tools({})
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "asm"
            out.mkdir()
            (out / "final.contigs.fa").write_text(">c len=4\nACGT\n")
            _, ran = assemble_extracted_reads([Path(tmp) / "r1.fq.gz", Path(tmp) / "r2.fq.gz"], out)
            assert not mock_run.called and ran is False
            (out / "final.contigs.fa").unlink()
            with caplog.at_level("WARNING"):
                _, ran = assemble_extracted_reads([Path(tmp) / "r1.fq.gz"], out)
            assert "holds no final.contigs.fa (interrupted run?); assembling again" in caplog.text
            assert mock_run.called and ran is True
            assert (out / "final.contigs.fa").exists()
            assert [p.name for p in Path(tmp).iterdir()] == ["asm"]  # no staging folder left behind

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_inferred_record_without_parameters_is_skipped(self, mock_run):
        """A record bootstrapped from disk has no parameters; each missing one is a wildcard."""
        state = {}
        mock_run.side_effect = _fake_tools(state)
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            first = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
            )
            inferred = {
                "SRR1": {
                    "mapped_reads": 5,
                    "files": [str(p) for p in first["SRR1"].files],
                    "unequal_mates": False,
                }
            }
            calls_before = len(state["calls"])
            again = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
                already_done=inferred,
            )
            assert again["SRR1"].skipped is True and len(state["calls"]) == calls_before

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_inferred_zero_mapped_record_is_skipped(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            results = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
                already_done={"SRR1": {"mapped_reads": 0, "files": [], "unequal_mates": False}},
            )
        assert results["SRR1"].skipped is True
        mock_run.assert_not_called()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_redo_reason_names_missing_output(self, mock_run, caplog):
        """The recorded parameters still match, but a recorded output file was removed from
        disk; the redo log must name that reason, not blame the parameters."""
        mock_run.side_effect = _fake_tools({})
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            record = {
                "SRR1": {
                    "genome_fasta": str(genome),
                    "preset": "sr",
                    "threshold": 0.5,
                    "mapped_reads": 5,
                    "unequal_mates": False,
                    "files": [str(root / "targeted" / "SRR1" / "GCF_1_1.fastq.gz")],
                }
            }
            with caplog.at_level("INFO"):
                extract_target_reads(
                    parsed_containment=table,
                    genome_id="GCF_1",
                    genome_fasta=genome,
                    fastq_folder=root / "fastq",
                    output_folder=root / "targeted",
                    threshold=0.5,
                    already_done=record,
                )
        assert "recorded output file missing" in caplog.text
        assert "recorded parameters differ" not in caplog.text

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_relative_genome_path_matches_recorded_path(self, mock_run, monkeypatch):
        """./genomes/x.fna and genomes/x.fna are the same file, so the record still matches."""
        mock_run.side_effect = _fake_tools({})
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            monkeypatch.chdir(root)
            record = {
                "SRR1": {
                    "genome_fasta": "./GCF_1.fna",
                    "preset": "sr",
                    "threshold": 0.5,
                    "mapped_reads": 0,
                    "files": [],
                    "unequal_mates": False,
                }
            }
            results = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta="GCF_1.fna",
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
                already_done=record,
            )
        assert results["SRR1"].skipped is True
        mock_run.assert_not_called()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_resolved_record_still_matches_after_the_project_directory_is_renamed(self, mock_run, tmp_path):
        """A registry record is written under the project's original location, the project
        directory is then renamed, and the registry (now loaded from the new location) is
        resolved through resolve_project_path exactly as the CLI does; the resolved record
        must still match so extract_target_reads skips instead of remapping."""
        mock_run.side_effect = _fake_tools({})
        old_root = tmp_path / "a"
        old_root.mkdir()
        _make_tree(old_root, paired=True)
        genome = old_root / "GCF_1.fna"
        extracted_dir = old_root / "targeted" / "SRR1"
        extracted_dir.mkdir(parents=True)
        extracted_files = [extracted_dir / "GCF_1_1.fastq.gz", extracted_dir / "GCF_1_2.fastq.gz"]
        for f in extracted_files:
            f.write_text("@r\nACGT\n+\nIIII\n")

        registry_file = old_root / "prepare_registry.json"
        registry = load_registry(registry_file)
        record_extraction(
            registry,
            "SRR1",
            "GCF_1",
            extracted_files,
            42,
            False,
            {"genome_fasta": genome, "preset": "sr", "threshold": 0.5},
        )
        save_registry(registry)

        new_root = tmp_path / "b"
        old_root.rename(new_root)

        reloaded = load_registry(new_root / "prepare_registry.json")
        stored = reloaded.datasets["SRR1"]["extractions"]["GCF_1"]
        resolved_record = {
            "SRR1": {
                **stored,
                "genome_fasta": str(resolve_project_path(reloaded, stored["genome_fasta"])),
                "files": [str(resolve_project_path(reloaded, p)) for p in stored["files"]],
            }
        }

        results = extract_target_reads(
            parsed_containment=new_root / "parsed_containment.txt",
            genome_id="GCF_1",
            genome_fasta=new_root / "GCF_1.fna",
            fastq_folder=new_root / "fastq",
            output_folder=new_root / "targeted",
            threshold=0.5,
            already_done=resolved_record,
        )
        assert results["SRR1"].skipped is True
        mock_run.assert_not_called()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_dry_run_reports_would_be_skips(self, mock_run, caplog):
        mock_run.side_effect = _fake_tools({})
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            record = {
                "SRR1": {
                    "genome_fasta": str(genome),
                    "preset": "sr",
                    "threshold": 0.5,
                    "mapped_reads": 42,
                    "files": [],
                    "unequal_mates": False,
                }
            }
            with caplog.at_level("INFO"):
                results = extract_target_reads(
                    parsed_containment=table,
                    genome_id="GCF_1",
                    genome_fasta=genome,
                    fastq_folder=root / "fastq",
                    output_folder=root / "targeted",
                    threshold=0.5,
                    dry_run=True,
                    already_done=record,
                )
        assert results["SRR1"].skipped is True and results["SRR1"].mapped_records == 42
        assert "would skip SRR1 (already extracted, 42 mapped reads); use --force to redo" in caplog.text
        mock_run.assert_not_called()

    def test_dry_run_summarises_missing_fastq(self, tmp_path, caplog):
        """With ``available`` given, samples without FASTQ are reported once as a count,
        not with a warning per sample (an audit found one WARNING per non-downloaded
        sample on a 16603-sample dry run)."""
        root, table, genome = _make_tree(tmp_path)
        (root / "fastq" / "SRR2").rename(root / "fastq" / "SRR2_gone")
        with caplog.at_level(logging.INFO):
            extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.0,
                dry_run=True,
                available={"SRR1", "SRR3"},
            )
        assert caplog.text.count("No FASTQ files found for") == 0
        assert "1 of the 3 selected sample(s) have no FASTQ" in caplog.text

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_available_filter_keeps_already_done_samples_without_fastq(self, mock_run, tmp_path, caplog):
        """A sample whose extraction is already recorded must not be dropped as 'no FASTQ'
        just because its input FASTQ folder is gone (e.g. after store_gc/store_unlink runs
        once a project's extraction against a genome is complete); it is still reported as
        already extracted, not as missing."""
        mock_run.side_effect = _fake_tools({})
        root, table, genome = _make_tree(tmp_path, paired=True)
        extracted_dir = root / "targeted" / "SRR1"
        extracted_dir.mkdir(parents=True)
        extracted_files = [extracted_dir / "GCF_1_1.fastq.gz", extracted_dir / "GCF_1_2.fastq.gz"]
        for f in extracted_files:
            f.write_text("@r\nACGT\n+\nIIII\n")
        record = {
            "SRR1": {
                "genome_fasta": str(genome),
                "preset": "sr",
                "threshold": 0.5,
                "mapped_reads": 42,
                "unequal_mates": False,
                "files": [str(p) for p in extracted_files],
            }
        }
        # The input FASTQ folder is gone, so a caller computing "available" from the
        # registry's download state and the FASTQ folder on disk would not include SRR1.
        shutil.rmtree(root / "fastq" / "SRR1")

        with caplog.at_level("INFO"):
            results = extract_target_reads(
                parsed_containment=table,
                genome_id="GCF_1",
                genome_fasta=genome,
                fastq_folder=root / "fastq",
                output_folder=root / "targeted",
                threshold=0.5,
                already_done=record,
                available={"SRR2"},
            )
        assert results["SRR1"].skipped is True and results["SRR1"].mapped_records == 42
        assert "already extracted against GCF_1" in caplog.text
        assert "have no FASTQ" not in caplog.text
        mock_run.assert_not_called()


class TestAssembleExtractedReads:
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_paired_uses_1_2(self, mock_run, tmp_path):
        mock_run.side_effect = _fake_tools({})
        out, ran = assemble_extracted_reads([Path("a_1.fastq.gz"), Path("a_2.fastq.gz")], tmp_path / "asm")
        assert out == tmp_path / "asm" and ran is True
        args = mock_run.call_args.args[1]
        assert "-1" in args and "-2" in args and "-r" not in args

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_single_uses_r(self, mock_run, tmp_path):
        mock_run.side_effect = _fake_tools({})
        assemble_extracted_reads([Path("a.fastq.gz")], tmp_path / "asm", min_contig_len=500)
        args = mock_run.call_args.args[1]
        assert "-r" in args and "--min-contig-len" in args

    def test_bad_read_count_raises(self):
        with pytest.raises(ProcessingError):
            assemble_extracted_reads([], "asm")

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_default_preset_added_to_the_megahit_args(self, mock_run, tmp_path):
        mock_run.side_effect = _fake_tools({})
        assemble_extracted_reads([Path("a.fastq.gz")], tmp_path / "asm")
        args = mock_run.call_args.args[1]
        assert args[args.index("--presets") + 1] == "meta-sensitive"

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_explicit_preset_added_to_the_megahit_args(self, mock_run, tmp_path):
        mock_run.side_effect = _fake_tools({})
        assemble_extracted_reads([Path("a.fastq.gz")], tmp_path / "asm", preset="meta-large")
        args = mock_run.call_args.args[1]
        assert args[args.index("--presets") + 1] == "meta-large"

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_default_and_none_preset_omit_the_flag(self, mock_run, tmp_path):
        mock_run.side_effect = _fake_tools({})
        assemble_extracted_reads([Path("a.fastq.gz")], tmp_path / "asm", preset="default")
        assert "--presets" not in mock_run.call_args.args[1]
        assemble_extracted_reads([Path("a.fastq.gz")], tmp_path / "asm", preset=None, force=True)
        assert "--presets" not in mock_run.call_args.args[1]

    def test_preset_with_explicit_k_values_raises(self):
        with pytest.raises(ProcessingError, match="cannot be combined"):
            assemble_extracted_reads([Path("a.fastq.gz")], "asm", preset="meta-sensitive", k_flags={"k-min": 21})

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_explicit_k_values_used_without_a_preset(self, mock_run, tmp_path):
        mock_run.side_effect = _fake_tools({})
        assemble_extracted_reads(
            [Path("a.fastq.gz")], tmp_path / "asm", preset=None, k_flags={"k-min": 21, "k-max": 141}
        )
        args = mock_run.call_args.args[1]
        assert args[args.index("--k-min") + 1] == "21"
        assert args[args.index("--k-max") + 1] == "141"
        assert "--presets" not in args

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_intermediate_contigs_removed_by_default(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "asm"
            assemble_extracted_reads([Path(tmp) / "r1.fq.gz"], out)
            assert (out / "final.contigs.fa").exists()
            assert not (out / "intermediate_contigs").exists()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_intermediate_contigs_kept_when_requested(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "asm"
            assemble_extracted_reads([Path(tmp) / "r1.fq.gz"], out, keep_intermediate=True)
            assert (out / "intermediate_contigs").exists()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_a_reused_assembly_keeps_its_intermediate_contigs(self, mock_run):
        """The reuse path returns before the cleanup, so a kept folder is not removed by a
        later run that does not reassemble."""
        mock_run.side_effect = _fake_tools({})
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "asm"
            assemble_extracted_reads([Path(tmp) / "r1.fq.gz"], out, keep_intermediate=True)
            assert (out / "intermediate_contigs").exists()

            _, ran = assemble_extracted_reads([Path(tmp) / "r1.fq.gz"], out)

            assert ran is False
            assert (out / "intermediate_contigs").exists()

    def test_megahit_failure_reports_stderr(self, tmp_path):
        """A megahit failure (e.g. on a filesystem without FIFOs) must surface its stderr,
        not just the bare CalledProcessError."""
        reads = [tmp_path / "r.fastq.gz"]
        reads[0].write_bytes(b"x")
        err = subprocess.CalledProcessError(
            1, ["megahit"], output="", stderr="line1\nOSError: [Errno 45] Operation not supported\n"
        )
        with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=err):
            with pytest.raises(ProcessingError) as excinfo:
                assemble_extracted_reads(reads, tmp_path / "asm", threads=1, tmp_dir=tmp_path / "scratch")
        assert "Operation not supported" in str(excinfo.value)

    def test_megahit_receives_tmp_dir(self, tmp_path):
        state = {"calls": []}
        reads = [tmp_path / "r.fastq.gz"]
        reads[0].write_bytes(b"x")
        with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
            assemble_extracted_reads(reads, tmp_path / "asm", threads=1, tmp_dir=tmp_path / "scratch")
        megahit_call = next(c for c in state["calls"] if c[0] == "megahit")
        assert "--tmp-dir" in megahit_call[1]

    def test_tmp_dir_inside_output_dir_is_rejected(self, tmp_path):
        """megahit requires its -o directory not to already exist. A tmp_dir nested under
        output_dir would create output_dir itself (via tmp_dir.mkdir(parents=True)) before
        megahit ever runs, so it must be rejected up front rather than passed through."""
        reads = [tmp_path / "r.fastq.gz"]
        reads[0].write_bytes(b"x")
        out_dir = tmp_path / "asm"
        with patch("prepare.data.read_extraction.SecureSubprocess.run_secure") as mock_run:
            with pytest.raises(ProcessingError, match="must not"):
                assemble_extracted_reads(reads, out_dir, threads=1, tmp_dir=out_dir / "scratch")
        mock_run.assert_not_called()
        assert not out_dir.exists()


class TestResolveAssemblyThreads:
    def test_explicit_request_wins(self, monkeypatch):
        monkeypatch.setattr("prepare.data.read_extraction.platform.system", lambda: "Darwin")
        assert resolve_assembly_threads(8, 4) == 8

    def test_macos_defaults_to_single_thread(self, monkeypatch):
        monkeypatch.setattr("prepare.data.read_extraction.platform.system", lambda: "Darwin")
        assert resolve_assembly_threads(None, 4) == 1

    def test_non_macos_uses_fallback(self, monkeypatch):
        monkeypatch.setattr("prepare.data.read_extraction.platform.system", lambda: "Linux")
        assert resolve_assembly_threads(None, 4) == 4


class TestRunMinimap2Retry:
    """_run_minimap2 retries once against the FASTA when the prebuilt index fails."""

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_index_failure_retries_against_the_fasta_and_warns(self, mock_run, tmp_path, caplog):
        genome = tmp_path / "g.fna"
        genome.write_text(">s\nACGT\n")
        index_path = tmp_path / "g.sr.mmi"
        sam_path = tmp_path / "out.sam"
        reads = [tmp_path / "r1.fastq.gz"]
        ok_result = MagicMock(returncode=0, stdout="", stderr="")
        # What a minimap2 that rejects the index actually raises: a non-zero exit.
        mock_run.side_effect = [subprocess.CalledProcessError(1, "minimap2"), ok_result]

        with caplog.at_level("WARNING"):
            result = _run_minimap2("SRR1", "sr", 4, sam_path, index_path, genome, reads)

        assert result is ok_result
        assert mock_run.call_count == 2
        first_call, second_call = mock_run.call_args_list
        assert first_call.args[0] == "minimap2" and str(index_path) in first_call.args[1]
        assert second_call.args[0] == "minimap2" and str(genome) in second_call.args[1]
        assert str(index_path) not in second_call.args[1]
        assert "retrying against the FASTA directly" in caplog.text

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_a_rejected_command_also_falls_back(self, mock_run, tmp_path):
        ok_result = MagicMock(returncode=0, stdout="", stderr="")
        mock_run.side_effect = [SecurityError("command timed out"), ok_result]

        result = _run_minimap2(
            "SRR1", "sr", 4, tmp_path / "out.sam", tmp_path / "g.sr.mmi", tmp_path / "g.fna", [tmp_path / "r1.fastq"]
        )

        assert result is ok_result
        assert mock_run.call_count == 2

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_an_interrupt_is_not_swallowed_by_the_fallback(self, mock_run, tmp_path):
        """Ctrl-C is not a reason to run minimap2 a second time."""
        mock_run.side_effect = KeyboardInterrupt()

        with pytest.raises(KeyboardInterrupt):
            _run_minimap2(
                "SRR1",
                "sr",
                4,
                tmp_path / "out.sam",
                tmp_path / "g.sr.mmi",
                tmp_path / "g.fna",
                [tmp_path / "r1.fastq"],
            )

        assert mock_run.call_count == 1


class TestSummariseContigsRicherStats:
    def test_n90_gc_and_size_bucket_from_real_sequence(self, tmp_path):
        contigs = tmp_path / "contigs.fa"
        # One 1000 bp contig, all G/C, and one 100 bp contig, all A/T.
        contigs.write_text(">c1\n" + "GC" * 500 + "\n>c2\n" + "AT" * 50 + "\n")

        stats = summarise_contigs(contigs)

        assert stats["contigs"] == 2
        assert stats["total_bp"] == 1100
        assert stats["largest"] == 1000
        assert stats["n50"] == 1000
        assert stats["n90"] == 1000
        assert stats["gc"] == round(1000 / 1100, 4)
        assert stats["contigs_ge_1kb"] == 1

    def test_missing_file_returns_a_zeroed_dict_with_the_new_keys(self, tmp_path):
        stats = summarise_contigs(tmp_path / "missing.fa")
        assert stats == {
            "contigs": 0,
            "total_bp": 0,
            "n50": 0,
            "n90": 0,
            "largest": 0,
            "gc": 0.0,
            "contigs_ge_1kb": 0,
        }

    def test_header_length_is_still_honoured_for_total_bp(self, tmp_path):
        """megahit's len= header value still wins for length even though gc is read from the
        (here deliberately short) sequence lines."""
        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1 len=100\nACGT\n>c2 len=50\nACGT\n")
        stats = summarise_contigs(contigs)
        assert stats["total_bp"] == 150
        assert stats["largest"] == 100


class TestFastaLength:
    def test_counts_sequence_bases_only(self, tmp_path):
        fasta = tmp_path / "g.fna"
        fasta.write_text(">chr1 some description\nACGTACGT\nACGT\n>chr2\nTTTT\n")
        assert fasta_length(fasta) == 16


class TestMegahitInputs:
    def test_megahit_gets_singletons_as_r(self, tmp_path):
        from prepare.data.assembly import _megahit_args

        mates = [tmp_path / "G_1.fastq.gz", tmp_path / "G_2.fastq.gz"]
        singles = tmp_path / "G_s.fastq.gz"
        args = _megahit_args(mates + [singles], tmp_path / "asm", 1, None, None, None, None)
        assert args[:6] == ["-1", str(mates[0]), "-2", str(mates[1]), "-r", str(singles)]
        assert "-r" not in _megahit_args(mates, tmp_path / "asm", 1, None, None, None, None)

    def test_three_files_without_a_singleton_file_are_refused(self, tmp_path):
        from prepare.data.assembly import _megahit_args

        reads = [tmp_path / f"G_{i}.fastq.gz" for i in (1, 2, 3)]
        with pytest.raises(ProcessingError):
            _megahit_args(reads, tmp_path / "asm", 1, None, None, None, None)

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_coverage_maps_the_mates_only_when_singletons_are_present(self, mock_run, tmp_path):
        """minimap2 pairs its query files positionally, so the singleton file is left out of the back-mapping."""
        mock_run.side_effect = _fake_tools({"mapped": 4})
        contigs = tmp_path / "contigs.fa"
        contigs.write_text(">c1\n" + "A" * 100 + "\n")
        reads = [tmp_path / "G_1.fastq.gz", tmp_path / "G_2.fastq.gz", tmp_path / "G_s.fastq.gz"]
        for path in reads:
            with gzip.open(path, "wt") as handle:
                handle.write("@r\nACGT\n+\nIIII\n")
        assembly_coverage(contigs, reads, "sr", 1, tmp_path, mapped_reads=4)
        minimap_args = next(c.args[1] for c in mock_run.call_args_list if c.args[0] == "minimap2")
        assert str(reads[2]) not in minimap_args and str(reads[1]) in minimap_args


class TestAssemblyCoverage:
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_coverage_numbers_from_faked_counts(self, mock_run, tmp_path):
        mock_run.side_effect = _fake_tools({"coverage_mapped": 8})
        contigs = tmp_path / "final.contigs.fa"
        contigs.write_text(">c1 len=100\nACGT\n>c2 len=50\nACGT\n")
        reads_file = tmp_path / "r_1.fastq.gz"
        with gzip.open(reads_file, "wt") as handle:
            handle.write("@r\nACGTACGTAC\n+\nIIIIIIIIII\n" * 5)

        result = assembly_coverage(contigs, [reads_file], "sr", 4, tmp_path, mapped_reads=40)

        assert result["reads_mapped"] == 8
        assert result["mapping_rate"] == pytest.approx(0.2)
        assert result["mean_depth_estimate"] == pytest.approx(8 * 10 / 150)
        assert not (tmp_path / "coverage.sam").exists()
        assert not (tmp_path / "coverage.bam").exists()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_mapping_rate_is_none_when_mapped_reads_is_unknown_or_zero(self, mock_run, tmp_path):
        mock_run.side_effect = _fake_tools({})
        contigs = tmp_path / "final.contigs.fa"
        contigs.write_text(">c1 len=10\nACGT\n")
        reads_file = tmp_path / "r.fastq.gz"
        with gzip.open(reads_file, "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n")

        assert assembly_coverage(contigs, [reads_file], "sr", 2, tmp_path, mapped_reads=None)["mapping_rate"] is None
        assert assembly_coverage(contigs, [reads_file], "sr", 2, tmp_path, mapped_reads=0)["mapping_rate"] is None

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_mean_depth_is_none_for_an_assembly_with_no_contigs(self, mock_run, tmp_path):
        """No contigs means no depth to report, which 0.0 would misstate as a measurement."""
        mock_run.side_effect = _fake_tools({})
        contigs = tmp_path / "final.contigs.fa"
        contigs.write_text("")
        reads_file = tmp_path / "r.fastq.gz"
        with gzip.open(reads_file, "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n")

        result = assembly_coverage(contigs, [reads_file], "sr", 2, tmp_path, mapped_reads=10)

        assert result["mean_depth_estimate"] is None

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_minimap2_and_samtools_args_shape(self, mock_run, tmp_path):
        mock_run.side_effect = _fake_tools({})
        contigs = tmp_path / "final.contigs.fa"
        contigs.write_text(">c1 len=10\nACGT\n")
        reads_file = tmp_path / "r.fastq.gz"
        with gzip.open(reads_file, "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n")

        assembly_coverage(contigs, [reads_file], "map-ont", 3, tmp_path, mapped_reads=10)

        minimap2_call = next(c for c in mock_run.call_args_list if c.args[0] == "minimap2")
        args = minimap2_call.args[1]
        assert args[:3] == ["-a", "-x", "map-ont"]
        assert args[args.index("-t") + 1] == "3"
        assert str(contigs) in args and str(reads_file) in args

        filter_call = next(c for c in mock_run.call_args_list if c.args[0] == "samtools" and c.args[1][0] == "view")
        fargs = filter_call.args[1]
        assert fargs[fargs.index("-F") + 1] == "0x904"
        assert fargs[fargs.index("-@") + 1] == "3"


class TestReferenceCoverage:
    """Breadth and mean depth against the reference genome, from the filtered BAM."""

    def _extract(self, root, table, genome, **kwargs):
        return extract_target_reads(
            parsed_containment=table,
            genome_id="GCF_1",
            genome_fasta=genome,
            fastq_folder=root / "fastq",
            output_folder=root / "targeted",
            threshold=0.5,
            **kwargs,
        )

    def test_summary_aggregates_contigs_by_length(self, tmp_path):
        tsv = tmp_path / "cov.tsv"
        tsv.write_text(coverage_table([("c1", 100, 50, 2.0), ("c2", 300, 300, 10.0)]))
        summary = summarise_coverage_table(tsv)
        assert summary["breadth"] == 0.875
        assert summary["mean_depth"] == 8.0
        assert summary["covered_bases"] == 350
        assert summary["reference_bp"] == 400

    def test_summary_rounds_to_four_decimals(self, tmp_path):
        tsv = tmp_path / "cov.tsv"
        tsv.write_text(coverage_table([("c1", 3, 1, 1.0)]))
        summary = summarise_coverage_table(tsv)
        assert summary["breadth"] == 0.3333

    def test_summary_of_an_empty_table_reports_none(self, tmp_path):
        tsv = tmp_path / "cov.tsv"
        tsv.write_text(coverage_table([]))
        summary = summarise_coverage_table(tsv)
        assert summary["reference_bp"] == 0
        assert summary["breadth"] is None and summary["mean_depth"] is None

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_reference_coverage_returns_the_summary_and_table_path(self, mock_run, tmp_path):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        bam = tmp_path / "GCF_1.mapped.bam"
        bam.write_bytes(b"")
        (tmp_path / "out").mkdir()
        (tmp_path / "sam").mkdir()
        coverage = reference_coverage("SRR1", "GCF_1", bam, tmp_path / "out", tmp_path / "sam", 2)
        sorted_bam = tmp_path / "sam" / "GCF_1.mapped.sorted.bam"
        tsv = tmp_path / "out" / "GCF_1_coverage.tsv"
        assert coverage["breadth"] == 0.875 and coverage["mean_depth"] == 8.0
        assert coverage["coverage_tsv"] == tsv
        assert tsv.exists()
        assert not sorted_bam.exists()
        assert state["calls"][0][1] == ["sort", "-@", "2", "-o", str(sorted_bam), str(bam)]
        assert state["calls"][1][1] == ["coverage", "-o", str(tsv), str(sorted_bam)]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_extraction_writes_the_table_and_records_breadth(self, mock_run, tmp_path):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        root, table, genome = _make_tree(tmp_path, paired=True)
        results = self._extract(root, table, genome)
        sample_dir = root / "targeted" / "SRR1"
        assert (sample_dir / "GCF_1_coverage.tsv").exists()
        coverage = results["SRR1"].coverage
        assert coverage["breadth"] == 0.875 and coverage["mean_depth"] == 8.0
        assert coverage["covered_bases"] == 350 and coverage["reference_bp"] == 400
        assert Path(coverage["coverage_tsv"]) == sample_dir / "GCF_1_coverage.tsv"
        assert list(root.rglob("*.sorted.bam")) == []
        assert not (sample_dir / "GCF_1.mapped.bam").exists()
        subcommands = [args[0] for tool, args in state["calls"] if tool == "samtools"]
        assert subcommands.index("fastq") < subcommands.index("sort") < subcommands.index("coverage")

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_zero_mapped_sample_runs_neither_sort_nor_coverage(self, mock_run, tmp_path):
        state = {"mapped": 0, "mapped_total": 0}
        mock_run.side_effect = _fake_tools(state)
        root, table, genome = _make_tree(tmp_path, paired=True)
        results = self._extract(root, table, genome)
        assert results["SRR1"].coverage is None
        subcommands = [args[0] for tool, args in state["calls"] if tool == "samtools"]
        assert "sort" not in subcommands and "coverage" not in subcommands

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_zero_kept_after_filter_runs_neither_sort_nor_coverage(self, mock_run, tmp_path):
        state = {"mapped": 0, "mapped_total": 5}
        mock_run.side_effect = _fake_tools(state)
        root, table, genome = _make_tree(tmp_path, paired=True)
        results = self._extract(root, table, genome)
        assert results["SRR1"].coverage is None
        subcommands = [args[0] for tool, args in state["calls"] if tool == "samtools"]
        assert "sort" not in subcommands and "coverage" not in subcommands

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_coverage_failure_is_a_warning_and_the_fastq_still_lands(self, mock_run, tmp_path, caplog):
        state = {"coverage_fail": True}
        mock_run.side_effect = _fake_tools(state)
        root, table, genome = _make_tree(tmp_path, paired=True)
        with caplog.at_level(logging.WARNING, logger="prepare.data.read_extraction"):
            results = self._extract(root, table, genome)
        sample_dir = root / "targeted" / "SRR1"
        assert results["SRR1"].coverage is None
        assert [p.name for p in results["SRR1"].files] == ["GCF_1_1.fastq.gz", "GCF_1_2.fastq.gz"]
        assert all(p.exists() for p in results["SRR1"].files)
        assert not (sample_dir / "GCF_1_coverage.tsv").exists()
        assert list(root.rglob("*.sorted.bam")) == []
        assert any(
            r.levelno == logging.WARNING and "SRR1: reference coverage skipped" in r.getMessage()
            for r in caplog.records
        )

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_sorted_bam_is_placed_under_the_temp_folder(self, mock_run, tmp_path):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        root, table, genome = _make_tree(tmp_path, paired=True)
        temp_folder = root / "scratch"
        self._extract(root, table, genome, temp_folder=temp_folder)
        sort_args = next(args for tool, args in state["calls"] if tool == "samtools" and args[0] == "sort")
        assert Path(sort_args[sort_args.index("-o") + 1]) == temp_folder / "GCF_1.mapped.sorted.bam"
        assert list(temp_folder.glob("*.sorted.bam")) == []

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_forced_rerun_with_zero_mapped_removes_an_earlier_table(self, mock_run, tmp_path):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        root, table, genome = _make_tree(tmp_path, paired=True)
        tsv = root / "targeted" / "SRR1" / "GCF_1_coverage.tsv"
        self._extract(root, table, genome)
        assert tsv.exists()
        state.update(mapped=0, mapped_total=0)
        results = self._extract(root, table, genome)
        assert results["SRR1"].coverage is None
        assert not tsv.exists()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_forced_rerun_with_zero_kept_removes_an_earlier_table(self, mock_run, tmp_path):
        state = {}
        mock_run.side_effect = _fake_tools(state)
        root, table, genome = _make_tree(tmp_path, paired=True)
        tsv = root / "targeted" / "SRR1" / "GCF_1_coverage.tsv"
        self._extract(root, table, genome)
        assert tsv.exists()
        state.update(mapped=0, mapped_total=5)
        results = self._extract(root, table, genome)
        assert results["SRR1"].coverage is None
        assert not tsv.exists()


class TestOldMinimap2WithoutSamHitOnly:
    """A minimap2 older than 2.17 rejects --sam-hit-only; the error names the minimum version."""

    OLD_STDERR = '[ERROR] unknown option in "--sam-hit-only"\n'

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_unknown_option_is_a_configuration_error(self, mock_run, tmp_path):
        from prepare.core.exceptions import ConfigurationError

        state = {"minimap2_fail_stderr": self.OLD_STDERR}
        mock_run.side_effect = _fake_tools(state)
        with pytest.raises(ConfigurationError, match="minimap2 2.17 or later") as caught:
            _run_minimap2(
                "SRR1", "sr", 4, tmp_path / "out.sam", tmp_path / "g.sr.mmi", tmp_path / "g.fna", [tmp_path / "r.fq"]
            )
        assert "unknown option" in str(caught.value)
        # No retry against the FASTA: the same option would be rejected again.
        assert [call[0] for call in state["calls"]] == ["minimap2"]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_other_failures_still_fall_back_to_the_fasta(self, mock_run, tmp_path):
        state = {"minimap2_fail_stderr": "[ERROR] failed to open file 'g.sr.mmi'\n"}
        mock_run.side_effect = _fake_tools(state)
        with pytest.raises(subprocess.CalledProcessError):
            _run_minimap2(
                "SRR1", "sr", 4, tmp_path / "out.sam", tmp_path / "g.sr.mmi", tmp_path / "g.fna", [tmp_path / "r.fq"]
            )
        assert len(state["calls"]) == 2


class TestMegahitVersion:
    """The version probe always passes VERSION_PROBE_TIMEOUT, not the run's subprocess_timeout."""

    @patch("prepare.utils.tools.shutil.which", return_value="/usr/bin/megahit")
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_passes_the_fixed_probe_timeout(self, mock_run, _which):
        from prepare.core.constants import VERSION_PROBE_TIMEOUT

        mock_run.return_value = MagicMock(stdout="MEGAHIT v1.2.9\n", stderr="", returncode=0)
        assert megahit_version() == "MEGAHIT v1.2.9"
        mock_run.assert_called_once_with("megahit", ["--version"], timeout=VERSION_PROBE_TIMEOUT, check=False)

    @patch("prepare.utils.tools.shutil.which", return_value="/usr/bin/megahit")
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_returns_empty_when_the_tool_cannot_be_run(self, mock_run, _which):
        mock_run.side_effect = SecurityError("not installed")
        assert megahit_version() == ""
