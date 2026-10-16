"""Tests for the extract_target_reads CLI command."""

import argparse
import functools
import gzip
import json
import logging
import tempfile
import threading
from pathlib import Path
from unittest.mock import patch

import pytest

from prepare.cli.commands.extraction_assembly import assemble_samples
from prepare.cli.commands.read_extraction import ExtractTargetReadsCommand
from prepare.cli.commands.status import StatusCommand
from prepare.core.exceptions import ProcessingError, SecurityError
from prepare.data.read_extraction import ExtractionResult
from prepare.data.registry import load_registry, record_download, record_extraction, save_registry
from helpers_extraction import _fake_tools, tools_present  # noqa: F401 (autouse fixture)
from prepare.utils import tools


def _args(tmp, **kwargs):
    base = dict(
        parsed_containment="parsed_containment.txt",
        genome_id="GCF_1",
        genome_fasta="GCF_1.fna",
        fastq_folder="fastq",
        output_folder="targeted",
        threshold=0.1,
        preset="sr",
        threads=4,
        min_mapq=0,
        temp_folder=None,
        allow_truncated=False,
        debug_keep_sam=False,
        assemble=False,
        assembly_threads=None,
        min_contig_len=None,
        assembly_preset="meta-sensitive",
        keep_intermediate=False,
        no_coverage=False,
        dry_run=False,
        force=False,
        registry=str(Path(tmp) / "registry.json"),
        data_root=None,
    )
    base.update(kwargs)
    return argparse.Namespace(**base)


def _tree(tmp):
    root = Path(tmp)
    table = root / "parsed_containment.txt"
    table.write_text("\tGCF_1\nSRR1\t0.9\nSRR2\t0.05\n")
    d = root / "fastq" / "SRR1"
    d.mkdir(parents=True)
    (d / "SRR1_1.fastq.gz").write_text("x")
    (d / "SRR1_2.fastq.gz").write_text("x")
    genome = root / "GCF_1.fna"
    genome.write_text(">s\nACGT\n")
    return root, table, genome


def _bootstrapped_tree(tmp):
    """A project whose extraction already happened on disk, with no registry yet."""
    root = Path(tmp)
    table = root / "parsed_containment.txt"
    table.write_text("\tGCF_1\nSRR1\t0.9\n")
    reads = root / "fastq" / "SRR1"
    reads.mkdir(parents=True)
    (reads / "SRR1_1.fastq.gz").write_text("x")
    (reads / "SRR1_2.fastq.gz").write_text("x")
    genome = root / "genomes" / "GCF_1.fna"
    genome.parent.mkdir(parents=True)
    genome.write_text(">s\nACGT\n")
    extracted = root / "targeted" / "SRR1"
    extracted.mkdir(parents=True)
    for name in ("GCF_1_1.fastq.gz", "GCF_1_2.fastq.gz"):
        with gzip.open(extracted / name, "wt") as handle:
            handle.write("@r1\nACGT\n+\nIIII\n")
    return root, table, genome


def _status_args(root, **overrides):
    base = dict(
        fastq_folder=str(root / "fastq"),
        metadata_folder=str(root / "metadata"),
        genomes_folder=str(root / "genomes"),
        targeted_folder=str(root / "targeted"),
        matches_folder=str(root / "matches"),
        registry=str(root / "prepare_registry.json"),
        data_root=None,
        accessions_file=None,
        parsed_containment=None,
        stage=None,
        genome=None,
        init=True,
        reconcile=False,
        export_tsv=None,
        next=False,
        list_missing=False,
        json=True,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _gzip_tree(tmp, mate1_reads=3, mate2_reads=3):
    """Like _tree, but with real gzip FASTQ content so count_fastq_reads can read it."""
    root = Path(tmp)
    table = root / "parsed_containment.txt"
    table.write_text("\tGCF_1\nSRR1\t0.9\nSRR2\t0.05\n")
    d = root / "fastq" / "SRR1"
    d.mkdir(parents=True)
    for name, n in (("SRR1_1.fastq.gz", mate1_reads), ("SRR1_2.fastq.gz", mate2_reads)):
        with gzip.open(d / name, "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n" * n)
    genome = root / "GCF_1.fna"
    genome.write_text(">s\nACGT\n")
    return root, table, genome


def _two_sample_tree(tmp):
    """Like _tree, but with two samples above the threshold."""
    root = Path(tmp)
    table = root / "parsed_containment.txt"
    table.write_text("\tGCF_1\nSRR1\t0.9\nSRR2\t0.8\n")
    for acc in ("SRR1", "SRR2"):
        d = root / "fastq" / acc
        d.mkdir(parents=True)
        (d / f"{acc}_1.fastq.gz").write_text("x")
        (d / f"{acc}_2.fastq.gz").write_text("x")
    genome = root / "GCF_1.fna"
    genome.write_text(">s\nACGT\n")
    return root, table, genome


class TestExtractTargetReadsCommand:
    def test_command_properties(self):
        cmd = ExtractTargetReadsCommand()
        assert cmd.name == "extract_target_reads"
        assert "target" in cmd.help.lower()

    def test_configure_parser_timeout_flag_defaults_to_none_and_parses(self):
        """--timeout defaults to None (run_secure then falls back to the active setting)."""
        cmd = ExtractTargetReadsCommand()
        parser = argparse.ArgumentParser()
        cmd.configure_parser(parser)
        required = ["--parsed-containment", "p.txt", "--genome-id", "GCF_1", "--genome-fasta", "g.fna"]

        args = parser.parse_args(required)
        assert args.timeout is None

        args = parser.parse_args(required + ["--timeout", "120"])
        assert args.timeout == 120.0

    @patch("prepare.utils.tools.shutil.which", return_value=None)
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_missing_minimap2_or_samtools_exits_3_before_any_work(self, mock_run, _which):
        """A missing minimap2/samtools is reported once, with an install hint, before any
        tool ever runs -- not as a raw subprocess error partway through extraction."""
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                )
            )
        assert rc == 3
        mock_run.assert_not_called()

    @patch("prepare.utils.tools.shutil.which")
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_missing_minimap2_logs_install_hint(self, mock_run, mock_which, caplog):
        mock_which.side_effect = lambda name: None if name == "minimap2" else f"/usr/bin/{name}"
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            with caplog.at_level("ERROR"):
                rc = cmd.execute(
                    _args(
                        tmp,
                        parsed_containment=str(table),
                        genome_fasta=str(genome),
                        fastq_folder=str(root / "fastq"),
                    )
                )
        assert rc == 3
        assert "minimap2 not found on PATH" in caplog.text
        assert "conda install -c conda-forge -c bioconda 'minimap2>=2.17'" in caplog.text

    @patch("prepare.utils.tools.shutil.which")
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_missing_megahit_only_checked_when_assembling(self, mock_run, mock_which):
        """megahit is only required (and checked) when --assemble is set."""
        mock_run.side_effect = _fake_tools({})
        mock_which.side_effect = lambda name: None if name == "megahit" else f"/usr/bin/{name}"
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    assemble=False,
                )
            )
        assert rc == 0
        mock_run.assert_called()

    @patch("prepare.utils.tools.shutil.which")
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_missing_megahit_exits_3_when_assembling(self, mock_run, mock_which):
        mock_which.side_effect = lambda name: None if name == "megahit" else f"/usr/bin/{name}"
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    assemble=True,
                )
            )
        assert rc == 3
        mock_run.assert_not_called()

    @patch("prepare.utils.tools.shutil.which", return_value=None)
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_dry_run_skips_the_tool_check(self, mock_run, _which):
        """--dry-run never touches any tool, so the pre-flight check is skipped entirely."""
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    dry_run=True,
                )
            )
        assert rc == 0
        mock_run.assert_not_called()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_extracts(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                )
            )
        assert rc == 0
        assert mock_run.called

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_returns_1_when_no_sample_yields_reads(self, mock_run, caplog):
        mock_run.side_effect = _fake_tools({"mapped": 0})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            with caplog.at_level("ERROR"):
                rc = cmd.execute(
                    _args(
                        tmp,
                        parsed_containment=str(table),
                        genome_fasta=str(genome),
                        fastq_folder=str(root / "fastq"),
                        output_folder=str(root / "targeted"),
                        threshold=0.5,
                    )
                )
        assert rc == 1
        assert "No reads mapped to GCF_1" in caplog.text

    def test_execute_returns_1_when_no_sample_meets_threshold(self, caplog):
        """threshold above every sample's containment -> nothing selected."""
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            with caplog.at_level("ERROR"):
                rc = cmd.execute(
                    _args(
                        tmp,
                        parsed_containment=str(table),
                        genome_fasta=str(genome),
                        fastq_folder=str(root / "fastq"),
                        output_folder=str(root / "targeted"),
                        threshold=2.0,
                    )
                )
        assert rc == 1
        assert "No sample meets containment >= 2.0 for GCF_1" in caplog.text

    def test_execute_returns_1_when_no_fastq_for_selected_samples(self, caplog):
        """A sample is selected but --fastq-folder does not contain its reads."""
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            with caplog.at_level("ERROR"):
                rc = cmd.execute(
                    _args(
                        tmp,
                        parsed_containment=str(table),
                        genome_fasta=str(genome),
                        fastq_folder=str(root / "no-such-fastq"),
                        output_folder=str(root / "targeted"),
                        threshold=0.5,
                    )
                )
        assert rc == 1
        assert "No FASTQ files found for the 1 selected sample(s) under" in caplog.text

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_dry_run(self, mock_run):
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    dry_run=True,
                )
            )
        assert rc == 0
        mock_run.assert_not_called()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_with_assembly(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    assemble=True,
                )
            )
        assert rc == 0
        tools = [c.args[0] for c in mock_run.call_args_list]
        assert "megahit" in tools

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_assembly_single_thread_on_macos(self, mock_run, monkeypatch):
        mock_run.side_effect = _fake_tools({})
        monkeypatch.setattr("prepare.data.read_extraction.platform.system", lambda: "Darwin")
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    threads=4,
                    assemble=True,
                )
            )
        assert rc == 0
        megahit_call = next(c for c in mock_run.call_args_list if c.args[0] == "megahit" and "-o" in c.args[1])
        args = megahit_call.args[1]
        assert args[args.index("--num-cpu-threads") + 1] == "1"

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_assembly_threads_override_on_macos(self, mock_run, monkeypatch):
        mock_run.side_effect = _fake_tools({})
        monkeypatch.setattr("prepare.data.read_extraction.platform.system", lambda: "Darwin")
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    threads=4,
                    assemble=True,
                    assembly_threads=6,
                )
            )
        assert rc == 0
        megahit_call = next(c for c in mock_run.call_args_list if c.args[0] == "megahit" and "-o" in c.args[1])
        args = megahit_call.args[1]
        assert args[args.index("--num-cpu-threads") + 1] == "6"

    def test_execute_missing_genome_column_returns_1(self):
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_id="GCF_absent",
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                )
            )
        assert rc == 1

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_records_extraction_and_assembly(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    assemble=True,
                    registry=str(registry_file),
                )
            )
            assert rc == 0
            data = json.loads(registry_file.read_text())
        extraction = data["datasets"]["SRR1"]["extractions"]["GCF_1"]
        assert extraction["mapped_reads"] > 0
        assert extraction["assembly"]["contigs"] == 2

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_records_extraction_and_assembly_timing(self, mock_run, monkeypatch):
        """Each extraction and each megahit run is timed; a rerun that skips both keeps the times."""
        import itertools
        import re

        ticks = itertools.count(0.0, 1.5)
        monkeypatch.setattr("prepare.data.registry_timing.monotonic", lambda: next(ticks))
        mock_run.side_effect = _fake_tools({})
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                assemble=True,
                registry=str(registry_file),
            )
            assert ExtractTargetReadsCommand().execute(args) == 0
            first = json.loads(registry_file.read_text())["datasets"]["SRR1"]["extractions"]["GCF_1"]
            assert ExtractTargetReadsCommand().execute(args) == 0
            second = json.loads(registry_file.read_text())["datasets"]["SRR1"]["extractions"]["GCF_1"]
        iso = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\+00:00$")
        assert iso.match(first["started"]) and first["seconds"] == 1.5
        assert iso.match(first["assembly"]["started"]) and first["assembly"]["seconds"] == 1.5
        assert (second["started"], second["seconds"]) == (first["started"], first["seconds"])
        assert second["assembly"]["seconds"] == first["assembly"]["seconds"]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_records_usage_for_extraction_and_assembly(self, mock_run):
        """Extraction and assembly are recorded as store catalogue usage for this genome."""
        from prepare.data.registry import load_registry as _load, save_registry as _save
        from prepare.store.catalog import Catalog
        from prepare.store.layout import init_store, store_paths

        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            store_root = root / "store"
            init_store(store_root)
            registry_file = root / "registry.json"

            registry = _load(registry_file)
            registry.project = {"id": "proj1", "name": "demo", "path": str(root), "created": "now"}
            _save(registry)

            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    assemble=True,
                    registry=str(registry_file),
                    data_root=str(store_root),
                )
            )
            assert rc == 0

            with Catalog(store_paths(store_root)) as catalog:
                catalog.migrate()
                rows = {
                    (r["accession"], r["genome_id"], r["stage"])
                    for r in catalog.conn.execute("SELECT accession, genome_id, stage FROM usage").fetchall()
                }
        assert ("SRR1", "GCF_1", "extracted") in rows
        assert ("SRR1", "GCF_1", "assembled") in rows

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_catalog_failure_leaves_extraction_outcome_unchanged(self, mock_run):
        """A broken catalogue write never changes the extraction's exit code or registry record."""
        from prepare.data.registry import load_registry as _load, save_registry as _save
        from prepare.store.layout import init_store

        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            store_root = root / "store"
            init_store(store_root)
            registry_file = root / "registry.json"

            registry = _load(registry_file)
            registry.project = {"id": "proj1", "name": "demo", "path": str(root), "created": "now"}
            _save(registry)

            with patch("prepare.store.usage.catalog_write", side_effect=RuntimeError("locked")):
                rc = cmd.execute(
                    _args(
                        tmp,
                        parsed_containment=str(table),
                        genome_fasta=str(genome),
                        fastq_folder=str(root / "fastq"),
                        output_folder=str(root / "targeted"),
                        threshold=0.5,
                        assemble=True,
                        registry=str(registry_file),
                        data_root=str(store_root),
                    )
                )
            assert rc == 0
            data = json.loads(registry_file.read_text())
        extraction = data["datasets"]["SRR1"]["extractions"]["GCF_1"]
        assert extraction["mapped_reads"] > 0
        assert extraction["assembly"]["contigs"] == 2

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_zero_mapped_sample_is_never_assembled(self, mock_run):
        """Leftover files of a zero-mapped sample must not reach megahit."""
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            leftovers = root / "targeted" / "SRR1"
            leftovers.mkdir(parents=True)
            names = ["GCF_1_1.fastq.gz", "GCF_1_2.fastq.gz", "GCF_1_s.fastq.gz"]
            for name in names:
                (leftovers / name).write_text("")
            seeded = load_registry(registry_file)
            record_extraction(
                seeded,
                "SRR1",
                "GCF_1",
                [leftovers / name for name in names],
                0,
                False,
                {"genome_fasta": str(genome), "preset": "sr", "threshold": 0.5},
            )
            save_registry(seeded)

            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    assemble=True,
                    registry=str(registry_file),
                )
            )
        assert rc == 1
        assert "megahit" not in [c.args[0] for c in mock_run.call_args_list]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_dry_run_reports_would_be_skips(self, mock_run, caplog):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                registry=str(registry_file),
            )
            assert cmd.execute(args) == 0
            with caplog.at_level("INFO"):
                rc = cmd.execute(_args(tmp, **{**vars(args), "dry_run": True}))
        assert rc == 0
        assert "already extracted, would be skipped: SRR1" in caplog.text

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_rerun_keeps_the_recorded_assembly(self, mock_run):
        """A rerun that does not run megahit again leaves the recorded assembly untouched."""
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                assemble=True,
                registry=str(registry_file),
            )
            assert cmd.execute(args) == 0
            # Stand in for the megahit version that produced the assembly on the first run.
            seeded = load_registry(registry_file)
            seeded.datasets["SRR1"]["extractions"]["GCF_1"]["assembly"]["version"] = "v1.2.9"
            save_registry(seeded)

            assert cmd.execute(args) == 0
            second = json.loads(registry_file.read_text())["datasets"]["SRR1"]["extractions"]["GCF_1"]["assembly"]
        assert second["version"] == "v1.2.9"

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_bootstrapped_registry_skips_the_extraction(self, mock_run):
        """status --init infers records without parameters; a rerun must still skip the mapping."""
        mock_run.side_effect = _fake_tools({})
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _bootstrapped_tree(tmp)
            registry_file = root / "prepare_registry.json"
            assert StatusCommand().execute(_status_args(root)) == 0
            assert registry_file.exists()
            calls_before = len(mock_run.call_args_list)

            rc = ExtractTargetReadsCommand().execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    registry=str(registry_file),
                )
            )
        assert rc == 0
        assert "minimap2" not in [c.args[0] for c in mock_run.call_args_list[calls_before:]]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_records_are_checkpointed_when_a_later_assembly_fails(self, mock_run):
        """Each extraction and assembly is written as it completes, so a later failure keeps them."""
        mock_run.side_effect = _fake_tools({})

        def fail_on_second(reads, out_dir, **kwargs):
            if "SRR2" in str(out_dir):
                raise ProcessingError("megahit crashed")
            Path(out_dir).mkdir(parents=True, exist_ok=True)
            (Path(out_dir) / "final.contigs.fa").write_text(">c1 len=100\nACGT\n")
            return Path(out_dir), True

        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _two_sample_tree(tmp)
            registry_file = root / "registry.json"
            with patch("prepare.cli.commands.extraction_assembly.assemble_extracted_reads", side_effect=fail_on_second):
                rc = cmd.execute(
                    _args(
                        tmp,
                        parsed_containment=str(table),
                        genome_fasta=str(genome),
                        fastq_folder=str(root / "fastq"),
                        output_folder=str(root / "targeted"),
                        threshold=0.5,
                        assemble=True,
                        registry=str(registry_file),
                    )
                )
            assert rc == 1
            datasets = json.loads(registry_file.read_text())["datasets"]
        assert datasets["SRR1"]["extractions"]["GCF_1"]["mapped_reads"] > 0
        assert datasets["SRR2"]["extractions"]["GCF_1"]["mapped_reads"] > 0
        assert datasets["SRR1"]["extractions"]["GCF_1"]["assembly"]["contigs"] == 1
        assert datasets["SRR2"]["extractions"]["GCF_1"].get("assembly") is None

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_second_run_skips_and_returns_0(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                registry=str(registry_file),
            )
            rc = cmd.execute(args)
            assert rc == 0
            calls_before = len(mock_run.call_args_list)

            rc = cmd.execute(args)
            assert rc == 0
            tools_after_second_run = [c.args[0] for c in mock_run.call_args_list[calls_before:]]
            assert "minimap2" not in tools_after_second_run
            assert "samtools" not in tools_after_second_run

            forced_args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                registry=str(registry_file),
                force=True,
            )
            calls_before_forced = len(mock_run.call_args_list)
            rc = cmd.execute(forced_args)
            assert rc == 0
            tools_after_forced_run = [c.args[0] for c in mock_run.call_args_list[calls_before_forced:]]
            assert "minimap2" in tools_after_forced_run
            assert "samtools" in tools_after_forced_run

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_second_run_skips_after_the_project_directory_is_renamed(self, mock_run, tmp_path):
        """The registry stores paths relative to itself, so a project directory renamed
        between two runs still lets the second run recognise the earlier extraction."""
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        old_root = tmp_path / "a"
        old_root.mkdir()
        root, table, genome = _tree(old_root)
        registry_file = root / "registry.json"
        rc = cmd.execute(
            _args(
                str(old_root),
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                registry=str(registry_file),
            )
        )
        assert rc == 0
        calls_before = len(mock_run.call_args_list)

        new_root = tmp_path / "b"
        old_root.rename(new_root)
        registry_file = new_root / "registry.json"
        rc = cmd.execute(
            _args(
                str(new_root),
                parsed_containment=str(new_root / "parsed_containment.txt"),
                genome_fasta=str(new_root / "GCF_1.fna"),
                fastq_folder=str(new_root / "fastq"),
                output_folder=str(new_root / "targeted"),
                threshold=0.5,
                registry=str(registry_file),
            )
        )
        assert rc == 0
        tools_after_rename = [c.args[0] for c in mock_run.call_args_list[calls_before:]]
        assert "minimap2" not in tools_after_rename
        assert "samtools" not in tools_after_rename

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_records_filter_flags_min_mapq_and_index(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    min_mapq=20,
                    registry=str(registry_file),
                )
            )
            assert rc == 0
            data = json.loads(registry_file.read_text())
        extraction = data["datasets"]["SRR1"]["extractions"]["GCF_1"]
        assert extraction["filter_flags"] == "0x904"
        assert extraction["min_mapq"] == 20
        assert extraction["index"].endswith("GCF_1.sr.mmi")
        assert ".index" in extraction["index"]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_caches_mate_read_counts_in_the_registry(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _gzip_tree(tmp, mate1_reads=3, mate2_reads=3)
            registry_file = root / "registry.json"
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    registry=str(registry_file),
                )
            )
            assert rc == 0
            data = json.loads(registry_file.read_text())
        assert data["datasets"]["SRR1"]["download"]["mate_reads"] == [3, 3]
        signature = data["datasets"]["SRR1"]["download"]["mate_reads_signature"]
        assert [entry[0] for entry in signature] == ["SRR1_1.fastq.gz", "SRR1_2.fastq.gz"]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_recounts_mates_after_the_files_change(self, mock_run):
        """A pair counted from a truncated download must not survive the re-download that
        fixed it: an unequal stale pair would force single-end mapping of a complete pair."""
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _gzip_tree(tmp, mate1_reads=3, mate2_reads=3)
            registry_file = root / "registry.json"
            registry = load_registry(registry_file)
            # What an earlier run recorded from the truncated copy of these files.
            registry.datasets["SRR1"] = {
                "download": {
                    "attempts": 1,
                    "mate_reads": [10, 9],
                    "mate_reads_signature": [["SRR1_1.fastq.gz", 1, 1.0], ["SRR1_2.fastq.gz", 1, 1.0]],
                }
            }
            save_registry(registry)

            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    registry=str(registry_file),
                )
            )
            assert rc == 0
            data = json.loads(registry_file.read_text())

        assert data["datasets"]["SRR1"]["download"]["mate_reads"] == [3, 3]

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_reuses_the_store_statistics_record_for_the_mate_counts(self, mock_run):
        """A store-backed accession already has exact per-file read counts; the mate pre-count
        reuses them instead of streaming both files again."""
        from prepare.store.sidecar import Sidecar, write_sidecar
        from prepare.store.stats import compute_dataset_stats

        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _gzip_tree(tmp, mate1_reads=3, mate2_reads=3)
            store_acc_dir = root / "store" / "sra" / "SRR1"
            store_acc_dir.mkdir(parents=True)
            project_acc_dir = root / "fastq" / "SRR1"
            for path in sorted(project_acc_dir.iterdir()):
                path.rename(store_acc_dir / path.name)
            project_acc_dir.rmdir()
            project_acc_dir.symlink_to(store_acc_dir)

            record = compute_dataset_stats(sorted(store_acc_dir.iterdir()), use_seqkit=False)
            sidecar = Sidecar(accession="SRR1", state="complete", stats=record)
            write_sidecar(store_acc_dir / "SRR1.json", sidecar)

            registry_file = root / "registry.json"
            with patch("prepare.cli.commands.read_extraction.count_fastq_reads") as counted:
                rc = cmd.execute(
                    _args(
                        tmp,
                        parsed_containment=str(table),
                        genome_fasta=str(genome),
                        fastq_folder=str(root / "fastq"),
                        output_folder=str(root / "targeted"),
                        threshold=0.5,
                        registry=str(registry_file),
                    )
                )

        assert rc == 0
        counted.assert_not_called()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_does_not_count_mates_of_a_sample_it_will_skip(self, mock_run):
        """An already extracted sample is skipped, so reading both its mate files is pure cost."""
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _gzip_tree(tmp, mate1_reads=3, mate2_reads=3)
            registry_file = root / "registry.json"
            extracted = root / "targeted" / "SRR1"
            extracted.mkdir(parents=True)
            for name in ("GCF_1_1.fastq.gz", "GCF_1_2.fastq.gz"):
                with gzip.open(extracted / name, "wt") as handle:
                    handle.write("@r1\nACGT\n+\nIIII\n")
            registry = load_registry(registry_file)
            record_extraction(
                registry,
                "SRR1",
                "GCF_1",
                [extracted / "GCF_1_1.fastq.gz", extracted / "GCF_1_2.fastq.gz"],
                7,
                False,
                {"genome_fasta": str(genome), "preset": "sr", "threshold": 0.5},
            )
            save_registry(registry)

            with patch("prepare.cli.commands.read_extraction.count_fastq_reads") as counted:
                rc = cmd.execute(
                    _args(
                        tmp,
                        parsed_containment=str(table),
                        genome_fasta=str(genome),
                        fastq_folder=str(root / "fastq"),
                        output_folder=str(root / "targeted"),
                        threshold=0.5,
                        registry=str(registry_file),
                    )
                )

        assert rc == 0
        counted.assert_not_called()

    def test_min_mapq_rejects_a_negative_value(self):
        """A negative --min-mapq would reach samtools as an unknown flag, after minimap2 ran."""
        parser = argparse.ArgumentParser()
        ExtractTargetReadsCommand().configure_parser(parser)

        required = ["--parsed-containment", "t.txt", "--genome-id", "G", "--genome-fasta", "g.fna"]
        assert parser.parse_args(required + ["--min-mapq", "0"]).min_mapq == 0
        assert parser.parse_args(required + ["--min-mapq", "20"]).min_mapq == 20
        with pytest.raises(SystemExit):
            parser.parse_args(required + ["--min-mapq", "-5"])

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_skips_a_truncated_download_unless_allowed(self, mock_run, caplog):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            registry = load_registry(registry_file)
            registry.datasets["SRR1"] = {
                "download": {"complete": {"verdict": "truncated", "reads_r1": 5, "expected_spots": 20}}
            }
            save_registry(registry)

            with caplog.at_level("WARNING"):
                rc = cmd.execute(
                    _args(
                        tmp,
                        parsed_containment=str(table),
                        genome_fasta=str(genome),
                        fastq_folder=str(root / "fastq"),
                        output_folder=str(root / "targeted"),
                        threshold=0.5,
                        registry=str(registry_file),
                    )
                )
        assert rc == 1
        assert "skipped SRR1: download truncated (5 of 20 spots); use --allow-truncated" in caplog.text
        assert not mock_run.called

        mock_run.reset_mock()
        with tempfile.TemporaryDirectory() as tmp2:
            root, table, genome = _tree(tmp2)
            registry_file = root / "registry.json"
            registry = load_registry(registry_file)
            registry.datasets["SRR1"] = {
                "download": {"complete": {"verdict": "truncated", "reads_r1": 5, "expected_spots": 20}}
            }
            save_registry(registry)
            rc = cmd.execute(
                _args(
                    tmp2,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    registry=str(registry_file),
                    allow_truncated=True,
                )
            )
        assert rc == 0
        assert mock_run.called

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_debug_keep_sam_leaves_the_sam_on_disk(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    debug_keep_sam=True,
                )
            )
            assert rc == 0
            assert list((root / "targeted" / "SRR1").glob("*.sam"))

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_records_mapped_total(self, mock_run):
        mock_run.side_effect = _fake_tools({"mapped_total": 100, "mapped": 80})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    registry=str(registry_file),
                )
            )
            assert rc == 0
            data = json.loads(registry_file.read_text())
        assert data["datasets"]["SRR1"]["extractions"]["GCF_1"]["mapped_total"] == 100

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_records_reference_coverage(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    registry=str(registry_file),
                )
            )
            assert rc == 0
            data = json.loads(registry_file.read_text())
            entry = data["datasets"]["SRR1"]["extractions"]["GCF_1"]
            assert entry["breadth"] == 0.875 and entry["mean_depth"] == 8.0
            assert entry["coverage_tsv"] == "targeted/SRR1/GCF_1_coverage.tsv"
            assert (root / entry["coverage_tsv"]).exists()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_assembly_preset_recorded_and_passed_to_megahit(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    assemble=True,
                    assembly_preset="meta-large",
                    registry=str(registry_file),
                )
            )
            assert rc == 0
            data = json.loads(registry_file.read_text())
        megahit_call = next(c for c in mock_run.call_args_list if c.args[0] == "megahit" and "-o" in c.args[1])
        args = megahit_call.args[1]
        assert args[args.index("--presets") + 1] == "meta-large"
        assembly = data["datasets"]["SRR1"]["extractions"]["GCF_1"]["assembly"]
        assert assembly["params"]["preset"] == "meta-large"

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_intermediate_contigs_removed_by_default(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    assemble=True,
                )
            )
            assert rc == 0
            asm_dir = root / "targeted" / "SRR1" / "GCF_1_assembly"
            assert not (asm_dir / "intermediate_contigs").exists()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_keep_intermediate_flag_keeps_the_folder(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    assemble=True,
                    keep_intermediate=True,
                )
            )
            assert rc == 0
            asm_dir = root / "targeted" / "SRR1" / "GCF_1_assembly"
            assert (asm_dir / "intermediate_contigs").exists()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_execute_records_coverage_stats(self, mock_run):
        mock_run.side_effect = _fake_tools({"coverage_mapped": 7})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    assemble=True,
                    registry=str(registry_file),
                )
            )
            assert rc == 0
            data = json.loads(registry_file.read_text())
        assembly = data["datasets"]["SRR1"]["extractions"]["GCF_1"]["assembly"]
        assert assembly["reads_mapped"] == 7
        assert assembly["mapping_rate"] is not None
        assert "mean_depth_estimate" in assembly
        assert "genome_fraction_estimate" in assembly
        assert "n90" in assembly and "gc" in assembly and "contigs_ge_1kb" in assembly

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_no_coverage_skips_the_mapping_calls(self, mock_run):
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            rc = cmd.execute(
                _args(
                    tmp,
                    parsed_containment=str(table),
                    genome_fasta=str(genome),
                    fastq_folder=str(root / "fastq"),
                    output_folder=str(root / "targeted"),
                    threshold=0.5,
                    assemble=True,
                    no_coverage=True,
                    registry=str(registry_file),
                )
            )
            assert rc == 0
            data = json.loads(registry_file.read_text())
        coverage_calls = [
            c
            for c in mock_run.call_args_list
            if c.args[0] == "minimap2" and any(str(a).endswith("final.contigs.fa") for a in c.args[1])
        ]
        assert coverage_calls == []
        assembly = data["datasets"]["SRR1"]["extractions"]["GCF_1"]["assembly"]
        assert "reads_mapped" not in assembly

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_default_megahit_scratch_is_unique_per_run(self, mock_run):
        """Two runs sharing one output folder each get their own default megahit scratch
        directory, so a concurrent run cannot remove another run's still-in-use scratch."""
        state = {}
        mock_run.side_effect = _fake_tools(state)
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                registry=str(registry_file),
                assemble=True,
                force=True,
                no_coverage=True,
            )
            reads = [root / "fastq" / "SRR1" / "SRR1_1.fastq.gz", root / "fastq" / "SRR1" / "SRR1_2.fastq.gz"]
            with_reads = {"SRR1": reads}
            results = {"SRR1": ExtractionResult(reads, 10, False, mapped_total=10)}

            assemble_samples(cmd, args, with_reads, results)
            assemble_samples(cmd, args, with_reads, results)

            tmp_dirs = [
                Path(call_args[call_args.index("--tmp-dir") + 1])
                for exe, call_args in state["calls"]
                if exe == "megahit" and "--tmp-dir" in call_args
            ]
            expected_prefix = str(Path(args.output_folder) / ".megahit-tmp-")
        assert len(tmp_dirs) == 2
        assert tmp_dirs[0] != tmp_dirs[1]
        for tmp_dir in tmp_dirs:
            assert str(tmp_dir).startswith(expected_prefix)
            assert not tmp_dir.exists()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_dangling_links_are_reported_once(self, mock_run, caplog):
        """A fastq/<ACC> symlink whose target is missing (e.g. an unmounted store) is
        reported once, by name, before extraction runs."""
        mock_run.side_effect = _fake_tools({})
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = root / "parsed_containment.txt"
            table.write_text("\tGCF_1\n")
            fastq = root / "fastq"
            fastq.mkdir()
            (fastq / "SRR9").symlink_to(root / "nowhere" / "SRR9")
            genome = root / "GCF_1.fna"
            genome.write_text(">s\nACGT\n")
            registry_file = root / "registry.json"
            args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(fastq),
                output_folder=str(root / "targeted"),
                dry_run=True,
                registry=str(registry_file),
            )
            with caplog.at_level(logging.WARNING):
                rc = ExtractTargetReadsCommand().execute(args)
        assert rc == 0
        warnings = [r.getMessage() for r in caplog.records if r.levelname == "WARNING"]
        dangling_warnings = [w for w in warnings if "dangling" in w and "SRR9" in w]
        assert len(dangling_warnings) == 1

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_forced_assembly_failure_clears_the_stale_record(self, mock_run):
        """A forced --assemble redo clears the recorded assembly before megahit reruns. If
        megahit then fails, the previous folder stays on disk (the new one is published by
        rename only on success) and the registry no longer describes it."""
        state = {}
        mock_run.side_effect = _fake_tools(state)
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                assemble=True,
                registry=str(registry_file),
            )
            assert cmd.execute(args) == 0
            before = json.loads(registry_file.read_text())["datasets"]["SRR1"]["extractions"]["GCF_1"]["assembly"]
            assert before is not None and before["contigs"] > 0

            fake = _fake_tools(state)

            def raise_on_megahit_assembly_run(executable, run_args, **kwargs):
                if executable == "megahit" and "-o" in run_args:
                    raise SecurityError("Parameter '--tmp-dir' not allowed for megahit")
                return fake(executable, run_args, **kwargs)

            mock_run.side_effect = raise_on_megahit_assembly_run
            forced_args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                assemble=True,
                force=True,
                registry=str(registry_file),
            )
            rc = cmd.execute(forced_args)
            after = json.loads(registry_file.read_text())["datasets"]["SRR1"]["extractions"]["GCF_1"]["assembly"]
        assert rc == 1
        assert after is None

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_forced_assembly_failure_clears_a_record_whose_folder_is_gone(self, mock_run):
        """The assembly folder was removed by hand after a successful run, so a forced redo has
        nothing to remove on disk; the recorded assembly must still be cleared, so a failed
        rerun never leaves a record describing a missing folder."""
        import shutil

        state = {}
        mock_run.side_effect = _fake_tools(state)
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            common = dict(
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                assemble=True,
                registry=str(registry_file),
            )
            assert cmd.execute(_args(tmp, **common)) == 0
            before = json.loads(registry_file.read_text())["datasets"]["SRR1"]["extractions"]["GCF_1"]["assembly"]
            assert before is not None and before["contigs"] > 0
            shutil.rmtree(root / "targeted" / "SRR1" / "GCF_1_assembly")

            fake = _fake_tools(state)

            def raise_on_megahit_assembly_run(executable, run_args, **kwargs):
                if executable == "megahit" and "-o" in run_args:
                    raise SecurityError("megahit failed")
                return fake(executable, run_args, **kwargs)

            mock_run.side_effect = raise_on_megahit_assembly_run
            rc = cmd.execute(_args(tmp, force=True, **common))
            after = json.loads(registry_file.read_text())["datasets"]["SRR1"]["extractions"]["GCF_1"]["assembly"]
        assert rc == 1
        assert after is None

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_forced_extraction_without_assemble_drops_the_assembly_record(self, mock_run):
        """A --force redo that does not also pass --assemble leaves the assembly folder on
        disk untouched, but its registry record is dropped: it was built from the replaced reads."""
        state = {}
        mock_run.side_effect = _fake_tools(state)
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                assemble=True,
                registry=str(registry_file),
            )
            assert cmd.execute(args) == 0
            before = json.loads(registry_file.read_text())["datasets"]["SRR1"]["extractions"]["GCF_1"]["assembly"]
            assert before is not None and before["contigs"] > 0

            forced_no_assemble = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                assemble=False,
                force=True,
                registry=str(registry_file),
            )
            assert cmd.execute(forced_no_assemble) == 0
            after = json.loads(registry_file.read_text())["datasets"]["SRR1"]["extractions"]["GCF_1"]["assembly"]
            assert (root / "targeted" / "SRR1" / "GCF_1_assembly" / "final.contigs.fa").is_file()
        assert after is None

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_successful_forced_assembly_ends_with_a_fresh_record(self, mock_run):
        """A forced --assemble redo that succeeds must not leave the record cleared; it
        ends with the new assembly block, not None."""
        state = {}
        mock_run.side_effect = _fake_tools(state)
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                assemble=True,
                registry=str(registry_file),
            )
            assert cmd.execute(args) == 0

            forced_args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                assemble=True,
                force=True,
                registry=str(registry_file),
            )
            assert cmd.execute(forced_args) == 0
            after = json.loads(registry_file.read_text())["datasets"]["SRR1"]["extractions"]["GCF_1"]["assembly"]
        assert after is not None
        assert after["contigs"] > 0

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_signal_between_samples_stops_before_the_next_one(self, mock_run):
        """A stop noticed after SRR1 is recorded ends the run before SRR2 is ever mapped."""
        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _two_sample_tree(tmp)
            args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
            )
            seen = []
            original = cmd._record_result

            def spy(args_, accession, outcome, *rest):
                original(args_, accession, outcome, *rest)
                seen.append(accession)
                if len(seen) == 1:
                    args_._termination.stop.set()

            with patch.object(cmd, "_record_result", side_effect=spy):
                rc = cmd.run(args)

            assert rc == 130
            assert seen == ["SRR1"]
            registry = load_registry(args.registry)
            assert "SRR1" in registry.datasets
            assert "SRR2" not in registry.datasets
            assert not (root / "targeted" / "SRR2").exists()

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_concurrent_runs_skip_the_loser_and_record_only_the_winner(self, mock_run):
        """Two commands extracting the same sample at once (overlapping SLURM shards): the
        one that loses the per-sample lock is skipped, not recorded, and the registry ends up
        with exactly the winner's extraction."""
        state = {}
        base_run = _fake_tools(state)
        holding = threading.Event()
        release = threading.Event()

        def run(executable, args, **kwargs):
            if executable == "minimap2" and "-o" in args:
                holding.set()
                release.wait(timeout=5)
            return base_run(executable, args, **kwargs)

        mock_run.side_effect = run

        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _tree(tmp)
            registry_file = root / "registry.json"
            kwargs = dict(
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                registry=str(registry_file),
            )
            return_codes = {}

            def worker(name):
                return_codes[name] = ExtractTargetReadsCommand().run(_args(tmp, **kwargs))

            first = threading.Thread(target=worker, args=("first",))
            first.start()
            assert holding.wait(timeout=5)

            second = threading.Thread(target=worker, args=("second",))
            second.start()
            second.join(timeout=5)
            assert not second.is_alive()  # the per-sample lock is non-blocking

            release.set()
            first.join(timeout=5)

            assert sorted(return_codes.values()) == [0, 1]
            data = json.loads(registry_file.read_text())

        extraction = data["datasets"]["SRR1"]["extractions"]["GCF_1"]
        assert extraction["mapped_reads"] > 0
        map_calls = [c for c in state["calls"] if c[0] == "minimap2" and "-o" in c[1]]
        assert len(map_calls) == 1

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_signal_between_assembly_samples_stops_before_the_next_one(self, mock_run):
        """A stop noticed after SRR1's assembly is recorded ends the run before SRR2's
        megahit ever starts."""
        import prepare.cli.commands.extraction_assembly as cli_module

        mock_run.side_effect = _fake_tools({})
        cmd = ExtractTargetReadsCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _two_sample_tree(tmp)
            args = _args(
                tmp,
                parsed_containment=str(table),
                genome_fasta=str(genome),
                fastq_folder=str(root / "fastq"),
                output_folder=str(root / "targeted"),
                threshold=0.5,
                assemble=True,
            )
            seen = []
            real_record_assembly = cli_module.record_assembly

            def spy(*a, **kw):
                result = real_record_assembly(*a, **kw)
                seen.append(a[1])
                if len(seen) == 1:
                    args._termination.stop.set()
                return result

            with patch.object(cli_module, "record_assembly", side_effect=spy):
                rc = cmd.run(args)

            assert rc == 130
            assert seen == ["SRR1"]
            registry = load_registry(args.registry)
            assert registry.datasets["SRR1"]["extractions"]["GCF_1"]["assembly"] is not None
            # Both samples were already mapped (the mapping loop ran to completion); only the
            # assembly loop was stopped, so SRR2's extraction record carries no assembly.
            assert registry.datasets["SRR2"]["extractions"]["GCF_1"]["assembly"] is None
            assert not (root / "targeted" / "SRR2" / "GCF_1_assembly").exists()


def test_assemble_loads_the_registry_at_most_once_per_sample(tmp_path):
    """A sample whose assembly is kept costs one registry load (for its extraction date) and no write."""
    import prepare.cli.commands.extraction_assembly as cli_mod
    from prepare.data.registry import Registry, record_assembly

    registry_file = tmp_path / "registry.json"
    registry = Registry(path=registry_file)
    accessions = [f"SRR{i}" for i in range(1, 6)]
    for accession in accessions:
        record_extraction(registry, accession, "GCF_1", [], 10, False, {})
        record_assembly(registry, accession, "GCF_1", tmp_path / accession, {"contigs": 1}, "1.2.9", {})
    save_registry(registry, registry_file)
    for accession in accessions:
        # A marker on disk: the kept folder needs no marker written and no inputs recorded.
        marked = tmp_path / "targeted" / accession / "GCF_1_assembly"
        marked.mkdir(parents=True)
        (marked / ".prepare-assembly.json").write_text("{}")
    args = _args(tmp_path, registry=str(registry_file), assemble=True, output_folder=str(tmp_path / "targeted"))
    with_reads = {accession: [tmp_path / f"{accession}.fq"] for accession in accessions}
    results = {accession: ExtractionResult(with_reads[accession], 10, False) for accession in accessions}

    with (
        patch.object(cli_mod, "assemble_extracted_reads", return_value=(None, False)),
        patch.object(cli_mod, "assembly_state", return_value=("current", "marker matches")),
        patch.object(cli_mod, "megahit_version", return_value="1.2.9"),
        patch.object(cli_mod, "fasta_length", return_value=100),
        patch.object(cli_mod, "registry_update") as writes,
        patch.object(cli_mod, "load_registry", wraps=cli_mod.load_registry) as spy,
    ):
        outcome = assemble_samples(ExtractTargetReadsCommand(), args, with_reads, results)
    assert outcome.reused == accessions
    # One load per sample, for its extraction date; a kept folder adds none.
    assert spy.call_count == len(accessions)
    writes.assert_not_called()


@pytest.mark.parametrize("every", [1, 0])
def test_extraction_reports_progress_and_moves_per_sample_lines_to_the_item_level(tmp_path, caplog, monkeypatch, every):
    """With summaries on, per-sample lines (the command's own and the data module's) are DEBUG."""
    monkeypatch.setenv("PREPARE_PROGRESS_EVERY", str(every))
    root, table, genome = _two_sample_tree(tmp_path)
    args = _args(
        tmp_path,
        parsed_containment=str(table),
        genome_fasta=str(genome),
        fastq_folder=str(root / "fastq"),
        output_folder=str(root / "targeted"),
        threshold=0.5,
    )
    without_versions = functools.partial(tools.require_tools, check_versions=False)
    with (
        patch("prepare.utils.tools.shutil.which", return_value="/usr/bin/tool"),
        patch("prepare.cli.commands.read_extraction.require_tools", without_versions),
        patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=_fake_tools({})),
        caplog.at_level(logging.DEBUG),
    ):
        assert ExtractTargetReadsCommand().execute(args) == 0

    summaries = [r for r in caplog.records if r.getMessage().startswith("extract_target_reads: ")]
    assert all(r.levelno == logging.INFO for r in summaries)
    per_sample = [r for r in caplog.records if r.getMessage().startswith(("SRR1: ", "SRR2: ", "Extracted 10 mapped"))]
    assert per_sample
    if every:
        assert [r.getMessage().split(" (")[0] for r in summaries] == [
            "extract_target_reads: 1/2 done",
            "extract_target_reads: finished 2/2",
        ]
        assert all(r.levelno == logging.DEBUG for r in per_sample)
    else:
        assert [r.getMessage().split(" (")[0] for r in summaries] == ["extract_target_reads: finished 2/2"]
        assert all(r.levelno == logging.INFO for r in per_sample)
    # The filter is removed after the run.
    assert not logging.getLogger("prepare.data.read_extraction").filters


def test_sample_line_templates_exist_in_the_data_module():
    """Each template the command moves to the item level is still logged by the data module."""
    import prepare.data.read_extraction as data_mod
    from prepare.cli.commands.read_extraction import SAMPLE_LINE_TEMPLATES

    source = Path(data_mod.__file__).read_text()
    for template in SAMPLE_LINE_TEMPLATES:
        assert f'"{template}"' in source, template


@pytest.fixture
def isolated_store(tmp_path, monkeypatch):
    """A tmp_path store, with the environment kept away from the developer's own store and config."""
    from prepare.store.layout import init_store

    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    monkeypatch.delenv("PREPARE_DATA", raising=False)
    return init_store(tmp_path / "store")


def _store_linked_sample(root, store, sidecar_fields, verdict="unverified"):
    """Move SRR1's reads into ``store``, link them back into fastq/ and write the store's sidecar.

    The registry records the download with ``verdict``, as a project whose record was written
    before the store copy was found short would.
    """
    from prepare.store.layout import sidecar_path, sra_dir
    from prepare.store.sidecar import Sidecar, write_sidecar

    project_dir = root / "fastq" / "SRR1"
    store_dir = sra_dir(store, "SRR1")
    store_dir.mkdir(parents=True)
    for path in sorted(project_dir.iterdir()):
        path.rename(store_dir / path.name)
    project_dir.rmdir()
    project_dir.symlink_to(store_dir)
    write_sidecar(sidecar_path(store, "SRR1"), Sidecar(accession="SRR1", layout="PAIRED", **sidecar_fields))
    registry_file = root / "registry.json"
    registry = load_registry(registry_file)
    complete = {"verdict": verdict, "reads_r1": 5}
    record_download(
        registry, "SRR1", "downloaded", root / "fastq", complete=complete, source="store", store_name="SRR1"
    )
    save_registry(registry)
    return registry_file


def _run_extraction(root, table, genome, registry_file, store, **overrides):
    """Run the command over the one-sample tree, returning the exit code and whether a tool ran."""
    args = _args(
        root,
        parsed_containment=str(table),
        genome_fasta=str(genome),
        fastq_folder=str(root / "fastq"),
        output_folder=str(root / "targeted"),
        threshold=0.5,
        registry=str(registry_file),
        data_root=str(store.root) if store is not None else None,
        **overrides,
    )
    with patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=_fake_tools({})) as run:
        rc = ExtractTargetReadsCommand().execute(args)
    return rc, run.called


_PARTIAL = {"state": "partial", "reads_per_mate": 5, "ncbi": {"spots": 20}}


def test_partial_store_copy_with_an_unverified_registry_verdict_is_skipped(tmp_path, isolated_store, caplog):
    """The registry says unverified, but the store's sidecar records a partial copy: skipped."""
    root, table, genome = _gzip_tree(tmp_path)
    registry_file = _store_linked_sample(root, isolated_store, _PARTIAL)
    with caplog.at_level(logging.WARNING):
        rc, ran = _run_extraction(root, table, genome, registry_file, isolated_store)
    assert rc == 1
    assert not ran
    assert "skipped SRR1: store copy partial (5 of 20 spots); use --allow-truncated" in caplog.text


def test_allow_truncated_extracts_a_partial_store_copy(tmp_path, isolated_store):
    """--allow-truncated covers a partial store copy as it covers a truncated download."""
    root, table, genome = _gzip_tree(tmp_path)
    registry_file = _store_linked_sample(root, isolated_store, _PARTIAL)
    rc, ran = _run_extraction(root, table, genome, registry_file, isolated_store, allow_truncated=True)
    assert rc == 0
    assert ran
    assert (root / "targeted" / "SRR1" / "GCF_1_1.fastq.gz").exists()


def test_failed_store_copy_is_skipped_with_its_verification_error(tmp_path, isolated_store, caplog):
    """A sidecar left in the failed state names its error in the skip line."""
    root, table, genome = _gzip_tree(tmp_path)
    failed = {"state": "failed", "error": "SRR1_2.fastq.gz: truncated gzip stream"}
    registry_file = _store_linked_sample(root, isolated_store, failed)
    with caplog.at_level(logging.WARNING):
        rc, ran = _run_extraction(root, table, genome, registry_file, isolated_store)
    assert rc == 1
    assert not ran
    expected = "skipped SRR1: store copy failed verification: SRR1_2.fastq.gz: truncated gzip stream"
    assert f"{expected}; use --allow-truncated" in caplog.text


def test_unverified_download_in_a_plain_project_is_not_skipped(tmp_path, isolated_store):
    """Without a store sidecar nothing proves the copy short without reading it, so it is extracted."""
    root, table, genome = _gzip_tree(tmp_path)
    registry_file = root / "registry.json"
    registry = load_registry(registry_file)
    record_download(registry, "SRR1", "downloaded", root / "fastq", complete={"verdict": "unverified", "reads_r1": 3})
    save_registry(registry)
    rc, ran = _run_extraction(root, table, genome, registry_file, None)
    assert rc == 0
    assert ran


def test_partial_copy_mode_folder_is_skipped_through_its_own_sidecar(tmp_path, isolated_store, caplog):
    """A copy-mode link is a real folder holding ``<ACC>.json``; that sidecar is read in place."""
    from prepare.store.sidecar import Sidecar, write_sidecar

    root, table, genome = _gzip_tree(tmp_path)
    write_sidecar(root / "fastq" / "SRR1" / "SRR1.json", Sidecar(accession="SRR1", layout="PAIRED", **_PARTIAL))
    registry_file = root / "registry.json"
    registry = load_registry(registry_file)
    record_download(registry, "SRR1", "downloaded", root / "fastq", complete={"verdict": "unverified", "reads_r1": 5})
    save_registry(registry)
    with caplog.at_level(logging.WARNING):
        rc, ran = _run_extraction(root, table, genome, registry_file, None)
    assert rc == 1
    assert not ran
    assert "skipped SRR1: store copy partial (5 of 20 spots); use --allow-truncated" in caplog.text


def test_a_malformed_sidecar_of_an_unselected_sample_does_not_stop_the_run(tmp_path, isolated_store):
    """Only the selected samples' sidecars are read: SRR2 (containment 0.05) is not selected."""
    root, table, genome = _gzip_tree(tmp_path)
    unselected = root / "fastq" / "SRR2"
    unselected.mkdir()
    (unselected / "SRR2.json").write_text("[1, 2]\n")
    rc, ran = _run_extraction(root, table, genome, root / "registry.json", None)
    assert rc == 0
    assert ran


@patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
def test_extraction_from_a_subfolder_writes_under_the_project_targeted(mock_run, tmp_path, monkeypatch):
    """Relative --fastq-folder and --output-folder resolve beside the registry (audit H4 follow-up)."""
    mock_run.side_effect = _fake_tools({})
    (tmp_path / "proj").mkdir()
    root, table, genome = _tree(tmp_path / "proj")
    save_registry(load_registry(root / "prepare_registry.json"))
    (root / "sub").mkdir()
    monkeypatch.chdir(root / "sub")
    args = _args(root, parsed_containment=str(table), genome_fasta=str(genome), threshold=0.5, registry=None)

    assert ExtractTargetReadsCommand().execute(args) == 0

    assert (root / "targeted" / "SRR1").is_dir()
    assert not (root / "sub" / "targeted").exists()
    registry = load_registry(root / "prepare_registry.json")
    files = registry.datasets["SRR1"]["extractions"]["GCF_1"]["files"]
    assert files and all(f.startswith("targeted/SRR1/") for f in files)
