"""Tests for the `status` (local inventory) CLI command."""

import argparse
import gzip
import json
import shlex
import tempfile
from pathlib import Path

import pytest

from prepare.cli.commands.status import StatusCommand
from prepare.cli.commands.status.render_text import selection_detail
from prepare.cli.commands.status.suggest import download_next_steps, next_steps, reselect_command
from prepare.data.registry import (
    Registry,
    SCHEMA_VERSION,
    load_registry,
    record_exclusion,
    record_extraction,
    record_genome,
    record_selection,
    registry_transaction,
    save_registry,
    upsert_dataset,
)
from prepare.processing.status_report import inventory_report

# A minimal single-run NCBI metadata XML, modeled on the fixture in
# tests/test_registry_reconcile_reverify.py: {accession} and {spots} are filled in per test.
_XML = """<?xml version="1.0"?>
<EXPERIMENT_PACKAGE_SET>
    <EXPERIMENT_PACKAGE>
        <EXPERIMENT>
            <IDENTIFIERS><PRIMARY_ID>EXP1</PRIMARY_ID></IDENTIFIERS>
            <LIBRARY_DESCRIPTOR>
                <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
                <LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT>
            </LIBRARY_DESCRIPTOR>
            <PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina HiSeq 2500</INSTRUMENT_MODEL></ILLUMINA></PLATFORM>
        </EXPERIMENT>
        <RUN_SET>
            <RUN accession="{accession}" total_spots="{spots}" total_bases="200" size="300">
                <IDENTIFIERS><PRIMARY_ID>{accession}</PRIMARY_ID></IDENTIFIERS>
            </RUN>
        </RUN_SET>
    </EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>"""


def _args(registry, **kwargs):
    base = dict(
        registry=str(registry),
        data_root=None,
        fastq_folder="fastq",
        metadata_folder="metadata",
        genomes_folder="genomes",
        accessions_file=None,
        parsed_containment=None,
        list_missing=False,
        json=False,
        targeted_folder="targeted",
        matches_folder="matches",
        stage=None,
        genome=None,
        init=False,
        reconcile=False,
        export_tsv=None,
        next=False,
    )
    base.update(kwargs)
    return argparse.Namespace(**base)


def _make_tree(tmp):
    """Create a small fixture tree: one accession fully present, one absent."""
    root = Path(tmp)
    (root / "fastq" / "SRR1").mkdir(parents=True)
    (root / "fastq" / "SRR1" / "SRR1.fastq.gz").write_text("@r\nACGT\n+\nIIII\n")
    (root / "metadata").mkdir()
    (root / "metadata" / "SRR1_metadata.xml").write_text("<xml/>")
    (root / "genomes").mkdir()
    (root / "genomes" / "GCF_000006945.2.fna").write_text(">s\nACGT\n")
    return root


class TestStatusCommand:
    def test_command_properties(self):
        cmd = StatusCommand()
        assert cmd.name == "status"
        assert "local" in cmd.help.lower()

    def test_on_disk_inventory(self, capsys):
        cmd = StatusCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root = _make_tree(tmp)
            result = cmd.execute(
                _args(
                    root / "prepare_registry.json",
                    fastq_folder=str(root / "fastq"),
                    metadata_folder=str(root / "metadata"),
                    genomes_folder=str(root / "genomes"),
                )
            )
        assert result == 0
        out = capsys.readouterr().out
        assert "FASTQ accessions on disk : 1" in out
        assert "Metadata XML on disk     : 1" in out
        assert "Genome FASTA on disk     : 1" in out

    def test_reconcile_present_and_missing(self, capsys):
        cmd = StatusCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root = _make_tree(tmp)
            accs = root / "accs.txt"
            accs.write_text("SRR1\nSRR2\n")
            result = cmd.execute(
                _args(
                    root / "prepare_registry.json",
                    fastq_folder=str(root / "fastq"),
                    metadata_folder=str(root / "metadata"),
                    genomes_folder=str(root / "genomes"),
                    accessions_file=str(accs),
                    list_missing=True,
                )
            )
        assert result == 0
        out = capsys.readouterr().out
        assert "FASTQ    : 1 present, 1 missing" in out
        assert "Metadata : 1 present, 1 missing" in out
        assert "SRR2" in out

    def test_json_output(self, capsys):
        cmd = StatusCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root = _make_tree(tmp)
            accs = root / "accs.txt"
            accs.write_text("SRR1\nSRR2\n")
            result = cmd.execute(
                _args(
                    root / "prepare_registry.json",
                    fastq_folder=str(root / "fastq"),
                    metadata_folder=str(root / "metadata"),
                    genomes_folder=str(root / "genomes"),
                    accessions_file=str(accs),
                    json=True,
                )
            )
        assert result == 0
        report = json.loads(capsys.readouterr().out)
        assert report["on_disk"]["fastq_accessions"] == 1
        assert report["wanted"]["total"] == 2
        assert report["wanted"]["fastq_missing"] == ["SRR2"]

    def test_on_disk_inventory_ignores_transient_temp_folder(self, capsys):
        """A <acc>_temp folder holding a partial FASTQ must not inflate the on-disk count."""
        cmd = StatusCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root = _make_tree(tmp)
            (root / "fastq" / "SRR2_temp").mkdir(parents=True)
            (root / "fastq" / "SRR2_temp" / "SRR2_temp_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            result = cmd.execute(
                _args(
                    root / "prepare_registry.json",
                    fastq_folder=str(root / "fastq"),
                    metadata_folder=str(root / "metadata"),
                    genomes_folder=str(root / "genomes"),
                    json=True,
                )
            )
        assert result == 0
        report = json.loads(capsys.readouterr().out)
        assert report["on_disk"]["fastq_accessions"] == 1

    def test_parsed_containment_supplies_wanted_list(self, capsys):
        cmd = StatusCommand()
        with tempfile.TemporaryDirectory() as tmp:
            root = _make_tree(tmp)
            table = root / "parsed_containment.txt"
            table.write_text("accession\tGCF_x\nSRR1\t0.9\nSRR2\t0.4\n")
            result = cmd.execute(
                _args(
                    root / "prepare_registry.json",
                    fastq_folder=str(root / "fastq"),
                    metadata_folder=str(root / "metadata"),
                    genomes_folder=str(root / "genomes"),
                    parsed_containment=str(table),
                    json=True,
                )
            )
        assert result == 0
        report = json.loads(capsys.readouterr().out)
        assert report["wanted"]["total"] == 2

    def test_missing_accessions_file_errors(self, tmp_path):
        cmd = StatusCommand()
        result = cmd.execute(_args(tmp_path / "prepare_registry.json", accessions_file="/nonexistent/accs.txt"))
        assert result == 1


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
        init=False,
        reconcile=False,
        export_tsv=None,
        next=False,
        list_missing=False,
        json=True,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def test_inventory_ignores_appledouble(tmp_path):
    root = _make_tree(tmp_path)
    (root / "metadata" / "._SRR1_metadata.xml").write_bytes(b"\x00\x05")
    (root / "genomes" / "._g.fna").write_bytes(b"\x00\x05")
    # A hidden accession directory with a real, non-empty FASTQ file inside: only the
    # directory's own dotted name should exclude it, not an empty-folder accident.
    hidden_acc = root / "fastq" / "._SRR1"
    hidden_acc.mkdir()
    (hidden_acc / "reads.fastq.gz").write_bytes(b"x" * 10)
    report = inventory_report(_status_args(root), Registry())
    assert report["on_disk"]["metadata_xml"] == 1
    assert report["on_disk"]["genome_fasta"] == 1
    assert report["on_disk"]["fastq_accessions"] == 1


def _project_tree(root):
    for acc in ("SRR1", "SRR2"):
        d = root / "fastq" / acc
        d.mkdir(parents=True)
        (d / f"{acc}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
    (root / "matches").mkdir()
    (root / "matches" / "GCF_1.csv").write_text("acc,containment,cANI\nSRR1,0.9,0.99\nSRR2,0.4,0.9\nSRR3,0.2,0.8\n")
    (root / "accessions.txt").write_text("SRR1\nSRR2\nSRR3\n")


def test_reconcile_from_a_subfolder_uses_the_project_above(tmp_path, monkeypatch, capsys):
    """Folders are taken relative to the registry's project, not the working directory (audit H4)."""
    root = tmp_path / "proj"
    root.mkdir()
    _project_tree(root)
    assert StatusCommand().execute(_status_args(root, init=True)) == 0
    sub = root / "sub"
    sub.mkdir()
    monkeypatch.chdir(sub)
    args = _args("unused", reconcile=True, json=True)
    args.registry = None
    capsys.readouterr()
    assert StatusCommand().execute(args) == 0
    assert json.loads(capsys.readouterr().out)["on_disk"]["fastq_accessions"] == 2
    reg = load_registry(root / "prepare_registry.json")
    assert reg.datasets["SRR1"]["download"]["state"] == "downloaded"
    assert not (sub / "fastq").exists()


class TestStatusWithRegistry:
    def test_without_registry_bootstraps_in_memory_and_hints(self, tmp_path, capsys):
        _project_tree(tmp_path)
        rc = StatusCommand().execute(_status_args(tmp_path, accessions_file=str(tmp_path / "accessions.txt")))
        out = json.loads(capsys.readouterr().out)
        assert rc == 0
        assert out["registry"]["exists"] is False
        assert out["stages"]["screened"]["count"] == 3 and out["stages"]["downloaded"]["count"] == 2
        assert out["wanted"]["total"] == 3 and out["on_disk"]["fastq_accessions"] == 2
        assert not (tmp_path / "prepare_registry.json").exists()

    def test_init_persists_bootstrap(self, tmp_path, capsys):
        _project_tree(tmp_path)
        rc = StatusCommand().execute(
            _status_args(tmp_path, init=True, accessions_file=str(tmp_path / "accessions.txt"))
        )
        assert rc == 0 and (tmp_path / "prepare_registry.json").exists()
        data = json.loads((tmp_path / "prepare_registry.json").read_text())
        assert data["datasets"]["SRR1"]["download"]["inferred"] is True
        # Bootstrapped screening records the match CSV relative to the project root.
        assert data["datasets"]["SRR1"]["screening"]["genomes"]["GCF_1"]["csv"] == "matches/GCF_1.csv"

    def test_second_init_refuses_to_overwrite(self, tmp_path, capsys):
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True, accessions_file=str(tmp_path / "accessions.txt")))
        capsys.readouterr()
        registry_path = tmp_path / "prepare_registry.json"
        before = registry_path.read_text()
        rc = StatusCommand().execute(_status_args(tmp_path, init=True))
        assert rc == 1
        assert registry_path.read_text() == before

    def test_stage_and_genome_filters_list_accessions(self, tmp_path, capsys):
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        capsys.readouterr()
        rc = StatusCommand().execute(_status_args(tmp_path, stage="screened", genome=["GCF_1"], json=False))
        out = capsys.readouterr().out
        assert rc == 0 and "SRR1" in out and "SRR3" in out

    def test_next_suggests_download_and_extraction(self, tmp_path, capsys):
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True, accessions_file=str(tmp_path / "accessions.txt")))
        capsys.readouterr()
        StatusCommand().execute(_status_args(tmp_path, next=True))
        out = json.loads(capsys.readouterr().out)
        commands = [n["command"] for n in out["next"]]
        assert any(c.startswith("prepare download_sra") for c in commands)  # SRR3 selected, not downloaded
        assert any(
            "extract_target_reads" in c and "GCF_1" in c for c in commands
        )  # SRR1/SRR2 downloaded, not extracted

    def test_next_does_not_suggest_a_no_skip_excluded_list(self, tmp_path, monkeypatch):
        """A selection recorded with --no-skip-excluded may still list an excluded accession,
        so --next must not suggest downloading that selection's output file directly; it
        should instead point at re-running select_datasets with --skip-excluded.

        SRR1 and SRR2 already have FASTQ on disk in `_project_tree`, so `--init` marks them
        downloaded and `download_next_steps` would drop them regardless of selection
        criteria; SRR3 has no FASTQ on disk, so it is the accession left to download and the
        one whose selection criteria this test exercises."""
        root = tmp_path
        _project_tree(root)
        monkeypatch.chdir(root)
        StatusCommand().execute(_status_args(root, init=True))
        with registry_transaction(str(root / "prepare_registry.json")) as reg:
            record_exclusion(reg, "SRR2", "isolate")
            record_selection(
                reg,
                ["SRR1", "SRR2", "SRR3"],
                {"skip_excluded": False, "column": "GCF_A", "threshold": 0.5},
                "sel_noskip.txt",
            )
        steps = download_next_steps(load_registry(str(root / "prepare_registry.json")))
        commands = [s["command"] for s in steps]
        download_commands = [c for c in commands if c.startswith("prepare download_sra")]
        # sel_noskip.txt is never downloaded directly; it is only the --output of the
        # runnable reselect command (checked in test_reselect_suggestion_is_runnable).
        assert not any("sel_noskip.txt" in c for c in download_commands)
        assert any("select_datasets" in c and "--skip-excluded" in c for c in commands)

    def test_mixed_skip_and_no_skip_selections_yield_both_a_download_and_a_reselect(self, tmp_path, monkeypatch):
        """A registry holding both an ordinary (--skip-excluded) selection for one accession
        and a --no-skip-excluded selection for another, each still to download, must surface
        both kinds of next step at once: a direct download command for the ordinary
        selection's output file, and a reselect suggestion for the --no-skip-excluded one.
        Neither must swallow or replace the other.

        SRR4's selection is written directly rather than through a second `record_selection`
        call: that helper unselects anything selected earlier but absent from its own
        `accessions` argument, which would otherwise clear SRR3's selection made just above."""
        root = tmp_path
        _project_tree(root)
        monkeypatch.chdir(root)
        StatusCommand().execute(_status_args(root, init=True))
        with registry_transaction(str(root / "prepare_registry.json")) as reg:
            record_selection(
                reg, ["SRR3"], {"skip_excluded": True, "column": "GCF_A", "threshold": 0.5}, "sel_normal.txt"
            )
            upsert_dataset(reg, "SRR4")["selection"] = {
                "selected": True,
                "date": "2026-01-01T00:00:00+00:00",
                "criteria": {"skip_excluded": False, "column": "GCF_B", "threshold": 0.3},
                "output": "sel_noskip.txt",
            }
        steps = download_next_steps(load_registry(str(root / "prepare_registry.json")))
        commands = [s["command"] for s in steps]

        download_commands = [c for c in commands if c.startswith("prepare download_sra")]
        reselect_commands = [c for c in commands if c.startswith("prepare select_datasets")]

        assert any("sel_normal.txt" in c for c in download_commands)
        assert len(reselect_commands) == 1
        assert "sel_noskip.txt" in reselect_commands[0] and "--skip-excluded" in reselect_commands[0]

    def test_reselect_suggestion_is_runnable(self, tmp_path, monkeypatch):
        """The reselect suggestion for a --no-skip-excluded selection is not a placeholder: it
        is built from the recorded criteria and can be run as-is to redo the selection with
        --skip-excluded, writing back to the same output file."""
        root = tmp_path
        _project_tree(root)
        monkeypatch.chdir(root)
        StatusCommand().execute(_status_args(root, init=True))
        with registry_transaction(str(root / "prepare_registry.json")) as reg:
            record_exclusion(reg, "SRR2", "isolate")
            record_selection(
                reg,
                ["SRR1", "SRR2", "SRR3"],
                {"skip_excluded": False, "column": "GCF_A", "threshold": 0.5},
                "sel_noskip.txt",
            )
        steps = download_next_steps(load_registry(str(root / "prepare_registry.json")))
        commands = [s["command"] for s in steps]
        reselect = next(c for c in commands if c.startswith("prepare select_datasets"))
        assert reselect == (
            "prepare select_datasets --genome-id GCF_A --threshold 0.5 --skip-excluded --output sel_noskip.txt"
        )

    def test_reselect_suggestion_uses_genome_ids_when_recorded(self, tmp_path, monkeypatch):
        """A --no-skip-excluded selection made with --genome-ids reselects the same way."""
        root = tmp_path
        _project_tree(root)
        monkeypatch.chdir(root)
        StatusCommand().execute(_status_args(root, init=True))
        with registry_transaction(str(root / "prepare_registry.json")) as reg:
            record_selection(
                reg,
                ["SRR1", "SRR2", "SRR3"],
                {"skip_excluded": False, "genome_ids": ["GCF_A", "GCF_B"], "require": "all", "threshold": 0.3},
                "sel_noskip.txt",
            )
        steps = download_next_steps(load_registry(str(root / "prepare_registry.json")))
        commands = [s["command"] for s in steps]
        reselect = next(c for c in commands if c.startswith("prepare select_datasets"))
        assert reselect == (
            "prepare select_datasets --genome-ids GCF_A GCF_B --require all "
            "--threshold 0.3 --skip-excluded --output sel_noskip.txt"
        )

    def test_reselect_suggestion_omits_malformed_numbers_and_require(self, tmp_path, monkeypatch):
        """Values read back from a hand-edited registry are coerced (threshold to float, top-N
        to a positive int, require to any/all); a value that does not coerce is omitted rather
        than pasted into the command unquoted."""
        root = tmp_path
        _project_tree(root)
        monkeypatch.chdir(root)
        StatusCommand().execute(_status_args(root, init=True))
        with registry_transaction(str(root / "prepare_registry.json")) as reg:
            record_selection(
                reg,
                ["SRR1", "SRR2", "SRR3"],
                {
                    "skip_excluded": False,
                    "genome_ids": ["GCF_A", "GCF_B"],
                    "require": "all; touch pwned",
                    "threshold": "0.5; touch pwned",
                    "top_n": "3 && touch pwned",
                },
                "sel_noskip.txt",
            )
        steps = download_next_steps(load_registry(str(root / "prepare_registry.json")))
        reselect = next(s["command"] for s in steps if s["command"].startswith("prepare select_datasets"))
        assert "pwned" not in reselect
        assert reselect == ("prepare select_datasets --genome-ids GCF_A GCF_B --skip-excluded --output sel_noskip.txt")

    def test_reselect_command_warns_on_malformed_threshold(self, caplog):
        """A recorded threshold that does not coerce to a float is dropped from the suggested
        command (so select_datasets falls back to its own default on rerun) and now also logs
        a warning naming the bad value, instead of failing silently."""
        with caplog.at_level("WARNING"):
            command = reselect_command({"column": "GCF_A", "threshold": "abc"}, "out.txt")
        assert "threshold" in caplog.text and "not a number" in caplog.text
        assert "'abc'" in caplog.text or '"abc"' in caplog.text
        assert "--threshold" not in command

    def test_reselect_command_warns_on_unknown_require(self, caplog):
        """A recorded require value outside any/all is dropped from the suggested command and
        logs a warning with its own wording (not "is not a number", which does not describe an
        enum value), naming the bad value, instead of failing silently."""
        with caplog.at_level("WARNING"):
            command = reselect_command({"genome_ids": ["GCF_A", "GCF_B"], "require": "maybe"}, "out.txt")
        assert "Recorded require 'maybe' is not 'any' or 'all'; the suggested command omits --require" in caplog.text
        assert "--require" not in command

    def test_reselect_suggestion_includes_metadata_file(self, tmp_path, monkeypatch):
        """A selection recorded with --metadata-file must reproduce that flag on rerun, next
        to the metadata column/value it already reproduces, so the reselect actually rereads
        the same metadata table rather than falling back to metadata table autodetection. The
        path holds a space, so this also exercises that it is shell-quoted like the metadata
        value already is (test_reselect_suggestion_quotes_values_with_spaces)."""
        root = tmp_path
        _project_tree(root)
        monkeypatch.chdir(root)
        StatusCommand().execute(_status_args(root, init=True))
        meta_path = str(root / "my meta.txt")
        with registry_transaction(str(root / "prepare_registry.json")) as reg:
            record_selection(
                reg,
                ["SRR1", "SRR2", "SRR3"],
                {
                    "skip_excluded": False,
                    "column": "GCF_A",
                    "threshold": 0.5,
                    "metadata_file": meta_path,
                    "metadata_column": "country",
                    "metadata_value": "Sweden",
                },
                "sel_noskip.txt",
            )
        steps = download_next_steps(load_registry(str(root / "prepare_registry.json")))
        commands = [s["command"] for s in steps]
        reselect = next(c for c in commands if c.startswith("prepare select_datasets"))
        assert f"--metadata-file {shlex.quote(meta_path)}" in reselect
        assert "--metadata-column country --metadata-value Sweden" in reselect

    def test_reselect_suggestion_coerces_numeric_strings(self):
        command = reselect_command({"column": "GCF_A", "threshold": "0.25", "top_n": "7", "require": "any"}, "out.txt")
        assert command == (
            "prepare select_datasets --genome-id GCF_A --threshold 0.25 --skip-excluded --output out.txt --top-n 7"
        )

    def test_reselect_suggestion_reproduces_run_filters(self):
        """A selection recorded with the run size, spot count and platform filters reselects with
        the same four flags: the size as an integer byte count, the platform shell-quoted."""
        command = reselect_command(
            {
                "column": "GCF_A",
                "threshold": 0.5,
                "max_run_size": 500000000,
                "min_spots": "1000",
                "max_spots": 2000000,
                "platform": "OXFORD NANOPORE",
            },
            "out.txt",
        )
        assert command == (
            "prepare select_datasets --genome-id GCF_A --threshold 0.5 --skip-excluded --output out.txt"
            " --max-run-size 500000000 --min-spots 1000 --max-spots 2000000 --platform 'OXFORD NANOPORE'"
        )
        assert shlex.split(command)[-2:] == ["--platform", "OXFORD NANOPORE"]

    def test_reselect_suggestion_keeps_min_spots_zero(self):
        """--min-spots 0 is a valid flag (it drops runs with no spot count), so it is reproduced."""
        command = reselect_command({"column": "GCF_A", "min_spots": 0}, "out.txt")
        assert command.endswith(" --min-spots 0")

    def test_reselect_command_warns_on_malformed_max_run_size(self, caplog):
        with caplog.at_level("WARNING"):
            command = reselect_command({"column": "GCF_A", "max_run_size": "big; rm -rf ~", "max_spots": -3}, "out.txt")
        assert "Recorded max_run_size 'big; rm -rf ~'" in caplog.text
        assert "Recorded max_spots -3" in caplog.text
        assert "--max-run-size" not in command and "--max-spots" not in command and "rm -rf" not in command

    def test_reselect_suggestion_reproduces_metadata_filter_top_n_and_table(self, tmp_path, monkeypatch):
        """A selection made with a metadata filter, a top-N cap and a non-default containment
        table must reselect the same way: dropping any of those would reproduce a different,
        looser selection or read the wrong file."""
        root = tmp_path
        _project_tree(root)
        monkeypatch.chdir(root)
        StatusCommand().execute(_status_args(root, init=True))
        with registry_transaction(str(root / "prepare_registry.json")) as reg:
            record_selection(
                reg,
                ["SRR1", "SRR2", "SRR3"],
                {
                    "skip_excluded": False,
                    "column": "GCF_A",
                    "threshold": 0.5,
                    "metadata_column": "country",
                    "metadata_value": "Sweden",
                    "top_n": 10,
                    "table": "custom_containment.txt",
                },
                "sel_noskip.txt",
            )
        steps = download_next_steps(load_registry(str(root / "prepare_registry.json")))
        commands = [s["command"] for s in steps]
        reselect = next(c for c in commands if c.startswith("prepare select_datasets"))
        assert reselect == (
            "prepare select_datasets --genome-id GCF_A --threshold 0.5 --skip-excluded "
            "--output sel_noskip.txt --metadata-column country --metadata-value Sweden "
            "--top-n 10 --parsed-containment custom_containment.txt"
        )

    def test_reselect_suggestion_omits_table_flag_when_it_is_the_default(self, tmp_path, monkeypatch):
        """A recorded table equal to select_datasets' own default must not be echoed back as
        --parsed-containment; only a non-default table needs to be named explicitly."""
        root = tmp_path
        _project_tree(root)
        monkeypatch.chdir(root)
        StatusCommand().execute(_status_args(root, init=True))
        with registry_transaction(str(root / "prepare_registry.json")) as reg:
            record_selection(
                reg,
                ["SRR1", "SRR2", "SRR3"],
                {"skip_excluded": False, "column": "GCF_A", "threshold": 0.5, "table": "parsed_containment.txt"},
                "sel_noskip.txt",
            )
        steps = download_next_steps(load_registry(str(root / "prepare_registry.json")))
        commands = [s["command"] for s in steps]
        reselect = next(c for c in commands if c.startswith("prepare select_datasets"))
        assert reselect == (
            "prepare select_datasets --genome-id GCF_A --threshold 0.5 --skip-excluded --output sel_noskip.txt"
        )

    def test_reselect_suggestion_quotes_values_with_spaces(self, tmp_path, monkeypatch):
        """A metadata value containing a space must be shell-quoted, so the suggested command
        is actually safe to paste and run rather than breaking the shell's argument split."""
        root = tmp_path
        _project_tree(root)
        monkeypatch.chdir(root)
        StatusCommand().execute(_status_args(root, init=True))
        with registry_transaction(str(root / "prepare_registry.json")) as reg:
            record_selection(
                reg,
                ["SRR1", "SRR2", "SRR3"],
                {
                    "skip_excluded": False,
                    "column": "GCF_A",
                    "threshold": 0.5,
                    "metadata_column": "city",
                    "metadata_value": "New York",
                },
                "sel_noskip.txt",
            )
        steps = download_next_steps(load_registry(str(root / "prepare_registry.json")))
        commands = [s["command"] for s in steps]
        reselect = next(c for c in commands if c.startswith("prepare select_datasets"))

        assert "--metadata-value 'New York'" in reselect
        assert shlex.split(reselect) == [
            "prepare",
            "select_datasets",
            "--genome-id",
            "GCF_A",
            "--threshold",
            "0.5",
            "--skip-excluded",
            "--output",
            "sel_noskip.txt",
            "--metadata-column",
            "city",
            "--metadata-value",
            "New York",
        ]

    def test_next_extraction_command_is_runnable(self, tmp_path, capsys):
        """The extract suggestion carries the table it was selected from and a FASTA that exists."""
        _project_tree(tmp_path)
        (tmp_path / "genomes").mkdir()
        (tmp_path / "genomes" / "GCF_1.fasta").write_text(">s\nACGT\n")
        table = tmp_path / "tables" / "containment.txt"
        table.parent.mkdir()
        table.write_text("\tGCF_1\nSRR1\t0.9\nSRR2\t0.4\n")
        StatusCommand().execute(_status_args(tmp_path, init=True))
        seeded = load_registry(tmp_path / "prepare_registry.json")
        record_selection(seeded, ["SRR1", "SRR2"], {"column": "GCF_1", "table": str(table)}, tmp_path / "acc.txt")
        save_registry(seeded)
        capsys.readouterr()

        StatusCommand().execute(_status_args(tmp_path, next=True))
        out = json.loads(capsys.readouterr().out)
        extract = next(s for s in out["next"] if "extract_target_reads" in s["command"])
        assert f"--parsed-containment {table}" in extract["command"]
        assert f"--genome-fasta {tmp_path / 'genomes' / 'GCF_1.fasta'}" in extract["command"]

    def test_next_uses_the_recorded_genome_fasta(self, tmp_path, capsys):
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        seeded = load_registry(tmp_path / "prepare_registry.json")
        record_genome(seeded, "GCF_1", tmp_path / "refs" / "wMel.fna", tmp_path / "manifest.csv")
        save_registry(seeded)
        capsys.readouterr()

        StatusCommand().execute(_status_args(tmp_path, next=True))
        out = json.loads(capsys.readouterr().out)
        extract = next(s for s in out["next"] if "extract_target_reads" in s["command"])
        assert f"--genome-fasta {tmp_path / 'refs' / 'wMel.fna'}" in extract["command"]

    def test_a_moved_project_still_reports_the_recorded_genome_fasta(self, tmp_path, capsys):
        """A registry recorded under the old project location still resolves its genome FASTA
        after the project directory is renamed, since the path is stored relative to the
        registry file rather than to the working directory at record time."""
        old_root = tmp_path / "a"
        _project_tree(old_root)
        StatusCommand().execute(_status_args(old_root, init=True))
        seeded = load_registry(old_root / "prepare_registry.json")
        record_genome(seeded, "GCF_1", old_root / "refs" / "wMel.fna", old_root / "manifest.csv")
        save_registry(seeded)
        capsys.readouterr()

        new_root = tmp_path / "b"
        old_root.rename(new_root)

        StatusCommand().execute(_status_args(new_root, next=True))
        out = json.loads(capsys.readouterr().out)
        extract = next(s for s in out["next"] if "extract_target_reads" in s["command"])
        assert f"--genome-fasta {new_root / 'refs' / 'wMel.fna'}" in extract["command"]

    def test_next_drops_excluded_and_already_extracted_accessions(self, tmp_path, capsys):
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True, accessions_file=str(tmp_path / "accessions.txt")))
        seeded = load_registry(tmp_path / "prepare_registry.json")
        record_exclusion(seeded, "SRR3", "16S amplicon")  # selected, not downloaded
        record_exclusion(seeded, "SRR2", "16S amplicon")  # downloaded
        record_extraction(seeded, "SRR1", "GCF_1", [], 0, False, {})  # zero mapped, but recorded
        save_registry(seeded)
        capsys.readouterr()

        StatusCommand().execute(_status_args(tmp_path, next=True))
        out = json.loads(capsys.readouterr().out)
        download = [s for s in out["next"] if "download_sra" in s["command"]]
        extract = [
            s for s in out["next"] if "extract_target_reads" in s["command"] and "--assemble" not in s["command"]
        ]
        assert download == []  # SRR3 was the only one left to download
        assert extract == []  # SRR1 has a record, SRR2 is excluded

    def test_reconcile_marks_missing_and_untracked(self, tmp_path, capsys):
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        capsys.readouterr()
        (tmp_path / "fastq" / "SRR2" / "SRR2_1.fastq").unlink()
        (tmp_path / "fastq" / "SRR8").mkdir()
        (tmp_path / "fastq" / "SRR8" / "SRR8_1.fastq").write_text("@r\nA\n+\nI\n")
        StatusCommand().execute(_status_args(tmp_path, reconcile=True))
        out = json.loads(capsys.readouterr().out)
        assert out["drift"]["recorded_missing"] == ["SRR2"] and out["drift"]["untracked_fastq"] == ["SRR8"]
        data = json.loads((tmp_path / "prepare_registry.json").read_text())
        assert data["datasets"]["SRR2"]["download"]["state"] == "missing"

    def test_reconcile_registers_untracked_work(self, tmp_path, capsys):
        """Untracked FASTQ and extractions are recorded as inferred, and still reported."""
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        capsys.readouterr()
        (tmp_path / "fastq" / "SRR8").mkdir()
        (tmp_path / "fastq" / "SRR8" / "SRR8_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        extracted = tmp_path / "targeted" / "SRR1"
        extracted.mkdir(parents=True)
        for name in ("GCF_1_1.fastq.gz", "GCF_1_2.fastq.gz"):
            with gzip.open(extracted / name, "wt") as handle:
                handle.write("@r1\nACGT\n+\nIIII\n")

        StatusCommand().execute(_status_args(tmp_path, reconcile=True))
        out = json.loads(capsys.readouterr().out)
        assert out["drift"]["untracked_fastq"] == ["SRR8"]
        assert out["drift"]["untracked_extractions"] == [["SRR1", "GCF_1"]]

        datasets = json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]
        assert datasets["SRR8"]["download"]["state"] == "downloaded"
        assert datasets["SRR8"]["download"]["inferred"] is True
        assert datasets["SRR8"]["download"]["attempts"] == 0
        extraction = datasets["SRR1"]["extractions"]["GCF_1"]
        assert extraction["mapped_reads"] > 0 and extraction["inferred"] is True

    def test_reconcile_without_a_registry_refuses(self, tmp_path, caplog):
        _project_tree(tmp_path)
        with caplog.at_level("ERROR"):
            rc = StatusCommand().execute(_status_args(tmp_path, reconcile=True))
        assert rc == 1
        assert not (tmp_path / "prepare_registry.json").exists()
        assert "status --init" in caplog.text

    def test_report_carries_the_schema_version(self, tmp_path, capsys):
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        out = json.loads(capsys.readouterr().out)
        assert out["registry"]["version"] == SCHEMA_VERSION

    def test_text_report_shows_selection_criteria_and_exclusion_reasons(self, tmp_path, capsys):
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        seeded = load_registry(tmp_path / "prepare_registry.json")
        record_selection(
            seeded,
            ["SRR1"],
            {"column": "GCF_1", "threshold": 0.5, "metadata_column": "organism", "metadata_value": "soil"},
            tmp_path / "accessions.txt",
        )
        record_exclusion(seeded, "SRR2", "16S amplicon")
        record_exclusion(seeded, "SRR3", "16S amplicon")
        save_registry(seeded)
        capsys.readouterr()

        StatusCommand().execute(_status_args(tmp_path, json=False))
        out = capsys.readouterr().out
        assert "column GCF_1" in out and "threshold 0.5" in out and "organism = soil" in out
        assert "16S amplicon: 2" in out

    def test_text_report_shows_the_selection_run_filters(self, tmp_path, capsys):
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        seeded = load_registry(tmp_path / "prepare_registry.json")
        record_selection(
            seeded,
            ["SRR1"],
            {
                "column": "GCF_1",
                "threshold": 0.5,
                "max_run_size": 600_000_000,
                "min_spots": 3_000_000,
                "max_spots": 9_000_000,
                "platform": "illumina",
            },
            tmp_path / "accessions.txt",
        )
        save_registry(seeded)
        capsys.readouterr()

        StatusCommand().execute(_status_args(tmp_path, json=False))
        out = capsys.readouterr().out
        assert (
            "column GCF_1, threshold 0.5, max size 600 MB, spots >= 3000000, spots <= 9000000, platform illumina" in out
        )

    @pytest.mark.parametrize(
        "size,text",
        [(2_000_000_000, "max size 2 GB"), (1_500_000, "max size 1.5 MB"), (512, "max size 512 bytes")],
    )
    def test_selection_detail_formats_the_run_size(self, tmp_path, size, text):
        registry = load_registry(tmp_path / "prepare_registry.json")
        record_selection(registry, ["SRR1"], {"column": "GCF_1", "max_run_size": size}, tmp_path / "a.txt")
        assert text in selection_detail(registry)

    def test_selection_detail_leaves_out_a_malformed_run_filter(self, tmp_path):
        registry = load_registry(tmp_path / "prepare_registry.json")
        record_selection(
            registry, ["SRR1"], {"column": "GCF_1", "max_run_size": "lots", "min_spots": -3}, tmp_path / "a.txt"
        )
        detail = selection_detail(registry)
        assert "max size" not in detail and "spots" not in detail

    def test_export_tsv(self, tmp_path, capsys, caplog):
        _project_tree(tmp_path)
        extracted = tmp_path / "targeted" / "SRR1"
        extracted.mkdir(parents=True)
        for name in ("GCF_1_1.fastq.gz", "GCF_1_2.fastq.gz"):
            with gzip.open(extracted / name, "wt") as handle:
                handle.write("@r1\nACGT\n+\nIIII\n")

        with caplog.at_level("INFO"):
            StatusCommand().execute(_status_args(tmp_path, init=True, export_tsv=str(tmp_path / "registry")))

        assert (tmp_path / "registry_datasets.tsv").exists() and (tmp_path / "registry_extractions.tsv").exists()
        ext = (tmp_path / "registry_extractions.tsv").read_text().splitlines()[0]
        assert ext.startswith("accession\t")
        assert ext.split("\t")[-4:] == ["coverage_tsv", "n90", "largest", "assembly_dir"]
        assert "registry_datasets.tsv" in caplog.text and "registry_extractions.tsv" in caplog.text

    def test_text_report_shows_stage_matrix(self, tmp_path, capsys):
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        capsys.readouterr()
        StatusCommand().execute(_status_args(tmp_path, json=False))
        out = capsys.readouterr().out
        assert "Local inventory" in out and "Stages" in out and "screened" in out and "GCF_1" in out

    def test_downloads_report_lists_truncated_and_unverified(self, tmp_path, capsys):
        """report['downloads'] surfaces registry verdicts; stages.downloaded is unaffected."""
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        seeded = load_registry(tmp_path / "prepare_registry.json")
        seeded.datasets["SRR1"]["download"]["complete"] = {"verdict": "truncated"}
        seeded.datasets["SRR2"]["download"]["complete"] = {"verdict": "unverified"}
        save_registry(seeded)
        capsys.readouterr()

        rc = StatusCommand().execute(_status_args(tmp_path))
        out = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert out["downloads"]["truncated"] == ["SRR1"]
        assert out["downloads"]["unverified"] == ["SRR2"]
        assert out["stages"]["downloaded"]["count"] == 2

    def test_text_report_shows_truncated_downloads_line(self, tmp_path, capsys):
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        seeded = load_registry(tmp_path / "prepare_registry.json")
        seeded.datasets["SRR1"]["download"]["complete"] = {"verdict": "truncated"}
        save_registry(seeded)
        capsys.readouterr()

        StatusCommand().execute(_status_args(tmp_path, json=False))
        out = capsys.readouterr().out
        assert "truncated downloads: 1 (SRR1)" in out


class TestStatusStorePlumbing:
    def test_no_store_configured_omits_store_section(self, tmp_path, capsys):
        _project_tree(tmp_path)
        rc = StatusCommand().execute(_status_args(tmp_path))
        out = json.loads(capsys.readouterr().out)
        assert rc == 0
        assert "store" not in out

    def test_no_store_configured_text_report_has_no_store_block(self, tmp_path, capsys):
        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, json=False))
        out = capsys.readouterr().out
        assert "Store" not in out

    def test_data_root_flag_adds_store_section(self, tmp_path, capsys):
        from prepare.store.catalog import catalog_write
        from prepare.store.layout import init_store

        _project_tree(tmp_path)
        store_root = tmp_path / "store"
        paths = init_store(store_root)
        with catalog_write(paths) as cat:
            cat.upsert_project("proj1", "P1", "/p1", "reg1")

        rc = StatusCommand().execute(_status_args(tmp_path, data_root=str(store_root)))
        out = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert out["store"]["root"] == str(store_root.resolve())
        assert out["store"]["datasets"] == {}

    def test_data_root_flag_adds_text_store_block_after_local_inventory(self, tmp_path, capsys):
        from prepare.store.layout import init_store

        _project_tree(tmp_path)
        store_root = tmp_path / "store"
        init_store(store_root)

        StatusCommand().execute(_status_args(tmp_path, data_root=str(store_root), json=False))
        out = capsys.readouterr().out

        assert "Local inventory" in out
        assert "Store" in out
        assert out.index("Local inventory") < out.index("Store")

    def test_registry_recorded_store_root_is_used_when_no_flag(self, tmp_path, capsys):
        from prepare.store.layout import init_store

        _project_tree(tmp_path)
        store_root = tmp_path / "store"
        init_store(store_root)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        seeded = load_registry(tmp_path / "prepare_registry.json")
        seeded.store = {"root": str(store_root), "mode": "symlink", "linked": []}
        save_registry(seeded)
        capsys.readouterr()

        rc = StatusCommand().execute(_status_args(tmp_path))
        out = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert out["store"]["root"] == str(store_root.resolve())


class TestStatusIncompleteStoreLinks:
    def test_status_names_links_to_incomplete_store_datasets(self, tmp_path, capsys):
        """A fastq/<ACC> symlink into a store dataset whose sidecar state is not ready (here
        'failed') is missing, but status must say why rather than reporting it as an ordinary
        missing accession."""
        from prepare.store.layout import init_store, sidecar_path, sra_dir
        from prepare.store.link import link_dataset
        from prepare.store.sidecar import Sidecar, write_sidecar

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        dataset_dir = sra_dir(paths, "SRR1")
        dataset_dir.mkdir(parents=True)
        (dataset_dir / "SRR1.fastq.gz").write_bytes(b"x" * 10)
        write_sidecar(sidecar_path(paths, "SRR1"), Sidecar(accession="SRR1", state="failed"))

        fastq_dir = tmp_path / "fastq"
        fastq_dir.mkdir()
        link_dataset(fastq_dir, "SRR1", paths, mode="absolute")

        (tmp_path / "matches").mkdir()
        (tmp_path / "matches" / "GCF_1.csv").write_text("acc,containment,cANI\nSRR1,0.9,0.99\n")
        (tmp_path / "accessions.txt").write_text("SRR1\n")

        rc = StatusCommand().execute(
            _status_args(
                tmp_path,
                accessions_file=str(tmp_path / "accessions.txt"),
                data_root=str(store_root),
                json=False,
            )
        )
        out = capsys.readouterr().out

        assert rc == 0
        assert "linked to a store dataset that is not complete" in out

    def test_json_report_lists_incomplete_store_links(self, tmp_path, capsys):
        from prepare.store.layout import init_store, sidecar_path, sra_dir
        from prepare.store.link import link_dataset
        from prepare.store.sidecar import Sidecar, write_sidecar

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        dataset_dir = sra_dir(paths, "SRR1")
        dataset_dir.mkdir(parents=True)
        (dataset_dir / "SRR1.fastq.gz").write_bytes(b"x" * 10)
        write_sidecar(sidecar_path(paths, "SRR1"), Sidecar(accession="SRR1", state="partial"))

        fastq_dir = tmp_path / "fastq"
        fastq_dir.mkdir()
        link_dataset(fastq_dir, "SRR1", paths, mode="absolute")

        (tmp_path / "matches").mkdir()
        (tmp_path / "matches" / "GCF_1.csv").write_text("acc,containment,cANI\nSRR1,0.9,0.99\n")
        (tmp_path / "accessions.txt").write_text("SRR1\n")

        rc = StatusCommand().execute(
            _status_args(
                tmp_path,
                accessions_file=str(tmp_path / "accessions.txt"),
                data_root=str(store_root),
            )
        )
        out = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert out["wanted"]["fastq_incomplete_store_links"] == ["SRR1"]
        assert "SRR1" in out["wanted"]["fastq_missing"]

    def test_a_complete_store_link_is_not_reported_as_incomplete(self, tmp_path, capsys):
        from prepare.store.layout import init_store, sidecar_path, sra_dir
        from prepare.store.link import link_dataset
        from prepare.store.sidecar import Sidecar, write_sidecar

        store_root = tmp_path / "store"
        paths = init_store(store_root)
        dataset_dir = sra_dir(paths, "SRR1")
        dataset_dir.mkdir(parents=True)
        (dataset_dir / "SRR1.fastq.gz").write_bytes(b"x" * 10)
        write_sidecar(sidecar_path(paths, "SRR1"), Sidecar(accession="SRR1", state="complete"))

        fastq_dir = tmp_path / "fastq"
        fastq_dir.mkdir()
        link_dataset(fastq_dir, "SRR1", paths, mode="absolute")

        (tmp_path / "matches").mkdir()
        (tmp_path / "matches" / "GCF_1.csv").write_text("acc,containment,cANI\nSRR1,0.9,0.99\n")
        (tmp_path / "accessions.txt").write_text("SRR1\n")

        rc = StatusCommand().execute(
            _status_args(
                tmp_path,
                accessions_file=str(tmp_path / "accessions.txt"),
                data_root=str(store_root),
                json=False,
            )
        )
        out = capsys.readouterr().out

        assert rc == 0
        assert "linked to a store dataset that is not complete" not in out


class TestStatusReconcileVerdict:
    def test_reconcile_fills_missing_download_verdict(self, tmp_path, capsys):
        """A download recorded before completeness verification existed (no `download.complete`)
        gets a verdict computed from what is on disk, once the registry records NCBI's spot count."""
        from prepare.data.registry import record_metadata

        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        capsys.readouterr()

        registry = load_registry(tmp_path / "prepare_registry.json")
        record_metadata(registry, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {"run_total_spots": 100})
        save_registry(registry)

        StatusCommand().execute(_status_args(tmp_path, reconcile=True))
        capsys.readouterr()

        data = json.loads((tmp_path / "prepare_registry.json").read_text())
        complete = data["datasets"]["SRR1"]["download"]["complete"]
        assert complete["verdict"] == "truncated"

    def test_reconcile_backfilling_a_verdict_does_not_disturb_the_recorded_date(self, tmp_path, capsys):
        """Filling in a missing verdict uses set_download_verdict, not record_download, so the
        download's recorded date (and files/bytes_total/message) survive untouched."""
        from prepare.data.registry import record_metadata

        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        capsys.readouterr()

        registry = load_registry(tmp_path / "prepare_registry.json")
        record_metadata(registry, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {"run_total_spots": 100})
        original_date = registry.datasets["SRR1"]["download"]["date"]
        save_registry(registry)

        StatusCommand().execute(_status_args(tmp_path, reconcile=True))
        capsys.readouterr()

        data = json.loads((tmp_path / "prepare_registry.json").read_text())
        assert data["datasets"]["SRR1"]["download"]["date"] == original_date
        assert data["datasets"]["SRR1"]["download"]["complete"]["verdict"] == "truncated"

    def test_reconcile_leaves_existing_verdict_alone(self, tmp_path, capsys):
        """A verdict already on file (e.g. from a fresh download run) is never recomputed."""
        from prepare.data.registry import record_download, record_metadata

        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        capsys.readouterr()

        registry = load_registry(tmp_path / "prepare_registry.json")
        record_metadata(registry, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {"run_total_spots": 100})
        record_download(
            registry,
            "SRR1",
            "downloaded",
            tmp_path / "fastq",
            attempt=False,
            complete={"method": "spots", "ratio": 1.0, "verdict": "complete"},
        )
        save_registry(registry)

        StatusCommand().execute(_status_args(tmp_path, reconcile=True))
        capsys.readouterr()

        data = json.loads((tmp_path / "prepare_registry.json").read_text())
        assert data["datasets"]["SRR1"]["download"]["complete"]["verdict"] == "complete"

    def test_reconcile_skips_store_sourced_downloads(self, tmp_path, capsys):
        """A download whose reads came from the shared store already has its own verdict
        pipeline; reconcile must not overwrite it."""
        from prepare.data.registry import record_download, record_metadata

        _project_tree(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        capsys.readouterr()

        registry = load_registry(tmp_path / "prepare_registry.json")
        record_metadata(registry, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {"run_total_spots": 100})
        record_download(
            registry, "SRR1", "downloaded", tmp_path / "fastq", attempt=False, source="store", store_name="SRR1"
        )
        save_registry(registry)

        StatusCommand().execute(_status_args(tmp_path, reconcile=True))
        capsys.readouterr()

        data = json.loads((tmp_path / "prepare_registry.json").read_text())
        assert "complete" not in data["datasets"]["SRR1"]["download"]


class TestStatusWithoutAReachableStore:
    """An unmounted or moved store must not stop a report about the project."""

    def test_status_reports_the_store_as_unavailable_and_still_lists_dangling_links(self, tmp_path, capsys):
        _project_tree(tmp_path)
        # The link points into a store that is not there: exactly what an unmounted volume
        # leaves behind, and exactly the thing status exists to report.
        gone = tmp_path / "unmounted"
        (tmp_path / "fastq" / "SRR3").symlink_to(gone / "sra" / "SRR3")

        assert StatusCommand().execute(_status_args(tmp_path, init=True)) == 0
        capsys.readouterr()

        rc = StatusCommand().execute(_status_args(tmp_path, data_root=str(gone), reconcile=True))
        out = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert out["store"]["available"] is False
        assert out["store"]["root"] == str(gone.resolve())
        assert "SRR3" in out["drift"]["dangling_links"]

    def test_text_report_says_the_store_could_not_be_read(self, tmp_path, capsys):
        _project_tree(tmp_path)
        gone = tmp_path / "unmounted"

        rc = StatusCommand().execute(_status_args(tmp_path, data_root=str(gone), json=False))
        out = capsys.readouterr().out

        assert rc == 0
        assert "Unavailable" in out

    def test_a_store_root_with_no_catalogue_is_reported_as_unavailable(self, tmp_path, capsys):
        from prepare.store.layout import init_store

        _project_tree(tmp_path)
        store_root = tmp_path / "store"
        init_store(store_root)
        # A store folder that has never been written to has no catalog.sqlite yet.

        rc = StatusCommand().execute(_status_args(tmp_path, data_root=str(store_root)))
        out = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert out["store"]["available"] is False


# ------------------------------------------- one pass over the registry and the folder listings


def _inventory_tree(root):
    """Folders covering each way an accession can read as present or missing on disk."""
    fastq = root / "fastq"
    for acc in ("SRR1", "SRR2", "SRR3", "SRR4", "SRR6"):
        (fastq / acc).mkdir(parents=True)
    (fastq / "SRR1" / "SRR1.fastq.gz").write_bytes(b"x" * 10)
    (fastq / "SRR2" / "SRR2_1.fq").write_bytes(b"x" * 10)
    (fastq / "SRR3" / "SRR3.fastq").write_bytes(b"")  # zero bytes: not present
    # SRR4 is an empty folder; SRR6 holds reads but its sidecar says the download failed.
    (fastq / "SRR6" / "SRR6.fastq").write_bytes(b"x" * 10)
    (fastq / "SRR6" / "SRR6.json").write_text(json.dumps({"state": "failed"}))
    # A transient folder is left out of the listing, but a wanted name equal to it is still probed.
    (fastq / "SRR7_temp").mkdir()
    (fastq / "SRR7_temp" / "SRR7.fastq").write_bytes(b"x" * 10)
    # A hidden folder is left out of the listing too; a wanted hidden name is probed on disk.
    (fastq / ".SRR8").mkdir()
    (fastq / ".SRR8" / "SRR8.fastq").write_bytes(b"x" * 10)
    # A fastq/<ACC> symlink to a folder outside fastq/ (a store link) is listed through the link.
    outside = root / "elsewhere" / "SRR10"
    outside.mkdir(parents=True)
    (outside / "SRR10.fastq.gz").write_bytes(b"x" * 10)
    (fastq / "SRR10").symlink_to(outside, target_is_directory=True)
    (root / "metadata").mkdir()
    for acc in ("SRR1", "SRR4"):
        (root / "metadata" / f"{acc}_metadata.xml").write_text("<xml/>")
    wanted = ["SRR1", "SRR2", "SRR3", "SRR4", "SRR5", "SRR6", "SRR7_temp", "SRR1", ".SRR8", "SRR10"]
    (root / "accessions.txt").write_text("\n".join(wanted) + "\n")


def _reference_inventory(root, wanted):
    """The wanted block as the per-accession probes of the earlier implementation computed it."""
    from prepare.data.sra import accession_has_fastq

    fastq_present = [a for a in wanted if accession_has_fastq(root / "fastq" / a)]
    meta_present = [a for a in wanted if (root / "metadata" / f"{a}_metadata.xml").exists()]
    return {
        "total": len(wanted),
        "fastq_present": len(fastq_present),
        "fastq_missing": [a for a in wanted if a not in fastq_present],
        "fastq_incomplete_store_links": [],
        "metadata_present": len(meta_present),
        "metadata_missing": [a for a in wanted if a not in meta_present],
    }


def test_inventory_wanted_block_matches_per_accession_probes(tmp_path):
    _inventory_tree(tmp_path)
    args = _status_args(tmp_path, accessions_file=str(tmp_path / "accessions.txt"))
    report = inventory_report(args, Registry())
    wanted = ["SRR1", "SRR2", "SRR3", "SRR4", "SRR5", "SRR6", "SRR7_temp", ".SRR8", "SRR10"]
    assert report["wanted"] == _reference_inventory(tmp_path, wanted)
    assert report["wanted"]["fastq_missing"] == ["SRR3", "SRR4", "SRR5", "SRR6"]
    assert report["on_disk"]["fastq_accessions"] == 3


def test_inventory_does_not_probe_the_disk_per_wanted_accession(tmp_path, monkeypatch):
    import pathlib

    from prepare.processing import status_report

    _inventory_tree(tmp_path)
    with open(tmp_path / "accessions.txt", "a") as handle:
        handle.write("\n".join(f"SRR{n}" for n in range(100, 400)) + "\n")
    probed = []
    real_has_fastq = status_report.accession_has_fastq

    def counting_has_fastq(acc_dir):
        probed.append(Path(acc_dir).name)
        return real_has_fastq(acc_dir)

    exists_calls = []
    real_exists = pathlib.Path.exists

    def counting_exists(self, *args, **kwargs):
        exists_calls.append(self)
        return real_exists(self, *args, **kwargs)

    monkeypatch.setattr(status_report, "accession_has_fastq", counting_has_fastq)
    monkeypatch.setattr(pathlib.Path, "exists", counting_exists)

    report = inventory_report(_status_args(tmp_path, accessions_file=str(tmp_path / "accessions.txt")), Registry())

    assert report["wanted"]["total"] == 309
    # Once per folder of the listing, plus the two wanted names the listing leaves out by design.
    assert sorted(probed) == [".SRR8", "SRR1", "SRR10", "SRR2", "SRR3", "SRR4", "SRR6", "SRR7_temp"]
    assert [p.name for p in exists_calls if p.name.endswith("_metadata.xml")] == [".SRR8_metadata.xml"]
    assert len(exists_calls) < 20


def test_wanted_from_parsed_containment_reads_the_index_column(tmp_path):
    table = tmp_path / "parsed.tsv"
    table.write_text("accession\tGCF_A\tmax_containment\nSRR2\t0.5\t0.5\nSRR1\tnot-a-number\t0.1\nSRR2\t0.1\t0.1\n")
    args = _status_args(tmp_path, parsed_containment=str(table))
    assert inventory_report(args, Registry())["wanted"]["fastq_missing"] == ["SRR2", "SRR1"]


def _mixed_registry(tmp_path):
    from prepare.data.registry import record_assembly, record_download, record_screening

    r = Registry(path=tmp_path / "prepare_registry.json")
    for n in range(1, 9):
        record_screening(r, f"SRR{n}", "G1", 0.1 * n, None, "matches", 0.0, None)
    record_screening(r, "SRR2", "G2", 0.3, None, "matches", 0.0, None)
    record_selection(r, ["SRR1", "SRR2", "SRR3", "SRR5"], {"min_containment": 0.1}, tmp_path / "sel.txt")
    record_exclusion(r, "SRR3", "amplicon")
    for acc in ("SRR1", "SRR2", "SRR4"):
        record_download(r, acc, "downloaded", tmp_path / "fastq", attempt=False)
    record_extraction(r, "SRR1", "G1", [], 50, False, {})
    record_extraction(r, "SRR2", "G1", [], 0, False, {})
    record_extraction(r, "SRR2", "G2", [], 9, False, {})
    record_extraction(r, "SRR4", "G3", [], 0, False, {})
    record_assembly(r, "SRR1", "G1", tmp_path / "asm", {"contigs": 4}, "v1", {})
    record_assembly(r, "SRR2", "G2", tmp_path / "asm2", {"contigs": 0}, "v1", {})
    r.datasets["SRR8"]["extractions"] = {"G1": None}
    return r


def test_stage_members_and_genome_counts_match_the_per_stage_queries(tmp_path):
    from prepare.data.registry import STAGES, query, stage_counts, stage_members
    from prepare.processing.status_report import genome_counts

    r = _mixed_registry(tmp_path)
    assert stage_members(r) == {s: query(r, s) for s in STAGES}
    assert genome_counts(r) == stage_counts(r)["genomes"]
    assert set(genome_counts(r)) == {"G1", "G2", "G3"}


def test_build_report_stages_match_the_per_stage_queries(tmp_path):
    from prepare.data.registry import ProjectPaths, STAGES, query
    from prepare.processing.status_report import build_report

    r = _mixed_registry(tmp_path)
    paths = ProjectPaths(fastq=tmp_path / "fastq", targeted=tmp_path / "targeted")
    report = build_report(r, _status_args(tmp_path), paths, tmp_path / "prepare_registry.json", True)
    assert report["stages"] == {s: {"count": len(query(r, s)), "accessions": query(r, s)} for s in STAGES}
    assert report["genomes"]["G1"]["zero_mapped"] == ["SRR2"]
    expected_g3 = {"extracted": 0, "assembled": 0, "zero_mapped": ["SRR4"], "empty_assembly_dirs": []}
    assert report["genomes"]["G3"] == expected_g3


def test_stage_filter_accessions_keeps_first_seen_order_without_duplicates(tmp_path):
    from prepare.processing.status_report import stage_filter_accessions

    r = _mixed_registry(tmp_path)
    screened = ["SRR2", "SRR1", "SRR3", "SRR4", "SRR5", "SRR6", "SRR7", "SRR8"]
    assert stage_filter_accessions(r, "screened", ["G2", "G1"]) == screened
    assert stage_filter_accessions(r, "extracted", ["G2", "G1"]) == ["SRR2", "SRR1"]


def test_download_next_steps_lists_selected_not_excluded_not_downloaded(tmp_path):
    r = _mixed_registry(tmp_path)
    assert [acc for step in download_next_steps(r) for acc in step["accessions"]] == ["SRR5"]


def _timed_project(tmp_path):
    """A project with an initialised registry holding a download, an extraction and an assembly time."""
    from prepare.data.registry import record_assembly
    from prepare.data.registry_timing import set_assembly_timing, set_download_timing, set_extraction_timing

    _project_tree(tmp_path)
    StatusCommand().execute(_status_args(tmp_path, init=True))
    seeded = load_registry(tmp_path / "prepare_registry.json")
    set_download_timing(seeded, "SRR1", "2026-10-01T10:00:00+00:00", 10.0)
    set_download_timing(seeded, "SRR2", "2026-10-01T10:00:10+00:00", 30.0)
    record_extraction(seeded, "SRR1", "GCF_1", [], 5, False, {})
    set_extraction_timing(seeded, "SRR1", "GCF_1", "2026-10-01T10:01:00+00:00", 2.5)
    record_assembly(seeded, "SRR1", "GCF_1", tmp_path / "asm", {"contigs": 1}, "v1", {})
    set_assembly_timing(seeded, "SRR1", "GCF_1", "2026-10-01T10:02:00+00:00", 45.0)
    save_registry(seeded)


def test_build_report_timing_on_a_fixture_registry(tmp_path):
    from prepare.data.registry import ProjectPaths
    from prepare.data.registry_timing import set_download_timing, set_extraction_timing
    from prepare.processing.status_report import build_report

    r = _mixed_registry(tmp_path)
    set_download_timing(r, "SRR1", "2026-10-01T10:00:00+00:00", 4.0)
    set_download_timing(r, "SRR2", "2026-10-01T10:00:00+00:00", 8.0)
    set_extraction_timing(r, "SRR2", "G2", "2026-10-01T10:00:00+00:00", 1.5)
    paths = ProjectPaths(fastq=tmp_path / "fastq", targeted=tmp_path / "targeted")
    report = build_report(r, _status_args(tmp_path), paths, tmp_path / "prepare_registry.json", True)
    assert report["timing"] == {
        "downloads_timed": 2,
        "download_seconds_total": 12.0,
        "download_seconds_median": 6.0,
        "extractions_timed": 1,
        "extraction_seconds_total": 1.5,
        "extraction_seconds_median": 1.5,
        "assemblies_timed": 0,
        "assembly_seconds_total": 0.0,
        "assembly_seconds_median": None,
    }


def test_status_json_and_text_report_timing(tmp_path, capsys):
    _timed_project(tmp_path)
    capsys.readouterr()
    StatusCommand().execute(_status_args(tmp_path))
    timing = json.loads(capsys.readouterr().out)["timing"]
    assert (timing["downloads_timed"], timing["download_seconds_total"], timing["download_seconds_median"]) == (
        2,
        40.0,
        20.0,
    )
    assert (timing["extractions_timed"], timing["assemblies_timed"]) == (1, 1)
    StatusCommand().execute(_status_args(tmp_path, json=False))
    out = capsys.readouterr().out
    assert (
        "  timing     : downloads 2 (40.0 s in total, median 20.0 s), extractions 1 (2.5 s in total, "
        "median 2.5 s), assemblies 1 (45.0 s in total, median 45.0 s)"
    ) in out.splitlines()


def test_text_report_has_no_timing_line_without_timing(tmp_path, capsys):
    _project_tree(tmp_path)
    StatusCommand().execute(_status_args(tmp_path, init=True))
    capsys.readouterr()
    StatusCommand().execute(_status_args(tmp_path, json=False))
    assert "timing" not in capsys.readouterr().out


# -------------------------------------------- Task 24: cross-stage funnel


def test_funnel_line_matches_the_specified_format():
    from prepare.cli.commands.status.render_text import funnel_line

    funnel = {
        "screened": {"accessions": 1200},
        "selected": {"accessions": 300, "excluded": 50},
        "downloaded": {"accessions": 280, "bytes": 1_200_000_000_000, "seconds": 147600.0, "failed": 5},
        "analysed": {"accessions": 200},
        "extracted": {"accessions": 250, "pairs": 260, "seconds": 1000.0},
        "assembled": {"accessions": 90, "pairs": 95, "total_bp": 123456, "seconds": 2000.0},
    }
    assert funnel_line(funnel) == (
        "funnel: 1200 screened, 300 selected, 280 downloaded (1.2 TB, 41 h), 250 extracted, 90 assembled"
    )


def _funnel_with_download(bytes_total, seconds):
    return {
        "screened": {"accessions": 10},
        "selected": {"accessions": 8, "excluded": 0},
        "downloaded": {"accessions": 7, "bytes": bytes_total, "seconds": seconds, "failed": 0},
        "analysed": {"accessions": 0},
        "extracted": {"accessions": 0, "pairs": 0, "seconds": None},
        "assembled": {"accessions": 0, "pairs": 0, "total_bp": 0, "seconds": None},
    }


def test_funnel_line_rounds_the_size_to_one_decimal_and_omits_an_unrecorded_time():
    # The crispatus project's downloads predate 0.7.0, so none carries a time: the line used to
    # read "(4.48731 GB, 0.0 h)".
    from prepare.cli.commands.status.render_text import funnel_line

    line = funnel_line(_funnel_with_download(4_487_310_000, None))
    assert "7 downloaded (4.5 GB), 0 extracted" in line


def test_funnel_line_keeps_a_recorded_time_even_with_a_size_in_terabytes():
    from prepare.cli.commands.status.render_text import funnel_line

    assert "7 downloaded (2.0 TB, 1.5 h)," in funnel_line(_funnel_with_download(2_040_000_000_000, 5400.0))
    assert "7 downloaded (0 bytes, 0.0 h)," in funnel_line(_funnel_with_download(0, 0.0))


def test_build_report_funnel_matches_project_funnel_directly(tmp_path):
    from prepare.data.registry import ProjectPaths
    from prepare.processing.project_funnel import funnel as project_funnel
    from prepare.processing.status_report import build_report

    r = _mixed_registry(tmp_path)
    paths = ProjectPaths(fastq=tmp_path / "fastq", targeted=tmp_path / "targeted")
    report = build_report(r, _status_args(tmp_path), paths, tmp_path / "prepare_registry.json", True)
    assert report["funnel"] == project_funnel(r)


def test_status_json_and_text_report_funnel(tmp_path, capsys):
    _timed_project(tmp_path)
    capsys.readouterr()
    StatusCommand().execute(_status_args(tmp_path))
    report = json.loads(capsys.readouterr().out)
    assert report["funnel"] == {
        "screened": {"accessions": 3},
        "selected": {"accessions": 0, "excluded": 0},
        "downloaded": {"accessions": 2, "bytes": 30, "seconds": 40.0, "failed": 0, "failed_seconds": None},
        "analysed": {"accessions": 0},
        "extracted": {"accessions": 1, "pairs": 1, "seconds": 2.5},
        "assembled": {"accessions": 1, "pairs": 1, "total_bp": 0, "seconds": 45.0},
    }
    StatusCommand().execute(_status_args(tmp_path, json=False))
    lines = capsys.readouterr().out.splitlines()
    assert "funnel: 3 screened, 0 selected, 2 downloaded (30 bytes, 0.0 h), 1 extracted, 1 assembled" in lines
    # The funnel line is additive; the pre-existing timing line is unchanged by this task.
    assert (
        "  timing     : downloads 2 (40.0 s in total, median 20.0 s), extractions 1 (2.5 s in total, "
        "median 2.5 s), assemblies 1 (45.0 s in total, median 45.0 s)"
    ) in lines


def test_export_tsv_carries_the_timing_columns(tmp_path, capsys):
    import pandas as pd

    _timed_project(tmp_path)
    prefix = tmp_path / "registry"
    StatusCommand().execute(_status_args(tmp_path, export_tsv=str(prefix)))
    datasets = pd.read_csv(f"{prefix}_datasets.tsv", sep="\t", index_col=0)
    extractions = pd.read_csv(f"{prefix}_extractions.tsv", sep="\t")
    assert datasets.loc["SRR1", "download_seconds"] == 10.0 and datasets.loc["SRR2", "download_seconds"] == 30.0
    row = extractions.set_index("accession").loc["SRR1"]
    assert (row["extraction_seconds"], row["assembly_seconds"]) == (2.5, 45.0)


def test_export_tsv_extractions_gain_coverage_and_assembly_extra_columns(tmp_path):
    """``to_dataframes``'s extractions frame appends coverage_tsv, n90, largest, assembly_dir,
    all empty for an extraction with no assembly and, for one with an assembly, project-relative."""
    import pandas as pd

    from prepare.data.registry import record_assembly
    from prepare.processing.status_report import to_dataframes

    r = load_registry(tmp_path / "prepare_registry.json")
    coverage = {"breadth": 0.8, "mean_depth": 12.5, "coverage_tsv": tmp_path / "targeted" / "SRR1" / "cov.tsv"}
    record_extraction(r, "SRR1", "GCF_1", [], 50, False, {}, coverage=coverage)
    stats = {"contigs": 5, "total_bp": 40000, "n50": 3000, "n90": 900, "largest": 9000, "gc": 0.4123}
    record_assembly(r, "SRR1", "GCF_1", tmp_path / "targeted" / "SRR1" / "GCF_1_assembly", stats, "v1", {})
    record_extraction(r, "SRR2", "GCF_1", [], 5, False, {})  # no assembly, no coverage
    _, extractions = to_dataframes(r)
    rows = extractions.set_index("accession")
    assembled = rows.loc["SRR1"]
    assert assembled["coverage_tsv"] == "targeted/SRR1/cov.tsv"
    assert (assembled["n90"], assembled["largest"]) == (900, 9000)
    assert assembled["assembly_dir"] == "targeted/SRR1/GCF_1_assembly"
    bare = rows.loc["SRR2"]
    assert bare["coverage_tsv"] is None or pd.isna(bare["coverage_tsv"])
    assert pd.isna(bare["n90"]) and pd.isna(bare["largest"])
    assert bare["assembly_dir"] is None or pd.isna(bare["assembly_dir"])


# -------------------------------------------- --init fills metadata; reconcile's new report fields


def test_init_fills_metadata_from_xml_with_a_spot_count(tmp_path, capsys):
    """--init bootstraps a metadata block with no spot count for every XML on disk; it must then
    be filled from that same XML in the same write, not left for a later --reconcile."""
    _project_tree(tmp_path)
    (tmp_path / "metadata").mkdir()
    (tmp_path / "metadata" / "SRR1_metadata.xml").write_text(_XML.format(accession="SRR1", spots=4))

    rc = StatusCommand().execute(_status_args(tmp_path, init=True, accessions_file=str(tmp_path / "accessions.txt")))

    assert rc == 0
    data = json.loads((tmp_path / "prepare_registry.json").read_text())
    metadata = data["datasets"]["SRR1"]["metadata"]
    assert metadata["run_total_spots"] == 4
    assert metadata["platform"] == "ILLUMINA"
    # Bootstrap marked the block inferred; filling it in keeps that mark.
    assert metadata["inferred"] is True


class TestStatusReportsStoreUnavailable:
    """A fastq/<ACC> link into a store that is not mounted at all (its sra/ folder is also gone)."""

    def _unmounted_link(self, tmp_path):
        _project_tree(tmp_path)
        gone = tmp_path / "unmounted"
        (tmp_path / "fastq" / "SRR3").symlink_to(gone / "sra" / "SRR3")
        return gone

    def test_reconcile_json_reports_store_unavailable_and_changes_no_records(self, tmp_path, capsys):
        gone = self._unmounted_link(tmp_path)
        assert StatusCommand().execute(_status_args(tmp_path, init=True)) == 0
        capsys.readouterr()
        before = json.loads((tmp_path / "prepare_registry.json").read_text())

        rc = StatusCommand().execute(_status_args(tmp_path, data_root=str(gone), reconcile=True, json=True))
        report = json.loads(capsys.readouterr().out)

        assert rc == 0
        assert report["drift"]["store_unavailable"] == ["SRR3"]
        assert report["drift"]["recorded_missing"] == []

        after = json.loads((tmp_path / "prepare_registry.json").read_text())
        assert after["datasets"] == before["datasets"]
        assert after["store"] == before["store"]

    def test_reconcile_text_warns_about_an_unmounted_store(self, tmp_path, capsys):
        gone = self._unmounted_link(tmp_path)
        StatusCommand().execute(_status_args(tmp_path, init=True))
        capsys.readouterr()

        rc = StatusCommand().execute(_status_args(tmp_path, data_root=str(gone), reconcile=True, json=False))
        out = capsys.readouterr().out

        assert rc == 0
        assert "WARNING" in out
        assert "SRR3" in out
        assert "not mounted" in out


def test_reconcile_json_drift_always_carries_the_new_keys(tmp_path, capsys):
    """Even when nothing is found, the JSON document names store_unavailable, metadata_filled,
    verdicts_rechecked and assemblies_dropped, so a reader does not have to special-case their absence."""
    _project_tree(tmp_path)
    StatusCommand().execute(_status_args(tmp_path, init=True))
    capsys.readouterr()

    rc = StatusCommand().execute(_status_args(tmp_path, reconcile=True, json=True))
    report = json.loads(capsys.readouterr().out)

    assert rc == 0
    assert report["drift"]["store_unavailable"] == []
    assert report["drift"]["metadata_filled"] == []
    assert report["drift"]["verdicts_rechecked"] == []
    assert report["drift"]["assemblies_dropped"] == []


def test_text_report_unchanged_when_the_three_new_lists_are_empty(tmp_path, capsys):
    """A reconcile that finds nothing new under the new reporting prints no WARNING line and no
    re-checked/filled line; the existing Drift section text is untouched."""
    _project_tree(tmp_path)
    StatusCommand().execute(_status_args(tmp_path, init=True))
    capsys.readouterr()

    rc = StatusCommand().execute(_status_args(tmp_path, reconcile=True, json=False))
    out = capsys.readouterr().out

    assert rc == 0
    assert "Drift against disk" in out
    assert "WARNING" not in out
    assert "Verdicts re-checked" not in out
    assert "Metadata filled from XML" not in out
    assert "Assemblies older than their extraction" not in out


def _setup_metadata_fill_and_verdict_recheck(tmp_path, capsys):
    """SRR1 downloaded with an unverified verdict and a metadata block with no spot count, then a
    real XML on disk with a spot count: one --reconcile should fill both. verdicts_rechecked only
    counts a verdict that was literally "unverified", not a merely missing one (see
    _fill_verdicts), so the download must be recorded with that verdict explicitly."""
    from prepare.data.registry import record_download, record_metadata

    _project_tree(tmp_path)
    StatusCommand().execute(_status_args(tmp_path, init=True))
    capsys.readouterr()

    registry = load_registry(tmp_path / "prepare_registry.json")
    record_metadata(registry, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {})
    record_download(
        registry,
        "SRR1",
        "downloaded",
        tmp_path / "fastq",
        attempt=False,
        complete={
            "method": "unverified",
            "ratio": None,
            "verdict": "unverified",
            "expected_spots": None,
            "reads_r1": None,
        },
    )
    save_registry(registry)
    (tmp_path / "metadata").mkdir()
    (tmp_path / "metadata" / "SRR1_metadata.xml").write_text(_XML.format(accession="SRR1", spots=1))


def test_reconcile_json_reports_metadata_filled_and_verdicts_rechecked(tmp_path, capsys):
    _setup_metadata_fill_and_verdict_recheck(tmp_path, capsys)

    rc = StatusCommand().execute(_status_args(tmp_path, reconcile=True, json=True))
    report = json.loads(capsys.readouterr().out)

    assert rc == 0
    assert report["drift"]["metadata_filled"] == ["SRR1"]
    assert report["drift"]["verdicts_rechecked"] == ["SRR1"]


def test_reconcile_text_reports_metadata_filled_and_verdicts_rechecked(tmp_path, capsys):
    _setup_metadata_fill_and_verdict_recheck(tmp_path, capsys)

    rc = StatusCommand().execute(_status_args(tmp_path, reconcile=True, json=False))
    out = capsys.readouterr().out

    assert rc == 0
    assert "Metadata filled from XML" in out
    assert "Verdicts re-checked" in out


def _seed_stale_assembly(tmp_path, capsys):
    """A 0.7.0-shaped registry: SRR1's assembly for GCF_1 is dated before its (re-)extraction."""
    _project_tree(tmp_path)
    StatusCommand().execute(_status_args(tmp_path, init=True))
    capsys.readouterr()
    registry = load_registry(tmp_path / "prepare_registry.json")
    registry.datasets["SRR1"]["extractions"] = {
        "GCF_1": {
            "date": "2026-09-02T10:00:00+00:00",
            "files": [],
            "reads_mapped": 10,
            "assembly": {"date": "2026-09-01T10:00:00+00:00", "contigs": 12, "total_bp": 9000, "n50": 900},
        },
        "GCF_2": {
            "date": "2026-09-02T10:00:00+00:00",
            "files": [],
            "reads_mapped": 10,
            "assembly": {"date": "2026-09-03T10:00:00+00:00", "contigs": 3, "total_bp": 3000, "n50": 1000},
        },
    }
    save_registry(registry)


def test_reconcile_drops_an_assembly_recorded_before_its_extraction(tmp_path, capsys, monkeypatch):
    from prepare.cli.commands.results import ResultsTableCommand

    _seed_stale_assembly(tmp_path, capsys)

    rc = StatusCommand().execute(_status_args(tmp_path, reconcile=True, json=True))
    first = json.loads(capsys.readouterr().out)
    assert rc == 0
    assert first["drift"]["assemblies_dropped"] == ["SRR1/GCF_1"]
    registry = load_registry(tmp_path / "prepare_registry.json")
    assert registry.datasets["SRR1"]["extractions"]["GCF_1"].get("assembly") is None
    assert registry.datasets["SRR1"]["extractions"]["GCF_1"]["reads_mapped"] == 10
    assert registry.datasets["SRR1"]["extractions"]["GCF_2"]["assembly"]["contigs"] == 3

    StatusCommand().execute(_status_args(tmp_path, reconcile=True, json=True))
    second = json.loads(capsys.readouterr().out)
    assert second["drift"]["assemblies_dropped"] == []
    assert load_registry(tmp_path / "prepare_registry.json").datasets == registry.datasets

    monkeypatch.chdir(tmp_path)
    args = argparse.Namespace(
        output=str(tmp_path / "results.tsv"),
        genome_id=None,
        parsed_containment=str(tmp_path / "missing.txt"),
        min_containment=0.0,
        registry=str(tmp_path / "prepare_registry.json"),
        no_record=True,
    )
    assert ResultsTableCommand().execute(args) == 0
    lines = (tmp_path / "results.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    rows = {tuple(line.split("\t")[:2]): dict(zip(header, line.split("\t"))) for line in lines[1:]}
    assert rows[("SRR1", "GCF_1")]["contigs"] == ""
    assert rows[("SRR1", "GCF_1")]["n50"] == ""
    assert rows[("SRR1", "GCF_2")]["contigs"] == "3"


def test_reconcile_text_names_the_dropped_assemblies(tmp_path, capsys):
    _seed_stale_assembly(tmp_path, capsys)

    rc = StatusCommand().execute(_status_args(tmp_path, reconcile=True, json=False))
    out = capsys.readouterr().out

    assert rc == 0
    assert "Assemblies older than their extraction  : 1 record(s) dropped: SRR1/GCF_1" in out


class TestNextStepsAreRunnableAndRelevant:
    """P4: extraction is suggested only for screened pairs, and every value is shell-quoted."""

    def _registry(self, tmp_path, genomes):
        registry = load_registry(tmp_path / "prepare_registry.json")
        for genome_id in genomes:
            record_genome(registry, genome_id, tmp_path / "genomes" / f"{genome_id}.fna", tmp_path / "m.json")
        return registry

    def test_extraction_is_suggested_only_for_screened_pairs(self, tmp_path):
        from prepare.data.registry import ProjectPaths, record_download, record_screening

        registry = self._registry(tmp_path, ["G1", "G2"])
        for acc in ("SRR1", "SRR3"):
            record_download(registry, acc, "downloaded", tmp_path / "fastq" / acc)
        record_screening(registry, "SRR1", "G1", 0.9, None, "matches", 0.0, None)
        record_screening(registry, "SRR1", "G2", 0.8, None, "matches", 0.0, None)
        record_screening(registry, "SRR3", "G1", 0.7, None, "matches", 0.0, None)
        steps = next_steps(registry, ProjectPaths.anchored(tmp_path))
        by_genome = {("G1" if "--genome-id G1" in s["command"] else "G2"): s["accessions"] for s in steps}
        assert by_genome == {"G1": ["SRR1", "SRR3"], "G2": ["SRR1"]}

    def test_extraction_is_suggested_only_for_what_the_command_would_extract(self, tmp_path):
        """Review: a suggested accession below the default threshold, excluded, or outside the
        selection would be listed on every --next and never extracted."""
        from prepare.data.registry import ProjectPaths, record_download, record_screening

        registry = self._registry(tmp_path, ["G1"])
        for acc, value in (("SRR1", 0.05), ("SRR2", 0.9), ("SRR3", 0.7), ("SRR4", 0.8)):
            record_download(registry, acc, "downloaded", tmp_path / "fastq" / acc)
            record_screening(registry, acc, "G1", value, None, "matches", 0.0, None)
        record_exclusion(registry, "SRR4", "isolate")
        record_selection(registry, ["SRR2"], {"column": "G1", "threshold": 0.5}, "accessions.txt")
        steps = next_steps(registry, ProjectPaths.anchored(tmp_path))
        extract = [s for s in steps if "extract_target_reads" in s["command"]]
        assert [s["accessions"] for s in extract] == [["SRR2"]]

    def test_an_unscreened_project_still_suggests_every_downloaded_accession(self, tmp_path):
        """A status --init project without screening blocks is not silenced."""
        from prepare.data.registry import ProjectPaths, record_download

        registry = self._registry(tmp_path, ["G1"])
        record_download(registry, "SRR1", "downloaded", tmp_path / "fastq" / "SRR1")
        steps = next_steps(registry, ProjectPaths.anchored(tmp_path))
        assert [s["accessions"] for s in steps] == [["SRR1"]]

    def test_suggestions_quote_paths_and_ids(self, tmp_path):
        from prepare.data.registry import ProjectPaths, record_download

        registry = self._registry(tmp_path, ["L crispatus"])
        record_download(registry, "SRR1", "downloaded", tmp_path / "fastq" / "SRR1")
        record_selection(registry, ["SRR1", "SRR2"], {"table": "my tables/parsed.txt"}, "my lists/acc.txt")
        steps = next_steps(registry, ProjectPaths.anchored(tmp_path))
        commands = [s["command"] for s in steps]
        download = next(c for c in commands if c.startswith("prepare download_sra"))
        assert shlex.split(download)[-1] == "my lists/acc.txt"
        extract = next(c for c in commands if "extract_target_reads" in c)
        tokens = shlex.split(extract)
        assert tokens[tokens.index("--genome-id") + 1] == "L crispatus"
        assert tokens[tokens.index("--parsed-containment") + 1] == "my tables/parsed.txt"
        assert tokens[tokens.index("--genome-fasta") + 1] == str(tmp_path / "genomes" / "L crispatus.fna")


class TestInitRegistersGenomesWithEvidence:
    """P13: a FASTA under genomes/ alone is a file; a genome needs screening or extraction evidence."""

    def _init(self, root):
        (root / "fastq").mkdir(exist_ok=True)
        assert StatusCommand().execute(_status_args(root, init=True)) == 0
        return json.loads((root / "prepare_registry.json").read_text())

    def test_a_fasta_alone_is_not_a_project_genome(self, tmp_path):
        alone = tmp_path / "alone"
        (alone / "genomes").mkdir(parents=True)
        (alone / "genomes" / "X.fna").write_text(">s\nACGT\n")
        assert "X" not in self._init(alone).get("genomes", {})

        screened = tmp_path / "screened"
        (screened / "genomes").mkdir(parents=True)
        (screened / "genomes" / "X.fna").write_text(">s\nACGT\n")
        (screened / "matches").mkdir()
        (screened / "matches" / "X.csv").write_text("acc,containment\nSRR1,0.9\n")
        assert "X" in self._init(screened)["genomes"]

    def test_extraction_folders_still_name_their_genome(self, tmp_path):
        root = tmp_path / "extracted"
        (root / "targeted" / "SRR1").mkdir(parents=True)
        (root / "targeted" / "SRR1" / "Y_assembly").mkdir()
        assert "Y" in self._init(root)["genomes"]
