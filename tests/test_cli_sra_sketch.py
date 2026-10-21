"""Tests for ``sra_sketch``: one sourmash signature per downloaded dataset, with local containment."""

import argparse
import gzip
import json
import logging
import os
import random
import sys
from unittest.mock import patch

import pytest

from prepare.cli.commands.sra_sketch import ANALYSIS_NAME, SRASketchCommand
from prepare.cli.main import main
from prepare.core.exceptions import ConfigurationError
from prepare.data import run_log
from prepare.utils.termination import Termination


def _reads(seed, count, length=1500):
    rng = random.Random(seed)
    return ["".join(rng.choice("ACGT") for _ in range(length)) for _ in range(count)]


def _write_fastq(path, reads, gz=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(reads))
    if gz:
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return path


def _dataset(tmp_path, accession="SRR1", seed=1):
    acc_dir = tmp_path / "fastq" / accession
    _write_fastq(acc_dir / f"{accession}_1.fastq", _reads(seed, 2))
    _write_fastq(acc_dir / f"{accession}_2.fastq.gz", _reads(seed + 10, 2), gz=True)
    return acc_dir


def _args(tmp_path, **overrides):
    parser = argparse.ArgumentParser()
    SRASketchCommand().configure_parser(parser)
    args = parser.parse_args([])
    args.fastq_folder = str(tmp_path / "fastq")
    args.output_dir = str(tmp_path / "sketches")
    args.matches_folder = str(tmp_path / "matches")
    args.registry = str(tmp_path / "prepare_registry.json")
    args._termination = Termination()
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def _registry(tmp_path):
    return json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]


def _genome(tmp_path, seed=1):
    reads = _reads(seed, 2)
    tail = "".join(random.Random(99).choice("ACGT") for _ in range(20000))
    path = tmp_path / "genome.fasta"
    path.write_text(f">g\n{reads[0]}{tail}\n")
    return path


@pytest.fixture(autouse=True)
def _sourmash():
    pytest.importorskip("sourmash")


def test_parser_defaults_and_exclusive_genome_flags():
    parser = argparse.ArgumentParser()
    SRASketchCommand().configure_parser(parser)
    args = parser.parse_args(["--accession", "SRR1", "--accession", "SRR2"])
    assert args.accession == ["SRR1", "SRR2"]
    assert args.fastq_folder == "fastq" and args.output_dir == "sketches" and args.matches_folder == "matches"
    assert args.max_workers is None and args.force is False and args.dry_run is False
    with pytest.raises(SystemExit):
        parser.parse_args(["--genome-fasta", "g.fa", "--signature", "g.sig"])
    with pytest.raises(SystemExit):
        parser.parse_args(["--max-workers", "0"])


def test_sketch_writes_signature_and_records(tmp_path):
    _dataset(tmp_path)

    assert SRASketchCommand().execute(_args(tmp_path)) == 0

    signature = tmp_path / "sketches" / "SRR1.sig.gz"
    assert signature.is_file()
    entry = _registry(tmp_path)["SRR1"]["analyses"][ANALYSIS_NAME]
    assert entry["output"] == "sketches/SRR1.sig.gz"
    assert set(entry["summary"]) == {
        "ksizes",
        "scaled",
        "track_abundance",
        "reads",
        "bases",
        "hashes",
        "files",
        "source",
    }
    assert entry["summary"]["reads"] == 4 and entry["summary"]["source"] == "computed"
    assert set(entry["summary"]["files"]) == {"SRR1_1.fastq", "SRR1_2.fastq.gz"}


def test_run_log_summary_and_rows(tmp_path, monkeypatch):
    _dataset(tmp_path)
    from prepare.core import settings
    from prepare.data.registry import load_registry, save_registry

    monkeypatch.setenv(settings.SETTINGS["run_log"].env, "true")
    settings.reset_for_tests()
    monkeypatch.chdir(tmp_path)
    save_registry(load_registry(tmp_path / "prepare_registry.json"))

    assert main(["sra_sketch", "--accession", "SRR1", "--genome-fasta", str(_genome(tmp_path))]) == 0

    records = run_log.read_runs(tmp_path)
    record = [r for r in records if r.command == "sra_sketch"][-1]
    assert record.summary["accessions"] == 1 and record.summary["sketched"] == 1
    assert record.summary["skipped"] == 0 and record.summary["failed"] == 0
    assert record.summary["genome_id"] == "genome" and record.summary["containment_rows"] == 1
    detail = run_log.read_detail(tmp_path, record)
    rows = detail["analyses"][ANALYSIS_NAME]
    assert set(rows) == {"SRR1"}
    assert all(not isinstance(v, (dict, list)) for v in rows["SRR1"].values())
    assert {"reads", "bases", "hashes_k21", "hashes_k31", "hashes_k51", "source", "containment", "cani"} <= set(
        rows["SRR1"]
    )


def test_second_run_skips_and_force_or_a_changed_file_resketches(tmp_path):
    acc_dir = _dataset(tmp_path)
    assert SRASketchCommand().execute(_args(tmp_path)) == 0
    signature = tmp_path / "sketches" / "SRR1.sig.gz"
    first = signature.stat().st_mtime_ns

    with patch("prepare.cli.commands.sra_sketch.sketch_dataset", side_effect=AssertionError("must skip")):
        assert SRASketchCommand().execute(_args(tmp_path)) == 0
    assert signature.stat().st_mtime_ns == first

    assert SRASketchCommand().execute(_args(tmp_path, force=True)) == 0
    assert signature.stat().st_mtime_ns != first

    second = signature.stat().st_mtime_ns
    os.utime(acc_dir / "SRR1_1.fastq", (1, 1))
    assert SRASketchCommand().execute(_args(tmp_path)) == 0
    assert signature.stat().st_mtime_ns != second


def test_dry_run_writes_nothing(tmp_path, capsys):
    _dataset(tmp_path)
    args = _args(tmp_path, dry_run=True)

    assert SRASketchCommand().execute(args) == 0

    assert "SRR1" in capsys.readouterr().out
    assert not (tmp_path / "sketches").exists()
    assert not (tmp_path / "prepare_registry.json").exists()
    assert SRASketchCommand().records_run(args) is False
    assert SRASketchCommand().records_run(_args(tmp_path)) is True


def test_corrupt_gzip_fails_that_accession_only(tmp_path, caplog):
    _dataset(tmp_path, "SRR1")
    _dataset(tmp_path, "SRR2", seed=5)
    (tmp_path / "fastq" / "SRR2" / "SRR2_2.fastq.gz").write_bytes(b"\x1f\x8b\x08\x00broken")

    with caplog.at_level(logging.WARNING):
        assert SRASketchCommand().execute(_args(tmp_path, max_workers=1)) == 1

    assert (tmp_path / "sketches" / "SRR1.sig.gz").is_file()
    assert not (tmp_path / "sketches" / "SRR2.sig.gz").exists()
    assert "SRR2" in caplog.text
    assert "SRR2" not in _registry(tmp_path)


def test_interrupt_records_finished_and_exits_130(tmp_path):
    from prepare.sra import sketch as sketch_module

    _dataset(tmp_path, "SRR1")
    _dataset(tmp_path, "SRR2", seed=5)
    args = _args(tmp_path, max_workers=1)
    real = sketch_module.sketch_dataset

    def stop_after_first(accession, files, stop=None):
        result = real(accession, files, stop=stop)
        args._termination.stop.set()
        return result

    with patch("prepare.cli.commands.sra_sketch.sketch_dataset", side_effect=stop_after_first):
        code = SRASketchCommand().run(args)

    assert code == 130
    written = sorted(p.name for p in (tmp_path / "sketches").iterdir())
    assert len(written) == 1 and written[0].endswith(".sig.gz")
    assert set(_registry(tmp_path)) == {written[0][: -len(".sig.gz")]}


def test_interrupt_counts_the_stopped_accession_as_neither(tmp_path, caplog):
    """The accession stopped mid-sketch is neither sketched nor failed, and the stop is reported once."""
    from prepare.sra import sketch as sketch_module

    _dataset(tmp_path, "SRR1")
    _dataset(tmp_path, "SRR2", seed=5)
    args = _args(tmp_path, max_workers=1)
    real = sketch_module.sketch_dataset
    summaries = []

    def stop_during_second(accession, files, stop=None):
        if accession == "SRR2":
            # A stop set by a worker, with no KeyboardInterrupt on the main thread: the
            # result is collected inside the loop, the path a signal handler never takes.
            args._termination.stop.set()
            raise InterruptedError("stopped")
        return real(accession, files, stop=stop)

    with (
        patch("prepare.cli.commands.sra_sketch.sketch_dataset", side_effect=stop_during_second),
        patch(
            "prepare.cli.commands.sra_sketch.run_log.note_run", side_effect=lambda a, summary: summaries.append(summary)
        ),
        caplog.at_level(logging.WARNING),
    ):
        code = SRASketchCommand().run(args)

    assert code == 130
    assert summaries[0]["sketched"] == 1 and summaries[0]["failed"] == 0
    assert set(_registry(tmp_path)) == {"SRR1"}
    warnings = [r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING]
    assert len(warnings) == 1 and warnings[0].startswith("Interrupted"), warnings


def test_a_single_accession_stopped_in_the_loop_is_not_a_failure(tmp_path):
    """One accession whose worker sets the stop and raises: exit 130, failed 0, nothing recorded."""
    _dataset(tmp_path, "SRR1")
    args = _args(tmp_path, max_workers=1)
    summaries = []

    def stop_at_once(accession, files, stop=None):
        args._termination.stop.set()
        raise InterruptedError("stopped")

    with (
        patch("prepare.cli.commands.sra_sketch.sketch_dataset", side_effect=stop_at_once),
        patch(
            "prepare.cli.commands.sra_sketch.run_log.note_run", side_effect=lambda a, summary: summaries.append(summary)
        ),
    ):
        code = SRASketchCommand().run(args)

    assert code == 130
    assert summaries[0]["failed"] == 0 and summaries[0]["sketched"] == 0
    assert not (tmp_path / "sketches").exists() or list((tmp_path / "sketches").iterdir()) == []


def test_worker_error_keeps_the_finished_records(tmp_path, caplog):
    """A bug in one worker fails that accession, names it, and keeps the others' records."""
    from prepare.sra import sketch as sketch_module

    _dataset(tmp_path, "SRR1")
    _dataset(tmp_path, "SRR2", seed=5)
    args = _args(tmp_path, max_workers=1)
    real = sketch_module.sketch_dataset
    summaries = []

    def break_on_second(accession, files, stop=None):
        if accession == "SRR2":
            raise RuntimeError("bug")
        return real(accession, files, stop=stop)

    with (
        patch("prepare.cli.commands.sra_sketch.sketch_dataset", side_effect=break_on_second),
        patch(
            "prepare.cli.commands.sra_sketch.run_log.note_run", side_effect=lambda a, summary: summaries.append(summary)
        ),
        caplog.at_level(logging.WARNING),
    ):
        code = SRASketchCommand().run(args)

    assert code == 1
    assert set(_registry(tmp_path)) == {"SRR1"}
    assert summaries[0]["sketched"] == 1 and summaries[0]["failed"] == 1
    assert any("SRR2" in r.getMessage() and "bug" in r.getMessage() for r in caplog.records)


def test_containment_csv_written_and_parsed(tmp_path):
    from prepare.data.branchwater import parse_containment_data
    from prepare.sra.sketch import load_dataset_minhash
    from prepare.sra.sketch_containment import genome_minhash
    from prepare.data.branchwater_search import sketch_fasta

    _dataset(tmp_path, "SRR1", seed=1)
    _dataset(tmp_path, "SRR2", seed=7)
    genome = _genome(tmp_path, seed=1)

    assert SRASketchCommand().execute(_args(tmp_path, genome_fasta=str(genome))) == 0

    csv_path = tmp_path / "matches" / "genome.csv"
    assert csv_path.is_file()
    genome_mh = genome_minhash(sketch_fasta(genome))
    expected = genome_mh.contained_by(load_dataset_minhash(tmp_path / "sketches" / "SRR1.sig.gz"))
    summary = parse_containment_data(tmp_path / "matches", tmp_path / "parsed.txt", tmp_path / "top.txt")
    assert summary.table.loc["SRR1", "genome"] == pytest.approx(expected, abs=1e-4)
    assert summary.table.index[0] == "SRR1"


def test_signature_input_and_genome_id(tmp_path):
    from prepare.data.branchwater_search import sketch_fasta

    _dataset(tmp_path)
    genome = _genome(tmp_path)
    sig = tmp_path / "wMel.sig"
    sig.write_text(json.dumps([sketch_fasta(genome)]))

    assert SRASketchCommand().execute(_args(tmp_path, signature=str(sig), genome_id="wMel-2026")) == 0

    assert (tmp_path / "matches" / "wMel-2026.csv").is_file()


def test_genome_id_of_a_gz_fasta_drops_both_suffixes(tmp_path):
    _dataset(tmp_path)
    genome = _genome(tmp_path)
    packed = tmp_path / "G.fna.gz"
    with gzip.open(packed, "wb") as handle:
        handle.write(genome.read_bytes())

    assert SRASketchCommand().execute(_args(tmp_path, genome_fasta=str(packed))) == 0

    assert (tmp_path / "matches" / "G.csv").is_file()


def test_store_cache_round_trip(tmp_path, monkeypatch):
    from prepare.store.catalog import Catalog
    from prepare.store.layout import init_store, store_paths
    from prepare.store.sidecar import Sidecar, read_sidecar, write_sidecar

    store_root = tmp_path / "store"
    monkeypatch.setenv("PREPARE_DATA", str(store_root))
    paths = init_store(store_root)
    store_acc = paths.sra / "SRR1"
    _write_fastq(store_acc / "SRR1_1.fastq", _reads(1, 2))
    write_sidecar(store_acc / "SRR1.json", Sidecar(accession="SRR1", state="complete"))

    def project(name):
        root = tmp_path / name
        (root / "fastq").mkdir(parents=True)
        (root / "fastq" / "SRR1").symlink_to(store_acc)
        from prepare.data.registry import load_registry, save_registry

        registry = load_registry(root / "prepare_registry.json")
        registry.project = {"id": name, "name": name, "path": str(root), "created": "now"}
        save_registry(registry)
        return root

    first = project("p1")
    assert SRASketchCommand().execute(_args(first, data_root=str(store_root))) == 0
    assert (store_acc / "SRR1.sig.gz").is_file()
    assert read_sidecar(store_acc / "SRR1.json").sketch["file"] == "SRR1.sig.gz"

    second = project("p2")
    with patch("prepare.cli.commands.sra_sketch.sketch_dataset", side_effect=AssertionError("must reuse")):
        assert SRASketchCommand().execute(_args(second, data_root=str(store_root))) == 0
    assert (second / "sketches" / "SRR1.sig.gz").is_file()
    assert _registry(second)["SRR1"]["analyses"][ANALYSIS_NAME]["summary"]["source"] == "store"

    with Catalog(store_paths(store_root)) as catalog:
        catalog.migrate()
        rows = catalog.conn.execute("SELECT stage, detail FROM usage WHERE accession = ?", ("SRR1",)).fetchall()
    assert {(r["stage"], r["detail"]) for r in rows} == {("analysed", ANALYSIS_NAME)}


def test_a_store_cache_write_error_only_costs_the_cache(tmp_path, caplog):
    _dataset(tmp_path)
    with caplog.at_level(logging.WARNING):
        with patch("prepare.cli.commands.sra_sketch.store_sketch", side_effect=PermissionError("read-only store")):
            assert SRASketchCommand().execute(_args(tmp_path)) == 0
    assert (tmp_path / "sketches" / "SRR1.sig.gz").is_file()
    assert ANALYSIS_NAME in _registry(tmp_path)["SRR1"]["analyses"]
    assert "read-only store" in caplog.text


def test_a_failed_accession_gets_no_containment_from_a_stale_signature(tmp_path):
    import csv

    _dataset(tmp_path, "SRR1", seed=1)
    _dataset(tmp_path, "SRR2", seed=7)
    genome = _genome(tmp_path, seed=1)
    assert SRASketchCommand().execute(_args(tmp_path, genome_fasta=str(genome))) == 0
    stale = tmp_path / "sketches" / "SRR2.sig.gz"
    assert stale.is_file()

    (tmp_path / "fastq" / "SRR2" / "SRR2_2.fastq.gz").write_bytes(b"\x1f\x8b\x08\x00broken")
    assert SRASketchCommand().execute(_args(tmp_path, genome_fasta=str(genome), max_workers=1)) == 1

    with open(tmp_path / "matches" / "genome.csv", newline="") as handle:
        assert [r["acc"] for r in csv.DictReader(handle)] == ["SRR1"]


def test_an_accession_without_files_keeps_its_branchwater_row(tmp_path):
    """Only an accession that failed while sketching loses its row; one with no files is unknown, not wrong."""
    import csv

    _dataset(tmp_path, "SRR1", seed=1)
    genome = _genome(tmp_path, seed=1)
    matches = tmp_path / "matches"
    matches.mkdir()
    (matches / "genome.csv").write_text(
        "acc,containment,cANI,biosample\nSRR1,0.5000,0.97,SAMN1\nSRR9,0.4000,0.96,SAMN9\n"
    )

    code = SRASketchCommand().execute(_args(tmp_path, accession=["SRR1", "SRR9"], genome_fasta=str(genome)))

    assert code == 1
    with open(matches / "genome.csv", newline="") as handle:
        rows = {r["acc"]: r for r in csv.DictReader(handle)}
    assert set(rows) == {"SRR1", "SRR9"}
    assert rows["SRR9"]["containment"] == "0.4000" and rows["SRR9"]["biosample"] == "SAMN9"


def test_missing_sourmash_is_exit_3(tmp_path, monkeypatch, caplog):
    _dataset(tmp_path)
    monkeypatch.setitem(sys.modules, "sourmash", None)

    with caplog.at_level(logging.ERROR):
        assert SRASketchCommand().execute(_args(tmp_path)) == 3

    assert "prepare[sourmash]" in caplog.text
    assert not (tmp_path / "sketches").exists()


def test_missing_sourmash_raises_configuration_error_before_work(tmp_path, monkeypatch):
    _dataset(tmp_path)
    monkeypatch.setitem(sys.modules, "sourmash", None)
    with pytest.raises(ConfigurationError):
        SRASketchCommand()._run(_args(tmp_path))


def test_folders_that_are_not_accessions_are_skipped_with_a_note(tmp_path, caplog):
    """Without --accession, a folder under fastq/ that is not an accession is skipped, not failed."""
    _dataset(tmp_path)
    (tmp_path / "fastq" / "notes").mkdir()
    (tmp_path / "fastq" / "notes" / "readme.txt").write_text("not reads\n")
    with caplog.at_level(logging.INFO):
        assert SRASketchCommand().execute(_args(tmp_path)) == 0
    assert (tmp_path / "sketches" / "SRR1.sig.gz").is_file()
    assert set(_registry(tmp_path)) == {"SRR1"}
    assert any("notes" in r.getMessage() and r.levelno == logging.INFO for r in caplog.records)


def test_an_accession_without_fastq_files_fails(tmp_path, caplog):
    _dataset(tmp_path)
    with caplog.at_level(logging.WARNING):
        assert SRASketchCommand().execute(_args(tmp_path, accession=["SRR1", "SRR9"])) == 1
    assert "SRR9" in caplog.text
    assert (tmp_path / "sketches" / "SRR1.sig.gz").is_file()
