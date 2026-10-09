"""Commands that record into the registry after a long step keep what other writers committed meanwhile.

Each test holds one command at its long step (ranking, parsing, profiling, validation, read
counting) on a ``threading.Event``, records a download for SRR9 from another thread through
``registry_transaction``, then lets the command finish, and checks that the registry file holds
both the command's change and SRR9. Before these commands recorded through ``registry_update``
or a registry batch, they wrote back the snapshot they had loaded before the long step, and SRR9
was lost. The writer must finish before the command is released, which shows the long step runs
without the registry lock; ``blacklist`` is the exception, since it does its whole read-modify-write
under the lock, so there the writer waits for the lock and still lands.
"""

import argparse
import json
import threading
import time
from pathlib import Path
from typing import Any, Callable, Dict

import pytest

from prepare.cli.commands import sra_profile as sra_profile_module
from prepare.cli.commands import sra_report as sra_report_module
from prepare.cli.commands.blacklist import BlacklistCommand
from prepare.cli.commands.metadata import DownloadMetadataCommand, ParseMetadataCommand
from prepare.cli.commands.select import SelectDatasetsCommand
from prepare.cli.commands.sra_enhanced import SRAValidateCommand
from prepare.cli.commands.sra_profile import SRAProfileCommand
from prepare.cli.commands.sra_report import SRAReportCommand
from prepare.cli.commands.status import StatusCommand
from prepare.cli.commands.status import command as status_module
from prepare.data.registry import load_registry, record_download, registry_transaction
from tests.perf_fixtures import write_metadata_folder

WAIT = 20.0
# Upper bound on how long a command may take to reach its gate or finish after it (coverage runs).
CAP = 300.0


class _Gate:
    """Holds the first call of a wrapped function until ``release`` is set."""

    def __init__(self) -> None:
        self.entered = threading.Event()
        self.release = threading.Event()

    def wrap(self, real: Callable[..., Any]) -> Callable[..., Any]:
        def gated(*args: Any, **kwargs: Any) -> Any:
            self.entered.set()
            assert self.release.wait(WAIT), "the test never released the gate"
            return real(*args, **kwargs)

        return gated


def _record_srr9(registry_file: Path, errors: list) -> None:
    try:
        with registry_transaction(registry_file) as registry:
            record_download(registry, "SRR9", "failed", registry_file.parent / "fastq", message="concurrent")
    except Exception as e:  # noqa: B902 - reported to the test thread
        errors.append(e)


def _race(gate: _Gate, run: Callable[[], int], registry_file: Path, command_holds_lock: bool = False) -> int:
    """Run ``run`` in a thread, record SRR9 while it is held at the gate, release it, return its code.

    With ``command_holds_lock`` False the writer must finish while the command is still held,
    which proves the command's long step runs without the registry lock. ``blacklist`` holds the
    lock at its gate, so there the writer is given a second and then waits for the lock.
    """
    outcome: Dict[str, Any] = {}

    def command() -> None:
        try:
            outcome["rc"] = run()
        except BaseException as e:  # noqa: B902 - reported to the test thread
            outcome["error"] = e

    worker = threading.Thread(target=command)
    worker.start()
    deadline = time.monotonic() + CAP
    while not gate.entered.wait(0.2):
        if not worker.is_alive() or time.monotonic() > deadline:
            break
    if not gate.entered.is_set():
        worker.join(CAP)
        if "error" in outcome:
            raise outcome["error"]
        pytest.fail(f"the command returned rc={outcome.get('rc')} before its long step")
    errors: list = []
    writer = threading.Thread(target=_record_srr9, args=(registry_file, errors))
    writer.start()
    if command_holds_lock:
        writer.join(timeout=1.0)
    else:
        writer.join(WAIT)
        assert not writer.is_alive(), "the writer waited for a registry lock the long step should not hold"
    gate.release.set()
    worker.join(CAP)
    writer.join(CAP)
    assert not worker.is_alive() and not writer.is_alive()
    assert not errors, errors
    if "error" in outcome:
        raise outcome["error"]
    return outcome["rc"]


def _datasets(registry_file: Path) -> Dict[str, Any]:
    return json.loads(registry_file.read_text())["datasets"]


def _assert_srr9_kept(registry_file: Path) -> None:
    srr9 = _datasets(registry_file)["SRR9"]["download"]
    assert srr9["state"] == "failed" and srr9["message"] == "concurrent"


def _write_fastq(path: Path, reads) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(reads)))


@pytest.fixture
def no_store(tmp_path, monkeypatch):
    """Keep every test away from a real store and the developer's config."""
    monkeypatch.delenv("PREPARE_DATA", raising=False)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    monkeypatch.setenv("HOME", str(tmp_path / "home"))


pytestmark = pytest.mark.usefixtures("no_store")


def test_select_records_its_selection_without_reverting_a_concurrent_download(tmp_path, monkeypatch):
    from prepare.cli.commands import select as select_module

    (tmp_path / "parsed_containment.txt").write_text("\tGCF_A\tmax_containment\nSRR1\t0.9\t0.9\nSRR2\t0.1\t0.1\n")
    registry_file = tmp_path / "prepare_registry.json"
    args = argparse.Namespace(
        parsed_containment=str(tmp_path / "parsed_containment.txt"),
        genome_id=None,
        genome_ids=None,
        require="any",
        threshold=0.5,
        top_n=None,
        metadata_file=None,
        metadata_column=None,
        metadata_value=None,
        output=str(tmp_path / "accessions.txt"),
        registry=str(registry_file),
        skip_excluded=True,
        skip_downloaded=False,
        no_record=False,
        max_run_size=None,
        min_spots=None,
        max_spots=None,
        platform=None,
    )
    gate = _Gate()
    monkeypatch.setattr(select_module, "select_accessions_ranked", gate.wrap(select_module.select_accessions_ranked))

    assert _race(gate, lambda: SelectDatasetsCommand().execute(args), registry_file) == 0

    assert _datasets(registry_file)["SRR1"]["selection"]["selected"] is True
    _assert_srr9_kept(registry_file)


def test_blacklist_runs_its_read_modify_write_under_the_registry_lock(tmp_path, monkeypatch):
    from prepare.cli.commands import blacklist as blacklist_module

    registry_file = tmp_path / "prepare_registry.json"
    args = argparse.Namespace(
        add=["SRR1"],
        remove=None,
        from_file=None,
        reason="contaminated",
        list=False,
        blacklist_file=str(tmp_path / "blacklist.txt"),
        registry=str(registry_file),
    )
    gate = _Gate()
    monkeypatch.setattr(blacklist_module, "read_blacklist_file", gate.wrap(blacklist_module.read_blacklist_file))

    assert _race(gate, lambda: BlacklistCommand().execute(args), registry_file, command_holds_lock=True) == 0

    assert _datasets(registry_file)["SRR1"]["exclusion"]["excluded"] is True
    assert (tmp_path / "blacklist.txt").read_text() == "SRR1  # contaminated\n"
    _assert_srr9_kept(registry_file)


def test_download_metadata_parses_outside_the_lock_and_keeps_a_concurrent_download(tmp_path, monkeypatch):
    from prepare.cli.commands import metadata as metadata_module

    folder = tmp_path / "metadata"
    paths = write_metadata_folder(folder, count=2, per_file=14, pool=30)
    registry_file = tmp_path / "prepare_registry.json"
    args = argparse.Namespace(
        email="a@b.c",
        matches_folder=str(tmp_path / "matches"),
        metadata_folder=str(folder),
        threshold=0.0,
        dry_run=False,
        accessions_file=None,
        api_key=None,
        batch_size=200,
        registry=str(registry_file),
        data_root=None,
    )
    downloaded = {path.name.split("_")[0]: path for path in paths}
    monkeypatch.setattr(metadata_module, "download_metadata", lambda **_kwargs: downloaded)
    gate = _Gate()
    monkeypatch.setattr(metadata_module, "parse_metadata_xml", gate.wrap(metadata_module.parse_metadata_xml))

    assert _race(gate, lambda: DownloadMetadataCommand().execute(args), registry_file) == 0

    datasets = _datasets(registry_file)
    assert all(datasets[acc]["metadata"]["run_total_spots"] for acc in downloaded)
    _assert_srr9_kept(registry_file)


def test_parse_metadata_records_the_table_without_reverting_a_concurrent_download(tmp_path, monkeypatch):
    from prepare.cli.commands import metadata as metadata_module

    folder = tmp_path / "metadata"
    write_metadata_folder(folder, count=3, per_file=14, pool=30)
    registry_file = tmp_path / "prepare_registry.json"
    args = argparse.Namespace(
        metadata_folder=str(folder),
        metadata_table_file=str(tmp_path / "metadata_table.txt"),
        registry=str(registry_file),
    )
    gate = _Gate()
    monkeypatch.setattr(metadata_module, "_metadata_fields", gate.wrap(metadata_module._metadata_fields))

    assert _race(gate, lambda: ParseMetadataCommand().execute(args), registry_file) == 0

    datasets = _datasets(registry_file)
    assert sum(1 for record in datasets.values() if "metadata" in record) == 3
    _assert_srr9_kept(registry_file)


def _assert_lock_released_for_usage(module, monkeypatch, registry_file: Path) -> list:
    """Wrap ``module.record_usage_many`` so it records whether the registry lock was held."""
    calls: list = []
    real = module.record_usage_many

    def checked(store, registry, rows):
        calls.append(Path(f"{registry_file}.lock").exists())
        return real(store, registry, rows)

    monkeypatch.setattr(module, "record_usage_many", checked)
    return calls


def test_sra_profile_records_analyses_in_one_batch_after_a_concurrent_download(tmp_path, monkeypatch):
    for mate in ("1", "2"):
        _write_fastq(tmp_path / "fastq" / "SRR1" / f"SRR1_{mate}.fastq", ["GGGGCCCCAT", "GCGCATATAT"])
    registry_file = tmp_path / "prepare_registry.json"
    parser = argparse.ArgumentParser()
    SRAProfileCommand().configure_parser(parser)
    args = parser.parse_args([])
    args.fastq_folder = str(tmp_path / "fastq")
    args.output_report = str(tmp_path / "sra_statistics.csv")
    args.output_dir = str(tmp_path / "profiles")
    args.registry = str(registry_file)
    gate = _Gate()
    module = sra_profile_module
    monkeypatch.setattr(module, "resolve_command_store", gate.wrap(module.resolve_command_store))
    usage_calls = _assert_lock_released_for_usage(module, monkeypatch, registry_file)

    assert _race(gate, lambda: SRAProfileCommand().execute(args), registry_file) == 0

    assert "profile" in _datasets(registry_file)["SRR1"]["analyses"]
    _assert_srr9_kept(registry_file)
    assert usage_calls == [False]


def test_sra_report_records_analyses_in_one_batch_after_a_concurrent_download(tmp_path, monkeypatch):
    for i, accession in enumerate(("SRR1", "SRR2")):
        _write_fastq(tmp_path / "fastq" / accession / f"{accession}.fastq", ["GGCCATATCG"[: 8 + i], "ATATATGCAA"])
    registry_file = tmp_path / "prepare_registry.json"
    parser = argparse.ArgumentParser()
    SRAReportCommand().configure_parser(parser)
    args = parser.parse_args([])
    args.fastq_folder = str(tmp_path / "fastq")
    args.output_dir = str(tmp_path / "reports")
    args.registry = str(registry_file)
    args.no_open = True
    args.no_report = True
    (tmp_path / "accessions.txt").write_text("SRR1\nSRR2\n")
    args.accessions_file = str(tmp_path / "accessions.txt")
    gate = _Gate()
    module = sra_report_module
    monkeypatch.setattr(module, "resolve_command_store", gate.wrap(module.resolve_command_store))
    usage_calls = _assert_lock_released_for_usage(module, monkeypatch, registry_file)

    assert _race(gate, lambda: SRAReportCommand().execute(args), registry_file) == 0

    datasets = _datasets(registry_file)
    assert "report" in datasets["SRR1"]["analyses"] and "report" in datasets["SRR2"]["analyses"]
    _assert_srr9_kept(registry_file)
    assert usage_calls == [False]


def test_sra_validate_records_results_after_a_concurrent_download(tmp_path, monkeypatch):
    from prepare.cli.commands import sra_enhanced as sra_enhanced_module

    _write_fastq(tmp_path / "fastq" / "SRR1" / "SRR1.fastq", ["ACGT"])
    registry_file = tmp_path / "prepare_registry.json"
    parser = argparse.ArgumentParser()
    SRAValidateCommand().configure_parser(parser)
    args = parser.parse_args([])
    args.fastq_folder = str(tmp_path / "fastq")
    args.registry = str(registry_file)
    gate = _Gate()
    real = SRAValidateCommand._validate_directory
    gated = gate.wrap(real)
    monkeypatch.setattr(SRAValidateCommand, "_validate_directory", lambda self, *a, **k: gated(self, *a, **k))
    usage_calls = _assert_lock_released_for_usage(sra_enhanced_module, monkeypatch, registry_file)

    assert _race(gate, lambda: SRAValidateCommand().execute(args), registry_file) == 0

    assert _datasets(registry_file)["SRR1"]["analyses"]["validate"]["summary"]["passed"] is True
    _assert_srr9_kept(registry_file)
    assert usage_calls == [False]


def _status_args(root: Path, **overrides) -> argparse.Namespace:
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
        list_missing=False,
        json=True,
        next=False,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def _status_tree(root: Path) -> None:
    for acc in ("SRR1", "SRR2"):
        _write_fastq(root / "fastq" / acc / f"{acc}_1.fastq", ["ACGT"])


def test_status_reconcile_scans_unlocked_and_keeps_a_concurrent_download(tmp_path, monkeypatch, capsys):
    _status_tree(tmp_path)
    registry_file = tmp_path / "prepare_registry.json"
    assert StatusCommand().execute(_status_args(tmp_path, init=True)) == 0
    (tmp_path / "fastq" / "SRR2" / "SRR2_1.fastq").unlink()
    gate = _Gate()
    monkeypatch.setattr(status_module, "scan_reconcile", gate.wrap(status_module.scan_reconcile))

    assert _race(gate, lambda: StatusCommand().execute(_status_args(tmp_path, reconcile=True)), registry_file) == 0

    datasets = _datasets(registry_file)
    assert datasets["SRR2"]["download"]["state"] == "missing"
    assert datasets["SRR1"]["download"]["state"] == "downloaded"
    _assert_srr9_kept(registry_file)
    capsys.readouterr()


def test_status_init_refuses_to_overwrite_a_registry_created_meanwhile(tmp_path, monkeypatch, capsys, caplog):
    _status_tree(tmp_path)
    registry_file = tmp_path / "prepare_registry.json"
    gate = _Gate()
    monkeypatch.setattr(status_module, "bootstrap_from_disk", gate.wrap(status_module.bootstrap_from_disk))

    assert _race(gate, lambda: StatusCommand().execute(_status_args(tmp_path, init=True)), registry_file) == 1

    assert set(_datasets(registry_file)) == {"SRR9"}
    _assert_srr9_kept(registry_file)
    assert "created by another process" in caplog.text
    capsys.readouterr()


def test_read_extraction_records_usage_after_releasing_the_registry_lock(tmp_path, monkeypatch):
    from unittest.mock import patch

    import prepare.store.usage as usage_module
    from helpers_extraction import _fake_tools
    from prepare.cli.commands.read_extraction import ExtractTargetReadsCommand
    from prepare.store.layout import init_store

    table = tmp_path / "parsed_containment.txt"
    table.write_text("\tGCF_1\nSRR1\t0.9\n")
    for mate in ("1", "2"):
        (tmp_path / "fastq" / "SRR1").mkdir(parents=True, exist_ok=True)
        (tmp_path / "fastq" / "SRR1" / f"SRR1_{mate}.fastq.gz").write_text("x")
    genome = tmp_path / "GCF_1.fna"
    genome.write_text(">s\nACGT\n")
    store_root = tmp_path / "store"
    init_store(store_root)
    registry_file = tmp_path / "registry.json"
    with registry_transaction(registry_file) as registry:
        registry.project = {"id": "proj1", "name": "demo", "path": str(tmp_path), "created": "now"}
    lock = Path(f"{registry_file}.lock")
    held_at_write: list = []
    real_write = usage_module.catalog_write

    def checked_write(paths):
        held_at_write.append(lock.exists())
        return real_write(paths)

    monkeypatch.setattr(usage_module, "catalog_write", checked_write)
    args = argparse.Namespace(
        parsed_containment=str(table),
        genome_id="GCF_1",
        genome_fasta=str(genome),
        fastq_folder=str(tmp_path / "fastq"),
        output_folder=str(tmp_path / "targeted"),
        threshold=0.5,
        preset="sr",
        threads=4,
        min_mapq=0,
        temp_folder=None,
        allow_truncated=False,
        debug_keep_sam=False,
        assemble=True,
        assembly_threads=None,
        min_contig_len=None,
        assembly_preset="meta-sensitive",
        keep_intermediate=False,
        no_coverage=False,
        dry_run=False,
        force=False,
        registry=str(registry_file),
        data_root=str(store_root),
    )
    # The external tools are faked, and the pre-flight check is told they are present, so the test
    # runs on a machine without minimap2, samtools or megahit.
    with (
        patch("prepare.utils.tools.shutil.which", return_value="/usr/bin/tool"),
        patch("prepare.data.read_extraction.SecureSubprocess.run_secure", side_effect=_fake_tools({})),
    ):
        assert ExtractTargetReadsCommand().execute(args) == 0

    # One write for the extraction and one for the assembly, neither under the registry lock.
    assert held_at_write == [False, False]
    assert load_registry(registry_file).datasets["SRR1"]["extractions"]["GCF_1"]["assembly"]["contigs"] == 2
