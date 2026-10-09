"""Tests for the completeness verdicts download_sra records (prepare/cli/commands/sra_verdicts.py).

An expected spot count falls back to the metadata XML when the registry has none; a store relink
never turns a recorded ``truncated`` verdict into ``unverified``; and an accession found on disk
without a ``downloaded`` record (a run killed after the files were in place) is verified against
its spot count instead of being recorded with no verdict. No real tool, store or network is used:
fake tools come from ``tests/helpers_processes.py`` and every store lives under ``tmp_path``.
"""

import argparse
import gzip
import json
from unittest.mock import patch

import pytest

from prepare.cli.commands import sra_verdicts
from prepare.cli.commands.sra import DownloadSraCommand
from prepare.data import registry_blocks as rb
from prepare.data.registry import Registry, load_registry, save_registry
from prepare.store.layout import init_store

from tests.helpers_processes import FAKE_READS, cli_env, install_fake_tools, run_cli
from tests.test_store_sidecar import NCBI_XML

ACC = "SRR1"
UNVERIFIED = {"verdict": "unverified", "ratio": None, "expected_spots": None, "reads_r1": None}
TRUNCATED = {"verdict": "truncated", "ratio": 0.5, "expected_spots": 10, "reads_r1": 5}


@pytest.fixture(autouse=True)
def isolated_home(tmp_path, monkeypatch):
    """No developer store or configuration is read: PREPARE_DATA unset, HOME under tmp_path."""
    home = tmp_path / "home"
    (home / ".config").mkdir(parents=True)
    monkeypatch.delenv("PREPARE_DATA", raising=False)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(home / ".config"))
    return home


def _write_xml(folder, total_spots, accession=ACC):
    """Write ``<accession>_metadata.xml`` with ``total_spots`` as its run's spot count into ``folder``."""
    folder.mkdir(parents=True, exist_ok=True)
    text = NCBI_XML.replace('total_spots="10"', f'total_spots="{total_spots}"').replace("SRR1", accession)
    (folder / f"{accession}_metadata.xml").write_text(text)


def _write_fastq(acc_dir, reads, accession=ACC):
    """Write a gzipped mate pair with ``reads`` records each into ``acc_dir``."""
    acc_dir.mkdir(parents=True, exist_ok=True)
    for mate in (1, 2):
        text = "".join(f"@{accession}.{i} {i}/{mate}\nACGT\n+\nIIII\n" for i in range(1, reads + 1))
        with gzip.open(acc_dir / f"{accession}_{mate}.fastq.gz", "wt") as handle:
            handle.write(text)


def _parse(tmp_path, *extra):
    """``download_sra`` arguments for a project under ``tmp_path``, parsed by the command's own parser."""
    parser = argparse.ArgumentParser()
    DownloadSraCommand().configure_parser(parser)
    accessions = tmp_path / "accessions.txt"
    if not accessions.exists():
        accessions.write_text(f"{ACC}\n")
    return parser.parse_args(
        [
            "--accessions-file",
            str(accessions),
            "--fastq-folder",
            str(tmp_path / "fastq"),
            "--registry",
            str(tmp_path / "prepare_registry.json"),
            "--min-free-gb",
            "0",
            *extra,
        ]
    )


def _write_registry(tmp_path, datasets):
    """Write a registry file holding ``datasets`` under ``tmp_path``."""
    path = tmp_path / "prepare_registry.json"
    registry = Registry(path=path)
    registry.datasets.update(json.loads(json.dumps(datasets)))
    return save_registry(registry, path)


# ------------------------------------------------------------------ the three behaviour changes


def test_present_accession_without_record_gets_truncated_verdict(tmp_path, monkeypatch):
    """Files left in place by a killed run (no ``downloaded`` record) are verified, not recorded blind."""
    install_fake_tools(tmp_path / "bin", tmp_path / "barrier")
    monkeypatch.setenv("PATH", str(tmp_path / "bin"))
    _write_fastq(tmp_path / "fastq" / ACC, reads=2)
    _write_xml(tmp_path / "metadata", total_spots=10)

    assert DownloadSraCommand().execute(_parse(tmp_path)) == 0

    block = rb.download_block(load_registry(tmp_path / "prepare_registry.json"), ACC)
    assert block.state == "downloaded"
    assert block.complete is not None
    assert block.complete.verdict == "truncated"
    assert block.complete.reads_r1 == 2
    assert block.complete.expected_spots == 10


def test_store_relink_never_turns_truncated_into_unverified(tmp_path, monkeypatch):
    """A link to a store copy whose sidecar holds no counts keeps the recorded ``truncated`` verdict."""
    install_fake_tools(tmp_path / "bin", tmp_path / "barrier")
    monkeypatch.setenv("PATH", str(tmp_path / "bin"))
    store = init_store(tmp_path / "store")
    truncated = {"verdict": "truncated", "ratio": 0.5, "expected_spots": 10, "reads_r1": 5}
    _write_registry(tmp_path, {ACC: {"download": {"attempts": 1, "state": "failed", "complete": truncated}}})
    unverified = {"verdict": "unverified", "ratio": None, "expected_spots": None, "reads_r1": None}

    def fake_download_sra(**kwargs):
        kwargs["on_result"](ACC, True, "linked from store, 2 files")
        return {"total": 1, "successful": 1, "failed": 0, "failed_accessions": [], "results": {}}

    with (
        patch("prepare.cli.commands.sra.download_sra", side_effect=fake_download_sra),
        patch.object(DownloadSraCommand, "_sidecar_completeness", staticmethod(lambda _s, _a: dict(unverified))),
    ):
        assert DownloadSraCommand().execute(_parse(tmp_path, "--data-root", str(store.root))) == 0

    block = rb.download_block(load_registry(tmp_path / "prepare_registry.json"), ACC)
    assert block.source == "store"
    assert block.complete.verdict == "truncated"
    assert block.complete.reads_r1 == 5


def test_xml_only_spots_verify_a_plain_download(tmp_path):
    """With no registry metadata, the project's metadata XML supplies the count a plain download is judged by."""
    install_fake_tools(tmp_path / "bin", tmp_path / "barrier")

    _write_xml(tmp_path / "metadata", total_spots=FAKE_READS)
    (tmp_path / "accessions.txt").write_text(f"{ACC}\n")
    env = cli_env(tmp_path, tmp_path / "bin")
    result = run_cli(
        [
            "download_sra",
            "--accessions-file",
            "accessions.txt",
            "--fastq-folder",
            "fastq",
            "--registry",
            "prepare_registry.json",
            "--min-free-gb",
            "0",
        ],
        cwd=tmp_path,
        env=env,
        timeout=60,
    )
    assert result.returncode == 0, result.stderr

    block = rb.download_block(load_registry(tmp_path / "prepare_registry.json"), ACC)
    assert block.state == "downloaded"
    assert block.complete is not None
    assert block.complete.verdict == "complete"
    assert block.complete.expected_spots == FAKE_READS


def test_present_verdicts_are_computed_outside_the_registry_lock(tmp_path):
    """``verify_download`` reads FASTQ files before the one registry write takes its lock."""
    _write_fastq(tmp_path / "fastq" / ACC, reads=3)
    args = _parse(tmp_path)
    lock = tmp_path / "prepare_registry.json.lock"
    held = []
    real = sra_verdicts.verify_download

    def spy(*a, **kw):
        held.append(lock.exists())
        return real(*a, **kw)

    with patch.object(sra_verdicts, "verify_download", side_effect=spy):
        DownloadSraCommand()._record_run_outcomes(
            args, {"already_downloaded_accessions": [ACC]}, tmp_path / "fastq", None, expected_spots={ACC: 3}
        )
    assert held == [False]
    block = rb.download_block(load_registry(args.registry), ACC)
    assert block.complete.verdict == "complete"
    assert block.complete.reads_r1 == 3


# ------------------------------------------------------------------ registry_inputs


def test_registry_inputs_reads_project_xml_before_store_xml(tmp_path):
    store = init_store(tmp_path / "store")
    _write_xml(tmp_path / "metadata", total_spots=7)
    _write_xml(store.metadata, total_spots=9)
    _write_xml(store.metadata, total_spots=11, accession="SRR2")
    (tmp_path / "accessions.txt").write_text(f"{ACC}\nSRR2\nSRR3\n")
    registry = load_registry(tmp_path / "prepare_registry.json")

    _, spots, _, _ = sra_verdicts.registry_inputs(_parse(tmp_path), registry, store)

    assert spots == {ACC: 7, "SRR2": 11}


def test_registry_inputs_prefers_registry_metadata_over_xml(tmp_path):
    _write_registry(tmp_path, {ACC: {"metadata": {"run_total_spots": 5, "run_size": 123}}})
    _write_xml(tmp_path / "metadata", total_spots=7)
    registry = load_registry(tmp_path / "prepare_registry.json")

    excluded, spots, truncated, sizes = sra_verdicts.registry_inputs(_parse(tmp_path), registry, None)

    assert spots == {ACC: 5}
    assert sizes == {ACC: 123}
    assert excluded == set() and truncated == set()


def test_registry_inputs_without_verification_collects_no_spots(tmp_path):
    _write_xml(tmp_path / "metadata", total_spots=7)
    registry = load_registry(tmp_path / "prepare_registry.json")

    _, spots, _, _ = sra_verdicts.registry_inputs(_parse(tmp_path, "--no-verify-downloads"), registry, None)

    assert spots == {}


def test_registry_inputs_tolerates_a_missing_accessions_file(tmp_path):
    _write_xml(tmp_path / "metadata", total_spots=7)
    args = _parse(tmp_path)
    (tmp_path / "accessions.txt").unlink()
    registry = load_registry(tmp_path / "prepare_registry.json")

    _, spots, _, _ = sra_verdicts.registry_inputs(args, registry, None)

    assert spots == {}


def test_registry_inputs_dry_run_is_empty(tmp_path):
    _write_xml(tmp_path / "metadata", total_spots=7)
    registry = load_registry(tmp_path / "prepare_registry.json")

    assert sra_verdicts.registry_inputs(_parse(tmp_path, "--dry-run"), registry, None) == (set(), {}, set(), {})


# ------------------------------------------------------------------ present_verdicts


def test_present_verdicts_skips_recorded_and_linked_accessions(tmp_path):
    fastq = tmp_path / "fastq"
    for acc in ("SRR1", "SRR2", "SRR3", "SRR4"):
        _write_fastq(fastq / acc, reads=2, accession=acc)
    (fastq / "SRR5").symlink_to(fastq / "SRR4")
    registry = Registry()
    registry.datasets["SRR2"] = {"download": {"attempts": 1, "state": "downloaded"}}
    registry.datasets["SRR3"] = {"download": {"attempts": 1, "state": "failed"}}
    expected = {"SRR1": 2, "SRR2": 2, "SRR3": 4, "SRR5": 2}

    verdicts = sra_verdicts.present_verdicts(fastq, ["SRR1", "SRR2", "SRR3", "SRR4", "SRR5"], registry, expected)

    assert set(verdicts) == {"SRR1", "SRR3", "SRR4"}
    assert verdicts["SRR4"] == UNVERIFIED
    assert verdicts["SRR1"]["verdict"] == "complete"
    assert verdicts["SRR3"]["verdict"] == "truncated"
    assert verdicts["SRR3"]["reads_r1"] == 2
    assert "bytes_total" not in verdicts["SRR3"]


def test_present_verdicts_logs_an_unreadable_file_and_moves_on(tmp_path, caplog):
    fastq = tmp_path / "fastq"
    (fastq / ACC).mkdir(parents=True)
    (fastq / ACC / f"{ACC}_1.fastq.gz").write_bytes(b"not gzip")

    verdicts = sra_verdicts.present_verdicts(fastq, [ACC], None, {ACC: 4})

    assert verdicts == {ACC: dict(UNVERIFIED, expected_spots=4)}
    assert "Could not verify" in caplog.text


# ------------------------------------------------------------------ linked_verdict


def test_linked_verdict_recomputes_from_sidecar_reads_and_known_spots():
    previous = {"verdict": "truncated", "ratio": 0.5, "expected_spots": 10, "reads_r1": 5}
    sidecar = {"verdict": "unverified", "ratio": None, "expected_spots": None, "reads_r1": 10}

    verdict = sra_verdicts.linked_verdict(previous, sidecar, 10)

    assert verdict["verdict"] == "complete"
    assert verdict["reads_r1"] == 10


def test_linked_verdict_uses_the_sidecar_spot_count_when_none_is_given():
    sidecar = {"verdict": "unverified", "ratio": None, "expected_spots": 20, "reads_r1": 10}

    assert sra_verdicts.linked_verdict(None, sidecar, None)["verdict"] == "truncated"


def test_linked_verdict_keeps_complete_and_truncated_over_unverified():
    unverified = {"verdict": "unverified", "ratio": None, "expected_spots": None, "reads_r1": None}
    for kept in ("complete", "truncated"):
        previous = rb.Verdict.from_dict({"verdict": kept, "ratio": 1.0, "expected_spots": 4, "reads_r1": 4})
        assert sra_verdicts.linked_verdict(previous, unverified, None)["verdict"] == kept


def test_linked_verdict_without_a_sidecar_keeps_the_previous_verdict():
    previous = {"verdict": "truncated", "ratio": 0.5, "expected_spots": 10, "reads_r1": 5}

    assert sra_verdicts.linked_verdict(previous, None, 10) == previous
    assert sra_verdicts.linked_verdict(None, None, 10) is None


# ------------------------------------------------------------------ fix round 1


def test_present_accession_without_a_spot_count_is_recorded_unverified(tmp_path):
    """No count anywhere: the verdict is ``unverified``, recorded without reading the files."""
    _write_fastq(tmp_path / "fastq" / ACC, reads=2)
    args = _parse(tmp_path)

    with patch.object(sra_verdicts, "verify_download", side_effect=AssertionError("no file read")):
        DownloadSraCommand()._record_run_outcomes(args, {"already_downloaded_accessions": [ACC]}, tmp_path / "fastq")

    block = rb.download_block(load_registry(args.registry), ACC)
    assert block.state == "downloaded"
    assert block.complete is not None
    assert block.complete.to_dict() == UNVERIFIED


def test_present_accession_without_a_spot_count_keeps_a_recorded_truncated_verdict(tmp_path):
    _write_fastq(tmp_path / "fastq" / ACC, reads=2)
    _write_registry(tmp_path, {ACC: {"download": {"attempts": 1, "state": "failed", "complete": TRUNCATED}}})
    args = _parse(tmp_path)

    DownloadSraCommand()._record_run_outcomes(args, {"already_downloaded_accessions": [ACC]}, tmp_path / "fastq")

    block = rb.download_block(load_registry(args.registry), ACC)
    assert block.state == "downloaded"
    assert block.complete.verdict == "truncated"
    assert block.complete.reads_r1 == 5


def test_present_store_link_keeps_a_recorded_truncated_verdict(tmp_path):
    """The already-present store path never turns ``truncated`` into ``unverified`` either."""
    store = init_store(tmp_path / "store")
    _write_registry(tmp_path, {ACC: {"download": {"attempts": 1, "state": "failed", "complete": TRUNCATED}}})
    args = _parse(tmp_path)

    with (
        patch("prepare.cli.commands.sra.is_store_link", return_value=True),
        patch.object(DownloadSraCommand, "_sidecar_completeness", staticmethod(lambda _s, _a: dict(UNVERIFIED))),
        patch("prepare.cli.commands.sra.record_usage_many"),
    ):
        DownloadSraCommand()._record_run_outcomes(
            args, {"already_downloaded_accessions": [ACC]}, tmp_path / "fastq", store=store
        )

    block = rb.download_block(load_registry(args.registry), ACC)
    assert block.source == "store"
    assert block.complete.verdict == "truncated"
    assert block.complete.reads_r1 == 5


def _redownload(tmp_path, monkeypatch, message, recorded=TRUNCATED):
    """Run download_sra over ``ACC`` recorded with ``recorded`` with a fake download reporting ``message``."""
    install_fake_tools(tmp_path / "bin", tmp_path / "barrier")
    monkeypatch.setenv("PATH", str(tmp_path / "bin"))
    _write_fastq(tmp_path / "fastq" / ACC, reads=10)
    _write_registry(tmp_path, {ACC: {"download": {"attempts": 1, "state": "downloaded", "complete": recorded}}})

    def fake_download_sra(**kwargs):
        kwargs["on_result"](ACC, True, message)
        return {"total": 1, "successful": 1, "failed": 0, "failed_accessions": [], "results": {}}

    with patch("prepare.cli.commands.sra.download_sra", side_effect=fake_download_sra):
        assert DownloadSraCommand().execute(_parse(tmp_path, "--redownload-truncated")) == 0
    return rb.download_block(load_registry(tmp_path / "prepare_registry.json"), ACC)


def test_unverified_redownload_keeps_a_recorded_truncated_verdict(tmp_path, monkeypatch):
    block = _redownload(tmp_path, monkeypatch, "Downloaded 2 files, unverified")

    assert block.attempts == 2
    assert block.complete.verdict == "truncated"
    assert block.complete.reads_r1 == 5


def test_unverified_redownload_keeps_a_recorded_complete_verdict(tmp_path, monkeypatch):
    """A plain re-download with no spot count is no evidence either way: ``complete`` stays."""
    recorded = {"verdict": "complete", "ratio": 1.0, "expected_spots": 10, "reads_r1": 10}
    block = _redownload(tmp_path, monkeypatch, "Downloaded 2 files, unverified", recorded=recorded)

    assert block.attempts == 2
    assert block.complete.verdict == "complete"
    assert block.complete.reads_r1 == 10
    assert block.complete.expected_spots == 10


def test_complete_redownload_replaces_a_recorded_truncated_verdict(tmp_path, monkeypatch):
    block = _redownload(tmp_path, monkeypatch, "Downloaded 2 files, complete (10 of 10 spots)")

    assert block.complete.verdict == "complete"
    assert block.complete.reads_r1 == 10
    assert block.complete.expected_spots == 10


def test_registry_inputs_falls_back_to_the_recorded_verdict_count_last(tmp_path):
    _write_registry(
        tmp_path,
        {
            ACC: {"download": {"attempts": 1, "state": "downloaded", "complete": TRUNCATED}},
            "SRR2": {"download": {"attempts": 1, "state": "downloaded", "complete": TRUNCATED}},
        },
    )
    _write_xml(tmp_path / "metadata", total_spots=12, accession="SRR2")
    (tmp_path / "accessions.txt").write_text(f"{ACC}\nSRR2\n")
    registry = load_registry(tmp_path / "prepare_registry.json")

    _, spots, _, _ = sra_verdicts.registry_inputs(_parse(tmp_path), registry, None)

    assert spots == {ACC: 10, "SRR2": 12}


def test_recorded_verdict_without_a_new_verdict_leaves_the_record_alone():
    registry = Registry()
    registry.datasets[ACC] = {"download": {"attempts": 1, "state": "failed", "complete": TRUNCATED}}

    assert sra_verdicts.recorded_verdict(registry, ACC, None, 10, linked=False) is None
    assert sra_verdicts.recorded_verdict(registry, ACC, None, 10, linked=True) is None
