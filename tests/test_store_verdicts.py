"""Tests for the one verdict a store accession records (``sra_verdicts.store_verdict``).

``store_link``, ``store_adopt`` and ``download_sra`` record the store sidecar's verdict judged
against the project's spot count and merged with the verdict on file, so a relink never turns a
recorded ``truncated`` verdict into ``unverified``, and a ready copy whose sidecar is ``unverified``
but short against that count is refused by ``store_link`` without ``--accept-partial``. A store
copy that ``download_sra`` keeps but does not link (an ``incomplete:`` result) records the
sidecar's verdict too, so a copy-mode project copy is skipped by ``extract_target_reads``. The spot
count follows ``expected_spots``'s order: registry metadata, store sidecar, metadata XML, previous
verdict. Every store lives under ``tmp_path``; no real tool, store or network is used.
"""

import argparse
import json
import threading
from pathlib import Path
from unittest.mock import patch

from prepare.cli.commands import sra_verdicts
from prepare.cli.commands.read_extraction import ExtractTargetReadsCommand
from prepare.cli.commands.sra import DownloadSraCommand
from prepare.cli.commands.store import StoreAdoptCommand, StoreLinkCommand
from prepare.data import registry_blocks as rb
from prepare.data.registry import Registry, load_registry, save_registry
from prepare.store.layout import init_store, sidecar_path, sra_dir
from prepare.store.link import link_dataset
from prepare.store.sidecar import Sidecar, read_sidecar, write_sidecar

from tests.helpers_processes import install_fake_tools
from tests.test_cli_sra_verdicts import _write_xml

ACC = "SRR1"
FASTQ_RECORD = "@r\nACGT\n+\nIIII\n"


def _registry(project, datasets):
    """Write a registry holding ``datasets`` into ``project``; return its path."""
    path = project / "prepare_registry.json"
    registry = Registry(path=path)
    registry.datasets.update(json.loads(json.dumps(datasets)))
    save_registry(registry, path)
    return path


def _store_copy(paths, reads, state="complete", verdict="unverified", spots=None):
    """A store copy of ``ACC`` holding ``reads`` records, with a sidecar recording ``reads_per_mate=reads``."""
    acc_dir = sra_dir(paths, ACC)
    acc_dir.mkdir(parents=True, exist_ok=True)
    (acc_dir / f"{ACC}_1.fastq").write_text(FASTQ_RECORD * max(reads or 0, 1))
    write_sidecar(
        sidecar_path(paths, ACC),
        Sidecar(
            accession=ACC,
            state=state,
            reads_per_mate=reads,
            ncbi={"spots": spots} if spots else {},
            completeness={"method": "unverified", "ratio": None, "verdict": verdict},
        ),
    )
    return acc_dir


def _link_args(project, root, accept_partial=False):
    return argparse.Namespace(
        accessions=[ACC],
        fastq_folder=str(project / "fastq"),
        registry=str(project / "prepare_registry.json"),
        data_root=str(root),
        link_mode="auto",
        accept_partial=accept_partial,
    )


def _verdict(registry_path):
    return rb.download_verdict(load_registry(registry_path), ACC)


TRUNCATED_10_OF_100 = {"verdict": "truncated", "ratio": 0.1, "expected_spots": 100, "reads_r1": 10}


# ------------------------------------------------------------------------------------ store_link


def test_store_link_keeps_a_recorded_truncated_verdict_over_an_unverified_sidecar(tmp_path):
    """A sidecar with no read count is no evidence: the recorded ``truncated`` verdict stays."""
    paths = init_store(tmp_path / "store")
    _store_copy(paths, reads=None)
    project = tmp_path / "project"
    project.mkdir()
    registry_path = _registry(
        project,
        {
            ACC: {
                "download": {
                    "attempts": 1,
                    "state": "missing",
                    "complete": {"verdict": "truncated", "ratio": 0.5, "expected_spots": 10, "reads_r1": 5},
                }
            }
        },
    )

    assert StoreLinkCommand().execute(_link_args(project, paths.root)) == 0

    verdict = _verdict(registry_path)
    assert verdict.verdict == "truncated"
    assert verdict.reads_r1 == 5


def test_store_link_refuses_a_short_unverified_copy_without_accept_partial(tmp_path):
    paths = init_store(tmp_path / "store")
    _store_copy(paths, reads=10)
    project = tmp_path / "project"
    project.mkdir()
    registry_path = _registry(project, {ACC: {"metadata": {"run_total_spots": 100}}})

    assert StoreLinkCommand().execute(_link_args(project, paths.root)) == 1

    assert not (project / "fastq" / ACC).exists()
    sidecar = read_sidecar(sidecar_path(paths, ACC))
    assert sidecar.state == "partial"
    assert sidecar.completeness["verdict"] == "truncated"
    assert sidecar.ncbi["spots"] == 100
    assert _verdict(registry_path) is None

    assert StoreLinkCommand().execute(_link_args(project, paths.root, accept_partial=True)) == 0

    assert (project / "fastq" / ACC).is_symlink()
    verdict = _verdict(registry_path)
    assert verdict.verdict == "truncated"
    assert verdict.reads_r1 == 10
    assert verdict.expected_spots == 100


def test_store_link_of_a_recorded_truncated_short_copy_is_refused_and_keeps_the_verdict(tmp_path):
    """The branch review's probe: registry truncated 10 of 100, sidecar complete/unverified with 10 reads."""
    paths = init_store(tmp_path / "store")
    _store_copy(paths, reads=10)
    project = tmp_path / "project"
    (project / "fastq").mkdir(parents=True)
    registry_path = _registry(
        project,
        {
            ACC: {
                "metadata": {"run_total_spots": 100},
                "download": {"attempts": 1, "state": "downloaded", "complete": TRUNCATED_10_OF_100},
            }
        },
    )

    assert StoreLinkCommand().execute(_link_args(project, paths.root)) == 1

    assert not (project / "fastq" / ACC).exists()
    verdict = _verdict(registry_path)
    assert verdict.verdict == "truncated"
    assert verdict.reads_r1 == 10
    assert verdict.expected_spots == 100


def test_store_link_judges_a_short_copy_by_the_project_metadata_xml(tmp_path):
    """No count in the registry or sidecar: the project's metadata XML supplies it."""
    paths = init_store(tmp_path / "store")
    _store_copy(paths, reads=2)
    project = tmp_path / "project"
    project.mkdir()
    _registry(project, {})
    _write_xml(project / "metadata", total_spots=10)

    assert StoreLinkCommand().execute(_link_args(project, paths.root)) == 1
    assert read_sidecar(sidecar_path(paths, ACC)).state == "partial"


def test_store_link_leaves_the_sidecar_alone_while_another_run_holds_the_lock(tmp_path, caplog):
    from prepare.store.locks import dataset_lock

    paths = init_store(tmp_path / "store")
    _store_copy(paths, reads=10)
    project = tmp_path / "project"
    project.mkdir()
    _registry(project, {ACC: {"metadata": {"run_total_spots": 100}}})

    held, release = threading.Event(), threading.Event()

    def _hold():
        with dataset_lock(paths, ACC):
            held.set()
            release.wait(10)

    holder = threading.Thread(target=_hold)
    holder.start()
    try:
        assert held.wait(10)
        assert StoreLinkCommand().execute(_link_args(project, paths.root)) == 1
    finally:
        release.set()
        holder.join()

    assert read_sidecar(sidecar_path(paths, ACC)).state == "complete"
    assert "another run holds the dataset lock" in caplog.text


# ----------------------------------------------------------------------------------- store_adopt


def _adopt(project, root):
    args = argparse.Namespace(
        fastq_folder="fastq",
        data_root=str(root),
        registry=str(project / "prepare_registry.json"),
        move=True,
        dry_run=False,
        compress=False,
        metadata_folder="metadata",
        lock_wait=0.0,
    )
    return StoreAdoptCommand().execute(args)


def _project_copy(project, reads):
    acc_dir = project / "fastq" / ACC
    acc_dir.mkdir(parents=True)
    (acc_dir / f"{ACC}.fastq").write_text(FASTQ_RECORD * reads)


def test_store_adopt_records_truncated_for_a_copy_short_against_the_registry_count(tmp_path, monkeypatch):
    """Adopt with no metadata XML: the registry's spot count still judges the adopted copy."""
    paths = init_store(tmp_path / "store")
    project = tmp_path / "project"
    _project_copy(project, reads=1)
    monkeypatch.chdir(project)
    registry_path = _registry(project, {ACC: {"metadata": {"run_total_spots": 100}}})

    assert _adopt(project, paths.root) == 0

    verdict = _verdict(registry_path)
    assert verdict.verdict == "truncated"
    assert verdict.reads_r1 == 1
    assert verdict.expected_spots == 100


def test_store_adopt_keeps_a_recorded_truncated_verdict_without_a_spot_count(tmp_path, monkeypatch):
    paths = init_store(tmp_path / "store")
    project = tmp_path / "project"
    _project_copy(project, reads=1)
    monkeypatch.chdir(project)
    recorded = {"verdict": "truncated", "ratio": None, "expected_spots": None, "reads_r1": 1}
    registry_path = _registry(
        project, {ACC: {"download": {"attempts": 1, "state": "downloaded", "complete": recorded}}}
    )

    assert _adopt(project, paths.root) == 0

    assert _verdict(registry_path).verdict == "truncated"


# ---------------------------------------------------------------------------------- download_sra


def _download_args(project, root, *extra):
    parser = argparse.ArgumentParser()
    DownloadSraCommand().configure_parser(parser)
    accessions = project / "accessions.txt"
    accessions.write_text(f"{ACC}\n")
    return parser.parse_args(
        [
            "--accessions-file",
            str(accessions),
            "--fastq-folder",
            str(project / "fastq"),
            "--registry",
            str(project / "prepare_registry.json"),
            "--data-root",
            str(root),
            "--min-free-gb",
            "0",
            "--max-retries",
            "1",
            *extra,
        ]
    )


def _fake_fetch(calls, reads):
    """A download_accession stand-in writing ``reads`` records under <output_folder>/<acc>."""

    def _download(accession, output_folder, *args, **kwargs):
        calls.append(accession)
        acc_dir = Path(output_folder) / accession
        acc_dir.mkdir(parents=True, exist_ok=True)
        (acc_dir / f"{accession}_1.fastq").write_text(FASTQ_RECORD * reads)
        return True, f"Downloaded 1 files, truncated ({reads} of 100 spots)"

    return _download


def test_copy_mode_project_records_truncated_for_an_incomplete_store_result(tmp_path, monkeypatch):
    """A forced refetch leaves the store copy short: the copy-mode project copy is then skipped."""
    install_fake_tools(tmp_path / "bin", tmp_path / "barrier")
    monkeypatch.setenv("PATH", str(tmp_path / "bin"))
    paths = init_store(tmp_path / "store")
    _store_copy(paths, reads=1)
    project = tmp_path / "project"
    link_dataset(project / "fastq", ACC, paths, mode="copy")
    assert not (project / "fastq" / ACC).is_symlink()
    unverified = {"verdict": "unverified", "ratio": None, "expected_spots": None, "reads_r1": 1}
    registry_path = _registry(
        project,
        {
            ACC: {
                "metadata": {"run_total_spots": 100},
                "download": {"attempts": 1, "state": "downloaded", "complete": unverified},
            }
        },
    )
    calls = []

    with patch("prepare.data.sra.accession.download_accession", side_effect=_fake_fetch(calls, 1)):
        rc = DownloadSraCommand().execute(_download_args(project, paths.root, "--link-mode", "copy", "--force"))

    assert rc == 1
    assert calls == [ACC]
    registry = load_registry(registry_path)
    block = rb.download_block(registry, ACC)
    assert block.state == "failed"
    assert block.message.startswith("incomplete:")
    assert block.complete.verdict == "truncated"
    assert block.complete.reads_r1 == 1
    assert block.complete.expected_spots == 100
    # The project copy still carries its old, complete-looking sidecar; the registry verdict skips it.
    assert read_sidecar(project / "fastq" / ACC / f"{ACC}.json").state == "complete"
    assert ACC in ExtractTargetReadsCommand._unusable_downloads(registry, project / "fastq", [ACC])


def test_a_partial_store_result_is_fetched_once_per_run(tmp_path, monkeypatch):
    """An ``incomplete:`` result is settled: the retry pass does not fetch it a second time."""
    install_fake_tools(tmp_path / "bin", tmp_path / "barrier")
    monkeypatch.setenv("PATH", str(tmp_path / "bin"))
    paths = init_store(tmp_path / "store")
    project = tmp_path / "project"
    project.mkdir()
    _registry(project, {ACC: {"metadata": {"run_total_spots": 100}}})
    calls = []

    with patch("prepare.data.sra.accession.download_accession", side_effect=_fake_fetch(calls, 5)):
        rc = DownloadSraCommand().execute(_download_args(project, paths.root))

    assert rc == 1
    assert calls == [ACC]
    sidecar = read_sidecar(sidecar_path(paths, ACC))
    assert sidecar.state == "partial"
    assert sidecar.refetch is None


# ------------------------------------------------------------------------- spot count lookup order


def test_registry_inputs_follows_the_documented_spot_count_order(tmp_path):
    """Metadata, then the store sidecar, then the XML, then the previous verdict."""
    paths = init_store(tmp_path / "store")
    _store_copy(paths, reads=20, spots=20)
    project = tmp_path / "project"
    project.mkdir()
    previous = {"verdict": "truncated", "ratio": 0.5, "expected_spots": 10, "reads_r1": 5}
    registry_path = _registry(
        project, {ACC: {"download": {"attempts": 1, "state": "downloaded", "complete": previous}}}
    )
    _write_xml(project / "metadata", total_spots=30)
    args = _download_args(project, paths.root)
    registry = load_registry(registry_path)

    _, store_run, _, _ = sra_verdicts.registry_inputs(args, registry, paths)
    _, plain_run, _, _ = sra_verdicts.registry_inputs(args, registry, None)

    assert store_run == {ACC: 20}
    assert plain_run == {ACC: 30}

    (project / "metadata" / f"{ACC}_metadata.xml").unlink()
    _, plain_without_xml, _, _ = sra_verdicts.registry_inputs(args, registry, None)
    assert plain_without_xml == {ACC: 10}

    registry.datasets[ACC]["metadata"] = {"run_total_spots": 40}
    _, with_metadata, _, _ = sra_verdicts.registry_inputs(args, registry, paths)
    assert with_metadata == {ACC: 40}


def test_store_spot_count_prefers_a_known_count_and_reads_the_project_xml(tmp_path):
    paths = init_store(tmp_path / "store")
    project = tmp_path / "project"
    project.mkdir()
    registry = load_registry(_registry(project, {}))
    _write_xml(project / "metadata", total_spots=30)

    assert sra_verdicts.store_spot_count(registry, ACC, paths, known=7) == 7
    assert sra_verdicts.store_spot_count(registry, ACC, paths) == 30
    assert sra_verdicts.store_spot_count(None, ACC, paths) is None
    assert sra_verdicts.store_spot_count(None, ACC, paths, metadata_folder=project / "metadata") == 30


def test_store_verdict_without_a_sidecar_leaves_the_record_alone():
    registry = Registry()
    registry.datasets[ACC] = {"download": {"attempts": 1, "state": "downloaded", "complete": TRUNCATED_10_OF_100}}

    assert sra_verdicts.store_verdict(registry, ACC, None, 100) is None
    assert sra_verdicts.store_sidecar_verdict(None, ACC) is None
