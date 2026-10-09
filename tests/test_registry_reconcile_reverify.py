"""Tests for the re-verification part of ``status --reconcile`` (prepare.data.registry_reconcile).

Covers re-checking ``unverified`` download verdicts once a spot count is known, filling registry
metadata from the XML folder before the verdicts, and leaving the records of a store that is not
mounted alone. Every store lives under ``tmp_path``; no real store, tool or network is used.
"""

import copy
import itertools
import json
import logging
import os
import shutil
from pathlib import Path
from typing import List

import pytest

import prepare.data.registry_reconcile as rr
from prepare.data import registry as reg
from prepare.data import registry_blocks as rb
from prepare.data.registry_batch import registry_update
from prepare.data.registry_reconcile import StoreReconcileReport, apply_reconcile, reconcile, scan_reconcile
from prepare.data.sra import fastq as fastq_module
from prepare.data.sra.spots import verdict_for_count

_UNVERIFIED = {"method": "unverified", "ratio": None, "verdict": "unverified", "expected_spots": None, "reads_r1": None}

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


@pytest.fixture(autouse=True)
def isolated_env(tmp_path, monkeypatch):
    """Keep every test away from the developer's store and configuration."""
    monkeypatch.setenv("PREPARE_DATA", str(tmp_path / "no-store"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setattr(reg, "_now", lambda: "2026-10-01T00:00:00+00:00")


def _fastq(path: Path, reads: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(reads)))
    return path


def _paths(tmp_path: Path) -> reg.ProjectPaths:
    root = tmp_path / "project"
    return reg.ProjectPaths(
        fastq=root / "fastq",
        metadata=root / "metadata",
        genomes=root / "genomes",
        targeted=root / "targeted",
        matches=root / "matches",
    )


def _registry(tmp_path: Path) -> reg.Registry:
    return reg.Registry(path=tmp_path / "project" / reg.REGISTRY_FILENAME)


def _paired_download(registry: reg.Registry, paths: reg.ProjectPaths, acc: str, spots) -> None:
    """A plain paired download (3 reads per mate, 1 orphan read) recorded as unverified."""
    _fastq(paths.fastq / acc / f"{acc}_1.fastq", 3)
    _fastq(paths.fastq / acc / f"{acc}_2.fastq", 3)
    _fastq(paths.fastq / acc / f"{acc}.fastq", 1)
    reg.record_download(registry, acc, "downloaded", paths.fastq, complete=dict(_UNVERIFIED))
    reg.record_metadata(registry, acc, paths.metadata / f"{acc}_metadata.xml", {"run_total_spots": spots})


def _store_with(tmp_path: Path, acc: str, reads_per_mate=None, spots=None, under: str = "vol"):
    """A store under ``tmp_path/<under>`` holding ``acc``, with a sidecar recording an unverified verdict."""
    from prepare.store.layout import init_store, sidecar_path
    from prepare.store.sidecar import Sidecar, write_sidecar

    store = init_store(tmp_path / under / "store")
    _fastq(store.sra / acc / f"{acc}_1.fastq", reads_per_mate or 1)
    write_sidecar(
        sidecar_path(store, acc),
        Sidecar(
            accession=acc,
            state="complete",
            reads_per_mate=reads_per_mate,
            ncbi={"spots": spots} if spots else {},
            completeness={"method": "unverified", "ratio": None, "verdict": "unverified"},
        ),
    )
    return store


def _linked_download(tmp_path: Path, acc: str, reads_per_mate=None, spots=None, under: str = "vol"):
    from prepare.store.link import link_dataset

    paths = _paths(tmp_path)
    store = _store_with(tmp_path, acc, reads_per_mate=reads_per_mate, spots=spots, under=under)
    link_dataset(paths.fastq, acc, store)
    registry = _registry(tmp_path)
    registry.store = {"root": str(store.root)}
    reg.record_download(
        registry, acc, "downloaded", paths.fastq, attempt=False, source="store", store_name=acc, complete=_UNVERIFIED
    )
    reg.update_linked(registry, acc, add=True)
    return registry, paths, store


def _refuse(*_args, **_kwargs):
    raise AssertionError("no FASTQ may be read here")


class TestUnverifiedPlainDownloadIsRechecked:
    def test_mate_one_and_orphan_are_counted_once_and_the_full_verdict_is_written(self, tmp_path, monkeypatch):
        paths = _paths(tmp_path)
        registry = _registry(tmp_path)
        _paired_download(registry, paths, "SRR1", 4)
        counted: List[str] = []
        real_count = fastq_module.count_fastq_reads

        def counting(path):
            counted.append(Path(path).name)
            return real_count(path)

        monkeypatch.setattr(fastq_module, "count_fastq_reads", counting)

        report = reconcile(registry, paths)

        assert sorted(counted) == ["SRR1.fastq", "SRR1_1.fastq"]
        assert rb.download_verdict(registry, "SRR1").to_dict() == verdict_for_count(4, 4)
        assert isinstance(report, StoreReconcileReport)
        assert report.verdicts_rechecked == ["SRR1"]

    def test_a_short_download_is_rechecked_as_truncated(self, tmp_path):
        paths = _paths(tmp_path)
        registry = _registry(tmp_path)
        _paired_download(registry, paths, "SRR1", 100)

        report = reconcile(registry, paths)

        assert rb.download_verdict(registry, "SRR1").to_dict() == verdict_for_count(4, 100)
        assert rb.download_verdict(registry, "SRR1").verdict == "truncated"
        assert report.verdicts_rechecked == ["SRR1"]

    def test_without_a_spot_count_it_stays_unverified_and_is_not_counted(self, tmp_path, monkeypatch):
        paths = _paths(tmp_path)
        registry = _registry(tmp_path)
        _paired_download(registry, paths, "SRR1", None)
        monkeypatch.setattr(fastq_module, "count_fastq_reads", _refuse)

        report = reconcile(registry, paths)

        assert rb.download_verdict(registry, "SRR1").verdict == "unverified"
        assert report.verdicts_rechecked == []

    def test_a_known_verdict_is_not_recounted(self, tmp_path, monkeypatch):
        paths = _paths(tmp_path)
        registry = _registry(tmp_path)
        _paired_download(registry, paths, "SRR1", 4)
        reg.set_download_verdict(registry, "SRR1", verdict_for_count(4, 4))
        monkeypatch.setattr(fastq_module, "count_fastq_reads", _refuse)

        report = reconcile(registry, paths)

        assert rb.download_verdict(registry, "SRR1").verdict == "complete"
        assert report.verdicts_rechecked == []


class TestUnverifiedStoreLinkIsRechecked:
    def test_the_sidecar_read_count_is_compared_without_reading_fastq(self, tmp_path, monkeypatch):
        registry, paths, _store = _linked_download(tmp_path, "SRR1", reads_per_mate=10)
        reg.record_metadata(registry, "SRR1", paths.metadata / "SRR1_metadata.xml", {"run_total_spots": 10})
        monkeypatch.setattr(fastq_module, "count_fastq_reads", _refuse)
        monkeypatch.setattr(rr, "count_fastq_reads", _refuse)
        monkeypatch.setattr(rr, "verify_download", _refuse)

        report = reconcile(registry, paths)

        assert rb.download_verdict(registry, "SRR1").to_dict() == verdict_for_count(10, 10)
        assert report.verdicts_rechecked == ["SRR1"]

    def test_the_sidecar_spot_count_is_used_when_the_registry_has_none(self, tmp_path, monkeypatch):
        registry, paths, _store = _linked_download(tmp_path, "SRR1", reads_per_mate=5, spots=10)
        monkeypatch.setattr(fastq_module, "count_fastq_reads", _refuse)

        report = reconcile(registry, paths)

        assert rb.download_verdict(registry, "SRR1").to_dict() == verdict_for_count(5, 10)
        assert report.verdicts_rechecked == ["SRR1"]

    def test_a_sidecar_without_a_read_count_leaves_it_unverified(self, tmp_path, monkeypatch):
        registry, paths, _store = _linked_download(tmp_path, "SRR1", reads_per_mate=None, spots=10)
        monkeypatch.setattr(fastq_module, "count_fastq_reads", _refuse)

        report = reconcile(registry, paths)

        assert rb.download_verdict(registry, "SRR1").verdict == "unverified"
        assert report.verdicts_rechecked == []


class TestMetadataFromXmlThenVerdict:
    def test_xml_fills_the_spot_count_and_the_verdict_in_one_run(self, tmp_path):
        paths = _paths(tmp_path)
        registry = _registry(tmp_path)
        _paired_download(registry, paths, "SRR1", None)
        paths.metadata.mkdir(parents=True, exist_ok=True)
        (paths.metadata / "SRR1_metadata.xml").write_text(_XML.format(accession="SRR1", spots=4))

        report = reconcile(registry, paths)

        assert rb.metadata_block(registry, "SRR1").run_total_spots == 4
        assert rb.metadata_block(registry, "SRR1").platform == "ILLUMINA"
        assert rb.download_verdict(registry, "SRR1").to_dict() == verdict_for_count(4, 4)
        assert report.metadata_filled == ["SRR1"]
        assert report.verdicts_rechecked == ["SRR1"]

    def test_the_inferred_mark_of_the_metadata_block_is_kept(self, tmp_path):
        paths = _paths(tmp_path)
        registry = _registry(tmp_path)
        _paired_download(registry, paths, "SRR1", None)
        rb.mark_inferred(registry, "SRR1", "metadata")
        paths.metadata.mkdir(parents=True, exist_ok=True)
        (paths.metadata / "SRR1_metadata.xml").write_text(_XML.format(accession="SRR1", spots=4))

        reconcile(registry, paths)

        assert rb.metadata_block(registry, "SRR1").inferred is True


class TestStoreNotMounted:
    def test_links_into_an_unmounted_store_are_not_marked_missing(self, tmp_path, caplog):
        registry, paths, _store = _linked_download(tmp_path, "SRR1", reads_per_mate=2, spots=2)
        reg.set_download_verdict(registry, "SRR1", verdict_for_count(2, 2))
        before = copy.deepcopy((registry.datasets, registry.store))
        # Unmounting the volume takes the whole store, sra/ folder included, out of reach.
        os.rename(tmp_path / "vol", tmp_path / "vol-unmounted")

        with caplog.at_level(logging.WARNING, logger=rr.__name__):
            report = reconcile(registry, paths)

        # Nothing about the dataset or the store link changes: state, dates, verdict, linked list.
        assert (registry.datasets, registry.store) == before
        assert report.recorded_missing == []
        assert report.store_unavailable == ["SRR1"]
        assert report.dangling_links == ["SRR1"]
        assert "not mounted" in caplog.text
        assert (paths.fastq / "SRR1").is_symlink()  # the link is never removed

    def test_a_dataset_gone_from_a_mounted_store_is_marked_missing(self, tmp_path):
        registry, paths, store = _linked_download(tmp_path, "SRR1", reads_per_mate=2, spots=2)
        shutil.rmtree(store.sra / "SRR1")

        report = reconcile(registry, paths)

        assert rb.download_block(registry, "SRR1").state == "missing"
        assert report.recorded_missing == ["SRR1"]
        assert report.store_unavailable == []
        assert report.dangling_links == ["SRR1"]
        assert (paths.fastq / "SRR1").is_symlink()


class TestReRecordKeepsTheStoreSource:
    def test_a_store_dataset_back_on_disk_is_rerecorded_as_linked_from_the_store(self, tmp_path):
        registry, paths, store = _linked_download(tmp_path, "SRR1", reads_per_mate=2, spots=2)
        # An earlier reconcile, run while the store was unmounted, marked it missing.
        download = rb.download_block(registry, "SRR1")
        download.state = "missing"
        rb.set_download_block(registry, "SRR1", download)
        reg.update_linked(registry, "SRR1", add=False)

        report = reconcile(registry, paths)

        assert report.untracked_fastq == ["SRR1"]
        block = rb.download_block(registry, "SRR1")
        assert block.state == "downloaded"
        assert block.source == "store"
        assert block.store_name == "SRR1"
        assert "SRR1" in rb.store_block(registry).linked

    def test_a_real_folder_replacing_a_store_link_is_recorded_as_a_plain_download(self, tmp_path):
        registry, paths, _store = _linked_download(tmp_path, "SRR1", reads_per_mate=2, spots=2)
        download = rb.download_block(registry, "SRR1")
        download.state = "missing"
        rb.set_download_block(registry, "SRR1", download)
        reg.update_linked(registry, "SRR1", add=False)
        # Outside PREPARE, the link was replaced by a folder of the project's own reads.
        (paths.fastq / "SRR1").unlink()
        _fastq(paths.fastq / "SRR1" / "SRR1_1.fastq", 2)

        report = reconcile(registry, paths)

        assert report.untracked_fastq == ["SRR1"]
        block = rb.download_block(registry, "SRR1")
        assert block.state == "downloaded"
        assert block.source is None
        assert block.store_name is None
        assert "SRR1" not in (rb.store_block(registry).linked or [])

    def test_a_plain_download_back_on_disk_has_no_source(self, tmp_path):
        paths = _paths(tmp_path)
        registry = _registry(tmp_path)
        _fastq(paths.fastq / "SRR1" / "SRR1_1.fastq", 2)
        reg.record_download(registry, "SRR1", "missing", paths.fastq)

        reconcile(registry, paths)

        block = rb.download_block(registry, "SRR1")
        assert block.state == "downloaded"
        assert block.source is None
        assert "SRR1" not in (rb.store_block(registry).linked or [])


class TestScanApplySplit:
    def _project(self, tmp_path: Path):
        registry, paths, store = _linked_download(tmp_path, "SRR5", reads_per_mate=10, spots=10)
        _paired_download(registry, paths, "SRR1", 4)
        _paired_download(registry, paths, "SRR2", None)
        paths.metadata.mkdir(parents=True, exist_ok=True)
        (paths.metadata / "SRR2_metadata.xml").write_text(_XML.format(accession="SRR2", spots=100))
        return registry, paths

    def test_apply_reads_no_fastq_sidecar_or_xml(self, tmp_path, monkeypatch):
        import prepare.data.metadata_fields as mf
        import prepare.data.sra.spots as spots_module

        registry, paths = self._project(tmp_path)
        expected = copy.deepcopy(registry)
        reconcile(expected, paths)
        plan = scan_reconcile(registry, paths)
        for module, name in (
            (fastq_module, "count_fastq_reads"),
            (rr, "count_fastq_reads"),
            (rr, "verify_download"),
            (mf, "parse_metadata_xml"),
            (spots_module, "_sidecar_spots"),
        ):
            monkeypatch.setattr(module, name, _refuse)

        report = apply_reconcile(registry, plan)

        assert registry.datasets == expected.datasets
        assert report.verdicts_rechecked == ["SRR1", "SRR2", "SRR5"]
        assert report.metadata_filled == ["SRR2"]
        assert rb.download_verdict(registry, "SRR2").verdict == "truncated"

    def test_scan_leaves_the_snapshot_unchanged(self, tmp_path):
        registry, paths = self._project(tmp_path)
        before = copy.deepcopy(registry.datasets)

        scan_reconcile(registry, paths)

        assert registry.datasets == before

    def test_a_verdict_recorded_after_the_scan_is_kept(self, tmp_path):
        registry, paths = self._project(tmp_path)
        plan = scan_reconcile(registry, paths)
        # Another process re-downloads SRR1 and records a verdict before the apply.
        recorded = verdict_for_count(3, 4)
        reg.set_download_verdict(registry, "SRR1", recorded)

        report = apply_reconcile(registry, plan)

        assert rb.download_verdict(registry, "SRR1").to_dict() == recorded
        assert "SRR1" not in report.verdicts_rechecked

    def test_metadata_recorded_after_the_scan_is_kept(self, tmp_path):
        registry, paths = self._project(tmp_path)
        plan = scan_reconcile(registry, paths)
        reg.record_metadata(registry, "SRR2", paths.metadata / "SRR2_metadata.xml", {"run_total_spots": 7})

        report = apply_reconcile(registry, plan)

        assert rb.metadata_block(registry, "SRR2").run_total_spots == 7
        assert report.metadata_filled == []
        # The scan counted SRR2 against 100 spots, not 7, so its verdict waits for the next run.
        assert rb.download_verdict(registry, "SRR2").verdict == "unverified"


class TestReconcileIsIdempotent:
    """A second ``status --reconcile`` on an unchanged project writes nothing new."""

    @staticmethod
    def _status_reconcile(registry_file: Path, paths: reg.ProjectPaths) -> StoreReconcileReport:
        # What the status command does: scan a snapshot without the lock, apply under it.
        plan = scan_reconcile(reg.load_registry(registry_file), paths)
        return registry_update(registry_file, lambda registry: apply_reconcile(registry, plan))

    @staticmethod
    def _content(registry_file: Path) -> dict:
        # Every save stamps "updated"; everything else in the file must stay as it was.
        data = json.loads(registry_file.read_text())
        data.pop("updated")
        return data

    def test_a_second_reconcile_changes_nothing(self, tmp_path, monkeypatch):
        ticks = itertools.count()
        monkeypatch.setattr(reg, "_now", lambda: f"2026-10-01T00:00:{next(ticks) % 60:02d}+00:00")
        registry, paths, _store = _linked_download(tmp_path, "SRR5", reads_per_mate=10, spots=10)
        _paired_download(registry, paths, "SRR1", 4)
        _paired_download(registry, paths, "SRR2", None)
        _paired_download(registry, paths, "SRR3", None)
        paths.metadata.mkdir(parents=True, exist_ok=True)
        (paths.metadata / "SRR2_metadata.xml").write_text(_XML.format(accession="SRR2", spots=100))
        # An XML for a run NCBI has not loaded: no spot count, so nothing to fill, now or later.
        no_spots = _XML.replace('total_spots="{spots}" ', "").format(accession="SRR3")
        (paths.metadata / "SRR3_metadata.xml").write_text(no_spots)
        _fastq(paths.fastq / "SRR7" / "SRR7_1.fastq", 2)
        registry_file = registry.path
        reg.save_registry(registry, registry_file)

        first = self._status_reconcile(registry_file, paths)
        after_first = self._content(registry_file)
        second = self._status_reconcile(registry_file, paths)

        assert first.metadata_filled == ["SRR2"]
        assert first.verdicts_rechecked == ["SRR1", "SRR2", "SRR5"]
        assert first.untracked_fastq == ["SRR7"]
        assert self._content(registry_file) == after_first
        assert second.metadata_filled == []
        assert second.verdicts_rechecked == []
        assert second.untracked_fastq == []
        assert second.recorded_missing == []
        assert rb.metadata_block(reg.load_registry(registry_file), "SRR3").run_total_spots is None
