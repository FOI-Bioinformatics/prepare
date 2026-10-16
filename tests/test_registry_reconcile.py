"""Tests for prepare.data.registry_reconcile and registry_update.

``reconcile`` was split into ``scan_reconcile`` (reads the disk, no lock) and ``apply_reconcile``
(fast, under the lock). ``_reference_reconcile`` below is the single-step function as it stood
before the split, kept so the two-step form can be checked against it. Two changes made since
then are mirrored in it: the report is a ``StoreReconcileReport`` (store, metadata and re-check
fields, all empty for this fixture), and a filled-in verdict is the full dict
``verdict_for_count`` returns rather than a hand-built subset.
"""

import copy
import gzip
import json
from pathlib import Path
from typing import Any, Dict, Optional

import pytest

from prepare.core.exceptions import DataAccessError
from prepare.data import registry as reg
from prepare.data import registry_blocks as rb
from prepare.data.registry_batch import registry_update
from prepare.data.registry_reconcile import StoreReconcileReport, apply_reconcile, reconcile, scan_reconcile
from prepare.data.sra import verify_download
from prepare.data.sra.spots import expected_spots, verdict_for_count

# ------------------------------------------------------------- reference (pre-split)


def _reference_reconcile(registry: reg.Registry, paths: reg.ProjectPaths) -> reg.ReconcileReport:
    from prepare.store.link import dangling_links

    report = StoreReconcileReport()
    on_disk = reg.scan_downloads(paths.fastq)
    for acc in registry.datasets:
        download = rb.download_block(registry, acc)
        if download is not None and download.state == "downloaded" and acc not in on_disk:
            download.state, download.date = "missing", reg._now()
            rb.set_download_block(registry, acc, download)
            report.recorded_missing.append(acc)
    tracked = set(reg.query(registry, "downloaded"))
    report.untracked_fastq = sorted(acc for acc in on_disk if acc not in tracked)
    for acc in report.untracked_fastq:
        reg.record_download(registry, acc, "downloaded", paths.fastq, attempt=False)
        rb.mark_inferred(registry, acc, "download", attempts=0)
    genome_ids = sorted(reg._genome_ids_on_disk(paths, registry))
    assemblies = reg.scan_assemblies(paths.targeted, genome_ids)
    for acc, per_genome_files in reg.scan_extractions(paths.targeted, genome_ids).items():
        for genome_id, files in per_genome_files.items():
            if rb.extraction_block(registry, acc, genome_id) is None:
                report.untracked_extractions.append((acc, genome_id))
                reg._infer_extraction(registry, acc, genome_id, files)
                asm_dir = assemblies.get(acc, {}).get(genome_id)
                if asm_dir is not None:
                    reg._infer_assembly(registry, acc, genome_id, asm_dir)
    report.empty_assembly_dirs = reg.empty_assembly_dirs(paths.targeted, genome_ids)
    report.dangling_links = dangling_links(paths.fastq)
    _reference_fill(registry, paths)
    return report


def _reference_fill(registry: reg.Registry, paths: reg.ProjectPaths) -> None:
    for acc in registry.datasets:
        download = rb.download_block(registry, acc) or rb.DownloadBlock()
        if download.state != "downloaded" or (download.complete is not None and download.complete.to_dict()):
            continue
        if download.source == "store":
            complete = _reference_store_verdict(registry, acc)
            if complete is not None:
                # The spot count in the order every caller uses: registry metadata, then the sidecar.
                spots = expected_spots(registry, acc, store=rb.store_block(registry).root)
                if complete["reads_r1"] is not None and spots:
                    complete = verdict_for_count(complete["reads_r1"], spots)
                reg.set_download_verdict(registry, acc, complete)
            continue
        spots = (rb.metadata_block(registry, acc) or rb.MetadataBlock()).run_total_spots
        if not spots:
            continue
        acc_dir = paths.fastq / acc
        if not acc_dir.is_dir():
            continue
        verify = verify_download(acc, acc_dir, spots)
        reg.set_download_verdict(registry, acc, verdict_for_count(verify["reads_r1"], spots))


def _reference_store_verdict(registry: reg.Registry, accession: str) -> Optional[Dict[str, Any]]:
    from prepare.store.layout import sidecar_path, store_paths
    from prepare.store.sidecar import sidecar_completeness

    root = rb.store_block(registry).root
    if not root:
        return None
    try:
        return sidecar_completeness(sidecar_path(store_paths(Path(root)), accession))
    except (OSError, DataAccessError):
        return None


# ------------------------------------------------------------------------- fixtures


def _fastq(path: Path, reads: int = 2, gz: bool = False) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(reads))
    if gz:
        with gzip.open(path, "wt") as handle:
            handle.write(body)
    else:
        path.write_text(body)
    return path


def _project(tmp_path: Path) -> reg.ProjectPaths:
    return reg.ProjectPaths(
        fastq=tmp_path / "fastq",
        metadata=tmp_path / "metadata",
        genomes=tmp_path / "genomes",
        targeted=tmp_path / "targeted",
        matches=tmp_path / "matches",
    )


@pytest.fixture
def fixed_clock(monkeypatch):
    monkeypatch.setattr(reg, "_now", lambda: "2026-09-30T00:00:00+00:00")


def _drifted_project(tmp_path: Path):
    """A registry and a disk that disagree in every way reconcile looks at."""
    from prepare.store.layout import init_store, sidecar_path
    from prepare.store.link import link_dataset
    from prepare.store.sidecar import Sidecar, write_sidecar

    paths = _project(tmp_path)
    for acc in ("SRR1", "SRR2", "SRR3", "SRR4"):
        _fastq(paths.fastq / acc / f"{acc}_1.fastq")
    paths.metadata.mkdir()
    paths.genomes.mkdir()
    (paths.genomes / "GCF_1.fna").write_text(">c\nACGT\n")
    _fastq(paths.targeted / "SRR1" / "GCF_1_1.fastq.gz", gz=True)
    _fastq(paths.targeted / "SRR1" / "GCF_1_2.fastq.gz", gz=True)
    asm = paths.targeted / "SRR1" / "GCF_1_assembly"
    asm.mkdir(parents=True)
    (asm / "final.contigs.fa").write_text(">k len=5\nACGTA\n")
    (paths.targeted / "SRR2" / "GCF_1_assembly").mkdir(parents=True)
    paths.matches.mkdir()
    (paths.matches / "GCF_1.csv").write_text("acc,containment,cANI\nSRR1,0.9,0.99\nSRR2,0.3,0.9\n")

    registry = reg.bootstrap_from_disk(paths, target_path=tmp_path / reg.REGISTRY_FILENAME)
    # Drift after the bootstrap: a download gone, a new one, an extraction nobody recorded.
    (paths.fastq / "SRR2" / "SRR2_1.fastq").unlink()
    _fastq(paths.fastq / "SRR7" / "SRR7_1.fastq")
    registry.datasets["SRR1"].pop("extractions")
    _fastq(paths.targeted / "SRR3" / "GCF_1_1.fastq", reads=5)
    # A download with NCBI's spot count and no verdict yet, one with a spot count that does not
    # match its reads, and one linked from a store whose sidecar holds a verdict.
    reg.record_metadata(registry, "SRR3", paths.metadata / "SRR3_metadata.xml", {"run_total_spots": 2})
    reg.record_metadata(registry, "SRR4", paths.metadata / "SRR4_metadata.xml", {"run_total_spots": 100})
    store = init_store(tmp_path / "store")
    _fastq(store.sra / "SRR5" / "SRR5_1.fastq")
    write_sidecar(
        sidecar_path(store, "SRR5"),
        Sidecar(
            accession="SRR5",
            state="complete",
            reads_per_mate=2,
            ncbi={"spots": 2},
            completeness={"method": "spots", "ratio": 1.0, "verdict": "complete"},
        ),
    )
    link_dataset(paths.fastq, "SRR5", store)
    registry.store = {"root": str(tmp_path / "store")}
    reg.record_download(registry, "SRR5", "downloaded", paths.fastq, attempt=False, source="store", store_name="SRR5")
    return registry, paths


def _dump(registry: reg.Registry) -> str:
    return json.dumps({"datasets": registry.datasets, "genomes": registry.genomes}, sort_keys=True)


# ---------------------------------------------------------------------------- tests


class TestScanApplyMatchesTheOldReconcile:
    def test_apply_of_scan_equals_the_single_step_reconcile(self, tmp_path, fixed_clock):
        registry, paths = _drifted_project(tmp_path)
        old, new = copy.deepcopy(registry), copy.deepcopy(registry)

        old_report = _reference_reconcile(old, paths)
        new_report = apply_reconcile(new, scan_reconcile(new, paths))

        assert new_report == old_report
        assert _dump(new) == _dump(old)
        # The fixture exercises every branch, so the comparison above is not vacuous.
        assert old_report.recorded_missing == ["SRR2"]
        assert old_report.untracked_fastq == ["SRR7"]
        assert old_report.untracked_extractions == [("SRR1", "GCF_1"), ("SRR3", "GCF_1")]
        assert old_report.empty_assembly_dirs == [("SRR2", "GCF_1")]
        assert rb.download_verdict(old, "SRR3").verdict == "complete"
        assert rb.download_verdict(old, "SRR4").verdict == "truncated"
        assert rb.download_verdict(old, "SRR5").verdict == "complete"
        assert rb.extraction_block(old, "SRR1", "GCF_1").assembly.contigs == 1

    def test_reconcile_from_the_registry_module_is_the_two_step_form(self, tmp_path, fixed_clock):
        registry, paths = _drifted_project(tmp_path)
        old, new = copy.deepcopy(registry), copy.deepcopy(registry)

        assert reg.reconcile(new, paths) == _reference_reconcile(old, paths)
        assert _dump(new) == _dump(old)
        assert reconcile is not reg.reconcile  # the registry module keeps a thin wrapper

    def test_scan_leaves_the_snapshot_unchanged(self, tmp_path, fixed_clock):
        registry, paths = _drifted_project(tmp_path)
        before = _dump(registry)

        scan_reconcile(registry, paths)

        assert _dump(registry) == before


class TestApplyRechecksAgainstTheRegistryItIsGiven:
    def test_apply_reads_no_fastq(self, tmp_path, fixed_clock, monkeypatch):
        import prepare.data.registry_reconcile as rr

        registry, paths = _drifted_project(tmp_path)
        plan = scan_reconcile(registry, paths)

        def refuse(*_args, **_kwargs):
            raise AssertionError("apply_reconcile must not read FASTQ")

        monkeypatch.setattr(rr, "count_fastq_reads", refuse)
        monkeypatch.setattr(rr, "verify_download", refuse)
        monkeypatch.setattr(rr, "summarise_contigs", refuse)
        apply_reconcile(registry, plan)

        assert rb.download_verdict(registry, "SRR4").verdict == "truncated"

    def test_changes_committed_after_the_scan_are_kept(self, tmp_path, fixed_clock):
        registry, paths = _drifted_project(tmp_path)
        plan = scan_reconcile(registry, paths)
        # Meanwhile another process records SRR9, re-records SRR2 and records SRR3's extraction.
        fresh = copy.deepcopy(registry)
        reg.record_download(fresh, "SRR9", "failed", paths.fastq, message="timeout")
        reg.record_download(fresh, "SRR2", "failed", paths.fastq)
        reg.record_extraction(fresh, "SRR3", "GCF_1", [], 42, False, {})

        report = apply_reconcile(fresh, plan)

        assert rb.download_block(fresh, "SRR9").state == "failed"
        assert rb.download_block(fresh, "SRR2").state == "failed"
        assert "SRR2" not in report.recorded_missing
        assert rb.extraction_block(fresh, "SRR3", "GCF_1").mapped_reads == 42
        assert ("SRR3", "GCF_1") not in report.untracked_extractions

    def test_a_download_finished_after_the_scan_is_not_marked_missing(self, tmp_path, fixed_clock):
        registry, paths = _drifted_project(tmp_path)
        plan = scan_reconcile(registry, paths)
        # Another process downloads SRR8 and records it once the scan has listed the folders.
        _fastq(paths.fastq / "SRR8" / "SRR8_1.fastq")
        fresh = copy.deepcopy(registry)
        reg.record_download(fresh, "SRR8", "downloaded", paths.fastq)

        report = apply_reconcile(fresh, plan)

        assert rb.download_block(fresh, "SRR8").state == "downloaded"
        assert report.recorded_missing == ["SRR2"]

    def test_an_accession_removed_after_the_scan_is_not_recorded_as_downloaded(self, tmp_path, fixed_clock):
        registry, paths = _drifted_project(tmp_path)
        plan = scan_reconcile(registry, paths)
        assert "SRR7" in plan.on_disk
        # store_unlink (or a user) removes the folder between the scan and the apply.
        (paths.fastq / "SRR7" / "SRR7_1.fastq").unlink()

        report = apply_reconcile(registry, plan)

        assert "SRR7" not in report.untracked_fastq
        assert rb.download_block(registry, "SRR7") is None

    def test_a_store_link_marked_missing_leaves_the_linked_list(self, tmp_path, fixed_clock):
        """An accession on the registry's list of store links whose project entry is gone is
        marked missing and taken off that list, so the store no longer counts it as linked."""
        paths = _project(tmp_path)
        paths.fastq.mkdir()
        registry = reg.Registry()
        registry.store = {"root": str(tmp_path / "store"), "linked": ["SRR1"]}
        reg.record_download(
            registry, "SRR1", "downloaded", paths.fastq, attempt=False, source="store", store_name="SRR1"
        )

        report = apply_reconcile(registry, scan_reconcile(registry, paths))

        assert report.recorded_missing == ["SRR1"]
        assert rb.download_block(registry, "SRR1").state == "missing"
        assert rb.store_block(registry).linked == []

    def test_a_verdict_is_used_only_for_the_spot_count_it_was_computed_for(self, tmp_path, fixed_clock):
        registry, paths = _drifted_project(tmp_path)
        plan = scan_reconcile(registry, paths)
        fresh = copy.deepcopy(registry)
        reg.record_metadata(fresh, "SRR4", paths.metadata / "SRR4_metadata.xml", {"run_total_spots": 2})

        apply_reconcile(fresh, plan)

        assert rb.download_verdict(fresh, "SRR4") is None
        assert rb.download_verdict(fresh, "SRR3").verdict == "complete"

    def test_a_verdict_is_not_used_for_a_download_redone_after_the_scan(self, tmp_path, monkeypatch):
        monkeypatch.setattr(reg, "_now", lambda: "2026-09-30T00:00:00+00:00")
        registry, paths = _drifted_project(tmp_path)
        plan = scan_reconcile(registry, paths)
        # SRR3 is downloaded again (same spot count, no verdict recorded) before the apply.
        monkeypatch.setattr(reg, "_now", lambda: "2026-09-30T01:00:00+00:00")
        reg.record_download(registry, "SRR3", "downloaded", paths.fastq, "Downloaded 1 files, unverified")

        apply_reconcile(registry, plan)

        assert rb.download_verdict(registry, "SRR3") is None
        assert rb.download_verdict(registry, "SRR4").verdict == "truncated"

    def test_an_extraction_that_appeared_after_the_scan_is_left_for_the_next_run(self, tmp_path, fixed_clock):
        registry, paths = _drifted_project(tmp_path)
        plan = scan_reconcile(registry, paths)
        registry.datasets["SRR3"]["extractions"] = {"GCF_1": {"mapped_reads": 5}}
        plan_with_it_recorded = scan_reconcile(registry, paths)
        registry.datasets["SRR3"].pop("extractions")

        report = apply_reconcile(registry, plan_with_it_recorded)

        assert ("SRR3", "GCF_1") not in report.untracked_extractions
        assert rb.extraction_block(registry, "SRR3", "GCF_1") is None
        assert ("SRR3", "GCF_1") in apply_reconcile(copy.deepcopy(registry), plan).untracked_extractions


class TestRegistryUpdate:
    def test_applies_the_mutation_under_the_lock_and_returns_its_result(self, tmp_path):
        path = tmp_path / reg.REGISTRY_FILENAME
        with reg.registry_transaction(path) as r:
            reg.record_exclusion(r, "SRR1", "contaminated")

        def mutation(r):
            assert (tmp_path / f"{reg.REGISTRY_FILENAME}.lock").exists()
            reg.record_exclusion(r, "SRR2", "low depth")
            return len(r.datasets)

        assert registry_update(path, mutation) == 2
        assert set(reg.query(reg.load_registry(path), "excluded")) == {"SRR1", "SRR2"}
        assert not (tmp_path / f"{reg.REGISTRY_FILENAME}.lock").exists()

    def test_a_mutation_that_raises_writes_nothing(self, tmp_path):
        path = tmp_path / reg.REGISTRY_FILENAME
        with reg.registry_transaction(path) as r:
            reg.record_exclusion(r, "SRR1", "contaminated")
        before = path.read_text()

        def mutation(r):
            reg.record_exclusion(r, "SRR2", "low depth")
            raise DataAccessError("refused")

        with pytest.raises(DataAccessError):
            registry_update(path, mutation)
        assert path.read_text() == before
