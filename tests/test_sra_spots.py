"""Tests for prepare.data.sra.spots: one lookup for an accession's expected spot count, the
verdict a read count gives against it, and the merge rule that never turns a recorded ``complete``
or ``truncated`` verdict into ``unverified``."""

import json

import pytest

from prepare.data import registry_blocks as rb
from prepare.data.registry import Registry
from prepare.data.sra import expected_spots, merged_verdict, spots_from_xml, verdict_for_count, verify_download
from prepare.store.layout import sidecar_path, store_paths

from tests.test_store_sidecar import NCBI_XML

ACC = "SRR1"


def _registry(metadata=None, verdict=None):
    """An in-memory registry with ``ACC``'s metadata block and download verdict as given."""
    registry = Registry()
    dataset = {}
    if metadata is not None:
        dataset["metadata"] = metadata
    if verdict is not None:
        dataset["download"] = {"attempts": 1, "state": "complete", "complete": verdict}
    if dataset:
        registry.datasets[ACC] = dataset
    return registry


def _store_with_spots(tmp_path, spots):
    """A store root on ``tmp_path`` holding a sidecar for ``ACC`` whose ``ncbi.spots`` is ``spots``."""
    paths = store_paths(tmp_path / "store")
    path = sidecar_path(paths, ACC)
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({"accession": ACC, "ncbi": {"spots": spots}, "schema": 1}))
    return paths


def _xml_folder(tmp_path, name="metadata", total_spots="10"):
    """A folder holding ``ACC``'s metadata XML with ``total_spots`` as its run's spot count."""
    folder = tmp_path / name
    folder.mkdir()
    (folder / f"{ACC}_metadata.xml").write_text(NCBI_XML.replace('total_spots="10"', f'total_spots="{total_spots}"'))
    return folder


# ------------------------------------------------------------------ verdict_for_count


def test_verdict_for_count_exact_threshold_is_complete():
    result = verdict_for_count(99, 100)
    assert result == {"method": "spots", "ratio": 0.99, "verdict": "complete", "expected_spots": 100, "reads_r1": 99}


def test_verdict_for_count_just_below_threshold_is_truncated():
    result = verdict_for_count(9899, 10000)
    assert result["verdict"] == "truncated"
    assert result["ratio"] == 0.9899


def test_verdict_for_count_without_spot_count_is_unverified():
    for expected in (None, 0):
        result = verdict_for_count(50, expected)
        assert result["verdict"] == "unverified"
        assert result["method"] == "unverified"
        assert result["ratio"] is None
        assert result["reads_r1"] == 50


def test_verdict_for_count_without_read_count_is_unverified():
    result = verdict_for_count(None, 100)
    assert result["verdict"] == "unverified"
    assert result["expected_spots"] == 100
    assert result["reads_r1"] is None


@pytest.mark.parametrize(
    "records, expected, verdict, ratio",
    [(99, 100, "complete", 0.99), (98, 100, "truncated", 0.98), (5, None, "unverified", None)],
)
def test_verify_download_output_unchanged(tmp_path, records, expected, verdict, ratio):
    acc_dir = tmp_path / ACC
    acc_dir.mkdir()
    (acc_dir / f"{ACC}.fastq").write_text("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(records)))
    result = verify_download(ACC, acc_dir, expected)
    assert list(result) == ["reads_r1", "expected_spots", "ratio", "verdict", "bytes_total"]
    assert result["reads_r1"] == records
    assert result["expected_spots"] == expected
    assert result["ratio"] == ratio
    assert result["verdict"] == verdict


# ------------------------------------------------------------------ spots_from_xml


def test_spots_from_xml_reads_run_total_spots(tmp_path):
    folder = _xml_folder(tmp_path, total_spots="1234")
    assert spots_from_xml(ACC, [folder]) == 1234


def test_spots_from_xml_skips_folders_without_the_file(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    folder = _xml_folder(tmp_path, total_spots="77")
    assert spots_from_xml(ACC, [empty, str(folder)]) == 77


def test_spots_from_xml_parse_error_returns_none(tmp_path):
    folder = tmp_path / "metadata"
    folder.mkdir()
    (folder / f"{ACC}_metadata.xml").write_text("<EXPERIMENT_PACKAGE_SET><unclosed>")
    assert spots_from_xml(ACC, [folder]) is None


def test_spots_from_xml_no_folders_returns_none(tmp_path):
    assert spots_from_xml(ACC, []) is None
    assert spots_from_xml(ACC, [tmp_path / "missing"]) is None


# ------------------------------------------------------------------ expected_spots


def test_expected_spots_prefers_registry_metadata(tmp_path):
    registry = _registry(metadata={"run_total_spots": 500}, verdict={"verdict": "truncated", "expected_spots": 400})
    store = _store_with_spots(tmp_path, 300)
    folder = _xml_folder(tmp_path, total_spots="200")
    assert expected_spots(registry, ACC, store=store, xml_folders=[folder]) == 500


def test_expected_spots_accepts_a_numeric_string_in_metadata():
    assert expected_spots(_registry(metadata={"run_total_spots": "600"}), ACC) == 600


def test_expected_spots_sidecar_wins_over_previous_verdict(tmp_path):
    registry = _registry(metadata={"run_total_spots": None}, verdict={"verdict": "truncated", "expected_spots": 400})
    store = _store_with_spots(tmp_path, 300)
    folder = _xml_folder(tmp_path, total_spots="200")
    assert expected_spots(registry, ACC, store=store, xml_folders=[folder]) == 300


def test_expected_spots_xml_wins_over_previous_verdict(tmp_path):
    registry = _registry(verdict={"verdict": "truncated", "expected_spots": 400})
    folder = _xml_folder(tmp_path, total_spots="200")
    assert expected_spots(registry, ACC, store=store_paths(tmp_path / "store"), xml_folders=[folder]) == 200


def test_expected_spots_previous_verdict_is_the_last_resort(tmp_path):
    registry = _registry(metadata={"run_total_spots": None}, verdict={"verdict": "truncated", "expected_spots": 400})
    empty = tmp_path / "empty"
    empty.mkdir()
    assert expected_spots(registry, ACC, store=store_paths(tmp_path / "store"), xml_folders=[empty]) == 400


def test_expected_spots_falls_back_to_sidecar(tmp_path):
    registry = _registry(verdict={"verdict": "unverified"})
    store = _store_with_spots(tmp_path, 300)
    folder = _xml_folder(tmp_path, total_spots="200")
    assert expected_spots(registry, ACC, store=store, xml_folders=[folder]) == 300


def test_expected_spots_accepts_a_store_root_path(tmp_path):
    store = _store_with_spots(tmp_path, 300)
    assert expected_spots(Registry(), ACC, store=store.root) == 300


def test_expected_spots_falls_back_to_xml(tmp_path):
    store = store_paths(tmp_path / "store")  # no sidecar written
    folder = _xml_folder(tmp_path, total_spots="200")
    assert expected_spots(Registry(), ACC, store=store, xml_folders=[folder]) == 200


def test_expected_spots_unknown_everywhere_is_none(tmp_path):
    assert expected_spots(Registry(), ACC) is None
    assert expected_spots(None, ACC, store=store_paths(tmp_path / "store"), xml_folders=[tmp_path]) is None


def test_expected_spots_ignores_an_unreadable_sidecar(tmp_path):
    paths = store_paths(tmp_path / "store")
    path = sidecar_path(paths, ACC)
    path.parent.mkdir(parents=True)
    path.write_text("[1, 2]")
    folder = _xml_folder(tmp_path, total_spots="200")
    assert expected_spots(Registry(), ACC, store=paths, xml_folders=[folder]) == 200


# ------------------------------------------------------------------ merged_verdict


TRUNCATED = {"method": "spots", "ratio": 0.5, "verdict": "truncated", "expected_spots": 100, "reads_r1": 50}
COMPLETE = {"method": "spots", "ratio": 1.0, "verdict": "complete", "expected_spots": 100, "reads_r1": 100}
UNVERIFIED = {"method": "unverified", "ratio": None, "verdict": "unverified", "expected_spots": None, "reads_r1": 50}


def test_merged_verdict_never_downgrades_truncated_to_unverified():
    assert merged_verdict(TRUNCATED, UNVERIFIED, None, None) == TRUNCATED
    assert merged_verdict(TRUNCATED, None, None, None) == TRUNCATED


def test_merged_verdict_recomputes_when_the_count_becomes_known():
    result = merged_verdict(TRUNCATED, UNVERIFIED, 100, 100)
    assert result["verdict"] == "complete"
    assert result["ratio"] == 1.0
    result = merged_verdict(None, UNVERIFIED, 50, 100)
    assert result["verdict"] == "truncated"


def test_merged_verdict_keeps_a_new_definite_verdict():
    assert merged_verdict(TRUNCATED, COMPLETE, None, None) == COMPLETE
    assert merged_verdict(COMPLETE, TRUNCATED, None, None) == TRUNCATED


def test_merged_verdict_never_downgrades_complete_to_unverified():
    assert merged_verdict(COMPLETE, UNVERIFIED, None, None) == COMPLETE
    assert merged_verdict(COMPLETE, None, None, None) == COMPLETE


def test_merged_verdict_relink_without_read_count_keeps_complete():
    """A relink to a store copy whose sidecar records no read count keeps the recorded verdict."""
    relinked = {"verdict": "unverified", "ratio": None, "expected_spots": None, "reads_r1": None}
    assert merged_verdict(COMPLETE, relinked, None, None) == COMPLETE


def test_merged_verdict_recompute_wins_over_a_recorded_complete():
    result = merged_verdict(COMPLETE, UNVERIFIED, 50, 100)
    assert result["verdict"] == "truncated"
    assert result["ratio"] == 0.5


def test_merged_verdict_unverified_over_unverified_takes_new():
    newer = dict(UNVERIFIED, reads_r1=60)
    assert merged_verdict(UNVERIFIED, newer, None, None) == newer


def test_merged_verdict_accepts_registry_verdict_blocks():
    previous = rb.Verdict(method="spots", ratio=0.5, verdict="truncated", expected_spots=100, reads_r1=50)
    result = merged_verdict(previous, rb.Verdict(verdict="unverified"), None, None)
    assert isinstance(result, dict)
    assert result == TRUNCATED
    assert merged_verdict(None, previous, None, None) == TRUNCATED


def test_merged_verdict_rejects_a_value_that_is_not_a_verdict():
    with pytest.raises(TypeError):
        merged_verdict("complete", None, None, None)


def test_merged_verdict_nothing_known_is_none():
    assert merged_verdict(None, None, None, None) is None
    assert merged_verdict(None, None, 10, None) is None


def test_merged_verdict_keeps_previous_when_nothing_new():
    assert merged_verdict(COMPLETE, None, None, None) == COMPLETE


def test_merged_verdict_returns_a_copy():
    result = merged_verdict(TRUNCATED, None, None, None)
    result["verdict"] = "changed"
    assert TRUNCATED["verdict"] == "truncated"
