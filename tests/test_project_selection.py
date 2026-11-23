"""Tests for prepare.data.project_selection: which recorded selection applies to a genome."""

from prepare.data import registry_blocks as rb
from prepare.data.project_selection import selection_for_genome
from prepare.data.registry import Registry


def _select(registry, accession, criteria):
    block = rb.SelectionBlock(selected=True, date="2026-11-23", criteria=criteria, output="accessions.txt")
    rb.set_selection_block(registry, accession, block)


def test_mixed_criteria_are_judged_per_accession(tmp_path):
    """Selections made for different genomes coexist; each genome gets only its own, on every run."""
    registry = Registry(path=tmp_path / "prepare_registry.json")
    _select(registry, "SRR1", {"column": "G1", "threshold": 0.5})
    _select(registry, "SRR3", {"column": "G2", "threshold": 0.5})
    _select(registry, "SRR4", {"column": "G2", "threshold": 0.5})

    assert selection_for_genome(registry, "G1") == {"SRR1"}
    assert selection_for_genome(registry, "G2") == {"SRR3", "SRR4"}


def test_max_containment_and_genome_ids_cover_a_genome(tmp_path):
    registry = Registry(path=tmp_path / "prepare_registry.json")
    _select(registry, "SRR1", {"column": "max_containment"})
    _select(registry, "SRR2", {"column": "A+B", "genome_ids": ["A", "B"]})
    _select(registry, "SRR3", {})

    assert selection_for_genome(registry, "B") == {"SRR1", "SRR2", "SRR3"}
    assert selection_for_genome(registry, "C") == {"SRR1", "SRR3"}


def test_no_selection_is_an_empty_set(tmp_path):
    assert selection_for_genome(Registry(path=tmp_path / "prepare_registry.json"), "G1") == set()


def test_selection_threshold_is_the_lowest_covering_threshold(tmp_path):
    from prepare.data.project_selection import selection_threshold

    registry = Registry(path=tmp_path / "prepare_registry.json")
    _select(registry, "SRR1", {"column": "G", "threshold": 0.05})
    _select(registry, "SRR2", {"column": "G", "threshold": 0.2})
    _select(registry, "SRR3", {"column": "H", "threshold": 0.01})

    assert selection_threshold(registry, "G") == 0.05


def test_threshold_without_a_recorded_selection_threshold_is_the_default(tmp_path):
    from prepare.data.project_selection import selection_threshold

    registry = Registry(path=tmp_path / "prepare_registry.json")
    _select(registry, "SRR1", {"column": "G"})

    assert selection_threshold(registry, "G") is None
    assert selection_threshold(Registry(path=tmp_path / "other.json"), "G") is None
