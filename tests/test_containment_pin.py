"""The containment summary and the screening registry, pinned on a 200 x 4 table.

The expected outputs in ``tests/fixtures/containment_pin`` were written by the row-by-row
implementation before the numpy one replaced it (2026-09-26). The table has ties in
``max_containment``, zeros, and one genome that is 0 everywhere; the registry it is recorded
into already holds screening entries, an inferred block, an exclusion-only record and an empty
record, so the merge with earlier records and the per-genome cap are both covered.
"""

import itertools
import json
from pathlib import Path

import pandas as pd
import pytest

from prepare.data import registry as reg
from prepare.data.branchwater import _generate_containment_summary
from tests.perf_fixtures import containment_data

PIN = Path(__file__).parent / "fixtures" / "containment_pin"


def seeded_registry(path):
    """A registry with the earlier records the run has to merge with, all dated ``EARLIER``."""
    r = reg.load_registry(path)
    reg.record_screening(r, "SRR1000003", "GCF_000", 0.99, 0.98, "branchwater", 0.1, None)
    reg.record_screening(r, "SRR9999999", "GCF_001", 0.4, None, "matches", 0.0, None)
    reg.record_screening(r, "SRR9999998", "GCF_OTHER", 0.7, None, "matches", 0.0, None)
    r.datasets["SRR9999998"]["screening"]["inferred"] = True
    reg.record_screening(r, "SRR1000010", "GCF_OTHER", 0.6, None, "matches", 0.0, None)
    r.datasets["SRR1000010"]["screening"]["inferred"] = True
    reg.record_exclusion(r, "SRR1000020", "16S amplicon")
    r.datasets["SRR1000030"] = {}
    for record in r.datasets.values():
        for block in ("screening", "exclusion"):
            if block in record:
                record[block]["date"] = "EARLIER"
    return r


def _summary_fields(summary):
    return {
        "thresholds": summary.thresholds,
        "counts": summary.counts,
        "max_containment": summary.max_containment,
        "genome_to_samples": summary.genome_to_samples,
        "sample_to_genomes": summary.sample_to_genomes,
    }


def _new_dates(datasets):
    """Every screening date the run wrote (the seeded ones read ``EARLIER``)."""
    return [
        record["screening"]["date"]
        for record in datasets.values()
        if "screening" in record and record["screening"]["date"] != "EARLIER"
    ]


def _without_new_dates(datasets):
    copy = json.loads(json.dumps(datasets))
    for record in copy.values():
        if "screening" in record and record["screening"]["date"] != "EARLIER":
            record["screening"]["date"] = "RUN"
    return copy


@pytest.fixture
def summary(tmp_path):
    """The summary of the 200 x 4 table at step 0.001; the parsed table is in ``tmp_path``."""
    return _generate_containment_summary(
        containment_data(200, 4), tmp_path / "parsed.txt", tmp_path / "summary.txt", 0.001
    )


@pytest.fixture
def ticking_clock(monkeypatch):
    """Make every registry timestamp distinct, so a run that stamps each cell shows it."""
    counter = itertools.count()
    monkeypatch.setattr(reg, "_now", lambda: f"T{next(counter):06d}")


def test_summary_matches_the_pinned_output(summary):
    expected = json.loads((PIN / "summary_200x4.json").read_text())
    assert _summary_fields(summary) == expected
    assert all(type(count) is int for count in summary.counts)


def test_summary_counts_at_a_coarse_step(tmp_path):
    """Step 0.1 gives 11 thresholds, 1.0 down to 0.0; counts are the rows at or above each."""
    data = containment_data(200, 4)
    result = _generate_containment_summary(data, tmp_path / "p.txt", tmp_path / "s.txt", 0.1)
    maxima = [max(row.values()) for row in data.values()]
    assert result.thresholds == [round(i * 0.1, 2) for i in range(10, -1, -1)]
    assert result.counts == [sum(m >= i * 0.1 for m in maxima) for i in range(10, -1, -1)]


@pytest.mark.parametrize("cap", ["5000", "40"])
@pytest.mark.parametrize("source", ["path", "table"])
def test_screening_matches_the_pinned_registry(summary, tmp_path, cap, source, ticking_clock):
    expected = json.loads((PIN / "screening_200x4.json").read_text())[cap]
    r = seeded_registry(tmp_path / "prepare_registry.json")
    table = tmp_path / "parsed.txt" if source == "path" else summary.table
    recorded = reg.record_screening_from_table(r, table, "matches", max_screened=int(cap))
    assert recorded == expected["recorded"]
    assert list(r.datasets) == expected["order"]
    assert list(r.genomes) == expected["genomes"]
    assert _without_new_dates(r.datasets) == _without_new_dates(expected["datasets"])
    dates = _new_dates(r.datasets)
    assert dates and len(set(dates)) == 1, "one timestamp per run"


def test_screening_from_the_in_memory_table_reads_no_file(summary, tmp_path, monkeypatch):
    def refuse(*args, **kwargs):
        raise AssertionError("record_screening_from_table re-read the parsed table")

    monkeypatch.setattr(pd, "read_csv", refuse)
    r = reg.load_registry(tmp_path / "prepare_registry.json")
    assert isinstance(summary.table, pd.DataFrame)
    assert reg.record_screening_from_table(r, summary.table, "matches") > 0


def test_screening_leaves_nothing_for_cap_screening_to_trim(summary, tmp_path):
    r = seeded_registry(tmp_path / "prepare_registry.json")
    reg.record_screening_from_table(r, summary.table, "matches", max_screened=40)
    for genome in ("GCF_000", "GCF_001", "GCF_002", "GCF_003"):
        assert reg.cap_screening(r, genome, 40) == 0


def test_parse_containment_command_records_the_in_memory_table(summary, tmp_path, monkeypatch):
    """The command hands the registry the table it just built instead of reading the file back."""
    import argparse

    from prepare.cli.commands import containment as command_module

    monkeypatch.setattr(command_module, "parse_containment_data", lambda *args, **kwargs: summary)
    seen = []
    real = command_module.record_screening_from_table
    monkeypatch.setattr(
        command_module,
        "record_screening_from_table",
        lambda registry, table, *args, **kwargs: seen.append(table) or real(registry, table, *args, **kwargs),
    )
    args = argparse.Namespace(
        matches_folder=str(tmp_path / "matches"),
        parsed_containment_file=str(tmp_path / "parsed.txt"),
        summary_containment_file=str(tmp_path / "summary.txt"),
        step_size=0.001,
        details_file=None,
        registry=str(tmp_path / "prepare_registry.json"),
        registry_max_screened=40,
    )
    assert command_module.ParseContainmentCommand().execute(args) == 0
    assert len(seen) == 1 and seen[0] is summary.table
    written = reg.load_registry(tmp_path / "prepare_registry.json")
    assert all(reg.cap_screening(written, genome, 40) == 0 for genome in ("GCF_000", "GCF_001", "GCF_002"))
