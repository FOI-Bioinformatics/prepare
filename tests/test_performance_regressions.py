"""Performance regression tests with time and memory bounds.

Every test here guards one of the performance fixes made across this project (registry batching
and serialisation, the FASTQ sampler, the containment summary and screening step, metadata XML
parsing, and the status/results report builders); each replaces a slower implementation that a
docstring below describes. Bounds are built the same way throughout: run the fixed code, take the
best of several timed runs on the development machine (Apple silicon laptop, macOS, 2026-09-26)
under ``make test``'s coverage instrumentation, and set the bound to about three times that
measurement, so normal machine and coverage variance does not make the test flaky while an actual
regression (an accidentally reintroduced quadratic pass, for example) still fails it. A failing
test reports the measured time in its assertion message, so a slower machine can be told apart
from a real regression by comparing that figure with the one in the docstring.

Fixtures live in ``tests/perf_fixtures.py``: a synthetic 20,000-dataset project registry built
through the real ``record_*`` writers, a 100,000-record gzip FASTQ file, 300 synthetic NCBI efetch
metadata XML files, and a 20,000 x 5 containment table.

Every time bound (not the memory bound) is multiplied by the ``PREPARE_PERF_SCALE`` environment
variable (default 1, at least 1; CI sets 4 for its shared runners), and each failure message prints the
scale in effect.

Run only these tests with ``python -m pytest -m perf`` (``make test-perf``).
"""

import argparse
import os
import shutil
import sys
import time
import tracemalloc
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from prepare.cli.commands.sra import DownloadSraCommand
from prepare.data import registry as registry_mod
from prepare.data import registry_blocks as rb
from prepare.data.branchwater import _generate_containment_summary
from prepare.data.metadata import parse_metadata
from prepare.data.sra import sample_records
from prepare.processing.results import results_rows
from prepare.processing.status_report import build_report
from tests.perf_fixtures import (
    build_registry,
    containment_data,
    write_illumina_like_fastq_gz,
    write_metadata_folder,
)


def _read_scale(text):
    """The perf time-bound factor parsed from ``PREPARE_PERF_SCALE``; values below 1 are rejected."""
    try:
        scale = float(text)
    except ValueError:
        raise ValueError(f"PREPARE_PERF_SCALE must be a number of at least 1, got {text!r}") from None
    if not scale >= 1:
        raise ValueError(f"PREPARE_PERF_SCALE must be at least 1, got {text!r}")
    return scale


PERF_SCALE = _read_scale(os.environ.get("PREPARE_PERF_SCALE", "1"))


def _scaled(bound):
    """``bound`` (seconds) multiplied by the module's perf scale factor."""
    return bound * PERF_SCALE


def _over(label, elapsed, bound):
    """Failure message naming the measured time, the scaled bound and the scale factor."""
    return f"{label} took {elapsed:.2f} s (bound {_scaled(bound):.2f} s, PREPARE_PERF_SCALE={PERF_SCALE:g})"


@pytest.fixture(scope="module")
def registry_20k(tmp_path_factory):
    """A 20,000-dataset registry file and the accession lists ``build_registry`` wrote."""
    root = tmp_path_factory.mktemp("perf_registry")
    path = root / registry_mod.REGISTRY_FILENAME
    names = build_registry(path)
    return path, names


@pytest.fixture(scope="module")
def containment_20000x5():
    """20,000 accessions by 5 genomes, about half the cells positive, the last genome all zeros."""
    return containment_data(20000, 5)


@pytest.fixture(scope="module")
def metadata_folder_300(tmp_path_factory):
    """300 synthetic metadata XML files, 40 sample attributes each, from overlapping tag sets."""
    folder = tmp_path_factory.mktemp("perf_metadata") / "metadata"
    write_metadata_folder(folder, count=300, per_file=40)
    return folder


@pytest.fixture(scope="module")
def fastq_100k(tmp_path_factory):
    """A 100,000-record gzip FASTQ file of 150 bp reads."""
    return write_illumina_like_fastq_gz(tmp_path_factory.mktemp("perf_fastq") / "SRR0_1.fastq.gz", 100_000)


def _best_of(call, repeats=3):
    """The shortest of ``repeats`` timed runs of ``call``, with the result of the last run."""
    timings = []
    result = None
    for _ in range(repeats):
        start = time.perf_counter()
        result = call()
        timings.append(time.perf_counter() - start)
    return min(timings), result


@pytest.mark.perf
def test_record_run_outcomes_with_5000_skipped_is_one_fast_transaction(registry_20k, monkeypatch, tmp_path):
    """5,000 skipped accessions are recorded in one registry write.

    Measured 2026-09-26 under coverage on an Apple-silicon laptop (macOS): 0.34 s after the fix
    (about 0.4 s per accession, over 30 minutes in total, before it, when every accession wrote its
    own transaction). Bound: 1.1 s, about three times the measured time, the best of three runs,
    each on its own copy of the registry.
    """
    source, names = registry_20k
    skipped = names["all"][-5000:]
    writes = []
    real_write = registry_mod._write_registry

    def counting_write(registry, target, lock):
        writes.append(target)
        return real_write(registry, target, lock)

    # Each run gets a fresh copy, so every run records the 5,000 outcomes rather than finding them done.
    copies = iter(range(3))

    def one_run():
        path = tmp_path / f"run_{next(copies)}" / registry_mod.REGISTRY_FILENAME
        path.parent.mkdir()
        shutil.copy(source, path)
        args = argparse.Namespace(registry=str(path))
        start = time.perf_counter()
        DownloadSraCommand()._record_run_outcomes(args, {"skipped_accessions": skipped}, path.parent / "fastq")
        return time.perf_counter() - start, path

    monkeypatch.setattr(registry_mod, "_write_registry", counting_write)
    runs = [one_run() for _ in range(3)]
    monkeypatch.undo()
    elapsed = min(t for t, _ in runs)
    assert len(writes) == 3
    assert rb.download_block(registry_mod.load_registry(runs[-1][1]), skipped[-1]).state == "skipped"
    assert elapsed < _scaled(1.1), _over("_record_run_outcomes for 5,000 skipped accessions", elapsed, 1.1)


@pytest.mark.perf
def test_write_registry_on_20000_datasets_is_fast(registry_20k, tmp_path):
    """Serialising and writing the 20,000-dataset registry.

    Measured 2026-09-26 under coverage on an Apple-silicon laptop (macOS): 0.08 s compact (0.33 s
    with the earlier indented output). Bound: 0.3 s, about three times the measured time, the best
    of three writes so a single slow disk flush does not fail the test.
    """
    # A copy of its own, so the write never touches the module-scoped registry the other tests read.
    path = Path(shutil.copy(registry_20k[0], tmp_path / registry_mod.REGISTRY_FILENAME))
    registry = registry_mod.load_registry(path)
    with registry_mod._acquire_lock(path.with_name(path.name + ".lock")) as lock:
        best, _ = _best_of(lambda: registry_mod._write_registry(registry, path, lock))
    assert "\n  " not in path.read_text()[:10000]
    assert best < _scaled(0.3), _over("_write_registry on 20,000 datasets", best, 0.3)


@pytest.mark.perf
def test_status_report_on_20000_datasets_is_fast(registry_20k):
    """``build_report``, the dict ``status --json`` prints, on the 20,000-dataset registry.

    Measured 2026-09-26 under coverage on an Apple-silicon laptop (macOS): 0.18 s after the fix
    (0.55 s before it, when the stages were counted with one pass over the registry per stage and
    genome and then listed again). Bound: 1.0 s, the best of three runs (the 3x figure would be
    0.6 s; held at the 1.0 s floor the task-9 ruling set for this test instead).
    """
    path, names = registry_20k
    root = path.parent
    registry = registry_mod.load_registry(path)
    args = argparse.Namespace(
        fastq_folder=str(root / "fastq"),
        metadata_folder=str(root / "metadata"),
        genomes_folder=str(root / "genomes"),
        accessions_file=None,
        parsed_containment=None,
        data_root=None,
        genome=None,
        init=False,
    )
    paths = registry_mod.ProjectPaths(fastq=root / "fastq", targeted=root / "targeted")
    best, report = _best_of(lambda: build_report(registry, args, paths, path, True))
    assert report["wanted"]["total"] == len(names["selected"])
    assert report["stages"]["screened"]["count"] == len(names["all"])
    assert best < _scaled(1.0), _over("build_report on 20,000 datasets", best, 1.0)


@pytest.mark.perf
def test_results_rows_on_20000_datasets_by_3_genomes_is_fast(registry_20k):
    """``results_rows`` for 60,000 (accession, genome) pairs, with a parsed containment table of 20,000 rows.

    Measured 2026-09-26 under coverage on an Apple-silicon laptop (macOS): 0.38 s after the fix
    (1.80 s before it, when every pair converted about six registry blocks and every table row
    became a pandas Series). Bound: 1.5 s, the best of three runs (the 3x figure would be 1.2 s;
    held at the 1.5 s floor the task-9 ruling set for this test instead).
    """
    path, names = registry_20k
    registry = registry_mod.load_registry(path)
    values = np.random.default_rng(0).random((len(names["all"]), 3))
    table = pd.DataFrame(values, index=names["all"], columns=["G1", "G2", "G3"])
    table["max_containment"] = table.max(axis=1)
    best, rows = _best_of(lambda: results_rows(registry, parsed_table=table))
    assert len(rows) == 3 * len(names["all"])
    assert best < _scaled(1.5), _over("results_rows for 60,000 pairs", best, 1.5)


@pytest.mark.perf
def test_stage_counts_on_5000_datasets_is_fast(tmp_path):
    """``stage_counts`` on 5,000 datasets, reading single registry fields rather than typed blocks.

    Measured 2026-09-26 under coverage on an Apple-silicon laptop (macOS): 0.05 s, the best of
    three runs (about 8x slower before the fix, when every dataset was converted into a typed
    ``RegistryBlock``). Bound: 0.5 s, well above 3x the measured time, to hold up against the
    slowdown a busier or shared machine can add (a fix-round-1 review measured this same call at a
    stable 0.13-0.14 s on such a machine; see the task-9 report's fix-round-1 section).
    """
    r = registry_mod.Registry()
    for i in range(5000):
        accession = f"SRR{i:07d}"
        registry_mod.record_screening(r, accession, "G1", 0.5, None, "matches", 0.0, None)
        r.datasets[accession]["download"] = {
            "attempts": 1,
            "state": "downloaded",
            "date": "d",
            "files": [{"path": "p", "bytes": 1, "mtime": "m"}, {"path": "q", "bytes": 1, "mtime": "m"}],
            "bytes_total": 2,
            "message": "",
            "complete": {"verdict": "complete", "ratio": 1.0},
        }
        if i % 10 == 0:
            r.datasets[accession]["extractions"] = {"G1": {"mapped_reads": i % 20, "files": [], "assembly": None}}
    conversions = []
    original = rb.RegistryBlock.from_dict.__func__

    def counting(cls, data):
        conversions.append(cls.__name__)
        return original(cls, data)

    with patch.object(rb.RegistryBlock, "from_dict", classmethod(counting)):
        best, counts = _best_of(lambda: registry_mod.stage_counts(r))
    assert counts["stages"]["downloaded"] == 5000 and counts["genomes"]["G1"]["extracted"] == 250
    assert len(counts["genomes"]["G1"]["zero_mapped"]) == 250
    assert conversions == []
    assert best < _scaled(0.5), _over("stage_counts on 5000 datasets", best, 0.5)


@pytest.mark.perf
def test_sample_records_100k_to_10k_is_fast(fastq_100k):
    """Sampling 10,000 of 100,000 gzip records, record count included.

    Measured 2026-09-26 under coverage on an Apple-silicon laptop (macOS): 0.05 s; the reservoir
    sampler this replaced took 0.11 s here (from the original Task 2 measurement) and scaled far
    worse on a real run (25 s for ``sra_profile --sample-size 10000`` on an 11.3M-record mate
    file). Bound: 0.2 s, about three times the measured time, the best of three runs.
    """
    best, got = _best_of(lambda: sample_records([fastq_100k], 10_000, seed=0))
    assert len(got) == 10_000
    assert best < _scaled(0.2), _over("sample_records 100k to 10k", best, 0.2)


@pytest.mark.perf
def test_summary_and_screening_on_20000_rows_is_fast(containment_20000x5, tmp_path):
    """Summary at step 0.001 plus screening into an empty registry, as ``parse_containment`` runs them.

    Measured 2026-09-26 under coverage on an Apple-silicon laptop (macOS): 0.37 s after the change
    (4.34 s before it, with one frame copy per threshold, one ``iterrows`` pass, a re-read of the
    written table and one screening block conversion and timestamp per positive cell). Bound: 1.2 s,
    about three times the measured time, the best of two runs.
    """
    timings = []
    for attempt in range(2):
        registry = registry_mod.load_registry(tmp_path / f"registry_{attempt}.json")
        start = time.perf_counter()
        summary = _generate_containment_summary(
            containment_20000x5, tmp_path / "parsed.txt", tmp_path / "summary.txt", 0.001
        )
        recorded = registry_mod.record_screening_from_table(registry, summary.table, tmp_path / "matches")
        timings.append(time.perf_counter() - start)
    best = min(timings)
    assert len(summary.sample_to_genomes) == 20000
    assert recorded > 0
    assert best < _scaled(1.2), _over("summary and screening on 20,000 x 5", best, 1.2)


@pytest.mark.perf
def test_parse_metadata_on_300_files_is_fast(metadata_folder_300, tmp_path):
    """Parsing 300 files into a 300-row table with every distinct attribute as a column.

    Measured 2026-09-26 under coverage on an Apple-silicon laptop (macOS): 0.11 s after the change
    (0.28 s before it, with two parses per file, list-based tag lookups, dense rows and pandas'
    default CSV chunking). Bound: 0.4 s, about three times the measured time, the best of three
    runs. The bound alone does not tell the two apart at this size; the one-parse-per-file spy test
    in test_data_metadata.py guards the main change.
    """
    best, table = _best_of(lambda: parse_metadata(metadata_folder_300, tmp_path / "metadata_table.txt"))
    assert table.shape[0] == 300
    assert table.shape[1] > 30 + 40
    assert best < _scaled(0.4), _over("parse_metadata on 300 files", best, 0.4)


@pytest.mark.perf
def test_load_registry_peak_memory_on_20000_datasets_is_bounded(registry_20k):
    """``tracemalloc`` peak while ``load_registry`` parses the 20,000-dataset registry file.

    Measured 2026-09-26 under coverage on an Apple-silicon laptop (macOS): about 60 MB, well under
    the 250 MB bound the task-9 brief set, so that figure is kept rather than tightened, to leave
    headroom for a slower machine's allocator overhead and coverage instrumentation.
    """
    path, _ = registry_20k
    tracemalloc.start()
    try:
        registry_mod.load_registry(path)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    peak_mb = peak / (1024 * 1024)
    assert peak_mb < 250, f"load_registry peak was {peak_mb:.1f} MB on 20,000 datasets"


def test_scaled_multiplies_bounds_by_the_perf_scale(monkeypatch):
    """``_scaled`` applies the module's scale factor, as CI's ``PREPARE_PERF_SCALE=4`` sets it."""
    monkeypatch.setattr(sys.modules[__name__], "PERF_SCALE", 4.0)
    assert _scaled(1.0) == 4.0
    assert "PREPARE_PERF_SCALE=4" in _over("x", 5.0, 1.0)


@pytest.mark.parametrize("text", ["0.5", "0", "-2", "fast", "nan"])
def test_perf_scale_below_one_or_not_a_number_is_rejected(text):
    """A scale that would tighten the bounds, or is not a number, is refused with a clear message."""
    with pytest.raises(ValueError, match="PREPARE_PERF_SCALE"):
        _read_scale(text)


def test_perf_scale_reads_one_and_larger_values():
    """The default and a CI-style value parse to floats."""
    assert _read_scale("1") == 1.0
    assert _read_scale("4") == 4.0
