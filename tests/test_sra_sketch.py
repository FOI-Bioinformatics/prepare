"""Tests for the dataset sketcher behind ``sra_sketch``.

The block reader joins every sequence of a FASTQ file with ``N`` so one ``add_sequence`` call
per block gives the same hashes as one call per read; the sketcher, writer and containment
helpers are checked against sourmash's own API on the same reads.
"""

import gzip
import threading

import pytest

from prepare.data.sra.sampling import iter_sequence_blocks
from prepare.data.branchwater_search import BRANCHWATER_COLUMNS, KSIZE, SCALED, _select_sketch, sketch_fasta
from prepare.sra.sketch import (
    KSIZES,
    load_dataset_minhash,
    signature_path,
    sketch_dataset,
    sketch_is_current,
    sketch_summary,
    write_signature,
)
from prepare.sra.sketch_containment import containment_rows, genome_minhash, merge_containment_csv


def _write_fastq(path, reads, gz=False, crlf=False, final_newline=True):
    """Write ``reads`` as a four-line FASTQ file (plain or gzip), with the line ends requested."""
    end = "\r\n" if crlf else "\n"
    text = "".join(f"@r{i}{end}{seq}{end}+{end}{'I' * len(seq)}{end}" for i, seq in enumerate(reads))
    if not final_newline:
        text = text.rstrip("\r\n")
    data = text.encode()
    if gz:
        with gzip.open(path, "wb") as handle:
            handle.write(data)
    else:
        path.write_bytes(data)
    return path


def _reads_from_blocks(blocks):
    """The sequences a sequence of ``(block, n)`` pairs holds, and the sum of the counts."""
    reads = []
    total = 0
    for block, n in blocks:
        if n:
            reads.extend(block.split(b"N"))
        total += n
    return reads, total


class TestIterSequenceBlocks:
    def test_blocks_join_sequences_with_n(self, tmp_path):
        path = _write_fastq(tmp_path / "a.fastq", ["ACGT", "GGCC", "TTAA"])
        blocks = list(iter_sequence_blocks(path))
        assert blocks == [(b"ACGTNGGCCNTTAA", 3)]

    def test_blocks_read_gzip(self, tmp_path):
        path = _write_fastq(tmp_path / "a.fastq.gz", ["ACGT", "GGCC"], gz=True)
        assert list(iter_sequence_blocks(path)) == [(b"ACGTNGGCC", 2)]

    def test_blocks_strip_cr(self, tmp_path):
        path = _write_fastq(tmp_path / "a.fastq", ["ACGT", "GGCC"], crlf=True)
        assert list(iter_sequence_blocks(path)) == [(b"ACGTNGGCC", 2)]

    def test_blocks_accept_missing_final_newline(self, tmp_path):
        path = _write_fastq(tmp_path / "a.fastq", ["ACGT", "GGCC"], final_newline=False)
        assert _reads_from_blocks(iter_sequence_blocks(path)) == ([b"ACGT", b"GGCC"], 2)

    def test_record_spanning_block_boundary(self, tmp_path, monkeypatch):
        reads = [
            f"ACGT{i:04d}".replace("0", "A").replace("1", "C").replace("2", "G").replace("3", "T") for i in range(10)
        ]
        reads = [r.translate(str.maketrans("456789", "ACGTAC")) for r in reads]
        path = _write_fastq(tmp_path / "a.fastq", reads)
        monkeypatch.setattr("prepare.data.sra.sampling.CHUNK_SIZE", 64)
        blocks = list(iter_sequence_blocks(path))
        assert len(blocks) > 1
        got, total = _reads_from_blocks(blocks)
        assert got == [r.encode() for r in reads]
        assert total == 10

    def test_truncated_record_raises(self, tmp_path):
        path = tmp_path / "a.fastq"
        path.write_bytes(b"@r0\nACGT\n+\nIIII\n@r1\nACGT\n")
        with pytest.raises(ValueError, match="Truncated FASTQ record"):
            list(iter_sequence_blocks(path))

    def test_not_fastq_raises(self, tmp_path):
        path = tmp_path / "a.fastq"
        path.write_bytes(b">r0\nACGT\n")
        with pytest.raises(ValueError, match="Not a FASTQ file"):
            list(iter_sequence_blocks(path))

    def test_empty_file_yields_nothing(self, tmp_path):
        path = tmp_path / "a.fastq"
        path.write_bytes(b"")
        assert list(iter_sequence_blocks(path)) == []

    def test_stop_token_raises_interrupted_error(self, tmp_path):
        path = _write_fastq(tmp_path / "a.fastq", ["ACGT"])
        stop = threading.Event()
        stop.set()
        with pytest.raises(InterruptedError):
            list(iter_sequence_blocks(path, stop=stop))


# ------------------------------------------------------------------------------- sketcher and writer

sourmash = pytest.importorskip("sourmash")


def _random_reads(seed, count, length=1500):
    """Long random reads, so that a scaled-1000 sketch keeps some of their k-mers."""
    import random

    rng = random.Random(seed)
    return ["".join(rng.choice("ACGT") for _ in range(length)) for _ in range(count)]


READS_1 = _random_reads(1, 3)
READS_1[1] = READS_1[1].lower()
READS_1[2] = READS_1[2][:700] + "NN" + READS_1[2][702:]
READS_2 = _random_reads(2, 2)
ORPHAN = _random_reads(3, 1)


def _dataset(tmp_path, gz=False, accession="SRR1"):
    acc_dir = tmp_path / "fastq" / accession
    acc_dir.mkdir(parents=True)
    suffix = ".fastq.gz" if gz else ".fastq"
    files = [
        _write_fastq(acc_dir / f"{accession}_1{suffix}", READS_1, gz=gz),
        _write_fastq(acc_dir / f"{accession}_2{suffix}", READS_2, gz=gz),
        _write_fastq(acc_dir / f"{accession}{suffix}", ORPHAN, gz=gz),
    ]
    return acc_dir, sorted(files)


def _reference_minhash(ksize, reads):
    mh = sourmash.MinHash(n=0, ksize=ksize, scaled=SCALED, track_abundance=True)
    for read in reads:
        mh.add_sequence(read, force=True)
    return mh


class TestSketchDataset:
    @pytest.mark.parametrize("gz", [False, True])
    def test_sketch_equals_per_read_add_sequence(self, tmp_path, gz):
        _acc_dir, files = _dataset(tmp_path, gz=gz)
        sketch = sketch_dataset("SRR1", files)
        all_reads = ORPHAN + READS_1 + READS_2  # name order: SRR1.fastq, SRR1_1, SRR1_2
        for k in KSIZES:
            assert sketch.minhashes[k].hashes == _reference_minhash(k, all_reads).hashes, k
            assert sketch.hashes[k] == len(sketch.minhashes[k])
        assert sketch.reads == 6
        assert sketch.bases == sum(len(r) for r in all_reads)
        assert set(sketch.signature) == {f.name for f in files}

    def test_sketch_logs_progress_within_a_dataset(self, tmp_path, caplog, monkeypatch):
        import logging

        _acc_dir, files = _dataset(tmp_path)
        monkeypatch.setattr("prepare.data.sra.sampling.CHUNK_SIZE", 512)
        monkeypatch.setattr("prepare.sra.sketch.PROGRESS_BLOCKS", 2)
        with caplog.at_level(logging.DEBUG, logger="prepare.sra.sketch"):
            sketch_dataset("SRR1", files)
        assert any("reads so far" in r.getMessage() for r in caplog.records)

    def test_sketch_stops_on_the_token(self, tmp_path):
        _acc_dir, files = _dataset(tmp_path)
        stop = threading.Event()
        stop.set()
        with pytest.raises(InterruptedError):
            sketch_dataset("SRR1", files, stop=stop)


class TestWriteSignature:
    def test_written_file_is_a_sourmash_sig_gz_with_k21_first(self, tmp_path):
        _acc_dir, files = _dataset(tmp_path)
        sketch = sketch_dataset("SRR1", files)
        path = write_signature(sketch, signature_path("SRR1", tmp_path / "sketches"))

        assert path == tmp_path / "sketches" / "SRR1.sig.gz"
        assert path.read_bytes()[:2] == b"\x1f\x8b"
        loaded = list(sourmash.load_file_as_signatures(str(path)))
        assert [s.minhash.ksize for s in loaded] == list(KSIZES)
        assert all(s.name == "SRR1" and s.minhash.scaled == SCALED and s.minhash.track_abundance for s in loaded)
        import json

        document = json.loads(gzip.decompress(path.read_bytes()))
        assert document[0]["signatures"][0]["ksize"] == KSIZE
        assert _select_sketch(document[0], path)["ksize"] == KSIZE

    def test_no_file_is_left_when_writing_fails(self, tmp_path, monkeypatch):
        _acc_dir, files = _dataset(tmp_path)
        sketch = sketch_dataset("SRR1", files)
        out_dir = tmp_path / "sketches"

        def boom(*args, **kwargs):
            raise OSError("disk full")

        monkeypatch.setattr("sourmash.signature.save_signatures_to_json", boom)
        with pytest.raises(OSError):
            write_signature(sketch, signature_path("SRR1", out_dir))
        assert not (out_dir / "SRR1.sig.gz").exists()
        assert [p.name for p in out_dir.iterdir()] == []

    def test_summary_and_currency(self, tmp_path):
        _acc_dir, files = _dataset(tmp_path)
        sketch = sketch_dataset("SRR1", files)
        path = write_signature(sketch, signature_path("SRR1", tmp_path / "sketches"))
        summary = sketch_summary(sketch)

        assert set(summary) == {"ksizes", "scaled", "track_abundance", "reads", "bases", "hashes", "files", "source"}
        assert summary["ksizes"] == [21, 31, 51] and summary["scaled"] == 1000 and summary["source"] == "computed"
        assert set(summary["hashes"]) == {"21", "31", "51"}
        assert sketch_is_current(summary, files, path)
        assert not sketch_is_current(None, files, path)
        assert not sketch_is_current({**summary, "ksizes": [21]}, files, path)
        path.unlink()
        assert not sketch_is_current(summary, files, path)

    def test_load_dataset_minhash_flattens_and_checks_scaled(self, tmp_path):
        _acc_dir, files = _dataset(tmp_path)
        sketch = sketch_dataset("SRR1", files)
        path = write_signature(sketch, signature_path("SRR1", tmp_path / "sketches"))

        loaded = load_dataset_minhash(path)
        assert loaded.ksize == KSIZE and not loaded.track_abundance
        assert set(loaded.hashes) == set(sketch.minhashes[KSIZE].hashes)

        other = sourmash.MinHash(n=0, ksize=21, scaled=2000)
        other.add_sequence("ACGT" * 50, force=True)
        with open(tmp_path / "wrong.sig", "w") as handle:
            sourmash.signature.save_signatures_to_json([sourmash.SourmashSignature(other, name="x")], handle)
        from prepare.core.exceptions import DataAccessError

        with pytest.raises(DataAccessError, match="scaled"):
            load_dataset_minhash(tmp_path / "wrong.sig")


# ----------------------------------------------------------------------------------- containment


def _genome(tmp_path):
    """A genome whose k-mers partly overlap the dataset reads."""
    import random

    rng = random.Random(4)
    tail = "".join(rng.choice("ACGT") for _ in range(20000))
    path = tmp_path / "genome.fasta"
    path.write_text(">g1\n" + READS_1[0] + tail + "\n>g2\n" + READS_2[0] + "\n")
    return path


class TestContainment:
    def test_containment_matches_contained_by(self, tmp_path):
        _acc_dir, files = _dataset(tmp_path)
        path = write_signature(sketch_dataset("SRR1", files), signature_path("SRR1", tmp_path / "sketches"))
        genome_mh = genome_minhash(sketch_fasta(_genome(tmp_path)))

        rows = containment_rows(genome_mh, {"SRR1": path, "SRR9": tmp_path / "missing.sig.gz"})

        expected = genome_mh.contained_by(load_dataset_minhash(path))
        assert rows == [("SRR1", expected)]
        assert 0 < expected < 1

    def test_merge_keeps_other_rows_and_metadata(self, tmp_path):
        import csv

        matches = tmp_path / "matches"
        matches.mkdir()
        target = matches / "genome.csv"
        with open(target, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=BRANCHWATER_COLUMNS)
            writer.writeheader()
            writer.writerow({"acc": "SRR9", "containment": "0.5", "cANI": "0.9675", "organism": "soil metagenome"})
            writer.writerow({"acc": "SRR1", "containment": "0.9", "cANI": "0.995"})

        merge_containment_csv(target, [("SRR1", 0.25), ("SRR2", 0.75)])

        with open(target, newline="") as handle:
            rows = list(csv.DictReader(handle))
        assert [r["acc"] for r in rows] == ["SRR2", "SRR9", "SRR1"]
        assert rows[1]["organism"] == "soil metagenome"
        by_acc = {r["acc"]: r for r in rows}
        assert by_acc["SRR1"]["containment"] == "0.2500"
        assert by_acc["SRR2"]["cANI"] == str(round(0.75 ** (1 / KSIZE), 4))
        assert list(rows[0].keys()) == BRANCHWATER_COLUMNS

    def test_merge_keeps_the_metadata_of_a_replaced_branchwater_row(self, tmp_path):
        """Replacing a Branchwater row updates containment and cANI only; its metadata survives."""
        import csv

        target = tmp_path / "matches" / "genome.csv"
        target.parent.mkdir()
        with open(target, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=BRANCHWATER_COLUMNS)
            writer.writeheader()
            writer.writerow(
                {"acc": "SRR1", "containment": "0.9", "cANI": "0.995", "organism": "gut metagenome", "lat_lon": "1 2"}
            )
            writer.writerow({"acc": "SRR1", "containment": "0.2", "cANI": "0.9", "organism": "duplicate"})

        merge_containment_csv(target, [("SRR1", 0.5)])

        with open(target, newline="") as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == 1
        assert rows[0]["containment"] == "0.5000" and rows[0]["organism"] == "gut metagenome"
        assert rows[0]["lat_lon"] == "1 2"

    def test_merge_drops_rows_of_accessions_that_failed(self, tmp_path):
        target = tmp_path / "matches" / "genome.csv"
        merge_containment_csv(target, [("SRR1", 0.5), ("SRR2", 0.4)])
        merge_containment_csv(target, [("SRR1", 0.6)], drop={"SRR2"})
        import csv

        with open(target, newline="") as handle:
            assert [r["acc"] for r in csv.DictReader(handle)] == ["SRR1"]

    def test_merged_csv_round_trips_through_parse_containment_data(self, tmp_path):
        from prepare.data.branchwater import parse_containment_data

        matches = tmp_path / "matches"
        matches.mkdir()
        merge_containment_csv(matches / "wMel.csv", [("SRR1", 0.25), ("SRR2", 0.0)])

        summary = parse_containment_data(matches, tmp_path / "parsed.txt", tmp_path / "top.txt")
        assert "wMel" in summary.table.columns
        assert summary.table.loc["SRR1", "wMel"] == pytest.approx(0.25)


# ------------------------------------------------------------------------------------ store cache


class TestStoreCache:
    """A signature cached beside the store copy, mirroring the statistics record."""

    def _store(self, tmp_path, monkeypatch):
        from prepare.store.layout import init_store, sidecar_path as dataset_sidecar_path, sra_dir
        from prepare.store.sidecar import Sidecar, write_sidecar

        monkeypatch.setenv("PREPARE_DATA", str(tmp_path / "store"))
        paths = init_store(tmp_path / "store")
        acc_dir = sra_dir(paths, "SRR1")
        acc_dir.mkdir(parents=True)
        files = [
            _write_fastq(acc_dir / "SRR1_1.fastq", READS_1),
            _write_fastq(acc_dir / "SRR1_2.fastq", READS_2),
        ]
        sidecar_path = dataset_sidecar_path(paths, "SRR1")
        write_sidecar(sidecar_path, Sidecar(accession="SRR1", state="complete"))
        return paths, acc_dir, sidecar_path, files

    def _signature(self, tmp_path, files):
        return write_signature(sketch_dataset("SRR1", files), signature_path("SRR1", tmp_path / "sketches"))

    def test_cached_sketch_requires_a_matching_block_and_file(self, tmp_path, monkeypatch):
        from prepare.store.sketches import cached_sketch, store_sketch
        from prepare.store.sidecar import read_sidecar

        paths, acc_dir, sidecar_path, files = self._store(tmp_path, monkeypatch)
        assert cached_sketch(acc_dir, sidecar_path, files) is None
        assert cached_sketch(acc_dir, None, files) is None

        signature = self._signature(tmp_path, files)
        summary = sketch_summary(sketch_dataset("SRR1", files))
        store_sketch(sidecar_path, signature, summary)

        cached = cached_sketch(acc_dir, sidecar_path, files)
        assert cached is not None
        path, block = cached
        assert path == acc_dir / "SRR1.sig.gz" and path.is_file()
        assert block["files"] == summary["files"] and block["file"] == "SRR1.sig.gz" and block["computed"]
        sidecar = read_sidecar(sidecar_path)
        assert sidecar.state == "complete" and sidecar.sketch == block

        path.unlink()
        assert cached_sketch(acc_dir, sidecar_path, files) is None

    def test_a_changed_file_invalidates_the_cache(self, tmp_path, monkeypatch):
        import os

        from prepare.store.sketches import cached_sketch, store_sketch

        _paths, acc_dir, sidecar_path, files = self._store(tmp_path, monkeypatch)
        store_sketch(sidecar_path, self._signature(tmp_path, files), sketch_summary(sketch_dataset("SRR1", files)))
        assert cached_sketch(acc_dir, sidecar_path, files) is not None

        os.utime(files[0], (1, 1))
        assert cached_sketch(acc_dir, sidecar_path, files) is None

    def test_store_sketch_takes_the_dataset_lock(self, tmp_path, monkeypatch):
        from prepare.store import sketches as sketches_module
        from prepare.store.sketches import store_sketch
        from prepare.store.stats import STATS_LOCK_WAIT_SECONDS
        from unittest.mock import patch
        from contextlib import nullcontext

        paths, _acc_dir, sidecar_path, files = self._store(tmp_path, monkeypatch)
        held = []

        def record_lock(lock_paths, accession, wait_seconds=0.0, should_stop=None):
            held.append((lock_paths.root, accession, wait_seconds))
            return nullcontext()

        with patch.object(sketches_module, "dataset_lock", side_effect=record_lock):
            store_sketch(sidecar_path, self._signature(tmp_path, files), sketch_summary(sketch_dataset("SRR1", files)))
        assert held == [(paths.root, "SRR1", STATS_LOCK_WAIT_SECONDS)]

    def test_a_held_lock_only_costs_the_cache(self, tmp_path, monkeypatch, caplog):
        import logging

        from prepare.store import locks as locks_module
        from prepare.store.layout import lock_path
        from prepare.store.sidecar import read_sidecar
        from prepare.store.sketches import store_sketch

        paths, acc_dir, sidecar_path, files = self._store(tmp_path, monkeypatch)
        lock_path(paths, "SRR1").write_text('{"pid": 999, "host": "other-host", "started": "now"}')
        monkeypatch.setattr(locks_module, "LOCK_POLL_SECONDS", 0.01)
        monkeypatch.setattr("prepare.store.sketches.STATS_LOCK_WAIT_SECONDS", 0.02)

        with caplog.at_level(logging.WARNING):
            store_sketch(sidecar_path, self._signature(tmp_path, files), sketch_summary(sketch_dataset("SRR1", files)))

        assert read_sidecar(sidecar_path).sketch is None
        assert not (acc_dir / "SRR1.sig.gz").exists()
        assert "Not caching sketch for SRR1" in caplog.text

    def test_store_sketch_stops_waiting_when_stop_is_set(self, tmp_path, monkeypatch, caplog):
        """A stop request ends the wait for a held lock at once; the cache is skipped with a warning."""
        import logging
        import time

        from prepare.store import locks as locks_module
        from prepare.store.layout import lock_path
        from prepare.store.sidecar import read_sidecar
        from prepare.store.sketches import store_sketch

        paths, acc_dir, sidecar_path, files = self._store(tmp_path, monkeypatch)
        lock_path(paths, "SRR1").write_text('{"pid": 999, "host": "other-host", "started": "now"}')
        monkeypatch.setattr(locks_module, "LOCK_POLL_SECONDS", 0.01)
        monkeypatch.setattr("prepare.store.sketches.STATS_LOCK_WAIT_SECONDS", 30.0)
        signature = self._signature(tmp_path, files)
        summary = sketch_summary(sketch_dataset("SRR1", files))

        started = time.monotonic()
        with caplog.at_level(logging.WARNING):
            store_sketch(sidecar_path, signature, summary, should_stop=lambda: True)

        assert time.monotonic() - started < 2.0
        assert read_sidecar(sidecar_path).sketch is None
        assert not (acc_dir / "SRR1.sig.gz").exists()
        assert "Not caching sketch for SRR1" in caplog.text

    def test_store_sketch_is_a_noop_without_a_sidecar(self, tmp_path, monkeypatch):
        from prepare.store.sketches import store_sketch

        _paths, _acc_dir, _sidecar_path, files = self._store(tmp_path, monkeypatch)
        signature = self._signature(tmp_path, files)
        store_sketch(None, signature, sketch_summary(sketch_dataset("SRR1", files)))
        store_sketch(tmp_path / "nowhere" / "SRR9.json", signature, sketch_summary(sketch_dataset("SRR1", files)))
        assert not (tmp_path / "nowhere").exists()
