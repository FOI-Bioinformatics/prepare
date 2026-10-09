"""
Tests for prepare.store.stats: one cached FASTQ statistics record per dataset, computed
once and shared by sra_profile and sra_report instead of each command re-parsing every read
on its own.
"""

import gzip
import threading
from contextlib import nullcontext
from pathlib import Path
from unittest.mock import Mock, patch

from prepare.store.sidecar import Sidecar, read_sidecar, write_sidecar
from prepare.store.stats import cached_stats, compute_dataset_stats, store_stats


def _write_fastq(path: Path, reads, gz: bool = False) -> None:
    """Write a synthetic FASTQ file: one record per string in ``reads``, quality all 'I'."""
    lines = []
    for i, seq in enumerate(reads):
        lines.append(f"@read{i}")
        lines.append(seq)
        lines.append("+")
        lines.append("I" * len(seq))
    content = "\n".join(lines) + "\n"
    if gz:
        with gzip.open(path, "wt") as f:
            f.write(content)
    else:
        path.write_text(content)


class TestComputeDatasetStatsCounts:
    """Read counts must be correct for both plain and gzip-compressed FASTQ files."""

    def test_counts_plain_file(self, tmp_path):
        fastq = tmp_path / "SRR1.fastq"
        _write_fastq(fastq, ["ACGTACGTAC"] * 20)

        stats = compute_dataset_stats([fastq], sample_size=1000, use_seqkit=False)

        assert stats["reads_per_file"][fastq.name] == 20
        assert stats["reads_total"] == 20
        assert stats["sample_size"] == 1000
        assert stats["sampled"] is False
        assert "computed" in stats and stats["computed"]

    def test_counts_gz_file(self, tmp_path):
        fastq = tmp_path / "SRR1.fastq.gz"
        _write_fastq(fastq, ["ACGTACGTAC"] * 20, gz=True)

        stats = compute_dataset_stats([fastq], sample_size=1000, use_seqkit=False)

        assert stats["reads_per_file"][fastq.name] == 20
        assert stats["reads_total"] == 20

    def test_counts_a_zero_length_record(self, tmp_path):
        """A read trimmed to length zero is a legal record, not a truncated file."""
        fastq = tmp_path / "SRR1.fastq"
        fastq.write_text("@r1\nACGT\n+\nIIII\n@r2\n\n+\n\n@r3\nACGT\n+\nIIII\n")

        stats = compute_dataset_stats([fastq], sample_size=1000, use_seqkit=False)

        assert stats["reads_total"] == 3
        assert stats["min_read_length"] == 0

    def test_signature_tracks_size_and_mtime(self, tmp_path):
        fastq = tmp_path / "SRR1.fastq"
        _write_fastq(fastq, ["ACGT"] * 5)

        stats = compute_dataset_stats([fastq], use_seqkit=False)

        st = fastq.stat()
        assert stats["signature"][fastq.name] == [st.st_size, st.st_mtime]

    def test_multiple_files_aggregate(self, tmp_path):
        f1 = tmp_path / "SRR1_1.fastq"
        f2 = tmp_path / "SRR1_2.fastq"
        _write_fastq(f1, ["ACGT"] * 7)
        _write_fastq(f2, ["ACGT"] * 3)

        stats = compute_dataset_stats([f1, f2], use_seqkit=False)

        assert stats["reads_per_file"][f1.name] == 7
        assert stats["reads_per_file"][f2.name] == 3
        assert stats["reads_total"] == 10


class TestComputeDatasetStatsFigures:
    """Length, base and GC figures of the record (moved here from the removed
    calculate_read_statistics tests, the second statistics path before 0.5.0)."""

    def test_lengths_bases_and_gc(self, tmp_path):
        fastq = tmp_path / "SRR1.fastq"
        _write_fastq(fastq, ["ATCGATCGATCGATCG", "GCTAGCTAGCTAGCTA", "A" * 16, "G" * 16])

        stats = compute_dataset_stats([fastq], use_seqkit=False)

        assert stats["reads_total"] == 4
        assert stats["bases_total"] == 64
        assert stats["avg_read_length"] == 16.0
        assert (stats["min_read_length"], stats["max_read_length"]) == (16, 16)
        assert stats["gc_content"] == 0.5  # a 0-1 fraction; sra_profile reports it times 100

    def test_n50(self, tmp_path):
        fastq = tmp_path / "SRR1.fastq"
        _write_fastq(fastq, ["A" * 100, "T" * 200, "C" * 300, "G" * 400])

        # 1000 bases in total; the cumulative sum from the longest read reaches 500 at 300 bp.
        assert compute_dataset_stats([fastq], use_seqkit=False)["n50"] == 300

    def test_small_file_is_not_sampled(self, tmp_path):
        fastq = tmp_path / "SRR1.fastq"
        _write_fastq(fastq, ["ACGT"] * 5)

        stats = compute_dataset_stats([fastq], use_seqkit=False)

        assert stats["reads_total"] == 5 and stats["sampled"] is False


class TestUniformSampling:
    """The reservoir sample must be uniform over the whole file, not biased to its head."""

    def test_sampled_flag_set_when_file_exceeds_sample_size(self, tmp_path):
        fastq = tmp_path / "SRR1.fastq"
        _write_fastq(fastq, ["ACGTACGT"] * 2000)

        stats = compute_dataset_stats([fastq], sample_size=500, use_seqkit=False)

        assert stats["reads_total"] == 2000
        assert stats["sampled"] is True
        assert stats["sample_size"] == 500

    def test_sample_reflects_tail_composition_of_a_50000_read_file(self, tmp_path):
        """40,000 head reads are all-AT (GC=0); 10,000 tail reads are all-GC (GC=1).

        A sampler biased toward the head of the file would see gc_content == 0.0; a
        uniform reservoir sample of 1,000 reads drawn from 50,000 must pick up some of
        the tail's GC-rich reads too (expected ~20% of the sample, well above noise).
        """
        head_reads = ["AAAAAAAA"] * 40000
        tail_reads = ["CCCCCCCC"] * 10000
        fastq = tmp_path / "SRR1.fastq"
        _write_fastq(fastq, head_reads + tail_reads)

        stats = compute_dataset_stats([fastq], sample_size=1000, use_seqkit=False)

        assert stats["reads_total"] == 50000
        assert stats["sampled"] is True
        assert stats["gc_content"] > 0.05


class TestCachedStats:
    """``cached_stats`` returns the sidecar's stats only while the signature matches."""

    def _seed(self, acc_dir: Path, reads):
        acc_dir.mkdir(parents=True, exist_ok=True)
        fastq = acc_dir / f"{acc_dir.name}.fastq"
        _write_fastq(fastq, reads)
        stats = compute_dataset_stats([fastq], use_seqkit=False)
        sidecar_path = acc_dir / f"{acc_dir.name}.json"
        write_sidecar(sidecar_path, Sidecar(accession=acc_dir.name, stats=stats))
        return fastq, sidecar_path, stats

    def test_cache_hit_when_signature_matches(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        _fastq, sidecar_path, stats = self._seed(acc_dir, ["ACGT"] * 10)

        hit = cached_stats(acc_dir, sidecar_path)

        assert hit == stats

    def test_cache_invalidated_by_size_change(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        fastq, sidecar_path, _stats = self._seed(acc_dir, ["ACGT"] * 10)

        # Rewrite the file with different content -> different size (and mtime).
        _write_fastq(fastq, ["ACGTACGTAC"] * 20)

        assert cached_stats(acc_dir, sidecar_path) is None

    def test_no_sidecar_path_returns_none(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()

        assert cached_stats(acc_dir, None) is None

    def test_missing_sidecar_file_returns_none(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()

        assert cached_stats(acc_dir, acc_dir / "missing.json") is None

    def test_sidecar_without_stats_returns_none(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        sidecar_path = acc_dir / "SRR1.json"
        write_sidecar(sidecar_path, Sidecar(accession="SRR1"))

        assert cached_stats(acc_dir, sidecar_path) is None


class TestStoreStats:
    """``store_stats`` writes through the sidecar, and is a no-op without one."""

    def test_writes_stats_and_timestamp_into_sidecar(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        sidecar_path = acc_dir / "SRR1.json"
        write_sidecar(sidecar_path, Sidecar(accession="SRR1"))

        stats = {"reads_total": 10, "signature": {}}
        store_stats(sidecar_path, stats)

        sidecar = read_sidecar(sidecar_path)
        assert sidecar.stats == stats
        assert sidecar.stats_computed is not None

    def test_noop_without_sidecar_path(self):
        # Must not raise: a project without a shared store has no sidecar to update.
        store_stats(None, {"reads_total": 1})

    def test_noop_when_sidecar_file_missing(self, tmp_path):
        store_stats(tmp_path / "missing.json", {"reads_total": 1})
        assert not (tmp_path / "missing.json").exists()


class TestStoreStatsConcurrency:
    """A statistics write must not clobber what another project published meanwhile."""

    def _store_sidecar(self, tmp_path, state="partial"):
        """A store rooted at ``tmp_path/store`` holding one accession's sidecar."""
        from prepare.store.layout import init_store, sidecar_path as dataset_sidecar_path, sra_dir

        paths = init_store(tmp_path / "store")
        acc_dir = sra_dir(paths, "SRR1")
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        path = dataset_sidecar_path(paths, "SRR1")
        write_sidecar(path, Sidecar(accession="SRR1", state=state))
        return paths, path

    def test_a_concurrent_state_change_survives_the_statistics_write(self, tmp_path):
        """A statistics write only replaces the two statistics fields.

        Computing the record takes seconds to minutes on a real dataset. Whatever another
        project published for the same accession in that time -- a completed re-download,
        with its new state and file list -- is still there afterwards.
        """
        _, sidecar_path = self._store_sidecar(tmp_path, state="partial")

        published = read_sidecar(sidecar_path)  # another project re-downloads and publishes
        published.state = "complete"
        published.files = [{"name": "SRR1_1.fastq", "bytes": 1, "md5": "x", "reads": 48000000}]
        write_sidecar(sidecar_path, published)

        store_stats(sidecar_path, {"reads_total": 1, "signature": {}})

        after = read_sidecar(sidecar_path)
        assert after.state == "complete"
        assert [f["name"] for f in after.files] == ["SRR1_1.fastq"]
        assert after.stats["reads_total"] == 1

    def test_the_lock_is_held_across_the_read_and_the_write(self, tmp_path, monkeypatch):
        """Another sidecar writer is excluded for the whole read-modify-write.

        That exclusion is what stops a concurrent publish being overwritten: without it, a
        download or ``store_verify --fix-state`` running under the same lock could write
        between this function's read and its write.
        """
        from prepare.core.exceptions import DataAccessError
        from prepare.store import locks as locks_module
        from prepare.store import stats as stats_module
        from prepare.store.locks import dataset_lock

        paths, sidecar_path = self._store_sidecar(tmp_path)
        monkeypatch.setattr(locks_module, "LOCK_POLL_SECONDS", 0.001)
        observed = {}
        real_write = stats_module._write_stats

        def other_writer():
            try:
                with dataset_lock(paths, "SRR1", wait_seconds=0.2):
                    observed["other_writer"] = "acquired the lock"
            except DataAccessError:
                observed["other_writer"] = "blocked"

        def write_while_another_writer_tries(path, stats):
            # Another writer, in a thread of its own, waits for the lock and gives up.
            thread = threading.Thread(target=other_writer)
            thread.start()
            thread.join()
            real_write(path, stats)

        with patch.object(stats_module, "_write_stats", side_effect=write_while_another_writer_tries):
            store_stats(sidecar_path, {"reads_total": 7, "signature": {}})

        assert observed["other_writer"] == "blocked"
        assert read_sidecar(sidecar_path).stats["reads_total"] == 7

    def test_the_write_takes_the_dataset_lock(self, tmp_path):
        from prepare.store import stats as stats_module

        paths, sidecar_path = self._store_sidecar(tmp_path)
        held = []

        def record_lock(lock_paths, accession, wait_seconds=0.0):
            held.append((lock_paths.root, accession, wait_seconds))
            return nullcontext(tmp_path)

        with patch.object(stats_module, "dataset_lock", side_effect=record_lock):
            store_stats(sidecar_path, {"reads_total": 7, "signature": {}})

        assert held == [(paths.root, "SRR1", stats_module.STATS_LOCK_WAIT_SECONDS)]
        assert read_sidecar(sidecar_path).stats["reads_total"] == 7

    def test_a_locked_dataset_skips_the_cache_with_a_warning(self, tmp_path, caplog, monkeypatch):
        """A download holding the lock for hours must not stall the analysis that asked for
        the statistics: the record is simply not cached this time."""
        import logging

        from prepare.store import locks as locks_module
        from prepare.store.layout import lock_path

        paths, sidecar_path = self._store_sidecar(tmp_path)
        lock_path(paths, "SRR1").write_text('{"pid": 999, "host": "other-host", "started": "now"}')

        monkeypatch.setattr(locks_module, "LOCK_POLL_SECONDS", 0.01)
        monkeypatch.setattr("prepare.store.stats.STATS_LOCK_WAIT_SECONDS", 0.02)

        with caplog.at_level(logging.WARNING):
            store_stats(sidecar_path, {"reads_total": 7, "signature": {}})

        assert read_sidecar(sidecar_path).stats == {}
        assert "Not caching statistics for SRR1" in caplog.text
        assert "other-host" in caplog.text

    def test_a_sidecar_outside_a_store_is_written_without_a_lock(self, tmp_path):
        """A sidecar beside a plain project folder has no store lock directory to use."""
        from prepare.store import stats as stats_module

        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        sidecar_path = acc_dir / "SRR1.json"
        write_sidecar(sidecar_path, Sidecar(accession="SRR1"))

        with patch.object(stats_module, "dataset_lock") as lock:
            store_stats(sidecar_path, {"reads_total": 3, "signature": {}})

        lock.assert_not_called()
        assert read_sidecar(sidecar_path).stats["reads_total"] == 3


class TestSeqkitPath:
    """When seqkit is installed, exact per-file counts come from ``seqkit stats -T``."""

    def test_uses_seqkit_when_available(self, tmp_path):
        fastq = tmp_path / "SRR1.fastq"
        _write_fastq(fastq, ["ACGTACGTAC"] * 20)  # 20 reads x 10bp = 200 bases

        table = (
            "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n"
            f"{fastq}\tFASTQ\tDNA\t20\t200\t10\t10.0\t10\n"
        )
        fake_result = Mock(stdout=table, returncode=0)

        with patch("prepare.store.stats.shutil.which", return_value="/usr/bin/seqkit"):
            with patch("prepare.store.stats.SecureSubprocess.run_secure", return_value=fake_result) as mock_run:
                stats = compute_dataset_stats([fastq], sample_size=1000, use_seqkit=True)

        mock_run.assert_called_once()
        call_args = mock_run.call_args[0]
        assert call_args[0] == "seqkit"
        assert call_args[1][:3] == ["stats", "-T", "-j"]

        assert stats["reads_total"] == 20
        assert stats["bases_total"] == 200
        assert stats["avg_read_length"] == 10.0
        assert stats["min_read_length"] == 10
        assert stats["max_read_length"] == 10

    def test_rows_are_matched_to_files_by_name_not_output_order(self, tmp_path):
        """``seqkit stats -j`` makes no promise about row order, and swapping two mates'
        counts would misreport each file."""
        mate1 = tmp_path / "SRR1_1.fastq"
        mate2 = tmp_path / "SRR1_2.fastq"
        _write_fastq(mate1, ["ACGTACGTAC"] * 20)
        _write_fastq(mate2, ["ACGTACGTAC"] * 30)

        table = (
            "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n"
            f"{mate2}\tFASTQ\tDNA\t30\t300\t10\t10.0\t10\n"  # second file reported first
            f"{mate1}\tFASTQ\tDNA\t20\t200\t10\t10.0\t10\n"
        )

        with patch("prepare.store.stats.shutil.which", return_value="/usr/bin/seqkit"):
            with patch("prepare.store.stats.SecureSubprocess.run_secure", return_value=Mock(stdout=table)):
                stats = compute_dataset_stats([mate1, mate2], sample_size=1000, use_seqkit=True)

        assert stats["reads_per_file"] == {"SRR1_1.fastq": 20, "SRR1_2.fastq": 30}
        assert stats["reads_total"] == 50

    def test_a_row_naming_an_unrequested_file_falls_back(self, tmp_path):
        fastq = tmp_path / "SRR1.fastq"
        _write_fastq(fastq, ["ACGTACGTAC"] * 5)
        table = (
            "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n"
            "somewhere/else.fastq\tFASTQ\tDNA\t999\t9990\t10\t10.0\t10\n"
        )

        with patch("prepare.store.stats.shutil.which", return_value="/usr/bin/seqkit"):
            with patch("prepare.store.stats.SecureSubprocess.run_secure", return_value=Mock(stdout=table)):
                stats = compute_dataset_stats([fastq], sample_size=1000, use_seqkit=True)

        assert stats["reads_total"] == 5  # the streaming count, not seqkit's row

    def test_falls_back_when_seqkit_not_installed(self, tmp_path):
        fastq = tmp_path / "SRR1.fastq"
        _write_fastq(fastq, ["ACGTACGTAC"] * 5)

        with patch("prepare.store.stats.shutil.which", return_value=None):
            with patch("prepare.store.stats.SecureSubprocess.run_secure") as mock_run:
                stats = compute_dataset_stats([fastq], use_seqkit=True)

        mock_run.assert_not_called()
        assert stats["reads_total"] == 5

    def test_use_seqkit_false_skips_seqkit_even_if_installed(self, tmp_path):
        fastq = tmp_path / "SRR1.fastq"
        _write_fastq(fastq, ["ACGT"] * 3)

        with patch("prepare.store.stats.shutil.which", return_value="/usr/bin/seqkit"):
            with patch("prepare.store.stats.SecureSubprocess.run_secure") as mock_run:
                stats = compute_dataset_stats([fastq], use_seqkit=False)

        mock_run.assert_not_called()
        assert stats["reads_total"] == 3

    def test_seqkit_failure_falls_back_to_streaming_counts(self, tmp_path):
        fastq = tmp_path / "SRR1.fastq"
        _write_fastq(fastq, ["ACGT"] * 4)

        with patch("prepare.store.stats.shutil.which", return_value="/usr/bin/seqkit"):
            with patch("prepare.store.stats.SecureSubprocess.run_secure", side_effect=RuntimeError("boom")):
                stats = compute_dataset_stats([fastq], use_seqkit=True)

        assert stats["reads_total"] == 4
