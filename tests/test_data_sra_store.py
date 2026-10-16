"""Tests for prepare.data.sra: downloads that use the shared data store, its locks and atomic publishing."""

import gzip
import json
import threading
import pytest
from pathlib import Path
from unittest.mock import patch
from prepare.data.sra.fastq import accession_has_fastq
from prepare.data.sra.accession import classify_download_error
from prepare.data.sra.download import download_sra
from prepare.data.sra import store_handoff as store_handoff_mod


class TestDownloadSraStore:
    """download_sra with a shared store: link what is already there, download what is not."""

    @staticmethod
    def _store(tmp_path):
        from prepare.store.layout import init_store

        return init_store(tmp_path / "store")

    @staticmethod
    def _accessions(tmp_path, *accessions):
        acc = tmp_path / "acc.txt"
        acc.write_text("".join(f"{a}\n" for a in accessions))
        return acc

    @staticmethod
    def _store_dataset(paths, accession="SRR1", state="complete"):
        from prepare.store.sidecar import Sidecar, write_sidecar

        acc_dir = paths.sra / accession
        acc_dir.mkdir(parents=True, exist_ok=True)
        (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        write_sidecar(acc_dir / f"{accession}.json", Sidecar(accession=accession, state=state))
        return acc_dir

    @staticmethod
    def _fake_download(calls):
        """A download_accession stand-in that writes one FASTQ under <output_folder>/<acc>."""

        def _download(accession, output_folder, *args, **kwargs):
            calls.append((accession, Path(output_folder), kwargs))
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True, exist_ok=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            return True, "Downloaded 1 files, unverified"

        return _download

    def test_complete_in_store_links_without_downloading(self, tmp_path):
        paths = self._store(tmp_path)
        self._store_dataset(paths)
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert calls == []
        assert stats["successful"] == 1
        assert stats["results"]["SRR1"].startswith("linked from store")
        assert (fastq_folder / "SRR1").is_symlink()
        assert (fastq_folder / "SRR1").resolve() == (paths.sra / "SRR1").resolve()

    def test_summary_counts_links_separately(self, tmp_path, caplog):
        """The final summary must not count a store link as a newly downloaded dataset
        (audit deferred S5-3)."""
        paths = self._store(tmp_path)
        self._store_dataset(paths, accession="SRR1")
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            with caplog.at_level("INFO", logger="prepare.data.sra"):
                stats = download_sra(
                    fastq_folder, self._accessions(tmp_path, "SRR1", "SRR2"), store=paths, max_retries=0
                )

        assert stats["results"]["SRR1"].startswith("linked from store")
        assert stats["results"]["SRR2"].endswith("; stored")
        assert "Newly downloaded: 1" in caplog.text
        assert "Linked from store: 1" in caplog.text

    def test_unverified_sidecar_also_links(self, tmp_path):
        paths = self._store(tmp_path)
        self._store_dataset(paths, state="unverified")
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert calls == []
        assert stats["successful"] == 1
        assert (fastq_folder / "SRR1").is_symlink()

    def test_a_fresh_store_download_links_under_its_own_dataset_lock(self, tmp_path):
        """The hand-off already holds the dataset lock when it links a fresh copy, so it asks
        link_dataset not to take the lock a second time."""
        from prepare.store import link as link_mod
        from prepare.store.layout import lock_path

        paths = self._store(tmp_path)
        fastq_folder = tmp_path / "project" / "fastq"
        original = link_mod.link_dataset
        seen = []

        def _spy(project_fastq, accession, store, mode="auto", locked=False):
            seen.append((locked, lock_path(store, accession).exists()))
            return original(project_fastq, accession, store, mode=mode, locked=locked)

        calls = []
        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            with patch("prepare.store.link.link_dataset", side_effect=_spy):
                stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert stats["successful"] == 1
        assert seen == [(True, True)]
        assert (fastq_folder / "SRR1").is_symlink()

    def test_a_ready_copy_whose_lock_is_briefly_held_waits_and_then_links(self, tmp_path):
        """A ready copy found without the lock is linked under it; when another run holds the
        lock at that moment the hand-off waits for it instead of failing."""
        from prepare.store.locks import dataset_lock

        paths = self._store(tmp_path)
        self._store_dataset(paths)
        fastq_folder = tmp_path / "project" / "fastq"
        held, release = threading.Event(), threading.Event()

        def _holder():
            with dataset_lock(paths, "SRR1"):
                held.set()
                release.wait(5)

        thread = threading.Thread(target=_holder)
        thread.start()
        held.wait(5)
        timer = threading.Timer(0.3, release.set)
        timer.start()
        try:
            result = store_handoff_mod._store_download("SRR1", fastq_folder, paths)
        finally:
            release.set()
            thread.join(5)
            timer.cancel()

        assert result[0] is True, result
        assert result[1].startswith("linked from store")
        assert (fastq_folder / "SRR1").is_symlink()

    def test_absent_downloads_into_the_store_and_links(self, tmp_path):
        from prepare.store.catalog import Catalog
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert len(calls) == 1
        accession, output_folder, kwargs = calls[0]
        assert accession == "SRR1"
        # The download lands in the store's tmp folder and is only published into sra/ once
        # it is verified, compressed and described.
        assert output_folder == paths.tmp
        assert Path(kwargs["staging_folder"]) == paths.tmp
        assert stats["successful"] == 1
        assert stats["results"]["SRR1"].endswith("; stored")
        assert "Downloaded 1 files, unverified" in stats["results"]["SRR1"]

        sidecar = read_sidecar(paths.sra / "SRR1" / "SRR1.json")
        assert sidecar is not None and sidecar.accession == "SRR1"
        assert sidecar.compression == "none"
        with Catalog(paths) as catalog:
            assert catalog.get_dataset("SRR1") is not None
        assert (fastq_folder / "SRR1").is_symlink()
        assert not (paths.locks / "SRR1.lock").exists()

    def test_catalogue_busy_during_publish_still_links_and_reindex_repairs_it(self, tmp_path, monkeypatch, caplog):
        """A catalogue lock held by another writer through the whole publish window must not
        turn a completed, published download into a reported failure: the dataset is linked
        into the project and the message says the catalogue write is pending, and a later
        store_reindex (called directly here) picks up the sidecar already on disk."""
        import threading

        from prepare.store.catalog import Catalog, catalog_write
        from prepare.store.sidecar import read_sidecar

        monkeypatch.setattr("prepare.store.catalog.CATALOG_LOCK_WAIT_SECONDS", 0.2)

        paths = self._store(tmp_path)
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        barrier_entered = threading.Event()
        release = threading.Event()

        def hold_catalogue_lock():
            with catalog_write(paths) as catalog:
                catalog.migrate()
                barrier_entered.set()
                release.wait(timeout=5)

        holder = threading.Thread(target=hold_catalogue_lock)
        holder.start()
        assert barrier_entered.wait(timeout=2)

        try:
            with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
                with caplog.at_level("WARNING", logger="prepare.data.sra.store_handoff"):
                    stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)
        finally:
            release.set()
            holder.join(timeout=5)

        assert stats["successful"] == 1
        result = stats["results"]["SRR1"]
        assert result.endswith("; stored")
        assert "catalogue pending" in result
        assert (fastq_folder / "SRR1").is_symlink()
        assert "store_reindex" in caplog.text
        assert str(paths.root) in caplog.text

        with Catalog(paths) as catalog:
            assert catalog.get_dataset("SRR1") is None

        sidecar = read_sidecar(paths.sra / "SRR1" / "SRR1.json")
        assert sidecar is not None
        with catalog_write(paths) as catalog:
            catalog.reindex([sidecar])
        with Catalog(paths) as catalog:
            assert catalog.get_dataset("SRR1") is not None

    def test_gzipped_store_download_records_gzip_compression(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        fastq_folder = tmp_path / "project" / "fastq"

        def _download(accession, output_folder, *args, **kwargs):
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True, exist_ok=True)
            with gzip.open(acc_dir / f"{accession}_1.fastq.gz", "wt") as handle:
                handle.write("@r\nACGT\n+\nIIII\n")
            return True, "Downloaded 1 files, unverified"

        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert read_sidecar(paths.sra / "SRR1" / "SRR1.json").compression == "gzip"

    def test_store_download_uses_metadata_xml_for_expected_spots(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        fastq_folder = tmp_path / "project" / "fastq"
        project_metadata = tmp_path / "project" / "metadata"
        project_metadata.mkdir(parents=True)
        (project_metadata / "SRR1_metadata.xml").write_text(
            '<RunSet><RUN total_spots="1" total_bases="4" size="10"/></RunSet>'
        )
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            download_sra(
                fastq_folder,
                self._accessions(tmp_path, "SRR1"),
                store=paths,
                store_metadata=project_metadata,
                max_retries=0,
            )

        sidecar = read_sidecar(paths.sra / "SRR1" / "SRR1.json")
        assert sidecar.ncbi["spots"] == 1
        assert sidecar.completeness["verdict"] == "complete"
        assert sidecar.state == "complete"

    def test_partial_in_store_is_redownloaded_when_resume_partial(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        self._store_dataset(paths, state="partial")
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert len(calls) == 1
        # The replacement is built beside the store, then moved over the partial copy.
        assert calls[0][1] == paths.tmp
        assert stats["successful"] == 1
        assert read_sidecar(paths.sra / "SRR1" / "SRR1.json").state == "complete"
        assert (fastq_folder / "SRR1").is_symlink()
        assert not (paths.tmp / "SRR1").exists()

    def test_failed_in_store_is_redownloaded_and_clears_the_error(self, tmp_path):
        """A dataset whose sidecar recorded a prior failure (state="failed", with an error
        message from that failed attempt) is treated as incomplete just like "partial" and
        redownloaded. build_sidecar always constructs a fresh Sidecar for the new attempt, so
        a successful redownload must not carry the stale error forward: the published sidecar
        must land on state="complete" with error=None."""
        from prepare.store.sidecar import Sidecar, read_sidecar, write_sidecar

        paths = self._store(tmp_path)
        acc_dir = self._store_dataset(paths, state="failed")
        write_sidecar(
            acc_dir / "SRR1.json",
            Sidecar(accession="SRR1", state="failed", error="fasterq-dump exit 3"),
        )
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert len(calls) == 1
        assert calls[0][1] == paths.tmp
        assert stats["successful"] == 1
        sidecar = read_sidecar(paths.sra / "SRR1" / "SRR1.json")
        assert sidecar.state == "complete"
        assert sidecar.error is None
        assert (fastq_folder / "SRR1").is_symlink()

    def test_files_without_a_sidecar_are_redownloaded(self, tmp_path):
        paths = self._store(tmp_path)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert len(calls) == 1

    def test_partial_in_store_is_refused_without_resume_partial(self, tmp_path):
        paths = self._store(tmp_path)
        self._store_dataset(paths, state="partial")
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            stats = download_sra(
                fastq_folder,
                self._accessions(tmp_path, "SRR1"),
                store=paths,
                resume_partial=False,
                max_retries=0,
            )

        assert calls == []
        assert stats["failed"] == 1
        assert "rerun with --resume-partial" in stats["results"]["SRR1"]
        assert not (fastq_folder / "SRR1").exists()

    def test_accept_partial_links_the_partial_copy(self, tmp_path):
        paths = self._store(tmp_path)
        self._store_dataset(paths, state="partial")
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            stats = download_sra(
                fastq_folder,
                self._accessions(tmp_path, "SRR1"),
                store=paths,
                resume_partial=False,
                accept_partial=True,
                max_retries=0,
            )

        assert calls == []
        assert stats["successful"] == 1
        assert stats["results"]["SRR1"].startswith("linked from store")
        assert (fastq_folder / "SRR1").is_symlink()

    def test_failed_store_download_keeps_the_partial_copy(self, tmp_path):
        paths = self._store(tmp_path)
        self._store_dataset(paths, state="partial")
        fastq_folder = tmp_path / "project" / "fastq"

        with patch(
            "prepare.data.sra.accession.download_accession", return_value=(False, "network: Download failed: timeout")
        ):
            stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert stats["failed"] == 1
        assert (paths.sra / "SRR1" / "SRR1_1.fastq").exists()
        assert not (fastq_folder / "SRR1").exists()

    @staticmethod
    def _fake_reads_download(calls, reads, expected=None):
        """A download_accession stand-in writing ``reads`` FASTQ records under <output_folder>/<acc>."""

        def _download(accession, output_folder, *args, **kwargs):
            calls.append((accession, Path(output_folder), kwargs))
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True, exist_ok=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n" * reads)
            if expected is None:
                return True, f"Downloaded 1 files, unverified ({reads} reads)"
            verdict = "complete" if reads >= expected else "truncated"
            return True, f"Downloaded 1 files, {verdict} ({reads} of {expected} spots)"

        return _download

    INCOMPLETE_SRR1 = (
        "incomplete: SRR1 store copy holds 1 of 100 spots; kept for --resume-partial; "
        "rerun with --accept-partial to use it"
    )

    def test_partial_store_download_is_not_linked_by_default(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch(
            "prepare.data.sra.accession.download_accession", side_effect=self._fake_reads_download(calls, 1, 100)
        ):
            stats = download_sra(
                fastq_folder,
                self._accessions(tmp_path, "SRR1"),
                store=paths,
                expected_spots={"SRR1": 100},
                max_retries=0,
            )

        assert len(calls) == 1
        assert stats["failed"] == 1
        assert stats["results"]["SRR1"] == self.INCOMPLETE_SRR1
        # Published for a later --resume-partial, but never presented to the project.
        sidecar = read_sidecar(paths.sra / "SRR1" / "SRR1.json")
        assert sidecar.state == "partial"
        assert sidecar.reads_per_mate == 1
        assert not (fastq_folder / "SRR1").exists()
        assert not (fastq_folder / "SRR1").is_symlink()

    def test_partial_store_download_is_linked_with_accept_partial(self, tmp_path):
        paths = self._store(tmp_path)
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch(
            "prepare.data.sra.accession.download_accession", side_effect=self._fake_reads_download(calls, 1, 100)
        ):
            stats = download_sra(
                fastq_folder,
                self._accessions(tmp_path, "SRR1"),
                store=paths,
                expected_spots={"SRR1": 100},
                accept_partial=True,
                max_retries=0,
            )

        assert stats["successful"] == 1
        assert stats["results"]["SRR1"] == "Downloaded 1 files, truncated (1 of 100 spots); stored"
        assert (fastq_folder / "SRR1").is_symlink()

    def test_partial_refetch_removes_an_existing_project_link(self, tmp_path):
        from prepare.store.link import link_dataset

        paths = self._store(tmp_path)
        self._store_dataset(paths, state="partial")
        fastq_folder = tmp_path / "project" / "fastq"
        link_dataset(fastq_folder, "SRR1", paths)
        calls = []

        with patch(
            "prepare.data.sra.accession.download_accession", side_effect=self._fake_reads_download(calls, 1, 100)
        ):
            stats = download_sra(
                fastq_folder,
                self._accessions(tmp_path, "SRR1"),
                store=paths,
                expected_spots={"SRR1": 100},
                max_retries=0,
            )

        assert len(calls) == 1
        assert stats["results"]["SRR1"] == self.INCOMPLETE_SRR1
        assert not (fastq_folder / "SRR1").is_symlink()
        assert (paths.sra / "SRR1" / "SRR1_1.fastq").exists()

    def test_partial_store_fetch_never_removes_a_real_project_folder(self, tmp_path):
        from prepare.store.locks import dataset_lock

        paths = self._store(tmp_path)
        for directory in (paths.sra, paths.tmp, paths.locks):
            directory.mkdir(parents=True, exist_ok=True)
        fastq_folder = tmp_path / "project" / "fastq"
        own = fastq_folder / "SRR1"
        own.mkdir(parents=True)
        (own / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        calls = []

        with patch(
            "prepare.data.sra.accession.download_accession", side_effect=self._fake_reads_download(calls, 1, 100)
        ):
            # _store_download holds the dataset lock around every _store_fetch call.
            with dataset_lock(paths, "SRR1"):
                success, message = store_handoff_mod._store_fetch(
                    "SRR1", fastq_folder, paths, "auto", None, expected_spots=100
                )

        assert success is False
        assert message == self.INCOMPLETE_SRR1
        assert own.is_dir() and not own.is_symlink()
        assert (own / "SRR1_1.fastq").exists()

    def test_complete_previous_copy_is_kept_when_the_refetch_is_short(self, tmp_path):
        from prepare.store.sidecar import Sidecar, read_sidecar, write_sidecar

        paths = self._store(tmp_path)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n" * 2)
        write_sidecar(
            acc_dir / "SRR1.json",
            Sidecar(
                accession="SRR1",
                state="complete",
                reads_per_mate=2,
                ncbi={"spots": 2},
                completeness={"method": "spots", "ratio": 1.0, "verdict": "complete"},
            ),
        )
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_reads_download(calls, 1, 2)):
            stats = download_sra(
                fastq_folder,
                self._accessions(tmp_path, "SRR1"),
                store=paths,
                expected_spots={"SRR1": 2},
                force=True,
                max_retries=0,
            )

        assert len(calls) == 1
        assert stats["successful"] == 1
        assert stats["results"]["SRR1"].startswith("linked from store")
        sidecar = read_sidecar(acc_dir / "SRR1.json")
        assert sidecar.state == "complete"
        assert sidecar.reads_per_mate == 2
        assert (acc_dir / "SRR1_1.fastq").read_text().count("@r") == 2
        assert (fastq_folder / "SRR1").is_symlink()
        assert not (paths.tmp / "SRR1").exists()

    def test_store_sidecar_takes_registry_spots_when_the_xml_has_none(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_reads_download(calls, 1, 1)):
            stats = download_sra(
                fastq_folder,
                self._accessions(tmp_path, "SRR1"),
                store=paths,
                expected_spots={"SRR1": 1},
                max_retries=0,
            )

        assert stats["successful"] == 1
        sidecar = read_sidecar(paths.sra / "SRR1" / "SRR1.json")
        assert sidecar.ncbi["spots"] == 1
        assert sidecar.completeness["verdict"] == "complete"
        assert sidecar.state == "complete"

    @staticmethod
    def _short_unverified_dataset(paths):
        from prepare.store.sidecar import Sidecar, write_sidecar

        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        write_sidecar(
            acc_dir / "SRR1.json",
            Sidecar(
                accession="SRR1",
                state="unverified",
                reads_per_mate=1,
                ncbi={},
                completeness={"method": "unverified", "ratio": None, "verdict": "unverified"},
            ),
        )
        return acc_dir

    def test_short_unverified_store_copy_is_refetched_with_the_sidecar_rewritten(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        acc_dir = self._short_unverified_dataset(paths)
        fastq_folder = tmp_path / "project" / "fastq"
        seen = []
        fake = self._fake_reads_download([], 100, 100)

        def _download(accession, output_folder, *args, **kwargs):
            # The sidecar is rewritten under the lock before the refetch starts.
            seen.append(read_sidecar(acc_dir / "SRR1.json"))
            return fake(accession, output_folder, *args, **kwargs)

        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            stats = download_sra(
                fastq_folder,
                self._accessions(tmp_path, "SRR1"),
                store=paths,
                expected_spots={"SRR1": 100},
                max_retries=0,
            )

        assert len(seen) == 1
        assert seen[0].state == "partial"
        assert seen[0].ncbi["spots"] == 100
        assert seen[0].completeness == {"method": "spots", "ratio": 0.01, "verdict": "truncated"}
        assert stats["successful"] == 1
        assert read_sidecar(acc_dir / "SRR1.json").state == "complete"
        assert (fastq_folder / "SRR1").is_symlink()

    def test_short_unverified_store_copy_is_not_linked_without_resume(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        acc_dir = self._short_unverified_dataset(paths)
        fastq_folder = tmp_path / "project" / "fastq"

        with patch("prepare.data.sra.accession.download_accession", side_effect=AssertionError("no download expected")):
            stats = download_sra(
                fastq_folder,
                self._accessions(tmp_path, "SRR1"),
                store=paths,
                expected_spots={"SRR1": 100},
                resume_partial=False,
                max_retries=0,
            )

        assert stats["failed"] == 1
        assert "rerun with --resume-partial" in stats["results"]["SRR1"]
        assert read_sidecar(acc_dir / "SRR1.json").state == "partial"
        assert not (fastq_folder / "SRR1").exists()

    def test_store_state_counts_a_short_unverified_copy_as_incomplete(self, tmp_path):
        paths = self._store(tmp_path)
        self._short_unverified_dataset(paths)

        assert store_handoff_mod._store_state(paths, "SRR1") == "ready"
        assert store_handoff_mod._store_state(paths, "SRR1", expected_spots=1) == "ready"
        assert store_handoff_mod._store_state(paths, "SRR1", expected_spots=100) == "incomplete"

    def test_incomplete_message_classifies_as_unknown(self):
        assert classify_download_error(self.INCOMPLETE_SRR1) == "unknown"

    @pytest.mark.parametrize(
        "previous, new, expected",
        [
            (None, ("partial", 1), "publish"),
            (("partial", 1), ("partial", 5), "publish"),
            (("partial", 5), ("partial", 1), "keep_previous"),
            (("complete", 100), ("partial", 50), "keep_previous"),
            (("complete", 100), ("complete", 100), "publish"),
            (("partial", 5), ("failed", None), "keep_previous"),
            (("complete", 100), ("failed", None), "keep_previous"),
            (("failed", None), ("partial", 1), "publish"),
            (("failed", None), ("failed", None), "publish"),
            (("partial", 50), ("complete", 99), "publish"),
        ],
    )
    def test_publish_decision(self, previous, new, expected):
        from prepare.store.sidecar import Sidecar

        def _sidecar(spec):
            if spec is None:
                return None
            state, reads = spec
            return Sidecar(accession="SRR1", state=state, reads_per_mate=reads, ncbi={"spots": 100})

        assert store_handoff_mod._publish_decision(_sidecar(previous), _sidecar(new)) == expected

    # Task 6: the refetch limit on a partial store copy.

    EXHAUSTED_SRR1 = (
        "incomplete: 2 refetches of SRR1 gained no reads (1 of 100 spots); NCBI's count may not be reachable; "
        "use --accept-partial, or --force to fetch again"
    )

    @staticmethod
    def _partial_dataset(paths, reads, refetch=None, state="partial"):
        from prepare.store.sidecar import Sidecar, write_sidecar

        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True, exist_ok=True)
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n" * reads)
        write_sidecar(
            acc_dir / "SRR1.json",
            Sidecar(
                accession="SRR1",
                state=state,
                reads_per_mate=reads,
                ncbi={"spots": 100},
                completeness={"method": "spots", "ratio": reads / 100, "verdict": "truncated"},
                refetch=refetch,
            ),
        )
        return acc_dir

    def _run_partial(self, tmp_path, paths, calls, reads, **kwargs):
        with patch(
            "prepare.data.sra.accession.download_accession", side_effect=self._fake_reads_download(calls, reads, 100)
        ):
            return download_sra(
                tmp_path / "project" / "fastq",
                self._accessions(tmp_path, "SRR1"),
                store=paths,
                expected_spots={"SRR1": 100},
                max_retries=0,
                **kwargs,
            )

    def test_third_identical_partial_refetch_does_not_download(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        calls = []
        # The first download and two refetches that each gain nothing.
        for _ in range(3):
            stats = self._run_partial(tmp_path, paths, calls, 1)
            assert stats["results"]["SRR1"] == self.INCOMPLETE_SRR1
        assert len(calls) == 3
        record = read_sidecar(paths.sra / "SRR1" / "SRR1.json").refetch
        assert record["unchanged"] == 2
        assert record["reads_per_mate"] == 1
        assert record["expected_spots"] == 100
        assert record["last"]

        stats = self._run_partial(tmp_path, paths, calls, 1)

        assert len(calls) == 3
        assert stats["failed"] == 1
        assert stats["results"]["SRR1"] == self.EXHAUSTED_SRR1
        assert classify_download_error(self.EXHAUSTED_SRR1) == "unknown"
        assert not (tmp_path / "project" / "fastq" / "SRR1").exists()

    def test_an_exhausted_copy_is_linked_as_partial_with_accept_partial(self, tmp_path):
        paths = self._store(tmp_path)
        self._partial_dataset(paths, 1, refetch={"unchanged": 2, "reads_per_mate": 1, "expected_spots": 100})
        calls = []

        stats = self._run_partial(tmp_path, paths, calls, 1, accept_partial=True)

        assert calls == []
        assert stats["successful"] == 1
        assert stats["results"]["SRR1"] == "linked from store (partial), 1 files"
        assert (tmp_path / "project" / "fastq" / "SRR1").is_symlink()

    def test_force_refetches_an_exhausted_copy(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        self._partial_dataset(paths, 1, refetch={"unchanged": 2, "reads_per_mate": 1, "expected_spots": 100})
        calls = []

        stats = self._run_partial(tmp_path, paths, calls, 1, force=True)

        assert len(calls) == 1
        assert stats["results"]["SRR1"] == self.INCOMPLETE_SRR1
        assert read_sidecar(paths.sra / "SRR1" / "SRR1.json").refetch["unchanged"] == 3

    def test_a_refetch_that_gains_reads_resets_the_count(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        self._partial_dataset(paths, 1, refetch={"unchanged": 1, "reads_per_mate": 1, "expected_spots": 100})
        calls = []

        stats = self._run_partial(tmp_path, paths, calls, 2)

        assert len(calls) == 1
        assert stats["failed"] == 1
        sidecar = read_sidecar(paths.sra / "SRR1" / "SRR1.json")
        assert sidecar.reads_per_mate == 2
        assert sidecar.refetch["unchanged"] == 0
        assert sidecar.refetch["reads_per_mate"] == 2
        assert not store_handoff_mod._refetch_exhausted(sidecar)

    def test_a_complete_refetch_clears_the_record(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        acc_dir = self._partial_dataset(paths, 1, refetch={"unchanged": 2, "reads_per_mate": 1, "expected_spots": 100})
        calls = []

        stats = self._run_partial(tmp_path, paths, calls, 100, force=True)

        assert len(calls) == 1
        assert stats["successful"] == 1
        sidecar = read_sidecar(acc_dir / "SRR1.json")
        assert sidecar.state == "complete"
        assert sidecar.refetch is None
        assert "refetch" not in json.loads((acc_dir / "SRR1.json").read_text())

    def test_a_worse_refetch_counts_against_the_kept_copy(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        acc_dir = self._partial_dataset(paths, 5)
        calls = []

        stats = self._run_partial(tmp_path, paths, calls, 1)

        assert len(calls) == 1
        assert stats["results"]["SRR1"].startswith("incomplete: SRR1 store copy holds 5 of 100 spots")
        sidecar = read_sidecar(acc_dir / "SRR1.json")
        assert sidecar.reads_per_mate == 5
        assert sidecar.refetch["unchanged"] == 1
        assert sidecar.refetch["reads_per_mate"] == 5
        assert (acc_dir / "SRR1_1.fastq").read_text().count("@r") == 5

    @pytest.mark.parametrize(
        "previous, new, expected",
        [
            (None, ("partial", 1), None),
            (("partial", 1, None), ("partial", 1), 1),
            (("partial", 1, 1), ("partial", 1), 2),
            (("partial", 5, 1), ("partial", 1), 2),
            (("partial", 1, 3), ("partial", 2), 0),
            (("partial", 1, 3), ("failed", None), 4),
            (("partial", 1, 3), ("complete", 100), None),
            # A count of 0 is known: a refetch that again returns no reads is no gain.
            (("partial", 0, 1), ("partial", 0), 2),
            (("partial", 0, None), ("partial", 0), 1),
            (("partial", 0, 1), ("partial", 1), 0),
        ],
    )
    def test_refetch_record(self, previous, new, expected):
        from prepare.store.sidecar import Sidecar

        prior = None
        if previous is not None:
            state, reads, unchanged = previous
            record = None if unchanged is None else {"unchanged": unchanged, "reads_per_mate": reads}
            prior = Sidecar(accession="SRR1", state=state, reads_per_mate=reads, ncbi={"spots": 100}, refetch=record)
        fresh = Sidecar(accession="SRR1", state=new[0], reads_per_mate=new[1], ncbi={"spots": 100})

        record = store_handoff_mod._refetch_record(prior, fresh)

        if expected is None:
            assert record is None
        else:
            assert record["unchanged"] == expected
            assert record["expected_spots"] == 100
            assert record["last"]

    def test_refetch_exhausted(self):
        from prepare.store.sidecar import Sidecar

        limit = store_handoff_mod.STORE_PARTIAL_REFETCH_LIMIT
        assert limit == 2
        assert not store_handoff_mod._refetch_exhausted(None)
        assert not store_handoff_mod._refetch_exhausted(Sidecar(state="partial"))
        assert not store_handoff_mod._refetch_exhausted(Sidecar(state="partial", refetch={"unchanged": limit - 1}))
        assert store_handoff_mod._refetch_exhausted(Sidecar(state="partial", refetch={"unchanged": limit}))

    # Task 5 review follow-ups.

    def test_forced_worse_refetch_never_links_a_short_unverified_copy(self, tmp_path):
        from prepare.store.sidecar import Sidecar, read_sidecar, write_sidecar

        paths = self._store(tmp_path)
        acc_dir = paths.sra / "SRR1"
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n" * 5)
        write_sidecar(
            acc_dir / "SRR1.json",
            Sidecar(
                accession="SRR1",
                state="unverified",
                reads_per_mate=5,
                ncbi={},
                completeness={"method": "unverified", "ratio": None, "verdict": "unverified"},
            ),
        )
        calls = []

        stats = self._run_partial(tmp_path, paths, calls, 1, force=True)

        assert len(calls) == 1
        assert stats["failed"] == 1
        assert stats["results"]["SRR1"].startswith("incomplete: SRR1 store copy holds 5 of 100 spots")
        sidecar = read_sidecar(acc_dir / "SRR1.json")
        assert sidecar.state == "partial"
        assert sidecar.reads_per_mate == 5
        assert not (tmp_path / "project" / "fastq" / "SRR1").exists()

    def test_precheck_refusal_removes_a_store_link(self, tmp_path):
        from prepare.store.link import link_dataset

        paths = self._store(tmp_path)
        self._store_dataset(paths, state="partial")
        fastq_folder = tmp_path / "project" / "fastq"
        link_dataset(fastq_folder, "SRR1", paths)

        with patch("prepare.data.sra.accession.download_accession", side_effect=AssertionError("no download expected")):
            stats = download_sra(
                fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, resume_partial=False, max_retries=0
            )

        assert "rerun with --resume-partial" in stats["results"]["SRR1"]
        assert not (fastq_folder / "SRR1").is_symlink()
        assert (paths.sra / "SRR1" / "SRR1_1.fastq").exists()

    def test_precheck_refusal_leaves_a_link_outside_the_store(self, tmp_path):
        paths = self._store(tmp_path)
        self._store_dataset(paths, state="partial")
        fastq_folder = tmp_path / "project" / "fastq"
        fastq_folder.mkdir(parents=True)
        elsewhere = tmp_path / "elsewhere" / "SRR1"
        elsewhere.mkdir(parents=True)
        (elsewhere / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (fastq_folder / "SRR1").symlink_to(elsewhere)

        with patch("prepare.data.sra.accession.download_accession", side_effect=AssertionError("no download expected")):
            download_sra(
                fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, resume_partial=False, max_retries=0
            )

        assert (fastq_folder / "SRR1").is_symlink()
        assert (fastq_folder / "SRR1").resolve() == elsewhere.resolve()

    def test_space_guard_counts_a_short_unverified_copy_due_for_refetch(self, tmp_path):
        from prepare.data.sra import download as download_mod

        paths = self._store(tmp_path)
        self._short_unverified_dataset(paths)
        args = (tmp_path / "project" / "fastq", None, None, paths, True, {"SRR1": 1000}, 1.0, False, ["SRR1"])

        assert download_mod._space_guard(*args).needs("SRR1") == {}
        assert download_mod._space_guard(*args, expected_spots={"SRR1": 100}).needs("SRR1") != {}

    def test_space_guard_exempts_an_exhausted_copy(self, tmp_path):
        from prepare.data.sra import download as download_mod

        paths = self._store(tmp_path)
        self._partial_dataset(paths, 1, refetch={"unchanged": 2, "reads_per_mate": 1, "expected_spots": 100})
        args = (tmp_path / "project" / "fastq", None, None, paths, True, {"SRR1": 1000}, 1.0)

        assert download_mod._space_guard(*args, False, ["SRR1"]).needs("SRR1") == {}
        assert download_mod._space_guard(*args, True, ["SRR1"]).needs("SRR1") != {}

    def test_second_project_links_without_downloading_again(self, tmp_path):
        paths = self._store(tmp_path)
        first = tmp_path / "project_a" / "fastq"
        second = tmp_path / "project_b" / "fastq"
        acc = self._accessions(tmp_path, "SRR1")
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            download_sra(first, acc, store=paths, max_retries=0)
            download_sra(second, acc, store=paths, max_retries=0)

        assert len(calls) == 1
        assert (first / "SRR1").resolve() == (second / "SRR1").resolve() == (paths.sra / "SRR1").resolve()

    def test_link_mode_copy_copies_instead_of_linking(self, tmp_path):
        paths = self._store(tmp_path)
        self._store_dataset(paths)
        fastq_folder = tmp_path / "project" / "fastq"

        with patch("prepare.data.sra.accession.download_accession", side_effect=AssertionError("no download expected")):
            download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, link_mode="copy", max_retries=0)

        assert (fastq_folder / "SRR1").is_dir() and not (fastq_folder / "SRR1").is_symlink()

    def test_an_existing_project_link_counts_as_already_downloaded(self, tmp_path):
        paths = self._store(tmp_path)
        self._store_dataset(paths)
        fastq_folder = tmp_path / "project" / "fastq"
        from prepare.store.link import link_dataset

        link_dataset(fastq_folder, "SRR1", paths)

        with patch("prepare.data.sra.accession.download_accession", side_effect=AssertionError("no download expected")):
            stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert stats["already_downloaded_accessions"] == ["SRR1"]

    def test_a_link_to_a_partial_dataset_is_not_already_downloaded(self, tmp_path):
        paths = self._store(tmp_path)
        self._store_dataset(paths, state="partial")
        fastq_folder = tmp_path / "project" / "fastq"
        from prepare.store.link import link_dataset

        link_dataset(fastq_folder, "SRR1", paths)
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert stats["already_downloaded_accessions"] == []
        assert len(calls) == 1

    def test_without_a_store_the_project_folder_is_used(self, tmp_path):
        fastq_folder = tmp_path / "fastq"
        calls = []

        with patch("prepare.data.sra.accession.download_accession", side_effect=self._fake_download(calls)):
            stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), max_retries=0)

        # Built in the project's own staging folder, then published into fastq/SRR1.
        assert calls[0][1] == fastq_folder / ".prepare-tmp"
        assert accession_has_fastq(fastq_folder / "SRR1")
        assert not (fastq_folder / "SRR1").is_symlink()
        assert stats["results"]["SRR1"] == "Downloaded 1 files, unverified"

    @staticmethod
    def _fake_download_using_scratch(calls):
        """A download_accession stand-in that actually creates kwargs["temp_folder"] and
        drops a marker file in it before returning, so a test asserting the folder is gone
        afterwards is only satisfied if something really removed it.
        """

        def _download(accession, output_folder, *args, **kwargs):
            calls.append((accession, Path(output_folder), kwargs))
            scratch = Path(kwargs["temp_folder"])
            scratch.mkdir(parents=True, exist_ok=True)
            (scratch / "marker").write_text("fasterq-dump scratch")
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True, exist_ok=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            return True, "Downloaded 1 files, unverified"

        return _download

    @staticmethod
    def _fake_download_using_scratch_then_failing(calls):
        """Like ``_fake_download_using_scratch``, but raises after creating the scratch folder,
        so the removal path exercised is _store_fetch's ``finally`` block, not the success path.
        """

        def _download(accession, output_folder, *args, **kwargs):
            calls.append((accession, Path(output_folder), kwargs))
            scratch = Path(kwargs["temp_folder"])
            scratch.mkdir(parents=True, exist_ok=True)
            (scratch / "marker").write_text("fasterq-dump scratch")
            raise OSError("fasterq-dump exploded")

        return _download

    def test_store_download_uses_store_tmp_for_fasterq_scratch(self, tmp_path):
        """fasterq-dump's scratch folder is created under the store's tmp, not the system temp,
        and is removed again once the download finishes, since download_accession only cleans
        up a temp_folder it created itself.
        """
        paths = self._store(tmp_path)
        calls: list = []

        with patch(
            "prepare.data.sra.accession.download_accession", side_effect=self._fake_download_using_scratch(calls)
        ):
            download_sra(tmp_path / "fastq", self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        kwargs = calls[0][2]
        assert Path(kwargs["temp_folder"]).parent == paths.tmp
        assert Path(kwargs["temp_folder"]).name == "SRR1_fqtmp"
        assert not Path(kwargs["temp_folder"]).exists()

    def test_store_download_removes_fasterq_scratch_after_a_failure(self, tmp_path):
        """A download_accession that creates the scratch folder and then fails must not leak
        it: _store_fetch's ``finally`` block removes it regardless of success or failure.
        """
        paths = self._store(tmp_path)
        calls: list = []

        with patch(
            "prepare.data.sra.accession.download_accession",
            side_effect=self._fake_download_using_scratch_then_failing(calls),
        ):
            stats = download_sra(tmp_path / "fastq", self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        kwargs = calls[0][2]
        assert stats["failed"] == 1
        assert not Path(kwargs["temp_folder"]).exists()


class TestStoreDownloadLockWait:
    """download_sra's --lock-wait: give up on an accession another project is downloading."""

    @staticmethod
    def _store(tmp_path):
        from prepare.store.layout import init_store

        return init_store(tmp_path / "store")

    def test_lock_wait_gives_up_naming_the_accession_and_holder(self, tmp_path):
        import json as json_module

        from prepare.store.layout import lock_path

        paths = self._store(tmp_path)
        lock = lock_path(paths, "SRR1")
        lock.parent.mkdir(parents=True, exist_ok=True)
        lock.write_text(json_module.dumps({"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00"}))
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        fastq_folder = tmp_path / "project" / "fastq"

        with patch("prepare.data.sra.accession.download_accession") as mock_download:
            stats = download_sra(
                fastq_folder, acc, store=paths, max_retries=0, lock_wait=0.1, num_threads=1, max_workers=1
            )

        mock_download.assert_not_called()
        assert stats["failed"] == 1
        message = stats["results"]["SRR1"]
        assert "SRR1" in message
        assert "4242" in message and "otherhost" in message
        # The other project's lock is left exactly as it was.
        assert json_module.loads(lock.read_text())["pid"] == 4242


class TestStoreDownloadPublishesAtomically:
    """sra/<ACC> only ever appears complete: verified, compressed and described."""

    @staticmethod
    def _store(tmp_path):
        from prepare.store.layout import init_store

        return init_store(tmp_path / "store")

    @staticmethod
    def _accessions(tmp_path, *accessions):
        acc = tmp_path / "acc.txt"
        acc.write_text("".join(f"{a}\n" for a in accessions))
        return acc

    def test_the_store_folder_is_absent_until_the_dataset_is_described(self, tmp_path):
        from prepare.store.sidecar import read_sidecar

        paths = self._store(tmp_path)
        fastq_folder = tmp_path / "project" / "fastq"
        seen = {}

        def _download(accession, output_folder, *args, **kwargs):
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True, exist_ok=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            # The files exist now, unverified and undescribed: the published location must
            # still be empty, so no other project can see or claim a half-finished dataset.
            seen["published_during_download"] = (paths.sra / accession).exists()
            seen["output_folder"] = Path(output_folder)
            return True, "Downloaded 1 files, unverified"

        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert seen["published_during_download"] is False
        assert seen["output_folder"] == paths.tmp
        assert stats["successful"] == 1
        # The folder is published with its sidecar already inside it.
        assert sorted(p.name for p in (paths.sra / "SRR1").iterdir()) == ["SRR1.json", "SRR1_1.fastq"]
        assert read_sidecar(paths.sra / "SRR1" / "SRR1.json") is not None
        assert not (paths.tmp / "SRR1").exists()

    def test_a_failure_before_publish_never_creates_the_store_folder(self, tmp_path):
        paths = self._store(tmp_path)
        fastq_folder = tmp_path / "project" / "fastq"

        def _download(accession, output_folder, *args, **kwargs):
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True, exist_ok=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            return True, "Downloaded 1 files, unverified"

        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            with patch("prepare.store.sidecar.build_sidecar", side_effect=OSError("disk went away")):
                stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert stats["failed"] == 1
        assert not (paths.sra / "SRR1").exists()
        # The half-finished work stays in the store's tmp folder, out of the published tree.
        assert (paths.tmp / "SRR1").is_dir()
        assert not (fastq_folder / "SRR1").exists()

    def test_a_lock_taken_over_during_the_download_publishes_nothing(self, tmp_path):
        import json as json_module

        from prepare.store.layout import lock_path

        paths = self._store(tmp_path)
        fastq_folder = tmp_path / "project" / "fastq"
        other = {"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00", "token": "feedbeef"}

        def _download(accession, output_folder, *args, **kwargs):
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True, exist_ok=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            # Another project reclaimed the lock while this holder was stalled.
            lock_path(paths, accession).write_text(json_module.dumps(other))
            return True, "Downloaded 1 files, unverified"

        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            stats = download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert stats["failed"] == 1
        assert stats["results"]["SRR1"].startswith("lock lost: ")
        assert not (paths.sra / "SRR1").exists()
        assert not (fastq_folder / "SRR1").exists()
        # The staged folder and the lock now belong to the other holder and are left alone.
        assert (paths.tmp / "SRR1").is_dir()
        assert json_module.loads(lock_path(paths, "SRR1").read_text())["token"] == "feedbeef"

    def test_a_failed_replacement_leaves_the_existing_copy_in_place(self, tmp_path):
        from prepare.store.sidecar import Sidecar, write_sidecar

        paths = self._store(tmp_path)
        existing = paths.sra / "SRR1"
        existing.mkdir(parents=True)
        (existing / "SRR1_1.fastq").write_text("@old\nACGT\n+\nIIII\n")
        write_sidecar(existing / "SRR1.json", Sidecar(accession="SRR1", state="partial"))
        fastq_folder = tmp_path / "project" / "fastq"

        def _download(accession, output_folder, *args, **kwargs):
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True, exist_ok=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@new\nACGT\n+\nIIII\n")
            return True, "Downloaded 1 files, unverified"

        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            with patch("prepare.store.sidecar.build_sidecar", side_effect=OSError("disk went away")):
                download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0)

        assert (existing / "SRR1_1.fastq").read_text().startswith("@old")
        assert (existing / "SRR1.json").is_file()

    def test_a_forced_refetch_of_an_absent_dataset_stays_forced(self, tmp_path):
        paths = self._store(tmp_path)
        fastq_folder = tmp_path / "project" / "fastq"
        calls = []

        def _download(accession, output_folder, *args, **kwargs):
            calls.append(kwargs)
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True, exist_ok=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            return True, "Downloaded 1 files, unverified"

        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            download_sra(fastq_folder, self._accessions(tmp_path, "SRR1"), store=paths, max_retries=0, force=True)

        assert calls[0]["force"] is True


def test_summary_counts_a_store_link_made_by_the_retry_pass_as_a_link(tmp_path, caplog):
    """A retry result ("Retry 1: linked from store ...") is a store link, not a new download."""
    from prepare.data.sra.download import _log_download_run_summary

    with caplog.at_level("INFO", logger="prepare.data.sra"):
        _log_download_run_summary(
            all_accessions=["SRR1"],
            already_downloaded=[],
            blacklisted=[],
            successful_count=1,
            failed_count=0,
            download_results={"SRR1": "Retry 1: linked from store (complete)"},
            abort_reason=None,
            failed_accessions=[],
            fastq_path=tmp_path,
        )

    assert "Newly downloaded: 0" in caplog.text
    assert "Linked from store: 1" in caplog.text
