"""Tests for prepare.data.sra: reading, verifying and downloading accessions, and project-level downloads."""

import gzip
import inspect
import json
import logging
import os
import shutil
import subprocess
import threading
import time
import pytest
from pathlib import Path
from unittest.mock import Mock, patch
from prepare.core.exceptions import DataAccessError, SecurityError, ValidationError
from prepare.data.sra.fastq import (
    accession_has_fastq,
    compress_fastq,
    count_fastq_reads,
    fastq_digest,
    fastq_files,
    iter_fastq_records,
    orphan_fastq,
    parse_verdict_message,
    primary_fastq,
    verify_download,
)
from prepare.data.sra.cleanup import (
    _ignore_missing,
    _prepare_temp_folder,
    _safe_rmtree,
    is_transient_folder,
    transient_bytes,
)
from prepare.data.sra.accession import (
    _cached_sra_archive,
    _check_existing_download,
    _handle_download_output,
    classify_download_error,
    download_accession,
    fasterq_dump_version,
)
from prepare.data.sra.download import (
    _check_existing_downloads,
    _read_blacklist_files,
    default_max_workers,
    download_sra,
)
from prepare.data.sra import accession as accession_mod
from prepare.data.sra import fastq as fastq_mod
from prepare.data.sra import retry as retry_mod
from prepare.data.sra import store_handoff as store_handoff_mod
from prepare.utils.security import SecureSubprocess
from helpers_extraction import _fake_tools
from _sra_fixtures import (
    _staging_download,
    _slow_download_tools,
    _no_prefetch_with_pigz,
    _with_prefetch,
    _scratch_recording_tools,
)


class TestReadBlacklistFiles:
    """Test _read_blacklist_files function."""

    def test_read_blacklist_files_success(self, tmp_path):
        """Test successful blacklist file reading."""
        blacklist1 = tmp_path / "blacklist1.txt"
        blacklist2 = tmp_path / "blacklist2.txt"

        blacklist1.write_text("SRR123\nSRR456\n\n")  # Include empty line
        blacklist2.write_text("SRR789\nSRR123\n")  # Include duplicate

        result = _read_blacklist_files([blacklist1, blacklist2])

        assert result == {"SRR123", "SRR456", "SRR789"}

    def test_read_blacklist_files_empty_list(self):
        """Test with empty blacklist files list."""
        result = _read_blacklist_files([])
        assert result == set()

    def test_read_blacklist_files_none(self):
        """Test with None blacklist files."""
        result = _read_blacklist_files(None)
        assert result == set()

    def test_read_blacklist_files_nonexistent_file(self, tmp_path):
        """Test handling nonexistent blacklist file."""
        nonexistent = tmp_path / "nonexistent.txt"

        with patch("prepare.data.sra.download.logger") as mock_logger:
            result = _read_blacklist_files([nonexistent])

        assert result == set()
        mock_logger.warning.assert_called_once()

    def test_read_blacklist_files_empty_file(self, tmp_path):
        """Test handling empty blacklist file."""
        empty_file = tmp_path / "empty.txt"
        empty_file.write_text("")

        result = _read_blacklist_files([empty_file])
        assert result == set()

    def test_read_blacklist_files_whitespace_handling(self, tmp_path):
        """Test proper whitespace handling."""
        blacklist = tmp_path / "blacklist.txt"
        blacklist.write_text("  SRR123  \n\t SRR456 \t\n   \n")

        result = _read_blacklist_files([blacklist])
        assert result == {"SRR123", "SRR456"}

    def test_blacklist_reader_ignores_inline_reasons_and_comment_lines(self, tmp_path):
        from prepare.data.sra.download import _read_blacklist_files

        bl = tmp_path / "blacklist.txt"
        bl.write_text("# written by prepare blacklist\nSRR1  # amplicon\nSRR2\n\n")
        assert _read_blacklist_files([bl]) == {"SRR1", "SRR2"}


class TestPrepareTempFolder:
    """Test _prepare_temp_folder function."""

    def test_prepare_temp_folder_create_new(self, tmp_path):
        """Test creating new temp folder."""
        temp_base = tmp_path / "temp"

        result = _prepare_temp_folder(temp_base)

        assert result.exists()
        assert result.is_dir()

    def test_prepare_temp_folder_existing_folder(self, tmp_path):
        """Test with existing temp folder."""
        temp_base = tmp_path / "temp"
        temp_base.mkdir()

        result = _prepare_temp_folder(temp_base)

        assert result.exists()
        assert result.is_dir()

    def test_prepare_temp_folder_none_input(self):
        """Test with None temp folder."""

        with patch("tempfile.mkdtemp") as mock_mkdtemp:
            mock_mkdtemp.return_value = "/tmp/test_temp"

            result = _prepare_temp_folder(None)

            mock_mkdtemp.assert_called_once()
            assert str(result) == "/tmp/test_temp"


class TestCheckExistingDownload:
    """Test _check_existing_download function."""

    def test_check_existing_download_no_files(self, tmp_path):
        """Test when no files exist."""
        output_path = tmp_path / "downloads" / "SRR123"

        result = _check_existing_download(output_path, force=False)
        assert result is False

    def test_check_existing_download_with_files_no_force(self, tmp_path):
        """Test when files exist and force=False."""
        output_path = tmp_path / "downloads" / "SRR123"
        output_path.mkdir(parents=True)
        (output_path / "SRR123.fastq").write_text("test data")

        result = _check_existing_download(output_path, force=False)
        assert result is True

    def test_check_existing_download_with_files_force(self, tmp_path):
        """Test when files exist and force=True."""
        output_path = tmp_path / "downloads" / "SRR123"
        output_path.mkdir(parents=True)
        (output_path / "SRR123.fastq").write_text("test data")

        with patch("shutil.rmtree") as mock_rmtree:
            result = _check_existing_download(output_path, force=True)

        assert result is False
        mock_rmtree.assert_called_once_with(output_path)

    def test_check_existing_download_empty_directory(self, tmp_path):
        """Test when directory exists but is empty."""
        output_path = tmp_path / "downloads" / "SRR123"
        output_path.mkdir(parents=True)

        with patch("prepare.data.sra.accession.logger"):
            result = _check_existing_download(output_path, force=False)

        assert result is False
        # Directory should be removed
        assert not output_path.exists()

    def test_check_existing_download_rmdir_error(self, tmp_path):
        """Test handling rmdir error."""
        output_path = tmp_path / "downloads" / "SRR123"
        output_path.mkdir(parents=True)

        with patch.object(Path, "rmdir", side_effect=OSError("Permission denied")):
            with patch("prepare.data.sra.accession.logger") as mock_logger:
                result = _check_existing_download(output_path, force=False)

        assert result is False
        mock_logger.warning.assert_called()

    def test_check_existing_download_dangling_symlink(self, tmp_path):
        """A broken symlink left by an interrupted run must be cleared, never raise."""
        output_path = tmp_path / "downloads" / "SRR123"
        output_path.parent.mkdir(parents=True)
        output_path.symlink_to(tmp_path / "gone-target")

        result = _check_existing_download(output_path, force=False)

        assert result is False
        assert not output_path.is_symlink()

    def test_check_existing_download_symlinked_empty_directory_is_unlinked(self, tmp_path):
        """A symlink to an empty directory is unlinked, not rmdir'd (rmdir rejects symlinks)."""
        real_dir = tmp_path / "real_empty"
        real_dir.mkdir()
        output_path = tmp_path / "downloads" / "SRR123"
        output_path.parent.mkdir(parents=True)
        output_path.symlink_to(real_dir)

        result = _check_existing_download(output_path, force=False)

        assert result is False
        assert not output_path.is_symlink()

    def test_check_existing_download_redownload_truncated_forces_fresh(self, tmp_path):
        """force=True (standing in for a truncated accession) removes the directory."""
        output_path = tmp_path / "downloads" / "SRR123"
        output_path.mkdir(parents=True)
        (output_path / "SRR123.fastq").write_text("partial")

        result = _check_existing_download(output_path, force=True)

        assert result is False
        assert not output_path.exists()


class TestFastqFiles:
    """Test fastq_files, the single source of truth for what counts as a FASTQ file on disk."""

    def test_returns_sorted_nonempty_matches(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        (acc_dir / "SRR1_2.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (acc_dir / "SRR1_empty.fastq").write_text("")
        (acc_dir / "notes.txt").write_text("not fastq")

        result = fastq_files(acc_dir)

        assert [p.name for p in result] == ["SRR1_1.fastq", "SRR1_2.fastq"]

    def test_matches_fq_gz(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        with gzip.open(acc_dir / "SRR1.fq.gz", "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n")

        assert [p.name for p in fastq_files(acc_dir)] == ["SRR1.fq.gz"]

    def test_missing_directory_returns_empty(self, tmp_path):
        assert fastq_files(tmp_path / "missing") == []

    def test_follows_symlinked_directory(self, tmp_path):
        real_dir = tmp_path / "real"
        real_dir.mkdir()
        (real_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        link = tmp_path / "SRR1"
        link.symlink_to(real_dir)

        assert [p.name for p in fastq_files(link)] == ["SRR1_1.fastq"]


class TestFastqFilesIgnoresAppleDouble:
    def test_appledouble_not_listed(self, tmp_path):
        acc = tmp_path / "SRR1"
        acc.mkdir()
        (acc / "SRR1_1.fastq.gz").write_bytes(b"x" * 10)
        (acc / "._SRR1_1.fastq.gz").write_bytes(b"\x00\x05\x16\x07")
        assert [p.name for p in fastq_files(acc)] == ["SRR1_1.fastq.gz"]

    def test_folder_with_only_appledouble_has_no_fastq(self, tmp_path):
        acc = tmp_path / "SRR2"
        acc.mkdir()
        (acc / "._SRR2_1.fastq.gz").write_bytes(b"\x00\x05\x16\x07")
        assert fastq_files(acc) == []
        assert accession_has_fastq(acc) is False


class TestAccessionHasFastq:
    """Test accession_has_fastq, the single source of truth for "already downloaded"."""

    def test_false_for_zero_byte_file(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        (acc_dir / "SRR1.fastq").write_text("")
        assert accession_has_fastq(acc_dir) is False

    def test_true_for_fq_gz(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        with gzip.open(acc_dir / "SRR1.fq.gz", "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n")
        assert accession_has_fastq(acc_dir) is True

    def test_false_when_sidecar_marks_partial(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        (acc_dir / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (acc_dir / "SRR1.json").write_text(json.dumps({"state": "partial"}))
        assert accession_has_fastq(acc_dir) is False

    def test_false_when_sidecar_marks_failed_or_downloading(self, tmp_path):
        for state in ("failed", "downloading"):
            acc_dir = tmp_path / f"SRR_{state}"
            acc_dir.mkdir()
            (acc_dir / f"{acc_dir.name}.fastq").write_text("@r\nACGT\n+\nIIII\n")
            (acc_dir / f"{acc_dir.name}.json").write_text(json.dumps({"state": state}))
            assert accession_has_fastq(acc_dir) is False

    def test_true_when_sidecar_marks_complete(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        (acc_dir / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (acc_dir / "SRR1.json").write_text(json.dumps({"state": "complete"}))
        assert accession_has_fastq(acc_dir) is True

    def test_true_when_sidecar_is_unreadable(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        (acc_dir / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (acc_dir / "SRR1.json").write_text("not valid json{")
        assert accession_has_fastq(acc_dir) is True

    def test_false_for_missing_directory(self, tmp_path):
        assert accession_has_fastq(tmp_path / "missing") is False


class TestCountFastqReads:
    """count_fastq_reads: a chunked binary newline count, // 4 per FASTQ record."""

    def test_plain_and_gz(self, tmp_path):
        plain = tmp_path / "a.fastq"
        plain.write_text("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(3)))
        assert count_fastq_reads(plain) == 3

        gz = tmp_path / "b.fastq.gz"
        with gzip.open(gz, "wt") as handle:
            handle.write("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(5)))
        assert count_fastq_reads(gz) == 5

    def test_multiple_internal_blocks(self, tmp_path):
        """A file over 1 MiB must still be counted correctly across chunk boundaries."""
        record = "@r{0}\n" + "A" * 32 + "\n+\n" + "I" * 32 + "\n"
        path = tmp_path / "big.fastq"
        with open(path, "w") as handle:
            for i in range(20000):
                handle.write(record.format(i))
        assert path.stat().st_size > 1024 * 1024
        assert count_fastq_reads(path) == 20000

    def test_missing_trailing_newline_plain(self, tmp_path):
        """A file whose last record has no trailing newline must not be undercounted."""
        content = "".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(2))
        content = content[:-1]  # drop the final newline
        path = tmp_path / "a.fastq"
        path.write_text(content)
        assert count_fastq_reads(path) == 2

    def test_missing_trailing_newline_gz(self, tmp_path):
        content = "".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(2))
        content = content[:-1]
        path = tmp_path / "a.fastq.gz"
        with gzip.open(path, "wt") as handle:
            handle.write(content)
        assert count_fastq_reads(path) == 2

    def test_empty_file(self, tmp_path):
        path = tmp_path / "empty.fastq"
        path.write_text("")
        assert count_fastq_reads(path) == 0


class TestFastqDigest:
    """fastq_digest: one pass over the stored bytes yields the record count, md5 and size."""

    @staticmethod
    def _expected(path):
        from prepare.store.sidecar import md5_file

        return (count_fastq_reads(path), md5_file(path), path.stat().st_size)

    @staticmethod
    def _as_tuple(digest):
        return (digest.records, digest.md5, digest.size)

    def test_plain(self, tmp_path):
        path = tmp_path / "a.fastq"
        path.write_text("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(3)))
        assert self._as_tuple(fastq_digest(path)) == self._expected(path)
        assert fastq_digest(path).records == 3

    def test_gz_single_member(self, tmp_path):
        path = tmp_path / "a.fastq.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(5)))
        assert self._as_tuple(fastq_digest(path)) == self._expected(path)
        assert fastq_digest(path).records == 5

    def test_gz_multi_member(self, tmp_path):
        """Concatenated gzip members (as pigz or ``cat a.gz b.gz`` produce) are all counted."""
        path = tmp_path / "multi.fastq.gz"
        members = [gzip.compress("".join(f"@m{m}r{i}\nACGT\n+\nIIII\n" for i in range(4)).encode()) for m in range(3)]
        path.write_bytes(b"".join(members))
        assert self._as_tuple(fastq_digest(path)) == self._expected(path)
        assert fastq_digest(path).records == 12

    def test_large_gz_spans_many_blocks(self, tmp_path):
        """Random sequence keeps the compressed file over 1 MiB, so blocks and members split mid-record."""
        import random

        rng = random.Random(7)
        path = tmp_path / "big.fastq.gz"
        members = []
        for member in range(2):
            records = []
            for i in range(30000):
                seq = "".join(rng.choices("ACGT", k=80))
                records.append(f"@m{member}r{i}\n{seq}\n+\n{'I' * 80}\n")
            members.append(gzip.compress("".join(records).encode()))
        path.write_bytes(b"".join(members))
        assert path.stat().st_size > 1024 * 1024
        assert self._as_tuple(fastq_digest(path)) == self._expected(path)
        assert fastq_digest(path).records == 60000

    def test_empty_plain(self, tmp_path):
        path = tmp_path / "empty.fastq"
        path.write_text("")
        assert self._as_tuple(fastq_digest(path)) == self._expected(path)

    def test_empty_gz(self, tmp_path):
        path = tmp_path / "empty.fastq.gz"
        with gzip.open(path, "wt"):
            pass
        assert self._as_tuple(fastq_digest(path)) == self._expected(path)

    def test_missing_final_newline_plain_and_gz(self, tmp_path):
        content = "".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(2))[:-1]
        plain = tmp_path / "a.fastq"
        plain.write_text(content)
        gz = tmp_path / "a.fastq.gz"
        with gzip.open(gz, "wt") as handle:
            handle.write(content)
        for path in (plain, gz):
            assert self._as_tuple(fastq_digest(path)) == self._expected(path)
            assert fastq_digest(path).records == 2

    def test_truncated_gz_raises_like_count_fastq_reads(self, tmp_path):
        path = tmp_path / "cut.fastq.gz"
        data = gzip.compress(b"@r\nACGT\n+\nIIII\n" * 1000)
        path.write_bytes(data[: len(data) // 2])
        with pytest.raises(EOFError):
            count_fastq_reads(path)
        with pytest.raises(EOFError):
            fastq_digest(path)

    def test_non_gzip_file_raises_like_gzip_open(self, tmp_path):
        """Same class (an OSError, as store adoption expects) and wording as count_fastq_reads."""
        path = tmp_path / "plain_named.fastq.gz"
        path.write_bytes(b"@r\nACGT\n+\nIIII\n")
        with pytest.raises(gzip.BadGzipFile) as expected:
            count_fastq_reads(path)
        with pytest.raises(gzip.BadGzipFile) as got:
            fastq_digest(path)
        assert str(got.value) == str(expected.value)

    def test_corrupt_deflate_stream_raises_bad_gzip_file(self, tmp_path):
        path = tmp_path / "corrupt.fastq.gz"
        data = bytearray(gzip.compress(b"".join(b"@r%d\nACGT\n+\nIIII\n" % i for i in range(200))))
        data[20] ^= 0xFF
        path.write_bytes(bytes(data))
        with pytest.raises(gzip.BadGzipFile):
            fastq_digest(path)

    def test_content_md5_is_md5_of_decompressed_bytes(self, tmp_path):
        import hashlib

        text = b"@r\nACGT\n+\nIIII\n" * 3
        gz = tmp_path / "a.fastq.gz"
        gz.write_bytes(gzip.compress(text))
        plain = tmp_path / "b.fastq"
        plain.write_bytes(text)
        expected = hashlib.md5(text).hexdigest()
        assert fastq_digest(gz, content_md5=True).content_md5 == expected
        assert fastq_digest(plain, content_md5=True).content_md5 == expected
        assert fastq_digest(gz).content_md5 is None


class TestIterFastqRecords:
    """iter_fastq_records: the shared raw four-line reader."""

    def test_reads_every_record(self, tmp_path):
        path = tmp_path / "a.fastq"
        path.write_text("@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n")
        assert list(iter_fastq_records(path)) == [("ACGT", "IIII"), ("TTTT", "JJJJ")]

    def test_zero_length_record_in_the_middle_is_a_record(self, tmp_path):
        """An empty sequence with an empty quality line is legal FASTQ, not truncation.

        fasterq-dump writes such a record for the empty mate of a half-empty spot, and a
        trimmer can reduce a read to length zero.
        """
        path = tmp_path / "zero.fastq"
        path.write_text("@r1\nACGT\n+\nIIII\n@r2\n\n+\n\n@r3\nACGT\n+\nIIII\n")
        assert list(iter_fastq_records(path)) == [("ACGT", "IIII"), ("", ""), ("ACGT", "IIII")]

    def test_zero_length_record_in_a_gz_file(self, tmp_path):
        path = tmp_path / "zero.fastq.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("@r1\n\n+\n\n@r2\nACGT\n+\nIIII\n")
        assert list(iter_fastq_records(path)) == [("", ""), ("ACGT", "IIII")]

    def test_file_cut_after_the_sequence_line_raises(self, tmp_path):
        path = tmp_path / "cut.fastq"
        path.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGT\n")
        with pytest.raises(ValueError, match="Truncated FASTQ record"):
            list(iter_fastq_records(path))

    def test_file_cut_after_the_header_raises(self, tmp_path):
        path = tmp_path / "cut.fastq"
        path.write_text("@r1\n")
        with pytest.raises(ValueError, match="Truncated FASTQ record"):
            list(iter_fastq_records(path))

    def test_a_third_line_that_is_not_a_separator_raises(self, tmp_path):
        path = tmp_path / "bad.fastq"
        path.write_text("@r1\nACGT\nIIII\nACGT\n")
        with pytest.raises(ValueError, match="Malformed FASTQ record"):
            list(iter_fastq_records(path))


class TestPrimaryAndOrphanFastq:
    """primary_fastq / orphan_fastq: which file in an accession folder holds the spots."""

    def _write(self, acc_dir, *names):
        acc_dir.mkdir(parents=True, exist_ok=True)
        for name in names:
            (acc_dir / name).write_text("@r\nACGT\n+\nIIII\n")

    def test_primary_prefers_mate_one_over_bare_file(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        self._write(acc_dir, "SRR1.fastq", "SRR1_1.fastq", "SRR1_2.fastq")

        assert primary_fastq(acc_dir) == acc_dir / "SRR1_1.fastq"

    def test_primary_accepts_r1_naming(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        self._write(acc_dir, "SRR1_R1.fastq.gz", "SRR1_R2.fastq.gz")

        assert primary_fastq(acc_dir) == acc_dir / "SRR1_R1.fastq.gz"

    def test_primary_is_bare_file_for_single_end(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        self._write(acc_dir, "SRR1.fastq")

        assert primary_fastq(acc_dir) == acc_dir / "SRR1.fastq"

    def test_primary_none_for_empty_folder(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()

        assert primary_fastq(acc_dir) is None
        assert orphan_fastq(acc_dir) is None

    def test_orphan_only_when_mates_present(self, tmp_path):
        paired = tmp_path / "SRR1"
        self._write(paired, "SRR1.fastq", "SRR1_1.fastq", "SRR1_2.fastq")
        single = tmp_path / "SRR2"
        self._write(single, "SRR2.fastq")

        assert orphan_fastq(paired) == paired / "SRR1.fastq"
        assert orphan_fastq(single) is None


class TestVerifyDownload:
    """verify_download: reads_r1 against NCBI's expected spot count."""

    def _acc_dir(self, tmp_path, name="SRR1"):
        acc_dir = tmp_path / name
        acc_dir.mkdir()
        (acc_dir / f"{name}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        return acc_dir

    def test_truncated(self, tmp_path, monkeypatch):
        acc_dir = self._acc_dir(tmp_path)
        monkeypatch.setattr("prepare.data.sra.fastq.count_fastq_reads", lambda path: 300000)

        result = verify_download("SRR1", acc_dir, expected_spots=48000000)

        assert result["verdict"] == "truncated"
        assert result["reads_r1"] == 300000
        assert result["expected_spots"] == 48000000
        assert result["ratio"] == round(300000 / 48000000, 4)

    def test_complete(self, tmp_path, monkeypatch):
        acc_dir = self._acc_dir(tmp_path)
        monkeypatch.setattr("prepare.data.sra.fastq.count_fastq_reads", lambda path: 999)

        result = verify_download("SRR1", acc_dir, expected_spots=1000)

        assert result["verdict"] == "complete"
        assert result["ratio"] == 0.999

    def test_complete_at_exact_threshold(self, tmp_path, monkeypatch):
        acc_dir = self._acc_dir(tmp_path)
        monkeypatch.setattr("prepare.data.sra.fastq.count_fastq_reads", lambda path: 990)

        result = verify_download("SRR1", acc_dir, expected_spots=1000)

        assert result["verdict"] == "complete"

    def test_truncated_just_below_threshold(self, tmp_path, monkeypatch):
        acc_dir = self._acc_dir(tmp_path)
        monkeypatch.setattr("prepare.data.sra.fastq.count_fastq_reads", lambda path: 989)

        result = verify_download("SRR1", acc_dir, expected_spots=1000)

        assert result["verdict"] == "truncated"

    def test_unverified_when_spots_unknown(self, tmp_path):
        acc_dir = self._acc_dir(tmp_path)

        result = verify_download("SRR1", acc_dir, expected_spots=None)

        assert result["verdict"] == "unverified"
        assert result["ratio"] is None
        assert result["expected_spots"] is None

    def test_counts_mate_one_plus_orphan_file(self, tmp_path):
        """--split-3 writes unpaired spots to a bare <acc>.fastq; those spots count too."""
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        record = "@r\nACGT\n+\nIIII\n"
        (acc_dir / "SRR1_1.fastq").write_text(record * 1000)
        (acc_dir / "SRR1_2.fastq").write_text(record * 1000)
        (acc_dir / "SRR1.fastq").write_text(record * 7)

        result = verify_download("SRR1", acc_dir, expected_spots=1007)

        assert result["reads_r1"] == 1007
        assert result["ratio"] == 1.0
        assert result["verdict"] == "complete"

    def test_single_end_counts_the_bare_file(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        (acc_dir / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n" * 500)

        result = verify_download("SRR1", acc_dir, expected_spots=500)

        assert result["reads_r1"] == 500
        assert result["verdict"] == "complete"

    def test_no_expected_bytes_parameter(self):
        """NCBI's recorded size is the .sra archive size, so it never constrained the FASTQ bytes."""
        assert "expected_bytes" not in inspect.signature(verify_download).parameters

    def test_bytes_total_sums_fastq_files(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        (acc_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (acc_dir / "SRR1_2.fastq").write_text("@r\nACGT\n+\nIIII\n")
        expected_bytes = sum(p.stat().st_size for p in acc_dir.glob("*.fastq"))

        result = verify_download("SRR1", acc_dir, expected_spots=None)

        assert result["bytes_total"] == expected_bytes


class TestVerifyDownloadPrecomputedCounts:
    """verify_download takes read counts a caller already has and then reads no FASTQ file."""

    def test_precomputed_counts_open_no_file(self, tmp_path, monkeypatch):
        import builtins

        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        record = "@r\nACGT\n+\nIIII\n"
        (acc_dir / "SRR1_1.fastq").write_text(record * 1000)
        (acc_dir / "SRR1_2.fastq").write_text(record * 1000)
        (acc_dir / "SRR1.fastq").write_text(record * 7)

        def no_open(*args, **kwargs):
            raise AssertionError(f"verify_download opened {args[0]}")

        monkeypatch.setattr(builtins, "open", no_open)
        monkeypatch.setattr(gzip, "open", no_open)
        result = verify_download("SRR1", acc_dir, expected_spots=1007, reads_r1=1000, reads_orphan=7)
        monkeypatch.undo()

        assert result["reads_r1"] == 1007
        assert result["verdict"] == "complete"

    def test_orphan_is_counted_when_only_mate_one_is_given(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        record = "@r\nACGT\n+\nIIII\n"
        (acc_dir / "SRR1_1.fastq").write_text(record * 1000)
        (acc_dir / "SRR1_2.fastq").write_text(record * 1000)
        (acc_dir / "SRR1.fastq").write_text(record * 7)

        result = verify_download("SRR1", acc_dir, expected_spots=None, reads_r1=1000)

        assert result["reads_r1"] == 1007


class TestCachedSraArchive:
    """_cached_sra_archive: which prefetch output fasterq-dump is pointed at."""

    def test_prefers_sra_over_sralite(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        (acc_dir / "SRR1.sralite").write_bytes(b"")
        (acc_dir / "SRR1.sra").write_bytes(b"")

        assert _cached_sra_archive(acc_dir, "SRR1") == acc_dir / "SRR1.sra"

    def test_finds_a_sralite_only_download(self, tmp_path):
        acc_dir = tmp_path / "SRR1"
        acc_dir.mkdir()
        (acc_dir / "SRR1.sralite").write_bytes(b"")

        assert _cached_sra_archive(acc_dir, "SRR1") == acc_dir / "SRR1.sralite"

    def test_none_when_nothing_was_written(self, tmp_path):
        """An empty cache folder yields None, which download_accession reports as a network failure."""
        assert _cached_sra_archive(tmp_path / "SRR1", "SRR1") is None


class TestParseVerdictMessage:
    def test_complete(self):
        result = parse_verdict_message("Downloaded 2 files, complete (300000 of 300000 spots)")
        assert result == {"verdict": "complete", "reads_r1": 300000, "expected_spots": 300000, "ratio": 1.0}

    def test_truncated(self):
        result = parse_verdict_message("Downloaded 1 files, truncated (300000 of 48000000 spots)")
        assert result["verdict"] == "truncated"
        assert result["reads_r1"] == 300000
        assert result["expected_spots"] == 48000000

    def test_truncated_with_retry_prefix(self):
        result = parse_verdict_message("Retry 1: Downloaded 1 files, truncated (300000 of 48000000 spots)")
        assert result["verdict"] == "truncated"

    def test_unverified(self):
        assert parse_verdict_message("Downloaded 2 files, unverified") == {"verdict": "unverified"}

    def test_unrelated_message_returns_none(self):
        assert parse_verdict_message("Download failed: timeout") is None
        assert parse_verdict_message("already exists") is None

    def test_verdict_words_outside_the_message_format_are_ignored(self):
        """Only the ", <verdict>" text _handle_download_output writes counts as a verdict."""
        assert parse_verdict_message("Download failed: unverified checksum for SRR1") is None
        assert parse_verdict_message("not-found: could not complete (0 of 5 spots) for SRR1") is None


class TestHandleDownloadOutput:
    """Test _handle_download_output function."""

    def test_handle_download_output_success(self, tmp_path):
        """Test successful file moving."""
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output"
        temp_path.mkdir()

        # Create test FASTQ files
        (temp_path / "SRR123.fastq").write_text("@seq1\nACGT\n+\nIIII\n")
        (temp_path / "SRR123_2.fastq").write_text("@seq2\nTGCA\n+\nIIII\n")

        success, message = _handle_download_output(temp_path, output_path)

        assert success is True
        assert "2 files" in message
        assert (output_path / "SRR123.fastq").exists()
        assert (output_path / "SRR123_2.fastq").exists()
        assert not temp_path.exists()

    def test_handle_download_output_no_fastq_files(self, tmp_path):
        """No FASTQ files found: message is classified and the temp folder is kept for inspection."""
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output"
        temp_path.mkdir()

        # Create non-FASTQ file
        (temp_path / "other.txt").write_text("not a fastq")

        with patch("shutil.rmtree") as mock_rmtree:
            with patch("prepare.data.sra.accession.logger") as mock_logger:
                success, message = _handle_download_output(temp_path, output_path)

        assert success is False
        assert message == "unknown: No FASTQ files created"
        mock_logger.error.assert_called()
        mock_rmtree.assert_not_called()
        assert temp_path.exists()

    def test_handle_download_output_rmdir_error(self, tmp_path):
        """Test handling rmtree error during cleanup."""
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output"
        temp_path.mkdir()

        (temp_path / "SRR123.fastq").write_text("@seq1\nACGT\n+\nIIII\n")

        with patch("shutil.rmtree", side_effect=OSError("Permission denied")):
            with patch("prepare.data.sra.accession.logger") as mock_logger:
                success, message = _handle_download_output(temp_path, output_path)

        assert success is True
        mock_logger.warning.assert_called()

    def test_handle_download_output_complete_verdict(self, tmp_path, monkeypatch):
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output" / "SRR123"
        temp_path.mkdir()
        (temp_path / "SRR123_1.fastq").write_text("@seq1\nACGT\n+\nIIII\n")
        (temp_path / "SRR123_2.fastq").write_text("@seq2\nTGCA\n+\nIIII\n")
        monkeypatch.setattr("prepare.data.sra.fastq.count_fastq_reads", lambda path: 300000)

        success, message = _handle_download_output(temp_path, output_path, expected_spots=300000)

        assert success is True
        assert message == "Downloaded 2 files, complete (300000 of 300000 spots)"

    def test_handle_download_output_truncated_verdict_logs_warning(self, tmp_path, monkeypatch):
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output" / "SRR123"
        temp_path.mkdir()
        (temp_path / "SRR123_1.fastq").write_text("@seq1\nACGT\n+\nIIII\n")
        monkeypatch.setattr("prepare.data.sra.fastq.count_fastq_reads", lambda path: 300000)

        with patch("prepare.data.sra.accession.logger") as mock_logger:
            success, message = _handle_download_output(temp_path, output_path, expected_spots=48000000)

        assert success is False
        assert message == (
            "Downloaded 1 files, truncated (300000 of 48000000 spots); rerun with --accept-partial to keep it"
        )
        assert parse_verdict_message(message)["verdict"] == "truncated"
        mock_logger.warning.assert_called()

    def test_handle_download_output_truncated_is_accepted_with_accept_partial(self, tmp_path, monkeypatch):
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output" / "SRR123"
        temp_path.mkdir()
        (temp_path / "SRR123_1.fastq").write_text("@seq1\nACGT\n+\nIIII\n")
        monkeypatch.setattr("prepare.data.sra.fastq.count_fastq_reads", lambda path: 300000)

        success, message = _handle_download_output(temp_path, output_path, expected_spots=48000000, accept_partial=True)

        assert success is True
        assert message == "Downloaded 1 files, truncated (300000 of 48000000 spots)"

    def test_handle_download_output_unverified_without_expected_spots(self, tmp_path):
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output" / "SRR123"
        temp_path.mkdir()
        (temp_path / "SRR123_1.fastq").write_text("@seq1\nACGT\n+\nIIII\n")

        success, message = _handle_download_output(temp_path, output_path)

        assert success is True
        assert message == "Downloaded 1 files, unverified"

    def test_handle_download_output_reports_compression_failure_in_message(self, tmp_path):
        """A per-file compression failure must not fail the download, but must be visible."""
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output" / "SRR123"
        temp_path.mkdir()
        (temp_path / "SRR123_1.fastq").write_text("@seq1\nACGT\n+\nIIII\n")

        with patch("prepare.data.sra.fastq.compress_fastq", side_effect=OSError("disk full")):
            with patch("prepare.data.sra.accession.logger") as mock_logger:
                success, message = _handle_download_output(temp_path, output_path, compress=True)

        assert success is True
        assert "compression failed for SRR123_1.fastq" in message
        # The verdict text itself (parseable by parse_verdict_message) must still be present.
        assert message.startswith("Downloaded 1 files, unverified")
        mock_logger.warning.assert_called()

    def test_handle_download_output_short_gz_after_compression_is_a_failure(self, tmp_path):
        """A compressed file holding fewer records than the plain one fails the download."""
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output" / "SRR123"
        temp_path.mkdir()
        records = "".join(f"@seq{i}\nACGT\n+\nIIII\n" for i in range(4))
        (temp_path / "SRR123_1.fastq").write_text(records)

        def short_compress(path, threads, stop=None):
            target = path.with_suffix(path.suffix + ".gz")
            with gzip.open(target, "wt") as handle:
                handle.write("".join(f"@seq{i}\nACGT\n+\nIIII\n" for i in range(2)))
            path.unlink()
            return target

        with patch("prepare.data.sra.fastq.compress_fastq", side_effect=short_compress):
            success, message = _handle_download_output(temp_path, output_path, compress=True)

        assert success is False
        assert message == "compression corrupted SRR123_1.fastq.gz: 2 of 4 reads"

    def test_handle_download_output_unreadable_gz_is_a_failure(self, tmp_path):
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output" / "SRR123"
        temp_path.mkdir()
        (temp_path / "SRR123_1.fastq").write_text("@seq1\nACGT\n+\nIIII\n")

        def bad_compress(path, threads, stop=None):
            target = path.with_suffix(path.suffix + ".gz")
            target.write_bytes(b"\x1f\x8b\x08\x00garbage")
            path.unlink()
            return target

        with patch("prepare.data.sra.fastq.compress_fastq", side_effect=bad_compress):
            success, message = _handle_download_output(temp_path, output_path, compress=True)

        assert success is False
        assert message.startswith("compression corrupted SRR123_1.fastq.gz: unreadable")

    def test_handle_download_output_counts_match_after_compression(self, tmp_path):
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output" / "SRR123"
        temp_path.mkdir()
        (temp_path / "SRR123_1.fastq").write_text("".join(f"@seq{i}\nACGT\n+\nIIII\n" for i in range(4)))

        with patch("prepare.data.sra.accession.shutil.which", return_value=None):
            success, message = _handle_download_output(temp_path, output_path, compress=True)

        assert success is True
        assert message == "Downloaded 1 files, unverified"
        assert (output_path / "SRR123_1.fastq.gz").exists()

    def test_compress_files_counts_each_plain_file_once(self, tmp_path):
        """The verdict and the compression read-back share one count per plain file."""
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output" / "SRR123"
        temp_path.mkdir()
        for name in ("SRR123_1.fastq", "SRR123_2.fastq", "SRR123.fastq"):
            (temp_path / name).write_text("".join(f"@seq{i}\nACGT\n+\nIIII\n" for i in range(4)))
        counted = []
        real_count = fastq_mod.count_fastq_reads

        def counting(path):
            counted.append(Path(path).name)
            return real_count(path)

        with (
            patch("prepare.data.sra.accession.shutil.which", return_value=None),
            patch("prepare.data.sra.fastq.count_fastq_reads", side_effect=counting),
        ):
            success, message = _handle_download_output(temp_path, output_path, expected_spots=8, compress=True)

        assert success is True, message
        assert message == "Downloaded 3 files, complete (8 of 8 spots)"
        plain = [name for name in counted if not name.endswith(".gz")]
        assert sorted(plain) == ["SRR123.fastq", "SRR123_1.fastq", "SRR123_2.fastq"], counted

    def test_no_compress_counts_only_the_primary_and_orphan_files(self, tmp_path):
        """Without compression the mate-2 file is not read: only the verdict's files are counted."""
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output" / "SRR123"
        temp_path.mkdir()
        for name in ("SRR123_1.fastq", "SRR123_2.fastq"):
            (temp_path / name).write_text("".join(f"@seq{i}\nACGT\n+\nIIII\n" for i in range(4)))
        counted = []
        real_count = fastq_mod.count_fastq_reads

        def counting(path):
            counted.append(Path(path).name)
            return real_count(path)

        with patch("prepare.data.sra.fastq.count_fastq_reads", side_effect=counting):
            success, message = _handle_download_output(temp_path, output_path, expected_spots=4, compress=False)

        assert success is True and message == "Downloaded 2 files, complete (4 of 4 spots)"
        assert counted == ["SRR123_1.fastq"]

    def test_a_stale_file_in_the_output_folder_does_not_break_the_count(self, tmp_path):
        """A primary file left in the folder by an earlier attempt is counted, not a KeyError."""
        temp_path = tmp_path / "temp"
        output_path = tmp_path / "output" / "SRR123"
        temp_path.mkdir()
        output_path.mkdir(parents=True)
        (output_path / "SRR123_1.fastq").write_text("@old\nACGT\n+\nIIII\n")
        (temp_path / "SRR123.fastq").write_text("".join(f"@seq{i}\nACGT\n+\nIIII\n" for i in range(4)))

        success, message = _handle_download_output(temp_path, output_path, expected_spots=None, compress=False)

        assert success is True and message == "Downloaded 1 files, unverified"


class TestDownloadAccession:
    """Test download_accession function."""

    def test_download_accession_existing_skip(self, tmp_path):
        """Test skipping existing download."""
        output_folder = tmp_path / "downloads"
        output_path = output_folder / "SRR123"
        output_path.mkdir(parents=True)
        (output_path / "SRR123.fastq").write_text("existing")

        success, message = download_accession("SRR123", output_folder, force=False)

        assert success is True
        assert "already exists" in message

    def test_download_accession_force_existing(self, tmp_path):
        """Test forcing download of existing accession."""
        output_folder = tmp_path / "downloads"
        output_path = output_folder / "SRR123"
        output_path.mkdir(parents=True)
        (output_path / "SRR123.fastq").write_text("existing")

        with patch("prepare.data.sra.cleanup._prepare_temp_folder") as mock_prep:
            with patch("prepare.utils.security.SecureSubprocess.run_secure") as mock_run:
                with patch("prepare.data.sra.accession._handle_download_output") as mock_handle:
                    with patch("prepare.data.sra.accession.shutil.which", return_value=None):
                        mock_prep.return_value = tmp_path / "temp"
                        mock_run.return_value = Mock(returncode=0, stdout="success", stderr="")
                        mock_handle.return_value = (True, "Downloaded 1 files")

                        success, message = download_accession("SRR123", output_folder, force=True)

        assert "Downloaded 1 files" in message

    def test_download_accession_security_error(self, tmp_path):
        """Test handling security error."""
        output_folder = tmp_path / "downloads"

        with patch("prepare.data.sra.cleanup._prepare_temp_folder") as mock_prep:
            with patch("prepare.utils.security.SecureSubprocess.run_secure") as mock_run:
                mock_prep.return_value = tmp_path / "temp"
                mock_run.side_effect = SecurityError("Command blocked")

                success, message = download_accession("SRR123", output_folder)

        assert success is False
        assert "Security error" in message

    def test_download_accession_command_failure(self, tmp_path):
        """Test handling command failure."""
        import subprocess

        output_folder = tmp_path / "downloads"

        with patch("prepare.data.sra.cleanup._prepare_temp_folder") as mock_prep:
            with patch("prepare.utils.security.SecureSubprocess.run_secure") as mock_run:
                mock_prep.return_value = tmp_path / "temp"
                mock_run.side_effect = subprocess.CalledProcessError(1, "fasterq-dump", stderr="Error occurred")

                success, message = download_accession("SRR123", output_folder)

        assert success is False
        assert "Download failed" in message

    def test_download_accession_success(self, tmp_path):
        """Test successful download."""
        output_folder = tmp_path / "downloads"
        temp_path = tmp_path / "temp"

        with patch("prepare.data.sra.accession.shutil.which", return_value=None):
            with patch("prepare.data.sra.cleanup._prepare_temp_folder") as mock_prep:
                with patch("prepare.utils.security.SecureSubprocess.run_secure") as mock_run:
                    with patch("prepare.data.sra.accession._handle_download_output") as mock_handle:
                        mock_prep.return_value = temp_path
                        mock_run.return_value = Mock(returncode=0, stdout="success", stderr="")
                        mock_handle.return_value = (True, "Downloaded 2 files")

                        success, message = download_accession("SRR123", output_folder, num_threads=8)

        assert success is True
        assert "Downloaded 2 files" in message
        # Verify command was called with correct arguments
        mock_run.assert_called_once()
        executable = mock_run.call_args[0][0]
        args = mock_run.call_args[0][1]
        assert executable == "fasterq-dump"
        assert "SRR123" in args
        assert "--threads" in args
        assert "8" in args

    def test_download_accession_command_passes_validation(self, tmp_path, monkeypatch):
        """The command download_accession builds must survive SecureSubprocess validation unmocked.

        prefetch is reported missing from PATH, so this exercises the direct
        fasterq-dump fallback with --split-3/--skip-technical appended.
        """
        monkeypatch.chdir(tmp_path)
        proc = Mock(returncode=0)
        proc.communicate.return_value = ("", "")
        with patch("prepare.utils.security.subprocess.Popen", return_value=proc) as mock_popen:
            with patch("prepare.data.sra.accession.shutil.which", return_value=None):
                with patch(
                    "prepare.data.sra.accession._handle_download_output", return_value=(True, "Downloaded 2 files")
                ):
                    success, message = download_accession(
                        "SRR2517620", tmp_path / "fastq", num_threads=4, temp_folder=tmp_path / "tmp"
                    )
        assert success is True, message
        cmd = mock_popen.call_args[0][0]
        assert cmd[0] == "fasterq-dump"
        assert cmd[1:4] == ["--threads", "4", "--progress"]
        assert cmd[4] == "SRR2517620"
        assert "--split-3" in cmd
        assert "--skip-technical" in cmd

    def test_download_accession_registers_output_and_temp_roots(self, tmp_path, monkeypatch):
        from prepare.utils.security import SecureSubprocess

        SecureSubprocess._extra_roots.clear()
        monkeypatch.chdir(tmp_path)
        proc = Mock(returncode=0)
        proc.communicate.return_value = ("", "")
        with patch("prepare.utils.security.subprocess.Popen", return_value=proc):
            with patch("prepare.data.sra.accession.shutil.which", return_value=None):
                with patch(
                    "prepare.data.sra.accession._handle_download_output", return_value=(True, "Downloaded 2 files")
                ):
                    download_accession("SRR2517620", tmp_path / "fastq", temp_folder=tmp_path / "scratch")
        assert (tmp_path / "fastq").resolve() in SecureSubprocess._extra_roots
        assert (tmp_path / "scratch").resolve() in SecureSubprocess._extra_roots
        SecureSubprocess._extra_roots.clear()

    def test_download_accession_passes_expected_spots_to_handle_output(self, tmp_path):
        """expected_spots reaches _handle_download_output so the verdict can be computed."""
        output_folder = tmp_path / "downloads"

        with patch("prepare.data.sra.cleanup._prepare_temp_folder") as mock_prep:
            with patch("prepare.utils.security.SecureSubprocess.run_secure") as mock_run:
                with patch("prepare.data.sra.accession._handle_download_output") as mock_handle:
                    with patch("prepare.data.sra.accession.shutil.which", return_value=None):
                        mock_prep.return_value = tmp_path / "temp"
                        mock_run.return_value = Mock(returncode=0, stdout="", stderr="")
                        mock_handle.return_value = (True, "Downloaded 1 files, complete (10 of 10 spots)")

                        download_accession("SRR123", output_folder, expected_spots=10)

        mock_handle.assert_called_once()
        assert mock_handle.call_args.kwargs.get("expected_spots") == 10 or 10 in mock_handle.call_args.args

    def test_download_accession_redownload_truncated_forces_fresh_download(self, tmp_path):
        """A prior partial download on disk must not short-circuit as 'already exists'."""
        output_folder = tmp_path / "downloads"
        output_path = output_folder / "SRR123"
        output_path.mkdir(parents=True)
        (output_path / "SRR123.fastq").write_text("partial")

        with patch("prepare.data.sra.accession.shutil.which", return_value=None):
            with patch("prepare.data.sra.cleanup._prepare_temp_folder") as mock_prep:
                with patch("prepare.utils.security.SecureSubprocess.run_secure") as mock_run:
                    with patch("prepare.data.sra.accession._handle_download_output") as mock_handle:
                        mock_prep.return_value = tmp_path / "temp"
                        mock_run.return_value = Mock(returncode=0, stdout="", stderr="")
                        mock_handle.return_value = (True, "Downloaded 1 files, complete (10 of 10 spots)")

                        success, message = download_accession(
                            "SRR123", output_folder, force=False, redownload_truncated=True
                        )

        assert "already exists" not in message
        mock_run.assert_called_once()

    def test_download_accession_failure_removes_temp_folder(self, tmp_path):
        """A failed download removes <acc>_temp; a later attempt starts clean anyway."""
        import subprocess

        output_folder = tmp_path / "downloads"

        with patch("prepare.data.sra.cleanup._prepare_temp_folder") as mock_prep:
            with patch("prepare.utils.security.SecureSubprocess.run_secure") as mock_run:
                mock_prep.return_value = tmp_path / "temp"
                mock_run.side_effect = subprocess.CalledProcessError(1, "fasterq-dump", stderr="Connection timed out")

                success, message = download_accession("SRR123", output_folder)

        assert success is False
        temp_path = output_folder / "SRR123_temp"
        assert not temp_path.exists()

    def test_download_accession_command_failure_classifies_network_error(self, tmp_path):
        import subprocess

        output_folder = tmp_path / "downloads"

        with patch("prepare.data.sra.cleanup._prepare_temp_folder") as mock_prep:
            with patch("prepare.utils.security.SecureSubprocess.run_secure") as mock_run:
                mock_prep.return_value = tmp_path / "temp"
                mock_run.side_effect = subprocess.CalledProcessError(1, "fasterq-dump", stderr="Connection timed out")

                success, message = download_accession("SRR123", output_folder)

        assert success is False
        assert message.startswith("network:")

    def test_download_accession_security_error_removes_temp_folder(self, tmp_path):
        output_folder = tmp_path / "downloads"

        with patch("prepare.data.sra.cleanup._prepare_temp_folder") as mock_prep:
            with patch("prepare.utils.security.SecureSubprocess.run_secure") as mock_run:
                mock_prep.return_value = tmp_path / "temp"
                mock_run.side_effect = SecurityError("Command blocked")

                success, message = download_accession("SRR123", output_folder)

        assert success is False
        temp_path = output_folder / "SRR123_temp"
        assert not temp_path.exists()
        assert message.startswith("unknown:")

    def test_download_accession_prefetch_then_fasterq_dump_sequence(self, tmp_path):
        """With prefetch on PATH: prefetch runs first, then fasterq-dump --split-3 on the .sra file."""
        output_folder = tmp_path / "downloads"
        state = {"reads": 4}

        with patch("prepare.data.sra.accession.shutil.which", side_effect=lambda tool: f"/usr/bin/{tool}"):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
                success, message = download_accession("SRR123", output_folder, compress=False)

        assert success is True, message
        executables = [c[0] for c in state["calls"]]
        assert executables == ["prefetch", "fasterq-dump"]

        prefetch_args = state["calls"][0][1]
        assert prefetch_args[:2] == ["-O", str(output_folder / ".sra-cache")]
        assert "--max-size" in prefetch_args and "100G" in prefetch_args
        assert "--progress" in prefetch_args
        assert prefetch_args[-1] == "SRR123"

        dump_args = state["calls"][1][1]
        assert dump_args[0] == "--split-3"
        assert dump_args[1] == "--skip-technical"
        assert dump_args[-1] == str(output_folder / ".sra-cache" / "SRR123" / "SRR123.sra")

        assert accession_has_fastq(output_folder / "SRR123")
        # The .sra archive is removed after a successful, verified download.
        assert not (output_folder / ".sra-cache" / "SRR123").exists()

    def test_download_accession_keep_sra_preserves_cache(self, tmp_path):
        """keep_sra=True leaves the downloaded .sra archive in place after success."""
        output_folder = tmp_path / "downloads"
        state = {"reads": 4}

        with patch("prepare.data.sra.accession.shutil.which", side_effect=lambda tool: f"/usr/bin/{tool}"):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
                success, message = download_accession("SRR123", output_folder, compress=False, keep_sra=True)

        assert success is True, message
        assert (output_folder / ".sra-cache" / "SRR123" / "SRR123.sra").exists()

    def test_download_accession_force_clears_cached_archive_before_prefetch(self, tmp_path):
        """A forced download refetches the .sra instead of reusing a possibly truncated one."""
        output_folder = tmp_path / "downloads"
        cache_dir = output_folder / ".sra-cache" / "SRR123"
        cache_dir.mkdir(parents=True)
        (cache_dir / "SRR123.sra").write_bytes(b"stale archive")
        acc_dir = output_folder / "SRR123"
        acc_dir.mkdir(parents=True)
        (acc_dir / "SRR123_1.fastq").write_text("@r\nACGT\n+\nIIII\n")

        fake = _fake_tools({"reads": 4})
        cache_present_at_prefetch = []

        def run(executable, args, **kwargs):
            if executable == "prefetch":
                cache_present_at_prefetch.append(cache_dir.exists())
            return fake(executable, args, **kwargs)

        with patch("prepare.data.sra.accession.shutil.which", side_effect=lambda tool: f"/usr/bin/{tool}"):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=run):
                success, message = download_accession("SRR123", output_folder, compress=False, force=True)

        assert success is True, message
        assert cache_present_at_prefetch == [False]

    def test_download_accession_truncated_removes_cached_archive(self, tmp_path):
        """A truncated verdict leaves no .sra behind, so the next attempt fetches it again."""
        output_folder = tmp_path / "downloads"
        state = {"reads": 4}

        with patch("prepare.data.sra.accession.shutil.which", side_effect=lambda tool: f"/usr/bin/{tool}"):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
                success, message = download_accession("SRR123", output_folder, compress=False, expected_spots=1000)

        assert success is False, message
        assert "truncated" in message
        assert not (output_folder / ".sra-cache" / "SRR123").exists()

    def test_cached_archive_ignores_prefetch_leftovers(self, tmp_path):
        """A .sra.lock, .sra.prf or .sra.tmp beside a .sralite is never taken for the archive."""
        from prepare.data.sra.accession import _cached_sra_archive

        acc_dir = tmp_path / "cache" / "SRR123"
        acc_dir.mkdir(parents=True)
        for leftover in ("SRR123.sra.lock", "SRR123.sra.prf", "SRR123.sra.tmp"):
            (acc_dir / leftover).write_bytes(b"")
        (acc_dir / "SRR123.sralite").write_bytes(b"")

        assert _cached_sra_archive(acc_dir, "SRR123") == acc_dir / "SRR123.sralite"
        (acc_dir / "SRR123.sralite").unlink()
        assert _cached_sra_archive(acc_dir, "SRR123") is None

    def test_prefetch_without_archive_is_a_network_failure(self, tmp_path):
        """prefetch exiting 0 without writing an archive fails before fasterq-dump is run."""
        output_folder = tmp_path / "downloads"
        state = {"reads": 4}
        real = _fake_tools(state)

        def tools(executable, args, **kwargs):
            if executable == "prefetch":
                return Mock(returncode=0, stdout="", stderr="")
            return real(executable, args, **kwargs)

        with patch("prepare.data.sra.accession.shutil.which", side_effect=lambda tool: f"/usr/bin/{tool}"):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=tools):
                success, message = download_accession("SRR123", output_folder, compress=False)

        assert success is False
        assert message.startswith("network: ")
        assert "prefetch wrote no archive for SRR123" in message
        assert [c[0] for c in state.get("calls", [])] == []

    def test_failed_attempt_removes_its_temp_folder(self, tmp_path):
        output_folder = tmp_path / "downloads"

        def failing(executable, args, **kwargs):
            if executable == "fasterq-dump":
                out_dir = Path(args[args.index("-O") + 1])
                out_dir.mkdir(parents=True, exist_ok=True)
                (out_dir / "SRR123_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
                raise subprocess.CalledProcessError(1, ["fasterq-dump"], stderr="err: connection reset")
            return Mock(returncode=0, stdout="", stderr="")

        with patch("prepare.data.sra.accession.shutil.which", return_value=None):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=failing):
                success, _ = download_accession("SRR123", output_folder, compress=False)

        assert success is False
        assert not (output_folder / "SRR123_temp").exists()

    def test_keep_temp_keeps_the_temp_folder_and_logs_its_path(self, tmp_path, caplog):
        output_folder = tmp_path / "downloads"

        def failing(executable, args, **kwargs):
            out_dir = Path(args[args.index("-O") + 1])
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "SRR123_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            raise subprocess.CalledProcessError(1, ["fasterq-dump"], stderr="err: connection reset")

        with caplog.at_level(logging.INFO, logger="prepare.data.sra"):
            with patch("prepare.data.sra.accession.shutil.which", return_value=None):
                with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=failing):
                    success, _ = download_accession("SRR123", output_folder, compress=False, keep_temp=True)

        assert success is False
        assert (output_folder / "SRR123_temp" / "SRR123_1.fastq").exists()
        assert any(str(output_folder / "SRR123_temp") in r.getMessage() for r in caplog.records)

    def test_download_accession_falls_back_without_prefetch_binary(self, tmp_path):
        """shutil.which('prefetch') is None: fall back to the direct fasterq-dump call."""
        output_folder = tmp_path / "downloads"
        state = {"reads": 4}

        with patch("prepare.data.sra.accession.shutil.which", return_value=None):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
                success, message = download_accession("SRR123", output_folder, compress=False)

        assert success is True, message
        executables = [c[0] for c in state["calls"]]
        assert executables == ["fasterq-dump"]

        dump_args = state["calls"][0][1]
        assert "SRR123" in dump_args
        assert "--split-3" in dump_args
        assert "--skip-technical" in dump_args
        # No cache directory was ever created.
        assert not (output_folder / ".sra-cache").exists()
        assert accession_has_fastq(output_folder / "SRR123")

    def test_download_accession_logs_the_direct_call_when_prefetch_is_missing(self, tmp_path, caplog, monkeypatch):
        """The fallback changes where the data comes from: INFO for the first accession, then item level."""
        import threading

        from prepare.data.sra import accession as accession_mod

        monkeypatch.setattr(accession_mod, "_PREFETCH_FALLBACK_LOGGED", threading.Event())
        state = {"reads": 4}

        with patch("prepare.data.sra.accession.shutil.which", return_value=None):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
                with caplog.at_level("DEBUG", logger="prepare.data.sra"):
                    for accession in ("SRR123", "SRR124"):
                        success, message = download_accession(accession, tmp_path / "downloads", compress=False)
                        assert success is True, message

        levels = [r.levelname for r in caplog.records if "prefetch not found on PATH" in r.getMessage()]
        assert levels == ["INFO", "DEBUG"]

    def test_download_accession_dumps_a_sralite_archive(self, tmp_path):
        """Some runs are served only as .sralite; fasterq-dump must be pointed at that file."""
        output_folder = tmp_path / "downloads"
        state = {"reads": 4, "sralite": True}

        with patch("prepare.data.sra.accession.shutil.which", side_effect=lambda tool: f"/usr/bin/{tool}"):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
                success, message = download_accession("SRR123", output_folder, compress=False)

        assert success is True, message
        dump_args = state["calls"][1][1]
        assert dump_args[-1] == str(output_folder / ".sra-cache" / "SRR123" / "SRR123.sralite")
        assert accession_has_fastq(output_folder / "SRR123")

    def test_download_accession_use_prefetch_false_skips_prefetch(self, tmp_path):
        """use_prefetch=False falls back to the direct call even when prefetch is on PATH."""
        output_folder = tmp_path / "downloads"
        state = {"reads": 4}

        with patch("prepare.data.sra.accession.shutil.which", side_effect=lambda tool: f"/usr/bin/{tool}"):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
                success, message = download_accession("SRR123", output_folder, compress=False, use_prefetch=False)

        assert success is True, message
        executables = [c[0] for c in state["calls"]]
        assert executables == ["fasterq-dump"]

    def test_download_accession_single_end_output(self, tmp_path):
        """fasterq-dump --split-3 on a single-end run writes <acc>.fastq; that must be recognized."""
        output_folder = tmp_path / "downloads"
        state = {"reads": 4, "single": True}

        with patch("prepare.data.sra.accession.shutil.which", return_value=None):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
                success, message = download_accession("SRR123", output_folder, compress=False)

        assert success is True, message
        assert (output_folder / "SRR123" / "SRR123.fastq").exists()
        assert accession_has_fastq(output_folder / "SRR123")

    def test_download_accession_compresses_with_python_gzip_fallback(self, tmp_path):
        """No pigz on PATH: compression falls back to Python gzip; read count survives."""
        output_folder = tmp_path / "downloads"
        state = {"reads": 4}

        with patch("prepare.data.sra.accession.shutil.which", return_value=None):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
                success, message = download_accession("SRR123", output_folder, compress=True)

        assert success is True, message
        acc_dir = output_folder / "SRR123"
        r1 = acc_dir / "SRR123_1.fastq.gz"
        r2 = acc_dir / "SRR123_2.fastq.gz"
        assert r1.exists() and r2.exists()
        assert not (acc_dir / "SRR123_1.fastq").exists()
        assert count_fastq_reads(r1) == 4
        assert accession_has_fastq(acc_dir)

    def test_download_accession_compresses_with_pigz_when_available(self, tmp_path):
        """pigz on PATH: compress_fastq shells out to pigz -p N -f <file> instead of Python gzip."""
        output_folder = tmp_path / "downloads"
        state = {"reads": 4}

        def which(tool):
            return None if tool == "prefetch" else f"/usr/bin/{tool}"

        with patch("prepare.data.sra.accession.shutil.which", side_effect=which):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
                success, message = download_accession("SRR123", output_folder, compress=True, num_threads=4)

        assert success is True, message
        pigz_calls = [c for c in state["calls"] if c[0] == "pigz"]
        assert len(pigz_calls) == 2
        for _, args in pigz_calls:
            assert args[0] == "-p"
            assert args[1] == "4"
            assert "-f" in args

        acc_dir = output_folder / "SRR123"
        assert (acc_dir / "SRR123_1.fastq.gz").exists()
        assert (acc_dir / "SRR123_2.fastq.gz").exists()


class TestCompressFastq:
    """Test the module-level compress_fastq helper directly."""

    def test_compress_fastq_python_fallback(self, tmp_path):
        path = tmp_path / "SRR1_1.fastq"
        path.write_text("@r\nACGT\n+\nIIII\n" * 3)

        with patch("prepare.data.sra.fastq.shutil.which", return_value=None):
            result = compress_fastq(path, threads=4)

        assert result == tmp_path / "SRR1_1.fastq.gz"
        assert result.exists()
        assert not path.exists()
        assert count_fastq_reads(result) == 3

    def test_compress_fastq_uses_pigz_when_available(self, tmp_path):
        path = tmp_path / "SRR1_1.fastq"
        path.write_text("@r\nACGT\n+\nIIII\n")
        state = {}

        with patch("prepare.data.sra.fastq.shutil.which", return_value="/usr/bin/pigz"):
            with patch(
                "prepare.data.sra.fastq.SecureSubprocess.run_secure", side_effect=_fake_tools(state)
            ) as mock_run:
                result = compress_fastq(path, threads=8)

        assert result == tmp_path / "SRR1_1.fastq.gz"
        mock_run.assert_called_once()
        executable, args = mock_run.call_args[0][0], mock_run.call_args[0][1]
        assert executable == "pigz"
        assert args == ["-p", "8", "-f", str(path)]

    def test_compress_fastq_python_fallback_atomic_on_mid_write_failure(self, tmp_path):
        """A write failure partway through must not leave a truncated .gz beside the plain file."""
        path = tmp_path / "SRR1_1.fastq"
        path.write_text("@r\nACGT\n+\nIIII\n" * 3)

        real_gzip_file = gzip.GzipFile
        write_calls = {"n": 0}

        def flaky_open(*args, **kwargs):
            handle = real_gzip_file(*args, **kwargs)
            real_write = handle.write

            def flaky_write(data):
                write_calls["n"] += 1
                if write_calls["n"] == 1:
                    raise OSError("disk full")
                return real_write(data)

            handle.write = flaky_write
            return handle

        with patch("prepare.data.sra.fastq.shutil.which", return_value=None):
            with patch("prepare.data.sra.fastq.gzip.GzipFile", side_effect=flaky_open):
                with pytest.raises(OSError):
                    compress_fastq(path, threads=4)

        # The plain file survives untouched, no partial .gz and no leftover temp file.
        assert path.exists()
        assert fastq_files(tmp_path) == [path]
        assert not (tmp_path / "SRR1_1.fastq.gz").exists()
        assert sorted(p.name for p in tmp_path.iterdir()) == ["SRR1_1.fastq"]

    def test_compress_fastq_pigz_failure_leaves_no_partial_gz(self, tmp_path):
        """If pigz (via run_secure) raises, any partial .gz it left behind is removed."""
        path = tmp_path / "SRR1_1.fastq"
        path.write_text("@r\nACGT\n+\nIIII\n")
        target = tmp_path / "SRR1_1.fastq.gz"

        def failing_run(executable, args, **kwargs):
            # Simulate pigz having written a partial .gz before dying mid-compression.
            target.write_bytes(b"partial-garbage")
            raise SecurityError("pigz failed")

        with patch("prepare.data.sra.fastq.shutil.which", return_value="/usr/bin/pigz"):
            with patch("prepare.data.sra.fastq.SecureSubprocess.run_secure", side_effect=failing_run):
                with pytest.raises(SecurityError):
                    compress_fastq(path, threads=4)

        assert path.exists()
        assert not target.exists()


class TestClassifyDownloadError:
    """Test classify_download_error pure helper."""

    @pytest.mark.parametrize(
        "text,expected_class",
        [
            ("Connection timed out", "network"),
            (
                "Command timed out after 1 s (set --timeout or PREPARE_TIMEOUT; 0 disables)",
                "network",
            ),
            ("curl: (7) Failed to connect to host", "network"),
            ("Could not resolve host ftp.ncbi.nlm.nih.gov", "network"),
            ("SSL handshake failed", "network"),
            ("network is unreachable", "network"),
            ("No space left on device", "disk-full"),
            ("write failed: ENOSPC", "disk-full"),
            ("Disk full while writing output", "disk-full"),
            ("SRR000000 not found", "not-found"),
            ("Invalid accession format", "not-found"),
            ("403 Forbidden", "not-found"),
            ("404 Not Found", "not-found"),
            ("accession does not exist", "not-found"),
            ("Some completely unrelated failure", "unknown"),
            ("", "unknown"),
            # Messages as sra-tools actually prints them.
            (
                "prefetch.3.1.1 err: name not found while resolving query within virtual file "
                "system module - failed to resolve accession 'SRR1' - no data ( 404 )",
                "not-found",
            ),
            (
                "fasterq-dump.3.1.1 err: storage exhausted while writing file within file system module",
                "disk-full",
            ),
            ("fasterq-dump.3.1.1 fatal: disk-limit exeeded!", "disk-full"),
            # A class prefix classifies back to its own class, so a message that has already
            # been prefixed once is not reclassified into a different class on the retry path.
            ("disk-full: not attempted", "disk-full"),
            ("not-found: Download failed: no data", "not-found"),
            ("network: Download failed: connection reset by peer", "network"),
            ("unknown: Download failed: something else entirely", "unknown"),
            # A lock message names a pid, a time and a path that can contain 403 or 404.
            (
                "locked: Accession SRR1 is locked by pid 40404 on host node403 since "
                "2026-09-30T10:04:03.404000+00:00: /data/fastq/.locks/SRR1.lock",
                "unknown",
            ),
            ("Retry 2: locked: Accession SRR1 is locked by pid 404 on host h since t: x.lock", "unknown"),
            ("lock lost: Lost the lock /p/404/SRR1.lock: now held by pid 403", "unknown"),
            ("store error: Store dataset SRR1 is locked by pid 404 on host h since t: l", "unknown"),
            # Error codes count only as whole numbers.
            ("Download failed after 14041 bytes", "unknown"),
        ],
    )
    def test_classify_download_error_table(self, text, expected_class):
        assert classify_download_error(text) == expected_class

    def test_classify_download_error_case_insensitive(self):
        assert classify_download_error("CONNECTION TIMED OUT") == "network"
        assert classify_download_error("NO SPACE LEFT") == "disk-full"
        assert classify_download_error("NOT FOUND") == "not-found"


class TestDefaultMaxWorkers:
    """Test default_max_workers pure helper."""

    def test_default_max_workers_scales_with_cpu_and_threads(self, monkeypatch):
        monkeypatch.setattr("prepare.utils.resources.available_cpus", lambda: 8)
        assert default_max_workers(4) == 2

    def test_default_max_workers_capped_by_default_max_workers_constant(self, monkeypatch):
        monkeypatch.setattr("prepare.utils.resources.available_cpus", lambda: 64)
        assert default_max_workers(1) == 4

    def test_default_max_workers_capped_by_max_concurrent_downloads(self, monkeypatch):
        monkeypatch.setattr("prepare.utils.resources.available_cpus", lambda: 1)
        assert default_max_workers(1) <= 10

    def test_default_max_workers_cap_comes_from_the_setting(self, monkeypatch):
        monkeypatch.setattr("prepare.utils.resources.available_cpus", lambda: 64)
        monkeypatch.setenv("PREPARE_MAX_WORKERS_CAP", "8")
        assert default_max_workers(1) == 8

    def test_default_max_workers_at_least_one(self, monkeypatch):
        monkeypatch.setattr("prepare.utils.resources.available_cpus", lambda: 2)
        assert default_max_workers(16) == 1


class TestIsTransientFolder:
    """Test is_transient_folder pure helper."""

    @pytest.mark.parametrize(
        "name,expected",
        [
            ("SRR123_temp", True),
            ("SRR000000_temp", True),
            (".sra-cache", True),
            ("SRR123", False),
            ("SRR123_temporary", False),
            ("temp", False),
            ("SRR123.sra-cache", False),
        ],
    )
    def test_is_transient_folder(self, name, expected):
        assert is_transient_folder(name) == expected


class TestTransientBytes:
    """Test transient_bytes: total size of .sra-cache and <ACC>_temp folders under one folder."""

    def test_missing_folder_is_zero(self, tmp_path):
        assert transient_bytes(tmp_path / "does-not-exist") == 0

    def test_empty_folder_is_zero(self, tmp_path):
        assert transient_bytes(tmp_path) == 0

    def test_sums_sra_cache_contents(self, tmp_path):
        cache = tmp_path / ".sra-cache"
        cache.mkdir()
        (cache / "SRR1.sra").write_bytes(b"x" * 100)
        (cache / "SRR2.sra").write_bytes(b"x" * 50)
        assert transient_bytes(tmp_path) == 150

    def test_sums_acc_temp_contents(self, tmp_path):
        temp = tmp_path / "SRR1_temp"
        temp.mkdir()
        (temp / "SRR1.fastq").write_bytes(b"x" * 200)
        assert transient_bytes(tmp_path) == 200

    def test_sums_across_multiple_transient_folders(self, tmp_path):
        cache = tmp_path / ".sra-cache"
        cache.mkdir()
        (cache / "SRR1.sra").write_bytes(b"x" * 100)
        temp = tmp_path / "SRR2_temp"
        temp.mkdir()
        (temp / "SRR2.fastq").write_bytes(b"x" * 200)
        assert transient_bytes(tmp_path) == 300

    def test_ignores_non_transient_entries(self, tmp_path):
        real = tmp_path / "SRR1"
        real.mkdir()
        (real / "SRR1_1.fastq").write_bytes(b"x" * 999)
        assert transient_bytes(tmp_path) == 0

    def test_recurses_into_subdirectories(self, tmp_path):
        cache = tmp_path / ".sra-cache"
        nested = cache / "nested"
        nested.mkdir(parents=True)
        (nested / "SRR1.sra").write_bytes(b"x" * 42)
        assert transient_bytes(tmp_path) == 42


class TestCheckExistingDownloads:
    """Test _check_existing_downloads function."""

    def test_check_existing_downloads_none_exist(self, tmp_path):
        """Test when no downloads exist."""
        output_folder = tmp_path / "downloads"
        accessions = ["SRR123", "SRR456"]

        already_downloaded, to_download, blacklisted = _check_existing_downloads(accessions, output_folder, force=False)

        assert to_download == accessions
        assert already_downloaded == []

    def test_check_existing_downloads_some_exist(self, tmp_path):
        """Test when some downloads exist."""
        output_folder = tmp_path / "downloads"
        (output_folder / "SRR123").mkdir(parents=True)
        (output_folder / "SRR123" / "SRR123.fastq").write_text("data")

        accessions = ["SRR123", "SRR456"]

        already_downloaded, to_download, blacklisted = _check_existing_downloads(accessions, output_folder, force=False)

        assert to_download == ["SRR456"]
        assert already_downloaded == ["SRR123"]

    def test_check_existing_downloads_force_all(self, tmp_path):
        """Test with force=True."""
        output_folder = tmp_path / "downloads"
        (output_folder / "SRR123").mkdir(parents=True)
        (output_folder / "SRR123" / "SRR123.fastq").write_text("data")

        accessions = ["SRR123", "SRR456"]

        with patch("shutil.rmtree"):
            already_downloaded, to_download, blacklisted = _check_existing_downloads(
                accessions, output_folder, force=True
            )

        assert to_download == accessions
        assert already_downloaded == []

    def test_check_existing_downloads_truncated_accession_is_not_already_downloaded(self, tmp_path):
        """An accession with a 'truncated' registry verdict must be redownloaded, not skipped."""
        output_folder = tmp_path / "downloads"
        (output_folder / "SRR123").mkdir(parents=True)
        (output_folder / "SRR123" / "SRR123.fastq").write_text("partial")

        accessions = ["SRR123", "SRR456"]

        already_downloaded, to_download, blacklisted = _check_existing_downloads(
            accessions, output_folder, force=False, truncated_accessions={"SRR123"}
        )

        assert to_download == ["SRR123", "SRR456"]
        assert already_downloaded == []


class TestDownloadSra:
    """Test download_sra function."""

    def test_download_sra_refuses_traversal_before_touching_disk(self, tmp_path):
        victim = tmp_path / "victim"
        victim.mkdir()
        (tmp_path / "acc.txt").write_text("../../victim\n")
        fastq = tmp_path / "proj" / "fastq"
        with pytest.raises(ValidationError):
            download_sra(fastq, tmp_path / "acc.txt")
        assert victim.is_dir()
        assert not (tmp_path / "proj").exists()

    def test_download_sra_from_file(self, tmp_path):
        """Test downloading SRA data from accessions file."""
        accessions_file = tmp_path / "accessions.txt"
        output_folder = tmp_path / "downloads"
        accessions_file.write_text("SRR123\nSRR456\n")

        with patch("prepare.data.sra.download._check_existing_downloads") as mock_check:
            with patch("concurrent.futures.ThreadPoolExecutor") as mock_executor:
                with patch("prepare.data.sra.retry._process_download_results") as mock_process:
                    with patch("prepare.data.sra.retry._handle_download_failure"):
                        mock_check.return_value = ([], ["SRR123", "SRR456"], [])
                        mock_executor_instance = Mock()
                        mock_executor.return_value.__enter__.return_value = mock_executor_instance
                        mock_executor_instance.submit.return_value = Mock()
                        mock_process.return_value = (2, 0)

                        result = download_sra(output_folder, accessions_file, num_threads=4, max_workers=2)

        assert isinstance(result, dict)
        mock_check.assert_called_once()

    def test_download_sra_propagates_abort_reason(self, tmp_path):
        """A disk-full abort reported by the retry pass surfaces as stats['aborted']."""
        accessions_file = tmp_path / "accessions.txt"
        output_folder = tmp_path / "downloads"
        accessions_file.write_text("SRR123\n")

        with patch("prepare.data.sra.download._check_existing_downloads") as mock_check:
            with patch("prepare.data.sra.retry._download_with_retries") as mock_retries:
                with patch("prepare.data.sra.retry._handle_download_failure"):
                    mock_check.return_value = ([], ["SRR123"], [])
                    mock_retries.return_value = (
                        0,
                        1,
                        ["SRR123"],
                        {"SRR123": "disk-full: not attempted"},
                        "disk-full",
                    )

                    result = download_sra(output_folder, accessions_file)

        assert result["aborted"] == "disk-full"
        assert result["failed"] == 1

    def test_download_sra_aborted_is_none_when_run_completes_normally(self, tmp_path):
        accessions_file = tmp_path / "accessions.txt"
        output_folder = tmp_path / "downloads"
        accessions_file.write_text("SRR123\n")

        with patch("prepare.data.sra.download._check_existing_downloads") as mock_check:
            with patch("prepare.data.sra.retry._download_with_retries") as mock_retries:
                with patch("prepare.data.sra.retry._handle_download_failure"):
                    mock_check.return_value = ([], ["SRR123"], [])
                    mock_retries.return_value = (1, 0, [], {"SRR123": "Downloaded 1 files"}, None)

                    result = download_sra(output_folder, accessions_file)

        assert result["aborted"] is None

    def test_download_sra_with_blacklist(self, tmp_path):
        """Test downloading with blacklist filtering."""
        accessions_file = tmp_path / "accessions.txt"
        blacklist_file = tmp_path / "blacklist.txt"
        output_folder = tmp_path / "downloads"

        accessions_file.write_text("SRR123\nSRR456\nSRR789\n")
        blacklist_file.write_text("SRR456\n")

        with patch("prepare.data.sra.download._check_existing_downloads") as mock_check:
            with patch("concurrent.futures.ThreadPoolExecutor") as mock_executor:
                with patch("prepare.data.sra.retry._process_download_results") as mock_process:
                    with patch("prepare.data.sra.retry._handle_download_failure"):
                        mock_check.return_value = ([], ["SRR123", "SRR789"], [])
                        mock_executor_instance = Mock()
                        mock_executor.return_value.__enter__.return_value = mock_executor_instance
                        mock_process.return_value = (2, 0)

                        download_sra(output_folder, accessions_file, blacklist=[blacklist_file])

        # Verify that blacklist filtering works correctly
        check_call_args = mock_check.call_args
        all_accessions = check_call_args[0][0]  # First positional argument
        blacklisted_accessions = check_call_args[0][3]  # Fourth positional argument (set of blacklisted)

        # All accessions should be passed to the function
        assert "SRR123" in all_accessions
        assert "SRR456" in all_accessions
        assert "SRR789" in all_accessions
        # But SRR456 should be in the blacklist set
        assert "SRR456" in blacklisted_accessions

    def test_download_sra_nonexistent_file(self, tmp_path):
        """Test error handling for nonexistent accessions file."""
        nonexistent_file = tmp_path / "nonexistent.txt"
        output_folder = tmp_path / "downloads"

        with pytest.raises(DataAccessError):
            download_sra(output_folder, nonexistent_file)

    def test_download_sra_dry_run(self, tmp_path):
        """Test dry run mode."""
        accessions_file = tmp_path / "accessions.txt"
        output_folder = tmp_path / "downloads"
        accessions_file.write_text("SRR123\nSRR456\n")

        with patch("prepare.data.sra.download._check_existing_downloads") as mock_check:
            with patch("prepare.data.sra.download.logger") as mock_logger:
                mock_check.return_value = ([], ["SRR123", "SRR456"], [])
                result = download_sra(output_folder, accessions_file, dry_run=True)

        expected_result = {
            "total": 2,
            "already_downloaded": 0,
            "blacklisted": 0,
            "to_download": 2,
            "successful": 0,
            "failed": 0,
            "already_downloaded_accessions": [],
            "blacklisted_accessions": [],
            "skipped_accessions": [],
            "aborted": None,
        }
        assert result == expected_result
        mock_logger.info.assert_called_with("Dry run: would download 2 accessions")

    def test_download_sra_returns_accession_lists(self, tmp_path):
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\nSRR2\n")
        (tmp_path / "fastq" / "SRR2").mkdir(parents=True)
        (tmp_path / "fastq" / "SRR2" / "SRR2_1.fastq").write_text("@r\nA\n+\nI\n")
        black = tmp_path / "black.txt"
        black.write_text("SRR1\n")
        stats = download_sra(tmp_path / "fastq", acc, dry_run=True, blacklist=[black])
        assert stats["already_downloaded_accessions"] == ["SRR2"]
        assert stats["blacklisted_accessions"] == ["SRR1"]

    def test_on_result_called_on_main_thread_per_accession(self, tmp_path, monkeypatch):
        import threading

        seen = []
        main = threading.get_ident()

        def on_result(acc, ok, message):
            seen.append((acc, ok, threading.get_ident() == main))

        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\nSRR2\n")
        with patch(
            "prepare.data.sra.accession.download_accession",
            side_effect=_staging_download((True, "Downloaded 2 files"), (False, "Download failed: x")),
        ):
            stats = download_sra(tmp_path / "fastq", acc, max_workers=2, max_retries=0, on_result=on_result)
        assert sorted(a for a, _, _ in seen) == ["SRR1", "SRR2"] and all(on_main for _, _, on_main in seen)
        assert stats["failed_accessions"] == ["SRR2"] or stats["failed_accessions"] == ["SRR1"]

    def test_on_result_failure_does_not_corrupt_the_tally(self, tmp_path):
        """A raising on_result (e.g. a registry lock timeout) must not be mistaken for a download failure."""
        calls = []

        def flaky_on_result(acc, ok, message):
            calls.append(acc)
            if len(calls) == 1:
                raise RuntimeError("registry locked")

        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\nSRR2\n")
        with patch(
            "prepare.data.sra.accession.download_accession",
            side_effect=_staging_download((True, "Downloaded 2 files")),
        ):
            stats = download_sra(tmp_path / "fastq", acc, max_retries=0, on_result=flaky_on_result)
        assert stats["successful"] == 2
        assert stats["failed"] == 0
        assert stats["failed_accessions"] == []
        assert len(calls) == 2

    def test_max_downloads_cutoffs_are_returned(self, tmp_path):
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\nSRR2\nSRR3\n")
        with patch("prepare.data.sra.accession.download_accession", return_value=(True, "Downloaded 2 files")):
            stats = download_sra(tmp_path / "fastq", acc, max_downloads=1)
        assert len(stats["skipped_accessions"]) == 2

    def test_registry_exclusions_join_the_blacklist(self, tmp_path):
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\nSRR2\n")
        stats = download_sra(tmp_path / "fastq", acc, dry_run=True, blacklist_accessions={"SRR2"})
        assert stats["blacklisted_accessions"] == ["SRR2"]

    def test_expected_spots_reaches_download_accession_per_accession(self, tmp_path):
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\nSRR2\n")
        calls = {}

        def fake_download_accession(accession, *args, **kwargs):
            calls[accession] = kwargs.get("expected_spots")
            return True, "Downloaded 1 files, unverified"

        with patch("prepare.data.sra.accession.download_accession", side_effect=fake_download_accession):
            download_sra(tmp_path / "fastq", acc, expected_spots={"SRR1": 1000})

        assert calls == {"SRR1": 1000, "SRR2": None}

    def test_download_sra_records_verdict_through_on_result(self, tmp_path, monkeypatch):
        """A real download_accession run (files created, verified) reports its verdict via on_result."""
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        fastq_folder = tmp_path / "fastq"

        def fake_run_secure(executable, args, **kwargs):
            if executable == "pigz":
                # The runner may have pigz installed (GitHub's Ubuntu images do); behave like it.
                plain = Path(args[-1])
                with open(plain, "rb") as src, gzip.open(str(plain) + ".gz", "wb") as dest:
                    shutil.copyfileobj(src, dest)
                plain.unlink()
                return Mock(returncode=0, stdout="", stderr="")
            out_dir = Path(args[args.index("-O") + 1])
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            return Mock(returncode=0, stdout="", stderr="")

        seen = []
        monkeypatch.setattr("prepare.utils.security.SecureSubprocess.run_secure", fake_run_secure)
        monkeypatch.setattr("prepare.data.sra.accession.shutil.which", lambda tool: None)
        monkeypatch.setattr("prepare.data.sra.fastq.count_fastq_reads", lambda path: 1)

        download_sra(
            fastq_folder,
            acc,
            expected_spots={"SRR1": 2},
            on_result=lambda a, ok, msg: seen.append((a, ok, msg)),
        )

        assert seen[0][:2] == ("SRR1", False)
        assert seen[0][2].startswith("Downloaded 1 files, truncated (1 of 2 spots)")
        assert "--accept-partial" in seen[0][2]
        # The retry pass ran: the second outcome carries the retry prefix.
        assert len(seen) == 2 and seen[1][2].startswith("Retry 1: ")
        assert seen[1][1] is False

    def test_truncated_download_is_a_success_with_accept_partial(self, tmp_path, monkeypatch):
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        fastq_folder = tmp_path / "fastq"

        def fake_run_secure(executable, args, **kwargs):
            out_dir = Path(args[args.index("-O") + 1])
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            return Mock(returncode=0, stdout="", stderr="")

        seen = []
        monkeypatch.setattr("prepare.utils.security.SecureSubprocess.run_secure", fake_run_secure)
        monkeypatch.setattr("prepare.data.sra.accession.shutil.which", lambda tool: None)
        monkeypatch.setattr("prepare.data.sra.fastq.count_fastq_reads", lambda path: 1)

        stats = download_sra(
            fastq_folder,
            acc,
            expected_spots={"SRR1": 2},
            compress=False,
            accept_partial=True,
            on_result=lambda a, ok, msg: seen.append((a, ok, msg)),
        )

        assert seen == [("SRR1", True, "Downloaded 1 files, truncated (1 of 2 spots)")]
        assert stats["successful"] == 1
        assert (fastq_folder / "SRR1" / "SRR1_1.fastq").exists()

    def test_truncated_accessions_are_redownloaded_not_skipped(self, tmp_path):
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        fastq_folder = tmp_path / "fastq"
        (fastq_folder / "SRR1").mkdir(parents=True)
        (fastq_folder / "SRR1" / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")

        with patch(
            "prepare.data.sra.accession.download_accession",
            side_effect=_staging_download((True, "Downloaded 1 files, complete (2 of 2 spots)")),
        ) as mock_download:
            stats = download_sra(fastq_folder, acc, truncated_accessions={"SRR1"})

        mock_download.assert_called_once()
        assert stats["already_downloaded"] == 0
        assert stats["successful"] == 1


class TestFasterqDumpVersion:
    """The tool version recorded in a store dataset's sidecar."""

    def test_returns_the_last_non_empty_line(self):
        from prepare.core.constants import VERSION_PROBE_TIMEOUT

        with (
            patch("prepare.utils.tools.shutil.which", return_value="/usr/bin/fasterq-dump"),
            patch(
                "prepare.utils.security.SecureSubprocess.run_secure",
                return_value=Mock(stdout="\nfasterq-dump : 3.0.10\n", stderr="", returncode=0),
            ) as mock_run,
        ):
            assert fasterq_dump_version() == "fasterq-dump : 3.0.10"
        mock_run.assert_called_once_with("fasterq-dump", ["--version"], timeout=VERSION_PROBE_TIMEOUT, check=False)

    def test_returns_empty_when_the_tool_cannot_be_run(self):
        with (
            patch("prepare.utils.tools.shutil.which", return_value="/usr/bin/fasterq-dump"),
            patch("prepare.utils.security.SecureSubprocess.run_secure", side_effect=SecurityError("not installed")),
        ):
            assert fasterq_dump_version() == ""

    def test_returns_empty_when_the_tool_is_not_on_path(self):
        with (
            patch("prepare.utils.tools.shutil.which", return_value=None),
            patch("prepare.utils.security.SecureSubprocess.run_secure") as mock_run,
        ):
            assert fasterq_dump_version() == ""
        mock_run.assert_not_called()


class TestSafeRmtreeIgnoresMissingFiles:
    """On a volume storing each file's AppleDouble sidecar (``._<name>``) next to it, macOS can
    delete ``._X`` together with ``X``; ``_safe_rmtree`` must not treat that race as a failure."""

    def test_a_file_removed_concurrently_during_cleanup_is_not_logged(self, tmp_path, caplog):
        target = tmp_path / "acc_dir"
        target.mkdir()
        (target / "keep.fastq").write_text("@r\nACGT\n+\nIIII\n")
        missing_name = "._keep.fastq"
        (target / missing_name).write_bytes(b"\x00\x05\x16\x07")

        real_unlink = os.unlink

        def _flaky_unlink(path, *args, **kwargs):
            # `os.unlink` is called with just the entry name (plus `dir_fd`) on a platform
            # using shutil's fd-based rmtree, not the full path, so the match is by basename.
            # The race being simulated is "the sidecar is already gone when rmtree reaches
            # it": remove the file for real, then report it missing. Leaving it on disk would
            # make the final rmdir fail with "directory not empty" on every platform whose
            # rmtree removes entries in scandir order (Linux CI), which is not the race.
            if os.path.basename(str(path)) == missing_name:
                real_unlink(path, *args, **kwargs)
                raise FileNotFoundError(2, "No such file or directory", str(path))
            return real_unlink(path, *args, **kwargs)

        with caplog.at_level(logging.WARNING):
            with patch("os.unlink", side_effect=_flaky_unlink):
                _safe_rmtree(target)

        assert not target.exists()
        assert not any("Could not remove directory" in r.message for r in caplog.records)

    def test_a_genuine_removal_failure_is_still_logged(self, tmp_path, caplog):
        target = tmp_path / "acc_dir"
        target.mkdir()
        (target / "keep.fastq").write_text("@r\nACGT\n+\nIIII\n")

        with caplog.at_level(logging.WARNING):
            with patch("os.unlink", side_effect=PermissionError(13, "Permission denied")):
                _safe_rmtree(target)

        assert any("Could not remove directory" in r.message for r in caplog.records)

    def test_ignore_missing_reraises_a_permission_error_given_directly(self):
        # shutil.rmtree's onexc callback (3.12+) passes the exception object itself, not a
        # sys.exc_info() tuple; only FileNotFoundError is swallowed, anything else re-raises.
        with pytest.raises(PermissionError):
            _ignore_missing(os.unlink, "some/path", PermissionError(13, "Permission denied"))


class TestDownloadInterrupt:
    @pytest.fixture(autouse=True)
    def _clear_stop(self):

        accession_mod.STOP.clear()
        SecureSubprocess.clear_stopping()
        yield
        accession_mod.STOP.clear()
        SecureSubprocess.clear_stopping()

    def test_keyboard_interrupt_during_submission_reaches_the_handler(self, tmp_path):

        def accessions():
            yield "SRR1"
            raise KeyboardInterrupt

        worker = Mock(return_value=(True, "ok"))
        stop = threading.Event()
        with patch.object(SecureSubprocess, "terminate_children", return_value=0) as term:
            with pytest.raises(KeyboardInterrupt):
                retry_mod._execute_parallel_downloads(
                    accessions(), tmp_path, 1, 1, False, None, {}, [], downloader=worker, stop=stop
                )
        term.assert_called_once_with(stop=stop)
        assert stop.is_set()
        assert not accession_mod.STOP.is_set()

    def test_a_new_run_leaves_the_process_wide_stopping_flag_set(self, tmp_path):

        SecureSubprocess.terminate_children(grace=0.0)
        worker = Mock(return_value=(True, "ok"))
        retry_mod._execute_parallel_downloads(["SRR1"], tmp_path, 1, 1, False, None, {}, [], downloader=worker)
        assert SecureSubprocess._stopping is True

    def test_tool_killed_by_the_interrupt_reports_interrupted(self, tmp_path):
        import subprocess as sp

        def killed(executable, args, **kwargs):
            accession_mod.STOP.set()
            raise sp.CalledProcessError(-15, [executable], output="", stderr="")

        with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=killed):
            with patch("prepare.data.sra.accession.shutil.which", return_value=None):
                success, message = accession_mod.download_accession("SRR2517620", tmp_path / "fastq")
        assert (success, message) == (False, "interrupted")

    def test_tool_failure_without_interrupt_is_still_an_error(self, tmp_path):
        import subprocess as sp

        err = sp.CalledProcessError(3, ["fasterq-dump"], output="", stderr="network timeout")
        with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=err):
            with patch("prepare.data.sra.accession.shutil.which", return_value=None):
                success, message = accession_mod.download_accession("SRR2517620", tmp_path / "fastq")
        assert success is False
        assert message.startswith("network:")

    def test_compression_is_skipped_once_stopped(self, tmp_path):

        temp = tmp_path / "SRR1_temp"
        temp.mkdir()
        (temp / "SRR1_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (temp / "SRR1_2.fastq").write_text("@r\nACGT\n+\nIIII\n")
        accession_mod.STOP.set()
        with patch("prepare.data.sra.fastq.compress_fastq") as mock_compress:
            success, message = accession_mod._handle_download_output(temp, tmp_path / "SRR1", compress=True)
        assert (success, message) == (False, "interrupted")
        mock_compress.assert_not_called()

    def test_interrupt_during_compression_publishes_nothing(self, tmp_path):
        """The project folder does not appear and the prefetch archive stays for the next run."""
        fastq = tmp_path / "fastq"
        token = threading.Event()
        state = {"reads": 4}
        archive = fastq / ".sra-cache" / "SRR123" / "SRR123.sra"

        def compress_then_stop(path, threads, stop=None):
            token.set()
            raise InterruptedError("stopped")

        with patch("prepare.data.sra.accession.shutil.which", side_effect=lambda tool: f"/usr/bin/{tool}"):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
                with patch("prepare.data.sra.fastq.compress_fastq", side_effect=compress_then_stop):
                    success, message = accession_mod._project_download("SRR123", fastq, stop=token)

        assert (success, message) == (False, "interrupted")
        assert not (fastq / "SRR123").exists()
        assert not (fastq / ".prepare-tmp" / "SRR123").exists()
        assert archive.exists()

    def test_store_download_waiting_on_a_held_lock_returns_interrupted(self, tmp_path):
        from prepare.store.layout import init_store, lock_path

        store = init_store(tmp_path / "store")
        # A live holder in another process (another host, so it is never judged dead).
        lock = lock_path(store, "SRR1")
        lock.parent.mkdir(parents=True, exist_ok=True)
        lock.write_text(json.dumps({"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00"}))
        accession_mod.STOP.set()
        with patch("prepare.data.sra.store_handoff._store_precheck", return_value=None):
            with patch("prepare.data.sra.store_handoff._store_fetch") as fetch:
                started = time.monotonic()
                result = store_handoff_mod._store_download("SRR1", tmp_path / "fastq", store, lock_wait=0.0)
        assert result == (False, "interrupted")
        assert json.loads(lock.read_text())["pid"] == 4242
        assert time.monotonic() - started < 1.0
        fetch.assert_not_called()

    def test_keyboard_interrupt_cancels_pending_and_terminates_children(self, tmp_path):

        calls = []

        def worker(acc, *a, **kw):
            calls.append(acc)
            if acc == "SRR1":
                raise KeyboardInterrupt
            return True, "ok"

        stop = threading.Event()
        with patch.object(SecureSubprocess, "terminate_children", return_value=0) as term:
            with pytest.raises(KeyboardInterrupt):
                retry_mod._execute_parallel_downloads(
                    ["SRR1", "SRR2", "SRR3", "SRR4"], tmp_path, 1, 1, False, None, {}, [], downloader=worker, stop=stop
                )
        term.assert_called_once_with(stop=stop)
        assert stop.is_set()
        assert not accession_mod.STOP.is_set()

    def test_keyboard_interrupt_collects_downloads_that_had_already_finished(self, tmp_path):
        # SRR1 is yielded and collected; SRR2 finished too but the interrupt arrives from
        # as_completed before it is yielded. Both must reach on_result, each exactly once.
        from concurrent.futures import wait

        def interrupted_as_completed(futures):
            wait(list(futures))
            ordered = sorted(futures, key=lambda f: futures[f])
            yield ordered[0]
            raise KeyboardInterrupt

        worker = Mock(return_value=(True, "ok"))
        seen = []
        with patch.object(SecureSubprocess, "terminate_children", return_value=0):
            with patch.object(retry_mod, "as_completed", side_effect=interrupted_as_completed):
                with pytest.raises(KeyboardInterrupt):
                    retry_mod._execute_parallel_downloads(
                        ["SRR1", "SRR2"],
                        tmp_path,
                        1,
                        2,
                        False,
                        None,
                        {},
                        [],
                        on_result=lambda acc, ok, msg: seen.append((acc, ok, msg)),
                        downloader=worker,
                    )
        assert sorted(seen) == [("SRR1", True, "ok"), ("SRR2", True, "ok")]

    def test_keyboard_interrupt_skips_a_finished_worker_that_raised_a_non_download_error(self, tmp_path):
        from concurrent.futures import wait

        def interrupted_as_completed(futures):
            wait(list(futures))
            raise KeyboardInterrupt
            yield  # pragma: no cover - makes this a generator

        def worker(acc, *a, **kw):
            if acc == "SRR1":
                raise RuntimeError("programming error")
            if acc == "SRR2":
                raise OSError("disk went away")
            return True, "ok"

        seen = []
        with patch.object(SecureSubprocess, "terminate_children", return_value=0):
            with patch.object(retry_mod, "as_completed", side_effect=interrupted_as_completed):
                with pytest.raises(KeyboardInterrupt):
                    retry_mod._execute_parallel_downloads(
                        ["SRR1", "SRR2", "SRR3"],
                        tmp_path,
                        1,
                        3,
                        False,
                        None,
                        {},
                        [],
                        on_result=lambda acc, ok, msg: seen.append((acc, ok)),
                        downloader=worker,
                    )
        assert sorted(seen) == [("SRR2", False), ("SRR3", True)]

    @staticmethod
    def _interrupt_while_running(futures):
        # The interrupt arrives while every worker is still running (published but cleaning up).
        time.sleep(0.05)
        raise KeyboardInterrupt
        yield  # pragma: no cover - makes this a generator

    def _run_interrupted(self, tmp_path, accessions, worker):
        seen = []
        with patch.object(SecureSubprocess, "terminate_children", return_value=0):
            with patch.object(retry_mod, "as_completed", side_effect=self._interrupt_while_running):
                with pytest.raises(KeyboardInterrupt):
                    retry_mod._execute_parallel_downloads(
                        accessions,
                        tmp_path,
                        1,
                        len(accessions),
                        False,
                        None,
                        {},
                        [],
                        on_result=lambda acc, ok, msg: seen.append((acc, ok, msg)),
                        downloader=worker,
                    )
        return sorted(seen)

    def test_keyboard_interrupt_collects_workers_that_finish_after_the_signal(self, tmp_path):
        def worker(acc, *a, **kw):
            time.sleep(0.3)
            return True, "ok"

        seen = self._run_interrupted(tmp_path, ["SRR1", "SRR2"], worker)
        assert seen == [("SRR1", True, "ok"), ("SRR2", True, "ok")]

    def test_keyboard_interrupt_leaves_out_a_worker_that_returned_interrupted(self, tmp_path):
        def worker(acc, *a, **kw):
            time.sleep(0.3)
            if acc == "SRR1":
                return False, "interrupted"
            if acc == "SRR2":
                return False, "fasterq-dump failed: network error"
            return True, "ok"

        seen = self._run_interrupted(tmp_path, ["SRR1", "SRR2", "SRR3"], worker)
        assert seen == [("SRR2", False, "fasterq-dump failed: network error"), ("SRR3", True, "ok")]

    def test_download_accession_runs_no_tool_once_stopped(self, tmp_path):

        accession_mod.STOP.set()
        with patch("prepare.data.sra.accession.SecureSubprocess.run_secure") as mock_run:
            with patch("prepare.data.sra.accession.shutil.which", return_value="/usr/bin/prefetch"):
                success, message = accession_mod.download_accession("SRR2517620", tmp_path / "fastq")
        assert (success, message) == (False, "interrupted")
        mock_run.assert_not_called()

    def test_download_accession_skips_fasterq_dump_when_stopped_during_prefetch(self, tmp_path):

        def prefetch(executable, args, **kwargs):
            accession_mod.STOP.set()

        with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=prefetch) as mock_run:
            with patch("prepare.data.sra.accession.shutil.which", return_value="/usr/bin/prefetch"):
                with patch("prepare.data.sra.accession._cached_sra_archive", return_value=tmp_path / "SRR2517620.sra"):
                    success, message = accession_mod.download_accession("SRR2517620", tmp_path / "fastq")
        assert (success, message) == (False, "interrupted")
        assert [c.args[0] for c in mock_run.call_args_list] == ["prefetch"]

    def test_retry_pass_returns_early_once_stopped(self, tmp_path):

        accession_mod.STOP.set()
        worker = Mock(return_value=(True, "ok"))
        retried, failed, abort = retry_mod._retry_failed_downloads(
            ["SRR1"], 2, tmp_path, 1, None, {"SRR1": "Download failed: timeout"}, downloader=worker
        )
        assert (retried, failed, abort) == (0, ["SRR1"], None)
        worker.assert_not_called()


def test_fasterq_dump_version_is_empty_on_a_tool_error_and_propagates_a_bug():
    """Kind (b): the narrow error yields the default; a TypeError propagates."""
    from prepare.data.sra import fasterq_dump_version

    with patch("prepare.utils.tools.shutil.which", return_value="/usr/bin/fasterq-dump"):
        with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=SecurityError("not found")):
            assert fasterq_dump_version() == ""
        fasterq_dump_version.cache_clear()  # probed once per process; the second probe must run
        with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=TypeError("bug")):
            with pytest.raises(TypeError):
                fasterq_dump_version()


def test_accession_has_fastq_ignores_a_sidecar_that_is_not_an_object(tmp_path):
    """Kind (b): a JSON sidecar holding a list is treated as carrying no state."""
    acc = tmp_path / "SRR1"
    acc.mkdir()
    (acc / "SRR1.fastq").write_text("@r\nA\n+\nI\n")
    (acc / "SRR1.json").write_text("[1, 2]")
    assert accession_has_fastq(acc) is True


def test_data_access_error_from_download_sra_chains_the_cause(tmp_path):
    """Kind (c): the wrapped error keeps its cause."""
    with pytest.raises(DataAccessError) as exc:
        download_sra(tmp_path / "fastq", tmp_path / "missing.txt")
    assert isinstance(exc.value.__cause__.__cause__, OSError)


def test_safe_rmtree_logs_an_os_error_and_propagates_a_bug(tmp_path, caplog):
    """Kind (e): a filesystem error is logged; a TypeError propagates."""
    target = tmp_path / "d"
    target.mkdir()
    with patch("prepare.data.sra.cleanup.shutil.rmtree", side_effect=PermissionError("denied")):
        with caplog.at_level(logging.WARNING):
            _safe_rmtree(target)
    assert "Could not remove directory" in caplog.text
    with patch("prepare.data.sra.cleanup.shutil.rmtree", side_effect=TypeError("bug")):
        with pytest.raises(TypeError):
            _safe_rmtree(target)


class TestProjectLockPath:
    def test_the_lock_lives_in_a_hidden_folder_of_the_fastq_folder(self, tmp_path):
        assert (
            accession_mod.project_lock_path(tmp_path / "fastq", "SRR1") == tmp_path / "fastq" / ".locks" / "SRR1.lock"
        )


class TestProjectDownload:
    """_project_download: one download per accession across processes, published with one rename."""

    def _run(self, fastq, state, **kwargs):
        with patch("prepare.data.sra.accession.shutil.which", side_effect=lambda tool: None):
            with patch("prepare.data.sra.fastq.shutil.which", side_effect=_no_prefetch_with_pigz):
                with patch(
                    "prepare.utils.security.SecureSubprocess.run_secure", side_effect=_slow_download_tools(state)
                ):
                    return accession_mod._project_download("SRR1", fastq, 1, compress=True, **kwargs)

    def test_publishes_a_complete_compressed_folder_and_leaves_no_staging(self, tmp_path):
        fastq = tmp_path / "fastq"
        state = {}
        success, message = self._run(fastq, state)
        assert success is True, message
        assert sorted(p.name for p in (fastq / "SRR1").iterdir()) == ["SRR1_1.fastq.gz", "SRR1_2.fastq.gz"]
        assert list((fastq / ".prepare-tmp").iterdir()) == []
        assert not (fastq / ".locks" / "SRR1.lock").exists()
        assert sorted(p.name for p in fastq.iterdir() if not p.name.startswith(".")) == ["SRR1"]

    def test_an_existing_download_found_inside_the_lock_is_not_fetched_again(self, tmp_path):
        fastq = tmp_path / "fastq"
        (fastq / "SRR1").mkdir(parents=True)
        (fastq / "SRR1" / "SRR1_1.fastq.gz").write_bytes(gzip.compress(b"@r\nACGT\n+\nIIII\n"))
        state = {}
        assert self._run(fastq, state) == (True, "already exists")
        assert state.get("calls") is None

    def test_the_prefetch_cache_stays_in_the_fastq_folder(self, tmp_path):
        fastq = tmp_path / "fastq"
        seen = {}

        def _download(accession, output_folder, *args, **kwargs):
            seen["output_folder"] = Path(output_folder)
            seen["sra_cache"] = kwargs.get("sra_cache")
            seen["staging_folder"] = kwargs.get("staging_folder")
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            return True, "Downloaded 1 files, unverified"

        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            assert accession_mod._project_download("SRR1", fastq, sra_cache=None)[0] is True
        assert seen["output_folder"] == fastq / ".prepare-tmp"
        assert seen["staging_folder"] == fastq / ".prepare-tmp"
        assert Path(seen["sra_cache"]) == fastq / ".sra-cache"

    def test_a_forced_download_replaces_the_existing_folder_only_once_the_new_one_is_ready(self, tmp_path):
        fastq = tmp_path / "fastq"
        (fastq / "SRR1").mkdir(parents=True)
        (fastq / "SRR1" / "SRR1_1.fastq").write_text("@old\nACGT\n+\nIIII\n")
        seen = {}

        def _download(accession, output_folder, *args, **kwargs):
            seen["old_still_there"] = (fastq / "SRR1" / "SRR1_1.fastq").read_text().startswith("@old")
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@new\nACGT\n+\nIIII\n")
            return True, "Downloaded 1 files, unverified"

        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            assert accession_mod._project_download("SRR1", fastq, force=True)[0] is True
        assert seen["old_still_there"] is True
        assert (fastq / "SRR1" / "SRR1_1.fastq").read_text().startswith("@new")
        assert list((fastq / ".prepare-tmp").iterdir()) == []

    def test_only_accessions_in_the_truncated_set_skip_the_existing_check(self, tmp_path):
        fastq = tmp_path / "fastq"
        for accession in ("SRR1", "SRR2"):
            (fastq / accession).mkdir(parents=True)
            (fastq / accession / f"{accession}_1.fastq").write_text("@old\nACGT\n+\nIIII\n")
        calls = []

        def _download(accession, output_folder, *args, **kwargs):
            calls.append((accession, kwargs.get("redownload_truncated")))
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@new\nACGT\n+\nIIII\n")
            return True, "Downloaded 1 files, unverified"

        truncated = frozenset({"SRR1"})
        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            first = accession_mod._project_download("SRR1", fastq, truncated=truncated, redownload_truncated=True)
            second = accession_mod._project_download("SRR2", fastq, truncated=truncated, redownload_truncated=True)
        assert first == (True, "Downloaded 1 files, unverified")
        assert second == (True, "already exists")
        assert calls == [("SRR1", True)]
        assert (fastq / "SRR2" / "SRR2_1.fastq").read_text().startswith("@old")

    def test_a_failed_download_publishes_nothing(self, tmp_path):
        fastq = tmp_path / "fastq"

        def _download(accession, output_folder, *args, **kwargs):
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            return False, "network: Connection timed out"

        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            assert accession_mod._project_download("SRR1", fastq) == (False, "network: Connection timed out")
        assert not (fastq / "SRR1").exists()
        assert not (fastq / ".prepare-tmp" / "SRR1").exists()

    def test_a_success_that_staged_no_folder_is_a_failure(self, tmp_path):
        fastq = tmp_path / "fastq"
        with patch("prepare.data.sra.accession.download_accession", return_value=(True, "Downloaded 1 files")):
            success, message = accession_mod._project_download("SRR1", fastq)
        assert success is False
        assert "SRR1" in message
        assert not (fastq / "SRR1").exists()

    def test_a_lost_lock_publishes_nothing_and_leaves_the_staged_folder_to_the_new_holder(self, tmp_path):
        from prepare.utils.lockfile import LockLost

        fastq = tmp_path / "fastq"

        def _download(accession, output_folder, *args, **kwargs):
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            return True, "Downloaded 1 files, unverified"

        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            with patch("prepare.data.sra.accession.verify_held", side_effect=LockLost("Lost the lock x: removed")):
                success, message = accession_mod._project_download("SRR1", fastq)
        assert success is False
        assert "Lost the lock" in message
        assert not (fastq / "SRR1").exists()
        assert (fastq / ".prepare-tmp" / "SRR1").is_dir()

    def test_two_threads_on_one_accession_download_it_once(self, tmp_path, fast_project_lock_poll):
        import threading as threading_module

        fastq = tmp_path / "fastq"
        state = {}
        results = []
        barrier = threading_module.Barrier(2)

        def _worker():
            barrier.wait()
            results.append(accession_mod._project_download("SRR1", fastq, 1, compress=True))

        with patch("prepare.data.sra.accession.shutil.which", side_effect=lambda tool: None):
            with patch("prepare.data.sra.fastq.shutil.which", side_effect=_no_prefetch_with_pigz):
                with patch(
                    "prepare.utils.security.SecureSubprocess.run_secure", side_effect=_slow_download_tools(state)
                ):
                    threads = [threading_module.Thread(target=_worker) for _ in range(2)]
                    for thread in threads:
                        thread.start()
                    for thread in threads:
                        thread.join(timeout=10)

        assert [c[0] for c in state["calls"]].count("fasterq-dump") == 1
        assert all(success for success, _ in results)
        messages = [message for _, message in results]
        assert messages.count("already exists") == 1
        assert sum(message.startswith("Downloaded 2 files") for message in messages) == 1

    def test_a_watcher_never_sees_an_incomplete_folder(self, tmp_path):
        import threading as threading_module

        fastq = tmp_path / "fastq"
        state = {}
        observations = []
        done = threading_module.Event()

        def _watch():
            while not done.is_set():
                folder = fastq / "SRR1"
                has_fastq = accession_has_fastq(folder)
                try:
                    names = tuple(sorted(p.name for p in folder.iterdir())) if folder.exists() else None
                except OSError:
                    names = None
                observations.append((names, has_fastq))
                time.sleep(0.002)

        watcher = threading_module.Thread(target=_watch)
        watcher.start()
        try:
            success, message = self._run(fastq, state)
        finally:
            done.set()
            watcher.join(timeout=5)

        assert success is True, message
        complete = ("SRR1_1.fastq.gz", "SRR1_2.fastq.gz")
        assert observations
        for names, has_fastq in observations:
            assert names is None or names == complete, names
            # The check ran before the listing, and a published folder never goes away again.
            assert not has_fastq or names == complete

    def test_lock_wait_gives_up_naming_the_holder(self, tmp_path, fast_project_lock_poll):
        fastq = tmp_path / "fastq"
        lock = accession_mod.project_lock_path(fastq, "SRR1")
        lock.parent.mkdir(parents=True)
        lock.write_text(json.dumps({"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00"}))
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")

        started = time.monotonic()
        with patch("prepare.data.sra.accession.download_accession") as mock_download:
            stats = download_sra(fastq, acc, max_retries=0, lock_wait=1.0, num_threads=1, max_workers=1)
        elapsed = time.monotonic() - started

        mock_download.assert_not_called()
        assert 1.0 <= elapsed < 2.0
        assert stats["failed"] == 1
        message = stats["results"]["SRR1"]
        assert "SRR1" in message and "4242" in message and "otherhost" in message
        assert json.loads(lock.read_text())["pid"] == 4242

    def test_a_stop_while_waiting_reports_interrupted(self, tmp_path, fast_project_lock_poll):
        import threading as threading_module

        fastq = tmp_path / "fastq"
        lock = accession_mod.project_lock_path(fastq, "SRR1")
        lock.parent.mkdir(parents=True)
        lock.write_text(json.dumps({"pid": 4242, "host": "otherhost", "started": "2026-01-01T00:00:00+00:00"}))
        stop = threading_module.Event()
        stop.set()
        with patch("prepare.data.sra.accession.download_accession") as mock_download:
            assert accession_mod._project_download("SRR1", fastq, stop=stop) == (False, "interrupted")
        mock_download.assert_not_called()

    def test_download_sra_without_a_store_goes_through_the_project_lock(self, tmp_path):
        acc = tmp_path / "acc.txt"
        acc.write_text("SRR1\n")
        fastq = tmp_path / "fastq"
        seen = {}

        def _download(accession, output_folder, *args, **kwargs):
            seen["lock_held"] = accession_mod.project_lock_path(fastq, accession).exists()
            acc_dir = Path(output_folder) / accession
            acc_dir.mkdir(parents=True)
            (acc_dir / f"{accession}_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
            return True, "Downloaded 1 files, unverified"

        with patch("prepare.data.sra.accession.download_accession", side_effect=_download):
            stats = download_sra(fastq, acc, max_retries=0)
        assert stats["successful"] == 1
        assert seen["lock_held"] is True
        assert accession_has_fastq(fastq / "SRR1")


class TestProjectHiddenFolders:
    """``.locks`` and ``.prepare-tmp`` never count as downloaded accessions."""

    @staticmethod
    def _tree(tmp_path):
        fastq = tmp_path / "fastq"
        (fastq / ".locks").mkdir(parents=True)
        (fastq / ".locks" / "SRR2.lock").write_text("{}")
        staged = fastq / ".prepare-tmp" / "SRR2"
        staged.mkdir(parents=True)
        (staged / "SRR2_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (fastq / ".prepare-tmp" / "SRR2_temp").mkdir()
        (fastq / ".prepare-tmp" / "SRR2_temp" / "SRR2_1.fastq").write_text("@r\nACGT\n")
        return fastq

    def test_scan_downloads_ignores_them(self, tmp_path):
        from prepare.data.registry import scan_downloads

        assert scan_downloads(self._tree(tmp_path)) == {}

    def test_status_inventory_ignores_them(self, tmp_path):
        import argparse

        from prepare.data.registry import Registry
        from prepare.processing.status_report import inventory_report

        fastq = self._tree(tmp_path)
        (tmp_path / "accessions.txt").write_text("SRR2\n.prepare-tmp\n.locks\n")
        args = argparse.Namespace(
            fastq_folder=str(fastq),
            metadata_folder=str(tmp_path / "metadata"),
            genomes_folder=str(tmp_path / "genomes"),
            targeted_folder=str(tmp_path / "targeted"),
            matches_folder=str(tmp_path / "matches"),
            registry=str(tmp_path / "prepare_registry.json"),
            data_root=None,
            accessions_file=str(tmp_path / "accessions.txt"),
            parsed_containment=None,
            stage=None,
            genome=None,
            init=False,
            reconcile=False,
            export_tsv=None,
            next=False,
            list_missing=False,
            json=True,
        )
        report = inventory_report(args, Registry())
        assert report["on_disk"]["fastq_accessions"] == 0
        assert report["wanted"]["fastq_present"] == 0

    def test_staging_folders_are_transient_and_counted(self, tmp_path):
        fastq = self._tree(tmp_path)
        (fastq / "SRR3_temp").mkdir()
        (fastq / "SRR3_temp" / "part").write_bytes(b"x" * 7)
        assert is_transient_folder(".prepare-tmp") is True
        staged_bytes = len("@r\nACGT\n+\nIIII\n") + len("@r\nACGT\n")
        assert transient_bytes(fastq) == staged_bytes + 7


class TestPublishFolder:
    def test_replaces_an_existing_folder(self, tmp_path):
        from prepare.data.sra.cleanup import publish_folder

        tmp = tmp_path / "tmp"
        staged = tmp / "SRR1"
        staged.mkdir(parents=True)
        (staged / "new").write_text("new")
        target = tmp_path / "fastq" / "SRR1"
        target.mkdir(parents=True)
        (target / "old").write_text("old")
        publish_folder(staged, target, tmp)
        assert [p.name for p in target.iterdir()] == ["new"]
        assert list(tmp.iterdir()) == []

    def test_a_failed_rename_restores_the_old_folder(self, tmp_path):
        from prepare.data.sra import cleanup as cleanup_module

        tmp = tmp_path / "tmp"
        tmp.mkdir()
        target = tmp_path / "fastq" / "SRR1"
        target.mkdir(parents=True)
        (target / "old").write_text("old")
        real_replace = os.replace

        def _replace(src, dst):
            if Path(src).name == "SRR1" and Path(src).parent == tmp:
                raise OSError("cross-device")
            return real_replace(src, dst)

        with patch("prepare.data.sra.cleanup.os.replace", side_effect=_replace):
            with pytest.raises(OSError):
                cleanup_module.publish_folder(tmp / "SRR1", target, tmp)
        assert (target / "old").read_text() == "old"

    def test_replaces_a_dangling_symlink_and_a_link_to_a_folder(self, tmp_path):
        from prepare.data.sra.cleanup import publish_folder

        tmp = tmp_path / "tmp"
        for index, link_target in enumerate((tmp_path / "gone", tmp_path / "elsewhere")):
            if index == 1:
                link_target.mkdir()
                (link_target / "keep").write_text("keep")
            staged = tmp / "SRR1"
            staged.mkdir(parents=True)
            (staged / "new").write_text("new")
            target = tmp_path / "fastq" / "SRR1"
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.is_symlink() or target.exists():
                shutil.rmtree(target)
            target.symlink_to(link_target, target_is_directory=True)
            publish_folder(staged, target, tmp)
            assert not target.is_symlink()
            assert [p.name for p in target.iterdir()] == ["new"]
        # The folder a replaced link pointed at is not touched.
        assert (tmp_path / "elsewhere" / "keep").read_text() == "keep"


class TestCachedArchiveWithKeepSra:
    """A redownload always refetches the archive; ``keep_sra`` keeps only a complete or unverified one."""

    def _download(self, output_folder, run=None, state=None, **kwargs):
        fake = run or _fake_tools(state if state is not None else {"reads": 4})
        with patch("prepare.data.sra.accession.shutil.which", side_effect=_with_prefetch):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=fake):
                return download_accession("SRR123", output_folder, compress=False, **kwargs)

    @pytest.mark.parametrize("flag", ["force", "redownload_truncated"])
    def test_a_redownload_wipes_the_cached_archive_even_with_keep_sra(self, tmp_path, flag):
        output_folder = tmp_path / "downloads"
        cache_dir = output_folder / ".sra-cache" / "SRR123"
        cache_dir.mkdir(parents=True)
        (cache_dir / "SRR123.sra").write_bytes(b"short archive")
        fake = _fake_tools({"reads": 4})
        seen = []

        def run(executable, args, **kwargs):
            if executable == "prefetch":
                seen.append((cache_dir / "SRR123.sra").exists())
            return fake(executable, args, **kwargs)

        success, message = self._download(output_folder, run=run, keep_sra=True, **{flag: True})
        assert success is True, message
        assert seen == [False]

    def test_keep_sra_keeps_a_complete_archive(self, tmp_path):
        output_folder = tmp_path / "downloads"
        success, message = self._download(output_folder, keep_sra=True, expected_spots=4)
        assert success is True and "complete" in message, message
        assert (output_folder / ".sra-cache" / "SRR123" / "SRR123.sra").exists()

    def test_keep_sra_keeps_an_unverified_archive(self, tmp_path):
        output_folder = tmp_path / "downloads"
        success, message = self._download(output_folder, keep_sra=True)
        assert success is True and "unverified" in message, message
        assert (output_folder / ".sra-cache" / "SRR123" / "SRR123.sra").exists()

    def test_keep_sra_does_not_keep_a_truncated_archive(self, tmp_path):
        output_folder = tmp_path / "downloads"
        success, message = self._download(output_folder, keep_sra=True, expected_spots=1000)
        assert success is False and "truncated" in message, message
        assert not (output_folder / ".sra-cache" / "SRR123").exists()

    def test_discard_keeps_on_request_only_a_complete_or_unverified_archive(self, tmp_path):
        cases = {
            "Downloaded 2 files, complete (4 of 4 spots)": True,
            "Downloaded 2 files, unverified": True,
            "Downloaded 2 files, truncated (4 of 1000 spots)": False,
        }
        for message, kept in cases.items():
            archive = tmp_path / "SRR1" / "SRR1.sra"
            archive.parent.mkdir(parents=True, exist_ok=True)
            archive.write_bytes(b"x")
            accession_mod._discard_cached_archive(tmp_path, "SRR1", message, keep=True)
            assert archive.exists() is kept, message
            accession_mod._discard_cached_archive(tmp_path, "SRR1", message)
            assert not archive.exists(), message


class TestProjectScratchFolder:
    """Without ``--temp-folder``, fasterq-dump's scratch is ``<fastq>/.prepare-tmp/<ACC>_fqtmp``."""

    def _run(self, fastq, state, **kwargs):
        with patch("prepare.data.sra.accession.shutil.which", side_effect=lambda tool: None):
            with patch("prepare.data.sra.fastq.shutil.which", side_effect=_no_prefetch_with_pigz):
                with patch(
                    "prepare.utils.security.SecureSubprocess.run_secure", side_effect=_scratch_recording_tools(state)
                ):
                    with patch("tempfile.mkdtemp", side_effect=AssertionError("mkdtemp must not be called")):
                        return accession_mod._project_download("SRR1", fastq, 1, compress=False, **kwargs)

    def test_scratch_is_under_the_staging_folder_and_removed_afterwards(self, tmp_path):
        fastq = tmp_path / "fastq"
        scratch = fastq / ".prepare-tmp" / "SRR1_fqtmp"
        scratch.mkdir(parents=True)
        (scratch / "left-over").write_text("from an interrupted run")
        state = {}
        success, message = self._run(fastq, state)
        assert success is True, message
        assert state["temp"] == str(scratch.absolute())
        assert state["temp_existed"] is True
        assert state["left_over_present"] is False
        assert not scratch.exists()
        assert list((fastq / ".prepare-tmp").iterdir()) == []

    def test_scratch_is_removed_after_a_failed_download(self, tmp_path):
        fastq = tmp_path / "fastq"
        state = {"fail": True}
        success, message = self._run(fastq, state)
        assert success is False, message
        assert state["temp"] == str((fastq / ".prepare-tmp" / "SRR1_fqtmp").absolute())
        assert not (fastq / ".prepare-tmp" / "SRR1_fqtmp").exists()

    def test_scratch_is_left_to_the_new_holder_after_a_lost_lock(self, tmp_path):
        """A process that took the lock over uses the same scratch name, so it is not removed."""
        from prepare.utils.lockfile import LockLost

        fastq = tmp_path / "fastq"
        scratch = fastq / ".prepare-tmp" / "SRR1_fqtmp"
        state = {}
        with patch.object(accession_mod, "verify_held", side_effect=LockLost("taken over by another host")):
            success, message = self._run(fastq, state)
        assert success is False
        assert message.startswith("lock lost:")
        assert state["temp"] == str(scratch.absolute())
        assert scratch.is_dir()

    def test_a_given_temp_folder_is_used_and_left_in_place(self, tmp_path):
        fastq = tmp_path / "fastq"
        own = tmp_path / "local-scratch"
        state = {}
        success, message = self._run(fastq, state, temp_folder=own)
        assert success is True, message
        assert state["temp"] == str(own.absolute())
        assert own.is_dir()
        assert not (fastq / ".prepare-tmp" / "SRR1_fqtmp").exists()

    def test_the_scratch_name_is_a_transient_folder(self):
        assert is_transient_folder("SRR1_fqtmp")

    def test_the_space_guard_counts_the_scratch_on_the_project_filesystem(self, tmp_path):
        from prepare.data.sra import space as space_mod

        locations = space_mod.download_locations(tmp_path / "fastq")
        assert locations[space_mod.TEMP] == tmp_path / "fastq" / ".prepare-tmp"
        given = space_mod.download_locations(tmp_path / "fastq", tmp_path / "local")
        assert given[space_mod.TEMP] == tmp_path / "local"


class TestPrefetchMaxSize:
    """prefetch's ``--max-size`` comes from the ``prefetch_max_size`` setting."""

    def _prefetch_args(self, tmp_path):
        state = {"reads": 4}
        with patch("prepare.data.sra.accession.shutil.which", side_effect=_with_prefetch):
            with patch("prepare.data.sra.accession.SecureSubprocess.run_secure", side_effect=_fake_tools(state)):
                success, message = download_accession("SRR123", tmp_path / "downloads", compress=False)
        assert success is True, message
        args = state["calls"][0][1]
        assert state["calls"][0][0] == "prefetch"
        return args[args.index("--max-size") + 1]

    def test_the_default_is_100g(self, tmp_path):
        assert self._prefetch_args(tmp_path) == "100G"

    def test_the_environment_variable_reaches_prefetch(self, tmp_path, monkeypatch):
        monkeypatch.setenv("PREPARE_PREFETCH_MAX_SIZE", "250G")
        assert self._prefetch_args(tmp_path) == "250G"

    def test_the_config_value_reaches_prefetch(self, tmp_path):
        config = tmp_path / "config" / "prepare" / "config.toml"
        config.parent.mkdir(parents=True, exist_ok=True)
        config.write_text('[runtime]\nprefetch_max_size = "500g"\n')
        assert self._prefetch_args(tmp_path) == "500g"

    def test_an_invalid_value_stops_the_download_with_the_settings_error(self, tmp_path, monkeypatch):
        from prepare.core.exceptions import ConfigurationError

        monkeypatch.setenv("PREPARE_PREFETCH_MAX_SIZE", "lots")
        with pytest.raises(ConfigurationError, match="PREPARE_PREFETCH_MAX_SIZE"):
            self._prefetch_args(tmp_path)
