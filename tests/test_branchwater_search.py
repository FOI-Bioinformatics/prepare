"""Tests for prepare.data.branchwater_search (no network; sourmash optional)."""

import json
import random
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import pytest
import requests

from prepare.core.exceptions import ConfigurationError, DataAccessError
from prepare.data import branchwater_search
from prepare.data.branchwater_search import (
    BRANCHWATER_COLUMNS,
    DEFAULT_SERVER,
    load_signature,
    search_index,
    sketch_fasta,
    write_branchwater_csv,
)

SIG_OBJECT = {
    "class": "sourmash_signature",
    "name": "wmel",
    "signatures": [
        {"num": 0, "ksize": 21, "seed": 42, "max_hash": 18446744073709552, "mins": [1, 2, 3], "molecule": "DNA"}
    ],
}


def _write_fasta(path, seed=7, length=5000):
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    path.write_text(f">contig1 test\n{seq[:2500]}\n>contig2\n{seq[2500:]}\n")


class TestSketchFasta:
    def test_builds_k21_scaled1000_signature(self, tmp_path):
        pytest.importorskip("sourmash")
        fasta = tmp_path / "GCF_000000001.1.fna"
        _write_fasta(fasta)
        signature = sketch_fasta(fasta)
        sketch = signature["signatures"][0]
        assert signature["name"] == "GCF_000000001.1"
        assert sketch["ksize"] == 21
        assert round(2**64 / sketch["max_hash"]) == 1000
        assert len(sketch["mins"]) > 0

    def test_missing_fasta_raises(self, tmp_path):
        pytest.importorskip("sourmash")
        with pytest.raises(DataAccessError, match="not found"):
            sketch_fasta(tmp_path / "missing.fna")

    def test_empty_fasta_raises(self, tmp_path):
        pytest.importorskip("sourmash")
        fasta = tmp_path / "empty.fna"
        fasta.write_text("")
        with pytest.raises(DataAccessError, match="No sequences"):
            sketch_fasta(fasta)

    def test_missing_sourmash_gives_install_hint(self, tmp_path, monkeypatch):
        monkeypatch.setitem(sys.modules, "sourmash", None)
        fasta = tmp_path / "x.fna"
        _write_fasta(fasta)
        with pytest.raises(ConfigurationError, match=r"prepare\[sourmash\]"):
            sketch_fasta(fasta)

    def test_missing_sourmash_hint_names_the_interpreter(self, tmp_path, monkeypatch):
        """A user with several Python environments needs to know which one to install into."""
        monkeypatch.setitem(sys.modules, "sourmash", None)
        fasta = tmp_path / "x.fna"
        _write_fasta(fasta)
        with pytest.raises(ConfigurationError) as exc_info:
            sketch_fasta(fasta)
        assert sys.executable in str(exc_info.value)
        assert "-m pip install 'prepare[sourmash]'" in str(exc_info.value)
        assert str(exc_info.value).startswith("Sketching a genome needs the 'sourmash' package.")


class TestLoadSignature:
    def test_list_form(self, tmp_path):
        path = tmp_path / "wmel.sig"
        path.write_text(json.dumps([SIG_OBJECT]))
        assert load_signature(path)["name"] == "wmel"

    def test_object_form(self, tmp_path):
        path = tmp_path / "wmel.sig"
        path.write_text(json.dumps(SIG_OBJECT))
        assert load_signature(path)["signatures"][0]["ksize"] == 21

    def test_wrong_parameters_raise(self, tmp_path):
        bad = json.loads(json.dumps(SIG_OBJECT))
        bad["signatures"][0]["ksize"] = 31
        path = tmp_path / "k31.sig"
        path.write_text(json.dumps([bad]))
        with pytest.raises(DataAccessError, match="k=21, scaled=1000"):
            load_signature(path)

    def test_missing_or_invalid(self, tmp_path):
        with pytest.raises(DataAccessError, match="not found"):
            load_signature(tmp_path / "none.sig")
        path = tmp_path / "bad.sig"
        path.write_text("{not json")
        with pytest.raises(DataAccessError, match="not valid JSON"):
            load_signature(path)

    def test_multi_sketch_selects_k21(self, tmp_path):
        k31_sketch = {
            "num": 0,
            "ksize": 31,
            "seed": 42,
            "max_hash": 18446744073709552,
            "mins": [4, 5, 6],
            "molecule": "DNA",
        }
        k21_sketch = SIG_OBJECT["signatures"][0]
        multi = {**SIG_OBJECT, "signatures": [k31_sketch, k21_sketch]}
        path = tmp_path / "multi.sig"
        path.write_text(json.dumps([multi]))

        signature = load_signature(path)

        assert len(signature["signatures"]) == 1
        assert signature["signatures"][0]["ksize"] == 21


def _mock_session(status_code=200, text=""):
    """Build a MagicMock standing in for a requests.Session, its .post returning a canned response."""
    session = MagicMock()
    response = Mock(status_code=status_code, text=text)
    response.iter_lines.side_effect = lambda *args, **kwargs: iter(text.splitlines())
    session.post.return_value = response
    return session


class TestSessionFactory:
    def test_retry_adapter_is_mounted_on_both_schemes(self):
        session = branchwater_search._session()
        for scheme in ("https://", "http://"):
            adapter = session.adapters[scheme]
            assert adapter.max_retries.total == 4
            assert adapter.max_retries.backoff_factor == 2
            assert set(adapter.max_retries.status_forcelist) == {429, 500, 502, 503, 504}
            assert "POST" in adapter.max_retries.allowed_methods


class TestSearchIndex:
    @patch("prepare.data.branchwater_search._session")
    def test_posts_signature_and_parses_csv(self, mock_session_factory):
        session = _mock_session(text="SRA accession,containment\nSRR1,0.5\nSRR2,0.9\n")
        mock_session_factory.return_value = session
        matches = search_index(SIG_OBJECT, 0.1)
        assert matches == [("SRR2", 0.9), ("SRR1", 0.5)]
        args, kwargs = session.post.call_args
        assert args[0] == f"{DEFAULT_SERVER}/search"
        assert kwargs["json"] == {"threshold": 0.1, "signature": SIG_OBJECT}
        assert kwargs["timeout"] == 600
        assert kwargs["stream"] is True

    @patch("prepare.data.branchwater_search._session")
    def test_a_response_without_a_charset_is_still_parsed(self, mock_session_factory):
        """requests leaves ``encoding`` unset when the response declares no charset, and
        ``iter_lines(decode_unicode=True)`` then yields bytes, which csv.reader rejects with
        an error that is not a PrepareError."""
        text = "SRA accession,containment\nSRR1,0.5\n"
        session = MagicMock()
        response = Mock(status_code=200, text=text, encoding=None)

        def iter_lines(*args, **kwargs):
            # What requests does: decode only when an encoding is known.
            if response.encoding is None:
                return iter([line.encode() for line in text.splitlines()])
            return iter(text.splitlines())

        response.iter_lines.side_effect = iter_lines
        session.post.return_value = response
        mock_session_factory.return_value = session

        assert search_index(SIG_OBJECT, 0.1) == [("SRR1", 0.5)]
        assert response.encoding == "utf-8"

    @patch("prepare.data.branchwater_search._session")
    def test_http_error_raises(self, mock_session_factory):
        mock_session_factory.return_value = _mock_session(status_code=500, text="boom")
        with pytest.raises(DataAccessError, match="HTTP 500"):
            search_index(SIG_OBJECT, 0.1)

    @patch("prepare.data.branchwater_search._session")
    def test_404_raises_without_retry(self, mock_session_factory):
        session = _mock_session(status_code=404, text="not found")
        mock_session_factory.return_value = session
        with pytest.raises(DataAccessError, match="HTTP 404"):
            search_index(SIG_OBJECT, 0.1)
        assert session.post.call_count == 1

    @patch("prepare.data.branchwater_search._session")
    def test_connection_error_raises(self, mock_session_factory):
        session = MagicMock()
        session.post.side_effect = requests.exceptions.ConnectionError("down")
        mock_session_factory.return_value = session
        with pytest.raises(DataAccessError, match="search failed"):
            search_index(SIG_OBJECT, 0.1, server="https://example.org/")

    @patch("prepare.data.branchwater_search._session")
    def test_unexpected_header_raises(self, mock_session_factory):
        mock_session_factory.return_value = _mock_session(text="<html>redirect</html>")
        with pytest.raises(DataAccessError, match="Unexpected Branchwater response"):
            search_index(SIG_OBJECT, 0.1)

    @patch("prepare.data.branchwater_search._session")
    def test_empty_result(self, mock_session_factory):
        mock_session_factory.return_value = _mock_session(text="SRA accession,containment\n")
        assert search_index(SIG_OBJECT, 0.1) == []

    @patch("prepare.data.branchwater_search._session")
    def test_rows_below_threshold_are_dropped(self, mock_session_factory):
        mock_session_factory.return_value = _mock_session(
            text="SRA accession,containment\nSRR1,0.0009\nSRR2,0.5\nSRR3,0.1\n"
        )
        matches = search_index(SIG_OBJECT, 0.1)
        assert matches == [("SRR2", 0.5), ("SRR3", 0.1)]

    @patch("prepare.data.branchwater_search._session")
    def test_all_rows_below_threshold_gives_empty(self, mock_session_factory):
        mock_session_factory.return_value = _mock_session(text="SRA accession,containment\nSRR1,0.0009\nSRR2,0.001\n")
        assert search_index(SIG_OBJECT, 0.1) == []


class TestSearchIndexCache:
    @patch("prepare.data.branchwater_search._session")
    def test_cache_write_then_hit_without_request(self, mock_session_factory, tmp_path):
        mock_session_factory.return_value = _mock_session(text="SRA accession,containment\nSRR1,0.5\n")
        cache_dir = tmp_path / "cache"

        first = search_index(SIG_OBJECT, 0.1, cache_dir=cache_dir)
        assert first == [("SRR1", 0.5)]
        assert mock_session_factory.call_count == 1

        second = search_index(SIG_OBJECT, 0.1, cache_dir=cache_dir)
        assert second == [("SRR1", 0.5)]
        assert mock_session_factory.call_count == 1

    @patch("prepare.data.branchwater_search._session")
    def test_zero_row_response_is_not_cached(self, mock_session_factory, tmp_path, caplog):
        mock_session_factory.return_value = _mock_session(text="SRA accession,containment\n")
        cache_dir = tmp_path / "cache"

        with caplog.at_level("INFO", logger="prepare.data.branchwater_search"):
            assert search_index(SIG_OBJECT, 0.1, cache_dir=cache_dir) == []
            assert search_index(SIG_OBJECT, 0.1, cache_dir=cache_dir) == []

        assert mock_session_factory.call_count == 2
        assert "empty result not cached" in caplog.text

    @patch("prepare.data.branchwater_search._session")
    def test_cache_hit_logs_message(self, mock_session_factory, tmp_path, caplog):
        mock_session_factory.return_value = _mock_session(text="SRA accession,containment\nSRR1,0.5\n")
        cache_dir = tmp_path / "cache"
        search_index(SIG_OBJECT, 0.1, cache_dir=cache_dir)
        with caplog.at_level("INFO"):
            search_index(SIG_OBJECT, 0.1, cache_dir=cache_dir)
        assert "using cached Branchwater result from" in caplog.text

    @patch("prepare.data.branchwater_search._session")
    def test_refresh_bypasses_cache(self, mock_session_factory, tmp_path):
        mock_session_factory.return_value = _mock_session(text="SRA accession,containment\nSRR1,0.5\n")
        cache_dir = tmp_path / "cache"
        search_index(SIG_OBJECT, 0.1, cache_dir=cache_dir)
        search_index(SIG_OBJECT, 0.1, cache_dir=cache_dir, refresh=True)
        assert mock_session_factory.call_count == 2

    @patch("prepare.data.branchwater_search._session")
    def test_stale_cache_by_age_refetches(self, mock_session_factory, tmp_path):
        mock_session_factory.return_value = _mock_session(text="SRA accession,containment\nSRR1,0.5\n")
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()

        key = branchwater_search._cache_key(SIG_OBJECT, 0.1, DEFAULT_SERVER)
        (cache_dir / f"{key}.csv").write_text("SRA accession,containment\nSRR1,0.5\n")
        old_fetched = (datetime.now(timezone.utc) - timedelta(days=10)).isoformat()
        (cache_dir / f"{key}.json").write_text(
            json.dumps({"fetched": old_fetched, "server": DEFAULT_SERVER, "threshold": 0.1, "rows": 1})
        )

        search_index(SIG_OBJECT, 0.1, cache_dir=cache_dir, max_cache_age_days=1)
        assert mock_session_factory.call_count == 1

    @patch("prepare.data.branchwater_search._session")
    def test_fresh_cache_within_max_age_is_used(self, mock_session_factory, tmp_path):
        mock_session_factory.return_value = _mock_session(text="SRA accession,containment\nSRR1,0.5\n")
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()

        key = branchwater_search._cache_key(SIG_OBJECT, 0.1, DEFAULT_SERVER)
        (cache_dir / f"{key}.csv").write_text("SRA accession,containment\nSRR1,0.5\n")
        recent_fetched = (datetime.now(timezone.utc) - timedelta(hours=1)).isoformat()
        (cache_dir / f"{key}.json").write_text(
            json.dumps({"fetched": recent_fetched, "server": DEFAULT_SERVER, "threshold": 0.1, "rows": 1})
        )

        matches = search_index(SIG_OBJECT, 0.1, cache_dir=cache_dir, max_cache_age_days=1)
        assert matches == [("SRR1", 0.5)]
        assert mock_session_factory.call_count == 0

    @patch("prepare.data.branchwater_search._session")
    def test_cache_write_failure_is_logged_and_ignored(self, mock_session_factory, tmp_path, caplog, monkeypatch):
        mock_session_factory.return_value = _mock_session(text="SRA accession,containment\nSRR1,0.5\n")
        cache_dir = tmp_path / "cache"

        def fail_mkdir(*args, **kwargs):
            raise OSError("no space left on device")

        monkeypatch.setattr(Path, "mkdir", fail_mkdir)
        with caplog.at_level("WARNING"):
            result = search_index(SIG_OBJECT, 0.1, cache_dir=cache_dir)
        assert result == [("SRR1", 0.5)]
        assert "cache" in caplog.text.lower()


class TestWriteBranchwaterCsv:
    def test_header_and_cani(self, tmp_path):
        out = write_branchwater_csv([("SRR11011981", 0.9461), ("SRR2", 0.0)], tmp_path / "bw" / "wmel.csv")
        lines = out.read_text().splitlines()
        assert lines[0] == ",".join(BRANCHWATER_COLUMNS)
        assert lines[1] == "SRR11011981,0.9461,0.9974,,,,,,,"
        assert lines[2] == "SRR2,0.0000,0.0000,,,,,,,"

    def test_empty_matches_write_header_only(self, tmp_path):
        out = write_branchwater_csv([], tmp_path / "empty.csv")
        assert out.read_text().splitlines() == [",".join(BRANCHWATER_COLUMNS)]
