"""
Tests for prepare.data.metadata module.
"""

import os

import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock, call
from urllib.error import HTTPError, URLError
import xml.etree.ElementTree as ET

from prepare.core.exceptions import DataAccessError, ValidationError
from prepare.data.metadata import (
    _get_unique_accessions,
    _download_single_metadata,
    download_metadata,
    _download_accessions_metadata,
    _download_batch_metadata,
    _pace_requests,
    _split_efetch_packages,
    _extract_metadata_fields,
    _extract_sample_attributes,
    parse_metadata,
    parse_metadata_xml,
    check_metadata_attributes,
)


def _http_error(code, reason="Error"):
    """Build a urllib HTTPError with the given status code, as Entrez.efetch would raise."""
    return HTTPError("https://eutils.ncbi.nlm.nih.gov/", code, reason, hdrs=None, fp=None)


def _package_xml(accessions_and_attrs):
    """Build a minimal efetch response with one EXPERIMENT_PACKAGE per accession.

    ``accessions_and_attrs`` is a list of (accession, {attr: value}) pairs.
    """
    root = ET.Element("EXPERIMENT_PACKAGE_SET")
    for accession, attrs in accessions_and_attrs:
        package = ET.SubElement(root, "EXPERIMENT_PACKAGE")
        run_set = ET.SubElement(package, "RUN_SET")
        run = ET.SubElement(run_set, "RUN", accession=accession, **attrs)
        identifiers = ET.SubElement(run, "IDENTIFIERS")
        ET.SubElement(identifiers, "PRIMARY_ID").text = accession
    return ET.tostring(root, encoding="unicode")


class TestGetUniqueAccessions:
    """Test _get_unique_accessions function."""

    def test_get_unique_accessions_branchwater_format(self, tmp_path):
        """Test extracting unique accessions from branchwater format."""
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()

        # Create test CSV with branchwater format
        csv_content = "acc,containment,organism\nSRR123,0.85,E. coli\nSRR456,0.92,S. aureus\nSRR123,0.75,E. coli"
        (matches_dir / "genome1.csv").write_text(csv_content)

        result = _get_unique_accessions(matches_dir, threshold=0.8)
        assert result == {"SRR123", "SRR456"}

    def test_get_unique_accessions_sra_format(self, tmp_path):
        """Test extracting unique accessions from SRA format."""
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()

        # Create test CSV with SRA format
        csv_content = "SRA accession,containment,species\nSRR789,0.88,B. subtilis\nSRR101,0.95,P. aeruginosa"
        (matches_dir / "genome2.csv").write_text(csv_content)

        result = _get_unique_accessions(matches_dir, threshold=0.9)
        assert result == {"SRR101"}

    def test_get_unique_accessions_threshold_filtering(self, tmp_path):
        """Test threshold filtering works correctly."""
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()

        csv_content = "acc,containment,organism\nSRR123,0.85,E. coli\nSRR456,0.75,S. aureus\nSRR789,0.95,B. subtilis"
        (matches_dir / "genome1.csv").write_text(csv_content)

        result = _get_unique_accessions(matches_dir, threshold=0.8)
        assert result == {"SRR123", "SRR789"}

    def test_get_unique_accessions_no_csv_files(self, tmp_path):
        """Test handling when no CSV files are found."""
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()

        with pytest.raises(DataAccessError, match="No match CSV in"):
            _get_unique_accessions(matches_dir, threshold=0.0)

    def test_get_unique_accessions_unknown_format(self, tmp_path):
        """Test handling unknown file format."""
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()

        # Create CSV with unknown format
        csv_content = "unknown_col,other_col\nvalue1,value2"
        (matches_dir / "unknown.csv").write_text(csv_content)

        with pytest.raises(DataAccessError, match="no 'acc' or 'SRA accession' column"):
            _get_unique_accessions(matches_dir, threshold=0.0)

    def test_get_unique_accessions_read_error(self, tmp_path):
        """Test handling CSV read errors."""
        matches_dir = tmp_path / "matches"
        matches_dir.mkdir()

        # Create malformed CSV
        (matches_dir / "malformed.csv").write_text("invalid,csv\ncontent")

        with patch("pandas.read_csv", side_effect=pd.errors.ParserError("Read error")):
            with pytest.raises(DataAccessError, match="malformed.csv: Read error"):
                _get_unique_accessions(matches_dir, threshold=0.0)


class TestDownloadSingleMetadata:
    """Test _download_single_metadata function."""

    def test_download_single_metadata_success(self, tmp_path):
        """Test successful metadata download."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()

        mock_response = _package_xml([("SRR123", {})]).encode()

        with patch("prepare.data.metadata.Entrez.efetch") as mock_efetch:
            mock_handle = MagicMock()
            mock_handle.read.return_value = mock_response
            mock_efetch.return_value = mock_handle

            with patch("prepare.data.metadata_fetch._pace_requests"):
                success, result = _download_single_metadata("SRR123", metadata_path, "test@example.com")

        assert success is True
        assert isinstance(result, Path)
        assert result.name == "SRR123_metadata.xml"
        assert ET.fromstring(result.read_text()).find(".//RUN").get("accession") == "SRR123"

    def test_download_single_metadata_404_is_final(self, tmp_path):
        """A 404 on a single accession is a final failure, never retried."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()

        with patch("prepare.data.metadata.Entrez.efetch", side_effect=_http_error(404)) as mock_efetch:
            with patch("prepare.data.metadata_fetch._pace_requests"):
                with patch("time.sleep") as mock_sleep:
                    success, result = _download_single_metadata("SRR123", metadata_path, "test@example.com")

        assert success is False
        assert result == "HTTP 404: not found at NCBI"
        mock_efetch.assert_called_once()
        mock_sleep.assert_not_called()

    def test_download_single_metadata_503_retried_then_fails(self, tmp_path):
        """A 503 is retried up to 3 times with exponential backoff, then fails."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()

        with patch("prepare.data.metadata.Entrez.efetch", side_effect=_http_error(503)):
            with patch("prepare.data.metadata_fetch._pace_requests"):
                with patch("time.sleep") as mock_sleep:
                    success, result = _download_single_metadata("SRR123", metadata_path, "test@example.com")

        assert success is False
        assert result == "network: HTTP 503 after 3 attempts"
        assert mock_sleep.call_args_list == [call(2), call(4), call(8)]

    def test_download_single_metadata_partial_retry(self, tmp_path):
        """Test successful download after one retry."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()

        mock_response = _package_xml([("SRR123", {})]).encode()

        with patch("prepare.data.metadata.Entrez.efetch") as mock_efetch:
            # First call fails, second succeeds
            mock_handle = MagicMock()
            mock_handle.read.return_value = mock_response
            mock_efetch.side_effect = [URLError("Network error"), mock_handle]

            with patch("prepare.data.metadata_fetch._pace_requests"):
                with patch("time.sleep"):
                    success, result = _download_single_metadata("SRR123", metadata_path, "test@example.com")

        assert success is True
        assert isinstance(result, Path)

    def test_download_single_metadata_api_key_set_on_entrez(self, tmp_path):
        """Entrez.api_key equals the given key during the call; None leaves it None."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()
        mock_response = b"<?xml version='1.0'?><root>test metadata</root>"

        seen_keys = []

        def fake_efetch(**kwargs):
            from prepare.data.metadata import Entrez

            seen_keys.append(Entrez.api_key)
            handle = MagicMock()
            handle.read.return_value = mock_response
            return handle

        with patch("prepare.data.metadata.Entrez.efetch", side_effect=fake_efetch):
            with patch("prepare.data.metadata_fetch._pace_requests"):
                _download_single_metadata("SRR123", metadata_path, "test@example.com", api_key="my-key")
                _download_single_metadata("SRR124", metadata_path, "test@example.com", api_key=None)

        assert seen_keys == ["my-key", None]


class TestDownloadMetadata:
    """Test download_metadata function."""

    def test_download_metadata_success(self, tmp_path):
        """Test successful metadata download workflow."""
        matches_dir = tmp_path / "matches"
        metadata_dir = tmp_path / "metadata"
        matches_dir.mkdir()

        # Create test match file
        csv_content = "acc,containment,organism\nSRR123,0.95,E. coli"
        (matches_dir / "genome1.csv").write_text(csv_content)

        with patch("prepare.data.metadata._download_accessions_metadata") as mock_download:
            mock_download.return_value = {"SRR123": metadata_dir / "SRR123_metadata.xml"}

            result = download_metadata("test@example.com", matches_dir, metadata_dir, threshold=0.9)

        assert isinstance(result, dict)
        mock_download.assert_called_once()

    def test_download_metadata_dry_run(self, tmp_path):
        """Test dry run mode."""
        matches_dir = tmp_path / "matches"
        metadata_dir = tmp_path / "metadata"
        matches_dir.mkdir()

        csv_content = "acc,containment,organism\nSRR123,0.95,E. coli"
        (matches_dir / "genome1.csv").write_text(csv_content)

        with patch("prepare.data.metadata.logger") as mock_logger:
            result = download_metadata("test@example.com", matches_dir, metadata_dir, dry_run=True)

        assert result == {}
        mock_logger.info.assert_called_with("Dry run, not downloading metadata")

    def test_download_metadata_dry_run_creates_no_folder(self, tmp_path):
        matches_dir = tmp_path / "matches"
        metadata_dir = tmp_path / "metadata"
        matches_dir.mkdir()
        (matches_dir / "genome1.csv").write_text("acc,containment,organism\nSRR123,0.95,E. coli")

        result = download_metadata("test@example.com", matches_dir, metadata_dir, dry_run=True)

        assert result == {}
        assert not metadata_dir.exists()

    def test_download_metadata_skip_existing(self, tmp_path):
        """Test skipping existing metadata files."""
        matches_dir = tmp_path / "matches"
        metadata_dir = tmp_path / "metadata"
        matches_dir.mkdir()
        metadata_dir.mkdir()

        # Create existing metadata file
        (metadata_dir / "SRR123_metadata.xml").write_text("existing")

        csv_content = "acc,containment,organism\nSRR123,0.95,E. coli\nSRR456,0.85,S. aureus"
        (matches_dir / "genome1.csv").write_text(csv_content)

        with patch("prepare.data.metadata._download_accessions_metadata") as mock_download:
            mock_download.return_value = {"SRR456": metadata_dir / "SRR456_metadata.xml"}

            download_metadata("test@example.com", matches_dir, metadata_dir)

        # Should only download SRR456, not SRR123
        download_args = mock_download.call_args[0]
        accessions_to_download = download_args[0]
        assert "SRR456" in accessions_to_download
        assert "SRR123" not in accessions_to_download

    def test_download_metadata_validation_error(self, tmp_path):
        """Test handling validation errors."""
        with patch("prepare.data.metadata.validate_folder", side_effect=ValidationError("Invalid folder")):
            with pytest.raises(DataAccessError):
                download_metadata("test@example.com", "invalid", tmp_path)

    def test_download_metadata_accessions_file_bypasses_matches_folder(self, tmp_path):
        """With --accessions-file the wanted set comes from the file; the matches folder is not read."""
        metadata_dir = tmp_path / "metadata"
        accessions_file = tmp_path / "accessions.txt"
        accessions_file.write_text("SRR1\nSRR2\n# comment\n\nSRR3\n")

        calls = []

        def fake_batch(batch, metadata_path, email, api_key):
            calls.extend(batch)
            successes = {}
            for accession in batch:
                xml_path = metadata_path / f"{accession}_metadata.xml"
                xml_path.write_text("<root/>")
                successes[accession] = xml_path
            return successes, {}

        with patch("prepare.data.metadata.validate_folder") as mock_validate_folder:
            with patch("prepare.data.metadata._download_batch_metadata", side_effect=fake_batch):
                result = download_metadata(
                    "test@example.com",
                    "does-not-exist",
                    metadata_dir,
                    accessions_file=accessions_file,
                )

        mock_validate_folder.assert_not_called()
        assert sorted(calls) == ["SRR1", "SRR2", "SRR3"]
        assert set(result) == {"SRR1", "SRR2", "SRR3"}

    def test_download_metadata_batch_size_out_of_range_raises(self, tmp_path):
        """batch_size outside 1-500 raises ValueError, not DataAccessError."""
        with pytest.raises(ValueError):
            download_metadata("test@example.com", tmp_path, tmp_path, batch_size=0)
        with pytest.raises(ValueError):
            download_metadata("test@example.com", tmp_path, tmp_path, batch_size=501)

    def test_download_metadata_passes_api_key_and_batch_size(self, tmp_path):
        """api_key and batch_size reach _download_accessions_metadata."""
        matches_dir = tmp_path / "matches"
        metadata_dir = tmp_path / "metadata"
        matches_dir.mkdir()
        (matches_dir / "genome1.csv").write_text("acc,containment,organism\nSRR123,0.95,E. coli")

        with patch("prepare.data.metadata._download_accessions_metadata") as mock_download:
            mock_download.return_value = {}
            download_metadata("test@example.com", matches_dir, metadata_dir, api_key="my-key", batch_size=5)

        assert mock_download.call_args.kwargs["api_key"] == "my-key"
        assert mock_download.call_args.kwargs["batch_size"] == 5


class TestDownloadAccessionsMetadata:
    """Test _download_accessions_metadata function."""

    def test_download_accessions_metadata_success(self, tmp_path):
        """Test successful accessions download."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()

        accessions = ["SRR123", "SRR456"]

        with patch("prepare.data.metadata._download_batch_metadata") as mock_download:
            mock_download.return_value = (
                {
                    "SRR123": metadata_path / "SRR123_metadata.xml",
                    "SRR456": metadata_path / "SRR456_metadata.xml",
                },
                {},
            )

            result = _download_accessions_metadata(accessions, metadata_path, "test@example.com", 2)

        assert len(result) == 2
        assert "SRR123" in result
        assert "SRR456" in result
        mock_download.assert_called_once_with(accessions, metadata_path, "test@example.com", None)

    def test_download_accessions_metadata_partial_failure(self, tmp_path):
        """Test handling partial download failures."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()

        accessions = ["SRR123", "SRR456"]

        with patch("prepare.data.metadata._download_batch_metadata") as mock_download:
            mock_download.return_value = ({"SRR123": metadata_path / "SRR123_metadata.xml"}, {"SRR456": "not found"})

            with patch("prepare.data.metadata.logger") as mock_logger:
                result = _download_accessions_metadata(accessions, metadata_path, "test@example.com", 2)

        assert len(result) == 1
        assert "SRR123" in result
        assert "SRR456" not in result
        mock_logger.error.assert_called()

    def test_download_accessions_metadata_passes_email_and_api_key_to_each_batch(self, tmp_path):
        """email and api_key reach every batch call, so each batch sets its own Entrez credentials.

        _download_accessions_metadata no longer sets Entrez.email/api_key itself (that global
        state is set only for the duration of one efetch call, under _ENTREZ_LOCK, so concurrent
        downloads with different credentials do not race on it; see _entrez_credentials). Instead
        it must pass email/api_key through to _download_batch_metadata for every batch.
        """
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()

        with patch("prepare.data.metadata._download_batch_metadata", return_value=({}, {})) as mock_batch:
            _download_accessions_metadata(["SRR1"], metadata_path, "test@example.com", 1, api_key="my-key")
            assert mock_batch.call_args == call(["SRR1"], metadata_path, "test@example.com", "my-key")

            _download_accessions_metadata(["SRR1"], metadata_path, "test@example.com", 1, api_key=None)
            assert mock_batch.call_args == call(["SRR1"], metadata_path, "test@example.com", None)

    def test_download_accessions_metadata_batches_by_batch_size(self, tmp_path):
        """5 accessions with batch_size=2 make 3 efetch calls (batches of 2, 2, 1)."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()
        accessions = [f"SRR{i}" for i in range(5)]

        with patch("prepare.data.metadata._download_batch_metadata", return_value=({}, {})) as mock_download:
            _download_accessions_metadata(accessions, metadata_path, "test@example.com", 5, batch_size=2)

        assert mock_download.call_count == 3
        called_batches = [c.args[0] for c in mock_download.call_args_list]
        assert called_batches == [accessions[0:2], accessions[2:4], accessions[4:5]]


class TestDownloadBatchMetadata:
    """Test _download_batch_metadata: one efetch call per batch, split into per-accession files."""

    def test_batch_of_three_with_shared_experiment_package(self, tmp_path):
        """Three accessions in one response, two sharing a package with two RUNs: three files."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()

        root = ET.Element("EXPERIMENT_PACKAGE_SET")
        shared_package = ET.SubElement(root, "EXPERIMENT_PACKAGE")
        run_set = ET.SubElement(shared_package, "RUN_SET")
        for accession, spots, size in (("SRR1", "100", "1000"), ("SRR2", "150", "1500")):
            run = ET.SubElement(run_set, "RUN", accession=accession, total_spots=spots, size=size)
            ET.SubElement(ET.SubElement(run, "IDENTIFIERS"), "PRIMARY_ID").text = accession

        solo_package = ET.SubElement(root, "EXPERIMENT_PACKAGE")
        solo_run_set = ET.SubElement(solo_package, "RUN_SET")
        solo_run = ET.SubElement(solo_run_set, "RUN", accession="SRR3", total_spots="300", size="3000")
        ET.SubElement(ET.SubElement(solo_run, "IDENTIFIERS"), "PRIMARY_ID").text = "SRR3"

        xml_text = ET.tostring(root, encoding="unicode")
        mock_handle = MagicMock()
        mock_handle.read.return_value = xml_text.encode()

        with patch("prepare.data.metadata.Entrez.efetch", return_value=mock_handle) as mock_efetch:
            with patch("prepare.data.metadata_fetch._pace_requests"):
                successes, failures = _download_batch_metadata(
                    ["SRR1", "SRR2", "SRR3"], metadata_path, "test@example.com", None
                )

        assert mock_efetch.call_count == 1
        assert failures == {}
        assert set(successes) == {"SRR1", "SRR2", "SRR3"}
        for accession, spots, size in (("SRR1", "100", "1000"), ("SRR2", "150", "1500"), ("SRR3", "300", "3000")):
            parsed = parse_metadata_xml(successes[accession])
            assert parsed["Run_ID"] == accession
            assert parsed["Run_Total_Spots"] == spots
            assert parsed["Run_Size"] == size

    def test_accession_missing_from_response_is_reported_failed(self, tmp_path):
        """An accession NCBI doesn't return is a failure, and efetch was called only once."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()

        xml_text = _package_xml([("SRR1", {"total_spots": "10"}), ("SRR2", {"total_spots": "20"})])
        mock_handle = MagicMock()
        mock_handle.read.return_value = xml_text.encode()

        with patch("prepare.data.metadata.Entrez.efetch", return_value=mock_handle) as mock_efetch:
            with patch("prepare.data.metadata_fetch._pace_requests"):
                successes, failures = _download_batch_metadata(
                    ["SRR1", "SRR2", "SRR4"], metadata_path, "test@example.com", None
                )

        assert mock_efetch.call_count == 1
        assert set(successes) == {"SRR1", "SRR2"}
        assert failures == {"SRR4": "not in the NCBI response"}

    def test_404_on_batch_falls_back_to_single_calls(self, tmp_path):
        """A 404 on a batch of more than one accession falls back to per-accession fetches."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()
        batch = ["SRR1", "SRR2", "SRR3"]

        def fake_efetch(**kwargs):
            if "," in kwargs["id"]:
                raise _http_error(404)
            accession = kwargs["id"]
            handle = MagicMock()
            handle.read.return_value = _package_xml([(accession, {})]).encode()
            return handle

        with patch("prepare.data.metadata.Entrez.efetch", side_effect=fake_efetch) as mock_efetch:
            with patch("prepare.data.metadata.Entrez.email", "test@example.com"):
                with patch("prepare.data.metadata_fetch._pace_requests"):
                    successes, failures = _download_batch_metadata(batch, metadata_path, "test@example.com", None)

        assert failures == {}
        assert set(successes) == set(batch)
        # One batch call plus one call per accession in the fallback.
        assert mock_efetch.call_count == 1 + len(batch)

    def test_404_on_single_accession_batch_is_final(self, tmp_path):
        """A 404 on a batch of exactly one accession is a final failure, not retried."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()

        with patch("prepare.data.metadata.Entrez.efetch", side_effect=_http_error(404)) as mock_efetch:
            with patch("prepare.data.metadata_fetch._pace_requests"):
                with patch("time.sleep") as mock_sleep:
                    successes, failures = _download_batch_metadata(["SRR1"], metadata_path, "test@example.com", None)

        assert successes == {}
        assert failures == {"SRR1": "HTTP 404: not found at NCBI"}
        mock_efetch.assert_called_once()
        mock_sleep.assert_not_called()

    def test_503_retried_three_times_then_fails(self, tmp_path):
        """A 503 is retried with sleeps 2, 4, 8 and then fails the whole batch."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()
        batch = ["SRR1", "SRR2"]

        with patch("prepare.data.metadata.Entrez.efetch", side_effect=_http_error(503)):
            with patch("prepare.data.metadata_fetch._pace_requests"):
                with patch("time.sleep") as mock_sleep:
                    successes, failures = _download_batch_metadata(batch, metadata_path, "test@example.com", None)

        assert successes == {}
        assert failures == {accession: "network: HTTP 503 after 3 attempts" for accession in batch}
        assert mock_sleep.call_args_list == [call(2), call(4), call(8)]

    def test_success_path_sleeps_only_through_pace_requests(self, tmp_path):
        """The only sleep on the success path is the one _pace_requests issues."""
        metadata_path = tmp_path / "metadata"
        metadata_path.mkdir()

        xml_text = _package_xml([("SRR1", {"total_spots": "10"})])
        mock_handle = MagicMock()
        mock_handle.read.return_value = xml_text.encode()

        with patch("prepare.data.metadata.Entrez.efetch", return_value=mock_handle):
            with patch("prepare.data.metadata_fetch._pace_requests") as mock_pace:
                with patch("time.sleep") as mock_sleep:
                    successes, failures = _download_batch_metadata(
                        ["SRR1"], metadata_path, "test@example.com", "my-key"
                    )

        assert failures == {}
        mock_pace.assert_called_once_with("my-key")
        mock_sleep.assert_not_called()


class TestSplitEfetchPackages:
    """Test _split_efetch_packages directly."""

    def test_splits_shared_package_into_independent_documents(self):
        xml_text = _package_xml([("SRR1", {"total_spots": "10"}), ("SRR2", {"total_spots": "20"})])
        result = _split_efetch_packages(xml_text, {"SRR1", "SRR2"})

        assert set(result) == {"SRR1", "SRR2"}
        for accession in ("SRR1", "SRR2"):
            root = ET.fromstring(result[accession])
            assert root.tag == "EXPERIMENT_PACKAGE_SET"
            runs = root.findall(".//RUN")
            assert len(runs) == 1
            assert runs[0].get("accession") == accession

    def test_ignores_accessions_not_in_wanted_set(self):
        xml_text = _package_xml([("SRR1", {}), ("SRR2", {})])
        result = _split_efetch_packages(xml_text, {"SRR1"})
        assert set(result) == {"SRR1"}

    def test_missing_accession_absent_from_result(self):
        xml_text = _package_xml([("SRR1", {})])
        result = _split_efetch_packages(xml_text, {"SRR1", "SRR9"})
        assert set(result) == {"SRR1"}

    def test_includes_xml_declaration(self):
        xml_text = _package_xml([("SRR1", {})])
        result = _split_efetch_packages(xml_text, {"SRR1"})
        assert result["SRR1"].startswith("<?xml")


class TestPaceRequests:
    """Test _pace_requests: one RequestPacer per rate, shared across batches in one process."""

    def test_reuses_one_pacer_per_rate(self, monkeypatch):
        import prepare.data.metadata_fetch as metadata_module

        monkeypatch.setattr(metadata_module, "_PACERS", {})
        with patch("prepare.data.ncbi_pacing.RequestPacer.wait", autospec=True) as mock_wait:
            _pace_requests(None)
            _pace_requests(None)
            _pace_requests("my-key")

        pacers = [call_.args[0] for call_ in mock_wait.call_args_list]
        assert pacers[0] is pacers[1]
        assert pacers[0].delay == pytest.approx(1 / 3)
        assert pacers[2].delay == pytest.approx(1 / 10)

    def test_successive_requests_are_spaced_by_the_delay(self, monkeypatch):
        import prepare.data.metadata_fetch as metadata_module

        monkeypatch.setattr(metadata_module, "_PACERS", {})
        clock = {"now": 100.0}

        def fake_sleep(seconds):
            clock["now"] += seconds

        with patch("prepare.data.ncbi_pacing.time.monotonic", side_effect=lambda: clock["now"]):
            with patch("prepare.data.ncbi_pacing.time.sleep", side_effect=fake_sleep) as mock_sleep:
                _pace_requests(None)
                _pace_requests(None)

        mock_sleep.assert_called_once()
        assert mock_sleep.call_args.args[0] == pytest.approx(1 / 3)

    def test_concurrent_first_calls_share_one_pacer(self, monkeypatch):
        """Two threads racing on the first call for a key build one pacer, not two."""
        import threading
        import time as time_module

        import prepare.data.metadata_fetch as metadata_module

        monkeypatch.setattr(metadata_module, "_PACERS", {})
        built = []
        waited = []

        class SlowPacer:
            def __init__(self, api_key):
                built.append(self)
                time_module.sleep(0.05)  # widen the check-then-set window

            def wait(self):
                waited.append(self)

        monkeypatch.setattr(metadata_module, "RequestPacer", SlowPacer)
        barrier = threading.Barrier(2)

        def call():
            barrier.wait()
            _pace_requests(None)

        threads = [threading.Thread(target=call) for _ in range(2)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        assert len(built) == 1
        assert len(waited) == 2 and waited[0] is waited[1] is built[0]


class TestExtractMetadataFields:
    """Test _extract_metadata_fields function."""

    def test_extract_metadata_fields_success(self):
        """Test successful metadata field extraction."""
        xml_content = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE_SET>
            <EXPERIMENT_PACKAGE>
                <EXPERIMENT>
                    <IDENTIFIERS>
                        <PRIMARY_ID>EXP123</PRIMARY_ID>
                    </IDENTIFIERS>
                    <TITLE>Test experiment</TITLE>
                    <LIBRARY_DESCRIPTOR>
                        <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
                        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
                    </LIBRARY_DESCRIPTOR>
                </EXPERIMENT>
                <SAMPLE>
                    <IDENTIFIERS>
                        <PRIMARY_ID>SAMN123</PRIMARY_ID>
                    </IDENTIFIERS>
                    <SAMPLE_NAME>
                        <SCIENTIFIC_NAME>Escherichia coli</SCIENTIFIC_NAME>
                    </SAMPLE_NAME>
                </SAMPLE>
                <RUN_SET>
                    <RUN>
                        <IDENTIFIERS>
                            <PRIMARY_ID>SRR123</PRIMARY_ID>
                        </IDENTIFIERS>
                    </RUN>
                </RUN_SET>
            </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""

        tree = ET.fromstring(xml_content)

        result = _extract_metadata_fields(tree, "test.xml")

        assert result["Run_ID"] == "SRR123"
        assert result["Sample_ID"] == "SAMN123"
        assert result["Experiment_Title"] == "Test experiment"
        assert result["Experiment_ID"] == "EXP123"
        assert result["Sample_Scientific_Name"] == "Escherichia coli"
        assert result["Experiment_Library_Strategy"] == "WGS"

    def test_extract_metadata_fields_missing_elements(self):
        """Test handling missing XML elements."""
        xml_content = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE_SET>
            <EXPERIMENT_PACKAGE>
            </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""

        tree = ET.fromstring(xml_content)

        result = _extract_metadata_fields(tree, "test.xml")

        # Should have None values for missing elements
        assert result["Run_ID"] is None
        assert result["Sample_ID"] is None
        assert result["Experiment_Title"] is None

    def test_extract_metadata_fields_complex_structure(self):
        """Test extraction from complex XML structure."""
        xml_content = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE_SET>
            <EXPERIMENT_PACKAGE>
                <EXPERIMENT>
                    <IDENTIFIERS>
                        <PRIMARY_ID>EXP456</PRIMARY_ID>
                    </IDENTIFIERS>
                    <TITLE>Complex experiment</TITLE>
                    <LIBRARY_DESCRIPTOR>
                        <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
                        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
                        <LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION>
                    </LIBRARY_DESCRIPTOR>
                </EXPERIMENT>
                <SAMPLE>
                    <IDENTIFIERS>
                        <PRIMARY_ID>SAMN456</PRIMARY_ID>
                    </IDENTIFIERS>
                    <TITLE>Test sample</TITLE>
                </SAMPLE>
                <RUN_SET>
                    <RUN>
                        <IDENTIFIERS>
                            <PRIMARY_ID>SRR456</PRIMARY_ID>
                        </IDENTIFIERS>
                        <Total_spots>1000000</Total_spots>
                        <Total_bases>150000000</Total_bases>
                    </RUN>
                </RUN_SET>
            </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""

        tree = ET.fromstring(xml_content)

        result = _extract_metadata_fields(tree, "test.xml")

        assert result["Run_ID"] == "SRR456"
        assert result["Sample_ID"] == "SAMN456"
        assert result["Experiment_Library_Strategy"] == "WGS"
        assert result["Experiment_Library_Source"] == "GENOMIC"
        assert result["Experiment_Library_Selection"] == "RANDOM"
        assert result["Sample_Title"] == "Test sample"
        assert result["Run_Total_Spots"] == "1000000"
        assert result["Run_Total_Bases"] == "150000000"

    def test_extract_metadata_fields_ncbi_attributes(self):
        """Real NCBI efetch XML holds spots, bases, size and md5 as attributes, not child text."""
        xml_content = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE_SET>
            <EXPERIMENT_PACKAGE>
                <EXPERIMENT>
                    <IDENTIFIERS>
                        <PRIMARY_ID>EXP1</PRIMARY_ID>
                    </IDENTIFIERS>
                    <LIBRARY_DESCRIPTOR>
                        <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
                        <LIBRARY_LAYOUT>
                            <PAIRED/>
                        </LIBRARY_LAYOUT>
                    </LIBRARY_DESCRIPTOR>
                    <PLATFORM>
                        <ILLUMINA>
                            <INSTRUMENT_MODEL>Illumina HiSeq 2500</INSTRUMENT_MODEL>
                        </ILLUMINA>
                    </PLATFORM>
                </EXPERIMENT>
                <RUN_SET>
                    <RUN accession="SRR1" total_spots="47964651" total_bases="14389395300" size="4744553813">
                        <IDENTIFIERS>
                            <PRIMARY_ID>SRR1</PRIMARY_ID>
                        </IDENTIFIERS>
                        <SRAFiles>
                            <SRAFile filename="SRR1" md5="abc" semantic_name="run"/>
                        </SRAFiles>
                    </RUN>
                </RUN_SET>
            </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""

        tree = ET.fromstring(xml_content)

        result = _extract_metadata_fields(tree, "test.xml")

        assert result["Run_Total_Spots"] == "47964651"
        assert result["Run_Total_Bases"] == "14389395300"
        assert result["Run_Size"] == "4744553813"
        assert result["Run_MD5"] == "abc"
        assert result["Run_Filename"] == "SRR1"
        assert result["Experiment_Library_Layout"] == "PAIRED"
        assert result["Platform"] == "ILLUMINA"

    def test_extract_metadata_fields_srafile_prefers_run_semantic_name(self):
        """When multiple SRAFile entries exist, prefer the one marked semantic_name=run."""
        xml_content = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE_SET>
            <EXPERIMENT_PACKAGE>
                <RUN_SET>
                    <RUN accession="SRR2">
                        <SRAFiles>
                            <SRAFile filename="other" md5="wrong" semantic_name="other"/>
                            <SRAFile filename="SRR2" md5="right" semantic_name="run"/>
                        </SRAFiles>
                    </RUN>
                </RUN_SET>
            </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""

        tree = ET.fromstring(xml_content)

        result = _extract_metadata_fields(tree, "test.xml")

        assert result["Run_MD5"] == "right"
        assert result["Run_Filename"] == "SRR2"

    def test_extract_metadata_fields_normalized_url_is_taken_by_semantic_name(self):
        """The normalized URL comes from the SRA Normalized element, wherever it is listed."""
        xml_content = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE_SET>
            <EXPERIMENT_PACKAGE>
                <RUN_SET>
                    <RUN accession="SRR3">
                        <SRAFiles>
                            <SRAFile filename="SRR3" semantic_name="run" url="https://example.org/run"/>
                            <SRAFile filename="SRR3.lite" semantic_name="SRA Lite" url="https://example.org/lite"/>
                            <SRAFile filename="SRR3.nm" semantic_name="SRA Normalized" url="https://example.org/norm"/>
                            <SRAFile filename="SRR3.x" semantic_name="Other" url="https://example.org/x"/>
                        </SRAFiles>
                    </RUN>
                </RUN_SET>
            </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""

        result = _extract_metadata_fields(ET.fromstring(xml_content), "test.xml")

        assert result["SRA_Normalized_URL"] == "https://example.org/norm"


class TestExtractSampleAttributes:
    """Test _extract_sample_attributes function."""

    def test_extract_sample_attributes_success(self):
        """Test successful sample attributes extraction."""
        xml_content = """<?xml version="1.0"?>
        <root>
            <SAMPLE_ATTRIBUTES>
                <SAMPLE_ATTRIBUTE>
                    <TAG>organism</TAG>
                    <VALUE>Escherichia coli</VALUE>
                </SAMPLE_ATTRIBUTE>
                <SAMPLE_ATTRIBUTE>
                    <TAG>isolation_source</TAG>
                    <VALUE>clinical isolate</VALUE>
                </SAMPLE_ATTRIBUTE>
            </SAMPLE_ATTRIBUTES>
        </root>"""

        tree = ET.fromstring(xml_content)
        # Pre-populate unique_attributes since the function only extracts existing attributes
        unique_attributes = {"organism", "isolation_source"}

        result = _extract_sample_attributes(tree, unique_attributes)

        assert result["organism"] == "Escherichia coli"
        assert result["isolation_source"] == "clinical isolate"

    def test_extract_sample_attributes_no_attributes(self):
        """Test when no sample attributes are present."""
        xml_content = """<?xml version="1.0"?><root></root>"""

        tree = ET.fromstring(xml_content)
        unique_attributes = set()

        result = _extract_sample_attributes(tree, unique_attributes)

        assert result == {}
        assert len(unique_attributes) == 0

    def test_extract_sample_attributes_empty_values(self):
        """Test handling empty attribute values."""
        xml_content = """<?xml version="1.0"?>
        <root>
            <SAMPLE_ATTRIBUTES>
                <SAMPLE_ATTRIBUTE>
                    <TAG>organism</TAG>
                    <VALUE></VALUE>
                </SAMPLE_ATTRIBUTE>
                <SAMPLE_ATTRIBUTE>
                    <TAG>missing_value</TAG>
                </SAMPLE_ATTRIBUTE>
            </SAMPLE_ATTRIBUTES>
        </root>"""

        tree = ET.fromstring(xml_content)
        unique_attributes = set()

        result = _extract_sample_attributes(tree, unique_attributes)

        # Empty values should not be included
        assert "organism" not in result
        assert "missing_value" not in result


class TestParseMetadata:
    """Test parse_metadata function."""

    def test_parse_metadata_success(self, tmp_path):
        """Test successful metadata parsing."""
        metadata_dir = tmp_path / "metadata"
        metadata_dir.mkdir()
        output_file = tmp_path / "parsed_metadata.tsv"

        # Create test XML files
        xml_content = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE_SET>
            <EXPERIMENT_PACKAGE>
                <SAMPLE>
                    <IDENTIFIERS>
                        <PRIMARY_ID>SAMN123</PRIMARY_ID>
                    </IDENTIFIERS>
                </SAMPLE>
                <RUN_SET>
                    <RUN>
                        <IDENTIFIERS>
                            <PRIMARY_ID>SRR123</PRIMARY_ID>
                        </IDENTIFIERS>
                    </RUN>
                </RUN_SET>
            </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""

        (metadata_dir / "SRR123_metadata.xml").write_text(xml_content)

        result = parse_metadata(metadata_dir, output_file)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1
        assert result.iloc[0]["Run_ID"] == "SRR123"
        assert result.iloc[0]["Sample_ID"] == "SAMN123"
        assert output_file.exists()

    def test_parse_metadata_no_xml_files(self, tmp_path):
        """Test when no XML files are found."""
        metadata_dir = tmp_path / "metadata"
        metadata_dir.mkdir()
        output_file = tmp_path / "parsed_metadata.tsv"

        with pytest.raises(DataAccessError, match="No metadata XML in"):
            parse_metadata(metadata_dir, output_file)
        assert not output_file.exists()

    def test_parse_metadata_invalid_xml(self, tmp_path):
        """Test handling invalid XML files."""
        metadata_dir = tmp_path / "metadata"
        metadata_dir.mkdir()
        output_file = tmp_path / "parsed_metadata.tsv"

        # Create invalid XML
        (metadata_dir / "SRR123_metadata.xml").write_text("invalid xml content")

        errors: list = []
        with patch("prepare.data.metadata.logger") as mock_logger:
            result = parse_metadata(metadata_dir, output_file, errors=errors)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0
        assert errors == ["SRR123_metadata.xml"]
        mock_logger.error.assert_called()

    def test_parse_metadata_nonexistent_folder(self):
        """Test handling nonexistent metadata folder."""
        with pytest.raises(ValidationError):
            parse_metadata("/nonexistent/folder", "/tmp/output.tsv")


class TestCheckMetadataAttributes:
    """Test check_metadata_attributes function."""

    def test_check_metadata_attributes_success(self, tmp_path):
        """Test successful metadata attributes checking."""
        input_file = tmp_path / "metadata.tsv"
        output_file = tmp_path / "attributes.txt"

        # Create test metadata file
        df = pd.DataFrame(
            {
                "Run_ID": ["SRR123", "SRR456"],
                "organism": ["E. coli", "S. aureus"],
                "isolation_source": ["clinical", "environmental"],
            }
        )
        df.to_csv(input_file, sep="\t", index=False)

        result = check_metadata_attributes(input_file, output_file)

        assert isinstance(result, dict)
        assert "organism" in result
        assert "isolation_source" in result
        assert result["organism"] == 2
        assert result["isolation_source"] == 2
        assert output_file.exists()

    def test_check_metadata_attributes_empty_file(self, tmp_path):
        """Test handling empty metadata file."""
        input_file = tmp_path / "empty.tsv"
        output_file = tmp_path / "attributes.txt"

        # Create file with just headers
        input_file.write_text("Run_ID\tSample_ID\n")

        result = check_metadata_attributes(input_file, output_file)

        assert isinstance(result, dict)
        assert len(result) == 0

    def test_check_metadata_attributes_nonexistent_file(self):
        """Test handling nonexistent input file."""
        with pytest.raises(DataAccessError):
            check_metadata_attributes("/nonexistent/file.tsv", "/tmp/output.txt")


class TestParseMetadataXml:
    """Test parse_metadata_xml: single-file counterpart to parse_metadata used by the sidecar."""

    NCBI_XML = """<?xml version="1.0"?>
        <EXPERIMENT_PACKAGE_SET>
            <EXPERIMENT_PACKAGE>
                <EXPERIMENT>
                    <IDENTIFIERS>
                        <PRIMARY_ID>EXP1</PRIMARY_ID>
                    </IDENTIFIERS>
                    <LIBRARY_DESCRIPTOR>
                        <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
                        <LIBRARY_LAYOUT>
                            <PAIRED/>
                        </LIBRARY_LAYOUT>
                    </LIBRARY_DESCRIPTOR>
                    <PLATFORM>
                        <ILLUMINA>
                            <INSTRUMENT_MODEL>Illumina HiSeq 2500</INSTRUMENT_MODEL>
                        </ILLUMINA>
                    </PLATFORM>
                </EXPERIMENT>
                <RUN_SET>
                    <RUN accession="SRR1" total_spots="47964651" total_bases="14389395300" size="4744553813">
                        <IDENTIFIERS>
                            <PRIMARY_ID>SRR1</PRIMARY_ID>
                        </IDENTIFIERS>
                        <SRAFiles>
                            <SRAFile filename="SRR1" md5="abc" semantic_name="run"/>
                        </SRAFiles>
                    </RUN>
                </RUN_SET>
            </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""

    def test_parse_metadata_xml_reads_ncbi_attributes(self, tmp_path):
        """A single metadata XML file parses to the same fields _extract_metadata_fields gives."""
        xml_file = tmp_path / "SRR1.xml"
        xml_file.write_text(self.NCBI_XML)

        result = parse_metadata_xml(xml_file)

        assert result["Run_Total_Spots"] == "47964651"
        assert result["Run_Total_Bases"] == "14389395300"
        assert result["Run_Size"] == "4744553813"
        assert result["Run_MD5"] == "abc"
        assert result["Run_Filename"] == "SRR1"
        assert result["Experiment_Library_Layout"] == "PAIRED"

    def test_parse_metadata_xml_missing_file_raises_oserror(self):
        """A missing metadata file raises OSError."""
        with pytest.raises(OSError, match="Metadata XML file not found"):
            parse_metadata_xml("/nonexistent/SRR1.xml")

    def test_parse_metadata_xml_unparsable_file_raises_parseerror(self, tmp_path):
        """A malformed XML file raises ET.ParseError."""
        import xml.etree.ElementTree as ET

        xml_file = tmp_path / "bad.xml"
        xml_file.write_text("not xml at all <<<")

        with pytest.raises(ET.ParseError):
            parse_metadata_xml(xml_file)

    def test_library_fields_under_design(self, tmp_path):
        xml = (
            '<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><EXPERIMENT accession="SRX1"><DESIGN>'
            "<LIBRARY_DESCRIPTOR><LIBRARY_NAME>lib</LIBRARY_NAME><LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>"
            "<LIBRARY_SOURCE>METAGENOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION>"
            "<LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT></LIBRARY_DESCRIPTOR></DESIGN>"
            "<PLATFORM><ILLUMINA/></PLATFORM></EXPERIMENT>"
            '<RUN_SET><RUN accession="SRR1" total_spots="5"/></RUN_SET></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>'
        )
        path = tmp_path / "SRR1_metadata.xml"
        path.write_text(xml)
        fields = parse_metadata_xml(path)
        assert fields["Experiment_Library_Name"] == "lib"
        assert fields["Experiment_Library_Strategy"] == "WGS"
        assert fields["Experiment_Library_Source"] == "METAGENOMIC"
        assert fields["Experiment_Library_Selection"] == "RANDOM"


if __name__ == "__main__":
    pytest.main([__file__])


# --- Narrow exception handling (tech debt): a programming error must not be swallowed ---


def test_data_access_error_from_download_metadata_chains_the_cause(tmp_path):
    """Kind (c): the wrapped error keeps its cause."""
    with pytest.raises(DataAccessError) as exc:
        download_metadata(
            email="a@b.c",
            matches_folder=tmp_path,
            metadata_folder=tmp_path / "metadata",
            accessions_file=tmp_path / "missing.txt",
        )
    assert isinstance(exc.value.__cause__.__cause__, OSError)


def test_unexpected_error_while_parsing_one_metadata_file_propagates(tmp_path):
    """Kind (a): a bug in field extraction is no longer logged and skipped."""
    (tmp_path / "SRR1_metadata.xml").write_text("<EXPERIMENT_PACKAGE_SET/>")
    with patch("prepare.data.metadata._extract_metadata_fields", side_effect=TypeError("bug")):
        with pytest.raises(TypeError):
            parse_metadata(tmp_path, tmp_path / "out.tsv")


def test_malformed_metadata_file_is_still_skipped(tmp_path):
    """Kind (a): an XML syntax error skips that file only."""
    (tmp_path / "SRR1_metadata.xml").write_text("<not xml")
    assert parse_metadata(tmp_path, tmp_path / "out.tsv").empty


def test_incomplete_read_in_a_batch_is_retried_then_reported_as_failed(tmp_path):
    """Fix round 1: a truncated NCBI reply takes the retry path instead of escaping raw."""
    from http.client import IncompleteRead

    handle = MagicMock()
    handle.read.side_effect = IncompleteRead(b"<EXPERIMENT")
    with (
        patch("prepare.data.metadata.Entrez.efetch", return_value=handle),
        patch("prepare.data.metadata_fetch._pace_requests"),
        patch("time.sleep") as mock_sleep,
    ):
        successes, failures = _download_batch_metadata(["SRR1", "SRR2"], tmp_path, "a@b.c", None)
    assert successes == {}
    assert set(failures) == {"SRR1", "SRR2"}
    assert mock_sleep.call_count == 3


def test_incomplete_read_for_one_accession_is_retried_then_reported_as_failed(tmp_path):
    """Fix round 1: the single-accession path retries a truncated reply the same way."""
    from http.client import IncompleteRead

    handle = MagicMock()
    handle.read.side_effect = IncompleteRead(b"<EXPERIMENT")
    with (
        patch("prepare.data.metadata.Entrez.efetch", return_value=handle),
        patch("prepare.data.metadata_fetch._pace_requests"),
        patch("time.sleep") as mock_sleep,
    ):
        success, message = _download_single_metadata("SRR1", tmp_path, "a@b.c")
    assert success is False
    assert "IncompleteRead" in message
    assert mock_sleep.call_count == 3


def _single_download(tmp_path, body, accession="SRR2"):
    handle = MagicMock()
    handle.read.return_value = body.encode()
    with (
        patch("prepare.data.metadata.Entrez.efetch", return_value=handle),
        patch("prepare.data.metadata_fetch._pace_requests"),
        patch("time.sleep"),
    ):
        return _download_single_metadata(accession, tmp_path, "a@b.c")


def test_single_download_keeps_only_the_requested_run(tmp_path):
    """Audit M10: a package listing SRR1 and SRR2 is pruned to SRR2 before it is written."""
    root = ET.fromstring(_package_xml([("SRR1", {})]))
    second = ET.SubElement(root.find(".//RUN_SET"), "RUN", accession="SRR2")
    ET.SubElement(ET.SubElement(second, "IDENTIFIERS"), "PRIMARY_ID").text = "SRR2"
    success, result = _single_download(tmp_path, ET.tostring(root, encoding="unicode"))
    assert success is True
    assert result.name == "SRR2_metadata.xml"
    fields = _extract_metadata_fields(ET.parse(result), result)
    assert fields["Run_ID"] == "SRR2"


@pytest.mark.parametrize(
    "body, expected",
    [
        ("<html>error</html>", "no package for SRR2"),
        ("<html><body>Error", "non-XML"),
        (_package_xml([("SRR9", {})]), "no package for SRR2"),
    ],
)
def test_single_download_refuses_bodies_without_the_package(tmp_path, body, expected):
    """An HTTP 200 error page or a package set without the accession writes nothing."""
    success, message = _single_download(tmp_path, body)
    assert success is False
    assert expected in message
    assert list(tmp_path.iterdir()) == []


# --- Single-pass metadata parsing (performance plan, task 6) ---


def _reference_fields(tree):
    """``_extract_metadata_fields`` as it was before the single-pass change, copied as the reference."""

    def run_attr(element, name):
        return None if element is None else element.get(name)

    def first_srafile(elements):
        for element in elements:
            if element.get("semantic_name") == "run":
                return element
        return elements[0] if elements else None

    def first_child_tag(element):
        if element is None:
            return None
        children = list(element)
        return children[0].tag if children else None

    def spot_length():
        averages = []
        for read in tree.findall(".//RUN/Statistics/Read"):
            average = read.get("average")
            if average is None:
                continue
            try:
                averages.append(float(average))
            except ValueError:
                continue
        if not averages:
            return None
        total = sum(averages)
        return str(int(total)) if total.is_integer() else str(total)

    t = tree.findtext
    run_element = tree.find(".//RUN")
    srafiles = tree.findall(".//RUN/SRAFiles/SRAFile")
    srafile = first_srafile(srafiles)
    return {
        "Run_ID": t(".//RUN/IDENTIFIERS/PRIMARY_ID"),
        "Run_Total_Spots": run_attr(run_element, "total_spots") or t(".//RUN/Total_spots"),
        "Run_Total_Bases": run_attr(run_element, "total_bases") or t(".//RUN/Total_bases"),
        "Run_Size": run_attr(run_element, "size") or t(".//RUN/size"),
        "Run_Download_Path": t(".//RUN/download_path"),
        "Run_MD5": (srafile.get("md5") if srafile is not None else None) or t(".//RUN/md5"),
        "Run_Filename": (srafile.get("filename") if srafile is not None else None) or t(".//RUN/filename"),
        "Run_Spot_Length": spot_length() or t(".//RUN/spot_length"),
        "Run_Reads": t(".//RUN/reads"),
        "Run_FTP": t(".//RUN/ftp"),
        "Run_Aspera": t(".//RUN/aspera"),
        "Run_Galaxy": t(".//RUN/galaxy"),
        "Project_ID": t(".//STUDY/IDENTIFIERS/PRIMARY_ID"),
        "Project_Title": t(".//STUDY/DESCRIPTOR/STUDY_TITLE"),
        "Project_Abstract": t(".//STUDY/DESCRIPTOR/STUDY_ABSTRACT"),
        "Sample_ID": t(".//SAMPLE/IDENTIFIERS/PRIMARY_ID"),
        "Sample_External_ID": t(".//SAMPLE/IDENTIFIERS/EXTERNAL_ID"),
        "Sample_Name": t(".//SAMPLE/SAMPLE_NAME/TAXON_ID"),
        "Sample_Scientific_Name": t(".//SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME"),
        "Sample_Title": t(".//SAMPLE/TITLE"),
        "Experiment_ID": t(".//EXPERIMENT/IDENTIFIERS/PRIMARY_ID"),
        "Experiment_Title": t(".//EXPERIMENT/TITLE"),
        "Experiment_Design": t(".//EXPERIMENT/DESIGN/DESIGN_DESCRIPTION"),
        "Experiment_Library_Name": t(".//EXPERIMENT//LIBRARY_DESCRIPTOR/LIBRARY_NAME"),
        "Experiment_Library_Strategy": t(".//EXPERIMENT//LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY"),
        "Experiment_Library_Source": t(".//EXPERIMENT//LIBRARY_DESCRIPTOR/LIBRARY_SOURCE"),
        "Experiment_Library_Selection": t(".//EXPERIMENT//LIBRARY_DESCRIPTOR/LIBRARY_SELECTION"),
        "Experiment_Library_Layout": first_child_tag(tree.find(".//LIBRARY_LAYOUT")),
        "Platform": first_child_tag(tree.find(".//PLATFORM")),
        "SRA_Normalized_URL": srafiles[1].get("url") if len(srafiles) > 1 else None,
    }


def _reference_parse_metadata(folder):
    """``parse_metadata`` as it was before the single-pass change (two parses per file, dense rows)."""
    from lxml import etree

    files = sorted(p for p in Path(folder).glob("*.xml") if not p.name.startswith("."))
    unique = set()
    for xml_file in files:
        try:
            tree = etree.parse(str(xml_file))
        except etree.XMLSyntaxError:
            continue
        for tag in tree.findall(".//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/TAG"):
            unique.add(tag.text)
    unique_attributes = sorted(unique)
    records = []
    for xml_file in files:
        try:
            tree = etree.parse(str(xml_file))
        except etree.XMLSyntaxError:
            continue
        record = _reference_fields(tree)
        attributes = {}
        for attribute in tree.findall(".//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
            tag = attribute.findtext("TAG")
            if tag in unique_attributes:
                attributes[tag] = attribute.findtext("VALUE")
        for name in unique_attributes:
            record[name] = attributes.get(name, None)
        records.append(record)
    return pd.DataFrame(records)


_EDGE_CASE_XML = {
    # Only a run identifier: every other fixed field is None.
    "SRR900001": "<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><RUN_SET><RUN><IDENTIFIERS>"
    "<PRIMARY_ID>SRR900001</PRIMARY_ID></IDENTIFIERS></RUN></RUN_SET></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>",
    # Empty elements read as "" (not None); child-element fallbacks for size/md5/spots; a duplicated tag
    # (the later value wins); an attribute without VALUE; a tag equal to a fixed column name.
    "SRR900002": "<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE><EXPERIMENT><TITLE/><PLATFORM><OXFORD_NANOPORE/>"
    "</PLATFORM></EXPERIMENT><SAMPLE><TITLE></TITLE><SAMPLE_ATTRIBUTES>"
    "<SAMPLE_ATTRIBUTE><TAG>host</TAG><VALUE>first</VALUE></SAMPLE_ATTRIBUTE>"
    "<SAMPLE_ATTRIBUTE><TAG>host</TAG><VALUE>second</VALUE></SAMPLE_ATTRIBUTE>"
    "<SAMPLE_ATTRIBUTE><TAG>no_value</TAG></SAMPLE_ATTRIBUTE>"
    "<SAMPLE_ATTRIBUTE><TAG>Platform</TAG><VALUE>from attribute</VALUE></SAMPLE_ATTRIBUTE>"
    "<SAMPLE_ATTRIBUTE><TAG>empty_value</TAG><VALUE/></SAMPLE_ATTRIBUTE>"
    "</SAMPLE_ATTRIBUTES></SAMPLE><RUN_SET><RUN><IDENTIFIERS><PRIMARY_ID>SRR900002</PRIMARY_ID></IDENTIFIERS>"
    "<size>77</size><md5>m5</md5><Total_spots>9</Total_spots><spot_length>150</spot_length>"
    "<Statistics><Read average='75.5'/><Read average='x'/><Read/></Statistics>"
    "</RUN></RUN_SET></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>",
    # Two packages: every field comes from the first match in document order.
    "SRR900003": "<EXPERIMENT_PACKAGE_SET>"
    "<EXPERIMENT_PACKAGE><STUDY><DESCRIPTOR><STUDY_TITLE>first study</STUDY_TITLE></DESCRIPTOR></STUDY>"
    "<RUN_SET><RUN accession='SRR900003' total_spots='5'><IDENTIFIERS><PRIMARY_ID>SRR900003</PRIMARY_ID>"
    "</IDENTIFIERS><SRAFiles><SRAFile filename='a' md5='1' url='u1'/><SRAFile filename='b' md5='2' url='u2'/>"
    "</SRAFiles></RUN></RUN_SET></EXPERIMENT_PACKAGE>"
    "<EXPERIMENT_PACKAGE><STUDY><IDENTIFIERS><PRIMARY_ID>SRP2</PRIMARY_ID></IDENTIFIERS><DESCRIPTOR>"
    "<STUDY_TITLE>second study</STUDY_TITLE></DESCRIPTOR></STUDY><SAMPLE><SAMPLE_ATTRIBUTES><SAMPLE_ATTRIBUTE>"
    "<TAG>strain</TAG><VALUE>K-12</VALUE></SAMPLE_ATTRIBUTE></SAMPLE_ATTRIBUTES></SAMPLE>"
    "<RUN_SET><RUN accession='SRR900004'><IDENTIFIERS><PRIMARY_ID>SRR900004</PRIMARY_ID></IDENTIFIERS></RUN>"
    "</RUN_SET></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>",
    # Malformed: skipped by both.
    "SRR900005": "<EXPERIMENT_PACKAGE_SET><unclosed>",
}


@pytest.fixture
def mixed_metadata_folder(tmp_path):
    """Synthetic NCBI-like files with overlapping and disjoint attribute sets, plus the edge cases above."""
    from tests.perf_fixtures import write_metadata_folder

    folder = tmp_path / "metadata"
    write_metadata_folder(folder, count=25, per_file=14, pool=30)
    for accession, xml in _EDGE_CASE_XML.items():
        (folder / f"{accession}_metadata.xml").write_text(xml)
    return folder


def _rows(df):
    """The table's cells as plain lists, so None and "" are compared exactly."""
    return [list(row) for row in df.itertuples(index=False, name=None)]


def test_parse_metadata_output_matches_the_two_pass_reference(mixed_metadata_folder, tmp_path):
    """The single-pass table has the same columns, order, dtypes and cell values as the earlier two-pass one."""
    expected = _reference_parse_metadata(mixed_metadata_folder)
    result = parse_metadata(mixed_metadata_folder, tmp_path / "out.tsv")

    assert list(result.columns) == list(expected.columns)
    assert list(result.dtypes) == list(expected.dtypes)
    assert _rows(result) == _rows(expected)
    expected.to_csv(tmp_path / "expected.tsv", sep="\t", index=False)
    assert (tmp_path / "out.tsv").read_text() == (tmp_path / "expected.tsv").read_text()
    # Checks on the fixture itself, so the comparison above covers what it is meant to.
    assert len(result) == 28
    assert result.columns[:2].tolist() == ["Run_ID", "Run_Total_Spots"]
    assert result.columns[29] == "SRA_Normalized_URL"
    assert list(result.columns[30:]) == sorted(result.columns[30:])
    row = result.set_index("Run_ID").loc["SRR900002"]
    assert row["Sample_Title"] == "" and row["Experiment_Title"] == ""
    assert row["host"] == "second" and row["no_value"] is None and row["Platform"] == "from attribute"


def test_parse_metadata_parses_each_file_once(mixed_metadata_folder, tmp_path):
    """Each XML file is parsed by lxml exactly once (the attribute scan no longer re-reads the folder)."""
    import prepare.data.metadata as metadata_module

    real_parse = metadata_module.etree.parse
    with patch.object(metadata_module.etree, "parse", side_effect=real_parse) as spy:
        parse_metadata(mixed_metadata_folder, tmp_path / "out.tsv")

    parsed = sorted(Path(call_args.args[0]).name for call_args in spy.call_args_list)
    assert parsed == sorted(p.name for p in mixed_metadata_folder.glob("*.xml"))


def test_parse_metadata_xml_matches_the_folder_row(tmp_path):
    """The single-file parse returns the folder table's row for that file, minus attributes it lacks.

    The folder holds no attribute named like a fixed column (the SRR900002 case), since such a tag
    replaces that column for every row of the folder table but only for its own file here.
    """
    from tests.perf_fixtures import write_metadata_folder

    def missing_as_none(value):
        return None if not isinstance(value, str) and pd.isna(value) else value

    folder = tmp_path / "metadata"
    paths = write_metadata_folder(folder, count=8, per_file=14, pool=30)
    for accession in ("SRR900001", "SRR900003"):
        paths.append(folder / f"{accession}_metadata.xml")
        paths[-1].write_text(_EDGE_CASE_XML[accession])
    table = parse_metadata(folder, tmp_path / "out.tsv").set_index("Run_ID", drop=False)
    for path in paths:
        single = parse_metadata_xml(path)
        row = {key: missing_as_none(value) for key, value in table.loc[single["Run_ID"]].items()}
        assert single == {key: row[key] for key in single}
        assert all(value is None for key, value in row.items() if key not in single)


def test_download_metadata_lists_the_metadata_folder_once(tmp_path):
    """Already-downloaded accessions are found with one directory listing, not one stat per accession."""
    import prepare.data.metadata as metadata_module

    metadata_dir = tmp_path / "metadata"
    metadata_dir.mkdir()
    (metadata_dir / "SRR2_metadata.xml").write_text("<x/>")
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("SRR1\nSRR2\nSRR3\n")

    real_listdir = os.listdir
    with (
        patch.object(metadata_module.os, "listdir", side_effect=real_listdir) as listdir_spy,
        patch.object(metadata_module, "_download_accessions_metadata", return_value={}) as download,
    ):
        download_metadata(
            email="a@b.c", matches_folder=tmp_path, metadata_folder=metadata_dir, accessions_file=accessions
        )

    assert sorted(download.call_args.args[0]) == ["SRR1", "SRR3"]
    assert [Path(c.args[0]) for c in listdir_spy.call_args_list] == [metadata_dir]
