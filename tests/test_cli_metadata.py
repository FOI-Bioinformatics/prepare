"""Tests for the metadata CLI commands' handling of parsed metadata (single-pass parsing, task 6)."""

import argparse
import json
from unittest.mock import patch

import pandas as pd

import prepare.data.metadata as metadata_module
from prepare.cli.commands.metadata import DownloadMetadataCommand, ParseMetadataCommand
from tests.perf_fixtures import write_metadata_folder


def _parse_args(tmp_path, folder):
    return argparse.Namespace(
        metadata_folder=str(folder),
        metadata_table_file=str(tmp_path / "metadata_table.txt"),
        registry=str(tmp_path / "prepare_registry.json"),
    )


def test_parse_metadata_command_records_every_run_without_iterating_rows_as_series(tmp_path):
    """Every run is recorded from the table without building one pandas Series per row."""
    folder = tmp_path / "metadata"
    write_metadata_folder(folder, count=6, per_file=14, pool=30)

    with patch.object(pd.DataFrame, "iterrows", side_effect=AssertionError("iterrows used")):
        assert ParseMetadataCommand().execute(_parse_args(tmp_path, folder)) == 0

    datasets = json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]
    assert sorted(datasets) == [f"SRR{1000000 + i}" for i in range(6)]
    first = datasets["SRR1000000"]["metadata"]
    assert first["run_total_spots"] == 1000000
    assert first["run_md5"] == f"{1:032x}"
    assert first["organism"] == "soil metagenome"
    assert first["collection_date"] == "value 0 of collection_date"
    assert first["library_layout"] == "SINGLE"
    assert first["platform"] == "ILLUMINA"
    assert first["xml"] == "metadata/SRR1000000_metadata.xml"


def test_parse_metadata_summary_counts_a_dropped_record(tmp_path, caplog):
    """A record that fails is dropped; the exit stays 0 and the summary line says how many."""
    import prepare.cli.commands.metadata as command_module

    folder = tmp_path / "metadata"
    write_metadata_folder(folder, count=3, per_file=14, pool=30)
    real_record = command_module.record_metadata

    def flaky_record(registry, accession, **kwargs):
        if accession == "SRR1000001":
            raise ValueError("bad record")
        return real_record(registry, accession, **kwargs)

    with patch.object(command_module, "record_metadata", side_effect=flaky_record):
        with caplog.at_level("INFO", logger="prepare.cli.commands.metadata"):
            assert ParseMetadataCommand().execute(_parse_args(tmp_path, folder)) == 0

    datasets = json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]
    assert sorted(datasets) == ["SRR1000000", "SRR1000002"]
    assert "Recorded metadata for 2 accession(s) in the registry; 1 dropped" in caplog.text


def test_parse_metadata_command_records_nothing_for_an_empty_folder(tmp_path):
    """An empty folder is an error: exit 1 and no registry written."""
    folder = tmp_path / "metadata"
    folder.mkdir()
    assert ParseMetadataCommand().execute(_parse_args(tmp_path, folder)) == 1
    assert not (tmp_path / "prepare_registry.json").exists()


def test_download_metadata_command_parses_each_downloaded_file_once(tmp_path):
    """Recording the fields of freshly downloaded files costs one lxml parse per file."""
    folder = tmp_path / "metadata"
    paths = write_metadata_folder(folder, count=3, per_file=14, pool=30)
    downloaded = {path.name.split("_")[0]: path for path in paths}
    args = argparse.Namespace(
        email="a@b.c",
        matches_folder=str(tmp_path / "matches"),
        metadata_folder=str(folder),
        threshold=0.0,
        dry_run=False,
        accessions_file=None,
        api_key=None,
        batch_size=200,
        registry=str(tmp_path / "prepare_registry.json"),
        data_root=None,
    )

    real_parse = metadata_module.etree.parse
    with (
        patch("prepare.cli.commands.metadata.download_metadata", return_value=downloaded),
        patch.object(metadata_module.etree, "parse", side_effect=real_parse) as spy,
    ):
        assert DownloadMetadataCommand().execute(args) == 0

    assert sorted(str(c.args[0]) for c in spy.call_args_list) == sorted(str(p) for p in paths)
    datasets = json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]
    assert datasets["SRR1000002"]["metadata"]["run_total_spots"] == 1000026


def test_download_metadata_signal_between_accessions_stops_before_the_next_one(tmp_path):
    """A stop noticed after the first downloaded file is parsed ends the run before the
    second downloaded file's XML is ever parsed, even though both were already fetched."""
    folder = tmp_path / "metadata"
    paths = write_metadata_folder(folder, count=2, per_file=14, pool=30)
    downloaded = {path.name.split("_")[0]: path for path in paths}
    args = argparse.Namespace(
        email="a@b.c",
        matches_folder=str(tmp_path / "matches"),
        metadata_folder=str(folder),
        threshold=0.0,
        dry_run=False,
        accessions_file=None,
        api_key=None,
        batch_size=200,
        registry=str(tmp_path / "prepare_registry.json"),
        data_root=None,
    )

    real_parse = metadata_module.parse_metadata_xml
    seen = []

    def spy(xml_path):
        result = real_parse(xml_path)
        seen.append(xml_path)
        if len(seen) == 1:
            args._termination.stop.set()
        return result

    with (
        patch("prepare.cli.commands.metadata.download_metadata", return_value=downloaded),
        patch("prepare.cli.commands.metadata.parse_metadata_xml", side_effect=spy),
    ):
        rc = DownloadMetadataCommand().run(args)

    assert rc == 130
    assert len(seen) == 1
    datasets = json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]
    assert len(datasets) == 1


def test_parse_metadata_signal_between_rows_stops_before_the_next_one(tmp_path):
    """A stop noticed after the first row is parsed ends the run before the second row is
    ever added to what gets recorded."""
    folder = tmp_path / "metadata"
    write_metadata_folder(folder, count=2, per_file=14, pool=30)
    cmd = ParseMetadataCommand()
    args = _parse_args(tmp_path, folder)

    real_row_record = cmd._row_record
    seen = []

    def spy(metadata_folder, row):
        result = real_row_record(metadata_folder, row)
        seen.append(row.get("Run_ID"))
        if len(seen) == 1:
            args._termination.stop.set()
        return result

    with patch.object(cmd, "_row_record", side_effect=spy):
        rc = cmd.run(args)

    assert rc == 130
    assert len(seen) == 1
    datasets = json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]
    assert len(datasets) == 1


def test_download_metadata_exits_4_when_ncbi_cannot_be_reached(tmp_path):
    """Every request failing for a network reason (after its retries) is a retryable failure."""
    from urllib.error import URLError

    accessions = tmp_path / "accessions.txt"
    accessions.write_text("SRR1\nSRR2\nSRR3\n")
    args = argparse.Namespace(
        email="a@b.c",
        matches_folder=str(tmp_path / "matches"),
        metadata_folder=str(tmp_path / "metadata"),
        threshold=0.0,
        dry_run=False,
        accessions_file=str(accessions),
        api_key=None,
        batch_size=2,
        registry=str(tmp_path / "prepare_registry.json"),
        data_root=None,
    )
    with (
        patch("prepare.data.metadata.Entrez.efetch", side_effect=URLError("connection refused")),
        patch("prepare.data.metadata_fetch._pace_requests"),
        patch("prepare.data.metadata.time.sleep"),
    ):
        assert DownloadMetadataCommand().execute(args) == 4


def test_download_metadata_does_not_raise_when_a_failure_is_not_a_network_one(tmp_path):
    """An accession NCBI answered without is not a network failure, so the run is not retryable."""
    batches = [({}, {"SRR1": "network: connection refused"}), ({}, {"SRR2": "not in the NCBI response"})]
    with patch("prepare.data.metadata._download_batch_metadata", side_effect=batches):
        result = metadata_module._download_accessions_metadata(["SRR1", "SRR2"], tmp_path, "a@b.c", 2, batch_size=1)
    assert result == {}


def _download_args(tmp_path, **overrides):
    args = dict(
        email="a@b.c",
        matches_folder=str(tmp_path / "matches"),
        metadata_folder=str(tmp_path / "metadata"),
        threshold=0.0,
        dry_run=False,
        accessions_file=None,
        api_key=None,
        batch_size=200,
        registry=str(tmp_path / "prepare_registry.json"),
        data_root=None,
    )
    args.update(overrides)
    return argparse.Namespace(**args)


def test_download_metadata_exit_1_when_an_accession_has_no_metadata(tmp_path, caplog):
    """NCBI answering for one of two accessions: the other is named and the exit is 1."""
    import logging

    accessions = tmp_path / "accessions.txt"
    accessions.write_text("SRR1000000\nSRR1000001\n")
    folder = tmp_path / "metadata"
    (good,) = write_metadata_folder(tmp_path / "fetched", count=1, per_file=4, pool=10)
    batches = [({"SRR1000000": good}, {"SRR1000001": "not in the NCBI response"})]
    with (
        patch("prepare.data.metadata._download_batch_metadata", side_effect=batches),
        caplog.at_level(logging.ERROR),
    ):
        assert DownloadMetadataCommand().execute(_download_args(tmp_path, accessions_file=str(accessions))) == 1
    assert "SRR1000001" in caplog.text
    datasets = json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]
    assert "SRR1000000" in datasets  # the fetched one is still recorded
    assert folder.is_dir()


def test_download_metadata_exits_4_when_some_accessions_failed_for_a_network_reason(tmp_path, caplog):
    """Minor 6: the fetched one is recorded, the missing one is named, and the exit says retry."""
    import logging

    accessions = tmp_path / "accessions.txt"
    accessions.write_text("SRR1000000\nSRR1000001\n")
    (good,) = write_metadata_folder(tmp_path / "fetched", count=1, per_file=4, pool=10)
    batches = [({"SRR1000000": good}, {"SRR1000001": "network: connection refused"})]
    with (
        patch("prepare.data.metadata._download_batch_metadata", side_effect=batches),
        caplog.at_level(logging.ERROR),
    ):
        assert DownloadMetadataCommand().execute(_download_args(tmp_path, accessions_file=str(accessions))) == 4
    assert "SRR1000001" in caplog.text
    assert "SRR1000000" in json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]


def test_download_metadata_failure_line_names_at_most_ten_accessions(tmp_path, caplog):
    """Minor 6: twelve missing accessions give a line with ten names and a count of the rest."""
    import logging

    missing = [f"SRR{2000000 + i}" for i in range(12)]
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("".join(f"{acc}\n" for acc in missing))
    batches = [({}, {acc: "not in the NCBI response" for acc in missing})]
    with (
        patch("prepare.data.metadata._download_batch_metadata", side_effect=batches),
        caplog.at_level(logging.ERROR),
    ):
        assert DownloadMetadataCommand().execute(_download_args(tmp_path, accessions_file=str(accessions))) == 1
    line = next(
        r.getMessage() for r in caplog.records if r.levelno == logging.ERROR and "No metadata" in r.getMessage()
    )
    assert line.count("SRR") == 10 and line.endswith("and 2 more")


def test_download_metadata_unreadable_matches_csv_is_exit_1(tmp_path):
    """A match file in another format stops the run with exit 1."""
    matches = tmp_path / "matches"
    matches.mkdir()
    (matches / "g.csv").write_text("run;score\nSRR1;0.9\n")
    with patch("prepare.data.metadata._download_batch_metadata") as fetch:
        assert DownloadMetadataCommand().execute(_download_args(tmp_path)) == 1
    fetch.assert_not_called()


def test_parse_metadata_exit_1_on_a_broken_xml(tmp_path):
    """One of two XML files fails to parse: the table holds the good one and the exit is 1."""
    folder = tmp_path / "metadata"
    write_metadata_folder(folder, count=1, per_file=4, pool=10)
    (folder / "SRR9_metadata.xml").write_text("<EXPERIMENT_PACKAGE_SET><RUN")
    args = argparse.Namespace(
        metadata_folder=str(folder),
        metadata_table_file=str(tmp_path / "metadata_table.txt"),
        registry=str(tmp_path / "prepare_registry.json"),
    )
    assert ParseMetadataCommand().execute(args) == 1
    table = pd.read_csv(tmp_path / "metadata_table.txt", sep="\t")
    assert list(table["Run_ID"]) == ["SRR1000000"]


def test_parse_metadata_with_no_xml_writes_no_table_and_exits_1(tmp_path, caplog):
    """An empty metadata folder leaves an older table as it is and exits 1 with a message."""
    import logging

    folder = tmp_path / "metadata"
    folder.mkdir()
    table = tmp_path / "metadata_table.txt"
    table.write_text("Run_ID\nSRR_OLD\n")
    args = argparse.Namespace(
        metadata_folder=str(folder), metadata_table_file=str(table), registry=str(tmp_path / "prepare_registry.json")
    )
    with caplog.at_level(logging.ERROR):
        assert ParseMetadataCommand().execute(args) == 1
    assert "No metadata XML in" in caplog.text
    assert table.read_text() == "Run_ID\nSRR_OLD\n"
    assert not (tmp_path / "prepare_registry.json").exists()
