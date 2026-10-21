"""
Tests for prepare.store.sidecar: the per-dataset JSON sidecar
``<root>/sra/<ACC>/<ACC>.json`` that makes the store self-describing.
"""

import hashlib
import json
import os
import logging
from pathlib import Path

import pytest

from prepare.core.exceptions import DataAccessError
from prepare.store.sidecar import (
    SIDECAR_SCHEMA,
    Sidecar,
    build_sidecar,
    ncbi_from_metadata_xml,
    read_sidecar,
    write_sidecar,
)

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
            <RUN accession="SRR1" total_spots="10" total_bases="1400" size="4744553813">
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


def _write_fastq(path, n_records):
    """Write ``n_records`` plain FASTQ records (4 lines each) to ``path``."""
    lines = []
    for i in range(n_records):
        lines.append(f"@read{i}\nACGT\n+\nIIII\n")
    path.write_text("".join(lines))


def _make_sidecar(**overrides):
    fields = dict(
        accession="SRR1",
        state="complete",
        layout="PAIRED",
        downloaded="2026-09-06T00:00:00+00:00",
        tool="fasterq-dump",
        tool_version="3.0.0",
        compression="gzip",
        files=[{"name": "SRR1_1.fastq.gz", "bytes": 100, "md5": "deadbeef", "reads": 10}],
        reads_per_mate=10,
        bases_total=None,
        ncbi={"spots": 10},
        completeness={"method": "spots", "ratio": 1.0, "verdict": "complete"},
        stats={},
        stats_computed=None,
    )
    fields.update(overrides)
    return Sidecar(**fields)


def test_sidecar_round_trip(tmp_path):
    sidecar = _make_sidecar()
    path = tmp_path / "SRR1.json"

    written_path = write_sidecar(path, sidecar)
    assert written_path == path
    assert path.is_file()

    loaded = read_sidecar(path)
    assert loaded == sidecar


def test_write_sidecar_is_atomic_and_sorted(tmp_path):
    sidecar = _make_sidecar()
    path = tmp_path / "SRR1.json"
    write_sidecar(path, sidecar)

    raw = path.read_text()
    data = json.loads(raw)
    assert data["accession"] == "SRR1"
    assert data["schema"] == SIDECAR_SCHEMA
    # sorted keys means "accession" (first alphabetically among our keys) precedes "tool"
    assert raw.index('"accession"') < raw.index('"tool"')
    # no leftover temp file
    leftovers = [p for p in tmp_path.iterdir() if p.name != "SRR1.json"]
    assert leftovers == []


def test_read_sidecar_missing_file_returns_none(tmp_path, caplog):
    path = tmp_path / "missing.json"
    with caplog.at_level(logging.WARNING):
        result = read_sidecar(path)
    assert result is None


def test_read_sidecar_invalid_json_returns_none(tmp_path, caplog):
    path = tmp_path / "bad.json"
    path.write_text("{not valid json")
    with caplog.at_level(logging.WARNING):
        result = read_sidecar(path)
    assert result is None
    assert "bad.json" in caplog.text


def test_read_sidecar_json_list_raises_data_access_error(tmp_path):
    """A sidecar holding a JSON list is reported as a store error, not an AttributeError."""
    path = tmp_path / "SRR1.json"
    path.write_text(json.dumps([{"accession": "SRR1"}]))
    with pytest.raises(DataAccessError, match=r"SRR1\.json: sidecar is not a JSON object"):
        read_sidecar(path)


def test_sidecar_from_dict_tolerant_of_missing_keys():
    sidecar = Sidecar.from_dict({})
    assert sidecar.accession == ""
    assert sidecar.schema == SIDECAR_SCHEMA
    assert sidecar.files == []


def test_sidecar_from_dict_ignores_unknown_keys():
    sidecar = Sidecar.from_dict({"accession": "SRR1", "made_up_field": "x"})
    assert sidecar.accession == "SRR1"
    assert not hasattr(sidecar, "made_up_field")


def test_build_sidecar_paired_layout_and_md5(tmp_path):
    acc_dir = tmp_path / "SRR1"
    acc_dir.mkdir()
    mate1 = acc_dir / "SRR1_1.fastq"
    mate2 = acc_dir / "SRR1_2.fastq"
    _write_fastq(mate1, 10)
    _write_fastq(mate2, 10)

    sidecar = build_sidecar("SRR1", acc_dir, {"spots": 10}, "3.0.0", "none")

    assert sidecar.layout == "PAIRED"
    assert sidecar.accession == "SRR1"
    assert sidecar.tool_version == "3.0.0"
    assert sidecar.compression == "none"
    assert len(sidecar.files) == 2

    by_name = {f["name"]: f for f in sidecar.files}
    assert by_name["SRR1_1.fastq"]["md5"] == hashlib.md5(mate1.read_bytes()).hexdigest()
    assert by_name["SRR1_2.fastq"]["md5"] == hashlib.md5(mate2.read_bytes()).hexdigest()
    assert by_name["SRR1_1.fastq"]["bytes"] == mate1.stat().st_size
    assert by_name["SRR1_1.fastq"]["reads"] == 10


def test_build_sidecar_verdict_complete(tmp_path):
    acc_dir = tmp_path / "SRR1"
    acc_dir.mkdir()
    _write_fastq(acc_dir / "SRR1.fastq", 10)

    sidecar = build_sidecar("SRR1", acc_dir, {"spots": 10}, "3.0.0", "none")

    assert sidecar.layout == "SINGLE"
    assert sidecar.state == "complete"
    assert sidecar.completeness["method"] == "spots"
    assert sidecar.completeness["verdict"] == "complete"
    assert sidecar.reads_per_mate == 10


def test_build_sidecar_verdict_partial(tmp_path):
    acc_dir = tmp_path / "SRR1"
    acc_dir.mkdir()
    _write_fastq(acc_dir / "SRR1.fastq", 10)

    sidecar = build_sidecar("SRR1", acc_dir, {"spots": 1000}, "3.0.0", "none")

    assert sidecar.state == "partial"
    assert sidecar.completeness["method"] == "spots"
    assert sidecar.completeness["verdict"] == "truncated"


def test_build_sidecar_verdict_unverified(tmp_path):
    acc_dir = tmp_path / "SRR1"
    acc_dir.mkdir()
    _write_fastq(acc_dir / "SRR1.fastq", 10)

    sidecar = build_sidecar("SRR1", acc_dir, {}, "3.0.0", "none")

    assert sidecar.state == "complete"
    assert sidecar.completeness["method"] == "unverified"
    assert sidecar.completeness["verdict"] == "unverified"


def test_build_sidecar_verify_off_keeps_the_count_but_does_not_judge(tmp_path):
    """With ``verify=False`` the NCBI block keeps ``spots`` and the verdict is ``unverified``."""
    acc_dir = tmp_path / "SRR1"
    acc_dir.mkdir()
    _write_fastq(acc_dir / "SRR1.fastq", 10)

    sidecar = build_sidecar("SRR1", acc_dir, {"spots": 100}, "3.0.0", "none", verify=False)

    assert sidecar.state == "complete"
    assert sidecar.ncbi["spots"] == 100
    assert sidecar.completeness == {"method": "unverified", "ratio": None, "verdict": "unverified"}


def test_build_sidecar_truncated_gzip_is_failed_state(tmp_path):
    acc_dir = tmp_path / "SRR1"
    acc_dir.mkdir()
    (acc_dir / "SRR1.fastq.gz").write_bytes(b"not a real gzip stream")

    sidecar = build_sidecar("SRR1", acc_dir, {"spots": 10}, "3.0.0", "gzip")

    assert sidecar.state == "failed"
    assert sidecar.error is not None
    assert sidecar.completeness["verdict"] == "unverified"
    assert sidecar.reads_per_mate is None


def test_ncbi_from_metadata_xml(tmp_path):
    xml_path = tmp_path / "SRR1.xml"
    xml_path.write_text(NCBI_XML)

    ncbi = ncbi_from_metadata_xml(xml_path)

    assert ncbi["spots"] == 10
    assert ncbi["bases"] == 1400
    assert ncbi["size"] == 4744553813
    assert ncbi["layout"] == "PAIRED"
    assert ncbi["files"] == [{"name": "SRR1", "md5": "abc"}]


def test_ncbi_from_metadata_xml_missing_file_returns_empty(tmp_path, caplog):
    xml_path = tmp_path / "missing.xml"

    with caplog.at_level(logging.WARNING):
        assert ncbi_from_metadata_xml(xml_path) == {}

    assert "Could not read NCBI metadata" in caplog.text
    assert "missing.xml" in caplog.text


def test_ncbi_from_metadata_xml_corrupt_xml_returns_empty(tmp_path, caplog):
    xml_path = tmp_path / "corrupt.xml"
    xml_path.write_text("<EXPERIMENT_PACKAGE_SET><unclosed>")

    with caplog.at_level(logging.WARNING):
        assert ncbi_from_metadata_xml(xml_path) == {}

    assert "Could not read NCBI metadata" in caplog.text
    assert "corrupt.xml" in caplog.text


def _write_corrupt_gzip_fastq(path):
    """Write a gzip FASTQ whose deflate stream has one flipped byte, so reading it raises zlib.error."""
    import gzip as _gzip

    records = "".join(f"@r{i}\n{'ACGT' * 25}\n+\n{'I' * 100}\n" for i in range(200))
    data = bytearray(_gzip.compress(records.encode()))
    data[20] ^= 0xFF
    path.write_bytes(bytes(data))
    return path


def test_build_sidecar_records_a_corrupt_gzip_as_failed(tmp_path):
    """Fix round 1: a corrupt gzip file gives a failed sidecar with the error, not an exception."""
    acc_dir = tmp_path / "SRR1"
    acc_dir.mkdir()
    _write_corrupt_gzip_fastq(acc_dir / "SRR1.fastq.gz")
    sidecar = build_sidecar("SRR1", acc_dir, {"spots": 10}, "3.0.0", "gzip")
    assert sidecar.state == "failed"
    assert "SRR1.fastq.gz" in sidecar.error


def _count_fastq_opens(monkeypatch):
    """Record every FASTQ path opened through ``builtins.open`` (``gzip.open`` goes through it too)."""
    import builtins
    from collections import Counter

    opened = Counter()
    real_open = builtins.open

    def spy(file, *args, **kwargs):
        name = str(file)
        if ".fastq" in name or ".fq" in name:
            opened[os.path.basename(name)] += 1
        return real_open(file, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", spy)
    return opened


@pytest.mark.parametrize("gz", [False, True])
def test_build_sidecar_reads_each_fastq_exactly_once(tmp_path, monkeypatch, gz):
    """Size, md5, read count and the completeness verdict all come from one pass per file."""
    import gzip

    acc_dir = tmp_path / "SRR1"
    acc_dir.mkdir()
    suffix = ".fastq.gz" if gz else ".fastq"
    for name, n_records in (("SRR1_1", 10), ("SRR1_2", 10), ("SRR1", 3)):
        text = "".join(f"@read{i}\nACGT\n+\nIIII\n" for i in range(n_records))
        path = acc_dir / f"{name}{suffix}"
        if gz:
            path.write_bytes(gzip.compress(text.encode()))
        else:
            path.write_text(text)

    opened = _count_fastq_opens(monkeypatch)
    sidecar = build_sidecar("SRR1", acc_dir, {"spots": 13}, "3.0.0", "gzip" if gz else "none")
    monkeypatch.undo()

    assert dict(opened) == {f"SRR1_1{suffix}": 1, f"SRR1_2{suffix}": 1, f"SRR1{suffix}": 1}
    assert sidecar.reads_per_mate == 13
    assert sidecar.state == "complete"
    by_name = {f["name"]: f for f in sidecar.files}
    mate1 = acc_dir / f"SRR1_1{suffix}"
    assert by_name[mate1.name]["md5"] == hashlib.md5(mate1.read_bytes()).hexdigest()
    assert by_name[mate1.name]["bytes"] == mate1.stat().st_size
    assert by_name[mate1.name]["reads"] == 10


def test_build_sidecar_corrupt_gzip_still_records_the_stored_md5(tmp_path):
    """A file that cannot be decompressed keeps the md5 of its stored bytes, as before."""
    acc_dir = tmp_path / "SRR1"
    acc_dir.mkdir()
    bad = acc_dir / "SRR1.fastq.gz"
    bad.write_bytes(b"not a real gzip stream")

    sidecar = build_sidecar("SRR1", acc_dir, {"spots": 10}, "3.0.0", "gzip")

    assert sidecar.state == "failed"
    assert sidecar.files[0]["md5"] == hashlib.md5(bad.read_bytes()).hexdigest()
    assert sidecar.files[0]["reads"] is None


def test_sidecar_refetch_round_trips(tmp_path):
    refetch = {"reason": "truncated", "attempts": 1, "date": "2026-10-01T00:00:00+00:00"}
    sidecar = _make_sidecar(refetch=refetch)
    path = tmp_path / "SRR1.json"
    write_sidecar(path, sidecar)

    assert json.loads(path.read_text())["refetch"] == refetch
    loaded = read_sidecar(path)
    assert loaded == sidecar
    assert loaded.refetch == refetch


def test_sidecar_without_refetch_writes_no_refetch_key():
    sidecar = _make_sidecar()
    assert sidecar.refetch is None
    assert "refetch" not in sidecar.to_dict()


def test_sidecar_sketch_round_trips(tmp_path):
    sketch = {"ksizes": [21, 31, 51], "scaled": 1000, "file": "SRR1.sig.gz", "files": {"SRR1_1.fastq": [4, 1.0]}}
    sidecar = _make_sidecar(sketch=sketch)
    path = tmp_path / "SRR1.json"
    write_sidecar(path, sidecar)

    assert json.loads(path.read_text())["sketch"] == sketch
    assert read_sidecar(path).sketch == sketch


def test_sidecar_without_sketch_writes_no_sketch_key():
    sidecar = _make_sidecar()
    assert sidecar.sketch is None
    assert "sketch" not in sidecar.to_dict()


def test_old_sidecar_json_round_trips_byte_identical(tmp_path):
    """A sidecar written before ``refetch`` existed is rewritten byte for byte as it was read."""
    path = tmp_path / "SRR1.json"
    old = {key: value for key, value in _make_sidecar().to_dict().items() if key != "refetch"}
    original = json.dumps(old, indent=2, sort_keys=True) + "\n"
    path.write_text(original)

    write_sidecar(path, read_sidecar(path))

    assert path.read_text() == original
    assert SIDECAR_SCHEMA == 1


def test_v070_sidecar_fixture_round_trips_byte_identical(tmp_path):
    """A sidecar written by PREPARE 0.7.0 (``tests/fixtures/sidecar_v0.7.0.json``, produced by that
    release's ``write_sidecar``) is read and rewritten byte for byte, with no ``refetch`` key added."""
    fixture = Path(__file__).parent / "fixtures" / "sidecar_v0.7.0.json"
    original = fixture.read_bytes()
    path = tmp_path / "SRR1.json"
    path.write_bytes(original)

    sidecar = read_sidecar(path)
    write_sidecar(path, sidecar)

    assert path.read_bytes() == original
    assert sidecar.state == "partial"
    assert sidecar.reads_per_mate == 40
    assert sidecar.refetch is None
    assert "refetch" not in json.loads(original)
