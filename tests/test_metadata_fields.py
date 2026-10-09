"""Tests for metadata_fields: the row-to-registry-field mapping and the XML-folder fill-in."""

import logging

import pytest

from prepare.data import registry as reg
from prepare.data import registry_blocks as rb
from prepare.data.metadata_fields import (
    FIELD_COLUMNS,
    fill_metadata_from_xml,
    metadata_fields,
    metadata_fields_from_xml,
)

# A minimal single-run NCBI metadata XML, modeled on tests/test_data_metadata.py's NCBI_XML fixture.
_VALID_XML = """<?xml version="1.0"?>
<EXPERIMENT_PACKAGE_SET>
    <EXPERIMENT_PACKAGE>
        <EXPERIMENT>
            <IDENTIFIERS><PRIMARY_ID>EXP1</PRIMARY_ID></IDENTIFIERS>
            <LIBRARY_DESCRIPTOR>
                <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
                <LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT>
            </LIBRARY_DESCRIPTOR>
            <PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina HiSeq 2500</INSTRUMENT_MODEL></ILLUMINA></PLATFORM>
        </EXPERIMENT>
        <RUN_SET>
            <RUN accession="{accession}" total_spots="{spots}" total_bases="200" size="300">
                <IDENTIFIERS><PRIMARY_ID>{accession}</PRIMARY_ID></IDENTIFIERS>
                <SRAFiles><SRAFile filename="{accession}" md5="abc" semantic_name="run"/></SRAFiles>
            </RUN>
        </RUN_SET>
    </EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>"""


def _write_xml(folder, accession, spots="100"):
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{accession}_metadata.xml"
    path.write_text(_VALID_XML.format(accession=accession, spots=spots))
    return path


class TestMetadataFields:
    """metadata_fields: maps a parsed row (dict or pandas row) to registry field names."""

    def test_maps_known_columns_and_skips_none(self):
        row = {
            "Run_Size": "300",
            "Run_MD5": "abc",
            "Run_Total_Spots": "100",
            "Run_Total_Bases": "200",
            "Experiment_Library_Strategy": "WGS",
            "Sample_Scientific_Name": None,
        }
        fields = metadata_fields(row)
        assert fields["run_size"] == "300"
        assert fields["run_md5"] == "abc"
        assert fields["run_total_spots"] == "100"
        assert fields["run_total_bases"] == "200"
        assert fields["assay_type"] == "WGS"
        assert fields["library_strategy"] == "WGS"
        assert "organism" not in fields

    def test_field_columns_cover_every_record_metadata_key(self):
        field_names = {name for name, _ in FIELD_COLUMNS}
        assert field_names == {
            "run_size",
            "run_md5",
            "run_total_spots",
            "run_total_bases",
            "assay_type",
            "organism",
            "collection_date",
            "library_layout",
            "platform",
            "library_strategy",
        }


class TestMetadataFieldsFromXml:
    """metadata_fields_from_xml: parse_metadata_xml plus the row mapping, {} with a warning on error."""

    def test_parses_a_valid_file(self, tmp_path):
        path = _write_xml(tmp_path, "SRR1")
        fields = metadata_fields_from_xml(path)
        assert fields["run_total_spots"] == "100"
        assert fields["run_total_bases"] == "200"
        assert fields["run_size"] == "300"
        assert fields["run_md5"] == "abc"
        assert fields["library_layout"] == "PAIRED"
        assert fields["platform"] == "ILLUMINA"
        assert fields["library_strategy"] == "WGS"

    def test_missing_file_returns_empty_with_warning(self, tmp_path, caplog):
        with caplog.at_level(logging.WARNING):
            fields = metadata_fields_from_xml(tmp_path / "missing_metadata.xml")
        assert fields == {}
        assert "missing_metadata.xml" in caplog.text

    def test_malformed_xml_returns_empty_with_warning(self, tmp_path, caplog):
        path = tmp_path / "bad_metadata.xml"
        path.write_text("not xml at all <<<")
        with caplog.at_level(logging.WARNING):
            fields = metadata_fields_from_xml(path)
        assert fields == {}
        assert "bad_metadata.xml" in caplog.text

    def test_extraction_error_returns_empty_with_warning(self, tmp_path, caplog, monkeypatch):
        """An XML that parses but fails extraction (``ValueError``) is skipped like the other two cases."""
        path = tmp_path / "odd_metadata.xml"
        path.write_text(_VALID_XML)

        def _failing_parse(_path):
            raise ValueError("no RUN element")

        monkeypatch.setattr("prepare.data.metadata_fields.parse_metadata_xml", _failing_parse)
        with caplog.at_level(logging.WARNING):
            fields = metadata_fields_from_xml(path)
        assert fields == {}
        assert "odd_metadata.xml" in caplog.text
        assert "no RUN element" in caplog.text


class TestFillMetadataFromXml:
    """fill_metadata_from_xml: fills registry metadata blocks missing a spot count from the XML folder."""

    def test_fills_a_block_missing_spots(self, tmp_path):
        registry = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_metadata(registry, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {})
        assert rb.metadata_block(registry, "SRR1").run_total_spots is None
        folder = tmp_path / "metadata"
        _write_xml(folder, "SRR1", spots="555")

        filled = fill_metadata_from_xml(registry, folder)

        assert filled == ["SRR1"]
        block = rb.metadata_block(registry, "SRR1")
        assert block.run_total_spots == 555
        assert block.run_md5 == "abc"
        assert block.platform == "ILLUMINA"

    def test_block_already_holding_spots_is_left_untouched(self, tmp_path):
        registry = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_metadata(registry, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {"run_total_spots": "999"})
        folder = tmp_path / "metadata"
        _write_xml(folder, "SRR1", spots="1")

        filled = fill_metadata_from_xml(registry, folder)

        assert filled == []
        assert rb.metadata_block(registry, "SRR1").run_total_spots == 999

    @pytest.mark.parametrize("spots_attribute", ["", 'total_spots="0"', 'total_spots="unknown"'])
    def test_xml_without_a_positive_spot_count_changes_nothing(self, tmp_path, monkeypatch, spots_attribute):
        registry = reg.load_registry(tmp_path / "prepare_registry.json")
        monkeypatch.setattr(reg, "_now", lambda: "2026-10-01T00:00:00+00:00")
        reg.record_metadata(registry, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {})
        before = dict(registry.datasets["SRR1"]["metadata"])
        # A rewrite of the block would carry this later date.
        monkeypatch.setattr(reg, "_now", lambda: "2026-10-02T00:00:00+00:00")
        folder = tmp_path / "metadata"
        folder.mkdir(parents=True, exist_ok=True)
        xml = _VALID_XML.replace('total_spots="{spots}"', spots_attribute).format(accession="SRR1")
        (folder / "SRR1_metadata.xml").write_text(xml)

        filled = fill_metadata_from_xml(registry, folder)

        assert filled == []
        assert registry.datasets["SRR1"]["metadata"] == before

    def test_accession_without_xml_is_skipped(self, tmp_path):
        registry = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_metadata(registry, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {})
        folder = tmp_path / "metadata"
        folder.mkdir(parents=True, exist_ok=True)

        filled = fill_metadata_from_xml(registry, folder)

        assert filled == []
        assert rb.metadata_block(registry, "SRR1").run_total_spots is None

    def test_unreadable_xml_is_skipped_with_warning(self, tmp_path, caplog):
        registry = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_metadata(registry, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {})
        folder = tmp_path / "metadata"
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "SRR1_metadata.xml").write_text("not xml at all <<<")

        with caplog.at_level(logging.WARNING):
            filled = fill_metadata_from_xml(registry, folder)

        assert filled == []
        assert rb.metadata_block(registry, "SRR1").run_total_spots is None
        assert "SRR1_metadata.xml" in caplog.text

    def test_accession_without_metadata_block_is_skipped(self, tmp_path):
        registry = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_screening(registry, "SRR1", "GCF_1", 0.9, None, "branchwater", 0.01, None)
        folder = tmp_path / "metadata"
        _write_xml(folder, "SRR1")

        filled = fill_metadata_from_xml(registry, folder)

        assert filled == []
        assert rb.metadata_block(registry, "SRR1") is None

    def test_inferred_mark_survives_the_fill(self, tmp_path):
        registry = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_metadata(registry, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {})
        rb.mark_inferred(registry, "SRR1", "metadata")
        assert rb.metadata_block(registry, "SRR1").inferred is True
        folder = tmp_path / "metadata"
        _write_xml(folder, "SRR1", spots="42")

        filled = fill_metadata_from_xml(registry, folder)

        assert filled == ["SRR1"]
        block = rb.metadata_block(registry, "SRR1")
        assert block.run_total_spots == 42
        assert block.inferred is True

    def test_returns_only_filled_accessions_in_registry_order(self, tmp_path):
        registry = reg.load_registry(tmp_path / "prepare_registry.json")
        reg.record_metadata(registry, "SRR1", tmp_path / "metadata" / "SRR1_metadata.xml", {})
        reg.record_metadata(registry, "SRR2", tmp_path / "metadata" / "SRR2_metadata.xml", {"run_total_spots": "1"})
        reg.record_metadata(registry, "SRR3", tmp_path / "metadata" / "SRR3_metadata.xml", {})
        folder = tmp_path / "metadata"
        _write_xml(folder, "SRR1", spots="10")
        _write_xml(folder, "SRR3", spots="30")

        filled = fill_metadata_from_xml(registry, folder)

        assert filled == ["SRR1", "SRR3"]


if __name__ == "__main__":
    pytest.main([__file__])
