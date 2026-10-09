"""Tests for default-file resolution shared by CLI commands."""

from pathlib import Path

import pytest

from prepare.core.exceptions import DataAccessError
from prepare.data.defaults import read_matrix, read_table, resolve_metadata_table


class TestResolveMetadataTable:
    def test_explicit_path_wins(self, tmp_path):
        table = tmp_path / "mine.txt"
        table.write_text("Run_ID\tx\n")
        assert resolve_metadata_table(str(table)) == table

    def test_ncbi_table_preferred_over_branchwater(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "metadata_table.txt").write_text("Run_ID\tx\n")
        (tmp_path / "metadata").mkdir()
        (tmp_path / "metadata" / "branchwater_metadata.txt").write_text("Run_ID\tx\n")
        assert resolve_metadata_table(None) == Path("metadata_table.txt")

    def test_falls_back_to_branchwater_metadata(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "metadata").mkdir()
        (tmp_path / "metadata" / "branchwater_metadata.txt").write_text("Run_ID\tx\n")
        assert resolve_metadata_table(None) == Path("metadata/branchwater_metadata.txt")

    def test_missing_lists_candidates(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        with pytest.raises(DataAccessError, match="metadata_table.txt.*branchwater_metadata.txt"):
            resolve_metadata_table(None)

    def test_explicit_missing_raises(self, tmp_path):
        with pytest.raises(DataAccessError, match="nope.txt"):
            resolve_metadata_table(str(tmp_path / "nope.txt"))


class TestReadMatrix:
    def test_tsv_from_parse_containment_drops_text_columns(self, tmp_path):
        f = tmp_path / "parsed_containment.txt"
        f.write_text("\tGCF_A\tGCF_B\tmax_containment\tmax_containment_annotation\nSRR1\t0.9\t0.1\t0.9\tGCF_A\n")
        df = read_matrix(f)
        assert list(df.columns) == ["GCF_A", "GCF_B", "max_containment"]
        assert list(df.index) == ["SRR1"]

    def test_csv(self, tmp_path):
        f = tmp_path / "abundance.csv"
        f.write_text(",sp1,sp2\nS1,1,2\n")
        df = read_matrix(f)
        assert df.loc["S1", "sp2"] == 2


class TestReadTable:
    def test_keeps_text_columns(self, tmp_path):
        f = tmp_path / "metadata_table.txt"
        f.write_text("Run_ID\ttreatment\tabundance\nSRR1\tcontrol\t1.5\nSRR2\ttreated\t2.5\n")
        df = read_table(f)
        assert list(df.columns) == ["treatment", "abundance"]
        assert list(df.index) == ["SRR1", "SRR2"]
        assert df.loc["SRR1", "treatment"] == "control"

    def test_csv(self, tmp_path):
        f = tmp_path / "metadata.csv"
        f.write_text("Run_ID,treatment,abundance\nSRR1,control,1.5\nSRR2,treated,2.5\n")
        df = read_table(f)
        assert list(df.columns) == ["treatment", "abundance"]
        assert list(df.index) == ["SRR1", "SRR2"]
        assert df.loc["SRR2", "treatment"] == "treated"


class TestReadRecords:
    def test_tsv_keeps_first_column(self, tmp_path):
        from prepare.data.defaults import read_records

        f = tmp_path / "taxonomy.tsv"
        f.write_text("genome_id\tgenus\nGCF_1\tWolbachia\n")
        df = read_records(f)
        assert list(df.columns) == ["genome_id", "genus"]
        assert df.loc[0, "genome_id"] == "GCF_1"

    def test_csv(self, tmp_path):
        from prepare.data.defaults import read_records

        f = tmp_path / "v.csv"
        f.write_text("original_name,is_valid\nE. coli,True\n")
        assert list(read_records(f).columns) == ["original_name", "is_valid"]
