"""
Test data.file_io module functionality.

Tests for file I/O utilities including directory management, file operations, and data processing.
"""

from pathlib import Path
from unittest.mock import patch
import pandas as pd
import pytest

from prepare.data.file_io import (
    ensure_directory,
    is_hidden_name,
    list_files,
    copy_file,
    read_csv,
    visible_files,
    write_csv,
)
from prepare.core.exceptions import DataAccessError


class TestEnsureDirectory:
    """Test ensure_directory function."""

    def test_ensure_directory_success(self, tmp_path):
        """Test successful directory creation."""
        new_dir = tmp_path / "new_directory"

        result = ensure_directory(new_dir)

        assert result == new_dir
        assert new_dir.exists()
        assert new_dir.is_dir()

    def test_ensure_directory_existing(self, tmp_path):
        """Test with existing directory."""
        existing_dir = tmp_path / "existing"
        existing_dir.mkdir()

        result = ensure_directory(existing_dir)

        assert result == existing_dir
        assert existing_dir.exists()

    def test_ensure_directory_nested(self, tmp_path):
        """Test creating nested directories."""
        nested_dir = tmp_path / "level1" / "level2" / "level3"

        result = ensure_directory(nested_dir)

        assert result == nested_dir
        assert nested_dir.exists()
        assert nested_dir.is_dir()

    def test_ensure_directory_string_path(self, tmp_path):
        """Test with string path."""
        new_dir = str(tmp_path / "string_path")

        result = ensure_directory(new_dir)

        assert result == Path(new_dir)
        assert Path(new_dir).exists()

    def test_ensure_directory_failure(self):
        """Test failure in directory creation."""
        with patch("pathlib.Path.mkdir", side_effect=OSError("Permission denied")):
            with pytest.raises(DataAccessError, match="Failed to create directory"):
                ensure_directory("/invalid/path")


class TestListFiles:
    """Test list_files function."""

    def test_list_files_basic(self, tmp_path):
        """Test basic file listing."""
        # Create test files
        (tmp_path / "file1.txt").touch()
        (tmp_path / "file2.txt").touch()
        (tmp_path / "file3.csv").touch()

        result = list_files(tmp_path, "*.txt")

        assert len(result) == 2
        file_names = [f.name for f in result]
        assert "file1.txt" in file_names
        assert "file2.txt" in file_names
        assert "file3.csv" not in file_names

    def test_list_files_all_files(self, tmp_path):
        """Test listing all files with default pattern."""
        # Create test files
        (tmp_path / "file1.txt").touch()
        (tmp_path / "file2.csv").touch()
        (tmp_path / "file3.json").touch()

        result = list_files(tmp_path)

        assert len(result) == 3

    def test_list_files_no_matches(self, tmp_path):
        """Test with no matching files."""
        (tmp_path / "file1.txt").touch()

        result = list_files(tmp_path, "*.csv")

        assert result == []

    def test_list_files_empty_directory(self, tmp_path):
        """Test with empty directory."""
        result = list_files(tmp_path)

        assert result == []

    def test_list_files_nonexistent_directory(self):
        """Test with nonexistent directory."""
        result = list_files("/nonexistent/path")

        assert result == []

    def test_list_files_string_path(self, tmp_path):
        """Test with string path."""
        (tmp_path / "test.txt").touch()

        result = list_files(str(tmp_path), "*.txt")

        assert len(result) == 1
        assert result[0].name == "test.txt"


class TestCopyFile:
    """Test copy_file function."""

    def test_copy_file_success(self, tmp_path):
        """Test successful file copy."""
        source = tmp_path / "source.txt"
        dest = tmp_path / "dest.txt"

        source.write_text("test content")

        result = copy_file(source, dest)

        assert result == dest
        assert dest.exists()
        assert dest.read_text() == "test content"

    def test_copy_file_create_dest_directory(self, tmp_path):
        """Test copy with destination directory creation."""
        source = tmp_path / "source.txt"
        dest = tmp_path / "new_dir" / "dest.txt"

        source.write_text("test content")

        result = copy_file(source, dest)

        assert result == dest
        assert dest.exists()
        assert dest.parent.exists()
        assert dest.read_text() == "test content"

    def test_copy_file_string_paths(self, tmp_path):
        """Test copy with string paths."""
        source = tmp_path / "source.txt"
        dest = tmp_path / "dest.txt"

        source.write_text("test content")

        result = copy_file(str(source), str(dest))

        assert result == dest
        assert dest.exists()

    def test_copy_file_nonexistent_source(self, tmp_path):
        """Test copy with nonexistent source file."""
        source = tmp_path / "nonexistent.txt"
        dest = tmp_path / "dest.txt"

        with pytest.raises(DataAccessError, match="Failed to copy"):
            copy_file(source, dest)

    def test_copy_file_permission_error(self, tmp_path):
        """Test copy with permission error."""
        source = tmp_path / "source.txt"
        dest = tmp_path / "dest.txt"

        source.write_text("test content")

        with patch("shutil.copy2", side_effect=OSError("Permission denied")):
            with pytest.raises(DataAccessError, match="Failed to copy"):
                copy_file(source, dest)


class TestReadCsv:
    """Test read_csv function."""

    def test_read_csv_success(self, tmp_path):
        """Test successful CSV reading."""
        csv_file = tmp_path / "test.csv"
        csv_content = "name,age,city\nAlice,25,NYC\nBob,30,LA"
        csv_file.write_text(csv_content)

        result = read_csv(csv_file)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2
        assert list(result.columns) == ["name", "age", "city"]
        assert result.iloc[0]["name"] == "Alice"

    def test_read_csv_with_kwargs(self, tmp_path):
        """Test CSV reading with additional arguments."""
        csv_file = tmp_path / "test.csv"
        csv_content = "name;age;city\nAlice;25;NYC\nBob;30;LA"
        csv_file.write_text(csv_content)

        result = read_csv(csv_file, sep=";")

        assert len(result) == 2
        assert list(result.columns) == ["name", "age", "city"]

    def test_read_csv_string_path(self, tmp_path):
        """Test CSV reading with string path."""
        csv_file = tmp_path / "test.csv"
        csv_content = "name,age\nAlice,25"
        csv_file.write_text(csv_content)

        result = read_csv(str(csv_file))

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1

    def test_read_csv_nonexistent_file(self):
        """Test CSV reading with nonexistent file."""
        with pytest.raises(DataAccessError, match="Failed to read CSV file"):
            read_csv("/nonexistent/file.csv")

    def test_read_csv_invalid_content(self, tmp_path):
        """Test CSV reading with invalid content."""
        csv_file = tmp_path / "invalid.csv"
        csv_file.write_text("invalid csv content\nwith\nmismatched\ncolumns,too,many")

        # Should still work as pandas is flexible, but test the error path
        with patch("pandas.read_csv", side_effect=Exception("Parse error")):
            with pytest.raises(DataAccessError, match="Failed to read CSV file"):
                read_csv(csv_file)


class TestWriteCsv:
    """Test write_csv function."""

    def test_write_csv_success(self, tmp_path):
        """Test successful CSV writing."""
        df = pd.DataFrame({"name": ["Alice", "Bob"], "age": [25, 30], "city": ["NYC", "LA"]})
        csv_file = tmp_path / "output.csv"

        write_csv(df, csv_file)

        assert csv_file.exists()

        # Read back and verify
        result = pd.read_csv(csv_file, index_col=0)
        assert len(result) == 2
        assert list(result.columns) == ["name", "age", "city"]

    def test_write_csv_create_directory(self, tmp_path):
        """Test CSV writing with directory creation."""
        df = pd.DataFrame({"col1": [1, 2], "col2": [3, 4]})
        csv_file = tmp_path / "new_dir" / "output.csv"

        write_csv(df, csv_file)

        assert csv_file.exists()
        assert csv_file.parent.exists()

    def test_write_csv_with_kwargs(self, tmp_path):
        """Test CSV writing with additional arguments."""
        df = pd.DataFrame({"name": ["Alice"], "age": [25]})
        csv_file = tmp_path / "output.csv"

        write_csv(df, csv_file, sep=";", index=False)

        assert csv_file.exists()
        content = csv_file.read_text()
        assert ";" in content
        assert "0" not in content  # No index

    def test_write_csv_string_path(self, tmp_path):
        """Test CSV writing with string path."""
        df = pd.DataFrame({"col": [1, 2]})
        csv_file = str(tmp_path / "output.csv")

        write_csv(df, csv_file)

        assert Path(csv_file).exists()

    def test_write_csv_permission_error(self, tmp_path):
        """Test CSV writing with permission error."""
        df = pd.DataFrame({"col": [1]})
        csv_file = tmp_path / "output.csv"

        with patch("pandas.DataFrame.to_csv", side_effect=OSError("Permission denied")):
            with pytest.raises(DataAccessError, match="Failed to write CSV file"):
                write_csv(df, csv_file)


class TestFileIoIntegration:
    """Integration tests for file I/O functionality."""

    def test_csv_roundtrip(self, tmp_path):
        """Test complete CSV read/write roundtrip."""
        original_df = pd.DataFrame(
            {"name": ["Alice", "Bob", "Charlie"], "age": [25, 30, 35], "score": [95.5, 87.2, 92.1]}
        )
        csv_file = tmp_path / "roundtrip.csv"

        # Write then read
        write_csv(original_df, csv_file, index=False)
        result_df = read_csv(csv_file)

        # Compare (pandas may change dtypes slightly)
        assert len(result_df) == len(original_df)
        assert list(result_df.columns) == list(original_df.columns)
        assert result_df["name"].tolist() == original_df["name"].tolist()


class TestVisibleFiles:
    def test_hidden_names(self):
        assert is_hidden_name("._SRR1_1.fastq.gz")
        assert is_hidden_name(".DS_Store")
        assert not is_hidden_name("SRR1_1.fastq.gz")

    def test_visible_files_drops_appledouble_and_sorts(self, tmp_path):
        (tmp_path / "b.fna").write_text(">b\nA\n")
        (tmp_path / "a.fna").write_text(">a\nA\n")
        (tmp_path / "._a.fna").write_bytes(b"\x00\x05\x16\x07")
        (tmp_path / ".DS_Store").write_bytes(b"\x00")
        (tmp_path / "sub").mkdir()
        assert [p.name for p in visible_files(tmp_path, "*.fna")] == ["a.fna", "b.fna"]
        assert [p.name for p in visible_files(tmp_path, "*", dirs=True)] == ["sub"]
        assert visible_files(tmp_path / "missing", "*") == []

    def test_visible_files_two_patterns_no_duplicates(self, tmp_path):
        (tmp_path / "x.fastq.gz").write_text("@r\nA\n+\nI\n")
        assert [p.name for p in visible_files(tmp_path, "*.fastq.gz", "*.gz")] == ["x.fastq.gz"]

    def test_list_files_hides_dotfiles_by_default(self, tmp_path):
        (tmp_path / "g.csv").write_text("a\n")
        (tmp_path / "._g.csv").write_bytes(b"\x00")
        assert [p.name for p in list_files(tmp_path, "*.csv")] == ["g.csv"]
        assert len(list_files(tmp_path, "*.csv", include_hidden=True)) == 2


if __name__ == "__main__":
    pytest.main([__file__])
