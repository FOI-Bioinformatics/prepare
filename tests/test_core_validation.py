"""Tests for prepare.core.validation.validate_folder."""

import pytest

from prepare.core.exceptions import DataAccessError, ValidationError, exit_code_for
from prepare.core.validation import validate_folder


def test_a_folder_that_cannot_be_created_is_a_data_access_error(tmp_path):
    blocker = tmp_path / "plain_file"
    blocker.write_text("x")

    with pytest.raises(DataAccessError, match="Failed to create directory") as caught:
        validate_folder(blocker / "sub", create=True)

    assert not isinstance(caught.value, ValidationError)
    assert exit_code_for(caught.value) == 1


def test_a_missing_folder_without_create_is_invalid_input(tmp_path):
    with pytest.raises(ValidationError, match="Directory does not exist"):
        validate_folder(tmp_path / "absent")


class TestGenomeName:
    """One rule for the genome id of a FASTA file, shared by every command that names one (P11)."""

    @pytest.mark.parametrize(
        "name, expected",
        [
            ("G.fna", "G"),
            ("G.fna.gz", "G"),
            ("GCF_000008025.1.fasta.gz", "GCF_000008025.1"),
            ("wMel.fa", "wMel"),
            ("wMel.sig", "wMel"),
            ("wMel.sig.gz", "wMel"),
            ("notes.txt", "notes"),
        ],
    )
    def test_genome_name_drops_the_fasta_suffix_and_gz(self, name, expected):
        from pathlib import Path

        from prepare.core.validation import genome_name

        assert genome_name(Path("/x") / name) == expected
        assert genome_name(name) == expected

    @pytest.mark.parametrize("name", [".fna", ".fna.gz", ".sig"])
    def test_a_file_that_is_only_a_suffix_keeps_its_name(self, name):
        """Minor 4: a dot-file named like a suffix is not the empty genome id."""
        from prepare.core.validation import genome_name

        assert genome_name(name) == name
