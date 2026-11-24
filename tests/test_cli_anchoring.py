"""Project files are resolved from the project root, as project folders are (audit P15).

A command run from a subfolder of a project reads and writes the files the pipeline hands from one step
to the next in the project root; a list the user writes by hand stays relative to the working directory.
"""

import pytest

from prepare.cli.main import main
from prepare.core import settings
from prepare.data.registry import REGISTRY_FILENAME, load_registry, save_registry
from prepare.store.resolve import STORE_ENV


@pytest.fixture
def sub(tmp_path, monkeypatch):
    """A project with an empty registry and match CSVs; the working directory is its subfolder ``sub``."""
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    monkeypatch.delenv(STORE_ENV, raising=False)
    settings.reset_for_tests()
    project = tmp_path / "project"
    (project / "matches").mkdir(parents=True)
    (project / "matches" / "G1.csv").write_text("acc,containment\nSRR1,0.9\nSRR2,0.5\n")
    save_registry(load_registry(project / REGISTRY_FILENAME))
    (project / "sub").mkdir()
    monkeypatch.chdir(project / "sub")
    yield project
    settings.reset_for_tests()


def test_parse_containment_writes_its_tables_in_the_project_root(sub):
    assert main(["parse_containment"]) == 0
    assert (sub / "parsed_containment.txt").is_file()
    assert (sub / "top_containments.txt").is_file()
    assert not (sub / "sub" / "parsed_containment.txt").exists()


def test_select_datasets_reads_and_writes_in_the_project_root(sub):
    (sub / "parsed_containment.txt").write_text("accession\tG1\nSRR1\t0.9\nSRR2\t0.5\n")
    assert main(["select_datasets", "--genome-id", "G1", "--threshold", "0.1"]) == 0
    assert (sub / "accessions.txt").read_text().split() == ["SRR1", "SRR2"]
    assert not (sub / "sub" / "accessions.txt").exists()


def test_blacklist_writes_its_file_in_the_project_root(sub):
    assert main(["blacklist", "--add", "SRR9", "--reason", "host reads"]) == 0
    assert "SRR9" in (sub / "blacklist.txt").read_text()
    assert not (sub / "sub" / "blacklist.txt").exists()


def test_sra_profile_writes_its_folder_in_the_project_root(sub):
    (sub / "fastq" / "SRR1").mkdir(parents=True)
    (sub / "fastq" / "SRR1" / "SRR1.fastq").write_text("@r\nACGT\n+\nIIII\n")
    assert main(["sra_profile", "--accession", "SRR1"]) == 0
    assert (sub / "sra_quality_profiles").is_dir()
    assert not (sub / "sub" / "sra_quality_profiles").exists()


def test_an_absolute_path_is_used_unchanged(sub, tmp_path):
    (sub / "parsed_containment.txt").write_text("accession\tG1\nSRR1\t0.9\n")
    target = tmp_path / "elsewhere" / "picked.txt"
    target.parent.mkdir()
    assert main(["select_datasets", "--genome-id", "G1", "--threshold", "0.1", "--output", str(target)]) == 0
    assert target.read_text().split() == ["SRR1"]


def test_an_accessions_file_stays_relative_to_the_working_directory(sub, capsys):
    (sub / "sub" / "list.txt").write_text("SRR1\n")
    assert main(["status", "--accessions-file", "list.txt", "--json"]) == 0
    assert "SRR1" in capsys.readouterr().out
