"""The external tool table: version parsing, probing fake tools on PATH, and ``require_tools``.

Every tool here is a ``tests/helpers_tools.fake_tool`` script on a ``PATH`` that holds nothing
else, so no real tool is ever run.
"""

import pytest

from helpers_tools import fake_tool
from prepare.cli.main import main
from prepare.core import settings
from prepare.core.exceptions import ConfigurationError
from prepare.utils import tools
from prepare.utils.tools import TOOLS, parse_version, probe_tool, require_tools, version_at_least


@pytest.fixture(autouse=True)
def _isolated(tmp_path, monkeypatch):
    """No user config, no store variable, and an empty ``PATH`` unless a test adds fake tools."""
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    monkeypatch.delenv("PREPARE_DATA", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path / "empty-bin"))
    monkeypatch.chdir(tmp_path)
    settings.reset_for_tests()
    yield
    settings.reset_for_tests()


@pytest.mark.parametrize(
    "text, expected",
    [
        ("2.28-r1209", (2, 28)),
        ("samtools 1.21\nUsing htslib 1.21", (1, 21)),
        ("MEGAHIT v1.2.9", (1, 2, 9)),
        ("fasterq-dump : 3.1.1", (3, 1, 1)),
        ("datasets version: 16.27.0", (16, 27, 0)),
        ("pigz 2.8", (2, 8)),
        ("\n\nprefetch : 3.0.10\n", (3, 0, 10)),
        ("no version here", None),
        ("", None),
        ("build 7", None),
    ],
)
def test_parse_version(text, expected):
    assert parse_version(text) == expected


@pytest.mark.parametrize(
    "version, floor, expected",
    [
        ((2, 28), "2.17", True),
        ((2, 16), "2.17", False),
        ((1, 2, 9), "1.2.9", True),
        ((3,), "3.0", True),
        ((3, 0), "3", True),
        ((2, 11, 3), "3.0", False),
        ((1, 10), "1.9", True),
    ],
)
def test_version_at_least_pads_and_compares_numerically(version, floor, expected):
    assert version_at_least(version, floor) is expected


def test_the_table_names_every_tool_with_a_conda_package():
    assert set(TOOLS) >= {"fasterq-dump", "prefetch", "minimap2", "samtools", "megahit", "pigz", "datasets", "seqkit"}
    assert TOOLS["fasterq-dump"].conda_package == "sra-tools"
    assert TOOLS["datasets"].conda_package == "ncbi-datasets-cli"
    assert TOOLS["minimap2"].min_version == "2.17"
    assert TOOLS["samtools"].min_version == "1.10"
    assert TOOLS["megahit"].min_version == "1.2.9"
    assert TOOLS["prefetch"].min_version == TOOLS["fasterq-dump"].min_version == "3.0"
    assert TOOLS["seqkit"].optional and TOOLS["seqkit"].version_args == ("version",)
    for spec in TOOLS.values():
        assert spec.name in spec.install_hint()


def test_probe_reads_the_version_from_stdout(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(fake_tool(tmp_path, "samtools", "samtools 1.21\nUsing htslib 1.21")))
    status = probe_tool("samtools")
    assert status.found and status.path.endswith("/samtools")
    assert status.version == (1, 21)
    assert status.version_text == "samtools 1.21"
    assert status.meets_floor and status.problem() is None


def test_probe_reads_the_version_from_stderr(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(fake_tool(tmp_path, "pigz", "", stderr="pigz 2.8")))
    status = probe_tool("pigz")
    assert status.version == (2, 8)
    assert status.version_text == "pigz 2.8"


def test_probe_of_a_missing_tool_runs_nothing():
    status = probe_tool("minimap2")
    assert not status.found and status.path is None and status.version is None
    assert "minimap2 not found on PATH" in status.problem()
    assert "conda install" in status.problem() and "minimap2" in status.problem()


def test_probe_of_a_tool_that_fails_without_a_version_records_the_error(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(fake_tool(tmp_path, "megahit", "", stderr="boom", rc=2)))
    status = probe_tool("megahit")
    assert status.found and status.version is None
    assert "exited with code 2" in status.error
    # A tool that cannot even print its version is refused, naming the error and a reinstall hint.
    assert "is not runnable, exited with code 2: boom" in status.problem()
    assert "reinstall with: conda install" in status.problem()


def test_a_version_that_cannot_be_read_otherwise_is_not_a_refusal(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(fake_tool(tmp_path, "megahit", "no digits here")))
    status = probe_tool("megahit")
    assert status.found and status.version is None and status.error == "printed no version number"
    assert status.problem() is None


def test_a_tool_that_exits_non_zero_has_no_version_even_if_its_error_names_one(tmp_path, monkeypatch):
    """A broken conda samtools names libcrypto.so.1.0.0; that is not samtools 1.0.0."""
    loader = "samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file"
    monkeypatch.setenv("PATH", str(fake_tool(tmp_path, "samtools", "", stderr=loader, rc=127)))
    status = probe_tool("samtools")
    assert status.found and status.version is None
    assert status.error == f"not runnable, exited with code 127: {loader}"


def test_require_tools_refuses_a_tool_that_cannot_run_with_exit_3(tmp_path, monkeypatch):
    import argparse

    from prepare.cli.commands.read_extraction import ExtractTargetReadsCommand

    loader = "samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file"
    fake_tool(tmp_path, "minimap2", "2.28-r1209")
    monkeypatch.setenv("PATH", str(fake_tool(tmp_path, "samtools", "", stderr=loader, rc=1)))
    with pytest.raises(ConfigurationError, match="not runnable") as excinfo:
        require_tools(["minimap2", "samtools"])
    message = str(excinfo.value)
    assert "missing, too old or not runnable" in message and "libcrypto.so.1.0.0" in message
    assert "older than" not in message
    args = argparse.Namespace(dry_run=False, assemble=False)
    assert ExtractTargetReadsCommand().execute(args) == 3


def test_probe_below_the_floor_is_a_problem_with_the_floor_and_hint(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(fake_tool(tmp_path, "minimap2", "2.16-r922")))
    status = probe_tool("minimap2")
    assert status.version == (2, 16) and not status.meets_floor
    problem = status.problem()
    assert "minimap2 2.16" in problem and "2.17" in problem and "conda install" in problem


def test_require_tools_passes_when_every_tool_is_present_and_new_enough(tmp_path, monkeypatch):
    fake_tool(tmp_path, "minimap2", "2.28-r1209")
    monkeypatch.setenv("PATH", str(fake_tool(tmp_path, "samtools", "samtools 1.21")))
    require_tools(["minimap2", "samtools"])


def test_require_tools_lists_every_problem_at_once(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(fake_tool(tmp_path, "minimap2", "2.16-r922")))
    with pytest.raises(ConfigurationError) as caught:
        require_tools(["minimap2", "samtools", "megahit"])
    message = str(caught.value)
    assert "minimap2 2.16" in message and "2.17" in message
    assert "samtools not found on PATH" in message
    assert "megahit not found on PATH" in message
    assert message.count("conda install") == 3


def test_require_tools_without_version_checks_runs_no_tool(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(fake_tool(tmp_path, "minimap2", "2.16-r922")))
    called = []
    monkeypatch.setattr(tools, "probe_tool", lambda *a, **k: called.append(a))
    require_tools(["minimap2"], check_versions=False)
    assert called == []


def test_require_tools_rejects_a_name_outside_the_table():
    with pytest.raises(KeyError):
        require_tools(["not-a-tool"])


# --- commands refuse a missing or old tool with exit code 3 ---------------------


def _extraction_project(tmp_path):
    table = tmp_path / "parsed_containment.txt"
    table.write_text("\tGCF_1\nSRR1\t0.9\n")
    reads = tmp_path / "fastq" / "SRR1"
    reads.mkdir(parents=True)
    (reads / "SRR1_1.fastq.gz").write_text("x")
    genome = tmp_path / "GCF_1.fna"
    genome.write_text(">s\nACGT\n")
    return table, genome


def test_extract_target_reads_refuses_minimap2_2_16_with_exit_3(tmp_path, monkeypatch, caplog):
    fake_tool(tmp_path, "samtools", "samtools 1.21")
    monkeypatch.setenv("PATH", str(fake_tool(tmp_path, "minimap2", "2.16-r922")))
    table, genome = _extraction_project(tmp_path)
    rc = main(
        [
            "extract_target_reads",
            "--parsed-containment",
            str(table),
            "--genome-id",
            "GCF_1",
            "--genome-fasta",
            str(genome),
            "--fastq-folder",
            str(tmp_path / "fastq"),
            "--output-folder",
            str(tmp_path / "targeted"),
        ]
    )
    assert rc == 3
    assert "minimap2 2.16" in caplog.text and "2.17" in caplog.text
    assert not (tmp_path / "targeted").exists()


def test_download_sra_without_fasterq_dump_exits_3(tmp_path, caplog):
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("SRR1\n")
    rc = main(["download_sra", "--accessions-file", str(accessions), "--fastq-folder", str(tmp_path / "fastq")])
    assert rc == 3
    assert "fasterq-dump not found on PATH" in caplog.text
    assert "sra-tools" in caplog.text


def test_download_sra_dry_run_needs_no_tool(tmp_path):
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("SRR1\n")
    rc = main(
        ["download_sra", "--accessions-file", str(accessions), "--fastq-folder", str(tmp_path / "fastq"), "--dry-run"]
    )
    assert rc == 0


def test_genome_download_without_datasets_exits_3(tmp_path, caplog):
    rc = main(["genome_download", "--accessions", "GCF_000005845.2", "--output-dir", str(tmp_path / "genomes")])
    assert rc == 3
    assert "datasets not found on PATH" in caplog.text
    assert "ncbi-datasets-cli" in caplog.text


def test_genome_prepare_without_datasets_exits_3(tmp_path, caplog):
    accessions = tmp_path / "genomes.txt"
    accessions.write_text("GCF_000005845.2\n")
    rc = main(["genome_prepare", "--accession-file", str(accessions), "--output-dir", str(tmp_path / "genomes")])
    assert rc == 3
    assert "datasets not found on PATH" in caplog.text
