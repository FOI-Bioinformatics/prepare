"""Opt-in smoke tests that talk to NCBI SRA (run with: make test-network).

SRR2517620 is a 425-spot MiSeq mosquito metagenome; fasterq-dump fetches it in
a few seconds, which is enough to prove the download command line and the
output layout end to end.
"""

import gzip
import shutil
import subprocess
from pathlib import Path

import pytest

from prepare.data.sra import download_accession

pytestmark = pytest.mark.network

TINY_RUN = "SRR2517620"
TINY_RUN_SPOTS = 425

REPO_ROOT = Path(__file__).resolve().parents[1]
SMOKE_CHAIN_SCRIPT = REPO_ROOT / "scripts" / "smoke_chain.sh"

needs_fasterq_dump = pytest.mark.skipif(shutil.which("fasterq-dump") is None, reason="fasterq-dump not on PATH")
needs_cli = pytest.mark.skipif(shutil.which("prepare") is None, reason="prepare CLI not on PATH")
needs_minimap2 = pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 not on PATH")
needs_samtools = pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")


def _read_count(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return sum(1 for _ in handle) // 4


@needs_fasterq_dump
def test_download_accession_writes_paired_fastq(tmp_path, monkeypatch):
    """Default settings (prefetch when available, compression on): files land as .fastq.gz."""
    monkeypatch.chdir(tmp_path)
    ok, message = download_accession(TINY_RUN, tmp_path / "fastq", num_threads=2)
    assert ok, message

    files = sorted(p.name for p in (tmp_path / "fastq" / TINY_RUN).glob("*.fastq*"))
    assert files == [f"{TINY_RUN}_1.fastq.gz", f"{TINY_RUN}_2.fastq.gz"]
    assert _read_count(tmp_path / "fastq" / TINY_RUN / f"{TINY_RUN}_1.fastq.gz") == TINY_RUN_SPOTS


@needs_fasterq_dump
def test_download_accession_writes_paired_fastq_uncompressed(tmp_path, monkeypatch):
    """compress=False leaves the plain FASTQ files fasterq-dump wrote, uncompressed."""
    monkeypatch.chdir(tmp_path)
    ok, message = download_accession(TINY_RUN, tmp_path / "fastq", num_threads=2, compress=False)
    assert ok, message

    files = sorted(p.name for p in (tmp_path / "fastq" / TINY_RUN).glob("*.fastq*"))
    assert files == [f"{TINY_RUN}_1.fastq", f"{TINY_RUN}_2.fastq"]
    assert _read_count(tmp_path / "fastq" / TINY_RUN / f"{TINY_RUN}_1.fastq") == TINY_RUN_SPOTS


@needs_fasterq_dump
@needs_cli
def test_download_sra_cli_round_trip(tmp_path):
    (tmp_path / "accessions.txt").write_text(f"{TINY_RUN}\n")
    result = subprocess.run(
        ["prepare", "download_sra", "--accessions-file", "accessions.txt", "--fastq-folder", "fastq"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert result.returncode == 0, result.stderr
    # Compression defaults to on, so the CLI writes the gzipped file, not the plain one.
    assert (tmp_path / "fastq" / TINY_RUN / f"{TINY_RUN}_1.fastq.gz").exists()
    assert "Newly downloaded: 1" in result.stderr + result.stdout


@needs_fasterq_dump
@needs_minimap2
@needs_samtools
def test_smoke_chain_script(tmp_path):
    """scripts/smoke_chain.sh end to end: the same chain the nightly workflow runs
    (download_test_genome, download_sra, sra_validate, sra_profile,
    extract_target_reads, results_table) against the real SRR2517620 run, so
    `make test-network` exercises exactly what the nightly job does.

    A mosquito metagenome maps almost nothing onto the bacterial reference genome
    used here, so this only checks that every step exits 0 and results.tsv gets a
    header plus one data row, not how many reads mapped.
    """
    result = subprocess.run(
        ["bash", str(SMOKE_CHAIN_SCRIPT), str(tmp_path)],
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    results_path = tmp_path / "results.tsv"
    assert results_path.exists(), result.stdout + result.stderr
    lines = results_path.read_text().splitlines()
    assert len(lines) == 2, f"expected a header plus one data row, got: {lines}"
    assert lines[0].split("\t")[:2] == ["accession", "genome_id"]
    assert lines[1].startswith(f"{TINY_RUN}\t")
