"""A small project on disk for the ``project_report`` tests: a registry file, FASTQ files and a run log.

Every stage the report covers holds something: screening against two genomes, a selection, an
exclusion, downloads with each completeness verdict, a failed download with a retried network
message, extractions with coverage, an assembly, timings, a profile analysis, a results_table export
and three run-log lines. The FASTQ files exist so a test can check that the report never opens them.
"""

import json
from pathlib import Path

from prepare.data import run_log
from prepare.data.registry import (
    REGISTRY_FILENAME,
    Registry,
    record_analysis,
    record_assembly,
    record_download,
    record_exclusion,
    record_export,
    record_extraction,
    record_screening,
    record_selection,
    save_registry,
    set_download_verdict,
)
from prepare.data.registry_timing import set_assembly_timing, set_download_timing, set_extraction_timing
from prepare.data.run_log import RunRecord

FAILED_MESSAGE = "Retry 2: Connection reset by peer"


def _fastq(root: Path, accession: str) -> None:
    folder = root / "fastq" / accession
    folder.mkdir(parents=True, exist_ok=True)
    for mate in ("1", "2"):
        (folder / f"{accession}_{mate}.fastq").write_text("@r\nACGT\n+\nIIII\n")


def _runs(root: Path) -> None:
    folder = run_log.runs_dir(root)
    folder.mkdir(parents=True, exist_ok=True)
    lines = []
    for n, command in enumerate(("download_sra", "extract_target_reads", "results_table"), start=1):
        record = RunRecord(
            run_id=f"20261001T12000{n}Z-{command}-000{n}",
            command=command,
            started=f"2026-10-01T12:00:0{n}Z",
            finished=f"2026-10-01T12:00:0{n + 1}Z",
            seconds=float(n),
            exit_code=0,
            argv=[command],
            args={},
            version="0.8.0",
            host="host",
            pid=100 + n,
            summary={"step": n},
        )
        lines.append(json.dumps(record.to_dict()))
    (folder / run_log.RUNS_FILE).write_text("\n".join(lines) + "\n")


def build_project(root: Path, with_runs: bool = True) -> Path:
    """Write the project under ``root`` and return the registry file's path."""
    registry = Registry(path=root / REGISTRY_FILENAME)
    for n in range(1, 7):
        record_screening(registry, f"SRR{n}", "GCF_A", 0.1 * n, None, "matches", 0.0, None)
    record_screening(registry, "SRR2", "GCF_B", 0.3, None, "matches", 0.0, None)
    record_selection(registry, ["SRR1", "SRR2", "SRR3", "SRR4", "SRR5"], {"min_containment": 0.1}, root / "sel.txt")
    record_exclusion(registry, "SRR5", "amplicon | mislabelled")
    for accession in ("SRR1", "SRR2", "SRR3"):
        _fastq(root, accession)
        record_download(registry, accession, "downloaded", root / "fastq", attempt=True)
    set_download_verdict(registry, "SRR1", {"method": "spots", "verdict": "complete"})
    set_download_verdict(registry, "SRR2", {"method": "spots", "verdict": "truncated", "ratio": 0.5})
    set_download_verdict(registry, "SRR3", {"method": "spots", "verdict": "unverified"})
    record_download(registry, "SRR4", "failed", root / "fastq", message=FAILED_MESSAGE)
    record_download(registry, "SRR4", "failed", root / "fastq", message=FAILED_MESSAGE)
    set_download_timing(registry, "SRR1", "2026-10-01T10:00:00+00:00", 10.0)
    set_download_timing(registry, "SRR2", "2026-10-01T10:00:10+00:00", 20.0)
    set_download_timing(registry, "SRR3", "2026-10-01T10:00:30+00:00", 30.0)
    set_download_timing(registry, "SRR4", "2026-10-01T10:01:00+00:00", 40.0)
    coverage = {"breadth": 0.9, "mean_depth": 12.0, "coverage_tsv": root / "targeted" / "cov.tsv"}
    record_extraction(registry, "SRR1", "GCF_A", [], 500, False, {}, coverage=coverage)
    record_extraction(registry, "SRR2", "GCF_A", [], 0, False, {})
    record_extraction(registry, "SRR2", "GCF_B", [], 40, False, {}, coverage={"breadth": 0.5, "mean_depth": 2.0})
    record_extraction(registry, "SRR3", "GCF_A", [], 70, False, {}, coverage={"breadth": 0.7, "mean_depth": 4.0})
    set_extraction_timing(registry, "SRR1", "GCF_A", "2026-10-01T10:02:00+00:00", 5.0)
    record_assembly(registry, "SRR1", "GCF_A", root / "targeted" / "asm", {"contigs": 4, "total_bp": 40000}, "v1", {})
    set_assembly_timing(registry, "SRR1", "GCF_A", "2026-10-01T10:03:00+00:00", 60.0)
    record_analysis(registry, "SRR1", "profile", root / "profiles" / "SRR1.json", {"gc_percent": 41.0})
    record_analysis(registry, "SRR2", "profile", root / "profiles" / "SRR2.json", {"gc_percent": 43.0})
    record_export(registry, "results_table", root / "results.tsv", {"rows": 4})
    save_registry(registry)
    if with_runs:
        _runs(root)
    return root / REGISTRY_FILENAME
