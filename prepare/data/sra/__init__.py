"""
SRA data handling for PREPARE.

This package downloads SRA accessions and checks the FASTQ files they produce. It is split by
concern:

- ``fastq``: finding, counting, compressing and verifying FASTQ files already on disk
- ``spots``: an accession's expected spot count and the completeness verdict a read count gives
- ``sampling``: drawing a uniform sample of records from a dataset's FASTQ files
- ``cleanup``: recognising, sizing, preparing and removing temporary download folders
- ``accession``: downloading one accession with prefetch and fasterq-dump
- ``store_handoff``: routing a download through the shared data store
- ``retry``: running many downloads in parallel and retrying the failures
- ``download``: ``download_sra``, the entry point for a whole accession list
- ``run_report``: a run's failure reasons and attempt counts, its report CSV and ``download_run.json``
  (imported from the submodule, not re-exported here)

The public names other packages import are re-exported here. Private names are imported from
their submodule, and each submodule refers to a sibling's names through the sibling module
(``accession_mod.download_accession``), so a test patches a name on the module that owns it.
"""

from prepare.data.sra.fastq import (
    COMPLETE_RATIO_THRESHOLD,
    MATE1_SUFFIXES,
    MATE_SUFFIXES,
    FastqDigest,
    accession_has_fastq,
    compress_fastq,
    count_fastq_reads,
    fastq_digest,
    fastq_files,
    fastq_stem,
    iter_fastq_records,
    orphan_fastq,
    parse_verdict_message,
    primary_fastq,
    verify_download,
)
from prepare.data.sra.spots import expected_spots, merged_verdict, spots_from_xml, verdict_for_count
from prepare.data.sra.sampling import sample_records
from prepare.data.sra.cleanup import is_transient_folder, transient_bytes
from prepare.data.sra.accession import (
    ALREADY_EXISTS,
    STOP,
    classify_download_error,
    download_accession,
    fasterq_dump_version,
)
from prepare.data.sra.store_handoff import SETTLED_PREFIXES, STORE_LINKED_PREFIX, STORE_READY_STATES
from prepare.data.sra.download import default_max_workers, download_sra

__all__ = [
    "ALREADY_EXISTS",
    "COMPLETE_RATIO_THRESHOLD",
    "FastqDigest",
    "MATE1_SUFFIXES",
    "MATE_SUFFIXES",
    "SETTLED_PREFIXES",
    "STOP",
    "STORE_LINKED_PREFIX",
    "STORE_READY_STATES",
    "accession_has_fastq",
    "classify_download_error",
    "compress_fastq",
    "count_fastq_reads",
    "default_max_workers",
    "download_accession",
    "download_sra",
    "expected_spots",
    "fasterq_dump_version",
    "fastq_digest",
    "fastq_files",
    "fastq_stem",
    "iter_fastq_records",
    "is_transient_folder",
    "merged_verdict",
    "orphan_fastq",
    "parse_verdict_message",
    "primary_fastq",
    "sample_records",
    "spots_from_xml",
    "transient_bytes",
    "verdict_for_count",
    "verify_download",
]
