"""Synthetic data builders for the performance regression tests.

Kept in one module so every fixture used by a timing or memory bound (the tests marked ``perf``) is
easy to find: a 20,000-dataset project registry built through the real ``record_*`` writers
(screening, selection, download, extraction and assembly), a gzip FASTQ writer for the sampler test,
a folder of NCBI efetch metadata XML files for the metadata parser test, and a containment table for
the summary and screening test. Ordinary (non-performance) tests also use some of these functions as
a quick way to write realistic-looking fixture data.
"""

import gzip
import random
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Union
from xml.sax.saxutils import escape

import numpy as np

from prepare.data.registry import (
    Registry,
    record_assembly,
    record_download,
    record_extraction,
    record_screening,
    record_selection,
    save_registry,
)

# ---------------------------------------------------------------------------
# Project registry (build_report, results_rows, stage_counts, write_registry, record_run_outcomes)
# ---------------------------------------------------------------------------


def accession_names(n_datasets: int) -> List[str]:
    """The accessions ``build_registry`` writes, in order."""
    return [f"SRR{index:08d}" for index in range(1, n_datasets + 1)]


def build_registry(
    path: Union[str, Path],
    n_datasets: int = 20000,
    genomes: Sequence[str] = ("G1", "G2", "G3"),
    seed: int = 0,
    n_selected: int = 2000,
    n_downloaded: int = 200,
    n_extracted: int = 5,
) -> Dict[str, List[str]]:
    """Write a registry of ``n_datasets`` screened accessions to ``path`` and return its accession lists.

    Every accession gets one screening entry per genome. The first ``n_selected`` are selected,
    the first ``n_downloaded`` of those are recorded as downloaded with one small FASTQ file each
    under ``<path's folder>/fastq/<acc>/``, and the first ``n_extracted`` are recorded as extracted
    and assembled for the first genome. Returns ``{"all", "selected", "downloaded", "extracted"}``.
    """
    target = Path(path)
    root = target.parent
    rng = random.Random(seed)
    registry = Registry(path=target)
    accessions = accession_names(n_datasets)
    csv_path = root / "matches" / "containment.csv"
    for accession in accessions:
        for genome_id in genomes:
            containment = rng.random()
            record_screening(
                registry, accession, genome_id, containment, containment**0.5, "branchwater", 0.1, csv_path
            )
    selected = accessions[:n_selected]
    record_selection(registry, selected, {"min_containment": 0.1}, root / "selected.txt")
    fastq_dir = root / "fastq"
    downloaded = selected[:n_downloaded]
    for accession in downloaded:
        folder = fastq_dir / accession
        folder.mkdir(parents=True, exist_ok=True)
        (folder / f"{accession}.fastq.gz").write_bytes(b"\x1f\x8b" + bytes(64))
        record_download(registry, accession, "downloaded", fastq_dir, "Downloaded 1 files")
    extracted = downloaded[:n_extracted]
    genome_id = genomes[0]
    for accession in extracted:
        out = root / "targeted" / accession
        record_extraction(
            registry,
            accession,
            genome_id,
            [out / f"{genome_id}_R1.fastq.gz"],
            100,
            False,
            {"preset": "sr", "threshold": 0.1, "genome_fasta": root / "genomes" / f"{genome_id}.fna"},
        )
        record_assembly(
            registry,
            accession,
            genome_id,
            out / f"{genome_id}_assembly",
            {"contigs": 3, "total_bp": 3000, "n50": 1000, "largest": 1200},
            "1.2.9",
            {"threads": 1},
        )
    save_registry(registry, target)
    return {"all": accessions, "selected": selected, "downloaded": downloaded, "extracted": extracted}


# ---------------------------------------------------------------------------
# FASTQ files (sample_records)
# ---------------------------------------------------------------------------


def write_fastq(
    path: Path,
    n: int,
    seq_for: Optional[Callable[[int], str]] = None,
    qual_for: Optional[Callable[[int], str]] = None,
) -> Path:
    """Write ``n`` four-line records to ``path`` (gzip-compressed when it ends in ``.gz``).

    Record ``i`` is named ``r<i>``; its sequence is ``seq_for(i)`` (default ``f"A{i}"``) and
    its quality ``qual_for(i)`` (default ``"I"`` repeated to the sequence length).
    """
    seq_for = seq_for or (lambda i: f"A{i}")
    parts = []
    for i in range(n):
        seq = seq_for(i)
        qual = qual_for(i) if qual_for else "I" * len(seq)
        parts.append(f"@r{i}\n{seq}\n+\n{qual}\n")
    data = "".join(parts).encode("ascii")
    if str(path).endswith(".gz"):
        with gzip.open(path, "wb", compresslevel=1) as handle:
            handle.write(data)
    else:
        path.write_bytes(data)
    return path


def write_illumina_like_fastq_gz(path: Path, n: int, read_length: int = 150) -> Path:
    """Write ``n`` gzip-compressed records of ``read_length`` bases with varied bases and qualities.

    Sequences and qualities are shifted windows over fixed patterns, which is quick to
    generate; the content does not matter for timing.
    """
    bases = "ACGTTGCAAGCTTCGA" * (read_length // 16 + 2)
    quals = "".join(chr(33 + (i * 7) % 41) for i in range(read_length + 16))
    parts = []
    for i in range(n):
        offset = i % 16
        parts.append(f"@SRR0.{i} {i} length={read_length}\n")
        parts.append(bases[offset : offset + read_length] + "\n+\n")
        parts.append(quals[offset : offset + read_length] + "\n")
    with gzip.open(path, "wb", compresslevel=1) as handle:
        handle.write("".join(parts).encode("ascii"))
    return path


# ---------------------------------------------------------------------------
# NCBI efetch metadata XML files (parse_metadata)
# ---------------------------------------------------------------------------

# Each generated file follows the layout of a real single-run ``EXPERIMENT_PACKAGE_SET`` response
# (experiment, submission, organization, study, sample with SAMPLE_ATTRIBUTES, pool, run with
# SRAFiles and read statistics), at roughly the size of a real file. Attribute sets overlap between
# files: every file carries a shared core of tags plus a window of tags that shifts with the file
# index, so the folder as a whole holds many more distinct tags than any one file.

CORE_TAGS = [
    "collection_date",
    "geo_loc_name",
    "lat_lon",
    "env_broad_scale",
    "env_local_scale",
    "env_medium",
    "host",
    "isolation_source",
    "sample_type",
    "strain",
]

_METADATA_TEMPLATE = """<?xml version="1.0" encoding="UTF-8" ?>
<EXPERIMENT_PACKAGE_SET>
<EXPERIMENT_PACKAGE>
<EXPERIMENT accession="SRX{n}" alias="exp_{n}">
<IDENTIFIERS><PRIMARY_ID>SRX{n}</PRIMARY_ID></IDENTIFIERS>
<TITLE>Illumina sequencing of sample {n}</TITLE>
<STUDY_REF accession="SRP{p}"><IDENTIFIERS><PRIMARY_ID>SRP{p}</PRIMARY_ID></IDENTIFIERS></STUDY_REF>
<DESIGN>
<DESIGN_DESCRIPTION>Shotgun metagenome of sample {n}</DESIGN_DESCRIPTION>
<SAMPLE_DESCRIPTOR accession="SRS{n}"><IDENTIFIERS><PRIMARY_ID>SRS{n}</PRIMARY_ID></IDENTIFIERS></SAMPLE_DESCRIPTOR>
<LIBRARY_DESCRIPTOR>
<LIBRARY_NAME>lib_{n}</LIBRARY_NAME>
<LIBRARY_STRATEGY>{strategy}</LIBRARY_STRATEGY>
<LIBRARY_SOURCE>METAGENOMIC</LIBRARY_SOURCE>
<LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION>
<LIBRARY_LAYOUT>{layout}</LIBRARY_LAYOUT>
</LIBRARY_DESCRIPTOR>
</DESIGN>
<PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina NovaSeq 6000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM>
</EXPERIMENT>
<SUBMISSION accession="SRA{p}" lab_name="Lab {p}">\
<IDENTIFIERS><PRIMARY_ID>SRA{p}</PRIMARY_ID></IDENTIFIERS></SUBMISSION>
<Organization type="institute"><Name>Institute {p}</Name><Contact email="lab{p}@example.org"><Name><First>A</First>\
<Last>B</Last></Name></Contact></Organization>
<STUDY accession="SRP{p}" alias="study_{p}">
<IDENTIFIERS><PRIMARY_ID>SRP{p}</PRIMARY_ID><EXTERNAL_ID namespace="BioProject">PRJNA{p}</EXTERNAL_ID></IDENTIFIERS>
<DESCRIPTOR>
<STUDY_TITLE>Metagenomes of project {p}</STUDY_TITLE>
<STUDY_TYPE existing_study_type="Metagenomics"/>
<STUDY_ABSTRACT>{abstract}</STUDY_ABSTRACT>
</DESCRIPTOR>
</STUDY>
<SAMPLE accession="SRS{n}" alias="sample_{n}">
<IDENTIFIERS><PRIMARY_ID>SRS{n}</PRIMARY_ID><EXTERNAL_ID namespace="BioSample">SAMN{n}</EXTERNAL_ID></IDENTIFIERS>
<TITLE>Sample {n}</TITLE>
<SAMPLE_NAME><TAXON_ID>{taxon}</TAXON_ID><SCIENTIFIC_NAME>{organism}</SCIENTIFIC_NAME></SAMPLE_NAME>
<SAMPLE_LINKS><SAMPLE_LINK><XREF_LINK><DB>bioproject</DB><ID>{p}</ID></XREF_LINK></SAMPLE_LINK></SAMPLE_LINKS>
<SAMPLE_ATTRIBUTES>
{attributes}
</SAMPLE_ATTRIBUTES>
</SAMPLE>
<Pool><Member member_name="" accession="SRS{n}" sample_name="sample_{n}" spots="{spots}" bases="{bases}"\
 tax_id="{taxon}"\
 organism="{organism}"><IDENTIFIERS><PRIMARY_ID>SRS{n}</PRIMARY_ID></IDENTIFIERS></Member></Pool>
<RUN_SET runs="1" bases="{bases}" spots="{spots}" bytes="{size}">
<RUN accession="SRR{n}" alias="run_{n}" total_spots="{spots}" total_bases="{bases}" size="{size}" load_done="true"\
 published="2024-01-01" is_public="true" cluster_name="public" static_data_available="1">
<IDENTIFIERS><PRIMARY_ID>SRR{n}</PRIMARY_ID></IDENTIFIERS>
<EXPERIMENT_REF accession="SRX{n}"/>
<Pool><Member member_name="" accession="SRS{n}" sample_name="sample_{n}" spots="{spots}" bases="{bases}"\
 tax_id="{taxon}" organism="{organism}"><IDENTIFIERS><PRIMARY_ID>SRS{n}</PRIMARY_ID></IDENTIFIERS></Member></Pool>
<SRAFiles>
<SRAFile cluster="public" filename="sample_{n}_R1.fastq.gz" url="https://sra-pub-src-1.s3.amazonaws.com/SRR{n}/r1"\
 size="{half}" date="2024-01-01" md5="{md5a}" semantic_name="fastq" supertype="Original" sratoolkit="0"/>
<SRAFile cluster="public" filename="SRR{n}"\
 url="https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/SRR{n}/SRR{n}.lite.1"\
 size="{size}" date="2024-01-02" md5="{md5b}" version="1" semantic_name="run" supertype="Primary ETL" sratoolkit="1"/>
</SRAFiles>
<CloudFiles><CloudFile filetype="run" provider="gs" location="gs.us-east1"/>\
<CloudFile filetype="run" provider="s3" location="s3.us-east-1"/></CloudFiles>
<Statistics nreads="2" nspots="{spots}"><Read index="0" count="{spots}" average="151" stdev="0"/>\
<Read index="1" count="{spots}" average="151" stdev="0"/></Statistics>
<Bases cs_native="false" count="{bases}"><Base value="A" count="1"/><Base value="C" count="1"/>\
<Base value="G" count="1"/><Base value="T" count="1"/><Base value="N" count="0"/></Bases>
</RUN>
</RUN_SET>
</EXPERIMENT_PACKAGE>
</EXPERIMENT_PACKAGE_SET>
"""


def attribute_tags(index: int, per_file: int = 40, pool: int = 1000) -> List[str]:
    """The SAMPLE_ATTRIBUTE tags of file ``index``: the shared core plus a shifting window of the pool."""
    window = per_file - len(CORE_TAGS)
    start = (index * 7) % pool
    return CORE_TAGS + [f"attr_{(start + k) % pool:03d}" for k in range(window)]


def metadata_xml(index: int, per_file: int = 40, pool: int = 1000) -> str:
    """One synthetic single-run metadata XML document with ``per_file`` sample attributes."""
    n = 1000000 + index
    spots = 1000000 + index * 13
    bases = spots * 302
    size = bases // 3
    attributes = "\n".join(
        f"<SAMPLE_ATTRIBUTE><TAG>{tag}</TAG><VALUE>{escape(f'value {index} of {tag}')}</VALUE></SAMPLE_ATTRIBUTE>"
        for tag in attribute_tags(index, per_file, pool)
    )
    return _METADATA_TEMPLATE.format(
        n=n,
        p=500000 + index // 25,
        strategy="WGS" if index % 5 else "AMPLICON",
        layout="<PAIRED/>" if index % 3 else "<SINGLE/>",
        abstract=escape("Shotgun sequencing of environmental samples. " * 20),
        taxon=256318 + index % 4,
        organism="metagenome" if index % 4 else "soil metagenome",
        attributes=attributes,
        spots=spots,
        bases=bases,
        size=size,
        half=size // 2,
        md5a=f"{index:032x}",
        md5b=f"{index + 1:032x}",
    )


def write_metadata_folder(folder: Path, count: int = 300, per_file: int = 40, pool: int = 1000) -> List[Path]:
    """Write ``count`` synthetic metadata XML files named ``SRR<n>_metadata.xml`` into ``folder``."""
    folder.mkdir(parents=True, exist_ok=True)
    paths = []
    for index in range(count):
        path = folder / f"SRR{1000000 + index}_metadata.xml"
        path.write_text(metadata_xml(index, per_file, pool))
        paths.append(path)
    return paths


# ---------------------------------------------------------------------------
# Containment table (parse_containment's summary and screening step)
# ---------------------------------------------------------------------------


def containment_data(rows: int, genomes: int, seed: int = 7, zero_column: bool = True) -> Dict[str, Dict[str, float]]:
    """``rows`` accessions by ``genomes`` genomes of containment values, with ties and zeros.

    Accessions are ``SRR<n>`` and genomes ``GCF_<n>``. With ``zero_column`` the last genome is 0 for
    every accession. The same arguments always give the same values.
    """
    rng = np.random.default_rng(seed)
    levels = np.round(np.linspace(0.05, 1.0, 20), 2)
    values = rng.choice(levels, size=(rows, genomes))
    values[rng.random((rows, genomes)) < 0.5] = 0.0
    if zero_column:
        values[:, -1] = 0.0
    names = [f"GCF_{j:03d}" for j in range(genomes)]
    return {f"SRR{1000000 + i}": {names[j]: float(values[i, j]) for j in range(genomes)} for i in range(rows)}
