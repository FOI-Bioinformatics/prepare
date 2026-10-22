"""Disk scanners and the registry bootstrap: what is on disk, and a registry rebuilt from it.

The scanners list downloads, metadata files, extracted reads and assembly folders under a
project. ``bootstrap_from_disk`` turns that listing into a registry whose blocks are all marked
as inferred. This module is imported by ``prepare.data.registry``, which re-exports its public
names, so it reaches the registry's record functions through a call-time lookup of that module.
"""

import logging
from pathlib import Path
from types import ModuleType
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Sequence, Set, Tuple, Union

from prepare.data import registry_blocks as rb
from prepare.data.file_io import visible_files
from prepare.data.read_extraction import coverage_table_path, summarise_contigs, summarise_coverage_table
from prepare.data.sra import accession_has_fastq, count_fastq_reads, fastq_files, is_transient_folder

if TYPE_CHECKING:
    from prepare.data.registry import ProjectPaths, Registry

logger = logging.getLogger(__name__)

_MATE_SUFFIXES = ("_1", "_2", "_s", "_0")
_ASSEMBLY_SUFFIX = "_assembly"
_CONTIGS_NAME = "final.contigs.fa"


def _registry_module() -> ModuleType:
    """The ``prepare.data.registry`` module, looked up at call time because that module imports this one."""
    from prepare.data import registry

    return registry


def split_extract_filename(name: str, genome_ids: Sequence[str]) -> Optional[Tuple[str, str]]:
    """Split ``<genome><mate>.fastq(.gz)`` into (genome_id, mate suffix); known genome ids win, longest first."""
    stem = name
    for suffix in (".fastq.gz", ".fastq", ".fq.gz", ".fq"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    else:
        return None
    for genome_id in sorted(genome_ids, key=len, reverse=True):
        if stem == genome_id:
            return genome_id, ""
        for mate in _MATE_SUFFIXES:
            if stem == genome_id + mate:
                return genome_id, mate
    for mate in _MATE_SUFFIXES:
        if stem.endswith(mate):
            return stem[: -len(mate)], mate
    return stem, ""


def scan_downloads(fastq_folder: Path) -> Dict[str, Tuple[int, int]]:
    """Accession -> (number of FASTQ files, total bytes) for every per-accession folder with reads.

    Transient folders (``<acc>_temp``, ``.sra-cache``) are skipped even when they hold a
    partial FASTQ file, so an in-progress or failed download is never mistaken for a
    downloaded accession.
    """
    found: Dict[str, Tuple[int, int]] = {}
    for folder in visible_files(fastq_folder, dirs=True):
        if is_transient_folder(folder.name):
            continue
        if accession_has_fastq(folder):
            # fastq_files, not a raw glob: a zero-byte file or the hidden temporary leftover of
            # an interrupted compression is not a downloaded read file.
            files = fastq_files(folder)
            found[folder.name] = (len(files), sum(p.stat().st_size for p in files))
    return found


def scan_metadata(metadata_folder: Path) -> Set[str]:
    """Return the accessions with an ``<accession>_metadata.xml`` file directly under ``metadata_folder``."""
    return {p.name[: -len("_metadata.xml")] for p in visible_files(metadata_folder, "*_metadata.xml")}


def _genome_ids_on_disk(paths: "ProjectPaths", registry: Optional["Registry"]) -> Set[str]:
    """The project's genome ids: recorded ones, match CSVs and extraction folders.

    A FASTA under ``genomes/`` alone is a file, not a project genome: it is listed only when
    its id also appears in the matches folder, the targeted folder or a registry record.
    """
    _r = _registry_module()
    ids: Set[str] = set(_r.known_genome_ids(registry)) if registry else set()
    ids.update(p.stem for p in visible_files(paths.matches, "*.csv"))
    for acc_dir in visible_files(paths.targeted, dirs=True):
        for asm in visible_files(acc_dir, f"*{_ASSEMBLY_SUFFIX}", dirs=True):
            ids.add(asm.name[: -len(_ASSEMBLY_SUFFIX)])
    return ids


def scan_extractions(targeted_folder: Path, genome_ids: Sequence[str]) -> Dict[str, Dict[str, List[Path]]]:
    """Return the extracted read files found on disk, as accession -> genome id -> list of file paths.

    Walks each per-accession folder under ``targeted_folder`` and classifies its files with
    ``split_extract_filename`` against ``genome_ids``; files that do not match a known genome id
    or mate suffix pattern are skipped. Returns an empty mapping if ``targeted_folder`` does not
    exist as a directory.
    """
    found: Dict[str, Dict[str, List[Path]]] = {}
    if not targeted_folder.is_dir():
        return found
    for acc_dir in visible_files(targeted_folder, dirs=True):
        for path in visible_files(acc_dir):
            split = split_extract_filename(path.name, genome_ids)
            if split:
                found.setdefault(acc_dir.name, {}).setdefault(split[0], []).append(path)
    return found


def scan_assemblies(targeted_folder: Path, genome_ids: Sequence[str]) -> Dict[str, Dict[str, Path]]:
    """Return the assembly directories found on disk, as accession -> genome id -> assembly directory path.

    Walks each per-accession folder under ``targeted_folder`` for subfolders named
    ``<genome_id>_assembly``; ``genome_ids`` is accepted for symmetry with ``scan_extractions``
    but not otherwise used, since the assembly suffix alone identifies the genome id. Returns an
    empty mapping if ``targeted_folder`` does not exist as a directory.
    """
    found: Dict[str, Dict[str, Path]] = {}
    if not targeted_folder.is_dir():
        return found
    for acc_dir in visible_files(targeted_folder, dirs=True):
        for asm in visible_files(acc_dir, f"*{_ASSEMBLY_SUFFIX}", dirs=True):
            found.setdefault(acc_dir.name, {})[asm.name[: -len(_ASSEMBLY_SUFFIX)]] = asm
    return found


def empty_assembly_dirs(targeted_folder: Path, genome_ids: Sequence[str]) -> List[Tuple[str, str]]:
    """(accession, genome_id) pairs whose assembly directory holds zero contigs, sorted."""
    pairs: List[Tuple[str, str]] = []
    for acc, per_genome_asm in scan_assemblies(targeted_folder, genome_ids).items():
        for genome_id, asm_dir in per_genome_asm.items():
            if summarise_contigs(asm_dir / _CONTIGS_NAME)["contigs"] == 0:
                pairs.append((acc, genome_id))
    return sorted(pairs)


def _read_accession_list(path: Optional[Union[str, Path]]) -> List[str]:
    if not path or not Path(path).exists():
        return []
    return [ln.strip() for ln in Path(path).read_text().splitlines() if ln.strip() and not ln.startswith("#")]


def _screening_from_matches(registry: "Registry", matches_folder: Path) -> None:
    _r = _registry_module()
    import csv

    for csv_path in visible_files(matches_folder, "*.csv"):
        with open(csv_path, newline="") as handle:
            for row in csv.DictReader(handle):
                acc = (row.get("acc") or row.get("SRA accession") or "").strip()
                if not acc:
                    continue
                try:
                    containment = float(row.get("containment", ""))
                except ValueError:
                    continue
                cani_raw = row.get("cANI") or ""
                try:
                    cani = float(cani_raw) if cani_raw.strip() else None
                except ValueError:
                    cani = None
                _r.record_screening(registry, acc, csv_path.stem, containment, cani, "matches", 0.0, csv_path)
                rb.mark_inferred(registry, acc, "screening")


def _bootstrap_selection(
    registry: "Registry", accessions_file: Optional[Union[str, Path]], parsed_containment: Optional[Union[str, Path]]
) -> None:
    _r = _registry_module()
    wanted = _read_accession_list(accessions_file)
    if parsed_containment and Path(parsed_containment).exists():
        with open(parsed_containment) as handle:
            next(handle, None)
            wanted += [ln.split("\t", 1)[0].strip() for ln in handle if ln.strip()]
    if wanted:
        _r.record_selection(
            registry, list(dict.fromkeys(wanted)), {"source": "bootstrap"}, accessions_file or parsed_containment or ""
        )
        for acc in wanted:
            rb.mark_inferred(registry, acc, "selection")


def _bootstrap_downloads_and_metadata(registry: "Registry", paths: "ProjectPaths") -> None:
    _r = _registry_module()
    for acc in scan_downloads(paths.fastq):
        _r.record_download(registry, acc, "downloaded", paths.fastq)
        rb.mark_inferred(registry, acc, "download", attempts=0)
    root = _r.project_root(registry)
    for acc in sorted(scan_metadata(paths.metadata)):
        _r.record_metadata(registry, acc, paths.metadata / f"{acc}_metadata.xml", {}, root=root)
        rb.mark_inferred(registry, acc, "metadata")


def _infer_extraction(registry: "Registry", acc: str, genome_id: str, files: Sequence[Path]) -> None:
    """Record one (accession, genome) extraction found on disk, marked as inferred.

    Reads are counted in every file of the pair (both mates, singles and unpaired), so the
    inferred count approximates the number of BAM records a real extraction records. When the
    extraction's ``<genome>_coverage.tsv`` is on disk, breadth and mean depth are read from it;
    otherwise (or when it cannot be read) they are None.
    """
    _r = _registry_module()
    reads = sum(count_fastq_reads(f) for f in files)
    coverage = _infer_coverage(coverage_table_path(files[0].parent, genome_id)) if files else None
    _r.record_extraction(registry, acc, genome_id, files, reads, False, {}, coverage=coverage)
    rb.mark_inferred(registry, acc, "extractions", genome_id=genome_id)


def _infer_coverage(tsv: Path) -> Optional[Dict[str, Any]]:
    """Breadth, mean depth and path of an existing ``samtools coverage`` table, else None."""
    if not tsv.is_file():
        return None
    try:
        summary = summarise_coverage_table(tsv)
    except (OSError, ValueError, KeyError) as e:
        logger.warning("Cannot read coverage table %s (%s); recording no coverage for it", tsv, e)
        return None
    return {"breadth": summary["breadth"], "mean_depth": summary["mean_depth"], "coverage_tsv": tsv}


def _infer_assembly(registry: "Registry", acc: str, genome_id: str, asm_dir: Path) -> None:
    """Record the assembly in ``asm_dir``, marked as inferred, when it holds contigs."""
    _r = _registry_module()
    stats = summarise_contigs(asm_dir / _CONTIGS_NAME)
    if stats["contigs"] > 0:
        _r.record_assembly(registry, acc, genome_id, asm_dir, stats, "", {})
        rb.mark_inferred(registry, acc, "extractions", genome_id=genome_id)


def _bootstrap_extractions(registry: "Registry", paths: "ProjectPaths", genome_ids: List[str]) -> None:
    for acc, per_genome_files in scan_extractions(paths.targeted, genome_ids).items():
        for genome_id, files in per_genome_files.items():
            _infer_extraction(registry, acc, genome_id, files)
    for acc, per_genome_asm in scan_assemblies(paths.targeted, genome_ids).items():
        for genome_id, asm_dir in per_genome_asm.items():
            _infer_assembly(registry, acc, genome_id, asm_dir)


def bootstrap_from_disk(
    paths: "ProjectPaths",
    accessions_file: Optional[Union[str, Path]] = None,
    parsed_containment: Optional[Union[str, Path]] = None,
    target_path: Optional[Path] = None,
) -> "Registry":
    """Rebuild a registry from what is on disk; every reconstructed block carries ``"inferred": true``.

    ``target_path`` is the file the caller intends to save this registry to (``status --init``
    passes ``registry_file``); it is bound to the registry before any writer below runs, so every
    path recorded during bootstrap is already relative to the right project root, even when
    ``--registry`` points somewhere other than the working directory. Left as ``None`` (a caller
    with no file in mind, e.g. the in-memory report built when no registry exists yet), the
    writers fall back to the working directory as the project root, same as an unbound
    ``Registry()``.
    """
    _r = _registry_module()
    registry = _r.Registry(path=target_path)
    if paths.matches.is_dir():
        _screening_from_matches(registry, paths.matches)
    _bootstrap_selection(registry, accessions_file, parsed_containment)
    _bootstrap_downloads_and_metadata(registry, paths)
    genome_ids = sorted(_genome_ids_on_disk(paths, registry))
    for genome_id in genome_ids:
        registry.genomes.setdefault(genome_id, {})
    _bootstrap_extractions(registry, paths, genome_ids)
    return registry
