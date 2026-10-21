# sra_sketch: sourmash sketches of downloaded datasets, with local containment

Date: 2026-10-21. Status: approved design. Implementation plan:
`docs/superpowers/plans/2026-10-21-local-sketch.md`.

## Purpose

PREPARE finds candidate metagenomes through Branchwater, whose index holds FracMinHash sketches
(k = 21, scaled = 1000) of about 1.16 million SRA runs. Two situations are not covered: a run that
Branchwater has not indexed, and a Branchwater hit that should be confirmed on the reads actually
downloaded. Both need a sketch of the local dataset built the same way, and a containment value
computed from it. The command `sra_sketch` builds one sourmash signature per downloaded
accession and, when a genome is given, writes its containment in each dataset to the matches
folder in the layout `parse_containment` and `select_datasets` already read.

Decisions taken with the user on 2026-10-21: sketches come from FASTQ files already in the
project or the store, not from a stream; containment is included; the earlier draft of a
streaming command (first N spots through `fastq-dump`) stays a later addition that will reuse
the signature writer from this design.

## Behaviour

`prepare sra_sketch [--fastq-folder fastq] [--accession ACC ...| --accessions-file F]
[--output-dir sketches] [--genome-fasta G | --signature S] [--genome-id ID]
[--matches-folder matches] [--max-workers N] [--force] [--dry-run] [--registry R] [--data-root D]`

For each accession the command finds the mate and orphan FASTQ files
(`prepare.sra.dataset_stats.accession_fastq_files`), reads every record of every file, and
builds three DNA FracMinHash sketches, k = 21, 31 and 51, scaled 1000, abundance tracked. The
three are written as one gzip-compressed JSON signature file `<output-dir>/<ACC>.sig.gz`, in the
format `sourmash sketch -o x.sig.gz` writes, named by the accession, k = 21 first so that
`prepare.data.branchwater_search.load_signature` finds the Branchwater-compatible sketch in it.

An accession is skipped when the registry's `analyses.sketch.summary.files` equals the current
sizes and modification times of its files and the signature file exists; `--force` sketches it
again. When the accession folder is a symlink into a shared store and the store holds a
signature for the same files, the file is copied from the store instead of sketched; after a
fresh sketch the signature is copied into the store beside the dataset and its description is
written to the sidecar, so another project linking the same dataset reuses it.

With `--genome-fasta` the genome is sketched at k = 21, scaled 1000 (`sketch_fasta`); with
`--signature` a plain `.sig` JSON file is read and its k = 21 sketch selected (`load_signature`).
The containment of the genome in each dataset is `genome.contained_by(dataset)` on flattened
sketches, which is the fraction of the genome's hashes present in the dataset. Rows
`acc,containment,cANI` (cANI = containment to the power 1/21, four decimals, as
`write_branchwater_csv`) are merged into `<matches-folder>/<genome_id>.csv`: existing rows for
other accessions and their metadata columns are kept, rows for the accessions of this run are
replaced, and the file is sorted by containment. The stem stays the genome id because
`parse_containment_data` keys genomes on it and keeps the maximum of duplicate rows.

`--dry-run` lists each accession with the action (sketch, reuse from store, skip) and writes
nothing. The command needs no external tool; sourmash is imported at the point of use and a
missing package is a `ConfigurationError` naming the `sourmash` extra (exit 3) before any work.

## Records

Registry, per accession sketched or reused: `record_analysis(acc, "sketch", output, summary)`
with `summary = {ksizes: [21, 31, 51], scaled: 1000, track_abundance: true, reads, bases,
hashes: {"21": n, "31": n, "51": n}, files: {name: [bytes, mtime]}, source: "computed" | "store"}`,
and a store usage row `(acc, "", "analysed", "sketch")`. Run log summary: accessions, sketched,
reused, skipped, failed, reads, bases, max_workers, genome_id, containment_rows; rows under
`analyses/sketch` with `reads, bases, hashes_k21, hashes_k31, hashes_k51, source` and, with a
genome, `containment, cani`. `status --export-tsv` lists `sketch` in the analyses column.

Store sidecar: an optional `sketch` block with the summary fields plus `file` and `computed`,
absent from the JSON when unset, so the sidecar schema is unchanged and older readers ignore it.

## Construction

- `prepare/data/sra/sampling.iter_sequence_blocks(path, stop)`: binary 1 MiB blocks through
  the existing gzip-aware opener, newline index with numpy, the sequence line of every complete
  record sliced and joined with `N`; the partial record carried to the next block. One
  `add_sequence(block, force=True)` call per sketch per block: a k-mer spanning a junction
  contains the `N` and is skipped, so the hashes equal those of per-read calls (verified with
  sourmash 4.9.4 at all three k). A non-FASTQ first byte or a truncated last record is a
  `ValueError`; a set stop token ends the read with `InterruptedError`.
- `prepare/sra/sketch.py`: `sketch_dataset`, `write_signature` (atomic, `save_signatures_to_json`
  with compression), `sketch_summary`, `sketch_is_current`, `load_dataset_minhash`.
- `prepare/sra/sketch_containment.py`: `genome_minhash`, `containment_rows`,
  `merge_containment_csv`.
- `prepare/store/sketches.py`: `cached_sketch`, `store_sketch`, mirroring `store/stats.py`
  (file signature for invalidation, the dataset lock with the stats wait, a held lock is a
  warning and never a failure). The signature is copied with the atomic helper, so a dataset
  folder never holds a partial `.sig.gz`.
- `prepare/cli/commands/sra_sketch.py`: planning, a thread pool bounded by `--max-workers` or
  `min(available_cpus(), max_workers_cap)` (sourmash releases the GIL in `add_sequence`),
  interrupt handling as in the download layer (stop, cancel pending, record what finished,
  exit 130), recording through a `RegistryBatch` after the work.

## Limits

Three abundance-tracking sketches of a deep metagenome hold millions of hashes and can take a
few hundred megabytes per worker; `--max-workers` bounds it. A project linked with
`--link-mode copy` gets no store cache, as for the statistics record. `--signature` must be a
plain JSON `.sig`; a `.sig.gz` from this command is a dataset signature, not a genome. The
first reads of a run are not sketched preferentially; the whole run is read, so a 10 GB file
takes minutes, not seconds. The store's byte totals include the signature files.
