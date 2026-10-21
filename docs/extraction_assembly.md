# Targeted read extraction and assembly

This guide covers `extract_target_reads`, which maps reads to a target genome and optionally assembles the mapped reads of each sample. Return to the [documentation index](README.md).

A targeted assembly uses only the reads relevant to one genome, which gives a small assembly instead of a whole-metagenome assembly. For every sample whose containment for the target genome meets the threshold, the command maps the reads with minimap2, keeps the mapped reads with samtools and writes them to `targeted/<ACC>/`. An accession excluded in the registry (`blacklist`) is always skipped. When the registry holds a selection (`select_datasets`), only the selected accessions that meet the threshold are extracted; `--ignore-selection` extracts every sample at the threshold instead. Without a selection, the table alone decides.

```bash
prepare extract_target_reads \
  --parsed-containment parsed_containment.txt \
  --genome-id GCF_000006945.2 \
  --genome-fasta genomes/GCF_000006945.2.fna \
  --fastq-folder fastq --output-folder targeted \
  --threshold 0.5 --preset sr
```

## Options

- `--dry-run` lists the qualifying samples without running any tool, and never checks that the tools are installed. A sample that was selected by containment but is not yet downloaded is common on a broad search, so the dry run prints one summary line for all such samples. A `fastq/<ACC>` symlink whose target is missing, typically an unmounted shared store, is logged as one warning that names every dangling link.
- `--preset` matches the read type: `sr` for Illumina, and `map-ont`, `map-pb` or `map-hifi` for long reads.
- `--threads` sets the threads for minimap2 and samtools (default 4).
- `--assemble` runs megahit on each sample's extracted reads.
- `--data-root` names a shared data store to resolve `--fastq-folder` against (see [shared_store.md](shared_store.md)).
- `--timeout SECONDS` stops minimap2, samtools or megahit after that time. The default, 0, means no limit.
- `--force` redoes extraction and assembly even when the registry says they are done.

The command needs `minimap2` and `samtools`, and `megahit` with `--assemble`, on the `PATH`.

## Mapping and coverage

Mapped reads always drop unmapped, secondary and supplementary alignments. `--min-mapq` additionally discards records below a mapping-quality threshold (default 0, which keeps every mapped record). A value of 20 is reasonable for a close relative of the target genome. A divergent strain can map with a low MAPQ, so a higher threshold can discard real matches.

For each sample with kept reads, `samtools coverage` runs on the kept alignments. It skips duplicate and QC-fail reads by default and writes `targeted/<ACC>/<genome>_coverage.tsv`. The registry records the breadth of the reference covered at 1x or more and the length-weighted mean depth. A failure of this step is logged as a warning and leaves the extracted reads in place. `--debug-keep-sam` keeps the intermediate SAM files.

## Skipping and redoing

A sample already extracted with the same genome FASTA, preset and `--min-mapq` is skipped on a rerun, including samples that mapped zero reads. The containment threshold only decides which samples are tried; it is recorded with the extraction but a different threshold does not redo a sample. `--force` redoes extraction and assembly. Re-extracting a sample drops its recorded assembly from the registry, since that assembly was built from the reads being replaced. The assembly folder stays on disk.

## Incomplete downloads and --allow-truncated

The command does not extract a sample whose download is known to be incomplete. This covers a sample whose registry download verdict is `truncated` and a sample whose store copy is recorded by its sidecar as `partial` or `failed`, whatever the registry verdict says. The skip line gives the reason, for example:

```text
skipped SRR1: store copy partial (5 of 20 spots); use --allow-truncated
```

`--allow-truncated` extracts such a sample anyway. A project download outside the store that is recorded `unverified` is still extracted. Only the sidecars of the selected samples are read for this check. Verdicts are described in [reads.md](reads.md).

## Assembly

With `--assemble`, each selected sample is assembled on its own with megahit. The assembly goes to `<genome_id>_assembly` in the sample's folder.

- `--assembly-preset` selects megahit's `--presets` value: `meta-sensitive` (the default, suited to small targeted read sets), `meta-large`, or `default` (no `--presets` flag).
- `--min-contig-len` sets the minimum contig length, and `--keep-intermediate` keeps megahit's `intermediate_contigs/` folder.
- Unless `--no-coverage` is given, the extracted reads are mapped back onto the assembled contigs to report a mapping rate and an estimated mean depth next to the other assembly statistics.
- A megahit failure is reported with the tool's own error message (the last lines of its stderr), not only the exit code.

### Threads and memory on macOS and elsewhere

On macOS the assembly defaults to one thread, because the parallel k-mer sorting step of megahit 1.2.9 is unstable on recent macOS releases. Mapping with minimap2 and samtools still uses `--threads`. `--assembly-threads` overrides the assembly thread count if a megahit build handles more.

`--assembly-memory` sets megahit's `--memory`. The value `auto` (the default; also `PREPARE_ASSEMBLY_MEMORY` or `assembly_memory` under `[runtime]`) gives 90% of the memory limit detected for the job, and leaves megahit's default when no limit is found, as on macOS. A size such as `32G` or `32000M` is passed in bytes. A fraction such as `0.5` is passed as it is, and megahit applies it to the whole node's memory, so give a size under a scheduler. Memory sizing on a cluster is described in [hpc.md](hpc.md).

### Scratch space

megahit needs FIFOs for its scratch files, which some filesystems do not provide (ExFAT, some network shares). `--temp-folder DIR` points the scratch space at a local POSIX filesystem. Without it, each run creates its own folder directly under `--output-folder`, named `.megahit-tmp-<random suffix>`. The folder is a sibling of the per-accession assembly directories and never inside one, because megahit refuses to run when its output directory exists. It is removed afterwards, also on failure. The random suffix lets two concurrent runs share one output folder. The `._*` AppleDouble files that macOS leaves on ExFAT and SMB volumes are ignored wherever PREPARE lists a folder, so they never look like FASTQ or genome files.

### Assembly markers and reuse

Each assembly folder carries a marker, `.prepare-assembly.json`, that records what the assembly was built from:

- the extracted read files (names and sizes) and their extraction date,
- the megahit preset, the minimum contig length and the k values,
- the megahit version, which is kept but not compared.

The registry's assembly record holds the same `inputs`. On a rerun, an assembly whose marker matches is reused and reported as reused. An assembly whose inputs differ is assembled again, and the reason is logged. Deleting the extracted reads therefore leads to a re-extraction, a new extraction date and a new assembly. A folder without `final.contigs.fa`, as an interrupted megahit leaves it, is assembled again with a warning.

An assembly folder without a marker is accepted, and given a marker whose megahit version is recorded as unknown, only when its registry record matches the preset and minimum contig length and is not older than the extraction. A value missing from the record matches any setting. A recorded preset, `default` included, matches only the same `--assembly-preset`.

### Replacing an assembly safely

megahit writes into a hidden staging folder beside `<genome_id>_assembly`, and the result replaces the previous folder by rename. A failed or interrupted assembly therefore leaves the earlier one in place. A staging folder left by a killed run is removed at the start of the next one. A run killed between moving the earlier folder aside and renaming the new one in is completed by the next run. That run puts the new folder in place when megahit had finished and its marker was written, and puts the earlier folder back otherwise. megahit exiting without writing `final.contigs.fa` is an error. A forced or stale redo that fails keeps the previous folder on disk, but without a registry record. A later run records it again from its marker if the marker still matches.

### Failures and exit status

A failure of megahit, or of the coverage mapping onto its contigs, is logged for that sample and the remaining samples continue. The run ends with one error line that lists the failed samples and exits with status 1. A sample that another process is extracting or assembling is skipped and reported, and does not fail the run. The closing line counts the samples assembled, reused, failed and busy elsewhere. Each sample's assembly step loads the project registry once, also for a sample whose assembly is reused, which adds up on a large registry when every assembly is reused.

Results appear in `status` and `results_table` (see [project_state.md](project_state.md)).
