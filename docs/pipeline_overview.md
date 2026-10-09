# Pipeline Overview

PREPARE moves a population of public SRA datasets through six stages: screen, select, download,
analyse, extract, assemble. Each stage has one or two commands. Every command records what it decided
in the project registry (`prepare_registry.json` in the project root), so `status` can show where each
dataset stands at any time. This page lists the stages in order, the commands that belong to each, the
files they produce and what the registry keeps. The [README](../README.md) documents every flag; the
[Branchwater workflow](branchwater_workflow.md) covers the screening step in more depth.

## Stage summary

| Stage | Commands | Main output | Recorded in the registry |
|---|---|---|---|
| 1. Screen | `branchwater_search`, `use_branchwater`, `parse_containment` | `matches/*.csv`, `parsed_containment.txt` | containment and cANI per accession and genome |
| 2. Select | `select_datasets`, `blacklist` | `accessions.txt`, `blacklist.txt` | selection criteria and date; exclusions with reasons |
| 3. Download | `download_sra` | `fastq/<accession>/` | outcome per accession: downloaded, failed, skipped |
| 4. Analyse | `sra_validate`, `sra_profile`, `sra_report` | CSV, JSON and HTML reports | one analysis entry per accession |
| 5. Extract | `extract_target_reads` | `targeted/<accession>/<genome>_*.fastq.gz` | mapped read count, parameters, files |
| 6. Assemble | `extract_target_reads --assemble` | `targeted/<accession>/<genome>_assembly/final.contigs.fa` | contig count, total length, N50, megahit version |

Metadata commands run alongside stages 1 and 2 and are listed after the stages below. A genome FASTA is
needed before stage 1 can run at all: `genome_prepare` (search GTDB, download, write a manifest and a
registry entry) or `genome_download` (download by accession/species/genus, no registry entry) get one;
see the [README](../README.md#0-getting-a-target-genome-fasta).

## 1. Screen

Find the metagenomes that contain the target genome.

```bash
# Query the Branchwater index with a genome FASTA (writes branchwater/<name>.csv)
prepare branchwater_search --genome-fasta genomes/GCF_000008025.1.fna

# Or place CSVs downloaded from the Branchwater web site in matches/ and convert them
prepare use_branchwater --branchwater-folder branchwater --matches-folder matches

# Build the sample by genome containment table and its summary
prepare parse_containment --matches-folder matches
```

`parse_containment` writes `parsed_containment.txt` (one row per SRA accession, one column per genome,
plus `max_containment`) and, by default, `top_containments.txt` (how many samples pass each containment
step, at `--step-size` intervals, default 0.1; the name and step size shown in the README's example are
passed explicitly, not the defaults), as well as `parsed_containment_details.tsv` with cANI and the
sample metadata columns kept alongside the containment values. `plot_containment` writes a plot image
(png by default) rather than only displaying one. The registry keeps every screened accession with its
containment per genome, capped per genome by `--registry-max-screened` (default 5000). An unreadable
match file makes `use_branchwater`, `parse_containment`, and `extract_branchwater_metadata` exit with
status 1 rather than skipping it silently.

`branchwater_search` retries a failing request automatically and caches a response under
`.branchwater-cache/` next to its output, keyed by the genome, thresholds and server; `--no-cache`,
`--refresh` and `--max-cache-age-days` control the cache. The server is sent `--threshold` but is not
always relied on to honour it; the client re-applies the threshold to the response and drops anything
below it.

## 2. Select

Decide which datasets are worth downloading, and which ones to leave out on purpose.

```bash
# Accessions whose containment for the genome is at or above the threshold
prepare select_datasets --genome-id GCF_000008025.1 --threshold 0.5 --output accessions.txt

# The same, restricted to a metadata value
prepare select_datasets --genome-id GCF_000008025.1 --threshold 0.5 \
    --metadata-column Sample_Scientific_Name --metadata-value "Drosophila melanogaster"

# Only Illumina runs of at most 2 GB (needs the NCBI metadata table from parse_metadata)
prepare select_datasets --genome-id GCF_000008025.1 --threshold 0.5 --max-run-size 2G --platform ILLUMINA

# Exclude a dataset and say why
prepare blacklist --add SRR2517418 --reason "16S amplicon mislabelled as WGS"
```

The registry records the column, threshold, metadata filter, date and output file of the latest
selection; running `select_datasets` again replaces the previous selection. Exclusions are kept with
their reason in the registry and in `blacklist.txt`, and `download_sra` honours them without further
flags. `--top-n` keeps only the highest-containment accessions after every other filter; excluded
accessions are dropped by default (`--skip-excluded`) and `--skip-downloaded` also drops accessions the
registry already records as downloaded. `sra_info` (see [README](../README.md#downloading-reads))
reports sizes of any accession list before downloading. It filters per
experiment package, not per run: it lists every run of each experiment package that a requested run,
experiment, sample, study, BioProject or BioSample accession matches (including sibling lanes or
replicates of the same experiment as a requested run), and drops runs of packages that match none of the
requested accessions; a reply that would be filtered down to nothing is listed in full with a warning.
`--no-record` logs the counts and writes its list without touching the registry, for an exploratory
run that should not replace the recorded selection. It refuses to overwrite a file that a recorded
selection names (the default `accessions.txt` once a selection has been recorded), so `--output` should
name a scratch file.

`--max-run-size BYTES` (suffixes K, M, G, T as powers of 10), `--min-spots N`, `--max-spots N` and
`--platform NAME` filter on `Run_Size`, `Run_Total_Spots` and `Platform` in the NCBI metadata table,
after the metadata filter and before `--top-n`. A run absent from the table or without a value is dropped
and counted in the log. A table without the column a requested filter reads (every Branchwater-derived
table) is an error: the command exits with status 1, writes no output and records nothing. The registry
records these four values with the selection.

## 3. Download

```bash
prepare download_sra --accessions-file accessions.txt --fastq-folder fastq --max-workers 4
```

Reads land in `fastq/<accession>/` by default (as `<accession>_1.fastq.gz`, `_2` for paired runs):
`prefetch` fetches the archive, then `fasterq-dump` converts it (`--no-prefetch` calls fasterq-dump
directly), and the resulting FASTQ files are gzip-compressed unless `--no-compress` is given. As each
download finishes the registry records it as downloaded (with file sizes and a completeness verdict from
comparing the read count against NCBI's spot count -- this verdict is only available when
`download_metadata` (see "Metadata" below) has already fetched that accession's expected spot count, so
run `download_metadata` before `download_sra` if the verdict matters), failed (with the attempt count and
message) or skipped (blacklisted or cut by `--max-downloads`). Failed accessions are also written to
`fastq/failed_accessions.txt` for a retry, and `--report-file` writes one row per accession. A rerun
skips accessions already on disk; `--force` downloads them again. `--dry-run` reports the plan and
touches neither the disk nor the registry. `--temp-folder DIR` sets where `fasterq-dump` writes scratch
files; without it, a plain per-project download uses `fastq/.prepare-tmp/<accession>_fqtmp`, while a
store-backed one (below) defaults to `<data-root>/tmp/<accession>_fqtmp`; either is removed once the
accession's download ends. Ctrl-C cancels downloads not yet started and
stops the running `prefetch`/`fasterq-dump` processes for the ones in progress.

When a shared data store is configured (`--data-root`, `PREPARE_DATA`, or the project's recorded
store), `download_sra` looks up each accession in the store first. A complete or unverified dataset is
linked into `fastq/<accession>` without a new download; a partial dataset is re-downloaded unless
`--accept-partial`; an accession absent from the store is downloaded once, compressed, and added to the
store for every project to reuse. A per-accession lock (`--lock-wait`) lets two projects on the same
store cooperate on the same accession without downloading it twice. See the store commands
(`store_init`, `store_status`, `store_adopt`, `store_usage`, `store_gc`, and related commands) in the
[README](../README.md#shared-data-store) and the store package section of
[ARCHITECTURE.md](ARCHITECTURE.md) for the full design.

## 4. Analyse

```bash
prepare sra_validate --fastq-folder fastq
prepare sra_profile --fastq-folder fastq --accessions-file accessions.txt
prepare sra_report --fastq-folder fastq --groups-file groups.json --quality-profiles sra_quality_profiles
```

`sra_validate` checks the FASTQ files, and with `--check-pairs` compares mate read counts, and with
`--md5` compares each file's checksum against the value recorded for it. `sra_profile` writes one
statistics table (`sra_statistics.csv`: read and base totals, read lengths, GC content, mean base
quality, grade) and one profile JSON per accession into `sra_quality_profiles/`. `sra_report` writes
one HTML report on the quality of the datasets and, with `--groups-file`, a comparison of the groups
with t-test or ANOVA results; with `--quality-profiles` it reuses the profiles `sra_profile` saved
instead of profiling the reads again. Each of them records one analysis entry per accession in the
registry (`validate`, `profile`, `report`). All three take `--accessions-file`, and `sra_profile` and
`sra_validate` also take repeated `--accession` flags; without either, every accession folder is used.

Each dataset is profiled by one path. Read and base totals, mean read length and GC content come from
the dataset's statistics record: exact read counts, and GC from a uniform sample of the first mate
file's reads. Per-read quality, complexity and adapter figures come from a sample of `--sample-size`
reads per dataset (default 10000) drawn from every mate file (uniformly across each file by default,
or `--sampler head`). For a dataset held in the shared data store, the statistics record is cached in
the store sidecar (the `stats` block, invalidated when a FASTQ file's size or modification time
changes), so `sra_profile`, `sra_report` and `sra_validate --check-pairs` compute it once and reuse it
rather than re-reading the files each time. GC is reported in percent (0 to 100) under `gc_percent` in
the table, the profile JSON and the registry.

`sra_profile` and `sra_report` label every printed read total "(mates counted)", since a paired-end
run's two mate files are counted separately. Before 0.5.0 these steps were four commands,
`sra_stats`, `sra_profile_quality`, `sra_dashboard` and `sra_compare`; those names now print where
the command went and exit with status 2.

## 5. Extract

Map each selected sample against the target genome and keep only the reads that map.

```bash
prepare extract_target_reads --parsed-containment parsed_containment.txt \
    --genome-id GCF_000008025.1 --genome-fasta genomes/GCF_000008025.1.fna --threshold 0.1
```

Samples at or above the threshold in the containment table are mapped with minimap2 (preset `sr` by
default, `--threads` default 4) and the mapped reads are written with samtools to `targeted/<accession>/`
as `<genome>_*.fastq.gz`. The registry records the number of mapped reads, the genome FASTA, preset,
threshold and the files written, including samples that mapped no reads. For samples with mapped reads it
also records the breadth of the reference genome covered at 1x or more and the length-weighted mean depth,
computed with `samtools coverage` over the kept alignments (duplicate and QC-fail reads are skipped by its
defaults) and written per sequence to `targeted/<accession>/<genome>_coverage.tsv`.
A rerun with the same genome, preset and threshold skips samples already recorded; `--dry-run` lists what
would be extracted and what would be skipped without running any tool (and without requiring
minimap2/samtools/megahit to be on PATH), printing one summary line for samples not yet downloaded rather
than one warning per sample; `--force` redoes them.

The genome is indexed once per preset under `<output-folder>/.index/` and the index is reused across
samples rather than rebuilt each time. Unmapped, secondary and supplementary alignments are always
dropped; `--min-mapq` (default 0) additionally drops mapped records below a mapping-quality threshold.
`--fastq-folder` accepts a folder made of store symlinks the same way it accepts plain per-accession
folders, so extraction reads directly from the shared store when a project links into one; `--data-root`
resolves it against a shared store the same way the other store-aware commands do.

## 6. Assemble

```bash
prepare extract_target_reads --parsed-containment parsed_containment.txt \
    --genome-id GCF_000008025.1 --genome-fasta genomes/GCF_000008025.1.fna --assemble
```

With `--assemble` the extracted reads of each sample go through megahit into
`targeted/<accession>/<genome>_assembly/`. The registry records contig count, total length, N50, largest
contig, the megahit version and the parameters. An existing `final.contigs.fa` is kept unless `--force`;
an assembly folder without contigs is reported as interrupted, and `--force` removes it and runs again.
A megahit failure surfaces the tool's own stderr (last few lines), not just its exit code. On macOS
megahit runs single-threaded by default (`--assembly-threads` overrides).

megahit's scratch files need a filesystem with FIFOs, which ExFAT and some network shares do not
provide; `--temp-folder DIR` points that scratch at a local POSIX filesystem instead of the default.
Without `--temp-folder`, each run creates its own scratch folder directly under `--output-folder`,
named `.megahit-tmp-<random suffix>`, and removes it once the run finishes, even on failure; the random
suffix keeps two concurrent runs sharing one output folder from clashing. On such a volume, stray `._*`
AppleDouble sidecar files are ignored wherever PREPARE lists a folder's contents. A `fastq/<ACC>`
symlink whose target is missing is logged as one WARNING naming every such dangling link found.

`--assembly-preset` sets megahit's `--presets` value: `meta-sensitive` (the default), `meta-large`, or
`default` for no preset. `intermediate_contigs/` is removed after a successful assembly unless
`--keep-intermediate`. Unless `--no-coverage`, the extracted reads are mapped back onto the assembled
contigs to report `reads_mapped`, `mapping_rate` and an estimated `mean_depth` in the assembly stats.

## Metadata (alongside stages 1 and 2)

Run `download_metadata` before `download_sra` if the download's spot-count verdict matters (see stage 3
above); nothing in the command ordering below enforces this, it is on the user to sequence it that way.

```bash
# extract_branchwater_metadata only has real metadata to extract when the CSVs came from the
# Branchwater web site; a branchwater_search (API) CSV carries only accession, containment and cANI
prepare extract_branchwater_metadata --branchwater-folder branchwater --metadata-folder metadata

# Richer, and the only route that works when branchwater_search was used: NCBI metadata per
# accession, then a consolidated table
prepare download_metadata --email you@example.org --matches-folder matches --metadata-folder metadata
prepare parse_metadata --metadata-folder metadata --metadata-table-file metadata_table.txt

# Summaries
prepare count_metadata --metadata-column Sample_Scientific_Name --threshold 0.5
prepare plot_metadata_counts --file-path metadata_counts.txt --save-format png
```

`count_metadata` writes `metadata_counts.txt` by default and defaults `--threshold` to 0.5.
`plot_metadata_counts`, unlike `plot_containment`, only writes a file when `--save-format` is given.

`download_metadata --accessions-file accessions.txt` fetches metadata for the selected datasets only.
The registry keeps the metadata file, run size, run checksum, library strategy and organism per
accession. Accessions are fetched from NCBI in batches (`--batch-size`, default 200, one XML per
accession written from each batch response) rather than one request per accession; `--api-key` (or the
`NCBI_API_KEY` environment variable) raises the request rate limit. With `--data-root`, each fetched
XML is also copied into `<data-root>/metadata/` so every project sharing the store has it, and
`store_verify --spots` can find it there later even for a dataset downloaded before the metadata was
fetched.

## Project state at any point

```bash
prepare status                              # inventory and the accession by stage matrix, per genome
prepare status --stage extracted --genome GCF_000008025.1
prepare status --next                       # commands that would advance the most datasets
prepare status --init                       # build the registry for a project created before it existed
prepare status --reconcile                  # record files removed by hand, register untracked work
prepare status --export-tsv registry        # registry_datasets.tsv and registry_extractions.tsv
prepare status --json                       # the same report for scripts
prepare results_table --output results.tsv  # one row per accession and genome, the study's results
```

`results_table` joins the registry with the parsed containment table into one row per screened
accession and genome, from containment through read mapping, reference coverage and assembly.

Presence is always re-checked on disk, so a folder deleted by hand shows up as missing rather than done.
Decisions and provenance (why a dataset was selected or excluded, when it was downloaded, how many
reads mapped) come from the registry. Commit `prepare_registry.json` with the project if the decisions
should travel with the results.

`status --data-root` reports whether the project uses a shared data store, and if so the store's root,
whether each dataset is store-backed or local, and any dataset that is partial or whose store symlink is
dangling. A wanted accession whose `fastq/<ACC>` link points at a store dataset that is not yet
`complete` or `unverified` (still `downloading`, `failed`, or `partial`) is reported as linked to a store
dataset that is not complete, rather than simply listed as missing. `status --reconcile` records a
dangling link the same way it records other files removed by hand.

`status --next` also turns a `select_datasets --no-skip-excluded` selection's still-excluded accessions
into a runnable `select_datasets ... --skip-excluded` command, reproducing that run's genome, threshold,
metadata, top-N and run filter criteria, instead of suggesting `download_sra` directly on them.

## The walkthrough

`make pipeline` runs `local_test.sh` on the bundled sample data: `use_branchwater`, `parse_containment`,
`extract_branchwater_metadata`, `count_metadata`, `plot_metadata_counts`, `plot_containment`,
`select_datasets`, `status --init`, `blacklist`, `status --stage excluded`, `status`, a dry-run
`extract_target_reads`, and a `status --json` check of the registry counts. It also walks through the
shared data store on the same sample data: `store_init`, a `download_sra --dry-run` against the store
(so no network access is needed), `store_status --json`, and `select_datasets --top-n`. It needs no
network access and no external tools; the store steps redirect `HOME` and unset `PREPARE_DATA` so
they never touch a real store or the user's configuration file.

## External tools by stage

| Stage | Tool | Purpose |
|---|---|---|
| Screen | sourmash (Python package) | sketching the genome for `branchwater_search` |
| Download | sra-tools (`prefetch`, `fasterq-dump`) | fetching reads |
| Download | pigz (optional) | parallel gzip compression of downloaded FASTQ files |
| Extract | minimap2, samtools | mapping and filtering reads |
| Assemble | megahit | assembling the extracted reads |
| Analyse | seqkit (optional) | faster read statistics for `sra_profile` and `sra_report` |

Python packages beyond the core install sit behind extras: `analysis` (scikit-learn, scipy) for
`diversity_analysis` and the group statistics in `sra_report --groups-file`, `interactive` (plotly,
jinja2) for the HTML outputs (`explore_containment`, `sra_report`, `interactive_plot`), `maps`
(cartopy) and `sourmash`. A command that needs a missing extra exits with an error naming it; the
download, parse, metadata, status and extraction stages above need none of them.

The README's installation section lists them; `environment.yml` installs sra-tools, pigz,
ncbi-datasets-cli, minimap2, samtools, megahit and sourmash by default (pigz is a real dependency, not
optional there), comments out `seqkit` as one additional, optional line, and installs PREPARE
with every extra (`-e ".[all]"`).
