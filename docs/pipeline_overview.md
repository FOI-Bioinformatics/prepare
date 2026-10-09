# Pipeline overview

PREPARE moves a population of public SRA datasets through six stages: screen, select, download,
analyse, extract, assemble. Each stage has one or two commands. Every command records what it decided in
the project registry (`prepare_registry.json` in the project root), so `status` can show where each
dataset stands at any time. This page lists the stages in order, with the commands, the files they produce
and what the registry keeps. Depth on a topic is in the page named in the relevant section; the flags of
every command are in `prepare <command> --help`.

## Stage summary

| Stage | Commands | Main output | Recorded in the registry |
|---|---|---|---|
| 1. Screen | `branchwater_search`, `use_branchwater`, `parse_containment` | `matches/*.csv`, `parsed_containment.txt` | containment and cANI per accession and genome |
| 2. Select | `select_datasets`, `blacklist` | `accessions.txt`, `blacklist.txt` | selection criteria and date; exclusions with reasons |
| 3. Download | `download_sra` | `fastq/<accession>/` | outcome per accession: downloaded, failed, skipped |
| 4. Analyse | `sra_validate`, `sra_profile`, `sra_report` | CSV, JSON and HTML reports | one analysis entry per accession |
| 5. Extract | `extract_target_reads` | `targeted/<accession>/<genome>_*.fastq.gz` | mapped read count, parameters, files |
| 6. Assemble | `extract_target_reads --assemble` | `targeted/<accession>/<genome>_assembly/final.contigs.fa` | contig count, total length, N50, megahit version |

Metadata commands run alongside stages 1 and 2 and are described after the stages. The Environment
commands `doctor`, `runs` and `project_report` are described at the end.

## 0. A target genome

Stage 1 needs a genome FASTA. `genome_prepare` searches GTDB by species or genus, downloads the matching
assemblies and writes a manifest CSV and a registry entry for each genome. `genome_download` fetches
assemblies by accession, species or genus, leaves the downloaded `genomes/ncbi_dataset.zip` where it
lands and records nothing in the registry. `genome_download` skips an accession whose FASTA already
exists; `--force` downloads it again and `--dry-run` reports which genomes are present and which are
missing without downloading. `genome_search` lists GTDB accessions without downloading (`--format list`,
the default, or `tsv`; see [analysis_commands.md](analysis_commands.md#genome_search)).

```bash
prepare genome_prepare --species "Lactobacillus crispatus" --output-dir genomes
prepare genome_download --accessions GCF_000006945.2 --output-dir genomes
prepare genome_download --accessions GCF_000006945.2 --dry-run
```

## 1. Screen

Find the metagenomes that contain the target genome.

```bash
# Query the Branchwater index with a genome FASTA (writes branchwater/<name>.csv)
prepare branchwater_search --genome-fasta genomes/GCF_000008025.1.fna

# Or place CSVs downloaded from the Branchwater web site in branchwater/ and convert them
prepare use_branchwater --branchwater-folder branchwater --matches-folder matches

# Build the sample by genome containment table and its summary
prepare parse_containment --matches-folder matches
```

`branchwater_search` needs the `sourmash` extra and network access. It sketches the FASTA with sourmash
(or takes a k=21, scaled=1000 signature with `--signature`), queries the public search API and writes
`branchwater/<genome>.csv`. The server is sent `--threshold` but is not always relied on to apply it, so
the client drops any returned match below the threshold as well. A search that returns zero matches is
reported as a warning, not an error. Check the index with a genome known to be abundant in metagenomes
before trusting an empty result.

A request that fails with a network error or a 429 or 5xx response retries up to 4 times with growing
backoff. Other 4xx client errors fail at once. A successful response is cached under
`.branchwater-cache/` next to the output CSV, keyed by the genome, containment threshold, sourmash
parameters and server URL, so a repeated search reads the cache. `--no-cache` neither reads nor writes
it, `--refresh` queries again and updates it, and `--max-cache-age-days N` treats an entry older than N
days as absent. Without these options a cache entry is reused until the cache file is deleted.

`use_branchwater` and `parse_containment` read CSVs from the Branchwater web site
(<https://branchwater.sourmash.bio/>) and from `branchwater_search`. A CSV needs only the columns `acc`
and `containment`. An unreadable file makes `use_branchwater`, `parse_containment` and
`extract_branchwater_metadata` exit with status 1 instead of skipping it.

`parse_containment` records every sample and does not filter. It writes:

- `parsed_containment.txt`: one row per SRA accession, one column per genome, plus `max_containment`;
- `top_containments.txt` by default (`--summary-containment-file`): how many samples reach each
  containment step, at `--step-size` intervals (default 0.1);
- `parsed_containment_details.tsv`: one row per accession and genome pair that appears in a match file,
  with the columns `parsed_containment.txt` does not carry.

| Column of the details table | Meaning |
|---|---|
| `accession` | SRA run accession |
| `genome_id` | target genome the row's containment refers to |
| `containment` | containment value for this accession and genome |
| `cANI` | Branchwater's estimated average nucleotide identity for the match |
| `biosample`, `bioproject` | NCBI sample and project identifiers, when Branchwater reports them |
| `assay_type`, `organism` | sequencing assay and the organism label from the metagenome record |
| `geo_loc_name`, `lat_lon` | geographic label and coordinates, when present |

Use the details table when a later step needs cANI or sample provenance. `parsed_containment.txt` stays
the samples-by-genomes matrix that the rest of the pipeline reads. `plot_containment` writes a plot image
(png by default). The registry keeps the screened accessions with their containment per genome, up to a
cap described in [project_state.md](project_state.md#the-registry).

## 2. Select

Decide which datasets are worth downloading, and which to leave out on purpose.

```bash
# Accessions whose containment for the genome is at or above the threshold
prepare select_datasets --genome-id GCF_000008025.1 --threshold 0.5 --output accessions.txt

# Only Illumina runs of at most 2 GB (needs the NCBI metadata table from parse_metadata)
prepare select_datasets --genome-id GCF_000008025.1 --threshold 0.5 --max-run-size 2G --platform ILLUMINA

# Exclude a dataset and say why
prepare blacklist --add SRR2517418 --reason "16S amplicon mislabelled as WGS"
```

The registry records the column, threshold, metadata filter, date and output file of the latest
selection, and a new run replaces it. Exclusions are kept with their reason in the registry and in
`blacklist.txt`, and `download_sra` honours them without further flags. `--no-record`, the filters (`--metadata-column`, `--metadata-value`,
`--max-run-size`, `--min-spots`, `--max-spots`, `--platform`, `--top-n`, `--skip-excluded`,
`--skip-downloaded`) and `sra_info`, which reports the sizes of an accession list before downloading, are
described in [reads.md](reads.md).

## 3. Download

```bash
prepare download_sra --accessions-file accessions.txt --fastq-folder fastq --max-workers 4
```

Reads land in `fastq/<accession>/` (as `<accession>_1.fastq.gz`, with `_2` for paired runs). `prefetch`
fetches the archive and `fasterq-dump` converts it; with `--no-prefetch`, `fasterq-dump` is called
directly. The FASTQ files are gzip-compressed unless `--no-compress` is given. As each download
finishes, the registry records it as downloaded, failed (with the attempt count and message) or skipped
(blacklisted or cut by `--max-downloads`). Failed accessions are also written to
`fastq/failed_accessions.txt`, and `--report-file` writes one row per accession. A rerun skips accessions
already on disk, `--force` downloads them again and `--dry-run` reports the plan without touching the
disk or the registry. Ctrl-C cancels downloads not yet started and stops the tools of those in progress.

Each download gets a completeness verdict from comparing its read count with NCBI's spot count:
`complete`, `truncated` or `unverified` when no spot count is known. The four places the spot count is
looked up, the verdict rules and `--redownload-truncated` are in [reads.md](reads.md).

With a shared data store configured (`--data-root`, `PREPARE_DATA` or the project's recorded store),
`download_sra` looks each accession up in the store first, links a usable copy into `fastq/<accession>`
and otherwise downloads once and adds the dataset to the store. See [shared_store.md](shared_store.md)
for the commands and [ARCHITECTURE.md](ARCHITECTURE.md) for the mechanics.

## 4. Analyse

```bash
prepare sra_validate --fastq-folder fastq
prepare sra_profile --fastq-folder fastq --accessions-file accessions.txt
prepare sra_report --fastq-folder fastq --groups-file groups.json --quality-profiles sra_quality_profiles
```

`sra_validate` checks the FASTQ files. With `--check-pairs` it compares mate read counts and with `--md5`
it compares each file's checksum with the recorded value. `sra_profile` writes one statistics table
(`sra_statistics.csv`: read and base totals, read lengths, GC content in percent, mean base quality,
grade) and one profile JSON per accession into `sra_quality_profiles/`. `sra_report` writes one HTML
report and, with `--groups-file`, a comparison of the groups with t-test or ANOVA results. Each command
records one analysis entry per accession in the registry (`validate`, `profile`, `report`). The sampling
rules, the statistics cache in the store sidecar and the report contents are in [reads.md](reads.md).

## 5. Extract

Map each selected sample against the target genome and keep only the reads that map.

```bash
prepare extract_target_reads --parsed-containment parsed_containment.txt \
    --genome-id GCF_000008025.1 --genome-fasta genomes/GCF_000008025.1.fna --threshold 0.1
```

Samples at or above the threshold in the containment table are mapped with minimap2 (preset `sr` by
default, `--threads` default 4) and the mapped reads are written with samtools to
`targeted/<accession>/` as `<genome>_*.fastq.gz`. The registry records the mapped read count, the genome
FASTA, preset, threshold and files, including for samples that mapped no reads. For samples with mapped
reads it also records the breadth of the reference covered at 1x or more and the length-weighted mean
depth, from `samtools coverage`, and writes per-sequence values to `<genome>_coverage.tsv`. A rerun with
the same genome, preset and threshold skips recorded samples, `--force` redoes them and `--dry-run` lists
the plan without running a tool. A sample whose download is `truncated`, or whose store copy is `partial` or
`failed`, is skipped unless `--allow-truncated` is given (see
[extraction_assembly.md](extraction_assembly.md)).

## 6. Assemble

```bash
prepare extract_target_reads --parsed-containment parsed_containment.txt \
    --genome-id GCF_000008025.1 --genome-fasta genomes/GCF_000008025.1.fna --assemble
```

With `--assemble` the extracted reads of each sample go through megahit into
`targeted/<accession>/<genome>_assembly/`. The registry records the contig count, total length, N50,
largest contig, the megahit version and the parameters. An assembly whose marker matches its inputs is
reused; one whose reads or settings changed is built again. The index, filtering options, assembly markers, scratch folders and the macOS thread
default are described in [extraction_assembly.md](extraction_assembly.md); memory sizing for cluster jobs
is in [hpc.md](hpc.md).

## Metadata

Two routes give per-sample metadata. They can be combined.

| Route | Command | Source | What it gives |
|---|---|---|---|
| Embedded fields | `extract_branchwater_metadata` | Columns of a Branchwater web-site CSV | Organism, biosample, bioproject, assay type, collection date, location and cANI, where the CSV has them. No NCBI request. |
| NCBI records | `download_metadata`, `parse_metadata` | NCBI, one XML per accession | Full sample and run records, including run size, checksum, library strategy and spot count. Needs an email address and is rate limited. |

A CSV written by `branchwater_search` carries only accession, containment and cANI, so
`extract_branchwater_metadata` has little to extract from it and `download_metadata` is the route that
works. Run `download_metadata` before `download_sra` when the completeness verdict or the size and
platform filters of `select_datasets` are needed.

`extract_branchwater_metadata` renames a few Branchwater columns to canonical names; in particular
`organism` becomes `Sample_Scientific_Name`. Use the output names (for example `--metadata-column
Sample_Scientific_Name`) in `count_metadata` and `single_sample`. A CSV from `branchwater_search` leaves the
`organism` column empty, so on that route `Sample_Scientific_Name` has values only in the NCBI table from
`parse_metadata`.

```bash
prepare extract_branchwater_metadata --branchwater-folder branchwater --metadata-folder metadata

prepare download_metadata --email you@example.org --matches-folder matches --metadata-folder metadata
prepare parse_metadata --metadata-folder metadata --metadata-table-file metadata_table.txt
```

`download_metadata --accessions-file accessions.txt` fetches metadata for the selected datasets only.
Accessions are requested from NCBI in batches (`--batch-size`, default 200), and `--api-key` or the
`NCBI_API_KEY` environment variable raises the rate limit. With `--data-root`, each XML is also copied
into `<data-root>/metadata/`, where `store_verify --spots` finds it later. The registry keeps the
metadata file, run size, run checksum, library strategy and organism per accession.

### Metadata table commands

`check_metadata_attributes`, `count_metadata` and `single_sample` read a metadata table. Without
`--file-path` (for `check_metadata_attributes`) or `--metadata-file` (for the other two) they use the first
of these that exists:

1. `metadata_table.txt`
2. `metadata/metadata_table.txt`
3. `metadata/branchwater_metadata.txt` (from `extract_branchwater_metadata`)

Either metadata route therefore works with the defaults. A table that `parse_metadata` writes under
another name has to be passed explicitly to each of these commands.

`check_metadata_attributes` counts the samples with a value for each attribute column of the table
(columns starting with `Run_`, `Project_`, `Sample_` or `Experiment_` are left out), which shows which
attributes are filled often enough to count. It writes one line per attribute, most frequent first
(default output `metadata_attribute_counts.txt`):

```bash
prepare check_metadata_attributes --file-path metadata_table.txt --output-file metadata_table_overview.txt
```

`count_metadata` counts the values of one metadata column among the samples at or above `--threshold`
(default 0.5) for each genome of the containment table (`--summary-file`, default
`parsed_containment.txt`). It writes `metadata_counts.txt` (one row per value, one column per genome) and,
next to it, the column totals in `metadata_counts_stats.txt` (`--stat-file` names it otherwise).
`plot_metadata_counts` plots either file and writes an image only when `--save-format` is given.

```bash
prepare count_metadata --metadata-column Sample_Scientific_Name --threshold 0.9
prepare count_metadata --summary-file parsed_containment.txt --metadata-file metadata_table.txt \
    --metadata-column Sample_Scientific_Name --threshold 0.95 --output-file genome_counts.txt
prepare plot_metadata_counts --file-path metadata_counts.txt --save-format png
```

`single_sample` counts the values of one metadata column among the samples at or above `--threshold`
(default 0.1) for one genome column of the containment table (`--summary-column`). It logs the `--top-n`
most frequent values (default 100) and writes no file:

```bash
prepare single_sample --summary-column GCF_000008025.1 --metadata-column Sample_Scientific_Name --threshold 0.95
```

## Project state

```bash
prepare status                              # inventory and the accession by stage matrix, per genome
prepare status --next                       # commands that would advance the most datasets
prepare results_table --output results.tsv  # one row per accession and genome
```

Presence is always re-checked on disk, so a folder deleted by hand shows up as missing. Decisions and
provenance come from the registry. The registry, the status options (`--stage`, `--genome`, `--init`,
`--reconcile`, `--export-tsv`, `--json`) and the `results_table` columns are described in
[project_state.md](project_state.md).

## Environment commands

These three commands read the environment or the project and are not stages.

`doctor` reports whether the environment can run PREPARE and changes nothing. Its flags, checks and
exit codes are in [configuration.md](configuration.md#checking-the-environment).

`runs` lists the commands recorded in the project's run log (`<project>/.prepare/runs/`), shows one run
(`--show`), compares two (`--diff`) or follows one accession through the runs (`--accession`).
`project_report` writes one report of the project from its registry as Markdown, JSON and, with the
`interactive` extra, HTML. Both are described in [project_state.md](project_state.md).

```bash
prepare runs --diff previous latest
prepare project_report --output-dir project_report
```

## The walkthrough

`make pipeline` runs `local_test.sh` on the bundled sample data. It runs `use_branchwater`,
`parse_containment`, `extract_branchwater_metadata`, `count_metadata`, `plot_metadata_counts`,
`plot_containment`, `select_datasets`, `status --init`, `blacklist`, `status --stage excluded`, `status`,
a dry-run `extract_target_reads` and a `status --json` check of the registry counts. It also steps
through the shared store: `store_init`, a `download_sra --dry-run` against the store, `store_status
--json` and `select_datasets --top-n`. It needs no network access and no external tools, and the store
steps redirect `HOME` and unset `PREPARE_DATA` so they never touch a real store or configuration file.

## External tools by stage

| Stage | Tool | Purpose |
|---|---|---|
| Genome | datasets (`ncbi-datasets-cli`) | `genome_download` and `genome_prepare` |
| Screen | sourmash (Python package) | sketching the genome for `branchwater_search` |
| Download | sra-tools (`fasterq-dump`; `prefetch` optional) | fetching reads |
| Download | pigz (optional) | parallel gzip compression; gzip is used without it |
| Extract | minimap2, samtools | mapping and filtering reads |
| Assemble | megahit (optional) | assembling the extracted reads; needed only for `--assemble` |
| Analyse | seqkit (optional) | faster read statistics for `sra_profile` and `sra_report` |

Python packages beyond the core install sit behind extras: `analysis` (scikit-learn, scipy) for
`diversity_analysis`, the PCA, t-SNE, diversity and clustered heatmap plots of `interactive_plot` and the
group statistics in `sra_report --groups-file`, `interactive` (plotly,
jinja2) for the HTML outputs, `maps` (cartopy) and `sourmash`. A command that needs a missing extra exits
with an error naming it. The download, parse, metadata, status and extraction stages need none of them.

`environment.yml` installs sra-tools, pigz, ncbi-datasets-cli, minimap2, samtools, megahit and sourmash,
leaves `seqkit` as a commented optional line, and installs PREPARE with every extra (`-e ".[all]"`).
`prepare doctor` reports which of these are present. The tool table and the extras are listed in
[README.md](../README.md).
