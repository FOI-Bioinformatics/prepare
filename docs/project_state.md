# Project state, results and run history

This guide describes the project registry and the commands that read it: `status`, `results_table`, `blacklist`, `runs` and `project_report`. Return to the [documentation index](README.md).

## The registry

PREPARE keeps a journal of every dataset a project touches in `prepare_registry.json` in the project root. It records which genomes a dataset was screened against and with what containment, whether it was selected (by which threshold and filter), whether it is excluded and why, the download outcome with file sizes and dates, which analyses ran, and for each target genome the number of mapped reads and the assembly statistics. Commands update it as they finish. `status` reads it and always re-checks the disk, so a deleted folder shows up as missing rather than done. Commit `prepare_registry.json` with the project if the decisions should travel with the results. The design of the file is described in [ARCHITECTURE.md](ARCHITECTURE.md).

The registry keeps the latest outcome of each step. Two limits apply:

- Each `select_datasets` run replaces the previous selection. Only the accessions of the latest run count as selected.
- To keep the registry small on a broad search, `branchwater_search` and `parse_containment` record at most 5000 accessions per genome, the ones with the highest containment. `--registry-max-screened` of those two commands changes the limit. The match CSVs keep every hit.

## status

`status` prints an accession by stage matrix per genome from the registry and the disk. The stages are screened, selected, excluded, downloaded, analysed, extracted and assembled.

```bash
prepare status --init                      # create the registry from an existing project
prepare status                             # accession by stage matrix, per genome
prepare status --stage extracted --genome GCF_000008025.1
prepare status --next                      # which commands would advance the most datasets
prepare status --reconcile                 # record files removed by hand, register untracked work
prepare status --export-tsv registry       # registry_datasets.tsv and registry_extractions.tsv
prepare status --json                      # the report as one JSON document
```

### Creating and reconciling the registry

- `--init` creates the registry from what is on disk. It fills every metadata block it creates from an `<ACCESSION>_metadata.xml` that records a spot count (spot count, size, md5, assay type, organism, collection date, library layout and strategy, platform), in the same write that creates the registry.
- `--reconcile` compares the registry with the disk, records files that were removed by hand and registers untracked work. It re-checks `unverified` download verdicts and fills metadata in the same way as `--init` (see [reads.md](reads.md) for the spot-count sources). A store link costs nothing for this check, since its read count comes from the sidecar, but a plain project download is counted again from its mate-1 file and its file of unpaired reads. The first reconcile after spot counts become known therefore reads, once, the mate-1 file of every `unverified` plain download that has a spot count; later runs read nothing again. On a large project, run that first reconcile as a cluster job rather than on a login node.

The reconcile report, and the `drift` object of `status --json`, add four lists:

| List | Content |
| --- | --- |
| `store_unavailable` | Project links into a data store that is not mounted. They are left alone rather than marked missing and also appear under `dangling_links`. |
| `metadata_filled` | Accessions whose metadata block was filled from an XML file. |
| `verdicts_rechecked` | Accessions whose `unverified` verdict became `complete` or `truncated`. |
| `assemblies_dropped` | `ACC/GENOME` assembly records dated before their extraction record. Reconcile removes the record and leaves the folder on disk, and the next `extract_target_reads --assemble` builds the assembly again. |

The four keys are always present in `--json`, as empty lists when there is nothing to report. The text report prints a warning for the first and a line for each of the other three, only when the list is not empty.

### Suggested next commands

`status --next` lists a runnable `download_sra` command for accessions that are ready to download. For accessions that are still selected under a `select_datasets --no-skip-excluded` run, it prints a `select_datasets ... --skip-excluded` command instead. The command reproduces that run's genome, threshold, metadata and top-N and run filter criteria, so that rerunning it drops the excluded accessions from the selection file first.

### Wanted lists and folder locations

Given a wanted list, `status` also reports which of those accessions are present on disk and which are missing, so nothing is fetched twice. The list is an accession file or the accessions of a containment table:

```bash
prepare status --accessions-file accessions.txt --list-missing
prepare status --parsed-containment parsed_containment.txt --json
```

`--list-missing` also prints the missing accessions. Presence is checked against the standard layout: `fastq/<ACC>/`, `metadata/<ACC>_metadata.xml` and `genomes/<GCF>.fna`. `--fastq-folder`, `--metadata-folder` and `--genomes-folder` override these locations.

### Funnel and timing

The text report has a funnel line across the stages, for example:

```text
funnel: 120 screened, 40 selected, 38 downloaded (1.2 TB, 41 h), 30 extracted, 12 assembled
```

Each number is the count of datasets in that stage. The downloaded stage adds the size of the downloaded data and the recorded download time. The time is left out when no download was timed, for example for data linked or adopted from a store.

The `funnel` key of `status --json` gives each stage in full:

- `selected` adds the excluded count.
- `downloaded` adds the bytes and seconds of the downloaded datasets, the failed count and the seconds of the failed attempts (`failed_seconds`).
- `analysed` gives the analysed count.
- `extracted` and `assembled` count accession and genome pairs and give their seconds. Assemblies also give the total assembled bases.

The extraction and assembly figures cover only the pairs that reached the stage (mapped reads above zero, contigs above zero). They can therefore be lower than the totals of the timing line, which also count attempts that did not. A stage with no recorded time has `seconds` set to null.

`status --json` summarises the recorded times under `timing` (counts, totals and medians). The text report adds one timing line when anything was timed, and `status --export-tsv` adds the same columns.

### Exports

`status --export-tsv registry` writes `registry_datasets.tsv` and `registry_extractions.tsv`. The extraction table has the timing columns followed by `coverage_tsv`, `n90`, `largest` and `assembly_dir`.

`extract_target_reads` skips samples already extracted or assembled with the same genome, preset and threshold, and `--force` redoes them (see [extraction_assembly.md](extraction_assembly.md)).

## results_table

`results_table` writes one tab-separated row per screened accession and genome. It joins containment, selection, exclusion, download, run size, the dataset profile, mapped reads, reference coverage and assembly statistics.

```bash
prepare results_table --output results.tsv
```

The export is recorded in the registry when one exists, and no registry is created when there is none. `--no-record` leaves the registry and the run log unchanged. The columns are:

| Column | Content |
| --- | --- |
| `accession`, `genome_id` | The accession and the screened genome. |
| `containment` | Containment from the parsed table, unrounded. |
| `selected`, `excluded`, `exclusion_reason` | Selection and exclusion state. |
| `download_state`, `download_verdict` | Download state, and the recorded completeness verdict (`complete`, `truncated` or `unverified`; empty where none was recorded). |
| `run_total_spots`, `run_size` | Run metadata from NCBI. |
| `total_reads`, `gc_percent`, `quality_grade` | Dataset profile. GC is in percent. |
| `quality_source` | The analysis that supplied the profile columns: `profile`, `report`, `legacy` (an older analysis record), or empty when none was recorded. |
| `mapped_reads`, `mapping_rate_to_reference`, `breadth`, `mean_depth` | Read extraction and reference coverage. The mapping rate is mapped reads over the run's spot count. |
| `contigs`, `total_bp`, `n50` | Assembly statistics. |
| `assembly_largest`, `assembly_n90`, `assembly_gc_percent`, `assembly_contigs_ge_1kb` | Further contig statistics. `assembly_gc_percent` has two decimals. |
| `genome_fraction_estimate`, `assembly_mapping_rate` | Assembly figures as recorded by the assembly step. |
| `assembly_mean_depth_estimate` | Filled only when the coverage mapping onto the contigs ran. |
| `download_seconds`, `extraction_seconds`, `assembly_seconds` | Recorded run times, empty where a step was not timed or the dataset was linked from the store. |
| `assembly_dir`, `coverage_tsv` | The assembly folder and the coverage table, as paths relative to the project. |

Columns are empty where the registry has nothing to report. The profile columns come from the newer of the `profile` analysis of `sra_profile` and the `report` analysis of `sra_report`, and a value that the newer one lacks is filled from the other. A project that only ran `sra_report` therefore has them too. The extraction time of a sample is measured from the end of the previous sample, so for the first sample it includes reading the containment table and building the minimap2 index. The assembly time is the megahit run alone, without the contig summary and the coverage mapping.

## blacklist

`blacklist` records an excluded accession with a reason in the registry and keeps `blacklist.txt` in step. `download_sra` skips excluded accessions automatically, and `select_datasets` skips them by default.

```bash
prepare blacklist --add SRR2517418 --reason "16S amplicon mislabelled as WGS"
prepare blacklist --list
prepare blacklist --remove SRR2517418
```

## runs

A command that keeps a run log adds one record per run to `.prepare/runs/runs.jsonl` beside the registry. The record holds the command, its arguments, start and finish times in UTC, run time, exit code, PREPARE version, host, process ID and a short summary. A run with larger results also writes `.prepare/runs/<run_id>.json`. The last 10 of these files are kept per command, and an older run keeps its line without one. The NCBI API key is replaced by `***` in what is recorded. A line that cannot be read, because an append was interrupted or it holds a byte that is not UTF-8, is skipped with a warning.

Run log rules:

- Nothing is written in a folder without a registry. A run log that cannot be written is reported as a warning and does not change the exit code.
- `PREPARE_RUN_LOG=false`, or `run_log = false` under `[runtime]` in the config file, turns the run log off.
- The commands that keep a run log are `download_sra`, `sra_profile`, `sra_report`, `results_table`, `extract_target_reads`, `store_verify`, `project_report`, `select_datasets`, `blacklist`, and `status` with `--init` or `--reconcile`.
- A dry run, `blacklist --list`, a plain `status` (also with `--json`, `--next` or `--export-tsv`), `doctor` and `runs` record nothing. For `results_table`, `select_datasets` and `project_report`, `--no-record` also leaves the run out of the run log.
- The summary of a run holds counts and totals. Per-accession results go to the run's detail file, one row per accession (per accession and genome for `extract_target_reads`). Examples are failed and downloaded accessions with reason and attempts for `download_sra`, the figures of each accession for `sra_profile` and `sra_report`, mapped reads, breadth and mean depth for `extract_target_reads` (for a sample skipped as already extracted, the values recorded earlier with `skipped: true`), the verdict and fix for `store_verify`, and the rank of each accession for `select_datasets`.
- PREPARE does not add `.prepare/` to a project's `.gitignore`. Add it there if the run log should not be committed with the registry.

`runs` reads the run log and records nothing itself. It finds the project the same way as `status`: the nearest `prepare_registry.json`, or `--registry`.

```bash
prepare runs                               # the last 20 runs, newest first (--limit 0 for all)
prepare runs --command sra_profile         # only the runs of one command
prepare runs --show latest                 # one run: arguments, times, summary and detail
prepare runs --diff previous latest        # compare two runs
prepare runs --accession SRR11011981       # the values recorded for one accession, run by run
prepare runs --diff previous latest --json # one JSON document, for scripts
```

A run is selected by its ID (for example `20261001T120501Z-sra_profile-3fa2`), by a prefix that matches only one run, or as `latest` or `previous`. With `--command`, these refer to that command's runs only.

- **`--diff`** lists every summary value of the two runs with the difference for numbers. It then compares the details row by row (one row per accession, or per accession and genome) and lists the rows found in only one run and the fields that changed. An accession found in more than one section of a detail has each field named after its section, for example `analyses.profile.gc_percent`. A run has a detail only when its command noted one and it is among the last 10 per command. A run without one is reported as not kept and is compared on its summary alone.
- **`--accession`** lists the runs whose kept detail holds values for that accession. It adds runs without a kept detail whose command line names the accession, as an argument or as `--option=ACC`. An accessions file named on the command line is not read. The number of other runs that could not be searched is given at the end.
- **Exit status.** `runs` exits 1 when the project has no run log, and also when a selector matches no run or more than one.

## project_report

`project_report` summarises the whole project in one report, read from the registry and the run log. It finds the project the same way as `status` and reads no FASTQ file.

```bash
prepare project_report                         # project_report/project_report.md, .json and .html
prepare project_report --output-dir reports --max-rows 50
prepare project_report --html never --no-environment
```

The report has ten sections:

1. The project: registry path, last update, and the registry and PREPARE versions.
2. The funnel of `status`.
3. Per genome, the datasets extracted, assembled and with no mapped reads, with the median breadth and depth.
4. The extraction rows of `results_table`, sorted by genome and then by mapped reads.
5. The download states and completeness verdicts of the downloaded datasets, with the truncated and unverified ones among them.
6. The failed downloads with their reason, the number of attempts over all runs (`attempts_total`) and the date. The report CSV and `download_run.json` of `download_sra` count the attempts of one run.
7. Timing: the totals and medians of `status` plus the minimum, quartiles, 90th percentile and maximum of each step.
8. The environment checks of `doctor`, without network access. `--no-environment` leaves them out.
9. The recorded exports and analysis outputs.
10. The last 10 runs of the run log.

Tables and accession lists are cut to `--max-rows` rows (200 by default, 0 for all) and say how many rows there were. `results_table` writes every extraction row.

`project_report.md` and `project_report.json` are always written. `project_report.html` needs the `interactive` extra:

- `--html auto` (the default) reports a missing extra in one line and writes only the other two files.
- `--html always` exits with status 3 and writes nothing when the extra is missing.
- `--html never` skips the page.

The export is recorded in the registry as `project_report`, and the run in the run log, unless `--no-record` is given. Without a project registry the command exits with status 1 and writes nothing. Create a registry with `prepare status --init`. `--json` prints the written paths as one JSON document.
