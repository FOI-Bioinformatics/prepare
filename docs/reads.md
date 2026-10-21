# Choosing, downloading and profiling reads

This guide covers the commands that pick datasets from the screening results, fetch their reads and check their quality. Return to the [documentation index](README.md).

The commands run in this order: `select_datasets` writes an accession list, `sra_info` shows sizes before a download, `download_sra` fetches the reads, `sra_validate` checks the files, and `sra_profile` and `sra_report` describe the data. Downloaded reads live in `fastq/<accession>/<accession>_1.fastq.gz`, with `_2` for paired runs. This is the layout that `status`, `sra_profile`, `sra_report`, `sra_validate` and `extract_target_reads` read.

## select_datasets

`select_datasets` ranks the screened accessions by containment and writes the ones that pass every filter to a text file with one accession per line.

```bash
prepare select_datasets --threshold 0.9 --output accessions.txt
prepare select_datasets --genome-id GCF_000008025.1 --threshold 0.5 \
    --metadata-column geo_loc_name_country_calc --metadata-value France --output accessions.txt
prepare select_datasets --threshold 0.9 --top-n 20 --output accessions.txt
prepare select_datasets --genome-ids GCF_000008025.1 GCF_000006945.2 --require all --threshold 0.5 \
    --output accessions.txt
prepare select_datasets --threshold 0.5 --max-run-size 2G --min-spots 1000000 --platform ILLUMINA \
    --output accessions.txt
```

Filter rules:

- **Input.** `--parsed-containment` names the containment table (default `parsed_containment.txt`).
- **Genome columns.** `--genome-id` ranks on one genome column, and `max_containment` is used when the flag is omitted. `--genome-ids` ranks on several columns together. With `--require any` or `--require all`, one or every listed column must meet the threshold.
- **Top N.** `--top-n N` keeps the N accessions with the highest containment after every other filter.
- **Exclusions.** Excluded accessions are skipped by default (`--skip-excluded`; `--no-skip-excluded` keeps them). `--skip-downloaded` also drops accessions that the registry records as downloaded, which is useful when selection is repeated on an expanded search.
- **Run filters.** `--max-run-size BYTES` (suffixes K, M, G and T are powers of 10, so `500M` is 500000000 bytes), `--min-spots`, `--max-spots` and `--platform` (case-insensitive, for example `ILLUMINA`) filter on the `Run_Size`, `Run_Total_Spots` and `Platform` columns of the NCBI metadata table. The table is `metadata_table.txt` from `parse_metadata`, found automatically, or `--metadata-file`. These filters apply after the threshold, the exclusions and the metadata filter, and before `--top-n`.
- **Missing values.** A run absent from the table, or without a value in the filtered column, is dropped because the bound cannot be checked, and the log counts such runs. A requested filter whose column is missing from the table is an error: the command exits with status 1 and writes no output and records no selection. Branchwater-derived tables carry none of these columns, so run `download_metadata` and `parse_metadata` for the candidate list first. When the column exists but every remaining candidate lacks a value, a warning gives the same hint. The log also reports the summed size of the selected runs.
- **Recording.** Each recorded run replaces the previous selection. Only the accessions of the latest run count as selected, and `status --stage selected` lists them.
- **Exploratory runs.** `--no-record` writes the output file and logs the counts but leaves the registry and the run log unchanged. Because `status --next` points `download_sra` at the file of the recorded selection, `--no-record` refuses to overwrite a file that a recorded selection names, including the default `accessions.txt`. Give it a scratch `--output` instead.

## sra_info

`sra_info` fetches NCBI metadata for the accessions and prints sizes, layouts, technologies and an estimated download time. It downloads no reads.

```bash
prepare sra_info --accessions-file accessions.txt --email you@example.org --output-report sra_info_report.csv
```

NCBI asks for an email address. It comes from `--email`, from the `PREPARE_NCBI_EMAIL` variable or from `ncbi_email` under `[runtime]` in the config file. `--api-key` raises the rate limit and `--bandwidth-mbps` sets the bandwidth used for the download estimate (default 100).

`sra_info` filters per experiment package, not per run. It lists every run of each experiment package that a requested run, experiment, sample, study, BioProject or BioSample accession matches. This includes sibling lanes or replicates of the same experiment as a requested run. Runs of packages that match none of the requested accessions are dropped. If a reply would be filtered down to nothing, the command lists every run returned and logs a warning instead.

## download_sra

`download_sra` runs `fasterq-dump` for each accession in parallel, skips accessions whose FASTQ files already exist (unless `--force`), retries failures and writes `fastq/failed_accessions.txt` for reruns. Ctrl-C cancels downloads that have not started, stops the running `prefetch` and `fasterq-dump` processes and exits without leaving orphaned tool processes.

```bash
prepare download_sra --accessions-file accessions.txt --fastq-folder fastq --max-workers 4 --num-threads 4
prepare download_sra --accessions-file accessions.txt --dry-run
prepare download_sra --accessions-file accessions.txt --report-file download_report.csv
```

Accessions excluded with `blacklist` are skipped on a real run, without `--blacklist blacklist.txt` on every call (that flag still works). `--dry-run` neither reads nor writes the registry, so it applies no blacklist exclusions and records nothing. See [project_state.md](project_state.md) for `blacklist`.

### Report files

`--report-file` writes one row per accession. The status is `downloaded`, `failed`, `already_present`, `blacklisted`, or `skipped` (skipped by `--max-downloads`). The other columns are:

| Column | Meaning |
| --- | --- |
| `seconds` | Time this accession took in this run. Empty for a dataset linked from the store and for accessions that were not downloaded. |
| `reason` | For a failed accession only: `network`, `not-found`, `disk-full`, `insufficient-space`, `locked`, `interrupted` or `unknown`. |
| `attempts` | Attempts in this run that started a download. A retry counts again. An accession that was refused for lack of space, cancelled by a disk-full abort, stopped by an interrupt before it started, given up after waiting for another run's lock, refused because the store holds a partial copy under `--no-resume-partial` or an incomplete copy whose refetches are used up, already present, or linked from the store counts 0. |

A run that reaches the download stage, or that stops in its preflight checks when the FASTQ folder already exists, also writes `fastq/download_run.json`, with or without `--report-file`. It holds the start and finish times, the totals (accessions, downloaded, failed, already present, blacklisted, skipped, attempts), the failures by reason, each failed accession with its reason, attempts and message, whether the run was aborted (`disk-full` or `interrupted`), the exit code, the host, the settings and the paths used. Both files are written after an interrupt (exit 130) and after a run whose final registry write failed. After an interrupt they list only the accessions that reported a result. A dry run writes neither file, and a failure to write them is logged as a warning without changing the exit code.

### Prefetch, compression and locking

- `prefetch` runs before `fasterq-dump`. `--no-prefetch` calls `fasterq-dump` directly. The `--max-size` given to `prefetch` is the `prefetch_max_size` setting (default `100G`; see [configuration.md](configuration.md)), and a run whose `.sra` archive is larger is not fetched.
- FASTQ files are gzip-compressed. `--no-compress` leaves them plain. `pigz` is used when installed, otherwise Python's gzip module.
- `--sra-cache DIR` sets where `prefetch` keeps its `.sra` archives (default `<fastq-folder>/.sra-cache`). `--keep-sra` keeps the archive after conversion. Only the archive of a `complete` or `unverified` download is kept, a `truncated` one is removed, and a redownload (`--force` or `--redownload-truncated`) always fetches the archive again.
- A project without a shared store takes a per-accession lock (`<fastq-folder>/.locks/<ACCESSION>.lock`) with the same heartbeat and 10-minute takeover as the store, bounded by `--lock-wait`. Two runs on one project therefore download each accession once. A `--sra-cache` folder shared by two projects without a store is not locked.
- Each download is built, verified and compressed under `<fastq-folder>/.prepare-tmp/` and then moved to `<fastq-folder>/<ACCESSION>` with one rename, so another run never sees a partly written folder. `status` and the download inventory ignore both hidden folders.
- A failed attempt removes its `<fastq-folder>/.prepare-tmp/<ACCESSION>_temp` build folder. `--keep-temp` keeps it for inspection and logs its path.
- `--temp-folder DIR` sets where `fasterq-dump` writes its scratch files. Without it, a per-project download uses `<fastq-folder>/.prepare-tmp/<ACCESSION>_fqtmp`, which is removed when the accession finishes, and a download into a shared store uses `<store>/tmp/<ACCESSION>_fqtmp`. The free-space check and `doctor` measure the same location. When the project folder is on a network filesystem, set `--temp-folder` or `PREPARE_TEMP_FOLDER` to local scratch space on the node. `extract_target_reads` accepts the same flag for megahit.
- `--timeout SECONDS` stops `prefetch` or `fasterq-dump` after that time. The default, 0, means no limit (also `PREPARE_TIMEOUT`, or `timeout` under `[runtime]`). A stopped tool is reported as a network failure, so `download_sra` retries it, and when every failed accession failed for that reason the exit status is 4. A `--version` probe used to record a tool version always has a fixed 30-second limit.
- `--max-workers N` sets how many accessions download at once, at most 10. Without it, the number is the CPUs available to the job divided by `--num-threads`, at least 1 and at most 4. The CPU count is the affinity mask, which reflects a SLURM allocation, then `SLURM_CPUS_PER_TASK`, then the machine's CPU count. `PREPARE_MAX_WORKERS_CAP` (or `max_workers_cap` under `[runtime]`) changes the cap of 4.

### Free space

Before each accession starts, `download_sra` checks the filesystems it writes to: the FASTQ folder (or the store's `tmp` folder), the `fasterq-dump` temporary folder and, with `prefetch`, the `.sra` cache. An accession whose run size is in the registry (from `download_metadata`) needs about 8 times that size for temporary files and 10 times for the FASTQ folder, plus the size itself for the cache. Locations on one filesystem add up, and downloads already running are counted. An accession of unknown size needs `--min-free-gb` (default 10; also `PREPARE_MIN_FREE_GB` or `min_free_gb` under `[runtime]`) on each filesystem.

A download that does not fit while others are running waits until they finish. One that would not fit even with nothing else running is not started. It fails alone with an `insufficient-space: not enough free space on ...` message, and the other accessions continue. `--min-free-gb 0` turns the check off, and a filesystem whose free space cannot be read is assumed to have room.

A disk that fills during a download stops the pass. The downloads not yet started are recorded as `disk-full: not attempted`, the running ones finish, no retry pass runs and the command exits with status 1.

### Completeness verdicts

Each download's read count is compared with NCBI's recorded spot count and a verdict is recorded in the registry; `--no-verify-downloads` turns the check off. The check also runs in a project with no store. Compressed files are read back after compression, and a file whose record count differs from the plain file's fails the download.

| Verdict | Meaning |
| --- | --- |
| `complete` | The read count per mate is at least 0.99 of the spot count. |
| `truncated` | Fewer reads than expected. The download fails: the staged files and the prefetch archive are discarded, the accession is retried and listed in `failed_accessions.txt`, and a run that NCBI serves short is therefore fetched again on every attempt. `--accept-partial` keeps it as a success instead. |
| `unverified` | No expected spot count is known, for example because `download_metadata` never ran for the accession. |

A dataset held in the shared store carries the store's own verdict instead: `complete`, `partial` or `unverified` (see [shared_store.md](shared_store.md)).

`--redownload-truncated` fetches an accession whose registry verdict is `truncated` instead of skipping it on a rerun. A recorded verdict is never downgraded to `unverified`. If a redownload keeps coming back `unverified`, run `download_metadata` for the accession so that an expected spot count exists. Otherwise a `truncated` record without a spot count is fetched again on every run.

An accession found on disk without a download record, as a killed run leaves it, is verified against its spot count before it is recorded. It is recorded `unverified` when no count is known, and a truncated one is reported with a pointer to `--redownload-truncated`. `--no-verify-downloads` skips the check of the project's own downloads and of a download into the shared store, whose copy is then recorded as unverified. A partial store copy is still refused without `--accept-partial`.

`extract_target_reads` does not extract a sample whose download is known to be incomplete unless `--allow-truncated` is given (see [extraction_assembly.md](extraction_assembly.md)).

### Where the spot count comes from

`download_sra`, the store hand-off, `store_link`, `store_adopt` and `status --reconcile` take the expected spot count from one lookup. The sources are tried in this order:

1. The registry's metadata block (`Run_Total_Spots`, from `download_metadata`).
2. The store sidecar's NCBI count.
3. The `<ACCESSION>_metadata.xml` file of the project's `metadata/` folder, then of the store's `metadata/` folder.
4. The count recorded with the previous download verdict.

Verdict rules that follow from this lookup:

- Linking a dataset from the store (by `download_sra`, `store_link` or `store_adopt`), or reconciling, never replaces a recorded `complete` or `truncated` verdict with `unverified` when no count is known.
- When the sidecar's read count and a spot count are both known, the verdict is computed again. A relink to a short copy therefore records `truncated`.
- A plain redownload that comes back `unverified` keeps a recorded `truncated` verdict. One that comes back `complete` replaces it.
- In a project linked with `--link-mode copy`, a store copy that a later `download_sra` run finds short is recorded with the store's verdict, so `extract_target_reads` skips the project's copy even though that copy's own sidecar still reads complete.
- `status --reconcile` re-checks every `unverified` verdict once a spot count is known. A plain project download is counted again from its mate-1 file and its file of unpaired reads. A link into the store is compared with the read count its sidecar records, without reading any FASTQ file. Before that, reconcile fills a metadata block that lacks a spot count from `metadata/<ACCESSION>_metadata.xml`, so one run fills both. An XML file without a spot count leaves the block unchanged.

### Resume and retries

A rerun skips accessions whose FASTQ files exist, so an interrupted run resumes by running the same command again. An interrupt during compression publishes nothing, so that accession is dumped again from the kept prefetch archive. `--max-retries` sets the number of retry attempts for failed downloads (default 1), and `fastq/failed_accessions.txt` lists the accessions that still failed (a run with no failures removes an earlier file). A project download that was accepted as `truncated` is fetched again only with `--redownload-truncated` or `--force`. For a store copy that is incomplete, `--no-resume-partial` stops the automatic refetch and `--accept-partial` uses the copy as it is (see [shared_store.md](shared_store.md)).

## sra_validate

`sra_validate` checks that each accession folder under `--fastq-folder` holds readable FASTQ files. `--check-pairs` also verifies that paired files have matching read counts, and `--md5` re-hashes each FASTQ file and compares it with the md5 the store recorded when the file was stored. It detects a file changed since then and has no effect outside a store. `--accessions-file` and repeated `--accession` flags restrict the check to the accessions named.

```bash
prepare sra_validate --fastq-folder fastq --check-pairs --accession SRR123456
```

## sra_profile

`sra_profile` computes statistics and a quality profile for each downloaded dataset. By default it profiles every accession folder under `--fastq-folder`. `--accessions-file` and repeated `--accession` flags restrict the run to the accessions named.

```bash
# Every accession folder under fastq/
prepare sra_profile --fastq-folder fastq

# Selected accessions, printing the path of each profile JSON as it is written
prepare sra_profile \
    --accessions-file accessions.txt \
    --accession SRR123456 \
    --output-dir quality_profiles \
    --detailed-reports
```

Each dataset is profiled once. Read and base totals, mean read length, N50 and GC content come from the dataset's statistics record. Read counts are exact, from `seqkit stats` when it is installed and otherwise from a streaming count. GC content comes from a uniform sample of the first mate file's reads. The record is cached in the store sidecar and reused by `sra_profile`, `sra_report` and `sra_validate --check-pairs` until a file's size or modification time changes. Per-read quality, complexity, duplication and adapter figures come from a sample of `--sample-size` reads per dataset (default 10000) drawn from every mate file. The sample is uniform across each file by default, or taken from the start of each file with `--sampler head`.

The command writes three things:

- The table `--output-report` (default `sra_statistics.csv`), with one row per accession.
- One `<accession>_quality_profile.json` per accession in `--output-dir` (default `sra_quality_profiles`), together with `quality_summary.json`.
- A `profile` analysis per accession in the project registry.

GC content is given in percent (0 to 100) under `gc_percent` in all three, and the profile JSON carries the complexity score under `complexity_score`. Printed read totals are labelled "(mates counted)". A paired-end run's two mate files are counted separately, so the figure is twice the spot count NCBI reports for the run. The command exits with status 1 when any accession has no readable FASTQ files. The other accessions are still profiled and written.

## sra_sketch

`sra_sketch` builds one sourmash signature per downloaded dataset and, when a genome is given, writes the containment of the genome in each dataset to the matches folder. It needs the `sourmash` extra and no external tool.

```bash
# Every accession folder under fastq/
prepare sra_sketch --fastq-folder fastq

# Selected accessions, compared with a genome; the result goes to matches/wMel.csv
prepare sra_sketch --accession SRR11011981 --genome-fasta wMel.fasta --genome-id wMel
```

Every record of every mate and orphan file of an accession is read once and hashed into three FracMinHash sketches, at k = 21, 31 and 51, scaled 1000, with abundance tracked. The three are written as `<output-dir>/<ACCESSION>.sig.gz` (default `sketches/`), the gzip JSON format `sourmash sketch -o x.sig.gz` writes, named by the accession and with the k = 21 sketch first. The k = 21 sketch is built exactly as `branchwater_search` builds its query, so `sourmash search`, `gather` and `compare` work on the files directly. An accession whose files have not changed in size or modification time since the recorded sketch is skipped; `--force` sketches it again. Without `--accession` or `--accessions-file`, every accession folder under the FASTQ folder is sketched and folders whose names are not accessions are listed once and skipped. When the accession folder is a link into a shared store, a signature cached beside the store copy (`<store>/sra/<ACCESSION>/<ACCESSION>.sig.gz`, described in the sidecar) is copied into the project instead of recomputed, and a fresh sketch is cached there for the next project. A project linked with `--link-mode copy` gets no cache.

With `--genome-fasta` the genome is sketched at k = 21, scaled 1000; `--signature` takes a `.sig` or `.sig.gz` of the genome instead (a `.sig.gz` written by this command describes a dataset, not a genome, and is not a sensible input). The containment of the genome in each dataset is the fraction of the genome's hashes present in the dataset sketch. The rows `acc`, `containment` and `cANI` are merged into `<matches-folder>/<genome_id>.csv` (default `matches/`, the genome id defaulting to the input file's stem): rows of other accessions and their metadata columns are kept, rows of the accessions sketched in this run are replaced, the row of an accession that failed while sketching is removed, an accession with no FASTQ files keeps its row, and the file is sorted by containment. A later `use_branchwater` replaces the file with the Branchwater CSV again when the two differ. `parse_containment` and `select_datasets` read it unchanged.

Accessions are sketched in parallel on `--max-workers` threads, by default the CPUs available to the job bounded by the `max_workers_cap` setting. Each worker holds three abundance sketches of its dataset in memory, a few hundred megabytes for a deep metagenome. The whole run is read, and the time grows faster than the data: sourmash keeps a sketch's hashes in a sorted array, so adding hashes to a sketch that already holds millions of them (a deep or error-rich metagenome at scaled 1000) is slow, and one such dataset can take an hour or more for the three k-mer sizes. Progress within a dataset is logged at the per-item level (`--progress-every 0` shows it). `--dry-run` lists what would be sketched, reused or skipped and writes nothing.

The command records a `sketch` analysis per accession in the project registry (k-mer sizes, scaled, read and base totals, hashes per k, the file snapshot, and whether the sketch was computed or reused from the store) and one run-log row per accession, with the containment when a genome was given. It exits with status 1 when any accession has no readable FASTQ files or cannot be read; the other accessions are still written. An interrupt records the accessions finished so far and exits with status 130; the accession being sketched leaves no partial file.

## sra_report

`sra_report` writes one HTML report, `sra_report.html` in `--output-dir` (default `sra_reports`), with the figures behind it in `sra_report.json`. Without groups, the report covers the quality of the datasets named by `--accessions-file`. With `--groups-file`, it adds a comparison of the groups and their statistical tests: a t-test for two groups and one-way ANOVA for more.

```bash
# Quality of a set of datasets
prepare sra_report --accessions-file accessions.txt --title "Project SRA Analysis"

# Treatment against control, reusing the profiles an earlier sra_profile run saved
prepare sra_report \
    --groups-file comparison_groups.json \
    --quality-profiles sra_quality_profiles
```

The groups file is a JSON object that maps group names to accession lists:

```json
{
  "Treatment_Group": ["SRR123456", "SRR123457"],
  "Control_Group": ["SRR789012", "SRR789013"]
}
```

Each accession is profiled once, by the same path as `sra_profile`, and that profile serves both the quality section and the comparison. An accession found in `--quality-profiles` is not profiled again. With neither `--accessions-file` nor `--groups-file`, every profile in that folder is reported. An accession without readable FASTQ files is left out with a warning.

The report opens in a browser unless `--no-open` is given. `--no-report` writes only `sra_report.json` and prints the results, and needs neither plotly nor jinja2. With `--groups-file`, scipy (the `analysis` extra) is checked before any dataset is profiled, also under `--no-report`. Without it the command exits with status 3 and writes nothing. A `report` analysis is recorded per accession with the total reads and bases, the GC content and the quality grade, so the profile columns of `results_table` are filled also for a project that never ran `sra_profile`.
