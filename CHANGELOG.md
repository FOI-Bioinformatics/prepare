# Changelog

All notable changes to PREPARE are documented in this file. Dates are in YYYY-MM-DD format.

## [0.10.2] - 2026-10-21

This release fixes the findings of the audit of the SRA download path of 2026-10-21
(`docs/superpowers/audits/2026-10-21-sra-download-audit.md`, H1-H4, M1-M8, L1-L7).

### Changed

- A truncated download in a project without a store is a failure: the files stay in place for
  inspection, the accession is retried and listed in `failed_accessions.txt`, the registry records
  `failed` with the `truncated` verdict, and the exit code is 1. `--accept-partial` keeps the earlier
  behaviour (the folder is published and recorded as downloaded). The shared store already treated
  the same data as incomplete, so both paths now agree (H1).
- Each compressed FASTQ file is read back after compression; a file with a different record count,
  or one that cannot be read, fails the download instead of being published (H4).
- An interrupt during compression leaves nothing published; the next run redoes the dump from the
  kept prefetch archive (M1).
- `--no-verify-downloads` also applies to a download into the shared store, whose copy is then
  recorded as unverified instead of being judged against the sidecar or the metadata XML (H3).
- A failed attempt removes its `<ACC>_temp` build folder; `--keep-temp` keeps it and logs its path (M5).
- An explicit `--max-workers` is capped at 10 with a warning; `--num-threads` and `--max-workers`
  must be positive, and so must fasterq-dump's `--threads` and pigz's `-p` in the subprocess
  allowlist (M8).

### Added

- `download_sra --keep-temp`.

### Removed

- `download_sra --verify-downloads`, which did nothing (verification is the default; use
  `--no-verify-downloads` to turn it off) (L3).
- The unused `DEFAULT_MAX_RETRIES`, `MAX_SUBPROCESS_TIMEOUT`, `FASTERQ_DUMP_BOOLEAN_FLAGS` and
  `FASTERQ_DUMP_INTEGER_FLAGS` constants (L5).

### Fixed

- A download worker that raised was reported as "cannot unpack non-iterable NoneType object"; the
  exception text is now the result, so its class (network or not), the exit code and the retry
  split are right (H2).
- The retry pass records no attempt for a worker stopped by an interrupt, and its back-off between
  rounds ends as soon as the run is stopped (M4).
- A project folder `fastq/<ACC>` that is not a store link made a store download fail only after the
  fetch, and the failure was retried; it is now refused before the fetch with a message naming
  `store_adopt` (M7).
- A prefetch leftover (`.sra.lock`, `.sra.prf`, `.sra.tmp`) beside a `.sralite` archive was handed to
  fasterq-dump, and a prefetch that exited 0 without an archive was reported as a fasterq-dump
  failure; both are handled (M2, M3).
- A failed forced redownload no longer records an empty file list while the earlier folder is still
  on disk (M6).
- A clean run removes an earlier `failed_accessions.txt`; the retry hint is logged once; with a store,
  `--redownload-truncated` refetches only the accessions the registry calls truncated (L1, L2, L4).
- The fasterq-dump version is probed once per process, and one `read_count` serves the store
  hand-off and the verdict merge (L6, L7).
- The fake `fasterq-dump` and `prefetch` of the process tests gained failure modes, and five
  process tests cover a truncated run, a tool failure fixed by the retry pass, an empty prefetch and
  an unpaired `--split-3` file.

## [0.10.1] - 2026-10-15

This release fixes the findings of the code audit of 2026-10-12 (H1-H5, M1-M16, L1-L11 and the
likely bugs it lists). PERMANOVA results from earlier versions should be recomputed (see Changed).

### Added

- `plot_map` (`maps` extra) draws sample locations from a metadata table.

### Fixed

Data safety and the shared store:

- A malformed line in an accessions file (`..`, `../x` or an absolute path) could become a folder name
  and be deleted by `download_sra`; accessions are now checked before any path is built (H1).
- `store_link`, `store_unlink` and `store_verify` refuse malformed accessions instead of acting on
  paths outside `fastq/` or the store (M1).
- A project bound to a store with a relative `--registry` looked stale from other folders, and
  `store_gc --include-stale` could remove a dataset it still linked to; registry paths are now stored
  resolved and a live link always protects its dataset (H2).
- `status --reconcile`, `store_link` and `store_unlink` run from a subfolder marked downloads missing
  or created a stray `fastq/` folder; they now act on the project that owns the registry (H4).
- `store_link` takes the dataset lock, so `store_gc` cannot remove a dataset while it is being linked
  or copied.
- `store_adopt --move` removes the project folder and links the store copy under one dataset lock; a
  held lock leaves the project folder in place.
- `store_gc` drops the catalogue row of a dataset whose folder is already gone instead of listing it
  as a candidate on every run.
- `status --reconcile` takes an accession it marks missing off the store's linked list.
- A registry written by a newer version is refused instead of being saved back as version 2, and
  unknown top-level keys are kept; a non-integer version is reported as a data error.
- `store_adopt` and `store_verify` skip a folder that is not named like an SRA accession (notes
  kept in `fastq/` or `<store>/sra/`) with one warning, instead of stopping before the others.
- Journal rewrites keep every line they do not parse byte-identical and flush the file to disk.
- Match CSV paths in the registry are recorded relative to the project root, so a registry moves
  with its project.
- `branchwater_search` anchors its output folder on the project root, like the other commands.

Statistics and results:

- PERMANOVA in `diversity_analysis` and `prepare_diversity` overstated F and computed p without the
  +1 correction, and samples without a group label were counted; see Changed (H3).
- PCA and t-SNE plots of `interactive_plot` colour and size by metadata again; before, every value
  was read as missing (H5).
- Branchwater containment thresholds are exact decimals, so maxima of exactly 0.3, 0.6 or 0.7 are
  counted at their threshold (M6).
- A Branchwater file with an unexpected header is refused instead of being copied as a match file (M7).
- `sra_profile --sampler head` and `sra_report` no longer stop on a zero-length read (M5).
- A homopolymer run at the 3' end of a read is counted (L1).
- `taxonomic_summary --min-abundance` is applied (L3).
- Selection ties are broken by run name, so the order does not depend on input order.

Downloads and retries:

- The retry pass of `download_sra` keeps `--force`, so a failed forced re-download is not reported as
  already present (M2).
- Results of the retry pass are classified like first-pass results: a store link found on retry is
  recorded as a link with its usage row, without download timing, and settled incomplete results
  are not fetched again (M3).
- A duplicated accession no longer holds a free-space reservation for the rest of the run (M4).
- A forced re-extraction that maps no reads removes the earlier FASTQ files and record (M9).

NCBI, GTDB and Branchwater clients:

- The NCBI taxonomy client paced 10 requests per second without an API key and 3 with one; the
  limits are now the other way round (M8).
- A single-run metadata download kept the whole efetch response, so a package listing another run
  first was parsed as that run; only the requested run is kept, and an unusable response writes no
  file (M10).
- An NCBI reply that is not XML is no longer kept as a metadata file.
- The SRA metadata search uses `[ACCN]` terms, allows two records per accession and warns about
  requested runs that are not returned.
- `SRA_Normalized_URL` comes from the SRA Normalized file, not the second file listed.
- The taxonomy cache keeps tax IDs as text, so an empty entry no longer turns IDs into `562.0` (L2).
  Validating an empty species list returns an empty table with the cache columns.
- A Branchwater CSV with a byte-order mark is read correctly (L4).
- A Branchwater search that returns no rows is not cached, so a later search queries the server
  again.
- GTDB taxon searches request further pages instead of stopping at the first 100 rows, and stop when
  a page adds nothing new; the search fallback ignores rows of a different accession.
- Genome downloads extract into a per-process folder, so concurrent runs do not collide, and report
  requested genomes that are absent from the NCBI zip.
- A `.tsv.gz` table is read as tab-separated.

Visualization and reports:

- The containment explorer escapes taxonomy family names in its HTML (M11) and fills its Location and
  Date columns from the metadata columns prepare writes (M12).
- The map plugin reads NCBI `lat_lon` values such as `38.98 N 77.10 W` and numeric values (M13).
- `sourmash scripts prepare_parse` reports unreadable match files with a non-zero exit and records
  screening like `parse_containment` (M14); the sourmash plugins print their errors to stderr (M15).
- `plot_metadata_counts` saves beside its input, warns when no `--save-format` is given and so
  nothing is written, keeps `limit` in the radar fallback and keeps a given title (L6); a single
  `--colors` value is accepted (L5). Figures are closed after saving.
- `project_report` shows the analysed stage. The CSV export of HTML report tables guards cells that
  start with a formula character, and table sorting is consistent for mixed numbers and text.

Command line:

- A broken pipe (`find_by_taxonomy ... | head`) exits 0 without an error line (M16).
- `find_by_taxonomy` output has no unnamed index column (L7).
- The global `-v` works for every command, including `store_status` (L8).
- `download_metadata --dry-run` creates no folder (L9).
- `--log-level` is case-insensitive; a config error no longer drops `--log-file`, `-q` or `-v` for
  `doctor`; a startup error under `--json` prints `{"error": ...}`.
- `PREPARE_ASSEMBLY_MEMORY` is checked at startup with the parser the assembly step uses, so a value
  such as `0G` is refused before any work. A `store` value of the wrong type in the config file is a
  configuration error, not a traceback.
- Numeric options refuse `nan` and `inf`.
- The `make` hint names `plotly_js_script` (L10), and the wheel carries `py.typed` (L11).
- The plugin registry refuses an abstract plugin class.

### Changed

- Accessions are checked at every command-line entry point, including `--accession` values; a
  malformed accession is refused with exit code 2, and for an accessions file the error names the
  file and line. Accessions files skip `#` comments and repeats.
- Invalid input (a `ValidationError`) exits with code 2, the usage-error code, instead of 1.
- Relative project folders are taken against the folder that holds the project registry, not the
  working directory, so a command run from a subfolder acts on the project above it.
- A stale project's live symlink now keeps its dataset out of `store_gc`; the reason line names the
  project with " (stale)".
- PERMANOVA uses the standard sums of squares (Anderson 2001) and the p-value
  (1 + #{F_perm >= F}) / (1 + permutations), so F and p differ from 0.10.0.
- Containment threshold labels are the half-up decimal rendering of k * step; some labels at
  half-way values differ from 0.10.0, the counts do not.
- Linking from a store needs a writable `<store>/locks`; a read-only store cannot be linked from.
- `store_status --details` replaces `--verbose`, which is kept as a hidden, deprecated alias.
- The CSS classes of the HTML reports use the `pre-` prefix.
- The default cache of `taxonomy_validation` is `taxonomy_validation_cache.csv`; an existing
  `taxonomy_cache.csv` is not read by default.
- The SRA metadata client paces 3 requests per second without an NCBI API key and 10 with one.
- Size accounting in `store_gc`, `store_link` and cleanup no longer follows symlinks.
- Genome accession files are de-duplicated.
- A stop during the retry pass of `download_sra` leaves the accessions not yet attempted unrecorded
  instead of recording them as failed.
- A forced re-extraction that maps no reads removes the FASTQ files of the earlier extraction.
- `--log-level` together with `-v` or `-q` is a usage error.
- `MissingGenomesError`, a `DataAccessError`, reports requested genomes absent from an NCBI zip,
  after the extracted genomes are registered.
- The containment explorer shows at most 5000 table rows, highest containment first, and notes how
  many were left out.
- `--help` no longer prints "(default: False)" for on/off flags.

### Internal

- One thread-safe `RequestPacer` (`data/ncbi_pacing.py`) for the three NCBI clients.
- Plugin base: abstract `FormatPlugin` and `VisualizerPlugin` interfaces, idempotent registration,
  built-in plugins loaded on first lookup, and third-party plugins discovered through the
  `prepare_plugins` entry-point group.
- Module splits: `registry_io` and `registry_scan` from `registry`, `read_extraction_export` from
  `read_extraction`, `sra_record` from the `download_sra` command, `metadata_fetch` from `metadata`;
  the SRA download tests are split into three files.
- The store journal is pruned by `store_gc` and compacted by `store_reindex`; `store_gc` loads each
  project registry once; a catalogue written by a newer version is refused.
- Shared helpers for folder sizes (`utils/sizes.py`), matrix preparation and figure writing.
- CI runs `make check`; the release workflow checks the tag against the package version
  (`scripts/check_version.py`).

## [0.10.0] - 2026-10-09

### Changed

- The tool is renamed from MetaQuest to PREPARE, because the old name clashes with an unrelated
  consumer product. The package, the command, the distribution and the repository
  (`github.com/FOI-Bioinformatics/prepare`) are now `prepare`; the exception base class is
  `PrepareError`; the sourmash plugin entry points are `prepare_parse`, `prepare_plot`,
  `prepare_diversity` and `prepare_taxonomy`.
- Every name the tool reads or writes follows the rename, with no fallback to the old name. A
  project created with 0.9.0 or earlier is used with 0.10.0 by renaming, by hand:
  `metaquest_registry.json` to `prepare_registry.json`; `.metaquest/` to `.prepare/`;
  each `.metaquest-assembly.json` marker to `.prepare-assembly.json`; `fastq/.metaquest-tmp/`
  to `fastq/.prepare-tmp/` (or delete it); in a shared store, `metaquest_store.json` to
  `prepare_store.json` in the store root only (the per-dataset `<ACC>.json` sidecars and `catalog.sqlite`
  keep their names); `~/.config/metaquest/config.toml`
  to `~/.config/prepare/config.toml`; and every `METAQUEST_<NAME>` environment variable to
  `PREPARE_<NAME>` (`METAQUEST_DATA` to `PREPARE_DATA`, `METAQUEST_PERF_SCALE` to
  `PREPARE_PERF_SCALE`). Rename these files before running any `prepare` command in the project:
  in a folder that holds only the old registry, `status` and `doctor` treat the folder as not yet a
  project, and `status --init` would build a new registry beside the old one.
- In each project that uses a store, after renaming its registry file, run
  `prepare store_init --data-root <store>` once, with the same `--project-name` if one was given
  before. This records the renamed registry in the store catalogue. Until then `store_status` lists
  the project as stale, and `store_gc --include-stale` treats its datasets as unused. Without
  `--project-name`, `store_init` sets the project name to the folder name.
- The `metaquest_version` key in the data of `project_report` is now `prepare_version`, and the
  `metaquest` key in the data of the python check of `doctor --json` is now `prepare`. The registry and
  the store sidecars do not hold a version key under either name.
- `make env` now creates a conda environment named `prepare`. An earlier `metaquest` environment keeps a
  `metaquest` console script that no longer imports; it can be removed.
- The hidden pointers for the commands renamed in 0.5.0 (`sra_stats`, `sra_profile_quality`,
  `sra_dashboard`, `sra_compare`) are unchanged.

## [0.9.0] - 2026-10-01

### Added

- A per-project run log. A command that keeps one appends a record of each run to
  `<project>/.metaquest/runs/runs.jsonl`: the command, its arguments with the NCBI API key masked, start and finish
  times in UTC, run time, exit code, MetaQuest version, host, process ID and a summary of counts and totals. Larger
  per-accession results go to `<run_id>.json` beside it, and the last 10 of these files are kept per command; an
  older run keeps its line without one. The record is written only where a project registry exists, also after an
  interrupt (exit code 130) or a `SystemExit` (its own code), a line that cannot be read (cut short, or holding a
  byte that is not UTF-8) is skipped with a warning, and a run log that cannot be written is a warning that never
  changes the command's exit code (`data/run_log.py`).
- The `run_log` setting (`METAQUEST_RUN_LOG`, `[runtime] run_log`, default true) turns the run log off.
- `download_sra`, `sra_profile`, `sra_report`, `results_table`, `extract_target_reads`, `store_verify`,
  `project_report`, `select_datasets`, `blacklist` and `status --init`/`--reconcile` add a record of each run to the
  run log. Per-accession results are kept in the run's detail file, one row per accession (per accession and genome
  for `extract_target_reads`), so `runs --diff` and `runs --accession` can compare them: download failures with
  reason and attempts and downloaded accessions with time and attempts, profile and report figures, mapped reads
  with breadth and mean depth (a sample skipped as already extracted with its recorded values and `skipped: true`),
  store verdicts and fixes, and selection ranks. Dry runs, `blacklist --list`, plain `status`, `doctor` and `runs`
  record nothing, and for `results_table`, `select_datasets` and `project_report`, `--no-record` also leaves the run
  out of the run log.
- `runs` (Environment group) reads the run log. It lists recent runs (`--limit`, default 20, 0 for all;
  `--command`), shows one run (`--show RUN`), compares two runs (`--diff RUN_A RUN_B`: every summary value with
  the difference for numbers, then the detail rows added, removed and changed) and follows one accession across
  runs (`--accession ACC`). A run is selected by its ID, a unique ID prefix, `latest` or `previous`; `--json`
  prints one document. A run without a detail file, because its command noted none or because it was pruned,
  is reported as not kept. `runs` records nothing and exits 1 without a run log or for a selector that matches
  no run or more than one.
- `project_report` (Environment group) writes one report of the project from its registry and run log:
  `project_report.md` and `project_report.json` always, and `project_report.html` when the `interactive` extra is
  installed (`--html auto|always|never`; `always` exits 3 without the extra and writes nothing). Sections: project,
  funnel, genomes, extractions, downloads, failed downloads (with reason and `attempts_total`, the attempts over all
  runs), timing (totals, medians, quartiles, 90th percentile, extremes), environment (the `doctor` checks without
  network; `--no-environment` leaves them out), outputs and the last 10 runs. Tables are cut to `--max-rows`
  (default 200, 0 for all) with the total shown. The export is recorded in the registry as `project_report` unless
  `--no-record`. It reads no FASTQ file and exits 1 without a registry; `--json` prints the written paths.
- `status` reports a cross-stage funnel: one text line with the number of datasets screened, selected, downloaded
  (with the size of the downloaded data and, when any download was timed, their recorded download time), extracted
  and assembled, and a `funnel` key in `status --json` that adds the excluded and failed counts, the recorded time
  of failed downloads (`failed_seconds`, apart from the downloaded datasets' `seconds`), the analysed stage, and the
  accession and genome pairs, time and assembled bases of the extraction and assembly stages. The extraction and
  assembly figures cover only the pairs that reached the stage (mapped reads or contigs above zero), unlike the
  totals of the timing line (`processing/project_funnel.py`).
- `results_table` has eight more columns after `download_verdict`: `quality_source` (which analysis supplied
  the quality columns: `profile`, `report`, `legacy`, or empty), the assembly's `assembly_largest`,
  `assembly_n90`, `assembly_gc_percent` (GC content in percent, two decimals), `assembly_contigs_ge_1kb` and
  `assembly_mean_depth_estimate` (empty unless the coverage mapping onto the contigs ran), and `assembly_dir`
  and `coverage_tsv` (paths relative to the project). All eight are empty where there is no assembly or no
  quality analysis was recorded.
- `status --export-tsv`'s `registry_extractions.tsv` has `coverage_tsv`, `n90`, `largest` and `assembly_dir`
  after its existing columns.
- `store_verify --json` prints one document: `root`, `checks` (which of `md5`, `spots`, `rescan` and
  `fix_state` ran), `counts` (datasets per verdict that occurred), `datasets` (one row per accession, with its
  state before the check and, under `--fix-state`, a `fix` entry with the action taken, `updated`,
  `unchanged`, `skipped-in-use` or `skipped-changed`) and `fixed` (the datasets `--fix-state` rewrote, with
  their state before and after). With no store configured it prints `{"error": ...}` and exits 1. In text
  mode, `--fix-state` prints one `<accession>: fixed (<before> -> <after>)` line per rewritten dataset.
- `download_sra --report-file` has two columns after the unchanged `accession,status,message,seconds`: `reason` (for
  a failed accession: `network`, `not-found`, `disk-full`, `insufficient-space`, `locked`, `interrupted` or
  `unknown`) and `attempts` (attempts in this run that started a download; a retry counts again, while an accession
  refused for space, cancelled by a disk-full abort, stopped by an interrupt before it started, given up after a
  lock wait, refused as a partial store copy or as an incomplete one whose refetches are used up, already present or
  linked from the store counts 0).
- `download_sra` writes `<fastq-folder>/download_run.json` after a run that reaches the download stage, with or
  without `--report-file` (a dry run writes none, and a run stopped by its preflight checks, such as a missing
  `fasterq-dump`, writes none unless the folder already exists): start and finish times, totals, failures by reason,
  each failed accession with its reason, attempts and message, whether the run was aborted (`disk-full` or
  `interrupted`), the exit code, host, settings and paths (`data/sra/run_report.py`).

### Changed

- `sra_report --groups-file` checks for scipy (the `analysis` extra) before any work, also under `--no-report`.
  A run without scipy used to profile every accession and then fail when the statistical tests ran or, under
  `--no-report`, succeed with every statistical test skipped; it now fails at once with exit code 3 and writes
  nothing.
- `sra_report`'s `report` analysis in the registry also carries `total_reads` and `total_bases`.
- The `download_sra` report CSV and `download_run.json` are also written for an interrupted run (`aborted:
  "interrupted"`, exit 130), listing the accessions that reported a result, and for a run whose final registry
  write failed. `failed_accessions.txt` is unchanged.

### Fixed

- `results_table` and `status --export-tsv` fill the quality columns (`total_reads`, `gc_percent`,
  `quality_grade`) from the `report` analysis of `sra_report` when `sra_profile` was never run for an accession;
  they used to stay empty. When both analyses exist, the newer one supplies the values and the other fills any
  it lacks; a value both lack is taken from the pre-0.5.0 `sra_stats` and `quality` analyses when the registry
  holds them.

### Removed

- The unused report generator (`metaquest/visualization/reporting.py` and
  `metaquest/visualization/templates/report_template.html`) and its three test files. No command called it;
  `project_report` replaces it.

## [0.8.0] - 2026-10-01

### Added

- One lookup for an accession's expected spot count, used by `download_sra` (with or without a store),
  `store_link`, `store_adopt` and `status --reconcile`: the registry's metadata block, then the store sidecar,
  then the metadata XML, then the count recorded with the previous download verdict (`data/sra/spots.py`,
  `expected_spots`).
- The `prefetch_max_size` setting (`METAQUEST_PREFETCH_MAX_SIZE`, `[runtime] prefetch_max_size`, default
  `100G`) sets `prefetch --max-size`, which was fixed at 100G.
- `results_table` has a last column, `download_verdict`: the accession's recorded download completeness
  verdict (`complete`, `truncated` or `unverified`), empty where none was recorded.
- `status --reconcile` names, inside the `drift` section of the report and of `--json`, every accession whose
  project link points into a data store that is not mounted (`store_unavailable`), every metadata block filled
  from the project's metadata XML (`metadata_filled`), every download verdict re-checked from `unverified`
  to a known outcome (`verdicts_rechecked`) and every assembly record dropped because it was dated before its
  extraction (`assemblies_dropped`, as `ACC/GENOME`). The four keys are always present in the `--json` drift
  object, as empty lists when there is nothing to report, so a reader need not test for them. The text report
  adds a warning naming the first five unmounted-store accessions, and a line each for the re-check, fill and
  dropped-assembly counts, only when there is something to report.
- Store sidecars can carry an optional `refetch` record (see the refetch limit under Fixed); a sidecar without
  it is read and written unchanged (schema 1).

### Changed

- `extract_target_reads --assemble` redoes an assembly when the extracted reads (file names and sizes), their
  extraction date, the megahit preset, the minimum contig length or the k values differ from those it was built
  from. These are recorded in a marker in the assembly folder, `.metaquest-assembly.json` (which also keeps the
  megahit version, not compared), and in a new optional `inputs` entry of the registry's assembly record; a
  registry without it is read and rewritten unchanged. An assembly that still matches is reused and reported
  as reused. An assembly folder written by an earlier version has no marker; it is accepted, and given a marker
  whose megahit version is recorded as unknown, only when its registry record matches the preset and minimum
  contig length and is not older than the extraction. A value missing from the record matches any setting, but
  a recorded preset, `default` included, matches only the same `--assembly-preset`.
- Re-running a targeted read extraction for a sample drops that sample's recorded assembly from the registry,
  because the assembly was built from the reads the extraction replaces; the assembly folder on disk is left
  in place.
- megahit writes into a hidden staging folder beside `<genome_id>_assembly`, and the result replaces the
  previous folder by rename, so a failed or interrupted assembly leaves the earlier one in place. `--force` no
  longer removes the folder before megahit runs; a failed forced or stale redo keeps the previous folder on
  disk without a registry record, and a later run records it again from its marker when the marker still
  matches. megahit exiting without writing `final.contigs.fa` is reported as an error.
- A sample whose assembly fails (megahit, or the coverage mapping onto its contigs) no longer stops the
  remaining samples. The failed samples are listed in one error line and `extract_target_reads` exits 1. A
  sample being extracted or assembled by another process is skipped and reported, and does not fail the run.
  Each sample's assembly step loads the project registry once, also when its assembly is reused (about one
  second per sample on a registry of about 15,000 datasets).
- `download_sra --keep-sra` keeps a cached `.sra` archive only after a complete or unverified download; a
  truncated archive is removed, and `--force` or `--redownload-truncated` always fetches the archive again
  instead of reusing a cached one.
- In a project without a store, `fasterq-dump` scratch files go to `fastq/.metaquest-tmp/<accession>_fqtmp`
  (removed when the accession finishes) instead of the system temporary directory; the free-space check and
  `doctor` measure that location. Set `--temp-folder` or `METAQUEST_TEMP_FOLDER` to use local scratch space,
  for example on a cluster node whose project folder is on a network filesystem.
- `download_sra --no-verify-downloads` skips the check of the project's own downloads. A copy in the shared
  store is always judged against its recorded spot count, and a partial copy is still refused without
  `--accept-partial`.
- `status --init` fills every metadata block it bootstraps from an `<accession>_metadata.xml` on disk with a
  positive spot count (spot count, size, md5, assay type, organism, collection date, library layout and
  strategy, platform) in the same write that creates the registry, rather than leaving every field empty until
  the next `status --reconcile`.
- `status --reconcile` fills a metadata block recorded without a spot count from the project's
  `metadata/<accession>_metadata.xml` before checking verdicts, so one run fills both. An XML that records no
  spot count leaves the block unchanged, so repeated runs write nothing new.
- `status --reconcile` re-checks an `unverified` download verdict once NCBI's spot count is known. A plain
  project download is re-counted from its mate-1 file and its file of unpaired spots. A link into the shared
  store is compared with the read count its sidecar records, and no FASTQ is read. The first reconcile after
  upgrading reads the mate-1 file of every `unverified` plain download that has a spot count once; later runs
  read nothing again.

### Fixed

- A partial store download is no longer linked into the project without `--accept-partial`. The copy is
  published to the store for a later `--resume-partial` refetch, any project link into the store for that
  accession is removed (a real project folder is never touched), and the accession is reported as
  `incomplete: ...`, which counts as a failure (exit 1).
- A store download whose metadata XML records no spot count is verified against the registry's count, and a
  store copy recorded as unverified but short against a known count is treated as partial and refetched.
- A store refetch that comes back less complete no longer replaces a better published copy, and a forced
  refetch that comes back worse no longer links an earlier copy that is unverified but short against the
  expected spot count without `--accept-partial`.
- A partial store copy is no longer downloaded again on every run. After two refetches in a row that gain no
  reads (`STORE_PARTIAL_REFETCH_LIMIT`; a refetch that again returns no reads at all gains none), it is linked
  with `--accept-partial` or refused with a message naming `--force`, the only way to fetch it again. The
  counter is kept in the store sidecar, so it is shared by every project that uses the store.
- A store copy that `download_sra` keeps but does not link (a result starting `incomplete:` or
  `partial in store;`) is a settled outcome and is no longer retried by the retry pass of the same run; before,
  `--max-retries` fetched a partial copy a second time and used up one of its refetches.
- The free-space check counts a store copy that is unverified but short and due for refetch, and exempts a
  copy that will not be fetched again.
- Linking a dataset from the shared store, and `status --reconcile`, no longer replace a recorded `complete` or
  `truncated` download verdict with `unverified` when no spot count is known; the verdict is recomputed when the
  read count and the spot count are both known. A plain redownload that comes back `unverified` no longer
  replaces a recorded `truncated` verdict; a redownload that comes back `complete` still replaces it. A
  `truncated` record without a spot count is therefore fetched again by every `--redownload-truncated` run;
  `download_metadata` for that accession gives it a count to be judged against.
- `store_link` and `store_adopt` record the same verdict as `download_sra` for a store dataset: the sidecar's
  read count judged against the project's spot count (registry metadata, sidecar, the project's and then the
  store's metadata XML, previous verdict), merged with the verdict on file. Before, they recorded the sidecar's
  verdict as it stood, so a relink turned a recorded `truncated` verdict into `unverified` and the copy was
  then extracted. `store_link` also refuses, without `--accept-partial`, a copy whose sidecar is `unverified`
  but short against that spot count, and records its sidecar as `partial`.
- In a project linked with `--link-mode copy`, a store copy that a later `download_sra` run finds short (a
  result starting `incomplete:` or `partial in store;`) is recorded with the store sidecar's verdict, so its
  project copy, whose own sidecar still reads complete, is skipped by `extract_target_reads` through the
  registry's `truncated` verdict.
- `download_sra` takes an accession's expected spot count, when the registry's metadata records none, from the
  store sidecar (with a store), then the metadata XML (the project's `metadata` folder, then the store's), then
  the previous verdict, so a plain download is verified rather than recorded `unverified`. Before, the previous
  verdict's count was used ahead of the store sidecar.
- `status --reconcile` drops an assembly record dated before its extraction record, which 0.7.0 carried forward
  when a sample was extracted again, so `status` and `results_table` no longer report contig statistics of an
  assembly built from other reads. The assembly folder is left on disk, and the next
  `extract_target_reads --assemble` builds the assembly again.
- `results_table` no longer stops on a hand-edited registry whose download verdict is not a mapping; the
  verdict reads as empty.
- An accession found on disk without a download record, as a run killed after its files were in place leaves
  it, is verified against its spot count before it is recorded, and is recorded `unverified` when no count is
  known, rather than with no verdict. A truncated one is reported with a pointer to `--redownload-truncated`.
- `extract_target_reads` skips a sample whose store copy is recorded as partial or failed, whatever its
  registry download verdict says. Before, it skipped only a registry verdict of `truncated`, so a store copy
  found short after its download was recorded as `unverified` was extracted and assembled. The skip line gives
  the reason, for example `store copy partial (R of S spots)` or `store copy failed verification: <error>`.
  `--allow-truncated` extracts such a copy anyway. A project download outside the store recorded as
  `unverified` is still extracted, because only a store sidecar records a copy as short without reading its
  files. Only the selected samples' sidecars are read, so a malformed sidecar of another sample does not stop
  the run.
- A staging folder left by a killed assembly, or an assembly folder an interrupted megahit left without
  `final.contigs.fa`, no longer blocks the next run: the staging folder is removed and the sample is assembled
  again (with a warning), where a contig-less folder used to stop the run with an error. An earlier folder
  moved aside by a publish that was cut short is put back (or the finished new one, when its marker was
  written); an assembly folder removed after a finished publish is not brought back.
- `status --reconcile` no longer marks a download missing when its project link points into a data store that
  is not mounted. It reports such accessions as store unavailable and logs one warning. A dataset removed from
  a mounted store is still marked missing, and links are never removed.
- A store-linked dataset that an earlier reconcile marked missing is re-recorded as linked from the store,
  keeping `source` and `store_name`, once its project entry is again a symlink into the store; a real folder in
  its place is recorded as a plain download.
- `store_gc` reports, and with `--yes` removes, a bare `tmp/<ACC>` staging folder a killed download left
  behind, once that accession's dataset lock is no longer held.
- `store_link --link-mode copy` (and any copy-mode link) can refresh an earlier copy of the same accession: the
  earlier copy is moved aside only after the new copy is complete, and put back if the swap does not happen,
  also on an interrupt. A project folder without the store's `<ACC>.json` is still refused. A copy-mode link
  first checks that the project's filesystem has room for the copy, and a shortfall or a failed copy is
  reported as an error naming the accession.
- `store_adopt` writes a dataset's sidecar into the staging folder before publishing it, so an adoption killed
  at any point no longer leaves a store folder without a sidecar.
- `store_reindex` no longer refuses to run because a store folder has no sidecar. It rebuilds one from the file
  names and sizes, with state `failed` and an error that says to run `store_verify --rescan --fix-state <ACC>`,
  and catalogues it; such a dataset is not linked without `--accept-partial`, and when no project uses it,
  `store_gc` offers it for removal. A folder whose lock is held is skipped with a warning, and so is a folder in
  `sra/` whose name is not an SRA accession. An unreadable sidecar still stops the reindex.

## [0.7.0] - 2026-10-01

### Added

- Runtime settings with one order of precedence: a flag, then a `METAQUEST_<NAME>` environment variable,
  then the `[runtime]` table of `~/.config/metaquest/config.toml` (or
  `$XDG_CONFIG_HOME/metaquest/config.toml`), then a built-in default. They cover the external tool
  timeout, the download worker cap, the lock waits and stale limits, the NCBI email address and API key,
  the temporary folder, logging, the free-space floor and megahit's memory; `docs/configuration.md`
  lists each one with its type, default, variable and flag. A value that does not parse stops the
  command with exit code 3 and names the variable or config key; an unknown `[runtime]` key is a
  warning. At DEBUG a run logs every setting with its value and source.
- Logging options on every command, before or after the command name: `--log-file PATH` appends every
  line at INFO or above to a file, each stamped with host name and process ID and with the full traceback
  of a failure; `-q/--quiet` (warnings and errors only on the console) and `-v/--verbose` (DEBUG), which
  cannot be combined and override `--log-level`. At DEBUG a run logs the version, its command line (the
  `--api-key` value hidden), host, process ID and SLURM job and array task IDs. `METAQUEST_LOG_HOST=true`
  adds host and process ID to console lines. See "Logging" in the README.
- Progress summaries for `download_sra`, `download_metadata` and `extract_target_reads`: one line every
  `--progress-every N` items (default 50; also `METAQUEST_PROGRESS_EVERY` or `progress_every` in
  `[runtime]`) and at least every 5 minutes, such as `download_sra: 150/2000 done (148 ok, 2 failed),
  3.1/min, about 9 h 57 min left`, and a closing line with the totals and the time taken (for
  `download_sra`, after the retry pass, so it counts the accessions a retry fetched).
- Free-space check for `download_sra`: before each accession starts, the filesystems it writes to (FASTQ
  or store `tmp` folder, `fasterq-dump` temporary folder, `.sra` cache) must have room for it, counting
  the downloads already running. An accession with a registry run size needs 8 times that size for its
  temporary files (assumed equal to the uncompressed output, not measured), 10 times for the FASTQ folder
  (uncompressed and gzip files together), plus the size for the cache; the factors come from seven runs
  of the crispatus test store, whose uncompressed FASTQ was 6.98 to 7.73 times the `.sra` size and
  gzip-compressed FASTQ 1.49 to 1.66 times. An accession of unknown size needs `--min-free-gb` (default
  10; `METAQUEST_MIN_FREE_GB`, `[runtime] min_free_gb`); 0 turns the check off. An accession that does
  not fit while others run waits for them to release their space; one that would not fit even alone
  fails with `insufficient-space: not enough free space on <mount>: ...`, and the others continue; the
  retry pass does not try it again. Only a tool's own out-of-space error stops the pass (see Fixed).
  Before the first download a warning names a filesystem that the downloads of known size may not fit
  on together; it does not stop the run.
- `--assembly-memory` on `extract_target_reads` (also `METAQUEST_ASSEMBLY_MEMORY`, `[runtime]
  assembly_memory`): `auto`, the default, passes megahit `--memory` as 90% of the memory limit detected
  for the job (cgroup v2 or v1, else `SLURM_MEM_PER_NODE` or `SLURM_MEM_PER_CPU`) and omits it when
  none is found; a size such as `32G` is passed in bytes; a fraction is passed as it is and applies to
  the node's whole memory. A whole number of bytes below 1M is refused as a size missing its unit.
- The default number of parallel downloads is the CPUs available to the job (the affinity mask, which a
  SLURM cpuset limits, else `SLURM_CPUS_PER_TASK`, else the CPU count) divided by `--num-threads`, at
  most 4; `METAQUEST_MAX_WORKERS_CAP` (or `[runtime] max_workers_cap`) changes the cap, which the
  `--max-workers` help now names. `seqkit stats` in the store statistics uses at most as many threads.
- `metaquest doctor` (new Environment group; `make doctor`): checks the Python and MetaQuest versions,
  every external tool with its path, version and oldest supported version, that the config file parses
  (with each runtime setting and its source), the shared data store (marker, write access, free space),
  free space at the project, temporary folder and `.sra` cache against `min_free_gb`, CPUs, memory
  limit and SLURM variables, and the nearest project registry; `--network` adds NCBI and Branchwater
  (10 s each), `--for COMMAND` turns a missing tool that command needs into a failure, `--json` writes
  one JSON document. Exit code 0 without a failed check, 3 with one. It still runs, and reports the
  error once (as the config check), when `config.toml` does not parse. A tool that exits non-zero on
  its version flag (a broken install, such as a missing `libcrypto`) is not runnable: its first error
  line is reported, not a version read from that error; it is a failed check for a tool `--for` needs
  and a warning otherwise. `--for` takes the command names `metaquest --help` lists. `SRAMetadataClient` has
  `close()` and works as a context manager. See "Checking the environment" in the README.
- Timing: the registry records when each download, extraction and assembly started and how many
  seconds it took (`started`, `seconds`; absent in registries written earlier). `download_sra
  --report-file` has a trailing `seconds` column, `results_table` ends with `download_seconds`,
  `extraction_seconds` and `assembly_seconds`, `status --json` has a `timing` block (counts, totals
  and medians), the text report one timing line, and `status --export-tsv` the same columns. An
  extraction's time runs from the previous sample's checkpoint, so the first sample's includes reading
  the containment table and building the minimap2 index; an assembly's time is the megahit run alone,
  without the contig summary and the coverage mapping. A failed download attempt's time is recorded
  too and counted in the totals; a time left on a `missing` or `skipped` download block is not, and
  `store_link`, `store_unlink` and a skip clear it.
- `pypi`, a second job in the release workflow (`.github/workflows/release.yml`), publishes the
  distribution the `build` job already produces to PyPI via PyPI's trusted-publisher mechanism (no
  token in the repository). It only runs once the repository variable `PYPI_PUBLISH` is set to
  `true`; see `docs/packaging.md` for the checklist before that and a bioconda recipe template.

### Changed

- Exit codes: 3 for a configuration problem (`ConfigurationError`: a missing optional package, an
  external tool that is missing or too old, a missing NCBI email address, a malformed `config.toml` or a
  setting that does not parse, a log file that cannot be opened), 4 for a retryable failure (a wait for
  the registry or catalogue lock that reached its limit; an NCBI request of `sra_info` that could not
  connect, timed out or got HTTP 429 or 5xx after its retries; a `download_sra` run whose every failure
  was a network one; a `download_metadata` run that could not reach NCBI for any accession), 130 for an
  interrupt, 1 for any other failure. `download_metadata` still logs each accession NCBI did not return
  and exits with 0 when it fetched any or a failure was not a network one; `sra_validate` keeps 4 for a
  registry lock wait that gave up. See "Exit codes" in the README.
- A command failure is one line on the console unless it is at DEBUG; the traceback now always goes to
  the log file when there is one. The hint "Use --log-level DEBUG for full traceback" follows a failure
  only when there is no log file and the console is not at DEBUG.
- Error messages of the store commands start with the command name (`store_status: ...`).
- The per-accession INFO lines of `download_sra` (starting, skipping, the temporary folder, the files
  downloaded, the store link or copy, a download or a retry that succeeded), of `download_metadata` for
  each NCBI request, and of `extract_target_reads` for each sample, are now logged at DEBUG, replaced at
  INFO by the progress summaries. `--progress-every 0` logs them at INFO again. Still at INFO: the
  first "prefetch not found on PATH" line of a run (later accessions log it at DEBUG), and the lines
  about an unusual event for one accession: a download interrupted or redone with `--force`, a truncated
  archive removed for the next attempt, a wait for running downloads to free space, a wait for a lock
  another run holds.
- External tool timeout is now configurable and, by default, unlimited: `run_secure` used to give every
  external tool (`fasterq-dump`, `prefetch`, `minimap2`, `samtools`, `megahit`) a fixed one-hour limit
  (`MAX_SUBPROCESS_TIMEOUT`) even when a caller passed no timeout at all, and `timeout=0` was silently
  treated as "use the one-hour limit" instead of "no limit". `--timeout SECONDS` on `download_sra` and
  `extract_target_reads` (also `METAQUEST_TIMEOUT` or `[runtime] timeout` in `config.toml`) now sets it;
  0, the default, means no limit. A `--version` probe always uses a fixed 30-second timeout regardless of
  this setting, so a hung tool cannot stall version detection.
- Two runtime settings are checked against each other: a `lock_heartbeat` not below
  `dataset_lock_stale`, or a `registry_lock_stale` of 5 s or less, stops the command with exit code 3.
- One table of external tools (`metaquest/utils/tools.py`) with the oldest supported versions: sra-tools
  3.0, minimap2 2.17, samtools 1.10 (for `samtools coverage`), megahit 1.2.9, pigz 2.4 and
  ncbi-datasets-cli 16. `download_sra`, `extract_target_reads`, `genome_download` and `genome_prepare`
  check the tools they need before any work starts and exit with 3 (was 1 for `download_sra` and
  `extract_target_reads`; `genome_download` and `genome_prepare` had no check and failed inside the
  first `datasets` call) when one is missing, older than its floor or not runnable (its version probe
  exits non-zero), listing every problem with its `conda install` command. A version that cannot be read
  for another reason (a probe timeout, no version number printed) is a warning. The optional tools
  (`prefetch`, `pigz`, `seqkit`) are checked only by `doctor`. The check runs each required tool's
  version probe at the start (up to 30 s per tool).
- `environment.yml` pins the same floors as that table (`sra-tools>=3.0`, `minimap2>=2.17`,
  `samtools>=1.10`, `megahit>=1.2.9`, `pigz>=2.4`, `ncbi-datasets-cli>=16`; checked by
  `tests/test_environment_pins.py`); `seqkit` stays commented out and has no floor.
- `--email` of `download_metadata`, `sra_info` and `validate_taxonomy` is no longer required when an
  address is set in `METAQUEST_NCBI_EMAIL` or `ncbi_email` in `[runtime]`; `--api-key` falls back to
  `METAQUEST_NCBI_API_KEY`, then `NCBI_API_KEY`, then `ncbi_api_key`. Without any email address the
  command stops with exit code 3 and says how to give one.
- `--lock-wait`, `--temp-folder` and `--log-level` fall back to their environment variable and
  `[runtime]` key when not given; the lock limits of the registry, the store catalogue and the dataset
  locks can be set the same way (see `docs/configuration.md`). The built-in values are unchanged.
- A flag that sets a runtime setting is checked like its environment variable: `--email nope` or a
  negative `--min-free-gb` (which used to turn the free-space check off) stops the command with exit
  code 3. The error for a malformed `config.toml` quotes the offending line; an invalid API key value is
  shown as `(hidden)`. The sourmash plugin's `metaquest_taxonomy` takes the email address and API key
  from the settings too, so `--email` is no longer required there.
- The NCBI taxonomy client, the GTDB client and `SRAMetadataClient` (`sra_info`) share one retrying
  HTTP session (`metaquest.utils.http.retrying_session()`); `SRAMetadataClient` used to call
  `requests.get` without retries, and now retries a connection failure, HTTP 429 or 5xx with backoff
  before it raises `NetworkError`.
- `DEFAULT_MEMORY_LIMIT_GB`, `MAX_FILE_SIZE_MB`, `DEFAULT_PLUGIN_TIMEOUT`, `ERROR_MESSAGES`,
  `SUCCESS_MESSAGES` and `DEFAULT_MAX_WORKERS` are removed from `metaquest/core/constants.py`: none was
  read anywhere outside its own definition. The console log level choices and default come from
  `constants.LOG_LEVELS` and `DEFAULT_LOG_LEVEL` alone.

### Security

- NCBI and store metadata XML is now parsed with a hardened lxml parser
  (`metaquest/utils/xml.py`'s `SAFE_PARSER`/`parse_xml_file`, used by `data/metadata.py`'s
  per-file parse): entity resolution, network access and external DTD loading are refused, closing an
  XXE (external entity reading a local file or the network) risk that lxml's defaults leave open; a
  "billion laughs" document is rejected by libxml2's entity-amplification limit. The stdlib
  `xml.etree.ElementTree` parse sites elsewhere (`data/metadata.py`'s batch parse,
  `data/sra_metadata.py`, `data/taxonomy.py`), which parse small responses fetched directly from NCBI,
  are unchanged: ElementTree does not resolve external entities (an undefined-entity error) or fetch a
  DTD, so XXE does not apply to them, and expat 2.4.1 and later rejects entity amplification such as
  billion laughs (tested).
- The `--api-key` value is replaced by `***` in the command line logged at DEBUG, also when the flag is
  abbreviated (`--api-k`), and shown as `(set)` in the settings list and in `doctor`.

### Fixed

- `sra_info` exits with 4, not 1, when NCBI cannot be reached after retries. Two faults hid this: its
  `except Exception` turned every failure, including a retryable `NetworkError`, into exit code 1; and
  `SRAMetadataClient.get_sra_metadata` caught `NetworkError` per batch (it is a `DataAccessError`),
  logged it and returned the partial results, so the error never reached the command. It now re-raises
  `NetworkError` before its per-batch handling.
- megahit under cgroups (a SLURM job, a container) sized its memory from the node's total, not from the
  job's limit, and could be killed for exceeding it; `--assembly-memory auto` now passes the job's limit.
  On a cgroup v1 host the limit is read from the job's own group (`/proc/self/cgroup`), not only the
  controller root, and a `--mem-per-cpu` job is sized from `SLURM_MEM_PER_CPU` times its CPUs.
- `--timeout` and a termination signal did not stop a tool that runs its work in a child process holding
  the output pipes (the conda `megahit` wrapper runs `megahit_core`): the wrapper was killed and the
  child ran on until it finished. Each tool now runs in a process group of its own, and the whole group
  is signalled (Linux and macOS).
- A disk that filled up during the first download pass did not stop it: every remaining accession was
  still started and failed in turn, and the retry pass tried again. The first result in which a tool
  reports that it ran out of space now cancels the downloads not yet started (recorded as `disk-full:
  not attempted`, as the retry pass already did), lets the running ones finish, and skips the retry pass.
- `timeout=0` passed to `run_secure` meant a one-hour limit instead of no limit (see Changed).

### Documentation

- `docs/hpc.md`: running on a SLURM cluster, with a pre-flight `doctor` check, an array job template for
  `download_sra`, what the walltime signal flushes and how a rerun resumes, exit codes with
  `--dependency=afterok` and a resubmission loop on 4, a shared store on a group filesystem, NFS caveats
  (lock files, clock skew, temporary files, log files), megahit memory and one extraction job per genome.
- `docs/configuration.md`: every runtime setting with its type, default, environment variable, flag and
  config key; the `[store] data_root` key; the two orders of precedence and why they differ.
- `docs/packaging.md`: the PyPI and bioconda release order, the checklist before turning on
  `PYPI_PUBLISH`, and a bioconda recipe template.
- README: "Configuration", "Running on a cluster", one exit-code table, "Logging", "Checking the
  environment", and the external tool table with the oldest supported versions.
- `docs/ARCHITECTURE.md`: the configuration section now describes `core/settings.py` and
  `store/resolve.py`; new entries for the tool table, resources, progress, XML and HTTP helpers, the
  free-space guard and the timing module; the exit-code and logging policies. The catalogue's journal mode
  is corrected to `DELETE`.

### Upgrade notes

- A script that treated every non-zero exit as the same failure keeps working. A script that tested for
  exactly 1 should also accept 3 and 4; 4 means the same command may succeed if run again later, and 3
  means the environment or configuration needs fixing first.
- An external tool is no longer stopped after one hour. A job that relied on that limit to end a stuck
  tool should set `--timeout` (or `METAQUEST_TIMEOUT`); under a scheduler the walltime is the limit.
- A script that searched the INFO log of `download_sra` or `download_metadata` for one line per
  accession no longer finds most of those lines: run with `--progress-every 0` to keep them at INFO, or
  with `--log-level DEBUG`. Warnings and errors naming an accession are unchanged.
- `--email` can be dropped from scripts once `METAQUEST_NCBI_EMAIL` or `ncbi_email` is set; a script
  that passes it keeps working.
- Older tools are refused with exit code 3: sra-tools before 3.0 by `download_sra`, minimap2 before 2.17
  and samtools before 1.10 by `extract_target_reads` (megahit before 1.2.9 with `--assemble`), and
  ncbi-datasets-cli before 16 by `genome_download` and `genome_prepare`. `environment.yml` installs
  versions that pass; run `metaquest doctor` to see what is installed.
- `~/.config/metaquest/config.toml` is now read by every command, not only for the store root: a file
  that is not valid TOML, or a `[runtime]` value that does not parse, stops every command except
  `doctor` with exit code 3. Run `metaquest doctor` to see which key is at fault.
- For library users: `MAX_SUBPROCESS_TIMEOUT` is kept as an alias of `DEFAULT_SUBPROCESS_TIMEOUT` (0);
  `metaquest.utils.security.missing_tools` is removed (use `metaquest.utils.tools.require_tools` or
  `probe_tool`); the five constants listed under Changed are removed.

## [0.6.0] - 2026-09-30

### Added

- `registry_update(path, mutation)`: loads the registry, applies a mutation to it and writes it back, all
  under one lock, returning the mutation's result. Several commands were migrated onto it; see Changed.
- Process-level concurrency tests (`tests/test_concurrency_processes.py`, `tests/helpers_processes.py`):
  spawn the real CLI as a subprocess with fake tools on `PATH` and exercise same-accession contention for
  a plain project and for a store, registry contention during a download, two `SIGTERM`s during a
  download, a tool that ignores `SIGTERM`, a `SIGKILL`ed lock holder, and `store_gc` against a dataset
  another run is downloading.
- An atomic-writes gate in `make check` (`scripts/check_atomic_writes.sh`): fails on a direct
  `.write_text`, `.write_bytes`, `.to_csv`, `copy2`, `copyfile`, `write_html` or `savefig` call, or an
  `open(path, "w...")`/`"x..."`, `path.open("w...")` or `open(path, mode)` call, outside the small allowlist
  in `scripts/atomic_writes_allowlist.txt`. Figure writers (matplotlib `savefig`, plotly `write_html`) are
  on the allowlist: a figure is regenerable and overwritten whole, and no later step reads it.
- A `multiprocess` pytest marker for tests that start real subprocesses; it runs as part of the default
  suite.
- Locks for a shared minimap2 index build and for per-sample read extraction
  (`<index>.lock`; `<output-folder>/.locks/<ACCESSION>.<GENOME_ID>.lock`), so two overlapping processes
  building the same index or extracting the same sample no longer duplicate the work; a sample already
  being extracted elsewhere is reported skipped rather than mapped twice.

### Changed

- One lock mechanism (`metaquest/utils/lockfile.py`) now backs the registry lock, the store catalogue
  lock, every per-accession dataset lock, the new plain-project download lock, and the new index-build
  and per-sample-extraction locks. Every lock file carries a holder record (pid, host, start time, a
  random token) refreshed by a heartbeat thread. A lock judged stale is reclaimed only under a second,
  short-lived guard file, so two waiters judging one lock stale at the same instant still produce exactly
  one winner, and a holder that stalls past its stale window (a slow write, a sleeping laptop) is not
  destroyed as long as its heartbeat keeps running. A holder whose process has since died on the same
  host is taken over at once instead of waiting out the stale window.
- The registry lock's stale window is now 120 seconds (was 30), refreshed by the heartbeat every 5
  seconds; the wait before a writer gives up is unchanged at 30 seconds. The store catalogue lock waits
  up to 60 seconds and is itself judged stale after 120 seconds. Every lock's "gave up waiting" message is
  now one format across every lock kind, naming what is locked, the holder's pid and host, and since when.
- A plain project's downloads (no shared store) now take a per-accession lock the same way a store does,
  so two runs on one project no longer both write into `fastq/<ACCESSION>_temp` and delete each other's
  in-flight output. Each download is built, verified and compressed under
  `fastq/.metaquest-tmp/<ACCESSION>` and published into `fastq/<ACCESSION>` with one rename.
  `--lock-wait` now bounds this wait too, the same as it already bounded the store's lock. The new
  in-lock "already exists" check this adds (so a second process finding a finished download does not
  re-fetch it) honours `--redownload-truncated` the same way the existing pre-lock check already did: an
  accession the registry's last verification marked `truncated` is still re-downloaded once the lock is
  held when that flag was given, rather than being reported already present.
- `store_gc` treats a dataset used (linked or downloaded) within the last day as in use even when nothing
  currently holds its lock, and re-checks each candidate's lock, usage, and link state again immediately
  before removing it, not only at classification time; a dataset kept for either reason is listed under
  `in_use` with the specific reason. `store_verify --fix-state` re-reads the sidecar under the dataset
  lock immediately before writing it back, and skips the write if another process changed it meanwhile or
  if the dataset is locked.
- A store catalogue write that fails after a dataset has already been published, but before it is
  linked into the project, no longer fails the download: the result is reported as a success with
  "; catalogue pending; stored" in its message, the dataset is still linked, and `store_reindex` repairs
  the missing catalogue row from the sidecar already on disk.
- Every command now handles `SIGINT`, `SIGTERM`, and `SIGHUP`. The first turns into a `KeyboardInterrupt`
  so a `finally` block or a registry batch can flush what has been done; a second signal is logged rather
  than raised, so it cannot cut a final write short; a third abandons the write, terminates any running
  tool, and exits with status 130. A `store_adopt --dry-run` scan stopped this way also exits 130, the
  same as a real run, rather than completing its report.
- `import metaquest` no longer calls `setup_logging()`; a library host now configures its own logging, or
  calls `metaquest.utils.logging.setup_logging()` itself for console or file output. The CLI is
  unaffected, since `metaquest.cli.main.main()` already calls it explicitly. A second call to
  `setup_logging()` (for example a host application reconfiguring it) now removes only the handlers it
  installed itself, leaving a host's own handlers, or pytest's `caplog` handler, untouched.
- Seven commands that used to load the registry, do their work, and write it back under a lock covering
  only the write (`select`, `blacklist`, `download_metadata`, `parse_metadata`, `sra_profile`,
  `sra_report`, `sra_validate`, and `status --init`/`--reconcile`) now read a snapshot for their own work
  and record their outcome through `registry_update` or a batch, so the file they write is re-read under
  the lock and never overwrites a concurrent `download_sra`'s committed outcome. `status --reconcile` now
  splits its filesystem scan (no lock held) from applying the resulting plan (under the lock), so read
  counting no longer runs while the registry lock is held.
- As part of that same migration, `download_metadata` and `parse_metadata` now log and drop a single
  accession's failing `record_metadata` call instead of letting the exception end the command before
  anything is saved; the command still finishes and exits 0, with every other accession's metadata
  recorded, and its summary line ("Recorded metadata for N accession(s) in the registry; M dropped")
  counts the dropped records. A user relying on either command failing loudly on a bad record should
  watch the log instead.
- `status --reconcile` checks an untracked accession's folder again before recording it as downloaded,
  and uses a spot-count verdict only for the download it was computed on, so a folder removed, or a
  download redone, between the scan and the apply is not recorded from the stale scan.
- `parse_containment` reads a parsed table from disk before taking the registry lock, not under it.
- Every tabular and text output file metaquest writes is now written atomically, through a uniquely named
  temporary file replaced into place with `os.replace`: CSV and TSV tables, the registry, store sidecars,
  the minimap2 index and its record, per-accession metadata XML (and its copy into a store's `metadata/`
  folder), extracted genome FASTA files, the SRA HTML report, the containment HTML report and the
  explorer page; the registry and sidecar writes also call `fsync`. Figures (PNG, PDF, SVG and the
  interactive plotly HTML plots) are still written directly. Temporary file names are dot-prefixed and
  carry the hostname, pid, and a random token, so they no longer collide between two SLURM nodes writing
  into the same folder, and stay invisible to every folder listing. Text metaquest writes is now explicit
  UTF-8 rather than the platform's default locale encoding.
- `--lock-wait` and `--sra-cache`'s help text now cover a plain project's own per-accession lock, not only
  the store's.
- Five more commands stop between units on a signal, record what finished, and exit with status 130:
  `extract_target_reads` (in both the mapping and the assembly loops), `store_adopt` (including
  `--dry-run`, which used to exit 0), `sra_profile`, `download_metadata` and `parse_metadata`.
- New download result strings reach `failed_accessions.txt` and the registry: `interrupted` (the run was
  stopped before or during this accession), `locked: ...` (another process held the accession's lock
  for longer than `--lock-wait`) and `lock lost: ...` (another process took the lock over while this one
  was stalled, so nothing was published). `locked:` is classified `unknown`, so the retry pass tries the
  accession again and waits up to `--lock-wait` for the lock once more.
- New files on disk: `fastq/.locks/` (a plain project's per-accession locks), `<store>/locks/<ACC>.used`
  (its modification time records when a project was last handed the dataset, read by `store_gc`),
  `<lock>.reclaim` guards held for a few file operations during a stale-lock takeover or a release, and,
  rarely, a `<lock>.reclaimed.<hex>` file left when a process is killed in the middle of a takeover.
- `store_gc --json` dataset entries gained a `downloaded` field (the catalogue's download time).
- A write failure at a former `DataFrame.to_csv` site now raises `DataAccessError` with a "Failed to write
  CSV file" message instead of the underlying `OSError`.
- `Bio.Entrez.email` and `Bio.Entrez.api_key` are now set only around each NCBI efetch call and restored
  afterwards, so a library host that sets its own values keeps them.
- A second `extract_target_reads` that needs the minimap2 index another process is building waits, with
  no time limit, for that build to finish and then reuses the index.

### Fixed

- A concurrent write from `download_sra` could be lost when `select`, `blacklist`, `download_metadata`,
  `parse_metadata`, `sra_profile`, `sra_report`, `sra_validate`, or `status --init`/`--reconcile` ran at
  the same time and wrote back a registry snapshot loaded before the download's outcome was recorded.
- Two downloads of the same accession into one plain project could stage into the same
  `<ACCESSION>_temp` folder and delete each other's in-flight output.
- `store_gc --yes` and `store_verify --fix-state` could remove or rewrite a dataset another process had
  since started using, because each checked the lock only once, before doing its own work, rather than
  again immediately before writing.
- A second `Ctrl-C` (or `SIGTERM`/`SIGHUP`) during a command's final registry flush could abandon the
  write outright; it is now logged and the flush is allowed to finish.
- The store catalogue could be written while a process held the project registry lock, so a slow or
  contended catalogue write held up every other command waiting on the registry; `sra_profile`,
  `sra_report`, and `sra_validate` now record store usage after the registry batch that recorded their
  analysis results has already released the lock.
- `terminate_children()`, called after an interrupted download, left the whole process's stop flag set
  for the rest of that process's life, so an in-process caller that ran a second download after an
  interrupted one found its tools refused to start. Each download run now carries its own stop token;
  interrupting one run's tools no longer affects a concurrent run in the same process.
- `pigz`, started to compress a finished download, was not stopped by an interrupt and could keep
  running after its parent command exited; it now receives the same per-run stop token as the download
  tools and is terminated with them, and a `pigz` stopped this way is reported as "compression skipped
  (interrupted)" rather than "compression failed". Without `pigz`, the Python gzip fallback now checks
  the stop token between 1 MiB blocks, so an interrupted download no longer waits for a whole file to be
  compressed.
- A download into a shared store, or a `store_adopt`, whose holder stalled past the lock's stale window
  (a sleeping laptop, a suspended job) while another project took the lock over could still publish its
  staged folder into the store. Both now confirm the lock is still theirs immediately before the
  publishing rename and report `lock lost: ...` with nothing published otherwise. A registry write and a
  store catalogue commit make the same check and raise `LockLost`, leaving a registry batch's queue
  intact for the next flush.
- A per-sample read extraction interrupted while `samtools fastq` was writing left partial FASTQ files
  under their final names, which `status --reconcile` and `status --init` then recorded as a finished
  extraction, so the next `extract_target_reads` skipped the sample. The export now writes into a
  dot-prefixed folder inside the sample's output folder, and each file is renamed into place only once
  the sample is complete. The names published are recorded in a hidden `.<GENOME_ID>.extracted.json`
  beside them; a rerun that writes fewer files removes only the files that genome's previous run
  recorded, never a file another genome's record lists (a genome called `G1_1` next to a paired `G1`),
  and nothing beyond the files it overwrites in a folder written by an earlier version.
- A run that waited for another run's download of the same accession and then found the files counted
  a download attempt and replaced the other run's result message with "already exists"; it now records
  no attempt and leaves an existing record as it was.
- A "locked: ..." or "lock lost: ..." result whose holder pid, time or lock path contained "403" or "404"
  was classified not-found and skipped by the retry pass; lock messages are now classified before the
  error codes, and 403 and 404 count only as whole numbers.
- `store_gc --yes` removing several datasets could delete the catalogue row of one that a download
  published again after gc had moved the old copy aside; such a dataset now keeps its new copy and row
  and is reported `in_use` with "published again during removal".
- A signal arriving while a command was installing or restoring its signal handlers escaped as a
  traceback and could leave a handler installed; the command now logs "Interrupted" and exits 130, and
  every earlier handler is recorded before the new one goes in, so it is always restored.
- An atomic write onto a read-only file (mode 0o444) failed with `PermissionError`; the target's mode is
  now copied onto the temporary file only just before the rename, so such a file is replaced and stays
  read-only.
- A minimap2 that exited 0 without writing its index published a zero-length index that every later run
  reused; an empty index is now refused when built and rebuilt when found.
- Two first-time `store_init` runs on a new store could each write a marker with its own id, the last one
  replacing the first; the marker is now published with a hard link that fails when it exists, so the
  second run keeps the first run's id (a filesystem without hard links falls back to the old rename).
- The taxonomy validation caches are read and appended as UTF-8, matching the UTF-8 written when a cache
  is reset.
- `store_adopt` no longer copies leftover temporary files (`.<name>.<host>.<pid>.<token>.tmp`, or the
  `<name>.tmp.<pid>` of earlier versions) from a project folder into the store, nor counts them in its
  free-space estimate.
- The warning for a catalogue write that failed after a store publish names the store root and the
  `store_reindex --data-root` command to repair it.
- The copy made by `store_link --mode copy` (and by `download_sra` with a copy link mode) is built under a
  dot-prefixed name and renamed into place, so an interrupted copy no longer leaves a partial dataset
  folder in the project.

### Testing

- `make test`: 2707 passed, 4 deselected, 95% coverage (was 2449 passed at 0.5.1); the same 2707 pass
  with no bioinformatics tool on `PATH`.
- `make check` and `make pipeline` pass on the final state of the branch.

### Upgrade notes

- Lock and log wording changed. A waiter now logs "waiting for <what>: held by <holder>" while it waits;
  a dead-holder takeover logs "Took over the lock on <what>: its holder is no longer running (pid N on
  host H)"; a repeated signal logs "Received SIGTERM again; N more will abandon the write" (or `SIGINT`
  or `SIGHUP` in place of `SIGTERM`). A script that greps metaquest's logs for the old "Gave up waiting
  for ..." wording needs to match the new "<what> is locked by pid N on host H since T: <lock>" message
  instead.
- `fastq/<ACCESSION>_temp` leftovers from a version before this one are still recognized, counted as
  transient bytes, and cleaned up on the next download; new downloads stage under
  `fastq/.metaquest-tmp/` instead.
- A registry lock orphaned by a holder killed on a different host is now judged stale only after 120
  seconds (up from 30), so a writer can be blocked behind an orphaned lock for up to 120 seconds; each
  waiter still gives up and reports a failure after its own 30 second wait, so a busy project can see a
  burst of "locked" results that a retry resolves once the 120 seconds pass. A holder that died on the
  same host is still taken over at once.
- A `--sra-cache` folder shared by two projects without a store is still not locked; avoiding concurrent
  prefetch runs into a shared cache directory remains the caller's responsibility, as before this release.

## [0.5.1] - 2026-09-26

### Performance

Measurements below are from a 14,146-dataset project (the crispatus Lactobacillus study) on an
external USB volume on an Apple silicon laptop, unless marked synthetic.

- `download_sra` now records its outcome in one registry transaction per run instead of one write
  per accession: 549 accessions recorded in 184.7 s before, 0.8 s after. The registry file is also
  written as compact JSON instead of indented JSON, with the same keys and values, so the file
  shrank from 16.0 MB to 9.7 MB; a tool that reads the registry as JSON is unaffected.
- `sra_profile --sample-size 10000` on one 11.3-million-spot paired run: 25.1 s and 279 MB peak RSS
  before, 20.9 s and 219 MB after. The sampler now indexes records and reads the gzip stream in 1
  MiB blocks, and keeps quality scores as histograms rather than per-read lists; the remaining time
  is reading the compressed files from disk.
- `extract_target_reads` now runs minimap2 with `--sam-hit-only`, so only mapped records reach the
  SAM file: on SRR23946447 the SAM shrank from 4.09 GB to 742 MB and wall time fell from 79 s to 30
  s, with the same `mapped_reads` count as before.
- `parse_containment`: 2.29 s and 283 MB before, 0.92 s and 210 MB after on the same project; a
  synthetic 20,000-by-5 table went from 4.34 s to 0.39 s. The summary and screening tables are now
  built with numpy in one pass, and a screening run records one timestamp per run instead of one
  per row; output tables are byte-identical to before.
- Download and store adoption now compute md5 and read counts in one streaming pass
  (`fastq_digest`) instead of reading each mate file several times; the store path now reads mate 1
  twice instead of four times.
- `parse_metadata`: each XML file is now parsed once into a sparse attribute table; 2,201 files went
  from 0.74 s to 0.41 s with a warm cache, and RSS from 185 MB to 136 MB. The TSV output is
  byte-identical to before.
- `status` and `results_table` no longer make per-accession filesystem probes or quadratic passes
  over the registry: `status --json` went from 1.27 s to 0.63 s, and `results_table` (29,190 rows)
  from 1.15 s to 0.75 s. Both outputs are byte-identical to before.
- The NCBI and GTDB taxonomy clients now reuse one HTTP session with retries, and write their caches
  incrementally instead of once at the end of a run.

### Changed

- `find_by_taxonomy`'s containment/taxonomy annotation is now built with a `melt` instead of a
  per-row loop, and drops rows whose containment is zero (or missing); the previous behaviour kept
  those as zero-valued rows. Rows keep the previous order (sample by sample, genomes in column
  order), so a table without zero or missing cells is the same as before. `--output-format summary`
  now lists only the samples with a positive containment for the requested taxon; a sample whose
  containment is zero for every genome of that taxon was previously listed with 0.0.

### Removed

- `metaquest.data.metadata.get_unique_sample_attributes`; `parse_metadata` collects the attribute
  names in its single pass and no command called the function.

### Fixed

- A registry batch mutation that raises partway through is now rolled back in full: the registry
  is reloaded from the file and the mutations that succeeded before it are replayed, so no
  accession is left with some of its fields updated.
- During `download_sra`, once a SIGTERM or SIGHUP has been turned into an interrupt, both signals
  are ignored until the queued outcomes have been written, so a repeated `kill` no longer loses
  them.
- `scripts/check_ascii.sh` checks every matching file on disk when its root is not a git
  checkout (a `git archive` export, an unpacked sdist) and fails when it finds no source file at
  all, instead of passing after checking nothing.
- `parse_containment` sorts the parsed table by `max_containment` with a stable sort, so accessions
  with the same maximum keep their input order. The default sort ordered ties differently on
  different CPUs, so the same matches gave a different table, screening order and pinned test
  output on Linux and macOS.
- The CI lint step no longer passes `--extend-ignore` on the command line, which replaced the
  setup.cfg list and re-enabled the pydocstyle style codes the project switches off.

### Testing

- Added `tests/test_performance_regressions.py` with timing and memory bounds covering the registry,
  the FASTQ sampler, containment processing, metadata parsing, and `status`/`results_table`.

## [0.5.0] - 2026-09-25

### Breaking changes

- Optional packages moved behind extras: `analysis` (scikit-learn, scipy), `interactive` (plotly,
  jinja2), `maps` (cartopy), `sourmash`, and `all` (every extra). A plain `pip install .` no longer
  installs any of them; a command that needs a missing extra stops with a `ConfigurationError`
  naming the extra and the interpreter to install it into, for example
  `/path/to/python -m pip install 'metaquest[analysis]'`. `statsmodels`, `umap-learn`, `networkx`,
  `upsetplot`, and `seaborn` are no longer dependencies at all. A command that previously produced a
  degraded plain-text or plain-HTML output when a package was absent now stops with that same error
  instead; the correlation heatmap in `interactive_plot` is drawn with matplotlib, which is a core
  dependency, so it no longer needs an extra.
- `sra_stats`, `sra_profile_quality`, `sra_dashboard`, and `sra_compare` are replaced by two
  commands: `sra_profile` (statistics table and a quality profile JSON per accession) and
  `sra_report` (one HTML report, with a group comparison and statistical tests when `--groups-file`
  is given). The four old names still parse; each prints the new command name and exits with status
  2 rather than running. They are hidden from `metaquest --help`. See "Renamed SRA commands" in
  README.md for the full flag mapping.
- `sra_validate --accessions` is replaced by `--accessions-file` and repeatable `--accession`,
  matching `sra_profile` and `sra_report`. `--include-contamination` is removed; the adapter figures
  are always computed.
- GC content is reported in percent (`gc_percent`, 0-100) everywhere: the `sra_profile` CSV, the
  per-accession quality profile JSON (`complexity_score` also replaces the old ad hoc key), and the
  `results_table` output. The previous `gc_content` key (a 0-1 fraction in the JSON, already percent
  in the CSV) is gone. `avg_quality` in the CSV is now a per-read mean over a sample of reads from
  every mate file, rather than mate-1 records only.
- `results_table` gains `total_reads`, `gc_percent`, and `quality_grade` columns, inserted after
  `run_size` and before `mapped_reads`. Every later column shifts position accordingly; a script
  that reads `results_table`'s output by column index rather than by header name must be updated.
- `import metaquest.data.sra` now imports a package (`metaquest/data/sra/`) rather than a single
  module. The public functions re-exported from `metaquest.data.sra` are unchanged; private helper
  names moved to their new submodules (`fastq`, `cleanup`, `accession`, `store_handoff`, `retry`,
  `download`) and are no longer importable from the old single-file path.
- Public Python API changes outside `metaquest.data.sra`:
  `metaquest.data.sra_metadata.calculate_read_statistics` is removed (statistics come from
  `metaquest.sra.dataset_stats` and `metaquest.sra.profiles`), and
  `metaquest.data.sra_metadata.generate_statistics_report` now takes a sequence of statistics rows
  and an output path instead of a FASTQ folder and a sample size. `SequenceQualityAnalyzer` moved to
  `metaquest/sra/quality.py`; it is still importable from `metaquest.sra` and
  `metaquest.sra.analytics`. `metaquest.data.branchwater_search.sourmash_hint` and
  `SOURMASH_HINT` are removed; a missing sourmash is reported through
  `metaquest.core.optional.require` like every other extra.
- Error lines that end a command with exit status 1 (for example "No accessions found in file" and
  "FASTQ folder ... does not exist") are now written to the log on stderr instead of stdout. A
  script that parsed those lines from stdout must read stderr or check the exit status instead. The
  `--json` error document for a missing store stays on stdout.
- Registry analysis keys changed: `sra_profile` records its outcome under `profile` and `sra_report`
  under `report`, replacing the `sra_stats` and `quality` keys. Only `results_table` and `status`
  still fall back to the old keys when reading a registry written by an earlier version.

### Exception handling

- Added a `flake8-blind-except` lint gate (`B902`, bare `except Exception:`) with a per-module
  exemption list in `setup.cfg` that only shrinks as modules are narrowed.
- Narrowed broad exception handling in six modules (`metaquest/data/sra.py` and its later split,
  `sra_metadata.py`, `metadata.py`, `branchwater.py`, `visualization/reporting.py`, and
  `sra_intelligent.py`): a caught exception now names the specific error types the surrounding code
  can actually raise, so a programming error (for example a `TypeError` or `AttributeError`) is no
  longer swallowed and reported as a data or network failure. Corrupt gzip streams (`zlib.error`)
  and truncated NCBI replies are now handled explicitly rather than falling through a catch-all.

### Output

- Every command now writes to stdout through exactly one channel: `BaseCommand.emit`,
  `emit_raw`, and `emit_json` (plus the module-level `emit_error_json` for callers with no command
  instance). Library modules that used to print now log or return their lines instead, and the
  calling command writes them. The four `--json` commands (`status`, `store_status`, `store_usage`,
  `store_gc`) each emit exactly one JSON document.
- Added a `make check` gate (`scripts/check_no_print.sh`) that fails on any `print(` or
  `sys.stdout.write` outside `metaquest/cli/base.py`.

### Dependency audit

- Added `pip-audit --strict --desc .` to the CI lint job, and a new weekly `audit.yml` workflow that
  runs the same check on a schedule.
- The weekly Branchwater workflow's verify step now makes real assertions about the pipeline's
  output instead of only echoing that it ran.

### Documentation

- Added a `flake8-docstrings` gate (`D101`/`D102`/`D103`: a docstring on every public class, method,
  and function), enforced on a defined set of modules: the registry, the `metaquest/data/sra/` and
  `metaquest/store/` packages, the `metaquest/cli/commands/store/` and
  `metaquest/cli/commands/status/` packages, `processing/selection.py`, `processing/results.py`,
  `processing/status_report.py`, `data/registry_blocks.py`, `core/optional.py`, `sra/quality.py`,
  and the SRA analysis command and profiling modules. Every other
  module has the check switched off explicitly in `setup.cfg` rather than by a wildcard, so the
  enforced set stays visible and can grow one file at a time.

### Module layout

- Split `metaquest/data/sra.py` into the package `metaquest/data/sra/` (`fastq`, `cleanup`,
  `accession`, `store_handoff`, `retry`, `download`), each module under 600 lines.
- Split `metaquest/cli/commands/store.py` into `metaquest/cli/commands/store/` (one module per
  command: `init`, `status`, `reindex`, `adopt`, `verify`, `link`, `usage`, `gc`, plus a shared
  `_shared.py`), and `metaquest/cli/commands/status.py` into
  `metaquest/processing/status_report.py` (report building) and `metaquest/cli/commands/status/`
  (`command.py`, `suggest.py`, `render_text.py`).
- Added typed registry blocks (`metaquest/data/registry_blocks.py`): every block in
  `metaquest_registry.json` (screening, selection, exclusion, download, metadata, analysis,
  extraction, assembly, project, store) is now a dataclass with `from_dict`/`to_dict`, round-trips
  any key it does not itself declare, and keeps the on-disk JSON layout unchanged. `update_linked`
  moved from the store commands into `metaquest/data/registry.py`, next to the other registry
  writers.
- Added a module size and maintainability gate: 800 lines and a maintainability index of 20 or
  higher per module under `metaquest/`, with a shrinking `KNOWN_EXCEPTIONS` list in
  `tests/test_module_sizes.py` for modules that were already over a ceiling when the gate was
  added.
- Added an ASCII-only gate for `metaquest/`, `tests/`, `scripts/`, `Makefile`, `setup.cfg`, and
  `pyproject.toml`. A single line can be exempted with an `# ascii-ok` marker and a reason; the
  package author's name and the test fixtures that check Unicode handling use it.
- Added a documented-commands gate: every command visible in `metaquest --help` must be mentioned
  in README.md, and README.md must not invoke a command the CLI registry does not have.

### Continuous integration

- Added a nightly smoke workflow (`.github/workflows/smoke.yml`) that runs
  `scripts/smoke_chain.sh` against the real SRR2517620 run on NCBI/SRA, using conda-installed
  tools (fasterq-dump, prefetch, minimap2, samtools); `make test-network` runs the same chain
  locally.

## [0.4.0] - 2026-09-24

### Breaking changes

- Python 3.11 is no longer supported; MetaQuest now requires Python 3.12 or later.
- `download_sra` no longer prints its own download summary from the CLI layer; the summary line
  now comes from the store layer once, so a run that only links existing store datasets is no
  longer counted as a download.
- `select_datasets --no-record` now refuses to overwrite a recorded selection's output file; give
  it a different `--output` path or run without `--no-record`.

### Python floor and packaging

- Raised the minimum supported Python version to 3.12 and the floor for MetaQuest's own
  dependencies to versions tested on 3.12 and 3.13.
- Documented the Python 3.12 floor in the README and dropped stale test-count claims.

### Continuous integration

- The test job now runs on a matrix of Ubuntu and macOS with Python 3.12 and 3.13, with a mypy
  type-check step enabled.
- `make clean` no longer deletes `test_data`, so a local pipeline run's sample data survives a
  clean.
- Tests that exercise the SRA download path no longer depend on `sra-tools` being on `PATH`.

### Store

- The shared store now keeps an append-only journal of project and usage records; `store_reindex`
  replays it after rebuilding the catalogue from sidecar files, so a reindex no longer loses which
  projects used which datasets.
- `store_gc` refuses to run against a catalogue that has datasets but no project records, and
  reports its `--json` refusal the same way it reports other refusals.
- `store_adopt` refuses an empty accession folder, verifies the store's copy before removing
  anything from the project, and `--copy` no longer records a link when the accession is already a
  duplicate in the store.
- `store_verify --fix-state` only promotes a dataset to `complete` once its files have actually
  been read through, never on size or md5 alone; `--rescan` reports a folder with no FASTQ files as
  missing instead of silently skipping it; spot counts are also read from a project's own
  `metadata/` folder when the store does not have them.
- Folder scans, copies, and removals throughout the store and CLI now ignore macOS AppleDouble
  sidecar files (`._*`) and dotfiles, so they never look like real FASTQ, FASTA, CSV, XML, or
  profile files.
- Journal replay and catalogue backfill now keep the earliest recorded `first_used` and the latest
  `last_used` date for a dataset, and preserve the hostname and time already on a catalogue row
  instead of overwriting them.

### Downloads

- Ctrl-C during `download_sra` now cancels pending downloads and terminates any running `prefetch`
  or `fasterq-dump` process immediately, without waiting on a store lock or on a tool started after
  the interrupt.
- `fasterq-dump` scratch files are written under the store's own tmp folder and cleaned up
  afterwards.
- A relink of an already-downloaded accession keeps its recorded read count, and the download
  summary counts relinks separately from new downloads.
- The hint to install `sourmash` now names the Python interpreter that is actually running, not a
  generic `pip install` line.

### Selection

- `select_datasets` gained `--max-run-size`, `--min-spots`, and `--max-spots` (with a K/M/G/T
  suffix on run size) and `--platform` filters, and logs the excluded count and total volume of
  what a filtered selection kept. A requested filter whose column is missing from the metadata table
  (as in a Branchwater-derived table) is an error, and nothing is written or recorded.
- `select_datasets --no-record` skips writing the selection into the project registry, and
  `status --next` prints a runnable reselect command that reproduces the metadata filter, top-N,
  source table, and the new size/spot/platform filters, shell-quoted where needed.
- `select_datasets` now records the metadata table it actually used, whether autodetected or
  passed explicitly. The `status --next` reselect hint warns once on a malformed recorded filter
  value and leaves that flag out of the suggested command (or uses its default) instead of failing
  silently.

### Extraction and assembly

- `extract_target_reads` now records breadth of coverage (fraction of the reference genome covered
  at 1x or more) and length-weighted mean depth from `samtools coverage`, alongside the existing
  read counts.
- Each extraction and assembly run gets its own megahit scratch folder under the output directory,
  named with a random suffix, so concurrent runs sharing one output folder no longer collide; a
  dangling `fastq/<ACC>` symlink is reported as one warning naming every broken link.
- `--force` on read extraction or assembly now clears a stale registry record even when its output
  folder is already gone, and keeps an earlier recorded genome-download manifest instead of
  discarding it.
- Extraction and assembly dry runs summarise what would run; megahit errors surface the tool's own
  message, and tmp-dir and parameter changes are logged explicitly.
- A forced `extract_target_reads` rerun that keeps no reads now removes the earlier
  `<genome>_coverage.tsv`.

### Results

- Added the `results_table` command, which writes one row per screened (accession, genome) pair
  from the registry and the parsed containment table, joined with the extraction and coverage
  columns for that pair.

### SRA info and metadata

- `sra_info` now reports run accessions, sizes, base counts, and dates from the `RUN` element, and
  its accession filter also matches BioProject and BioSample identifiers, at the package level
  rather than the run level.
- Metadata parsing now fills in library strategy, source, and selection from the
  `EXPERIMENT`/`DESIGN` section when present.
- `sra_compare` writes JSON-safe results (no raw numpy types) and works from a single accession or
  a quality-profiles directory without requiring `--accessions-file`.
- A package whose inspection fails no longer makes `sra_info` fall back to listing every run
  returned; that fallback now applies only when every package was inspected without error, and a
  separate warning names how many packages could not be inspected.

### Documentation

- Reworded the README's testing summary from "All modules now thoroughly covered" to "Critical
  modules now thoroughly covered", matching the modules the README and CLAUDE.md elsewhere still
  list as uncovered.
- Reflowed a wrapped paragraph in the README's read-extraction section.
- Dropped stale per-module coverage percentages from CLAUDE.md's list of well-tested reference
  files, keeping the file list itself.
- Added a Releases section to the README describing the tag-triggered release workflow.

## [0.3.1] and earlier

See the git history for changes before 0.4.0; this file starts tracking releases from 0.4.0.
