# Code audit, SWOT and gap analysis, 2026-10-12

Scope: the whole `prepare` package at commit 16a9ac7 (v0.10.0 plus two commits), 38k lines
across about 140 modules, 3247 tests, 96% line coverage. Five reviewers each took one slice
(SRA download layer; store and registry; data and processing; CLI, settings, packaging and CI;
visualization and plugins). Every "confirmed" bug was either reproduced in a scratch folder or
checked against the source by two readers. Nothing in the repository was changed.

## 1. Confirmed bugs

Severity: H = data loss, security, or wrong scientific result; M = wrong behaviour without an
error, or an error on legal input; L = cosmetic or edge case.

### High

| # | Where | Defect | Failure scenario |
|---|-------|--------|------------------|
| H1 | `prepare/data/sra/download.py:323`, `data/sra/accession.py:404,424,548`, `data/sra/store_handoff.py:469` | Accessions read from the file are joined into paths and `rmtree`d before the format check in `run_secure` runs. | A line `../../victim` in the accessions file deletes the sibling folder `victim`; a line `..` deletes the whole FASTQ folder; an absolute path deletes that path. Then the run reports "Invalid SRA accession format". Reproduced. |
| H2 | `prepare/store/usage.py:31,184`, `data/registry.py:162` | The catalogue stores the registry path as given (`str(registry.path)`), unresolved. A relative `--registry` makes the project look stale from any other cwd, and `linked_by` skips stale projects, so the symlink does not protect the dataset. | `store_link --registry prepare_registry.json SRR2`, then after the 24 h grace `store_gc --include-stale --yes` from another folder removed `sra/SRR2` and left the project link dangling. Reproduced. |
| H3 | `prepare/processing/diversity.py:322,337` | PERMANOVA subtracts the mean squared distance from both sums of squares, which is not part of the formula; NaN group labels count in the total but in no group; the p-value lacks the +1 correction. | Two groups of 5 gave F = 56.9 where the reference formula gives 33.7. P-values from `diversity_analysis --permanova-formula` are therefore wrong. Reproduced. The test at `tests/test_diversity.py:95` only checks that keys exist. |
| H4 | `prepare/cli/commands/status/command.py:166-172` | `status --reconcile` resolves `ProjectPaths` against the cwd but finds the registry by walking up. | From `proj/sub`, every downloaded accession is marked "missing" and that is saved to `proj/prepare_registry.json`. Reproduced. Same cwd/registry mismatch in `store_link`, `store_unlink` and `ensure_project_identity` (project path recorded as cwd). |
| H5 | `prepare/visualization/interactive.py:53,67-75,238` | The metadata Series is indexed by sample name but assigned into a RangeIndex frame, so every value is NaN. | `interactive_plot --plot-type pca --color-by grp` draws one trace named "nan"; `--size-by` is dropped with a misleading warning. Colouring by metadata has never worked for PCA or t-SNE. Reproduced. |

### Medium

| # | Where | Defect | Failure scenario |
|---|-------|--------|------------------|
| M1 | `prepare/cli/commands/store/link.py:53,150,223`, `store/layout.py:118` | No accession validation in `store_link`, `store_unlink`, `store_verify`. | `store_unlink ../important_link` deleted a symlink outside `fastq/`; `store_link --accept-partial ../tmp` linked the store's tmp folder into the project and wrote `../tmp` into registry, catalogue and journal. Reproduced. |
| M2 | `prepare/data/sra/retry.py:249` | The retry pass always passes `force=False`. | A failed `--force` re-download is reported as "already exists" or "linked from store" on retry, exit 0, old data kept, registry overwritten to "downloaded". `tests/test_data_sra.py:1713` pins the wrong behaviour with an obsolete reason. |
| M3 | `prepare/data/sra/retry.py:160,266`, `cli/commands/sra.py:516,524` | Retry results carry a "Retry N:" prefix but callers compare unprefixed strings. | "Retry 1: linked from store" is counted as a fresh download: no store source, no `update_linked`, no catalogue usage row, a timing recorded; settled "incomplete:" results are refetched again with `--max-retries 2+`. |
| M4 | `prepare/data/sra/space.py:250` | The free-space guard keys reservations by accession, so a duplicate line leaks one reservation forever. | With a duplicated accession and tight disk, a later dataset waits in `reserve` until the run is killed. Reproduced. |
| M5 | `prepare/sra/quality.py:181` | `--sampler head` divides by `len(sequence)`. | A zero-length read (legal FASTQ, written by fasterq-dump for empty mates) aborts `sra_profile` and `sra_report` with ZeroDivisionError. Reproduced. |
| M6 | `prepare/data/branchwater.py:378` | Thresholds built as `k * step` give 0.30000000000000004 and the count uses `>=`. | Maxima of exactly 0.3, 0.6, 0.7 are not counted at their threshold; the summary table undercounts. Reproduced. |
| M7 | `prepare/data/branchwater.py` `process_branchwater_files`, `_process_genome_containments` | `validate_header()` return value is ignored. | A CSV with columns `foo,bar` is copied into the target as a Branchwater file. The test mocks the validator to raise, hiding the bug. Reproduced. |
| M8 | `prepare/data/taxonomy.py` `NCBITaxonomyClient.__init__` | Rate limit inverted: 0.1 s without an API key, 0.33 s with one. | Without a key, about 10 requests per second against NCBI's limit of 3, so 429s and retries. `tests/test_taxonomy_extended.py:159` asserts the inverted values. |
| M9 | `prepare/data/read_extraction.py` `_map_and_extract` | A forced re-extraction that maps nothing removes only the coverage table. | Old `targeted/<ACC>/<G>_[12].fastq.gz` and `.extracted.json` stay; registry says 0 reads, then `status --init`/`--reconcile` rediscovers the stale files. |
| M10 | `prepare/data/metadata.py` `_extract_metadata_fields`, `_download_single_metadata` | Fields come from the first `RUN` in the package; the single-accession fallback (after an HTTP 400 batch) writes the raw unfiltered response. | `SRR2_metadata.xml` whose package lists SRR1 first parses as SRR1 with SRR1's spot count; registry rows and truncation checks inherit it. Reproduced with synthetic XML. |
| M11 | `prepare/visualization/explorer.py:298` | Template has `autoescape=False` and renders `fam.name` raw. | A taxonomy family string containing `<img onerror=...>` becomes live HTML in the explorer header. Reproduced. |
| M12 | `prepare/visualization/explorer.py:120-121` | Reads `location` and `collection_date`; the metadata prepare writes uses `geo_loc_name_country_calc` and `collection_date_sam`. | Location and Date columns in the explorer are always empty. Reproduced. |
| M13 | `prepare/plugins/visualizers/map.py:83,104-108` | NCBI `lat_lon` "38.98 N 77.10 W" splits into four tokens and is dropped; a numeric value raises AttributeError. | Standard BioSample coordinates give an empty map. No CLI command reaches the plugin yet. Reproduced. |
| M14 | `prepare/plugins/sourmash_plugin.py:60-65` | `prepare_parse` never passes `errors=` and never updates the registry. | Unreadable match CSVs give a partial table and exit 0. |
| M15 | `prepare/plugins/sourmash_plugin.py:67,131,225` | Errors are logged but nothing configures logging under sourmash, and the package NullHandler silences the fallback. | `sourmash scripts prepare_plot --file-path /nonexistent` prints only the banner and exits 1. Reproduced. |
| M16 | `prepare/cli/main.py:252` | `BrokenPipeError` falls into the generic handler. | `find_by_taxonomy ... \| head -2` logs an ERROR and exits 1. Reproduced. |

### Low

| # | Where | Defect |
|---|-------|--------|
| L1 | `prepare/sra/quality.py:284` | Homopolymer run at the 3' end of a read is never counted (`ACGTAAAA` gives rate 0.0). Reproduced. |
| L2 | `prepare/data/taxonomy.py` `validate_taxonomic_assignments` | Cache read with plain `read_csv`; once an empty `tax_id` is cached, IDs come back as 562.0. Reproduced. |
| L3 | `prepare/cli/commands/advanced_analysis.py:387` | `taxonomic_summary --min-abundance` is parsed and never passed on. |
| L4 | `prepare/plugins/formats/branchwater.py:101` | Opened without `utf-8-sig`; a BOM makes `acc` "missing" and the genome gets an empty summary. Reproduced. |
| L5 | `prepare/visualization/plots.py:277` | `--colors red` passed as a string; pie iterates per character ("Invalid RGBA argument: 'e'"). Reproduced. |
| L6 | `prepare/visualization/plots.py:348,372`, `cli/commands/metadata.py:375` | `plot_metadata_counts` saves to cwd not beside the input; without `--save-format` writes nothing and exits 0; radar fallback drops `limit`. |
| L7 | `prepare/cli/commands/explore.py:299` | `find_by_taxonomy` TSV keeps the pandas index as a junk first column. Reproduced. |
| L8 | `prepare/cli/base.py:138` | Global `-v` is dropped for `store_status`, which defines `--verbose`. Reproduced. |
| L9 | `prepare/cli/commands/metadata.py:114` | `download_metadata --dry-run` creates `metadata/`. Reproduced. |
| L10 | `Makefile:58` | Hint names `plotly_cdn_script()`, which does not exist (`plotly_js_script`). |
| L11 | `pyproject.toml:95` | `py.typed` is declared in package data but the file does not exist, so the wheel carries no type marker. |

## 2. Likely bugs (read, not reproduced)

- Store: `link_dataset` takes no dataset lock, so a link landing between gc's `.used` check and
  its `os.replace` is dangling while the registry says complete; copy-mode `copytree` is also
  unlocked (`store/link.py:201`, `store/gc.py:513`).
- Store: a registry with a newer `version` is silently downgraded to 2 and unknown top-level keys
  are dropped; a non-int version gives a traceback instead of `DataAccessError`
  (`data/registry.py:188,216`).
- Store: a catalogue row whose folder is missing is a gc candidate forever; an interrupt between
  rename-aside and commit leaves `tmp/<ACC>_gc` with a live row (`gc.py:511-519`).
- Store: reconcile `_mark_missing` does not clear `store.linked` (`registry_reconcile.py:240`).
- SRA: accessions-file comment lines are not skipped and duplicates are not removed
  (`download.py:324`), unlike `cli/base.py` `read_accessions_file`; a stop token set mid-round
  makes the serial retry loop record every remaining accession as failed.
- Data: an HTTP 200 NCBI error body becomes a permanent metadata file; `gtdb.py
  get_accessions_for_genus` caps at 100 with no paging; the GTDB search fallback takes `rows[0]`
  unchecked; `sra_metadata._fetch_batch_metadata` uses `retmax=len(accessions)` on an all-fields
  search; a zero-row Branchwater response is cached and reused; `genome_download
  extract_and_organize` never compares requested accessions with the zip and uses fixed temp
  paths that collide between concurrent runs; `_separator_for` reads `.tsv.gz` as comma.
- CLI: `PREPARE_ASSEMBLY_MEMORY` accepts `0G` at startup and fails at assembly; `--log-level`
  is case-sensitive only on the command line; `--log-level ERROR -v` silently becomes DEBUG;
  `doctor` drops `--log-file`, `-q`, `-v` after a config error; startup errors under `--json`
  give an empty stdout while later errors give `{"error": ...}`; `store = "x"` in the config
  gives an AttributeError (`store/resolve.py:92`).
- Viz: no `plt.close` anywhere, so library loops leak figures; plugins register as an import
  side effect and re-registration raises; `exportCSV` does not guard formula-leading cells;
  the report table comparator mixes numeric and text compare and is non-transitive
  (`utils/html.py:310,342`).

## 3. SWOT

### Strengths
- Storage discipline. Every output goes through an atomic write with fsync; downloads are staged
  and published with one rename; the store marker is published with `os.link`; the registry is
  load-mutate-save under a lock with `verify_held`; gc is staged with per-accession locks,
  rename-aside, batched catalogue delete and restore on failure.
- Lock file implementation is rigorous (O_EXCL, token compare under reclaim guard, heartbeat
  through the open descriptor, signal mask during creation).
- Subprocess safety. No shell strings anywhere; allowlisted flags; own process group per tool;
  timeouts and interrupts kill the group.
- CLI contract. Exit codes are consistent (2 usage, 3 configuration), `--json` output stayed
  clean under `-v` and under config errors, settings keep their source and mask secrets, logging
  removes only its own handlers.
- Test discipline. 3247 tests, 96% coverage, conftest strips every `PREPARE_*` variable and the
  home folder; gates for print calls, ASCII, atomic writes, module size and README command list.
- Offline-safe reports. plotly.js is embedded; two of three HTML templates autoescape.
- Clear error classification for downloads (pure function, fixed order) and clean separation
  of transient from permanent NCBI errors.

### Weaknesses
- Input validation is applied late or not at all at three entry points (download accessions
  file, store link/unlink/verify, explorer family names). H1, M1, M11.
- Path handling assumes cwd is the project root in several commands (H2, H4), while the
  registry lookup walks upward. The two models disagree.
- Scientific numerics are under-tested: PERMANOVA (H3) and the containment threshold grid
  (M6) have tests that check shape, not values. No statistic is compared against a reference
  implementation.
- Several tests pin bugs (inverted rate limit, mocked header validator, `force=False` on
  retry).
- String-typed results. Download outcomes are prefixed strings parsed by `startswith` in at
  least five places; the "Retry N:" prefix already broke three of them (M3).
- Duplication: three NCBI rate limiters, two Branchwater header validators, two GTDB cache
  writers, two accession-file readers, three folder-size walkers, PCA/t-SNE preparation twice,
  four `write_html` blocks.
- Ten modules over 600 lines (`registry.py` 1001, `read_extraction.py` 960, `metadata.py`
  797, `analytics.py` 791, `sra.py` 775, `adopt.py` 718, `registry_blocks.py` 679,
  `lockfile.py` 672, `accession.py` 650, `store_handoff.py` 637); `test_data_sra.py` is
  4358 lines.
- Plugin system is half-built: no abstract methods on the base class, no entry-point
  discovery (`PLUGIN_REGISTRY_NAME` is unused), registration by import side effect, map
  plugin unreachable from the CLI while the `maps` extra installs cartopy.
- CI runs fewer gates than `make check` (no ASCII, no docs-command, no radon, black only on
  `prepare/`), the release workflow does not check tag against version and runs no tests.

### Opportunities
- The PyPI name `prepare` is free (404 on 2026-10-12). Publishing claims it before someone
  else does; the release workflow is already wired for trusted publishing behind
  `vars.PYPI_PUBLISH`.
- A bioconda recipe (template exists in `docs/packaging.md`) would reach the target audience,
  since every external tool used (sra-tools, minimap2, megahit, seqkit, samtools, sourmash) is
  on bioconda.
- A typed result object for download outcomes would remove the M2/M3 class entirely.
- A single `validate_accession()` at the CLI boundary closes H1 and M1 in one change.
- Branchwater and sourmash are both in motion; the plugin layer, once finished, is the natural
  place for new containment sources.
- The HTML reports are close to a good product; fixing the explorer columns and escaping, and
  capping the results table, would make them usable on large runs.

### Threats
- A user-supplied accessions file is the first input of the pipeline; H1 means a malformed or
  hostile line removes data outside the project. Shared stores on HPC amplify H2.
- Published PERMANOVA p-values are wrong; anyone who has used the feature has unreliable
  results.
- NCBI rate limiting (M8) risks the shared institutional IP being throttled.
- The name `prepare` is generic; squatting on PyPI and discoverability are both risks.
- Dependence on sra-tools and NCBI/ENA endpoints: outages (Branchwater zero-result outage on
  2026-09-04) are cached and reused (likely bug).
- Bus factor: a single author, 2 stars, no open issues, no external contributors.

## 4. Gap analysis

| Area | Current state | Target | Gap |
|------|---------------|--------|-----|
| Input validation | Accession format checked inside `run_secure` only | Validate at every CLI and file boundary with `SRA_ACCESSION_PATTERN` | H1, M1; one shared accession-file reader (skip comments, dedupe, validate) |
| Project root | Mixed cwd/registry-walk | One `ProjectPaths` derived from the registry location | H2, H4 plus `store_link`, `store_unlink`, `ensure_project_identity` |
| Statistics correctness | Shape-only tests | Reference comparisons (scikit-bio for PERMANOVA, exact decimals for thresholds) | H3, M6, and the pinned-bug tests |
| Download result model | Prefixed strings | Small dataclass with kind, source, message | M2, M3, five `startswith` sites |
| Retry semantics | `force=False` on retry | Carry the caller's `force`; drop the obsolete test | M2 |
| Store concurrency | gc locks, link does not | Link and copy take the dataset lock | Likely race |
| Registry forward compat | Downgrade on save | Refuse to write a newer-version registry | Likely bug |
| gc scaling | O(datasets x projects) registry loads | Load each project registry once per run | Performance on large stores |
| Journal | Unbounded | Prune on gc, compact on reindex | Growth |
| HTML safety | 2 of 3 templates autoescape | All templates autoescape; CSV export guards formulas | M11 |
| Explorer | Wrong metadata columns, no row cap, Python loop | Use the written column names, cap rows, vectorise | M12, scaling |
| Interactive plots | Metadata join broken, no `plt.close` | Index-aligned join, close figures, write atomically | H5 |
| Plugins | Registry without discovery or abstract API | Abstract base, entry-point loading, map command or drop `maps` extra | Half-built feature |
| sourmash plugin | Silent failures, partial tables exit 0 | Configure a stderr handler, pass `errors=`, update registry | M14, M15 |
| NCBI clients | Three rate limiters, inverted delay, unchecked bodies | One client with correct limits and body checks | M8, M10, likely bugs |
| CI parity | `make check` stronger than CI | CI runs `make check`; release checks tag == version and runs tests | Packaging gaps, `py.typed` |
| Gate scripts | README check matches words in prose; print gate misses `stdout.buffer` | Match command table rows; broaden the print pattern | False negatives |
| Module size | 10 modules over 600 lines | Split by responsibility (registry, read_extraction, sra command) | Maintainability |
| Rename leftovers | `mq-` CSS prefix (120 hits), historical plans, Makefile hint, unused constant | Rename CSS prefix or leave documented | Cosmetic |
| Distribution | Not on PyPI or bioconda | Publish 0.10.x, add bioconda recipe | Reach |
| Known 0.9.1 items | Second interrupt skips collection; scratch removal after lost lock | Fix | Open |

## 5. Suggested order of work

1. Security and data loss: H1, M1 (shared validated accession reader), H2 (resolve registry
   path, check stale projects' links), H4 (project root from registry). Add a test for each.
2. Wrong results: H3 (PERMANOVA against scikit-bio), M6 (decimal thresholds), M8 (rate limit),
   M10 (filter the single-run fallback), M9 (remove stale targeted reads), M7 (honour
   `validate_header`). Replace the three tests that pin bugs.
3. Download result model: typed outcome, carry `force` through retry (M2, M3, M4, M5, L1).
4. Reports: H5, M11, M12, M13, M14, M15, L5, L6, L7, and `plt.close`.
5. CLI polish and CI parity: M16, L8, L9, L10, L11; CI runs `make check`; release checks
   version; `py.typed` added.
6. Structural: result dataclass, single NCBI client, plugin base, module splits, `mq-` prefix.
