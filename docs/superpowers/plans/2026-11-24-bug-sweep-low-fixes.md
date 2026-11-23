# Bug sweep 2026-11-23: Low findings and deferred minors (0.14.0) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the remaining findings of the 2026-11-23 bug sweep (the Low list at the end of
`2026-11-23-bug-sweep-medium-fixes.md`) and the minors deferred by the final reviews of 0.12.3 and 0.13.0.

**Architecture:** Ten independent tasks grouped by layer, then a release task. Each finding gets a test that fails
on `main` (448c337, v0.13.0) before its fix.

**Tech Stack:** Python 3.12, pytest, numpy/pandas/scipy (analysis extra); conda environment `metaquest`.

**Spec:** The sweep findings and the two final reviews (no separate document). Each task names the defect.

**Why 0.14.0:** Task 8 records `branchwater_search --output` screening under the output's stem and rejects
containment thresholds outside 0-1, and Task 5 refuses to replace a symlink that does not point into the store.

## Global Constraints

- Branch `fix/bug-sweep-low` off `main` (448c337). Stage files by name; never `git add -A`.
- Run as `cd /Users/andreassjodin/Code/prepare && conda run -n metaquest <cmd>`; the cwd resets.
- Test first, watched failing, then the smallest fix (CLAUDE.md).
- Module ceiling 800 lines, 700 for `SPLIT_MODULES` in `tests/test_module_sizes.py`; no new radon rank C
  (`MAX_RANK_C` must not grow), no rank D.
- Errors through `self.fail`; exit codes 0, 1, 2, 3, 4, 130 only. Atomic writes, ASCII only, no `print(`.
- A test needing an extra (scipy for Task 2's statistics) skips when it is absent.
- Store tests use a `tmp_path` store with `PREPARE_DATA`/`HOME` monkeypatched; no network.
- Modest scientific language. Update README/docs/CHANGELOG with each behaviour change, in the same task.
- `make check` and `make test` before each commit.

## Review Focus

1. Task 6: a `catalog_write` session that raises after `upsert_project` must leave no journal line, and one that
   commits must leave exactly one (test `test_a_rolled_back_session_writes_no_journal_line`).
2. Task 5: `store_gc --yes` with `tmp/<ACC>_old` present and `sra/<ACC>` missing must keep the staged copy
   (test `test_gc_keeps_staged_copies_when_the_published_folder_is_missing`).
3. Task 2: a FASTQ whose quality characters are all `@`..`h` (Phred+64) gives qualities 31 lower than read as +33;
   a normal +33 file is unchanged (tests `test_phred64_is_detected`, `test_phred33_is_unchanged`).
4. Task 8: `branchwater_search` without `--output` still records under the genome name
   (test `test_default_output_records_under_the_genome_name`).
5. Task 7: a NaN reaching `record_run` is written as `null`, and the file parses with `json.loads(..., parse_constant=...)`
   raising on NaN (test `test_nan_in_the_run_log_is_written_as_null`).

---

### Task 1: Statistics library: enrichment, distances, diversity

- `calculate_enrichment` (`prepare/processing/statistics.py:~63-104`): `fold_enrichment` becomes the ratio of
  proportions, `(observed / total_observed) / (expected / total_expected)`, and the table is sorted by it.
  Test: observed {A: 5, B: 5}, expected {A: 50, B: 950} -> fold(A) == 10.0, fold(B) == pytest.approx(10/19), A first.
- `_calculate_jaccard_distance`/`_calculate_dice_distance`/`_calculate_cosine_distance` (`statistics.py:~126-175`):
  two all-zero samples have distance 0.0 (as `diversity._calculate_jaccard_distance`). Test per metric.
- `diversity._calculate_shannon_diversity` (`diversity.py:~240-246`): compute over the non-zero proportions only
  (no pseudo-count); `pielou_evenness` returns 0.0 (no RuntimeWarning) for a sample with fewer than 2 species.
  Tests: row [10, 0, 0] -> 0.0 exactly; [5, 5] -> ln 2; pielou of [10, 0] -> 0.0 with `warnings` turned into errors.
- Commit `fix: enrichment uses proportions; empty samples are identical; Shannon without a pseudo-count`.

### Task 2: SRA quality figures

- GC and N content in `store/stats.compute_dataset_stats` (`store/stats.py:~254-264`) and `sra/quality._summarise`
  are base-weighted (GC bases / all bases), not a mean of per-read fractions. The statistics record keeps its field
  names. Test: one 1000 bp all-GC read and ten 10 bp all-A reads -> gc_content == pytest.approx(1000/1100).
  Bump nothing in the cache key unless the record is cached across versions; if it is (`store/stats.py` cache),
  add the bump in this task so old records are recomputed, and say so in the CHANGELOG.
- Phred offset (`sra/quality.py`): detect +64 from the sample: when the lowest quality character seen is `@` (64)
  or above and the highest is above `J` (74), decode with 64, else 33. Remove the unused `quality_encodings`
  attribute or use it. Tests `test_phred64_is_detected`, `test_phred33_is_unchanged`.
- `analytics._perform_statistical_tests` (`~601-615`): a comparison whose groups have no within-group variance, or
  fewer than 2 values in a group, is reported as skipped (`"statistic": None, "p_value": None,
  "significant": False, "note": "no within-group variance"`), not t=inf/p=0; remove the dead `except RuntimeWarning`;
  a NaN p-value from scipy is treated as skipped; `holm_adjusted` never sees NaN. Test with groups [100, 100] vs
  [200, 200] (skip if scipy is absent).
- `format_statistics_summary` (`data/sra_metadata.py`): label the read-length and GC lines "mean of datasets".
- `sra_profile --sample-size` help (`cli/commands/sra_profile.py:~54`): GC is sampled from every mate file.
- Commit `fix: base-weighted GC, Phred+64 detection, skipped tests without variance`.

### Task 3: `validate_taxonomy` with no species

- `cli/commands/advanced_analysis.py:~339-341`: an empty species list logs "No species to validate" and returns 0
  (no ZeroDivisionError). Test with an empty `--species-file`.
- Commit `fix: validate_taxonomy accepts an empty species list`.

### Task 4: Download and metadata errors

- `data/genome_download.py:~109-112, 136-138`: a `CalledProcessError` from `datasets` is reported with its stderr
  (`e.stderr` stripped, last 500 characters), raised `from e`; a stderr naming a network failure (reuse
  `classify_download_error(...) == "network"`) raises `NetworkError` (exit 4). Tests for both.
- `data/metadata_fetch.py:~147-159, 230-236`: HTTP 400 is reported as `"HTTP 400: bad request (check the accession
  and --api-key)"`, not "not found"; 404 stays "not found at NCBI". 4xx other than 429 is not retried. No sleep
  after the last attempt (single and batch, every exception branch). Tests: 401 -> one call, no sleep; 503 x3 ->
  two sleeps; 400 message.
- `data/sra/accession.py:~267-270, 319-339`: a compression failure of any file returns `(False, "compression failed
  for <names>")` and publishes nothing, as an interrupted compression already does. Update the existing test that
  expects success (`test_handle_download_output_reports_compression_failure_in_message`).
- `data/accessions.py`/`core/constants.py:9` vs `core/validation.py:16-41`: the accessions-file reader uses
  `validate_accession` (SRR/ERR/DRR runs), so `SRX123` is refused up front with the existing invalid-accession
  error. Test: a file with `SRX1` -> `ValidationError` naming it.
- Commit `fix: report datasets errors, retry only transient NCBI errors, never publish a half-compressed folder`.

### Task 5: Store safety

- `cli/commands/store/gc.py` `_leftover_candidates` (~330-350): `tmp/<ACC>_old` and bare `tmp/<ACC>` are leftovers
  only when `sra/<ACC>` exists; otherwise they are kept and reported once as "staged copy of <ACC> with no published
  folder; move it to sra/<ACC> or run store_verify" (warning). Test
  `test_gc_keeps_staged_copies_when_the_published_folder_is_missing`.
- `store/link.py:~105-121` (`_clear_existing`): an existing symlink that does not point into this store is not
  replaced; `link_dataset` raises `DataAccessError("<path> links outside the store (<target>); remove it first")`.
  A dangling link into the store is still replaced. Tests for both.
- `store_link`/`store_unlink` (`cli/commands/store/link.py`): positional accessions are de-duplicated in order.
- docs/shared_store.md: one sentence for the foreign-link refusal.
- Commit `fix: store_gc keeps orphaned staged copies; store_link refuses foreign symlinks`.

### Task 6: Store journal order and timestamps

- Journal lines are written after the SQLite commit (store rule 4 still holds): `Catalog` collects pending journal
  appends during a `catalog_write` session (`_pending_journal: Optional[List[Callable[[], None]]]`, None outside a
  session, which keeps today's immediate write for `journal_enabled` callers outside a session);
  `catalog_write` runs them after `conn.commit()`, `commit_now` runs those queued so far after its commit, and a
  session that raises drops them. Tests `test_a_rolled_back_session_writes_no_journal_line`,
  `test_a_committed_session_writes_each_line_once`, `test_commit_now_flushes_the_lines_queued_so_far`.
- Replay keeps project timestamps: `append_project` already stamps `at`; `upsert_project(..., at: Optional[str] = None)`
  uses it for `created` (first insert) and `last_seen`, and `replay` passes `record.get("at")`.
  Test: build, delete the catalogue, `store_reindex`; `projects.created` and `last_seen` equal the originals.
- Commit `fix: journal lines follow the catalogue commit; replay keeps project timestamps`.

### Task 7: Status, exports and the run log

- `status --export-tsv` (`processing/status_report.py:~463`, `cli/commands/status/command.py:~146-153`): an empty
  table is written with its header (the column lists the non-empty export uses). Test: no extractions ->
  `pd.read_csv(path, sep="\t")` has the extraction columns and 0 rows.
- `data/run_log.py:~134-145` `_json_safe`: NaN and +-inf become None; the registry writer (`registry_io.py:~124`)
  passes `allow_nan=False` after the same conversion. Test `test_nan_in_the_run_log_is_written_as_null`, and one for
  the registry.
- `status --next` (`cli/commands/status/suggest.py:~278-283`): an extraction with an assembly record (any contigs,
  including 0) is not suggested for `--assemble` again. Test.
- Commit `fix: header for empty exports, no NaN in JSON, no repeated --assemble suggestion`.

### Task 8: Command-line consistency

- `branchwater_search --output` (`cli/commands/branchwater_search.py:76-77, 109-116`): the registry genome id is the
  output file's stem (what `parse_containment` reads); without `--output` it is the genome name, as now. When the two
  differ, one INFO line says which is used. Tests: custom `--output foo.csv` records `foo`;
  `test_default_output_records_under_the_genome_name`.
- `genome_download`/`genome_prepare` (`cli/commands/genome.py:171-187, 333`): accessions from all sources are
  de-duplicated in order. Test.
- `extract_target_reads` (`cli/commands/read_extraction.py`): assembly-only flags given without `--assemble`
  (`--assembly-threads`, `--assembly-memory`, `--min-contig-len`, a non-default `--assembly-preset`,
  `--keep-intermediate`, `--no-coverage`) log one WARNING naming them. Test.
- Containment thresholds: new `unit_interval` argparse type in `cli/base.py` (finite, 0 to 1) for every flag that
  compares a containment (find them with `grep -rn '"--threshold"\|--min-containment' prepare/cli prepare/plugins`;
  leave flags that are not containments, such as quality thresholds, alone). Test parametrised over them with 1.5
  and -0.1 -> exit 2.
- `--max-cache-age-days` -> `non_negative_int`. Test.
- Commit `fix: consistent genome ids, de-duplicated accessions, warnings for ignored flags, bounded thresholds`.

### Task 9: Usage errors through `self.fail`

- Replace each `self.logger.error(...); return 1` usage-error path with `return self.fail(ValidationError(msg), self.name)`
  (same exit code 1): `explore.py:58,144,153,166,255,260,265`, `store/gc.py:~373` (`--dry-run` with `--yes`),
  `store/usage.py` (ambiguous project name), `sra_enhanced.py:145,157,452`, `sra_profile.py:~275`, `metadata.py:~282`.
  Line numbers are from 0.12.2; find them with `grep -n "logger.error" -A1 <file> | grep -B1 "return 1"`.
  Test: one parametrised test that each listed path still exits 1 and logs one line (use two or three
  representative commands; the rest is mechanical).
- Commit `refactor: report usage errors through BaseCommand.fail`.

### Task 10: Deferred minors from the 0.12.3 and 0.13.0 reviews

- `data/registry.py` `record_selection`: `set(accessions)` is built once, not per kept accession. Test: 10,000 kept
  and 10,000 accessions record in under 1 s (perf marker, scaled by `PREPARE_PERF_SCALE`).
- `processing/selection.py:~193`: the "dropped all N remaining candidate(s); none has a <col> value" warning fires
  when every row that reaches the run filter lacks a value, counting only rows not excluded/downloaded. Test.
- `download_sra`/`store_adopt` after an interrupt still log the unrecorded-usage rerun hint before exiting 130.
  Test for `download_sra`.
- `DownloadSraCommand._unrecorded_usage` is reset at the start of each `_run`. Test: two `execute` calls on one
  instance; the second reports nothing from the first.
- `extract_target_reads._resolve_threshold`: the warning counts only selected samples that are not excluded, and
  the table is read once. Test.
- `project_funnel._tally_extractions`: pairs and seconds count only extractions whose genome is still on disk for
  that accession (use `_stage_on_disk` per genome, or pass the per-genome sets `status_report` already computes).
  Test with two genomes, one deleted.
- Commit `fix: deferred minors from the 0.12.3 and 0.13.0 reviews`.

### Task 11: Gates, changelog and 0.14.0

- `make check`, `make test`; CHANGELOG `## [0.14.0] - <date>` with `### Changed` (Task 5 foreign links, Task 8
  genome id and thresholds) and `### Fixed`; bump `pyproject.toml`, `prepare/__init__.py`, `CLAUDE.md:7`.
- Commit `release: 0.14.0`. Push, PR and tag after the user's go-ahead.

## Not in this plan

- `prefetch --max-size` skipping a run with exit 0 ("prefetch wrote no archive", classed network): cannot be
  verified offline; left as it is.
- `store_verify` recording a run in read-only mode: its docstring says every run is recorded, by design.
