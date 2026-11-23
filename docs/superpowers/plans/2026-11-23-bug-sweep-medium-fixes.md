# Bug sweep 2026-11-23: Medium-severity fixes (0.13.0) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the Medium findings of the 2026-11-23 six-agent bug sweep and release them as 0.13.0.

**Architecture:** Thirteen independent fixes, grouped by layer: the download and metadata layer (Tasks 1-3),
the hand-offs between pipeline steps (Tasks 4-6), `status` and `project_report` (Tasks 7-8), the
shared store (Task 9), command-line validation (Tasks 10-11) and the SRA quality analysis (Tasks 12-13).
Each task starts with a test that fails on the branch it starts from.

**Tech Stack:** Python 3.12, pytest, pandas, numpy; conda environment `metaquest`.

**Spec:** The sweep findings (no separate spec document); the High fixes are
`docs/superpowers/plans/2026-11-23-bug-sweep-high-fixes.md` (0.12.3). Each task quotes the defect it fixes.

**Why 0.13.0, not 0.12.4:** Task 4 changes the default of `extract_target_reads --threshold`, Task 5 removes
registry entries that earlier versions kept, and Task 10 rejects flag values that were accepted before.

## Global Constraints

- Branch `fix/bug-sweep-medium` off `main` after 0.12.3 is merged. Stage files by name; never `git add -A`.
- Run commands as `cd /Users/andreassjodin/Code/prepare && conda run -n metaquest <cmd>` (the cwd resets).
- A bug fix starts with a test that reproduces the bug and fails, then the smallest fix (CLAUDE.md).
- Module ceiling is **700 lines** (`tests/test_module_sizes.py`; CLAUDE.md's 800 is out of date),
  maintainability index >= 20, no new radon rank C or worse (`MAX_RANK_C` only shrinks). `sra.py` (698)
  and `registry.py` (700) have no room: put new code in other modules.
- Errors through `self.fail`; only exit codes 0, 1, 2, 3, 4, 130. Usage errors from argparse types exit 2.
- Atomic writes, no `print(`, ASCII only, per-item logging at `item_level(...)`.
- Store tests use a `tmp_path` store with `PREPARE_DATA` and `HOME` monkeypatched; no network.
- Modest scientific language in docs and comments. When a flag or behaviour changes, update README.md,
  the matching `docs/` file and CHANGELOG.md in the same task.
- `make check` and `make test` pass before each commit.

## Review Focus

1. A project whose registry holds a selection made before 0.13.0 with no `threshold` in its criteria: Task 4
   must fall back to 0.1, not fail (test `test_threshold_without_a_recorded_selection_threshold_is_the_default`).
2. Re-screening with `parse_containment` on a matches folder that holds only some genomes' CSVs: Task 5 must
   drop vanished entries only for genome columns present in the new table
   (test `test_rescreening_leaves_genomes_absent_from_the_table_alone`).
3. `status --stage downloaded --genome G`: Task 7 must still remove accessions gone from disk when a genome
   is named (test `test_stage_filter_with_a_genome_leaves_out_accessions_gone_from_disk`).
4. `store_unlink` of an accession the project never recorded as used: Task 9 must not fail
   (test `test_unlink_without_a_usage_row_succeeds`).
5. `--older-than 0` on `store_gc`: still accepted (0 means any age) (test in Task 10).

---

### Task 1: Entrez fetches time out

Defect: `Entrez.efetch` in `prepare/data/metadata_fetch.py:137` and `:214` opens its connection without a
timeout (Biopython's `urlopen(request)`), so a stalled NCBI connection hangs `download_metadata` forever while
`_ENTREZ_LOCK` is held; the retry loop and Ctrl-C handling never run.

Behaviour: every `Entrez.efetch` call runs with a socket timeout of `ENTREZ_TIMEOUT_SECONDS = 120`, set with
`socket.setdefaulttimeout` inside `_entrez_credentials` (which already holds `_ENTREZ_LOCK`) and restored on
exit. A timeout raises `TimeoutError`/`socket.timeout` (an `OSError`), which the existing handlers retry as
`network: ...`.

**Files:**
- Modify: `prepare/data/metadata_fetch.py` (`_entrez_credentials`, new constant)
- Test: `tests/test_metadata_fetch.py` (find the module that tests `_download_single_metadata`:
  `grep -ln "_download_single_metadata" tests/`)

- [ ] **Step 1: Write the failing tests**
  - `test_efetch_runs_with_a_socket_timeout`: patch `Entrez.efetch` with a function that records
    `socket.getdefaulttimeout()` and returns a handle with one package; assert the recorded value is
    `ENTREZ_TIMEOUT_SECONDS` and that `socket.getdefaulttimeout()` is back to its previous value afterwards.
  - `test_a_timed_out_fetch_is_retried_as_network`: `Entrez.efetch` raises `TimeoutError("timed out")` on every
    call (patch `time.sleep`); assert the result is `(False, msg)` with `msg.startswith("network: ")`.
- [ ] **Step 2:** Run them; expect the first to FAIL (timeout is None).
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run `conda run -n metaquest pytest <that module> tests/test_cli_commands.py -k metadata -q`; expect PASS.
- [ ] **Step 5:** Commit `fix: time out stalled NCBI efetch connections`.

---

### Task 2: A number such as 1,404,213 no longer makes a failure "not found"

Defect: `_NOT_FOUND_ERROR_RE` (`prepare/data/sra/accession.py:40-42`) matches a bare `\b403\b|\b404\b`,
and `classify_download_error` tests not-found before disk-full. `"spots read : 1,404,213 ... storage exhausted"`
was classified `not-found`, so the retry pass skipped it and the disk-full abort never ran.

Behaviour: disk-full is tested before not-found. A 403/404 counts only in the forms tools print it:
`( 404 )` (prefetch), `HTTP 404`, `status 404`, `code 404`, `404 Not Found`, `403 Forbidden`
(regex `\(\s*40[34]\s*\)|\b(?:http|status|code)[ :]*40[34]\b|\b404 not found\b|\b403 forbidden\b`, case-insensitive),
plus the existing word patterns.

**Files:** Modify `prepare/data/sra/accession.py:40-42,103-125`; Test `tests/test_data_sra_download.py`
(the class that tests `classify_download_error`).

- [ ] **Step 1: Write the failing tests** (parametrised `test_classify_download_error_cases`):
  - `"Download failed: spots read : 1,404,213\nerr: storage exhausted while writing file"` -> `"disk-full"`
  - `"Download failed: wrote 404 files to /data/run403/"` -> `"unknown"`
  - `"failed to resolve accession SRR1 - no data ( 404 )"` -> `"not-found"` (stays)
  - `"HTTP 403 Forbidden"` -> `"not-found"` (stays)
- [ ] **Step 2:** Run; expect the first two to FAIL.
- [ ] **Step 3:** Implement; update the docstring's order note.
- [ ] **Step 4:** Run `conda run -n metaquest pytest tests/test_data_sra_download.py tests/test_cli_commands.py -q`; expect PASS.
- [ ] **Step 5:** Commit `fix: classify disk-full before not-found and match 403/404 only as status codes`.

---

### Task 3: Ctrl-C in the retry pass removes the build folder

Defect: the retry pass runs `download_accession` on the main thread, where Ctrl-C raises `KeyboardInterrupt`.
Its `except` clauses (`prepare/data/sra/accession.py:~555-585`) do not catch it, so `<fastq>/.prepare-tmp/<ACC>_temp`
(often many GB) stays until that accession is downloaded again. Reproduced by the download agent.

Behaviour: `except KeyboardInterrupt:` calls `_settle_temp_folder(temp_path, keep_temp)` and re-raises.

**Files:** Modify `prepare/data/sra/accession.py` (`download_accession`); Test `tests/test_data_sra_download.py`.

- [ ] **Step 1:** `test_keyboard_interrupt_removes_the_build_folder_and_propagates`: patch `_run_download_tool`
  to create `temp_path/SRR1_1.fastq` and raise `KeyboardInterrupt`; call `download_accession` through the same
  entry the module's interrupt tests use; assert `pytest.raises(KeyboardInterrupt)` and the build folder is gone.
  A second case with `keep_temp=True` asserts the folder is kept.
- [ ] **Step 2:** Run; expect FAIL (folder remains).
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run `conda run -n metaquest pytest tests/test_data_sra_download.py tests/test_concurrency_processes.py -q`; expect PASS.
- [ ] **Step 5:** Commit `fix: remove the build folder when the retry pass is interrupted`.

---

### Task 4: Extraction uses the selection's threshold by default

Defect: `extract_target_reads --threshold` defaults to 0.1 (`cli/commands/read_extraction.py:117-121`) whatever
threshold `select_datasets` used, so a sample selected at 0.05 is skipped with only an INFO line; `status --next`
(`cli/commands/status/suggest.py:237-253`, `_would_extract`) repeats the gap.

Behaviour:
- New `selection_threshold(registry, genome_id) -> Optional[float]` in `prepare/data/project_selection.py`: the
  lowest `criteria["threshold"]` among selected accessions whose criteria cover `genome_id` (`_covers`), None
  when there is none.
- `--threshold` default becomes None. Resolved once at the start of `execute`: the explicit value, else
  `selection_threshold(...)` unless `--ignore-selection`, else `DEFAULT_CONTAINMENT_THRESHOLD`. The resolved
  value is logged at INFO and is the one recorded in the extraction record (`"threshold"`).
- `_would_extract` and the suggested command use the same resolution; the command line includes
  `--threshold X` when X differs from the default.
- When an explicit `--threshold` leaves out selected accessions, one WARNING names how many.

**Files:**
- Modify: `prepare/data/project_selection.py`, `prepare/cli/commands/read_extraction.py`,
  `prepare/cli/commands/status/suggest.py`
- Modify: README.md and `docs/reads.md` (the `--threshold` description), CHANGELOG.md
- Test: `tests/test_project_selection.py`, `tests/test_cli_read_extraction.py`, `tests/test_cli_status.py`

**Interfaces:**
- Produces: `selection_threshold(registry: Registry, genome_id: str) -> Optional[float]`.

- [ ] **Step 1: Write the failing tests**
  - `test_selection_threshold_is_the_lowest_covering_threshold`: SRR1 `{column: G, threshold: 0.05}`,
    SRR2 `{column: G, threshold: 0.2}`, SRR3 `{column: H, threshold: 0.01}` -> `selection_threshold(reg, "G") == 0.05`.
  - `test_threshold_without_a_recorded_selection_threshold_is_the_default`: criteria `{column: G}` -> None.
  - CLI: a table with SRR9 at 0.07 for G, selected at threshold 0.05; `extract_target_reads` without
    `--threshold` (tools mocked as in the module) extracts SRR9 and records `threshold == 0.05`.
  - CLI: with `--threshold 0.1` it skips SRR9 and logs a WARNING containing `1 selected`.
  - status: `status --next` for the same project suggests a command containing `--threshold 0.05` and lists SRR9.
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run `conda run -n metaquest pytest tests/test_project_selection.py tests/test_cli_read_extraction.py tests/test_read_extraction.py tests/test_cli_status.py -q`; expect PASS.
- [ ] **Step 5:** Commit `fix: extract_target_reads uses the selection's threshold unless --threshold is given`.

---

### Task 5: Re-screening drops entries for accessions that left the matches

Defect: `record_screening_table` (`prepare/data/screening_table.py`) only drops cells still in the table
(`_drop_zero_cells`), and `branchwater_search` (`cli/commands/branchwater_search.py:109-116`) only adds. An
accession that is no longer in the new match set keeps its old screening entry, so `results_table`, `status`
and `_would_extract` still use it. Reproduced by the pipeline agent.

Behaviour: new `drop_vanished_screening(registry, genome_id, present: AbstractSet[str]) -> List[str]` in
`screening_table.py` removes the `genome_id` screening entry of every accession not in `present` (through
`_drop_screening`) and returns them. `record_screening_table` calls it once per genome column of the new
table with the table's index; `branchwater_search` calls it with the accessions it found. One INFO line per
genome gives the count. Other blocks of the dataset (download, extraction) are untouched.

**Files:** Modify `prepare/data/screening_table.py`, `prepare/cli/commands/branchwater_search.py`;
docs: `docs/project_state.md` (one sentence on re-screening); Test `tests/test_screening_table.py`
(`grep -ln record_screening_table tests/`), `tests/test_cli_branchwater_search.py` (or the module testing it).

- [ ] **Step 1: Write the failing tests**
  - `test_rescreening_drops_accessions_missing_from_the_new_table`: record a table G = {SRR1: .9, SRR2: .8}, then
    one with G = {SRR1: .9}; `query(registry, "screened", "G") == ["SRR1"]`.
  - `test_rescreening_leaves_genomes_absent_from_the_table_alone`: first table has G and H; second has only G
    (SRR1); SRR2's H entry stays.
  - `test_rescreening_keeps_the_download_record`: SRR2 downloaded; after it vanishes, its download block stays.
  - `test_branchwater_search_drops_vanished_matches`: two runs with mocked matches [SRR1, SRR2] then [SRR1].
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the two modules plus `tests/test_data_registry.py tests/test_processing_results.py`; expect PASS.
- [ ] **Step 5:** Commit `fix: re-screening removes entries for accessions no longer matched`.

---

### Task 6: `sra_profile`/`sra_report` skip excluded and incomplete folders when discovering

Defect: with no `--accession`/`--accessions-file`, `_resolve_accessions` (`cli/commands/sra_profile.py:144-158`)
profiles every accession-named folder, including blacklisted accessions and folders whose sidecar state is
`partial`, `failed` or `downloading`; they enter totals, grades and registry analyses. Reproduced.

Behaviour: in the discovery branch only, leave out folders for which `accession_has_fastq(folder / name)` is
False and accessions the registry excludes (`query(load_registry(args.registry), "excluded")`, when a registry
is found). One INFO line gives how many of each were left out. Named accessions are profiled as given.
Check whether `sra_report` uses the same discovery; if it has its own, apply the same rule there.

**Files:** Modify `prepare/cli/commands/sra_profile.py` (and `sra_report.py` if needed); `docs/reads.md`;
Test `tests/test_cli_sra_profile.py` (`grep -ln "SraProfileCommand\|sra_profile" tests/`).

- [ ] **Step 1:** `test_discovery_skips_blacklisted_and_partial_folders`: folders SRR1 (excluded in the registry),
  SRR2 (sidecar state partial), SRR3 (complete); profiling is mocked; assert only SRR3 is profiled.
  `test_named_accessions_are_profiled_even_if_excluded`: `--accession SRR1` profiles SRR1.
- [ ] **Step 2:** Run; expect the first to FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the module and `tests/test_cli_sra_report.py`; expect PASS.
- [ ] **Step 5:** Commit `fix: sra_profile discovery leaves out excluded accessions and incomplete downloads`.

---

### Task 7: `status` leaves excluded accessions out of "wanted" and applies the disk check to `--stage` and Gaps

Defects:
- `_resolve_wanted` (`processing/status_report.py:69-73`) returns `query(registry, "selected")` without removing
  excluded accessions, so an excluded accession is listed under "Missing FASTQ". Reproduced.
- `_print_stage_filter` and `_print_gaps` (`cli/commands/status/render_text.py:216-236`) query the raw
  registry, so `--stage downloaded` lists an accession that the Stages table reports as gone from disk. Reproduced.

Behaviour:
- `_resolve_wanted` without explicit inputs returns selected minus excluded. Explicit `--accessions-file`/
  `--parsed-containment` lists are left as given.
- `_print_stage_filter` and `_print_gaps` take the report and leave out the accessions in
  `report["stages"][stage]["missing_on_disk"]` (for Gaps: downloaded minus `missing_on_disk` of `downloaded`)
  and excluded accessions in "Downloaded but not extracted". `--json` output of `--stage` (if any) follows
  the same rule.

**Files:** Modify `prepare/processing/status_report.py`, `prepare/cli/commands/status/render_text.py`;
Test `tests/test_cli_status.py`.

- [ ] **Step 1: Write the failing tests**
  - `test_excluded_selected_accession_is_not_missing`: SRR2 selected then excluded, no reads; `status --json`
    -> `wanted.total == 1` and SRR2 not in `fastq_missing`.
  - `test_stage_filter_leaves_out_accessions_gone_from_disk`: three downloads, SRR2's folder removed;
    `status --stage downloaded` prints SRR1 and SRR3 only.
  - `test_stage_filter_with_a_genome_leaves_out_accessions_gone_from_disk`: same with `--genome G` and stage `extracted`.
  - `test_gaps_leave_out_excluded_accessions`.
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run `conda run -n metaquest pytest tests/test_cli_status.py tests/test_project_funnel.py -q`; expect PASS.
- [ ] **Step 5:** Commit `fix: status leaves excluded and deleted accessions out of wanted, --stage and Gaps`.

---

### Task 8: Funnel figures and `project_report` use the same disk-checked sets as `status`

Defects: `funnel` (`processing/project_funnel.py:79-97`) counts accessions from the disk-checked `members` but
sums bytes, pairs, seconds and total_bp over every registry record ("2 downloaded (3 bytes)" after one was
deleted). `project_report` (`processing/project_report.py:151`, `:~341`) uses no disk check at all, so its counts
differ from `status` once files are removed.

Behaviour:
- `funnel` sums download bytes and seconds only for accessions in `members["downloaded"]` (failed downloads
  keep being counted from every record), extraction pairs/seconds only for accessions in `members["extracted"]`,
  assembly pairs/seconds/total_bp only for accessions in `members["assembled"]`.
- New `disk_checked_members(registry, paths) -> Dict[str, List[str]]` in `processing/status_report.py`, used by
  `build_report` and by `project_report` (with `genome_counts(registry, on_disk=True)`), so both report the
  same counts. `project_report` obtains `ProjectPaths` the way `status` does (check `cli/commands/status/command.py`).
  It still opens no FASTQ file (the checks are existence checks; say so in its docstring).

**Files:** Modify `prepare/processing/project_funnel.py`, `prepare/processing/status_report.py`,
`prepare/processing/project_report.py`, `prepare/cli/commands/project_report.py` (if it builds paths);
`docs/project_state.md` or the project_report doc; Test `tests/test_project_funnel.py`, `tests/test_project_report*.py`.

**Interfaces:**
- Produces: `disk_checked_members(registry: Registry, paths: ProjectPaths) -> Dict[str, List[str]]`.

- [ ] **Step 1: Write the failing tests**
  - `test_funnel_bytes_follow_the_downloaded_members`: three 1-byte downloads, members list two -> bytes 2.
  - `test_funnel_pairs_follow_the_extracted_members`.
  - `test_project_report_counts_match_status_after_a_deletion`: build a project, delete SRR2's folder;
    `project_report`'s funnel downloaded count equals `status`'s.
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run `conda run -n metaquest pytest tests/test_project_funnel.py tests/test_cli_status.py $(ls tests/test_*project_report*.py) -q`; expect PASS.
- [ ] **Step 5:** Commit `fix: funnel figures and project_report use the disk-checked stage sets`.

---

### Task 9: `store_unlink` removes the project's usage rows

Defect: `StoreUnlinkCommand.execute` (`cli/commands/store/link.py:~215-255`) removes the link and updates the
registry but leaves the project's `usage` rows, so `store_gc` keeps the dataset "in use by a live project"
for as long as the project exists. Reproduced by the store agent.

Behaviour: after the registry update, one `catalog_write` session calls a new
`Catalog.forget_usage(accession, project_id) -> int` for each removed accession, and the journal lines for those
(accession, project) pairs are pruned with an extended `journal.prune_usage(paths, accessions, project_id=None)`
(per-project when `project_id` is given), so `store_reindex` does not bring the rows back (store rule 4).
A failure is reported the way Task 4 of 0.12.3 does: link removed, registry updated, exit 4 with a message
naming `store_unlink` to run again. A dataset the project never used is not an error.

**Files:** Modify `prepare/store/catalog.py`, `prepare/store/journal.py`, `prepare/store/usage.py`
(a `forget_usage_many(paths, registry, accessions) -> bool` wrapper beside `record_usage_many`),
`prepare/cli/commands/store/link.py`; `docs/shared_store.md`; Test `tests/test_store_usage.py`,
`tests/test_store_journal.py`, `tests/test_cli_store.py`, `tests/test_store_gc.py`.

**Interfaces:**
- Produces: `Catalog.forget_usage(accession: str, project_id: str) -> int` (rows removed);
  `prune_usage(paths, accessions, project_id: Optional[str] = None) -> int`;
  `forget_usage_many(paths: Optional[StorePaths], registry: Registry, accessions: Sequence[str]) -> bool`.

- [ ] **Step 1: Write the failing tests**
  - `test_unlink_then_gc_offers_the_dataset`: one project links SRR1 (usage recorded), `store_unlink SRR1`, then
    `store_gc --json` lists SRR1 as a candidate (set the `.used` marker age past the grace period as the gc
    tests do).
  - `test_reindex_after_unlink_does_not_restore_the_usage`.
  - `test_unlink_keeps_another_projects_usage`: two projects use SRR1; one unlinks; SRR1 stays in use.
  - `test_unlink_without_a_usage_row_succeeds`.
- [ ] **Step 2:** Run; expect the first three to FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run `conda run -n metaquest pytest tests/test_store_usage.py tests/test_store_journal.py tests/test_cli_store.py tests/test_store_gc.py -q`; expect PASS.
- [ ] **Step 5:** Commit `fix: store_unlink forgets the project's use of the dataset`.

---

### Task 10: Numeric flags reject meaningless values

Defects (CLI agent): `store_gc --older-than -30` silently turns off the age filter (`cli/commands/store/gc.py:110`);
`extract_target_reads --threads`, `--assembly-threads` and `--min-contig-len` accept 0 and negatives
(`read_extraction.py:130,168,185`); `single_sample --top-n -5` returns all rows but the last five (`samples.py:55`);
`--registry-max-screened -1` caps nothing (`branchwater_search.py:60`, `containment.py:63`, `plugins/sourmash_plugin.py:71`);
`download_metadata --batch-size 0`/`600` fails as a raw `ValueError` with exit 1 (`metadata.py:96`).

Behaviour: argparse types, so a bad value is a usage error (exit 2) with a clear message:
- `positive_int` (exists in `cli/base.py`) for `--threads`, `--assembly-threads`, `--min-contig-len`,
  `single_sample --top-n`, `--registry-max-screened`.
- new `non_negative_int` in `cli/base.py` for `store_gc --older-than` (0 = any age).
- new `int_in_range(low, high)` factory in `cli/base.py` for `--batch-size` (1-500, the bound already enforced
  in `data/metadata.py:168`).
- Replace `select.py`'s private `_positive_int` with `positive_int`.

**Files:** Modify `prepare/cli/base.py`, the command modules named above; Test `tests/test_cli_arguments.py`
(create if no module tests parser validation: `grep -ln "positive_int" tests/`).

- [ ] **Step 1:** One parametrised test `test_bad_numeric_values_are_usage_errors` that builds the real parser
  (`cli/main.py`) and asserts `SystemExit` with code 2 for each of: `store_gc --older-than -1`,
  `extract_target_reads --threads 0 ...`, `--assembly-threads -2`, `--min-contig-len 0`, `single_sample --top-n -5`,
  `branchwater_search --registry-max-screened -1`, `parse_containment --registry-max-screened 0`,
  `download_metadata --batch-size 0`, `--batch-size 600`. Plus `test_older_than_zero_is_accepted`.
  Fill each command's other required arguments with dummies.
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the new module and `tests/test_cli_commands.py tests/test_cli_select.py tests/test_store_gc.py -q`; expect PASS.
- [ ] **Step 5:** Commit `fix: reject zero and negative values for counts, thread numbers and limits`.

---

### Task 11: `sra_report --groups-file` validates its accessions; `genome_download --dry-run` creates nothing

Defects:
- `load_groups` (`cli/commands/sra_report.py:46-66`) checks only types, so `{"g": ["../outside/x"]}` is read
  from outside `--fastq-folder` and recorded in the registry under that name. Reproduced.
- `genome_download --dry-run` creates `--output-dir` (`cli/commands/genome.py:~203`) before its dry-run return. Reproduced.

Behaviour:
- `load_groups` strips each entry and raises `ValidationError` naming every entry that fails
  `validate_accession`.
- `genome_download` creates the output folder only after the dry-run return; `partition_present_genomes`
  already treats a missing folder as "nothing present".

**Files:** Modify `prepare/cli/commands/sra_report.py`, `prepare/cli/commands/genome.py`;
Test `tests/test_cli_sra_report.py`, `tests/test_genome_download.py` (or the CLI module that tests `genome_download`).

- [ ] **Step 1:** `test_groups_file_with_a_path_is_refused` (exit 1 via `self.fail`, message names `../outside/x`),
  `test_groups_file_entries_are_stripped`, `test_genome_download_dry_run_creates_no_folder`.
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the two modules; expect PASS.
- [ ] **Step 5:** Commit `fix: validate sra_report group accessions; genome_download --dry-run creates no folder`.

---

### Task 12: GC-only anomalies are flagged; one adapter threshold

Defects (analysis agent):
- `_detect_profile_anomalies` gives `unusual_gc` weight 0.2 and `detect_dataset_anomalies` flags a dataset only
  when `severity > 0.2` (`prepare/sra/analytics.py:~164,~496`), so a dataset whose only anomaly is GC outside
  20-80 % appears under `anomaly_types` but never in `anomalous_datasets`.
- The anomaly flag uses adapter > 0.1 while the report's "Contamination issues" tile counts > 0.05
  (`prepare/sra/reporting.py:~341`).

Behaviour: the flag test becomes `severity >= 0.2` (every single anomaly flags its dataset). One constant
`HIGH_ADAPTER_FRACTION = 0.05` in `analytics.py`, used by the anomaly flag and the report tile;
`recommend_processing_params`' own 0.01/0.1 stringency steps are unchanged.

**Files:** Modify `prepare/sra/analytics.py`, `prepare/sra/reporting.py`; Test `tests/test_sra_analytics*.py`
(`grep -ln detect_dataset_anomalies tests/`).

- [ ] **Step 1:** `test_gc_only_anomaly_is_flagged` (GC 85 %, nothing else -> in `anomalous_datasets`),
  `test_adapter_at_seven_percent_is_flagged_and_counted` (anomaly entry exists and the summary tile count is 1).
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the analytics and reporting tests; expect PASS.
- [ ] **Step 5:** Commit `fix: flag single anomalies and use one adapter-contamination threshold`.

---

### Task 13: `sra_profile` totals leave out sampled profiles; metadata read length per read

Defects:
- `sra_profile._summary_stats` (`cli/commands/sra_profile.py:~218-225`) sums `total_reads`/`total_bases` over all
  profiles, including `sampled=True` ones whose "total" is the sample size; `sra_report` already leaves them out
  and counts them as `sampled_datasets`.
- `SRAMetadataClient` sets `avg_length = bases / spots` (`data/sra_metadata.py:~396`), which for a PAIRED run is
  twice the read length (300 for 2x150) and disagrees with `sra_profile`'s `avg_read_length`.

Behaviour:
- `_summary_stats` sums only profiles with `sampled` False and adds `"sampled_datasets": <count>`; the text
  summary prints "N dataset(s) sampled only, not counted" when N > 0. The run-log summary gains the same key.
- `avg_length = bases / (spots * 2)` when `layout` is PAIRED, else `bases / spots`; the docstring of
  `SRADatasetInfo.avg_length` says "mean length of one read".

**Files:** Modify `prepare/cli/commands/sra_profile.py`, `prepare/data/sra_metadata.py`; CHANGELOG;
Test `tests/test_cli_sra_profile.py`, `tests/test_sra_metadata_client.py`.

- [ ] **Step 1:** `test_summary_leaves_out_sampled_profiles` (one exact 1000-read profile + one sampled 10000 ->
  total_reads 1000, sampled_datasets 1); `test_avg_length_of_a_paired_run_is_per_read` (spots 100, bases 30000,
  PAIRED -> 150.0); `test_avg_length_of_a_single_run_is_unchanged`.
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the two modules and `tests/test_sra_metadata_extended.py`; expect PASS.
- [ ] **Step 5:** Commit `fix: sra_profile totals leave out sampled datasets; metadata read length is per read`.

---

### Task 14: Gates, changelog and version 0.13.0

**Files:** `CHANGELOG.md` (`## [0.13.0] - <date>`, `### Changed` for Task 4, 5 and 10's behaviour changes,
`### Fixed` for the rest), `pyproject.toml`, `prepare/__init__.py`, `CLAUDE.md:7` (and correct the module
ceiling in CLAUDE.md's "Quality gates" from 800 to 700 lines).

- [ ] **Step 1:** `conda run -n metaquest make check`; expect exit 0.
- [ ] **Step 2:** `conda run -n metaquest make test`; expect all pass.
- [ ] **Step 3:** Changelog, version bump, CLAUDE.md corrections.
- [ ] **Step 4:** `conda run -n metaquest pytest tests/test_module_sizes.py tests/test_lint_exemptions.py tests/test_environment_pins.py -q`; expect PASS.
- [ ] **Step 5:** Commit `release: 0.13.0`. Push, PR and tag wait for the user.

---

## Not in this plan (Low findings, unplanned)

Analysis: `calculate_enrichment` and the `statistics.py` distance functions (no CLI caller), per-read GC
averaging, `validate_taxonomy` on an empty list, t-test with zero variance, unused Phred+64, Shannon pseudo-count.
Downloads: `genome_download` drops stderr, HTTP 400 treated as not-found, permanent HTTP errors retried,
mixed plain/gzip after a compression failure, SRX/SRS accepted by the accessions file. Store: `publish_folder`
two-rename window, replay timestamps, journal-before-commit, replacing a foreign symlink. Other: header-less
`--export-tsv`, NaN in the run log, `status --next` for 0-contig assemblies, `branchwater_search --output` genome
id, `genome_download` duplicate accessions, assembly flags ignored without `--assemble`, `logger.error`+`return 1`
usage errors.
