# SRA download fixes 0.10.2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix every finding in `docs/superpowers/audits/2026-10-21-sra-download-audit.md` and close
the listed test gaps, released as 0.10.2.

**Architecture:** Four waves in the audit's order. Wave 1 gives the plain project and the store one
completeness rule and makes the verdict hold for the files that are published. Wave 2 carries the
real error text through the result path and makes the retry pass behave like the first pass on an
interrupt. Wave 3 removes the misleading failures and the leftovers. Wave 4 deletes dead code and
gives the fake tools failure modes so the process tests cover the paths above. Each task is one
commit on a branch `fix/sra-download-audit`.

**Tech Stack:** Python 3.12, pytest. Run everything as
`cd /Users/andreassjodin/Code/prepare && conda run -n metaquest <cmd>`.

**Spec:** `docs/superpowers/audits/2026-10-21-sra-download-audit.md`. Finding IDs (H1-H4, M1-M8,
L1-L9) refer to that file.

## Global constraints

- ASCII only in `prepare/`, `tests/`, `scripts/`; no `print(` outside `cli/base.py`; every output
  file through the atomic helpers; modules under 800 lines; radon rank below D; black at 120.
- Errors through `self.fail` and the exception classes; no new exit codes.
- Commit trailer `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`; never `git add -A`.
- Modest scientific language in docstrings and docs.
- Gate before every commit: `conda run -n metaquest make check` and the owned tests; full
  `make test` at the end of each wave.
- Behaviour changes are listed in CHANGELOG under 0.10.2. One flag is removed (`--verify-downloads`,
  a no-op) and one added (`--keep-temp`).

## Review focus

1. A truncated download in a plain project must give exit 1, a `failed` registry state with the
   verdict kept, and a retry; with `--accept-partial` it must behave as today (Task 1.1).
2. The message a worker raised must reach `download_results`, the report CSV, `failure_reason`
   and the exit code unchanged (Task 2.1).
3. An interrupt anywhere (download, compression, retry pass, back-off) must publish nothing new,
   record no attempt that did not run, and exit 130 (Tasks 1.4, 2.2).
4. Nothing in the store path changes its sidecar states or the hand-off order (every task).

---

## Wave 1: one completeness rule

### Task 1.1: truncated is a failure unless accepted (H1)

**Files:** `prepare/data/sra/accession.py` (`_handle_download_output`, `download_accession`,
`_project_download`), `prepare/data/sra/download.py` (the `_project_download` partial),
`prepare/data/sra/store_handoff.py` (`_store_fetch`), `prepare/cli/commands/sra.py`
(`--accept-partial` help), `tests/test_data_sra_download.py`.

- [ ] **Step 1: Failing tests.** In `tests/test_data_sra_download.py`:
  `test_truncated_download_is_a_failure_without_accept_partial` (`_handle_download_output` with
  `expected_spots` four times the reads returns `(False, msg)`, `msg` starts with the
  `"Downloaded N files, truncated (n of m spots)"` text so `parse_verdict_message` still reads it,
  and the files are left in `output_path`); `test_truncated_download_accepted_with_accept_partial`
  (`accept_partial=True` returns `True`); `test_project_download_truncated_records_failed_and_retries`
  through `download_sra` with a fake downloader: the accession is in `failed_accessions`, the
  retry pass is entered, and `on_result` received `success=False` with the verdict in the message.
- [ ] **Step 2: Run** `pytest tests/test_data_sra_download.py -k truncated -q`; FAIL.
- [ ] **Step 3: Implement.** `_handle_download_output(..., accept_partial: bool = False)`: after the
  verdict, when `verdict["verdict"] == "truncated" and not accept_partial`, skip compression, log
  the warning as today, and `return False, message`. `download_accession` takes
  `accept_partial: bool = False` and passes it on; `_discard_cached_archive` is not called on
  failure, so the archive is kept for the retry, which is wanted (prefetch resumes it).
  `download_sra` adds `accept_partial=accept_partial` to the `_project_download` partial.
  `_store_fetch` passes `accept_partial=True` to `download_accession`, since the sidecar states
  decide there (no behaviour change in the store). The `--accept-partial` help becomes
  "Keep a truncated download (plain project) or use an incomplete store dataset instead of
  failing". `record_download` already keeps the verdict for a `failed` state through
  `merged_verdict`; add an assertion to the test.
- [ ] **Step 4: Run** `pytest tests/test_data_sra_download.py tests/test_data_sra_retry.py tests/test_cli_sra.py -q`; PASS.
- [ ] **Step 5: Commit** `fix: a truncated download fails unless --accept-partial (audit H1)`.

### Task 1.2: `--no-verify-downloads` reaches the store path (H3)

**Files:** `prepare/data/sra/store_handoff.py` (`_resolve_expected_spots`, `_store_download`,
`_store_downloader`), `prepare/data/sra/download.py`, `prepare/cli/commands/sra.py`
(`_store_options`), `tests/test_data_sra_store.py`.

- [ ] **Step 1: Failing test.** `test_store_download_unverified_when_verification_is_off`: a store
  with a sidecar `ncbi.spots` and an XML, `_store_download(..., verify=False)` with a fake
  `download_accession` records the sidecar state `unverified` (`expected_spots` passed to the
  downloader is `None`) and links the dataset.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** `_resolve_expected_spots(accession, store, store_metadata, given, verify=True)`
  returns `None` when `verify` is false. `_store_download` and `_store_downloader` take
  `verify: bool = True`; `download_sra` takes `verify_downloads: bool = True` and passes it to
  `_store_downloader`; `_store_options` adds `"verify_downloads": getattr(args, "verify_downloads", True)`.
- [ ] **Step 4: Run** `pytest tests/test_data_sra_store.py tests/test_cli_sra.py -q`; PASS.
- [ ] **Step 5: Commit** `fix: --no-verify-downloads also skips the store verdict (audit H3)`.

### Task 1.3: verify the compressed files (H4)

**Files:** `prepare/data/sra/accession.py` (`_handle_download_output`),
`prepare/data/sra/fastq.py`, `tests/test_data_sra_download.py`.

- [ ] **Step 1: Failing test.** `test_short_gz_after_compression_is_a_failure`: patch
  `compress_fastq` to write a `.gz` holding half the records; `_handle_download_output` returns
  `(False, "compression corrupted <file>: n of m reads")` and the folder is not reported complete.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Count reads per file once before compression
  (`fastq_mod.count_fastq_reads`, already called by `verify_download` for the mate-1 file; keep
  one call per file by passing the counts into `verify_download` or counting once in
  `_handle_download_output`). After each successful `compress_fastq`, count the `.gz` file; a
  different count returns `(False, ...)`. The cost is one gzip read per file, the same as the
  store's `build_sidecar`; note in the docstring.
- [ ] **Step 4: Run** `pytest tests/test_data_sra_download.py -q`; PASS.
- [ ] **Step 5: Commit** `fix: the plain project verifies the compressed files (audit H4)`.

### Task 1.4: an interrupt during compression publishes nothing (M1)

**Files:** `prepare/data/sra/accession.py`, `tests/test_data_sra_download.py`.

- [ ] **Step 1: Failing test.** `test_interrupt_during_compression_is_not_published`: stop set after
  the first file compresses; `_handle_download_output` returns `(False, "interrupted")`; through
  `_project_download` with a fake tool, `fastq/<ACC>` does not exist afterwards and the prefetch
  archive is still in the cache.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Replace the `compression_skipped` branches with
  `return False, "interrupted"`. `download_accession` already keeps the archive on failure.
  Remove the `"compression skipped (interrupted)"` message and its tests. The store's
  `compression` label in `_store_fetch` can then stay as is (a published folder is uniform).
- [ ] **Step 4: Run** `pytest tests/test_data_sra_download.py tests/test_interrupt_e2e.py -q`; PASS.
- [ ] **Step 5: Commit** `fix: an interrupt during compression leaves the download unpublished (audit M1)`.
- [ ] **Wave gate:** `make check` and `make test`.

## Wave 2: error reporting

### Task 2.1: the worker's error text is the result (H2)

**Files:** `prepare/data/sra/retry.py` (`_collect_result`, `_process_download_results`),
`tests/test_data_sra_retry.py`.

- [ ] **Step 1: Failing test.** Replace the assertion at `tests/test_data_sra_retry.py:71` with
  `"SRR456": "Error: boom"` style expectations, and add
  `test_raised_network_error_is_classified_network`: a downloader raising `NetworkError("timed out")`
  ends with `download_results[acc]` containing the text, `classify_download_error` giving
  `network`, and the accession in the retry batch.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** `_collect_result` appends `(acc, (False, f"Error: {e}"))`. Delete the
  `except (TypeError, ValueError)` branch in `_process_download_results` and the docstring lines
  about `None`.
- [ ] **Step 4: Run** `pytest tests/test_data_sra_retry.py tests/test_cli_sra.py -q`; PASS.
- [ ] **Step 5: Commit** `fix: a worker's exception text is the recorded result (audit H2)`.

### Task 2.2: the retry pass handles an interrupt like the first pass (M4)

**Files:** `prepare/data/sra/retry.py` (`_retry_failed_downloads`), `tests/test_data_sra_retry.py`.

- [ ] **Step 1: Failing tests.** `test_retry_pass_drops_interrupted_results`: a downloader that sets
  `stop` and returns `(False, "interrupted")` on the first accession; the second is not attempted;
  neither is recorded through `on_result` as failed, and `retry_count` is 1.
  `test_retry_backoff_stops_on_interrupt`: with `max_retries=3`, set `stop` from a timer during
  the back-off; the function returns well before `2**retry` seconds (use a monkeypatched wait
  and assert it was called with the stop event).
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** In the loop, when `message == "interrupted"` do not append to
  `failed_accessions`, do not notify `on_result`, and `break` (the caller reports the rest as not
  completed, as the `stop_requested` branch does). Replace `time.sleep(2**retry)` with
  `stop.wait(2**retry)` (`stop` is a `threading.Event`, never None here).
- [ ] **Step 4: Run** `pytest tests/test_data_sra_retry.py -q`; PASS.
- [ ] **Step 5: Commit** `fix: the retry pass stops cleanly on an interrupt (audit M4)`.

### Task 2.3: a real project folder is refused before the store fetch (M7)

**Files:** `prepare/data/sra/store_handoff.py` (`_store_download`), `tests/test_data_sra_store.py`.

- [ ] **Step 1: Failing test.** `test_store_download_refuses_real_project_folder`: `fastq/<ACC>` is
  a real folder with files; `_store_download` returns `(False, "project folder exists: ...")`
  without calling the downloader, and `_split_not_found` keeps the accession out of the retry batch.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Add `STORE_PROJECT_FOLDER_PREFIX = "project folder exists:"` to
  `SETTLED_PREFIXES`; at the top of `_store_download`, when `project_path / accession` exists and
  is not a symlink and `force` is false, return the settled message naming `store_adopt`.
- [ ] **Step 4: Run** `pytest tests/test_data_sra_store.py tests/test_data_sra_retry.py -q`; PASS.
- [ ] **Step 5: Commit** `fix: a real project folder settles a store download before the fetch (audit M7)`.
- [ ] **Wave gate:** `make check` and `make test`.

## Wave 3: robustness

### Task 3.1: exact archive names and a prefetch that wrote nothing (M2, M3)

**Files:** `prepare/data/sra/accession.py` (`_cached_sra_archive`, `download_accession`),
`tests/test_data_sra_download.py`.

- [ ] **Step 1: Failing tests.** `test_cached_archive_ignores_prefetch_leftovers` (`.sra.lock`,
  `.sra.prf`, `.sra.tmp` beside `.sralite` gives the `.sralite`);
  `test_prefetch_without_archive_is_a_network_failure` (fake prefetch writes nothing; the result
  is `network: prefetch wrote no archive for <ACC>` and fasterq-dump is not run).
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** `_cached_sra_archive` checks `<ACC>.sra` then `<ACC>.sralite` only
  and returns `None` otherwise; `download_accession` raises `DownloadError` (or the existing
  network-classified exception type in `_DOWNLOAD_ERRORS`) when it is `None`. Add
  `"prefetch wrote no archive"` to the network patterns of `classify_download_error`.
- [ ] **Step 4: Run** `pytest tests/test_data_sra_download.py -q`; PASS.
- [ ] **Step 5: Commit** `fix: prefetch leftovers and empty prefetch results are handled (audit M2, M3)`.

### Task 3.2: scratch folders of failed attempts (M5)

**Files:** `prepare/data/sra/accession.py`, `prepare/data/sra/download.py`,
`prepare/cli/commands/sra.py` (`--keep-temp`), `tests/test_data_sra_download.py`, `tests/test_cli_sra.py`.

- [ ] **Step 1: Failing tests.** `test_failed_attempt_removes_temp_folder` and
  `test_keep_temp_keeps_the_temp_folder_and_logs_its_path`.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** `download_accession(..., keep_temp: bool = False)`: in the three
  failure branches, `_safe_rmtree(temp_path)` unless `keep_temp`, in which case log the path at
  INFO. Thread `keep_temp` through `download_sra` and the CLI flag `--keep-temp`
  ("Keep fasterq-dump's build folder of a failed attempt for inspection").
- [ ] **Step 4: Run** `pytest tests/test_data_sra_download.py tests/test_cli_sra.py -q`; PASS.
- [ ] **Step 5: Commit** `fix: a failed attempt removes its build folder unless --keep-temp (audit M5)`.

### Task 3.3: a failed forced redownload keeps the file record (M6)

**Files:** `prepare/data/registry.py` (`record_download`), `tests/test_data_registry.py`.

- [ ] **Step 1: Failing test.** `test_failed_record_keeps_files_still_on_disk`: a `downloaded`
  record, then `record_download(..., "failed")` with the folder still present keeps `files` and
  `bytes_total`; with the folder gone they are cleared.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** For `state == "failed"`, when `fastq_files(folder)` is non-empty,
  keep the existing `files`/`bytes_total`; document in the docstring.
- [ ] **Step 4: Run** `pytest tests/test_data_registry.py tests/test_cli_status.py -q`; PASS.
- [ ] **Step 5: Commit** `fix: a failed redownload keeps the record of the files on disk (audit M6)`.

### Task 3.4: worker and thread bounds (M8)

**Files:** `prepare/cli/commands/sra.py` (`_resolve_max_workers`, `--num-threads` type),
`prepare/utils/security.py`, `tests/test_cli_sra.py`, `tests/test_security_comprehensive.py`.

- [ ] **Step 1: Failing tests.** `--max-workers 64` resolves to `MAX_CONCURRENT_DOWNLOADS` with a
  warning; `--num-threads 0` is a parser error; `run_secure` rejects `--threads 0` for fasterq-dump.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Cap in `_resolve_max_workers`; a `positive_int` argparse type for
  `--num-threads` and `--max-workers` (check `cli/base.py` for an existing helper); the integer
  check in `security.py` requires a value above zero for `--threads` and `-p`.
- [ ] **Step 4: Run** the two test files; PASS.
- [ ] **Step 5: Commit** `fix: bound --max-workers and reject zero threads (audit M8)`.

### Task 3.5: small behaviours (L1-L4)

**Files:** `prepare/data/sra/retry.py`, `prepare/cli/commands/sra.py`,
`prepare/data/sra/store_handoff.py`, `prepare/data/sra/download.py`, tests as named.

- [ ] **Step 1: Failing tests.** `test_clean_run_removes_stale_failed_file` (L1);
  `test_retry_hint_logged_once` with `caplog` (L2); `--verify-downloads` is unknown to the parser
  (L3); `test_store_redownload_truncated_only_for_truncated_set` (L4: the store downloader gets
  `redownload_truncated` only for accessions in the truncated set).
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** `_handle_download_failure` unlinks the file when there are no
  failures; drop the hint from `retry.py` (the CLI one names the folder the user gave); remove
  the `--verify-downloads` argument; `download_sra` passes `truncated` to the store downloader
  the way it does to `_project_download` (apply the set in `_store_downloader._download`).
- [ ] **Step 4: Run** `pytest tests/test_data_sra_retry.py tests/test_cli_sra.py tests/test_data_sra_store.py -q`; PASS.
- [ ] **Step 5: Commit** `fix: stale failed file, duplicated hint, no-op flag, store truncated set (audit L1-L4)`.
- [ ] **Wave gate:** `make check` and `make test`.

## Wave 4: dead code and the missing process tests

### Task 4.1: dead and duplicated code (L5-L7)

**Files:** `prepare/core/constants.py`, `prepare/utils/security.py`,
`prepare/data/sra/store_handoff.py`, `prepare/data/sra/retry.py`,
`prepare/cli/commands/sra_verdicts.py`, `prepare/data/sra/accession.py`.

- [ ] **Step 1:** Remove `DEFAULT_MAX_RETRIES`, `FASTERQ_DUMP_BOOLEAN_FLAGS`,
  `MAX_SUBPROCESS_TIMEOUT`, `_publish_store_dataset`; make `downloader` a required argument of
  `_download_with_retries` and `_retry_failed_downloads` and update the tests that relied on the
  fallback. Keep one `_read_count` in `prepare/data/sra/spots.py` (the stricter type rule) and
  import it from both callers. Cache `fasterq_dump_version` with `functools.lru_cache` (clear it
  in the tests that patch `probe_tool`).
- [ ] **Step 2: Run** `pytest tests/test_data_sra_retry.py tests/test_data_sra_store.py tests/test_cli_sra.py tests/test_security_comprehensive.py -q`; PASS. `grep -rn` for each removed name returns nothing.
- [ ] **Step 3: Commit** `refactor: remove dead download helpers and unify the read-count parser (audit L5-L7)`.

### Task 4.2: fake tools with failure modes and the process tests

**Files:** `tests/helpers_processes.py`, `tests/test_concurrency_processes.py`.

- [ ] **Step 1:** Extend the fake tools: `FAKE_FASTERQ_EXIT` (exit code and a stderr line),
  `FAKE_FASTERQ_SPOTS` (records per mate, default 4), `FAKE_FASTERQ_UNPAIRED` (also write
  `<ACC>.fastq`), `FAKE_PREFETCH_NOOP` (exit 0, write nothing), `FAKE_FASTERQ_FAIL_ONCE`
  (fail the first call per accession, recorded in a marker file, succeed on the retry).
- [ ] **Step 2:** Add tests, each a real `prepare download_sra` process with the isolated
  environment: `test_truncated_download_fails_in_plain_project` (spots 4, registry
  `run_total_spots` 100; exit 1, registry `failed` with verdict `truncated`,
  `failed_accessions.txt` lists it); `test_truncated_download_accepted`; `test_tool_failure_then_retry`
  (exit 0 on the second attempt, `attempts` 2, the first attempt's stderr text in the run log
  detail); `test_prefetch_noop_is_network_failure` (exit 4); `test_unpaired_file_is_recorded`
  (`files` lists three entries).
- [ ] **Step 3: Run** `pytest tests/test_concurrency_processes.py -q`; PASS.
- [ ] **Step 4: Commit** `test: fake sra-tools with failure modes and process tests for the download verdicts`.

### Task 4.3: docs and release

- [ ] `CHANGELOG.md` 0.10.2: Changed (truncated fails unless `--accept-partial`; interrupt during
  compression leaves nothing published; compressed files verified; `--no-verify-downloads` with a
  store), Added (`--keep-temp`), Removed (`--verify-downloads`), Fixed (the rest by finding ID).
- [ ] `docs/reads.md` download section: the completeness rule, `--accept-partial`, `--keep-temp`;
  `docs/ARCHITECTURE.md` "Exit-code policy" if it lists the truncated case; README flag mentions
  (`python scripts/check_docs_commands.py`).
- [ ] Version 0.10.2 in `pyproject.toml` and `prepare/__init__.py` (see `scripts/check_version.py`).
- [ ] `make check`, `make test`, then `superpowers:finishing-a-development-branch`.
