# Tech debt quick wins 0.11.1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the quick wins of `docs/superpowers/audits/2026-10-21-tech-debt.md`: the twelve
deferred review minors (D5), the store command preamble (D4), a removal version for the renamed
commands (D8), an index of the planning documents (D7) and three guards against regrowth
(audit section 4). Released as 0.11.1. Findings D1, D2, D3, D6 and D9 are left for the
medium-term plan the audit describes.

**Architecture:** Three waves. Wave 1 fixes the minors of the download layer, wave 2 the minors
of `sra_sketch`, wave 3 the structural items and the guards. Every fix starts with a failing
test. Each task is one commit on a branch `fix/tech-debt-quick-wins` from main.

**Tech Stack:** Python 3.12, pytest. Run everything as
`cd /Users/andreassjodin/Code/prepare && conda run -n metaquest <cmd>`.

**Spec:** `docs/superpowers/audits/2026-10-21-tech-debt.md`. Finding IDs (D1-D9) and minor
numbers (1-12) refer to that file.

## Global constraints

- ASCII only in `prepare/`, `tests/`, `scripts/`; no `print(` outside `cli/base.py`; every output
  file through the atomic helpers; modules under 800 lines (`registry.py`, `sra.py` and the other
  `SPLIT_MODULES` under 700); radon rank below D; black at 120.
- Errors through `self.fail` and the exception classes; no new exit codes.
- Tests never touch the network, a real tool, or a real store (tmp_path store, `PREPARE_DATA`
  and `XDG_CONFIG_HOME`/`HOME` monkeypatched).
- Commit trailer `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`; never `git add -A`.
- Modest scientific language in docstrings and docs.
- Gate before every commit: `conda run -n metaquest make check` and the owned tests; full
  `make test` at the end of each wave.
- No flag changes. Behaviour changes are listed in CHANGELOG under 0.11.1.

## Review focus

1. An interrupted `sra_sketch` run must still exit 130, record every accession that finished,
   and count the interrupted one neither as failed nor as sketched (Task 2.1).
2. The store sketch cache must be bypassed, not trusted, when the cached signature was built
   with other k sizes (Task 2.4).
3. The seven store commands must refuse a missing store with the same log line and the same
   JSON document as before (Task 3.1; `tests/test_cli_store.py` pins both).
4. The three guards must pass on the branch as it is, so the counts they freeze are the
   measured ones, not aspirations (Task 3.3).

---

## Wave 1: download layer minors

### Task 1.1: count plain files once; classify a raised NetworkError by its class (minors 1, 3)

**Files:** `prepare/data/sra/accession.py` (`_compress_files`), `prepare/data/sra/retry.py`
(`_collect_result`), `tests/test_data_sra_download.py`, `tests/test_data_sra_retry.py`.

- [ ] **Step 1: Failing tests.** `test_compress_files_counts_each_plain_file_once`: patch
  `fastq_mod.count_fastq_reads` with a counting wrapper, run `_compress_files` on two plain files
  with compression off, and assert each path was counted at most once (the verdict already
  counted it). `test_raised_network_error_without_keyword_is_classified_network`: a fake
  downloader raising `NetworkError("boom")` through `_execute_parallel_downloads`; the stored
  result message starts with `"network:"` and `classify_download_error` of it is `"network"`.
- [ ] **Step 2: Run** both with `-q`; FAIL.
- [ ] **Step 3: Implement.** In `_compress_files` take the pre-compression counts from the
  verdict's per-file figures when present instead of recounting; count only files that are
  actually compressed. In `_collect_result` prefix the message with the class of the exception
  when it is a `NetworkError` (`"network: ..."`) or the disk-full `OSError` (`"disk-full: ..."`),
  so `classify_download_error` sees a keyword whatever the text.
- [ ] **Step 4: Run** the two files; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: count plain files once; a raised NetworkError keeps its class`.

### Task 1.2: `--no-verify-downloads` leaves a shared sidecar's spot count alone (minor 2)

**Files:** `prepare/data/sra/store_handoff.py` (`_sidecar_ncbi`), `tests/test_data_sra_store.py`.

- [ ] **Step 1: Failing test.** `test_no_verify_keeps_the_xml_spot_count_in_the_sidecar`: a
  store download with `verify=False` where the metadata XML records `spots`; the written sidecar
  keeps `ncbi.spots` from the XML and its verdict is `unverified`.
- [ ] **Step 2: Run**; FAIL (`spots` dropped).
- [ ] **Step 3: Implement.** `_sidecar_ncbi` with `verify=False` keeps whatever the XML says and
  only stops filling `spots` from `expected_spots`; the verdict stays `unverified` because
  `build_sidecar` is given no expected count. Update the docstring: a count on file is a fact
  from NCBI, the verdict is what `--no-verify-downloads` turns off.
- [ ] **Step 4: Run** `tests/test_data_sra_store.py tests/test_store_verdicts.py`; PASS.
  `make check`.
- [ ] **Step 5: Commit** `fix: --no-verify-downloads keeps the XML spot count in the sidecar`.

### Task 1.3: `failed_accessions.txt` is removed under the run's lock (minor 4)

**Files:** `prepare/data/sra/retry.py` (the unlink near line 348), `tests/test_data_sra_retry.py`.

- [ ] **Step 1: Failing test.** `test_clean_run_removes_failed_file_only_when_it_wrote_none`:
  two runs in one project folder; the file written by run A's failures must survive run B's clean
  finish. Model with the existing fake downloader and a pre-existing file whose content run B did
  not write.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Remove the file only when this run wrote it (track the path in the
  run's state) or when its content equals this run's failure list; otherwise log at DEBUG that
  another run's file was left.
- [ ] **Step 4: Run** `tests/test_data_sra_retry.py`; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: a clean run removes only its own failed_accessions.txt`.

### Task 1.4: one process test with compression on (minor 5)

**Files:** `tests/helpers_processes.py`, `tests/test_concurrency_processes.py`.

- [ ] **Step 1: Failing test.** `test_download_with_compression_publishes_gz_files`: the real CLI
  as a process with the fake `pigz` (which writes a real gzip) and without `--no-compress`; the
  project folder holds `.fastq.gz` files, the registry record is `downloaded`, the verdict is
  `complete`.
- [ ] **Step 2: Run** `pytest tests/test_concurrency_processes.py -k compression -q`; it must
  fail first (write the assertion against a name the fake does not yet produce, then fix the
  fake), or document in the ledger why the fake already passes and keep the test.
- [ ] **Step 3: Implement.** Extend the fake `pigz` if needed so it honours `-p N` and writes
  `<file>.gz`.
- [ ] **Step 4: Run** the file; PASS.
- [ ] **Step 5: Commit** `test: a process download with compression on`.

**Wave gate:** `make test` green.

## Wave 2: `sra_sketch` minors

### Task 2.1: an interrupt is neither a failure nor a double warning; finished records survive a worker error (minors 6, 7, 9)

**Files:** `prepare/cli/commands/sra_sketch.py` (`_sketch_all`, `_collect`, `_run`),
`tests/test_cli_sra_sketch.py`.

- [ ] **Step 1: Failing tests.** `test_interrupt_counts_the_stopped_accession_as_neither`: stop
  after the first of three accessions; exit 130, the run log summary has `sketched 1`,
  `failed 0`, and the interrupt warning appears once in caplog.
  `test_worker_error_keeps_the_finished_records`: the second worker raises `RuntimeError`; the
  command exits 1, the first accession is recorded in the registry and the run log, and the
  error names the accession.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** In `_sketch_all` do not append the interrupted accession to
  `failed`; log the interrupt once (drop the duplicate in `_run`). In `_collect` catch
  `Exception` for the worker (re-raising `KeyboardInterrupt`), record the accession as failed
  with the message, and let `_run` record the finished ones before returning 1. Keep `_run`
  under rank D; extract `_finish` if needed.
- [ ] **Step 4: Run** `tests/test_cli_sra_sketch.py`; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: sra_sketch counts an interrupt once and keeps finished records`.

### Task 2.2: the store sketch lock wait honours the stop token (minor 8)

**Files:** `prepare/store/sketches.py` (`store_sketch`), `prepare/cli/commands/sra_sketch.py`
(`_sketch_one`), `tests/test_sra_sketch.py`.

- [ ] **Step 1: Failing test.** `test_store_sketch_stops_waiting_when_stop_is_set`: hold the
  dataset lock in the test, call `store_sketch` with a `should_stop` that returns True; it returns
  within a second with a warning rather than waiting `STATS_LOCK_WAIT_SECONDS`.
- [ ] **Step 2: Run**; FAIL (waits 30 s or `TypeError` on the argument).
- [ ] **Step 3: Implement.** `store_sketch(..., should_stop=None)` passes it to `dataset_lock`
  and treats `LockWaitStopped` like a held lock (warning, no cache). `_sketch_one` passes
  `stop.is_set`.
- [ ] **Step 4: Run** `tests/test_sra_sketch.py tests/test_cli_sra_sketch.py`; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: the store sketch lock wait stops with the run`.

### Task 2.3: `--signature` accepts `.sig.gz` or says why not (minor 10)

**Files:** `prepare/data/branchwater_search.py` (`load_signature`), `tests/test_branchwater_search.py`.

- [ ] **Step 1: Failing test.** `test_load_signature_reads_a_gzipped_file`: a `.sig.gz` written
  with `save_signatures_to_json(..., compression=1)` loads and yields the k=21 sketch (skip
  without sourmash).
- [ ] **Step 2: Run**; FAIL (traceback from `json.load` on gzip bytes).
- [ ] **Step 3: Implement.** Open with `gzip.open` when the first two bytes are the gzip magic;
  a file that is neither JSON nor gzip raises `DataAccessError("not a sourmash signature: ...")`.
- [ ] **Step 4: Run** the file and `tests/test_cli_sra_sketch.py`; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: --signature reads a .sig.gz file`.

### Task 2.4: the store sketch cache checks the k sizes (minor 11)

**Files:** `prepare/store/sketches.py` (`cached_sketch`), `tests/test_sra_sketch.py`.

- [ ] **Step 1: Failing test.** `test_cached_sketch_ignores_a_block_with_other_ksizes`: a sidecar
  `sketch` block whose `ksizes` is `[21]` with the files unchanged; `cached_sketch` returns None.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Compare `block.get("ksizes")` and `block.get("scaled")` with
  `KSIZES` and `SCALED` (imported from `prepare.sra.sketch` lazily to avoid a cycle, or pass
  them in).
- [ ] **Step 4: Run** `tests/test_sra_sketch.py`; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: the store sketch cache requires the same k sizes`.

### Task 2.5: non-accession folders under `fastq/` are skipped (minor 12)

**Files:** `prepare/cli/commands/sra_sketch.py` (the folder listing near line 138),
`tests/test_cli_sra_sketch.py`.

- [ ] **Step 1: Failing test.** `test_folders_that_are_not_accessions_are_skipped_with_a_note`:
  `fastq/` holds `SRR000001/` and `notes/`; without `--accession` the command sketches
  `SRR000001` only and logs `notes` at INFO as skipped.
- [ ] **Step 2: Run**; FAIL (exit 1 on `notes`).
- [ ] **Step 3: Implement.** Filter the listing with `validate_accession`
  (`core/validation.py`); log the skipped names once, joined.
- [ ] **Step 4: Run** `tests/test_cli_sra_sketch.py`; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: sra_sketch skips folders that are not accessions`.

**Wave gate:** `make test` green.

## Wave 3: structure, documents and guards

### Task 3.1: one store preamble for the seven store commands (D4)

**Files:** `prepare/cli/commands/store/_shared.py`, `prepare/cli/commands/store/{gc,status,usage,reindex,link,verify,adopt}.py`, `tests/test_cli_store.py`.

- [ ] **Step 1: Failing test.** `test_required_store_returns_root_or_the_refusal`: a new helper
  `_required_store(args, command) -> Tuple[Optional[Path], Optional[int]]` in `_shared.py`
  returns `(root, None)` with a store and `(None, 1)` plus the `_no_store_hint` log line (and the
  JSON document with `json=True`) without one; a `DataAccessError` from `load_registry` returns
  `(None, command.fail(...))`.
- [ ] **Step 2: Run**; FAIL (no helper).
- [ ] **Step 3: Implement.** Write the helper; replace the preamble in the seven commands with
  `root, code = _required_store(args, self); if root is None: return code`. No change to any
  log line or JSON shape. Check each module's line count falls.
- [ ] **Step 4: Run** `tests/test_cli_store.py tests/test_store_verdicts.py`; PASS. `make check`.
- [ ] **Step 5: Commit** `refactor: one store preamble for the store commands`.

### Task 3.2: a removal version for the renamed commands and an index of the planning documents (D8, D7)

**Files:** `prepare/cli/commands/renamed.py`, `CHANGELOG.md`, `docs/superpowers/README.md` (new),
`docs/README.md`.

- [ ] **Step 1: Failing test.** `test_renamed_commands_name_their_removal_version` in
  `tests/test_cli_main.py`: the module exposes `REMOVAL_VERSION = "0.13.0"` and each pointer's
  log line ends with `"; this pointer is removed in 0.13.0"`.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Add the constant and the sentence. Write
  `docs/superpowers/README.md`: one table of specs, plans and audits with the date, the release
  it shipped in (from CHANGELOG) and a status column (`done`, `open`), the two open items being
  this plan and the medium-term items of the tech debt audit. Link it from `docs/README.md` under
  a "Design records" line. No plan or spec file is moved.
- [ ] **Step 4: Run** `tests/test_cli_main.py` and `python scripts/check_docs_commands.py`; PASS.
- [ ] **Step 5: Commit** `docs: index of the design records; removal version for renamed commands`.

### Task 3.3: three guards against regrowth (audit section 4)

**Files:** `tests/test_module_sizes.py`, `tests/test_lint_exemptions.py` (new),
`tests/test_layer_boundaries.py` (new).

- [ ] **Step 1: Write the three tests.**
  `test_lint_exemption_counts_do_not_grow`: parse `setup.cfg`'s `per-file-ignores`, count the
  entries carrying `B902` and the entries carrying `D101`, and assert each is at most the
  constant measured on this branch (`B902_EXEMPT`, `DOCSTRING_EXEMPT`); the docstring says the
  constants only go down.
  `test_only_the_cli_imports_the_cli`: walk `prepare/` with `ast`, and assert no module outside
  `prepare/cli/` imports `prepare.cli`, except the entries of `KNOWN_CLI_IMPORTERS = {"prepare/plugins/sourmash_plugin.py"}`
  (removed by the medium-term plan's D2 task).
  `test_rank_c_count_does_not_grow` in `tests/test_module_sizes.py`: radon's `cc_visit` over the
  package, count blocks with complexity 11-20, assert at most `MAX_RANK_C` (the measured 105);
  the constant only goes down.
- [ ] **Step 2: Run** the three; PASS at the measured values. Then raise each threshold by one
  temporarily and watch the test still pass, lower it by one and watch it FAIL, to prove the
  counts are read; restore the measured values.
- [ ] **Step 3: Record** the three constants in `CLAUDE.md` under "Quality gates" in one
  sentence each.
- [ ] **Step 4: Run** `make check`; PASS.
- [ ] **Step 5: Commit** `test: guards for lint exemptions, layer imports and rank-C count`.

### Task 3.4: docs and release 0.11.1

**Files:** `CHANGELOG.md`, `docs/reads.md` (sra_sketch section: `.sig.gz` accepted, folders
skipped), `docs/shared_store.md` (sketch cache and k sizes), `pyproject.toml`, `prepare/__init__.py`,
`CLAUDE.md` (version).

- [ ] **Step 1: Write** the 0.11.1 section: Fixed (minors 1-4, 6-12 in user terms), Changed
  (removal version for the renamed pointers), Internal (store preamble, guards, design record
  index). Move the CI note already under Unreleased into 0.11.1.
- [ ] **Step 2: Run** `python scripts/check_version.py`, `make check`, `make test`; all green.
- [ ] **Step 3: Commit** `release: 0.11.1`.

The branch is then ready for a fresh whole-branch review, merge into main and tag `v0.11.1` on
the user's instruction.
