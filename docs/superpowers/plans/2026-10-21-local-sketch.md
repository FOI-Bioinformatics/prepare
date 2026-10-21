# sra_sketch Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A new Reads command `sra_sketch` that builds a sourmash signature (k = 21, 31, 51,
scaled 1000, abundance) per downloaded accession, caches it in the shared store, records it in the
registry and run log, and writes the containment of a genome into the matches folder.

**Architecture:** Four layers, each one task: a block reader in `data/sra/sampling.py`; the
sketcher and writer in `sra/sketch.py` with the containment and CSV merge in
`sra/sketch_containment.py`; the store cache in `store/sketches.py` with one optional sidecar
field; the command in `cli/commands/sra_sketch.py` copied from `sra_profile`. Branch
`feat/sra-sketch`, one commit per task, released as 0.11.0.

**Tech Stack:** Python 3.12, numpy, sourmash 4.8 or later (the `sourmash` extra), pytest. Run
everything as `cd /Users/andreassjodin/Code/prepare && conda run -n metaquest <cmd>`.

**Spec:** `docs/superpowers/specs/2026-10-21-local-sketch-design.md`.

## Global constraints

- Sourmash only through `prepare.core.optional.require("sourmash", "sourmash", purpose)` at the
  point of use; `tests/test_optional_imports.py` must still build the parser with it blocked.
- Every file through `open_atomic` / `write_text_atomic` / `copy_file` (`data/file_io.py`); keep
  `open_atomic(path, "wb")` on one line for `scripts/check_atomic_writes.sh`.
- ASCII only; no `print(`; modules under 800 lines; radon rank below D; black at 120.
- Store tests use a `tmp_path` store with `PREPARE_DATA` and `XDG_CONFIG_HOME`/`HOME`
  monkeypatched (store rule 2). The store cache never touches the catalogue tables (rule 4).
- All sourmash tests `pytest.importorskip("sourmash")` at module level or in the test.
- Commit trailer `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`; never `git add -A`.

## Review focus

1. For every k, the hashes and abundances of a dataset sketch equal those of per-read
   `add_sequence` calls over the same reads, for paired plus orphan files, plain and gzip,
   lowercase and `N` (Task 2).
2. No `.sig.gz` exists after a failure or an interrupt; a partial file never appears under its
   final name (Tasks 2, 4).
3. The merged matches CSV keeps other accessions' rows and metadata and `parse_containment_data`
   reads it unchanged (Task 2).
4. The store cache is written only under the dataset lock, is invalidated by a changed file, and
   a held lock costs only the cache (Task 3).

## File structure

New: `prepare/sra/sketch.py`, `prepare/sra/sketch_containment.py`, `prepare/store/sketches.py`,
`prepare/cli/commands/sra_sketch.py`, `tests/test_sra_sketch.py`, `tests/test_cli_sra_sketch.py`.
Modified: `prepare/data/sra/sampling.py`, `prepare/data/sra/__init__.py`,
`prepare/store/sidecar.py`, `prepare/cli/main.py`, `tests/test_optional_imports.py`,
`tests/test_store_sidecar.py`, docs listed in Task 5.

---

### Task 1: block reader `iter_sequence_blocks`

**Files:** `prepare/data/sra/sampling.py`, `prepare/data/sra/__init__.py`, `tests/test_sra_sketch.py`.

- [ ] **Step 1: Failing tests** (no sourmash needed): `test_blocks_join_sequences_with_n` (plain
  file of 3 reads gives one block `b"ACGT...NACGT...N..."` with `n == 3`);
  `test_blocks_read_gzip`; `test_blocks_strip_cr`; `test_blocks_accept_missing_final_newline`;
  `test_record_spanning_block_boundary` (monkeypatch `CHUNK_SIZE` to 64 on a 10-read file; the
  concatenation of blocks split on `N` equals the reads, and the counts sum to 10);
  `test_truncated_record_raises`; `test_not_fastq_raises`; `test_stop_token_raises_interrupted_error`.
- [ ] **Step 2: Run** `pytest tests/test_sra_sketch.py -k blocks -q`; FAIL.
- [ ] **Step 3: Implement.**
  ```python
  def iter_sequence_blocks(path, stop=None) -> Iterator[Tuple[bytes, int]]:
      """Yield (sequences joined by b"N", record count) per block of ``path``; see module docstring."""
  ```
  Open through `_open_binary`; read `CHUNK_SIZE`; prepend the carry; the first byte of the file
  must be `@`; `newline_at = np.flatnonzero(np.frombuffer(buf, np.uint8) == 10)`; complete
  records are `len(newline_at) // 4`; for record `i` the sequence is
  `buf[newline_at[4*i] + 1 : newline_at[4*i + 1]]` (drop a trailing `\r` when the first record
  has one); yield `b"N".join(seqs), n`; carry `buf[newline_at[4*n - 1] + 1:]`. At EOF, a carry
  that is only whitespace is fine; a carry with a header and four lines without the final
  newline is one more record; anything else raises `ValueError("Truncated FASTQ record ...")`.
  Check `stop.is_set()` before each read and raise `InterruptedError(f"interrupted while reading {path}")`.
  Export from `prepare/data/sra/__init__.py`.
- [ ] **Step 4: Run** the tests; PASS. `make check`.
- [ ] **Step 5: Commit** `feat: block reader that joins FASTQ sequences for sketching`.

### Task 2: sketcher, writer, containment

**Files:** `prepare/sra/sketch.py`, `prepare/sra/sketch_containment.py`, `tests/test_sra_sketch.py`.

- [ ] **Step 1: Failing tests** (importorskip sourmash):
  `test_sketch_equals_per_read_add_sequence` (paired plus orphan, plain and gzip, lowercase and
  `N`; for k in 21, 31, 51 compare `mh.hashes` dicts with a reference built by
  `MinHash(n=0, ksize=k, scaled=1000, track_abundance=True).add_sequence(read, force=True)` per
  read; `reads` and `bases` exact); `test_write_signature_is_sourmash_sig_gz` (gzip magic,
  `load_file_as_signatures` yields three signatures named by the accession, ksizes 21/31/51 in
  that order, scaled 1000, abundance; `gzip.decompress` then `json.loads` then
  `branchwater_search._select_sketch` finds k = 21); `test_write_signature_leaves_no_file_on_error`
  (patch `save_signatures_to_json` to raise; no `.sig.gz`, no temp file in the folder);
  `test_sketch_is_current_rules`; `test_load_dataset_minhash_flattens_and_checks_scaled`;
  `test_containment_matches_contained_by` (genome FASTA whose k-mers partly overlap the reads;
  `containment_rows` equals `genome_mh.contained_by(dataset_mh)`);
  `test_merge_keeps_other_rows_and_metadata` (an existing CSV with `BRANCHWATER_COLUMNS` and a
  filled `organism` for SRR9 and a stale row for SRR1; after the merge SRR9 is unchanged, SRR1 is
  replaced, the file is sorted, `cANI == round(c ** (1/21), 4)`);
  `test_merged_csv_round_trips_through_parse_containment_data`.
- [ ] **Step 2: Run** `pytest tests/test_sra_sketch.py -q`; FAIL.
- [ ] **Step 3: Implement `prepare/sra/sketch.py`.**
  ```python
  KSIZES = (21, 31, 51)
  SCALED = branchwater_search.SCALED            # assert branchwater_search.KSIZE in KSIZES
  SIGNATURE_SUFFIX = ".sig.gz"

  @dataclass
  class DatasetSketch:
      accession: str; files: List[Path]; reads: int; bases: int
      hashes: Dict[int, int]; signature: Dict[str, List[Union[int, float]]]; minhashes: Dict[int, Any]

  def signature_path(accession, output_dir) -> Path
  def sketch_dataset(accession, files, stop=None) -> DatasetSketch
  def write_signature(sketch, path) -> Path
  def sketch_summary(sketch, source="computed") -> Dict[str, Any]
  def sketch_is_current(summary, files, path) -> bool
  def load_dataset_minhash(path, ksize=branchwater_search.KSIZE)
  ```
  `sketch_dataset`: `sourmash = require("sourmash", "sourmash", "Sketching reads")`; one
  `MinHash(n=0, ksize=k, scaled=SCALED, track_abundance=True)` per k; for each file, for
  `(block, n) in iter_sequence_blocks(file, stop)`: `reads += n`, `bases += len(block) - (n - 1)`
  when `n`, `mh.add_sequence(block, force=True)` per k; log per file at `item_level`.
  `write_signature`: `sigs = [SourmashSignature(sketch.minhashes[k], name=sketch.accession) for k in KSIZES]`;
  `with open_atomic(path, "wb") as fp: save_signatures_to_json(sigs, fp, compression=1)`
  (`from sourmash.signature import save_signatures_to_json`, through `require`).
  `sketch_summary`: `{"ksizes": list(KSIZES), "scaled": SCALED, "track_abundance": True,
  "reads", "bases", "hashes": {str(k): n}, "files": signature, "source": source}`.
  `sketch_is_current`: `path.is_file()` and the summary's `files`, `ksizes`, `scaled` match
  (`_file_signature` imported from `prepare.store.stats`).
  `load_dataset_minhash`: `next(iter(sourmash.load_file_as_signatures(str(path), ksize=ksize)))`,
  check `.minhash.scaled == SCALED` else `DataAccessError`, return `.minhash.flatten()`.
- [ ] **Step 4: Implement `prepare/sra/sketch_containment.py`.**
  `genome_minhash(signature: dict)`: `load_one_signature_from_json(json.dumps(signature)).minhash`
  (the dict `sketch_fasta`/`load_signature` return; wrap in a list if the loader wants one).
  `containment_rows(genome_mh, paths: Dict[str, Path]) -> List[Tuple[str, float]]`: per
  accession `genome_mh.contained_by(load_dataset_minhash(path))`; a `DataAccessError`/`OSError`
  on one file is logged and the row left out. `merge_containment_csv(path, rows, ksize=KSIZE)`:
  read existing rows with `csv.DictReader` (utf-8-sig) when the file exists; header = existing
  fieldnames extended with any missing `BRANCHWATER_COLUMNS` (or exactly `BRANCHWATER_COLUMNS`
  for a new file); drop existing rows whose `acc` is in `rows`; add new rows with
  `containment=f"{c:.4f}"`, `cANI=round(c ** (1 / ksize), 4) if c > 0 else 0.0`, other columns
  empty; sort descending by float containment; `with open_atomic(path, "w", newline="") as fp:
  csv.DictWriter(...)`. Check `open_atomic`'s signature in `data/file_io.py` for the text mode
  and `newline` argument first.
- [ ] **Step 5: Run** the tests; PASS. `make check`.
- [ ] **Step 6: Commit** `feat: dataset sketcher, signature writer and local containment`.

### Task 3: store cache

**Files:** `prepare/store/sidecar.py`, `prepare/store/sketches.py`, `tests/test_store_sidecar.py`,
`tests/test_sra_sketch.py`.

- [ ] **Step 1: Failing tests.** In `tests/test_store_sidecar.py`: `sketch` is absent from
  `to_dict()` while None and round-trips when set (copy the `refetch` tests). In
  `tests/test_sra_sketch.py` with a tmp_path store (`store_paths`, `write_sidecar`, the
  monkeypatched environment from `tests/test_store_stats.py`): `test_cached_sketch_requires_matching_signature`
  (no block, changed file, missing `.sig.gz` all give None; a match gives the path and block);
  `test_store_sketch_copies_file_and_writes_block_under_lock` (the file appears at
  `<store>/sra/<ACC>/<ACC>.sig.gz`, the sidecar gains `sketch` with `file`, `computed` and the
  summary fields, nothing else in the sidecar changes); `test_store_sketch_warns_when_lock_held`
  (hold `dataset_lock` in the test with `STATS_LOCK_WAIT_SECONDS` patched to 0.1; a warning,
  no exception, no file); `test_store_sketch_is_a_noop_without_sidecar`.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** `Sidecar.sketch: Optional[Dict[str, Any]] = None` after `refetch`;
  `to_dict` pops it while None. `prepare/store/sketches.py`:
  ```python
  def store_signature_path(sidecar_path: Path) -> Path      # <acc_dir>/<ACC>.sig.gz
  def cached_sketch(acc_dir, sidecar_path, files) -> Optional[Tuple[Path, Dict[str, Any]]]
  def store_sketch(sidecar_path, signature_file: Path, summary: Dict[str, Any]) -> None
  ```
  `cached_sketch`: `read_sidecar`, `sidecar.sketch`, `block["files"] == _file_signature(files)`,
  ksizes and scaled match, the file exists. `store_sketch`: None sidecar path is a no-op;
  `_store_location` (import from `store.stats`) decides whether `dataset_lock(paths, acc,
  wait_seconds=STATS_LOCK_WAIT_SECONDS)` is taken; inside, `copy_file(signature_file, target)`,
  re-read the sidecar, set `sketch = {**summary, "file": target.name, "computed": now}`,
  `write_sidecar`; `DataAccessError` from the lock logs
  `"Not caching sketch for %s: %s"` and returns.
- [ ] **Step 4: Run** `pytest tests/test_store_sidecar.py tests/test_sra_sketch.py tests/test_store_stats.py -q`; PASS. `make check`.
- [ ] **Step 5: Commit** `feat: cache dataset sketches beside the store copy`.

### Task 4: the command

**Files:** `prepare/cli/commands/sra_sketch.py`, `prepare/cli/main.py`, `tests/test_cli_sra_sketch.py`,
`tests/test_optional_imports.py`, `tests/test_cli_run_log_wiring.py`.

- [ ] **Step 1: Failing tests** (fixtures copied from `tests/test_cli_sra_profile.py`: `_write_fastq`,
  `_args`, `_paired_dataset`, `_store_linked`, `_project_with_store`):
  `test_sketch_writes_signature_and_records` (exit 0; `sketches/SRR1.sig.gz`; registry
  `analyses.sketch` with `output == "sketches/SRR1.sig.gz"` and the summary key set
  `{ksizes, scaled, track_abundance, reads, bases, hashes, files, source}`);
  `test_run_log_summary_and_rows` through `main([...])` as `tests/test_cli_run_log_wiring.py` does;
  `test_second_run_skips_and_force_resketches` (mtime unchanged after the second run, summary
  `skipped == 1`; `os.utime` on a file or `--force` re-sketches);
  `test_dry_run_writes_nothing` (no folder, registry untouched, `records_run` False);
  `test_corrupt_gzip_fails_that_accession_only` (exit 1, no file for it, the other written);
  `test_interrupt_records_finished_and_exits_130` (patch `sketch_dataset` to set
  `args._termination.stop` after the first accession and raise `KeyboardInterrupt` on the
  second; `BaseCommand.run` returns 130; first recorded; no `.sig.gz` and no temp file for the second);
  `test_containment_csv_written_and_parsed` (`--genome-fasta`; `matches/<stem>.csv` equals
  `contained_by`; `parse_containment_data` ranks the better dataset first);
  `test_signature_input` (`--signature` with a plain `.sig` from `sketch_fasta`);
  `test_store_cache_round_trip` (first project writes the store file and block; a second project
  linking the same dataset reuses it with `sketch_dataset` patched to raise; `source == "store"`;
  usage row `analysed/sketch` in the catalogue);
  `test_missing_sourmash_is_exit_3` (`monkeypatch.setitem(sys.modules, "sourmash", None)`;
  `execute` returns 3; log names `prepare[sourmash]`; nothing written);
  `test_parser_defaults_and_exclusive_genome_flags`.
  `tests/test_optional_imports.py`: add `sra_sketch --accession SRR1` to the parse list.
- [ ] **Step 2: Run** `pytest tests/test_cli_sra_sketch.py tests/test_optional_imports.py -q`; FAIL.
- [ ] **Step 3: Implement `SRASketchCommand(BaseCommand)`.** `name = "sra_sketch"`,
  `group = "Reads"`, help "Sourmash sketches (k=21,31,51, scaled 1000) of downloaded datasets,
  with local containment of a genome". Flags as in the spec; `--max-workers` type `positive_int`
  (reuse the helper in `sra_profile.py` or `cli/base.py`), default None; `--genome-fasta` and
  `--signature` in a mutually exclusive group. `records_run`: `not dry_run`.
  Structure (copy `sra_profile`): `_resolve_accessions`; `_plan(args, accessions, registry)` ->
  list of `Planned(accession, files, acc_dir, signature, action)` with action `skip`, `reuse`
  or `sketch` (a missing folder or no files is `failed` with a warning); `_dry_run` emits
  `ACC: sketch (3 files)` lines and totals; `_sketch_all(args, planned)` with a
  `ThreadPoolExecutor(max_workers)`, `as_completed`, a `ProgressReporter`, per-accession lines
  at `progress.item_level`, `DATASET_READ_ERRORS` and `DataAccessError` failing one accession,
  and the interrupt branch (set `term.stop`, `executor.shutdown(wait=False, cancel_futures=True)`,
  remember the interrupt, continue to record); `_sketch_one` = `sketch_dataset` +
  `write_signature` + `store_sketch(_resolved_sidecar_path(acc_dir), path, summary)`;
  `_reuse_one` = `copy_file(store_path, project_path)`; `_containment(args, genome, done)`
  only with a genome, over sketched, reused and skipped accessions whose file exists;
  `_record` and `_note_run` as in the spec; `_run` ends with
  `if interrupted: raise KeyboardInterrupt("sra_sketch stopped")`, `1` if any failed, else `0`.
  `execute`: anchor `--fastq-folder`, `--output-dir` and `--matches-folder` with `anchored_folder`;
  unless dry run, `require("sourmash", "sourmash", "Sketching reads")` first; then the genome
  (`sketch_fasta` or `load_signature`); the same `except (PrepareError, OSError) -> self.fail`
  and broad-except pattern as `sra_profile`. Worker count:
  `args.max_workers or max(1, min(resources.available_cpus(), settings.active().max_workers_cap))`.
  Register in `cli/main.py` after `SRAReportCommand()`.
- [ ] **Step 4: Run** `pytest tests/test_cli_sra_sketch.py tests/test_optional_imports.py tests/test_cli_run_log_wiring.py tests/test_plugins.py -q`; PASS. `make check` (module size, radon, atomic writes, docs check will fail until Task 5).
- [ ] **Step 5: Commit** `feat: sra_sketch command`.

### Task 5: docs and release

**Files:** `README.md`, `docs/reads.md`, `docs/pipeline_overview.md`, `docs/ARCHITECTURE.md`,
`docs/shared_store.md`, `docs/configuration.md`, `CHANGELOG.md`, `pyproject.toml`, `prepare/__init__.py`.

- [ ] README: Reads group row for `sra_sketch` (`python scripts/check_docs_commands.py`); the
  `sourmash` extra row names `sra_sketch`.
- [ ] `docs/reads.md`: a `sra_sketch` section after `sra_profile`: what is sketched, the file
  format, skip and `--force`, the store cache, `--genome-fasta`/`--signature` and the merged
  matches CSV, `--max-workers` and memory, exit codes, the limits from the spec.
- [ ] `docs/pipeline_overview.md`: the analyse stage mentions `sra_sketch` and that its matches
  CSV feeds `parse_containment` and `select_datasets`.
- [ ] `docs/ARCHITECTURE.md`: `sra/` and `store/` module lists, the sidecar `sketch` field;
  `docs/shared_store.md`: a dataset folder may hold `<ACC>.sig.gz`; `docs/configuration.md`:
  `max_workers_cap` also bounds `sra_sketch`.
- [ ] CHANGELOG `## [Unreleased]` -> 0.11.0 Added; version bump; `make check`, `make test`.
- [ ] Live check on the Wolbachia project: sketch the downloaded runs with the wMel genome,
  `sourmash sig describe` on one file, `parse_containment` on the matches folder ranks
  SRR11011981 first, `status --export-tsv` lists `sketch`, `runs` shows the record, Ctrl-C
  mid-run leaves no `.sig.gz` for the interrupted accession.
- [ ] `superpowers:finishing-a-development-branch`.
