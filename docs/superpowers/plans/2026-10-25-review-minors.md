# Review Minors (0.12.1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the seven Minor findings the 0.12.0 whole-branch review deferred, as patch release 0.12.1.

**Architecture:** Seven independent fixes in the data, processing and CLI layers, each one a failing test
first and the smallest change after. No new module; one new optional out-parameter on
`extract_target_reads` (the `failures` pattern of `download_metadata`) so the command, not the
library loop, decides what a recording failure does to the run.

**Tech Stack:** Python 3.12, pytest, the existing fake tools in `tests/helpers_extraction.py`.

**Spec:** No spec document. The authority is the deferred-minors list of the 2026-10-21 review of
`fix/pipeline-gui-mediums`, repeated under Context below, read against the rules in `CLAUDE.md`
(exit-code policy, `self.fail`, atomic writes, budgets that only shrink).

## Context

The review deferred, in this order:

1. A non-transient `PrepareError` from the extraction recording callback is a warning and the run exits 0.
2. A deferred lock timeout makes `extract_target_reads` raise after the loop, which skips `--assemble`
   and the run-log summary; and when the callback raises, the stop check after it is skipped, so the
   next sample starts before the interrupt is honoured.
3. `holm_adjusted` sorts a NaN p-value first and reports it as 0.0 (significant).
4. `genome_name(".fna")` is `""`.
5. A truncated or corrupt `.fna.gz` raises `EOFError`/`gzip.BadGzipFile` uncaught through `sketch_fasta`
   (`branchwater_search` and `sra_sketch`).
6. `download_metadata` with some accessions failing for a network reason exits 1 (not 4), and the error
   line names every failed accession.
7. `assembly_coverage` leaves the singleton file out of the back-mapping, so the mapping rate omits
   those reads.

## Global Constraints

- Branch `fix/review-minors` from main (v0.12.0, 8aa37f9); release 0.12.1; one commit per task.
- Every fix starts with a failing test; `conda run -n metaquest make check` and `make test` green before each commit.
- Exit codes: `PrepareError` 1, `ValidationError` 2, `ConfigurationError` 3, `TransientError` 4, interrupt 130; errors leave a command through `self.fail`.
- `prepare/data/read_extraction.py` stays at or under 700 lines, `prepare/sra/analytics.py` under 800; `MAX_RANK_C` (104) and the B902 count (25) only shrink.
- ASCII only in `prepare/`, `tests/`, `docs/` is free text; modest scientific language.
- Commit trailer `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`.

## Review Focus

1. A recording failure on the last sample: the summary, the assembly and the exit code must all still happen (Task 1 tests the first sample; `test_a_recording_failure_on_the_last_sample_still_assembles` is added there).
2. Both a lock timeout and a plain `PrepareError` in one run: exit 4 wins (Task 1, `test_a_lock_timeout_outranks_a_recording_error`).
3. All p-values NaN: `holm_adjusted` returns all NaN and `_run_statistical_tests` reports nothing significant (Task 2, second assertion).
4. A gzipped FASTA whose first two bytes are the gzip magic but whose body is plain text: `DataAccessError`, not a traceback (Task 4, `test_sketch_fasta_rejects_a_fake_gzip_header`).
5. Every accession failing for a network reason with `failures` given: still `NetworkError` from the data layer, exit 4 once (Task 5 keeps `test_download_metadata_exits_4_when_ncbi_cannot_be_reached`).

---

### Task 1: Recording failures reach the exit code; the stop check always runs

**Files:**
- Modify: `prepare/data/read_extraction.py` (`_notify_result`, `extract_target_reads`)
- Modify: `prepare/cli/commands/read_extraction.py` (`_record_result_and_check_stop`, `execute`)
- Test: `tests/test_read_extraction.py`, `tests/test_cli_read_extraction.py`

**Interfaces:**
- Produces: `extract_target_reads(..., record_failures: Optional[List[PrepareError]] = None)`. When a list is
  given, every `PrepareError` the callback raises is appended to it (a warning is logged, the run goes on)
  and nothing is raised after the loop. When it is `None`, the behaviour of 0.12.0 stands: a
  `TransientError` is raised after the loop, any other `PrepareError` is only logged.
- Produces: `_notify_result(on_result, accession, result, deferred: List[PrepareError]) -> None` appends every `PrepareError`.

- [ ] **Step 1: Write the failing library test** in `tests/test_read_extraction.py`, class holding `test_a_lock_timeout_while_recording_is_a_transient_error`:

```python
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_record_failures_collects_every_recording_error_and_raises_nothing(self, mock_run):
        from prepare.core.exceptions import DataAccessError, LockTimeoutError

        mock_run.side_effect = _fake_tools({})
        errors = {"SRR1": LockTimeoutError("registry lock held"), "SRR2": DataAccessError("disk full")}

        def on_result(accession, result):
            if accession in errors:
                raise errors[accession]

        failures = []
        with tempfile.TemporaryDirectory() as tmp:
            root, table, genome = _make_tree(tmp, paired=True)
            results = extract_target_reads(
                parsed_containment=table, genome_id="GCF_1", genome_fasta=genome,
                fastq_folder=root / "fastq", output_folder=root / "targeted", threshold=0.3,
                on_result=on_result, record_failures=failures,
            )
        assert set(results) == {"SRR1", "SRR2"}
        assert [type(e).__name__ for e in failures] == ["LockTimeoutError", "DataAccessError"]
```

- [ ] **Step 2: Run it**

Run: `conda run -n metaquest pytest tests/test_read_extraction.py -k record_failures_collects -q`
Expected: FAIL with `TypeError: ... unexpected keyword argument 'record_failures'`

- [ ] **Step 3: Implement** in `prepare/data/read_extraction.py`: `_notify_result` catches `PrepareError` once, logs the warning and appends to `deferred`; `extract_target_reads` gains `record_failures: Optional[List[PrepareError]] = None`, passes `record_failures if record_failures is not None else (deferred := [])` as the list, and after the loop, only when `record_failures is None`, raises the first `TransientError` in `deferred` if any. Update the docstring of `_notify_result` and the parameter list of `extract_target_reads`.

- [ ] **Step 4: Run the file**

Run: `conda run -n metaquest pytest tests/test_read_extraction.py -q`
Expected: all pass, including `test_a_lock_timeout_while_recording_is_a_transient_error` and `test_an_unexpected_callback_error_propagates`

- [ ] **Step 5: Write the failing CLI tests** in `tests/test_cli_read_extraction.py`, class `TestExtractionFollowsTheProject` (use `_project_tree`, `_args`, `_fake_tools` as the class's `test_a_real_run_skips_the_excluded_accession` does; patch `ReadExtractionCommand._record_result` with a side effect that raises for one accession and otherwise calls the real method):

```python
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_a_recording_error_is_exit_1_after_the_run(self, mock_run, tmp_path, caplog):
        """Minor 1: a DataAccessError from the registry write is reported and fails the run."""

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_a_lock_timeout_still_assembles_and_notes_the_run_then_exits_4(self, mock_run, tmp_path):
        """Minor 2: with --assemble, the lock timeout on SRR1 does not skip megahit for SRR2;
        the run-log summary is written; the exit is 4."""

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_a_recording_failure_on_the_last_sample_still_assembles(self, mock_run, tmp_path):
        """Review focus 1: the error on the last sample changes the exit code, nothing else."""

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_a_lock_timeout_outranks_a_recording_error(self, mock_run, tmp_path):
        """Review focus 2: one TransientError among the failures makes the exit 4."""

    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_a_stop_during_a_failing_record_is_honoured_before_the_next_sample(self, mock_run, tmp_path):
        """Minor 2: the stop flag set while recording SRR1 (which then raises) stops the run
        before SRR2 is started; exit 130 through run()."""
```

Assertions: exit codes as named; for the assemble test, `megahit` appears in `state["calls"]` for SRR2 and
`run_log` summary is noted (patch `prepare.cli.commands.read_extraction.run_log.note_run` and assert a call
with `with_reads` in its summary); for the stop test, `SRR2` has no `targeted/SRR2` folder and the
callback saw only `SRR1`. The error line: `"Recording failed for %d sample(s): %s"` at ERROR with the
accessions and the first line of each reason.

- [ ] **Step 6: Run them**

Run: `conda run -n metaquest pytest tests/test_cli_read_extraction.py -k "recording_error or lock_timeout or failing_record or recording_failure" -q`
Expected: 5 FAIL (exit 0 where 1/4 expected; assembly skipped; stop one sample late)

- [ ] **Step 7: Implement** in `prepare/cli/commands/read_extraction.py`:
  - `_record_result_and_check_stop`: wrap the `_record_result` call and the progress/clock lines in
    `try: ... finally:` so the stop check at the end runs whether or not the record raised; a set stop
    raises `KeyboardInterrupt` from the `finally` (it replaces the record error, 130 outranks it).
  - `execute`: `record_failures: List[PrepareError] = []` passed to `extract_target_reads`; after the
    assembly block (and after `_report_no_reads`'s own `return 1` path), a new helper
    `_fail_for_recording(self, failures: List[PrepareError]) -> int` logs the ERROR line and returns
    `self.fail(first_transient_or_first, "Error recording extraction results")`. Call it wherever
    `execute` would return 0 or 1 with a non-empty list, so a recording failure never drops the exit
    code below what the run earned (an assembly failure's 1 stays 1 unless a transient failure makes it 4).

- [ ] **Step 8: Run the two files and the gates**

Run: `conda run -n metaquest pytest tests/test_read_extraction.py tests/test_cli_read_extraction.py -q && conda run -n metaquest make check`
Expected: all pass; `make check` green; `read_extraction.py` at or under 700 lines (trim docstrings if not)

- [ ] **Step 9: Commit**

```bash
git add prepare/data/read_extraction.py prepare/cli/commands/read_extraction.py tests/test_read_extraction.py tests/test_cli_read_extraction.py
git commit -m "fix: a failed extraction record reaches the exit code after the summary and assembly; the stop check always runs"
```

---

### Task 2: `holm_adjusted` leaves NaN out of the correction

**Files:**
- Modify: `prepare/processing/statistics.py` (`holm_adjusted`)
- Test: `tests/test_sra_analytics.py` (next to `test_holm_adjusted_helper`)

**Interfaces:**
- Produces: `holm_adjusted(p_values: Sequence[float]) -> List[float]`: a NaN input gives NaN at its position and is not counted among the tests; the finite values are corrected as before.

- [ ] **Step 1: Write the failing test**

```python
    def test_holm_adjusted_leaves_nan_out_of_the_correction(self):
        import math

        from prepare.processing.statistics import holm_adjusted

        adjusted = holm_adjusted([float("nan"), 0.01, 0.03])
        assert math.isnan(adjusted[0])
        assert [round(p, 6) for p in adjusted[1:]] == [0.02, 0.03]  # m = 2, not 3
        assert all(math.isnan(p) for p in holm_adjusted([float("nan"), float("nan")]))
```

- [ ] **Step 2: Run it**

Run: `conda run -n metaquest pytest tests/test_sra_analytics.py -k nan_out -q`
Expected: FAIL (`adjusted[0] == 0.0`, and the finite values corrected with m = 3)

- [ ] **Step 3: Implement**: build `order` over the indices whose value is finite (`math.isfinite`), use their count as m, preset `adjusted` to `float("nan")` for every index.

- [ ] **Step 4: Run the file**

Run: `conda run -n metaquest pytest tests/test_sra_analytics.py -q`
Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add prepare/processing/statistics.py tests/test_sra_analytics.py
git commit -m "fix: holm_adjusted leaves a NaN p-value out of the correction instead of reporting it as 0.0"
```

---

### Task 3: `genome_name` of a dot-file keeps the name

**Files:**
- Modify: `prepare/core/validation.py` (`genome_name`)
- Test: `tests/test_core_validation.py` (class `TestGenomeName`)

- [ ] **Step 1: Write the failing test**

```python
    @pytest.mark.parametrize("name", [".fna", ".fna.gz", ".sig"])
    def test_a_file_that_is_only_a_suffix_keeps_its_name(self, name):
        from prepare.core.validation import genome_name

        assert genome_name(name) == name
```

- [ ] **Step 2: Run it**

Run: `conda run -n metaquest pytest tests/test_core_validation.py -k only_a_suffix -q`
Expected: 3 FAIL (`'' == '.fna'`)

- [ ] **Step 3: Implement**: strip a suffix only when something remains (`len(name) > len(suffix)`); otherwise return `name`. One line in the docstring.

- [ ] **Step 4: Run the file**

Run: `conda run -n metaquest pytest tests/test_core_validation.py -q`
Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add prepare/core/validation.py tests/test_core_validation.py
git commit -m "fix: genome_name keeps the name of a file that is only a suffix"
```

---

### Task 4: A corrupt gzipped FASTA is a `DataAccessError`

**Files:**
- Modify: `prepare/data/branchwater_search.py` (`sketch_fasta`)
- Test: `tests/test_branchwater_search.py` (class holding `test_sketch_fasta_reads_gzip`), `tests/test_cli_sra_sketch.py`

**Interfaces:**
- Produces: `sketch_fasta` raises `DataAccessError("Cannot read <path>: <reason>")` for `OSError` (which includes `gzip.BadGzipFile`), `EOFError` and `zlib.error` while reading; the `ConfigurationError` for a missing sourmash and the two existing `DataAccessError`s are unchanged.

- [ ] **Step 1: Write the failing library tests**

```python
    def test_sketch_fasta_reports_a_truncated_gzip(self, tmp_path):
        """A download cut short: the gzip stream ends without its trailer."""
        import gzip

        from prepare.core.exceptions import DataAccessError

        whole = gzip.compress(b">g\n" + b"ACGT" * 500 + b"\n")
        path = tmp_path / "g.fna.gz"
        path.write_bytes(whole[: len(whole) // 2])
        with pytest.raises(DataAccessError, match="Cannot read"):
            sketch_fasta(path)

    def test_sketch_fasta_rejects_a_fake_gzip_header(self, tmp_path):
        """Review focus 4: the gzip magic followed by plain text."""
        from prepare.core.exceptions import DataAccessError

        path = tmp_path / "g.fna.gz"
        path.write_bytes(b"\x1f\x8b>g\nACGT\n")
        with pytest.raises(DataAccessError, match="Cannot read"):
            sketch_fasta(path)
```

- [ ] **Step 2: Run them**

Run: `conda run -n metaquest pytest tests/test_branchwater_search.py -k "truncated_gzip or fake_gzip" -q`
Expected: 2 FAIL with `EOFError` / `gzip.BadGzipFile` (skipped if sourmash is absent; it is present in `metaquest`)

- [ ] **Step 3: Implement** in `sketch_fasta`: wrap the `with _open_text(path) as handle:` block in `try/except (OSError, EOFError, zlib.error) as e: raise DataAccessError(f"Cannot read {path}: {e}") from e` (`zlib` is already imported).

- [ ] **Step 4: Write the failing CLI test** in `tests/test_cli_sra_sketch.py` next to `test_genome_id_of_a_gz_fasta_drops_both_suffixes`, using that test's fixture shape with the genome file replaced by the truncated gzip of Step 1:

```python
def test_a_corrupt_genome_gzip_is_exit_1_not_a_traceback(tmp_path, caplog):
    ...
    assert rc == 1
    assert "Cannot read" in caplog.text
```

- [ ] **Step 5: Run the two files**

Run: `conda run -n metaquest pytest tests/test_branchwater_search.py tests/test_cli_sra_sketch.py tests/test_cli_branchwater_search.py -q`
Expected: all pass

- [ ] **Step 6: Commit**

```bash
git add prepare/data/branchwater_search.py tests/test_branchwater_search.py tests/test_cli_sra_sketch.py
git commit -m "fix: a truncated or corrupt gzipped genome FASTA is a DataAccessError in branchwater_search and sra_sketch"
```

---

### Task 5: A partial network failure in `download_metadata` exits 4; the failure line is short

**Files:**
- Modify: `prepare/cli/commands/metadata.py` (`DownloadMetadataCommand.execute`, new `_fail_for_missing`)
- Test: `tests/test_cli_metadata.py`

**Interfaces:**
- Produces: `DownloadMetadataCommand._fail_for_missing(self, failures: Dict[str, str]) -> int`: logs
  `"No metadata for %d accession(s): %s"` naming at most 10 accessions (sorted) and `"... and N more"`;
  returns `self.fail(NetworkError(...), "Error downloading metadata")` (exit 4) when any reason starts with
  `"network: "`, else `self.fail(DataAccessError(...), ...)` (exit 1). Called in place of the `if failures:` block.

- [ ] **Step 1: Write the failing tests** (next to `test_download_metadata_exit_1_when_an_accession_has_no_metadata`, same fixture style):

```python
def test_download_metadata_exits_4_when_some_accessions_failed_for_a_network_reason(tmp_path, caplog):
    """Minor 6: the fetched one is recorded, the missing one is named, and the exit says retry."""
    # batches = [({"SRR1000000": good}, {"SRR1000001": "network: connection refused"})]
    assert DownloadMetadataCommand().execute(_download_args(tmp_path, accessions_file=str(accessions))) == 4
    assert "SRR1000001" in caplog.text
    assert "SRR1000000" in json.loads((tmp_path / "prepare_registry.json").read_text())["datasets"]


def test_download_metadata_failure_line_names_at_most_ten_accessions(tmp_path, caplog):
    # 12 failures "not in the NCBI response"; exit 1
    line = next(r.getMessage() for r in caplog.records if r.levelno == logging.ERROR and "No metadata" in r.getMessage())
    assert line.count("SRR") == 10 and line.endswith("and 2 more")
```

- [ ] **Step 2: Run them**

Run: `conda run -n metaquest pytest tests/test_cli_metadata.py -k "network_reason or at_most_ten" -q`
Expected: 2 FAIL (exit 1 instead of 4; twelve accessions on the line)

- [ ] **Step 3: Implement** `_fail_for_missing` and call it; import `NetworkError`, `DataAccessError` from `prepare.core.exceptions` if not already.

- [ ] **Step 4: Run the file and the metadata data tests**

Run: `conda run -n metaquest pytest tests/test_cli_metadata.py tests/test_data_metadata.py -q`
Expected: all pass (`test_download_metadata_exits_4_when_ncbi_cannot_be_reached` and `..._exit_1_when_an_accession_has_no_metadata` unchanged)

- [ ] **Step 5: Commit**

```bash
git add prepare/cli/commands/metadata.py tests/test_cli_metadata.py
git commit -m "fix: download_metadata exits 4 when any accession failed for a network reason; the failure line names at most ten"
```

---

### Task 6: `assembly_coverage` maps the singleton file too

**Files:**
- Modify: `prepare/data/read_extraction_export.py` (`assembly_coverage`, new `_count_back_mapped`)
- Test: `tests/test_read_extraction.py` (class holding `test_assembly_coverage_mapping_rate_is_pinned`)

**Interfaces:**
- Produces: `_count_back_mapped(contigs, queries: List[Path], preset, threads, work_root: Path) -> int`: one
  minimap2 run over `queries` into `work_root / "coverage.sam"`, one samtools filter into
  `coverage.bam`, the count, both files removed in `finally`. `assembly_coverage` calls it once per
  query group: the mate files, then the singleton file when there is one, and sums the counts.

- [ ] **Step 1: Write the failing test**

```python
    @patch("prepare.data.read_extraction.SecureSubprocess.run_secure")
    def test_assembly_coverage_counts_the_singleton_file_in_a_second_run(self, mock_run, tmp_path):
        """Minor 7: mates in one minimap2 run, the singleton in another; the counts add up."""
        state = {"sam_records": self.RECORDS}
        mock_run.side_effect = _fake_tools(state)
        contigs = tmp_path / "final.contigs.fa"
        contigs.write_text(">c1 len=10\nACGT\n")
        reads = [tmp_path / f"GCF_1_{tag}.fastq.gz" for tag in ("1", "2", "s")]
        for path in reads:
            with gzip.open(path, "wt") as handle:
                handle.write("@r\nACGT\n+\nIIII\n")

        result = assembly_coverage(contigs, reads, "sr", 1, tmp_path, mapped_reads=4)

        aligns = [args for exe, args in state["calls"] if exe == "minimap2" and "-o" in args]
        assert [a[-2:] for a in aligns] == [[str(reads[0]), str(reads[1])], [str(reads[2])]]
        assert result["reads_mapped"] == 4  # 2 from each run of the fake
        assert result["mapping_rate"] == 1.0
        assert not (tmp_path / "coverage.sam").exists() and not (tmp_path / "coverage.bam").exists()
```

- [ ] **Step 2: Run it**

Run: `conda run -n metaquest pytest tests/test_read_extraction.py -k singleton_file_in_a_second_run -q`
Expected: FAIL (one minimap2 run; `reads_mapped == 2`)

- [ ] **Step 3: Implement**: extract the try/finally of `assembly_coverage` into `_count_back_mapped`; compute `groups = [mates, singles]` (drop an empty group; a `reads` of one or two files is one group) and `reads_mapped = sum(_count_back_mapped(...) for group in groups)`. Update the comment that explains why the singleton is mapped on its own.

- [ ] **Step 4: Run the file and the gates**

Run: `conda run -n metaquest pytest tests/test_read_extraction.py tests/test_extraction_staging.py tests/test_cli_read_extraction.py -q && conda run -n metaquest make check`
Expected: all pass; `make check` green

- [ ] **Step 5: Commit**

```bash
git add prepare/data/read_extraction_export.py tests/test_read_extraction.py
git commit -m "fix: assembly_coverage maps the singleton reads in a second minimap2 run so the mapping rate counts them"
```

---

### Task 7: Release 0.12.1

**Files:**
- Modify: `CHANGELOG.md`, `pyproject.toml`, `prepare/__init__.py`, `CLAUDE.md` (version line), `docs/packaging.md` (two version strings), `docs/extraction_assembly.md` (coverage sentence: singleton reads are mapped too), `docs/ARCHITECTURE.md` (the `download_metadata` sentence in the exit-code section: exit 4 when any accession failed for a network reason), `docs/superpowers/README.md` (row for this plan)

- [ ] **Step 1: CHANGELOG** section `## [0.12.1] - 2026-10-25` with one `### Fixed` list, one bullet per task in Task order, each naming the command or function and the behaviour change in one sentence.

- [ ] **Step 2: Version** `0.12.1` in the four files; grep confirms no `0.12.0` remains outside CHANGELOG and the README row for the 0.12.0 plan.

Run: `grep -rn "0\.12\.0" pyproject.toml prepare/__init__.py CLAUDE.md docs/packaging.md`
Expected: no output

- [ ] **Step 3: Docs** as listed; `python scripts/check_docs_commands.py` passes.

- [ ] **Step 4: Full gates**

Run: `conda run -n metaquest make check && conda run -n metaquest make test 2>&1 | tail -3`
Expected: `make check` green; `N passed` with N at least 3985 + the tests added (about 15), no failures

- [ ] **Step 5: Commit**

```bash
git add CHANGELOG.md pyproject.toml prepare/__init__.py CLAUDE.md docs/packaging.md docs/extraction_assembly.md docs/ARCHITECTURE.md docs/superpowers/README.md
git commit -m "release: 0.12.1"
```
