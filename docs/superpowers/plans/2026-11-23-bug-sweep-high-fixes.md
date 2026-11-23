# Bug sweep 2026-11-23: High-severity fixes (0.12.3) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the seven High findings of the 2026-11-23 six-agent bug sweep and release them as 0.12.3.

**Architecture:** Each fix is local to one layer: project selection (`data/project_selection.py`,
`data/registry.py`, `processing/selection.py`), download verification (`data/sra/fastq.py`,
`data/sra/accession.py`, `store/sidecar.py`), the shared store commands (`cli/commands/store/`,
`cli/commands/sra.py`), and the download/metadata data layer (`data/sra/download.py`,
`data/sra_metadata.py`, `utils/http.py`). Every task starts with a test that fails on current `main`.

**Tech Stack:** Python 3.12, pytest, pandas; conda environment `metaquest`.

**Spec:** The sweep findings (no separate spec document). Each task quotes the defect it fixes, the
reproduction the auditing agent ran, and the agreed behaviour.

## Global Constraints

- Run commands as `cd /Users/andreassjodin/Code/prepare && conda run -n metaquest <cmd>`; the shell cwd resets.
- Work on branch `fix/bug-sweep-high` off `main` (55efb3e). Stage files by name; never `git add -A`.
- A bug fix starts with a test that reproduces the bug and fails, then the smallest fix (CLAUDE.md).
- Errors through `self.fail`; raise the matching class from `core/exceptions.py` (`LockTimeoutError` -> exit 4,
  `DataAccessError` -> exit 1, `ConfigurationError` -> exit 3). No new exit codes.
- Output files through the atomic helpers in `data/file_io.py`; no `print(`; ASCII only in `prepare/`, `tests/`.
- Tests never touch the network, a real store or `~/.config/prepare`; store tests use `tmp_path` and
  monkeypatch `PREPARE_DATA` and `HOME`.
- Modules stay under 800 lines (`python tests/test_module_sizes.py --check`); no new radon rank D.
- Modest scientific language in docs and comments.
- `make check` and `make test` pass before each commit.

## Review Focus

1. A registry written by 0.12.2 that already holds selections with mixed criteria: `selection_for_genome`
   must give the same answer on every run (Task 1 test with mixed blocks covers it).
2. `select_datasets --skip-downloaded --top-n N`: a kept downloaded accession must not use up one of the N
   places (Task 2 test `test_skip_downloaded_kept_accessions_do_not_count_against_top_n`).
3. Paired data where `--split-3` also wrote an orphan `<acc>.fastq`: the orphan count must not be compared
   with mate 2 (Task 3 test `test_verify_download_orphan_file_does_not_count_against_mate_two`).
4. A single-end run (`<acc>.fastq` only): no mate-2 check, verdict unchanged (Task 3 test
   `test_verify_download_single_end_is_unchanged`).
5. A blacklist file with a UTF-8 BOM saved by Excel: it must still be read, not refused (Task 6 test
   `test_blacklist_with_bom_is_read`).

---

### Task 1: `selection_for_genome` decides per accession

Defect: `prepare/data/project_selection.py:30-35` takes the criteria of whichever selected accession the
set iteration yields first. With mixed criteria (left by `--skip-downloaded` keeps), the answer for a genome
depends on `PYTHONHASHSEED`; two agents reproduced it.

Behaviour: each selected accession is judged by its own selection block. An accession is in the result for
`genome_id` when its block's criteria are empty, or `column` is None, `genome_id` or `"max_containment"`,
or `genome_id` is in `criteria["genome_ids"]`. Accessions left out are counted in one INFO line.

**Files:**
- Modify: `prepare/data/project_selection.py:19-41`
- Test: `tests/test_project_selection.py` (create; no test module for this file exists)

**Interfaces:**
- Produces: `selection_for_genome(registry: Registry, genome_id: str) -> Set[str]` (unchanged signature;
  now a subset of the selected accessions instead of all-or-nothing), and a private
  `_covers(criteria: Dict[str, Any], genome_id: str) -> bool`.

- [ ] **Step 1: Write the failing tests**

```python
from prepare.data import registry_blocks as rb
from prepare.data.project_selection import selection_for_genome
from prepare.data.registry import Registry


def _select(reg, acc, criteria):
    rb.set_selection_block(reg, acc, rb.SelectionBlock(selected=True, date="2026-11-23", criteria=criteria, output="a.txt"))


def test_mixed_criteria_are_judged_per_accession(tmp_path):
    reg = Registry(path=tmp_path / "prepare_registry.json")
    _select(reg, "SRR1", {"column": "G1", "threshold": 0.5})
    _select(reg, "SRR3", {"column": "G2", "threshold": 0.5})
    _select(reg, "SRR4", {"column": "G2", "threshold": 0.5})

    assert selection_for_genome(reg, "G1") == {"SRR1"}
    assert selection_for_genome(reg, "G2") == {"SRR3", "SRR4"}


def test_max_containment_and_genome_ids_cover_a_genome(tmp_path):
    reg = Registry(path=tmp_path / "prepare_registry.json")
    _select(reg, "SRR1", {"column": "max_containment"})
    _select(reg, "SRR2", {"column": "A+B", "genome_ids": ["A", "B"]})
    _select(reg, "SRR3", {})

    assert selection_for_genome(reg, "B") == {"SRR1", "SRR2", "SRR3"}
    assert selection_for_genome(reg, "C") == {"SRR1", "SRR3"}


def test_no_selection_is_an_empty_set(tmp_path):
    assert selection_for_genome(Registry(path=tmp_path / "r.json"), "G1") == set()
```

- [ ] **Step 2: Run and confirm failure**

Run: `conda run -n metaquest pytest tests/test_project_selection.py -v`
Expected: the first two tests FAIL (current code returns all selected accessions or none).

- [ ] **Step 3: Implement** the per-accession filter in `selection_for_genome`. Update its docstring and the
  module docstring (the selection can now cover part of the accessions). Keep the lazy `query` import.

- [ ] **Step 4: Run** `conda run -n metaquest pytest tests/test_project_selection.py tests/test_read_extraction.py tests/test_cli_read_extraction.py tests/test_cli_status.py -v`
Expected: PASS. If an existing test asserted the all-or-nothing behaviour, update it to the per-accession
rule and say so in the commit message.

- [ ] **Step 5: Commit**

```bash
git add prepare/data/project_selection.py tests/test_project_selection.py
git commit -m "fix: judge each selected accession by its own criteria in selection_for_genome"
```

---

### Task 2: `--skip-downloaded` keeps downloaded accessions that meet the criteria selected

Defect: `record_selection(keep=...)` (`prepare/data/registry.py:218-248`) only avoids unselecting accessions
that were already selected; a downloaded accession that meets the new criteria but was never selected stays
unselected, so extraction (which narrows to the selection) drops it. Reproduced by the pipeline agent
(SRR2 at containment 0.8, threshold 0.6, missing from the extraction set). Also, `keep` is computed from
accessions dropped at the threshold stage, before the metadata, run and `--top-n` filters.

Behaviour:
- In `select_accessions_ranked` (`prepare/processing/selection.py:~360`) move `_drop_excluded` after the
  metadata and run filters and before `--top-n`. `dropped` then holds accessions that pass every criterion
  except the exclusion, and `--top-n` counts only new accessions.
- In `record_selection`, every accession in `keep` gets a fresh selection block with this run's criteria
  (`selected=True`, no `ranked` entry), whether or not it was selected before. All selected accessions then
  carry the same criteria after any `select_datasets` run.
- The output accessions file is unchanged: it lists only the new accessions (the ones to download).

**Files:**
- Modify: `prepare/processing/selection.py` (`select_accessions_ranked`, order of `_drop_excluded`)
- Modify: `prepare/data/registry.py:218-248` (`record_selection`)
- Modify: `prepare/cli/commands/select.py:293` (comment only: kept accessions are re-selected with these criteria)
- Modify: `docs/reads.md:27` (one sentence: downloaded accessions that meet the criteria stay selected)
- Test: `tests/test_cli_select.py`, `tests/test_data_registry.py`

**Interfaces:**
- Consumes: Task 1's per-accession `selection_for_genome`.
- Produces: `record_selection(registry, accessions, criteria, output, ranked=None, keep=())` (signature
  unchanged; `keep` accessions are now always selected with `criteria`).

- [ ] **Step 1: Write the failing tests**

In `tests/test_data_registry.py`:

```python
def test_record_selection_selects_kept_accessions_with_the_new_criteria(tmp_path):
    r = reg.Registry(path=tmp_path / "prepare_registry.json")
    reg.record_download(r, "SRR2", "downloaded", tmp_path / "fastq")  # never selected
    reg.record_selection(r, ["SRR1", "SRR3"], {"column": "G", "threshold": 0.6}, tmp_path / "a.txt", keep=["SRR2"])

    assert set(reg.query(r, "selected")) == {"SRR1", "SRR2", "SRR3"}
    assert rb.selection_block(r, "SRR2").criteria["column"] == "G"
```

(Use the module's existing import names and the registry's real download-recording helper; check the top of
the test module.)

In `tests/test_cli_select.py`, using the module's existing helpers for a containment table and registry:
- `test_skip_downloaded_keeps_a_never_selected_download_selected`: table G = {SRR1: 0.9, SRR2: 0.8,
  SRR3: 0.7}; SRR2 recorded downloaded, never selected; run `--genome-id G --threshold 0.6 --skip-downloaded`.
  Assert the output file lists `SRR1`, `SRR3` only, and `query(registry, "selected")` is
  `{"SRR1", "SRR2", "SRR3"}`.
- `test_skip_downloaded_kept_accessions_do_not_count_against_top_n`: same table plus SRR4: 0.65, SRR1
  downloaded; `--top-n 2 --skip-downloaded`. Assert the output file is `SRR2`, `SRR3`, and SRR1 is selected.
- `test_skip_downloaded_does_not_keep_a_download_below_the_threshold`: SRR5: 0.2 downloaded; threshold 0.6.
  Assert SRR5 is not selected.

- [ ] **Step 2: Run and confirm failure**

Run: `conda run -n metaquest pytest tests/test_data_registry.py tests/test_cli_select.py -k "kept or skip_downloaded" -v`
Expected: the first three FAIL; the below-threshold test may already pass.

- [ ] **Step 3: Implement** both changes. Keep the "N already downloaded, M excluded, K selected" log line;
  its counts now refer to accessions that passed every filter.

- [ ] **Step 4: Run** `conda run -n metaquest pytest tests/test_cli_select.py tests/test_data_registry.py tests/test_processing_selection.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add prepare/processing/selection.py prepare/data/registry.py prepare/cli/commands/select.py docs/reads.md tests/test_cli_select.py tests/test_data_registry.py
git commit -m "fix: select_datasets --skip-downloaded keeps downloads that meet the criteria selected"
```

---

### Task 3: A mate-2 file shorter than mate 1 is a truncated download

Defect: `verify_download` (`prepare/data/sra/fastq.py:334-370`) counts only the mate-1 (or single-end) file
plus the orphan file. A `_2` file with 40 reads beside a `_1` file with 100 was reported
`complete (100 of 100 spots)`; the store sidecar (`prepare/store/sidecar.py:165-186`) records `state=complete`.

Behaviour:
- When a mate-1 file exists and a mate-2 file (`_2`/`_R2`) exists, count mate 2 and compare it with the
  mate-1 count alone (not including the orphan file). If they differ, the verdict is `"truncated"`
  whatever `expected_spots` is (also when it is None). A mate-1 file with no mate-2 file is left as today
  (no mate check), since some runs yield only `_1`; this is noted in the docstring.
- The returned dict gains `"reads_r2": Optional[int]` (None for single-end) and `"mates_match": bool`
  (True for single-end).
- `_verdict_message` (`accession.py:~297`) writes `Downloaded N files, truncated (mate files differ: A and B reads)`
  when the mates differ; `parse_verdict_message` (`fastq.py:~373`) reads it back as
  `{"verdict": "truncated", "reads_r1": A, "reads_r2": B, "expected_spots": None, "ratio": None}`.
- Callers that already hold counts pass them: `_plain_verdict` (`accession.py:~275`) and `build_sidecar`
  (`sidecar.py:173`) pass `reads_r2=` from their per-file counts.

**Files:**
- Modify: `prepare/data/sra/fastq.py` (new `mate2_fastq`, `verify_download`, `parse_verdict_message`)
- Modify: `prepare/data/sra/accession.py` (`_plain_verdict`, `_verdict_message`)
- Modify: `prepare/store/sidecar.py:165-186`
- Test: `tests/test_data_sra_download.py`, `tests/test_store_sidecar.py`

**Interfaces:**
- Produces: `mate2_fastq(acc_dir: Union[str, Path]) -> Optional[Path]` (the `_2`/`_R2` file, None when absent);
  `verify_download(accession, acc_dir, expected_spots, reads_r1=None, reads_orphan=None, reads_r2=None) -> Dict[str, Any]`
  with the new keys `reads_r2` and `mates_match`.

- [ ] **Step 1: Write the failing tests** (in `tests/test_data_sra_download.py`, with a small helper that
  writes `n` four-line FASTQ records)

```python
def test_verify_download_short_mate_two_is_truncated(tmp_path):
    acc = tmp_path / "SRR1"; acc.mkdir()
    _fastq(acc / "SRR1_1.fastq", 100); _fastq(acc / "SRR1_2.fastq", 40)
    result = verify_download("SRR1", acc, 100)
    assert result["verdict"] == "truncated"
    assert (result["reads_r2"], result["mates_match"]) == (40, False)


def test_verify_download_short_mate_two_is_truncated_without_expected_spots(tmp_path):
    acc = tmp_path / "SRR1"; acc.mkdir()
    _fastq(acc / "SRR1_1.fastq", 100); _fastq(acc / "SRR1_2.fastq", 40)
    assert verify_download("SRR1", acc, None)["verdict"] == "truncated"


def test_verify_download_orphan_file_does_not_count_against_mate_two(tmp_path):
    acc = tmp_path / "SRR1"; acc.mkdir()
    _fastq(acc / "SRR1_1.fastq", 90); _fastq(acc / "SRR1_2.fastq", 90); _fastq(acc / "SRR1.fastq", 10)
    result = verify_download("SRR1", acc, 100)
    assert (result["verdict"], result["mates_match"]) == ("complete", True)


def test_verify_download_single_end_is_unchanged(tmp_path):
    acc = tmp_path / "SRR1"; acc.mkdir()
    _fastq(acc / "SRR1.fastq", 100)
    result = verify_download("SRR1", acc, 100)
    assert (result["verdict"], result["reads_r2"], result["mates_match"]) == ("complete", None, True)


def test_mate_mismatch_message_round_trips():
    parsed = parse_verdict_message("Downloaded 2 files, truncated (mate files differ: 100 and 40 reads)")
    assert parsed["verdict"] == "truncated" and parsed["reads_r2"] == 40
```

Plus `test_download_accession_with_short_mate_two_fails_without_accept_partial`: mock fasterq-dump the way
the module's existing `download_accession` tests do so it writes 100/40 mate files; assert the result is
`(False, ...)` containing `mate files differ`. In `tests/test_store_sidecar.py`:
`test_build_sidecar_short_mate_two_is_partial` asserting `state == "partial"`.

- [ ] **Step 2: Run and confirm failure**

Run: `conda run -n metaquest pytest tests/test_data_sra_download.py tests/test_store_sidecar.py -k "mate" -v`
Expected: FAIL (`KeyError: 'reads_r2'` or verdict `complete`).

- [ ] **Step 3: Implement** as described. Without compression `_plain_verdict` passes no counts, so
  `verify_download` now reads mate 2 as well; note this extra pass in the docstring.

- [ ] **Step 4: Run** `conda run -n metaquest pytest tests/test_data_sra_download.py tests/test_store_sidecar.py tests/test_sra_spots.py tests/test_cli_sra_verdicts.py tests/test_cli_store.py tests/test_registry_reconcile.py -v`
Expected: PASS. `store_verify` now reports a mate mismatch as truncated; adjust any fixture that wrote
unequal mates by accident.

- [ ] **Step 5: Commit**

```bash
git add prepare/data/sra/fastq.py prepare/data/sra/accession.py prepare/store/sidecar.py tests/test_data_sra_download.py tests/test_store_sidecar.py
git commit -m "fix: treat a download whose mate files differ in read count as truncated"
```

---

### Task 4: A failed usage write after linking is an error, not a warning

Defect: `record_usage_many` returns False when the catalogue write fails (for example the catalogue lock
wait runs out), and `store_link` (`cli/commands/store/link.py:178`), `store_adopt` (`adopt.py:254,273`)
and `download_sra` (`cli/commands/sra.py:401`) ignore it and exit 0. The project is then absent from the
catalogue, `linked_by` cannot see its symlink, and after the 24 h `.used` grace `store_gc --yes` removed the
dataset (reproduced by the store agent).

Behaviour: when `record_usage_many` returns False in those four places, the command finishes the rest of its
work and returns exit 4 through `self.fail(LockTimeoutError(msg), self.name)`, with
`msg = "Linked <n> accession(s) but could not record their use in the store catalogue; run store_link <accessions> again before store_gc"`
(for `download_sra` and `--copy`, the same message with "Downloaded"/"Copied"). The accession list in the
message is capped at 10 names followed by "and N more".

**Files:**
- Modify: `prepare/cli/commands/store/link.py:178`, `prepare/cli/commands/store/adopt.py:254,273`,
  `prepare/cli/commands/sra.py:399-401`
- Modify: `docs/shared_store.md` (one sentence under the store_link/adopt description: exit 4 means rerun)
- Test: `tests/test_cli_store.py`, plus the existing `download_sra` store test module (find it with
  `grep -ln "source.*store" tests/test_cli_*.py`)

**Interfaces:**
- Produces: a module-level helper in `prepare/store/usage.py`:
  `usage_not_recorded(verb: str, accessions: Sequence[str]) -> LockTimeoutError` building the message above,
  used by all four call sites.

- [ ] **Step 1: Write the failing tests**

```python
def test_store_link_exits_4_when_usage_cannot_be_recorded(self, tmp_path, monkeypatch, caplog):
    # same setup as test_links_accession_and_updates_registry
    monkeypatch.setattr("prepare.cli.commands.store.link.record_usage_many", lambda *a, **k: False)
    rc = StoreLinkCommand().execute(_link_args(["SRR1"], data_root=str(root), registry=str(registry_path),
                                               fastq_folder=str(project_dir / "fastq")))
    assert rc == 4
    assert (project_dir / "fastq" / "SRR1").is_symlink()
    assert "store_link SRR1" in caplog.text
```

Equivalent tests for `store_adopt` (move and `--copy`), patching `prepare.cli.commands.store.adopt.record_usage_many`,
and for `download_sra` with an already-downloaded store accession, patching
`prepare.cli.commands.sra.record_usage_many`. A unit test for `usage_not_recorded` with 12 accessions asserts
"and 2 more".

- [ ] **Step 2: Run and confirm failure** — `conda run -n metaquest pytest tests/test_cli_store.py -k "usage_cannot" -v`; expected rc 0 != 4.

- [ ] **Step 3: Implement.** In `adopt.py` the two recording helpers return the bool so `execute` can call
  `self.fail` once after the report is printed.

- [ ] **Step 4: Run** `conda run -n metaquest pytest tests/test_cli_store.py tests/test_store_usage.py tests/test_cli_commands.py tests/test_registry_concurrent_commands.py -v`; expected PASS.

- [ ] **Step 5: Commit**

```bash
git add prepare/store/usage.py prepare/cli/commands/store/link.py prepare/cli/commands/store/adopt.py prepare/cli/commands/sra.py docs/shared_store.md tests/test_cli_store.py tests/test_store_usage.py <download test module>
git commit -m "fix: exit 4 when a store link or download cannot record its use in the catalogue"
```

---

### Task 5: `store_gc --dry-run` never clears the rebuilt flag

Defect: `cli/commands/store/gc.py:394-397` deletes `REBUILT_WITHOUT_PROJECTS` whenever `--accept-rebuilt`
is given and a project is recorded, including under an explicit `--dry-run`. Clearing the flag on a plain
`store_gc --accept-rebuilt` (no `--yes`, no `--dry-run`) is the documented acknowledgement and stays, as
`test_accept_rebuilt_clears_the_flag_and_gc_proceeds` asserts.

Behaviour: with `--dry-run`, the flag is kept and the command logs
`"--dry-run: the flag set when the catalogue was rebuilt without projects on <date> is kept"`.

**Files:**
- Modify: `prepare/cli/commands/store/gc.py:394-400`
- Test: `tests/test_store_gc.py` (the class holding the rebuilt-flag tests, ~line 700)

- [ ] **Step 1: Write the failing test** `test_dry_run_with_accept_rebuilt_keeps_the_flag`: same setup as
  `test_accept_rebuilt_clears_the_flag_and_gc_proceeds`, run with `accept_rebuilt=True, dry_run=True`;
  assert `rc == 0`, `self._flag(paths)` is still set, and `sra_dir(paths, "SRR1").is_dir()`.
- [ ] **Step 2: Run** `conda run -n metaquest pytest tests/test_store_gc.py -k rebuilt -v`; expected the new test FAILS.
- [ ] **Step 3: Implement** the `--dry-run` guard; update the `--accept-rebuilt` help text and `docs/shared_store.md`
  to say `--dry-run` keeps the flag.
- [ ] **Step 4: Run** `conda run -n metaquest pytest tests/test_store_gc.py -v`; expected PASS.
- [ ] **Step 5: Commit**

```bash
git add prepare/cli/commands/store/gc.py docs/shared_store.md tests/test_store_gc.py
git commit -m "fix: store_gc --dry-run keeps the rebuilt-without-projects flag"
```

---

### Task 6: A blacklist file that cannot be read stops the download

Defect: `_read_blacklist_files` (`prepare/data/sra/download.py:43-70`) catches `OSError` and
`UnicodeDecodeError`, logs a warning and carries on, so a mistyped `--blacklist` path or a UTF-16 file
downloads the accessions the user meant to exclude.

Behaviour: the file is opened with `encoding="utf-8-sig"` (so a UTF-8 BOM is accepted), and an
`OSError` or `UnicodeDecodeError` raises `DataAccessError(f"Cannot read blacklist file {path}: {e}")`
before any download starts. `download_sra` reports it through `self.fail` (exit 1).

**Files:**
- Modify: `prepare/data/sra/download.py:43-70`
- Test: `tests/test_data_sra_download.py`

- [ ] **Step 1: Write the failing tests**
  - `test_missing_blacklist_file_raises`: `pytest.raises(DataAccessError, match="blacklist")` on
    `_read_blacklist_files([tmp_path / "nope.txt"])`.
  - `test_utf16_blacklist_file_raises`: file written with `encoding="utf-16"`; same assertion.
  - `test_blacklist_with_bom_is_read`: file written as `"﻿SRR1\nSRR2\n"` with `encoding="utf-8"`;
    result `== {"SRR1", "SRR2"}`. (Write the BOM as the escape `"﻿"` so the test file stays ASCII.)
  - A command-level test that `download_sra --blacklist missing.txt` returns 1 and calls no download
    (patch the download function the module's existing CLI tests patch).
- [ ] **Step 2: Run** `conda run -n metaquest pytest tests/test_data_sra_download.py -k blacklist -v`; expected the raise tests and the BOM test FAIL.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run** `conda run -n metaquest pytest tests/test_data_sra_download.py tests/test_cli_blacklist.py -v`; expected PASS.
- [ ] **Step 5: Commit**

```bash
git add prepare/data/sra/download.py tests/test_data_sra_download.py
git commit -m "fix: stop download_sra when a blacklist file cannot be read"
```

---

### Task 7: Keep the NCBI API key out of logs and error messages

Defect: `SRAMetadataClient._make_request` (`prepare/data/sra_metadata.py:98-110`) sends `api_key` as a query
parameter and embeds `str(e)` of a `requests` exception, which contains the full URL, in a log line and in
the raised `NetworkError`/`DataAccessError`. Reproduced with a mocked 503.

Behaviour: a helper replaces the value of any `api_key=` query parameter with `***` in a string; the log
line and both exception messages pass through it. `raise ... from e` is kept (the chained exception is not
shown to the user by `BaseCommand.fail`; confirm this in `cli/base.py`, and if the traceback goes to the log
file, log the redacted message only and raise `from None`).

**Files:**
- Modify: `prepare/utils/http.py` (new helper), `prepare/data/sra_metadata.py:98-110`
- Test: `tests/test_sra_metadata_client.py`, `tests/test_utils_http.py` (create if absent; check `ls tests | grep http`)

**Interfaces:**
- Produces: `redact_secrets(text: str) -> str` in `prepare/utils/http.py`; replaces
  `api_key=<value>` (value up to `&`, whitespace, quote or end) with `api_key=***`.

- [ ] **Step 1: Write the failing tests**

```python
def test_redact_secrets_hides_the_api_key():
    url = "https://eutils.ncbi.nlm.nih.gov/x?db=sra&api_key=SECRETKEY123&id=1"
    assert redact_secrets(f"503 Server Error for url: {url}") == (
        "503 Server Error for url: https://eutils.ncbi.nlm.nih.gov/x?db=sra&api_key=***&id=1"
    )


def test_make_request_failure_does_not_leak_the_api_key(caplog):
    client = SRAMetadataClient("a@b.c", api_key="SECRETKEY123")
    # mock client.session.get to raise requests.HTTPError("503 Server Error for url: ...&api_key=SECRETKEY123")
    with pytest.raises(NetworkError) as info:
        client._make_request("https://eutils.ncbi.nlm.nih.gov/x", {"db": "sra"})
    assert "SECRETKEY123" not in str(info.value)
    assert "SECRETKEY123" not in caplog.text
```

(Make the mocked `HTTPError` carry a `response` with `status_code=503` so `_is_transient` returns True.)

- [ ] **Step 2: Run** `conda run -n metaquest pytest tests/test_sra_metadata_client.py tests/test_utils_http.py -k "redact or leak" -v`; expected FAIL.
- [ ] **Step 3: Implement.** Grep `prepare/` for other places that format a `requests` exception built with
  `api_key` in `params` (`grep -rn "params\[\"api_key\"\]\|api_key=" prepare`); route any found through
  `redact_secrets` in this commit.
- [ ] **Step 4: Run** `conda run -n metaquest pytest tests/test_sra_metadata_client.py tests/test_sra_metadata_extended.py tests/test_utils_http.py -v`; expected PASS.
- [ ] **Step 5: Commit**

```bash
git add prepare/utils/http.py prepare/data/sra_metadata.py tests/test_sra_metadata_client.py tests/test_utils_http.py
git commit -m "fix: redact the NCBI API key from request error messages and logs"
```

---

### Task 8: Gates, changelog and version 0.12.3

**Files:**
- Modify: `CHANGELOG.md` (new `## [0.12.3] - <date>` section under `[Unreleased]`, one `### Fixed` bullet per
  task in the style of 0.12.2), `pyproject.toml:7`, `prepare/__init__.py:10`, `CLAUDE.md:7` (0.12.2 -> 0.12.3)

- [ ] **Step 1:** Run `conda run -n metaquest make check`; expected exit 0. Fix any black/flake8/mypy/size finding in the touched files.
- [ ] **Step 2:** Run `conda run -n metaquest make test`; expected all pass, none skipped that passed on `main`.
- [ ] **Step 3:** Write the changelog entries and bump the version in the three files.
- [ ] **Step 4:** Run `conda run -n metaquest pytest tests/test_environment_pins.py tests/test_module_sizes.py tests/test_lint_exemptions.py -v`; expected PASS.
- [ ] **Step 5: Commit**

```bash
git add CHANGELOG.md pyproject.toml prepare/__init__.py CLAUDE.md
git commit -m "release: 0.12.3"
```

Pushing, the PR and the tag wait for the user's go-ahead.

---

## Follow-up backlog (not in this plan)

Medium and Low findings of the same sweep, to be planned after 0.12.3. Mediums:
- Downloads: no socket timeout on `Entrez.efetch` (`metadata_fetch.py:137,214`); bare 403/404 regex classes a
  disk-full failure as not-found (`accession.py:41`); Ctrl-C in the retry pass leaves `<ACC>_temp`.
- Pipeline: extraction ignores the selection threshold; re-screening keeps vanished accessions;
  `sra_profile`/`sra_report` profile blacklisted and partial folders.
- Status: excluded accessions counted as wanted; funnel figures and `--stage` lists ignore the disk check;
  `project_report` counts differ from `status`.
- Store: `store_unlink` leaves usage rows, so the dataset is never reclaimed.
- CLI: negative `store_gc --older-than`, `--threads`/`--assembly-threads`, `single_sample --top-n`,
  `--registry-max-screened`; `sra_report --groups-file` accessions unvalidated; `genome_download --dry-run`
  creates its folder.
- Analysis: GC-only anomaly never flagged; adapter thresholds 0.05 vs 0.1; `sra_profile` totals include
  sampled datasets; metadata `avg_length` per spot not per read.
