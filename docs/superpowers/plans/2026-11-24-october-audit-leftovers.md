# October audit leftovers (0.15.0) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the Low findings of the 2026-10-21 pipeline and page audit that are still present on `main`
(17333fd), and the two minors deferred by the 0.12.1 review.

**Architecture:** Eleven tasks: five on the pipeline commands (path safety, cANI, path anchoring, exit codes,
the two deferred minors), one on the documentation, four on the HTML pages, then the release. Each defect gets
a test that fails on `main` before its fix.

**Tech Stack:** Python 3.12, pytest, plotly and jinja2 (`interactive` extra); conda environment `metaquest`.

**Spec:** `docs/superpowers/audits/2026-10-21-pipeline-gui-audit.md` (P14-P20, G14-G21) and the re-check of
2026-11-24 (in the conversation; summarised per task). P18 and P20 are fixed and not in this plan.

**Why 0.15.0:** Task 3 changes where project files are read and written when a command runs from a subfolder,
and Task 4 makes `sra_report` exit 1 when an accession is left out. The fixes merged in PR #53 (now under
`[Unreleased]`) ship in the same release.

## Global Constraints

- Branch `fix/october-audit-leftovers` off `main` (17333fd). Stage files by name; never `git add -A`.
- Run as `cd /Users/andreassjodin/Code/prepare && conda run -n metaquest <cmd>` (plain `make test` picks the
  base Python 3.11 and fails collection); the cwd resets between calls.
- Test first, watched failing, then the smallest fix (CLAUDE.md).
- Module ceiling 800 lines, 700 for `SPLIT_MODULES`; `cli/commands/read_extraction.py` is at 768, so code
  added there goes into a helper module if it needs more than ~25 lines. No new radon rank C
  (`MAX_RANK_C` = 104 must not grow), no rank D.
- Errors through `self.fail`; exit codes 0, 1, 2, 3, 4, 130 only. Atomic writes, ASCII only, no `print(`.
- Tests that need plotly/jinja2 skip when the `interactive` extra is absent.
- Modest scientific language. Update README/docs/CHANGELOG with each behaviour change, in the same task.
- `make check` and `make test` before each commit.

## Review Focus

1. An absolute path given to any flag of Task 3 is used unchanged; only relative paths move.
2. A user-written input list (`--accessions-file`, `--groups-file`, `blacklist --from-file`) given relative
   to the working directory is still found from a subfolder (it is not anchored).
3. `sra_sketch --genome-id` values that are plain names with dots or dashes (`E.coli_K-12`, `GCF_000005845.2`)
   are still accepted.
4. A rescreen whose containment for an accession changed must not keep the old cANI.
5. The pages still render with no registry (the explorer and `sra_report` run outside a project).

Tests for each line are in the owning task (1 and 2 in Task 3, 3 in Task 1, 4 in Task 2, 5 in Task 7).

---

### Task 1: `sra_sketch --genome-id` is a plain name (P16)

**Files:** Modify `prepare/core/constants.py`, `prepare/cli/commands/sra_sketch.py`; Test
`tests/test_cli_sra_sketch.py` (or the existing `sra_sketch` CLI test file).

**Interfaces:** Produces `GENOME_ID_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9._-]*$"` in `core/constants.py`.

- [ ] **Step 1:** Tests: `--genome-id a/b` and `--genome-id ../x` exit 2 with an error naming the value, and no
  file is written outside `matches/`; `--genome-id E.coli_K-12` and `GCF_000005845.2` write
  `matches/<id>.csv` (Review Focus 3).
- [ ] **Step 2:** Run them; expected: the two refusals FAIL (a subfolder or an escaped file is written).
- [ ] **Step 3:** In `execute`, before any work, a `--genome-id` not matching `GENOME_ID_PATTERN` is
  `return self.fail(ValidationError(f"--genome-id must be a plain name (letters, digits, '.', '_', '-'): {value!r}"), self.name)`.
  The default from the genome file's name is not checked (it comes from a file name). Help text names the rule.
- [ ] **Step 4:** Run the tests; expected PASS. `make check`, `make test`.
- [ ] **Step 5:** Commit `fix: sra_sketch --genome-id must be a plain name`.

### Task 2: `parse_containment` keeps cANI (P14)

**Files:** Modify `prepare/data/screening_table.py` (`_write_row`); Test `tests/test_screening_table.py`
(or the existing test module for `record_screening_table`).

- [ ] **Step 1:** Tests: an accession whose registry entry for genome G has `containment=0.5, cani=0.98`
  (as `branchwater_search` records it) keeps `cani == 0.98` after `record_screening_table` records 0.5 for
  it again; with a table value of 0.6 the entry has `cani is None` (Review Focus 4).
- [ ] **Step 2:** Run them; expected: the first FAILS (`cani` is None).
- [ ] **Step 3:** In `_write_row`, keep the existing entry's `cani` when its `containment` equals the new
  rounded value; otherwise `None`. No other field changes.
- [ ] **Step 4:** Run the tests; PASS. Gates.
- [ ] **Step 5:** Commit `fix: parse_containment keeps the cANI of an unchanged screening entry`.

### Task 3: Project files resolved from the project root (P15)

**Files:** Modify `prepare/cli/base.py`, and the commands below; README.md and `docs/` (the paragraph on
folders in `docs/project_state.md` and each command's guide); Test `tests/test_cli_anchoring.py` (new).

**Interfaces:** Produces `anchor_args(args: argparse.Namespace, *dests: str) -> None` in `cli/base.py`:
each named attribute that is not None is replaced by `str(anchored_folder(value, getattr(args, "registry", None)))`.

Rule (ruling for the reviewer): a flag whose value is a file the pipeline itself writes or reads back (its
default is a project file name) is anchored, like the folder flags; a list the user writes by hand stays
relative to the working directory. Anchored, by command:

| Command | Flags |
|---|---|
| `parse_containment` | `--parsed-containment-file`, `--summary-containment-file`, `--details-file` |
| `select_datasets` | `--parsed-containment`, `--metadata-file`, `--output` |
| `extract_target_reads` | `--parsed-containment` |
| `results_table` | `--parsed-containment`, `--output` |
| `parse_metadata` | `--metadata-table-file` |
| `check_metadata_attributes`, `count_metadata`, `plot_metadata_counts`, `plot_map` | their table and output file flags (`--file-path`, `--summary-file`, `--metadata-file`, `--stat-file`, `--output-file`, `--output`) |
| `sra_profile` | `--output-report`, `--output-dir` |
| `sra_report` | `--output-dir`, `--quality-profiles` |
| `blacklist` | `--blacklist-file` |
| `status` | `--parsed-containment` |
| `enrich_taxonomy`, `explore_containment`, `find_by_taxonomy` | `--parsed-containment`, `--output`, `--cache` |

Not anchored: `--accessions-file`, `--groups-file`, `blacklist --from-file`.

- [ ] **Step 1:** Tests, each in a `tmp_path` project with a registry and run with `monkeypatch.chdir(project / "sub")`:
  `parse_containment` writes `project/parsed_containment.txt`; `select_datasets` reads it and writes
  `project/accessions.txt`; `sra_profile` writes `project/sra_quality_profiles/`; `blacklist` writes
  `project/blacklist.txt`; an absolute `--output` is used unchanged (Review Focus 1); `status
  --accessions-file list.txt` with `sub/list.txt` reads that file (Review Focus 2).
- [ ] **Step 2:** Run them; expected: the anchoring tests FAIL (files land in `sub/`), the two others PASS.
- [ ] **Step 3:** Add `anchor_args`; call it at the top of each command's `execute`, inside the existing
  `try` that already anchors folders where there is one. Help texts of anchored flags end with
  "(relative to the project root)".
- [ ] **Step 4:** Run the tests; PASS. Gates.
- [ ] **Step 5:** README, `docs/project_state.md` (one paragraph: which paths are resolved from the project
  root, which from the working directory) and CHANGELOG `### Changed`.
- [ ] **Step 6:** Commit `fix: project files are resolved from the project root, as project folders are`.

### Task 4: Exit codes (P17)

**Files:** Modify `prepare/cli/commands/sra_report.py`, `sra_profile.py`, `metadata.py` (`parse_metadata`),
`sra_enhanced.py` (`sra_info`); `docs/reads.md`; Tests in the existing CLI test modules of those commands.

- [ ] **Step 1:** Tests: `sra_report` with one readable and one unreadable accession writes the report and exits
  1, and the error names the left-out accession; `parse_metadata` with one corrupt XML exits 1 through
  `self.fail` (one ERROR line, `caplog` record has `exc_info`); `sra_info` with no metadata fetched exits 1 the
  same way.
- [ ] **Step 2:** Run them; expected: the `sra_report` test FAILS (exit 0); the other two FAIL on `exc_info`.
- [ ] **Step 3:** `sra_report`: after writing the report, `return self.fail(DataAccessError(f"{len(failed)} of
  {len(accessions)} accession(s) left out of the report (no readable FASTQ files): {names}"), self.name)`,
  names capped at ten as elsewhere. `sra_profile` reports its partial failure the same way (exit 1 as now).
  `parse_metadata` and `sra_info` raise `DataAccessError` with their current messages into `self.fail`.
- [ ] **Step 4:** Run the tests; PASS. Gates.
- [ ] **Step 5:** `docs/reads.md` and CHANGELOG `### Changed` (sra_report exit 1).
- [ ] **Step 6:** Commit `fix: sra_report exits 1 when an accession is left out; errors through self.fail`.

### Task 5: The two minors deferred in 0.12.1

**Files:** Modify `prepare/cli/commands/read_extraction.py`, `prepare/cli/commands/metadata.py`
(`_fail_for_missing`); Tests `tests/test_cli_read_extraction.py`, the `download_metadata` CLI tests.

- [ ] **Step 1:** Tests: (a) `extract_target_reads` where recording sample 1 raises `DataAccessError` while the
  stop is set exits 130 and `caplog` holds an ERROR naming sample 1 and the error; (b) the same with an earlier
  sample's recording failed and the stop set after a later sample: the ERROR line
  "Recording failed for 1 sample(s): ..." is logged before exit 130; (c) `download_metadata` with failures
  `{"SRR1": "network: timeout", "SRR2": "no record"}` exits 4 with an error saying NCBI could not be reached for
  1 of 2 accession(s).
- [ ] **Step 2:** Run them; expected: all three FAIL.
- [ ] **Step 3:** `_record_result_and_check_stop`: in the `finally`, when the stop is set and a `PrepareError`
  is in flight (`sys.exc_info()[1]`), log `"Recording the extraction result for %s failed: %s"` at ERROR before
  raising `KeyboardInterrupt`. `execute`: `except KeyboardInterrupt:` around the extraction and assembly, log the
  `record_failures` line (reuse the text of `_fail_for_recording`, extracted into a small function) and
  re-raise. `_fail_for_missing`: count the network reasons, "NCBI could not be reached for {n} of {total}
  accession(s); try again later".
- [ ] **Step 4:** Run the tests; PASS. Gates (watch the 800-line ceiling).
- [ ] **Step 5:** Commit `fix: an interrupt keeps the record-failure lines; count only network failures`.

### Task 6: Documentation (P19, the parts still wrong)

**Files:** `docs/extraction_assembly.md:22`, `docs/project_state.md:60, 136`, `docs/pipeline_overview.md`
(cANI, after Task 2).

- [ ] **Step 1:** `extraction_assembly.md`: `--data-root` names the store used to record usage and to check a
  store copy's completeness; reads are always taken from `--fastq-folder`.
- [ ] **Step 2:** `project_state.md`: add `sra_sketch` to the run-log command list; the presence check covers
  `fastq/<ACC>/` and `metadata/<ACC>_metadata.xml`, and genome FASTA files are counted.
- [ ] **Step 3:** `pipeline_overview.md`: cANI is kept from the match CSVs while an accession's containment is
  unchanged.
- [ ] **Step 4:** `python scripts/check_docs_commands.py`; `make check`.
- [ ] **Step 5:** Commit `docs: extraction store use, run-log commands, presence check, cANI`.

### Task 7: One page head and footer (G15, G16, G19 templates)

**Files:** Modify `prepare/utils/html.py`, `prepare/visualization/explorer.py`, `prepare/sra/reporting.py`,
`prepare/visualization/project_report.py`, `prepare/cli/commands/explore.py`; Tests `tests/test_html_assets.py`,
`tests/test_explorer.py`, the `sra_report` and `project_report` page tests.

**Interfaces:** Produces in `utils/html.py`:
- `utc_stamp(iso: Optional[str] = None) -> str`: "YYYY-MM-DD HH:MM UTC" for `iso` (any offset) or now.
- `page_head(title: str) -> str`: charset, viewport, `<title>` (escaped), plotly script, `REPORT_CSS`.
- `page_footer() -> str`: `<footer class="pre-footer">Generated {utc_stamp()} by PREPARE {__version__}</footer>`.

- [ ] **Step 1:** Tests: each of the three pages starts `<!DOCTYPE html>\n<html lang="en">`, has one
  `pre-skip` link, and a footer matching `Generated \d{4}-\d\d-\d\d \d\d:\d\d UTC by PREPARE \d+\.\d+\.\d+`;
  `sra_report`'s readout shows a UTC stamp; `project_report`'s `<title>` is "Project report: <name>";
  `explore_containment --title "X"` gives `<title>X</title>`; the explorer and `sra_report` pages render with
  no registry (Review Focus 5).
- [ ] **Step 2:** Run them; expected FAIL on `sra_report` (no `lang`, no skip link), the footers, the titles and
  `--title` (unknown flag, exit 2).
- [ ] **Step 3:** Add the three helpers; the three templates take `head_html`/`footer_html` (`|safe`) in place of
  their own head and footer. `sra_report` gets `lang="en"` and a skip link to its `<main id="report">`.
  `project_report`'s `created`/`updated` cells go through `utc_stamp`. `explore_containment --title` (default
  the current "PREPARE Containment Explorer") is passed to the explorer.
- [ ] **Step 4:** Run the tests; PASS. Gates.
- [ ] **Step 5:** README (`explore_containment --title`), CHANGELOG `### Fixed`.
- [ ] **Step 6:** Commit `fix: one page head and footer with a UTC stamp and version; sra_report lang and skip link`.

### Task 8: Labels, chart text and the interactive figures (G20, G21, G14)

**Files:** Modify `prepare/sra/analytics.py` or `reporting.py`, `prepare/utils/html.py`,
`prepare/visualization/explorer.py`, `prepare/visualization/interactive.py`; Tests as in Task 7 plus the
interactive page tests.

**Interfaces:** Produces `METRIC_LABELS: Dict[str, str]` beside `NUMERIC_COLUMNS`: `avg_read_length` "Mean
read length (bp)", `gc_percent` "GC content (%)", `total_reads` "Reads", `complexity_score` "Complexity
score"; and `chart_panel(html: str, label: str) -> str` in `utils/html.py`, wrapping a figure in
`<div class="pre-panel" role="figure" aria-label="{label}">`.

- [ ] **Step 1:** Tests: the group box-plot title is "GC content (%) by group" and the test heading "GC content
  (%)"; every chart panel of the three report pages has `role="figure"` and a non-empty `aria-label`; the
  explorer's sortable `th` cells carry `scope="col"`; an interactive figure (pca) has no fixed `width` and uses
  `plotly_layout()`'s fonts.
- [ ] **Step 2:** Run them; expected FAIL.
- [ ] **Step 3:** Use `METRIC_LABELS` (falling back to the column name) for titles, axis labels and headings;
  wrap the panels in `chart_panel`; add `scope="col"`; in `interactive.py` replace `width=800` with
  `autosize=True` and apply `plotly_layout()` before the figure's own settings.
- [ ] **Step 4:** Run the tests; PASS. Gates.
- [ ] **Step 5:** Commit `fix: labels with units, chart text alternatives, header scope, responsive figures`.

Ruling for the reviewer: plotly figures keep a light style in a dark colour scheme (figures cannot follow the
CSS media query without script); not in this plan.

### Task 9: Explorer table and links from the project report (G18, G17)

**Files:** Modify `prepare/visualization/explorer.py`, `prepare/utils/html.py`,
`prepare/visualization/project_report.py`, `prepare/processing/project_report_markdown.py`; Tests
`tests/test_explorer.py`, `tests/test_html_assets.py`, the project report tests.

- [ ] **Step 1:** Tests: a containment of 1.2 shows "1.200" in the cell (the bar stays at full width); the
  export downloads `<slug of the page title>.csv` and the table note says the export holds the rows shown;
  the sort regex accepts `1e-05`; the project report lists every output of an analysis (not only the first),
  in HTML and Markdown, and an output ending in `.html` that exists is a link relative to the report.
- [ ] **Step 2:** Run them; expected FAIL.
- [ ] **Step 3:** Clamp only the bar's `--v`; the JS export builds the name from `document.title`
  (lowercase, runs of other characters to `_`); numeric regex `/^-?\d+(\.\d+)?([eE][-+]?\d+)?$/`;
  `_outputs` lists all outputs, the HTML renderer links `.html` ones with `os.path.relpath` from the report's
  folder.
- [ ] **Step 4:** Run the tests; PASS. Gates.
- [ ] **Step 5:** Commit `fix: explorer shows recorded containment, names its export; project report links outputs`.

### Task 10: Remove the unused report functions (G19)

**Files:** Modify `prepare/sra/reporting.py`; `tests/test_sra_reporting_extended.py`,
`tests/test_sra_reporting_starter.py`, `tests/test_optional_imports.py`.

- [ ] **Step 1:** Point `test_sra_dashboard_without_plotly_is_an_error_not_a_degraded_file` at
  `generate_report` (the path `sra_report` uses) and run it: PASS before and after the removal.
- [ ] **Step 2:** Remove `generate_quality_dashboard`, `create_comparative_analysis`, `_generate_quality_html`,
  `_generate_comparative_html` and their tests; keep every helper `generate_report` still calls.
- [ ] **Step 3:** `grep -rn` for the four names in `prepare/`, `tests/`, `docs/`: none. Gates; lower
  `MAX_RANK_C` if the count dropped.
- [ ] **Step 4:** Commit `refactor: remove the report functions only tests called`.

### Task 11: Gates, changelog and 0.15.0

- [ ] `make check`, `make test`; CHANGELOG: rename `[Unreleased]` to `## [0.15.0] - <date>` with `### Changed`
  (Tasks 3, 4) above `### Fixed` (PR #53 and Tasks 1, 2, 5-10); bump `pyproject.toml`, `prepare/__init__.py`,
  `CLAUDE.md:7`.
- [ ] Commit `release: 0.15.0`. Push, PR, merge and tag after the user's go-ahead.

## Not in this plan

- Dark styling of plotly figures (Task 8 ruling).
- G17's links from `sra_report` and the explorer back to the project report: neither page knows where the
  project report is written.
