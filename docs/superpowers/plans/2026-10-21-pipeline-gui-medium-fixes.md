# Pipeline and page Medium fixes 0.12.0 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the twenty-one Medium findings of
`docs/superpowers/audits/2026-10-21-pipeline-gui-audit.md` (P2-P13 on the command chain, G5-G13
on the HTML pages). Released as 0.12.0, since P2 and P3 change documented extraction
behaviour. The High findings shipped in 0.11.2; the Low findings stay open.

**Architecture:** Five waves. Wave 1 makes the screen and metadata stages write only what they
fully read and fail through `self.fail`. Wave 2 makes extraction follow the selection and the
exclusions, drops the threshold from the extraction identity, keeps singleton reads, and makes
`status --next` print runnable commands. Wave 3 closes the validate, genome and registry gaps.
Wave 4 fixes the pages. Wave 5 is docs and release. Each task is one commit on a branch
`fix/pipeline-gui-mediums` from main, test first.

**Tech Stack:** Python 3.12, pytest, pandas, plotly (the `interactive` extra; page tests skip
without it). Run everything as `cd /Users/andreassjodin/Code/prepare && conda run -n metaquest <cmd>`.

**Spec:** `docs/superpowers/audits/2026-10-21-pipeline-gui-audit.md`. Finding ids (P2-P13,
G5-G13) refer to that file.

## Global constraints

- ASCII only in `prepare/`, `tests/`, `scripts/`; no `print(` outside `cli/base.py`; every output
  file through the atomic helpers; modules under 800 lines (`read_extraction.py` and
  `registry.py` are `SPLIT_MODULES`, under 700); radon rank below D and the rank-C budget
  `MAX_RANK_C` in `tests/test_module_sizes.py` (105; extract a helper rather than raise it);
  the B902 and docstring exemption counts in `tests/test_lint_exemptions.py` do not grow;
  black at 120.
- Errors through `self.fail` and the exception classes; exit codes 0, 1, 2, 3, 4, 130 only.
- Tests never touch the network, a real tool or a real store; store tests use a `tmp_path`
  store with `PREPARE_DATA` and `HOME` monkeypatched; a page test uses `pytest.importorskip`.
- Commit trailer `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`; never `git add -A`.
- Modest scientific language in docstrings and docs.
- Gate before every commit: `conda run -n metaquest make check` and the owned tests; full
  `make test` at the end of each wave.
- Behaviour changes are listed in CHANGELOG under 0.12.0. New flags: `extract_target_reads
  --ignore-selection`; `use_branchwater` gains nothing (0.11.2 already refreshes).

## Rulings taken in this plan

- **P2**: `extract_target_reads` always skips excluded accessions, and when the registry holds a
  selection it extracts the selected accessions that meet `--threshold`; without a selection it
  uses the table alone as today. `--ignore-selection` restores the table-only choice. The
  selection wins because `select_datasets` is the documented step that decides what is
  downloaded and analysed; the threshold stays as a floor so a lowered selection threshold does
  not silently pull in weak hits. Cost if wrong: a project that selected at 0.95 and wants to
  extract at 0.5 passes `--ignore-selection` or reselects.
- **P3**: the threshold leaves the extraction identity (`_record_matches` no longer compares
  it) but is still recorded for information. Mapping does not depend on it.
- **P10**: singleton reads (`<genome>_s.fastq.gz`) are kept in the record and given to megahit
  as unpaired reads (`-r`). They mapped; leaving them out of the assembly was an omission, and
  `status --init` already lists the file. Cost if wrong: assemblies change slightly for
  reruns, which the assembly marker records.
- **P13**: a FASTA under `genomes/` becomes a project genome only when its id also appears in
  the matches folder, the targeted folder or a `genome_prepare` record. A file alone is a
  file. Cost if wrong: a genome downloaded by hand and never screened is not listed by
  `project_report`, which is the intended reading of "project genome".
- **G5**: the four comparison tests get a Holm correction; the page shows the raw and the
  adjusted p-value and "significant" follows the adjusted one. Four tests at 0.05 each is the
  textbook case for a correction and the cost is one line per test.
- **G13**: the funnel chart becomes a bar chart "Datasets per stage"; the stages are not nested
  (analysed, extracted and assembled do not require selected), so a funnel misreads.

## Review focus

1. Extraction never maps an excluded accession, with or without a selection, with or without
   `--ignore-selection` (Task 2.1).
2. A project extracted before 0.12.0 at any threshold skips every sample on a rerun with the
   same genome, preset and mapq (Task 2.2); a recorded threshold is ignored, never rejected.
3. `parse_containment` with one unreadable match file leaves `parsed_containment.txt` and the
   registry exactly as they were (Task 1.1).
4. `sra_report` with the corrected p-values still renders the pages of
   `tests/test_cli_sra_report.py` and the registry summary keeps its field names (Task 4.1).
5. The explorer without a taxonomy map shows no "Unknown" section and names `enrich_taxonomy`
   (Task 4.2).

---

## Wave 1: the screen and metadata stages write only what they read

### Task 1.1: `parse_containment` fails before writing on a partial read (P5, P6)

**Files:** `prepare/data/branchwater.py` (`parse_containment_data`, `_process_genome_containments`),
`prepare/cli/commands/containment.py` (`ParseContainmentCommand.execute`),
`tests/test_data_branchwater.py`, `tests/test_cli_commands.py`.

- [ ] **Step 1: Failing tests.** `test_parse_containment_with_an_unreadable_file_writes_nothing`:
  two match files, one with the wrong header; after the call `parsed_containment.txt` does not
  exist (or an earlier copy is byte-identical) and `errors == ["bad.csv"]`.
  `test_parse_containment_command_exit_1_leaves_registry_unchanged`: same through the command;
  the registry's screening blocks are unchanged and the exit is 1.
  `test_parse_containment_on_an_empty_matches_folder_is_exit_1`: an empty folder, or one with no
  `*.csv`, gives exit 1 with the message naming the folder, and the registry is not touched.
  `test_a_semicolon_separated_csv_names_its_separator`: a `;`-separated file is reported as
  "not comma separated (found ';')", not as missing columns.
- [ ] **Step 2: Run** `pytest tests/test_data_branchwater.py tests/test_cli_commands.py -k parse_containment -q`; FAIL.
- [ ] **Step 3: Implement.** In `parse_containment_data` collect every file's rows first and
  raise `DataAccessError` listing the unreadable files before `_write_details_table` or the
  summary is written; an empty folder raises `DataAccessError("No match CSV in <folder>")`.
  In `_process_genome_containments` sniff the header line: when it holds no comma but a `;`
  or TAB, raise with the separator named. In the command, catch the error with `self.fail`
  before the registry transaction; delete the `errors` list plumbing that is now unused.
- [ ] **Step 4: Run** the two files; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: parse_containment writes nothing on a partial or empty read`.

### Task 1.2: metadata commands report failures in the exit code (P9)

**Files:** `prepare/data/metadata.py` (`download_metadata` failures, the matches CSV reader near
line 89, `parse_metadata` near line 497), `prepare/cli/commands/metadata.py`,
`tests/test_data_metadata.py`, `tests/test_cli_metadata.py` (or where the commands are tested).

- [ ] **Step 1: Failing tests.** `test_download_metadata_exit_1_when_an_accession_has_no_metadata`:
  with the HTTP session mocked to answer for one of two accessions, the command exits 1 and the
  summary names the missing accession; `test_download_metadata_unreadable_matches_csv_is_exit_1`
  (the three sibling readers already do this); `test_parse_metadata_exit_1_on_a_broken_xml`:
  one of two XML files fails to parse, the table is written from the good one, exit 1;
  `test_parse_metadata_with_no_xml_writes_no_table_and_exits_1`: an older
  `metadata_table.txt` is left as it is and the exit is 1 with a message.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** `download_metadata` returns the failures to the command (it
  already has the dict); the command returns 1 when any remains after the network case (exit
  4 stays for an all-network failure). The matches reader raises `DataAccessError` like
  `parse_containment`. `parse_metadata` counts parse errors into the return value and the
  command exits 1 when any; no XML is a `DataAccessError`.
- [ ] **Step 4: Run** the files; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: metadata failures reach the exit code`.

### Task 1.3: the screen commands fail through `self.fail` (P17, the part that touches this wave)

**Files:** `prepare/cli/commands/branchwater.py`, `containment.py`, `select.py`,
`prepare/data/branchwater.py` (`extract_branchwater_metadata`, P18), their tests.

- [ ] **Step 1: Failing tests.** For `use_branchwater`, `extract_branchwater_metadata` and
  `select_datasets --no-record`: the error path logs through `fail` (caplog shows the
  traceback attached at DEBUG, one ERROR line) and returns the exception's code.
  `test_extract_branchwater_metadata_unreadable_file_is_an_error`: a file in another format
  is in `errors` and the exit is 1.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Replace `logger.error(...); return 1` with raising the matching
  `PrepareError` subclass and `return self.fail(e, context)`; `extract_branchwater_metadata`
  reads the header with the CSV module and appends to `errors`.
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: screen commands fail through the shared path`.

**Wave gate:** `make test` green.

## Wave 2: extraction follows the project

### Task 2.1: extraction respects the selection and the exclusions (P2)

**Files:** `prepare/data/read_extraction.py` (`select_samples_for_genome`, `extract_target_reads`),
`prepare/cli/commands/read_extraction.py` (`selected_samples`, `--ignore-selection`),
`prepare/data/registry.py` (`query(registry, "excluded")`, `query(registry, "selected")` already
exist), `tests/test_read_extraction.py`, `tests/test_cli_read_extraction.py`.

- [ ] **Step 1: Failing tests.** `test_excluded_accessions_are_never_extracted`: a project with
  SRR1 (0.9, blacklisted) and SRR2 (0.9): the dry run lists SRR2 only, and so does a real run
  with the tools faked. `test_a_recorded_selection_bounds_the_extraction`: selection at 0.5
  recorded for SRR2 and SRR3, SRR4 at 0.3 downloaded: the run extracts SRR2 and SRR3.
  `test_ignore_selection_uses_the_table`: the same project with `--ignore-selection` extracts
  SRR4 too. `test_without_a_selection_the_table_decides`: no selection block, table at the
  threshold as today.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** `selected_samples(table, genome_id, threshold, registry,
  ignore_selection)`: the table's samples at the threshold, minus `query(registry,
  "excluded")`, intersected with `query(registry, "selected")` when that is non-empty and
  `ignore_selection` is off. Log one INFO line with the three counts (table, excluded,
  selected). The flag help: "Extract every sample at --threshold, not only the selected ones;
  excluded samples are always skipped". Keep `extract_target_reads` under rank D by passing the
  final list in.
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: extraction follows the selection and the exclusions`.

### Task 2.2: the threshold leaves the extraction identity (P3)

**Files:** `prepare/data/read_extraction.py` (`_record_matches`, `_record_status`),
`prepare/data/registry_blocks.py` (`ExtractionBlock.threshold` stays),
`tests/test_read_extraction.py`, `docs/project_state.md`, `docs/extraction_assembly.md`.

- [ ] **Step 1: Failing test.** `test_a_different_threshold_does_not_remap`: a record with
  `threshold: 0.5` and the same genome, preset and mapq is `already_done` at 0.1.
  `test_a_different_mapq_still_remaps` pins the rest of the key.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Remove the threshold comparison; the docstring says the key is
  genome FASTA, preset and mapq, and that the threshold is recorded for information. Update the
  two docs' sentences on the skip key (they also omit `min_mapq`).
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: the containment threshold is not part of the extraction identity`.

### Task 2.3: `status --next` prints runnable, relevant commands (P4)

**Files:** `prepare/cli/commands/status/suggest.py` (`_extraction_next_steps`, the download
suggestion), `tests/test_cli_status.py`.

- [ ] **Step 1: Failing tests.** `test_extraction_is_suggested_only_for_screened_pairs`: with
  genomes G1 and G2 and SRR3 screened against G1 only, the G2 command lists no SRR3;
  `test_suggestions_quote_paths_and_ids`: a genome id `L crispatus` and a table path with a
  space appear quoted (`shlex.split` of the command gives the original tokens).
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Filter `to_extract` to accessions whose screening block holds the
  genome (or, without a screening block, any downloaded accession, so a `status --init`
  project is not silenced); `shlex.quote` every interpolated value in both suggestion
  builders.
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: status --next suggests screened pairs and quotes its commands`.

### Task 2.4: singleton reads are recorded and assembled (P10)

**Files:** `prepare/data/read_extraction_export.py` (`_export_mapped_fastq`),
`prepare/data/assembly.py` (`_megahit_args`), `prepare/data/registry_scan.py` (unchanged),
`tests/test_read_extraction.py`, `tests/test_cli_extraction_assembly.py`, `docs/extraction_assembly.md`.

- [ ] **Step 1: Failing tests.** `test_paired_export_keeps_the_singleton_file`: with `_1`, `_2`
  and a non-empty `_s`, the returned files are the three; `test_megahit_gets_singletons_as_r`:
  the argument list holds `-1`, `-2` and `-r <_s>`; without `_s` no `-r`.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Return `[out1, out2, singles]` when `singles` exists; the assembly
  builder maps a file ending `_s.fastq.gz` to `-r`. Describe `_s` and `_0` in the doc.
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: singleton reads are recorded and given to megahit`.

### Task 2.5: a failed extraction record is not swallowed (P12)

**Files:** `prepare/data/read_extraction.py` (`_notify_result`), `tests/test_read_extraction.py`.

- [ ] **Step 1: Failing test.** `test_a_lock_timeout_while_recording_is_a_transient_error`: the
  callback raises `LockTimeoutError`; the run finishes the remaining samples, then
  `extract_target_reads` raises the first `TransientError` so the command exits 4; a plain
  `DataAccessError` from the callback stays a warning (the existing behaviour, pinned).
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** `_notify_result` keeps a list of `TransientError`s and the caller
  re-raises the first after the loop; other `PrepareError`s stay warnings; anything else is not
  caught (a bug propagates).
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: a lock timeout while recording an extraction exits 4`.

**Wave gate:** `make test` green.

## Wave 3: validate, genomes and the registry

### Task 3.1: `sra_validate` names what it skips and records only accessions (P8)

**Files:** `prepare/cli/commands/sra_enhanced.py` (`_find_accession_dirs`, the result loop,
the generic handler), `prepare/data/registry.py:608` (the `analysed` stage), `tests/test_cli_commands.py`
(or the validate tests), `tests/test_data_registry.py`.

- [ ] **Step 1: Failing tests.** `test_named_accessions_without_a_folder_fail_the_run`: two
  named, one folder: the missing one is logged as an error, counted failed, exit 1;
  `test_folder_scan_skips_non_accession_and_transient_folders`: `fastq/old_runs` and
  `fastq/.prepare-tmp` are neither validated nor recorded; `test_a_failed_validation_is_not_analysed`:
  a dataset whose only analysis is a FAILED validate is not in `query(registry, "analysed")`.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Reuse `validate_accession` (`core/validation.py`) and
  `is_transient_folder` (`data/sra/cleanup.py`) in the scan, as `sra_sketch` does; report named accessions with no folder; `query(..., "analysed")`
  counts an analysis only when its record has no `failed`/`status == "FAILED"` marker (check the
  validate block's field name in `registry_blocks.py`). Replace the generic `except Exception`
  with `self.fail` for `PrepareError` and `OSError` and let the rest propagate, then remove the
  module from the B902 list (lower `B902_EXEMPT`).
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: sra_validate reports missing accessions, skips stray folders; a failed validation is not analysed`.

### Task 3.2: a gzipped genome FASTA works everywhere (P11)

**Files:** `prepare/data/branchwater_search.py` (`sketch_fasta`), `prepare/cli/commands/genome.py`
(`_genome_name`, move it to `prepare/data/genome_names.py` or `core/validation.py`),
`prepare/cli/commands/branchwater_search.py`, `prepare/cli/commands/sra_sketch.py` (`--genome-id`
default), `tests/test_branchwater_search.py`, `tests/test_cli_sra_sketch.py`.

- [ ] **Step 1: Failing tests.** `test_sketch_fasta_reads_gzip` (skip without sourmash);
  `test_genome_id_of_a_gz_fasta_drops_both_suffixes` for `branchwater_search` and `sra_sketch`
  (`G.fna.gz` -> `G`).
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Open with `gzip.open(..., "rt")` on the gzip magic; one shared
  `genome_name(path)` used by `genome_prepare`, `branchwater_search`, `sra_sketch` and
  `_genome_ids_on_disk`.
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: a gzipped genome FASTA is read and named like a plain one`.

### Task 3.3: `status --init` registers genomes with evidence (P13)

**Files:** `prepare/data/registry_scan.py` (`_genome_ids_on_disk`, the loop at line 267),
`tests/test_cli_status.py` (where `status --init` is tested).

- [ ] **Step 1: Failing test.** `test_a_fasta_alone_is_not_a_project_genome`: `genomes/X.fna`
  with no matches, targeted or record: `registry.genomes` lacks `X`; with `matches/X.csv` it has
  it; `test_extraction_folders_still_name_their_genome` pins the targeted path.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** The FASTA scan only confirms ids found elsewhere (it still resolves
  the FASTA path for `_genome_fasta`).
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: status --init registers a genome only with screening or extraction evidence`.

**Wave gate:** `make test` green.

## Wave 4: the pages

### Task 4.1: honest p-values (G5; G4's note already renders)

**Files:** `prepare/sra/analytics.py` (`compare_datasets`, `_generate_comparative_recommendations`),
`prepare/sra/reporting.py` (the tests section), `prepare/cli/commands/sra_report.py`
(`_print_comparison`), `tests/test_sra_analytics.py`, `tests/test_cli_sra_report.py`.

- [ ] **Step 1: Failing tests.** `test_holm_adjusted_p_values_decide_significance`: four tests
  with raw p 0.01, 0.02, 0.03, 0.2: adjusted 0.04, 0.06, 0.06, 0.2 and only the first
  significant; `test_report_shows_raw_and_adjusted_p`: the HTML carries "p 0.0100, adjusted
  0.0400" for that test; the registry summary field names are unchanged.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** Add `p_adjusted` (Holm over the tests that ran) and set
  `significant` from it; the template and the terminal print both; the recommendation text
  says "after Holm correction".
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: sra_report corrects its four tests for multiple testing`.

### Task 4.2: the explorer without taxonomy says so (G6, G10)

**Files:** `prepare/visualization/explorer.py` (`_build_summary_data`, `_assemble_html`, the
template), `prepare/cli/commands/explore.py`, `tests/test_explorer.py`.

- [ ] **Step 1: Failing tests.** `test_no_taxonomy_hides_the_taxonomy_sections`: with every
  genome "Unknown" the HTML has no sunburst, heatmap or distributions, the cards show samples
  and genomes only, and a note names `enrich_taxonomy`; `test_an_empty_result_explains_itself`:
  no row above `--min-containment`: the page says so and `total_genomes` is the number of
  genome columns.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** `has_taxonomy = (long_df["family"] != "Unknown").any()`; the
  template wraps the three sections and the two cards in it and otherwise renders the note;
  the empty branch keeps `total_genomes` and renders "No sample reaches containment X; lower
  --min-containment".
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: the explorer says when taxonomy or results are absent`.

### Task 4.3: the explorer's cut-offs and aggregates are labelled, placeholders blanked (G7, G8, G9)

**Files:** `prepare/visualization/explorer.py` (`_build_distributions`, `_build_heatmap`,
`_metadata_text`), `tests/test_explorer.py`.

- [ ] **Step 1: Failing tests.** `test_genus_bar_title_names_the_cut`: with 31 genera the title
  reads "Samples per genus (top 30 of 31)"; `test_heatmap_title_says_maximum_per_family`;
  `test_metadata_placeholders_are_blank`: "NP", "uncalculated", "not provided", "missing",
  "not applicable", "unknown" (case-insensitive) render as empty cells.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** One `_PLACEHOLDERS` set in `_metadata_text`; titles carry the
  rule.
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: explorer titles state their cut-offs; metadata placeholders are blank`.

### Task 4.4: `interactive_plot` writes a file or fails (G11, G12)

**Files:** `prepare/cli/commands/advanced_analysis.py` (`InteractivePlotCommand`),
`prepare/visualization/interactive.py` (`create_interactive_heatmap`, `_write_figure`),
`tests/test_interactive_plots.py`, `tests/test_cli_commands_advanced_analysis.py`.

- [ ] **Step 1: Failing tests.** `test_no_output_file_is_a_validation_error`: without
  `--output-file` the command exits 1 naming the flag (`--no-show` becomes the default: the
  page is written, never shown, and `--show` opens it); `test_heatmap_uses_a_sequential_scale_for_nonnegative_data`
  (`Viridis` when `z.min() >= 0`, `RdBu_r` otherwise, centred at 0);
  `test_sample_metadata_colours_the_rows` or the parameter is removed from the signature and
  the CLI (choose removal unless a caller needs it: it is unused);
  `test_heatmap_size_is_capped` at `EXPLORER_HEATMAP_MAX_HEIGHT`; `test_written_page_has_a_title_and_viewport`
  (the pages get `REPORT_CSS`'s head through `_write_figure`, G14 comes for free). Errors go
  through `self.fail`.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement** as the tests say; remove the `advanced_analysis.py` entry from the
  B902 list if the broad handler goes.
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: interactive_plot always writes its page; sequential scale and a size cap for the heatmap`.

### Task 4.5: the project report's stage chart is a bar chart (G13)

**Files:** `prepare/visualization/project_report.py` (`_funnel_figure`, `_genome_figure`),
`tests/test_project_report.py`.

- [ ] **Step 1: Failing tests.** `test_stage_chart_is_a_bar_chart`: the figure's trace type is
  `bar` and its title "Datasets per stage"; `test_genome_chart_axis_starts_at_zero`: with all
  zeros the y axis range is `[0, 1]`.
- [ ] **Step 2: Run**; FAIL.
- [ ] **Step 3: Implement.** `go.Bar` with the stage order kept; `yaxis_rangemode="tozero"` on
  both charts.
- [ ] **Step 4: Run**; PASS. `make check`.
- [ ] **Step 5: Commit** `fix: project_report draws stages as bars, axes from zero`.

**Wave gate:** `make test` green.

## Wave 5: docs and release

### Task 5.1: docs, changelog and 0.12.0

**Files:** `CHANGELOG.md`, `docs/pipeline_overview.md`, `docs/extraction_assembly.md`,
`docs/project_state.md`, `docs/reads.md`, `docs/analysis_commands.md`, `README.md` (the
`extract_target_reads` row if it names the threshold rule), `docs/superpowers/README.md`,
`pyproject.toml`, `prepare/__init__.py`, `CLAUDE.md`, `docs/packaging.md`.

- [ ] **Step 1: Write** the 0.12.0 section: Breaking (P2 extraction follows the selection;
  P3 identity; P10 singletons in the assembly; `interactive_plot` requires `--output-file`),
  Fixed (the rest), Internal (B902 list shrinks by the modules freed). Update each doc's
  sentence named in the audit's P19 that these tasks touched; mark the audit's Medium findings
  shipped in the design-record index.
- [ ] **Step 2: Run** `python scripts/check_version.py`, `python scripts/check_docs_commands.py`,
  `make check`, `make test`; all green.
- [ ] **Step 3: Commit** `release: 0.12.0`.

The branch is then ready for a fresh whole-branch review, a pull request, and a merge on green
on the user's instruction. The Low findings (P14-P20, G14-G21) remain for a later round.
