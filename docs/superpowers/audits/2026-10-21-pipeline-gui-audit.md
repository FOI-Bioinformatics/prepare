# Audit of the pipeline steps and the HTML pages, 2026-10-21

Scope: the six-stage command chain of `docs/pipeline_overview.md` (screen, select, download,
analyse, extract, assemble) with the metadata and environment commands, and the HTML pages
PREPARE writes (`explore_containment`, `sra_report`, `project_report`, and the plotly pages of
`interactive_plot`), at commit fb922e2 (v0.11.1). The SRA download internals were audited on
2026-10-21 (`2026-10-21-sra-download-audit.md`) and are out of scope.

Method: `make pipeline` (the bundled Salmonella walkthrough: 10 checks, all passed); two scratch
projects (the synthetic two-dataset project of the `sra_sketch` live check, and a three-genome
project built for this audit) exercised with the CLI and no store; the pages generated from those
projects opened in Chrome over a local HTTP server (filter, sort, console); two independent
read-only code reviews, one of the stage hand-offs and one of the page writers, with synthetic
injection payloads rendered in headless Chrome. "Reproduced" below means the CLI or the function
was run; "read" means source evidence only. Nothing in the repository was changed.

Severity: H = wrong result recorded as correct, data loss, a crash on realistic input, or a
chain that cannot run as documented; M = wrong behaviour without an error, a misleading number
or label, or a documented option that is broken; L = cosmetic, drift, duplication.

## 1. Pipeline findings

### High

| # | Where | Defect | Failure scenario |
|---|---|---|---|
| P1 | `cli/commands/sra_sketch.py:296-312`, `sra/sketch_containment.py:74-78` | With `--genome-fasta`, every accession `_plan` marks failed is passed as `drop=` to `merge_containment_csv`, and an accession with no FASTQ files on disk is "failed". Its row is deleted from `matches/<genome>.csv`. The default genome id is the FASTA stem, the same stem `branchwater_search` writes, so the file edited is the Branchwater result, and the Branchwater containment of the sketched accessions is replaced by the local value as well. | `sra_sketch --accessions-file accessions.txt --genome-fasta genomes/G.fna` before every accession is downloaded: the Branchwater rows of the undownloaded accessions vanish; the next `parse_containment` and `select_datasets` lose them with exit 0. `docs/reads.md` says rows of other accessions are kept. Reproduced. |

### Medium

| # | Where | Defect | Failure scenario |
|---|---|---|---|
| P2 | `data/read_extraction.py:592`, `cli/commands/read_extraction.py:570` | `extract_target_reads` chooses samples from the containment table at `--threshold` (default 0.1) and reads neither the selection nor the exclusions. | After `select_datasets` at 0.5 and a `blacklist`, extraction maps and assembles the blacklisted dataset and every downloaded dataset between 0.1 and 0.5, all recorded as normal extractions. Reproduced with a dry run. |
| P3 | `data/read_extraction.py:227`; `cli/commands/status/suggest.py:212` | The containment threshold is part of the extraction identity although it only chooses samples, and `status --next` prints its `--assemble` suggestion without `--threshold`. | A project extracted at 0.5 that runs the suggested command re-maps every sample before assembling. Read. |
| P4 | `cli/commands/status/suggest.py:175, 212-218` | `status --next` suggests extraction for every genome and every downloaded accession, whether or not the accession was screened against that genome, and pastes the genome id, FASTA path and table without `shlex.quote`. | A genome with containment 0 is suggested; the run exits 1 with "No FASTQ files found". A genome id with a space (a web CSV named `L crispatus.csv`) gives a broken command. Reproduced. |
| P5 | `cli/commands/containment.py:86-94`, `data/branchwater.py:528-534` | With one unreadable match file, `parse_containment` still writes the table from the rest, records screening from it, and only then exits 1. A `;` or TAB separated `.csv` reports "Missing required columns: acc, containment". | The table loses that genome's column; `select_datasets` on `max_containment` drops those accessions silently while the registry still holds their screening. Reproduced. |
| P6 | `data/branchwater.py:508`, `cli/commands/containment.py:86` | `parse_containment` on an empty or misspelt matches folder writes nothing, exits 0, then records screening from the old table with `csv` provenance pointing at the empty folder. | A typo in `--matches-folder` leaves provenance for files that do not exist. Reproduced. |
| P7 | `data/branchwater.py:125` | `use_branchwater` skips any existing `matches/<genome>.csv` at INFO, exit 0; there is no `--force`. | After `branchwater_search --refresh` or a newer web download, `matches/` keeps the old result; after P1 has edited the file it is never restored. Reproduced. |
| P8 | `cli/commands/sra_enhanced.py:217-219, 468`; `data/registry.py:609` | `sra_validate`: accessions named by `--accessions-file` with no folder are dropped silently with exit 0; a folder scan records stray folders (`fastq/old_runs` became a dataset); a FAILED validation counts as "analysed" in `status` and the funnel. | 1 of 2 named accessions validated, exit 0; `old_runs` listed under `status --stage analysed`. Reproduced. |
| P9 | `data/metadata.py:89, 212-216, 497, 506` | `download_metadata` logs each accession "not in the NCBI response" as ERROR and exits 0 even when every one fails; an unreadable matches CSV is a warning there but exit 1 in the three other readers; `parse_metadata` exits 0 when XML files fail to parse and leaves an older table in use. | An accession with no metadata never changes the exit status. Read (network). |
| P10 | `data/read_extraction_export.py:89-93` | Paired extraction returns only `_1` and `_2`; `<genome>_s.fastq.gz` (reads whose mate was filtered) is on disk but absent from `files` and never given to megahit, while `mapped_reads` counts those reads and `status --init` lists `_s`. | Mapped reads in divergent regions are missing from the assembly and the file record. Read; confirm whether intended. |
| P11 | `data/branchwater_search.py:73, 80`; `cli/commands/branchwater_search.py:75`; `cli/commands/genome.py:29` | `sketch_fasta` opens the FASTA as plain text, and the genome id is `Path.stem`, so `G.fna.gz` gives `G.fna` where `genome_prepare` strips `.fna.gz`. `GENOME_FASTA_GLOBS` and `status` accept `.fna.gz`. | `branchwater_search` and `sra_sketch` exit 1 with a bare `UnicodeDecodeError`; with a working read the registry would hold two ids for one genome. Reproduced. |
| P12 | `data/read_extraction.py:104` | `_notify_result` catches every exception, including `LockTimeoutError`, as a warning; the run exits 0 where a `TransientError` should give 4. | Reads on disk with no registry record; `status` shows the sample as not extracted and the next run maps it again. Read. |
| P13 | `data/registry_scan.py:258-270`, `cli/commands/project_report.py` | `status --init` registers every FASTA under `genomes/` as a project genome. A genome that no stage used appears in `project_report` ("Datasets per genome", the genome table) beside the real one. | The walkthrough project reports two genomes (`GCF_000008985.1` from the folder, `salmonella_subset` from the matches). Reproduced. |

### Low

| # | Where | Defect |
|---|---|---|
| P14 | `data/screening_table.py:137` | `parse_containment` records `cani=None`, overwriting the value `branchwater_search` stored; `status --init` records it from the CSV. Nothing reads the field; `pipeline_overview.md` says the registry keeps cANI. Reproduced. |
| P15 | several commands | Path anchoring differs between producers and consumers: `fastq`, `matches`, `metadata`, `targeted` and sketch folders are anchored to the project root; `--parsed-containment-file`, `--summary-containment-file`, `select_datasets --parsed-containment/--output`, `extract_target_reads --parsed-containment`, `results_table`, `parse_metadata --metadata-table-file`, `--branchwater-folder`, the `sra_profile`/`sra_report` output folders, `--blacklist-file` and `status --accessions-file/--parsed-containment` are relative to the working directory. From `proj/sub`, `parse_containment` reads `proj/matches` and writes `sub/parsed_containment.txt`. Read. |
| P16 | `cli/commands/sra_sketch.py:291, 311` | `--genome-id` is joined into `matches/<id>.csv` unvalidated: `a/b` writes a subfolder `parse_containment` never lists; `../x` escapes the folder. Read. |
| P17 | `use_branchwater`, `extract_branchwater_metadata`, `parse_containment`, `select_datasets --no-record`, `sra_validate`, `status`, `explore.py:146-166`, `advanced_analysis.py:254` | `logger.error(...); return 1` instead of `self.fail`; `sra_validate`'s generic `except Exception` logs no traceback; `sra_report` exits 0 when some accessions fail where `sra_profile` exits 1. Read. |
| P18 | `data/branchwater.py:197-205` | `extract_branchwater_metadata` splits the header on `,` by hand; a file in another format is skipped with a warning, not an error, exit 0. Read. |
| P19 | docs | `extraction_assembly.md:22` says `--data-root` resolves `--fastq-folder` (the store is only used for usage rows); `project_state.md:133` omits `sra_sketch` from the run-log commands; `project_state.md:57` says presence is checked against `genomes/<GCF>.fna` (only FASTQ and metadata are checked); the extraction skip key also compares `min_mapq` and the resolved FASTA path; `status --next` prints extraction and assembly suggestions the docs do not mention; `reads.md:184` omits that `--quality-profiles` alone works; the `_0` and `_s` extraction outputs are undescribed; `sra_report --groups-file` must be JSON, which the flag help does not say (a TSV is refused with a clear message). |
| P20 | `cli/base.py` `read_accessions_file` | `--accessions-file <(printf ...)` (a process substitution) fails; a regular file is required. Cosmetic. |

## 2. HTML page findings

Injection: none. All three jinja2 environments autoescape (`explorer.py:391`, `reporting.py:59`,
`visualization/project_report.py:85`); table cells use `html.escape`; plotly escapes `</` in
figure JSON and renders labels as SVG text. Payloads placed in accessions, genome and sample
ids, taxon names, metadata, titles, group names, the project name, a failed-download message and
PCA labels executed on none of six rendered pages. No page references a CDN; plotly.js is
embedded (4.8 MB per page), so every page works offline. No console errors in Chrome.

### High

| # | Where | Defect | Failure scenario |
|---|---|---|---|
| G1 | `cli/commands/sra_report.py:313`, `sra/analytics.py:682-697`, `sra/reporting.py:443` | A group with no profiled accession is kept as an empty list; `compare_datasets` then sets every test to `{}`; the template formats `test_result.p_value` and raises `UndefinedError`. `sra_report.json` is written, the HTML is not, exit 1 from the generic handler. | Groups file `{"A": [...], "B": ["SRRx"]}` where `SRRx` has no FASTQ. Reproduced: "Report failed: 'dict object' has no attribute 'p_value'". |
| G2 | `sra/reporting.py:326, 389`; `store/stats.py:187`; `sra/analytics.py:138-140` | The "Total reads" card sums `profile.total_reads`: both mate files are counted (the terminal and `reads.md` say "mates counted", the page does not), and a profile that only has sampled reads contributes its sample size with no mark. The "Total Reads by group" box plot has the same issue. | Profiles 10,000 (sampled) + 2M + 3M + 4M render "Total reads 9,010,000". Reproduced. |
| G3 | `visualization/explorer.py:188-200` | The sunburst's `values="sample_count"` is the number of samples per species; plotly sums children into parents, so a sample holding two species of a family counts twice at family level, contradicting the "Top families by sample count" card on the same page. | Two samples, two genomes in one family: card "F (2)", sunburst F = 4. Reproduced. |

### Medium

| # | Where | Defect |
|---|---|---|
| G4 | `sra/reporting.py:441-444`, `sra/analytics.py:701-716` | Skipped tests print "p-value nan · not significant"; "skipped (identical values)" prints "p-value 1.0000", never computed; the explanatory `note` is not rendered. Seen on the synthetic project. |
| G5 | `sra/analytics.py:725, 771`, `sra/reporting.py:440` | Four tests at p < 0.05 each with no multiple-testing correction, and the recommendations say "Significant differences detected in ...". |
| G6 | `visualization/explorer.py:150-245` | Without a taxonomy map the page is empty of information: cards "families 0, genera 0" beside "Top families: Unknown (83)"; a sunburst of one disc "Unknown / Unknown / Unknown"; a heatmap of 83 rows by one "Unknown" column; "Containment by family: Unknown". The page should say that no taxonomy is present (and name `enrich_taxonomy`) and hide the taxonomy sections. Seen on both projects. |
| G7 | `visualization/explorer.py:241` | "Samples per genus" keeps the top 30 genera; the title does not say so (the table cut-off has a note). |
| G8 | `visualization/explorer.py:211` | The heatmap aggregates with `max` per family but is titled "Sample-by-family containment": a family with several genomes shows only its best genome. Reproduced. |
| G9 | `visualization/explorer.py:255-277` | Metadata placeholders from Branchwater are shown raw in the Location column: "NP" and "uncalculated". Seen on the walkthrough project. |
| G10 | `visualization/explorer.py:158-171` | With no row above the threshold the page says "0 samples 0 genomes ... containment 0.000-0.000" with no explanation; `total_genomes` is reset to 0 in the empty branch while the non-empty branch counts genome columns with no hit. |
| G11 | `cli/commands/advanced_analysis.py:194-254`, `visualization/interactive.py:209-461` | `interactive_plot` without `--output-file` and with `--no-show` writes nothing and prints "Interactive plot created successfully!" with exit 0; without `--no-show` it calls `fig.show()`, which needs a browser; errors use `logger.error; return 1`. |
| G12 | `visualization/interactive.py:294, 364, 374-375` | The interactive heatmap uses the diverging `RdBu_r` scale for non-negative data; `sample_metadata` is accepted and ignored; the figure grows 15 px per sample and 20 px per feature with no cap. |
| G13 | `visualization/project_report.py:35-56`, `processing/project_report.py` | The funnel chart draws the stages as nested when they are not: the synthetic project shows screened 1, analysed 2. The "Datasets per genome" bar chart with all zeros renders an axis from -1 to 1. Seen on both projects. |

### Low

| # | Where | Defect |
|---|---|---|
| G14 | `visualization/interactive.py:56` | The pca/tsne/heatmap/diversity pages have no `<title>`, viewport or `lang`, a fixed 800 px width and the default plotly template, unlike the three report pages. |
| G15 | `sra/reporting.py:370` | `<html>` without `lang`, no skip link; the explorer and project report have both. |
| G16 | `explorer.py:376`, `reporting.py:227, 470`, `project_report_markdown.py:204`, `visualization/project_report.py:89` | Stamps differ: the explorer has no date or version; `sra_report` has local time without zone and no version; `project_report` has UTC plus version but local `created`/`updated` cells; its `<title>` is the constant "Project report"; `explore_containment` has no `--title`. |
| G17 | `project_report_markdown.py:163` | No links between pages; output paths are plain text and only the first output per analysis is listed. |
| G18 | `explorer.py:269`, `utils/html.py:315, 346-367` | The table clamps containment to [0, 1] without a mark; CSV export covers only the displayed rows (at most 5000) under a fixed name; the numeric sort regex ignores scientific notation (harmless at `%.3f`). |
| G19 | `explorer.py:286-381`, `reporting.py:368-474`, `visualization/project_report.py:100-155` | Three templates repeat doctype, head, plotly, CSS, header and footer; `project_report` adds an inline style to undo table-sort styling. `SRAReportGenerator.generate_quality_dashboard` and `create_comparative_analysis` (`reporting.py:70-184`) have no caller outside tests. |
| G20 | `reporting.py:313-314, 441` | Box-plot titles and test headings use `str.title()` ("Gc Percent By Group", "Avg Read Length") with no units. |
| G21 | accessibility | Charts have no text alternative; explorer header cells are `th role=button` without `scope`; interactive pages have no dark handling. |

## 3. Test gaps

Covered: explorer title and family escaping, table cap and note, heatmap height cap
(`tests/test_explorer.py`); `project_report` escapes a failure message; inline plotly without a
CDN and the CSV formula guard (`tests/test_html_assets.py`); every flag named in the four user
guides exists in the parsers (checked programmatically); older registries load.

Not covered: P1 (a failed accession's Branchwater row), P2 (extraction against the selection and
exclusions), P5 and P6 (partial or empty `parse_containment`), P7 (a changed Branchwater CSV),
P8 (named accessions dropped by `sra_validate`), P11 (a gzipped FASTA), P12 (a lock timeout in
`_notify_result`); escaping in `sra_report`; G1 (the empty group); G2 (sampled or mates-counted
totals); G3 (sunburst aggregation); G4 (p-value display); G7 (genus cut); G10 (empty explorer);
the interactive pages (their tests mock `write_html`, so the written page is never checked).

## 4. Checked and sound

- Matches CSVs with extra columns or a UTF-8 BOM parse; duplicate accessions keep their maximum.
- The metadata table has `Run_ID` first, tab separated, with the columns `select_datasets` needs;
  a missing filter column is an error.
- `download_sra` does not overwrite an earlier `downloaded` state for blacklisted or capped
  accessions; `results_table` and `status --export-tsv` read quality through one function.
- A newer registry schema is refused; missing extraction parameters act as wildcards; unequal
  mate counts switch extraction to per-mate single-end mapping; an extraction rerun with
  unchanged parameters skips samples and still assembles them; a failed assembly exits 1.
- The registry lock is never held during ranking or profiling.
- Explorer: NaN containment is dropped, `--min-containment` is inclusive, metadata is matched by
  sample name with duplicates tolerated, non-ASCII and 300-character names render, a 3000 by 20
  table renders in about 3 s with the "showing top 5000 of 60000 rows" note; filter and column
  sorting work in Chrome. GC is in percent throughout `sra_report`. Every writer uses the atomic
  helpers.

## 5. Suggested order

1. P1 and P7 together (drop only accessions that failed while sketching, never those without
   files; a `--force` for `use_branchwater`): one small plan, since P1 can delete screening data
   and P7 prevents its recovery.
2. G1, G2, G3, G4 and G6: the `sra_report` crash and the three numbers a reader would act on,
   plus the empty-taxonomy explorer. Each is a failing test and a small fix.
3. P2, P3, P4: extraction respects the selection and exclusions, the threshold leaves the
   extraction identity, `status --next` quotes and filters its suggestions. One plan; P2 changes
   documented behaviour and needs a decision on whether the selection or the threshold wins.
4. P5, P6, P8, P9, P12, P17: exit codes and partial writes in the screen, metadata and validate
   stages, which share one shape (write only when every input was read; fail through
   `self.fail`).
5. The rest as a documentation and polish round (P13-P16, P18-P20, G7-G21), including one shared
   page template.
