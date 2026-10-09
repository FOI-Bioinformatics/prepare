# SRA information, profiles, reports and validation

Four commands support the `download_sra` step. They read an accession list (`--accessions-file`,
one per line) or, for `sra_profile` and `sra_validate`, repeated `--accession` flags; none of them
download reads.

## sra_info

Fetches NCBI metadata for the accessions and prints sizes, layouts, technologies and an estimated
download time. Requires `--email` (NCBI asks for it); `--api-key` raises the rate limit.

```bash
prepare sra_info --accessions-file accessions.txt --email you@example.org --output-report sra_info_report.csv
```

## sra_profile

Reads the downloaded FASTQ files under `--fastq-folder` (default `fastq`) and writes, per
accession, read and base totals, read-length statistics, N50, GC content (in percent, under
`gc_percent`), mean base quality, complexity, duplication and a quality grade: one row each in
`--output-report` (default `sra_statistics.csv`) and one `<accession>_quality_profile.json` each in
`--output-dir` (default `sra_quality_profiles`). Totals and GC come from the dataset's statistics
record, cached in the store sidecar for a store-linked dataset; the per-read figures come from a
sample of `--sample-size` reads drawn from every mate file. GC content is computed from a sample of the first mate file; per-read quality from a sample of all mates.

```bash
prepare sra_profile --fastq-folder fastq --accession SRR123456 --accession SRR123457
```

## sra_report

Writes one HTML report, `sra_report.html` in `--output-dir` (default `sra_reports`), and its figures
in `sra_report.json`. With `--groups-file` (a JSON object mapping group names to accession lists) the
report adds a comparison of the groups with a t-test (two groups) or one-way ANOVA (more). Each
accession is profiled once, or read back from `--quality-profiles`; `--no-report` skips the HTML.

```bash
prepare sra_report --groups-file groups.json --quality-profiles sra_quality_profiles
```

## sra_validate

Checks that each accession folder holds readable FASTQ files; `--check-pairs` also verifies that
paired files have matching read counts. `--accessions-file` and repeated `--accession` flags restrict
the check to the accessions named (the `--accessions A B` form was removed in 0.5.0).

```bash
prepare sra_validate --fastq-folder fastq --check-pairs --accession SRR123456
```

## Renamed in 0.5.0

`sra_stats` and `sra_profile_quality` became `sra_profile`; `sra_dashboard` and `sra_compare` became
`sra_report`. The old names print the new one and exit with status 2.

Downloads themselves are documented in the README under "Downloading reads".
