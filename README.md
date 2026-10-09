# PREPARE

PREPARE is a command-line toolkit for selecting and analysing public metagenomic datasets by genome
containment. It searches the [Branchwater](https://branchwater.sourmash.bio/) index (a sourmash
project) for SRA runs that contain a target genome, collects their NCBI metadata, downloads the reads,
and extracts and assembles the reads that map to the target. Each step records its outcome in a
project registry, so a project can be inspected, resumed and reported on at any point.

The pipeline has six stages: screen, select, download, analyse, extract and assemble.
[docs/pipeline_overview.md](docs/pipeline_overview.md) describes each stage and the commands it uses.

## Installation

Requires Python 3.12 or newer.

### Quick start

```bash
git clone https://github.com/FOI-Bioinformatics/prepare.git
cd prepare
make dev-install  # every extra plus the development tools
```

### Core install and optional extras

A plain install brings only the core packages (pandas, numpy, matplotlib, biopython, lxml,
requests). That is enough to run the download, containment, metadata and matplotlib plotting steps.
Other features need an extra:

| Extra | Packages | Needed for |
|---|---|---|
| `analysis` | scikit-learn, scipy | `diversity_analysis`, the PCA, t-SNE, diversity and clustered heatmap plots of `interactive_plot`, group statistics in `sra_report --groups-file` |
| `interactive` | plotly, jinja2 | `interactive_plot`, `explore_containment`, the `sra_report` HTML report, the `project_report` HTML page |
| `maps` | cartopy | geographic sample maps |
| `sourmash` | sourmash | sketching a genome in `branchwater_search`, and the `sourmash scripts prepare_*` plugin |
| `all` | all of the above | everything |

```bash
pip install .                      # core only
pip install '.[analysis,interactive]'
pip install '.[all]'               # every extra
```

A command that needs a missing extra stops with exit code 3 and an error that names the extra and
the interpreter to install it into, for example:

```text
Diversity analysis needs the 'scikit-learn' package. Install it into this interpreter with: /path/to/python -m pip install 'prepare[analysis]'
```

`requirements.txt` lists the core packages and `requirements-all.txt` the core
plus every extra.

### Conda environment and external tools

`environment.yml` installs PREPARE with every extra and the external tools below, except the
optional `seqkit`:

```bash
conda env create -f environment.yml
conda activate prepare
```

Download, extraction and assembly steps call these command-line tools:

| Tool (conda package) | Oldest version | Used by |
|---|---|---|
| `fasterq-dump` (sra-tools) | 3.0 | `download_sra` |
| `prefetch` (sra-tools, optional) | 3.0 | `download_sra`; without it fasterq-dump reads directly from NCBI |
| `pigz` (optional) | 2.4 | `download_sra`, `store_adopt`; without it files are compressed with gzip |
| `datasets` (ncbi-datasets-cli) | 16 | `genome_download`, `genome_prepare` |
| `minimap2` | 2.17 (first with `--sam-hit-only`) | `extract_target_reads` |
| `samtools` | 1.10 (first with `samtools coverage`) | `extract_target_reads` |
| `megahit` (optional) | 1.2.9 | `extract_target_reads --assemble` |
| `seqkit` (optional) | any | `sra_profile`, `sra_report`; without it read statistics are computed in Python |

A command checks the tools it needs before it starts. A missing tool, one that cannot run, or one
older than the version above stops the command with exit code 3 and the `conda install` command that
provides it. The optional tools (`prefetch`, `pigz`, `seqkit`) are not checked by the commands, which
work without them. `download_test_genome` fetches its genome over HTTPS and needs none of these
tools. The table is defined in `prepare/utils/tools.py`, and `environment.yml` pins the same
versions.

## Checking the environment

`prepare doctor` reports whether this environment can run PREPARE, without changing anything; its
flags and checks are in [docs/configuration.md](docs/configuration.md#checking-the-environment).

## Configuration

Run-time settings are taken from a flag, then a `PREPARE_<NAME>` environment variable, then the
`[runtime]` table of `~/.config/prepare/config.toml`, then a built-in default. The shared data store
root, logging and every setting with its type, default, variable and flag are described in
[docs/configuration.md](docs/configuration.md).

## Commands

`prepare --help` lists the commands in seven groups, and `prepare COMMAND --help` lists the options
of one command.

| Group | Command | Purpose |
|---|---|---|
| Containment | `branchwater_search` | Search the Branchwater index with a genome and write a Branchwater CSV |
| | `use_branchwater` | Process Branchwater CSV files downloaded beforehand |
| | `parse_containment` | Parse containment data from match files |
| | `plot_containment` | Plot containment data |
| | `explore_containment` | Write an interactive HTML explorer for containment data |
| | `enrich_taxonomy` | Add GTDB taxonomy to genome accessions |
| | `find_by_taxonomy` | Filter containment results by taxonomic classification |
| Metadata | `extract_branchwater_metadata` | Extract metadata from Branchwater files |
| | `download_metadata` | Download metadata for SRA accessions |
| | `parse_metadata` | Parse downloaded metadata files into one table |
| | `check_metadata_attributes` | Count how often each metadata attribute is populated |
| | `count_metadata` | Count metadata values by genome |
| | `single_sample` | Analyse a single sample |
| | `plot_metadata_counts` | Plot metadata counts |
| Genomes | `genome_search` | Search GTDB for genome accessions by species or genus |
| | `genome_download` | Download genome assemblies from NCBI by accession, species or genus |
| | `genome_prepare` | Search GTDB, download genomes and write a manifest CSV |
| | `download_test_genome` | Download the test genome |
| Reads | `select_datasets` | Write an accessions file from parsed containment |
| | `blacklist` | Record datasets to exclude, with a reason |
| | `download_sra` | Download SRA datasets |
| | `status` | Report which reads, metadata and genomes are available locally |
| | `results_table` | Write one results row per accession and genome |
| | `sra_info` | Report information on SRA datasets before downloading |
| | `sra_validate` | Validate downloaded SRA datasets |
| | `sra_profile` | Statistics and a quality profile for each downloaded dataset |
| | `sra_report` | One HTML report on dataset quality, with an optional group comparison |
| | `extract_target_reads` | Keep the reads that map to a target genome, optionally assemble them |
| Store | `store_init` | Create the shared data store and record this project's use of it |
| | `store_status` | Report dataset counts, bytes and stale projects |
| | `store_reindex` | Rebuild the store catalogue from the dataset sidecar files |
| | `store_adopt` | Move or copy project FASTQ folders into the store and link them back |
| | `store_verify` | Verify store datasets against their sidecars |
| | `store_link` | Link project accessions to the store's copies |
| | `store_unlink` | Remove a project's link to the store |
| | `store_usage` | Report which projects and genomes used a dataset |
| | `store_gc` | Report, and with `--yes` remove, unused datasets and leftover temporary files |
| Analysis | `diversity_analysis` | Alpha and beta diversity metrics |
| | `interactive_plot` | Interactive visualizations |
| | `validate_taxonomy` | Validate species names against NCBI taxonomy |
| | `taxonomic_summary` | Taxonomic summaries at several levels |
| Environment | `doctor` | Check tools, configuration, store, disk space and resources |
| | `runs` | List, show and compare recorded runs of this project |
| | `project_report` | Write a Markdown, JSON and HTML report of the project |

One example per group follows. The linked documents describe every option and output.

**Containment** ([docs/pipeline_overview.md](docs/pipeline_overview.md)):

```bash
prepare branchwater_search --genome-fasta genomes/GCF_000008025.1.fna --threshold 0.1
prepare use_branchwater --branchwater-folder branchwater --matches-folder matches
prepare parse_containment --matches-folder matches --parsed-containment-file parsed_containment.txt
```

**Metadata** ([docs/pipeline_overview.md](docs/pipeline_overview.md)):

```bash
prepare download_metadata --matches-folder matches --metadata-folder metadata --threshold 0.95 --email [EMAIL]
prepare parse_metadata --metadata-folder metadata --metadata-table-file metadata_table.txt
```

**Genomes** ([docs/pipeline_overview.md](docs/pipeline_overview.md)):

```bash
prepare genome_prepare --species "Lactobacillus crispatus" --output-dir genomes
```

**Reads** ([docs/reads.md](docs/reads.md), [docs/project_state.md](docs/project_state.md),
[docs/extraction_assembly.md](docs/extraction_assembly.md)):

```bash
prepare select_datasets --threshold 0.9 --output accessions.txt
prepare download_sra --accessions-file accessions.txt --fastq-folder fastq --max-workers 4 --num-threads 4
prepare status
prepare extract_target_reads --parsed-containment parsed_containment.txt --genome-id GCF_000006945.2 \
  --genome-fasta genomes/GCF_000006945.2.fna --fastq-folder fastq --output-folder targeted --threshold 0.5
```

**Store** ([docs/shared_store.md](docs/shared_store.md)):

```bash
prepare store_init --data-root /data/prepare_store --set-default
prepare store_adopt --fastq-folder fastq --move
```

**Analysis** ([docs/analysis_commands.md](docs/analysis_commands.md)):

```bash
prepare enrich_taxonomy --parsed-containment parsed_containment.txt --output taxonomy.tsv
prepare taxonomic_summary --abundance-file parsed_containment.txt --taxonomy-file taxonomy.tsv \
  --levels phylum class order family genus --output-dir taxonomic_summaries
```

**Environment** ([docs/project_state.md](docs/project_state.md)):

```bash
prepare runs --diff previous latest
prepare project_report
```

## Exit codes

Every command exits with 0, 1, 2, 3, 4 or 130, so a batch script can decide whether to resubmit.
[docs/hpc.md](docs/hpc.md#exit-codes-and-job-dependencies) gives the meaning of each code.

## Documentation

[docs/README.md](docs/README.md) is the index of the documentation.

## Contributing and releases

See [CONTRIBUTING.md](CONTRIBUTING.md) for the development setup, quality checks and pull requests.
[CHANGELOG.md](CHANGELOG.md) lists the changes in each release, and
[docs/packaging.md](docs/packaging.md) describes the release workflow.

The nightly smoke test against the real NCBI and SRA services is described in
[CONTRIBUTING.md](CONTRIBUTING.md#network-smoke-test).

## License

PREPARE is distributed under the GNU General Public License v3.0; see [LICENSE](LICENSE).
