# Advanced Branchwater Workflows

This guide covers advanced workflows and decision points for using PREPARE with Branchwater containment data. For basic usage, see the main [README](../README.md).

## Workflow Decision Points

### Metadata Strategy Selection

Choose your metadata approach based on analysis requirements:

#### Option A: Embedded Metadata (Fast)
Use metadata already present in Branchwater CSV files. Suitable for rapid analysis with basic sample information:

```bash
prepare extract_branchwater_metadata --branchwater-folder branchwater_files --metadata-folder metadata
```

**Advantages:**
- No NCBI API calls required
- Immediate availability
- Includes: organism, bioproject, biosample, geographic data

**Limitations:** 
- Limited depth compared to full NCBI records
- May lack study-specific metadata

#### Option B: Full NCBI Metadata (Comprehensive)
Download complete metadata records from NCBI for detailed analysis:

```bash
prepare download_metadata --matches-folder matches --metadata-folder metadata --threshold 0.95 --email your@email.com
prepare parse_metadata --metadata-folder metadata --metadata-table-file parsed_metadata.txt
```

**Advantages:**
- Complete sample descriptions
- Study protocols and methods
- Enhanced publication metadata

**Considerations:**
- Requires NCBI API access
- Longer processing time
- Rate-limited downloads

## Searching from a FASTA

`branchwater_search` replaces the manual download. It needs the sourmash extra and network access:

```bash
prepare branchwater_search --genome-fasta genomes/GCF_000008025.1.fna --threshold 0.1 --branchwater-folder branchwater
prepare use_branchwater --branchwater-folder branchwater --matches-folder matches
```

A run that returns zero matches is reported as a warning, not an error; verify the index with a genome
known to be abundant in metagenomes before trusting an empty result.

### Retries and the response cache

A `branchwater_search` request that fails with a network or server error retries automatically, up to 4
attempts with a growing backoff between them; a 4xx client error other than a rate limit fails at once
instead of retrying. A successful response is cached under `.branchwater-cache/` next to the output CSV,
keyed by the genome, containment threshold, sourmash parameters and server URL, so repeating the same
search reads the cache instead of querying again:

```bash
prepare branchwater_search --genome-fasta genomes/GCF_000008025.1.fna --threshold 0.1   # caches
prepare branchwater_search --genome-fasta genomes/GCF_000008025.1.fna --threshold 0.1   # reads the cache
prepare branchwater_search --genome-fasta genomes/GCF_000008025.1.fna --threshold 0.1 --refresh  # forces a new query
prepare branchwater_search --genome-fasta genomes/GCF_000008025.1.fna --threshold 0.1 --no-cache  # never reads or writes the cache
```

`--max-cache-age-days N` treats a cached entry older than N days as absent. There is no default
expiry, so a cache entry is otherwise reused until `--refresh` or `--no-cache` is used, or the cache
file is deleted by hand.

## Advanced Filtering and Thresholds

### Containment steps

`parse_containment` does not filter; it records every sample and summarizes how many samples reach
each containment step. Choose the step size for the summary and apply thresholds downstream:

```bash
prepare parse_containment --matches-folder matches --step-size 0.05
prepare select_datasets --threshold 0.9 --output accessions.txt
prepare count_metadata --metadata-column Sample_Scientific_Name --threshold 0.9
```

### The details table

`parse_containment` also writes `parsed_containment_details.tsv`, one row per (accession, genome) pair
that appears in a match file, with the columns that `parsed_containment.txt` does not carry:

| Column | Meaning |
|---|---|
| `accession` | SRA run accession |
| `genome_id` | target genome the row's containment refers to |
| `containment` | containment value for this accession and genome |
| `cANI` | Branchwater's estimated average nucleotide identity for the match |
| `biosample`, `bioproject` | NCBI sample and project identifiers, when Branchwater reports them |
| `assay_type`, `organism` | sequencing assay and the organism label from the metagenome record |
| `geo_loc_name`, `lat_lon` | geographic label and coordinates, when present |

Use this table when a downstream step needs cANI or sample provenance; `parsed_containment.txt` stays
the samples-by-genomes matrix the rest of the pipeline reads.

## What the registry records

Most commands in this workflow update the project registry (`prepare_registry.json`) as they run;
`status` reads it back. One line per command:

- `branchwater_search` / `parse_containment`: screening, the containment per genome
- `select_datasets`: selection, the threshold and metadata filter used
- `blacklist`: exclusions, with a reason
- `download_sra`: download outcomes, with file sizes and dates
- `download_metadata` / `parse_metadata`: metadata
- `sra_validate` / `sra_profile` / `sra_report`: analyses (`validate`, `profile`, `report`)
- `extract_target_reads`: extraction and assembly, per target genome

See the README's "Project state" section and `prepare status --help` for the full set of flags.

## Troubleshooting Common Issues

### Format Validation Errors
If you encounter format errors, verify your CSV structure:

```bash
# Check file headers
head -1 your_branchwater_file.csv

# Expected format: acc,containment,organism,[additional_columns]
```

---

*For basic Branchwater usage, see the main [README](../README.md). For technical details, refer to [ARCHITECTURE.md](ARCHITECTURE.md).*