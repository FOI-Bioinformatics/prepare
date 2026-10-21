# Analysis and plotting commands

This guide describes the commands that summarise and plot containment and abundance tables, and the input format that each one expects. Return to the [documentation index](README.md).

## Input formats

- **Parsed containment table.** The output of `parse_containment` is a tab-separated table with accessions as rows and one containment column per genome. The first column is the accession and becomes the index. Commands that read it also accept a `max_containment` column and ignore non-numeric columns such as `max_containment_annotation`.
- **Abundance matrix.** `diversity_analysis`, `interactive_plot` and `taxonomic_summary` expect a numeric table with samples as rows and species or genomes as columns. The first column holds the sample names. A file with a `.txt` or `.tsv` suffix is read as tab-separated, any other as comma-separated. Non-numeric columns are dropped. The parsed containment table has this shape, so it can be used as an abundance matrix. No earlier pipeline step writes a file named `abundance_matrix.csv`; supply your own table in this layout.
- **Sample metadata table.** A table with the same first column as the matrix (the sample names) and one column per attribute, for example `treatment` or `site`. Text columns are kept. Samples in the matrix must be present in the metadata.

## Plotting containment

`plot_containment` plots the distribution of containment scores and always writes an image file next to the input table, named `<input-stem>_<plot-type>_<column>.<format>`. The default format is png, and `--save-format` selects png, jpg, pdf or svg.

```bash
prepare plot_containment --file-path parsed_containment.txt --column max_containment --plot-type rank --save-format png --threshold 0.05
```

This writes `parsed_containment_rank_max_containment.png`. The plot types are `rank`, `histogram`, `box` and `violin`. The sourmash plugin command `sourmash scripts prepare_plot`, registered by the `sourmash` extra, calls the same function and writes files with the same names.

## Plotting metadata counts

`plot_metadata_counts` plots the distribution of metadata attributes. Unlike `plot_containment`, it writes a file only when `--save-format` is given. Without it the plot is built and not saved.

```bash
prepare plot_metadata_counts --file-path metadata_counts.txt --plot-type bar --save-format png
```

This writes `metadata_counts_bar.png` in the current working directory. The plot types are `bar`, `pie` and `radar`. The input is a count table from `count_metadata`.

## Plotting sample locations

`plot_map` draws the sample locations of a metadata table on a world map. It reads the coordinates from one column, `lat_lon` by default, and accepts values such as `38.98 N 77.10 W` or `12.5,-3.1`. Rows without parseable coordinates are not drawn. It needs the `maps` extra.

```bash
prepare plot_map --metadata-file metadata_table.txt --output sample_map.png --lat-lon-column lat_lon
```

The image format follows the extension of `--output`: png, jpg, pdf or svg.

## diversity_analysis

`diversity_analysis` calculates alpha and beta diversity metrics for an abundance matrix and, with a metadata table, a PERMANOVA test. It needs the `analysis` extra.

```bash
prepare diversity_analysis \
    --abundance-file abundance_matrix.csv \
    --metadata-file sample_metadata.csv \
    --alpha-metrics shannon simpson chao1 \
    --beta-metric bray_curtis \
    --permanova-formula "treatment + site"
```

The alpha metrics default to `shannon`, `simpson`, `chao1` and `observed_species`. `--beta-metric` is one of `bray_curtis` (the default), `jaccard`, `euclidean` and `manhattan`. The terms of `--permanova-formula` are columns of the metadata table. Results go to `--output-dir` (default `diversity_results`).

Each term of `--permanova-formula` is tested on its own as a one-factor PERMANOVA (Anderson 2001) on the beta-diversity distance matrix, with 999 permutations. Samples with a missing value for the term are left out of that test. The pseudo-F statistic uses the standard sums of squares, SS_T = sum of d_ij^2 / N over all pairs and SS_W = sum over groups of the within-group pairs' d_ij^2 / n_g, and the p-value is (1 + the number of permutations with F >= the observed F) / (1 + 999). F and p are reported as nan when fewer than two groups remain. Since 0.10.1 these values differ from those of 0.10.0 and earlier, which used a different sum of squares and no +1 correction; results from earlier versions should be recomputed. The command has no seed option, so the p-value can vary slightly between runs.

## interactive_plot

`interactive_plot` writes a browser-based plot from a data matrix. It needs the `interactive` extra and the `analysis` extra (scikit-learn for `pca`, `tsne` and `diversity`, scipy for `heatmap`).

```bash
# Interactive PCA plot
prepare interactive_plot \
    --data-file abundance_matrix.csv \
    --metadata-file metadata.csv \
    --plot-type pca \
    --color-by treatment \
    --output-file pca_plot.html

# Interactive heatmap
prepare interactive_plot \
    --data-file abundance_matrix.csv \
    --plot-type heatmap \
    --title "Species Abundance Heatmap"

# Diversity comparison plots
prepare interactive_plot \
    --data-file abundance_matrix.csv \
    --metadata-file metadata.csv \
    --plot-type diversity \
    --color-by treatment_group
```

The plot types are `pca`, `tsne`, `heatmap` and `diversity`. `--color-by` and `--size-by` name metadata columns and need `--metadata-file`. `--no-show` writes the file without opening a browser.

## validate_taxonomy

`validate_taxonomy` checks species names against the NCBI taxonomy database. The input is a text file with one name per line, or a `.csv` file. For a CSV, `--species-column` names the column, and the first column is used when it is omitted. The email address comes from `--email`, `PREPARE_NCBI_EMAIL` or `ncbi_email` under `[runtime]`.

```bash
# Validate species from a text file
prepare validate_taxonomy \
    --species-file species_list.txt \
    --email your.email@domain.com \
    --output-file validation_results.csv

# Validate from a CSV with a named column
prepare validate_taxonomy \
    --species-file data.csv \
    --species-column organism_name \
    --email your.email@domain.com
```

`--cache-file` names the NCBI lookup cache (default `taxonomy_validation_cache.csv`), and `--api-key` raises the NCBI rate limit. This cache is separate from the GTDB cache `taxonomy_cache.tsv` that `enrich_taxonomy` and `explore_containment` use (see below); the two have different formats and are not interchangeable.

## enrich_taxonomy and taxonomic_summary

`enrich_taxonomy` looks up the GTDB taxonomy of every genome in a parsed containment table and writes a taxonomy map with the genome ids and their ranks (default `taxonomy_map.tsv`). `--cache` names the GTDB lookup cache that it reads and appends to (default `taxonomy_cache.tsv`).

`taxonomic_summary` summarises an abundance or containment table per taxonomic rank. The taxonomy table is either the map from `enrich_taxonomy` (the usual route after `parse_containment`) or the validation CSV from `validate_taxonomy`.

```bash
prepare enrich_taxonomy --parsed-containment parsed_containment.txt --output taxonomy.tsv
prepare taxonomic_summary --abundance-file parsed_containment.txt --taxonomy-file taxonomy.tsv \
    --levels phylum class order family genus --output-dir taxonomic_summaries

prepare taxonomic_summary --abundance-file abundance_matrix.csv --taxonomy-file validation_results.csv
```

`--levels` defaults to phylum, class, order, family and genus. `--min-abundance` (default 0.001) drops values below the threshold. Output goes to `--output-dir` (default `taxonomic_summaries`).

## genome_search

`genome_search` looks up GTDB genome accessions by name without downloading anything. Give `--species` or `--genus`. By default only GTDB representative genomes are listed, and `--all` includes every genome. `--format tsv` pairs each accession with the queried name, and `--output` writes to a file instead of standard output.

```bash
prepare genome_search --genus Francisella --format tsv
```

## explore_containment

`explore_containment` writes an interactive HTML browser of the containment table (default `containment_explorer.html`). `--min-containment` drops weak hits (a page with no hit left says so), `--metadata` adds an optional sample metadata table, and `--open` opens the page in a browser. Without a taxonomy for the genomes the page shows the results table only and points at `enrich_taxonomy`. It also writes a GTDB taxonomy cache (default `taxonomy_cache.tsv`) next to it. The taxonomy is computed on the fly unless `--taxonomy-map` names an existing map.

```bash
prepare explore_containment --parsed-containment parsed_containment.txt --min-containment 0.1
```

## find_by_taxonomy

`find_by_taxonomy` filters the containment table to one taxonomic group, given by `--genus`, `--family` or `--species`, using a taxonomy map from `enrich_taxonomy`. `--min-containment` adds a containment bound. `--format summary` prints a count per rank instead of the full table, and `--output` writes a TSV file instead of standard output.

```bash
prepare find_by_taxonomy --parsed-containment parsed_containment.txt --taxonomy-map taxonomy.tsv --genus Francisella
```
