# Heatmap and PCA plotting (R Markdown)

This folder contains the R Markdown script [Exploratory_analysis.Rmd](./scripts/Exploratory_analysis.Rmd) within the scripts directory that generates a heatmap and a PCA plot from RNA-seq read-count data to assess sample relationships and quality.

#### Inputs:

```bash
# A gene-by-sample raw-read-count table (rows = genes, columns = samples).
read_counts_table_all_samples.txt

# A sample metadata file (tab-delimited).
metadata.txt
```

#### Usage

- Edit the R Markdown parameters under `# INPUT files` to point to your count matrix and metadata files. Manuscript [metadata](./data/metadata.txt) and [read-counts](./data/read_counts_table_all_samples.txt) tables are located within the [data](./data/) directory.

#### Required R packages:

- DESeq2 (used for normalization, v1.46.0)
- pheatmap (for heatmaps, v1.0.12)
- ggplot2 (for plotting, v3.5.1)
- ggsci (for colors, v3.2.0)
- dplyr (for data manipulation, v1.1.4)

#### What the script does

- Normalizes counts (default method configurable in the Rmd).
- Performs variance-stabilizing transformation
- Produces a clustered heatmap of sample-to-sample distances.
- Computes PCA on transformed counts and plots PC1 vs PC2 with samples colored by treatment + genotype.

#### Outputs

- A heatmap image (PNG).
- A PCA scatter plot.
- Optional data tables of transformed counts and sample distances saved to disk (if enabled in the script).