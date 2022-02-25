# networker: Reproducible PPI Network Creation and Visualization in R

A suite of functions to make it simple to construct PPI networks inside of R,
with an emphasis on usability and reproducibility.

## Installation
You can install `networker` with the following (assuming you have `devtools`
already installed):
```r
devtools::install_github("https://github.com/travis-m-blimkie/networker")
```

If you run into installation issues relating to the package `biomaRt` (which is
a dependency of `networker`), you can install it with the following:
```r
if (!require("BiocManager")) { # Optionally install BiocManager if needed
  install.packages("BiocManager")
}
BiocManager::install("biomaRt")
```

## Example
```r
> library(tidyverse)
# ── Attaching packages ────────────────────────────────── tidyverse 1.3.1 ──
# ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
# ✓ tibble  3.1.6     ✓ dplyr   1.0.8
# ✓ tidyr   1.2.0     ✓ stringr 1.4.0
# ✓ readr   2.1.2     ✓ forcats 0.5.1
# ── Conflicts ───────────────────────────────────── tidyverse_conflicts() ──
# x dplyr::filter() masks stats::filter()
# x dplyr::lag()    masks stats::lag()

> library(networker)
# Thanks for using networker! If you encounter any bugs or problems, please 
# submit an issue at the Github page: 
# https://github.com/travis-m-blimkie/networker/issues

> fetch_ppi_data()
# Downloading data from InnateDB...
# Retrieving gene mapping data from biomaRt...                                                      
# Done, saving InnateDB PPI data to:
# /home/networker/innateDB_ppi_exp_data_mapped.Rds

> innatedb <- readRDS("innateDB_ppi_exp_data_mapped.Rds")
> ex_genes <- read_csv("ex_de_genes.csv")

> glimpse(ex_genes)
# Rows: 335
# Columns: 9
# $ ensembl_gene_id <chr> "ENSG00000003147", "ENSG00000005381", "ENSG0000000…
# $ hgnc_symbol     <chr> "ICA1", "MPO", "PGLYRP1", "MMP25", "ANLN", "LTF", …
# $ entrez_gene_id  <dbl> 3382, 4353, 8993, 64386, 54443, 4057, 6556, 3082, …
# $ baseMean        <dbl> 29.63805, 332.51245, 310.38400, 3654.71160, 22.318…
# $ log2FoldChange  <dbl> -0.6631473, -1.0411931, -0.6977906, -0.5962315, -0…
# $ lfcSE           <dbl> 0.12867106, 0.28542343, 0.18071716, 0.12173940, 0.…
# $ stat            <dbl> -5.153818, -3.647889, -3.861231, -4.897605, -3.332…
# $ pvalue          <dbl> 2.552353e-07, 2.644034e-04, 1.128174e-04, 9.701195…
# $ padj            <dbl> 1.010786e-05, 2.199332e-03, 1.115054e-03, 2.965279…

> ex_network <- build_network(
  df       = ex_genes,
  col      = "ensembl_gene_id",
  order    = "min_steiner",
  ppi_data = innatedb
)
# Finding interactions...
# Creating network...
# Performing 'Steiner' minimum network trimming...
# Done.

> plot_network(
  network      = ex_network,
  layout       = "force_atlas",
  fill_column  = log2FoldChange,
  edge_alpha   = 0.5,
  edge_colour  = "grey70",
  label        = TRUE,
  label_column = hgnc_symbol,
  label_filter = 40,
  label_size   = 5,
  label_colour = "dodgerblue4",
  label_face   = "bold"
)
# Calculating Force Atlas node positions...
# Warning message:
# Removed 221 rows containing missing values (geom_text_repel).
```

![](man/figures/network_example.png)

## Versioning
This package makes use of [SemVer](https://semver.org/).

## Authors
Travis Blimkie is the originator and principal contributor. You can check the
list of all contributors [here](https://github.com/travis-m-blimkie/networker/graphs/contributors).

## License
This project is written under the GPLv3 license, available
[here.](https://github.com/travis-m-blimkie/networker/blob/main/LICENSE.md)
