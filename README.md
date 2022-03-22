# networker: Reproducible PPI Network Creation and Visualization in R

A suite of functions to make it simple to construct PPI networks inside of R,
with an emphasis on usability and reproducibility.

## Installation
You can install `networker` with the following (assuming you have `devtools`
already installed):
```r
devtools::install_github("https://github.com/travis-m-blimkie/networker")
```

If you run into installation issues relating to the package `biomaRt` or `ReactomePA` (which are
dependencies of `networker`), you can install them with the following:
```r
if (!require("BiocManager")) { # Optionally install BiocManager if needed
  install.packages("BiocManager")
}
BiocManager::install("biomaRt", "ReactomePA")
```

## Example
```r
> library(networker)
# Thanks for using networker! If you encounter any bugs or problems, please 
# submit an issue at the Github page: 
# https://github.com/travis-m-blimkie/networker/issues

> ex_genes <- read_csv("ex_de_genes.csv")
> glimpse(ex_genes)
# Rows: 158
# Columns: 3
# $ ensembl_gene_id <chr> "ENSG00000160654", "ENSG00000018280", "ENSG00000179...
# $ hgnc_symbol     <chr> "CD3G", "SLC11A1", "CIITA", "EIF4E3", "RNASE6", "AL...
# $ log2FoldChange  <dbl> 0.07662557, 1.02935351, -2.59604031, -0.47667505, 2...

> ex_network <- build_network(
    df    = ex_genes,
    col   = "ensembl_gene_id",
    order = "min_steiner"
  )
# ==> INFO: Found 3 duplicate IDs in the input column, which have been removed:
# ENSG00000172724, ENSG00000169397, ENSG00000169385
# 
# Finding interactions...
# Creating network...
# Performing 'Steiner' minimum network trimming...
# Done.

> plot_network(
    network      = ex_net,
    fill_column  = log2FoldChange,
    layout       = "force_atlas",
    label        = TRUE,
    label_column = gene_name,
    label_filter = 2,
    legend = TRUE
  )
# Calculating Force Atlas node positions...
# Warning message:
# Removed 166 rows containing missing values (geom_text_repel).
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
