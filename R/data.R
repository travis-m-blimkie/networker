#' Experimentally verified human PPI data downloaded from InnateDB
#'
#' A data frame containing human PPI data from InnateDB, from the entry
#' "All Experimentally Validated Interactions (updated weekly)" at
#' <https://innatedb.com/redirect.do?go=downloadImported>. A few important steps
#' have been taken to filter the data, namely the removal of duplicate
#' interactions, and removing interactions that have the same components but are
#' swapped between A and B
#'
#' @format A data frame with 152259 rows and 4 columns:
#' \describe{
#'   \item{ensembl_gene_A}{Ensembl gene ID for the first gene/protein in the
#'     interaction}
#'   \item{hgnc_symbol_A}{HGNC symbol for the first gene/protein in the
#'     interaction}
#'   \item{ensembl_gene_B}{Ensembl gene ID for the second gene/protein in the
#'     interaction}
#'   \item{hgnc_symbol_B}{HGNC symbol for the second gene/protein in the
#'     interaction}
#' }
"innatedb_exp"


#' Human gene ID mappings from biomaRt
#'
#' A tibble containing gene ID mapping information for Ensembl, Entrez, and
#' HGNC gene identifiers.
#'
#' @format A data frame with 68005 rows and 3 columns"
#' \describe{
#'   \item{ensembl_gene_id}{Ensembl gene IDs}
#'   \item{hgnc_symbol}{HGNC symbols}
#'   \item{entrez_gene_id}{Entrez (NCBI) gene IDs}
#' }
"biomart_id_mapping_human"
