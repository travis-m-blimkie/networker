#' Extract a subnetwork based on pathway genes
#'
#' @param network Input network object; output from `build_network`
#' @param enrich_result Pathway enrichment result, output from `enrich_network`
#' @param pathway_name Name of the pathway determining what genes (nodes) are
#'   pulled from the input network.
#'
#' @return `tidygraph` object for plotting or further analysis
#'
#' @export
#'
#' @import dplyr
#'
#' @seealso <https://www.github.com/travis-m-blimkie/networker>
#'
extract_subnetwork <- function(network, enrich_result, pathway_name) {

  stopifnot(pathway_name %in% enrich_result[["description"]])

  pathway_genes_entrez <- enrich_result %>%
    filter(description == pathway_name) %>%
    tidyr::separate_rows(gene_id, sep = "/") %>%
    pull(gene_id)

  pathway_genes_ensembl <-
    filter(biomart_id_mapping_human, entrez_gene_id %in% pathway_genes_entrez) %>%
    distinct(ensembl_gene_id) %>%
    drop_na() %>%
    pull(ensembl_gene_id)

  print(
    filter(biomart_id_mapping_human, ensembl_gene_id %in% pathway_genes_ensembl) %>%
      pull(hgnc_symbol)
  )

  pathway_genes_num <- as_tibble(network) %>%
    mutate(rn = row_number()) %>%
    filter(name %in% pathway_genes_ensembl) %>%
    pull(rn)

  convert(
    .data = network,
    .f    = to_local_neighborhood,
    node  = pathway_genes_num,
    order = 2)
}
