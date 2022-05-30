#' Test a tidygraph network object for enriched Reactome pathways
#'
#' @param network tidygraph object, i.e. one output by `build_network` with
#'   Ensembl IDs as node names
#' @param filter Desired p-value threshold for filtering output pathways.
#'   Defaults to 0.05.
#' @param background Optional background genes used when testing
#'
#' @return A tibble of enriched Reactome pathways
#'
#' @export
#'
#' @import dplyr
#'
#' @references See <https://www.bioconductor.org/packages/ReactomePA/> for
#'   details on ReactomePA's methods.
#'
#' @seealso <https://github.com/travis-m-blimkie/networker>
#'
enrich_network <- function(network, filter = 0.05, background = NULL) {

  message("Checking input...")
  if (!all(c("tbl_graph", "igraph") %in% class(network_1))) {
    stop("Argument 'network' must be a tidygraph object, produced by ",
         "'build_network()'")
  }

  input_entrez <- as_tibble(network) %>%
    dplyr::select(name) %>%
    left_join(biomart_id_mapping_human, by = c("name" = "ensembl_gene_id")) %>%
    pull(entrez_gene_id)

  message("Testing ", length(input_entrez), " genes for enriched pathways...")
  output <- ReactomePA::enrichPathway(
    gene     = na.omit(input_entrez),
    organism = "human",
    universe = background
  ) %>%
    purrr::pluck("result") %>%
    tibble::remove_rownames() %>%
    janitor::clean_names() %>%
    as_tibble() %>%
    filter(p_adjust < filter) %>%
    dplyr::select(id, description, "p_value" = pvalue, p_adjust, gene_id)

  message("Done.\n")
  return(output)
}
