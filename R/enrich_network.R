#' Test a tidygraph network object for enriched Reactome pathways
#'
#' @param network tidygraph object, i.e. one output by `build_network`
#' @param filter Desired p-value threshold for filtering output pathways.
#'   Defaults to 0.05.
#' @param background Optional background genes used when testing
#'
#' @return
#' @export
#'
#' @import ReactomePA
#' @import dplyr
#'
#'
enrich_network <- function(network, filter = 0.05, background = NULL) {
  node_table <- as_tibble(network)

  input_entrez <- clusterProfiler::bitr(
    geneID   = node_table[["name"]],
    fromType = "ENSEMBL",
    toType   = "ENTREZID",
    OrgDb    = "org.Hs.eg.db"
  ) %>% pull(ENTREZID)

  enrichPathway(
    gene     = na.omit(input_entrez),
    organism = "human",
    universe = background
  ) %>%
    pluck("result") %>%
    remove_rownames() %>%
    janitor::clean_names() %>%
    as_tibble() %>%
    filter(p_adjust < filter)
}
