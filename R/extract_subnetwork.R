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
#' @importFrom igraph as.igraph V induced_subgraph decompose.graph simplify delete.vertices get.shortest.paths
#'
#' @details Uses functions from the igraph package to extract a minimally
#'   connected "module" from the starting network, using genes from a given
#'   pathway as the basis. To see what genes were pulled out for the pathway,
#'   use `attr(output, "starters")`.
#'
#' @references Code for network module (subnetwork) extraction was based off of
#' that used in jboktor/NetworkAnalystR on Github.
#'
#' @seealso <https://www.github.com/travis-m-blimkie/networker>
#'
extract_subnetwork <- function(network, enrich_result, pathway_name) {

  stopifnot(pathway_name %in% enrich_result[["description"]])

  message("Pulling genes for given pathway...", appendLF = FALSE)
  pathway_genes_entrez <- enrich_result %>%
    filter(description == pathway_name) %>%
    pull(gene_id) %>%
    strsplit(., split = "/") %>%
    unlist()

  message("found ", length(pathway_genes_entrez), " genes.")

  pathway_genes_ensembl <- biomart_id_mapping_human %>%
    filter(entrez_gene_id %in% pathway_genes_entrez) %>%
    pull(ensembl_gene_id) %>%
    unique()

  pathway_node_ids <- as_tibble(network) %>%
    mutate(rn = row_number()) %>%
    filter(name %in% pathway_genes_ensembl) %>%
    pull(rn)

  # Get subgraphs, which will only contain the specified nodes
  message("Calculating subgraphs from specified nodes...")
  pathway_subgraphs <- induced_subgraph(
    graph = as.igraph(network),
    vids = pathway_node_ids
  )

  # Decompose each subgraph
  pathway_components <-
    decompose.graph(graph = pathway_subgraphs, min.vertices = 1)

  # If all the specified nodes form a single, connected network, pull that...
  if (length(pathway_components) == 1) {
    message("All nodes form a single network...")
    module_network <- pathway_components[[1]]

    # ...or we need to minimally connect each subgraph we've identified
  } else {
    message("Determining shortest paths between nodes...")

    module_shortest_paths <- list()

    for (i in 1:length(pathway_node_ids)) {
      module_shortest_paths[[i]] <- get.shortest.paths(
        as.igraph(network),
        pathway_node_ids[i],
        pathway_node_ids[-(1:i)]
      )$vpath
    }

    module_nodes <- unique(unlist(module_shortest_paths))
    nodes_to_remove <- V(network)$name[-module_nodes]
    module_network <- simplify(delete.vertices(network, nodes_to_remove))
  }

  module_network_tidygraph <- as_tbl_graph(module_network)
  module_network_tibble    <- as_tibble(module_network_tidygraph)

  message("Done, new subnetwork contains ",
          nrow(module_network_tibble),
          " nodes.\n")

  attr(module_network_tidygraph, "starters") <- pathway_genes_ensembl
  return(module_network_tidygraph)
}
