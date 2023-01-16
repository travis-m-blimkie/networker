#' Extract a subnetwork based on pathway genes
#'
#' @param network Input network object; output from `build_network`
#' @param genes List of Ensembl gene IDs to use as the starting point to extract
#'   a subnetwork from the initial network. You must provide either `genes` or
#'   `enrich_result` argument.
#' @param enrich_result Pathway enrichment result, output from `enrich_network`.
#'   You must provide either `genes` or `enrich_result` argument.
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
#'   connected subnetwork or module from the starting network, using either a
#'   list of Ensembl genes or genes from an enriched pathway as the basis. To
#'   see what genes were pulled out for the pathway, check `attr(x,
#'   "starters")`.
#'
#' @references Code for network module (subnetwork) extraction was based off of
#' that used in jboktor/NetworkAnalystR on Github.
#'
#' @seealso <https://www.github.com/travis-m-blimkie/networker>
#'
extract_subnetwork <- function(network, genes = NULL, enrich_result = NULL, pathway_name) {

  if (is.null(genes) & is.null(enrich_result)) {
    stop("You must specify either 'genes' or 'enrich_result' to provide genes ",
         "to extract from the initial network")
  }

  message("Checking inputs...")
  if (!is.null(genes)) {
    if ( !grepl(x = genes[[1]], pattern = "^ENSG") ) {
      stop("Argument 'genes' must be a list of Ensembl gene IDs")
    }
  }

  if (!is.null(enrich_result)) {
    if ( !all(c("description", "gene_id") %in% colnames(enrich_result)) ) {
      stop("Argument 'enrich_result' must contain the columns 'description' ",
           "and 'gene_id'")
    }

    if (!pathway_name %in% enrich_result[["description"]]) {
      stop("Argument 'pathway_name' must be present in the 'description' ",
           "column of the 'enrich_result' object")
    }

    if ( !grepl(enrich_result[["gene_id"]][1], pattern = "([0-9]{2,5}/)+") ) {
      stop(
        "The 'gene_id' column must contain Entrez gene IDs separated with a '/'"
      )
    }
  }


  if (!is.null(genes)) {
    message("Using provided list of Ensembl genes...", appendLF = FALSE)
    genes_to_extract <- genes

  } else if (!is.null(enrich_result)) {
    message("Pulling genes for given pathway...", appendLF = FALSE)
    pathway_genes_entrez <- enrich_result %>%
      filter(description == pathway_name) %>%
      pull(gene_id) %>%
      strsplit(., split = "/") %>%
      unlist()

    genes_to_extract <- biomart_id_mapping_human %>%
      filter(entrez_gene_id %in% pathway_genes_entrez) %>%
      pull(ensembl_gene_id) %>%
      unique()
  }

  message("found ", length(genes_to_extract), " genes.")


  gene_node_ids <- as_tibble(network) %>%
    mutate(rn = row_number()) %>%
    filter(name %in% genes_to_extract) %>%
    pull(rn)

  # Get subgraphs, which will only contain the specified nodes
  message("Calculating subgraphs from specified nodes...")
  all_subgraphs <- induced_subgraph(
    graph = as.igraph(network),
    vids = gene_node_ids
  )

  # Decompose each subgraph
  all_components <-
    decompose.graph(graph = all_subgraphs, min.vertices = 1)

  # If all the specified nodes form a single, connected network, pull that...
  if (length(all_components) == 1) {
    message("All nodes form a single network...")
    module_network <- all_components[[1]]

    # ...or we need to minimally connect each subgraph we've identified
  } else {
    message("Determining shortest paths between nodes...")

    module_shortest_paths <- list()

    for (i in 1:length(gene_node_ids)) {
      module_shortest_paths[[i]] <- get.shortest.paths(
        as.igraph(network),
        gene_node_ids[i],
        gene_node_ids[-(1:i)]
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

  attr(module_network_tidygraph, "starters") <- genes_to_extract
  return(module_network_tidygraph)
}
