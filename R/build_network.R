#' Construct a PPI network from input genes and InnateDB's database
#'
#' @param df Input data frame containing genes of interest
#' @param col Column of input genes as Ensembl IDs (character)
#' @param order Desired network order. Possible options are "zero" (default),
#'   "first," "min_simple," or "min_steiner."
#' @param ppi_data Data frame of InnateDB PPI data; minimally should contain
#'   rows of interactions as pairs of Ensembl gene IDs, e.g. "ensembl_gene_A"
#'   and "ensembl_gene_B"
#' @param seed Number used in call to `set.seed()` to allow for reproducible
#'   network generation
#'
#' @return `tidygraph` object for plotting or further analysis
#' @export
#'
#' @importFrom igraph V components induced_subgraph
#' @import dplyr
#' @import tidygraph
#'
#' @references None.
#'
#' @seealso <https://www.github.com/travis-m-blimkie/networker>
#'
build_network <- function(df, col, order, ppi_data, seed = 1) {

  remove_subnetworks <- function(graph) {
    V(graph)$comp <- components(graph)$membership
    induced_subgraph(graph, V(graph)$comp == 1)
  }

  gene_vector <- df[[col]]

  if (!grepl(x = gene_vector[1], pattern = "^ENSG")) {
    stop("Input genes must be human Ensembl IDs")
  }

  if ( !all(c("ensembl_gene_A", "ensembl_gene_B") %in% colnames(ppi_data)) ) {
    stop("Argument 'ppi_data' must be a data frame containing columns ",
         "'ensembl_gene_A' and 'ensembl_gene_B' (case sensitive)")
  }

  ppi_data_ensembl <- ppi_data %>%
    select(starts_with("ensembl"))

  message("Finding interactions...")
  if (order == "zero") {
    edge_table <- ppi_data_ensembl %>%
      filter(ensembl_gene_A %in% gene_vector & ensembl_gene_B %in% gene_vector)
  } else if (order %in% c("first", "min_simple", "min_steiner")) {
    edge_table <- ppi_data_ensembl %>%
      filter(ensembl_gene_A %in% gene_vector | ensembl_gene_B %in% gene_vector)
  } else {
    stop(
      "Argument 'order' must be one of: ",
      "'zero', 'first', 'min_simple', or 'min_steiner'"
    )
  }

  message("Creating network...")
  network <- edge_table %>%
    as_tbl_graph(directed = FALSE) %>%
    remove_subnetworks() %>%
    as_tbl_graph(directed = FALSE) %>%
    mutate(
      degree      = centrality_degree(),
      betweenness = centrality_betweenness(),
      seed        = (name %in% gene_vector)
    ) %>%
    dplyr::select(-comp)

  if (order == "min_simple") {

    message("Performing 'simple' minimum network trimming...")
    network <- filter(
      network,
      !(degree == 1 & !seed),
      !(betweenness == 0 & !seed)
    )
  } else if (order == "min_steiner") {

    message("Performing 'Steiner' minimum network trimming...")
    set.seed(seed)

    terminals <- network %>%
      activate(nodes) %>%
      pull(name) %>%
      intersect(gene_vector)

    network <- SteinerNet::steinertree(
      type      = "SP",
      terminals = terminals,
      graph     = network,
      color     = FALSE
    ) %>%
      magrittr::extract2(1) %>%
      as_tbl_graph()
  }

  message("Done.")
  network %>% left_join(df, by = c("name" = col))
}
