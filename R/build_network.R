#' Construct a PPI network from input genes and InnateDB's database
#'
#' @param df Input data frame containing genes of interest
#' @param col Column of input genes as Ensembl IDs (character)
#' @param order Desired network order. Possible options are "zero" (default),
#'   "first," "min_simple," or "min_steiner."
#' @param ppi_data Data frame of PPI data; must contain rows of interactions as
#'   pairs of Ensembl gene IDs, with columns named "ensembl_gene_A" and
#'   "ensembl_gene_B". Defaults to pre-packaged InnateDB PPI data.
#' @param seed Number used in call to `set.seed()` to allow for reproducible
#'   network generation
#'
#' @return `tidygraph` object for plotting or further analysis
#' @export
#'
#' @import dplyr
#' @import tidygraph
#'
#' @details The "min_steiner" method is implemented with the `SteinerNet`
#'   package.
#'
#' @references See <https://cran.r-project.org/web/packages/SteinerNet/index.html>
#'   for details on the Steiner network trimming.
#'
#' @seealso <https://www.github.com/travis-m-blimkie/networker>
#'
build_network <- function(df, col, order, ppi_data = innatedb_exp, seed = 1) {

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
  network_init <- edge_table %>%
    as_tbl_graph(directed = FALSE) %>%
    remove_subnetworks() %>%
    as_tbl_graph() %>%
    mutate(
      degree      = centrality_degree(),
      betweenness = centrality_betweenness(),
      seed        = (name %in% gene_vector)
    ) %>%
    select(-comp)


  # Perform node filtering/trimming for minimum order networks, and recalculate
  # the network statistics
  if (order == "min_simple") {

    message("Performing 'simple' minimum network trimming...")
    network_out <- network_init %>%
      filter(
        !(degree == 1 & !seed),
        !(betweenness == 0 & !seed)
      ) %>%
      mutate(
        degree      = centrality_degree(),
        betweenness = centrality_betweenness(),
      )
  } else if (order == "min_steiner") {

    message("Performing 'Steiner' minimum network trimming...")
    set.seed(seed)

    terminals <- network_init %>%
      activate(nodes) %>%
      pull(name) %>%
      intersect(gene_vector)

    network_out <- SteinerNet::steinertree(
      type      = "SP",
      terminals = terminals,
      graph     = network_init,
      color     = FALSE
    ) %>%
      magrittr::extract2(1) %>%
      as_tbl_graph(directed = FALSE) %>%
      mutate(
        degree      = centrality_degree(),
        betweenness = centrality_betweenness(),
      )
  } else {
    network_out <- network_init
  }

  if (nrow(as_tibble(network_out)) > 2000) {
    message(
      "\nWarning:\nYour network contains more than 2000 nodes, and will likely ",
      "be difficult to interpret when plotted.",
      appendLF = FALSE
    )

    if (order != "zero") {
      message(
        " Consider using a zero order network to reduce the number of nodes.\n"
      )
    }
  }

  message("Done.")
  network_out %>% left_join(., df, by = c("name" = col))
}
