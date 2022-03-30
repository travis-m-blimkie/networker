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

  # Check for and remove any duplicate IDs, which will cause problems later.
  # Make sure to warn the user about this.
  message("Cleaning input data...")
  df_clean <- distinct(df, !!sym(col), .keep_all = TRUE)
  gene_vector <- unique(df_clean[[col]])
  lost_ids <- df[[col]][duplicated(df[[col]])]

  if (length(gene_vector) < nrow(df)) {
    message(paste0(
      "  INFO: Found ",
      nrow(df) - length(gene_vector),
      " duplicate IDs in the input column, which have been removed:"
    ))
    message("  ", paste(lost_ids, collapse = ", "))
  }

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
      seed        = (name %in% gene_vector),
      hub_score   = centrality_hub()
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
        hub_score   = centrality_hub()
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
        hub_score   = centrality_hub()
      )
  } else {
    network_out <- network_init
  }

  if (nrow(as_tibble(network_out)) > 2000) {
    message(
      "Warning:\nYour network contains more than 2000 nodes, and will likely ",
      "be difficult to interpret when plotted."
    )
  } else if (nrow(as_tibble(network_out)) > 2000 & order != "zero") {
    message(
      "\nWarning:\nYour network contains more than 2000 nodes, and will ",
      "likely be difficult to interpret when plotted. Consider switching to a ",
      "zero order network to improve legibility.\n"
    )
  }

  message("Mapping input Ensembl IDs to HGNC symbols...")
  ensembl_to_hgnc <- biomart_id_mapping_human %>%
    select("name" = ensembl_gene_id, "gene_name" = hgnc_symbol)

  network_mapped <- left_join(
    network_out,
    ensembl_to_hgnc,
    by = "name"
  )

  network_final <- left_join(
    network_mapped,
    df_clean,
    by = c("name" = col)
  )

  attr(network_final, "order") <- order

  message("Done.\n")
  return(network_final)
}
