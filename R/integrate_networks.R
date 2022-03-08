#' Combine two (or more?) gene lists into a single integrated network
#'
#' @param df_a First data frame of input genes, with the first column containing
#'   Ensembl gene IDs.
#' @param df_b Second data frame. The name and contents of its first column
#'   should be the same as `df_a`
#' @param col Column containing Ensembl gene IDs; must be present in both `df_a`
#'   and `df_b`
#' @param order Desired type of output network ("zero", "first", "min_simple" or
#'   "min_steiner"); passed to `build_network()`
#' @param ppi_data Data frame of PPI data; must contain rows of interactions as
#'   pairs of Ensembl gene IDs, with columns named "ensembl_gene_A" and
#'   "ensembl_gene_B". Defaults to pre-packaged InnateDB PPI data.
#' @param seed Number used in call to `set.seed()` to allow for reproducible
#'   network generation
#'
#' @return `tidygraph` object for plotting or further analysis
#'
#' @export
#'
#' @import dplyr
#'
#' @details The "min_steiner" method is implemented with the `SteinerNet`
#'   package.
#'
#' @references See <https://cran.r-project.org/web/packages/SteinerNet/index.html>
#'   for details on the Steiner network trimming.
#'
#' @seealso <https://www.github.com/travis-m-blimkie/networker>
#'
integrate_networks <- function(df_a, df_b, col, order, ppi_data = innatedb_exp, seed = 1) {

  # Clean the input df_a and df_b
  input_dfs <- list("A" = df_a, "B" = df_b)

  input_dfs_clean <- input_dfs %>% purrr::imap(
    ~dplyr::select(.x, "name" = col) %>%
      mutate(source = .y)
  )

  # Build a network for df_a and df_b
  message("Build individual networks...")
  individual_networks <- input_dfs_clean %>% purrr::map(
    ~build_network(
      df = .x,
      col = "name",
      order = order,
      ppi_data = ppi_data,
      seed = seed
    )
  )

  # Build the combined/integrated network
  message("\nBuild combined network...")
  combined_seed_nodes <-
    full_join(input_dfs_clean[[1]], input_dfs_clean[[2]], by = "name") %>%
    tidyr::unite(contains("source"), col = "source_network", na.rm = TRUE)

  combined_network <- build_network(
    df       = combined_seed_nodes,
    col      = "name",
    order    = order,
    ppi_data = ppi_data,
    seed     = seed
  )

  # Identify novel nodes present only in the combined network
  message("\nIdentify novel nodes...")
  individual_network_nodes <- individual_networks %>%
    purrr::map(~as_tibble(.x)) %>%
    bind_rows(., .id = "source_network") %>%
    pull("name")

  combined_network_nodes <- as_tibble(combined_network) %>%
    pull("name")

  novel_nodes <- setdiff(combined_network_nodes, individual_network_nodes)

  combined_network_meta <- combined_network %>%
    mutate(novel_node = case_when(
      name %in% novel_nodes ~ TRUE,
      TRUE ~ FALSE
    ))


  # Return all three network objects
  message("Done.")
  return(
    list(
      "integrated_network"  = combined_network_meta,
      "individual_networks" = individual_networks
    )
  )
}
