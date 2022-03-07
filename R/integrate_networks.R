#' Combine two (or more?) gene lists into a single integrated network
#'
#' @param df_a First data frame of input genes, with the first column containing
#'   Ensembl gene IDs.
#' @param df_b Second data frame. The name and contents of its first column
#'   should be the same as `df_a`
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
integrate_networks <- function(df_a, df_b, order, ppi_data = innatedb_exp, seed = 1) {

  stopifnot(colnames(df_a)[1] == colnames(df_b)[1])

  input_dfs <- list("A" = df_a, "B" = df_b) %>% purrr::imap(
    ~dplyr::select(.x, "name" = 1) %>%
      mutate(source = .y)
  )

  combined_seed_nodes <-
    full_join(input_dfs[[1]], input_dfs[[2]], by = "name") %>%
    tidyr::unite(contains("source"), col = "source_network", na.rm = TRUE)

  # Build a network using the combined node list
  combined_network <- build_network(
    df       = combined_seed_nodes,
    col      = "name",
    order    = order,
    ppi_data = ppi_data,
    seed     = seed
  )

  return(combined_network)
}
