#' Combine two (or more?) gene lists into a single integrated network
#'
#' @param df_list List of data frames; data to be integrated into a single
#'   network. If the list items do not have names, they will be assigned names
#'   using `LETTERS`.
#' @param col Column containing Ensembl gene IDs; must be present in all
#'   elements of `df_list`.
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
#' @importFrom purrr map imap
#'
#' @details The "min_steiner" method is implemented with the `SteinerNet`
#'   package.
#'
#' @references See <https://cran.r-project.org/web/packages/SteinerNet/index.html>
#'   for details on the Steiner network trimming.
#'
#' @seealso <https://www.github.com/travis-m-blimkie/networker>
#'
integrate_networks <- function(df_list, col, order, ppi_data = innatedb_exp, seed = 1) {

  # Clean the input df_list, starting with giving them names if not present
  if (is.null(names(df_list))) {
    names(df_list) <- LETTERS[seq(1, length(df_list))]
  }

  input_dfs_clean <- df_list %>% purrr::imap(
    ~dplyr::select(.x, "name" = col) %>%
      mutate(source = .y)
  )

  # Build a network for each data frame in df_list
  message("==> Build individual networks...")
  individual_networks <- input_dfs_clean %>% map(
    ~build_network(
      df = .x,
      col = "name",
      order = order,
      ppi_data = ppi_data,
      seed = seed
    )
  )

  # Build the combined/integrated network by first joining the inputs, making
  # sure to remove duplicates before assigning the sources
  message("\n==> Build combined network...")
  combined_seed_nodes <- input_dfs_clean %>%
    map(~distinct(.x, name, .keep_all = TRUE)) %>%
    bind_rows() %>%
    group_by(name) %>%
    summarise(source_network = paste(source, collapse = "_"))

  combined_network <- build_network(
    df       = combined_seed_nodes,
    col      = "name",
    order    = order,
    ppi_data = ppi_data,
    seed     = seed
  )

  # Identify novel nodes present only in the combined network
  message("\n==> Identify novel nodes...")
  individual_network_nodes <- individual_networks %>%
    map(~as_tibble(.x)) %>%
    bind_rows(., .id = "source_network") %>%
    pull("name")

  combined_network_nodes <- as_tibble(combined_network) %>%
    pull("name")

  novel_nodes <- setdiff(combined_network_nodes, individual_network_nodes)

  combined_network_meta <- combined_network %>%
    mutate(
      novel_node = case_when(
        (name %in% novel_nodes & !seed) ~ TRUE,
        TRUE ~ FALSE
      ),
      source_novel = case_when(
        novel_node ~ "novel_node",
        (!novel_node & !is.na(source_network)) ~ source_network,
        TRUE ~ NA_character_
      )
    )

  # Return both network objects
  message("\n==> Done.\n")
  return(
    list(
      "integrated_network"  = combined_network_meta,
      "individual_networks" = individual_networks
    )
  )
}
