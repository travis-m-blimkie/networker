#' Plot an undirected PPI network using ggraph
#'
#' @param network `tidygraph` object, output from `build_network`
#' @param fill_column Tidy-select column for mapping node colour. Designed to
#'   handle continuous numeric mappings (either positive/negative only, or
#'   both), and categorical mappings. Two-sided numeric columns will be mapped
#'   to green and red for negative and positive values, respectively.
#' @param fill_type Type of fill mapping to perform for nodes, based on supplied
#'   `fill_column`. Options are: "fold_change", "two_sided", "one_sided", or
#'   "categorical".
#' @param layout Layout of nodes in the network. Supports all layouts from
#'   `ggraph`/`igraph`, as well as "force_atlas" (see Details)
#' @param legend Should a legend be included? Defaults to FALSE
#' @param fontfamily Font to use for labels and legend (if present). Defaults to
#'   "Helvetica".
#' @param edge_colour Edge colour, defaults to "grey50"
#' @param edge_alpha Transparency of edges, defaults to 0.2
#' @param node_size Numeric vector of length two, specifying size range of nodes
#' (maps to node degree). Default is `c(3, 9)`.
#' @param int_colour Fill colour for non-seed nodes, i.e. interactors. Defaults
#'   to "grey70"
#' @param label Boolean, whether labels should be added to nodes. Defaults to
#'   FALSE. Note when TRUE, only seed nodes receive labels.
#' @param label_column Tidy-select column of the network/data to be used in
#'   labeling nodes.
#' @param label_filter Degree filter used to determine if a given node should be
#'   labeled. Defaults to 0. This value can be tweaked to reduce the number of
#'   node labels, to prevent the network from being too crowded.
#' @param label_size Size of node labels, defaults to 5.
#' @param label_colour Colour of node labels, defaults to "black"
#' @param hub_colour Colour of node labels for hubs. The top 2% of nodes (based
#'   on calculated hub score) are highlighted with this colour, if `label =
#'   TRUE`.
#' @param label_face Font face for node labels, defaults to "bold"
#' @param label_padding Padding around the label, defaults to 0.25 lines.
#' @param min_seg_length Minimum length of lines to be drawn from labels to
#'   points. The default specified here is 0.25, half of the normal default
#'   value.
#' @param subnet Logical determining if networks produced by
#'   `extract_subnetwork` should be treated as such, or just as a normal network
#'   from `build_network`.
#' @param seed Number used in call to `set.seed()` to allow for reproducible
#'   network generation. Can be changed to get slightly different layouts.
#' @param ... Further parameters to be passed on to `ggplot2::theme()`, e.g.
#'   `legend.position`
#'
#' @return An object of class "gg"
#'
#' @export
#'
#' @import ggplot2
#' @import ggraph
#' @import dplyr
#'
#' @details Any layout supported by ggraph can be specified here - see
#'   `?layout_tbl_graph_igraph` for a list of options. Additionally, there is
#'   support for the "force_atlas" method, implemented via the ForceAtlas2
#'   package.
#'
#'   The `fill_type` argument will determine how the node colour is mapped to
#'   the desired column. "fold_change" represents a special case, where the fill
#'   column is numeric and whose values should be mapped to up (> 0) or down (<
#'   0). "two_sided" and "one_sided" are designed for numeric data that contains
#'   either positive and negative values, or only positive/negative values,
#'   respectively. "categorical" is designed for non-numeric colour mapping.
#'
#'   Node statistics (degree, betweenness, and hub score) and calculated using
#'   the respective functions from the `tidygraph` package.
#'
#'   If plotting a network created from the `extract_subnetwork` function, the
#'   genes belonging to the extracted pathway (i.e. the contents of the
#'   "gene_id" column in the enrichment results) will be highlighted, instead of
#'   hub nodes. The colour used here can be controlled via the `hub_colour`
#'   argument, and this behaviour can be turned off altogether by setting the
#'   `subnet` argument to FALSE.
#'
#' @references See <https://github.com/analyxcompany/ForceAtlas2> for details on
#'   this method.
#'
#' @seealso <https://github.com/travis-m-blimkie/networker>
#'
plot_network <- function(
  network,
  fill_column,
  fill_type,
  layout         = "kk",
  legend         = FALSE,
  fontfamily     = "Helvetica",
  edge_colour    = "grey50",
  edge_alpha     = 0.2,
  node_size      = c(3, 9),
  int_colour     = "grey70",
  label          = FALSE,
  label_column,
  label_filter   = 0,
  label_size     = 5,
  label_colour   = "black",
  hub_colour     = "blue2",
  label_face     = "bold",
  label_padding  = 0.25,
  min_seg_length = 0.25,
  subnet         = TRUE,
  seed           = 1,
  ...
) {

  # Set up fill scaling based on argument `fill_type`
  if (fill_type == "fold_change") {
    network <- network %>%
      mutate(
        new_fill_col = case_when(
          {{fill_column}} < 0 ~ "Down",
          {{fill_column}} > 0 ~ "Up",
          TRUE ~ NA_character_
        )
      )
    network_fill_geom <- scale_fill_manual(
      values   = c("Up" = "firebrick3", "Down" = "#188119"),
      na.value = int_colour
    )
    network_fill_guide <-
      guides(fill = guide_legend(override.aes = list(size = 5)))

  } else if (fill_type == "two_sided") {
    network <- network %>%
      mutate(new_fill_col = {{fill_column}})
    network_fill_geom <- scale_fill_distiller(
      palette  = "RdYlBu",
      na.value = int_colour,
      guide    = ifelse(legend, "colourbar", "none")
    )
    network_fill_guide <- NULL

  } else if (fill_type == "one_sided") {
    network <- mutate(network, new_fill_col = {{fill_column}})
    network_fill_geom <- scale_fill_viridis_c()
    network_fill_guide <- NULL

  } else if (fill_type == "categorical") {
    network <- network %>%
      mutate(new_fill_col = {{fill_column}})
    network_fill_geom <- scale_fill_brewer(
      palette  = "Set1",
      na.value = int_colour,
      guide    = ifelse(legend, "legend", "none")
    )
    network_fill_guide <-
      guides(fill = guide_legend(override.aes = list(size = 5)))
  } else {
    stop("Argument '' must be one of 'fold_change', 'two_sided', 'one_sided', ",
         "or 'categorical'")
  }

  # If we're using the Force Atlas layout, we need to pre-calculate the node
  # positions using the appropriate function from the ForceAtlas2 package
  if (layout == "force_atlas") {
    message("Calculating Force Atlas node positions...")
    set.seed(seed)
    layout_object <- ForceAtlas2::layout.forceatlas2(
      graph    = network,
      directed = FALSE,
      plotstep = 0
    )
  } else {
    layout_object <- layout
  }

  # Set a plain white background
  set_graph_style(foreground = "white")

  # Get fill_column as a string, so we can clean it up if the legend is included
  legend_name <- match.call()$fill_column

  if ( subnet & "starters" %in% names(attributes(network)) ) {

    message(
      "Detected this is a sub-network generated by `extract_subnetwork()`.\n",
      "Blue node labels indicate genes from the extracted pathway."
    )

    starter_nodes <- attr(network, "starters")

    network <- network %>%
      mutate(
        node_label = case_when(
          degree > label_filter ~ {{label_column}},
          TRUE ~ NA_character_
        ),
        is_starter = case_when(
          name %in% starter_nodes ~ "y",
          TRUE ~ "n"
        )
      )

    ggraph(network, layout = layout_object) +
      geom_edge_link(show.legend = FALSE, alpha = edge_alpha, colour = edge_colour) +
      geom_node_point(aes(size = degree, fill = new_fill_col), pch = 21, colour = "grey70") +
      network_fill_geom +
      geom_node_text(
        aes(label = node_label, colour = is_starter),
        size          = label_size,
        repel         = TRUE,
        family        = fontfamily,
        fontface      = label_face,
        check_overlap = TRUE,
        show.legend   = FALSE,
        box.padding   = label_padding,
        min.segment.length = min_seg_length
      ) +
      scale_size_continuous(range = node_size, guide = "none") +
      scale_colour_manual(values = c("y" = hub_colour, "n" = label_colour)) +
      labs(fill = NULL) +
      theme(
        text         = element_text(family = fontfamily),
        plot.margin  = unit(rep(0, 4), "cm"),
        legend.text  = element_text(size = 14),
        ...
      ) +
      network_fill_guide

  } else if (label) {
    hub_nodes <- as_tibble(network) %>%
      arrange(desc(hub_score)) %>%
      slice_head(n = 3 + ceiling(nrow(as_tibble(network)) * 0.01)) %>%
      pull(name)

    network <- network %>%
      mutate(
        node_label = case_when(
          degree > label_filter ~ {{label_column}},
          TRUE ~ NA_character_
        ),
        is_hub = case_when(
          name %in% hub_nodes ~ "y",
          TRUE ~ "n"
        )
      )

    ggraph(network, layout = layout_object) +
      geom_edge_link(show.legend = FALSE, alpha = edge_alpha, colour = edge_colour) +
      geom_node_point(aes(size = degree, fill = new_fill_col), pch = 21, colour = "grey70") +
      network_fill_geom +
      geom_node_text(
        aes(label = node_label, colour = is_hub),
        size          = label_size,
        repel         = TRUE,
        family        = fontfamily,
        fontface      = label_face,
        check_overlap = TRUE,
        show.legend   = FALSE,
        box.padding   = label_padding,
        min.segment.length = min_seg_length
      ) +
      scale_size_continuous(range = node_size, guide = "none") +
      scale_colour_manual(values = c("y" = hub_colour, "n" = label_colour)) +
      labs(fill = NULL) +
      theme(
        text         = element_text(family = fontfamily),
        plot.margin  = unit(rep(0, 4), "cm"),
        legend.text  = element_text(size = 14),
        ...
      ) +
      network_fill_guide

  } else {
    ggraph(network, layout = layout_object) +
      geom_edge_link(show.legend = FALSE, alpha = edge_alpha, colour = edge_colour) +
      geom_node_point(aes(size = degree, fill = new_fill_col), pch = 21) +
      network_fill_geom +
      scale_size_continuous(range = node_size, guide = "none") +
      labs(fill = NULL) +
      theme(
        text         = element_text(family = fontfamily),
        plot.margin  = unit(rep(0, 4), "cm"),
        legend.text  = element_text(size = 14),
        ...
      )
  }
}
