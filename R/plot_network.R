#' Plot an undirected network using ggraph
#'
#' @param network `tidygraph` object, output from `build_netork`
#' @param layout Layout of nodes in the network. Supports all layouts from
#'   `ggraph`/`igraph`, as well as "force_atlas" (see Details)
#' @param fill_column Tidy-select column for mapping node colour. The colour
#'   scale set up is green-white-red, designed to map to fold change values
#' @param edge_colour Edge colour, defaults to "grey50"
#' @param edge_alpha Transparency of edges, defaults to 0.2
#' @param label Boolean, whether labels should be added to nodes. Defaults to
#'   FALSE. Note when TRUE, only seed nodes receive labels.
#' @param label_filter Degree filter used to determine if a given node should be
#'   labeled. Defaults to 40. This value can be tweaked to reduce the number of
#'   node labels, to prevent the network from being too crowded.
#' @param label_column Tidy-select column of the network/data to be used in
#'   labeling nodes.
#' @param label_size Size of node labels, defaults to 4.
#' @param label_colour Colour of node labels, defaults to "black"
#' @param label_face Font face for node labels, defaults to "bold"
#'
#' @return
#' @export
#'
#' @import tidygraph
#' @import ggraph
#' @import dplyr
#'
#' @references None.
#'
#' @seealso <https://www.github.com/travis-m-blimkie/networker>
#'
plot_network <- function(
  network,
  fill_column,
  layout       = "kk",
  edge_colour  = "grey50",
  edge_alpha   = 0.2,
  label        = FALSE,
  label_column,
  label_filter = 40,
  label_size   = 4,
  label_colour = "black",
  label_face   = "bold"
) {

  if (layout == "force_atlas") {
    message("Calculating Force Atlas node positions...")
    set.seed(1)
    layout_object <- ForceAtlas2::layout.forceatlas2(
      graph    = network,
      directed = FALSE,
      plotstep = 0
    )
  } else {
    layout_object <- layout
  }

  set_graph_style(foreground = "white")

  if (label) {

    network <- network %>%
      mutate(node_label = case_when(
        seed & (degree > label_filter) ~ {{label_column}},
        TRUE ~ NA_character_
      ))

    ggraph(network, layout = layout_object) +
      geom_edge_link(show.legend = FALSE, alpha = edge_alpha, colour = edge_colour) +
      geom_node_point(aes(size = degree, fill = {{fill_column}}), pch = 21) +
      scale_fill_gradient2(
        low   = "springgreen4",
        mid   = "white",
        high  = "firebrick",
        guide = "none"
      ) +
      geom_node_text(
        aes(label = node_label),
        size          = label_size,
        repel         = TRUE,
        colour        = label_colour,
        fontface      = label_face,
        check_overlap = TRUE
      ) +
      scale_size_continuous(range = c(3, 9), guide = "none") +
      theme(plot.margin = unit(rep(0.05, 4), "cm"))
  } else {
    ggraph(network, layout = layout_object) +
      geom_edge_link(show.legend = FALSE, alpha = edge_alpha, colour = edge_colour) +
      geom_node_point(aes(size = degree, fill = {{fill_column}}), pch = 21) +
      scale_fill_gradient2(
        low   = "springgreen4",
        mid   = "white",
        high  = "firebrick",
        guide = "none"
      ) +
      scale_size_continuous(range = c(3, 9), guide = "none") +
      theme(plot.margin = unit(rep(0.05, 4), "cm"))
  }
}