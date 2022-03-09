#' Plot an undirected PPI network using ggraph
#'
#' @param network `tidygraph` object, output from `build_network`
#' @param fill_column Tidy-select column for mapping node colour. The colour
#'   scale set up is green-white-red, designed to map to fold change values
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
#' @param label_size Size of node labels, defaults to 4.
#' @param label_colour Colour of node labels, defaults to "black"
#' @param label_face Font face for node labels, defaults to "bold"
#' @param label_padding Padding around the label, defaults to 0.25 lines.
#' @param min_seg_length Minimum length of lines to be drawn from labels to
#'   points. The default specified here is 0.25, half of the normal default
#'   value.
#' @param seed Number used in call to `set.seed()` to allow for reproducible
#'   network generation. Can be changed to get slightly different layouts.
#'
#' @return An object of class "gg"
#' @export
#'
#' @import ggraph
#' @import dplyr
#'
#' @details Any layout supported by ggraph can be specified here - see
#'   `?layout_tbl_graph_igraph` for a list of options. Additionally, there is
#'   support for the "force_atlas" method, implemented via the ForceAtlas2
#'   package
#'
#' @references See <https://github.com/analyxcompany/ForceAtlas2> for details on
#'   this method.
#'
#' @seealso <https://github.com/travis-m-blimkie/networker>
#'
plot_network <- function(
  network,
  fill_column,
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
  label_face     = "bold",
  label_padding  = 0.25,
  min_seg_length = 0.25,
  seed           = 1
) {

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

  set_graph_style(foreground = "white")

  if (is.numeric(pull(network, {{fill_column}}))) {
    network_fill_geom <- scale_fill_gradient2(
      low      = "springgreen4",
      mid      = "white",
      high     = "firebrick",
      na.value = int_colour,
      guide    = ifelse(legend, "colourbar", "none")
    )
  } else if (is.character(pull(network, {{fill_column}}))) {
    network_fill_geom <- scale_fill_brewer(
      palette  = "Set1",
      na.value = int_colour,
      guide    = ifelse(legend, "legend", "none")
    )
  }

  if (label) {

    network <- network %>%
      mutate(node_label = case_when(
        seed & (degree > label_filter) ~ {{label_column}},
        TRUE ~ NA_character_
      ))

    ggraph(network, layout = layout_object) +
      geom_edge_link(show.legend = FALSE, alpha = edge_alpha, colour = edge_colour) +
      geom_node_point(aes(size = degree, fill = {{fill_column}}), pch = 21, colour = "grey70") +
      network_fill_geom +
      geom_node_text(
        aes(label = node_label),
        size          = label_size,
        repel         = TRUE,
        family        = fontfamily,
        colour        = label_colour,
        fontface      = label_face,
        check_overlap = TRUE,
        box.padding   = label_padding,
        min.segment.length = min_seg_length
      ) +
      scale_size_continuous(range = node_size, guide = "none") +
      theme(
        text         = element_text(family = fontfamily),
        plot.margin  = unit(rep(0.05, 4), "cm"),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14)
      )
  } else {
    ggraph(network, layout = layout_object) +
      geom_edge_link(show.legend = FALSE, alpha = edge_alpha, colour = edge_colour) +
      geom_node_point(aes(size = degree, fill = {{fill_column}}), pch = 21) +
      network_fill_geom +
      scale_size_continuous(range = node_size, guide = "none") +
      theme(
        text         = element_text(family = fontfamily),
        plot.margin  = unit(rep(0.05, 4), "cm"),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14)
      )
  }
}
