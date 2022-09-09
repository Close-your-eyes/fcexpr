#' Plot gating tree as a graph to explore how populations are gated in a simple way
#'
#'
#'
#' @param PopulationFullPath vector of all populations path (full) gated in a (flowjo) workspace;
#' optionally provide named vector (e.g. with short names) which will be plotted in the graph at correponding nodes
#' @param layout layout for plotting, see ?ggraph::ggraph or ?ggraph::layout_tbl_graph_igraph
#' @param find_short_gating_path if PopulationFullPath is a vector without names find the shortest unique
#' gating path from these full paths for plotting
#' @param ... additional argument to ggrepel, like max.overlaps
#'
#' @return list of plot, graph and data.frame
#' @export
#'
#' @examples
#' \dontrun{
#' # Get count data frame from flowjo wsp
#' ps_counts <- fcexpr::wsx_get_popstats(ws = ws)[["counts"]]
#' # check if gatingtree is unique
#' dplyr::group_by(ps_counts, PopulationFullPath) %>% dplyr::count()
#' # generate graph
#' graph_list <- fcexpr::gating_tree_plot(PopulationFullPath = unique(ps_counts$PopulationFullPath))
#' # manual plotting graph to modify (e.g. layout)
#' plot <- ggraph::ggraph(graph_list[["graph"]], layout = "tree") +
#' ggraph::geom_edge_link() +
#' ggraph::geom_node_point(size = 4) +
#' ggraph::geom_node_label(ggplot2::aes(label = Population), repel = T)
#'
#' ## color nodes by a lineage marker
#'igraph::V(graph_list[["graph"]])$subname <- dplyr::case_when(grepl("CD3\\+", igraph::V(graph_list[["graph"]])$name) &
#'                                                          !grepl("CD8\\+|CD4\\+", igraph::V(graph_list[["graph"]])$name) ~ "CD3+",
#'                                                          grepl("CD8\\+", igraph::V(graph_list[["graph"]])$name) ~ "CD8+",
#'                                                          grepl("CD4\\+", igraph::V(graph_list[["graph"]])$name) ~ "CD4+",
#'                                                          grepl("CD19\\+", igraph::V(graph_list[["graph"]])$name) ~ "CD19+",
#'                                                          grepl("CD66b\\+", igraph::V(graph_list[["graph"]])$name) ~ "CD66b+",
#'                                                          grepl("CD11c\\+", igraph::V(graph_list[["graph"]])$name) ~ "CD11c+",
#'                                                          grepl("CD56\\+|NK", igraph::V(graph_list[["graph"]])$name) ~ "NK")
#' }
gating_tree_plot <- function(PopulationFullPath,
                             layout = "tree",
                             find_short_gating_path = T,
                             ...) {

  if (!requireNamespace("ggraph", quietly = T)) {
    utils::install.packages("ggraph")
  }
  if (!requireNamespace("igraph", quietly = T)) {
    utils::install.packages("igraph")
  }

  # make unique but keep names if present
  if (is.null(names(PopulationFullPath))) {
    PopulationFullPath <- unique(PopulationFullPath)
  } else {
    if(!all(dplyr::pull(dplyr::count(dplyr::group_by(dplyr::distinct(utils::stack(PopulationFullPath)), values)), n) == 1)) {
      stop("Population (short names) for PopulationFullPaths are not unique. Please fix or do not provide names.")
    }
    PopulationFullPath_un <- unique(PopulationFullPath)
    names(PopulationFullPath_un) <- stats::setNames(names(PopulationFullPath), PopulationFullPath)[PopulationFullPath_un]
    PopulationFullPath <- PopulationFullPath_un
  }

  from_to_df <- dplyr::mutate(data.frame(path = gsub("root/root", "root", paste0("root/", PopulationFullPath))), level = stringr::str_count(path, "/"))
  all_lev <- unique(from_to_df$level)

  from_to_df2 <- dplyr::bind_rows(purrr::keep(purrr::map(all_lev[-length(all_lev)], function(x) {
    out <- purrr::map(dplyr::distinct(from_to_df[which(from_to_df$level == x),])[,"path",drop=T], function(y) {
      # avoid regex which will have a problem with "+" and other special characters; hence check with str_locate is match starts at position 1
      to_df <- dplyr::distinct(from_to_df[intersect(which(from_to_df$level == (x+1)), which(stringr::str_locate(from_to_df$path, stringr::fixed(y))[,"start"] == 1)),])
      if (nrow(to_df) > 0) {
        return(data.frame(from = y, to = to_df$path))
      } else {
        return(NULL)
      }
    })
    out <- purrr::keep(out, purrr::negate(is.null))
    if (length(out) > 0) {
      return(dplyr::bind_rows(out))
    } else {
      return(NULL)
    }
  }), purrr::negate(is.null)))

  from_to_df3 <- dplyr::mutate(from_to_df2, to =  gsub("root/", "", to), from =  gsub("root/", "", from))
  graph <- igraph::graph_from_data_frame(d = from_to_df3, directed = T)

  if (is.null(names(PopulationFullPath))) {
    if (find_short_gating_path) {
      igraph::V(graph)$Population <- fcexpr:::shortest_unique_path(igraph::V(graph)$name)
    } else {
      igraph::V(graph)$Population <- igraph::V(graph)$name
    }
  } else {
    igraph::V(graph)$Population <- stats::setNames(names(PopulationFullPath), PopulationFullPath)[igraph::V(graph)$name]
  }

  plot <- ggraph::ggraph(graph, layout = layout) +
    ggraph::geom_edge_link() +
    ggraph::geom_node_point(size = 4) +
    ggraph::geom_node_label(ggplot2::aes(label = Population), repel = T, ...)

  return(list(plot = plot, graph = graph, df = from_to_df3))
}
