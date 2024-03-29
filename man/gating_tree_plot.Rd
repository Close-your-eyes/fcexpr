% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gating_tree_plot.R
\name{gating_tree_plot}
\alias{gating_tree_plot}
\title{Plot gating tree as a graph to explore how populations are gated in a simple way}
\usage{
gating_tree_plot(
  PopulationFullPath,
  layout = "tree",
  find_short_gating_path = T,
  ...
)
}
\arguments{
\item{PopulationFullPath}{vector of all populations path (full) gated in a (flowjo) workspace;
optionally provide named vector (e.g. with short names) which will be plotted in the graph at correponding nodes}

\item{layout}{layout for plotting, see ?ggraph::ggraph or ?ggraph::layout_tbl_graph_igraph}

\item{find_short_gating_path}{if PopulationFullPath is a vector without names find the shortest unique
gating path from these full paths for plotting}

\item{...}{additional argument to ggrepel, like max.overlaps}
}
\value{
list of plot, graph and data.frame
}
\description{
Plot gating tree as a graph to explore how populations are gated in a simple way
}
\examples{
\dontrun{
# Get count data frame from flowjo wsp
ps_counts <- fcexpr::wsx_get_popstats(ws = ws)[["counts"]]
# check if gatingtree is unique
dplyr::group_by(ps_counts, PopulationFullPath) \%>\% dplyr::count()
# generate graph
graph_list <- fcexpr::gating_tree_plot(PopulationFullPath = unique(ps_counts$PopulationFullPath))
# manual plotting graph to modify (e.g. layout)
plot <- ggraph::ggraph(graph_list[["graph"]], layout = "tree") +
ggraph::geom_edge_link() +
ggraph::geom_node_point(size = 4) +
ggraph::geom_node_label(ggplot2::aes(label = Population), repel = T)

## color nodes by a lineage marker
igraph::V(graph_list[["graph"]])$subname <- dplyr::case_when(grepl("CD3\\\\+", igraph::V(graph_list[["graph"]])$name) &
                                                         !grepl("CD8\\\\+|CD4\\\\+", igraph::V(graph_list[["graph"]])$name) ~ "CD3+",
                                                         grepl("CD8\\\\+", igraph::V(graph_list[["graph"]])$name) ~ "CD8+",
                                                         grepl("CD4\\\\+", igraph::V(graph_list[["graph"]])$name) ~ "CD4+",
                                                         grepl("CD19\\\\+", igraph::V(graph_list[["graph"]])$name) ~ "CD19+",
                                                         grepl("CD66b\\\\+", igraph::V(graph_list[["graph"]])$name) ~ "CD66b+",
                                                         grepl("CD11c\\\\+", igraph::V(graph_list[["graph"]])$name) ~ "CD11c+",
                                                         grepl("CD56\\\\+|NK", igraph::V(graph_list[["graph"]])$name) ~ "NK")
}
}
