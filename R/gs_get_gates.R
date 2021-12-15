#' Gates from gatingset for plotting with ggcyto
#'
#' @param gs gatingset
#' @param bins number of bins for plotting
#' @param low.lim lower range limit
#' @param up.lim upper range limit
#' @param up.lim.scatter upper range limit for scatter channels
#' @param gate.pct.stat.x.pos x-position of percent labels for gates
#' @param gate.pct.stat.y.pos y-position of percent labels for gates
#' @param gate.pct.stat.size size of percent labels for gates
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' gates <- gs_get_gates(gs = my_gs)
#'
#' out <- purrr::flatten(lapply(unique(gates$gate.level), function (z) {
#'g <- gates[which(gates[,"gate.level"] == z),]
#'# split into groups of equal gate.level, subset, xy-plane
#'g <- split(g, paste(g$gate.level, g$subset, g$x, g$y, sep = "__"))
#'p <- lapply(g, function(gg) {
#'  my.filter <- if (gg[1,"marginalFilter"]) {marginalFilter} else {NULL}
#'  p <- ggcyto(gs, subset = gg[1,"subset"], filter = my.filter, aes(!!sym(gg[1,"x"]), !!sym(gg[1,"y"])), max_nrow_to_plot = 5e4) +
#'    geom_hex(bins = gg[1,"bins"]) +
#'    theme_bw() +
#'    xlab(gg[1,"x.lab"]) +
#'    ylab(gg[1,"y.lab"]) +
#'    ggcyto_par_set(limits = list(x = c(gg[1,"x.low.lim"], gg[1,"x.up.lim"]), y = c(gg[1,"y.low.lim"], gg[1,"y.up.lim"]))) +
#'    scale_fill_gradientn(colours = scexpr::col_pal("spectral"), trans = "pseudo_log") +
#'    theme(legend.position = "none", strip.background = element_rect(fill = "white"), text = element_text(family = "Courier"), panel.grid = element_blank(), axis.text = element_blank()) +
#'    facet_grid(cols = vars(Patient, TCR), rows = vars(stimulus)) # order inline: facet_grid(cols = vars(PBMC.donor.short, factor(IL15.pre.stim.conc.ug.ml, levels=c("0", "0.12", "0.37", "1.11", "3.33", "10"))), rows = vars(RPTECs, RPTEC.IFNg.pre.stim))
  ## loop through multiple lines of gg for gates
#'  for (i in 1:nrow(gg)) {
#'    p <-
#'      p +
#'      geom_gate(gg[i,"gate.path.full"], colour = "black")
#'    #geom_stats(gg[i,"gate.path.full"], type = "percent", size = gg[i,"gate.pct.stat.size"], color = "black", digits = 1, adjust = c(gg[i,"gate.pct.stat.x.pos"], gg[i,"gate.pct.stat.y.pos"]), fill = alpha(c("white"),0.5))
#'  }
#'  return(p)
#'})
#'return(p)
#'}))
#'
#'
#' }
gs_get_gates <- function(gs,
                         bins = 200,
                         low.lim = 0,
                         up.lim = 300,
                         up.lim.scatter = 250000,
                         gate.pct.stat.x.pos = 0.8,
                         gate.pct.stat.y.pos = 0.2,
                         gate.pct.stat.size = 4) {

  if (!requireNamespace("flowWorkspace", quietly = T)){
    utils::install.packages("flowWorkspace")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }

  gates <-
    data.frame(gate.path.full = flowWorkspace::gs_get_pop_paths(gs),
               gate.path.auto = flowWorkspace::gs_get_pop_paths(gs, path = "auto"),
               gate.level = nchar(flowWorkspace::gs_get_pop_paths(gs)) - nchar(gsub("/", "", flowWorkspace::gs_get_pop_paths(gs)))) %>%
    dplyr::filter(gate.level > 0) %>%
    dplyr::mutate(subset = gsub("^/", "", dirname(gate.path.full))) %>%
    dplyr::mutate(subset = ifelse(subset == "", "root", subset)) %>%
    dplyr::mutate(bins = bins, gate.pct.stat.x.pos = gate.pct.stat.x.pos, gate.pct.stat.y.pos = gate.pct.stat.y.pos, gate.pct.stat.size = gate.pct.stat.size)

  gates$dims <- sapply(gates$gate.path.full, function(x) {list(unname(flowCore::parameters(flowWorkspace::gs_pop_get_gate(gs[[1]], x)[[1]])))})
  gates$x <- unname(sapply(gates$dims, function(x) {unlist(x)[1]}))
  gates$y <- unname(sapply(gates$dims, function(x) {unlist(x)[2]}))
  gates$x.lab <- unname(sapply(gates$dims, function(x) {unlist(x)[1]}))
  gates$y.lab <- unname(sapply(gates$dims, function(x) {unlist(x)[2]}))

  gates <-
    gates %>%
    dplyr::mutate(marginalFilter = ifelse(grepl("fsc|ssc", x, ignore.case = T) & grepl("fsc|ssc", y, ignore.case = T), T, F)) %>%
    dplyr::mutate(x.low.lim = low.lim, x.up.lim = ifelse(grepl("fsc|ssc", x, ignore.case = T), up.lim.scatter, up.lim), y.low.lim = low.lim, y.up.lim = ifelse(grepl("fsc|ssc", y, ignore.case = T), up.lim.scatter, up.lim))

  gates$facet.strip <- c(T, rep(F, nrow(gates)-1))

  return(gates)
}



