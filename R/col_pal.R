#' Color palettes
#'
#' Common interface to get color palettes from different packages:
#' viridisLite, colorRamps,
#' \href{"https://github.com/BlakeRMills/MetBrewer"}{MetBrewer},
#' href{"https://github.com/karthik/wesanderson"}{wesanderson}
#' and RColorBrewer.
#' Moreover three custom palettes: "custom", "dutch" and "spanish".
#'
#' See available palettes:
#' c(rev(ls("package:viridisLite"))[-c(1,2)],
#' rev(ls("package:colorRamps"))[-c(2,3)],
#' rownames(RColorBrewer::brewer.pal.info),
#' names(wesanderson::wes_palettes), "custom", "dutch", "spanish",
#' MetBrewer::MetPalettes)
#' Use scales::show_col() to plot color grid
#'
#' @param name name of the palette
#' @param n number of colors to return; may not work for every palette
#' @param nbrew number of color from brewer palettes
#' @param reverse reverse the palette
#'
#' @return a color palette as character vector
#' @export
#'
#' @examples
#' \dontrun{
#' }
col_pal <- function(name = "custom",
                    n = 100,
                    nbrew = NULL,
                    reverse = F) {


  #library(paletteer)
  #paletteer::paletteer_d("MetBrewer::VanGogh2")

  #scales::show_col(palette.colors(n = 9, palette = "Okabe-Ito", recycle = FALSE))
  #scales::show_col(paletteer::paletteer_d("colorblindr::OkabeIto"))

  '
  knitr::kable(dplyr::bind_rows(paletteer::palettes_c_names %>% dplyr::mutate(type2 = "continuous"),
                                paletteer::palettes_d_names %>% dplyr::mutate(type2 = "discrete")) %>%
                 dplyr::mutate(command = paste0(package, "::", palette)) %>%
                 dplyr::select(-c(novelty)),
               format = "jira")'

  #getOption("max.print") # change temporarily
  #options(max.print=3000)
  # order and color by discrete or continuous; sequential or divergent (color with crayon)
  dplyr::bind_rows(paletteer::palettes_c_names %>% dplyr::mutate(type2 = "continuous"),
                   paletteer::palettes_d_names %>% dplyr::mutate(type2 = "discrete")) %>%
    dplyr::mutate(command = paste0(package, "::", palette)) %>%
    dplyr::pull(command)


  # automate that choice (c vs. d) and in case of d, automatically make continuous if n > max_n
  colpal <- paletteer::paletteer_c("scico::berlin", n = 100)
  colpal <- paletteer::paletteer_d("nord::frost", n = 100, type = "continuous")


  scales::show_col(colpal)

  test2 <- paletteer::palettes_d_names %>% dplyr::mutate(type2 = "discrete")
  test <- paletteer::palettes_c_names %>% dplyr::mutate(type2 = "continuous")
  unique(test$type)
  ?grDevices::palette.colors()

  scl <- NULL

  if (any(grepl(name, names(MetBrewer::MetPalettes)))) {
    scl <- grDevices::colorRampPalette(MetBrewer::MetPalettes[[name]][[1]], interpolate = "linear")(n)
  }

  if (any(grepl(name, c(ls(loadNamespace("viridisLite")))))) {
    scl <- viridis::viridis_pal(option = name)(n)
  }

  if (any(grepl(name, ls(loadNamespace("colorRamps")), ignore.case = T))) {
    if (!isNamespaceLoaded("colorRamps")) {
      attachNamespace("colorRamps")
    }
    colfun <- match.fun(grep(name, ls(loadNamespace("colorRamps")), ignore.case = T, value = T)[1])
    scl <- colfun(n)
  }

  if (any(grepl(name, names(wesanderson::wes_palettes), ignore.case = T))) {
    scl <- as.character(wesanderson::wes_palette(grep(name, names(wesanderson::wes_palettes), ignore.case = T, value = T), n, type = "continuous"))
  }

  if (any(grepl(name, rownames(RColorBrewer::brewer.pal.info), ignore.case = T))) {
    name <- grep(name, rownames(RColorBrewer::brewer.pal.info), ignore.case = T, value = T)
    if (grepl("^spectral$", name, ignore.case = T)) {
      # needed so often, counter-intuitive by default
      reverse <- !reverse
    }
    if (is.null(nbrew) || nbrew > RColorBrewer::brewer.pal.info[name, "maxcolors"]) {
      nbrew <- RColorBrewer::brewer.pal.info[name, "maxcolors"]
    }
    scl <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(nbrew, name), interpolate = "linear")(n)
  }

  if (name %in% c("custom", "dutch", "spanish")) {
    scl <- switch(name,
                  custom = c("grey65", "darkgoldenrod1", "cornflowerblue", "forestgreen", "tomato2", "mediumpurple1", "turquoise3", "lightgreen", "navy", "plum1",
                             "red4", "khaki1", "tan4", "cadetblue1", "olivedrab3", "darkorange2", "burlywood2", "violetred3", "aquamarine3",
                             "grey30", "lavender", "yellow", "grey10", "pink3", "turquoise4", "darkkhaki", "magenta", "blue", "green", "blueviolet", "red",
                             "darkolivegreen", "orchid1", "springgreen", "dodgerblue4", "deepskyblue", "palevioletred4", "gold4", "maroon1", "lightyellow", "greenyellow", "purple4")[1:n],
                  dutch = c(
                    "#FFC312","#C4E538","#12CBC4","#FDA7DF","#ED4C67",
                    "#F79F1F","#A3CB38","#1289A7","#D980FA","#B53471",
                    "#EE5A24","#009432","#0652DD","#9980FA","#833471",
                    "#EA2027","#006266","#1B1464","#5758BB","#6F1E51")[1:n],
                  spanish = c(
                    "#40407a","#706fd3","#f7f1e3","#34ace0","#33d9b2",
                    "#2c2c54","#474787","#aaa69d","#227093","#218c74",
                    "#ff5252","#ff793f","#d1ccc0","#ffb142","#ffda79",
                    "#b33939","#cd6133","#84817a","#cc8e35","#ccae62")[1:n])
    scl <- scl[which(!is.na(scl))]
  }

  if (!missing(n) && length(scl) < n) {
    scl <- scales::hue_pal()(n)
  }

  if (is.null(scl)) {
    print("Palette name not found. Switching to spectral.")
    scl <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"), interpolate = "linear")(n)
  }

  if (reverse) {
    scl <- rev(scl)
  }

  return(scl)
}

'
# https://www.datanovia.com/en/blog/awesome-list-of-657-r-color-names/
# https://stackoverflow.com/questions/5392061/algorithm-to-check-similarity-of-colors
library(scales)
library(tidyverse)


tt <- col_pal("custom")
show_col(tt)

library(schemr)

r_color <- colors()
out <- as.data.frame(schemr::rgb_to_lab(t(col2rgb(r_color))))
rownames(out) <- r_color
r_color <- r_color[which(!grepl("^grey|^gray|black|white", r_color))]


euclidean <- function(x, mat) {
  sqrt(sum((as.numeric(mat[x[1],]) - as.numeric(mat[x[2],]))^2))
}
euclidean2 <- function(x, mat) {
  sqrt(sum((as.numeric(mat[x[1,1],]) - as.numeric(mat[x[1,2],]))^2))
}

for (i in 1:20) {
  # find max distant color
  res <- dplyr::bind_rows(pbapply::pblapply(split(expand.grid(setdiff(r_color, tt), tt, stringsAsFactors=F), 1:nrow(expand.grid(setdiff(r_color, tt), tt))), function(x) {
    data.frame(dist = euclidean2(x, mat = out), col1 = x[1,1], col2 = x[1,2])
  }))
  # min dist - best measure?!
  res2 <- res %>% dplyr::group_by(col1) %>% summarise(mean_dist = mean(dist), min_dist = min(dist))
  tt <- c(tt, res2 %>% dplyr::filter(min_dist == max(min_dist)) %>% dplyr::slice(1) %>% dplyr::pull(col1))
  show_col(tt)
}
# distance matrix
table <- data.frame(dist = unlist(pbapply::pblapply(combn(r_color[1:100], 2, simplify = F), mat = out, euclidean)))
table$col1 <- sapply(combn(r_color[1:100], 2, simplify = F), "[", 1)
table$col2 <- sapply(combn(r_color[1:100], 2, simplify = F), "[", 2)
library(tidyverse)
library(viridis)

ggplot(table %>% dplyr::filter(col1 != col2), aes(x=col1,y=col2,fill=dist)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank())

table <- data.frame(dist = unlist(pbapply::pblapply(split(expand.grid(tt, r_color), 1:nrow(expand.grid(tt, r_color))), mat = out, euclidean)))
'
