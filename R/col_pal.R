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
#' @param direction reverse palette with -1
#'
#' @return a color palette as character vector
#' @export
#'
#' @examples
#' \dontrun{
#' }
col_pal <- function(name = NULL,
                    n = NULL,
                    direction = c(1,-1)) {

  if (!requireNamespace("paletteer", quietly = T)) {
    utils::install.packages("paletteer")
  }

  paletteers <- dplyr::bind_rows(paletteer::palettes_c_names %>% dplyr::mutate(type2 = "continuous"),
                                 paletteer::palettes_d_names %>% dplyr::mutate(type2 = "discrete")) %>%
    dplyr::mutate(command = paste0(package, "::", palette))

  if (is.null(name)) {
    message("Select one palette by palette or command. Additional ones are 'custom' and 'ggplot' or 'hue'.")
    return(paletteers)
  }

  direction <- match.arg(direction, choices = c(1,-1))

  if (name %in% c("ggplot", "ggplot2", "hue", "hue_pal", "huepal")) {
    pal_select <- prismatic::color(scales::hue_pal()(n))
    if (direction == -1) {
      pal_select <- rev(pal_select)
    }
  } else if (name == "custom") {
    pal_select <- c("grey65", "darkgoldenrod1", "cornflowerblue", "forestgreen", "tomato2", "mediumpurple1", "turquoise3", "lightgreen", "navy", "plum1",
                    "red4", "khaki1", "tan4", "cadetblue1", "olivedrab3", "darkorange2", "burlywood2", "violetred3", "aquamarine3",
                    "grey30", "lavender", "yellow", "grey10", "pink3", "turquoise4", "darkkhaki", "magenta", "blue", "green", "blueviolet", "red",
                    "darkolivegreen", "orchid1", "springgreen", "dodgerblue4", "deepskyblue", "palevioletred4", "gold4", "maroon1", "lightyellow", "greenyellow", "purple4")[1:n]
    if (direction == -1) {
      pal_select <- rev(pal_select)
    }
  } else {
    if (grepl("::", name)) {
      pal_select <- paletteers %>% dplyr::filter(tolower(command) == tolower(name))
    } else {
      pal_select <- paletteers %>% dplyr::filter(tolower(palette) == tolower(name))
    }

    if (nrow(pal_select) == 0) {
      stop("Palette not found.")
    } else if (nrow(pal_select) > 1) {
      stop("Name is ambiguous. Please specify by command.")
    }

    if (pal_select$type2 == "discrete") {
      type <- "discrete"
      if (is.null(n)) {
        n <- pal_select$length
      }
      if (n > pal_select$length) {
        message("n = ", n, " larger than number of discrete color in palette (", pal_select$length, "). Going to interpolate to provide ", n, " colors.")
        type <- "continuous"
      }
      pal_return <- paletteer::paletteer_d(pal_select$command, n = n, type = type, direction = direction)
    } else if (pal_select$type2 == "continuous") {
      if (is.null(n)) {
        n <- 100
      }
      pal_return <- paletteer::paletteer_c(pal_select$command, n = n, direction = direction)
    }
  }
  return(pal_return)
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
