#' Color palettes
#'
#' Wrapper function around paletteer and two additional color palettes.
#' See paletteer::palettes_c_names and paletteer::palettes_d_names.
#'
#' @param name name of the palette
#' @param n number of colors to return; may not work for every palette
#' @param direction reverse palette with -1
#' @param contrast_filter remove colors if contrast to bg_color is below contrast_ratio_min.
#' @param contrast_ratio_min minimum ration
#' @param bg_color background hex color or r color name
#'
#' @return a color palette as character vector
#' @export
#'
#'
#' @importFrom magrittr "%>%"
#'
#' @examples
#' \dontrun{
#' col_pal(name = "custom", n = 10)
#' }
col_pal <- function(name = NULL,
                    n = NULL,
                    direction = c(1,-1),
                    contrast_filter = F,
                    contrast_ratio_min = 1.5,
                    bg_color = "white") {

  if (!requireNamespace("paletteer", quietly = T)) {
    utils::install.packages("paletteer")
  }

  paletteers <- dplyr::bind_rows(paletteer::palettes_c_names %>% dplyr::mutate(type2 = "continuous"),
                                 paletteer::palettes_d_names %>% dplyr::mutate(type2 = "discrete")) %>%
    dplyr::mutate(command = paste0(package, "::", palette))

  if (is.null(name)) {
    message("Select one palette by palette or command. Additional ones are 'custom', 'ggplot', 'hue', 'material'.")
    print(paletteers, n = 50)
    return(paletteers)
  }

  direction <- as.numeric(match.arg(as.character(direction), choices = c("1","-1")))

  if (name %in% c("ggplot", "ggplot2", "hue", "hue_pal", "huepal")) {
    if (is.null(n)) {
      n <- 100
    }
    pal_return <- prismatic::color(scales::hue_pal()(n))
    if (direction == -1) {
      pal_return <- rev(pal_return)
    }
  } else if (name == "custom") {
    pal_return <- prismatic::color(c("grey65", "darkgoldenrod1", "cornflowerblue", "forestgreen", "tomato2", "mediumpurple1", "turquoise3", "lightgreen", "navy", "plum1",
                                     "red4", "khaki1", "tan4", "cadetblue1", "olivedrab3", "darkorange2", "burlywood2", "violetred3", "aquamarine3",
                                     "grey30", "lavender", "yellow", "grey10", "pink3", "turquoise4", "darkkhaki", "magenta", "blue", "green", "blueviolet", "red",
                                     "darkolivegreen", "orchid1", "springgreen", "dodgerblue4", "deepskyblue", "palevioletred4", "gold4", "maroon1", "lightyellow", "greenyellow", "purple4"))
    if (!is.null(n)) {
      pal_return <- pal_return[1:n]
    }
    if (direction == -1) {
      pal_return <- rev(pal_return)
    }
  } else if (name == "material") {
    # bg_col: '#273238'
    pal_return <- prismatic::color(c("#CC6666", "#DF935F", "#81A3BE","#B5BD68", "#707880", "#B394BB", "#4F5C69", "#CC0000", "#F1CC37",
                                     "#5787DA", "#25B876", "#5E3582", "#93C5DE", "#B7A7C5", "#2E90B0", "#C5C8C6", "#50C186", "#34312F"))

    if (!is.null(n)) {
      pal_return <- pal_return[1:n]
    }
    if (direction == -1) {
      pal_return <- rev(pal_return)
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
      print(pal_select)
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

  if (contrast_filter) {
    ratios <- purrr::map_dbl(setNames(pal_return, pal_return), contrast_ratio, color2 = bg_color)
    ratios <- ratios[which(ratios > contrast_ratio_min)]
    pal_return <- prismatic::color(names(ratios))
  }


  return(pal_return)
}


hex_to_linear_rgb <- function(hex) {
  rgb <- col2rgb(hex) / 255
  linear_rgb <- ifelse(rgb <= 0.03928,
                       rgb / 12.92,
                       ((rgb + 0.055) / 1.055) ^ 2.4)
  return(linear_rgb)
}


relative_luminance <- function(hex) {
  linear_rgb <- hex_to_linear_rgb(hex)
  # Apply luminance formula
  luminance <- 0.2126 * linear_rgb[1] +
    0.7152 * linear_rgb[2] +
    0.0722 * linear_rgb[3]
  return(luminance)
}


contrast_ratio <- function(color1, color2) {
  L1 <- relative_luminance(color1)
  L2 <- relative_luminance(color2)
  if (L1 < L2) { temp <- L1; L1 <- L2; L2 <- temp }
  ratio <- (L1 + 0.05) / (L2 + 0.05)
  return(ratio)
}

#purrr::map_dbl(setNames(col_pal("material"), col_pal("material")), contrast_ratio, color2 = "#273238")
#purrr::map_dbl(setNames(col_pal("custom"), col_pal("custom")), contrast_ratio, color2 = "white")

