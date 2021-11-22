#' Split a flowjo export of multiple plots into separate graphics for rearrangement
#'
#' @param img character path to an image, preferable a png exported from flowjos layout editor, has to be equally spaced in x and y direction
#' @param nPlots_x integer, how many plots in x direction
#' @param nPlots_y integer, how many plots in y direction
#' @param folder folder name of where to save split images, this folder is written to the directory of img
#' @param pptx_name character, name a pptx to write the split images to; if missing no pptx is generated
#' @param pptx_image_size numeric, size of the images in pptx
#' @param pptx_border_space numeric, space between images in pptx
#'
#' @return writes a folder with split images and optionally a pptx
#' @export
#'
#' @examples
#' \dontrun{
#' # copy an example file to your disk
#' file.copy(from = system.file("extdata", "img.png", package = "fcexpr"), to = "your_path")
#' fcexpr::split_flowjo_export_image(img = "your_path/img.png",
#' nPlots_x = 6, nPlots_y = 6, pptx = "my_img_split.pptx")
#' }
split_flowjo_export_image <- function(img,
                                      nPlots_x = 4,
                                      nPlots_y = 3,
                                      folder,
                                      pptx_name,
                                      pptx_image_size = 1,
                                      pptx_border_space = 0.2) {

  if (requireNamespace("magick", quietly = T)){
    utils::install.packages("magick")
  }

  if (missing(folder)) {
    folder <- format(Sys.time(), "%Y%m%d_%H%M%S")
  }

  parent_folder <- dirname(img)
  dir.create(file.path(parent_folder, folder), recursive = T, showWarnings = F)
  img <-  magick::image_trim(magick::image_read(img))
  x.incr <- ceiling(magick::image_info(img)$width/nPlots_x)
  y.incr <- ceiling(magick::image_info(img)$height/nPlots_y)


  for (x in seq(1,nPlots_x,1)) {
    for (y in seq(1,nPlots_y,1)) {
      geom <- paste0(x.incr, "x", y.incr, "+", (x-1)*x.incr, "+", (y-1)*y.incr)
      img.crop <- magick::image_crop(img, geometry = geom)
      tryCatch({
        img.crop.trim <- magick::image_trim(img.crop)
        magick::image_write(img.crop.trim, path = file.path(parent_folder, folder, paste0("x", x, "_", "y", y, ".png")), format = "png")
      }, error = function(e) {
        print(paste0("Whitespace found at position x = ", x, ", y = ", y, "."))
      })
    }
  }
  print(paste0("Images saved to ", file.path(parent_folder, folder)))


  if (!missing(pptx_name)) {
    if (!"officer" %in% rownames(utils::installed.packages())) {utils::install.packages("officer")}
    pptx <- officer::read_pptx()
    pptx <- officer::add_slide(pptx, layout = "Title and Content", master = "Office Theme")

    for (i in list.files(file.path(parent_folder, folder), full.names = T)) {

      x.coord <- as.numeric(stringr::str_extract(stringr::str_extract(i, "\\/x[:digit:]{1,}"), "[:digit:]{1,}"))
      y.coord <- as.numeric(stringr::str_extract(stringr::str_extract(i, "_y[:digit:]{1,}\\."), "[:digit:]{1,}"))

      pptx <- officer::ph_with(pptx, value = officer::external_img(i, width = 0.5, height = 0.5),
                               location = officer::ph_location(left = x.coord + pptx_border_space,
                                                               top = y.coord + pptx_border_space,
                                                               width = pptx_image_size,
                                                               height = pptx_image_size))
    }

    if (rev(strsplit(pptx_name, "\\.")[[1]]) != "pptx") {
      pptx_name <- paste0(pptx_name, ".pptx")
    }
    print(pptx, target = file.path(parent_folder, folder, pptx_name))
    print(paste0("pptx saved as ", file.path(parent_folder, folder, pptx_name)))

  }
}
