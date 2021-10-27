split_flowjo_export_image <- function(img,
                                      n.plots.x = 4,
                                      n.plots.y = 3,
                                      export.folder.name = "split_img",
                                      pptx.name,
                                      pptx.image.size = 1,
                                      pptx.border.space = 0.2) {


  if ("magick" %in% rownames(installed.packages())) {install.packages(x)}

  img <- magick::image_read(img)

  img <- magick::image_trim(img)
  x.incr <- ceiling(magick::image_info(img)$width/n.plots.x)
  y.incr <- ceiling(magick::image_info(img)$height/n.plots.y)

  dir.create(file.path(dirname(img), export.folder.name), recursive = T, showWarnings = F)

  out <- lapply(1:n.plots.x, function(x) {
    lapply(1:n.plots.y, function(y) {
      geom <- paste0(x.incr, "x", y.incr, "+", (x-1)*x.incr, "+", (y-1)*y.incr)
      img.crop <- magick::image_crop(img, geometry = geom)

      tryCatch({
        img.crop.trim <- magick::image_trim(img.crop)
        magick::image_write(img.crop.trim, path = file.path(wd, export.folder.name, paste0("x", x, "_", "y", y, ".png")), format = "png")
      }, error = function(e) {
        print(paste0("Whitespace found at position x = ", x, ", y = ", y, "."))
      })

      return(NULL)
    })
    return(NULL)
  })
  print(paste0("Images saved to ", file.path(wd, export.folder.name)))


  if (!missing(pptx.name)) {
    if ("officer" %in% rownames(installed.packages())) {install.packages(x)}
    pptx <-
      read_pptx() %>%
      add_slide(layout = "Title and Content", master = "Office Theme")

    for (i in list.files(export.folder.name, full.names = T)) {

      x.coord <- as.numeric(stringr::str_extract(stringr::str_extract(i, "\\/x[:digit:]{1,}"), "[:digit:]{1,}"))
      y.coord <- as.numeric(stringr::str_extract(stringr::str_extract(i, "_y[:digit:]{1,}\\."), "[:digit:]{1,}"))

      pptx <-
        pptx %>%
        ph_with(value = external_img(i, width = 0.5, height = 0.5),
                location = ph_location(left = x.coord + pptx.border.space,
                                       top = y.coord + pptx.border.space,
                                       width = pptx.image.size,
                                       height = pptx.image.size))
    }

    print(pptx, target = file.path(wd, export.folder.name, pptx.name))
    print(paste0("pptx saved as ", file.path(wd, export.folder.name, pptx.name)))

  }
}
