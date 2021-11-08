wsx_compMats_to_fcs <- function(ws) {

  ws <- fcexpr:::check_ws(ws)

  ss <- xml2::xml_find_all(xml2::xml_child(ws, "SampleList"), "Sample")
  compMats <- lapply(seq_along(ss), function(n) {
    sp <- xml2::xml_child(ss[[n]], "transforms:spilloverMatrix")
    if (!is.na(sp)) {
      compMat <- t(do.call(cbind, lapply(xml2::xml_children(sp)[2:length(xml2::xml_children(sp))], function(x) {
        mat <- as.matrix(stats::setNames(as.numeric(xml_attr(xml2::xml_children(x), "value")), xml_attr(xml2::xml_children(x), "parameter")), nrow = 1, byrow = T)
        colnames(mat) <- xml2::xml_attr(x, "parameter")
        return(mat)
      })))
    } else {
      compMat <- NULL
    }
    return(compMat)
  })
  names(compMats) <- sapply(seq_along(ss), function(n) {
    gsub("^file:", "", xml2::xml_attr(xml2::xml_child(ss[[n]], "DataSet"), "uri"))
  })
  compMats <- compMats[which(!sapply(compMats, is.null))]

  ## read FCS and write compMat to SPILL keyword

}
n<-43
ws <- "/Users/vonskopnik/Downloads/20-Oct-2021.wsp"

x<-xml2::xml_children(sp)[2:length(xml2::xml_children(sp))][1]
