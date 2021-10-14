compMat_to_FCS <- function(fcs.file.paths, compMat) {

  comp.mat <- read.csv(comp.mat.csv.path, header = T, row.names = 1)
  rownames(comp.mat) <- gsub("_", "/", rownames(comp.mat))
  colnames(comp.mat) <- rownames(comp.mat)

  for (i in fcs.file.paths) {
    print(i)
    ff <- flowCore::read.FCS(i, truncate_max_range = F, emptyValue = F)
    sp <- flowCore::keyword(ff)[["SPILL"]]
    rownames(sp) <- colnames(sp)

    if (!all(colnames(comp.mat) %in% colnames(sp))) {
      stop("Not all colnames of comp.mat match with those of the SPILLOVER keyword matrix from the FCS file.")
    }

    if (!all(rownames(comp.mat) %in% rownames(sp))) {
      stop("Not all rownames of comp.mat match with those of the SPILLOVER keyword matrix from the FCS file.")
    }

    for (r in rownames(comp.mat)) {
      for (c in colnames(comp.mat)) {
        sp[which(rownames(sp) == r), which(colnames(sp) == c)] <- comp.mat[which(rownames(comp.mat) == r), which(colnames(comp.mat) == c)]
      }
    }


    rownames(sp) <- NULL
    flowCore::keyword(ff)[["SPILL"]] <- sp
    flowCore::write.FCS(ff, i)
  }
}
