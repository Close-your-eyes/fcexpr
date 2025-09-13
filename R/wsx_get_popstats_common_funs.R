get_group_df <- function(ws, groups, invert_groups) {
  group_df <- wsx_get_groups(ws, collapse = NULL)
  if (is.null(groups)) {
    groups <- unique(group_df[,"FlowJoGroup", drop=T])
  } else {
    groups <- unique(groups)
  }
  if (invert_groups) {
    group_df <- group_df[which(!group_df[,"FlowJoGroup", drop=T] %in% groups),]
  } else {
    group_df <- group_df[which(group_df[,"FlowJoGroup", drop=T] %in% groups),]
  }
  if (nrow(group_df) == 0) {
    stop("Non of provided groups found.")
  }
  # if any sample is in at least two groups, the group column becomes a list
  if (any(table(group_df$sampleID) > 1)) {
    # or nest with tidyr?
    group_df <- dplyr::summarise(dplyr::group_by(group_df, sampleID), FlowJoGroup = list(FlowJoGroup))
  }
  return(group_df)
}

get_sample_nodes <- function(ws, group_df) {
  # each sample which may be in multiple groups is only considered once
  ids <- unique(group_df[,"sampleID",drop=T])
  samples <- xml2::xml_children(xml2::xml_child(ws, "SampleList"))
  sample_ids <- which(sapply(xml2::xml_attrs(xml2::xml_child(samples, "DataSet")), "[[", "sampleID") %in% ids)
  if (anyDuplicated(basename(sapply(xml2::xml_attrs(xml2::xml_child(samples, "DataSet")), "[[", "uri"))) != 0) {
    message("Duplicate filenames detected in workspace:")
    x <- basename(sapply(xml2::xml_attrs(xml2::xml_child(samples, "DataSet")), "[[", "uri"))
    message(paste(x[which(duplicated(x))], collapse = ", "))
  }
  samples <- samples[sample_ids]
  return(samples)
}


check_filenames_hdd <- function(pop_df) {
  gates_paths <-
    pop_df |>
    dplyr::select(FilePath, FileName) |>
    dplyr::mutate(FileDir = dirname(FilePath)) |>
    dplyr::mutate(dir_exist = dir.exists(FileDir)) |>
    dplyr::mutate(file_exist = file.exists(FilePath)) |>
    dplyr::filter(dir_exist & !file_exist)
  if (nrow(gates_paths) > 0) {
    message("Some FileNames as in FlowJo not found on disk. Did you change fcs filenames on disk but did not
            open FlowJo? Do so, to update FileNames in FlowJo to as they are on disk. Otherwise joining
            with samplediscription may fail. Or join by identity only: ")
    message("dplyr::left_join(counts, sd, by = c('identity'))")
  }
}

add_identity <- function(pop_df, kwl) {
  fcs_idents <-
    stack(fcexpr:::get_fcs_identities(kwl)) |>
    dplyr::rename(identity = "values", "FileName" = ind) |>
    dplyr::mutate(FileName = as.character(FileName))
  return(dplyr::left_join(pop_df, fcs_idents, by  = "FileName"))
}

