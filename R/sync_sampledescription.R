#' Synchronize a meta data table (sampledescription) with FCS files
#'
#' When running this function the rows in the 'file.name' (e.g. 'sampledescription.xlsx') and FCS files
#' in the FCS.files.folder are synchronized. When new FCS files are added the sampledescription is appended.
#' When entries in the FileName column of sampledescription are changed, FCS files are renamed accordingly.
#' Always, sampledescription is written to the parent folder of FCS.file.folder.
#' See Details below of the cases that can be handled.
#'
#' Five cases can be handled: (i) if no 'file.name'-file exists, one will be initiated based on
#' FCS files in 'FCS.file.folder'. (ii) When new FCS files are added to 'FCS.file.folder'
#' these are added in order of acquisition to the 'file.name'-file. (iii) When file names
#' in the FileName column of the 'file.name'-file are altered, FCS files in FCS.file.folder
#' are renamed accordingly. (iv) If FCS files are to be excluded or removed, the entry in the FileName column has to be left
#' blank and the function has to be called. (v) When the order of rows in the xlsx-file
#' is changed prefixes will be re-numbered.
#'
#' Preferentially, do not have the file.name'-file open in another program when calling the the function.
#' Never edit the 'identity' column in the 'file.name'-file manually.
#' Never mix up rows of FileName and identity.
#' The identity column contains a concatenated string of the $FIL keyword from FCS files, the number of events ($TOT)
#' and the acquisition date time of the FCS file.
#'
#' @param FCS.file.folder path to the folder which contains FCS files
#' @param file.name name of the sampledescription file, one of the following file types: .xlsx, .ods, .txt, .tsv, .csv
#' @param exclude.folders character vector of folders to exclude when checking for FCS files
#' @param write.log write a hidden (not hidden on windows) log file every time changes take place
#'
#' @return No return value. Instead sampledescription table and FCS files are synchronized.
#' @export
#'
#' @examples
#' \dontrun{
#' sync_sampledescription(FCS.file.folder = file.path(wd, 'FCS_files'))
#' }
sync_sampledescription <- function(FCS.file.folder,
                                   file.name = "sampledescription.xlsx",
                                   exclude.folders = c("compensation", "other_fcs_files", "experiment.file", "deleted_fcs_files",
                                                       "8_peak_bead", "rainbow_bead", "8_peak_beads", "rainbow_beads"),
                                   write.log = T) {

  file.suffix <- tools::file_ext(file.name)
  if (!file.suffix %in% c("xlsx")) { # c("xlsx", "ods", "tsv")
    stop("file.name is expected to have one of the following suffixes: .xlsx")
  }
  if (file.suffix == "ods") {
    stop("ods not handled, yet.")
  }
  if (!dir.exists(FCS.file.folder)) {
    stop(FCS.file.folder, " not found.")
  }

  wd <- dirname(FCS.file.folder)

  fcs.files <- get_fcs_identities(folder_path = FCS.file.folder,
                                  exclude_folders = exclude.folders,
                                  allow_duplicates = F)

  ### initiate
  if (!file.exists(file.path(wd, file.name))) {
    desc_file <- init_desc(wd, file.name = file.name, fcs.files = fcs.files)
    message(desc_file, "initiated.")
    invisible(desc_file)
  }

  desc <- read_desc(wd = wd, file.name = file.name, fcs.files = fcs.files)
  write_desc_log(wd = wd, file.name = file.name, desc = desc, write.log = write.log)

  ### find files for deletion
  change1 <- move_fcs_update_desc_on_disk(desc = desc,
                                          FCS.file.folder = FCS.file.folder,
                                          fcs.files = fcs.files,
                                          wd = wd,
                                          file.name = file.name,
                                          write.log = write.log)
  if (change1) {
    fcs.files <- get_fcs_identities(folder_path = FCS.file.folder,
                                    exclude_folders = exclude.folders,
                                    allow_duplicates = F)
    desc <- read_desc(wd = wd, file.name = file.name, fcs.files = fcs.files)
  }

  ### find new files for addition to desc
  change2 <- add_fcs_update_desc_on_disk(desc = desc,
                                         fcs.files = fcs.files,
                                         wd = wd,
                                         file.name = file.name,
                                         write.log = write.log)
  if (change2) {
    fcs.files <- get_fcs_identities(folder_path = FCS.file.folder,
                                    exclude_folders = exclude.folders,
                                    allow_duplicates = F)
    desc <- read_desc(wd = wd, file.name = file.name, fcs.files = fcs.files)
  }

  ### find files for renaming
  change3 <- rename_fcs_update_desc_on_disk(desc = desc,
                                            fcs.files = fcs.files,
                                            wd = wd,
                                            file.name = file.name,
                                            write.log = write.log)
  if (change3) {
    fcs.files <- get_fcs_identities(folder_path = FCS.file.folder,
                                    exclude_folders = exclude.folders,
                                    allow_duplicates = F)
    desc <- read_desc(wd = wd, file.name = file.name, fcs.files = fcs.files)
  }
  ### new order --> rename
  change4 <- reorder_fcs_update_desc_on_disk(desc = desc,
                                             fcs.files = fcs.files,
                                             wd = wd,
                                             file.name = file.name,
                                             write.log = write.log)

  print(tibble::tibble(operation = c("delete", "add", "rename", "reorder"),
                       done = c(change1, change2, change3, change4)))
}

write_desc_log <- function(wd, file.name, desc, write.log) {
  if (write.log) {
    if (Sys.info()[["sysname"]] %in% c("Linux", "Darwin")) {
      file <- file.path(wd, paste0(".log_", file.name))
    }
    if (Sys.info()[["sysname"]] == "Windows") {
      file <- file.path(wd, paste0("log_", file.name))
    }

    if (file.exists(file)) {
      log <- openxlsx::loadWorkbook(file)
    } else {
      log <- openxlsx::createWorkbook()
    }

    time <- format(Sys.time(), "%Y.%m.%d-%H.%M.%S")
    tryCatch({
      openxlsx::addWorksheet(log, time)
      openxlsx::writeData(log, time, desc)
      openxlsx::saveWorkbook(log, file, overwrite = T)
    }, error = function(e) {
      time <- paste0(time, sample(1:100, 1))
      openxlsx::addWorksheet(log, time)
      openxlsx::writeData(log, time, desc)
      openxlsx::saveWorkbook(log, file, overwrite = T)
    })

    '        if (Sys.info()[["sysname"]] == "Windows") {
            command <- paste0("attrib +h ", file)
            system(command)
        }'
  }
}



write_desc <- function(named.sheet.list, wd, file.name) {

  ext <- tolower(tools::file_ext(file.name))
  ## make repetitive elements more compact
  if (ext == "xlsx") {

    tryCatch({
      openxlsx::write.xlsx(named.sheet.list, file = file.path(wd, file.name), firstRow = T, colWidths = "auto", overwrite = T)
    }, error = function(e) {
      choice <- utils::menu(c("Yes", "No"), title = paste0("Error when writing sampledescription. Is the file still opened? If so, close it and give saving an updated version another try (1) or not (2)?"))
      if (choice == 1) {
        openxlsx::write.xlsx(named.sheet.list, file = file.path(wd, file.name), firstRow = T, colWidths = "auto", overwrite = T)
      }
      if (choice == 2) {
        message("Exiting.")
        return(NULL)
      }
    }, warning=function(w) {
      choice <- utils::menu(c("Yes", "No"), title = paste0("Error when writing sampledescription. Is the file still opened? If so, close it and give saving an updated version another try (1) or not (2)?"))
      if (choice == 1) {
        openxlsx::write.xlsx(named.sheet.list, file = file.path(wd, file.name), firstRow = T, colWidths = "auto", overwrite = T)
      }
      if (choice == 2) {
        message("Exiting.")
        return(NULL)
      }
    })

    ## read the putative updated file and check if FileNames are updated
    #browser()
    if (!identical(named.sheet.list[[1]][,"FileName",drop=T],openxlsx::read.xlsx(xlsxFile = file.path(wd, file.name))[,"FileName",drop=T])) {
      choice <- utils::menu(c("Yes", "No"), title = paste0("FileNames in sampledescription seem to not have changed. Is the file still opened? If so, close it and give saving an updated version another try (1) or not (2)?"))
      if (choice == 1) {
        openxlsx::write.xlsx(named.sheet.list, file = file.path(wd, file.name), firstRow = T, colWidths = "auto", overwrite = T)
      }
      if (choice == 2) {
        message("Exiting.")
        return(NULL)
      }
    }
  }
  if (ext == "tsv") {

    tryCatch({
      utils::write.table(x = named.sheet.list[[1]], file = file.path(wd, file.name), sep = "\t", row.names = F, na = "")
    }, error = function(e) {
      choice <- utils::menu(c("Yes", "No"), title = paste0("Error when writing sampledescription. Is the file still opened? If so, close it and give saving an updated version another try (1) or not (2)?"))
      if (choice == 1) {
        utils::write.table(x = named.sheet.list[[1]], file = file.path(wd, file.name), sep = "\t", row.names = F, na = "")
      }
      if (choice == 2) {
        message("Exiting.")
        return(NULL)
      }
    }, warning=function(w) {
      choice <- utils::menu(c("Yes", "No"), title = paste0("Error when writing sampledescription. Is the file still opened? If so, close it and give saving an updated version another try (1) or not (2)?"))
      if (choice == 1) {
        utils::write.table(x = named.sheet.list[[1]], file = file.path(wd, file.name), sep = "\t", row.names = F, na = "")
      }
      if (choice == 2) {
        message("Exiting.")
        return(NULL)
      }
    })

    ## read the putative updated file and check if FileNames are updated
    if (!identical(named.sheet.list[[1]][,"FileName",drop=T],utils::read.table(file = file.path(wd, file.name), header = T, sep = "\t", check.names = F)[,"FileName",drop=T])) {
      choice <- utils::menu(c("Yes", "No"), title = paste0("FileNames in sampledescription seem to not have changed. Is the file still opened? If so, close it and give saving an updated version another try (1) or not (2)?"))
      if (choice == 1) {
        utils::write.table(x = named.sheet.list[[1]], file = file.path(wd, file.name), sep = "\t", row.names = F, na = "")
      }
      if (choice == 2) {
        message("Exiting.")
        return(NULL)
      }
    }
  }
  if (tools::file_ext(file.name) %in% c("ods")) {
    #to do
  }

}



# check.FCS.files <- function(FCS.file.folder,
#                             exclude.folders = NULL,
#                             recursive = T) {
#
#   fcs.file.paths <- list.files(
#     path = FCS.file.folder,
#     pattern = "\\.fcs", #w/o $ is fine
#     full.names = T,
#     recursive = recursive,
#     ignore.case = T
#   )
#   if (length(fcs.file.paths) == 0) {
#     stop("No FCS files found.")
#   }
#
#   if (!is.null(exclude.folders)) {
#     exclude.pattern <- paste0(tolower(exclude.folders), collapse = "|")
#     fcs.file.paths <- fcs.file.paths[!grepl(exclude.pattern, tolower(fcs.file.paths))]
#   }
#
#   if (length(fcs.file.paths) == 0) {
#     stop("No FCS files left after filtering for exclusion folders.")
#   }
#
#   idents <- get_fcs_identities(kwl = flowCore::read.FCSheader(fcs.file.paths, emptyValue = F))
#
#   check_fcs_ident_uniqueness(idents)
#
#   return(idents)
# }

reorder_fcs_update_desc_on_disk <- function(desc,
                                            fcs.files,
                                            wd,
                                            file.name,
                                            write.log) {

  if (!identical(sort(desc[["FileName"]]), desc[["FileName"]])) {
    if (all(grepl("^[[:digit:]]{1,}_-_", desc[["FileName"]]))) {
      fcs.files <- stats::setNames(names(fcs.files), fcs.files)
      fcs.files <- fcs.files[desc[, "identity"]]
      desc[["FileName"]] <- gsub("^[[:digit:]]{1,}", "", desc[["FileName"]])
      desc[["FileName"]] <- paste0(sprintf("%04d", 1:nrow(desc)), desc[["FileName"]])
      desc[["FileName"]] <- normalize_filename(desc[["FileName"]])

      rows <- which(desc[["FileName"]] != basename(fcs.files))
      print(data.frame(FileName = desc[rows, "FileName"], PreviousFileName = basename(fcs.files)[rows]))

      choice <- 1
      if (interactive()) {
        choice <- utils::menu(c("Yes", "No"), title = "Rename FCS files as indicated?")
      }

      if (choice == 1) {
        file.rename(fcs.files, file.path(dirname(fcs.files), desc[["FileName"]]))
        write_desc(named.sheet.list = stats::setNames(list(desc), c("samples")), wd = wd, file.name = file.name)
        write_desc_log(wd = wd, file.name = file.name, desc = desc, write.log = write.log)
      } else {
        message("No files renamed.")
      }
    } else {
      message("No reodering as prefix-numbers were not detected accurately.")
    }
    return(TRUE)
  }
  return(FALSE)
}

rename_fcs_update_desc_on_disk <- function(desc,
                                           fcs.files,
                                           wd,
                                           file.name,
                                           write.log) {

  idx <- which(!desc[["FileName"]] %in% basename(names(fcs.files)))
  if (length(idx) > 0) {
    fcs.files <- stats::setNames(names(fcs.files), fcs.files)
    fcs.files <- fcs.files[desc[["identity"]]]
    desc[["FileName"]] <- ifelse(!grepl("^[[:digit:]]{1,}_-_", desc[["FileName"]]), paste0(sprintf("%04d", 1:nrow(desc)), "_-_", desc[["FileName"]]), desc[["FileName"]])
    desc[["FileName"]] <- normalize_filename(desc[["FileName"]])

    print(data.frame(FileName = desc[idx, "FileName"], PreviousFileName = basename(fcs.files[idx]), stringsAsFactors = F))

    choice <- 1
    if (interactive()) {
      choice <- utils::menu(c("Yes", "No"), title = "Rename FCS files as indicated?")
    }

    if (choice == 1) {
      file.rename(fcs.files[idx], file.path(dirname(fcs.files[idx]), desc[idx, "FileName"]))
      write_desc(named.sheet.list = stats::setNames(list(desc), c("samples")), wd = wd, file.name = file.name)
      write_desc_log(wd = wd, file.name = file.name, desc = desc, write.log = write.log)
      return(TRUE)
    }
  }
  return(FALSE)
}

add_fcs_update_desc_on_disk <- function(desc,
                                        fcs.files,
                                        wd,
                                        file.name,
                                        write.log) {

  new_fcs <- fcs.files[which(!fcs.files %in% desc[["identity"]])]
  if (length(new_fcs)) {
    new_fcs <- new_fcs[order(lubridate::parse_date_time(sapply(strsplit(new_fcs, "_-_"), "[", 3), orders = "%Y.%m.%d-%H.%M.%S"))]
    desc.diff <- data.frame(FileName = paste0(sprintf(paste0("%04d"), (nrow(desc) + 1):(nrow(desc) + length(new_fcs))), "_-_", basename(names(new_fcs))),
                            identity = new_fcs, stringsAsFactors = F)
    desc.diff[, c(names(desc)[which(!names(desc) %in% names(desc.diff))])] <- ""
    desc <- rbind(desc, desc.diff)

    file.rename(names(new_fcs), file.path(dirname(names(new_fcs)), desc.diff[["FileName"]]))
    write_desc(named.sheet.list = stats::setNames(list(desc), c("samples")), wd = wd, file.name = file.name)
    write_desc_log(wd = wd, file.name = file.name, desc = desc, write.log = write.log)
    message(nrow(desc.diff), " new files have been found and added to the sampledescription.")
    return(TRUE)
  }
  return(FALSE)
}

move_fcs_update_desc_on_disk <- function(desc,
                                         FCS.file.folder,
                                         fcs.files,
                                         wd,
                                         file.name,
                                         write.log,
                                         del_folder = "deleted_FCS_files") {

  idx <- intersect(which(is.na(desc[["FileName"]])), which(!is.na(desc[["identity"]])))
  if (length(idx) > 0) {
    del_files <- fcs.files[which(fcs.files %in% desc[idx, "identity"])]

    choice <- 1
    # if (interactive()) {
    #   choice <- utils::menu(c("Yes", "No"), title = paste0("Move these FCS files to deleted_FCS_files and exclude them from sampledescription: ", paste(names(del_files), collapse = ", ")))
    # }

    if (choice == 1) {
      del_folder <- file.path(FCS.file.folder, del_folder)
      dir.create(del_folder, showWarnings = F, recursive = T)

      if (dir.exists(del_folder)) {
        mv_files_safe(from = names(del_files),
                      to = file.path(del_folder, basename(names(del_files))))

        desc <- dplyr::filter(desc, !is.na(FileName))
        desc[["FileName"]] <- ifelse(grepl("^[[:digit:]]{1,}_-_", desc[["FileName"]]),
                                     paste0(sprintf("%04d", 1:nrow(desc)), "_-_", substr(desc[["FileName"]], 8, nchar(desc[["FileName"]]))),
                                     paste0(sprintf("%04d", 1:nrow(desc)), "_-_", desc[["FileName"]]))

        write_desc(
          named.sheet.list = stats::setNames(list(desc), c("samples")),
          wd = wd,
          file.name = file.name
        )
        write_desc_log(
          wd = wd,
          file.name = file.name,
          desc = desc,
          write.log = write.log
        )
        message("FCS files moved to ", del_folder, " and ", file.name, " updated.")
      } else {
        message(del_folder, " folder could not be created - no files were removed.")
      }
    }
    return(TRUE)
  }
  return(FALSE)
}


read_desc <- function(wd, file.name, fcs.files) {

  if (tools::file_ext(file.name) == "xlsx") {
    desc <- as.data.frame(openxlsx::read.xlsx(file.path(wd, file.name), sheet = 1, skipEmptyCols = F, detectDates = T), stringsAsFactors = F)
  }
  if (tools::file_ext(file.name) == "tsv") {
    desc <- as.data.frame(utils::read.table(file = file.path(wd, file.name), header = T, sep = "\t", check.names = F))
    desc[is.na(desc)] <- ""
    if (ncol(desc) == 1) {
      stop("sampledescription has only one column. Did you provide the wrong seperator (file.sep)?")
    }
  }
  if (tools::file_ext(file.name) == "ods") {
    #to do
  }

  desc[apply(desc,c(1,2),function(x) grepl("^ {1,}$", x))] <- NA # replace only-whitespace containing cells with NA
  desc <- desc[which(rowSums(is.na(desc)) < ncol(desc)), ] # rm pure NA rows
  desc <- desc[which(!is.na(desc$identity)),]

  if (any(!c("FileName", "identity") %in% names(desc))) {
    stop("Columns FileName and identity have to exist is the sampledescription file.")
  }

  if (nrow(desc) > length(fcs.files)) {
    print(desc[which(!desc[, "identity"] %in% fcs.files), which(names(desc) %in% c("FileName", "identity"))])
    stop("More rows in sampledescription than files in FCS.files.folder. For entries above no matching FCS files were found. Did you delete them manually? Please fix by deleting those rows manually in the xlsx-file. Then save it, close it and run sync_sampledescription again.")
  }

  validate_filenames(desc[["FileName"]])
  return(desc)
}

validate_filenames <- function(x) {
  invalid_pattern <- "[^A-Za-z0-9._-]"

  bad <- x[grepl(invalid_pattern, x)]

  if (length(bad)) {
    stop(
      "Invalid FileName(s):\n",
      paste(unique(bad), collapse = "\n"),
      "\n\nIllegal characters: anything but A-Za-z0-9._-"
    )
  }
}

mv_files_safe <- function(from, to) {

  if (any(file.exists(to))) {
    stop("Target file already exists:\n",
         paste(to[file.exists(to)], collapse = "\n"))
  }

  ok <- file.rename(from, to)

  if (!all(ok)) {
    stop("Failed to move some files:\n",
         paste(from[!ok], collapse = "\n"))
  }
}


init_desc <- function(wd, file.name, fcs.files) {

  other_putative_sd <- stats::na.omit(purrr::map_chr(list.files(wd, "\\.xlsx$|\\.tsv$", full.names = T), function(x){
    ext <- tolower(tools::file_ext(x))
    colnames <- switch(ext,
                       tsv  = names(utils::read.table(x, header = T, nrows = 2, sep = "\t")),
                       xlsx = names(openxlsx::read.xlsx(x, rows = c(1,2))))
    if (all(c("FileName", "identity") %in% colnames)) {
      return(x)
    }
    return(NA)
  }))

  if (length(other_putative_sd) > 0) {
    choice <- utils::menu(c("Yes", "No"), title = paste0("Other putative sampledescriptions found in the parent folder of FCS.file.folder: ", basename(other_putative_sd), ". Do you want to continue initiating another file (1)? If not (2), change the file.name argument and rerun sync_sampledescription."))
    if (choice == 2) {
      return(NULL)
    }
  }

  fcs.files <- fcs.files[order(lubridate::parse_date_time(
    sapply(strsplit(fcs.files, "_-_"), "[", 3),
    orders = "%Y.%m.%d-%H.%M.%S"
  ))]

  desc <- data.frame(FileName = paste0(
    sprintf(paste0("%04d"), seq_along(fcs.files)),
    "_-_",
    basename(names(fcs.files))
  ),
  identity = fcs.files,
  stringsAsFactors = F
  )
  init.columns <- c("AbCalcFile", "AbCalcSheet", "ExpProtocolFile", "ExpPart")
  desc[, init.columns] <- ""
  rownames(desc) <- NULL

  file.rename(names(fcs.files), file.path(dirname(names(fcs.files)), desc[["FileName"]]))
  write_desc(stats::setNames(list(desc), nm = c("samples")), wd = wd, file.name = file.name)
  write_desc_log(wd = wd, file.name = file.name, desc = desc, write.log = write.log)
  return(file.path(wd, file.name))
}

normalize_filename <- function(x) {
  x <- ifelse(!grepl("\\.fcs$", tolower(x)), paste0(x, ".fcs"), x)
  x <- sub("\\.FCS$", ".fcs", x)
  return(x)
}
