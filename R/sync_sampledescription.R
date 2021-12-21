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
#' @param init.columns additional columns to add to the initial file
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
                                   exclude.folders = c("compensation", "other_fcs_files", "experiment.file", "deleted_fcs_files"),
                                   init.columns = c("AbCalcFile", "AbCalcSheet", "ExpProtocolFile", "ExpPart"),
                                   write.log = T) {


  file.suffix <- rev(strsplit(file.name, "\\.")[[1]])[1]
  if (!file.suffix %in% c("xlsx", "ods", "txt", "tsv", "csv")) {
    stop("file.name is expected to have one of the following suffixes: .xlsx, .ods, .txt, .tsv, .csv.")
  }
  if (file.suffix == "ods") {
    stop("ods not handled, yet.")
  }

  if (file.suffix %in% c("txt", "tsv")) {
    file.sep <- "\t"
  }
  if (file.suffix == "csv") {
    file.sep <- ";"
  }

  if (!dir.exists(FCS.file.folder)) {
    stop(paste0(FCS.file.folder, " not found."))
  }
  wd <- dirname(FCS.file.folder)
  fcs.files <- .check.FCS.files(FCS.file.folder = FCS.file.folder, exclude.folders = exclude.folders)

  # initiate
  if (!file.exists(file.path(wd, file.name))) {

    # check existing files (put in separate function)
    for (x in list.files(wd, "\\.xlsx$", full.names = T)) {
      if(all(c("FileName", "identity") %in% names(openxlsx::read.xlsx(x, rows = c(1,2))))) {
        choice <- utils::menu(c("Yes", "No"), title = paste0("Another putative sampledescription file was found in the parent folder of FCS.file.folder: ", basename(x), ". Do you want to continue initiating another file (type 1)? If not change the file.name argument to ", basename(x), " and type 2."))
        if (choice == 2) {
          return(NULL)
        }
      }
    }
    for (x in list.files(wd, "\\.txt$|\\.tsv$|\\.csv$", full.names = T)) {
      if(all(any(grepl("FileName", names(utils::read.table(x, header = T, nrows = 2, sep = "\t")))), any(grepl("identity", names(utils::read.table(x, header = T, nrows = 2, sep = "\t")))))) {
        choice <- utils::menu(c("Yes", "No"), title = paste0("Another putative sampledescription file was found in the parent folder of FCS.file.folder: ", basename(x), ". Do you want to continue initiating another file (type 1)? If not change the file.name argument to ", basename(x), " and type 2."))
        if (choice == 2) {
          return(NULL)
        }
      }
    }
    for (x in list.files(wd, "\\.ods$", full.names = T)) {
      # to do
      # readODS::read_ods()
      #readODS::write_ods(sd, file.path(wd, "sampledescription.ods"), sheet = "samples")
    }

    fcs.files <- fcs.files[order(lubridate::parse_date_time(sapply(strsplit(fcs.files, "_-_"), "[", 3), orders = "%Y.%m.%d-%H.%M.%S"))]
    sd <- data.frame(FileName = paste0(sprintf(paste0("%04d"), seq_along(fcs.files)), "_-_", basename(names(fcs.files))), identity = fcs.files, stringsAsFactors = F)
    sd[, init.columns] <- ""
    rownames(sd) <- NULL

    .write.sd(stats::setNames(list(sd), nm = c("samples")), wd = wd, file.name = file.name, file.sep = file.sep)
    .write.sd.log(wd = wd, file.name = file.name, sd = sd, write.log = write.log, file.sep = file.sep)
    file.rename(names(fcs.files), file.path(dirname(names(fcs.files)), sd[, "FileName"]))
    return(paste0(file.name, " initiated."))
  }

  sd <- .read.and.check.sd(wd = wd, file.name = file.name, fcs.files = fcs.files, file.sep = file.sep)
  .write.sd.log(wd = wd, file.name = file.name, sd = sd, write.log = write.log, file.sep = file.sep)

  # find files for deletion
  sd.delete.ind <- intersect(which(is.na(sd[, "FileName"])), which(!is.na(sd[, "identity"])))
  if (length(sd.delete.ind) > 0) {
    fcs.files.del <- fcs.files[which(fcs.files %in% sd[sd.delete.ind, "identity"])]
    print(names(fcs.files.del))

    if (interactive()) {
      choice <- utils::menu(c("Yes", "No"), title = "Move these FCS files to deleted_FCS_files and exclude them from sampledescription?")
    } else {
      choice <- 1
    }

    if (choice == 1) {
      dir.create(file.path(FCS.file.folder, "deleted_FCS_files"), showWarnings = F, recursive = T)

      if (dir.exists(file.path(FCS.file.folder, "deleted_FCS_files"))) {
        file.copy(names(fcs.files.del), file.path(FCS.file.folder, "deleted_FCS_files", basename(names(fcs.files.del))))
        file.remove(names(fcs.files.del))
        sd <- sd[which(!is.na(sd[, "FileName"])), ]
        sd[, "FileName"] <- ifelse(grepl("^[[:digit:]]{1,}_-_", sd[, "FileName"]), paste0(sprintf("%04d", 1:nrow(sd)), "_-_", substr(sd[, "FileName"],
                                                                                                                                     8, nchar(sd[, "FileName"]))), paste0(sprintf("%04d", 1:nrow(sd)), "_-_", sd[, "FileName"]))
        .write.sd(named.sheet.list = stats::setNames(list(sd), c("samples")), wd = wd, file.name = file.name, file.sep = file.sep)
        .write.sd.log(wd = wd, file.name = file.name, sd = sd, write.log = write.log, file.sep = file.sep)
        print(paste0("FCS files moved to deleted_FCS_files and ", file.name, " updated."))
      } else {
        print("deleted_FCS_files folder could not be created - no files were removed.")
      }
    }
    if (choice == 2) {
      return("No files removed.")
    }
  }

  fcs.files <- .check.FCS.files(FCS.file.folder = FCS.file.folder, exclude.folders = exclude.folders)
  sd <- .read.and.check.sd(wd = wd, file.name = file.name, fcs.files = fcs.files, file.sep = file.sep)

  # find new files for addition to sd
  fcs.files.diff <- fcs.files[which(!fcs.files %in% sd[, "identity"])]
  if (length(fcs.files.diff) != 0) {
    fcs.files.diff <- fcs.files.diff[order(lubridate::parse_date_time(sapply(strsplit(fcs.files.diff, "_-_"), "[", 3), orders = "%Y.%m.%d-%H.%M.%S"))]
    sd.diff <- data.frame(FileName = paste0(sprintf(paste0("%04d"), (nrow(sd) + 1):(nrow(sd) + length(fcs.files.diff))), "_-_", basename(names(fcs.files.diff))),
                          identity = fcs.files.diff, stringsAsFactors = F)
    sd.diff[, c(names(sd)[which(!names(sd) %in% names(sd.diff))])] <- ""
    sd <- rbind(sd, sd.diff)

    .write.sd(named.sheet.list = stats::setNames(list(sd), c("samples")), wd = wd, file.name = file.name, file.sep = file.sep)
    .write.sd.log(wd = wd, file.name = file.name, sd = sd, write.log = write.log, file.sep = file.sep)
    print(paste0(nrow(sd.diff), " new files have been found and added to the sampledescription."))
    file.rename(names(fcs.files.diff), file.path(dirname(names(fcs.files.diff)), sd.diff[, "FileName"]))
  }

  fcs.files <- .check.FCS.files(FCS.file.folder = FCS.file.folder, exclude.folders = exclude.folders)
  sd <- .read.and.check.sd(wd = wd, file.name = file.name, fcs.files = fcs.files, file.sep = file.sep)

  # find files for renaming
  sd.rename.ind <- which(!sd[, "FileName"] %in% basename(names(fcs.files)))
  if (length(sd.rename.ind) > 0) {
    fcs.files <- stats::setNames(names(fcs.files), fcs.files)
    fcs.files <- fcs.files[sd[, "identity"]]
    sd[, "FileName"] <- ifelse(!grepl("^[[:digit:]]{1,}_-_", sd[, "FileName"]), paste0(sprintf("%04d", 1:nrow(sd)), "_-_", sd[, "FileName"]), sd[, "FileName"])
    sd[, "FileName"] <- ifelse(!grepl("\\.fcs$", tolower(sd[, "FileName"])), paste0(sd[, "FileName"], ".fcs"), sd[, "FileName"])
    sd[, "FileName"] <- sub("\\.FCS$", ".fcs", sd[, "FileName"])

    print(data.frame(FileName = sd[sd.rename.ind, "FileName"], PreviousFileName = basename(fcs.files[sd.rename.ind]), stringsAsFactors = F))
    if (interactive()) {
      choice <- utils::menu(c("Yes", "No"), title = "Rename FCS files as indicated?")
    } else {
      choice <- 1
    }

    if (choice == 1) {
      .write.sd(named.sheet.list = stats::setNames(list(sd), c("samples")), wd = wd, file.name = file.name, file.sep = file.sep)
      .write.sd.log(wd = wd, file.name = file.name, sd = sd, write.log = write.log, file.sep = file.sep)
      file.rename(fcs.files[sd.rename.ind], file.path(dirname(fcs.files[sd.rename.ind]), sd[sd.rename.ind, "FileName"]))
    }
    if (choice == 2) {
      return("No files renamed.")
    }
  }

  fcs.files <- .check.FCS.files(FCS.file.folder = FCS.file.folder, exclude.folders = exclude.folders)
  sd <- .read.and.check.sd(wd = wd, file.name = file.name, fcs.files = fcs.files, file.sep = file.sep)

  # new order
  if (!identical(sort(sd[, "FileName"]), sd[, "FileName"])) {
    if (all(grepl("^[[:digit:]]{1,}_-_", sd[, "FileName"]))) {
      fcs.files <- stats::setNames(names(fcs.files), fcs.files)
      fcs.files <- fcs.files[sd[, "identity"]]
      sd[, "FileName"] <- gsub("^[[:digit:]]{1,}", "", sd[, "FileName"])
      sd[, "FileName"] <- paste0(sprintf("%04d", 1:nrow(sd)), sd[, "FileName"])
      sd[, "FileName"] <- ifelse(!grepl("\\.fcs$", tolower(sd[, "FileName"])), paste0(sd[, "FileName"], ".fcs"), sd[, "FileName"])
      sd[, "FileName"] <- sub("\\.FCS$", ".fcs", sd[, "FileName"])

      rows <- which(sd[, "FileName"] != basename(fcs.files))
      print(data.frame(FileName = sd[rows, "FileName"], PreviousFileName = basename(fcs.files)[rows]))
      if (interactive()) {
        choice <- utils::menu(c("Yes", "No"), title = "Rename FCS files as indicated?")
      } else {
        choice <- 1
      }

      if (choice == 1) {
        .write.sd(named.sheet.list = stats::setNames(list(sd), c("samples")), wd = wd, file.name = file.name, file.sep = file.sep)
        .write.sd.log(wd = wd, file.name = file.name, sd = sd, write.log = write.log, file.sep = file.sep)
        file.rename(fcs.files, file.path(dirname(fcs.files), sd[, "FileName"]))
      }
      if (choice == 2) {
        return("No files renamed.")
      }
    } else {
      stop("No reodering as prefix-numbers were not detected accurately.")
    }
  }

}

.write.sd.log <- function(wd, file.name, sd, write.log, file.sep) {
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
      openxlsx::writeData(log, time, sd)
      openxlsx::saveWorkbook(log, file, overwrite = T)
    }, error = function(e) {
      time <- paste0(time, sample(1:100, 1))
      openxlsx::addWorksheet(log, time)
      openxlsx::writeData(log, time, sd)
      openxlsx::saveWorkbook(log, file, overwrite = T)
    })

    '        if (Sys.info()[["sysname"]] == "Windows") {
            command <- paste0("attrib +h ", file)
            system(command)
        }'
  }
}

.write.sd <- function(named.sheet.list, wd, file.name, file.sep) {
  if (rev(strsplit(file.name, "\\.")[[1]])[1] == "xlsx") {
    tryCatch({
      openxlsx::write.xlsx(named.sheet.list, file = file.path(wd, file.name), firstRow = T, colWidths = "auto", overwrite = T)
    }, error = function(e) {
      new <- file.path(wd, paste0(format(Sys.time(), "%Y.%m.%d-%H.%M.%S_"), file.name))
      print(paste0("Is ", file.name, " still open in Excel? Saving as updated file as ", new, ". Please delete the former one manually and remove the date-prefix of the new file."))
      openxlsx::write.xlsx(named.sheet.list, file = new, firstRow = T, colWidths = "auto")
    })
  }
  if (rev(strsplit(file.name, "\\.")[[1]])[1] %in% c("txt", "tsv", "csv")) {
    tryCatch({
      utils::write.table(x = named.sheet.list[[1]], file = file.path(wd, file.name), sep = file.sep, row.names = F)
    }, error = function(e) {
      new <- file.path(wd, paste0(format(Sys.time(), "%Y.%m.%d-%H.%M.%S_"), file.name))
      print(paste0("Is ", file.name, " still opened? Saving as updated file as ", new, ". Please delete the former one manually and remove the date-prefix of the new file."))
      utils::write.table(x = named.sheet.list[[1]], file = new, sep = file.sep)
    })
  }
  if (rev(strsplit(file.name, "\\.")[[1]])[1] %in% c("ods")) {
    #to do
  }

}

.check.FCS.files <- function(FCS.file.folder, exclude.folders = NULL) {

  if (!requireNamespace("BiocManager", quietly = T)){
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("flowCore", quietly = T)){
    BiocManager::install("flowCore")
  }

  fcs.file.paths <- list.files(path = FCS.file.folder, pattern = "\\.fcs", full.names = T, recursive = T, ignore.case = T)
  if (!is.null(exclude.folders)) {
    fcs.file.paths <- fcs.file.paths[which(!grepl(paste0(tolower(exclude.folders), collapse = "|"), tolower(fcs.file.paths)))]
  }

  if (length(fcs.file.paths) == 0) {
    stop("No FCS files found or left after filtering for exclusion folders.")
  }

  fcs.files <- sapply(fcs.file.paths, function(x) {
    ff <- flowCore::read.FCS(x, which.lines = 1, emptyValue = F, truncate_max_range = F)
    dd <- flowCore::keyword(ff)[["$DATE"]]
    tt <- flowCore::keyword(ff)[["$BTIM"]]
    if (nchar(tt) - nchar(gsub(":", "", tt)) > 2) {
      tt <- paste(rev(rev(strsplit(tt, ":")[[1]])[-1]), collapse = ":")
    }
    datetime <- paste0(dd, "-", tt)
    # if analysis starts at 23:5x and ends at 00:xx then date of the next day is assigned - this is problematic though for ordering of samples and has
    # to be corrected; subtract the number of seconds of one day (86400) to get the correct date for ordering samples
    if (grepl("^2[[:digit:]]", flowCore::keyword(ff)[["$BTIM"]]) & grepl("^0[[:digit:]]", flowCore::keyword(ff)[["$ETIM"]])) {
      sub <- 86400
    } else {
      sub <- 0
    }
    datetime <- format(lubridate::parse_date_time(datetime, orders = c("%Y-%b-%d-%H:%M:%S", "%Y-%B-%d-%H:%M:%S", "%Y-%m-%d-%H:%M:%S", "%d-%b-%Y-%H:%M:%S",
                                                                       "%d-%m-%Y-%H:%M:%S", "%d-%B-%Y-%H:%M:%S", "%d-%b-%Y-%H:%M:%S")) - sub, "%Y.%m.%d-%H.%M.%S")
    if (is.na(datetime)) {
      print(paste0("datetime ", paste0(dd, "-", tt), " could not be converted to a uniform format. Please, provide this to the package-maintainer."))
    }
    identity <- paste0(flowCore::keyword(ff)[["$FIL"]], "_-_", trimws(flowCore::keyword(ff)[["$TOT"]]), "_-_", datetime)
    return(identity)
  })

  if (length(unique(fcs.files)) != length(fcs.files)) {
    stop(paste0("Duplicate FCS files found. This is not allowed. Please, remove one of each duplicates. \n", paste(names(fcs.files[duplicated(fcs.files) |
                                                                                                                                     duplicated(fcs.files, fromLast = T)]), collapse = "\n")))
  }
  return(fcs.files)
}

.read.and.check.sd <- function(wd, file.name, fcs.files, file.sep) {
  if (rev(strsplit(file.name, "\\.")[[1]])[1] == "xlsx") {
    sd <- as.data.frame(openxlsx::read.xlsx(file.path(wd, file.name), sheet = 1, skipEmptyCols = F, detectDates = T), stringsAsFactors = F)
  }
  if (rev(strsplit(file.name, "\\.")[[1]])[1] %in% c("txt", "tsv", "csv")) {
    sd <- as.data.frame(utils::read.table(file = file.path(wd, file.name), header = T, sep = file.sep))
    if (ncol(sd) == 1) {
      stop("sampledescription has only one column. Did you provide the wrong seperator (file.sep)?")
    }
  }
  if (rev(strsplit(file.name, "\\.")[[1]])[1] %in% c("ods")) {
    #to do
  }
  sd <- sd[which(rowSums(is.na(sd)) < ncol(sd)), ]
  if (any(!c("FileName", "identity") %in% names(sd))) {
    stop("Columns FileName and identity have to exist is the sampledescription file.")
  }
  if (nrow(sd) > length(fcs.files)) {
    print(sd[which(!sd[, "identity"] %in% fcs.files), which(names(sd) %in% c("FileName", "identity"))])
    stop("More rows in sampledescription than files in FCS.files.folder. For entries above no matching FCS files were found. Did you delete them manually? Please fix by deleting those rows manually in the xlsx-file. Then save it, close it and run sync_sampledescription again.")
  }
  if (any(sapply(c("/", ":", "\\|", "\\?", "\\!", "\\*", "<", ">", "'", "\""), function(x) grepl(x, sd[, "FileName"])))) {
    stop("There is at least one FileName with one or more illegal character(s) which may cause problems in file-naming ( / : | ? ! * < > ' \")")
  }
  return(sd)
}
