#' Have an excel file (sampledescription) created or modified to store meta data of FCS files in
#'
#' Four cases can be handled: (i) if no xlsx-file named xlsx.file.name exists
#' in wd a xlsx-file is initiated based on the FCS files in FCS.file.folder. (ii)
#' When new FCS files are added to FCS.file.folder and the function is called
#' they are added in order of acquisition to the xlsx-file. (iii) When file names
#' in the FileName column of the xlsx-file are changed and the function is called
#' the FCS files in FCS.file.folder are renamed accordingly. (iv) If FCS files
#' are to be excluded or removed the entry in the FileName column has to be left
#' blank and the function has to be called.
#'
#' Attention: Close the xlsx-file before calling the function so the file can be edited.
#' Never edit the 'identity' columns in the xlsx-file manually.
#'
#'
#' @param FCS.file.folder path to the root folder which contains FCS files
#' @param xlsx.file.name name of the sampledescription file
#' @param exclude.folders character vector of folders to exclude when checking for FCS files
#'
#' @return No return value. Instead sampledescription.xlsx and FCS files are synced.
#' @export
#'
#' @examples
#' \dontrun{
#' # When the script is saved to R_scripts in the experiment folder, this will
#' write the path of the experiment folder into wd:
#' wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
#' # If the sampledescription is to be initiated call
#' sync_sampledescription(FCS.file.folder = file.path(wd, "FCS_files"))
#' }
sync_sampledescription <- function(FCS.file.folder, xlsx.file.name = "sampledescription.xlsx", exclude.folders = c("compensation",
    "other_fcs_files", "experiment.file", "deleted_fcs_files")) {

    if (!dir.exists(FCS.file.folder)) {
        stop(paste0(FCS.file.folder, " not found."))
    }
    wd <- dirname(FCS.file.folder)

    fcs.files <- .check.FCS.files(FCS.file.folder = FCS.file.folder, exclude.folders = exclude.folders)

    # initiate
    if (!file.exists(file.path(wd, xlsx.file.name))) {
        fcs.files <- fcs.files[order(lubridate::parse_date_time(sapply(strsplit(fcs.files, "_-_"), "[", 3), orders = "%Y.%m.%d-%H.%M.%S", locale = "en_GB.UTF-8"))]
        sd <- data.frame(FileName = paste0(sprintf(paste0("%04d"), seq_along(fcs.files)), "_-_", basename(names(fcs.files))), identity = fcs.files, stringsAsFactors = FALSE)
        sd[, c("AbCalcFile", "AbCalcSheet", "ExpProtocolFile", "ExpPart")] <- ""

        .write.sd(stats::setNames(list(sd), nm = c("samples")), wd = wd, xlsx.file.name = xlsx.file.name)

        if (Sys.info()[["sysname"]] == "Darwin") {
            .write.sd.log(wd = wd, xlsx.file.name = xlsx.file.name, sd = sd)
        }
        file.rename(names(fcs.files), file.path(dirname(names(fcs.files)), sd[, "FileName"]))
        return(paste0(xlsx.file.name, " initiated."))
    }

    sd <- .read.and.check.sd(wd = wd, xlsx.file.name = xlsx.file.name, fcs.files = fcs.files)

    if (Sys.info()[["sysname"]] == "Darwin") {
        .write.sd.log(wd = wd, xlsx.file.name = xlsx.file.name, sd = sd)
    }

    # find files for deletion
    sd.delete.ind <- intersect(which(is.na(sd[,"FileName"])), which(!is.na(sd[,"identity"])))
    if (length(sd.delete.ind) > 0) {
        fcs.files.del <- fcs.files[which(fcs.files %in% sd[sd.delete.ind,"identity"])]
        print(names(fcs.files.del))
        choice <- utils::menu(c("Yes", "No"), title = "Move these FCS files to deleted_FCS_files and exclude them from sampledescription?")

        if (choice == 1) {
            dir.create(file.path(FCS.file.folder, "deleted_FCS_files"), showWarnings = F, recursive = T)

            if (dir.exists(file.path(FCS.file.folder, "deleted_FCS_files"))) {
                file.copy(names(fcs.files.del), file.path(FCS.file.folder, "deleted_FCS_files", basename(names(fcs.files.del))))
                file.remove(names(fcs.files.del))

                sd <- sd[which(!is.na(sd[,"FileName"])),]
                sd[,"FileName"] <- ifelse(grepl("^[[:digit:]]{1,}_-_", sd[,"FileName"]), paste0(sprintf("%04d", 1:nrow(sd)), "_-_", substr(sd[,"FileName"], 8, nchar(sd[,"FileName"]))), paste0(sprintf("%04d", 1:nrow(sd)), "_-_", sd[,"FileName"]))
                .write.sd(named.sheet.list = stats::setNames(list(sd), c("samples")), wd = wd, xlsx.file.name = xlsx.file.name)
                print(paste0("FCS files moved and ", xlsx.file.name, " updated."))
            } else {
                print("deleted_FCS_files folder could not be created - no files were removed.")
            }
        }
        if (choice == 2) {
            return("No files removed.")
        }
    }

    fcs.files <- .check.FCS.files(FCS.file.folder = FCS.file.folder, exclude.folders = exclude.folders)
    sd <- .read.and.check.sd(wd = wd, xlsx.file.name = xlsx.file.name, fcs.files = fcs.files)

    # find new files for addition to sd
    fcs.files.diff <- fcs.files[which(!fcs.files %in% sd[,"identity"])]
    if (length(fcs.files.diff) != 0) {

        fcs.files.diff <- fcs.files.diff[order(lubridate::parse_date_time(sapply(strsplit(fcs.files.diff, "_-_"), "[", 3), orders = "%Y.%m.%d-%H.%M.%S", locale = "en_GB.UTF-8"))]

        sd.diff <- data.frame(FileName = paste0(sprintf(paste0("%04d"), (nrow(sd) + 1):(nrow(sd) + length(fcs.files.diff))), "_-_", basename(names(fcs.files.diff))),
                              identity = fcs.files.diff,
                              stringsAsFactors = FALSE)
        sd.diff[,c(names(sd)[which(!names(sd) %in% names(sd.diff))])] <- ""
        sd <- rbind(sd, sd.diff)

        .write.sd(named.sheet.list = stats::setNames(list(sd), c("samples")), wd = wd, xlsx.file.name = xlsx.file.name)
        print(paste0(nrow(sd.diff), " new files have been found and added to the sampledescription."))
        file.rename(names(fcs.files.diff), file.path(dirname(names(fcs.files.diff)), sd.diff[,"FileName"]))
    }

    fcs.files <- .check.FCS.files(FCS.file.folder = FCS.file.folder, exclude.folders = exclude.folders)
    sd <- .read.and.check.sd(wd = wd, xlsx.file.name = xlsx.file.name, fcs.files = fcs.files)

    # find files for renaming
    sd.rename.ind <- which(!sd[,"FileName"] %in% basename(names(fcs.files)))
    if (length(sd.rename.ind) > 0) {
        fcs.files <- stats::setNames(names(fcs.files), fcs.files)
        fcs.files <- fcs.files[sd[,"identity"]]
        sd[,"FileName"] <- ifelse(!grepl("^[[:digit:]]{1,}_-_", sd[,"FileName"]), paste0(sprintf("%04d", 1:nrow(sd)), "_-_", sd[,"FileName"]), sd[,"FileName"])
        sd[,"FileName"] <- ifelse(!grepl("\\.fcs$", tolower(sd[,"FileName"])), paste0(sd[,"FileName"], ".fcs"), sd[,"FileName"])
        sd[,"FileName"] <- sub("\\.FCS$", ".fcs", sd[,"FileName"])

        print(data.frame(FileName = sd[sd.rename.ind, "FileName"], PreviousFileName = basename(fcs.files[sd.rename.ind])))
        choice <- utils::menu(c("Yes", "No"), title = "Rename FCS files as indicated?")

        if (choice == 1) {
            .write.sd(named.sheet.list = stats::setNames(list(sd), c("samples")), wd = wd, xlsx.file.name = xlsx.file.name)
            file.rename(fcs.files[sd.rename.ind], file.path(dirname(fcs.files[sd.rename.ind]), sd[sd.rename.ind,"FileName"]))
        }
        if (choice == 2) {
            return("No files renamed.")
        }
    }

    fcs.files <- .check.FCS.files(FCS.file.folder = FCS.file.folder, exclude.folders = exclude.folders)
    sd <- .read.and.check.sd(wd = wd, xlsx.file.name = xlsx.file.name, fcs.files = fcs.files)

    # new order
    if (!identical(sort(sd[,"FileName"]), sd[,"FileName"])) {
        if (all(grepl("^[[:digit:]]{1,}_-_", sd[,"FileName"]))) {
            fcs.files <- stats::setNames(names(fcs.files), fcs.files)
            fcs.files <- fcs.files[sd[,"identity"]]
            sd[,"FileName"] <- gsub("^[[:digit:]]{1,}", "", sd[,"FileName"])
            sd[,"FileName"] <- paste0(sprintf("%04d", 1:nrow(sd)), sd[,"FileName"])
            sd[,"FileName"] <- ifelse(!grepl("\\.fcs$", tolower(sd[,"FileName"])), paste0(sd[,"FileName"], ".fcs"), sd[,"FileName"])
            sd[,"FileName"] <- sub("\\.FCS$", ".fcs", sd[,"FileName"])

            print(data.frame(FileName = sd[, "FileName"], PreviousFileName = basename(fcs.files)))
            choice <- utils::menu(c("Yes", "No"), title = "Rename FCS files as indicated?")

            if (choice == 1) {
                .write.sd(named.sheet.list = stats::setNames(list(sd), c("samples")), wd = wd, xlsx.file.name = xlsx.file.name)
                file.rename(fcs.files, file.path(dirname(fcs.files), sd[,"FileName"]))
            }
            if (choice == 2) {
                return("No files renamed.")
            }
        } else {
            stop("No reodering as prefix-numbers not detected accurately.")
        }
    }

}

.write.sd.log <- function(wd, xlsx.file.name, sd) {
    # find out how to do on windows or omit on windows
    file <- file.path(wd, paste0(".log_", xlsx.file.name))
    if (file.exists(file)) {
        log <- openxlsx::loadWorkbook(file)
    } else {
        log <- openxlsx::createWorkbook()
    }
    time <- format(Sys.time(), "%Y.%m.%d-%H.%M.%S")
    openxlsx::addWorksheet(log, time)
    openxlsx::writeData(log, time, sd)
    openxlsx::saveWorkbook(log, file, overwrite = T)
}

.write.sd <- function(named.sheet.list, wd, xlsx.file.name) {
    tryCatch(openxlsx::write.xlsx(named.sheet.list, file = file.path(wd, xlsx.file.name), firstRow = T, colWidths = "auto", overwrite = T), error = function(e) {
        print(paste0("Is ", xlsx.file.name, " still open in Excel? Saving as updated file as ", file.path(wd, paste0(format(Sys.time(), "%Y.%m.%d-%H.%M.%S_"),
            xlsx.file.name)), ". Please delete the former one and remove the date-prefix of the new file."))
        openxlsx::write.xlsx(named.sheet.list, file = file.path(wd, paste0(format(Sys.time(), "%Y.%m.%d-%H.%M.%S_"), xlsx.file.name)), firstRow = T, colWidths = "auto")
    })
}

.check.FCS.files <- function(FCS.file.folder, exclude.folders) {
    fcs.file.paths <- list.files(path = FCS.file.folder, pattern = "\\.fcs", full.names = T, recursive = T, ignore.case = T)
    fcs.file.paths <- fcs.file.paths[which(!grepl(paste0(tolower(exclude.folders), collapse = "|"), tolower(fcs.file.paths)))]

    if (length(fcs.file.paths) == 0) {
        stop("No FCS files found or left over after filtering for exclusion folders.")
    }

    fcs.files <- sapply(fcs.file.paths, function(x) {
        ff <- flowCore::read.FCS(x, which.lines = 1, emptyValue = F, truncate_max_range = F)
        datetime <- paste0(flowCore::keyword(ff)[["$DATE"]], "-", flowCore::keyword(ff)[["$BTIM"]])
        # if analysis starts at 23:5x and ends at 00:xx then date of the next day is assigned - this is problematic though for ordering of samples and has
        # to be corrected; subtract the number of seconds of one day (86400) to get the correct date for ordering samples
        if (grepl("^2[[:digit:]]", flowCore::keyword(ff)[["$BTIM"]]) & grepl("^0[[:digit:]]", flowCore::keyword(ff)[["$ETIM"]])) {
            sub <- 86400
        } else {
            sub <- 0
        }
        datetime <- format(lubridate::parse_date_time(datetime, orders = c("%Y-%b-%d-%H:%M:%S", "%Y-%B-%d-%H:%M:%S", "%Y-%m-%d-%H:%M:%S", "%d-%b-%Y-%H:%M:%S",
            "%d-%m-%Y-%H:%M:%S", "%d-%B-%Y-%H:%M:%S"), locale = "en_GB.UTF-8") - sub, "%Y.%m.%d-%H.%M.%S")
        identity <- paste0(flowCore::keyword(ff)[["$FIL"]], "_-_", trimws(flowCore::keyword(ff)[["$TOT"]]), "_-_", datetime)
        return(identity)
    })

    if (length(unique(fcs.files)) != length(fcs.files)) {
        stop(paste0("Duplicate FCS files found. This is not allowed. Please, remove one of each duplicates. \n", paste(names(fcs.files[duplicated(fcs.files) | duplicated(fcs.files, fromLast = T)]),
            collapse = "\n")))
    }
    return(fcs.files)
}

.read.and.check.sd <- function(wd, xlsx.file.name, fcs.files) {
    sd <- as.data.frame(openxlsx::read.xlsx(file.path(wd, xlsx.file.name), sheet = 1, skipEmptyCols = F, detectDates = T))
    sd <- sd[which(rowSums(is.na(sd)) < ncol(sd)), ]
    if (any(!c("FileName", "identity") %in% names(sd))) {
        stop("Columns FileName and identity have to exist is the sampledescription file.")
    }
    if (nrow(sd) > length(fcs.files)) {
        print(sd[which(!sd[,"identity"] %in% fcs.files),which(names(sd) %in% c("FileName", "identity"))])
        stop("More rows in sampledescription than files in FCS.files.folder. For entries above no matching FCS files were found. Did you delete them manually? Please fix by deleting those rows manually in the xlsx-file. Then save it, close it and run sync_sampledescription again.")
    }
    if (any(sapply(c("/", ":", "\\|", "\\?", "\\!", "\\*", "<", ">", "'", "\""), function(x) grepl(x, sd[,"FileName"])))) {stop("There is at least one FileName with one or more illegal character(s) which may cause problems in file-naming ( / : | ? ! * < > ' \")")}
    return(sd)
}
