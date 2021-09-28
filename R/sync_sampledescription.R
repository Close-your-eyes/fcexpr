#' Title
#'
#' @param wd character
#' @param FCS.file.folder character
#' @param xlsx.file.name character
#' @param sample.sheet.name character
#' @param exclude.folders character
#'
#' @return No return value. Instead a sampledescription and FCS files are synced.
#' @export
#'
#' @examples
sync_sampledescription <- function(wd, FCS.file.folder, xlsx.file.name = "sampledescription.xlsx", sample.sheet.name = "samples", exclude.folders = c("compensation",
    "other_fcs_files", "experiment.file", "deleted_fcs_files")) {

    if (missing(wd)) {
        wd <- getwd()
    }
    if (missing(FCS.file.folder)) {
        FCS.file.folder <- file.path(wd, "FCS_files")
    }
    if (!dir.exists(FCS.file.folder)) {
        stop(paste0(FCS.file.folder, " not found."))
    }
    fcs.files <- check.FCS.files(FCS.file.folder = FCS.file.folder, exclude.folders = exclude.folders)

    # initiate
    if (!file.exists(file.path(wd, xlsx.file.name))) {
        # order by acquisition datetime
        fcs.files <- fcs.files[order(lubridate::parse_date_time(sapply(strsplit(fcs.files, "_-_"), "[", 3), orders = "%Y.%m.%d-%H.%M.%S", locale = "en_GB.UTF-8"))]
        sd.init <- data.frame(FileName = paste0(sprintf(paste0("%04d"), seq_along(fcs.files)), "_-_", basename(names(fcs.files))), FilePath = names(fcs.files),
            identity = fcs.files, stringsAsFactors = FALSE)
        sd.init[, c("AbCalcFile", "AbCalcSheet", "ExpProtocolFile", "ExpPart")] <- ""

        write.sd(stats::setNames(list(sd.init, sd.init[, c("FilePath", "identity")]), nm = c(sample.sheet.name, "initial.file.names")), wd = wd, xlsx.file.name = xlsx.file.name)

        if (Sys.info()[["sysname"]] == "Darwin") {
            write.sd.log(wd = wd, xlsx.file.name = xlsx.file.name, sd = sd.init)
        }
        file.rename(sd.init[, "FilePath"], file.path(dirname(sd.init[, "FilePath"]), sd.init[, "FileName"]))
        return(paste0(xlsx.file.name, " initiated."))
    }

    sd <- read.and.check.sd(wd = wd, xlsx.file.name = xlsx.file.name, sample.sheet.name = sample.sheet.name)

    if (Sys.info()[["sysname"]] == "Darwin") {
        write.sd.log(wd = wd, xlsx.file.name = xlsx.file.name, sd = sd)
    }

    # find files for deletion
    sd.delete.ind <- intersect(which(is.na(sd[, "FileName"])), which(!is.na(sd[, "identity"])))
    if (length(sd.delete.ind) > 0) {
        print(names(fcs.files[which(fcs.files %in% sd[sd.delete.ind, "identity"])]))
        choice <- utils::menu(c("Yes", "No"), title = "Move these FCS files to deleted_FCS_files and exclude them from sampledescription?")

        if (choice == 1) {
            dir.create(file.path(FCS.file.folder, "deleted_FCS_files"), showWarnings = F, recursive = T)

            if (dir.exists(file.path(FCS.file.folder, "deleted_FCS_files"))) {
                file.copy(names(fcs.files[which(fcs.files %in% sd[sd.delete.ind, "identity"])]), file.path(FCS.file.folder, "deleted_FCS_files", basename(names(fcs.files[which(fcs.files %in%
                  sd[sd.delete.ind, "identity"])]))))
                file.remove(names(fcs.files[which(fcs.files %in% sd[sd.delete.ind, "identity"])]))

                sd <- sd[which(!is.na(sd[, "FileName"])), ]
                sd[, "FileName"] <- ifelse(grepl("^[[:digit:]]{1,}_-_", sd[, "FileName"]), paste0(sprintf("%04d", 1:nrow(sd)), "_-_", substr(sd[, "FileName"],
                  8, nchar(sd[, "FileName"]))), paste0(sprintf("%04d", 1:nrow(sd)), "_-_", sd[, "FileName"]))
                fcs.files <- fcs.files[which(fcs.files %in% sd[, "identity"])]


                if ("initial.file.names" %in% names(openxlsx::loadWorkbook(file.path(wd, xlsx.file.name)))) {
                  initial.file.names <- openxlsx::read.xlsx(file.path(wd, xlsx.file.name), sheet = "initial.file.names")
                  named.sheet.list <- stats::setNames(list(sd, initial.file.names), c(sample.sheet.name, "initial.file.names"))
                } else {
                  named.sheet.list <- stats::setNames(list(sd), c(sample.sheet.name))
                }
                write.sd(named.sheet.list = named.sheet.list, wd = wd, xlsx.file.name = xlsx.file.name)

                print(paste0("FCS files removed and ", xlsx.file.name, " updated."))
            } else {
                print("deleted_FCS_files folder could not be created - no files were removed.")
            }
        }
        if (choice == 2) {
            return("No files removed.")
        }
    }

    fcs.files <- check.FCS.files(FCS.file.folder = FCS.file.folder, exclude.folders = exclude.folders)
    sd <- read.and.check.sd(wd = wd, xlsx.file.name = xlsx.file.name, sample.sheet.name = sample.sheet.name)

    # find new files for addition to sd
    fcs.files.diff <- fcs.files[which(!fcs.files %in% sd[, "identity"])]
    if (length(fcs.files.diff) != 0) {

        fcs.files.diff <- fcs.files.diff[order(lubridate::parse_date_time(sapply(strsplit(fcs.files.diff, "_-_"), "[", 3), orders = "%Y.%m.%d-%H.%M.%S", locale = "en_GB.UTF-8"))]

        sd.diff <- data.frame(FileName = paste0(sprintf(paste0("%04d"), (nrow(sd) + 1):(nrow(sd) + length(fcs.files.diff))), "_-_", basename(names(fcs.files.diff))),
            FilePath = names(fcs.files.diff), identity = fcs.files.diff, stringsAsFactors = FALSE)
        sd.diff[, c(names(sd)[which(!names(sd) %in% names(sd.diff))])] <- ""
        sd <- rbind(sd, sd.diff)

        if ("initial.file.names" %in% names(openxlsx::loadWorkbook(file.path(wd, xlsx.file.name)))) {
            initial.file.names <- rbind(openxlsx::read.xlsx(file.path(wd, xlsx.file.name), sheet = "initial.file.names"), sd.diff[, c("FilePath", "identity")])
            named.sheet.list <- stats::setNames(list(sd, initial.file.names), c(sample.sheet.name, "initial.file.names"))
        } else {
            named.sheet.list <- stats::setNames(list(sd), c(sample.sheet.name))
        }
        write.sd(named.sheet.list = named.sheet.list, wd = wd, xlsx.file.name = xlsx.file.name)
        print(paste0("New files added and renamed: ", paste(sd.diff[, "FileName"], collapse = ",")))
        file.rename(sd.diff[, "FilePath"], file.path(dirname(sd.diff[, "FilePath"]), sd.diff[, "FileName"]))
    }

    fcs.files <- check.FCS.files(FCS.file.folder = FCS.file.folder, exclude.folders = exclude.folders)
    sd <- read.and.check.sd(wd = wd, xlsx.file.name = xlsx.file.name, sample.sheet.name = sample.sheet.name)

    # find files for renaming
    sd.rename.ind <- which(!sd[, "FileName"] %in% basename(names(fcs.files)))

    if (length(sd.rename.ind) > 0) {

        sd[, "FileName"] <- ifelse(!grepl("^[[:digit:]]{1,}_-_", sd[, "FileName"]), paste0(sprintf("%04d", 1:nrow(sd)), "_-_", sd[, "FileName"]), sd[, "FileName"])
        sd[, "FileName"] <- ifelse(!grepl("\\.fcs$", tolower(sd[, "FileName"])), paste0(sd[, "FileName"], ".fcs"), sd[, "FileName"])
        sd[, "FileName"] <- sub("\\.FCS$", ".fcs", sd[, "FileName"])
        sd[, "FilePath"] <- file.path(dirname(names(fcs.files)), sd[, "FileName"])  # if folder is moved to another top-folder or another machine, FilePaths are changed, here already new paths are written

        fcs.files.rev <- stats::setNames(names(fcs.files), fcs.files)
        FileName = sd[order(sd[, "FileName"]), "FileName"][sd.rename.ind]
        FilePath.towrite <- fcs.files.rev[fcs.files %in% sd[order(sd[, "FileName"]), "identity"][sd.rename.ind]]

        if (any(!names(FilePath.towrite) == sd[order(sd[, "FileName"]), "identity"][sd.rename.ind])) {
            stop("File Identites do not match. Contact the code maintainer.")
        }
        print(data.frame(FileName = FileName, PreviousFileName = basename(FilePath.towrite)))
        choice <- utils::menu(c("Yes", "No"), title = "Rename FCS files as indicated?")

        if (choice == 1) {
            if ("initial.file.names" %in% names(openxlsx::loadWorkbook(file.path(wd, xlsx.file.name)))) {
                initial.file.names <- openxlsx::read.xlsx(file.path(wd, xlsx.file.name), sheet = "initial.file.names")
                named.sheet.list <- stats::setNames(list(sd, initial.file.names), c(sample.sheet.name, "initial.file.names"))
            } else {
                named.sheet.list <- stats::setNames(list(sd), c(sample.sheet.name))
            }
            write.sd(named.sheet.list = named.sheet.list, wd = wd, xlsx.file.name = xlsx.file.name)

            file.rename(FilePath.towrite, file.path(dirname(FilePath.towrite), FileName))
        }
        if (choice == 2) {
            return("No files renamed.")
        }
    }
}

write.sd.log <- function(wd, xlsx.file.name, sd) {
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

write.sd <- function(named.sheet.list, wd, xlsx.file.name) {
    tryCatch(openxlsx::write.xlsx(named.sheet.list, file = file.path(wd, xlsx.file.name), firstRow = T, colWidths = "auto", overwrite = T), error = function(e) {
        print(paste0("Is ", xlsx.file.name, " still open in Excel? Saving as updated file as ", file.path(wd, paste0(format(Sys.time(), "%Y.%m.%d-%H.%M.%S_"),
            xlsx.file.name)), ". Please delete the former one and remove the date-prefix of the new file."))
        openxlsx::write.xlsx(named.sheet.list, file = file.path(wd, paste0(format(Sys.time(), "%Y.%m.%d-%H.%M.%S_"), xlsx.file.name)), firstRow = T, colWidths = "auto")
    })
}

check.FCS.files <- function(FCS.file.folder, exclude.folders) {
    fcs.file.paths <- list.files(path = FCS.file.folder, pattern = "\\.fcs", full.names = TRUE, recursive = TRUE, ignore.case = T)
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
        stop(paste0("Duplicate FCS files found. This is not allowed. \n", paste(names(fcs.files[duplicated(fcs.files) | duplicated(fcs.files, fromLast = T)]),
            collapse = "\n")))
    }
    return(fcs.files)
}

read.and.check.sd <- function(wd, xlsx.file.name, sample.sheet.name) {
    ##### to do: check sd for consistency (e.g. FileName, FilePath, identity etc.)
    if (!sample.sheet.name %in% names(openxlsx::loadWorkbook(file.path(wd, xlsx.file.name)))) {
        stop(paste0("Sheet '", sample.sheet.name, "' not found in ", file.path(wd, xlsx.file.name), "."))
    }
    sd <- as.data.frame(openxlsx::read.xlsx(file.path(wd, xlsx.file.name), sheet = sample.sheet.name, skipEmptyCols = F, detectDates = T))
    if (purrr::is_empty(sd)) {
        stop(paste0("Please open ", xlsx.file.name, " with excel once and save again. This error has to do with previous writing of xlxs files with the writexl package. The function causing the error is read.xlsx from openxlsx which cannot handle the files written with writexl. Read.xlsx is superior though as it correctly detects dates in an excel sheet."))
    }
    if (!"FileName" %in% names(sd)) {
        stop(paste0("FileName column not found in ", xlsx.file.name, "."))
    }
    if (any(sapply(c("/", ":", "\\|", "\\?", "\\!", "\\*", "<", ">", "'", "\""), function(x) grepl(x, sd[, "FileName"])))) {
        stop("There is at least one FileName with one or more illegal character(s) which may cause problems in file-naming ( / : | ? ! * < > ' \")")
    }
    return(sd)
}
