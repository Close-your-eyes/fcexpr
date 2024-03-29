#' Complement the antibody panel table with additional information from your antibody list
#'
#' This function matches the antibodies in a panel to those listed in an antibody list.
#' The antibody list may contain additional information like clone or catalog number which are to
#' be added to the panel table for proper annotation of the experiment.
#' The panel_file should be a xlsx-file as found in the Protocols-folder when running fcexpr::new_exp().
#' Necessary columns in panel_file and antibody_list used for matching are Antigen, Conjugate, Box, Lot.
#' If they do not exist this function cannot be used, yet. Columns added to the panel_file will appear just right to
#' the LiveDeadMarker column. Anything written there will be overwritten.
#'
#' The following typos in Antigen and Conjugate are rescued (matching may still work with the ab list):
#' Space, special characters (e.g. hyphen), case.
#'
#' @param panel_file character, path to the xlsx file in which the antibody panel was calculated
#' @param panel_sheet character, the name of the sheet in the panel_file
#' @param antibody_list character, path to the xlsx file representing the antibody list with additional information for the used antibodies
#' @param antibody_list_cols character vector, which columns from antibody_list to append to the panel_sheet in panel_file
#' @param antibody_list_sheet character or numeric, the name or index of the sheet in the antibody_list
#'
#' @return No return value, but an appended panel_file written to disk.
#' @export
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' ab_info_to_panel(panel_file = "my.path/20220202_Exp.part.1.xlsx",
#' antibody_list = "/Volumes/AG_Hiepe/_AG-HIEPE_Common/Antibody_List/20200705_antibody_list.xlsx")
#' }
ab_info_to_panel <- function(panel_file,
                             panel_sheet = "Panel",
                             antibody_list,
                             antibody_list_sheet = 1,
                             antibody_list_cols = c("Reactivity", "Isotype", "Clone", "Vendor", "Cat", "Expiry.date", "Concentration.ug.ml", "Recomm.dilution")) {

  if (missing(panel_file)) {
    stop("Please provide a panel_file.")
  }
  if (!grepl("xlsx$", panel_file)) {
    stop("panel_file has to be an xlsx file.")
  }
  if (basename(panel_file) == panel_file) {
    panel_file <- file.path(getwd(), panel_file)
    message(paste0("panel_file set to ", panel_file))
  }
  if (!file.exists(panel_file)) {
    stop("panel_file not found.")
  }
  panel <- openxlsx::read.xlsx(panel_file, sheet = panel_sheet, skipEmptyRows = F)
  if (any(!c("Antigen", "Conjugate", "Box", "Lot") %in% names(panel))) {
    stop("At least one of the required columns Antigen, Conjugate, Box, Lot is not found in panel_sheet of panel_file.")
  }
  panel <- panel[,c("Antigen", "Conjugate", "Box", "Lot")]

  if (missing(antibody_list)) {
    if (Sys.info()[["sysname"]] == "Darwin") {
      message("Trying to gues antibody_list location: '/Volumes/AG_Hiepe/_AG-HIEPE_Common/Antibody_List/20200705_antibody_list.xlsx'")
    }
    antibody_list <- "/Volumes/AG_Hiepe/_AG-HIEPE_Common/Antibody_List/20200705_antibody_list.xlsx"
  }

  if (missing(antibody_list)) {
    if (Sys.info()[["sysname"]] == "Darwin") {
      antibody_list <- "/Volumes/AG_Hiepe/_AG-HIEPE_Common/Antibody_List/20200705_antibody_list.xlsx"
    } else if (Sys.info()[["sysname"]] == "Windows") {
      antibody_list <- "Y:/AG_Hiepe/_AG-HIEPE_Common/Antibody_List/20200705_antibody_list.xlsx"
    }
    if (file.exists(antibody_list)) {
      ans <- readline(prompt = paste0("Use this antibody list as reference: ", antibody_list, "? y/n ?."))
      if (ans == "n") {
        stop("no antibody_list provided")
      } else if (ans == "y") {
        message("antibody list found.")
      } else {
        stop("Please type y for yes or n for n or provide the path to the antibody list manually.")
      }
    }
  }

  if (!grepl("xlsx$", antibody_list)) {
    stop("antibody_list has to be an xlsx.")
  }
  ab.list <- openxlsx::read.xlsx(antibody_list, sheet = antibody_list_sheet)
  if (any(!c("Antigen", "Conjugate", "Box", "Lot") %in% names(ab.list))) {
    stop("At least one of the required columns Antigen, Conjugate, Box, Lot is not found in the antibody_list_sheet of antibody_list.")
  }

  panel$row.num <- 1:nrow(panel)
  panel <- panel[which(!is.na(panel$Antigen)),]

  panel_add <-
    panel %>%
    dplyr::mutate(antigen_case = make.names(gsub(" ", "", tolower(Antigen))), conjugate_case = make.names(gsub(" ", "", tolower(Conjugate)))) %>%
    dplyr::left_join(ab.list %>% dplyr::mutate(antigen_case = make.names(gsub(" ", "", tolower(Antigen))), conjugate_case = make.names(gsub(" ", "", tolower(Conjugate)))), by = c("antigen_case", "conjugate_case")) %>%
    dplyr::select(-c(antigen_case, conjugate_case)) %>%
    #fuzzyjoin::regex_left_join(ab.list, by = c("Antigen", "Conjugate"), ignore_case = T) %>% # this would match CD45RO and CD4 - not useful
    dplyr::filter(!is.na(Antigen.y) & !is.na(Conjugate.y)) %>% # take names from Ab list
    dplyr::select(-c(Antigen.x, Conjugate.x)) %>%
    dplyr::rename("Antigen" = Antigen.y, "Conjugate" = Conjugate.y) %>%
    dplyr::filter(is.na(Box.x) | Box.x == Box.y) %>%
    dplyr::mutate(Box = ifelse(!is.na(Box.x), Box.x, Box.y)) %>%
    dplyr::filter(is.na(Lot.x) | Lot.x == Lot.y) %>%
    dplyr::mutate(Lot = ifelse(!is.na(Lot.x), Lot.x, Lot.y)) %>%
    dplyr::select(Antigen, Conjugate, Box, Lot, dplyr::all_of(antibody_list_cols), row.num) %>%
    dplyr::arrange(row.num) %>%
    dplyr::distinct()

  if (nrow(panel_add) == 0) {
    stop("No antibodies could be matched. Check for correct Antigen, Conjugate and Box as indicated in the ab_list.")
  }

  if (any(duplicated(panel_add$row.num))) {
    print(as.data.frame(panel_add))
    stop("Ambiguous entries found. Please check and make specific by providing a Box and/or Lot.")
  }

  if (length(setdiff(min(panel_add$row.num):max(panel_add$row.num), panel_add$row.num)) > 0) {
    panel_add <-
      panel_add %>%
      dplyr::bind_rows(data.frame(row.num = setdiff(min(panel_add$row.num):max(panel_add$row.num), panel_add$row.num))) %>%
      dplyr::arrange(row.num)
  }

  #panel_add <- dplyr::select(panel_add, -row.num)
  panel_add <- tidyr::drop_na(panel_add, Antigen)

  # avoid ablating formulas in cells - hence load wb and append
  wb <- openxlsx::loadWorkbook(panel_file)
  orig.sheet <- openxlsx::read.xlsx(panel_file, skipEmptyCols = F, skipEmptyRows = F, sheet = panel_sheet)

  for (y in seq_along(panel_add$row.num)) {
    openxlsx::writeData(wb, sheet = panel_sheet, xy = c(which(names(orig.sheet) == "Box"), panel_add$row.num[y]+1), x = panel_add[y,"Box", drop = T])
    openxlsx::writeData(wb, sheet = panel_sheet, xy = c(which(names(orig.sheet) == "Lot"), panel_add$row.num[y]+1), x = panel_add[y,"Lot", drop = T])
    for (z in seq_along(antibody_list_cols)) {
      openxlsx::writeData(wb, sheet = panel_sheet, xy = c(which(colnames(orig.sheet) == "LiveDeadMarker")+2+z,panel_add$row.num[y]+1), x = panel_add[y,antibody_list_cols[z]])
    }
  }
  for (z in seq_along(antibody_list_cols)) {
    openxlsx::writeData(wb, sheet = panel_sheet, xy = c(which(colnames(orig.sheet) == "LiveDeadMarker")+2+z,1), x = antibody_list_cols[z])
  }


  #openxlsx::writeData(wb, sheet = panel_sheet, xy = c(which(names(orig.sheet) == "Box"), 1), x = as.data.frame(panel_add[,"Box", drop = F]))
  #openxlsx::writeData(wb, sheet = panel_sheet, xy = c(which(names(orig.sheet) == "Lot"), 1), x = as.data.frame(panel_add[,"Lot", drop = F]))
  #openxlsx::writeData(wb, sheet = panel_sheet, xy = c(which(colnames(orig.sheet) == "LiveDeadMarker")+1,1), x = as.data.frame(panel_add[,antibody_list_cols]))
  openxlsx::saveWorkbook(wb, panel_file, overwrite = T)

  print(as.data.frame(panel_add))
  message(paste0("Saved as ", panel_file, "."))
}
