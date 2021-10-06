## ---- include = F-----------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(knitr.kable.NA = '', width = 60)

## ----setup, eval=F----------------------------------------
#  install.packages("BiocManager")
#  BiocManager::install("flowCore")
#  BiocManager::install("flowWorkspace")
#  BiocManager::install("CytoML")
#  install.packages("devtools")
#  devtools::install_github("Close-your-eyes/fcexpr")
#  library(fcexpr)

## ---------------------------------------------------------
dir <- base::system.file("extdata", package = "fcexpr")
fcexpr::new_exp(path = dir, name = "my_experiment")

## ---------------------------------------------------------
wd <- file.path(dir, paste0(base::gsub("-", "", base::Sys.Date()), "_my_experiment"))
list.files(wd)
list.files(wd, recursive = T)

## ---- echo=F----------------------------------------------
utils::untar(base::system.file("extdata", "Part_1.tgz", package = "fcexpr"), exdir = file.path(wd, "FCS_files"))

## ---------------------------------------------------------
list.files(file.path(wd, "FCS_files", "Part_1"))

## ---------------------------------------------------------
fcexpr::sync_sampledescription(FCS.file.folder = file.path(wd, "FCS_files"), init.columns = c())

## ---------------------------------------------------------
list.files(file.path(wd))

## ---- echo=F----------------------------------------------
knitr::kable(openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx")))

## ---------------------------------------------------------
list.files(file.path(wd, "FCS_files", "Part_1"))

## ---- echo=F,include=F------------------------------------
sd <- openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx"))
sd[2,1] <- "File_2"
sd[5,1] <- "File_5"
openxlsx::write.xlsx(sd, file.path(wd, "sampledescription.xlsx"), overwrite = T)

## ---- echo=F----------------------------------------------
knitr::kable(openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx")))

## ---------------------------------------------------------
fcexpr::sync_sampledescription(FCS.file.folder = file.path(wd, "FCS_files"), write.log = F)

## ---- echo=F----------------------------------------------
knitr::kable(openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx")))

## ---------------------------------------------------------
list.files(file.path(wd, "FCS_files", "Part_1"))

## ---- echo=F,include=F------------------------------------
sd <- openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx"))
sd$Patient_ID <- sort(rep(1:5, 2))
sd$Staining <- rep(c("FMO", "full"), 5)
sd$Patient_sex <- c(rep("M",4),rep("F",6))
openxlsx::write.xlsx(sd, file.path(wd, "sampledescription.xlsx"), overwrite = T)

## ---- echo=F----------------------------------------------
knitr::kable(openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx")))

## ---- echo=F,include=F------------------------------------
sd <- openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx"))
sd <- sd[c(1,2,9,10,5,6,3,4,7,8),]
openxlsx::write.xlsx(sd, file.path(wd, "sampledescription.xlsx"), overwrite = T)

## ---- echo=F----------------------------------------------
knitr::kable(openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx")))

## ---------------------------------------------------------
fcexpr::sync_sampledescription(FCS.file.folder = file.path(wd, "FCS_files"), write.log = F)

## ---------------------------------------------------------
list.files(file.path(wd, "FCS_files", "Part_1"))

## ---- echo=F----------------------------------------------
utils::untar(base::system.file("extdata", "Part_2.tgz", package = "fcexpr"), exdir = file.path(wd, "FCS_files"))

## ---------------------------------------------------------
list.files(file.path(wd, "FCS_files"))
fcexpr::sync_sampledescription(FCS.file.folder = file.path(wd, "FCS_files"), write.log = F)

## ---- echo=F----------------------------------------------
knitr::kable(openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx")))

## ---- echo=F,include=F------------------------------------
sd <- openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx"))
sd$FileName <- ifelse(grepl("comp", sd$FileName), NA, sd$FileName)
openxlsx::write.xlsx(sd, file.path(wd, "sampledescription.xlsx"), overwrite = T)

## ---- echo=F----------------------------------------------
knitr::kable(openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx")))

## ---------------------------------------------------------
fcexpr::sync_sampledescription(FCS.file.folder = file.path(wd, "FCS_files"), write.log = F)

## ---- echo=F----------------------------------------------
knitr::kable(openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx")))

## ---- echo=F,include=F------------------------------------
sd <- openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx"))
sd$Patient_ID[11:14] <- sort(rep(6:7, 2))
sd$Staining[11:14] <- rep(c("FMO", "full"), 2)
sd$Patient_sex[11:14] <- c(rep("M",2),rep("F",2))
sd$FileName <- paste0("Pat", sd$Patient_ID, "_", sd$Staining)
openxlsx::write.xlsx(sd, file.path(wd, "sampledescription.xlsx"), overwrite = T)

## ---- echo=F----------------------------------------------
knitr::kable(openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx")))

## ---------------------------------------------------------
fcexpr::sync_sampledescription(FCS.file.folder = file.path(wd, "FCS_files"), write.log = F)

## ---- echo=F----------------------------------------------
knitr::kable(openxlsx::read.xlsx(file.path(wd, "sampledescription.xlsx")))

## ---------------------------------------------------------
list.files(file.path(wd, "FCS_files"), recursive = T)

## ---- echo=F----------------------------------------------
unlink(wd, recursive = T)

