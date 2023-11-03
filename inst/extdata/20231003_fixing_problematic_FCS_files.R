## Problem with FCS files from MacsQuant at AG_Mashreghi
# Somehow a few files lead to crashing FlowJo 10.9.0. But only on my Mac, not on Windows.
# This could be solved by reading the FCS file with flowCore::read.FCS and just change the order of events in the exprs slot: fcs@exprs <- fcs@exprs[sample(1:nrow(fcs@exprs), nrow(fcs@exprs)),]
# But then writing the FCS file with flowCore::write.FCS failed because the spillover matrix in the $SPILLOVER keyword was a 0x0 matrix. This caused an error in a C++ script within flowCore::write.FCS.
# So, first a proper neutral matrix has to be written into this keyword with: fcexpr::compMat_neutral_to_fcs.


# FYI, for finding the error it was necessary to source all scripts from flowCore:
'r_files <- list.files("/Users/christopher.skopnik/Downloads/flowCore-devel/R", full.names = T)
for (i in r_files) {
  source(i)
}'
# load MatrixStats
'library(matrixStats)'
# and add flowCore:::to spill_to_string(mat, cols) in IO.R








'temp <- flowCore::read.FCS("/Volumes/CMS_SSD_2TB/2023_UriSeq/FCS_files/Session_07/0141_-_GS2329_07.16_pH5_PBS.fcs", emptyValue = F, truncate_max_range = F)
temp@exprs <- temp@exprs[sample(1:nrow(temp@exprs), nrow(temp@exprs)-2),]
r_files <- list.files("/Users/christopher.skopnik/Downloads/flowCore-devel/R", full.names = T)
for (i in r_files) {
  source(i)
}
library(matrixStats)
# add flowCore::: to spill_to_string(mat, cols) in IO.R
write.FCS(temp, "/Volumes/CMS_SSD_2TB/2023_UriSeq/FCS_files/Session_07/0141_-_GS2329_07.16_pH5_PBS_order_changed.fcs")

compMat_neutral_to_fcs(fcs_file_path = "/Users/christopher.skopnik/Desktop/0141_-_GS2329_07.16_pH5_PBS.fcs", compMat_keyword = "$SPILLOVER")

temp2 <- flowCore::read.FCS("/Volumes/CMS_SSD_2TB/2023_UriSeq/FCS_files/Session_06/0086_-_GS2329_06.3_PBS.fcs", emptyValue = F, truncate_max_range = F)
temp3 <- flowCore::read.FCS("/Users/christopher.skopnik/Desktop/0141_-_GS2329_07.16_pH5_PBS.fcs")
temp3@exprs <- temp3@exprs[sample(1:nrow(temp3@exprs), nrow(temp3@exprs)),]
flowCore::write.FCS(temp3, "/Users/christopher.skopnik/Desktop/0141_-_GS2329_07.16_pH5_PBS.fcs")
'


ls("package:flowCore")
lsf.str("package:flowCore")
ls(getNamespace("flowCore"))

install.packages("pacman")
pacman::p_funs(flowCore, TRUE)

# https://stackoverflow.com/questions/44696431/how-to-get-the-package-name-of-a-function-in-r

install.packages("sos")
sos::findFn("spill_to_string")
environmentName(environment("spill_to_string"))
environment(select)
find("p_del")
findFunction("p_del")

findAllFun <- function(f) {
  h <- help.search(paste0("^",f,"$"),agrep=FALSE)
  h$matches[,"Package"]
}

findAllFun(f = "string_to_spill")
findAllFun(f = "p_del")

