library(fcexpr)
path <- new_exp(path = "~/Documents", name = "fcexpr_test")

ff_list <- ff_simulate2(npop = 8, dims = 10, m = 10)
fcs_paths <- file.path(path, "FCS_files", paste0(names(ff_list), ".fcs"))
purrr::map2(.x = ff_list, .y = fcs_paths, ~flowCore::write.FCS(x = .x, filename = .y))


# file.edit(file.path(path, "R_scripts", "R_Template.R"))
sync_sampledescription(FCS.file.folder = file.path(wd, "FCS_files"))


xx <- wsx_get_popstats(ws = "/Volumes/CMS_SSD_2TB 1/example_workspaces/no_group_noOrAndGates.wsp",
                       strip_data = F)
yy <- wsx_get_popstats2(ws = "/Volumes/CMS_SSD_2TB 1/example_workspaces/no_group_noOrAndGates.wsp",
                        more_gate_data = T)
comcols <- intersect(colnames(xx[["counts"]]), colnames(yy[["counts"]]))

out <- waldo::compare(xx[["counts"]][,comcols],
                      yy[["counts"]][,comcols])
length(out) # len zero for equality

colnames(xx[["counts"]])
colnames(yy2 <- yy[["counts"]])

out <- outout <- waldo::compare(xx[["counts"]][,comcols],
                      as.data.frame(yy[["counts"]][,comcols]))

tt <- dplyr::anti_join(xx[["counts"]], yy[["counts"]])
tt2 <- dplyr::anti_join(df2_sub, df1_sub)


sd_path <- file.path(wd, "sampledescription.xlsx")
keywords <- c("$OP", "$TOT")
keywords_to_sampledescription(sd_path, keywords)
