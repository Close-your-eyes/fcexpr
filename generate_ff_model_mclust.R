library(fcexpr)
library(mclust)
ws <- "/Volumes/CMS_SSD_2TB/Experiment_data/20230222_TNF_IFNg_CD107a_induction_by_IL15_IFNb/FJ_workspaces/Exp_part_1_2_3_4_5_6.wsp"

# check samples, groups, populations
groups <- fcexpr::wsx_get_groups(ws)
paths <- fcexpr::wsx_get_fcs_paths(ws, filter_AllSamples = T)
samplenames <- basename(grep("AnKr", paths$Exp_part_1, value = T))
pops <- fcexpr::wsx_get_poppaths(ws, groups = "Exp_part_1")


# read flow frames
fflist <- fcexpr::wsp_get_ff(wsp = ws,
                             FCS.file.folder = "/Volumes/CMS_SSD_2TB/Experiment_data/20230222_TNF_IFNg_CD107a_induction_by_IL15_IFNb/FCS_files",
                             population = "Single Cells",
                             samples = samplenames[c(1,4)])
ff <- purrr::list_flatten(fflist[["flowframes"]], name_spec = "{outer}")
ff <- purrr::list_flatten(ff, name_spec = "{outer}")

# index 1: unstim, index 2: stim
ncol(ff@exprs)

# saveRDS(exprs_model[which(names(exprs_model) != "z")], file.path("/Users/vonskopnik/Documents/R_packages/fcexpr/inst", paste0(names(ff)[2], "_mclustmodel.rds")))

