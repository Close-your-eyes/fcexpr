# ws <- "/Users/vonskopnik/Desktop/Exp_part_20_21.wsp"
# ws <- "/Users/vonskopnik/Desktop/ExpPart_6_for_pub.wsp"
# ws <- "/Users/vonskopnik/Desktop/20231005_FJ_exp_wsp.wsp"

ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/no_group_noOrAndGates.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_more_gate_types.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_different_gating_trees.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_1D_gates_range_and_2Sector.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_1D_gates_range_and_2Sector_with_NotGate.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_1D_gates_range_and_2Sector_with_NotGate_differentGatingTrees.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/Multiple_OrNodes_AndNodes_sameDims_sameGatingTrees.wsp" # fail

system.time(df1 <- fcexpr::wsx_get_popstats(ws, strip_data = F)[["counts"]])

system.time(df2 <- wsx_get_popstats2(ws, more_gate_data = T)[["counts"]])



df1_sub <- df1[,which(names(df1) %in% c("FileName", "PopulationFullPath", "Population", "Count"))]
df2_sub <- df2[,which(names(df2) %in% c("FileName", "PopulationFullPath", "Population", "Count"))]

tt <- dplyr::anti_join(df1_sub, df2_sub)

df2$GateDef[which(is.null(df2$GateDef[[1]]))] <- NA
out <- df2$GateDef
df2$GateDef[which(sapply(df2$GateDef, is.null))] <- NA
df2$GateDef[[2]]
