# ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/Exp_part_20_21.wsp"
# ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/ExpPart_6_for_pub.wsp"
# ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/20231005_FJ_exp_wsp.wsp" # error: ids are wrong, multiple ids for AndOrNodes - fixed!

ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/no_group_noOrAndGates.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_more_gate_types.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_different_gating_trees.wsp"

ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_1D_gates_range_and_2Sector.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_1D_gates_range_and_2Sector_with_NotGate.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/noOrAndGates_1D_gates_range_and_2Sector_with_NotGate_differentGatingTrees.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/Multiple_OrNodes_AndNodes_sameDims_sameGatingTrees.wsp" # one full path missing - check why; same id assigned to different AndNodes?!; check what happens if one of the multiple nodes is removed
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/Multiple_OrNodes_AndNodes_sameDims_differentGatingTrees.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/Multiple_OrNodes_AndNodes_NotNode_on_OrAndNodes_sameDims_sameGatingTrees.wsp" # NotNodes of OrNodes or AndNodes are different! - check if code lines are compatible with ordinary NotNodes
ws <- "/Users/vonskopnik/Desktop/example_workspaces/Multiple_OrNodes_AndNodes_NotNode_on_OrAndNodes_sameDims_sameGatingTrees.wsp"
ws <- "/Users/vonskopnik/Desktop/example_workspaces/Multiple_OrNodes_AndNodes_NotNode_on_OrAndNodes_with_children_sameDims_sameGatingTrees.wsp"
ws <- "/Users/vonskopnik/Desktop/example_workspaces/OrAndNodes_from_different_Dims_sameGatingTrees.wsp" # some gate types missing
ws <- "/Users/vonskopnik/Desktop/example_workspaces/Exp_part_20_21.wsp"
ws <- "/Volumes/CMS_SSD_2TB/example_workspaces/Complicated_OrAndGates_OrGate_at_diff_hierachies_sameGatingTree.wsp" # error

system.time(df1 <- fcexpr::wsx_get_popstats(ws, strip_data = F)[["counts"]])
system.time(df2 <- fcexpr::wsx_get_popstats2(ws, more_gate_data = T)[["counts"]])

## compare even more columns!

df1_sub <- df1[,which(names(df1) %in% c("FileName", "PopulationFullPath", "Population", "Count"))]#, "gate_id", "parentgate_id"))]
names(df1_sub)[which(names(df1_sub) == "gate_id")] <- "id"
names(df1_sub)[which(names(df1_sub) == "parentgate_id")] <- "parent_id"
df2_sub <- df2[,which(names(df2) %in% c("FileName", "PopulationFullPath", "Population", "Count"))]# "id", "parent_id"))]

tt <- dplyr::anti_join(df1_sub, df2_sub)
tt2 <- dplyr::anti_join(df2_sub, df1_sub)



dup_id <- unique(pop_df$id[duplicated(df2$id)])
df22 <- df2[which(df2$id %in% dup_id),]
df2$GateDef[which(is.null(df2$GateDef[[1]]))] <- NA
out <- df2$GateDef
df2$GateDef[which(sapply(df2$GateDef, is.null))] <- NA
df2$GateDef[[2]]
