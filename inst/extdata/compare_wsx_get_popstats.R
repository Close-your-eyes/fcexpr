# ws <- "/Users/vonskopnik/Desktop/Exp_part_20_21.wsp"
# ws <- "/Users/vonskopnik/Desktop/ExpPart_6_for_pub.wsp"
# ws <- "/Users/vonskopnik/Desktop/20231005_FJ_exp_wsp.wsp"

system.time(df1 <- fcexpr::wsx_get_popstats(ws, strip_data = F)[["counts"]])

system.time(df2 <- wsx_get_popstats2(ws, more_gate_data = T)[["counts"]])



df1_sub <- df1[,which(names(df1) %in% c("FileName", "PopulationFullPath", "Population", "Count"))]
df2_sub <- df2[,which(names(df2) %in% c("FileName", "PopulationFullPath", "Population", "Count"))]

tt <- dplyr::anti_join(df1_sub, df2_sub)
