library(ordinal)
library(lme4)
library(ggplot2)
library(ggeffects)
library(dplyr)
df <- results

df <- results %>%
  mutate(tag = if_else(grepl("DBH",  filename),sub("^[^_]*_([0-9]+)_.*", "\\1",  filename),sub("^[^_]*_(\\d+)\\.txt$", "\\1", filename) ),
         tag = as.integer(tag))

df_field <- read.csv("/Users/sruthikp/Work/Analysis/BCI2019/2019_TLS_data/tls_trees_with_subplot_info.csv")
# Make sure LianaLoad is an ordered factor
df_field <- data.frame(subplot = df_field$Subplot.no.,
                       tag = df_field$Tag,
                       DBH = df_field$DBH,
                       Species = df_field$Species,
                       LianaLoad = df_field$Lianas,
                       RegQuality = df_field$Reg..quality,
                       ptQuality = df_field$PointCloudQuality)

final_df <- merge(x = df, y = df_field, by = "tag", all = TRUE)

final_df <- data_frame(tag = final_df$tag,
                     dbh = final_df$DBH,
                     lianaLoad = final_df$LianaLoad,
                     species = final_df$Species,
                     subplot = final_df$subplot,
                     height = final_df$height,
                     crownArea = final_df$crown_area,
                     crownDepth = final_df$crown_depth,
                     crownVol = final_df$crown_vol,
                     regQuality = final_df$RegQuality,
                     ptQuality = final_df$ptQuality)

write.csv(final_df, file = "/Users/sruthikp/Work/Analysis/BCI2019/2019_TLS_data/bci_trees_tls_metrics.csv", row.names = FALSE)

################################# Final QSM Metrics File ######################


pt_struct_metrics_2019 <- read.csv("bci_trees_tls_pointcloud_metrics.csv")
qsm_struct_metrics_2019 <- read.csv("QSMv2_4_opt_model_mean_std_struct_metrics_v2.csv")

pt_struct_metrics_2019 <- data.frame(tag = pt_struct_metrics_2019$tag,
                                     dbh = pt_struct_metrics_2019$dbh/10,
                                     lianaLoad = pt_struct_metrics_2019$lianaLoad,
                                     species = pt_struct_metrics_2019$species,
                                     subplot = pt_struct_metrics_2019$subplot,
                                     ptQuality = pt_struct_metrics_2019$ptQuality)

qsm_struct_metrics_2019 <- data.frame(tag = qsm_struct_metrics_2019$ID,
                                      totVol = qsm_struct_metrics_2019$MeanTotalVolume,
                                      totVolSD = qsm_struct_metrics_2019$StdTotalVolume,
                                      totArea = qsm_struct_metrics_2019$MeanTotalArea,
                                      totAreaSD = qsm_struct_metrics_2019$StdTotalArea,
                                      trunkArea = qsm_struct_metrics_2019$MeanTrunkArea,
                                      trunkAreaSD = qsm_struct_metrics_2019$StdTrunkArea,
                                      trunkVol = qsm_struct_metrics_2019$MeanTrunkVolume,
                                      trunkVolSD = qsm_struct_metrics_2019$StdTrunkVolume,
                                      branchVol12 = qsm_struct_metrics_2019$MeanBranch_1_2_Volume,
                                      branchVol12SD = qsm_struct_metrics_2019$StdBranch_1_2_Volume,
                                      branchArea12 = qsm_struct_metrics_2019$MeanBranch_1_2_Area,
                                      branchArea12SD = qsm_struct_metrics_2019$StdBranch_1_2_Area,
                                      branchLen12 = qsm_struct_metrics_2019$MeanBranch_1_2_Len,
                                      branchLen12SD = qsm_struct_metrics_2019$StdBranch_1_2_Len,
                                      branchNo12 = qsm_struct_metrics_2019$MeanBranch_1_2_No,
                                      branchNo12SD = qsm_struct_metrics_2019$StdBranch_1_2_No)

df_2019 <- merge(pt_struct_metrics_2019, qsm_struct_metrics_2019, by = "tag")

df_2019 <- df_2019 %>% filter(ptQuality == "ok")

write.csv(df_2019, file = "/Users/sruthikp/Work/Analysis/BCI2019/Datasets/bci_trees_tls_final_qsmv2_4_metrics.csv", row.names = FALSE)

################# Liana biomass data subplot-level ###################

all_subplot_trees_2019 <- read.csv("liana_census_all_trees_within_subplot_subqua_sorted.csv")
all_50ha_liana <- read.csv("2025-03-21_liana_per_quad_agb_2017.csv")

all_subplot_trees_2019 <- all_subplot_trees_2019 %>% select(Subplot.number, Quadrant)

quad_to_subplot <- all_subplot_trees_2019 %>%
  distinct(Quadrant, Subplot.number)


all_50ha_liana_with_subplot <- all_50ha_liana %>%
  left_join(quad_to_subplot, by = c("quad" = "Quadrant"))



all_50ha_liana_with_subplot <- all_50ha_liana_with_subplot %>% filter(!is.na(Subplot.number))

all_50ha_liana_with_subplot <- all_50ha_liana_with_subplot %>%
  group_by(Subplot.number) %>%
  summarise(total_liana_agb = sum(X2017_liana_agb_total))

all_50ha_liana_with_subplot <- data_frame(subplot = all_50ha_liana_with_subplot$Subplot.number,
                                     lianaAGB = all_50ha_liana_with_subplot$total_liana_agb)
write.csv(all_50ha_liana_with_subplot, file = "/Users/sruthikp/Work/Analysis/BCI2019/Datasets/liana_agb_subplot_level.csv", row.names = FALSE)


