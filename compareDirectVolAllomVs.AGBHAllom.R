library(ordinal)
library(lme4)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(brms)  
library(dplyr)
library(tidyr)

setwd("/Users/sruthikp/Work/Analysis/BCI2019/Datasets")

all_subplot_trees <- read.csv("liana_census_all_trees_within_subplot_subqua_sorted.csv")

all_subplot_trees <- data_frame(tag = all_subplot_trees$Tag,
                                species = all_subplot_trees$Species,
                                subplot = all_subplot_trees$Subplot.number,
                                dbh = all_subplot_trees$DBH/10,
                                lianaLoad = all_subplot_trees$Lianas)

all_subplot_trees <- all_subplot_trees %>%
  mutate(lianaLoad = case_when(
    lianaLoad == 0 ~ 0,
    lianaLoad <= 50 ~ 1,
    lianaLoad > 50 ~ 2,
  ))

all_subplot_trees <- all_subplot_trees %>% filter(dbh >= 40)


all_subplot_trees$species <- factor(all_subplot_trees$species)
all_subplot_trees$lianaLoad <- as.factor(all_subplot_trees$lianaLoad)
all_subplot_trees$subplot <- as.factor(all_subplot_trees$subplot)

sp_wood_density <- read.csv("BCI_50ha_dbh_sp_WD.csv")
sp_wood_density <- data_frame(species = sp_wood_density$species,
                              wd = sp_wood_density$wd)
sp_wood_density <- distinct(sp_wood_density)

all_subplot_trees <- merge(all_subplot_trees, sp_wood_density, by = "species")

df_actual <- all_subplot_trees
df_no <- all_subplot_trees %>%
  mutate(lianaLoad = "0")

df_actual$lianaLoad <- factor(df_actual$lianaLoad, levels=c("0","1","2"))
df_no$lianaLoad     <- factor(df_no$lianaLoad, levels=c("0","1","2"))


# SCENARIO A: actual liana
#   allow_new_levels=TRUE if you have new subplots or new species not in training
post_mat_direct_allom <- posterior_predict(
  mod_twoStage_without_subplot,
  newdata = df_actual,
  ndraws = 6000,           # request 2000 draws
  allow_new_levels = TRUE, # if subplots/species are new
  sample_new_levels = "gaussian"  # optional if you want random effect sampling for new levels
)

post_mat_h_allom <- posterior_predict(
  mod_height_sp_liana_pred,
  newdata = df_actual,
  ndraws = 6000,           # request 2000 draws
  allow_new_levels = TRUE, # if subplots/species are new
  sample_new_levels = "gaussian"  # optional if you want random effect sampling for new levels
)

post_mat_direct_allom <- exp(post_mat_direct_allom)
post_mat_h_allom     <- exp(post_mat_h_allom)

wd_vec <- df_actual$wd   # length = nrow(df_actual) = #trees
# Make a matrix with same dimension
wd_mat <- matrix(wd_vec, nrow=nrow(post_mat_direct_allom), ncol=ncol(post_mat_direct_allom), byrow=TRUE)

dbh_vec <- df_actual$dbh
dbh_mat <- matrix(dbh_vec, nrow=nrow(post_mat_direct_allom), ncol=ncol(post_mat_direct_allom), byrow=TRUE)

post_biomass_direct_allom <- post_mat_direct_allom * wd_mat
post_biomass_h_allom <- 0.0559 * (wd_mat * (dbh_mat * dbh_mat) * post_mat_h_allom)

library(dplyr)
library(tidyr)

df_draws_direct <- as.data.frame(post_biomass_direct_allom) %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = -".draw",
    names_to="tree_index",
    values_to="biomass_direct"
  ) %>%
  mutate(tree_index = as.integer(gsub("V","", tree_index)))

df_draws_height <- as.data.frame(post_biomass_h_allom) %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = -".draw",
    names_to="tree_index",
    values_to="biomass_height"
  ) %>%
  mutate(tree_index = as.integer(gsub("V","", tree_index)))

df_direct_mean <- df_draws_direct %>%
  group_by(tree_index) %>%
  summarise(biomass_direct = mean(biomass_direct), .groups="drop")

df_height_mean <- df_draws_height %>%
  group_by(tree_index) %>%
  summarise(biomass_height = mean(biomass_height), .groups="drop")

df_compare <- df_direct_mean %>%
  inner_join(df_height_mean, by="tree_index")


library(DescTools)

ccc_result <- CCC(
  x = df_compare$biomass_direct,
  y = df_compare$biomass_height
)

ccc_val <- ccc_result$rho.c

vol_comp_ccc_plot <- ggplot(df_compare, aes(x=biomass_direct, y=biomass_height)) +
  geom_point(alpha=0.7) +
  # 1:1 line:
  geom_abline(slope=1, intercept=0, linetype="dashed", color="darkblue") +
  labs(
    x="Mean Biomass (kg) \n (Direct Allometry)",
    y="Mean Biomass (kg) \n (Height Allometry)"
  ) +
  theme_minimal() +
  theme(
        axis.title = element_text(size=14),
        legend.text =  element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

ccc_label <- paste0("CCC = ", round(ccc_val$est,3))

vol_comp_ccc_plot <- vol_comp_ccc_plot + annotate(
  "text",
  x = max(df_compare$biomass_direct)*0.15, # near the right side
  y = max(df_compare$biomass_height)*0.80, # near the top
  label = ccc_label,
  color = "black",
  size = 5
) +
  annotate(
    "text", 
    x = max(df_compare$biomass_direct)*0.17, 
    y = max(df_compare$biomass_height)*0.95, 
    label = "1:1 reference line", 
    size = 5, fontface = "italic", color = "darkblue"
  )
