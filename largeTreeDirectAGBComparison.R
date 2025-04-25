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
post_mat_actual <- posterior_predict(
  mod_twoStage_without_subplot,
  newdata = df_actual,
  ndraws = 6000,           # request 2000 draws
  allow_new_levels = TRUE, # if subplots/species are new
  sample_new_levels = "gaussian"  # optional if you want random effect sampling for new levels
)

# post_mat_actual is a matrix (#draws x #trees)
dim(post_mat_actual)

# SCENARIO B: no-liana
post_mat_no <- posterior_predict(
  mod_twoStage_without_subplot,
  newdata = df_no,
  ndraws = 6000,
  allow_new_levels = TRUE,
  sample_new_levels = "gaussian"
)
dim(post_mat_no)

post_mat_actual <- exp(post_mat_actual)
post_mat_no     <- exp(post_mat_no)

wd_vec <- df_actual$wd   # length = nrow(df_actual) = #trees
# Make a matrix with same dimension
wd_mat <- matrix(wd_vec, nrow=nrow(post_mat_actual), ncol=ncol(post_mat_actual), byrow=TRUE)

post_biomass_actual <- post_mat_actual * wd_mat

# Similarly for no-liana:
wd_vec_no <- df_no$wd
wd_mat_no <- matrix(wd_vec_no, nrow=nrow(post_mat_no), ncol=ncol(post_mat_no), byrow=TRUE)
post_biomass_no <- post_mat_no * wd_mat_no

# 2) Build a map data frame with subplot/species for each tree
df_tree_map <- df_no %>%
  mutate(tree_index = row_number()) %>%
  select(tree_index, subplot, species)

df_draws_no <- as.data.frame(post_biomass_no) %>%
  mutate(.draw = row_number()) %>%  # label each row by draw index
  pivot_longer(
    cols = -".draw",
    names_to = "tree_index", 
    values_to = "biomass"
  ) %>%
  mutate(tree_index = as.integer(gsub("V","", tree_index))) %>%
  # optionally join your df_tree_map if you want to see which subplot/species
  left_join(df_tree_map, by="tree_index")
# => columns: .draw, tree_index, biomass, subplot, species, etc.

df_draws_actual <- as.data.frame(post_biomass_actual) %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = -".draw",
    names_to = "tree_index",
    values_to = "biomass"
  ) %>%
  mutate(tree_index = as.integer(gsub("V","", tree_index))) %>%
  left_join(df_tree_map, by="tree_index")

###############################################################################
########## Site-level comparison #################
###############################################################################

df_site_no <- df_draws_no %>%
  group_by(.draw) %>%
  summarise(
    site_biomass = sum(biomass, na.rm=TRUE), 
    agbd_density = sum(biomass)/1000/5.28,
    .groups="drop"
  ) %>%
  mutate(scenario = "Liana-free")

df_site_actual <- df_draws_actual %>%
  group_by(.draw) %>%
  summarise(
    site_biomass = sum(biomass,na.rm=TRUE),
    agbd_density = sum(biomass)/1000/5.28,
    .groups="drop"
  ) %>%
  mutate(scenario="Liana-laden \n (Trees only)")

df_site_both <- bind_rows(df_site_no, df_site_actual)  #

# Make scenario a factor with no_liana first:
df_site_both$scenario <- factor(df_site_both$scenario,
                                levels=c("Liana-free","Liana-laden \n (Trees only)")
)

# We'll pivot wide:
df_no_wide <- df_site_no %>%
  rename(agb_no = agbd_density) %>%
  select(.draw, agb_no)

df_act_wide <- df_site_actual %>%
  rename(agb_act = agbd_density) %>%
  select(.draw, agb_act)

df_site_compare <- df_no_wide %>%
  inner_join(df_act_wide, by=".draw") %>%
  mutate(
    diff = agb_act - agb_no,      # e.g. negative if 'actual' is smaller than 'no'
    perc_diff = (agb_act - agb_no)/agb_no * 100
  )

# Summaries
p_neg <- mean(df_site_compare$diff < 0)    # e.g. 0.89
# coverage that excludes zero: 79% => means the central 79% interval doesn't cross zero
# (You already derived that from quantiles)
ci_78 <- quantile(df_site_compare$diff, c(0.11, 0.89))

mean_diff <- mean(df_site_compare$diff)    # average difference in mg/ha
mean_perc_diff <- mean(df_site_compare$perc_diff)
ci_78 <- quantile(df_site_compare$perc_diff, c(0.11, 0.89))
lab_str <- glue(
  "Displaced AGB = {round(mean_perc_diff,1)} % [79% CI: {round(ci_78[1],2)} to {round(ci_78[2],2)}]\n",
)

ggplot(df_site_both, aes(x=scenario, y=agbd_density,fill=scenario)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#E69F00", "#999999")) +
  labs(
    y="AGB difference (%)"
  ) +
  ylim(115,180) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  annotate("text",
           x=1.5,
           y=130,
           label=lab_str,
           size=4, hjust=0.5
  )

###############################################################################
########## Subplot-level comparison #################
###############################################################################

df_sub_no <- df_draws_no %>%
  group_by(.draw, subplot) %>%
  summarise(tree_biomass_no = sum(biomass), .groups="drop")  # 'biomass' in kg

df_sub_actual <- df_draws_actual %>%
  group_by(.draw, subplot) %>%
  summarise(tree_biomass_act = sum(biomass), .groups="drop")

subplot_area_ha <- 0.16  # adapt to your real subplot area

df_sub_no <- df_sub_no %>%
  mutate(
    tree_biomass_no_MgHa = tree_biomass_no / 1000 / subplot_area_ha
    # divides by 1000 to get Mg, then divides by 0.16 ha => Mg/ha
  )

df_sub_actual <- df_sub_actual %>%
  mutate(
    tree_biomass_act_MgHa = tree_biomass_act / 1000 / subplot_area_ha
  )

library(dplyr)

df_sub_compare <- df_sub_no %>%
  select(.draw, subplot, tree_biomass_no_MgHa) %>%
  inner_join(
    df_sub_actual %>%
      select(.draw, subplot, tree_biomass_act_MgHa),
    by=c(".draw","subplot")
  ) %>%
  mutate(
    displaced_biomass = tree_biomass_no_MgHa - tree_biomass_act_MgHa
  )

df_fraction_liana <- all_subplot_trees %>%
  group_by(subplot) %>%
  summarise(
    fraction_liana12 = mean(lianaLoad %in% c("1","2"))  # if lianaLoad is factor
    # or if lianaLoad is numeric, do: mean(lianaLoad > 0)
  )

df_sub_final <- df_sub_compare %>%
  # also join fraction_liana if you have it
  left_join(df_fraction_liana, by="subplot")%>%
  mutate(
    displaced_biomass_perc = ((tree_biomass_no_MgHa - tree_biomass_act_MgHa)/tree_biomass_no_MgHa)*100
  )

df_sub_summary <- df_sub_final %>%
  group_by(subplot, fraction_liana12) %>%
  summarise(
    mean_displaced = mean(displaced_biomass_perc),
    
    .groups="drop"
  )

ggplot(df_sub_summary, aes(x=fraction_liana12, y=mean_displaced)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "pearson",  # or "spearman"
           label.x.npc = "left", # position the label
           label.y.npc = "top") +
  labs(x="Fraction of Trees with Lianas", y="Mean displaced AGB (%)") +
  theme_minimal()
