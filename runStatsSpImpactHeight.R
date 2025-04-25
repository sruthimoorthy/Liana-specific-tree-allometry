library(ordinal)
library(lme4)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(brms)  
library(dplyr)
library(tidyr)
library(ggpubr)

setwd("/Users/sruthikp/Work/Analysis/BCI2019/Datasets")

df_2019 <- read.csv("bci_trees_tls_metrics.csv")

df_2019 <- data_frame(tag = df_2019$tag,
                      lianaLoad = df_2019$lianaLoad,
                      species = df_2019$species,
                      dbh = df_2019$dbh/10,
                      height = df_2019$height,
                      source = "tls",
                      year = 2019)

df_2019 <- df_2019 %>%
  mutate(lianaLoad = case_when(
    lianaLoad == 0 ~ 0,
    lianaLoad <= 50 ~ 1,
    lianaLoad > 50 ~ 2,
  ))


speciesCodes <- read.csv("speciesCodeBook.csv")

height_data_2015 <- read.csv("new_height_liana_data_2015_no_outliers_q20_2lianalevels.csv")

df_2015 <- height_data_2015 %>%
  left_join(speciesCodes, by = c("Species"="Code"))

df_2015 <- data_frame(tag = df_2015$Tag,
                      lianaLoad = df_2015$Lianas,
                      species = df_2015$Full.Name,
                      dbh = df_2015$DBH_cm,
                      height = df_2015$Height,
                      source = "drone_photogrammetry",
                      year = 2015)

height_data_2011 <- read.csv("new_height_liana_data_2011_no_outliers_q20_2lianalevels.csv")
height_data_2011$Species <- casefold(height_data_2011$Species, upper = FALSE)

df_2011 <- height_data_2011 %>%
  left_join(speciesCodes, by = c("Species"="Code"))

df_2011 <- data_frame(tag = df_2011$Tag,
                      species = df_2011$Full.Name,
                      dbh = df_2011$DBH_cm,
                      lianaLoad = df_2011$Lianas,
                      height = df_2011$Height,
                      source = "laser_rangefinder",
                      year = 2011)

df_all <- rbind(df_2011,df_2015,df_2019)
df_all$species <- trimws(df_all$species)
df_all$lianaLoad <- as.integer(df_all$lianaLoad)

df_all <- df_all %>%
  mutate(combLianaCat = ifelse(lianaLoad > 0, 1, 0))


domSpeciesBCI <- read.csv("speciesTraitsBCIData.csv")


df_counts <- df_all %>%
  group_by(species, combLianaCat) %>%
  summarise(n = n(), .groups = "drop")

# 2) Pivot wider so each row is one species with columns for each category count
df_counts_wide <- df_counts %>%
  pivot_wider(
    names_from = combLianaCat, 
    values_from = n,
    values_fill = 0  # fill missing combos with 0
  )

# 3) Filter species with at least 5 observations in each category
df_counts_filtered <- df_counts_wide %>%
  filter(`0` >= 10, `1` >= 10)

# 4) Extract the species names that pass this criterion
sp_enough <- df_counts_filtered$species

df_all$combLianaCat <- factor(df_all$combLianaCat, ordered=FALSE, levels=c("0","1"))
df_all$species <- factor(df_all$species)



mod_height <- brm(
  formula = log(height) ~ log(dbh) + combLianaCat + source + (1 + combLianaCat |species) + (1|tag),
  data = df_all,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 2
)

summary(mod_height)

ce_dbh_liana <- conditional_effects(
  mod_height, 
  effects = "dbh:combLianaCat",   # or "lianaLoad:logDBH"
  re_formula = NA                 # ignore random effects for pop-level visualization
)
plot(ce_dbh_liana, points=TRUE)

ce_source <- conditional_effects(
  mod_height, 
  effects = "dbh:source",
  re_formula = NA
)
plot(ce_source)

re_height <- ranef(mod_height)$species  # array
dimnames(re_height)

df_sp_combLianaCat_slope <- data.frame(
  species = rownames(re_height),
  est = re_height[,"Estimate","combLianaCat1"],
  lower = re_height[,"Q2.5","combLianaCat1"],
  upper = re_height[,"Q97.5","combLianaCat1"]
)

df_sp_combLianaCat_slope_sub <- subset(df_sp_combLianaCat_slope, species %in% sp_enough)

ggplot(df_sp_combLianaCat_slope_sub, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Species-Specific Slopes for LianaLoad",
       x="Species", 
       y="Deviation from Population Slope")

df_sp_lianaload_slope_dom_sp <- subset(df_sp_combLianaCat_slope_sub , species %in% domSpeciesBCI$full.name)

df_sp_lianaload_slope_dom_sp_traits <- df_sp_lianaload_slope_dom_sp %>% left_join(domSpeciesBCI, by = c("species" = "full.name"))



df_sp_lianaload_slope_dom_sp_traits_long <- df_sp_lianaload_slope_dom_sp_traits %>%
  pivot_longer(
    cols = c("wsgbest","maxht", "RGR_100", "MORT_100", "lmadiscbest"),  # list your trait columns here
    names_to = "traitName",
    values_to = "traitValue"
  )

df_sp_lianaload_slope_dom_sp_traits_lim <- data.frame(name = df_sp_lianaload_slope_dom_sp_traits_long$species,
                                                       est = df_sp_lianaload_slope_dom_sp_traits_long$est,
                                                       lower95 = df_sp_lianaload_slope_dom_sp_traits_long$lower,
                                                       upper95 = df_sp_lianaload_slope_dom_sp_traits_long$upper,
                                                       traitName = df_sp_lianaload_slope_dom_sp_traits_long$traitName,
                                                       traitVal = df_sp_lianaload_slope_dom_sp_traits_long$traitValue)


ggplot(df_sp_lianaload_slope_dom_sp_traits_lim, aes(x=traitVal, y=est)) +
  geom_point()+
  geom_smooth(method = "lm") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "spearman",  # or "spearman"
           label.x.npc = "left", # position the label
           label.y.npc = "top") +
  facet_wrap(~traitName, scales="free_x") +
  theme_minimal()

