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

df_2019 <- read.csv("bci_trees_tls_pointcloud_metrics.csv")

df_2019_quadrant <- read.csv("liana_census_all_trees_within_subplot_subqua_sorted.csv")
df_2019_quadrant <- data_frame(tag = df_2019_quadrant$Tag,
                               quadrant = df_2019_quadrant$Quadrant)

df_2019_new <- merge(df_2019, df_2019_quadrant, by = "tag" )

df_2019_unique <- distinct(df_2019_new, tag, .keep_all = TRUE)


df_2019 <- data_frame(tag = df_2019_unique$tag,
                      quadrant = df_2019_unique$quadrant,
                      lianaLoad = df_2019_unique$lianaLoad,
                      species = df_2019_unique$species,
                      dbh = df_2019_unique$dbh/10,
                      height = df_2019_unique$height,
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
                      quadrant = df_2015$q20,
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
                      quadrant = df_2011$q20,
                      species = df_2011$Full.Name,
                      dbh = df_2011$DBH_cm,
                      lianaLoad = df_2011$Lianas,
                      height = df_2011$Height,
                      source = "laser_rangefinder",
                      year = 2011)

df_all <- rbind(df_2011,df_2015,df_2019)
df_all$species <- trimws(df_all$species)

domSpeciesBCI <- read.csv("speciesTraitsBCIData.csv")


df_counts <- df_all %>%
  group_by(species, lianaLoad) %>%
  summarise(n = n(), .groups = "drop")

# 2) Pivot wider so each row is one species with columns for each category count
df_counts_wide <- df_counts %>%
  pivot_wider(
    names_from = lianaLoad, 
    values_from = n,
    values_fill = 0  # fill missing combos with 0
  )

# 3) Filter species with at least 5 observations in each category
df_counts_filtered <- df_counts_wide %>%
  filter(`0` >= 5, `1` >= 5, `2` >= 5)

# 4) Extract the species names that pass this criterion
sp_enough <- df_counts_filtered$species

df_all$lianaLoad <- factor(df_all$lianaLoad, ordered=FALSE, levels=c("0","1","2"))
df_all$species <- factor(df_all$species)

mod_height_simple <- brm(
  formula = log(height) ~ log(dbh) + lianaLoad + source + (1 |species) + (1|tag) + (1|quadrant),
  data = df_all,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

summary(mod_height_simple)

mod_height <- brm(
  formula = log(height) ~ log(dbh) + lianaLoad + source + (1 + lianaLoad |species) + (1|tag) + (1|quadrant),
  data = df_all,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

summary(mod_height)

mod_height_interaction <- brm(
  formula = log(height) ~ lianaLoad*log(dbh) + lianaLoad + source + (1 + lianaLoad |species) + (1|tag) + (1|quadrant),
  data = df_all,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

summary(mod_height_interaction)

ce_source <- conditional_effects(
  mod_height_interaction, 
  effects = "dbh:source",
  re_formula = NA
)
ce_source <- plot(ce_source, plot=FALSE)

df_plot <- ce_source[[1]]$data

# Exponentiate the relevant columns if your model is log(...) ~ ...
df_plot$estimate__ <- exp(df_plot$estimate__)
df_plot$lower__    <- exp(df_plot$lower__)
df_plot$upper__    <- exp(df_plot$upper__)

# 5) Update the ggplot's data with the exponentiated values
ce_source[[1]]$data <- df_plot

plot_ce_source_custom <- ce_source[[1]] +
  labs(
    x = "DBH (cm)",
    y = "Height (m)"
  ) +
  scale_color_manual(
    name = "Sensor", 
    values = c("tls"="#4682B4","drone_photogrammetry"="#B47846","laser_rangefinder"="#B4464B"), # or your own color palette
    labels = c("tls"="TLS", "drone_photogrammetry"="Drone Photogrammetry", "laser_rangefinder"="Rangefinder")  
  ) + scale_fill_manual(
    values = c("tls"="#4682B4", "drone_photogrammetry"="#B47846", "laser_rangefinder"="#B4464B"),
    guide = "none"       # <- Hide the fill legend
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size=14),
    legend.text =  element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    plot.title = element_text(size=14, face="bold")
  )

ggsave("source_effect_on_height.png", width = 8, height = 5, plot = plot_ce_source_custom, dpi = 600)

ce_dbh_liana <- conditional_effects(
  mod_height_interaction, 
  effects = "dbh:lianaLoad",   # or "lianaLoad:logDBH"
  re_formula = NA,
  method = "fitted"# ignore random effects for pop-level visualization
)


ce_dbh_liana <- plot(ce_dbh_liana, plot=FALSE)

df_plot <- ce_dbh_liana[[1]]$data

# Exponentiate the relevant columns if your model is log(...) ~ ...
df_plot$estimate__ <- exp(df_plot$estimate__)
df_plot$lower__    <- exp(df_plot$lower__)
df_plot$upper__    <- exp(df_plot$upper__)

# 5) Update the ggplot's data with the exponentiated values
ce_dbh_liana[[1]]$data <- df_plot


plot_ce_dbh_liana_custom <- ce_dbh_liana[[1]] +
  labs(
 
  x = "DBH (cm)",
  y = "Height (m)"
) +
  scale_color_manual(
    name = "Liana Load", 
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"), # or your own color palette
    labels = c("0"="No Liana", "1"="Medium", "2"="Heavy")  
  ) + scale_fill_manual(
    values = c("0"="#009e73", "1"="#D55E00", "2"="#000000"),
    guide = "none"       # <- Hide the fill legend
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size=14),
    legend.text =  element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    plot.title = element_text(size=14, face="bold")
  )


re_height <- ranef(mod_height_interaction)$species  # array
dimnames(re_height)
# Convert it to a data frame
df_height_int <- data.frame(
  species = rownames(re_height),
  est     = re_height[, "Estimate","Intercept" ],
  lower   = re_height[, "Q2.5","Intercept" ],
  upper   = re_height[, "Q97.5","Intercept" ]
)

# Subset to only species with >=10 observations
df_height_int_sub <- subset(df_height_int, species %in% sp_enough)

ggplot(df_height_int_sub, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Species-Specific Random Intercepts",
       x="Species", y="Deviation from Overall Intercept (log-height)")

###################### Species specific impact for heavy liana load category ########################

df_height_int_dom_sp <- subset(df_height_int_sub, species %in% domSpeciesBCI$full.name)

df_height_int_dom_sp_traits <- df_height_int_dom_sp %>% left_join(domSpeciesBCI, by = c("species" = "full.name"))

ggplot(df_height_int_dom_sp_traits, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Species-Specific Random Intercept",
       x="Species", 
       y="Deviation from Overall Intercept") 

df_sp_lianaload2_slope <- data.frame(
  species = rownames(re_height),
  est = re_height[,"Estimate","lianaLoad2"],
  lower = re_height[,"Q2.5","lianaLoad2"],
  upper = re_height[,"Q97.5","lianaLoad2"]
)

df_sp_lianaload2_slope_sub <- subset(df_sp_lianaload2_slope, species %in% sp_enough)

ggplot(df_sp_lianaload2_slope_sub, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Species-Specific Slopes for LianaLoad=2",
       x="Species", 
       y="Deviation from Population Slope")

df_sp_lianaload2_slope_dom_sp <- subset(df_sp_lianaload2_slope_sub , species %in% domSpeciesBCI$full.name)

df_sp_lianaload2_slope_dom_sp_traits <- df_sp_lianaload2_slope_dom_sp %>% left_join(domSpeciesBCI, by = c("species" = "full.name"))


ggplot(df_sp_lianaload2_slope_dom_sp_traits, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Species-Specific Slopes for LianaLoad=2",
       x="Species", 
       y="Deviation from Population Slope") 

###################### Correlation between traits and liana effect on height for heavy category ########################

df_sp_lianaload2_slope_dom_sp_traits_long <- df_sp_lianaload2_slope_dom_sp_traits %>%
  pivot_longer(
    cols = c("wsgbest","maxht", "RGR_100", "MORT_100"),  # list your trait columns here
    names_to = "traitName",
    values_to = "traitValue"
  )

df_sp_lianaload2_slope_dom_sp_traits_long <- df_sp_lianaload2_slope_dom_sp_traits_long %>% left_join(df_counts_filtered, by = "species")

df_sp_lianaload2_slope_dom_sp_traits_lim <- data.frame(name = df_sp_lianaload2_slope_dom_sp_traits_long$species,
                                              est = df_sp_lianaload2_slope_dom_sp_traits_long$est,
                                              lower95 = df_sp_lianaload2_slope_dom_sp_traits_long$lower,
                                              upper95 = df_sp_lianaload2_slope_dom_sp_traits_long$upper,
                                              traitName = df_sp_lianaload2_slope_dom_sp_traits_long$traitName,
                                              traitVal = df_sp_lianaload2_slope_dom_sp_traits_long$traitValue,
                                              freq = df_sp_lianaload2_slope_dom_sp_traits_long$`2`)



trait_names <- list(
  'wsgbest'="Wood Density",
  'maxht'="Max. Height",
  'RGR_100'="Rel. Growth Rate",
  'MORT_100'="Mortality Rate"
)

trait_labeller <- function(variable,value){
  return(trait_names[value])
}

df_stats <- df_sp_lianaload2_slope_dom_sp_traits_lim %>%
  group_by(traitName) %>%
  summarise(
    p_value = cor.test(traitVal, est, method = "spearman")$p.value
  ) %>%
  mutate(
    lineType = if_else(p_value < 0.005, "solid", "dotted")
  )

# Now merge (join) back to the main data frame
df_plot <- df_sp_lianaload2_slope_dom_sp_traits_lim %>%
  left_join(df_stats, by="traitName")

sp_height_impact_plot <- ggplot(df_plot, aes(x=traitVal, y=est)) +
  geom_point(aes(size=freq)) +
  geom_smooth(
    method = "lm", 
    aes(weight=freq, linetype=lineType), 
    se = TRUE
  ) +
  # you can still do stat_cor, though it won't affect line type directly
  stat_cor(
    aes(label = paste(..r.label.., sep="~`,`~")), 
    method = "spearman",  
    label.x.npc = "left",
    label.y = max(df_plot$est)
  ) +
  labs(
    size = "No. of individuals",
    x = "Trait value",
    y = "Effect of heavy liana load"
  ) +
  facet_wrap(~traitName, scales="free_x", labeller=trait_labeller) +
  scale_linetype_manual(values=c("dotted","solid"), guide = "none") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size=14),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    legend.text  = element_text(size=14),
    strip.text   = element_text(size=14)
  )


ggsave("sp_height_impact_combined_plot.png", width = 8, height = 7, plot = sp_height_impact_plot, dpi = 600)

###################### Species specific impact for medium liana load category ########################

df_sp_lianaload1_slope <- data.frame(
  species = rownames(re_height),
  est = re_height[,"Estimate","lianaLoad1"],
  lower = re_height[,"Q2.5","lianaLoad1"],
  upper = re_height[,"Q97.5","lianaLoad1"]
)

df_sp_lianaload1_slope_sub <- subset(df_sp_lianaload1_slope, species %in% sp_enough)

ggplot(df_sp_lianaload1_slope_sub, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Species-Specific Slopes for LianaLoad=1",
       x="Species", 
       y="Deviation from Population Slope")

###################### Correlation between traits and liana effect on height for medium category ########################

df_sp_lianaload1_slope_dom_sp <- subset(df_sp_lianaload1_slope_sub , species %in% domSpeciesBCI$full.name)

df_sp_lianaload1_slope_dom_sp_traits <- df_sp_lianaload1_slope_dom_sp %>% left_join(domSpeciesBCI, by = c("species" = "full.name"))


ggplot(df_sp_lianaload1_slope_dom_sp_traits, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Species-Specific Slopes for LianaLoad=1",
       x="Species", 
       y="Deviation from Population Slope") 

df_sp_lianaload1_slope_dom_sp_traits_long <- df_sp_lianaload1_slope_dom_sp_traits %>%
  pivot_longer(
    cols = c("wsgbest","maxht", "RGR_100", "MORT_100", "lmadiscbest"),  # list your trait columns here
    names_to = "traitName",
    values_to = "traitValue"
  )

df_sp_lianaload1_slope_dom_sp_traits_long <- df_sp_lianaload1_slope_dom_sp_traits_long %>% left_join(df_counts_filtered, by = "species")

df_sp_lianaload1_slope_dom_sp_traits_lim <- data.frame(name = df_sp_lianaload1_slope_dom_sp_traits_long$species,
                                                       est = df_sp_lianaload1_slope_dom_sp_traits_long$est,
                                                       lower95 = df_sp_lianaload1_slope_dom_sp_traits_long$lower,
                                                       upper95 = df_sp_lianaload1_slope_dom_sp_traits_long$upper,
                                                       traitName = df_sp_lianaload1_slope_dom_sp_traits_long$traitName,
                                                       traitVal = df_sp_lianaload1_slope_dom_sp_traits_long$traitValue,
                                                       freq = df_sp_lianaload2_slope_dom_sp_traits_long$`1`)

ggplot(df_sp_lianaload1_slope_dom_sp_traits_lim, aes(x=traitVal, y=est)) +
  geom_point(aes(size=freq)) +
  geom_smooth(method = "lm",aes(weight=freq)) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "spearman",  # or "spearman"
           label.x.npc = "left", # position the label
           label.y.npc = "top") +
  labs(size = "No. of individuals",x = "Trait Value", y ="Deviation from Population Slopes (LianaLoad=1)") +
  facet_wrap(~traitName, scales="free_x",labeller=trait_labeller) +
  theme_minimal() 

####################################### Population vs. species effects ###############################################

pop_df <- ce_dbh_liana[["dbh:lianaLoad"]] 
pop_df <- pop_df$data %>%
  rename(
    dbh_val = dbh,      # if the column is exactly "log(dbh)"
    lianaLoad = lianaLoad,
    estimate = estimate__,     # might be "estimate__"
    lower = lower__, 
    upper = upper__
  ) %>%
  mutate(which = "Community Avg.")

# 2) Plot the population-level result:
p_pop <- ggplot(pop_df, aes(x=dbh_val, y=estimate, color=liana_cat)) +
  geom_line() +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=liana_cat), alpha=0.2, color=NA) +
  labs(
    title="Community-Level: LianaLoad * log(DBH)",
    x="log(DBH)", # or "DBH" if it shows actual DBH
    y="Predicted log(Height)"  # or "Predicted Height" if you used transform=exp
  ) +
  theme_minimal()

dominant_species <- df_height_int_dom_sp$species

dbh_seq <- seq(min(pop_df$dbh_val), max(pop_df$dbh_val), by=5)

species_df_list <- list()

for(sp in dominant_species){
  
  # For each species, create a grid of DBH x lianaLoad
  tmp_grid <- expand.grid(
    dbh = dbh_seq,     # numeric DBH
    lianaLoad = c("0","1","2"),
    species = sp,
    source = "tls", # or pick typical
    tag = NA,        # if needed, or pick a baseline tag?
    quadrant = NA    # or a baseline quadrant?
  )
  # If your model uses 'log(dbh)' internally, we might do 'log(dbh)' as well:
  tmp_grid$`log(dbh)` <- log(tmp_grid$dbh)
  
  # Use 'fitted()' or 'predict()':
  pred_sp <- fitted(
    mod_height_interaction,
    newdata = tmp_grid,
    re_formula = ~(1+lianaLoad|species),  # keep species random slopes
    summary = TRUE
  )
  
  tmp_out <- cbind(tmp_grid, pred_sp) %>%
    as_tibble() %>%
    mutate(which = "Species")
  
  species_df_list[[sp]] <- tmp_out
}

species_df <- bind_rows(species_df_list)

species_df$which <- "Species"

pop_df2 <- pop_df %>%
  rename(dbh = dbh_val) %>%
  select(dbh, lianaLoad, estimate, lower, upper, which) 


species_df <- species_df %>% rename(estimate=Estimate, lower=Q2.5, upper=Q97.5)
species_df$estimate <- exp(species_df$estimate)
species_df <- species_df %>% select(dbh, lianaLoad, species, estimate, lower, upper, which) 

# Option A: replicate each row for all 12 species:
unique_species <- unique(species_df$species)
pop_df2_big <- pop_df2 %>%
  tidyr::expand_grid(species = unique_species)

plot_df <- bind_rows(pop_df2_big, species_df) 

ggplot(plot_df, aes(x=dbh, y=estimate, color=lianaLoad, linetype=which)) +
  geom_line() +
  scale_color_manual(
    name = "Liana Load", 
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"), # or your own color palette
    labels = c("0"="No Liana", "1"="Medium", "2"="Heavy")  
  )+
  scale_linetype_manual(
    name = "Effect",
    values = c("Species" = "solid", "Community Avg." = "dashed")
  ) +
  #geom_ribbon(
  #  aes(ymin=lower, ymax=upper, fill=lianaLoad), alpha=0.2, color=NA
  #) +
  facet_wrap(~ species, scales="free_y") +  # or "fixed" if you prefer
  labs(
    title="Height ~ LianaLoad * DBH",
    x="DBH (cm)",
    y="Predicted Height"
  ) +
  # If you want distinct fill for each lianaLoad but don't want to pick colors manually:
  # scale_color_discrete()  # or scale_color_brewer/palette
  # scale_fill_discrete() 
  # You do *not* have to specify color= for each species => we only have 3 categories for lianaLoad
  theme_minimal()


