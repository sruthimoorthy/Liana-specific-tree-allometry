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

df_2019 <- read.csv("liana_census_all_trees_within_subplot_subqua_sorted.csv")

df_2019 <- data_frame(tag = df_2019$Tag,
                      lianaLoad = df_2019$Lianas,
                      species = df_2019$Species,
                      dbh = df_2019$DBH/10,
                      quadrant = df_2019$Quadrant,
                      year = 2019)

df_2019 <- df_2019 %>%
  mutate(lianaLoad = case_when(
    lianaLoad == 0 ~ 0,
    lianaLoad <= 50 ~ 1,
    lianaLoad > 50 ~ 2,
  ))

df_2019$lianaLoad <- factor(df_2019$lianaLoad, ordered=FALSE, levels=c("0","1","2"))

speciesCodes <- read.csv("speciesCodeBook.csv")

height_data_2015 <- read.csv("new_height_liana_data_2015_no_outliers_q20_2lianalevels.csv")

df_2015 <- height_data_2015 %>%
  left_join(speciesCodes, by = c("Species"="Code"))

df_2015 <- data_frame(tag = df_2015$Tag,
                      lianaLoad = df_2015$Lianas,
                      species = df_2015$Full.Name,
                      dbh = df_2015$DBH_cm,
                      quadrant = df_2015$q20,
                      year = 2015)

height_data_2011 <- read.csv("new_height_liana_data_2011_no_outliers_q20_2lianalevels.csv")
height_data_2011$Species <- casefold(height_data_2011$Species, upper = FALSE)

df_2011 <- height_data_2011 %>%
  left_join(speciesCodes, by = c("Species"="Code"))

df_2011 <- data_frame(tag = df_2011$Tag,
                      species = df_2011$Full.Name,
                      dbh = df_2011$DBH_cm,
                      lianaLoad = df_2011$Lianas,
                      quadrant = df_2011$q20,
                      year = 2011)

df_all <- rbind(df_2011,df_2015,df_2019)
df_all$species <- trimws(df_all$species)

domSpeciesBCI <- read.csv("speciesTraitsBCIData.csv")


counts <- table(df_all$species)
counts_df <- as.data.frame(counts)
colnames(counts_df) <- c("species", "Freq")

# Subset to species with >= 10 individuals
sp_enough <- subset(counts_df, Freq >= 10)$species

############################## Liana load presence analysis #########################

df_all$lianaLoad <- factor(df_all$lianaLoad, ordered=TRUE, levels=c("0","1","2"))
df_all$species <- factor(df_all$species)

mod_infest <- brm(
  formula = lianaLoad ~ log(dbh) + (1|species) + (1|tag) + (1|quadrant),
  data = df_all,
  family = cumulative("logit"),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, cores = 4
)

summary(mod_infest)

ce_dbh_lianaLoad <- conditional_effects(
  mod_infest, 
  effects = "dbh", 
  categorical = TRUE ,# or "lianaLoad:logDBH"
  re_formula = NA                 # ignore random effects for pop-level visualization
)

plot_ce_lianaload <- plot(ce_dbh_lianaLoad, plot=FALSE)

p_custom <- plot_ce_lianaload[[1]] +
  labs(
    x = "DBH (cm)",
    y = "Probability"
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
    plot.title = element_text(size=14, face="bold")
  )

ggsave("dbh_propensity_plot.png", width = 9, height = 5, plot = p_custom, dpi = 600)

ranef_infest <- ranef(mod_infest)$species

# Convert it to a data frame
df_infest_int <- data.frame(
  species = rownames(ranef_infest),
  est     = ranef_infest[, "Estimate",],
  lower   = ranef_infest[, "Q2.5",],
  upper   = ranef_infest[, "Q97.5",]
)

# Subset to only species with >=10 observations
df_infest_int_sub <- subset(df_infest_int, species %in% sp_enough)

ggplot(df_infest_int_sub, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Species Propensity for Higher Liana Load (Ordinal Model)",
    x = "Species",
    y = "Random Intercept (log-odds scale)"
  ) 

###################### Look at dominant species to link with traits ########################

df_infest_int_dom_sp <- subset(df_infest_int_sub, species %in% domSpeciesBCI$full.name)

df_infest_int_dom_sp_traits <- df_infest_int_dom_sp %>% left_join(domSpeciesBCI, by = c("species" = "full.name"))

df_infest_int_dom_sp_traits_long <- df_infest_int_dom_sp_traits %>%
  pivot_longer(
    cols = c("wsgbest","maxht", "RGR_100", "MORT_100"),  # list your trait columns here
    names_to = "traitName",
    values_to = "traitValue"
  )

df_infest_int_dom_sp_traits_long <- df_infest_int_dom_sp_traits_long %>% left_join(counts_df, by = "species")

df_infest_int_dom_sp_traits_lim <- data.frame(name = df_infest_int_dom_sp_traits_long$species,
                                              est = df_infest_int_dom_sp_traits_long$est,
                                              lower95 = df_infest_int_dom_sp_traits_long$lower,
                                              upper95 = df_infest_int_dom_sp_traits_long$upper,
                                              traitName = df_infest_int_dom_sp_traits_long$traitName,
                                              traitVal = df_infest_int_dom_sp_traits_long$traitValue,
                                              freq = df_infest_int_dom_sp_traits_long$Freq)

ggplot(df_infest_int_dom_sp_traits, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Species Propensity for Higher Liana Load (Ordinal Model)",
    x = "Species",
    y = "Random Intercept (log-odds scale)"
  ) 
  #+
  #scale_color_viridis_c(option="plasma") 

trait_names <- list(
  'wsgbest'="Wood Density",
  'maxht'="Max. Height",
  'RGR_100'="Rel. Growth Rate",
  'MORT_100'="Mortality Rate"
)

trait_labeller <- function(variable,value){
  return(trait_names[value])
}


propensity_sp_plot <- ggplot(df_infest_int_dom_sp_traits_lim, aes(x=traitVal, y=est)) +
  geom_point(aes(size=freq)) +
  geom_smooth(method = "lm",aes(weight=freq)) +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "spearman",  # or "spearman"
           label.x.npc = "left", # position the label
           label.y.npc = "top") +
  labs(size = "No. of individuals",x = "Trait value", y ="Random intercept (log-odds scale)") +
  facet_wrap(~traitName, scales="free_x",labeller=trait_labeller) +
  theme_minimal() 

ggsave("propensity_sp_plot.png", width = 9, height = 5, plot = propensity_sp_plot, dpi = 600)
######################################################################## Height analysis ####################################

