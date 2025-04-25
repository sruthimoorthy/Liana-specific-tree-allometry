library(ordinal)
library(lme4)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(brms)  
library(dplyr)

setwd("Datasets")

all_tree_data_2019 <- read.csv("liana_census_all_trees_within_subplot_subqua_sorted.csv")

df_2019 <- data_frame(tag = all_tree_data_2019$Tag,
                      lianaLoad = all_tree_data_2019$Lianas,
                      species = all_tree_data_2019$Species,
                      dbh = all_tree_data_2019$DBH,
                      year = 2019)

df_2019 <- df_2019 %>%
  mutate(lianaLoad = case_when(
    lianaLoad == 0 ~ 0,
    lianaLoad <= 50 ~ 1,
    lianaLoad > 50 ~ 2,
  ))

df_2019$lianaLoad <- factor(df_2019$lianaLoad, ordered=TRUE, levels=c("0","1","2"))
df_2019$species <- factor(df_2019$species)
df_2019$dbh <- df_2019$dbh/10

speciesCodes <- read.csv("speciesCodeBook.csv")

height_data_2015 <- read.csv("new_height_liana_data_2015_no_outliers_q20_2lianalevels.csv")

df_2015 <- height_data_2015 %>%
  left_join(speciesCodes, by = c("Species"="Code"))

df_2015 <- data_frame(tag = df_2015$Tag,
                      species = df_2015$Full.Name,
                      dbh = df_2015$DBH_cm,
                      lianaLoad = df_2015$Lianas,
                      year = 2015
                      )

height_data_2011 <- read.csv("new_height_liana_data_2011_no_outliers_q20_2lianalevels.csv")
height_data_2011$Species <- casefold(height_data_2011$Species, upper = FALSE)

df_2011 <- height_data_2011 %>%
  left_join(speciesCodes, by = c("Species"="Code"))

df_2011 <- data_frame(tag = df_2011$Tag,
                      species = df_2011$Full.Name,
                      dbh = df_2011$DBH_cm,
                      lianaLoad = df_2011$Lianas,
                      year = 2011
)

df_lianaload <- rbind(df_2011,df_2015,df_2019)
df_lianaload$species <- trimws(df_lianaload$species)


df_lianaload$lianaLoad <- factor(df_lianaload$lianaLoad, ordered=TRUE, levels=c("0","1","2"))
df_lianaload$species <- factor(df_lianaload$species)
df_lianaload <- df_lianaload %>%
  mutate(dbh = ifelse(year == 2019, dbh/10, dbh))

df_lianaload <- na.omit(df_lianaload)

df_filtered <- df_lianaload %>%
  group_by(species) %>%
  filter(n() >= 10) %>%
  ungroup()

df_filtered$species <- factor(df_filtered$species)

mod_brms <- brm(
  formula = lianaLoad ~ log(dbh) + (1 | species) + (1|tag), 
  data = df_filtered,
  family = cumulative(link = "logit"), 
  # MCMC settings:
  chains = 4, 
  cores = 4, 
  iter = 3000,
  # (Optional) set some mildly informative priors
  # prior = c(set_prior("normal(0,5)", class="b"))
)

summary(mod_brms)

ce <- conditional_effects(mod_brms, "dbh", categorical=TRUE)
plot(ce)


re_species <- ranef(mod_brms)$species

# Suppose we want the median (or mean) effect (column 1) and the lower & upper 95% CI (columns 3 and 4).
# The exact column indexing might differ by brms version. Check names with dimnames(re_species).

re_df <- data.frame(
  species = rownames(re_species[, , 1]),             # species names
  intercept_est = re_species[, 1, ],       # posterior mean/median
  intercept_lower = re_species[, 1, ],         # 2.5% quantile
  intercept_upper = re_species[, 1, ]         # 97.5% quantile
)

library(ggplot2)

ggplot(re_df, aes(x = reorder(species, intercept_est), y = intercept_est)) +
  geom_point() +
  geom_errorbar(aes(ymin = intercept_lower, ymax = intercept_upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Species Random Intercepts (Posterior Estimates)",
       x = "Species",
       y = "Random Intercept Estimate")


library(tidybayes)
post_re <- mod_brms %>%
  spread_draws(r_species[species,Intercept]) 
# This yields one row per MCMC draw per species with the random intercept.

head(post_re)
post_re_summary <- post_re %>%
  median_qi(.width = c(.95))  # 95% credible intervals

# Now plot
ggplot(post_re_summary, aes(x = reorder(species, r_species), y = r_species)) +
  geom_pointinterval(aes(ymin = .lower, ymax = .upper)) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Random Intercepts by Species",
       x = "Species", y = "Posterior Median (Random Intercept)")
