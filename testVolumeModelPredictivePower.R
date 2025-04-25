library(ordinal)
library(lme4)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(brms)  
library(dplyr)
library(tidyr)

setwd("Datasets")

df_2019 <- read.csv("bci_trees_tls_final_qsmv2_4_metrics.csv")

df_2019 <- df_2019 %>%
  mutate(lianaLoad = case_when(
    lianaLoad == 0 ~ 0,
    lianaLoad <= 50 ~ 1,
    lianaLoad > 50 ~ 2,
  ))

df_2019$lianaLoad <- factor(df_2019$lianaLoad, ordered=FALSE, levels=c("0","1","2"))
df_2019$species <- factor(df_2019$species)
df_2019$subplot <- factor(df_2019$subplot)

M <- 20
df_long <- df_2019 %>%
  rowwise() %>%
  mutate(
    volume_draws = list(
      rnorm(M, mean=totVol, sd=totVolSD)  # or a custom distribution
    )
  ) %>%
  unnest(cols=volume_draws) %>%
  rename(volume=volume_draws) %>%
  # if you want to keep track of repeated draws per tree:
  mutate(drawID = row_number())  # optional

df_long$tag <- factor(df_long$tag)

mod_without_liana_without_subplot <- brm(
  formula = bf(
    log(volume) ~ log(dbh)  + (1  | species) + (1|tag) 
  ),
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99),
  save_pars = save_pars(all=TRUE)
)

mod_without_liana <- brm(
  formula = bf(
    log(volume) ~ log(dbh)  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99),
  save_pars = save_pars(all=TRUE)
)

# (3) Fit the model:
mod_twoStage_with_subplot <- brm(
  formula = bf(
    log(volume) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99),
  save_pars = save_pars(all=TRUE)
)

# (3) Fit the model:
mod_twoStage_without_subplot <- brm(
  formula = bf(
    log(volume) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) 
  ),
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99),
  save_pars = save_pars(all=TRUE)
)

####### Analyze the predictive power ############
best_volume_model <- loo_compare(loo(mod_without_liana_without_subplot), loo(mod_without_liana), loo(mod_twoStage_without_subplot), loo(mod_twoStage_with_subplot))

best_volume_model_two <- loo_compare(loo(mod_twoStage_without_subplot,moment_match=TRUE), loo(mod_without_liana,moment_match=TRUE))
waic(mod_without_liana)
waic(mod_twoStage_without_subplot)
waic(mod_twoStage_with_subplot)
waic(mod_without_liana_without_subplot)
