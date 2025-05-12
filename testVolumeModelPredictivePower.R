library(ordinal)
library(lme4)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(brms)  
library(dplyr)
library(tidyr)

setwd("")

df_2019 <- read.csv("bci_trees_tls_all_model_qsmv2_4_metrics.csv")

df_2019 <- df_2019 %>% filter(dbh >= 37)

df_2019 <- df_2019 %>%
  mutate(lianaLoad = case_when(
    lianaLoad == 0 ~ 0,
    lianaLoad <= 50 ~ 1,
    lianaLoad > 50 ~ 2,
  ))

df_2019$lianaLoad <- factor(df_2019$lianaLoad, ordered=FALSE, levels=c("0","1","2"))
df_2019$species <- factor(df_2019$species)
df_2019$subplot <- factor(df_2019$subplot)


df_2019$tag <- factor(df_2019$tag)


mod_without_liana_without_subplot <- brm(
  formula = bf(
    log(totVol) ~ log(dbh)  + (1  | species) + (1|tag) 
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99),
  save_pars = save_pars(all=TRUE)
)

mod_without_liana <- brm(
  formula = bf(
    log(totVol) ~ log(dbh)  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99),
  save_pars = save_pars(all=TRUE)
)

# (3) Fit the model:
mod_twoStage_with_subplot <- brm(
  formula = bf(
    log(totVol) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99),
  save_pars = save_pars(all=TRUE)
)

# (3) Fit the model:
mod_twoStage_without_subplot <- brm(
  formula = bf(
    log(totVol) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) 
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99),
  save_pars = save_pars(all=TRUE)
)

####### Analyze the predictive power ############

vol_model_comparison <- loo_compare(loo(mod_without_liana_without_subplot), loo(mod_without_liana), loo(mod_twoStage_without_subplot), loo(mod_twoStage_with_subplot))
print(vol_model_comparison)
