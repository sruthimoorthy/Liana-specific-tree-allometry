library(ordinal)
library(lme4)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(brms)  
library(dplyr)
library(tidyr)

setwd("Datasets")


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

df_all <- rbind(df_2015,df_2019)
df_all$species <- trimws(df_all$species)


df_all$lianaLoad <- factor(df_all$lianaLoad, ordered=FALSE, levels=c("0","1","2"))
df_all$species <- factor(df_all$species)

########################### log-log linear form ######################

mod_height_base_pred <- brm(
  formula = log(height) ~ log(dbh)  + (1 |species),
  data = df_all,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

mod_height_liana_pred <- brm(
  formula = log(height) ~ log(dbh) + lianaLoad  + (1 |species),
  data = df_all,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

mod_height_liana_inter_pred <- brm(
  formula = log(height) ~ log(dbh) * lianaLoad  + (1 |species),
  data = df_all,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

mod_height_sp_liana_pred <- brm(
  formula = log(height) ~ log(dbh) + lianaLoad  + (1 + lianaLoad |species),
  data = df_all,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4,
  save_pars = save_pars(all=TRUE)
)

mod_height_sp_liana_inter_pred <- brm(
  formula = log(height) ~ log(dbh) * lianaLoad  + (1 + lianaLoad |species),
  data = df_all,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

######################### Michealis-menten form #################################

############## BRMS model formula and prior ##############

mm_form1 <- brmsformula(log(height) ~ (Vm + x * log(dbh))-log(K+(dbh^x)), Vm + x + K ~ 1 + (1|species), nl=TRUE)

mm_form2 <- brmsformula(log(height) ~ (Vm + x * log(dbh))-log(K+(dbh^x)), Vm ~ 1 + lianaLoad + (1|species), x + K ~ 1 + (1|species), nl=TRUE)

mm_form3 <- brmsformula(log(height) ~ (Vm + x * log(dbh))-log(K+(dbh^x)),  x ~ 1 + lianaLoad + (1|species), Vm + K ~ 1 + (1|species), nl=TRUE)

mm_form4 <- brmsformula(log(height) ~ (Vm + x * log(dbh))-log(K+(dbh^x)),  K ~ 1 + lianaLoad + (1|species), Vm + x ~ 1 + (1|species), nl=TRUE)

mm_form5 <- brmsformula(log(height) ~ (Vm + x * log(dbh))-log(K+(dbh^x)),  x + K ~ 1 + lianaLoad + (1|species), Vm  ~ 1 + (1|species), nl=TRUE)

mm_form6 <- brmsformula(log(height) ~ (Vm + x * log(dbh))-log(K+(dbh^x)),  Vm + K ~ 1 + lianaLoad + (1|species), x  ~ 1 + (1|species), nl=TRUE)

mm_form7 <- brmsformula(log(height) ~ (Vm + x * log(dbh))-log(K+(dbh^x)),  Vm + x ~ 1 + lianaLoad + (1|species), K  ~ 1 + (1|species), nl=TRUE)

mm_form8 <- brmsformula(log(height) ~ (Vm + x * log(dbh))-log(K+(dbh^x)),  Vm + x + K ~ 1 + lianaLoad + (1|species), nl=TRUE)


####### Prior set from Martinnez Cano et. al. (2019) paper #########

prior_1 <- c(set_prior("normal(4, 3.5)",  nlpar = "Vm"),
             set_prior("normal(1,1)", nlpar = "x"),
             set_prior("normal(35,25)",  nlpar = "K"))

prior_2 <- c(set_prior("cauchy(0,3.5)", class = "sd",nlpar = "Vm"),
             set_prior("cauchy(0,3.5)",class = "sd", nlpar = "x"),
             set_prior("cauchy(0,3.5)",class = "sd", nlpar = "K"))

prior_mm_fit1 <- rbind(prior_1, prior_2)

mm_fit1 = brm(mm_form1, data=df_all, chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.95,max_treedepth = 20), cores = 4)
mm_fit2 = brm(mm_form2, data=df_all,prior=prior_mm_fit1, chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.9), cores = 4)
mm_fit3 = brm(mm_form3, data=df_all,prior=prior_mm_fit1, chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.9), cores = 4)
mm_fit4 = brm(mm_form4, data=df_all,prior=prior_mm_fit1, chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.9), cores = 4)
mm_fit5 = brm(mm_form5, data=df_all,prior=prior_mm_fit1, chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.9), cores = 4)
mm_fit6 = brm(mm_form6, data=df_all,prior=prior_mm_fit1, chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.9), cores = 4)
mm_fit7 = brm(mm_form7, data=df_all,prior=prior_mm_fit1, chains = 4,  warmup = 2000,
              iter = 4000,control = list(adapt_delta = 0.9), cores = 4)
mm_fit8 = brm(mm_form8, data=df_all,prior=prior_mm_fit1, chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.95), cores = 4)

####### Analyze the predictive power ############
best_height_pred_model <- loo_compare(waic(mod_height_base_pred), waic(mod_height_liana_pred), waic(mod_height_liana_inter_pred), waic(mod_height_sp_liana_pred,moment_match = TRUE),waic(mod_height_sp_liana_inter_pred),
            waic(mm_fit1), waic(mm_fit2),waic(mm_fit3),waic(mm_fit4),waic(mm_fit5),waic(mm_fit6), waic(mm_fit7),waic(mm_fit8))

