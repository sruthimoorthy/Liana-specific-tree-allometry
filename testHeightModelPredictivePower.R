library(ordinal)
library(lme4)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(brms)  
library(dplyr)
library(tidyr)

setwd("")

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

df_2019 <- df_2019 %>% filter(dbh >= 37)

speciesCodes <- read.csv("speciesCodeBook.csv")

height_data_2015 <- read.csv("new_height_liana_data_2015_no_outliers_q20.csv")


df_2015 <- height_data_2015 %>%
  left_join(speciesCodes, by = c("spcode"="Code"))

df_2015 <- data_frame(tag = df_2015$tag,
                      quadrant = df_2015$q20,
                      lianaLoad = df_2015$liana_merge,
                      species = df_2015$Full.Name,
                      dbh = df_2015$dbh_use/10,
                      height = df_2015$MaxHt,
                      source = "drone_photogrammetry",
                      year = 2015)

df_2015 <- df_2015 %>%
  mutate(lianaLoad = case_when(
    lianaLoad == 0 ~ 0,
    lianaLoad <= 2 ~ 1,
    lianaLoad > 2 ~ 2,
  ))

df_2015$lianaLoad <- factor(df_2015$lianaLoad, ordered=FALSE, levels=c("0","1","2"))
df_2015$species <- factor(df_2015$species)

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



mm_form1 <- brmsformula(log(height) ~ Vm + x * log(dbh) - log(K+(dbh^x)), Vm + x + K ~ 1 + (1|species), nl=TRUE)

mm_form2 <- brmsformula(log(height) ~ Vm + x * log(dbh) - log(K+(dbh^x)), Vm ~ 1 + lianaLoad + (1|species), x + K ~ 1 + (1|species), nl=TRUE)

mm_form3 <- brmsformula(log(height) ~ Vm + x * log(dbh) - log(K+(dbh^x)),  x ~ 1 + lianaLoad + (1|species), Vm + K ~ 1 + (1|species), nl=TRUE)

mm_form4 <- brmsformula(log(height) ~ Vm + x * log(dbh) - log(K+(dbh^x)),  K ~ 1 + lianaLoad + (1|species), Vm + x ~ 1 + (1|species), nl=TRUE)

mm_form5 <- brmsformula(log(height) ~ Vm + x * log(dbh) - log(K+(dbh^x)),  x + K ~ 1 + lianaLoad + (1|species), Vm  ~ 1 + (1|species), nl=TRUE)

mm_form6 <- brmsformula(log(height) ~ Vm + x * log(dbh) - log(K+(dbh^x)),  Vm + K ~ 1 + lianaLoad + (1|species), x  ~ 1 + (1|species), nl=TRUE)

mm_form7 <- brmsformula(log(height) ~ Vm + x * log(dbh) - log(K+(dbh^x)),  Vm + x ~ 1 + lianaLoad + (1|species), K  ~ 1 + (1|species), nl=TRUE)

mm_form8 <- brmsformula(log(height) ~ Vm + x * log(dbh) - log(K+(dbh^x)),  Vm + x + K ~ 1 + lianaLoad + (1|species), nl=TRUE)




mm_fit1 = brm(mm_form1, data=df_all, chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.9), cores = 4)
mm_fit2 = brm(mm_form2, data=df_all,chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.9), cores = 4)
mm_fit3 = brm(mm_form3, data=df_all, chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.9), cores = 4)
mm_fit4 = brm(mm_form4, data=df_all, chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.9), cores = 4)
mm_fit5 = brm(mm_form5, data=df_all, chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.9), cores = 4)
mm_fit6 = brm(mm_form6, data=df_all,chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.9), cores = 4)
mm_fit7 = brm(mm_form7, data=df_all, chains = 4,  warmup = 2000,
              iter = 4000,control = list(adapt_delta = 0.9), cores = 4)
mm_fit8 = brm(mm_form8, data=df_all,chains = 4,  warmup = 2000,
              iter = 4000, control = list(adapt_delta = 0.9), cores = 4)

####### Analyze the predictive power ############
best_height_pred_model <- loo_compare(waic(mod_height_base_pred), waic(mod_height_liana_pred), waic(mod_height_liana_inter_pred), waic(mod_height_sp_liana_pred,moment_match = TRUE),waic(mod_height_sp_liana_inter_pred),
            waic(mm_fit1), waic(mm_fit2),waic(mm_fit3),waic(mm_fit4),waic(mm_fit5),waic(mm_fit6), waic(mm_fit7),waic(mm_fit8))

############# Plotting the residuals of the log-log and michalis menten model ##############

####################### Height residual plot ##############################################

log_log_fitted_vals <- fitted(mod_height_liana_inter_pred,
                      newdata = df_all,
                      re_formula = NULL, 
                      allow_new_levels = TRUE,# includes random effects
                      scale = "response", # return predictions on height (original scale)
                      summary = TRUE)

# 'fitted_vals' has columns: "Estimate", "Est.Error", "Q2.5", "Q97.5", etc.,
# with one row per row of df_all. Let's merge them into df_all:
df_all$log_log_pred_height <- exp(log_log_fitted_vals[, "Estimate"])  # mean predicted height
df_all$log_log_pred_lower  <- exp(log_log_fitted_vals[, "Q2.5"])
df_all$log_log_pred_upper  <- exp(log_log_fitted_vals[, "Q97.5"])

# 2) Compute residuals => observed minus predicted
#    Make sure df_all$height is on the same scale as the model (i.e. actual height in meters).
df_all$log_log_resid <- df_all$height - df_all$log_log_pred_height

# 3) Plot residuals
# Option A: residuals vs. fitted
p_height_log_log_resid_fit <-ggplot(df_all, aes(x=log_log_pred_height, y=log_log_resid, color=lianaLoad)) +
  # 1) Smaller points
  geom_point(alpha=0.3, size=0.5) +  
  # 2) Horizontal reference line
  geom_hline(yintercept=0, linetype="dashed") +
  # 3) Smoother, grouped by lianaLoad
  geom_smooth(aes(group=lianaLoad), 
              method="gam",  # or method="gam", method="lm", etc.
              se=FALSE) +      # typically no ribbons for residual plots
  scale_color_manual(
    name = "Liana Load", 
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"),
    labels = c("0"="No Liana", "1"="Light", "2"="Heavy")  
  ) +
  labs(
    x="Predicted Height (m)",
    y="Residual",
    title = "log-log model"
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title = element_text(size=14),
    axis.text  = element_text(size=12)
  )



# Option B: residuals vs. DBH
p_log_log_resid_dbh <- ggplot(df_all, aes(x=dbh, y=log_log_resid, color=lianaLoad)) +
  # 1) Smaller points
  geom_point(alpha=0.3, size=0.5) +  
  # 2) Horizontal reference line
  geom_hline(yintercept=0, linetype="dashed") +
  # 3) Smoother, grouped by lianaLoad
  geom_smooth(aes(group=lianaLoad), 
              method="gam",  # or method="gam", method="lm", etc.
              se=FALSE) +      # typically no ribbons for residual plots
  scale_color_manual(
    name = "Liana Load", 
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"),
    labels = c("0"="No Liana", "1"="Light", "2"="Heavy")  
  ) +
  labs(
    x="DBH (cm)",
    y="Residual",
    title = "log-log model"
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title = element_text(size=14),
    axis.text  = element_text(size=12)
  )

gmm_fitted_vals <- fitted(mm_fit4,
                              newdata = df_all,
                              re_formula = NULL, 
                              allow_new_levels = TRUE,# includes random effects
                              scale = "response", # return predictions on height (original scale)
                              summary = TRUE)

# 'fitted_vals' has columns: "Estimate", "Est.Error", "Q2.5", "Q97.5", etc.,
# with one row per row of df_all. Let's merge them into df_all:
df_all$gmm_pred_height <- exp(gmm_fitted_vals[, "Estimate"])  # mean predicted height
df_all$gmm_pred_lower  <- exp(gmm_fitted_vals[, "Q2.5"])
df_all$gmm_pred_upper  <- exp(gmm_fitted_vals[, "Q97.5"])

# 2) Compute residuals => observed minus predicted
#    Make sure df_all$height is on the same scale as the model (i.e. actual height in meters).
df_all$gmm_resid <- df_all$height - df_all$gmm_pred_height

# 3) Plot residuals
# Option A: residuals vs. fitted
p_height_gmm_resid_fit <-ggplot(df_all, aes(x=gmm_pred_height, y=gmm_resid, color=lianaLoad)) +
  # 1) Smaller points
  geom_point(alpha=0.3, size=0.5) +  
  # 2) Horizontal reference line
  geom_hline(yintercept=0, linetype="dashed") +
  # 3) Smoother, grouped by lianaLoad
  geom_smooth(aes(group=lianaLoad), 
              method="gam",  # or method="gam", method="lm", etc.
              se=FALSE) +      # typically no ribbons for residual plots
  scale_color_manual(
    name = "Liana Load", 
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"),
    labels = c("0"="No Liana", "1"="Light", "2"="Heavy")  
  ) +
  labs(
    x="Predicted Height (m)",
    y="Residual",
    title = "gMM model"
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title = element_text(size=14),
    axis.text  = element_text(size=12)
  )



# Option B: residuals vs. DBH
p_gmm_resid_dbh <- ggplot(df_all, aes(x=dbh, y=gmm_resid, color=lianaLoad)) +
  # 1) Smaller points
  geom_point(alpha=0.3, size=0.5) +  
  # 2) Horizontal reference line
  geom_hline(yintercept=0, linetype="dashed") +
  # 3) Smoother, grouped by lianaLoad
  geom_smooth(aes(group=lianaLoad), 
              method="gam",  # or method="gam", method="lm", etc.
              se=FALSE) +      # typically no ribbons for residual plots
  scale_color_manual(
    name = "Liana Load", 
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"),
    labels = c("0"="No Liana", "1"="Light", "2"="Heavy")  
  ) +
  labs(
    x="DBH (cm)",
    y="Residual",
    title = "gMM model"
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title = element_text(size=14),
    axis.text  = element_text(size=12)
  )


resid_comp <- ggarrange(plotlist = c(list(p_log_log_resid_dbh), list(p_gmm_resid_dbh), list(p_height_log_log_resid_fit), list(p_height_gmm_resid_fit)),
                        nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom", labels = c("a","b","c","d"))

ggsave("resid_comp_best_pred_models.png", plot = resid_comp, height = 9, width = 12, dpi = 600)

