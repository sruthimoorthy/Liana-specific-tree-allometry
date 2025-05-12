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

df_2019 <- data_frame(tag = df_2019$tag,
                      lianaLoad = df_2019$lianaLoad,
                      species = df_2019$species,
                      dbh = df_2019$dbh/10,
                      height = df_2019$height,
                      cpa = df_2019$crownArea,
                      cDepth = df_2019$crownDepth,
                      cVol = df_2019$crownVol,
                      source = "tls",
                      subplot = df_2019$subplot,
                      year = 2019)

df_2019 <- df_2019 %>%
  mutate(lianaLoad = case_when(
    lianaLoad == 0 ~ 0,
    lianaLoad <= 50 ~ 1,
    lianaLoad > 50 ~ 2,
  ))

df_2019 <- df_2019 %>%
  mutate(combLianaLoad = ifelse(lianaLoad > 0, 1, 0))

df_2019$combLianaLoad <- as.factor(df_2019$combLianaLoad)

df_2019 <- df_2019 %>% filter(dbh >= 37)

df_2019$lianaLoad <- factor(df_2019$lianaLoad, ordered=FALSE, levels=c("0","1","2"))
df_2019$species <- factor(df_2019$species)
df_2019$subplot <- factor(df_2019$subplot)

domSpeciesBCI <- read.csv("speciesTraitsBCIData.csv")

mod_cpa_simple <- brm(
  formula = log(cpa) ~ log(dbh) + lianaLoad  + (1 |species) + (1|subplot),
  data = df_2019,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

summary(mod_cpa_simple)

ce_cpa_dbh_liana <- conditional_effects(
  mod_cpa_simple, 
  effects = "dbh:lianaLoad",   # or "lianaLoad:logDBH"
  re_formula = NA   
 
)

ce_cpa_dbh_liana <- plot(ce_cpa_dbh_liana, plot=FALSE)

df_plot <- ce_cpa_dbh_liana[[1]]$data

# Exponentiate the relevant columns if your model is log(...) ~ ...
df_plot$estimate__ <- exp(df_plot$estimate__)
df_plot$lower__    <- exp(df_plot$lower__)
df_plot$upper__    <- exp(df_plot$upper__)

# 5) Update the ggplot's data with the exponentiated values
ce_cpa_dbh_liana[[1]]$data <- df_plot


plot_ce_cpa_dbh_liana_custom <- ce_cpa_dbh_liana[[1]] +
  labs(
    x = "DBH (cm)",
    y = expression(paste("Crown Area " ,(m^{2})))
  ) +
  scale_color_manual(
    name = "Liana Load", 
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"), # or your own color palette
    labels = c("0"="No Liana", "1"="Light", "2"="Heavy")  
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

mod_cd_simple <- brm(
  formula = log(cDepth) ~ log(dbh) + lianaLoad  + (1 |species)+ (1|subplot) ,
  data = df_2019,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

mod_cd_no_subplot <- brm(
  formula = log(cDepth) ~ log(dbh) + lianaLoad  + (1 |species) ,
  data = df_2019,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

summary(mod_cd_simple)
summary(mod_cd_no_subplot)

ce_cDepth_dbh_liana <- conditional_effects(
  mod_cd_simple, 
  effects = "dbh:lianaLoad",   # or "lianaLoad:logDBH"
  re_formula = NA                 # ignore random effects for pop-level visualization
)

ce_cDepth_dbh_liana <- plot(ce_cDepth_dbh_liana, plot=FALSE)

df_plot <- ce_cDepth_dbh_liana[[1]]$data

# Exponentiate the relevant columns if your model is log(...) ~ ...
df_plot$estimate__ <- exp(df_plot$estimate__)
df_plot$lower__    <- exp(df_plot$lower__)
df_plot$upper__    <- exp(df_plot$upper__)

# 5) Update the ggplot's data with the exponentiated values
ce_cDepth_dbh_liana[[1]]$data <- df_plot

plot_ce_cDepth_dbh_liana_custom <- ce_cDepth_dbh_liana[[1]] +
  labs(
    x = "DBH (cm)",
    y = expression(paste("Crown Depth " ,(m)))
  ) +
  scale_color_manual(
    name = "Liana Load", 
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"), # or your own color palette
    labels = c("0"="No Liana", "1"="Light", "2"="Heavy")  
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


mod_cvol_simple <- brm(
  formula = log(cVol) ~ log(dbh) + lianaLoad  + (1 |species) +(1|subplot),
  data = df_2019,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

summary(mod_cvol_simple)


ce_cvol_dbh_liana <- conditional_effects(
  mod_cvol_simple, 
  effects = "dbh:lianaLoad",   # or "lianaLoad:logDBH"
  re_formula = NA                 # ignore random effects for pop-level visualization
)

mod_cvol_simple_nosubplot <- brm(
  formula = log(cVol) ~ log(dbh) + lianaLoad  + (1 |species) ,
  data = df_2019,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

summary(mod_cvol_simple_nosubplot)

ce_cvol_dbh_liana <- plot(ce_cvol_dbh_liana, plot=FALSE)

df_plot <- ce_cvol_dbh_liana[[1]]$data

# Exponentiate the relevant columns if your model is log(...) ~ ...
df_plot$estimate__ <- exp(df_plot$estimate__)
df_plot$lower__    <- exp(df_plot$lower__)
df_plot$upper__    <- exp(df_plot$upper__)

# 5) Update the ggplot's data with the exponentiated values
ce_cvol_dbh_liana[[1]]$data <- df_plot


plot_ce_cvol_dbh_liana_custom <- ce_cvol_dbh_liana[[1]] +
  labs(
    x = "DBH (cm)",
    y = expression(paste("Crown Volume " ,(m^3)))
  ) +
  scale_color_manual(
    name = "Liana Load", 
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"), # or your own color palette
    labels = c("0"="No Liana", "1"="Light", "2"="Heavy")  
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

###################### Crown metrics residual plot #########################


fitted_vals <- fitted(mod_cpa_simple,
                      newdata = df_2019,
                      re_formula = NULL, 
                      allow_new_levels = TRUE,# includes random effects
                      scale = "response", # return predictions on height (original scale)
                      summary = TRUE)

# 'fitted_vals' has columns: "Estimate", "Est.Error", "Q2.5", "Q97.5", etc.,
# with one row per row of df_all. Let's merge them into df_all:
df_2019$resid_pred_cpa <- fitted_vals[, "Estimate"]  # mean predicted height
df_2019$resid_cpa_lower  <- fitted_vals[, "Q2.5"]
df_2019$resid_cpa_upper  <- fitted_vals[, "Q97.5"]

# 2) Compute residuals => observed minus predicted
#    Make sure df_all$height is on the same scale as the model (i.e. actual height in meters).
df_2019$cpa_resid <- log(df_2019$cpa) - df_2019$resid_pred_cpa

# 3) Plot residuals
# Option A: residuals vs. fitted
p_cpa_resid_fit <- ggplot(df_2019, aes(x=resid_pred_cpa, y=cpa_resid, color = lianaLoad)) +
  geom_point(alpha=0.6) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_smooth(aes(group=lianaLoad), 
              method="gam",  # or method="gam", method="lm", etc.
              se=FALSE) +
  scale_color_manual(
    name = "Liana Load", 
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"), # or your own color palette
    labels = c("0"="No Liana", "1"="Light", "2"="Heavy")  
  ) +
  labs(
    x=expression(paste("Predicted crown projected area ",(m^2))),
    y=expression(paste("Residual = Observed - Predicted ",(m^2)))
  ) +
  theme_minimal()

fitted_vals <- fitted(mod_cd_simple,
                      newdata = df_2019,
                      re_formula = NULL, 
                      allow_new_levels = TRUE,# includes random effects
                      scale = "response", # return predictions on height (original scale)
                      summary = TRUE)

# 'fitted_vals' has columns: "Estimate", "Est.Error", "Q2.5", "Q97.5", etc.,
# with one row per row of df_all. Let's merge them into df_all:
df_2019$resid_pred_cd <- fitted_vals[, "Estimate"]  # mean predicted height
df_2019$resid_cd_lower  <- fitted_vals[, "Q2.5"]
df_2019$resid_cd_upper  <- fitted_vals[, "Q97.5"]

# 2) Compute residuals => observed minus predicted
#    Make sure df_all$height is on the same scale as the model (i.e. actual height in meters).
df_2019$cd_resid <- log(df_2019$cDepth) - df_2019$resid_pred_cd

# 3) Plot residuals
# Option A: residuals vs. fitted
p_cd_resid_fit <- ggplot(df_2019, aes(x=resid_pred_cd, y=cd_resid, color = lianaLoad)) +
  geom_point(alpha=0.6) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_smooth(aes(group=lianaLoad), 
              method="gam",  # or method="gam", method="lm", etc.
              se=FALSE) +
  scale_color_manual(
    name = "Liana Load", 
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"), # or your own color palette
    labels = c("0"="No Liana", "1"="Light", "2"="Heavy")  
  ) +
  labs(
    x=expression(paste("Predicted crown depth ",(m))),
    y=expression(paste("Residual = Observed - Predicted ",(m)))
  ) +
  theme_minimal()

fitted_vals <- fitted(mod_cvol_simple,
                      newdata = df_2019,
                      re_formula = NULL, 
                      allow_new_levels = TRUE,# includes random effects
                      scale = "response", # return predictions on height (original scale)
                      summary = TRUE)

# 'fitted_vals' has columns: "Estimate", "Est.Error", "Q2.5", "Q97.5", etc.,
# with one row per row of df_all. Let's merge them into df_all:
df_2019$resid_pred_cvol <- fitted_vals[, "Estimate"]  # mean predicted height
df_2019$resid_cvol_lower  <- fitted_vals[, "Q2.5"]
df_2019$resid_cvol_upper  <- fitted_vals[, "Q97.5"]

# 2) Compute residuals => observed minus predicted
#    Make sure df_all$height is on the same scale as the model (i.e. actual height in meters).
df_2019$cvol_resid <- log(df_2019$cVol) - df_2019$resid_pred_cvol

# 3) Plot residuals
# Option A: residuals vs. fitted
p_cvol_resid_fit <- ggplot(df_2019, aes(x=resid_pred_cvol, y=cvol_resid, color = lianaLoad)) +
  geom_point(alpha=0.6) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_smooth(aes(group=lianaLoad), 
              method="gam",  # or method="gam", method="lm", etc.
              se=FALSE) +
  scale_color_manual(
    name = "Liana Load", 
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"), # or your own color palette
    labels = c("0"="No Liana", "1"="Light", "2"="Heavy")  
  ) +
  labs(
    x=expression(paste("Predicted crown volume ",(m^3))),
    y=expression(paste("Residual = Observed - Predicted ",(m^3)))
  ) +
  theme_minimal()


crown_metrics_residual_plot_list <- c(list(p_cd_resid_fit), list(p_cpa_resid_fit), list(p_cvol_resid_fit))
crown_metrics_resid_plot <- ggarrange(plotlist = crown_metrics_residual_plot_list,
                                      nrow=1, ncol=3, common.legend = TRUE,
                                      legend = "bottom", labels = c("a","b","c"))

ggsave("log_log_crown_model_residual_plot.png", width =9, height = 5, plot = crown_metrics_resid_plot, dpi = 600)



###################### Plot predicted data ##################################

# We'll create a helper function that:
#  1) Gets the fitted values on the *linear predictor* scale (which is log(totVol)).
#  2) Exponentiates them to get totVol on the original scale.
#  3) Merges into df_2019.

get_predictions_for_log_model <- function(model, df, response_var = "totVol") {
  
  # 1) Use fitted() with summary=TRUE to get mean + intervals for each row
  #    re_formula=NULL includes random effects if your model used them
  preds_log <- fitted(model, 
                      newdata = df, 
                      re_formula = NA, 
                      summary = TRUE)
  
  # 2) We assume your model is log(response_var) on the LHS.
  #    So exponentiate:
  df$pred_est   <- exp(preds_log[,"Estimate"])
  df$pred_lower <- exp(preds_log[,"Q2.5"])
  df$pred_upper <- exp(preds_log[,"Q97.5"])
  
  # optionally rename them for clarity
  names(df)[names(df) %in% c("pred_est","pred_lower","pred_upper")] <-
    paste0(response_var, c("_pred","_lower","_upper"))
  
  return(df)
}

# For totVolModel:
df_2019 <- get_predictions_for_log_model(mod_cpa_simple, df_2019, "cpa")
df_2019 <- get_predictions_for_log_model(mod_cd_simple, df_2019, "cDepth")
df_2019 <- get_predictions_for_log_model(mod_cvol_simple, df_2019, "cVol")

df_2019 <- df_2019 %>%
  group_by(lianaLoad) %>%
  arrange(dbh) %>%
  ungroup()


p_cpa <- ggplot(df_2019, aes(x=dbh)) +
  geom_point(aes(y=cpa, color=lianaLoad), alpha=0.6) +
  geom_line(aes(y=cpa_pred, color=lianaLoad), size=1) +
  geom_ribbon(
    aes(ymin=cpa_lower, ymax=cpa_upper, fill=lianaLoad),
    alpha=0.2
  ) +
  scale_color_manual(
    name="Liana Load",
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"),
    labels = c("0"="No Liana","1"="Light","2"="Heavy")
  ) +
  scale_fill_manual(
    name="Liana Load",
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"),
    labels = c("0"="No Liana","1"="Light","2"="Heavy")
  ) +
  scale_x_log10(labels=label_number()) +  # x-axis log
  scale_y_log10(labels=label_number()) +  # y-axis log
  labs(
    x="DBH (cm)",
    y=expression(paste("Crown projected area ", (m^2)))
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title=element_text(size=14),
    axis.text=element_text(size=12),
    plot.title=element_text(size=14, face="bold")
  )

p_cDepth <- ggplot(df_2019, aes(x=dbh)) +
  geom_point(aes(y=cDepth, color=lianaLoad), alpha=0.6) +
  geom_line(aes(y=cDepth_pred, color=lianaLoad), size=1) +
  geom_ribbon(
    aes(ymin=cDepth_lower, ymax=cDepth_upper, fill=lianaLoad),
    alpha=0.2
  ) +
  scale_color_manual(
    name="Liana Load",
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"),
    labels = c("0"="No Liana","1"="Light","2"="Heavy")
  ) +
  scale_fill_manual(
    name="Liana Load",
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"),
    labels = c("0"="No Liana","1"="Light","2"="Heavy")
  ) +
  scale_x_log10(labels=label_number()) +  # x-axis log
  scale_y_log10(labels=label_number()) +  # y-axis log
  labs(
    x="DBH (cm)",
    y=expression(paste("Crown Depth", " (m)"))
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title=element_text(size=14),
    axis.text=element_text(size=12),
    plot.title=element_text(size=14, face="bold")
  )

p_cVol <- ggplot(df_2019, aes(x=dbh)) +
  geom_point(aes(y=cVol, color=lianaLoad), alpha=0.6) +
  geom_line(aes(y=cVol_pred, color=lianaLoad), size=1) +
  geom_ribbon(
    aes(ymin=cVol_lower, ymax=cVol_upper, fill=lianaLoad),
    alpha=0.2
  ) +
  scale_color_manual(
    name="Liana Load",
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"),
    labels = c("0"="No Liana","1"="Light","2"="Heavy")
  ) +
  scale_fill_manual(
    name="Liana Load",
    values = c("0"="#009e73","1"="#D55E00","2"="#000000"),
    labels = c("0"="No Liana","1"="Light","2"="Heavy")
  ) +
  scale_x_log10(labels=label_number()) +  # x-axis log
  scale_y_log10(labels=label_number()) +  # y-axis log
  labs(
    x="DBH (cm)",
    y = expression(paste("Crown volume ", (m^3)))
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title=element_text(size=14),
    axis.text=element_text(size=12),
    plot.title=element_text(size=14, face="bold")
  )

crown_plot_list <- c(list(p_cDepth), list(p_cpa), list(p_cVol))
crown_metrics_plot <- ggarrange(plotlist = crown_plot_list, 
                                     ncol=3, nrow=1, 
                                     common.legend=TRUE, legend="bottom",  labels = c("a", "b", "c"))

ggsave("crown_combined_plot.png", width = 11, height = 6, plot = crown_metrics_plot, dpi = 600)



#################### Random effects of species on crown metrics #####################

##### Area ######

re_sp_cpa <- ranef(mod_cpa_simple)$species
dimnames(re_sp_cpa)

re_sp_cpa_estimates <- data.frame(
  species = rownames(re_sp_cpa),
  est = re_sp_cpa[,"Estimate","Intercept"],
  lower = re_sp_cpa[,"Q2.5","Intercept"],
  upper = re_sp_cpa[,"Q97.5","Intercept"]
)

re_sp_cpa_estimates_dom <-  re_sp_cpa_estimates %>% filter(species %in% domSpeciesBCI$full.name)

plot_re_cpa <- ggplot(re_sp_cpa_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
       x="", y="Effect on \nlog(crown area)")

##### Volume ######

re_sp_cvol <- ranef(mod_cvol_simple)$species
dimnames(re_sp_cvol)

re_sp_cvol_estimates <- data.frame(
  species = rownames(re_sp_cvol),
  est = re_sp_cvol[,"Estimate","Intercept"],
  lower = re_sp_cvol[,"Q2.5","Intercept"],
  upper = re_sp_cvol[,"Q97.5","Intercept"]
)

re_sp_cvol_estimates_dom <-  re_sp_cvol_estimates %>% filter(species %in% domSpeciesBCI$full.name)

plot_re_cvol <- ggplot(re_sp_cvol_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
       x="", y="Effect on \nlog(crown volume)")

##### Depth ######

re_sp_cdep <- ranef(mod_cd_simple)$species
dimnames(re_sp_cdep)

re_sp_cdep_estimates <- data.frame(
  species = rownames(re_sp_cdep),
  est = re_sp_cdep[,"Estimate","Intercept"],
  lower = re_sp_cdep[,"Q2.5","Intercept"],
  upper = re_sp_cdep[,"Q97.5","Intercept"]
)

re_sp_cdep_estimates_dom <-  re_sp_cdep_estimates %>% filter(species %in% domSpeciesBCI$full.name)

plot_re_cdep <- ggplot(re_sp_cdep_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
       x="", y="Effect on \nlog(crown depth)")

#################### Random effects of subplots on crown metrics #######################


##### Area ######

re_subplot_cpa <- ranef(mod_cpa_simple)$subplot
dimnames(re_subplot_cpa)

re_subplot_cpa_estimates <- data.frame(
  subplot = rownames(re_subplot_cpa),
  est = re_subplot_cpa[,"Estimate","Intercept"],
  lower = re_subplot_cpa[,"Q2.5","Intercept"],
  upper = re_subplot_cpa[,"Q97.5","Intercept"]
)


ggplot(re_subplot_cpa_estimates, aes(x=reorder(subplot, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Subplot-Specific Random Intercepts",
       x="Subplots", y="Effect on \n log(crown projected area)")

##### Volume ######

re_subplot_cvol <- ranef(mod_cvol_simple)$subplot
dimnames(re_subplot_cvol)

re_subplot_cvol_estimates <- data.frame(
  subplot = rownames(re_subplot_cvol),
  est = re_subplot_cvol[,"Estimate","Intercept"],
  lower = re_subplot_cvol[,"Q2.5","Intercept"],
  upper = re_subplot_cvol[,"Q97.5","Intercept"]
)


ggplot(re_subplot_cvol_estimates, aes(x=reorder(subplot, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Subplot-Specific Random Intercepts",
       x="Subplots", y="Effect on \nlog(crown volume)")

##### Depth ######

re_subplot_cdep <- ranef(mod_cd_simple)$subplot
dimnames(re_subplot_cdep)

re_subplot_cdep_estimates <- data.frame(
  subplot = rownames(re_subplot_cdep),
  est = re_subplot_cdep[,"Estimate","Intercept"],
  lower = re_subplot_cdep[,"Q2.5","Intercept"],
  upper = re_subplot_cdep[,"Q97.5","Intercept"]
)


ggplot(re_subplot_cdep_estimates, aes(x=reorder(subplot, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Subplot-Specific Random Intercepts",
       x="Subplots", y="Effect on \nlog(crown depth)")

#################### Subplot estimate and liana load relation #####################

df_fraction <- df_2019 %>%
  group_by(subplot) %>%
  summarise(
    fraction_liana12 = mean(lianaLoad %in% c("1","2"))  # if lianaLoad is factor
    # or if lianaLoad is numeric, do: mean(lianaLoad > 0)
  )

################ Volume ###########ß


df_sub_re_cvol_join <- left_join(re_subplot_cvol_estimates, df_fraction, by="subplot")

ggplot(df_sub_re_cvol_join, aes(x=fraction_liana12, y=est)) +
  geom_point() +
  geom_smooth(method="lm") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "pearson",  # or "spearman"
           label.x.npc = "left", # position the label
           label.y.npc = "top") +
  labs(x="Fraction of Trees with Lianas", y="Subplot Random Intercept") +
  theme_minimal()

############## Area ##############


df_sub_re_cvol_join <- left_join(re_subplot_cvol_estimates, df_fraction, by="subplot")

subplot_cvol_re_effect <- ggplot(df_sub_re_cvol_join, aes(x=fraction_liana12, y=est)) +
  geom_point() +
  geom_smooth(method="lm") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "pearson",  # or "spearman"
           label.x.npc = "left", # position the label
           label.y.npc = "top") +
  labs(x="Fraction of Trees with Lianas", y="Subplot Random Intercept", title = expression(paste("Crown volume ", (m^3) ))) +
  theme_minimal()
ggsave("subplot_cvol_re_effect.png", width = 6, height = 5, plot = subplot_cvol_re_effect, dpi = 600)
