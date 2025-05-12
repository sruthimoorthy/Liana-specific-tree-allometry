library(ordinal)
library(lme4)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(brms)  
library(dplyr)
library(tidyr)
library(ggpubr)
library(patchwork)
library(cowplot)
library(scales)

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

# (3) Fit the model:
totVolModel <- brm(
  formula = bf(
    log(totVol) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

totVolModelNoSubplot <- brm(
  formula = bf(
    log(totVol) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) 
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

trunkVolModel <- brm(
  formula = bf(
    log(trunkVol) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

brVolModel <- brm(
  formula = bf(
    log(branchVol12) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)




ce_totVol_model <- conditional_effects(
  totVolModel, 
  effects = "dbh:lianaLoad",   # or "lianaLoad:logDBH"
  re_formula = NA   
  
)

ce_trunkVol_model <- conditional_effects(
  trunkVolModel, 
  effects = "dbh:lianaLoad",   # or "lianaLoad:logDBH"
  re_formula = NA   
  
)

ce_brVol_model <- conditional_effects(
  brVolModel, 
  effects = "dbh:lianaLoad",   # or "lianaLoad:logDBH"
  re_formula = NA   
  
)

####### Total volume conditional effects plot ###########

ce_totVol_model <- plot(ce_totVol_model, plot=FALSE)

df_totVol_plot <- ce_totVol_model[[1]]$data

# Exponentiate the relevant columns if your model is log(...) ~ ...
df_totVol_plot$estimate__ <- exp(df_totVol_plot$estimate__)/1000
df_totVol_plot$lower__    <- exp(df_totVol_plot$lower__)/1000
df_totVol_plot$upper__    <- exp(df_totVol_plot$upper__)/1000

# 5) Update the ggplot's data with the exponentiated values
ce_totVol_model[[1]]$data <- df_totVol_plot


plot_ce_totVol_model <- ce_totVol_model[[1]] +
  geom_point(
    data = df_2019,
    aes(x = dbh, y = totVol/1000, color = factor(lianaLoad)),
    inherit.aes = FALSE,    # don't inherit the conditional_effects aes
    alpha = 0.5,            # slightly transparent points
    size = 2
  )+
  labs(
    x = "DBH (cm)",
    y = expression(paste("Total Volume"["0,1,2"], ' ' ,(m^{3})))
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

############# Trunk Volume conditional effects plot ################

ce_trunkVol_model <- plot(ce_trunkVol_model, plot=FALSE)

df_trunkVol_plot <- ce_trunkVol_model[[1]]$data

# Exponentiate the relevant columns if your model is log(...) ~ ...
df_trunkVol_plot$estimate__ <- exp(df_trunkVol_plot$estimate__)/1000
df_trunkVol_plot$lower__    <- exp(df_trunkVol_plot$lower__)/1000
df_trunkVol_plot$upper__    <- exp(df_trunkVol_plot$upper__)/1000

# 5) Update the ggplot's data with the exponentiated values
ce_trunkVol_model[[1]]$data <- df_trunkVol_plot


plot_ce_trunkVol_model <- ce_trunkVol_model[[1]] +
  labs(
    x = "DBH (cm)",
    y = expression(paste("Trunk Volume"["0"], ' ' ,(m^{3})))
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

############# Branch Volume conditional effects plot ################

ce_brVol_model <- plot(ce_brVol_model, plot=FALSE)

df_brVol_plot <- ce_brVol_model[[1]]$data

# Exponentiate the relevant columns if your model is log(...) ~ ...
df_brVol_plot$estimate__ <- exp(df_brVol_plot$estimate__)/1000
df_brVol_plot$lower__    <- exp(df_brVol_plot$lower__)/1000
df_brVol_plot$upper__    <- exp(df_brVol_plot$upper__)/1000

# 5) Update the ggplot's data with the exponentiated values
ce_brVol_model[[1]]$data <- df_brVol_plot


plot_ce_brVol_model <- ce_brVol_model[[1]] +
  geom_point(
    data = df_2019,
    aes(x = dbh, y = branchVol, color = factor(lianaLoad)),
    inherit.aes = FALSE,    # don't inherit the conditional_effects aes
    alpha = 0.5,            # slightly transparent points
    size = 2
  ) +
  labs(
    x = "DBH (cm)",
    y = expression(paste("Branch Volume"["1,2"], ' ' ,(m^{3})))
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


qsm_volume_plot_list <- c(list(plot_ce_totVol_model), list(plot_ce_trunkVol_model), list(plot_ce_brVol_model))
qsm_volume_metrics_plot <- ggarrange(plotlist = qsm_volume_plot_list, 
                                     ncol=3, nrow=1, 
                                     common.legend=TRUE, legend="bottom",  labels = c("a", "b", "c"))

ggsave("qsm_volume_combined_plot.png", width = 11, height = 6, plot = qsm_volume_metrics_plot, dpi = 600)

###################### Volume metrics residual plot #########################
df_2019_pred <- read.csv("bci_trees_tls_final_qsmv2_4_metrics.csv")

df_2019_pred <- df_2019_pred %>% filter(dbh >= 37)
df_2019_pred <- df_2019_pred %>%
  mutate(lianaLoad = case_when(
    lianaLoad == 0 ~ 0,
    lianaLoad <= 50 ~ 1,
    lianaLoad > 50 ~ 2,
  ))

df_2019_pred$lianaLoad <- factor(df_2019_pred$lianaLoad, ordered=FALSE, levels=c("0","1","2"))
df_2019_pred$species <- factor(df_2019_pred$species)
df_2019_pred$subplot <- factor(df_2019_pred$subplot)

df_2019_pred <- data.frame(dbh = df_2019_pred$dbh,
                           lianaLoad = df_2019_pred$lianaLoad,
                           species = df_2019_pred$species,
                           totVol = df_2019_pred$totVol,
                           trunkVol = df_2019_pred$trunkVol,
                           branchVol12 = df_2019_pred$branchVol12)

fitted_vals <- fitted(totVolModel,
                      newdata = df_2019_pred,
                      re_formula = NULL, 
                      allow_new_levels = TRUE,# includes random effects
                      scale = "response", # return predictions on height (original scale)
                      summary = TRUE)

# 'fitted_vals' has columns: "Estimate", "Est.Error", "Q2.5", "Q97.5", etc.,
# with one row per row of df_all. Let's merge them into df_all:
df_2019_pred$resid_pred_totvol <- fitted_vals[, "Estimate"]  # mean predicted height
df_2019_pred$resid_totvol_lower  <- fitted_vals[, "Q2.5"]
df_2019_pred$resid_totvol_upper  <- fitted_vals[, "Q97.5"]

# 2) Compute residuals => observed minus predicted
#    Make sure df_all$height is on the same scale as the model (i.e. actual height in meters).
df_2019_pred$totvol_resid <- log(df_2019_pred$totVol) - df_2019_pred$resid_pred_totvol

# 3) Plot residuals
# Option A: residuals vs. fitted
p_totvol_resid_fit <- ggplot(df_2019_pred, aes(x=resid_pred_totvol, y=totvol_resid, color = lianaLoad)) +
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
    x=expression(paste("Predicted total volume ",(L))),
    y=expression(paste("Residual = Observed - Predicted ",(L)))
  ) +
  theme_minimal()

fitted_vals <- fitted(trunkVolModel,
                      newdata = df_2019_pred,
                      re_formula = NULL, 
                      allow_new_levels = TRUE,# includes random effects
                      scale = "response", # return predictions on height (original scale)
                      summary = TRUE)

# 'fitted_vals' has columns: "Estimate", "Est.Error", "Q2.5", "Q97.5", etc.,
# with one row per row of df_all. Let's merge them into df_all:
df_2019_pred$resid_pred_trunkvol <- fitted_vals[, "Estimate"]  # mean predicted height
df_2019_pred$resid_trunkvol_lower  <- fitted_vals[, "Q2.5"]
df_2019_pred$resid_trunkvol_upper  <- fitted_vals[, "Q97.5"]

# 2) Compute residuals => observed minus predicted
#    Make sure df_all$height is on the same scale as the model (i.e. actual height in meters).
df_2019_pred$trunkvol_resid <- log(df_2019_pred$trunkVol) - df_2019_pred$resid_pred_trunkvol

# 3) Plot residuals
# Option A: residuals vs. fitted
p_trunkvol_resid_fit <- ggplot(df_2019_pred, aes(x=resid_pred_trunkvol, y=trunkvol_resid, color = lianaLoad)) +
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
    x=expression(paste("Predicted trunk volume ",(L))),
    y=expression(paste("Residual = Observed - Predicted ",(L)))
  ) +
  theme_minimal()

fitted_vals <- fitted(brVolModel,
                      newdata = df_2019_pred,
                      re_formula = NULL, 
                      allow_new_levels = TRUE,# includes random effects
                      scale = "response", # return predictions on height (original scale)
                      summary = TRUE)

# 'fitted_vals' has columns: "Estimate", "Est.Error", "Q2.5", "Q97.5", etc.,
# with one row per row of df_all. Let's merge them into df_all:
df_2019_pred$resid_pred_brvol <- fitted_vals[, "Estimate"]  # mean predicted height
df_2019_pred$resid_brvol_lower  <- fitted_vals[, "Q2.5"]
df_2019_pred$resid_brvol_upper  <- fitted_vals[, "Q97.5"]

# 2) Compute residuals => observed minus predicted
#    Make sure df_all$height is on the same scale as the model (i.e. actual height in meters).
df_2019_pred$brvol_resid <- log(df_2019_pred$branchVol12) - df_2019_pred$resid_pred_brvol

# 3) Plot residuals
# Option A: residuals vs. fitted
p_brvol_resid_fit <- ggplot(df_2019_pred, aes(x=resid_pred_brvol, y=brvol_resid, color = lianaLoad)) +
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
    x=expression(paste("Predicted branch volume ",(L))),
    y=expression(paste("Residual = Observed - Predicted ",(L)))
  ) +
  theme_minimal()


vol_metrics_residual_plot_list <- c(list(p_totvol_resid_fit), list(p_trunkvol_resid_fit), list(p_brvol_resid_fit))
vol_metrics_resid_plot <- ggarrange(plotlist = vol_metrics_residual_plot_list,
                                      nrow=1, ncol=3, common.legend = TRUE,
                                      legend = "bottom", labels = c("a","b","c"))

ggsave("log_log_volume_model_residual_plot.png", width =9, height = 5, plot = vol_metrics_resid_plot, dpi = 600)



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
df_2019_pred <- get_predictions_for_log_model(totVolModel, df_2019_pred, "totVol")
df_2019_pred <- get_predictions_for_log_model(trunkVolModel, df_2019_pred, "trunkVol")
df_2019_pred <- get_predictions_for_log_model(brVolModel, df_2019_pred, "branchVol12")

df_2019_pred <- df_2019_pred %>%
  group_by(lianaLoad) %>%
  arrange(dbh) %>%
  ungroup()


p_totVol <- ggplot(df_2019_pred, aes(x=dbh)) +
  geom_point(aes(y=totVol, color=lianaLoad), alpha=0.6) +
  geom_line(aes(y=totVol_pred, color=lianaLoad), size=1) +
  geom_ribbon(
    aes(ymin=totVol_lower, ymax=totVol_upper, fill=lianaLoad),
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
    y=expression(paste("Total Volume"["0,1,2"], " (L)"))
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title=element_text(size=14),
    axis.text=element_text(size=12),
    plot.title=element_text(size=14, face="bold")
  )

p_trunkVol <- ggplot(df_2019_pred, aes(x=dbh)) +
  geom_point(aes(y=trunkVol, color=lianaLoad), alpha=0.6) +
  geom_line(aes(y=trunkVol_pred, color=lianaLoad), size=1) +
  geom_ribbon(
    aes(ymin=trunkVol_lower, ymax=trunkVol_upper, fill=lianaLoad),
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
    y=expression(paste("Trunk Volume"["0"], " (L)"))
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title=element_text(size=14),
    axis.text=element_text(size=12),
    plot.title=element_text(size=14, face="bold")
  )

p_brVol <- ggplot(df_2019_pred, aes(x=dbh)) +
  geom_point(aes(y=branchVol12, color=lianaLoad), alpha=0.6) +
  geom_line(aes(y=branchVol12_pred, color=lianaLoad), size=1) +
  geom_ribbon(
    aes(ymin=branchVol12_lower, ymax=branchVol12_upper, fill=lianaLoad),
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
    y = expression(paste("Branch Volume"["1,2"], " (L)"))
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title=element_text(size=14),
    axis.text=element_text(size=12),
    plot.title=element_text(size=14, face="bold")
  )

qsm_volume_plot_list <- c(list(p_totVol), list(p_trunkVol), list(p_brVol))

top_vol_plot <- ggarrange(p_totVol, 
                             ncol=1, nrow=1, 
                             common.legend=TRUE, legend="none",  labels = c( "a"))

bottom_vol_plot <- ggarrange(p_trunkVol,p_brVol, 
          ncol=2, nrow=1, 
          common.legend=TRUE, legend="bottom",  labels = c( "b", "c"))


qsm_volume_metrics_plot <- ggarrange(top_vol_plot,bottom_vol_plot, 
                                     ncol=1, nrow=2, 
                                     common.legend=TRUE, legend="bottom")

ggsave("qsm_volume_combined_plot.png", width = 11, height = 9, plot = qsm_volume_metrics_plot, dpi = 600)


#################### Random effects of species on Volume #####################

domSpeciesBCI <- read.csv("speciesTraitsBCIData.csv")

######### Total Volume #########

re_sp_vol <- ranef(totVolModel)$species
dimnames(re_sp_vol)

re_sp_vol_estimates <- data.frame(
  species = rownames(re_sp_vol),
  est = re_sp_vol[,"Estimate","Intercept"],
  lower = re_sp_vol[,"Q2.5","Intercept"],
  upper = re_sp_vol[,"Q97.5","Intercept"]
)

re_sp_vol_estimates_dom <-  re_sp_vol_estimates %>% filter(species %in% domSpeciesBCI$full.name)

totvol_re_plot <- ggplot(re_sp_vol_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    x="", y="Effect on \n log(total volume)")

######### Trunk Volume #########

re_sp_vol <- ranef(trunkVolModel)$species
dimnames(re_sp_vol)

re_sp_vol_estimates <- data.frame(
  species = rownames(re_sp_vol),
  est = re_sp_vol[,"Estimate","Intercept"],
  lower = re_sp_vol[,"Q2.5","Intercept"],
  upper = re_sp_vol[,"Q97.5","Intercept"]
)

re_sp_vol_estimates_dom <-  re_sp_vol_estimates %>% filter(species %in% domSpeciesBCI$full.name)

trunkvol_re_plot <- ggplot(re_sp_vol_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    x="", y="Effect on \n log(trunk volume)")

######### Branch Volume #########

re_sp_vol <- ranef(brVolModel)$species
dimnames(re_sp_vol)

re_sp_vol_estimates <- data.frame(
  species = rownames(re_sp_vol),
  est = re_sp_vol[,"Estimate","Intercept"],
  lower = re_sp_vol[,"Q2.5","Intercept"],
  upper = re_sp_vol[,"Q97.5","Intercept"]
)

re_sp_vol_estimates_dom <-  re_sp_vol_estimates %>% filter(species %in% domSpeciesBCI$full.name)

brvol_re_plot <- ggplot(re_sp_vol_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    x="", y="Effect on \n log(branch volume)")

############## Species random effects conbined figure ##################

qsm_random_plots <- ggarrange(totvol_re_plot, trunkvol_re_plot, brvol_re_plot, 
                              ncol=3, nrow=1, 
                              labels = c("a", "b", "c"))

ggsave("qsm_volume_re_plots_combined.png", width = 12, height = 8, plot = qsm_random_plots, dpi = 600)

#################### Random effects of subplots on volume #######################


re_subplot_vol <- ranef(totVolModel)$subplot
dimnames(re_subplot_vol)

re_subplot_vol_estimates <- data.frame(
  subplot = rownames(re_subplot_vol),
  est = re_subplot_vol[,"Estimate","Intercept"],
  lower = re_subplot_vol[,"Q2.5","Intercept"],
  upper = re_subplot_vol[,"Q97.5","Intercept"]
)


ggplot(re_subplot_vol_estimates, aes(x=reorder(subplot, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Subplot-Specific Random Intercepts",
       x="Subplots", y="Effect on \n log(total volume)")


re_subplot_brvol <- ranef(brVolModel)$subplot
dimnames(re_subplot_brvol)

re_subplot_brvol_estimates <- data.frame(
  subplot = rownames(re_subplot_brvol),
  est = re_subplot_brvol[,"Estimate","Intercept"],
  lower = re_subplot_brvol[,"Q2.5","Intercept"],
  upper = re_subplot_brvol[,"Q97.5","Intercept"]
)

#################### Subplot estimate and liana load relation #####################

df_fraction <- df_2019 %>%
  group_by(subplot) %>%
  summarise(
    fraction_liana12 = mean(lianaLoad %in% c("1","2"))  # if lianaLoad is factor
    # or if lianaLoad is numeric, do: mean(lianaLoad > 0)
  )

################ Volume ###########ß


df_sub_re_vol_join <- left_join(re_subplot_vol_estimates, df_fraction, by="subplot")

subplot_totvol_re_effect <- ggplot(df_sub_re_vol_join, aes(x=fraction_liana12, y=est)) +
  geom_point() +
  geom_smooth(method="lm") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "pearson",  # or "spearman"
           label.x.npc = "left", # position the label
           label.y.npc = "top") +
  labs(x="Fraction of Trees with Lianas", y="Subplot Random Intercept", title = expression(paste("Total Volume"["0,1,2"]," ", (m^3) ))) +
  theme_minimal()

ggsave("subplot_totvol_re_effect.png", width = 6, height = 5, plot = subplot_totvol_re_effect, dpi = 600)

################Branch Volume ###########ß


df_sub_re_brvol_join <- left_join(re_subplot_brvol_estimates, df_fraction, by="subplot")

ggplot(df_sub_re_brvol_join, aes(x=fraction_liana12, y=est)) +
  geom_point() +
  geom_smooth(method="lm") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "pearson",  # or "spearman"
           label.x.npc = "left", # position the label
           label.y.npc = "top") +
  labs(x="Fraction of Trees with Lianas", y="Subplot Random Intercept") +
  theme_minimal()


