library(ordinal)
library(lme4)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(brms)  
library(dplyr)
library(tidyr)
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

df_2019$totAreaVolRatio <- df_2019$totArea/(df_2019$totVol/1000)
df_2019$trunkAreaVolRatio <- df_2019$trunkArea/(df_2019$trunkVol/1000)
df_2019$brAreaVolRatio <- df_2019$branchArea12/(df_2019$branchVol12/1000)

# (3) Fit the model:
mod_totSArea <- brm(
  formula = bf(
    log(totArea) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

mod_trunkSArea <- brm(
  formula = bf(
    log(trunkArea) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

# (3) Fit the model:
mod_branchSArea <- brm(
  formula = bf(
    log(branchArea12) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

# (3) Fit the model:
mod_totAreaVolRatio <- brm(
  formula = bf(
    log(totAreaVolRatio) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

mod_trunkAreaVolRatio <- brm(
  formula = bf(
    log(trunkAreaVolRatio) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

# (3) Fit the model:
mod_brAreaVolRatio <- brm(
  formula = bf(
    log(brAreaVolRatio) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

mod_brLen <- brm(
  formula = bf(
    log(branchLen12) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

# (3) Fit the model:
mod_brNo <- brm(
  formula = bf(
    log(branchNo12) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_2019,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)


# 1) Create a named list of your models
brmetrics_model_list <- list(
  totArea        = mod_totSArea,
  trunkArea      = mod_trunkSArea,
  branchArea12      = mod_branchSArea,
  totAreaVolRatio   = mod_totAreaVolRatio,
  trunkAreaVolRatio = mod_trunkAreaVolRatio,
  brAreaVolRatio  = mod_brAreaVolRatio,
  branchLen12 = mod_brLen,
  branchNo12 = mod_brNo 
)

#################### Plot fitted models with data points ####################
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
                           totArea = df_2019_pred$totArea,
                           trunkArea = df_2019_pred$trunkArea,
                           branchArea12 = df_2019_pred$branchArea12,
                           totAreaVolRatio = df_2019_pred$totArea/(df_2019_pred$totVol/1000),
                           trunkAreaVolRatio = df_2019_pred$trunkArea/(df_2019_pred$trunkVol/1000),
                           brAreaVolRatio = df_2019_pred$branchArea12/(df_2019_pred$branchVol12/1000),
                           branchLen12 = df_2019_pred$branchLen12,
                           branchNo12 = df_2019_pred$branchNo12
                           )

# This function will:
#   1) For each model in sarea_model_list, do partial-effect predictions on df
#      (assuming log-transformed response).
#   2) Merge predictions (Estimate + intervals) into df.
#   3) Plot the observed vs. predicted lines + ribbons, 
#      returning a list of ggplot objects.


plot_br_metrics <- function(model_list, df, 
                               color_map = c("0"="#009e73","1"="#D55E00","2"="#000000")) {
  
  # Y-axis labels: customize for each metric's meaning
  metric_labels <- c(
    totArea          = expression(paste("Total S. area "["0,1,2"], (m^2))),
    trunkArea        = expression(paste("Trunk S. area "["0"], (m^2))),
    branchArea12     = expression(paste("Branch S. area "["1,2"], (m^2))),
    totAreaVolRatio  = expression(paste("Total S. area / volume "["0,1,2"], (m^-1))),
    trunkAreaVolRatio= expression(paste("Trunk S. area / volume "["0"], (m^-1))),
    brAreaVolRatio   = expression(paste("Branch S. area / volume "["1,2"], (m^-1))),
    branchLen12   = expression(paste("Branch length "["1,2"], (m))),
    branchNo12   = expression(paste("No. of branches "["1,2"], (n)))
  )
  
  # Decide which metrics use log scale on y
  # e.g., surface-area metrics => TRUE, ratio metrics => FALSE
  use_log_scale <- c(
    totArea = TRUE,
    trunkArea = TRUE,
    branchArea12 = TRUE,
    totAreaVolRatio = TRUE,
    trunkAreaVolRatio = TRUE,
    brAreaVolRatio = TRUE,
    branchLen12 = TRUE,
    branchNo12 = TRUE
  )
  
  # Helper function for partial-effect predictions on log link
  get_predictions_for_log_model <- function(model, df, response_var) {
    preds <- fitted(model,
                    newdata = df,
                    re_formula = NA,    # ignore random effects
                    summary = TRUE,
                    scale = "response") # exponentiate if log link
    df[[paste0(response_var, "_pred")]]   <- exp(preds[,"Estimate"])
    df[[paste0(response_var, "_lower")]]  <- exp(preds[,"Q2.5"])
    df[[paste0(response_var, "_upper")]]  <- exp(preds[,"Q97.5"])
    return(df)
  }
  
  plot_list <- list()
  
  for (var_name in names(model_list)) {
    cat("Generating plot for model:", var_name, "\n")
    mod <- model_list[[var_name]]
    
    # 1) Merge predictions
    df2 <- get_predictions_for_log_model(mod, df, var_name)
    
    # 2) Sort by DBH & lianaLoad for a clean line
    df2 <- df2 %>%
      group_by(lianaLoad) %>%
      arrange(dbh) %>%
      ungroup()
    
    # 3) Build the base plot
    p <- ggplot(df2, aes(x=dbh)) +
      # Observed points => the actual variable = var_name in df2
      geom_point(aes_string(y=var_name, color="lianaLoad"), alpha=0.6) +
      # Fitted line => var_name_pred
      geom_line(aes_string(y=paste0(var_name, "_pred"), color="lianaLoad"), size=1) +
      # Ribbon => var_name_lower, var_name_upper
      geom_ribbon(aes_string(
        ymin=paste0(var_name, "_lower"),
        ymax=paste0(var_name, "_upper"),
        fill="lianaLoad"
      ), alpha=0.2) +
      # color/fill scale
      scale_color_manual(
        name="Liana Load",
        values=color_map,
        labels=c("0"="No Liana","1"="Light","2"="Heavy")
      ) +
      scale_fill_manual(
        name="Liana Load",
        values=color_map,
        labels=c("0"="No Liana","1"="Light","2"="Heavy")
      ) +
      scale_x_log10(labels=label_number()) +  # x-axis in log scale, numeric labels
      labs(
        x="DBH (cm)",
        y=metric_labels[var_name],
  
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.title=element_text(size=14),
        axis.text=element_text(size=12)
      )
    
    # 4) Conditionally apply y log scale if use_log_scale[var_name]==TRUE
    if (use_log_scale[var_name]) {
      p <- p + scale_y_log10(labels = label_number())
    }
    
    plot_list[[var_name]] <- p
  }
  
  return(plot_list)
}

# Usage:
# Suppose your data is in df_2019_pred, 
# and your 6 brms models are in sarea_model_list:

brmetrics_plots <- plot_br_metrics(brmetrics_model_list, df_2019_pred)
# Then display them:
# sarea_plots[["totArea"]]
# sarea_plots[["branchArea12"]]
# etc.

br_metrics_final_plot <- ggarrange(plotlist = brmetrics_plots, 
                                   ncol=3, nrow=3, 
                                   common.legend=TRUE, legend="bottom", labels = c("a","b","c","d","e","f", "g","h"))

# 5) Display
print(br_metrics_final_plot)
ggsave("qsm_brmetrics_combined_plot.png", width = 12, height = 12, plot = br_metrics_final_plot, dpi = 800)

############################ Branch metrics residual plot #################################

plot_residuals_brmetrics <- function(model_list, df) {
  metric_labels <- c(
    totArea          = expression(paste("Predicted total S. area "["0,1,2"], (m^2))),
    trunkArea        = expression(paste("Predicted trunk S. area "["0"], (m^2))),
    branchArea12     = expression(paste("Predicted branch S. area "["1,2"], (m^2))),
    totAreaVolRatio  = expression(paste("Predicted total S. area / volume "["0,1,2"], (m^-1))),
    trunkAreaVolRatio= expression(paste("Predicted trunk S. area / volume "["0"], (m^-1))),
    brAreaVolRatio   = expression(paste("Predicted branch S. area / volume "["1,2"], (m^-1))),
    branchLen12   = expression(paste("Predicted branch length "["1,2"], (m))),
    branchNo12   = expression(paste("Predicted no. of branches "["1,2"], (n)))
  )
  
  # We'll store each residual plot in a list:
  plot_list <- list()
  
  for (var_name in names(model_list)) {
    cat("Building residual plot for:", var_name, "\n")
    model <- model_list[[var_name]]
    
    # 1) Get predictions on *original* scale (exponentiated if LHS is log(...))
    preds <- fitted(model,
                    newdata          = df,
                    re_formula       = NULL,    # include random effects
                    scale            = "response",
                    allow_new_levels = TRUE,    # if new species/subplots
                    summary          = TRUE)
    
    # preds has columns: Estimate, Est.Error, Q2.5, Q97.5, ...
    # The 'Estimate' is the predicted value on the original scale (e.g. m^2).
    
    # 2) Merge predicted => df
    df2 <- df %>% mutate(
      pred_est   = preds[,"Estimate"],
      pred_lower = preds[,"Q2.5"],
      pred_upper = preds[,"Q97.5"]
    )
    
    # 3) Compute residual: log(observed) - log(predicted).
    #    But first check if df2[[var_name]] is strictly positive for log(...).
    #    If zero or negative, you'll need a different approach.
    
    # We'll store the residual in e.g. 'resid'
    df2 <- df2 %>%
      mutate(
        # log(observed). E.g. if var_name="totArea", then log(totArea).
        log_obs = log(.data[[var_name]]),
        # log(predicted):
        log_pred = pred_est,
        # residual = log(observed) - log(pred)
        resid = log_obs - log_pred
      )
    
    # 4) Build the residual plot => residual vs. predicted (on original scale or log scale)
    #    If you want the x-axis to be 'pred_est' in linear scale, do that. 
    #    Or if you prefer log scale, do log10 transformation.
    
    y_label <- "Residual"
    x_label <-  metric_labels[var_name]
    
    p_resid <- ggplot(df2, aes(x=pred_est, y=resid, color=lianaLoad)) +
      geom_point(alpha=0.6) +
      geom_hline(yintercept=0, linetype="dashed") +
      geom_smooth(aes(group=lianaLoad), 
                  method="gam",  # or method="gam", method="lm", etc.
                  se=FALSE) +
      scale_color_manual(
        name="Liana Load", 
        values=c("0"="#009e73","1"="#D55E00","2"="#000000"),
        labels=c("0"="No Liana","1"="Light","2"="Heavy")
      ) +
      labs(
        x=x_label,
        y=y_label
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.title=element_text(size=14),
        axis.text=element_text(size=12)
      )
    
    # Optionally log scale the x-axis if it's big:
    # p_resid <- p_resid + scale_x_log10(labels=label_number())
    
    # store in list
    plot_list[[var_name]] <- p_resid
  }
  
  return(plot_list)
}

brmetrics_residual_plots <- plot_residuals_brmetrics(brmetrics_model_list, df_2019_pred)

brmetrics_metrics_resid_plot <- ggarrange(plotlist = brmetrics_residual_plots,
                                    nrow=3, ncol=3,
                                    common.legend=TRUE, legend="bottom", labels = c("a","b","c","d","e","f", "g","h"))

ggsave("log_log_brmetrics_model_residual_plot.png", width =12, height = 12, plot = brmetrics_metrics_resid_plot, dpi = 600)



#################### Random effects of species on Volume #####################

domSpeciesBCI <- read.csv("speciesTraitsBCIData.csv")

############# Surface Area #################

re_sp_brarea <- ranef(mod_totSArea)$species
dimnames(re_sp_brarea)

re_sp_brarea_estimates <- data.frame(
  species = rownames(re_sp_brarea),
  est = re_sp_brarea[,"Estimate","Intercept"],
  lower = re_sp_brarea[,"Q2.5","Intercept"],
  upper = re_sp_brarea[,"Q97.5","Intercept"]
)

re_sp_brarea_estimates_dom <-  re_sp_brarea_estimates %>% filter(species %in% domSpeciesBCI$full.name)

re_sp_totsarea_plot <- ggplot(re_sp_brarea_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    x="", y="Effect on \n log(Total S.area)")

#----------------------------#
re_sp_brarea <- ranef(mod_trunkSArea)$species
dimnames(re_sp_brarea)

re_sp_brarea_estimates <- data.frame(
  species = rownames(re_sp_brarea),
  est = re_sp_brarea[,"Estimate","Intercept"],
  lower = re_sp_brarea[,"Q2.5","Intercept"],
  upper = re_sp_brarea[,"Q97.5","Intercept"]
)

re_sp_brarea_estimates_dom <-  re_sp_brarea_estimates %>% filter(species %in% domSpeciesBCI$full.name)

re_sp_trunksarea_plot <- ggplot(re_sp_brarea_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    x="", y="Effect on \n log(Trunk S.area)")

#----------------------------#

re_sp_brarea <- ranef(mod_branchSArea)$species
dimnames(re_sp_brarea)

re_sp_brarea_estimates <- data.frame(
  species = rownames(re_sp_brarea),
  est = re_sp_brarea[,"Estimate","Intercept"],
  lower = re_sp_brarea[,"Q2.5","Intercept"],
  upper = re_sp_brarea[,"Q97.5","Intercept"]
)

re_sp_brarea_estimates_dom <-  re_sp_brarea_estimates %>% filter(species %in% domSpeciesBCI$full.name)

re_sp_brsarea_plot <- ggplot(re_sp_brarea_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    x="", y="Effect on \n log(Branch S.area)")

############# Surface Area to Volume ratio #################

re_sp_brvol <- ranef(mod_totAreaVolRatio)$species
dimnames(re_sp_brvol)

re_sp_brvol_estimates <- data.frame(
  species = rownames(re_sp_brvol),
  est = re_sp_brvol[,"Estimate","Intercept"],
  lower = re_sp_brvol[,"Q2.5","Intercept"],
  upper = re_sp_brvol[,"Q97.5","Intercept"]
)

re_sp_brvol_estimates_dom <-  re_sp_brvol_estimates %>% filter(species %in% domSpeciesBCI$full.name)

re_sp_totsareavolratio_plot <- ggplot(re_sp_brvol_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    x="", y="Effect on \n log(Total S. area / volume)")


#------------------------------------#

re_sp_brvol <- ranef(mod_trunkAreaVolRatio)$species
dimnames(re_sp_brvol)

re_sp_brvol_estimates <- data.frame(
  species = rownames(re_sp_brvol),
  est = re_sp_brvol[,"Estimate","Intercept"],
  lower = re_sp_brvol[,"Q2.5","Intercept"],
  upper = re_sp_brvol[,"Q97.5","Intercept"]
)

re_sp_brvol_estimates_dom <-  re_sp_brvol_estimates %>% filter(species %in% domSpeciesBCI$full.name)

re_sp_trunksareavolratio_plot <- ggplot(re_sp_brvol_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    x="", y="Effect on \n log(Trunk S. area / volume)")


#------------------------------------#

re_sp_brArVolRatio <- ranef(mod_brAreaVolRatio)$species
dimnames(re_sp_brArVolRatio)

re_sp_brArVolRatio_estimates <- data.frame(
  species = rownames(re_sp_brArVolRatio),
  est = re_sp_brArVolRatio[,"Estimate","Intercept"],
  lower = re_sp_brArVolRatio[,"Q2.5","Intercept"],
  upper = re_sp_brArVolRatio[,"Q97.5","Intercept"]
)

re_sp_brArVolRatio_estimates_dom <-  re_sp_brArVolRatio_estimates %>% filter(species %in% domSpeciesBCI$full.name)

re_sp_brArVolRatio_plot <- ggplot(re_sp_brArVolRatio_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    x="", y="Effect on \n log(Branch S. area / volume)")


############# Length #################

re_sp_brlen <- ranef(mod_brLen)$species
dimnames(re_sp_brlen)

re_sp_brlen_estimates <- data.frame(
  species = rownames(re_sp_brlen),
  est = re_sp_brlen[,"Estimate","Intercept"],
  lower = re_sp_brlen[,"Q2.5","Intercept"],
  upper = re_sp_brlen[,"Q97.5","Intercept"]
)

re_sp_brlen_estimates_dom <-  re_sp_brlen_estimates %>% filter(species %in% domSpeciesBCI$full.name)

re_sp_brlen_plot <- ggplot(re_sp_brlen_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    x="", y="Effect on \n log(Branch length)")


############# Number #################

re_sp_brNo <- ranef(mod_brNo)$species
dimnames(re_sp_brNo)

re_sp_brNo_estimates <- data.frame(
  species = rownames(re_sp_brNo),
  est = re_sp_brNo[,"Estimate","Intercept"],
  lower = re_sp_brNo[,"Q2.5","Intercept"],
  upper = re_sp_brNo[,"Q97.5","Intercept"]
)

re_sp_brNo_estimates_dom <-  re_sp_brNo_estimates %>% filter(species %in% domSpeciesBCI$full.name)

re_sp_brNo_plot <- ggplot(re_sp_brNo_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    x="", y="Effect on \n log(Branch Number)")

########################## Combine #################################

sp_plot_list = c(list(re_sp_totsarea_plot), list(re_sp_trunksarea_plot), list(re_sp_brsarea_plot), list(re_sp_totsareavolratio_plot), list(re_sp_trunksareavolratio_plot), 
                 list(re_sp_brArVolRatio_plot),list(re_sp_brlen_plot), list(re_sp_brNo_plot) )
sp_final_plot <- ggarrange(plotlist = sp_plot_list, 
                        ncol=3, nrow=3, 
                        common.legend=TRUE, legend="bottom",labels = c("a","b","c","d","e","f", "g","h"))

# 5) Display
print(sp_final_plot)
ggsave("qsm_re_brmetrics_combined_plot.png", width = 14, height = 14, plot = sp_final_plot, dpi = 600)



