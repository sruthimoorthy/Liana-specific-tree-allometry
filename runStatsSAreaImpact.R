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
    # Generate a single tibble of size M with columns brVolume, brArea, brLen, brNo
    draws = list(
      tibble(
        totVolume = rnorm(M, mean=totVol, sd=totVolSD),
        trunkVolume = rnorm(M, mean=trunkVol, sd=trunkVolSD),
        branchVolume = rnorm(M, mean=branchVol12, sd=branchVol12SD),
        totSArea = rnorm(M, mean=totArea, sd=totAreaSD),
        trunkSArea = rnorm(M, mean=trunkArea, sd=trunkAreaSD),
        branchSArea = rnorm(M, mean=branchArea12, sd=branchArea12SD)
      )
    )
  ) %>%
  unnest(cols = draws) %>%  # unnest once, in parallel
  mutate(drawID = row_number())  # or a better grouping if needed

df_long$tag <- factor(df_long$tag)

df_long$totAreaVolRatio <- df_long$totSArea/(df_long$totVolume/1000)
df_long$trunkAreaVolRatio <- df_long$trunkSArea/(df_long$trunkVolume/1000)
df_long$brAreaVolRatio <- df_long$branchSArea/(df_long$branchVolume/1000)

# (3) Fit the model:
mod_totSArea <- brm(
  formula = bf(
    log(totSArea) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

mod_trunkSArea <- brm(
  formula = bf(
    log(trunkSArea) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

# (3) Fit the model:
mod_branchSArea <- brm(
  formula = bf(
    log(branchSArea) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_long,
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
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

mod_trunkAreaVolRatio <- brm(
  formula = bf(
    log(trunkAreaVolRatio) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_long,
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
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

# 1) Create a named list of your models
model_list <- list(
  totSArea        = mod_totSArea,
  trunkSArea      = mod_trunkSArea,
  branchSArea      = mod_branchSArea,
  totAreaVolRatio   = mod_totAreaVolRatio,
  trunkAreaVolRatio = mod_trunkAreaVolRatio,
  brAreaVolRatio  = mod_brAreaVolRatio
)
# 2) A small helper function that returns a customized ggplot 
#    for each model's dbh:lianaLoad effect.
make_cond_plot <- function(fit, metricName) {
  if (metricName == "totSArea") {
    labelName = expression(paste("Total Surface Volume"["0,1,2"], ' ' ,(m^{2})))
  }else if(metricName == "trunkSArea"){
    labelName = expression(paste("Trunk Surface Area"["0"], ' ' ,(m^{2})))
  }else if (metricName == "branchSArea") {
    labelName = expression(paste("Branch Surface Volume"["0,1,2"], ' ' ,(m^{2})))
  }else if(metricName == "totAreaVolRatio"){
    labelName = expression(paste("Total S. Area" ["0,1,2"],"/ Total Volume"["0,1,2"], ' ' ,(m^{-1})))
  }else if (metricName == "trunkAreaVolRatio") {
    labelName = expression(paste("Trunk S. Area" ["0"],"/ Trunk Volume"["0"], ' ' ,(m^{-1})))
  }else if(metricName == "brAreaVolRatio"){
    labelName = expression(paste("Br. S. Area" ["1,2"],"/ Br. Volume"["1,2"], ' ' ,(m^{-1})))
  }
  # compute conditional effects for dbh:lianaLoad, ignoring random effects
  ce_obj <- conditional_effects(
    fit, 
    effects = "dbh:lianaLoad",   # or "lianaLoad:logDBH"
    re_formula = NA 
  )
  
  # The "plot(..., plot=FALSE)" returns a list of ggplots, typically we want [[1]]
  ce_plot <- plot(ce_obj, plot=FALSE)[[1]]
  
  # 4) Exponentiate the y-data so we display on original scale
  #    Usually the data columns are named estimate__, lower__, upper__, 
  #    plus maybe effect1__, effect2__ for x-values, etc.
  df_plot <- ce_plot$data
  
    # Exponentiate the relevant columns if your model is log(...) ~ ...
  df_plot$estimate__ <- exp(df_plot$estimate__)
  df_plot$lower__    <- exp(df_plot$lower__)
  df_plot$upper__    <- exp(df_plot$upper__)
  
  
  # 5) Update the ggplot's data with the exponentiated values
  ce_plot$data <- df_plot
  
  # 6) Now we can add custom labels, color scales, etc.
  p_custom <- ce_plot +
    labs(
      x = "DBH (cm)",
      y = labelName
    ) +
    scale_color_manual(
      name   = "Liana Load",
      values = c("0"="#009e73","1"="#D55E00","2"="#000000"),
      labels = c("0"="No Liana", "1"="Medium", "2"="Heavy")
    ) +
    scale_fill_manual(
      values = c("0"="#009e73","1"="#D55E00","2"="#000000"),
      guide = "none"
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
  
  return(p_custom)
}

# 3) Loop over each model, generate the conditional-effects plot
plot_list <- lapply(names(model_list), function(mName) {
  make_cond_plot(model_list[[mName]], mName)
})
# 4) Combine them in a grid
# if you have 5 plots, 2 rows x 3 columns might work
# or just let ggarrange decide
final_plot <- ggarrange(plotlist = plot_list, 
                        ncol=3, nrow=2, 
                        common.legend=TRUE, legend="bottom")

# 5) Display
print(final_plot)
ggsave("qsm_surfArea_combined_plot.png", width = 11, height = 8, plot = final_plot, dpi = 600)

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
    x="", y="Deviation from Overall Intercept \n (log-Total S.area)")

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
    x="", y="Deviation from Overall Intercept \n (log-Trunk S.area)")

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
    x="", y="Deviation from Overall Intercept \n (log-Branch S.area)")

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
    x="", y="Deviation from Overall Intercept \n (log-Total S. area / volume)")


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
    x="", y="Deviation from Overall Intercept \n (log-Trunk S. area / volume)")


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
    x="", y="Deviation from Overall Intercept \n (log-Branch S. area / volume)")


sp_plot_list = c(list(re_sp_totsarea_plot), list(re_sp_trunksarea_plot), list(re_sp_brsarea_plot), list(re_sp_totsareavolratio_plot), list(re_sp_trunksareavolratio_plot), list(re_sp_brArVolRatio_plot))
sp_final_plot <- ggarrange(plotlist = sp_plot_list, 
                        ncol=3, nrow=2, 
                        common.legend=TRUE, legend="bottom")

# 5) Display
print(sp_final_plot)
ggsave("qsm_re_surfArea_combined_plot.png", width = 14, height = 11, plot = sp_final_plot, dpi = 600)



