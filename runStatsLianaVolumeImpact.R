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
       branchVolume = rnorm(M, mean=branchVol12, sd=branchVol12SD)
      )
    )
  ) %>%
  unnest(cols = draws) %>%  # unnest once, in parallel
  mutate(drawID = row_number())  # or a better grouping if needed



df_long$tag <- factor(df_long$tag)

# (3) Fit the model:
totVolModel <- brm(
  formula = bf(
    log(totVolume) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

trunkVolModel <- brm(
  formula = bf(
    log(trunkVolume) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

brVolModel <- brm(
  formula = bf(
    log(branchVolume) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)


domSpeciesBCI <- read.csv("speciesTraitsBCIData.csv")

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
  labs(
    x = "DBH (cm)",
    y = expression(paste("Total Volume"["0,1,2"], ' ' ,(m^{3})))
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
    labels = c("0"="No Liana", "1"="Medium", "2"="Heavy")  
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
  labs(
    x = "DBH (cm)",
    y = expression(paste("Branch Volume"["1,2"], ' ' ,(m^{3})))
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


#################### Random effects of species on Volume #####################

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
       x="", y="Deviation from Overall Intercept \n (log- total volume)")

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
    x="", y="Deviation from Overall Intercept \n (log- trunk volume)")

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
    x="", y="Deviation from Overall Intercept \n (log- branch volume)")

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
       x="Subplots", y="Deviation from Overall Intercept (log-volume)")



#################### Subplot estimate and liana load relation #####################

df_fraction <- df_long %>%
  group_by(subplot) %>%
  summarise(
    fraction_liana12 = mean(lianaLoad %in% c("1","2"))  # if lianaLoad is factor
    # or if lianaLoad is numeric, do: mean(lianaLoad > 0)
  )

################ Volume ###########ß


df_sub_re_vol_join <- left_join(re_subplot_vol_estimates, df_fraction, by="subplot")

ggplot(df_sub_re_vol_join, aes(x=fraction_liana12, y=est)) +
  geom_point() +
  geom_smooth(method="lm") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "pearson",  # or "spearman"
           label.x.npc = "left", # position the label
           label.y.npc = "top") +
  labs(x="Fraction of Trees with Lianas", y="Subplot Random Intercept") +
  theme_minimal()


