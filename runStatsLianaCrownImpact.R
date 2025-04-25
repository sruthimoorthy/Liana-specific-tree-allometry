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
       x="", y="Deviation from Overall Intercept \n (log-crown area)")

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
       x="", y="Deviation from Overall Intercept \n (log-crown volume)")

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
       x="", y="Deviation from Overall Intercept \n (log-crown depth)")

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
       x="Subplots", y="Deviation from Overall Intercept (log-cpa)")

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
       x="Subplots", y="Deviation from Overall Intercept (log-crown volume)")

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
       x="Subplots", y="Deviation from Overall Intercept (log-crown depth)")

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

ggplot(df_sub_re_cvol_join, aes(x=fraction_liana12, y=est)) +
  geom_point() +
  geom_smooth(method="lm") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "pearson",  # or "spearman"
           label.x.npc = "left", # position the label
           label.y.npc = "top") +
  labs(x="Fraction of Trees with Lianas", y="Subplot Random Intercept") +
  theme_minimal()
