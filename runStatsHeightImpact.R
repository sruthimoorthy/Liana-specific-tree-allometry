library(ordinal)
library(lme4)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(brms)  
library(dplyr)
library(tidyr)
library(ggpubr)
library(scales)

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


height_data_2011 <- read.csv("new_height_liana_data_2011_no_outliers_q20.csv")
height_data_2011$Species <- casefold(height_data_2011$Species, upper = FALSE)

df_2011 <- height_data_2011 %>%
  left_join(speciesCodes, by = c("Species"="Code"))

df_2011 <- data_frame(tag = df_2011$Tag,
                      quadrant = df_2011$q20,
                      species = df_2011$Full.Name,
                      dbh = df_2011$DBH_cm,
                      lianaLoad = df_2011$Lianas,
                      height = df_2011$Height,
                      source = "laser_rangefinder",
                      year = 2011)


df_2011 <- df_2011 %>%
  mutate(lianaLoad = case_when(
    lianaLoad == 0 ~ 0,
    lianaLoad <= 2 ~ 1,
    lianaLoad > 2 ~ 2,
  ))

df_all <- rbind(df_2011,df_2015,df_2019)
df_all$species <- trimws(df_all$species)

domSpeciesBCI <- read.csv("speciesTraitsBCIData.csv")


df_counts <- df_all %>%
  group_by(species, lianaLoad) %>%
  summarise(n = n(), .groups = "drop")

# 2) Pivot wider so each row is one species with columns for each category count
df_counts_wide <- df_counts %>%
  pivot_wider(
    names_from = lianaLoad, 
    values_from = n,
    values_fill = 0  # fill missing combos with 0
  )

# 3) Filter species with at least 5 observations in each category
df_counts_filtered <- df_counts_wide %>%
  filter(`0` >= 5, `1` >= 5, `2` >= 5)

# 4) Extract the species names that pass this criterion
sp_enough <- df_counts_filtered$species

df_all$lianaLoad <- factor(df_all$lianaLoad, ordered=FALSE, levels=c("0","1","2"))
df_all$species <- factor(df_all$species)

mod_height_simple <- brm(
  formula = log(height) ~ log(dbh) + lianaLoad + source + (1 |species) + (1|tag) + (1|quadrant),
  data = df_all,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

summary(mod_height_simple)

mod_height <- brm(
  formula = log(height) ~ log(dbh) + lianaLoad + source + (1 + lianaLoad |species) + (1|tag) + (1|quadrant),
  data = df_all,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

summary(mod_height)

mod_height_interaction <- brm(
  formula = log(height) ~ lianaLoad*log(dbh) + lianaLoad + source + (1 + lianaLoad |species) + (1|tag) + (1|quadrant),
  data = df_all,
  family = gaussian(),  # ordinal
  chains = 4, 
  warmup = 2000,
  iter = 4000, 
  control = list(adapt_delta = 0.95),
  cores = 4
)

summary(mod_height_interaction)

ce_source <- conditional_effects(
  mod_height_interaction, 
  effects = "dbh:source",
  re_formula = NA
)
ce_source <- plot(ce_source, plot=FALSE)

df_plot <- ce_source[[1]]$data

# Exponentiate the relevant columns if your model is log(...) ~ ...
df_plot$estimate__ <- exp(df_plot$estimate__)
df_plot$lower__    <- exp(df_plot$lower__)
df_plot$upper__    <- exp(df_plot$upper__)

# 5) Update the ggplot's data with the exponentiated values
ce_source[[1]]$data <- df_plot

plot_ce_source_custom <- ce_source[[1]] +
  labs(
    x = "DBH (cm)",
    y = "Height (m)"
  ) +
  scale_color_manual(
    name = "Sensor", 
    values = c("tls"="#4682B4","drone_photogrammetry"="#B47846","laser_rangefinder"="#B4464B"), # or your own color palette
    labels = c("tls"="TLS", "drone_photogrammetry"="Drone Photogrammetry", "laser_rangefinder"="Rangefinder")  
  ) + scale_fill_manual(
    values = c("tls"="#4682B4", "drone_photogrammetry"="#B47846", "laser_rangefinder"="#B4464B"),
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

ggsave("source_effect_on_height.png", width = 8, height = 5, plot = plot_ce_source_custom, dpi = 600)

ce_dbh_liana <- conditional_effects(
  mod_height_interaction, 
  effects = "dbh:lianaLoad",   # or "lianaLoad:logDBH"
  re_formula = NA,
  method = "fitted"# ignore random effects for pop-level visualization
)


ce_dbh_liana <- plot(ce_dbh_liana, plot=FALSE)

df_plot <- ce_dbh_liana[[1]]$data

# Exponentiate the relevant columns if your model is log(...) ~ ...
df_plot$estimate__ <- exp(df_plot$estimate__)
df_plot$lower__    <- exp(df_plot$lower__)
df_plot$upper__    <- exp(df_plot$upper__)

# 5) Update the ggplot's data with the exponentiated values
ce_dbh_liana[[1]]$data <- df_plot


plot_ce_dbh_liana_custom <- ce_dbh_liana[[1]] +
  labs(
    
    x = "DBH (cm)",
    y = "Height (m)"
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


####################### Height residual plot ##############################################

fitted_vals <- fitted(mod_height,
                      newdata = df_all,
                      re_formula = NULL, 
                      allow_new_levels = TRUE,# includes random effects
                      scale = "response", # return predictions on height (original scale)
                      summary = TRUE)

# 'fitted_vals' has columns: "Estimate", "Est.Error", "Q2.5", "Q97.5", etc.,
# with one row per row of df_all. Let's merge them into df_all:
df_all$pred_height <- exp(fitted_vals[, "Estimate"])  # mean predicted height
df_all$pred_lower  <- exp(fitted_vals[, "Q2.5"])
df_all$pred_upper  <- exp(fitted_vals[, "Q97.5"])

# 2) Compute residuals => observed minus predicted
#    Make sure df_all$height is on the same scale as the model (i.e. actual height in meters).
df_all$resid <- df_all$height - df_all$pred_height

# 3) Plot residuals
# Option A: residuals vs. fitted
p_height_resid_fit_no_inter <-ggplot(df_all, aes(x=pred_height, y=resid, color=lianaLoad)) +
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
    y="Residual = Observed - Predicted (m)"
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title = element_text(size=14),
    axis.text  = element_text(size=12)
  )



# Option B: residuals vs. DBH
p_resid_dbh_no_inter <- ggplot(df_all, aes(x=dbh, y=resid, color=lianaLoad)) +
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
    y="Residual = Observed - Predicted (m)"
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title = element_text(size=14),
    axis.text  = element_text(size=12)
  )

fitted_vals <- fitted(mod_height_interaction,
                      newdata = df_all,
                      re_formula = NULL, 
                      allow_new_levels = TRUE,# includes random effects
                      scale = "response", # return predictions on height (original scale)
                      summary = TRUE)

# 'fitted_vals' has columns: "Estimate", "Est.Error", "Q2.5", "Q97.5", etc.,
# with one row per row of df_all. Let's merge them into df_all:
df_all$pred_height <- exp(fitted_vals[, "Estimate"])  # mean predicted height
df_all$pred_lower  <- exp(fitted_vals[, "Q2.5"])
df_all$pred_upper  <- exp(fitted_vals[, "Q97.5"])

# 2) Compute residuals => observed minus predicted
#    Make sure df_all$height is on the same scale as the model (i.e. actual height in meters).
df_all$resid <- df_all$height - df_all$pred_height

# 3) Plot residuals
# Option A: residuals vs. fitted
p_height_resid_fit_w_inter <-ggplot(df_all, aes(x=pred_height, y=resid, color=lianaLoad)) +
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
    y="Residual = Observed - Predicted (m)"
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title = element_text(size=14),
    axis.text  = element_text(size=12)
  )



# Option B: residuals vs. DBH
p_resid_dbh_w_inter <- ggplot(df_all, aes(x=dbh, y=resid, color=lianaLoad)) +
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
    y="Residual = Observed - Predicted (m)"
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title = element_text(size=14),
    axis.text  = element_text(size=12)
  )


height_causal_resid_comp <- ggarrange(p_height_resid_fit_no_inter, p_resid_dbh_no_inter, p_height_resid_fit_w_inter, 
                                      p_resid_dbh_w_inter, ncol = 2, nrow = 2, common.legend=TRUE, legend="bottom", 
                                      labels = c("a", "b", "c", "d"))

ggsave("causal_height_model_residual_plots.png", width = 12, height = 12, plot = height_causal_resid_comp, dpi = 600)
# Print or save the plots
ggsave("p_height_resid_fit_w_inter.png", width = 6, height = 5, plot = p_height_resid_fit_w_inter, dpi = 600)

###################### Plot predicted data ##################################

# 1) Create a regular DBH grid and replicate for each liana load.
df_grid <- expand.grid(
  dbh = seq(20, 200, length.out=100),  # from 10 cm to 200 cm
  lianaLoad = c("0","1","2"),
  source = "drone_photogrammetry"  # or whichever baseline if your formula includes 'source'
)

# 2) Get partial-effect predictions on the *response* scale
#    re_formula=NA means ignore random effects
pred_log <- fitted(
  mod_height_interaction,
  newdata     = df_grid,
  re_formula  = NA,
  scale       = "response",  # => return predictions in meters, not log(m)
  summary     = TRUE
)
# 'pred_log' is a matrix with columns "Estimate","Est.Error","Q2.5","Q97.5", etc.

# Merge them back into df_grid
df_grid$height_est   <- exp(pred_log[,"Estimate"])
df_grid$height_lower <- exp(pred_log[,"Q2.5"])
df_grid$height_upper <- exp(pred_log[,"Q97.5"])

# 3) Plot your observed data (points) + smooth partial curves (lines + ribbons).
p_height <- ggplot() +
  # A) Observed data as points
  geom_point(
    data=df_all,
    aes(x=dbh, y=height, color=lianaLoad),
    alpha=0.3, size=0.5
  ) +
  # B) Model partial effect lines
  geom_line(
    data = df_grid %>% arrange(lianaLoad, dbh),
    aes(x=dbh, y=height_est, color=lianaLoad),
    size=1
  ) +
  # C) Optional ribbon for 95% CI
  geom_ribbon(
    data = df_grid %>% arrange(lianaLoad, dbh),
    aes(x=dbh, ymin=height_lower, ymax=height_upper, fill=lianaLoad),
    alpha=0.2
  ) +
  # Choose color/fill scale
  scale_color_manual(
    name="Liana Load",
    values=c("0"="#009e73","1"="#D55E00","2"="#000000"),
    labels=c("0"="No Liana","1"="Light","2"="Heavy")
  ) +
  scale_fill_manual(
    name="Liana Load",
    values=c("0"="#009e73","1"="#D55E00","2"="#000000"),
    labels=c("0"="No Liana","1"="Light","2"="Heavy")
  ) +
  # D) Transform axes to log, but keep numeric labels
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  labs(
    x="DBH (cm)",
    y="Height (m)"
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title      = element_text(size=14),
    axis.text       = element_text(size=12),
    plot.title      = element_text(size=14, face="bold")
  )

# Finally, print or save the plot
print(p_height)
ggsave("causal_best_height_plot.png", width = 6, height = 5, plot = p_height, dpi = 600)
#################### Plotting datapoints for source #################

# 1) Create a regular DBH grid and replicate for each liana load.
df_source_grid <- expand.grid(
  dbh = seq(20, 200, length.out=100),  # from 10 cm to 200 cm
  source = c("drone_photogrammetry","tls","laser_rangefinder"),
  lianaLoad = "0"  # or whichever baseline if your formula includes 'source'
)

# 2) Get partial-effect predictions on the *response* scale
#    re_formula=NA means ignore random effects
source_pred_log <- fitted(
  mod_height_interaction,
  newdata     = df_source_grid,
  re_formula  = NA,
  scale       = "response",  # => return predictions in meters, not log(m)
  summary     = TRUE
)
# 'pred_log' is a matrix with columns "Estimate","Est.Error","Q2.5","Q97.5", etc.

# Merge them back into df_grid
df_source_grid$height_est   <- exp(source_pred_log[,"Estimate"])
df_source_grid$height_lower <- exp(source_pred_log[,"Q2.5"])
df_source_grid$height_upper <- exp(source_pred_log[,"Q97.5"])

# 3) Plot your observed data (points) + smooth partial curves (lines + ribbons).
p_height_source_impact <- ggplot() +
  # A) Observed data as points
  geom_point(
    data=df_all,
    aes(x=dbh, y=height, color=source),
    alpha=0.4, size = 0.5
  ) +
  # B) Model partial effect lines
  geom_line(
    data = df_source_grid %>% arrange(source, dbh),
    aes(x=dbh, y=height_est, color=source),
    size=1
  ) +
  # C) Optional ribbon for 95% CI
  geom_ribbon(
    data = df_grid %>% arrange(source, dbh),
    aes(x=dbh, ymin=height_lower, ymax=height_upper, fill=source),
    alpha=0.2
  ) +
  # Choose color/fill scale
  scale_color_manual(
    name = "Sensor", 
    values = c("tls"="#4682B4","drone_photogrammetry"="darkgreen","laser_rangefinder"="#B4464B"), # or your own color palette
    labels = c("tls"="TLS", "drone_photogrammetry"="Drone Photogrammetry", "laser_rangefinder"="Rangefinder")  
  ) + scale_fill_manual(
    values = c("tls"="#4682B4", "drone_photogrammetry"="darkgreen", "laser_rangefinder"="#B4464B"),
    guide = "none"       # <- Hide the fill legend
  )  +
  # D) Transform axes to log, but keep numeric labels
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  labs(
    x="DBH (cm)",
    y="Height (m)"
  ) +
  theme_minimal() +
  theme(
    legend.position="bottom",
    axis.title      = element_text(size=14),
    axis.text       = element_text(size=12),
    plot.title      = element_text(size=14, face="bold")
  )

# Finally, print or save the plot
print(p_height_source_impact)
ggsave("source_vs_height_plot.png", width = 6, height = 5, plot = p_height_source_impact, dpi = 600)
###################### Random effects on height ########################

re_sp_height <- ranef(mod_height)$species
dimnames(re_sp_height)

re_sp_height_estimates <- data.frame(
  species = rownames(re_sp_height),
  est = re_sp_height[,"Estimate","Intercept"],
  lower = re_sp_height[,"Q2.5","Intercept"],
  upper = re_sp_height[,"Q97.5","Intercept"]
)

re_sp_height_estimates_dom <-  re_sp_height_estimates %>% filter(species %in% domSpeciesBCI$full.name)

height_re_plot <- ggplot(re_sp_height_estimates_dom, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    x="", y="Effect on \n log(height)")

ggsave("height_re_plot.png", width = 6, height = 5, plot = height_re_plot, dpi = 600)


###################### Species specific impact for heavy liana load category ########################

df_sp_lianaload2_slope <- data.frame(
  species = rownames(re_sp_height),
  est = re_sp_height[,"Estimate","lianaLoad2"],
  lower = re_sp_height[,"Q2.5","lianaLoad2"],
  upper = re_sp_height[,"Q97.5","lianaLoad2"]
)

df_sp_lianaload2_slope_sub <- subset(df_sp_lianaload2_slope, species %in% sp_enough)

height_re_lianaload2_plot <- ggplot(df_sp_lianaload2_slope_sub, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Species-Specific Slopes for LianaLoad=2",
       x="Species", 
       y="Effect of heavy liana load on \n log(height)")

ggsave("height_re_lianaload2_plot.png", width = 6, height = 5, plot = height_re_lianaload2_plot, dpi = 600)

###################### Species specific impact for light liana load category ########################

df_sp_lianaload1_slope <- data.frame(
  species = rownames(re_sp_height),
  est = re_sp_height[,"Estimate","lianaLoad1"],
  lower = re_sp_height[,"Q2.5","lianaLoad1"],
  upper = re_sp_height[,"Q97.5","lianaLoad1"]
)

df_sp_lianaload1_slope_sub <- subset(df_sp_lianaload1_slope, species %in% sp_enough)

height_re_lianaload1_plot <- ggplot(df_sp_lianaload1_slope_sub, aes(x=reorder(species, est), y=est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title="Species-Specific Slopes for LianaLoad=1",
       x="Species", 
       y="Effect of light liana load on \n log(height)")

ggsave("height_re_lianaload1_plot.png", width = 6, height = 5, plot = height_re_lianaload1_plot, dpi = 600)
