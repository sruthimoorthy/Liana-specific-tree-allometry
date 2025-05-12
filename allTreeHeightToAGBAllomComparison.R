library(ordinal)
library(lme4)
library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(brms)  
library(dplyr)
library(tidyr)
library(glue)

setwd("")

all_50ha_liana_subplot <- read.csv("liana_agb_subplot_level.csv")

all_subplot_trees <- read.csv("liana_census_all_trees_within_subplot_subqua_sorted.csv")

all_subplot_trees <- data_frame(tag = all_subplot_trees$Tag,
                                species = all_subplot_trees$Species,
                                subplot = all_subplot_trees$Subplot.number,
                                dbh = all_subplot_trees$DBH/10,
                                lianaLoad = all_subplot_trees$Lianas)

all_subplot_trees <- all_subplot_trees %>%
  mutate(lianaLoad = case_when(
    lianaLoad == 0 ~ 0,
    lianaLoad <= 50 ~ 1,
    lianaLoad > 50 ~ 2,
  ))

df_size <- all_subplot_trees %>%
  mutate(
    size_class = if_else(dbh < 40, "small", "large")  # 40 cm cutoff
  )


df_liana_summary <- df_size %>%
  group_by(subplot, size_class) %>%
  summarise(
    n_trees = n(),
    n_liana = sum(lianaLoad != "0"),  # or lianaLoad != 0 if numeric
    fraction_liana = mean(lianaLoad != "0"),  # fraction
    .groups="drop"
  )

df_size %>%
  mutate(size_class = if_else(dbh < 40, "small", "large")) %>%
  group_by(size_class) %>%
  summarise(
    n_trees = n(),
    n_liana = sum(lianaLoad != "0"),
    fraction_liana = mean(lianaLoad != "0"),
    .groups="drop"
  )


all_subplot_trees$species <- factor(all_subplot_trees$species)
all_subplot_trees$lianaLoad <- as.factor(all_subplot_trees$lianaLoad)
all_subplot_trees$subplot <- as.factor(all_subplot_trees$subplot)

sp_wood_density <- read.csv("BCI_50ha_dbh_sp_WD.csv")
sp_wood_density <- data_frame(species = sp_wood_density$species,
                              wd = sp_wood_density$wd)
sp_wood_density <- distinct(sp_wood_density)

all_subplot_trees <- merge(all_subplot_trees, sp_wood_density, by = "species")

df_actual <- all_subplot_trees
df_no <- all_subplot_trees %>%
  mutate(lianaLoad = "0")

df_actual$lianaLoad <- factor(df_actual$lianaLoad, levels=c("0","1","2"))
df_no$lianaLoad     <- factor(df_no$lianaLoad, levels=c("0","1","2"))

m_sum <- summary(mm_fit3)
resid_sd <- m_sum$spec_pars["sigma", "Estimate"]  # or check row/col indices

# correction factor
cf <- exp(0.5 * resid_sd^2)

# SCENARIO A: actual liana
#   allow_new_levels=TRUE if you have new subplots or new species not in training
post_mat_actual <- posterior_predict(
  mod_height_liana_inter_pred,
  newdata = df_actual,
  ndraws = 6000,           # request 2000 draws
  allow_new_levels = TRUE, # if subplots/species are new
  sample_new_levels = "gaussian"  # optional if you want random effect sampling for new levels
)

# post_mat_actual is a matrix (#draws x #trees)
dim(post_mat_actual)

# SCENARIO B: no-liana
post_mat_no <- posterior_predict(
  mod_height_liana_inter_pred,
  newdata = df_no,
  ndraws = 6000,
  allow_new_levels = TRUE,
  sample_new_levels = "gaussian"
)
dim(post_mat_no)

post_mat_actual <- exp(post_mat_actual) * cf
post_mat_no     <- exp(post_mat_no) * cf

wd_vec <- df_actual$wd   # length = nrow(df_actual) = #trees
# Make a matrix with same dimension
wd_mat <- matrix(wd_vec, nrow=nrow(post_mat_actual), ncol=ncol(post_mat_actual), byrow=TRUE)

dbh_vec <- df_actual$dbh
dbh_mat <- matrix(dbh_vec, nrow=nrow(post_mat_actual), ncol=ncol(post_mat_actual), byrow=TRUE)

post_biomass_actual <- 0.0559 * (wd_mat * (dbh_mat * dbh_mat) * post_mat_actual)

post_biomass_no <- 0.0559 * (wd_mat * (dbh_mat * dbh_mat) * post_mat_no)

# 2) Build a map data frame with subplot/species for each tree
df_tree_map <- df_no %>%
  mutate(tree_index = row_number()) %>%
  select(tree_index, subplot, species)

df_draws_no <- as.data.frame(post_biomass_no) %>%
  mutate(.draw = row_number()) %>%  # label each row by draw index
  pivot_longer(
    cols = -".draw",
    names_to = "tree_index", 
    values_to = "biomass"
  ) %>%
  mutate(tree_index = as.integer(gsub("V","", tree_index))) %>%
  # optionally join your df_tree_map if you want to see which subplot/species
  left_join(df_tree_map, by="tree_index")
# => columns: .draw, tree_index, biomass, subplot, species, etc.

df_draws_actual <- as.data.frame(post_biomass_actual) %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = -".draw",
    names_to = "tree_index",
    values_to = "biomass"
  ) %>%
  mutate(tree_index = as.integer(gsub("V","", tree_index))) %>%
  left_join(df_tree_map, by="tree_index")

###############################################################################
########## Site-level comparison #################
###############################################################################

df_site_no <- df_draws_no %>%
  group_by(.draw) %>%
  summarise(
    site_biomass = sum(biomass, na.rm=TRUE), 
    agbd_density = sum(biomass)/1000/5.28,
    .groups="drop"
  ) %>%
  mutate(scenario = "Liana-free")

df_site_actual <- df_draws_actual %>%
  group_by(.draw) %>%
  summarise(
    site_biomass = sum(biomass,na.rm=TRUE),
    agbd_density = sum(biomass)/1000/5.28,
    .groups="drop"
  ) %>%
  mutate(scenario="Liana-laden \n (Trees only)")


liana_agbd <- sum(all_50ha_liana_subplot$lianaAGB) / 1000 / 5.28

df_site_actual_plus <- df_draws_actual %>%
  group_by(.draw) %>%
  summarise(
    site_biomass = sum(biomass,na.rm=TRUE),
    agbd_density = (sum(biomass)/1000/5.28) + liana_agbd,
    .groups="drop"
  ) %>%
  mutate(scenario="Liana-laden \n (Trees + Lianas)")

df_site_both <- bind_rows(df_site_no, df_site_actual)  #

df_site_three <- bind_rows(df_site_no, df_site_actual,df_site_actual_plus ) 

# Make scenario a factor with no_liana first:
df_site_both$scenario <- factor(df_site_both$scenario,
                                levels=c("Liana-free","Liana-laden \n (Trees only)")
)

df_site_three$scenario <- factor(df_site_three$scenario,
                                 levels=c("Liana-free","Liana-laden \n (Trees only)","Liana-laden \n (Trees + Lianas)")
)

# We'll pivot wide:
df_no_wide <- df_site_no %>%
  rename(agb_no = agbd_density) %>%
  select(.draw, agb_no)

df_act_wide <- df_site_actual %>%
  rename(agb_act = agbd_density) %>%
  select(.draw, agb_act)

df_site_compare <- df_no_wide %>%
  inner_join(df_act_wide, by=".draw") %>%
  mutate(
    diff = agb_act - agb_no,      # e.g. negative if 'actual' is smaller than 'no'
    perc_diff = (agb_act - agb_no)/agb_no * 100
  )

# Summaries
p_neg <- mean(df_site_compare$diff < 0)    # e.g. 0.89
# coverage that excludes zero: 79% => means the central 79% interval doesn't cross zero
# (You already derived that from quantiles)
ci_90 <- quantile(df_site_compare$diff, c(0.05, 0.95))

mean_diff <- mean(df_site_compare$diff)    # average difference in mg/ha
mean_perc_diff <- mean(df_site_compare$perc_diff) 


lab_str <- glue(
  "Displaced AGB = {round(mean_diff,1)} Mg/ha [90% CI: {round(ci_90[1],2)} to {round(ci_90[2],2)}]\n",
  "Liana AGB  ~ {round(liana_agbd,1)} Mg/ha"
)

site_level_plot <- ggplot(df_site_three, aes(x=scenario, y=agbd_density,fill=scenario)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#E69F00", "#999999","darkred")) +
  labs(
    y="AGB (Mg/ha)",
    title = "50-ha site-level"
  ) +
  ylim(215, 270)+
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size=14),
        legend.text =  element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank()) +
  annotate("text",
           x=2,
           y=220,
           label=lab_str,
           size=4, hjust=0.5
  )

###############################################################################
########## Subplot-level comparison #################
###############################################################################


df_sub_no <- df_draws_no %>%
  group_by(.draw, subplot) %>%
  summarise(tree_biomass_no = sum(biomass), .groups="drop")  # 'biomass' in kg

df_sub_actual <- df_draws_actual %>%
  group_by(.draw, subplot) %>%
  summarise(tree_biomass_act = sum(biomass), .groups="drop")

subplot_area_ha <- 0.16  # adapt to your real subplot area

df_sub_no <- df_sub_no %>%
  mutate(
    tree_biomass_no_MgHa = tree_biomass_no / 1000 / subplot_area_ha
    # divides by 1000 to get Mg, then divides by 0.16 ha => Mg/ha
  )

df_sub_actual <- df_sub_actual %>%
  mutate(
    tree_biomass_act_MgHa = tree_biomass_act / 1000 / subplot_area_ha
  )


library(dplyr)

df_sub_compare <- df_sub_no %>%
  select(.draw, subplot, tree_biomass_no_MgHa) %>%
  inner_join(
    df_sub_actual %>%
      select(.draw, subplot, tree_biomass_act_MgHa),
    by=c(".draw","subplot")
  ) %>%
  mutate(
    displaced_biomass = tree_biomass_no_MgHa - tree_biomass_act_MgHa
  )

df_liana <- all_50ha_liana_subplot %>%
  mutate(
    liana_agb_MgHa = lianaAGB / 1000 / subplot_area_ha
  )
df_liana$subplot <- factor(df_liana$subplot)

df_fraction_liana <- all_subplot_trees %>%
  group_by(subplot) %>%
  summarise(
    fraction_liana = mean(lianaLoad %in% c("1","2"))  # if lianaLoad is factor
    # or if lianaLoad is numeric, do: mean(lianaLoad > 0)
  )



df_sub_final <- df_sub_compare %>%
  left_join(
    df_liana %>% select(subplot, liana_agb_MgHa),
    by="subplot"
  ) %>%
  # also join fraction_liana if you have it
  left_join(df_fraction_liana, by="subplot") %>%
  mutate(
    compensated_percent = (liana_agb_MgHa / displaced_biomass)*100
  )


df_sub_summary <- df_sub_final %>%
  group_by(subplot, fraction_liana) %>%
  summarise(
    mean_displaced = mean(displaced_biomass),
    lower_displaced = quantile(displaced_biomass, 0.025, na.rm=TRUE),
    upper_displaced = quantile(displaced_biomass, 0.975, na.rm=TRUE),
    mean_compensated = mean(liana_agb_MgHa),
    .groups="drop"
  )

df_long <- df_sub_summary %>%
  pivot_longer(
    cols=c("mean_displaced","mean_compensated"),
    names_to="metric", values_to="value"
  )

library(ggplot2)
library(ggpubr)

df_r2 <- df_long %>%
  group_by(metric) %>%
  do({
    fit <- lm(value ~ fraction_liana, data=.)
    data.frame(
      r2 = summary(fit)$r.squared,
      # store coordinates for annotation if you like
      x_pos=0.1,  # fraction_liana near the left
      y_pos=max(.$value)*0.9
    )
  })

my_colors <- c("mean_displaced"="#999999","mean_compensated"="darkred")
my_labels <- c("mean_displaced"="Liana displaced AGB","mean_compensated"="Liana AGB")


subplot_level_disp_vs_comp_plot <- ggplot(df_long, aes(x=fraction_liana, y=value, color=metric)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  
  # Now place R² text in matching color
  geom_text(
    data = df_r2,
    aes(
      x = x_pos, 
      y = y_pos, 
      label = paste0("R²=", round(r2,3)),
      color=metric
    ),
    size = 5,
    inherit.aes = FALSE,  # so we don't use fraction_liana/value from the main data
    show.legend=FALSE      # avoid separate legend for text
  ) +
  # Customize color scale
  scale_color_manual(
    name = "AGB Metric",    # legend title
    breaks=names(my_labels),
    labels=my_labels,
    values=my_colors
  ) +
  labs(
    x="Fraction of trees with lianas",
    y="AGB (Mg/ha)"      # depending on your actual units
  ) +
  theme_minimal() +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        axis.title = element_text(size=14),
        legend.text =  element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

subplot_level_compensation_plot <- ggplot(df_sub_summary, aes(x=mean_displaced, y=mean_compensated)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm", se=FALSE, color="black") +
  geom_smooth(
    data = df_sub_summary, 
    aes(x = mean_displaced, y = mean_displaced), 
    method = "lm", 
    se = FALSE, 
    linetype = "dashed", 
    color = "darkblue"
  ) +
  annotate(
    "text", 
    x = min(df_sub_summary$mean_displaced, na.rm = TRUE), 
    y = max(df_sub_summary$mean_displaced, na.rm = TRUE), 
    label = "1:1 reference line", 
    hjust = 0, vjust = 1,
    size = 5, fontface = "italic", color = "darkblue"
  ) +
  stat_cor(
    method = "pearson", 
    label.x.npc = "left",   # position correlation label
    label.y = max(df_sub_summary$mean_displaced-5, na.rm = TRUE),
    size = 5 # near top
  ) +
  labs(
    x="Mean liana displaced AGB (Mg/ha)",
    y="Liana AGB (Mg/ha)"
  ) +
  
  theme_minimal() +
  theme( axis.title = element_text(size=14),
         legend.text =  element_text(size=14),
         axis.text.x = element_text(size=12),
         axis.text.y = element_text(size=12))


################## Basal Area of infestation ####################

df_basal_infest <- all_subplot_trees %>%
  # keep only liana-laden trees
  filter(lianaLoad %in% c("1","2")) %>%
  # compute basal area in cm^2: π*(dbh/2)^2
  mutate(
    basal_area_m2 =( pi * (dbh / 2)^2)/10000
  ) %>%
  group_by(subplot) %>%
  summarise(
    # sum over all liana-infested trees in that subplot
    basal_area_liana = sum(basal_area_m2, na.rm=TRUE),
    .groups="drop"
  )

df_sub_summary2 <- df_sub_summary %>%
  left_join(df_basal_infest, by="subplot")

df_long_basal <- df_sub_summary2 %>%
  tidyr::pivot_longer(
    cols = c("mean_displaced","mean_compensated"),
    names_to = "metric",
    values_to = "value"
  )

df_r2_basal <- df_long_basal %>%
  group_by(metric) %>%
  do({
    fit <- lm(value ~ basal_area_liana, data=.)
    data.frame(
      r2 = summary(fit)$r.squared,
      # pick coordinates for annotation
      x_pos = max(.$basal_area_liana)*0.05,
      y_pos = max(.$value)*0.9
    )
  })


my_colors <- c("mean_displaced"="#999999","mean_compensated"="darkred")
my_labels <- c("mean_displaced"="Liana-displaced AGB","mean_compensated"="Liana AGB")

subplot_level_disp_vs_comp_basal_plot <- ggplot(df_long_basal, 
                                                aes(x=basal_area_liana, y=value, color=metric)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  
  # Annotate R² in matching color
  geom_text(
    data = df_r2_basal,
    aes(
      x = 0.3,
      y = y_pos,
      label = paste0("R²=", round(r2,3)),
      color = metric
    ),
    inherit.aes = FALSE,
    size = 5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name = "AGB Metric",
    breaks = names(my_labels),
    labels = my_labels,
    values = my_colors
  ) +
  labs(
    x="Basal area (m²) of liana-infested trees per subplot",
    y="AGB (Mg/ha)"
  ) +
  theme_minimal() +
  theme(
    legend.position="top",
    legend.title=element_blank(),
    axis.title=element_text(size=14),
    legend.text=element_text(size=14),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12)
  )

subplot_level_disp_vs_comp_basal_plot

#################### Only plotting displaced tree biomass #################

# 1) Fit a simple linear model
fit_fraction <- lm(mean_displaced ~ fraction_liana, data=df_sub_summary)
r2_fraction <- summary(fit_fraction)$r.squared

# 2) Plot
plot_fraction_disp <- ggplot(df_sub_summary, aes(x=fraction_liana, y=mean_displaced)) +
  geom_point(alpha=0.6) +
  geom_smooth(method="lm", se=FALSE, color="black") +
  annotate(
    "text",
    x = 0.05,  # near left side, adapt as needed
    y = max(df_sub_summary$mean_displaced, na.rm=TRUE)*0.95,
    label = paste0("R²=", round(r2_fraction, 3)),
    size=5,
    color="black",
    hjust=0
  ) +
  labs(
    x="Fraction of trees with lianas",
    y="Mean liana-displaced\n AGB (Mg/ha)"
  ) +
  theme_minimal() +
  theme(
    axis.title=element_text(size=14),
    axis.text=element_text(size=12)
  )


# 1) Fit a simple model
fit_basal <- lm(mean_displaced ~ basal_area_liana, data=df_sub_summary2)
r2_basal <- summary(fit_basal)$r.squared

# 2) Plot
plot_basal_disp <- ggplot(df_sub_summary2, aes(x=basal_area_liana, y=mean_displaced)) +
  geom_point(alpha=0.6) +
  geom_smooth(method="lm", se=FALSE, color="black") +
  annotate(
    "text",
    x = min(df_sub_summary2$basal_area_liana, na.rm=TRUE)*1.05,
    y = max(df_sub_summary2$mean_displaced, na.rm=TRUE)*0.95,
    label = paste0("R²=", round(r2_basal,3)),
    size=5,
    color="black",
    hjust=0
  ) +
  labs(
    x="Basal area liana-laden trees (m²)",
    y="Mean liana-displaced\n AGB (Mg/ha)"
  ) +
  theme_minimal() +
  theme(
    axis.title=element_text(size=14),
    axis.text=element_text(size=12)
  )

