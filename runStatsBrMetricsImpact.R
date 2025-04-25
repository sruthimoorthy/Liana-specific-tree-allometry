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
       
        brLen    = rnorm(M, mean = branchLen12, sd = branchLen12SD),
        brNo     = rnorm(M, mean = branchNo12, sd = branchNo12SD)
      )
    )
  ) %>%
  unnest(cols = draws) %>%  # unnest once, in parallel
  mutate(drawID = row_number())  # or a better grouping if needed

df_long$tag <- factor(df_long$tag)

# (3) Fit the model:
mod_brLen <- brm(
  formula = bf(
    log(brLen) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)

# (3) Fit the model:
mod_brNo <- brm(
  formula = bf(
    log(brNo) ~ log(dbh) + lianaLoad  + (1  | species) + (1|tag) + (1|subplot)
  ),
  data   = df_long,
  family = gaussian(),
  chains = 4, iter=4000, warmup=2000,
  cores  = 4,
  control= list(max_treedepth=15,adapt_delta=0.99)
)


# 1) Create a named list of your models
model_list <- list(
 
  brLength        = mod_brLen,
  brNumber        = mod_brNo
  
)
# 2) A small helper function that returns a customized ggplot 
#    for each model's dbh:lianaLoad effect.
make_cond_plot <- function(fit, metricName) {
  if (metricName == "brVolume") {
    labelName = expression(paste("Br. Volume"["1,2"], ' ' ,(m^{3})))
  }else if(metricName == "brSurfaceArea"){
    labelName = expression(paste("Br. Surface Area"["1,2"], ' ' ,(m^{2})))
  }else if(metricName == "brLength"){
    labelName = expression(paste("Br. Length"["1,2"], ' ' ,(m)))
  }else if(metricName == "brNumber"){
    labelName = expression(paste("No. Branches"["1,2"], ' ', (n)))
  }else{
    labelName = expression(paste("Br. Area" ["1,2"],"/ Br. Volume"["1,2"], ' ' ,(m^{-1})))
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
  if (metricName == "brVolume") {
    df_plot$estimate__ <- exp(df_plot$estimate__)/1000
    df_plot$lower__    <- exp(df_plot$lower__)/1000
    df_plot$upper__    <- exp(df_plot$upper__)/1000
  }else{
    # Exponentiate the relevant columns if your model is log(...) ~ ...
    df_plot$estimate__ <- exp(df_plot$estimate__)
    df_plot$lower__    <- exp(df_plot$lower__)
    df_plot$upper__    <- exp(df_plot$upper__)
}
  
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
                        ncol=2, nrow=1, 
                        common.legend=TRUE, legend="bottom")

# 5) Display
print(final_plot)
ggsave("qsm_brmetrics_combined_plot.png", width = 9, height = 6, plot = final_plot, dpi = 600)
#################### Random effects of species on Volume #####################

domSpeciesBCI <- read.csv("speciesTraitsBCIData.csv")


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
    x="", y="Deviation from Overall Intercept \n (log-Branch length)")


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
    x="", y="Deviation from Overall Intercept \n (log-Branch Number)")

sp_plot_list = c(list(re_sp_brlen_plot), list(re_sp_brNo_plot))
sp_final_plot <- ggarrange(plotlist = sp_plot_list, 
                           ncol=2, nrow=1, 
                           common.legend=TRUE, legend="bottom")

# 5) Display
print(sp_final_plot)
ggsave("qsm_re_brmetrics_combined_plot.png", width = 9, height =6, plot = sp_final_plot, dpi = 600)





