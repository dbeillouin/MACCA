###############################################################################
# Project   : MACCA Database – Carbon Stock Analysis in Agroforestry Systems
# Author    : Damien Beillouin
# Date      : 30 July 2024
# Update    : September 2025
# Purpose   : Clean and prepare SOC stock data for LMP experiments comparing
#             Full Sun vs. Agroforestry systems.
# Notes     : FAIR principles: documented, reproducible, structured, transparent
# Focus on sequestration rates
###############################################################################

# --- Load libraries ----------------------------------------------------------
library(dplyr)
library(forcats)
library(ggplot2)
library(readr)
library(metafor)

# ========================== 1. Load data ====================================
data_raw <- read_csv("~/Documents/MACCA/Data_for_analysis.csv", col_types = cols())

# Filter LMP experiments with Full Sun control
FULL_sun <- data_raw %>%
  filter(experiment_type == "LMP",
         control_type == "Full sun",
         !is.na(delta_stock_T) | Keep_rate == "YES")

# Convert numeric columns
FULL_sun <- FULL_sun %>%
  mutate(delta_stock_C.yr = as.numeric(gsub(",", ".", delta_stock_C.yr)),
         delta_stock_T   = as.numeric(delta_stock_T))

# Correct known data entry
FULL_sun[FULL_sun$ID == "442", "control_soc_mean_T_ha"] <- 40.5

# ========================== 2. Compute derived variables =====================
FULL_sun <- FULL_sun %>%
  mutate(
    # Adjust delta_stock_T
    delta_stock_T = if_else(!is.na(delta_stock_T) & !is.na(delta_stock_C.yr),
                            delta_stock_T - delta_stock_C.yr,
                            delta_stock_T),
    # Assign design for SOC change
    design = if_else(!is.na(delta_stock_T), "BACI", Design),
    # Sequestration rate
    seq_rate = if_else(is.na(delta_stock_T),
                       (treatment_soc_mean_T_ha - control_soc_mean_T_ha) / time_since_conversion,
                       delta_stock_T),
    seq_rate_sd = sqrt(treatment_soc_sd_T_ha^2 + control_soc_sd_T_ha^2) / time_since_conversion,
    vi = seq_rate_sd,
    seq_rate_vi = seq_rate_sd^2,
    # Fill history
    history_C = if_else(!is.na(delta_stock_T), history_T, history_C)
  )

# ========================== 3. Treatment categorization =====================
FULL_sun <- FULL_sun %>%
  mutate(
    NEW_treatment_type = fct_lump_min(NEW_treatment_type, 10),
    NEW_treatment_type2 = as.character(NEW_treatment_type),
    NEW_treatment_type2 = if_else(NEW_treatment_type2 %in% c("Alley cropping", "Hedgerow"),
                                  "Alley/Hedgerow", NEW_treatment_type2)
  )

# Meta-analysis: effect of treatment
m_treatment_orig <- rma.mv(seq_rate, seq_rate_vi,
                           mods = ~ -1 + NEW_treatment_type,
                           random = list(~ 1 | id_article, ~1 | id_expérimentation),
                           data = FULL_sun, method = "ML")
cat("AIC m_treatment_orig:", AIC(m_treatment_orig), "\n")

m_treatment_grouped <- rma.mv(seq_rate, seq_rate_vi,
                              mods = ~ -1 + NEW_treatment_type2,
                              random = list(~ 1 | id_article, ~1 | id_expérimentation),
                              data = FULL_sun, method = "ML")
cat("AIC m_treatment_grouped:", AIC(m_treatment_grouped), "\n")

# ========================== 4. Biodiversity & Woody diversity ==============
FULL_sun <- FULL_sun %>%
  mutate(
    Woody_category = case_when(
      treatment_N_Woody < 1.5 ~ "1 species",
      treatment_N_Woody >= 2 ~ "2+ species"
    ),
    diff_species = treatment_NB_species - as.numeric(control_NB_species),
    diff_species_class = cut(diff_species, breaks = c(-5, 1.1, 25),
                             labels = c("A:<0/+1","C+2+"), right = FALSE),
    diff_species_class = if_else(is.na(diff_species_class), "C+2+", diff_species_class),
    diff_Woody = treatment_N_Woody - control_N_Woody,
    diff_Woody_cat = case_when(
      diff_Woody >= 0 & diff_Woody <= 2 ~ "1",
      diff_Woody > 2 ~ "2+",
      TRUE ~ "0"
    )
  )

# Meta-analysis
m_woody <- rma.mv(seq_rate, seq_rate_vi, mods = ~ -1 + Woody_category,
                  random = list(~1|id_article, ~1|id_expérimentation),
                  data = FULL_sun, method = "ML")
cat("AIC m_woody:", AIC(m_woody), "\n")

m_diff_woody <- rma.mv(seq_rate, seq_rate_vi, mods = ~ -1 + diff_Woody_cat,
                       random = list(~1|id_article, ~1|id_expérimentation),
                       data = FULL_sun, method = "ML")
cat("AIC m_diff_woody:", AIC(m_diff_woody), "\n")

m_diff_species <- rma.mv(seq_rate, seq_rate_vi, mods = ~ -1 + diff_species_class,
                         random = list(~1|id_article, ~1|id_expérimentation),
                         data = FULL_sun, method = "ML")
cat("AIC m_diff_species:", AIC(m_diff_species), "\n")

# ========================== 5. Land-use history =============================
FULL_sun <- FULL_sun %>%
  mutate(Hist = paste(history_C, history_T),
         History_reclass = case_when(
           Hist %in% c("cropland forest", "forest cropland","NA NA",
                       "unknown pasture","unknown unknown") ~ "Unknown/mixed",
           Hist == "cropland cropland" ~ "Cropland",
           Hist == "forest forest" ~ "Forest",
           Hist %in% c("grassland grassland","pasture pasture") ~ "Grassland",
           TRUE ~ "NA"),
         Hist2 = paste(History_reclass, control_Land_Use),
         Hist2 = if_else(Hist2 %in% c("Forest Grassland","Grassland Grassland"),
                         "Forest/Grassland Grassland", Hist2))

# Meta-analysis
m_history <- rma.mv(seq_rate, seq_rate_vi, mods = ~ -1 + Hist2,
                    random = list(~1|id_article, ~1|id_expérimentation),
                    data = FULL_sun, method = "ML")
cat("AIC m_history:", AIC(m_history), "\n")

# ========================== 6. Climate variables ===========================
FULL_sun <- FULL_sun %>%
  mutate(FAO_Climate = fct_lump_min(FAO_Climate, 10),
         Temp_class = cut(temperature, breaks = c(0,20,25,37)))

m_climate <- rma.mv(seq_rate, seq_rate_vi, mods = ~ FAO_Climate,
                    random = list(~1|id_article, ~1|id_expérimentation),
                    data = FULL_sun, method = "ML")
cat("AIC m_climate:", AIC(m_climate), "\n")

m_temp <- rma.mv(seq_rate, seq_rate_vi, mods = ~ -1 + Temp_class,
                 random = list(~1|id_article, ~1|id_expérimentation),
                 data = FULL_sun, method = "ML")
cat("AIC m_temp:", AIC(m_temp), "\n")

m_precip <- rma.mv(seq_rate, seq_rate_vi, mods = ~ -1 + precipitation,
                   random = list(~1|id_article, ~1|id_expérimentation),
                   data = FULL_sun, method = "ML")
cat("AIC m_precip:", AIC(m_precip), "\n")

# ========================== 7. Main crop classification ====================
FULL_sun <- FULL_sun %>%
  mutate(main_culture2 = case_when(
    grepl("Cocoa|Coffee|cocoa", main_culture) ~ "Cocoa/Coffee",
    grepl("Beans|Maize", main_culture) ~ "Beans/Maize",
    grepl("Banana", main_culture) ~ "Others",
    TRUE ~ as.character(main_culture)
  ))

m_crop <- rma.mv(seq_rate, seq_rate_vi, mods = ~ -1 + main_culture2,
                 random = list(~1|id_article, ~1|id_expérimentation),
                 data = FULL_sun, method = "ML")
cat("AIC m_crop:", AIC(m_crop), "\n")

# ========================== 8. Soil depth ==================================
FULL_sun <- FULL_sun %>%
  mutate(MEAN_depthD = cut(MEAN_depth, breaks = c(0,20,40,60,80)))

m_depth <- rma.mv(seq_rate, seq_rate_vi, mods = ~ -1 + MEAN_depthD,
                  random = list(~1|id_article, ~1|id_expérimentation),
                  data = FULL_sun, method = "ML")
cat("AIC m_depth:", AIC(m_depth), "\n")

# ========================== 9. Time since conversion ======================
FULL_sun <- FULL_sun %>%
  mutate(time_since_conversionD = cut(time_since_conversion, breaks = c(0,11,20,80),
                                      include.lowest = TRUE),
         time_since_conversionD = as.character(time_since_conversionD),
         time_since_conversionD = if_else(is.na(time_since_conversionD), "unknown", time_since_conversionD))

m_time_cont <- rma.mv(seq_rate, seq_rate_vi, mods = ~ -1 + time_since_conversion,
                      random = list(~1|id_article, ~1|id_expérimentation),
                      data = FULL_sun, method = "ML")


# ========================== 10. Experimental design =========================
FULL_sun <- FULL_sun %>%
  mutate(Grouped_Design = case_when(
    is.na(Design) ~ "Non spécifié",
    Design %in% c("RCBD", "RCT") ~ "Randomized Designs",
    Design %in% c("BACI", "BA") ~ "Before-After Designs",
    Design == "CI" ~ "Control-Impact Designs",
    TRUE ~ "Autre"
  ))

# Random-effects only model
m_design <- rma.mv(seq_rate, seq_rate_vi,
                   random = list(~1|id_article, ~1|id_expérimentation),
                   data = FULL_sun, method = "ML")
cat("AIC m_design:", AIC(m_design), "\n")

# ========================== 11. Quality score ===============================
FULL_sun <- FULL_sun %>%
  mutate(quality_score = rowSums(cbind(
    (!is.na(history_C) & history_C != "unknown")*1,
    is.na(Confunding_factors)*1,
    (Grouped_Design %in% c("Before-After Designs","Randomized Designs"))*1
  )))

m_quality <- rma.mv(seq_rate, seq_rate_vi, mods = ~ -1 + quality_score,
                    random = list(~1|id_article, ~1|id_expérimentation),
                    data = FULL_sun, method = "ML")
cat("AIC m_quality:", AIC(m_quality), "\n")

# Null model
m_null <- rma.mv(seq_rate, seq_rate_vi,
                 random = list(~1|id_article, ~1|id_expérimentation),
                 data = FULL_sun, method = "ML")
cat("AIC m_null:", AIC(m_null), "\n")

# ========================== 12. Impute missing SDs ==========================
# Flag missing treatment SD
FULL_sun <- FULL_sun %>%
  mutate(estSD = if_else(is.na(treatment_soc_sd_T_ha), TRUE, FALSE))

# Impute treatment SD by CV*mean per Grouped_Design
TAB <- FULL_sun %>%
  group_by(Grouped_Design) %>%
  summarise(Mean_m = mean(treatment_soc_mean_T_ha, na.rm = TRUE),
            SD_m = mean(treatment_soc_sd_T_ha, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(CV_m = SD_m / Mean_m * 1.5)

FULL_sun <- FULL_sun %>%
  left_join(TAB, by = "Grouped_Design") %>%
  mutate(treatment_soc_sd_T_ha = if_else(is.na(treatment_soc_sd_T_ha),
                                         CV_m * treatment_soc_mean_T_ha,
                                         treatment_soc_sd_T_ha)) %>%
  dplyr::select(-Mean_m, -SD_m, -CV_m)

# Impute control SD by CV*mean per Grouped_Design
TAB <- FULL_sun %>%
  group_by(Grouped_Design) %>%
  summarise(Mean_m = mean(control_soc_mean_T_ha, na.rm = TRUE),
            SD_m = mean(control_soc_sd_T_ha, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(CV_m = SD_m / Mean_m * 1.5)

FULL_sun <- FULL_sun %>%
  left_join(TAB, by = "Grouped_Design") %>%
  mutate(control_soc_sd_T_ha = if_else(is.na(control_soc_sd_T_ha),
                                       CV_m * control_soc_mean_T_ha,
                                       control_soc_sd_T_ha)) %>%
  dplyr::select(-Mean_m, -SD_m, -CV_m)

# ========================== 13. Recalculate sequestration rates =============
FULL_sun <- FULL_sun %>%
  mutate(
    delta_stock_T = as.numeric(delta_stock_T),
    seq_rate = if_else(is.na(delta_stock_T),
                       (treatment_soc_mean_T_ha - control_soc_mean_T_ha)/time_since_conversion,
                       delta_stock_T),
    seq_rate_sd = sqrt(treatment_soc_sd_T_ha^2 + control_soc_sd_T_ha^2)/time_since_conversion,
    history_C = if_else(!is.na(delta_stock_T), history_T, history_C),
    vi = seq_rate_sd,
    seq_rate_vi = seq_rate_sd^2
  )


# ========================== 10. Experimental design =========================
FULL_sun <- FULL_sun %>%
  mutate(Grouped_Design = case_when(
    is.na(Design) ~ "Non spécifié",
    Design %in% c("RCBD", "RCT") ~ "Randomized Designs",
    Design %in% c("BACI", "BA") ~ "Before-After Designs",
    Design == "CI" ~ "Control-Impact Designs",
    TRUE ~ "Autre"
  ))

# Random-effects only model
m_design <- rma.mv(seq_rate, seq_rate_vi,
                   random = list(~1|id_article, ~1|id_expérimentation),
                   data = FULL_sun, method = "ML")
cat("AIC m_design:", AIC(m_design), "\n")

# ========================== 11. Quality score ===============================
FULL_sun <- FULL_sun %>%
  mutate(quality_score = rowSums(cbind(
    (!is.na(history_C) & history_C != "unknown")*1,
    is.na(Confunding_factors)*1,
    (Grouped_Design %in% c("Before-After Designs","Randomized Designs"))*1
  )))

m_quality <- rma.mv(seq_rate, seq_rate_vi, mods = ~ -1 + quality_score,
                    random = list(~1|id_article, ~1|id_expérimentation),
                    data = FULL_sun, method = "ML")
cat("AIC m_quality:", AIC(m_quality), "\n")

# Null model
m_null <- rma.mv(seq_rate, seq_rate_vi,
                 random = list(~1|id_article, ~1|id_expérimentation),
                 data = FULL_sun, method = "ML")
cat("AIC m_null:", AIC(m_null), "\n")

# ========================== 12. Impute missing SDs ==========================
# Flag missing treatment SD
FULL_sun <- FULL_sun %>%
  mutate(estSD = if_else(is.na(treatment_soc_sd_T_ha), TRUE, FALSE))

# Impute treatment SD by CV*mean per Grouped_Design
TAB <- FULL_sun %>%
  group_by(Grouped_Design) %>%
  summarise(Mean_m = mean(treatment_soc_mean_T_ha, na.rm = TRUE),
            SD_m = mean(treatment_soc_sd_T_ha, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(CV_m = SD_m / Mean_m * 1.5)

FULL_sun <- FULL_sun %>%
  left_join(TAB, by = "Grouped_Design") %>%
  mutate(treatment_soc_sd_T_ha = if_else(is.na(treatment_soc_sd_T_ha),
                                         CV_m * treatment_soc_mean_T_ha,
                                         treatment_soc_sd_T_ha)) %>%
  dplyr::select(-Mean_m, -SD_m, -CV_m)

# Impute control SD by CV*mean per Grouped_Design
TAB <- FULL_sun %>%
  group_by(Grouped_Design) %>%
  summarise(Mean_m = mean(control_soc_mean_T_ha, na.rm = TRUE),
            SD_m = mean(control_soc_sd_T_ha, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(CV_m = SD_m / Mean_m * 1.5)

FULL_sun <- FULL_sun %>%
  left_join(TAB, by = "Grouped_Design") %>%
  mutate(control_soc_sd_T_ha = if_else(is.na(control_soc_sd_T_ha),
                                       CV_m * control_soc_mean_T_ha,
                                       control_soc_sd_T_ha)) %>%
  dplyr::select(-Mean_m, -SD_m, -CV_m)

# ========================== 13. Recalculate sequestration rates =============
FULL_sun <- FULL_sun %>%
  mutate(
    delta_stock_T = as.numeric(delta_stock_T),
    seq_rate = if_else(is.na(delta_stock_T),
                       (treatment_soc_mean_T_ha - control_soc_mean_T_ha)/time_since_conversion,
                       delta_stock_T),
    seq_rate_sd = sqrt(treatment_soc_sd_T_ha^2 + control_soc_sd_T_ha^2)/time_since_conversion,
    history_C = if_else(!is.na(delta_stock_T), history_T, history_C),
    vi = seq_rate_sd,
    seq_rate_vi = seq_rate_sd^2
  )



# ========================== 14 Data checks / visualization =================
FULL_sun <- FULL_sun %>%
  mutate(yie = seq_rate,
         COL = "1")

qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)

library(ggExtra)
library(ggpubr)

scatter_plot <- ggplot(FULL_sun) +
  geom_point(aes(x = vi, y = yie, fill = Grouped_Design, color = Grouped_Design, size = 1/vi),
             shape = 21, color = "black") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
  labs(x = "Variance", y = "Sequestration rate (Mg.ha.an-1)") +
  scale_size(range = c(1, 6)) +
  theme_pubr(legend = "bottom")

scatter_plot_marginal <- ggMarginal(scatter_plot, type = "boxplot", groupFill = TRUE)
scatter_plot_marginal


# ========================== 15 Variance structure: scale models =================
# Goal: Test whether the variance of sequestration rates depends on
#       time since conversion

# --- 1. Scale model: variance as a function of time since conversion ---
m_scale <- rma(
  yi = yie,
  vi = vi,
  scale = ~ time_since_conversion,
  data = FULL_sun
)

# --- 2. Null model: constant variance ---
m_null <- rma(
  yi = yie,
  vi = vi,
  scale = ~1,
  data = FULL_sun
)

# --- 3. Compare models using likelihood ratio test ---
anova_result <- anova(m_scale, m_null)
print(anova_result)


# ========================== Interpretation =================
# Comparison of scale models: variance as a function of time since conversion
#
# Full model (m_scale): variance depends on time_since_conversion
# Reduced model (m_null): variance is constant
#
# Key points from anova_result:
# - LRT = 19.70, p < 0.0001 → variance of sequestration rates depends significantly on time since conversion
# - AIC improves by ~17.7 points → full model fits the data much better
# - Residual heterogeneity (QE = 332.56) remains high → other factors may still explain variability
#
# Practical takeaway:
# The effect of time since conversion on variance is real and statistically robust.
# Analyses and predictions should account for this scale effect.



# Extract coefficients from scale model
scale_coefs <- coef(m_scale)  # intercept and slope

# Extract numeric coefficients
beta0 <- m_scale$alpha["intrcpt", 1]                   # intercept
beta1 <- m_scale$alpha["time_since_conversion", 1]

# Predicted variance at time_since_conversion = 1 and 10
sigma2_1  <- exp(beta0 + beta1*1)
sigma2_10 <- exp(beta0 + beta1*10)

# Relative reduction (%)
reduction_pct <- 100 * (1 - sigma2_10 / sigma2_1)

cat("Predicted variance at time = 1:", round(sigma2_1,4), "\n")
cat("Predicted variance at time = 10:", round(sigma2_10,4), "\n")
cat("Relative reduction in heterogeneity:", round(reduction_pct,1), "%\n")



# Create a sequence of times
time_seq <- seq(0, 20, by = 0.1)  # from 0 to 20 years, adjust as needed

# Compute predicted variance
pred_var <- exp(beta0 + beta1 * time_seq)

# Combine into a data frame for plotting
df_var <- data.frame(
  time_since_conversion = time_seq,
  predicted_variance = pred_var
)

# Plot
p <- ggplot(df_var, aes(x = time_since_conversion, y = predicted_variance)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(data = data.frame(time_since_conversion = c(1,10),
                               predicted_variance = exp(beta0 + beta1 * c(1,10))),
             aes(x = time_since_conversion, y = predicted_variance),
             color = "red", size = 3) +
  labs(
    title = "Predicted Variance of Sequestration Rates over Time Since Conversion",
    x = "Time Since Conversion (years)",
    y = "Predicted Variance"
  ) +
  theme_minimal(base_size = 14)

# Add marginal density for time_since_conversion
ggMarginal(p, type = "density", margins = "x", fill="grey80", alpha=0.5)







# ========================== 16. Variance structure by design =================
# Goal: Show how predicted variance changes with time since conversion,
#       separately for each Grouped_Design


# --- 1. Scale model: variance as a function of time since conversion ---
m_scale <- rma(
  yi = yie,
  vi = vi,
  scale = ~ Grouped_Design+ time_since_conversion,
  data = FULL_sun
)


m_scale2 <- rma(
  yi = yie,
  vi = vi,
  scale = ~ Grouped_Design* time_since_conversion,
  data = FULL_sun
)
anova_result <- anova(m_scale, m_scale2)
print(anova_result)
AIC(m_scale)



# --- 2. Null model: constant variance ---
m_null <- rma(
  yi = yie,
  vi = vi,
  scale = ~1,
  data = FULL_sun
)

# --- 3. Compare models using likelihood ratio test ---
anova_result <- anova(m_scale, m_null)
print(anova_result)



# --- Extraire les coefficients de la partie scale ---
scale_coefs <- coef(m_scale)$alpha  # vecteur nommé

# --- Définir les niveaux du facteur de design ---
levels_design <- c("Reference", "Control-Impact", "Randomized")

# --- Calculer les log-variances par niveau ---
logvar <- c(
  scale_coefs["intrcpt"],                                            # Reference
  scale_coefs["intrcpt"] + scale_coefs["Grouped_DesignControl-Impact Designs"],
  scale_coefs["intrcpt"] + scale_coefs["Grouped_DesignRandomized Designs"]
)

# --- Convertir log-variance en variance ---
variance <- exp(logvar)

# --- Calculer le facteur multiplicatif par rapport à la référence ---
factor_mult <- variance / variance[1]

# --- Créer un tableau propre ---
df_variance <- data.frame(
  Design = levels_design,
  log_variance = round(logvar, 3),
  variance = round(variance, 4),
  factor_vs_reference = round(factor_mult, 2)
)

# --- Afficher le tableau ---
print(df_variance)




# ========================== 17. Random effects structure ===================


# Multi-level models
m_multi  <- rma.mv(seq_rate, seq_rate_vi, random = list(~1 | id_article/id_expérimentation), struct = "UN", data = FULL_sun)
m_multi1 <- rma.mv(seq_rate, seq_rate_vi, random = list(~1 | id_article, ~1 | id_expérimentation), struct = "UN", data = FULL_sun)

# Null models for variance partitioning
m_within_null  <- rma.mv(seq_rate, seq_rate_vi, random = list(~1 | id_article, ~1 | id_expérimentation), sigma2 = c(0, NA), data = FULL_sun)
m_between_null <- rma.mv(seq_rate, seq_rate_vi, random = list(~1 | id_article, ~1 | id_expérimentation), sigma2 = c(NA, 0), data = FULL_sun)
m_both_null    <- rma.mv(seq_rate, seq_rate_vi, random = list(~1 | id_article, ~1 | id_expérimentation), sigma2 = c(0, 0), data = FULL_sun)

# Compare models
#aov_scale     <- anova(m_multi, m_scale)
aov_within    <- anova(m_multi, m_within_null)
aov_between   <- anova(m_multi, m_between_null)
aov_bothnull  <- anova(m_multi, m_both_null)

# Summarize results
aov_table <- rbind(
  c(df = aov_between$p.f, aov_between$fit.stats.f[c(3:4, 1)], LRT = NA, p = NA),
  c(df = aov_within$p.r, aov_within$fit.stats.r[c(3:4, 1)], LRT = aov_within$LRT, p = aov_within$pval),
  c(df = aov_between$p.r, aov_between$fit.stats.r[c(3:4, 1)], LRT = aov_between$LRT, p = aov_between$pval),
  c(df = aov_bothnull$p.r, aov_bothnull$fit.stats.r[c(3:4, 1)], LRT = aov_bothnull$LRT, p = aov_bothnull$pval)
)
rownames(aov_table) <- c("Three-level model", "Within-studies variance constrained",
                         "Between-studies variance constrained", "Both variance components constrained")
aov_table


# 1. Three-level model (full model)
#    - AIC = 976.14, BIC = 987.33
#    - Baseline model with both within- and between-study variance freely estimated.

# 2. Within-studies variance constrained
#    - AIC = 1043.83 → much worse fit (increase of ~67 points)
#    - LRT = 69.69, p << 0.001 → constraining within-study variance is highly significant
#      ⇒ within-study variability is a meaningful contributor to overall heterogeneity.

# 3. Between-studies variance constrained
#    - AIC = 1043.39 → similar deterioration in fit
#    - LRT = 69.26, p << 0.001 → between-study variance is also a significant contributor.

# 4. Both variance components constrained
#    - AIC = 1538.81 → massive deterioration
#    - LRT = 566.68, p << 0.001 → clearly, both variance components are essential to model heterogeneity.

# ========================== Practical Takeaways =================
# - Both within- and between-study heterogeneity are significant and should be retained.
# - Ignoring either variance component would result in severely biased inference.
# - The three-level random-effects model is fully justified for these data.


# ========================== 18. Summary of three-level model ===================

summary(m_multi1)

# Multivariate meta-analysis (k = 303; method REML)
# - logLik = -467.93 ; AIC = 941.87 ; BIC = 952.99
# - Model includes two random effects:
#     • Between-articles variance (σ² = 0.2638, SD = 0.51, 44 levels)
#     • Within-articles / between-experiments variance (σ² = 0.0607, SD = 0.25, 303 levels)

# Heterogeneity test:
# - Q(302) = 852.89, p < 0.0001 → strong evidence of heterogeneity.

# Overall effect size (original scale):
# - Estimate = 0.2613 ± 0.1054 (SE), z = 2.48, p = 0.013
# - 95% CI = [0.0548 ; 0.4679]
# - No back-transformation needed for this outcome.

# ========================== Practical Takeaways =================
# - Both variance components are important, justifying the three-level structure.
# - The overall effect is positive and statistically significant.
# - Effect magnitude ≈ 0.26, with CI spanning 0.05 to 0.47.


# ========================== 19. Heterogeneity partition ===================

# Goal: Quantify the proportion of total variance explained by each random effect
#       in a three-level meta-analysis.

# Step 1: Construct the "hat" matrix for weighted least squares
W <- diag(1 / m_multi$vi)               # inverse-variance weights
X <- model.matrix(m_multi)              # design matrix for fixed effects
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W  # projection matrix

# Step 2: Compute total proportion of variance attributable to random effects
total_var_pct <- 100 * sum(m_multi$sigma2) /
  (sum(m_multi$sigma2) + (m_multi$k - m_multi$p)/sum(diag(P)))

cat("✅ Total proportion of variance due to random effects:", round(total_var_pct, 1), "%\n")

# Step 3: Compute proportion of variance by level
level_var_pct <- 100 * m_multi$sigma2 /
  (sum(m_multi$sigma2) + (m_multi$k - m_multi$p)/sum(diag(P)))

names(level_var_pct) <- c("Between-article", "Within-article")  # adjust names to your model
cat("Proportion of variance by level:\n")
print(round(level_var_pct, 1))




# ========================== 19. Moderator analysis: Language ===================
unique(FULL_sun$language)

FULL_sun$language <- ifelse(FULL_sun$language == "English\n\nSpanish",
                              "Spanish",
                              FULL_sun$language)

# Model with 'language' as moderator
m_lang <- rma.mv(seq_rate, seq_rate_vi,
                 mods = ~ language,
                 random = list(~1 | id_article/id_expérimentation),
                 struct = "UN",
                 data = FULL_sun)

# Compare model with and without moderator
aov_lang <- anova(m_multi, m_lang,refit=TRUE)

# Summarize results
summary(m_lang)
print(aov_lang)

# ========================== Interpretation ===================
# Intercept (English): 0.1926 [95% CI: -0.0573, 0.4425], p = 0.1309
#   → Overall effect size for English-language studies is positive but not statistically significant.
#
# Language (Spanish): 0.2228 [95% CI: -0.2260, 0.6715], p = 0.3305
#   → Effect sizes from Spanish-language studies are not significantly different
#     from those in English-language studies.
#
# Test of moderators: QM(1) = 0.947, p = 0.3305
#   → 'language' does not explain a significant share of heterogeneity.
#
# Model comparison (likelihood-ratio test): LRT = 0.974, p = 0.324
#   → Adding 'language' as moderator does not improve model fit (AIC/BIC unchanged).
#
# ========================== Practical takeaway =================
# There is no evidence that effect sizes differ systematically
# between English and Spanish publications once within- and between-study
# variance components are accounted for.

exclude <- c(128, 130, 198, 200, 201)
FULL_sun2_clean <- FULL_sun[-exclude, ]

dim(FULL_sun2_clean)
length(unique(FULL_sun2_clean$id_article))

