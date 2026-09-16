###############################################################################
# Project   : MACCA Database – Carbon Stock Analysis in Agroforestry Systems
# Author    : Damien Beillouin
# Date      : 30 July 2024
# Update    : September 2025
# Purpose   : Clean and prepare SOC stock data for LMP experiments comparing
#             Full Sun vs. Agroforestry systems.
# Notes     : FAIR principles: documented, reproducible, structured, transparent
# Focus on relative ratios!
###############################################################################

library(dplyr)
library(ggplot2)
library(ggExtra)
library(metafor)
library(forcats)
library(ggpubr)
library(vcd)

# ========================== Part I: Data Loading =========================
# Load the MACCA dataset
data_path <- "~/Documents/MACCA/Data_for_analysis.csv"
FULL_sunRR <- read.csv(data_path) %>%
  filter(experiment_type == "LMP",         # Keep LMP experiments only
         control_type == "Full sun",       # Full sun control
         Keep_ratio == "YES")              # Keep valid ratios

# Impute missing 'time_since_conversion' with the median
med_conversion <- median(FULL_sunRR$time_since_conversion, na.rm = TRUE)
FULL_sunRR$time_since_conversion[is.na(FULL_sunRR$time_since_conversion)] <- med_conversion

# ========================== Part II: Main Crops =========================
# Merge similar crop types
FULL_sunRR <- FULL_sunRR %>%
  mutate(main_culture2 = as.character(main_culture),
         main_culture2 = case_when(
           grepl("Cocoa|Coffee|cocoa", main_culture2) ~ "Cocoa/Coffee",
           grepl("Beans|Maize", main_culture2)       ~ "Beans/Maize",
           grepl("Banana", main_culture2)            ~ "Others",
           TRUE                                      ~ main_culture2
         ),
         main_culture = fct_lump_min(main_culture, 11))

# ========================== Part III: Agroforestry Systems =========================
FULL_sunRR <- FULL_sunRR %>%
  mutate(NEW_treatment_type = fct_lump_min(NEW_treatment_type, 10),
         NEW_treatment_type2 = as.character(NEW_treatment_type),
         NEW_treatment_type2 = ifelse(NEW_treatment_type2 %in% c("Alley cropping", "Hedgerow"),
                                      "Alley/Hedgerow", NEW_treatment_type2))

# ========================== Part IV: Nitrogen-Fixing Species =========================
FULL_sunRR <- FULL_sunRR %>%
  mutate(diffN2 = as.numeric(treatment_N2_fixing) - as.numeric(control_N2_fixing),
         diffN2Cat = case_when(
           diffN2 <= 0 ~ 0,
           diffN2 == 1 ~ 1,
           diffN2 > 1  ~ 2
         ))

# ========================== Part V: Species and Woody Plant Counts =========================
FULL_sunRR <- FULL_sunRR %>%
  mutate(diff_Woody = treatment_N_Woody - control_N_Woody,
         diff_Woody_cat = case_when(
           diff_Woody >= 0 & diff_Woody <= 2  ~ "1",
           diff_Woody > 2 & diff_Woody <= 99 ~ "2+",
           TRUE                              ~ "0"
         ),
         diff_species = treatment_NB_species - as.numeric(control_NB_species),
         diff_species_class = cut(diff_species, breaks = c(-5, 1.1, 25),
                                         labels = c("A:<0/+1", "C+2+"), right = FALSE),
                diff_species_class = factor(ifelse(is.na(diff_species_class), "C+2+", as.character(diff_species_class))))

# ========================== Part VI: Soil and Climate =========================
FULL_sunRR <- FULL_sunRR %>%
  mutate(soil_type = fct_lump_min(soil_type, 10),
         FAO_Climate = fct_lump_min(FAO_Climate, 10))

# ========================== Part VII: Experimental Design =========================
FULL_sunRR <- FULL_sunRR %>%
  mutate(Grouped_Design = case_when(
    is.na(Design) ~ "Not specified",
    Design %in% c("RCBD", "RCT") ~ "Randomized Designs",
    Design %in% c("BACI", "BA") ~ "Before-After Designs",
    Design == "CI" ~ "Control-Impact Designs",
    TRUE ~ "Other"
  ))

# ========================== Part VIII: Quality Score =========================
FULL_sunRR <- FULL_sunRR %>%
  mutate(quality_score = rowSums(cbind(
    (!is.na(treatment_soc_sd_T_ha)) * 1,
    (!is.na(history_C) & history_C != "unknown") * 1,
    (document_type == "article") * 1,
    is.na(Confunding_factors) * 1,
    Grouped_Design %in% c("Before-After Designs", "Randomized Designs") * 1
  )))

# ========================== Part IX: Impute Missing SDs =========================
impute_sd <- function(df, mean_var, sd_var, group_var) {
  TAB <- df %>% group_by(!!sym(group_var)) %>%
    summarise(Mean_m = mean(!!sym(mean_var), na.rm = TRUE),
              SD_m = mean(!!sym(sd_var), na.rm = TRUE),
              .groups = "drop")
  TAB$CV_m <- TAB$SD_m / TAB$Mean_m
  df <- df %>%
    left_join(TAB, by = group_var) %>%
    mutate(!!sd_var := ifelse(is.na(!!sym(sd_var)), CV_m * !!sym(mean_var), !!sym(sd_var))) %>%
    dplyr::select(-Mean_m, -SD_m, -CV_m)
  return(df)
}

FULL_sunRR <- FULL_sunRR %>%
  impute_sd("treatment_soc_mean_T_ha", "treatment_soc_sd_T_ha", "Grouped_Design") %>%
  impute_sd("control_soc_mean_T_ha", "control_soc_sd_T_ha", "Grouped_Design")

# ========================== Part X: Compute Effect Sizes =========================
lnRR <- escalc(measure = "ROM",
               n1i = FULL_sunRR$treatment_replicate_nb,
               n2i = FULL_sunRR$control_replicate_nb,
               m1i = FULL_sunRR$treatment_soc_mean_T_ha,
               m2i = FULL_sunRR$control_soc_mean_T_ha,
               sd1i = FULL_sunRR$treatment_soc_sd_T_ha,
               sd2i = FULL_sunRR$control_soc_sd_T_ha)

FULL_sunRR <- cbind(FULL_sunRR, lnRR)

# ========================== Part XI: Density =========================
transform_density <- function(value) {
  if (is.na(value)) return(NA)
  value <- gsub(",", ".", value)
  if (grepl("-", value)) {
    bounds <- as.numeric(unlist(strsplit(value, "-")))
    return(mean(bounds))
  } else if (grepl("<|>", value)) {
    return(as.numeric(gsub("<|>", "", value)))
  } else {
    return(as.numeric(value))
  }
}

FULL_sunRR$density_numeric <- sapply(FULL_sunRR$density, transform_density)

# ========================== Part XII: Reclass History =========================
FULL_sunRR <- FULL_sunRR %>%
  mutate(Hist = paste(history_C, history_T),
         History_reclass = case_when(
           Hist %in% c("cropland forest", "forest cropland", "NA NA",
                       "unknown pasture", "unknown unknown") ~ "Unknown/mixed",
           Hist == "cropland cropland" ~ "Cropland",
           Hist %in% c("forest forest", "pasture pasture") ~ "Forest",
           Hist %in% c("grassland grassland", "pasture pasture") ~ "Grassland",
           TRUE ~ "Unknown/mixed"
         ),
         Hist2 = paste(History_reclass, control_Land_Use),
         Hist2 = case_when(
           Hist2 %in% c("Forest Grassland", "Grassland Grassland") ~ "Forest/Grassland Grassland",
           TRUE ~ Hist2
         ))

# ========================== Part XIII: Visual Checks =========================
scatter_plot_color <- ggplot(FULL_sunRR) +
  geom_point(aes(y = yi, x = vi, fill = Grouped_Design, color = Grouped_Design, size = 1/vi),
             shape = 21, color = "black") +
  labs(x = "Variance", y = "Sequestration rate (Mg.ha^-1.yr^-1)") +
  theme_pubr(base_size = 14) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype = 2, col = 'grey70') +
  geom_hline(yintercept = 0, linetype = 2, col = 'grey70') +
  scale_size(range = c(1, 6))

ggMarginal(scatter_plot_color, type = "boxplot", groupFill = TRUE)

# ========================== Part XIV: Meta-analysis Models =========================
# Goal: Test whether the variance of sequestration rates depends on
#       time since conversion

# --- 1. Scale model: variance as a function of time since conversion ---
m_scale <- rma(
  yi = yi,
  vi = vi,
  scale = ~ time_since_conversion,
  data = FULL_sunRR
)

# --- 2. Null model: constant variance ---
m_null <- rma(
  yi = yi,
  vi = vi,
  scale = ~1,
  data = FULL_sunRR
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
# Key results from anova_result:
# - LRT = 11.01, p = 0.0009 → variance of sequestration rates significantly depends on time since conversion
# - AIC improves from 21.64 → 12.63 (~9 points) → full model provides a better fit
# - Residual heterogeneity (QE = 2639.57) remains very high → additional sources of variability exist
#
# Practical takeaway:
# Time since conversion exerts a real, statistically robust effect on the variance of sequestration rates.
# Predictions and meta-analytic models should explicitly account for this scale effect.


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
    title = "Predicted Variance of relative Ratios over Time Since Conversion",
    x = "Time Since Conversion (years)",
    y = "Predicted Variance"
  ) +
  theme_minimal(base_size = 14)

# Add marginal density for time_since_conversion
ggMarginal(p, type = "density", margins = "x", fill="grey80", alpha=0.5)



# Extraire la covariance des coefficients de scale en tant que matrice
vcov_alpha_mat <- as.matrix(vcov(m_scale)$alpha)

# Vérifier les noms
rownames(vcov_alpha_mat)
colnames(vcov_alpha_mat)

# Calculer l'erreur standard du prédicteur linéaire
se_eta <- sqrt(
  vcov_alpha_mat["intrcpt","intrcpt"] +
    (time_seq^2) * vcov_alpha_mat["time_since_conversion","time_since_conversion"] +
    2 * time_seq * vcov_alpha_mat["intrcpt","time_since_conversion"]
)


eta <- beta0 + beta1 * time_seq
# --- 4. 95% CI on log scale ---
ci_upper <- eta + 1.96 * se_eta
ci_lower <- eta - 1.96 * se_eta

# --- 5. Transform back to variance scale ---
df_var$lower_CI <- exp(ci_lower)
df_var$upper_CI <- exp(ci_upper)


p <- ggplot(df_var, aes(x = time_since_conversion, y = predicted_variance)) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "lightblue", alpha = 0.3) +  # CI band
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(data = data.frame(time_since_conversion = c(1,10),
                               predicted_variance = exp(beta0 + beta1 * c(1,10))),
             aes(x = time_since_conversion, y = predicted_variance),
             color = "red", size = 3) +
  labs(
    title = "Predicted Variance of relative Ratios over Time Since Conversion",
    x = "Time Since Conversion (years)",
    y = "Predicted Variance"
  ) +
  theme_minimal(base_size = 14)

# ========================== 16. Variance structure by design =================
# Goal: Show how predicted variance changes with time since conversion,
#       separately for each Grouped_Design

# --- 1. Scale model: variance as a function of time since conversion ---
m_scale <- rma(
  yi = yi,
  vi = vi,
  scale = ~ Grouped_Design*time_since_conversion,
  data = FULL_sunRR
)

m_scale2 <- rma(
  yi = yi,
  vi = vi,
  scale = ~ Grouped_Design+time_since_conversion,
  data = FULL_sunRR
)
anova_result <- anova(m_scale, m_scale2)
print(anova_result)
AIC(m_scale2)

# --- 2. Null model: constant variance ---
m_null <- rma(
  yi = yi,
  vi = vi,
  scale = ~1,
  data = FULL_sunRR
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
m_multi  <- rma.mv(yi, vi, random = list(~1 | id_article/id_expérimentation), struct = "UN", data = FULL_sunRR)
m_multi1 <- rma.mv(yi, vi, random = list(~1 | id_article, ~1 | id_expérimentation), struct = "UN", data = FULL_sunRR)

# Null models for variance partitioning
m_within_null  <- rma.mv(yi, vi, random = list(~1 | id_article, ~1 | id_expérimentation), sigma2 = c(0, NA), data = FULL_sunRR)
m_between_null <- rma.mv(yi, vi, random = list(~1 | id_article, ~1 | id_expérimentation), sigma2 = c(NA, 0), data = FULL_sunRR)
m_both_null    <- rma.mv(yi, vi, random = list(~1 | id_article, ~1 | id_expérimentation), sigma2 = c(0, 0), data = FULL_sunRR)

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
#    - AIC = -61.52, BIC = -49.99
#    - Baseline model with both within- and between-study variance freely estimated.

# 2. Within-studies variance constrained
#    - AIC = 20.75 → much worse fit (increase of ~82 points)
#    - LRT = 84.27, p = 4.32e-20 → constraining within-study variance is highly significant
#      ⇒ within-study variability is a meaningful contributor to overall heterogeneity.

# 3. Between-studies variance constrained
#    - AIC = 341.28 → dramatic deterioration in fit
#    - LRT = 404.80, p = 4.97e-90 → between-study variance is also a highly significant contributor.

# 4. Both variance components constrained
#    - AIC = 1908.10 → massive deterioration
#    - LRT = 1973.62, p ≈ 0 → clearly, both variance components are essential to model heterogeneity.

# ========================== Practical Takeaways =================
# - Both within- and between-study heterogeneity are significant and should be retained.
# - Ignoring either variance component would lead to severely biased inference.
# - The three-level random-effects model is fully justified for these data.


# ========================== 18. Summary of three-level model ===================

summary(m_multi1)

# Multivariate meta-analysis (k = 346; method REML)
# - logLik = 34.77 ; AIC = -63.53 ; BIC = -52.00
# - Model includes two random effects:
#     • Between-articles variance (σ² = 0.0224, SD = 0.15, 54 levels)
#     • Within-articles / between-experiments variance (σ² = 0.0111, SD = 0.11, 332 levels)

# Heterogeneity test:
# - Q(345) = 2641.85, p < 0.0001 → very strong evidence of heterogeneity.

# Overall effect size (log scale):
# - Estimate = 0.0686 ± 0.0239 (SE), z = 2.87, p = 0.0041
# - 95% CI = [0.0217 ; 0.1155]

# -------------------------------------------------------------------
# Back-transformation (log scale → natural ratio scale)
est_bt   <- exp(coef(m_multi1))   # median effect (ratio)
ci.lb_bt <- exp(0.0217)           # lower CI (ratio)
ci.ub_bt <- exp(0.1155)           # upper CI (ratio)

# Back-transformation (ratio → % change, median)
est_pct   <- 100 * (est_bt   - 1)
ci.lb_pct <- 100 * (ci.lb_bt - 1)
ci.ub_pct <- 100 * (ci.ub_bt - 1)

# Mean of the lognormal distribution (true expectation on original scale)
mu     <- coef(m_multi1)        # mean on log scale
sigma2 <- vcov(m_multi1)        # variance of log-scale estimate
mean_bt  <- exp(mu + 0.5 * sigma2)       # mean effect (ratio)
mean_pct <- 100 * (mean_bt - 1)          # mean effect (%)

# Display results
results <- data.frame(
  scale = c("Median (exp(mu))", "Mean (exp(mu+0.5σ²))"),
  ratio = c(est_bt, mean_bt),
  percent = c(est_pct, mean_pct),
  CI95.lb = c(ci.lb_pct, NA),   # CI only for median (from model)
  CI95.ub = c(ci.ub_pct, NA)
)
results

# Expected output (approx):
#   - Median ratio ≈ 1.071  → +6.9% (95% CI: +2.2% to +12.2%)
#   - Mean ratio   ≈ 1.072  → +7.0% (no CI, depends on lognormal variance)

# ========================== Practical Takeaways =================
# - Both article-level and experiment-level variance components
#   contribute significantly to heterogeneity → justifies 3-level model.
# - Back-transformation shows a median increase of ~7% (CI: +2% to +12%).
# - The lognormal mean is almost identical here (~7%), confirming robustness.



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
unique(FULL_sunRR$language)
unique(FULL_sunRR$language)

FULL_sunRR$language <- ifelse(FULL_sunRR$language == "English\n\nSpanish",
                              "Spanish",
                              FULL_sunRR$language)

# Model with 'language' as moderator
m_lang <- rma.mv(yi, vi,
                 mods = ~ language,
                 random = list(~1 | id_article/id_expérimentation),
                 struct = "UN",
                 data = FULL_sunRR)

# Compare model with and without moderator
aov_lang <- anova(m_multi, m_lang,refit=TRUE)

# Summarize results
summary(m_lang)
print(aov_lang)

# ========================== Interpretation ===================
# - Intercept (English): 0.0835 [95% CI: 0.0223, 0.1447], p = 0.0075 →
#   overall effect size is significantly positive for English-language studies.
#
# - Language (Spanish): -0.0365 [95% CI: -0.1322, 0.0592], p = 0.4546 →
#   effect sizes from Spanish-language studies are not significantly different
#   from those in English-language studies.
#
# - Test of moderators: QM(1) = 0.56, p = 0.45 →
#   'language' does not explain a significant share of heterogeneity.
#
# - Model comparison (likelihood-ratio test): LRT = 0.58, p = 0.45 →
#   adding 'language' as moderator does not improve model fit (AIC, BIC unchanged).
#
# ========================== Practical takeaway =================
# - There is no evidence that the effect sizes differ systematically
#   between English and Spanish publications once within- and between-study
#   variance components are accounted for.

# Indices des points 
Exclude <- c(107, 121, 155, 259, 262)

# Dataset sans les outliers
FULL_sunRR_clean <- FULL_sunRR[-Exclude, ]

 dim(FULL_sunRR_clean)
length(unique(FULL_sunRR_clean$Authors))
