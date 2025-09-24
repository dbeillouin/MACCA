##############################################
# Partial Dependence and XGBoost Analysis
# Author: Damien Beillouin
# Purpose: Analyze effect of initial SOC on SEQ rate
#          across soil depth groups using XGBoost and PDPs
##############################################

# -----------------------------
# 0. Load libraries
# -----------------------------
library(pdp)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mgcv)
library(gratia)

# -----------------------------
# 1. Filter dataset if needed
# -----------------------------
unique(FULL_sun2_clean$NEW_treatment_type2)
FULL_sun2_sub <- FULL_sun2_clean 

# -----------------------------
# 2. Build design matrix for XGBoost
# -----------------------------
X_full <- model.matrix(
  seq_rate ~ NEW_treatment_type2 + diff_species_class + History_reclass +
    control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
    MEAN_depth + time_since_conversion + Grouped_Design,
  data = FULL_sun2_sub
)[, -1]


y <- FULL_sun2_sub$seq_rate

params <- list(
  objective = "reg:squarederror",
  max_depth = 4,
  eta = 0.05,
  subsample = 0.9,
  colsample_bytree = 0.8
)

# -----------------------------
# 4. Cross-validation to find optimal rounds
# -----------------------------
cv_results <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 500,
  nfold = 5,
  verbose = 1,
  early_stopping_rounds = 20,
  metrics = "rmse"
)

best_nrounds <- cv_results$best_iteration

# -----------------------------
# 5. Train final XGBoost model
# -----------------------------
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_nrounds
)

best_nrounds <- cv_results$best_iteration
best_params <- params  # contient max_depth, eta, colsample_bytree, subsample...

xgb_model <- xgboost(
  data = X_full,
  label = y,
  nrounds = best_nrounds,
  objective = best_params$objective,
  max_depth = best_params$max_depth,
  eta = best_params$eta,
  subsample = best_params$subsample,
  colsample_bytree = best_params$colsample_bytree,
  verbose = 0
)


# -----------------------------
# 3. Define soil depth groups (10-cm bins)
# -----------------------------
FULL_sun2_sub$depth_group <- cut(
  FULL_sun2_sub$MEAN_depth,
  breaks = c(0, 15, 30, 45, 55, 75),
  include.lowest = TRUE,
  right = FALSE
)

# Check for NAs
FULL_sun2_sub[is.na(FULL_sun2_sub$depth_group), c("MEAN_depth", "depth_group")]

# -----------------------------
# 4. Compute Partial Dependence (ICE)
# -----------------------------
range_vals <- quantile(FULL_sun2_sub$temperature, probs = c(0.05, 0.99))

pdp_initial_stock <- pdp::partial(
  object = xgb_model,
  pred.var = "temperature",
  train = X_full,
  ice = TRUE,
  center = FALSE,
  pred.grid = data.frame(temperature = seq(range_vals[1], range_vals[2], length.out = 50))
)


# -----------------------------
# 5. Add depth group and other covariates to PDP
# -----------------------------
n_obs <- nrow(FULL_sun2_sub)
n_vals <- length(seq(range_vals[1], range_vals[2], length.out = 50))

pdp_initial_stock$depth_group <- rep(FULL_sun2_sub$depth_group, each = n_vals)

vars_to_add <- c("depth_group", "Grouped_Design", "NEW_treatment_type2",
                 "diff_species_class", "History_reclass", "main_culture2")

for (v in vars_to_add) {
  pdp_initial_stock[[v]] <- rep(FULL_sun2_sub[[v]], each = n_vals)
}

# -----------------------------
# 6. Define colors for depth groups
# -----------------------------
depth_colors <- c(
  "[0,15)"   = "#FFCC66",
  "[15,30)"  = "#FF9966",
  "[30,45)"  = "#FF6666",
  "[45,55)"  = "#CC3333",
  "[55,75]"  = "#990000"
)

# -----------------------------
# 7. PDP plot per depth
# -----------------------------
ggplot(pdp_initial_stock, aes(x = temperature, y = yhat, 
                              color = depth_group, fill = depth_group)) +
  
  # Horizontal reference line at y = 0
  geom_hline(aes(yintercept = 0), linetype = 2) +
  
  # GAM smoothing with confidence interval
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 6),
    se = TRUE,
    alpha = 0.2,
    size = 1.2
  ) +
  
  # Overlay observed points
  geom_point(
    data = FULL_sun2_sub,
    aes(x = temperature, y = seq_rate, color = depth_group),
    inherit.aes = FALSE,
    alpha = 0.3,
    size = scales::rescale(1 / FULL_sun2_sub$seq_rate_sd, to = c(1, 5))
  ) +
  scale_size_continuous(range = c(1, 3)) + # ajuster la plage des tailles
  
  # Color and fill scales
  scale_color_manual(values = depth_colors) +
  scale_fill_manual(values = depth_colors) +
  
  # Labels
  labs(
    x = "Initial SOC stock (t/ha)",
    y = "Partial dependence",
    color = "Soil depth group",
    fill = "Soil depth group"
  ) +
  
  # Theme adjustments
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80", size = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14)
  ) +
  
  # Rug plot for SOC distribution
  geom_rug(
    data = FULL_sun2_sub,
    aes(x = temperature, color = depth_group),
    inherit.aes = FALSE,
    sides = "b",
    length = unit(0.02, "npc"),
    alpha = 0.5
  )

# -----------------------------
# 8. Compute pseudo-R² by depth group
# -----------------------------
calc_r2_group <- function(group_name, data_full, model, X_full) {
  idx <- data_full$depth_group == group_name
  y_obs <- data_full$seq_rate[idx]
  y_pred <- predict(model, newdata = X_full[idx, ])
  cor(y_obs, y_pred)^2
}

groups <- unique(FULL_sun2_sub$depth_group)
r2_table <- sapply(groups, calc_r2_group, data_full = FULL_sun2_sub, model = xgb_model, X_full = X_full)

# -----------------------------
# 9. Add GAM-smoothed predictions to dataset
# -----------------------------
aic_results <- pdp_initial_stock %>%
  group_by(depth_group) %>%
  group_map(~ {
    tibble(
      k = 2:10,
      AIC = map_dbl(2:10, function(kk) {
        m <- gam(yhat ~ s(temperature, bs = "cs", k = kk), data = .x)
        AIC(m)
      }),
      depth_group = as.character(.y$depth_group)   # récupération correcte du nom du groupe
    )
  }) %>%
  bind_rows()

print(aic_results)

ggplot(aic_results, aes(x = k, y = AIC, color = depth_group)) +
  geom_line() +
  geom_point() +
  theme_minimal()

gam_models <- pdp_initial_stock %>%
  group_by(depth_group) %>%
  group_map(~ gam(yhat ~ s(temperature, bs = "cs", k = 4), data = .x), .keep = TRUE)

names(gam_models) <- levels(pdp_initial_stock$depth_group)

smooth_deriv <- lapply(seq_along(gam_models), function(i) {
  gm <- gam_models[[i]]
  depth <- names(gam_models)[i]
  d <- derivatives(gm, term = "s(temperature)", interval = "confidence", n = 100)
  d$depth_group <- depth
  d
}) %>% bind_rows()

smooth_deriv <- smooth_deriv %>%
  mutate(
    slope_signif = (.lower_ci > 0) | (.upper_ci < 0)
  )

smooth_deriv <- smooth_deriv %>%
  group_by(depth_group) %>%
  arrange(temperature) %>%
  mutate(
    slope_sign = sign(.derivative),
    sign_change = slope_sign != lag(slope_sign, default = first(slope_sign))
  ) %>%
  ungroup()

summary_deriv <- smooth_deriv %>%
  group_by(depth_group) %>%
  summarize(
    n_points = n(),  # nombre total d'observations dans le groupe
    mean_slope = mean(.derivative),
    overall_sign = case_when(
      mean_slope > 0 ~ "positive",
      mean_slope < 0 ~ "negative",
      TRUE ~ "neutral"
    ),
    first_significant_time = min(temperature[slope_signif], na.rm = TRUE),
    inflection_times = paste(round(temperature[sign_change], 1), collapse = ", ")
  )

ggplot(smooth_deriv, aes(x = temperature, y = .derivative, color = depth_group)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = depth_group), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(data = subset(smooth_deriv, sign_change == TRUE),
             aes(x = temperature, y = .derivative), color = "black", size = 3) +
  labs(y = "Slope (first derivative)", x = "Time since conversion") +
  theme_minimal()

# -----------------------------
# 10. Compute metrics by depth
# -----------------------------
p <- ggplot(pdp_initial_stock, aes(x = temperature, y = yhat, color = depth_group)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5), se = TRUE)
smooth_data <- ggplot_build(p)$data[[1]]

group_levels <- levels(pdp_initial_stock$depth_group)
smooth_data$depth_group <- factor(group_levels[smooth_data$group], levels = group_levels)

FULL_sun2_sub <- FULL_sun2_sub %>%
  rowwise() %>%
  mutate(
    y_gam = {
      smooth_g <- smooth_data %>%
        filter(depth_group == depth_group) %>%
        arrange(x) %>%
        distinct(x, .keep_all = TRUE)
      approx(x = smooth_g$x, y = smooth_g$y, xout = control_soc_mean_T_ha, rule = 2)$y
    }
  ) %>%
  ungroup()

metrics_by_depth <- FULL_sun2_sub %>%
  group_by(depth_group) %>%
  summarise(
    n_points = n(),
    pseudo_R2 = cor(seq_rate, y_gam)^2,
    RMSE = sqrt(mean((seq_rate - y_gam)^2, na.rm = TRUE)),
    MAE  = mean(abs(seq_rate - y_gam), na.rm = TRUE),
    Spearman = cor(seq_rate, y_gam, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )

metrics_by_depth$RGBOOST <- r2_table

# 4. GAM with depth interaction
gam_depth <- gam(yhat ~ s(temperature, by = depth_group, bs = "cs", k = 6) +
                   depth_group,
                 data = pdp_initial_stock)

summary(gam_depth)  # Check significance of terms

# 5. GAM without interaction (same smooth for all depths)
gam_noint <- gam(yhat ~ s(temperature, bs = "cs", k = 6) + depth_group,
                 data = pdp_initial_stock)

# Compare models to test if depth-specific smooths improve fit
anova(gam_noint, gam_depth, test = "Chisq")


# Moyenne marginale des prédictions pour chaque depth_group
emm <- emmeans(gam_depth, ~ depth_group)

# Comparaisons par paires
pairs(emm)

# Extraire les paires depuis emmeans
pair_df <- as.data.frame(pairs(emm)) %>%
  mutate(Significant = ifelse(p.value < 0.05, "Yes", "No")) %>%
  dplyr::select(contrast, estimate, Significant)

# Affichage
print(pair_df)

summary_pairs <- pair_df %>%
  filter(Significant == "Yes") %>%
  dplyr::rowwise() %>%
  mutate(
    group1 = strsplit(contrast, " - ")[[1]][1],
    group2 = strsplit(contrast, " - ")[[1]][2]
  ) %>%
  dplyr::select(group1, group2) %>%
  group_by(group1) %>%
  summarise(Significant_vs = paste(group2, collapse = ", "), .groups = "drop")

print(summary_pairs)


# 6. Predict at 40 T ha with confidence intervals
library(dplyr)

# Target SOC values
target_soc <- c(20, 28)
depths <- levels(pdp_initial_stock$depth_group)

# 1. Predictions from GAM
pred_df <- lapply(target_soc, function(soc) {
  newdata <- data.frame(
    temperature = soc,
    depth_group = depths
  )
  pred <- predict(gam_depth, newdata = newdata, se.fit = TRUE)
  data.frame(
    target_soc = soc,
    depth_group = depths,
    fit = pred$fit,
    lower = pred$fit - 1.96 * pred$se.fit,
    upper = pred$fit + 1.96 * pred$se.fit
  )
}) %>% bind_rows()


# pred_df = GAM predictions
pred_df




# -----------------------------
# 11. Spearman correlation test by depth
# -----------------------------

# Average PDP per depth/SOC
pdp_unique <- pdp_initial_stock %>%
  group_by(depth_group, temperature) %>%
  summarise(yhat = mean(yhat), .groups = "drop")

# Interpolate for each observation
FULL_sun2_sub <- FULL_sun2_sub %>%
  rowwise() %>%
  mutate(
    yhat_num = approx(
      x = pdp_unique$temperature[pdp_unique$depth_group == depth_group],
      y = pdp_unique$yhat[pdp_unique$depth_group == depth_group],
      xout = temperature,
      rule = 2
    )$y
  ) %>%
  ungroup() %>%
  mutate(yhat_num = as.numeric(yhat_num))

# Compute Spearman correlation & trend table
summary_table <- FULL_sun2_sub %>%
  group_by(depth_group) %>%
  summarise(
    n_points = n(),
    rho = cor(yhat_num, temperature, method = "spearman"),
    p_val = cor.test(yhat_num, temperature, method = "spearman")$p.value,
    trend = case_when(
      rho > 0 ~ "increase",
      rho < 0 ~ "decrease",
      TRUE    ~ "none"
    ),
    reliability = case_when(
      n_points < 5 ~ "low",
      p_val > 0.05 ~ "moderate",
      TRUE ~ "high"
    ),
    .groups = "drop"
  ) %>%
  arrange(depth_group)

print(summary_table)


##############################################
# Partial Dependence and XGBoost Analysis - RR
# Purpose: Analyze effect of initial SOC on RR
#          across soil depth groups using XGBoost and PDPs
##############################################
##############################################

# -----------------------------
# 0. Load libraries
# -----------------------------
library(pdp)
library(ggplot2)
library(dplyr)
library(tidyr)
library(xgboost)

# -----------------------------
# 1. Filter dataset if needed
# -----------------------------
unique(FULL_sun2_clean$Grouped_Design)
FULL_sun2_sub <- FULL_sunRR #%>% filter(control_soc_mean_T_ha<150)
#%>%
#  filter(Grouped_Design == "Randomized Designs")


# -----------------------------
# 2. Build design matrix for XGBoost
# -----------------------------
X_full <- model.matrix(
  yi ~ NEW_treatment_type2 + diff_species_class + History_reclass +
    control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
    MEAN_depth + time_since_conversion + Grouped_Design,
  data = FULL_sun2_sub
)[, -1]

y <- FULL_sun2_sub$yi



params <- list(
  objective = "reg:squarederror",
  max_depth = 4,
  eta = 0.05,
  subsample = 0.9,
  colsample_bytree = 0.8
)
dtrain <- xgb.DMatrix(data = X_full, label = y)

# -----------------------------
# 4. Cross-validation to find optimal rounds
# -----------------------------
cv_results <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 500,
  nfold = 5,
  verbose = 1,
  early_stopping_rounds = 20,
  metrics = "rmse"
)

best_nrounds <- cv_results$best_iteration

# -----------------------------
# 5. Train final XGBoost model
# -----------------------------
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_nrounds
)


best_nrounds <- cv_results$best_iteration
best_params <- params  # contient max_depth, eta, colsample_bytree, subsample...

xgb_model <- xgboost(
  data = X_full,
  label = y,
  nrounds = best_nrounds,
  objective = best_params$objective,
  max_depth = best_params$max_depth,
  eta = best_params$eta,
  subsample = best_params$subsample,
  colsample_bytree = best_params$colsample_bytree,
  verbose = 0
)

# -----------------------------
# 3. Define soil depth groups (10-cm bins)
# -----------------------------
FULL_sun2_sub$depth_group <- cut(
  FULL_sun2_sub$MEAN_depth,
  breaks = c(0, 15, 30, 45, 55, 75),
  include.lowest = TRUE,
  right = FALSE
)

# Check for NAs
FULL_sun2_sub[is.na(FULL_sun2_sub$depth_group), c("MEAN_depth", "depth_group")]

# -----------------------------
# 4. Compute Partial Dependence (ICE)
# -----------------------------
range_vals <- quantile(FULL_sun2_sub$temperature, probs = c(0.05, 0.99))


pdp_initial_stock <- pdp::partial(
  object = xgb_model,
  pred.var = "temperature",
  train = X_full,
  ice = TRUE,
  center = FALSE,
  pred.grid = data.frame(temperature = seq(range_vals[1], range_vals[2], length.out = 50))
)

# -----------------------------
# 5. Add depth group and other covariates to PDP
# -----------------------------
n_obs <- nrow(FULL_sun2_sub)
n_vals <- length(seq(range_vals[1], range_vals[2], length.out = 50))

pdp_initial_stock$depth_group <- rep(FULL_sun2_sub$depth_group, each = n_vals)

vars_to_add <- c("depth_group", "Grouped_Design", "NEW_treatment_type2",
                 "diff_species_class", "History_reclass", "main_culture2")

for (v in vars_to_add) {
  pdp_initial_stock[[v]] <- rep(FULL_sun2_sub[[v]], each = n_vals)
}

# -----------------------------
# 6. Define colors for depth groups
# -----------------------------
depth_colors <- c(
  "[0,15)"   = "#FFCC66",
  "[15,30)"  = "#FF9966",
  "[30,45)"  = "#FF6666",
  "[45,55)"  = "#CC3333",
  "[55,75]"  = "#990000"
)

# -----------------------------
# 7. PDP plot per depth
# -----------------------------
ggplot(pdp_initial_stock, aes(x = temperature, y = exp(yhat), 
                              color = depth_group, fill = depth_group)) +
  
  # Horizontal reference line at y = 1 (log(0) -> exp(0) = 1)
  geom_hline(aes(yintercept = 1), linetype = 2) +
  
  # GAM smoothing with confidence interval
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs", k =7),
    se = TRUE,
    alpha = 0.2,
    size = 1.2
  ) +
  
  # Overlay observed points (back-transformed)
  geom_point(
    data = FULL_sun2_sub,
    aes(x = temperature, y = exp(yi), color = depth_group),
    inherit.aes = FALSE,
    alpha = 0.3,
    size = 1
  ) +
  
  # Color and fill scales
  scale_color_manual(values = depth_colors) +
  scale_fill_manual(values = depth_colors) +
  
  # Labels
  labs(
    x = "Initial SOC stock (t/ha)",
    y = "RR (back-transformed)",
    color = "Soil depth group",
    fill = "Soil depth group"
  ) +
  
  # Theme adjustments
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey80", size = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14)
  ) +
  
  # Rug plot for SOC distribution
  geom_rug(
    data = FULL_sun2_sub,
    aes(x = temperature, color = depth_group),
    inherit.aes = FALSE,
    sides = "b",
    length = unit(0.02, "npc"),
    alpha = 0.5
  )+
  ylim(0.5,2.5)

# -----------------------------
# 6. Compute pseudo-R² by depth group for RR
# -----------------------------
r2_table_RR <- sapply(groups, function(g) {
  idx <- FULL_sun2_sub$depth_group == g
  cor(FULL_sun2_sub$yi[idx], predict(xgb_model, newdata = X_full[idx, ]))^2
})

# -----------------------------
# 7. Add GAM-smoothed predictions to dataset
# -----------------------------
gam_models <- pdp_initial_stock %>%
  group_by(depth_group) %>%
  group_map(~ gam(yhat ~ s(temperature, bs = "cs", k = 6), data = .x), .keep = TRUE)

names(gam_models) <- levels(pdp_initial_stock$depth_group)

smooth_deriv <- lapply(seq_along(gam_models), function(i) {
  gm <- gam_models[[i]]
  depth <- names(gam_models)[i]
  d <- derivatives(gm, term = "s(temperature)", interval = "confidence", n = 100)
  d$depth_group <- depth
  d
}) %>% bind_rows()

smooth_deriv <- smooth_deriv %>%
  mutate(
    slope_signif = (.lower_ci > 0) | (.upper_ci < 0)
  )

smooth_deriv <- smooth_deriv %>%
  group_by(depth_group) %>%
  arrange(temperature) %>%
  mutate(
    slope_sign = sign(.derivative),
    sign_change = slope_sign != lag(slope_sign, default = first(slope_sign))
  ) %>%
  ungroup()

summary_deriv <- smooth_deriv %>%
  group_by(depth_group) %>%
  summarize(
    n_points = n(),  # nombre total d'observations dans le groupe
    mean_slope = mean(.derivative),
    overall_sign = case_when(
      mean_slope > 0 ~ "positive",
      mean_slope < 0 ~ "negative",
      TRUE ~ "neutral"
    ),
    first_significant_time = min(temperature[slope_signif], na.rm = TRUE),
    inflection_times = paste(round(temperature[sign_change], 1), collapse = ", ")
  )

ggplot(smooth_deriv, aes(x = temperature, y = .derivative, color = depth_group)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = depth_group), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(data = subset(smooth_deriv, sign_change == TRUE),
             aes(x = temperature, y = .derivative), color = "black", size = 3) +
  labs(y = "Slope (first derivative)", x = "Time since conversion") +
  theme_minimal()

# -----------------------------
# 10. Compute metrics by depth
# -----------------------------

p <- ggplot(pdp_initial_stock, aes(x = temperature, y = yhat, color = depth_group)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5), se = TRUE)
smooth_data <- ggplot_build(p)$data[[1]]

group_levels <- levels(pdp_initial_stock$depth_group)
smooth_data$depth_group <- factor(group_levels[smooth_data$group], levels = group_levels)

FULL_sun2_sub <- FULL_sun2_sub %>%
  rowwise() %>%
  mutate(
    y_gam = {
      smooth_g <- smooth_data %>%
        filter(depth_group == depth_group) %>%
        arrange(x) %>%
        distinct(x, .keep_all = TRUE)
      approx(x = smooth_g$x, y = smooth_g$y, xout = control_soc_mean_T_ha, rule = 2)$y
    }
  ) %>%
  ungroup()


metrics_by_depth <- FULL_sun2_sub %>%
  group_by(depth_group) %>%
  summarise(
    n_points = n(),
    pseudo_R2 = cor(yi, y_gam)^2,
    RMSE = sqrt(mean((yi - y_gam)^2, na.rm = TRUE)),
    MAE  = mean(abs(yi - y_gam), na.rm = TRUE),
    Spearman = cor(yi, y_gam, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )

metrics_by_depth$RGBOOST <- r2_table


# 4. GAM with depth interaction
gam_depth <- gam(yhat ~ s(temperature, by = depth_group, bs = "cs", k = 6) +
                   depth_group,
                 data = pdp_initial_stock)

summary(gam_depth)  # Check significance of terms

# 5. GAM without interaction (same smooth for all depths)
gam_noint <- gam(yhat ~ s(temperature, bs = "cs", k = 6) + depth_group,
                 data = pdp_initial_stock)

# Compare models to test if depth-specific smooths improve fit
anova(gam_noint, gam_depth, test = "Chisq")

# Moyenne marginale des prédictions pour chaque depth_group
emm <- emmeans(gam_depth, ~ depth_group)

# Comparaisons par paires
pairs(emm)

# Extraire les paires depuis emmeans
pair_df <- as.data.frame(pairs(emm)) %>%
  mutate(Significant = ifelse(p.value < 0.05, "Yes", "No")) %>%
  dplyr::select(contrast, estimate, Significant)

# Affichage
print(pair_df)

summary_pairs <- pair_df %>%
  filter(Significant == "Yes") %>%
  dplyr::rowwise() %>%
  mutate(
    group1 = strsplit(contrast, " - ")[[1]][1],
    group2 = strsplit(contrast, " - ")[[1]][2]
  ) %>%
  dplyr::select(group1, group2) %>%
  group_by(group1) %>%
  summarise(Significant_vs = paste(group2, collapse = ", "), .groups = "drop")

print(summary_pairs)


# 6. Predict at 25 t Ha with confidence intervals
library(dplyr)

# Target SOC values
target_soc <- c(20, 28)
depths <- levels(pdp_initial_stock$depth_group)

# 1. Predictions from GAM
pred_df <- lapply(target_soc, function(soc) {
  newdata <- data.frame(
    temperature = soc,
    depth_group = depths
  )
  pred <- predict(gam_depth, newdata = newdata, se.fit = TRUE)
  data.frame(
    target_soc = soc,
    depth_group = depths,
    fit = pred$fit,
    lower = pred$fit - 1.96 * pred$se.fit,
    upper = pred$fit + 1.96 * pred$se.fit
  )
}) %>% bind_rows()


# pred_df = GAM predictions
pred_df



# -----------------------------
# 11. Spearman correlation test by depth
# -----------------------------
# Average PDP per depth/SOC
pdp_unique <- pdp_initial_stock %>%
  group_by(depth_group, temperature) %>%
  summarise(yhat = mean(yhat), .groups = "drop")

# Interpolate for each observation
FULL_sun2_sub <- FULL_sun2_sub %>%
  rowwise() %>%
  mutate(
    yhat_num = approx(
      x = pdp_unique$temperature[pdp_unique$depth_group == depth_group],
      y = pdp_unique$yhat[pdp_unique$depth_group == depth_group],
      xout = temperature,
      rule = 2
    )$y
  ) %>%
  ungroup() %>%
  mutate(yhat_num = as.numeric(yhat_num))

# Compute Spearman correlation & trend table
summary_table <- FULL_sun2_sub %>%
  group_by(depth_group) %>%
  summarise(
    n_points = n(),
    rho = cor(yhat_num, temperature, method = "spearman"),
    p_val = cor.test(yhat_num, temperature, method = "spearman")$p.value,
    trend = case_when(
      rho > 0 ~ "increase",
      rho < 0 ~ "decrease",
      TRUE    ~ "none"
    ),
    reliability = case_when(
      n_points < 5 ~ "low",
      p_val > 0.05 ~ "moderate",
      TRUE ~ "high"
    ),
    .groups = "drop"
  ) %>%
  arrange(depth_group)

print(summary_table)

