###############################################################################
# 06 – Table 2: standardized XGBoost predictions per management category
#
# Changes vs. former 6_PDP_categorial.R:
#  * the former loop `colnames(Xb)[sapply(Xb, is.numeric)]` iterated over matrix
#    CELLS, not columns, so most covariates were not set as described in the
#    caption; counterfactual data are now built on the data frame and converted
#    with the same design-matrix function as the model;
#  * precipitation and temperature are sampled jointly (same site), keeping
#    their correlation;
#  * intervals come from refitting the model on study-level bootstrap samples
#    (the caption said "bootstrapped models (1,000 iterations)" but the former
#    code used one model and resampled rows 500 times);
#  * N is reported as observations AND studies per level;
#  * a support flag shows when the standardized time lies outside the range
#    observed for that category.
###############################################################################
suppressPackageStartupMessages(library(xgboost))

STD_TIME <- c(RR = 25, storage_rate = 5)   # years after conversion (caption of Table 2)
N_DRAW   <- 2000
CAT_VARS <- c("NEW_treatment_type2", "History_reclass", "diff_species_class")

counterfactual <- function(d, var, level, time) {
  base <- d[sample(nrow(d), N_DRAW, replace = TRUE), ]
  same_level <- d[d[[var]] == level, ]
  topsoil <- d$MEAN_depth[d$MEAN_depth < 30]
  clim <- d[sample(nrow(d), N_DRAW, replace = TRUE), c("precipitation", "temperature")]
  base[[var]] <- factor(level, levels = levels(d[[var]]))
  base$time_since_conversion <- time
  base$MEAN_depth <- sample(topsoil, N_DRAW, replace = TRUE)
  base$control_soc_mean_T_ha <- sample(same_level$control_soc_mean_T_ha, N_DRAW, replace = TRUE)
  base$precipitation <- clim$precipitation
  base$temperature <- clim$temperature
  make_X(base)
}

rows <- list()
for (nm in c("RR", "storage_rate")) {
  obj <- readRDS(dir_out("models", paste0("xgb_", nm, ".rds")))
  d <- obj$data; X <- obj$X; y <- d[[obj$y_name]]
  back <- if (nm == "RR") exp else identity
  t_std <- STD_TIME[[nm]]

  boot_models <- lapply(seq_len(N_BOOT_PDP), function(b) {
    idx <- boot_rows_by_study(d$id_article)
    fit_xgb(X[idx, , drop = FALSE], y[idx], obj$nrounds)
  })

  for (v in CAT_VARS) for (L in levels(d[[v]])) {
    sub <- d[d[[v]] == L, ]
    if (nrow(sub) < 3) next
    Xc <- counterfactual(d, v, L, t_std)
    est <- mean(predict(obj$model, Xc))
    bs <- sapply(boot_models, function(m) mean(predict(m, counterfactual(d, v, L, t_std))))
    rows[[length(rows) + 1]] <- tibble(
      metric = nm, predictor = v, level = L,
      xgb_estimate = back(est), xgb_lb = back(quantile(bs, .025)), xgb_ub = back(quantile(bs, .975)),
      n_obs = nrow(sub), n_studies = n_distinct(sub$id_article),
      max_time_observed = max(sub$time_since_conversion),
      standardized_time = t_std,
      time_within_support = t_std <= max(sub$time_since_conversion))
  }
}

xgb_tab <- bind_rows(rows) %>%
  mutate(xgb_formatted = sprintf("%.2f [%.2f; %.2f]", xgb_estimate, xgb_lb, xgb_ub))
meta_tab <- read_csv(dir_out("tables", "table2_metafor_and_conditions.csv"), show_col_types = FALSE) %>%
  select(metric, predictor, level, metafor_formatted = formatted, depth, initial_soc, precipitation, temperature)
table2 <- xgb_tab %>% left_join(meta_tab, by = c("metric", "predictor", "level"))
write_csv(table2, dir_out("tables", "Table2_full.csv"))

if (any(!table2$time_within_support))
  warning("Some categories have no observation as old as the standardized time: see column time_within_support.")
