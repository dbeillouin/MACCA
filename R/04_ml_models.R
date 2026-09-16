###############################################################################
# 04 – Machine-learning models (XGBoost)
#
# Changes vs. former 3__Models.R / 4_Importance_plot_GB.R:
#  * outliers identified by observation ID (not row numbers) and documented;
#  * the R² reported as "out-of-sample" is now truly out-of-sample: the former
#    script computed cor(y, predict(model, X))² on the TRAINING data. Folds are
#    grouped by study so that no study contributes to both training and test;
#  * the same datasets (after outlier removal) are used for CV, importance,
#    interactions and partial dependence (the former importance script used
#    the storage-rate data WITHOUT outlier removal but the RR data WITH it);
#  * bootstrap resamples studies, not rows (observations are nested in studies);
#  * pairwise interactions computed from SHAP interaction values of the final
#    model, with code (Fig. 3b and the interaction percentages in 3.2).
###############################################################################
suppressPackageStartupMessages(library(xgboost))

rr <- read_csv(dir_out("data", "analysis_dataset_RR.csv"), show_col_types = FALSE)
sq <- read_csv(dir_out("data", "analysis_dataset_storage_rate.csv"), show_col_types = FALSE)

ml <- list(RR = list(d = rr, y = "yi"), storage_rate = list(d = sq, y = "seq_rate"))
ml_out <- list()

for (nm in names(ml)) {
  message("-- ", nm)
  y_name <- ml[[nm]]$y
  d <- ml[[nm]]$d %>%
    filter(if_all(all_of(c(PREDICTORS, y_name)), ~ !is.na(.x))) %>%
    mutate(ID = as.character(ID))

  # ---------------------------------------------------------------------------
  # 1. Outliers (|residual| > 3 SD), identified by ID
  # ---------------------------------------------------------------------------
  if ((!exists("USE_METAFOREST") || USE_METAFOREST) && requireNamespace("metaforest", quietly = TRUE)) {
    d$.y <- d[[y_name]]
    d$.v <- if (nm == "RR") d$vi else d$seq_rate_vi
    f <- as.formula(paste(".y ~", paste(PREDICTORS, collapse = " + ")))
    mf <- metaforest::MetaForest(f, data = as.data.frame(d[, c(".y", ".v", "id_experiment", PREDICTORS)]),
                                 vi = ".v", study = "id_experiment",
                                 whichweights = "random", num.trees = 8000)
    res <- d$.y - predict(mf, data = as.data.frame(d))$predictions
    rule <- "MetaForest residuals"
    record(paste0(nm, "_MetaForest_OOB_R2"), mf$forest$r.squared, "3.2")
  } else {
    warning("metaforest not installed: outliers flagged from XGBoost out-of-fold residuals instead")
    dd <- freeze_levels(d); X0 <- make_X(dd); fo <- study_folds(dd$id_article, 5, SEED)
    oof <- numeric(nrow(dd))
    for (k in 1:5) {
      m0 <- fit_xgb(X0[fo != k, , drop = FALSE], dd[[y_name]][fo != k], 200)
      oof[fo == k] <- predict(m0, X0[fo == k, , drop = FALSE])
    }
    res <- dd[[y_name]] - oof
    rule <- "XGBoost out-of-fold residuals"
  }
  stopifnot(!anyNA(d$ID), !anyDuplicated(d$ID))
  outliers <- d$ID[abs(res) > 3 * sd(res)]
  write_csv(tibble(ID = outliers, rule = rule, dataset = nm), dir_out("data", paste0("outliers_", nm, ".csv")))
  record(paste0(nm, "_n_outliers_removed"), length(outliers), "2.5 (new outlier sentence)", rule)

  d <- d %>% filter(!ID %in% outliers) %>% freeze_levels()
  X <- make_X(d); y <- d[[y_name]]

  # ---------------------------------------------------------------------------
  # 2. Number of boosting rounds: early stopping on study-grouped folds
  # ---------------------------------------------------------------------------
  fo <- study_folds(d$id_article, 5, SEED)
  cv <- xgb.cv(params = XGB_PARAMS, data = xgb.DMatrix(X, label = y), nrounds = 1000,
               folds = lapply(1:5, function(k) which(fo == k)),
               early_stopping_rounds = 20, verbose = 0)
  best_n <- if (!is.null(cv$best_iteration)) cv$best_iteration else cv$early_stop$best_iteration  # xgboost 1.x / 2.x
  if (is.null(best_n)) best_n <- which.min(cv$evaluation_log$test_rmse_mean)
  record(paste0(nm, "_xgb_nrounds"), best_n, "2.5")

  # ---------------------------------------------------------------------------
  # 3. Predictive performance
  # ---------------------------------------------------------------------------
  perf <- map_dfr(seq_len(N_CV_REPEATS), function(r) {
    fo_g <- study_folds(d$id_article, 5, SEED + r)             # grouped by study
    fo_p <- study_folds(make_profile(d), 5, SEED + r)          # grouped by soil profile
    set.seed(SEED + r); fo_r <- sample(rep_len(1:5, nrow(d)))   # random rows (former practice)
    oof <- function(fo) {
      p <- numeric(nrow(d))
      for (k in 1:5) p[fo == k] <- predict(fit_xgb(X[fo != k, , drop = FALSE], y[fo != k], best_n),
                                           X[fo == k, , drop = FALSE])
      p
    }
    pg <- oof(fo_g); pr <- oof(fo_r); pp <- oof(fo_p)
    tibble(repeat_id = r, R2_grouped = r2(y, pg), RMSE_grouped = sqrt(mean((y - pg)^2)),
           R2_profile = r2(y, pp), RMSE_profile = sqrt(mean((y - pp)^2)),
           R2_random = r2(y, pr), RMSE_random = sqrt(mean((y - pr)^2)))
  })
  final <- fit_xgb(X, y, best_n)
  r2_in <- r2(y, predict(final, X))
  write_csv(perf, dir_out("tables", paste0("cv_performance_", nm, ".csv")))
  record(paste0(nm, "_R2_out_of_sample_study_grouped"), mean(perf$R2_grouped), "3.2, response letter",
         "report this one as out-of-sample")
  record(paste0(nm, "_RMSE_out_of_sample_study_grouped"), mean(perf$RMSE_grouped), "3.2")
  record(paste0(nm, "_R2_random_row_CV"), mean(perf$R2_random), "3.2",
         "within-dataset: other layers of the same profile can be in the training folds")
  record(paste0(nm, "_R2_profile_grouped_CV"), mean(perf$R2_profile), "3.2",
         "new agroforestry/control comparison (all depths of a profile held out together)")
  record(paste0(nm, "_n_profiles"), n_distinct(make_profile(d)), "3.2")
  record(paste0(nm, "_R2_CV_sd_across_repeats"),
         sprintf("random %.2f; profile %.2f; study %.2f", sd(perf$R2_random), sd(perf$R2_profile), sd(perf$R2_grouped)), "3.2")
  record(paste0(nm, "_R2_in_sample"), r2_in, "3.2", "former scripts reported this as R²")

  # ---------------------------------------------------------------------------
  # 4. Variable importance (Gain), cluster bootstrap
  # ---------------------------------------------------------------------------
  imp <- map_dfr(seq_len(N_BOOT_IMP), function(b) {
    idx <- boot_rows_by_study(d$id_article)
    mb <- fit_xgb(X[idx, , drop = FALSE], y[idx], best_n)
    xgb.importance(model = mb) %>%
      as_tibble() %>%
      mutate(variable = feature_origin(Feature)) %>%
      group_by(variable) %>% summarise(gain = sum(Gain), .groups = "drop") %>%
      mutate(share = 100 * gain / sum(gain), run = b)
  }) %>%
    complete(variable = PREDICTORS, run, fill = list(gain = 0, share = 0))
  imp_sum <- imp %>% group_by(variable) %>%
    summarise(share_mean = mean(share), share_q17 = quantile(share, .17), share_q83 = quantile(share, .83),
              .groups = "drop") %>%
    arrange(desc(share_mean))
  write_csv(imp_sum, dir_out("tables", paste0("importance_", nm, ".csv")))
  for (i in seq_len(nrow(imp_sum)))
    record(paste0(nm, "_importance_pct_", imp_sum$variable[i]), imp_sum$share_mean[i], "3.2, Abstract, Conclusion, Fig. 3a")

  # ---------------------------------------------------------------------------
  # 5. Interactions: SHAP interaction values of the final model
  # ---------------------------------------------------------------------------
  shap <- predict(final, X, predinteraction = TRUE)        # n x (p+1) x (p+1), last = bias
  p <- ncol(X); fnames <- colnames(X)
  A <- apply(abs(shap[, 1:p, 1:p, drop = FALSE]), c(2, 3), mean)
  dimnames(A) <- list(fnames, fnames)
  g <- feature_origin(fnames)
  G <- sapply(PREDICTORS, function(cj) sapply(PREDICTORS, function(ci) sum(A[g == ci, g == cj, drop = FALSE])))
  # rows = focal variable, columns = partner
  main  <- diag(G)
  inter <- rowSums(G) - main
  inter_share <- tibble(variable = PREDICTORS, main = main, interaction = inter,
                        interaction_share = inter / (main + inter))
  pair <- G; diag(pair) <- NA
  pair_share <- sweep(pair, 1, rowSums(pair, na.rm = TRUE), "/") * 100   # rows sum to 100 %
  write_csv(inter_share, dir_out("tables", paste0("interaction_share_", nm, ".csv")))
  write.csv(round(pair_share, 1), dir_out("tables", paste0("interaction_pairs_rowpct_", nm, ".csv")))
  for (i in seq_len(nrow(inter_share)))
    record(paste0(nm, "_interaction_share_", inter_share$variable[i]), inter_share$interaction_share[i],
           "3.2, 4 (0.86 / 0.807 ...)")
  record(paste0(nm, "_depth_share_of_initialSOC_interactions_pct"),
         pair_share["control_soc_mean_T_ha", "MEAN_depth"], "3.4 (63%)")

  saveRDS(list(model = final, X = X, data = d, y_name = y_name, nrounds = best_n),
          dir_out("models", paste0("xgb_", nm, ".rds")))
  ml_out[[nm]] <- imp_sum
}
