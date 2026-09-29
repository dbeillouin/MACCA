suppressPackageStartupMessages(library(xgboost))

d <- read_csv(dir_out("data", "analysis_dataset_RR.csv"), show_col_types = FALSE) %>%
  filter(!is.na(treatment_soc_mean_T_ha), !is.na(control_soc_mean_T_ha),
         !is.na(treatment_soc_sd_T_ha),   !is.na(control_soc_sd_T_ha),
         !is.na(control_replicate_nb),    !is.na(treatment_replicate_nb),
         control_replicate_nb > 0, treatment_replicate_nb > 0, !is.na(yi)) %>%
  mutate(md   = treatment_soc_mean_T_ha - control_soc_mean_T_ha,
         v_md = treatment_soc_sd_T_ha^2 / treatment_replicate_nb +
                control_soc_sd_T_ha^2   / control_replicate_nb) %>%
  filter(v_md > 0) %>%
  filter(if_all(all_of(PREDICTORS), ~ !is.na(.x))) %>%
  freeze_levels()

record("abs_importance_n_obs_studies",
       sprintf("%d observations / %d studies, identical for both responses",
               nrow(d), n_distinct(d$id_article)), "Table S11")

X <- make_X(d)
studies <- unique(d$id_article)

xgb_rank <- function(y, lab) {
  cv <- xgb.cv(params = XGB_PARAMS, data = xgb.DMatrix(X, label = y), nrounds = 1000,
               nfold = 5, early_stopping_rounds = 25, verbose = 0)
  n <- if (!is.null(cv$best_iteration)) cv$best_iteration else cv$early_stop$best_iteration
  map_dfr(seq_len(N_BOOT_IMP), function(b) {
    s   <- sample(studies, length(studies), replace = TRUE)
    idx <- unlist(lapply(s, function(x) which(d$id_article == x)))
    xgb.importance(model = fit_xgb(X[idx, , drop = FALSE], y[idx], n)) %>%
      as_tibble() %>% mutate(variable = feature_origin(Feature)) %>%
      group_by(variable) %>% summarise(gain = sum(Gain), .groups = "drop") %>%
      mutate(share = 100 * gain / sum(gain))
  }) %>% group_by(variable) %>%
    summarise(m = round(mean(share), 1),
              a = round(quantile(share, .17), 1),
              b = round(quantile(share, .83), 1), .groups = "drop") %>%
    setNames(c("variable", paste0("xgb_", lab), paste0("xgb_", lab, "_q17"),
               paste0("xgb_", lab, "_q83")))
}

mf_rank <- function(yv, vv, lab) {
  if (!requireNamespace("metaforest", quietly = TRUE)) return(NULL)
  dd <- as.data.frame(d[, c("id_experiment", PREDICTORS)]); dd$.y <- yv; dd$.v <- vv
  m <- metaforest::MetaForest(as.formula(paste(".y ~", paste(PREDICTORS, collapse = " + "))),
         data = dd, vi = ".v", study = "id_experiment",
         whichweights = "random", num.trees = 8000)
  v <- m$forest$variable.importance
  record(paste0("abs_importance_MetaForest_OOB_R2_", lab), round(m$forest$r.squared, 3), "Table S11")
  setNames(tibble(variable = names(v), value = round(100 * v / sum(abs(v)), 1)),
           c("variable", paste0("mf_", lab)))
}

tab <- list(mf_rank(d$md, d$v_md, "absolute"), mf_rank(d$yi, d$vi, "ratio"),
            xgb_rank(d$md, "absolute"), xgb_rank(d$yi, "ratio")) %>%
  compact() %>% reduce(full_join, by = "variable")
if ("mf_absolute" %in% names(tab)) tab <- arrange(tab, desc(mf_absolute))

write_csv(tab, dir_out("tables", "importance_absolute_vs_RR.csv"))
for (v in c("control_soc_mean_T_ha", "time_since_conversion"))
  if (v %in% tab$variable && "mf_absolute" %in% names(tab))
    record(paste0("abs_importance_", v),
           sprintf("MetaForest: %.1f%% on the absolute difference vs %.1f%% on the response ratio",
                   tab$mf_absolute[tab$variable == v], tab$mf_ratio[tab$variable == v]),
           "3.2, 4.1, Table S11")
