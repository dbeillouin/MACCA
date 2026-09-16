###############################################################################
# 07 – Data support by soil depth and time since conversion (Reviewer #2)
# Provides the [n], [N], [k] placeholders of the Discussion caveat and a
# supplementary table that shows where predictions are supported by data.
###############################################################################

for (nm in c("RR", "storage_rate")) {
  obj <- readRDS(dir_out("models", paste0("xgb_", nm, ".rds")))
  d <- obj$data

  cov <- d %>%
    mutate(time_class = cut(time_since_conversion, c(0, 5, 10, 20, Inf),
                            labels = c("<=5 y", "5-10 y", "10-20 y", ">20 y"), include.lowest = TRUE)) %>%
    group_by(depth_group, time_class, .drop = FALSE) %>%
    summarise(n_obs = n(), n_studies = n_distinct(id_article), .groups = "drop")
  write_csv(cov, dir_out("tables", paste0("coverage_depth_by_time_", nm, ".csv")))

  old <- d %>% filter(time_since_conversion > 20)
  record(paste0(nm, "_n_obs_total_after_outliers"), nrow(d), "4 caveat [N]")
  record(paste0(nm, "_n_obs_older_than_20y"), nrow(old), "4 caveat [n]")
  record(paste0(nm, "_n_studies_older_than_20y"), n_distinct(old$id_article), "4 caveat [k]")
  record(paste0(nm, "_max_depth_midpoint_older_than_20y_cm"), if (nrow(old)) max(old$MEAN_depth) else NA, "4 caveat, response letter")
  record(paste0(nm, "_max_time_below_40cm_y"),
         if (any(d$MEAN_depth > 40)) max(d$time_since_conversion[d$MEAN_depth > 40]) else NA, "response letter")
  record(paste0(nm, "_n_obs_below_40cm_and_20y_or_more"),
         sum(d$MEAN_depth > 40 & d$time_since_conversion >= 20), "response letter")
}
