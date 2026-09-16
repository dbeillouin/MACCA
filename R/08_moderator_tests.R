###############################################################################
# 08 – Moderator tests on OBSERVED effect sizes (three-level meta-regressions)
# Complements the machine-learning analyses, whose out-of-sample performance on
# unseen studies is low: every claim about drivers should also be supported
# here. One model per moderator, ML estimation, LRT against the null model.
###############################################################################
suppressPackageStartupMessages({ library(metafor); library(splines) })

out <- list()
for (nm in c("RR", "storage_rate")) {
  d <- read_csv(dir_out("data", paste0("analysis_dataset_", nm, ".csv")), show_col_types = FALSE)
  d$.y <- if (nm == "RR") d$yi else d$seq_rate
  d$.v <- if (nm == "RR") d$vi else d$seq_rate_vi
  d <- d %>% filter(if_all(all_of(PREDICTORS), ~ !is.na(.x)))
  rnd <- list(~ 1 | id_article, ~ 1 | id_experiment)
  m0 <- rma.mv(.y, .v, random = rnd, data = d, method = "ML")
  spec <- list(
    initial_SOC_linear = ~ control_soc_mean_T_ha, initial_SOC_spline = ~ ns(control_soc_mean_T_ha, 3),
    time_linear = ~ time_since_conversion, time_spline = ~ ns(time_since_conversion, 3),
    precipitation_spline = ~ ns(precipitation, 3), temperature_spline = ~ ns(temperature, 3),
    depth_spline = ~ ns(MEAN_depth, 3), depth_linear = ~ MEAN_depth,
    agroforestry_type = ~ factor(NEW_treatment_type2), land_use_history = ~ factor(History_reclass),
    species_added = ~ factor(diff_species_class), crop = ~ factor(main_culture2), design = ~ factor(Grouped_Design),
    all_biophysical = ~ ns(control_soc_mean_T_ha, 3) + ns(time_since_conversion, 3) + ns(precipitation, 3) +
                        ns(temperature, 3) + MEAN_depth,
    all_management  = ~ factor(NEW_treatment_type2) + factor(History_reclass) + factor(diff_species_class) +
                        factor(main_culture2) + factor(Grouped_Design))
  # heterogeneity explained: proportional reduction of the variance components (REML fits)
  m0_reml <- rma.mv(.y, .v, random = rnd, data = d, method = "REML")
  r2_het <- function(mods) {
    m <- tryCatch(rma.mv(.y, .v, mods = mods, random = rnd, data = d, method = "REML"), error = function(e) NULL)
    if (is.null(m)) return(c(total = NA, between = NA, within = NA))
    setNames(pmax(0, c(1 - sum(m$sigma2) / sum(m0_reml$sigma2),
                       1 - m$sigma2[1] / m0_reml$sigma2[1],
                       1 - m$sigma2[2] / m0_reml$sigma2[2])) * 100,
             c("total", "between", "within"))
  }
  for (k in names(spec)) {
    m <- tryCatch(rma.mv(.y, .v, mods = spec[[k]], random = rnd, data = d, method = "ML"), error = function(e) NULL)
    if (is.null(m)) next
    a <- anova(m0, m)
    slope <- if (grepl("linear", k)) sprintf("%.4f [%.4f; %.4f]", m$b[2], m$ci.lb[2], m$ci.ub[2]) else NA_character_
    h <- r2_het(spec[[k]])
    out[[length(out) + 1]] <- tibble(metric = nm, moderator = k, n_obs = nrow(d), n_studies = n_distinct(d$id_article),
                                     df = a$parms.f - a$parms.r, LRT = a$LRT, p = a$pval, slope_log_or_rate = slope,
                                     heterogeneity_explained_total_pct = h["total"],
                                     heterogeneity_explained_between_study_pct = h["between"],
                                     heterogeneity_explained_within_study_pct = h["within"])
    if (k %in% c("all_biophysical", "all_management"))
      record(paste0(nm, "_heterogeneity_explained_", k),
             sprintf("total %.1f%%; between-study %.1f%%; within-study %.1f%%", h["total"], h["between"], h["within"]),
             "3.2, 4, response letter", "proportional reduction of variance components (REML)")
    record(paste0(nm, "_modtest_", k), sprintf("LRT=%.2f (df=%d), p=%.3g%s", a$LRT, a$parms.f - a$parms.r, a$pval,
                                               ifelse(is.na(slope), "", paste0("; slope ", slope))), "3.2-3.6, 4")
  }
}
write_csv(bind_rows(out), dir_out("tables", "moderator_tests_observed_data.csv"))
