###############################################################################
# 03 – Classical meta-analysis
#
# Changes vs. former scripts:
#  * confidence intervals taken from the fitted model (the former RR script
#    back-transformed hard-coded bounds: exp(0.0217), exp(0.1155));
#  * location-scale models use sampling VARIANCES (the former storage-rate
#    script passed standard deviations as `vi`);
#  * everything is fitted on the datasets written by step 02.
###############################################################################
suppressPackageStartupMessages(library(metafor))

rr <- read_csv(dir_out("data", "analysis_dataset_RR.csv"), show_col_types = FALSE)
sq <- read_csv(dir_out("data", "analysis_dataset_storage_rate.csv"), show_col_types = FALSE)

datasets <- list(
  RR           = list(d = rr, y = "yi",       v = "vi",          back = function(x) 100 * (exp(x) - 1), unit = "%"),
  storage_rate = list(d = sq, y = "seq_rate", v = "seq_rate_vi", back = identity,                     unit = "Mg C ha-1 yr-1")
)

fit3 <- function(d, y, v, method = "REML", mods = NULL, ...) {
  d$.y <- d[[y]]; d$.v <- d[[v]]
  rnd <- list(~ 1 | id_article, ~ 1 | id_experiment)
  if (is.null(mods)) rma.mv(.y, .v, random = rnd, data = d, method = method, ...)
  else rma.mv(.y, .v, mods = mods, random = rnd, data = d, method = method, ...)
}

# Multilevel I² (Nakagawa & Santos 2012)
i2_multilevel <- function(m) {
  W <- diag(1 / m$vi); X <- model.matrix(m)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  typical_v <- (m$k - m$p) / sum(diag(P))
  tot <- sum(m$sigma2) + typical_v
  c(total = 100 * sum(m$sigma2) / tot,
    between_study = 100 * m$sigma2[1] / tot,
    within_study  = 100 * m$sigma2[2] / tot)
}

meta_tables <- list()

for (nm in names(datasets)) {
  S <- datasets[[nm]]; d <- S$d
  message("-- ", nm)

  # ---- overall effect -------------------------------------------------------
  m <- fit3(d, S$y, S$v)
  est <- c(est = as.numeric(coef(m)), lb = m$ci.lb, ub = m$ci.ub)
  record(paste0(nm, "_overall_estimate_", S$unit), S$back(est["est"]), "3.1, Abstract, Conclusion")
  record(paste0(nm, "_overall_CI"), paste0("[", round(S$back(est["lb"]), 2), "; ", round(S$back(est["ub"]), 2), "]"), "3.1")
  i2 <- i2_multilevel(m)
  record(paste0(nm, "_I2_total"), i2["total"], "3.1, Abstract")
  record(paste0(nm, "_I2_between_study_share_of_heterogeneity"), 100 * i2["between_study"] / i2["total"], "3.1 (~60%)")

  # ---- is the three-level structure needed? ---------------------------------
  m_within0  <- fit3(d, S$y, S$v, sigma2 = c(NA, 0))
  m_between0 <- fit3(d, S$y, S$v, sigma2 = c(0, NA))
  record(paste0(nm, "_LRT_within_study_level"), paste0("LRT=", round(anova(m, m_within0)$LRT, 2), ", p=", signif(anova(m, m_within0)$pval, 3)), "2.4")
  record(paste0(nm, "_LRT_between_study_level"), paste0("LRT=", round(anova(m, m_between0)$LRT, 2), ", p=", signif(anova(m, m_between0)$pval, 3)), "2.4")

  # ---- topsoil (0-30 cm) ----------------------------------------------------
  if ("soil_depth_end" %in% names(d)) {
    dt <- d %>% filter(soil_depth_end <= 30)
    if (nrow(dt) > 10) {
      mt <- fit3(dt, S$y, S$v)
      record(paste0(nm, "_topsoil_0_30_estimate_", S$unit),
             paste0(round(S$back(coef(mt)), 2), " [", round(S$back(mt$ci.lb), 2), "; ", round(S$back(mt$ci.ub), 2), "]"),
             "3.1, Abstract", "subset of observations with lower sampling bound <= 30 cm; check this matches the former definition")
    }
  }

  # ---- language (ML fits for the LRT) ---------------------------------------
  dl <- d %>% filter(!is.na(language))
  m0l <- fit3(dl, S$y, S$v, method = "ML")
  m1l <- fit3(dl, S$y, S$v, method = "ML", mods = ~ language)
  a <- anova(m0l, m1l)
  record(paste0(nm, "_language_LRT"), paste0("LRT=", round(a$LRT, 2), ", p=", signif(a$pval, 3)), "3.1")

  # ---- small-study effects: multilevel Egger test ----------------------------
  d$.sei <- sqrt(d[[S$v]])
  me <- fit3(d, S$y, S$v, mods = ~ .sei)
  record(paste0(nm, "_Egger_p"), me$pval[2], "3.1")

  # ---- precision vs duration and design (location-scale models) -------------
  ds <- d %>% filter(!is.na(time_since_conversion))
  ms <- rma(yi = ds[[S$y]], vi = ds[[S$v]], scale = ~ time_since_conversion, data = ds)
  tau <- function(t) sqrt(exp(ms$alpha[1] + ms$alpha[2] * t))
  record(paste0(nm, "_tau_change_2_to_10_years_pct"), 100 * (tau(10) / tau(2) - 1), "3.1",
         "residual heterogeneity SD; the text currently says 'uncertainty (SD)'")
  msd <- rma(yi = ds[[S$y]], vi = ds[[S$v]], scale = ~ Grouped_Design + time_since_conversion, data = ds)
  al <- msd$alpha[, 1]
  rn <- rownames(msd$alpha)
  if (any(grepl("Randomized", rn))) {
    ratio <- sqrt(exp(al[grepl("Randomized", rn)]))  # SD ratio vs reference design (same time)
    record(paste0(nm, "_SD_randomized_vs_reference_design_pct"), 100 * (ratio - 1), "3.1 (+481%)",
           paste("reference level:", levels(factor(ds$Grouped_Design))[1]))
  }
  record(paste0(nm, "_scale_model_LRT_design_time"),
         paste0("p=", signif(anova(msd, rma(yi = ds[[S$y]], vi = ds[[S$v]], scale = ~ 1, data = ds))$pval, 3)), "3.1, Fig. S7")

  # ---- Table 2: subgroup estimates -------------------------------------------
  for (v in c("NEW_treatment_type2", "History_reclass", "diff_species_class")) {
    dv <- d %>% filter(!is.na(.data[[v]]))
    dv$.g <- factor(dv[[v]])
    mv <- tryCatch(fit3(dv, S$y, S$v, mods = ~ -1 + .g), error = function(e) NULL)
    if (is.null(mv)) next
    lev <- sub("^\\.g", "", rownames(mv$b))
    env <- dv %>% group_by(level = .g) %>%
      summarise(n_obs = n(), n_studies = n_distinct(id_article),
                depth = sprintf("%.1f ± %.1f", mean(MEAN_depth, na.rm = TRUE), sd(MEAN_depth, na.rm = TRUE)),
                initial_soc = sprintf("%.1f ± %.1f", mean(control_soc_mean_T_ha, na.rm = TRUE), sd(control_soc_mean_T_ha, na.rm = TRUE)),
                precipitation = sprintf("%.0f ± %.0f", mean(precipitation, na.rm = TRUE), sd(precipitation, na.rm = TRUE)),
                temperature = sprintf("%.1f ± %.1f", mean(temperature, na.rm = TRUE), sd(temperature, na.rm = TRUE)),
                .groups = "drop")
    tab <- tibble(metric = nm, predictor = v, level = lev,
                  estimate = if (nm == "RR") exp(as.numeric(mv$b)) else as.numeric(mv$b),
                  ci_lb = if (nm == "RR") exp(mv$ci.lb) else mv$ci.lb,
                  ci_ub = if (nm == "RR") exp(mv$ci.ub) else mv$ci.ub) %>%
      left_join(env %>% mutate(level = as.character(level)), by = "level")
    meta_tables[[paste(nm, v)]] <- tab
  }
}

table2_meta <- bind_rows(meta_tables) %>%
  mutate(across(c(estimate, ci_lb, ci_ub), ~ round(.x, 2)),
         formatted = sprintf("%.2f [%.2f; %.2f]", estimate, ci_lb, ci_ub))
write_csv(table2_meta, dir_out("tables", "table2_metafor_and_conditions.csv"))
message("Table 2 (metafor part) written; RR estimates are ratios, N = observations and studies per level.")
