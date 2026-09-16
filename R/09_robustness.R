###############################################################################
# 09 – Robustness of the retained findings (observed effect sizes only)
#
#  A. Initial SOC gradient            – predictions, slope, leave-one-study-out
#  B. Temperature                     – data density, threshold vs linear,
#                                       confounding with system type, LOSO
#  C. Initial SOC x temperature       – interaction, slopes at 22 and 28 °C, LOSO
#  D. Short-term studies              – time x design (inflated early gains?)
#  E. Depth x precipitation           – exploratory, LOSO, number of deep studies
#  F. Grey literature                 – document type as moderator
#
# Leave-one-study-out (LOSO): every model is refitted without each study in
# turn; a finding is flagged "study-dependent" if removing a single study
# makes p > 0.05 or changes the sign of the key coefficient.
###############################################################################
suppressPackageStartupMessages({ library(metafor); library(splines) })
if (!exists("N_LOSO_MAX")) N_LOSO_MAX <- Inf

rnd <- list(~ 1 | id_article, ~ 1 | id_experiment)
fitml <- function(d, mods = NULL) {
  if (is.null(mods)) rma.mv(.y, .v, random = rnd, data = d, method = "ML")
  else rma.mv(.y, .v, mods = mods, random = rnd, data = d, method = "ML")
}
lrt <- function(d, mods0, mods1) {
  m0 <- tryCatch(fitml(d, mods0), error = function(e) NULL)
  m1 <- tryCatch(fitml(d, mods1), error = function(e) NULL)
  if (is.null(m0) || is.null(m1)) return(list(LRT = NA, p = NA, m1 = NULL))
  a <- anova(m0, m1); list(LRT = a$LRT, p = a$pval, m1 = m1)
}
# LOSO: returns one row per removed study with p-value and a chosen coefficient
loso <- function(d, mods0, mods1, coef_name = NULL) {
  studies <- unique(d$id_article)
  if (is.finite(N_LOSO_MAX)) studies <- head(studies, N_LOSO_MAX)
  map_dfr(studies, function(s) {
    r <- lrt(d[d$id_article != s, ], mods0, mods1)
    b <- if (!is.null(coef_name) && !is.null(r$m1) && coef_name %in% rownames(r$m1$b)) r$m1$b[coef_name, 1] else NA
    tibble(removed_study = s, LRT = r$LRT, p = r$p, coef = b)
  })
}
loso_summary <- function(tab, full_coef = NA, full_p = 0) {
  if (!is.na(full_p) && full_p > 0.05) return(sprintf("not significant with all studies (p = %.3g): LOSO not informative", full_p))
  flips <- tab %>% filter(p > 0.05 | (!is.na(coef) & !is.na(full_coef) & sign(coef) != sign(full_coef)))
  sprintf("max p without one study = %.3g; studies whose removal makes p > 0.05 or flips sign: %s",
          max(tab$p, na.rm = TRUE), if (nrow(flips)) paste(flips$removed_study, collapse = ", ") else "none")
}
pred_ns <- function(m, x_obs, grid, back) {
  X <- unname(as.matrix(predict(ns(x_obs, 3), grid)))
  p <- predict(m, newmods = X)
  tibble(x = grid, pred = back(p$pred), lb = back(p$ci.lb), ub = back(p$ci.ub))
}
fmt_pred <- function(tb, unit = "") paste0(tb$x, unit, ": ", sprintf("%.2f [%.2f; %.2f]", tb$pred, tb$lb, tb$ub), collapse = " | ")

rob <- list()
for (nm in c("RR", "storage_rate")) {
  message("-- ", nm)
  d <- read_csv(dir_out("data", paste0("analysis_dataset_", nm, ".csv")), show_col_types = FALSE) %>%
    filter(!is.na(control_soc_mean_T_ha), !is.na(temperature), !is.na(time_since_conversion),
           !is.na(precipitation), !is.na(MEAN_depth), !is.na(NEW_treatment_type2))
  d$.y <- if (nm == "RR") d$yi else d$seq_rate
  d$.v <- if (nm == "RR") d$vi else d$seq_rate_vi
  back <- if (nm == "RR") function(x) 100 * (exp(x) - 1) else identity
  S <- function(k) paste0(nm, "_rob_", k)

  # ---------------------------------------------------------------- A. initial SOC
  soc_grid <- c(20, 40, 60, 90, 130); soc_grid <- soc_grid[soc_grid <= max(d$control_soc_mean_T_ha)]
  m_soc <- rma.mv(.y, .v, mods = ~ ns(control_soc_mean_T_ha, 3), random = rnd, data = d, method = "REML")
  record(S("A_initialSOC_predictions"), fmt_pred(pred_ns(m_soc, d$control_soc_mean_T_ha, soc_grid, back)), "3.x initial SOC, Abstract")
  bins <- d %>% mutate(bin = cut(control_soc_mean_T_ha, c(-Inf, 30, 50, 70, 100, Inf))) %>%
    group_by(bin) %>% summarise(n_obs = n(), n_studies = n_distinct(id_article), .groups = "drop")
  record(S("A_initialSOC_data_density"), paste0(bins$bin, ": ", bins$n_obs, " obs/", bins$n_studies, " studies", collapse = " | "), "3.x, caption")
  lin <- lrt(d, NULL, ~ I(control_soc_mean_T_ha / 10))
  b10 <- lin$m1$b[2, 1]
  record(S("A_initialSOC_slope_per_10MgC"),
         sprintf("%.3f [%.3f; %.3f] (%s), p = %.3g", back(b10) - back(0), back(lin$m1$ci.lb[2]) - back(0),
                 back(lin$m1$ci.ub[2]) - back(0), if (nm == "RR") "percentage points" else "Mg C ha-1 yr-1", lin$p), "3.x")
  L <- loso(d, NULL, ~ I(control_soc_mean_T_ha / 10), "I(control_soc_mean_T_ha/10)")
  record(S("A_initialSOC_LOSO"), loso_summary(L, b10, lin$p), "3.x, response letter")
  rob[[paste(nm, "A")]] <- L %>% mutate(metric = nm, test = "initial SOC linear")

  # ---------------------------------------------------------------- B. temperature
  tb <- d %>% mutate(bin = cut(temperature, c(-Inf, 20, 22, 24, 26, 28, Inf))) %>%
    group_by(bin) %>% summarise(n_obs = n(), n_studies = n_distinct(id_article), .groups = "drop")
  record(S("B_temperature_data_density"), paste0(tb$bin, ": ", tb$n_obs, " obs/", tb$n_studies, " studies", collapse = " | "), "3.x, caption")
  t_grid <- c(18, 21, 24, 27, 29); t_grid <- t_grid[t_grid >= min(d$temperature) & t_grid <= max(d$temperature)]
  m_t <- rma.mv(.y, .v, mods = ~ ns(temperature, 3), random = rnd, data = d, method = "REML")
  record(S("B_temperature_predictions"), fmt_pred(pred_ns(m_t, d$temperature, t_grid, back), " C"), "3.x")
  # threshold vs linear vs spline (AIC, ML)
  aic <- function(mods) tryCatch(AIC(fitml(d, mods)), error = function(e) NA)
  a_lin <- aic(~ temperature); a_spl <- aic(~ ns(temperature, 3))
  hinge <- map_dfr(22:28, function(k) tibble(k = k, AIC = aic(as.formula(paste0("~ temperature + pmax(temperature - ", k, ", 0)")))))
  best <- hinge[which.min(hinge$AIC), ]
  record(S("B_temperature_shape_AIC"),
         sprintf("linear %.1f; spline %.1f; best threshold at %d C: %.1f (dAIC vs linear = %.1f)",
                 a_lin, a_spl, best$k, best$AIC, best$AIC - a_lin), "3.x",
         "a threshold is supported only if dAIC < -2 AND the LOSO below is stable")
  # confounding with system type, initial SOC and precipitation
  c1 <- lrt(d, ~ factor(NEW_treatment_type2), ~ factor(NEW_treatment_type2) + ns(temperature, 3))
  c2 <- lrt(d, ~ factor(NEW_treatment_type2) + ns(control_soc_mean_T_ha, 3) + ns(precipitation, 3),
            ~ factor(NEW_treatment_type2) + ns(control_soc_mean_T_ha, 3) + ns(precipitation, 3) + ns(temperature, 3))
  record(S("B_temperature_given_system_type"), sprintf("LRT=%.2f, p=%.3g", c1$LRT, c1$p), "3.x, 4")
  record(S("B_temperature_given_type_SOC_precip"), sprintf("LRT=%.2f, p=%.3g", c2$LRT, c2$p), "3.x, 4")
  sys_t <- d %>% group_by(NEW_treatment_type2) %>%
    summarise(mean_T = mean(temperature), share_above_26C = mean(temperature > 26), n = n(), .groups = "drop")
  record(S("B_temperature_by_system"), paste0(sys_t$NEW_treatment_type2, ": ", round(sys_t$mean_T, 1), " C (",
                                              round(100 * sys_t$share_above_26C), "% > 26 C, n=", sys_t$n, ")", collapse = " | "), "4")
  if ("altitude" %in% names(d))
    record(S("B_corr_temperature_altitude"), cor(d$temperature, as.numeric(d$altitude), use = "complete.obs"), "4")
  tl <- lrt(d, NULL, ~ temperature)
  L <- loso(d, NULL, ~ temperature, "temperature")
  record(S("B_temperature_linear_LOSO"), loso_summary(L, tl$m1$b[2, 1], tl$p), "3.x, response letter")
  tsp <- lrt(d, NULL, ~ ns(temperature, 3))
  L2 <- loso(d, NULL, ~ ns(temperature, 3))
  record(S("B_temperature_spline_LOSO"), loso_summary(L2, NA, tsp$p), "3.x, response letter",
         "non-linear shape: check monotonicity in B_temperature_predictions before interpreting")
  rob[[paste(nm, "B2")]] <- L2 %>% mutate(metric = nm, test = "temperature spline")
  # sensitivity: temperatures flagged "to_verify" in data/data_corrections.csv
  if (!exists("CORRECTIONS_FILE")) CORRECTIONS_FILE <- file.path(PROJECT_DIR, "data", "data_corrections.csv")
  tv <- flag_to_verify(d, CORRECTIONS_FILE, "temperature")
  if (any(tv)) {
    dv <- d[!tv, ]
    s1 <- lrt(dv, NULL, ~ ns(temperature, 3))
    dv$soc_c <- (dv$control_soc_mean_T_ha - 50) / 10; dv$t_c <- dv$temperature - 25
    s2 <- lrt(dv, ~ soc_c + t_c, ~ soc_c * t_c)
    record(S("B_sensitivity_without_temperatures_to_verify"),
           sprintf("%d obs removed; temperature spline p=%.3g; SOC x temperature p=%.3g", sum(tv), s1$p, s2$p), "response letter")
  }
  # dependence on the hottest sites (few studies above 28 C)
  hot <- d %>% filter(temperature > 28)
  cool <- d %>% filter(temperature <= 28)
  th <- lrt(cool, NULL, ~ ns(temperature, 3))
  record(S("B_temperature_without_sites_above_28C"),
         sprintf("sites > 28 C: %d obs / %d studies; spline test without them: LRT=%.2f, p=%.3g",
                 nrow(hot), n_distinct(hot$id_article), th$LRT, th$p), "3.x, caption")
  # temperature or elevation? (r ~ -0.8 in these data)
  if ("altitude" %in% names(d)) {
    d$alt <- suppressWarnings(as.numeric(d$altitude))
    da <- d %>% filter(!is.na(alt))
    record(S("B_AIC_temperature_vs_altitude"),
           sprintf("spline temperature AIC %.1f vs spline altitude AIC %.1f (n=%d)",
                   tryCatch(AIC(fitml(da, ~ ns(temperature, 3))), error = function(e) NA),
                   tryCatch(AIC(fitml(da, ~ ns(alt, 3))), error = function(e) NA), nrow(da)), "4")
  }
  rob[[paste(nm, "B")]] <- L %>% mutate(metric = nm, test = "temperature linear")

  # ---------------------------------------------------------------- C. SOC x temperature
  d$soc_c <- (d$control_soc_mean_T_ha - 50) / 10       # per 10 Mg C ha-1, centred at 50
  d$t_c <- d$temperature - 25
  ci <- lrt(d, ~ soc_c + t_c, ~ soc_c * t_c)
  record(S("C_SOCxTemperature_LRT"), sprintf("LRT=%.2f, p=%.3g", ci$LRT, ci$p), "3.x, 4")
  if (!is.null(ci$m1)) {
    b <- ci$m1$b[, 1]; V <- ci$m1$vb
    slope_at <- function(T) {
      g <- c(0, 1, 0, T - 25); est <- sum(g * b); se <- sqrt(as.numeric(t(g) %*% V %*% g))
      sprintf("%d C: %.3f [%.3f; %.3f]", T, est, est - 1.96 * se, est + 1.96 * se)
    }
    record(S("C_SOC_slope_per_10MgC_at_22_and_28C_logscale"), paste(slope_at(22), slope_at(28), sep = " | "), "3.x",
           if (nm == "RR") "log response ratio per 10 Mg C ha-1" else "Mg C ha-1 yr-1 per 10 Mg C ha-1")
  }
  L <- loso(d, ~ soc_c + t_c, ~ soc_c * t_c, "soc_c:t_c")
  record(S("C_SOCxTemperature_LOSO"), loso_summary(L, if (!is.null(ci$m1)) ci$m1$b["soc_c:t_c", 1] else NA, ci$p), "3.x, response letter")
  ct <- lrt(d, ~ factor(NEW_treatment_type2) + soc_c + t_c, ~ factor(NEW_treatment_type2) + soc_c * t_c)
  record(S("C_SOCxTemperature_given_system_type"), sprintf("LRT=%.2f, p=%.3g", ct$LRT, ct$p), "3.x")
  rob[[paste(nm, "C")]] <- L %>% mutate(metric = nm, test = "SOC x temperature")

  # ---------------------------------------------------------------- D. short-term studies
  d$short <- factor(ifelse(d$time_since_conversion <= 3, "<=3 y", ">3 y"))
  tab_short <- d %>% group_by(short, Grouped_Design) %>% summarise(n_obs = n(), n_studies = n_distinct(id_article), .groups = "drop")
  record(S("D_short_term_by_design_counts"), paste0(tab_short$short, " ", tab_short$Grouped_Design, ": ", tab_short$n_obs, "/", tab_short$n_studies, collapse = " | "), "3.x")
  sd_ <- lrt(d, ~ short + factor(Grouped_Design), ~ short * factor(Grouped_Design))
  record(S("D_short_x_design_LRT"), sprintf("LRT=%.2f, p=%.3g", sd_$LRT, sd_$p), "3.x, 4")
  ms <- tryCatch(rma.mv(.y, .v, mods = ~ short - 1, random = rnd, data = d, method = "REML"), error = function(e) NULL)
  if (!is.null(ms)) record(S("D_effect_short_vs_longer"),
                           paste0(sub("short", "", rownames(ms$b)), ": ", sprintf("%.2f [%.2f; %.2f]", back(ms$b[, 1]), back(ms$ci.lb), back(ms$ci.ub)), collapse = " | "), "3.x")
  s_ctrl <- lrt(d, ~ short, ~ short + ns(control_soc_mean_T_ha, 3) + ns(temperature, 3))
  record(S("D_short_term_effect_persists_given_SOC_temp"),
         { m_a <- tryCatch(fitml(d, ~ ns(control_soc_mean_T_ha, 3) + ns(temperature, 3)), error = function(e) NULL)
           m_b <- tryCatch(fitml(d, ~ short + ns(control_soc_mean_T_ha, 3) + ns(temperature, 3)), error = function(e) NULL)
           if (is.null(m_a) || is.null(m_b)) NA else sprintf("p=%.3g", anova(m_a, m_b)$pval) }, "3.x")

  # ---------------------------------------------------------------- E. depth x precipitation (exploratory)
  dp <- lrt(d, ~ precipitation + MEAN_depth, ~ precipitation * MEAN_depth)
  deep <- d %>% filter(MEAN_depth >= 30)
  record(S("E_depthxprecip_LRT"), sprintf("LRT=%.2f, p=%.3g; deep (>=30 cm): %d obs / %d studies, time %s-%s y",
                                          dp$LRT, dp$p, nrow(deep), n_distinct(deep$id_article),
                                          round(min(deep$time_since_conversion), 1), round(max(deep$time_since_conversion), 1)), "3.x (exploratory)")
  L <- loso(d, ~ precipitation + MEAN_depth, ~ precipitation * MEAN_depth, "precipitation:MEAN_depth")
  record(S("E_depthxprecip_LOSO"), loso_summary(L, if (!is.null(dp$m1)) dp$m1$b["precipitation:MEAN_depth", 1] else NA, dp$p), "response letter")
  rob[[paste(nm, "E")]] <- L %>% mutate(metric = nm, test = "depth x precipitation")

  # ---------------------------------------------------------------- F. grey literature
  if ("document_type" %in% names(d)) {
    d$doc <- factor(ifelse(tolower(d$document_type) == "article", "peer-reviewed article", "grey literature"))
    g <- lrt(d, NULL, ~ doc)
    record(S("F_grey_literature"), sprintf("LRT=%.2f, p=%.3g; grey: %d obs / %d studies", g$LRT, g$p,
                                           sum(d$doc == "grey literature"), n_distinct(d$id_article[d$doc == "grey literature"])), "3.1, 4")
  }
}
write_csv(bind_rows(rob), dir_out("tables", "robustness_leave_one_study_out.csv"))
