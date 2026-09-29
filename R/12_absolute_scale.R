###############################################################################
# 12 – Absolute scale: is the initial-stock effect a normalisation artefact?
#
# Both response variables are normalised by the control stock: the response
# ratio divides by it, the storage rate subtracts it and then divides by the
# duration. Slessarev et al. (2023) show that this creates a negative
# relationship with the initial stock even when none exists, and recommend
# (i) using the absolute difference SOCtreatment - SOCcontrol, (ii) treating
# time as a covariate rather than a divisor, and (iii) correcting the slope
# for regression to the mean with Blomqvist's (1977) formula:
#
#        beta_corrected = (beta_observed + lambda) / (1 - lambda)
#        lambda = s2_u / s2_z
#          s2_u = error variance of the initial value = mean(SD^2 / n)
#          s2_z = variance of the initial values across observations
#
# Everything quoted in section 3.3 of the manuscript is produced here.
###############################################################################
suppressPackageStartupMessages({ library(metafor); library(splines) })

rnd <- list(~ 1 | id_article, ~ 1 | id_experiment)

d <- read_csv(dir_out("data", "analysis_dataset_RR.csv"), show_col_types = FALSE) %>%
  filter(!is.na(treatment_soc_mean_T_ha), !is.na(control_soc_mean_T_ha),
         !is.na(treatment_soc_sd_T_ha),   !is.na(control_soc_sd_T_ha),
         !is.na(control_replicate_nb),    !is.na(treatment_replicate_nb),
         control_replicate_nb > 0, treatment_replicate_nb > 0,
         !is.na(temperature), !is.na(time_since_conversion)) %>%
  mutate(md    = treatment_soc_mean_T_ha - control_soc_mean_T_ha,
         v_md  = treatment_soc_sd_T_ha^2 / treatment_replicate_nb +
                 control_soc_sd_T_ha^2   / control_replicate_nb,
         s2u   = control_soc_sd_T_ha^2   / control_replicate_nb,
         soc_c = control_soc_mean_T_ha - mean(control_soc_mean_T_ha),
         t_c   = temperature            - mean(temperature),
         yr_c  = time_since_conversion  - mean(time_since_conversion)) %>%
  filter(v_md > 0)

message(sprintf("[12] absolute scale: %d observations / %d studies",
                nrow(d), n_distinct(d$id_article)))
record("abs_n_obs_studies",
       sprintf("%d observations / %d studies", nrow(d), n_distinct(d$id_article)), "3.3")

fitml <- function(x, mods = NULL) {
  if (is.null(mods)) rma.mv(md, v_md, random = rnd, data = x, method = "ML")
  else               rma.mv(md, v_md, mods = mods, random = rnd, data = x, method = "ML")
}
lrt <- function(x, mods0, mods1) {
  m0 <- tryCatch(fitml(x, mods0), error = function(e) NULL)
  m1 <- tryCatch(fitml(x, mods1), error = function(e) NULL)
  if (is.null(m0) || is.null(m1)) return(list(LRT = NA, df = NA, p = NA, m1 = NULL))
  a <- anova(m0, m1)
  list(LRT = a$LRT, df = a$parms.f - a$parms.r, p = a$pval, m1 = m1)
}
lambda <- function(x) mean(x$s2u) / var(x$control_soc_mean_T_ha)
blomqvist <- function(b, lam) (b + lam) / (1 - lam)

rows <- list()
add  <- function(...) rows[[length(rows) + 1L]] <<- tibble(...)

# ------------------------------------------------------- 1. overall difference
m_overall <- rma.mv(md, v_md, random = rnd, data = d, method = "REML")
record("abs_overall_difference",
       sprintf("%.2f [%.2f; %.2f] Mg C/ha, p = %.3g",
               m_overall$b[1, 1], m_overall$ci.lb[1], m_overall$ci.ub[1], m_overall$pval[1]),
       "3.3")
add(quantity = "Mean absolute difference (Mg C/ha)",
    estimate = sprintf("%.2f [%.2f; %.2f]", m_overall$b[1, 1], m_overall$ci.lb[1], m_overall$ci.ub[1]),
    test = sprintf("p = %.3g", m_overall$pval[1]), lambda = NA_character_, corrected = NA_character_)

# --------------------------------------------- 2. initial stock, absolute scale
for (spec in list(list("linear", ~ soc_c), list("spline", ~ ns(soc_c, 3)))) {
  a <- lrt(d, NULL, spec[[2]])
  record(sprintf("abs_initialSOC_%s", spec[[1]]),
         sprintf("LRT=%.2f, df=%d, p=%.3g", a$LRT, a$df, a$p), "3.3")
  add(quantity = sprintf("Absolute difference ~ initial SOC (%s)", spec[[1]]),
      estimate = if (spec[[1]] == "linear" && !is.null(a$m1))
                   sprintf("%+.4f [%+.4f; %+.4f]", a$m1$b[2, 1], a$m1$ci.lb[2], a$m1$ci.ub[2]) else "-",
      test = sprintf("LRT = %.2f, df = %d, p = %.3g", a$LRT, a$df, a$p),
      lambda = NA_character_, corrected = NA_character_)
}
a_adj <- lrt(d, ~ yr_c + t_c, ~ yr_c + t_c + soc_c)
record("abs_initialSOC_given_time_temperature",
       sprintf("LRT=%.2f, df=%d, p=%.3g; slope=%+.4f [%+.4f; %+.4f]",
               a_adj$LRT, a_adj$df, a_adj$p,
               a_adj$m1$b["soc_c", 1],
               a_adj$m1$ci.lb[which(rownames(a_adj$m1$b) == "soc_c")],
               a_adj$m1$ci.ub[which(rownames(a_adj$m1$b) == "soc_c")]), "3.3")
add(quantity = "Absolute difference ~ initial SOC, given time and temperature",
    estimate = sprintf("%+.4f [%+.4f; %+.4f]", a_adj$m1$b["soc_c", 1],
                       a_adj$m1$ci.lb[which(rownames(a_adj$m1$b) == "soc_c")],
                       a_adj$m1$ci.ub[which(rownames(a_adj$m1$b) == "soc_c")]),
    test = sprintf("LRT = %.2f, df = %d, p = %.3g", a_adj$LRT, a_adj$df, a_adj$p),
    lambda = NA_character_, corrected = NA_character_)

# ------------------------------------ 3. interaction with temperature, absolute
a_int <- lrt(d, ~ soc_c + t_c, ~ soc_c * t_c)
record("abs_initialSOC_x_temperature",
       sprintf("LRT=%.2f, df=%d, p=%.3g", a_int$LRT, a_int$df, a_int$p), "3.3")
add(quantity = "Absolute difference ~ initial SOC x temperature",
    estimate = "-", test = sprintf("LRT = %.2f, df = %d, p = %.3g", a_int$LRT, a_int$df, a_int$p),
    lambda = NA_character_, corrected = NA_character_)

# ------------------------------------------- 4. time since conversion, absolute
a_yr <- lrt(d, ~ soc_c + t_c, ~ soc_c + t_c + ns(yr_c, 3))
record("abs_time_spline_given_stock_temperature",
       sprintf("LRT=%.2f, df=%d, p=%.3g", a_yr$LRT, a_yr$df, a_yr$p), "3.3")
a_yrl <- lrt(d, ~ soc_c + t_c, ~ soc_c + t_c + yr_c)
record("abs_time_linear_given_stock_temperature",
       sprintf("LRT=%.2f, df=%d, p=%.3g; slope=%+.4f Mg C/ha/yr",
               a_yrl$LRT, a_yrl$df, a_yrl$p, a_yrl$m1$b["yr_c", 1]), "3.3")

studies <- unique(d$id_article); pv <- rep(NA_real_, length(studies))
for (i in seq_along(studies)) {
  x <- d[d$id_article != studies[i], , drop = FALSE]
  pv[i] <- tryCatch(lrt(x, ~ soc_c + t_c, ~ soc_c + t_c + ns(yr_c, 3))$p,
                    error = function(e) NA_real_)
}
k <- which.max(pv)
record("abs_time_spline_LOSO",
       sprintf("max p without one study = %.3g (study %s); studies whose removal makes p > 0.05: %d",
               pv[k], studies[k], sum(pv > 0.05, na.rm = TRUE)), "3.3")
add(quantity = "Absolute difference ~ time since conversion (spline), given stock and temperature",
    estimate = sprintf("%+.4f Mg C/ha/yr (linear term)", a_yrl$m1$b["yr_c", 1]),
    test = sprintf("LRT = %.2f, df = %d, p = %.3g", a_yr$LRT, a_yr$df, a_yr$p),
    lambda = NA_character_,
    corrected = sprintf("leave-one-study-out maximum p = %.3g", max(pv, na.rm = TRUE)))

# -------------------------- 5. Blomqvist correction, overall and by sub-group
subsets <- list(
  "All observations"        = d,
  "Control-impact designs"  = d[d$Grouped_Design == "Control-Impact Designs", , drop = FALSE],
  "Before-after designs"    = d[d$Grouped_Design == "Before-After Designs",   , drop = FALSE],
  "Randomized designs"      = d[d$Grouped_Design == "Randomized Designs",     , drop = FALSE],
  "Surface layers (0-30 cm)"= d[d$soil_depth_start == 0 & d$soil_depth_end <= 30, , drop = FALSE],
  "At least 3 replicates"   = d[d$control_replicate_nb >= 3, , drop = FALSE])

for (nm in names(subsets)) {
  x <- subsets[[nm]]
  if (nrow(x) < 30 || n_distinct(x$id_article) < 5) next
  m <- tryCatch(rma.mv(md, v_md, mods = ~ control_soc_mean_T_ha, random = rnd,
                       data = x, method = "REML"), error = function(e) NULL)
  if (is.null(m)) next
  lam <- lambda(x); b <- m$b[2, 1]
  record(paste0("abs_blomqvist_", gsub("[^A-Za-z0-9]+", "_", nm)),
         sprintf("n=%d/%d studies; observed %+.4f [%+.4f; %+.4f], p=%.3g; lambda=%.4f; corrected %+.4f",
                 nrow(x), n_distinct(x$id_article), b, m$ci.lb[2], m$ci.ub[2], m$pval[2],
                 lam, blomqvist(b, lam)), "3.3")
  add(quantity = paste0("Slope on initial SOC: ", nm),
      estimate = sprintf("%+.4f [%+.4f; %+.4f] (n = %d, %d studies)",
                         b, m$ci.lb[2], m$ci.ub[2], nrow(x), n_distinct(x$id_article)),
      test = sprintf("p = %.3g", m$pval[2]),
      lambda = sprintf("%.4f", lam),
      corrected = sprintf("%+.4f", blomqvist(b, lam)))
}

write_csv(bind_rows(rows), dir_out("tables", "absolute_scale.csv"))
message("[12] written: ", dir_out("tables", "absolute_scale.csv"))
