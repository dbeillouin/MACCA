###############################################################################
# 05 – Partial dependence by soil depth class (Fig. 4, Figs. S12-S13, 3.3-3.5)
#
# Changes vs. former 5_PDP_*.R:
#  * NO EXTRAPOLATION: within each depth class, curves are computed only for
#    the observations of that class and only over the range of the moderator
#    observed in that class (formerly every observation's ICE curve was drawn
#    over the full range and then labelled with its depth class, which produced
#    e.g. 45-75 cm predictions at 20 years without any such data);
#  * uncertainty = cluster bootstrap of the XGBoost model (studies resampled,
#    model refitted), instead of GAM confidence bands fitted on thousands of
#    PDP/ICE pseudo-observations;
#  * depth x moderator tests are run on the OBSERVED effect sizes with
#    three-level meta-regressions (natural splines), because p-values of GAMs
#    fitted to model predictions (p < 2e-16, p = 1) are not valid inference.
###############################################################################
suppressPackageStartupMessages({ library(xgboost); library(metafor); library(splines); library(ggplot2) })

MODERATORS <- c(time_since_conversion = "Time since conversion (years)",
                control_soc_mean_T_ha = "Initial soil organic carbon stock (Mg C ha-1)",
                temperature           = "Mean annual temperature (°C)",
                precipitation         = "Mean annual precipitation (mm)")
LAYERS <- list(`0-30 cm` = c(0, 30), `45-75 cm` = c(45, 75))   # pooled layers quoted in the text
GRID_N <- 40

pdp_on_rows <- function(model, X, rows, var, grid) {
  sapply(grid, function(g) { Xg <- X[rows, , drop = FALSE]; Xg[, var] <- g; mean(predict(model, Xg)) })
}

groups_for <- function(d) {
  cls <- split(seq_len(nrow(d)), d$depth_group, drop = TRUE)
  lay <- lapply(LAYERS, function(b) which(d$MEAN_depth >= b[1] & d$MEAN_depth < b[2]))
  c(cls, lay)
}

all_pdp <- list(); tests <- list()

for (nm in c("RR", "storage_rate")) {
  obj <- readRDS(dir_out("models", paste0("xgb_", nm, ".rds")))
  d <- obj$data; X <- obj$X; y <- d[[obj$y_name]]
  back <- if (nm == "RR") function(x) 100 * (exp(x) - 1) else identity
  grp <- groups_for(d)
  grp <- grp[sapply(grp, length) >= 2]

  # grids restricted to the observed range within each group
  grids <- lapply(setNames(names(MODERATORS), names(MODERATORS)), function(v)
    lapply(grp, function(rows) { r <- range(d[[v]][rows]); if (diff(r) == 0) NULL else seq(r[1], r[2], length.out = GRID_N) }))

  point <- map_dfr(names(MODERATORS), function(v) map_dfr(names(grp), function(gname) {
    g <- grids[[v]][[gname]]; if (is.null(g)) return(NULL)
    tibble(metric = nm, moderator = v, group = gname, x = g, n_obs = length(grp[[gname]]),
           n_studies = n_distinct(d$id_article[grp[[gname]]]),
           fit = pdp_on_rows(obj$model, X, grp[[gname]], v, g))
  }))

  message("   bootstrap PDP (", N_BOOT_PDP, " refits) for ", nm)
  boot <- map_dfr(seq_len(N_BOOT_PDP), function(b) {
    idx <- boot_rows_by_study(d$id_article)
    mb <- fit_xgb(X[idx, , drop = FALSE], y[idx], obj$nrounds)
    map_dfr(names(MODERATORS), function(v) map_dfr(names(grp), function(gname) {
      g <- grids[[v]][[gname]]; if (is.null(g)) return(NULL)
      tibble(moderator = v, group = gname, x = g, fit_b = pdp_on_rows(mb, X, grp[[gname]], v, g))
    }))
  })
  ci <- boot %>% group_by(moderator, group, x) %>%
    summarise(lwr = quantile(fit_b, .025), upr = quantile(fit_b, .975), .groups = "drop")
  pd <- point %>% left_join(ci, by = c("moderator", "group", "x")) %>%
    mutate(across(c(fit, lwr, upr), back, .names = "{.col}_bt"))
  all_pdp[[nm]] <- pd

  # ---- values quoted in the text (NA when outside the observed range) --------
  at <- function(v, gname, x0) {
    s <- pd %>% filter(moderator == v, group == gname)
    if (nrow(s) == 0 || x0 < min(s$x) || x0 > max(s$x))
      return(paste0("outside observed range (", ifelse(nrow(s), paste0(round(min(s$x), 1), "-", round(max(s$x), 1)), "no data"), ")"))
    f <- function(col) approx(s$x, s[[col]], x0)$y
    sprintf("%.2f [%.2f; %.2f]", f("fit_bt"), f("lwr_bt"), f("upr_bt"))
  }
  for (gname in names(LAYERS)) {
    for (tt in c(2, 7, 20)) record(paste0(nm, "_PDP_time_", tt, "y_", gname), at("time_since_conversion", gname, tt), "3.3, 4, Abstract")
    record(paste0(nm, "_PDP_initialSOC_25_", gname), at("control_soc_mean_T_ha", gname, 25), "3.4")
    s <- pd %>% filter(moderator == "control_soc_mean_T_ha", group == gname)
    if (nrow(s) > 2) {
      record(paste0(nm, "_PDP_slope_per_10MgC_initialSOC_", gname), 10 * coef(lm(fit_bt ~ x, data = s))[2], "3.4",
             if (nm == "RR") "percentage points per 10 Mg C ha-1" else "Mg C ha-1 yr-1 per 10 Mg C ha-1")
      zero <- s$x[which(diff(sign(s$fit_bt)) != 0)]
      record(paste0(nm, "_PDP_initialSOC_zero_crossing_", gname), if (length(zero)) round(zero[1], 1) else "none", "3.4, Abstract (~80)")
    }
  }

  # ---- tests on observed data: does depth modify the response curve? ---------
  dt <- d; dt$.y <- y; dt$.v <- if (nm == "RR") d$vi else d$seq_rate_vi
  for (v in names(MODERATORS)) {
    dt$.m <- dt[[v]]
    base <- function(mods) rma.mv(.y, .v, mods = mods, random = list(~ 1 | id_article, ~ 1 | id_experiment),
                                  data = dt, method = "ML")
    m_add <- tryCatch(base(~ ns(.m, df = 3) + MEAN_depth), error = function(e) NULL)
    m_int <- tryCatch(base(~ ns(.m, df = 3) * MEAN_depth), error = function(e) NULL)
    m_nodepth <- tryCatch(base(~ ns(.m, df = 3)), error = function(e) NULL)
    if (is.null(m_add) || is.null(m_int) || is.null(m_nodepth)) next
    a1 <- anova(m_nodepth, m_add); a2 <- anova(m_add, m_int)
    tests[[paste(nm, v)]] <- tibble(metric = nm, moderator = v,
                                    depth_additive_LRT = a1$LRT, depth_additive_p = a1$pval,
                                    depth_interaction_LRT = a2$LRT, depth_interaction_df = a2$parms.f - a2$parms.r,
                                    depth_interaction_p = a2$pval)
    record(paste0(nm, "_test_depth_x_", v),
           sprintf("additive: LRT=%.2f, p=%.3g; interaction: LRT=%.2f (df=%d), p=%.3g",
                   a1$LRT, a1$pval, a2$LRT, a2$parms.f - a2$parms.r, a2$pval), "3.3-3.5, 4")
  }
}

pdp_tab <- bind_rows(all_pdp)
write_csv(pdp_tab, dir_out("tables", "partial_dependence_by_depth.csv"))
write_csv(bind_rows(tests), dir_out("tables", "tests_depth_interactions_observed_data.csv"))

# =============================================================================
# Fig. 4 draft following the ASD artwork rules: 4 panels (a-d), no frame, only
# x and y axes, no grid, letters inside panels, horizontal y label above axis,
# depth classes distinguished by grey level and line type (no colour needed).
# =============================================================================
if (requireNamespace("patchwork", quietly = TRUE)) {
  cls_levels <- levels(cut(0, DEPTH_BREAKS, include.lowest = TRUE, right = FALSE))
  asd_theme <- theme_classic(base_size = 9, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.5, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)),
          plot.title.position = "panel", legend.position = "bottom",
          legend.title = element_blank(), panel.grid = element_blank())
  panel <- function(nm, v, letter, ylab) {
    s <- pdp_tab %>% filter(metric == nm, moderator == v, group %in% cls_levels) %>%
      mutate(group = factor(group, levels = cls_levels))
    obj <- readRDS(dir_out("models", paste0("xgb_", nm, ".rds")))
    pts <- obj$data %>% mutate(yobs = if (nm == "RR") 100 * (exp(yi) - 1) else seq_rate,
                               w = if (nm == "RR") 1 / sqrt(vi) else 1 / seq_rate_sd,
                               group = factor(as.character(depth_group), levels = cls_levels))
    ggplot(s, aes(x, fit_bt, group = group)) +
      geom_hline(yintercept = 0, linewidth = 0.3, linetype = 3) +
      geom_point(data = pts, aes(.data[[v]], yobs, size = w, shape = group), inherit.aes = FALSE,
                 colour = "grey70", alpha = 0.35, show.legend = FALSE) +
      geom_ribbon(aes(ymin = lwr_bt, ymax = upr_bt, fill = group), alpha = 0.15, show.legend = FALSE) +
      geom_line(aes(linetype = group, colour = group), linewidth = 0.8) +
      geom_rug(data = pts, aes(.data[[v]], colour = group), sides = "b", inherit.aes = FALSE,
               length = unit(0.02, "npc"), show.legend = FALSE) +
      scale_colour_grey(start = 0.75, end = 0) + scale_fill_grey(start = 0.75, end = 0) +
      scale_size(range = c(0.3, 1.8)) +
      annotate("text", x = -Inf, y = Inf, label = paste0("(", letter, ")"), hjust = -0.4, vjust = 1.4, size = 3.2) +
      labs(x = MODERATORS[[v]], title = ylab) + asd_theme
  }
  fig4 <- (panel("RR", "time_since_conversion", "a", "Response ratio (%)") +
           panel("storage_rate", "time_since_conversion", "b", expression("Storage rate (Mg C ha"^-1*" yr"^-1*")"))) /
          (panel("RR", "control_soc_mean_T_ha", "c", "Response ratio (%)") +
           panel("storage_rate", "control_soc_mean_T_ha", "d", expression("Storage rate (Mg C ha"^-1*" yr"^-1*")"))) +
          patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
  ggsave(dir_out("figures", "Fig4_draft.pdf"), fig4, width = 174, height = 160, units = "mm")
  ggsave(dir_out("figures", "Fig4_draft.png"), fig4, width = 174, height = 160, units = "mm", dpi = 600)
}
