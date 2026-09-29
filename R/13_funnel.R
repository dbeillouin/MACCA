###############################################################################
# 13 – Funnel plots (Fig. S4)
#
# The funnel plots of the supplementary material were not produced by the
# workflow and therefore did not reflect the corrected sampling variances.
# They are now drawn here, from the same datasets and the same three-level
# models as the Egger-type tests, and follow the journal's figure rules
# (axes only, no frame, no grid, y-axis label horizontal above the axis,
# panel letters inside the panels).
###############################################################################
suppressPackageStartupMessages({
  library(metafor); library(ggplot2); library(patchwork)
})

rnd <- list(~ 1 | id_article, ~ 1 | id_experiment)

if (!exists("UNIT_RATE")) {
  UNIT_RATE <- expression("Storage rate (Mg C ha"^-1*" yr"^-1*")")
}

# theme and helper: defined in 10_figures.R, redefined here so that the script
# can also be run on its own with "only13"
if (!exists("theme_asd")) {
  theme_asd <- theme_classic(base_size = 8, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.4, colour = "black"),
          axis.ticks = element_line(linewidth = 0.4, colour = "black"),
          axis.text = element_text(colour = "black", size = 7),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 8, hjust = 0, margin = margin(b = 3)),
          plot.title.position = "panel",
          panel.grid = element_blank(), plot.margin = margin(4, 6, 4, 4))
}
if (!exists("save_fig")) {
  save_fig <- function(p, name, w, h) {
    ggsave(dir_out("figures", paste0(name, ".pdf")),  p, width = w, height = h, units = "mm")
    ggsave(dir_out("figures", paste0(name, ".tiff")), p, width = w, height = h, units = "mm",
           dpi = 600, compression = "lzw", bg = "white")
    ggsave(dir_out("figures", paste0(name, ".png")),  p, width = w, height = h, units = "mm",
           dpi = 300, bg = "white")
  }
}

panel <- function(nm, letter, xlab) {
  d <- read_csv(dir_out("data", paste0("analysis_dataset_", nm, ".csv")), show_col_types = FALSE)
  yy <- if (nm == "RR") d$yi else d$seq_rate
  vv <- if (nm == "RR") d$vi else d$seq_rate_vi
  keep <- !is.na(yy) & !is.na(vv) & vv > 0
  dd <- data.frame(y = yy[keep], v = vv[keep],
                   id_article = d$id_article[keep], id_experiment = d$id_experiment[keep])
  dd$se <- sqrt(dd$v)

  m <- rma.mv(y, v, random = rnd, data = dd, method = "REML")
  b <- as.numeric(m$b[1, 1])

  se_max <- max(dd$se)
  lim <- data.frame(se = seq(0, se_max * 1.02, length.out = 200))
  lim$lo <- b - 1.96 * lim$se
  lim$hi <- b + 1.96 * lim$se

  # panel letter: explicit coordinates, because the y axis is reversed
  x_rng <- range(c(dd$y, lim$lo, lim$hi), finite = TRUE)
  x_tag <- x_rng[1] + 0.02 * diff(x_rng)
  y_tag <- 0.02 * se_max

  ggplot() +
    geom_ribbon(data = lim, aes(y = se, xmin = lo, xmax = hi), fill = "grey92", colour = NA) +
    geom_vline(xintercept = b, linewidth = 0.4, colour = "grey40") +
    geom_point(data = dd, aes(x = y, y = se), size = 0.8, alpha = 0.55, colour = "#1B6CA8") +
    annotate("text", x = x_tag, y = y_tag, label = paste0("(", letter, ")"),
             hjust = 0, vjust = 0, size = 3, fontface = "bold") +
    scale_y_reverse(expand = expansion(mult = c(0.06, 0.04))) +
    coord_cartesian(xlim = x_rng) +
    labs(x = xlab, title = "Standard error") +
    theme_asd
}

pa <- panel("RR",           "a", "Log response ratio")
pb <- panel("storage_rate", "b", UNIT_RATE)

save_fig(pa | pb, "FigS_funnel", w = 174, h = 80)
message("[13] written: ", dir_out("figures", "FigS_funnel.png"))
