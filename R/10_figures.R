###############################################################################
# 10 – Figures (main text and supplementary)
#  * colour-blind-safe palettes (Okabe–Ito, viridis); captions never rely on
#    colour alone (line types, symbols and panel letters carry the information)
#  * ASD artwork rules for graphs: no frame, x and y axes only, no grid, panel
#    letters inside panels below the top of the y axis, horizontal y label above
#  * PDF (vector) + TIFF 600 dpi + PNG; widths 84 mm or 174 mm
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(metafor); library(splines); library(patchwork) })

OI <- c(blue = "#0072B2", orange = "#E69F00", green = "#009E73", vermillion = "#D55E00",
        purple = "#CC79A7", sky = "#56B4E9", yellow = "#F0E442", grey = "#7F7F7F")
COL_METRIC <- c("Response ratio" = OI[["blue"]], "Storage rate" = OI[["vermillion"]])
COL_SYSTEM <- c("Silvopasture" = OI[["green"]], "Shaded perennial crops" = OI[["orange"]],
                "Alley cropping/hedgerow" = OI[["blue"]], "Multistrata" = OI[["purple"]], "Other" = OI[["grey"]])
SYS_LAB <- c("SilvoPasture" = "Silvopasture", "Shaded perennial" = "Shaded perennial crops",
             "Alley/Hedgerow" = "Alley cropping/hedgerow", "Multistrata system" = "Multistrata", "Other" = "Other")
LAB <- c(control_soc_mean_T_ha = "Initial SOC stock", time_since_conversion = "Time since conversion",
         precipitation = "Precipitation", temperature = "Temperature", MEAN_depth = "Soil depth",
         NEW_treatment_type2 = "Agroforestry type", History_reclass = "Previous land use",
         diff_species_class = "Tree species added", main_culture2 = "Main crop", Grouped_Design = "Study design")
SHORT <- c(control_soc_mean_T_ha = "Initial SOC", time_since_conversion = "Time", precipitation = "Precipitation",
           temperature = "Temperature", MEAN_depth = "Depth", NEW_treatment_type2 = "AF type",
           History_reclass = "Land use", diff_species_class = "Species", main_culture2 = "Crop", Grouped_Design = "Design")
BIOPHYS <- c("control_soc_mean_T_ha", "time_since_conversion", "precipitation", "temperature", "MEAN_depth")
UNIT_RATE <- expression("Storage rate (Mg C ha"^-1*" yr"^-1*")")
UNIT_SOC  <- expression("Initial soil organic carbon stock (Mg C ha"^-1*")")

theme_asd <- theme_classic(base_size = 8, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        axis.text = element_text(colour = "black", size = 7),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 8, hjust = 0, margin = margin(b = 3)),
        plot.title.position = "panel",
        legend.title = element_text(size = 7), legend.text = element_text(size = 7),
        legend.key.height = unit(3.2, "mm"), legend.background = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 8, face = "bold"),
        panel.grid = element_blank(), plot.margin = margin(4, 6, 4, 4))
tag <- function(letter) annotate("text", x = -Inf, y = Inf, label = paste0("(", letter, ")"),
                                 hjust = -0.35, vjust = 1.3, size = 3, fontface = "bold")
# panel letter for faceted plots: drawn in the first facet only
tag_facet <- function(letter) geom_text(data = data.frame(grp = factor("Biophysical", levels = c("Biophysical", "System")),
                                                          lab = paste0("(", letter, ")")),
                                        aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
                                        hjust = -0.35, vjust = 1.3, size = 3, fontface = "bold")
save_fig <- function(p, name, w, h) {
  dev_pdf <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(dir_out("figures", paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = dev_pdf)
  ggsave(dir_out("figures", paste0(name, ".tiff")), p, width = w, height = h, units = "mm", dpi = 600, compression = "lzw", bg = "white")
  ggsave(dir_out("figures", paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
}

rr <- read_csv(dir_out("data", "analysis_dataset_RR.csv"), show_col_types = FALSE)
sq <- read_csv(dir_out("data", "analysis_dataset_storage_rate.csv"), show_col_types = FALSE)
rnd <- list(~ 1 | id_article, ~ 1 | id_experiment)
pct <- function(x) 100 * (exp(x) - 1)

# =============================================================================
# Fig. 2 – map: countries shaded by number of studies, sites sized by observations,
#          zoom on the Costa Rica–Nicaragua cluster
# =============================================================================
if (requireNamespace("maps", quietly = TRUE)) {
  getll <- function(d) d %>% mutate(
    lon_ = if ("lon" %in% names(d)) lon else rescale_coord(`X_(WGS84)`, -120, -30),
    lat_ = if ("lat" %in% names(d)) lat else rescale_coord(`Y_(WGS84)`, -60, 35))
  s_all <- bind_rows(getll(rr) %>% transmute(id_article, country, lon_, lat_),
                     getll(sq) %>% transmute(id_article, country, lon_, lat_))
  record("fig2_obs_without_coordinates", sum(is.na(s_all$lat_) | is.na(s_all$lon_)), "Fig. 2 caption")
  s <- s_all %>% filter(!is.na(lat_), !is.na(lon_)) %>%
    mutate(lon_ = round(lon_, 2), lat_ = round(lat_, 2)) %>%
    group_by(lon_, lat_) %>% summarise(n_obs = n(), .groups = "drop") %>% arrange(desc(n_obs))
  record("fig2_n_sites", nrow(s), "Fig. 2 caption")
  n_country <- s_all %>% distinct(country, id_article) %>% count(country, name = "n_studies")
  write_csv(n_country, dir_out("tables", "studies_by_country.csv"))

  # countries covered by the search string (Table S1); others are outside the search scope
  SCOPE <- c("Mexico", "Belize", "Guatemala", "Honduras", "El Salvador", "Nicaragua", "Costa Rica", "Panama",
             "Cuba", "Jamaica", "Bahamas", "Cayman Islands", "Turks and Caicos", "Haiti", "Dominican Republic", "Puerto Rico",
             "Virgin Islands", "Anguilla", "Saint-Martin", "Saint-Barthelemy", "Antigua", "Barbuda", "Aruba", "Barbados",
             "Bonaire", "Curacao", "Dominica", "Grenada", "Guadeloupe", "Martinique", "Montserrat", "Saint Eustatius",
             "Saint Kitts", "Nevis", "Saint Lucia", "Saint Vincent", "Trinidad", "Tobago",
             "Venezuela", "Colombia", "Ecuador", "Peru")
  hires <- requireNamespace("mapdata", quietly = TRUE)
  if (hires) suppressPackageStartupMessages(library(mapdata))
  world <- map_data(if (hires) "worldHires" else "world", xlim = c(-110, -50), ylim = c(-25, 35)) %>%
    left_join(n_country, by = c("region" = "country")) %>%
    mutate(n_studies = ifelse(is.na(n_studies), 0, n_studies),
           cls = case_when(!region %in% SCOPE ~ "out", n_studies == 0 ~ "c0", n_studies <= 2 ~ "c1",
                           n_studies <= 5 ~ "c3", n_studies <= 10 ~ "c6", TRUE ~ "c10"),
           cls = factor(cls, levels = c("c0", "c1", "c3", "c6", "c10", "out")))
  COL_N <- c(c0 = "#F7F4EA", c1 = "#D4E9C8", c3 = "#95CDA8", c6 = "#4DA88F", c10 = "#1F7872", out = "#DADADA")
  LAB_N <- c("0", "1\u20132", "3\u20135", "6\u201310", "> 10", "Outside search scope")
  SEA <- "#DCEAF3"; BORDER <- "#8F8F8F"
  ZOOM <- c(xmin = -86.4, xmax = -82.3, ymin = 8.3, ymax = 13.4)
  countries <- tibble::tribble(~lon, ~lat, ~lab,
    -100.5, 22.0, "MEXICO", -86.9, 15.4, "HONDURAS", -78.2, 21.7, "CUBA",
    -71.6, 3.0, "COLOMBIA", -65.5, 7.2, "VENEZUELA", -76.4, -2.9, "ECUADOR", -74.3, -11.0, "PERU", -58.5, -6.0, "BRAZIL")
  land <- function(lw) list(
    geom_polygon(data = world, aes(long, lat, group = group, fill = cls), colour = BORDER, linewidth = lw),
    scale_fill_manual(values = COL_N, labels = LAB_N, drop = FALSE, name = "Number of studies per country",
                      guide = guide_legend(order = 1, nrow = 1, override.aes = list(colour = "#8F8F8F", linewidth = 0.2))))
  pts <- list(
    geom_point(data = s, aes(lon_, lat_, size = n_obs), shape = 21, fill = "#2B2B2B", colour = "white", stroke = 0.45, alpha = 0.85),
    scale_size_area(max_size = 7, limits = c(1, max(s$n_obs)), breaks = c(2, 10, 25, 50), name = "Observations per site",
                    guide = guide_legend(order = 2, nrow = 1)))
  deg_x <- function(x) paste0(abs(x), "\u00b0W")
  deg_y <- function(y) ifelse(y < 0, paste0(abs(y), "\u00b0S"), ifelse(y > 0, paste0(y, "\u00b0N"), "0\u00b0"))
  map_theme <- theme_asd + theme(panel.background = element_rect(fill = SEA, colour = NA),
                                 panel.border = element_rect(fill = NA, colour = "#6E6E6E", linewidth = 0.35),
                                 panel.grid.major = element_line(colour = "white", linewidth = 0.25),
                                 axis.line = element_blank(), axis.text = element_text(size = 6.5, colour = "#3A3A3A"),
                                 legend.key = element_rect(fill = NA, colour = NA))
  lab_box <- function(x, y, l) annotate("label", x = x, y = y, label = l, size = 3, fontface = "bold", hjust = 0, vjust = 1,
                                        label.size = 0, fill = "white")
  main_map <- ggplot() + land(0.18) +
    geom_text(data = countries, aes(lon, lat, label = lab), size = 2.1, colour = "#4F4F4F", fontface = "italic") +
    annotate("rect", xmin = ZOOM[["xmin"]], xmax = ZOOM[["xmax"]], ymin = ZOOM[["ymin"]], ymax = ZOOM[["ymax"]],
             fill = NA, colour = "#2B2B2B", linewidth = 0.35, linetype = 2) +
    annotate("text", x = ZOOM[["xmax"]] + 0.5, y = ZOOM[["ymax"]] + 0.2, label = "(b)", size = 2.6, fontface = "bold", hjust = 0, vjust = 0) +
    pts +
    scale_x_continuous(breaks = seq(-100, -60, 10), labels = deg_x) + scale_y_continuous(breaks = seq(-10, 20, 10), labels = deg_y) +
    coord_quickmap(xlim = c(-104, -56), ylim = c(-14, 26), expand = FALSE) +
    labs(x = NULL, title = NULL) + lab_box(-103.3, 25.2, "(a)") + map_theme
  zoom_map <- ggplot() + land(0.25) +
    annotate("text", x = c(-84.35, -85.35), y = c(9.25, 12.35), label = c("COSTA RICA", "NICARAGUA"),
             size = 2.1, colour = c("#F7F7F7", "#3A3A3A"), fontface = "italic") +
    pts +
    scale_x_continuous(breaks = seq(-86, -83, 1), labels = deg_x) + scale_y_continuous(breaks = seq(9, 13, 1), labels = deg_y) +
    coord_quickmap(xlim = c(ZOOM[["xmin"]], ZOOM[["xmax"]]), ylim = c(ZOOM[["ymin"]], ZOOM[["ymax"]]), expand = FALSE) +
    labs(x = NULL, title = NULL) + lab_box(ZOOM[["xmin"]] + 0.08, ZOOM[["ymax"]] - 0.08, "(b)") +
    map_theme + theme(legend.position = "none")
  fig2 <- (main_map | zoom_map) + plot_layout(widths = c(1.6, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
          legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 7),
          legend.spacing.y = unit(1, "mm"), legend.margin = margin(0, 0, 0, 0))
  save_fig(fig2, "Fig2_map", 174, 122)
}

# =============================================================================
# Fig. 3 – importance (a), share of interactions (b), pairwise interactions (c, d)
# =============================================================================
imp <- bind_rows(read_csv(dir_out("tables", "importance_RR.csv"), show_col_types = FALSE) %>% mutate(metric = "Response ratio"),
                 read_csv(dir_out("tables", "importance_storage_rate.csv"), show_col_types = FALSE) %>% mutate(metric = "Storage rate"))
ord <- imp %>% group_by(variable) %>% summarise(m = mean(share_mean), .groups = "drop") %>%
  mutate(grp = ifelse(variable %in% BIOPHYS, "Biophysical", "System")) %>% arrange(grp, desc(m))
fac <- function(d) d %>% mutate(var = factor(LAB[variable], levels = rev(LAB[ord$variable])),
                                grp = factor(ifelse(variable %in% BIOPHYS, "Biophysical", "System"), levels = c("Biophysical", "System")),
                                metric = factor(metric, levels = names(COL_METRIC)))
imp <- fac(imp)
pd <- position_dodge(width = -0.55)
p3a <- ggplot(imp, aes(share_mean, var, colour = metric)) +
  geom_linerange(aes(xmin = share_q17, xmax = share_q83), position = pd, linewidth = 0.6) +
  geom_point(position = pd, size = 1.9) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = COL_METRIC, name = NULL) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Share of total importance (Gain, %)") + tag_facet("a") +
  theme_asd + theme(legend.position = "top", legend.justification = "left",
                    strip.text.y = element_blank(), panel.spacing = unit(3, "mm"))
ish <- fac(bind_rows(read_csv(dir_out("tables", "interaction_share_RR.csv"), show_col_types = FALSE) %>% mutate(metric = "Response ratio"),
                     read_csv(dir_out("tables", "interaction_share_storage_rate.csv"), show_col_types = FALSE) %>% mutate(metric = "Storage rate")))
p3b <- ggplot(ish, aes(100 * interaction_share, var, colour = metric)) +
  geom_vline(xintercept = 50, linetype = 3, linewidth = 0.3) +
  geom_point(position = pd, size = 1.9, show.legend = FALSE) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = COL_METRIC) + scale_x_continuous(limits = c(0, 100)) +
  labs(x = "Contribution due to interactions (%)") + tag_facet("b") +
  theme_asd + theme(axis.text.y = element_blank(), strip.text.y = element_text(angle = -90, size = 7.5),
                    panel.spacing = unit(3, "mm"))
heat <- function(file, letter, ttl, show_y = TRUE) {
  M <- read.csv(dir_out("tables", file), row.names = 1, check.names = FALSE)
  lev <- ord$variable
  df <- as.data.frame(as.table(as.matrix(M))) %>% setNames(c("focal", "partner", "pct")) %>%
    mutate(focal = factor(SHORT[as.character(focal)], levels = rev(SHORT[lev])),
           partner = factor(SHORT[as.character(partner)], levels = SHORT[lev]),
           dark = !is.na(pct) & pct > 35)
  ggplot(df, aes(partner, focal, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(pct), "", round(pct)), colour = dark), size = 2.1, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 60), oob = scales::squish, na.value = "white",
                         name = "Share of the moderator's interactions (%)") +
    scale_y_discrete(expand = expansion(add = c(0.5, 1.4))) +
    labs(x = "Partner moderator", title = ttl) + tag(letter) + theme_asd +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank(),
          legend.position = "bottom", legend.key.width = unit(12, "mm"), legend.key.height = unit(2.6, "mm"),
          legend.title = element_text(size = 7, vjust = 0.9),
          axis.text.y = if (show_y) element_text(colour = "black", size = 7) else element_blank())
}
fig3 <- (p3a | p3b) / ((heat("interaction_pairs_rowpct_RR.csv", "c", "Response ratio") |
                        heat("interaction_pairs_rowpct_storage_rate.csv", "d", "Storage rate", show_y = FALSE)) +
                       plot_layout(guides = "collect") & theme(legend.position = "bottom")) +
  plot_layout(heights = c(1, 1.15))
save_fig(fig3, "Fig3_importance_interactions", 174, 195)

# =============================================================================
# Fig. 4 – drivers estimated on observed effect sizes
# =============================================================================
prep <- function(d, y, v) d %>% filter(!is.na(control_soc_mean_T_ha), !is.na(temperature), !is.na(time_since_conversion)) %>%
  mutate(.y = .data[[y]], .v = .data[[v]], w = 1 / sqrt(.v))
dr <- prep(rr, "yi", "vi"); ds <- prep(sq, "seq_rate", "seq_rate_vi")
ylim_q <- function(x) quantile(x, c(0.01, 0.99), na.rm = TRUE)
PT <- "#8C8C8C"

m_a <- rma.mv(.y, .v, mods = ~ ns(control_soc_mean_T_ha, 3), random = rnd, data = dr, method = "REML")
ga <- seq(min(dr$control_soc_mean_T_ha), max(dr$control_soc_mean_T_ha), length.out = 100)
pa <- predict(m_a, newmods = unname(as.matrix(predict(ns(dr$control_soc_mean_T_ha, 3), ga))))
fa <- tibble(x = ga, fit = pct(pa$pred), lb = pct(pa$ci.lb), ub = pct(pa$ci.ub))
p4a <- ggplot() +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
  geom_point(data = dr, aes(control_soc_mean_T_ha, pct(.y), size = w), shape = 16, colour = PT, alpha = 0.35) +
  geom_ribbon(data = fa, aes(x, ymin = lb, ymax = ub), fill = OI[["blue"]], alpha = 0.22) +
  geom_line(data = fa, aes(x, fit), colour = OI[["blue"]], linewidth = 0.9) +
  geom_rug(data = dr, aes(control_soc_mean_T_ha), sides = "b", length = unit(1.2, "mm"), linewidth = 0.2, alpha = 0.5) +
  scale_size(range = c(0.3, 2.5), guide = "none") + coord_cartesian(ylim = ylim_q(pct(dr$.y))) +
  labs(x = UNIT_SOC, title = "Response ratio (%)") + tag("a") + theme_asd

dr <- dr %>% mutate(soc_c = (control_soc_mean_T_ha - 50) / 10, t_c = temperature - 25)
m_b <- rma.mv(.y, .v, mods = ~ soc_c * t_c, random = rnd, data = dr, method = "REML")
fb <- map_dfr(c(22, 28), function(T) {
  sub <- if (T < 25) dr$control_soc_mean_T_ha[dr$temperature <= 24] else dr$control_soc_mean_T_ha[dr$temperature >= 26]
  qq <- quantile(sub, c(0.10, 0.90)); g <- seq(qq[1], min(qq[2], 100), length.out = 80)
  p <- predict(m_b, newmods = cbind((g - 50) / 10, T - 25, (g - 50) / 10 * (T - 25)))
  tibble(x = g, fit = pct(p$pred), lb = pct(p$ci.lb), ub = pct(p$ci.ub), temp = paste0(T, " \u00b0C"))
})
COL_T <- setNames(c(OI[["blue"]], OI[["vermillion"]]), c("22 \u00b0C", "28 \u00b0C"))
p4b <- ggplot(fb, aes(x, fit, linetype = temp, colour = temp, fill = temp)) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.18, colour = NA) + geom_line(linewidth = 0.9) +
  scale_colour_manual(values = COL_T, name = NULL) + scale_fill_manual(values = COL_T, name = NULL) +
  scale_linetype_manual(values = c(2, 1), name = NULL) +
  labs(x = UNIT_SOC, title = "Predicted response ratio (%)") + tag("b") + theme_asd +
  theme(legend.position = c(0.8, 0.86), legend.key.width = unit(7, "mm"))

m_c <- rma.mv(.y, .v, mods = ~ ns(time_since_conversion, 3), random = rnd, data = ds, method = "REML")
gc <- seq(min(ds$time_since_conversion), max(ds$time_since_conversion), length.out = 120)
pc <- predict(m_c, newmods = unname(as.matrix(predict(ns(ds$time_since_conversion, 3), gc))))
fc <- tibble(x = gc, fit = pc$pred, lb = pc$ci.lb, ub = pc$ci.ub)
p4c <- ggplot() +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
  geom_point(data = ds, aes(time_since_conversion, .y, size = w), shape = 16, colour = PT, alpha = 0.35) +
  geom_ribbon(data = fc, aes(x, ymin = lb, ymax = ub), fill = OI[["vermillion"]], alpha = 0.22) +
  geom_line(data = fc, aes(x, fit), colour = OI[["vermillion"]], linewidth = 0.9) +
  geom_rug(data = ds, aes(time_since_conversion), sides = "b", length = unit(1.2, "mm"), linewidth = 0.2, alpha = 0.5) +
  scale_x_sqrt(breaks = c(1, 5, 10, 20, 40, 70)) + scale_size(range = c(0.3, 2.5), guide = "none") +
  coord_cartesian(ylim = ylim_q(ds$.y)) +
  labs(x = "Time since conversion (years, square-root scale)", title = UNIT_RATE) + tag("c") + theme_asd

mt <- read_csv(dir_out("tables", "moderator_tests_observed_data.csv"), show_col_types = FALSE) %>%
  filter(moderator %in% c("all_biophysical", "all_management")) %>%
  mutate(group = ifelse(moderator == "all_biophysical", "Biophysical", "System"),
         metric = ifelse(metric == "RR", "Response ratio", "Storage rate"),
         lab = paste0(sprintf("%.1f%%", heterogeneity_explained_total_pct), "\n",
                      ifelse(p < 0.001, "p < 0.001", ifelse(p < 0.01, sprintf("p = %.3f", p), sprintf("p = %.2f", p)))))
pdd <- position_dodge(width = 0.5)
p4d <- ggplot(mt, aes(metric, heterogeneity_explained_total_pct, colour = group)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "black") +
  geom_linerange(aes(ymin = 0, ymax = heterogeneity_explained_total_pct), position = pdd, linewidth = 1.1) +
  geom_point(position = pdd, size = 2.6) +
  geom_text(aes(label = lab, y = heterogeneity_explained_total_pct + 1.6), position = pdd, size = 2.3, vjust = 0, show.legend = FALSE) +
  scale_colour_manual(values = c(Biophysical = OI[["green"]], System = OI[["grey"]]), name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.18))) +
  labs(x = NULL, title = "Heterogeneity explained (%)") + tag("d") + theme_asd +
  theme(legend.position = c(0.28, 0.85))
save_fig((p4a | p4b) / (p4c | p4d), "Fig4_drivers_observed", 174, 150)

# =============================================================================
# Supplementary figures
# =============================================================================
cvp <- bind_rows(read_csv(dir_out("tables", "cv_performance_RR.csv"), show_col_types = FALSE) %>% mutate(metric = "Response ratio"),
                 read_csv(dir_out("tables", "cv_performance_storage_rate.csv"), show_col_types = FALSE) %>% mutate(metric = "Storage rate")) %>%
  select(metric, R2_random, R2_profile, R2_grouped) %>%
  pivot_longer(-metric, names_to = "scheme", values_to = "R2") %>%
  mutate(scheme = factor(scheme, levels = c("R2_random", "R2_profile", "R2_grouped"),
                         labels = c("Random observations", "Soil profiles", "Studies")))
pS5 <- ggplot(cvp, aes(scheme, R2, fill = metric)) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
  geom_boxplot(width = 0.6, outlier.size = 0.6, linewidth = 0.3, alpha = 0.85) +
  scale_fill_manual(values = COL_METRIC, name = NULL) + facet_wrap(~ metric) +
  labs(x = "Units withheld in cross-validation", title = expression("Out-of-fold R"^2)) +
  theme_asd + theme(legend.position = "none")
save_fig(pS5, "FigS_cross_validation", 174, 80)

if (file.exists(dir_out("tables", "robustness_leave_one_study_out.csv"))) {
  lo <- read_csv(dir_out("tables", "robustness_leave_one_study_out.csv"), show_col_types = FALSE) %>%
    mutate(metric = factor(ifelse(metric == "RR", "Response ratio", "Storage rate"), levels = names(COL_METRIC)),
           test = factor(test, levels = c("SOC x temperature", "depth x precipitation", "temperature linear", "temperature spline", "initial SOC linear"),
                         labels = c("Initial SOC \u00d7 temperature", "Depth \u00d7 precipitation", "Temperature (linear)",
                                    "Temperature (spline)", "Initial SOC (linear)")))
  pS6 <- ggplot(lo, aes(p, test, colour = metric)) +
    geom_vline(xintercept = 0.05, linetype = 2, linewidth = 0.3) +
    geom_jitter(height = 0.15, width = 0, size = 0.7, alpha = 0.6) +
    scale_x_log10(breaks = c(1e-4, 1e-3, 0.01, 0.05, 1), labels = c("0.0001", "0.001", "0.01", "0.05", "1")) +
    scale_colour_manual(values = COL_METRIC, guide = "none") + facet_wrap(~ metric) +
    labs(x = "p-value after removing one study (log scale)", title = "Test") + theme_asd
  save_fig(pS6, "FigS_leave_one_study_out", 174, 80)
}

cov <- bind_rows(read_csv(dir_out("tables", "coverage_depth_by_time_RR.csv"), show_col_types = FALSE) %>% mutate(metric = "Response ratio"),
                 read_csv(dir_out("tables", "coverage_depth_by_time_storage_rate.csv"), show_col_types = FALSE) %>% mutate(metric = "Storage rate")) %>%
  filter(!is.na(depth_group)) %>%
  mutate(depth_group = factor(depth_group, levels = levels(cut(0, DEPTH_BREAKS, include.lowest = TRUE, right = FALSE)),
                              labels = paste0(DEPTH_BREAKS[-length(DEPTH_BREAKS)], "\u2013", DEPTH_BREAKS[-1], " cm")),
         time_class = factor(time_class, levels = c("<=5 y", "5-10 y", "10-20 y", ">20 y"),
                             labels = c("\u22645 y", "5\u201310 y", "10\u201320 y", ">20 y")),
         lab = ifelse(n_obs == 0, "0", paste0(n_obs, "\n(", n_studies, ")")), dark = n_obs > 40)
pS7 <- ggplot(cov, aes(time_class, depth_group, fill = n_obs)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = lab, colour = dark), size = 2, lineheight = 0.85, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_y_discrete(limits = rev) + scale_fill_viridis_c(option = "mako", direction = -1, name = "Observations") +
  facet_wrap(~ metric) + labs(x = "Time since conversion", title = "Soil depth class") +
  theme_asd + theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "right")
save_fig(pS7, "FigS_data_coverage", 174, 80)

dsd <- ds %>% filter(!is.na(Grouped_Design), !Grouped_Design %in% c("Not specified", "Other"))
msd <- rma(yi = .y, vi = .v, scale = ~ Grouped_Design + time_since_conversion, data = dsd)
lev <- levels(factor(dsd$Grouped_Design)); tg <- seq(2, 20, 0.5); al <- msd$alpha[, 1]
sdg <- map_dfr(lev, function(L) {
  z <- al[1] + ifelse(L == lev[1], 0, al[paste0("Grouped_Design", L)]) + al["time_since_conversion"] * tg
  tibble(time = tg, tau = sqrt(exp(z)), design = c("Before-After Designs" = "Before\u2013after", "Control-Impact Designs" = "Control\u2013impact",
                                                  "Randomized Designs" = "Randomized")[[L]])
})
COL_D <- setNames(c(OI[["green"]], OI[["orange"]], OI[["purple"]]), sort(unique(sdg$design)))
pS8 <- ggplot(sdg, aes(time, tau, colour = design, linetype = design)) + geom_line(linewidth = 0.8) +
  scale_colour_manual(values = COL_D, name = NULL) + scale_linetype_manual(values = c(1, 2, 4), name = NULL) +
  labs(x = "Time since conversion (years)", title = expression("Heterogeneity SD of storage rates (Mg C ha"^-1*" yr"^-1*")")) +
  theme_asd + theme(legend.position = c(0.72, 0.85), legend.key.width = unit(7, "mm"))
save_fig(pS8, "FigS_precision_design_time", 84, 70)

# S3 – partial dependence by depth class (from step 05 outputs)
if (file.exists(dir_out("tables", "partial_dependence_by_depth.csv"))) {
  pdt <- read_csv(dir_out("tables", "partial_dependence_by_depth.csv"), show_col_types = FALSE)
  cls <- levels(cut(0, DEPTH_BREAKS, include.lowest = TRUE, right = FALSE))
  nice <- paste0(DEPTH_BREAKS[-length(DEPTH_BREAKS)], "\u2013", DEPTH_BREAKS[-1], " cm")
  pdt <- pdt %>% filter(group %in% cls) %>% mutate(group = factor(group, levels = cls, labels = nice))
  COL_DEP <- setNames(c("#440154", "#3B528B", "#21908C", "#5DC863", "#C9B400"), nice)
  pan <- function(metric_, mod, letter, ttl, xl) {
    s <- pdt %>% filter(metric == metric_, moderator == mod)
    ggplot(s, aes(x, fit_bt, colour = group, fill = group)) +
      geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
      geom_ribbon(aes(ymin = lwr_bt, ymax = upr_bt), alpha = 0.12, colour = NA) + geom_line(linewidth = 0.7) +
      scale_colour_manual(values = COL_DEP, name = "Soil depth", drop = FALSE) +
      scale_fill_manual(values = COL_DEP, name = "Soil depth", drop = FALSE) +
      labs(x = xl, title = ttl) + tag(letter) + theme_asd
  }
  pS3 <- (pan("RR", "time_since_conversion", "a", "Response ratio (%)", "Time since conversion (years)") |
          pan("storage_rate", "time_since_conversion", "b", UNIT_RATE, "Time since conversion (years)")) /
         (pan("RR", "control_soc_mean_T_ha", "c", "Response ratio (%)", UNIT_SOC) |
          pan("storage_rate", "control_soc_mean_T_ha", "d", UNIT_RATE, UNIT_SOC)) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  save_fig(pS3, "FigS_partial_dependence", 174, 150)
}
message("Figures written to ", dir_out("figures"))
