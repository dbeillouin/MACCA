##########################################
# Forest plot (modeled PDP only, Nature Comm style)
# Facet by custom categories
# Includes observed counts (SEQ / RR) in labels
# Author: Damien Beillouin
##########################################

library(dplyr)
library(ggplot2)

# -----------------------------
# 0. Prepare observed counts
# -----------------------------

# Function to compute counts for a dataset by three categories
compute_counts <- function(df, prefix) {
  # 1. History_reclass
  hist <- df %>%
    group_by(History_reclass) %>%
    summarise(n = n(),
              n_studies = n_distinct(id_article),
              .groups = "drop") %>%
    mutate(category = "History_reclass", level = History_reclass) %>%
    dplyr::select(category, level, n, n_studies)
  
  # 2. NEW_treatment_type2
  treat <- df %>%
    group_by(NEW_treatment_type2) %>%
    summarise(n = n(),
              n_studies = n_distinct(id_article),
              .groups = "drop") %>%
    mutate(category = "NEW_treatment_type2", level = NEW_treatment_type2) %>%
    dplyr::select(category, level, n,n_studies)
  
  # 3. diff_species_class
  species <- df %>%
    group_by(diff_species_class) %>%
    summarise(n = n(),
              n_studies = n_distinct(id_article),
              .groups = "drop") %>%
    mutate(category = "diff_species_class", level = diff_species_class) %>%
    dplyr::select(category, level, n,n_studies)
  
  # Stack all rows
  counts <- bind_rows(hist, treat, species)
  names(counts)[3] <- paste0("n_", prefix)
  names(counts)[4] <- paste0("n_studies_", prefix)
  return(counts)
}

# Compute counts for SEQ and RR datasets
obs_seq_counts <- compute_counts(FULL_sun2_clean, "seq")
obs_RR_counts  <- compute_counts(FULL_sunRR, "RR")

# Merge counts
obs_counts <- full_join(obs_seq_counts, obs_RR_counts,
                        by = c("category", "level")) %>%
  mutate(
    n_seq = ifelse(is.na(n_seq), 0, n_seq),
    n_RR  = ifelse(is.na(n_RR), 0, n_RR),
    count_label = paste0(level, " (", n_seq, "/", n_RR, ")")
  )

# -----------------------------
# 1. Assign facet groups
# -----------------------------
results_table <- results_table %>%
  mutate(facet_group = case_when(
    level %in% c("Other", "Multistrata system", "Alley/Hedgerow", "Shaded perennial", "SilvoPasture") ~ "Type of agroforestry",
    level %in% c("C+2+", "A:<0/+1") ~ "Number of species",
    level %in% c("Grassland", "Cropland", "Forest", "Unknown/mixed") ~ "Previous land-use",
    TRUE ~ "Other"
  ))

# -----------------------------
# 2. Merge counts into results_table
# -----------------------------
results_table <- results_table %>%
  left_join(obs_counts %>% dplyr::select(level, count_label),
            by = "level") %>%
  group_by(facet_group) %>%
  arrange(desc(pred), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(level_f = factor(count_label, levels = unique(count_label))) %>% 
  filter(!level=="Other")


# -----------------------------
# 3. Create forest plot
# -----------------------------

library(patchwork)
library(ggplot2)


# ---- Plot SEQ ----

# Colonnes nécessaires
cols <- names(results_table)

# Créer ligne vide pour Multistrata system
placeholder <- tibble::tibble(
  level      = "Multistrata system",
  pred       = NA_real_,
  se         = NA_real_,
  y_mean     = NA_real_,
  y_se       = NA_real_,
  y_lci      = NA_real_,
  y_uci      = NA_real_,
  var        = "seq",
  facet_group= "Type of agroforestry",  # à ajuster selon ton facet
  count_label= "Multistrata system (0/0)",
  level_f    = "Multistrata system (0/32)"
)

# Ajouter au tableau
results_table <- bind_rows(results_table, placeholder) %>%
  arrange(facet_group, var, desc(pred)) %>%
  mutate(level_f = factor(level_f, levels = unique(level_f)))

results_table<- results_table %>% filter(!level=="Unknown/mixed")


forest_seq <- ggplot(results_table %>% filter(var=="seq")) +
  geom_point(aes(x = pred, y = level_f), color = "steelblue", size = 4) +
  geom_errorbarh(aes(xmin = lci,
                     xmax = uci,
                     y = level_f),
                 color = "steelblue", height = 0.2) +
  facet_wrap(~facet_group, scales="free_y", ncol=1) +
  labs(x = "Partial Dependence (mean ± 95% CI)",
       y = "",
       title = "SEQ") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(face="bold", size=14),
    axis.text.y = element_text(face="italic", size=12),
    plot.title = element_text(face="bold", hjust=0.5)
  )

# ---- Plot RR ----
forest_RR <- ggplot(results_table %>% filter(var=="rr")) +
  geom_point(aes(x = exp(pred), y = level_f), color = "steelblue", size = 4) +
  geom_errorbarh(aes(xmin = exp(lci),
                     xmax = exp(uci),
                     y = level_f),
                 color = "steelblue", height = 0.2) +
  facet_wrap(~facet_group, scales="free_y", ncol=1) +
  labs(x = "Partial Dependence (mean ± 95% CI)",
       y = "",
       title = "RR") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(face="bold", size=14),
    plot.title = element_text(face="bold", hjust=0.5)
  )

# ---- Combine side by side with equal widths ----
combined_plot <- forest_seq | forest_RR + 
  plot_layout(widths = c(69999,1)) +
  theme(plot.margin = margin(5,5,5,5))

combined_plot



library(ggplot2)
library(dplyr)
library(ggplot2)
library(dplyr)
library(ggplot2)

dodge_width <- 0.6

# Échelle visuelle pour RR
rr_min <- 0.10
rr_max <- 0.30

rr_vals <- results_table %>% filter(var == "rr") %>% pull(pred)
rr_lci  <- results_table %>% filter(var == "rr") %>% pull(lci)
rr_uci  <- results_table %>% filter(var == "rr") %>% pull(uci)

results_table <- results_table %>%
  mutate(
    pred_plot = case_when(
      var == "rr"  ~ rr_min + (pred - min(rr_vals)) / (max(rr_vals) - min(rr_vals)) * (rr_max - rr_min),
      TRUE ~ pred
    ),
    lci_plot = case_when(
      var == "rr"  ~ rr_min + (lci - min(rr_vals)) / (max(rr_vals) - min(rr_vals)) * (rr_max - rr_min),
      TRUE ~ lci
    ),
    uci_plot = case_when(
      var == "rr"  ~ rr_min + (uci - min(rr_vals)) / (max(rr_vals) - min(rr_vals)) * (rr_max - rr_min),
      TRUE ~ uci
    )
  )


ggplot(results_table, aes(y = level_f, x = pred_plot, color = var, group = var)) +
  geom_point(size = 4, position = position_dodge(width = dodge_width)) +
  geom_errorbarh(aes(xmin = lci_plot, xmax = uci_plot),
                 height = 0.2,
                 position = position_dodge(width = dodge_width)) +
  scale_color_manual(values = c("steelblue", "darkorange")) +
  
  # Axes
  scale_x_continuous(
    name = "SEQ",
    sec.axis = sec_axis(~ ., name = "RR (log scale)")
  ) +
  
  facet_wrap(~facet_group, scales = "free_y", ncol = 1) +
  
  labs(y = "", title = "Partial Dependence") +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


library(ggh4x)

ggplot(results_table, aes(y = level_f, x = pred_plot, color = var, group = var)) +
  geom_point(size = 4, position = position_dodge(width = dodge_width)) +
  geom_errorbarh(aes(xmin = lci_plot, xmax = uci_plot),
                 height = 0.2,
                 position = position_dodge(width = dodge_width)) +
  scale_color_manual(values = c("steelblue", "darkorange")) +
  
  # Axes
  scale_x_continuous(
    name = "SEQ",
    sec.axis = sec_axis(~ ., name = "RR (log scale)")
  ) +
  
  # Grid à une colonne -> identique à wrap, mais avec espace proportionnel
  ggh4x::facet_grid2(
    rows = vars(facet_group),
    scales = "free_y",
    space = "free_y"
  ) +
  
  labs(y = "", title = "Partial Dependence") +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


