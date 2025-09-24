##########################################
# MetaForest Analysis and Visualization
# Author: Damien Beillouin
# Purpose: High-quality, reproducible workflow for publication
# Description: Preselection of important predictors and visualization of weighted means
#              for two datasets (SEQ and RR) using MetaForest and ggdist
##########################################

# Load required packages
library(caret)        # for model training and preprocessing
library(ggdist)       # for half-eye plots with uncertainty intervals
library(metaforest)   # for MetaForest meta-analytic random forest
library(dplyr)        # data manipulation
library(tidyr)        # data reshaping
library(ggplot2)      # plotting
library(ggpubr)       # for professional themes
library(patchwork)    # combine multiple plots

##############################
# Step 1: Prepare SEQ dataset
##############################

# Convert and recode variables
FULL_sun <- FULL_sun %>%
  dplyr::mutate(
    dplyr::across(
      c(Grouped_Design, NEW_treatment_type2, History_reclass, main_culture2),
      as.factor
    ),
    Crate = as.numeric(gsub(",", ".", Crate)),
    vi    = as.numeric(gsub(",", ".", vi)),
    vi    = vi * (1 + (4 - quality_score) / 10)
  )

# Filter out rows with missing target variable
FULL_sun2 <- FULL_sun %>% filter(!is.na(control_soc_mean_T_ha))

# Fit MetaForest model
mf_seq <- MetaForest(
  formula = seq_rate ~ NEW_treatment_type2 + diff_species_class + History_reclass +
    control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
    MEAN_depth + time_since_conversion + Grouped_Design,
  data = FULL_sun2,
  study = "id_article",
  whichweights = "random",
  num.trees = 8000
)

# Preselection of important predictors (recursive selection, stratified)
preselected_seq <- preselect(
  mf_seq,
  replications = 200,
  algorithm = "recursive"
)

# Normalize preselection frequencies
selected_seq <- preselected_seq$selected
selected_seq_norm <- t(apply(selected_seq, 1, function(x) x / max(x, na.rm = TRUE)))
selected_seq_long <- pivot_longer(
  as.data.frame(selected_seq_norm),
  cols = everything(),
  names_to = "variable",
  values_to = "value"
)
selected_seq_long$num <- rep(1:200, each = ncol(selected_seq_norm))
selected_seq_long$var <- "SEQ"

##############################
# Step 2: Prepare RR dataset
##############################

# Convert and recode variables
FULL_sunRR <- FULL_sunRR %>%
  dplyr::mutate(
    Grouped_Design     = as.factor(Grouped_Design),
    Ratio              = as.numeric(gsub(",", ".", Ratio)),
    NEW_treatment_type2 = dplyr::recode(
      as.character(NEW_treatment_type),
      "Alley cropping" = "Alley/Hedgerow",
      "Hedgerow"       = "Alley/Hedgerow",
      .default         = as.character(NEW_treatment_type)
    ),
    NEW_treatment_type2 = as.factor(NEW_treatment_type2),
    History_reclass    = as.factor(History_reclass),
    main_culture2      = as.factor(main_culture2),
    Crate_vi           = as.numeric(gsub(",", ".", Crate_vi)),
    vi                 = vi * (1 + (4 - quality_score) / 10)
  )



# Reclassify diff_species into categorical variable
FULL_sunRR$diff_species_class <- cut(
  FULL_sunRR$diff_species,
  breaks = c(-5, 1.1, 25),
  labels = c("A:<0/+1", "C+2+"),
  right = FALSE
)
FULL_sunRR$diff_species_class[is.na(FULL_sunRR$diff_species_class)] <- "C+2+"

# Impute missing values with median
FULL_sunRR$time_since_conversion[is.na(FULL_sunRR$time_since_conversion)] <- median(FULL_sunRR$time_since_conversion, na.rm = TRUE)

# Fit MetaForest model
mf_rr <- MetaForest(
  formula = yi ~ NEW_treatment_type2 + diff_species_class + History_reclass +
    control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
    MEAN_depth + time_since_conversion + Grouped_Design,
  data = FULL_sunRR,
  study = "id_article",
  whichweights = "random",
  num.trees = 8000
)

# Preselection of important predictors
preselected_rr <- preselect(
  mf_rr,
  replications = 100,
  algorithm = "recursive"
)

# Normalize preselection frequencies
selected_rr <- preselected_rr$selected
selected_rr_norm <- t(apply(selected_rr, 1, function(x) x / max(x, na.rm = TRUE)))
selected_rr_long <- pivot_longer(
  as.data.frame(selected_rr_norm),
  cols = everything(),
  names_to = "variable",
  values_to = "value"
)
selected_rr_long$num <- rep(1:100, each = ncol(selected_rr_norm))
selected_rr_long$var <- "RR"

##############################
# Step 3: Combine datasets for visualization
##############################

combined_data <- rbind(selected_seq_long, selected_rr_long)

# Common plotting theme for half-eye plots
halfeye_theme <- theme_pubr(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 14),
    plot.caption = element_text(size = 10, face = "italic"),
    panel.background = element_rect(fill = "white", color = "black"),
    legend.position = "bottom"
  )

# Main combined half-eye plot
global_plot <- ggplot(combined_data, aes(x = value, y = reorder(variable, value, FUN = mean, na.rm = TRUE), fill = var)) +
  stat_halfeye(
    position = position_dodge(width = 0.5),
    scale = 0.35, normalize = "xy",
    point_interval = median_qi,
    .width = c(0.50, 0.66),
    interval_colour = "black",
    point_colour = "black",
    fatten_point = 0.1
  ) +
  scale_fill_manual(values = c("#8A9B7D", "#ff7f0e")) +
  scale_color_manual(values = c("#8A9B7D", "#ff7f0e")) +
  halfeye_theme +
  geom_vline(aes(xintercept = 0.5), linetype = "dashed", color = "gray80", size = 0.7) +
  stat_summary(
    fun = median, geom = "point", shape = 21, size = 3,
    fill = "white", color = "black",
    position = position_dodge(width = 0.5)
  )

##############################
# Step 4: Compute percentage of non-NA values per variable
##############################

hist_data <- combined_data %>%
  group_by(variable) %>%
  summarize(non_na_count = sum(!is.na(value)), .groups = "drop") %>%
  mutate(perc = non_na_count / nrow(combined_data))

# Optional ordering for aesthetics
variable_order <- data.frame(
  variable = c("temperature", "control_soc_mean_T_ha", "time_since_conversion",
               "precipitation", "MEAN_depth", "main_culture2", "History_reclass", "Grouped_Design",
               "NEW_treatment_type2", "diff_species_class"),
  order = 1:10
)
hist_data <- left_join(hist_data, variable_order, by = "variable")

# Histogram of data coverage
hist_plot <- ggplot(hist_data, aes(x = reorder(variable, -order), y = perc)) +
  geom_bar(stat = "identity", fill = "gray80", color = "black", width = 0.7, size = 0.3) +
  labs(x = "", y = "") +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 16),
    plot.caption = element_text(size = 10, face = "italic"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "black")
  ) +
  coord_flip()

##############################
# Step 5: Combine plots for publication-quality figure
##############################

final_plot <- global_plot + hist_plot +
  plot_layout(ncol = 2, widths = c(2, 1), guides = "collect") +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    axis.title = element_blank(),
    axis.text = element_text(size = 12)
  )

# Display final figure
print(final_plot)



###
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# ---------------------------
# 1. Combine SEQ and RR preselection data
# ---------------------------

combined_data <- bind_rows(
  selected_seq_long %>% mutate(varset = "SEQ"),
  selected_rr_long  %>% mutate(varset = "RR")
)

# ---------------------------
# 2. Compute importance statistics
# ---------------------------

importance_summary <- combined_data %>%
  group_by(variable) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.17, na.rm = TRUE),   # 66% CI lower
    upper = quantile(value, 0.83, na.rm = TRUE),   # 66% CI upper
    sel_freq = mean(value, na.rm = TRUE),          # normalized selection frequency
    .groups = "drop"
  ) %>%
  arrange(desc(mean_value))

# ---------------------------
# 3. Discrete selection frequency classes
# ---------------------------

summary_data <- importance_summary %>%
  mutate(
    sel_class = case_when(
      sel_freq < 0.33 ~ "Low",
      sel_freq < 0.66 ~ "Medium",
      TRUE ~ "High"
    ),
    sel_class = factor(sel_class, levels = c("Low", "Medium", "High"))
  )

# ---------------------------
# 4. Publication-ready barplot
# ---------------------------

barplot_final <- ggplot(summary_data, aes(
  x = reorder(variable, mean_value),
  y = mean_value,
  fill = sel_class
)) +
  geom_col(width = 0.6, color = "black", size = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "black") +
  scale_fill_manual(values = c("Low" = "#d3e5ff", "Medium" = "#6baed6", "High" = "#08306b")) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(
    fill = "Selection Frequency",
    title = "Variable Importance and Robustness in MetaForest Models",
    subtitle = "Bar height = mean importance; error bars = 66% CI; color = selection class",
    caption = "Source: MetaForest preselection analysis"
  )

print(barplot_final)


#####
#2#
####


importance_stats <- combined_data %>%
  group_by(variable, varset) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.17, na.rm = TRUE),
    upper = quantile(value, 0.83, na.rm = TRUE),
    .groups = "drop"
  )

overall_summary <- importance_stats %>%
  group_by(variable) %>%
  summarize(
    mean_value = mean(mean_value),
    .groups = "drop"
  )

# Compute overall mean importance across datasets
overall_summary <- overall_summary %>%
  mutate(
    sel_class = case_when(
      mean_value < 0.33 ~ "Low",
      mean_value < 0.66 ~ "Medium",
      TRUE ~ "High"
    ),
    sel_class = factor(sel_class, levels = c("Low", "Medium", "High"))
  )



clean_variable_names <- function(var) {
  var <- case_when(
    var == "Grouped_Design"           ~ "Experimental design",
    var == "History_reclass"          ~ "Land-use history",
    var == "MEAN_depth"               ~ "Soil depth",
    var == "NEW_treatment_type2"      ~ "Agroforestry type",
    var == "control_soc_mean_T_ha"    ~ "Initial SOC stock",
    var == "diff_species_class"       ~ "Number of species added",
    var == "main_culture2"            ~ "Main crop type",
    var == "precipitation"            ~ "Precipitation",
    var == "temperature"              ~ "Temperature",
    var == "time_since_conversion"    ~ "Time since conversion",
    TRUE                              ~ str_replace_all(var, "_", " ")
  )
  (var)
}

# Apply to both datasets
overall_summary <- overall_summary %>%
  mutate(variable = clean_variable_names(variable))

importance_stats <- importance_stats %>%
  mutate(variable = clean_variable_names(variable))



barplot_simple <- ggplot(overall_summary, aes(
  x = reorder(variable, mean_value),
  y = mean_value,
  fill = sel_class
)) +
  geom_col(width = 0.7, color = "black", size = 0.5) +  # thicker bars & border
  geom_errorbar(
    data = importance_stats,
    aes(x = reorder(variable, mean_value), ymin = lower, ymax = upper, color = varset),
    width = 0.25,
    position = position_dodge(width = 0.7),
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Low"    = "#c6dbef",  # pastel bleu clair
      "Medium" = "#6baed6",  # bleu moyen
      "High"   = "#2171b5"   # bleu foncé
    )
  ) +
  scale_color_manual(
    values = c(
      "SEQ" = "#35978f",  # violet pastel
      "RR"  = "#bf812d"   # violet moyen
    )
  ) +
  coord_flip()+
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(
    fill = "Selection Frequency",
    color = "Dataset",
  )

print(barplot_simple)
