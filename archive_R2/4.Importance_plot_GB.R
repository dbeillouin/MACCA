##########################################
# XGBoost Recursive Bootstrap Importance
# Author: Damien Beillouin
# Purpose: Compute robust variable importance for SEQ and RR datasets
##########################################

library(dplyr)
library(ggplot2)
library(stringr)
library(xgboost)
library(tidyr)

#-----------------------------------
# 1. Parameters
#-----------------------------------
n_repeats <- 100   # number of bootstrap runs
threshold <- 0.01  # minimal importance to keep variable

#-----------------------------------
# 2. Function to run bootstrap XGBoost and extract importance
#-----------------------------------
compute_xgb_importance <- function(data, response, predictors, n_repeats = 100, threshold = 0.01) {
  
  importance_list <- vector("list", n_repeats)
  
  for(k in 1:n_repeats){
    # Bootstrap sample
    idx <- sample(1:nrow(data), replace = TRUE)
    data_boot <- data[idx, ]
    
    # Build design matrix
    X <- model.matrix(as.formula(paste(response, "~", paste(predictors, collapse = "+"))), data = data_boot)[,-1]
    y <- data_boot[[response]]
    dtrain <- xgb.DMatrix(data = X, label = y)
    
    # Train XGBoost
    model <- xgb.train(
      params = list(objective = "reg:squarederror",
                    max_depth = 4,
                    eta = 0.05,
                    subsample = 0.9,
                    colsample_bytree = 0.8),
      data = dtrain,
      nrounds = 500,
      verbose = 0
    )
    
    # Extract importance
    imp <- xgb.importance(model = model) %>%
      dplyr::select(Feature, Gain)
    imp$run <- k
    importance_list[[k]] <- imp
  }
  
  # Combine and group categorical variables
  importance_combined <- bind_rows(importance_list) %>%
    mutate(variable_origin = Feature,
           variable_origin = str_replace(variable_origin, "NEW_treatment_type2.*", "NEW_treatment_type2"),
           variable_origin = str_replace(variable_origin, "Grouped_Design.*", "Grouped_Design"),
           variable_origin = str_replace(variable_origin, "History_reclass.*", "History_reclass"),
           variable_origin = str_replace(variable_origin, "main_culture2.*", "main_culture2"),
           variable_origin = str_replace(variable_origin, "diff_species_class.*", "diff_species_class")) 
  
  # Compute summary statistics per variable
  importance_stats <- importance_combined %>%
    group_by(variable_origin) %>%
    summarize(
      mean_value = mean(Gain, na.rm = TRUE),
      lower = quantile(Gain, 0.17, na.rm = TRUE),
      upper = quantile(Gain, 0.83, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(mean_value > threshold)
  
  # Normalize for Low/Medium/High
  importance_stats <- importance_stats %>%
    mutate(mean_value_scaled = (mean_value - min(mean_value)) / (max(mean_value) - min(mean_value)))
  
  return(importance_stats)
}

#-----------------------------------
# 3. Define datasets and predictors
#-----------------------------------
predictors <- c("NEW_treatment_type2", "diff_species_class", "History_reclass",
                "control_soc_mean_T_ha", "precipitation", "temperature",
                "main_culture2", "MEAN_depth", "time_since_conversion",
                "Grouped_Design")

FULL_sun2_clean
# SEQ dataset
importance_seq <- compute_xgb_importance(FULL_sun2_clean, "seq_rate", predictors, n_repeats, threshold) %>%
  mutate(varset = "SEQ")

# RR dataset
FULL_sunRR_cleanDD<- FULL_sunRR_clean %>% filter(!id_expérimentation %in% c(934,939))
importance_rr <- compute_xgb_importance(FULL_sunRR_cleanDD, "yi", predictors, n_repeats, threshold) %>%
  mutate(varset = "RR")

#-----------------------------------
# 4. Clean variable names
#-----------------------------------
clean_variable_names <- function(var) {
  case_when(
    var == "Grouped_Design"        ~ "Experimental design",
    var == "History_reclass"       ~ "Land-use history",
    var == "MEAN_depth"            ~ "Soil depth",
    var == "NEW_treatment_type2"   ~ "Agroforestry type",
    var == "control_soc_mean_T_ha" ~ "Initial SOC stock",
    var == "diff_species_class"    ~ "Number of species added",
    var == "main_culture2"         ~ "Main crop type",
    var == "precipitation"         ~ "Precipitation",
    var == "temperature"           ~ "Temperature",
    var == "time_since_conversion" ~ "Time since conversion",
    TRUE                           ~ str_replace_all(var, "_", " ")
  )
}

importance_seq <- importance_seq %>%
  mutate(variable = clean_variable_names(variable_origin))
importance_rr <- importance_rr %>%
  mutate(variable = clean_variable_names(variable_origin))

#-----------------------------------
# 5. Combine datasets
#-----------------------------------
combined_data <- bind_rows(importance_seq, importance_rr)

#-----------------------------------
# 6. Compute overall summary for plotting
#-----------------------------------
importance_stats <- combined_data %>%
  group_by(variable, varset) %>%
  summarize(
    mean_value = mean(mean_value, na.rm = TRUE),
    lower = mean(lower, na.rm = TRUE),
    upper = mean(upper, na.rm = TRUE),
    .groups = "drop"
  )

overall_summary <- importance_stats %>%
  group_by(variable) %>%
  summarize(mean_value = mean(mean_value), .groups = "drop") %>%
  mutate(
    mean_value_scaled = (mean_value - min(mean_value)) / (max(mean_value) - min(mean_value)),
    sel_class = case_when(
      mean_value_scaled < 0.33 ~ "Low",
      mean_value_scaled < 0.66 ~ "Medium",
      TRUE ~ "High"
    ),
    sel_class = factor(sel_class, levels = c("Low", "Medium", "High"))
  )

#-----------------------------------
# 7. Publication-ready barplot
#-----------------------------------
barplot_final <- ggplot(overall_summary, aes(
  x = reorder(variable, mean_value),
  y = mean_value,
  fill = sel_class
)) +
  geom_col(width = 0.6, color = "black", size = 0.3) +
  geom_errorbar(
    data = importance_stats,
    aes(x = reorder(variable, mean_value), ymin = lower, ymax = upper, color = varset),
    width = 0.2,
    position = position_dodge(width = 0.7),
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("Low" = "#d3e5ff", "Medium" = "#6baed6", "High" = "#08306b")) +
  scale_color_manual(values = c("SEQ" = "#35978f", "RR" = "#bf812d")) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  theme(
  #  axis.title = element_blank(),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(
    fill = "Importance Class",
    color = "Dataset",
    y= "Relative mean importance",
    x="",
    title = "",
    subtitle = "",
    caption = ""
  )

print(barplot_final)

overall_summary$relative <- overall_summary$mean_value / sum(overall_summary$mean_value)

