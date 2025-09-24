# Fit MetaForest model
mf_rr <- MetaForest(
  formula = yi ~ NEW_treatment_type2 + diff_species_class + History_reclass +
    control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
    MEAN_depth + time_since_conversion + Grouped_Design,
  data = FULL_sunRR,
  study = "id_expérimentation",
  whichweights = "random",
  num.trees = 8000
)

plot(mf_rr)


pred_values <- predict(mf_rr)

# Compute residuals
residuals <- mf_rr$data$yi - pred_values$predictions
boxplot(residuals, main="Residuals")

which(abs(residuals) > 3*sd(residuals))


# Indices des outliers
outliers <- c(107, 121, 155, 259, 262)

# Dataset sans les outliers
FULL_sunRR_clean <- FULL_sunRR[-outliers, ]
unique(FULL_sunRR_clean$NEW_treatment_type2)


# # Refitting du modèle MetaForest
# mf_rr_clean <- MetaForest(
#   formula = yi ~ NEW_treatment_type2 + diff_species_class + History_reclass +
#     control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
#     MEAN_depth + time_since_conversion + Grouped_Design,
#   data = FULL_sunRR_clean,
#   study = "id_expérimentation",
#   whichweights = "random",
#   num.trees = 8000
# )
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# 
# 
# # Variables numériques
# num_vars <- c("control_soc_mean_T_ha", "precipitation", "temperature", "MEAN_depth", "time_since_conversion")
# 
# # Préparer les données
# num_data <- FULL_sunRR_clean %>%
#   select(all_of(num_vars)) %>%
#   pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
# 
# # Plot
# ggplot(num_data, aes(x = value)) +
#   geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
#   facet_wrap(~variable, scales = "free_x") +
#   theme_minimal() +
#   labs(title = "Distribution des variables numériques",
#        x = "", y = "")
# 
# 
# # Variables catégorielles
# cat_vars <- c("NEW_treatment_type2", "diff_species_class", "History_reclass", "main_culture2", "Grouped_Design")
# 
# # Préparer les données
# cat_data <- FULL_sunRR_clean %>%
#   select(all_of(cat_vars)) %>%
#   pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
#   mutate(value = as.factor(value))  # assure que c'est un factor
# 
# # Plot
# ggplot(cat_data, aes(x = value)) +
#   geom_bar(fill = "tomato", color = "black", alpha = 0.7) +
#   facet_wrap(~variable, scales = "free_x") +
#   theme_minimal() +
#   labs(title = "Distribution des variables catégorielles",
#        x = "", y = "")+ coord_flip()
# 
# 
# 


##########################################
# MetaForest & XGBoost Analysis
# Author: Damien Beillouin
# Purpose: Hyperparameter tuning, in-sample prediction, and variable importance
# Dataset: FULL_sunRR_clean
##########################################


library(caret)
library(xgboost)
library(Matrix)

#-----------------------------
# 1. Tuning MetaForest
#-----------------------------
# Define tuning grid
tune_grid <- expand.grid(
  mtry = c(2, 4, 6, 8),
  min.node.size = c(5, 10, 15),
  num.trees = c(2000, 5000),
  sample.fraction = c(0.6, 0.9)
)
tune_grid$R2_mean <- NA
tune_grid$R2_sd <- NA

best_r2 <- 0
best_model <- NULL
n_repeats <- 4   # number of runs per combination for stability

for(i in 1:nrow(tune_grid)) {
  r2_vals <- numeric(n_repeats)
  
  for(k in 1:n_repeats){
    set.seed(100 + k)
    
    mf <- MetaForest(
      formula = yi ~ NEW_treatment_type2 + diff_species_class + History_reclass +
        control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
        MEAN_depth + time_since_conversion + Grouped_Design,
      data = FULL_sunRR_clean,
      study = "id_expérimentation",
      whichweights = "random",
      mtry = tune_grid$mtry[i],
      min.node.size = tune_grid$min.node.size[i],
      sample.fraction = tune_grid$sample.fraction[i],
      num.trees = tune_grid$num.trees[i]
    )
    
    # In-sample R²
    pred <- mf$predicted$predicted
    obs <- mf$data$yi
    r2_vals[k] <- mf$forest$r.squared
  }
  
  # Store mean and SD
  tune_grid$R2_mean[i] <- mean(r2_vals, na.rm = TRUE)
  tune_grid$R2_sd[i] <- sd(r2_vals, na.rm = TRUE)
  
  # Update best model
  if(tune_grid$R2_mean[i] > best_r2){
    best_r2 <- tune_grid$R2_mean[i]
    best_model <- mf
  }
}

# Inspect tuning results
print(best_r2)
print(best_model)
head(tune_grid)

#-----------------------------
# 2. Final MetaForest on full dataset
#-----------------------------
set.seed(252)
mf_final <- MetaForest(
  formula = yi ~ NEW_treatment_type2 + diff_species_class + History_reclass +
    control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
    MEAN_depth + time_since_conversion + Grouped_Design,
  data = FULL_sunRR_clean,
  study = "id_expérimentation",
  whichweights = "random",
  mtry = 6,
  min.node.size = 5,
  sample.fraction = 0.9,
  num.trees = 10000
)

# Predictions & metrics
pred_mf <- predict(mf_final)$predictions
obs <- mf_final$data$yi
R2_mf <- 1 - sum((obs - pred_mf)^2) / sum((obs - mean(obs))^2)
RMSE_mf <- sqrt(mean((obs - pred_mf)^2))

# Residual plot
boxplot(obs - pred_mf, main = "MetaForest Residuals")

# Variable importance
var_imp <- sort(mf_final$variable.importance, decreasing = TRUE)
head(var_imp, 5)

cat("MetaForest R²:", round(R2_mf, 3), "\n")
cat("MetaForest RMSE:", round(RMSE_mf, 3), "\n")

#-----------------------------
# 3. Gradient Boosting (XGBoost)
#-----------------------------
# Prepare matrices
X <- model.matrix(yi ~ NEW_treatment_type2 + diff_species_class + History_reclass +
                    control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
                    MEAN_depth + time_since_conversion + Grouped_Design, 
                  data = FULL_sunRR_clean)[,-1]
y <- FULL_sunRR_clean$yi
dtrain <- xgb.DMatrix(data = X, label = y)

# Parameters
params <- list(
  objective = "reg:squarederror",
  max_depth = 4,
  eta = 0.05,
  subsample = 0.9,
  colsample_bytree = 0.8
)

# Cross-validation
cv_results <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 500,
  nfold = 5,
  verbose = 1,
  early_stopping_rounds = 20,
  metrics = "rmse"
)
best_nrounds <- cv_results$best_iteration

# Train final model
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_nrounds
)

# Predictions & metrics
pred_xgb <- predict(xgb_model, newdata = X)
R2_xgb <- cor(y, pred_xgb)^2
RMSE_xgb <- sqrt(mean((y - pred_xgb)^2))

# Variable importance
importance_matrix <- xgb.importance(model = xgb_model)
print(head(importance_matrix, 10))

# Plot predictions
plot(y, pred_xgb, xlab = "Observed", ylab = "Predicted", main = "XGBoost Predictions")
abline(0,1,col="red")

cat("XGBoost R²:", round(R2_xgb,3), "\n")
cat("XGBoost RMSE:", round(RMSE_xgb,3), "\n")


#####



##########################################
# MetaForest Analysis for seq_rate
# Author: Damien Beillouin
# Purpose: Hyperparameter tuning, in-sample prediction, and variable importance
# Dataset: FULL_sun2
##########################################
FULL_sun2<- FULL_sun %>% filter(!is.na(control_soc_mean_T_ha))
FULL_sun2<- FULL_sun2 %>% filter(!is.na(time_since_conversion))

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


plot(mf_seq) 


pred_values <- predict(mf_seq)

# Compute residuals
residuals <- mf_seq$data$seq_rate - pred_values$predictions
boxplot(residuals, main="Residuals")

outliers<-which(abs(residuals) > 3*sd(residuals)) 


FULL_sun2_clean <- FULL_sun2[-outliers, ]
unique(FULL_sun2_clean$NEW_treatment_type2)
unique(FULL_sun$NEW_treatment_type2)

library(caret)

#-----------------------------
# 1. Tuning MetaForest
#-----------------------------
# Define tuning grid
tune_grid_seq <- expand.grid(
  mtry = c(2, 4, 6, 8),
  min.node.size = c(5, 10, 15),
  num.trees = c(2000, 5000),
  sample.fraction = c(0.6, 0.9)
)
tune_grid_seq$R2_mean <- NA
tune_grid_seq$R2_sd <- NA

best_r2 <- 0
best_model <- NULL
n_repeats <- 5   # number of runs per combination

for(i in 1:nrow(tune_grid_seq)) {
  r2_vals <- numeric(n_repeats)
  
  for(k in 1:n_repeats){
    set.seed(100 + k)
    
    mf <- MetaForest(
      formula = seq_rate ~ NEW_treatment_type2 + diff_species_class + History_reclass +
        control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
        MEAN_depth + time_since_conversion + Grouped_Design,
      data = FULL_sun2_clean,
      study = "id_article",
      whichweights = "random",
      mtry = tune_grid_seq$mtry[i],
      min.node.size = tune_grid_seq$min.node.size[i],
      sample.fraction = tune_grid_seq$sample.fraction[i],
      num.trees = tune_grid_seq$num.trees[i]
    )
    
    # In-sample R²
    pred <- mf$predicted$predicted
    obs <- mf$data$seq_rate
    r2_vals[k] <- mf$forest$r.squared
  }
  
  # Store mean and SD
  tune_grid_seq$R2_mean[i] <- mean(r2_vals, na.rm = TRUE)
  tune_grid_seq$R2_sd[i] <- sd(r2_vals, na.rm = TRUE)
  
  # Update best model
  if(tune_grid_seq$R2_mean[i] > best_r2){
    best_r2 <- tune_grid_seq$R2_mean[i]
    best_model <- mf
  }
}

# Inspect tuning results
print(best_r2)
print(best_model)
head(tune_grid_seq)

#-----------------------------
# 2. Final MetaForest on full dataset
#-----------------------------
set.seed(252)
mf_seq_final <- MetaForest(
  formula = seq_rate ~ NEW_treatment_type2 + diff_species_class + History_reclass +
    control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
    MEAN_depth + time_since_conversion + Grouped_Design,
  data = FULL_sun2_clean,
  study = "id_article",
  whichweights = "random",
  mtry = 6,
  min.node.size = 5,
  sample.fraction = 0.9,
  num.trees = 10000
)

# Predictions & metrics
pred_mf <- predict(mf_seq_final)$predictions
obs <- mf_seq_final$data$seq_rate
R2_seq <- 1 - sum((obs - pred_mf$predictions)^2) / sum((obs - mean(obs))^2)
RMSE_seq <- sqrt(mean((obs - pred_mf$predictions)^2))

# Residual plot
boxplot(obs - pred_mf, main = "MetaForest Residuals (seq_rate)")

# Variable importance
var_imp_seq <- sort(mf_seq_final$variable.importance, decreasing = TRUE)
head(var_imp_seq, 5)

cat("MetaForest seq_rate R²:", round(R2_seq, 3), "\n")
cat("MetaForest seq_rate RMSE:", round(RMSE_seq, 3), "\n")





#-----------------------------
# 3. Gradient Boosting (XGBoost)
#-----------------------------
# Prepare matrices
X <- model.matrix(seq_rate ~ NEW_treatment_type2 + diff_species_class + History_reclass +
                    control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
                    MEAN_depth + time_since_conversion + Grouped_Design, 
                  data = FULL_sun2_clean)[,-1]
y <- FULL_sun2_clean$seq_rate
dtrain <- xgb.DMatrix(data = X, label = y)

# Parameters
params <- list(
  objective = "reg:squarederror",
  max_depth = 4,
  eta = 0.05,
  subsample = 0.9,
  colsample_bytree = 0.8
)

# Cross-validation
cv_results <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 500,
  nfold = 5,
  verbose = 1,
  early_stopping_rounds = 20,
  metrics = "rmse"
)
best_nrounds <- cv_results$best_iteration

# Train final model
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_nrounds
)

# Predictions & metrics
pred_xgb <- predict(xgb_model, newdata = X)
R2_xgb <- cor(y, pred_xgb)^2
RMSE_xgb <- sqrt(mean((y - pred_xgb)^2))

# Variable importance
importance_matrix <- xgb.importance(model = xgb_model)
print(head(importance_matrix, 10))

# Plot predictions
plot(y, pred_xgb, xlab = "Observed", ylab = "Predicted", main = "XGBoost Predictions")
abline(0,1,col="red")

cat("XGBoost R²:", round(R2_xgb,3), "\n")
cat("XGBoost RMSE:", round(RMSE_xgb,3), "\n")



