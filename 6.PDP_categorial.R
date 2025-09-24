  ##############################################
  # Partial Dependence for Categorical Variables
  # Clean, fully data-driven
  # Author: Damien Beillouin
  ##############################################
  
  library(xgboost)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  set.seed(594)
  
  # ---------------------------------------------
  # 0. Ensure factors are consistent (keep original levels)
  # ---------------------------------------------
  factor_vars <- c("NEW_treatment_type2", "History_reclass", "diff_species_class")
  FULL_sun2_sub_seq[factor_vars] <- lapply(FULL_sun2_sub_seq[factor_vars], factor)
  FULL_sun2_sub_rr[factor_vars]  <- lapply(FULL_sun2_sub_rr[factor_vars], factor)
  
  VERIF<- FULL_sun2_sub_seq %>% filter(NEW_treatment_type2=="Other")
  # ---------------------------------------------
  # 1. Prepare X matrix for XGBoost (contr.treatment)
  # ---------------------------------------------
  X_seq <- model.matrix(
    ~ NEW_treatment_type2 + diff_species_class + History_reclass +
      control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
      MEAN_depth + time_since_conversion + Grouped_Design,
    data = FULL_sun2_sub_seq,
    contrasts.arg = list(
      NEW_treatment_type2 = "contr.treatment",
      diff_species_class = "contr.treatment",
      History_reclass = "contr.treatment"
    )
  )[, -1]  # remove intercept
  
  y_seq <- FULL_sun2_sub_seq$seq_rate
  dtrain_seq <- xgb.DMatrix(data = X_seq, label = y_seq)
  
  # ---------------------------------------------
  # 2. Train XGBoost
  # ---------------------------------------------
  cv_seq <- xgb.cv(
    params = list(objective = "reg:squarederror", max_depth = 4, eta = 0.05,
                  subsample = 0.9, colsample_bytree = 0.8),
    data = dtrain_seq,
    nrounds = 500,
    nfold = 5,
    early_stopping_rounds = 20,
    metrics = "rmse",
    verbose = 0
  )
  
  best_nrounds_seq <- cv_seq$best_iteration
  
  xgb_seq <- xgboost(
    data = X_seq,
    label = y_seq,
    nrounds = best_nrounds_seq,
    objective = "reg:squarederror",
    max_depth = 4,
    eta = 0.05,
    subsample = 0.9,
    colsample_bytree = 0.8,
    verbose = 0
  )
  
  # ---------------------------------------------
  # 3. Generic bootstrap PDP function
  # ---------------------------------------------
  bootstrap_pdp_cat <- function(cat_var, xgb_model, X_full, n_boot = 500,
                                continuous_refs = list(),
                                categorical_refs = list(),
                                depth_vals = NULL,
                                cond_sample_vars = NULL) {
    
    # Colonnes dummy pour la variable catégorielle
    dummy_cols <- grep(cat_var, colnames(X_full), value = TRUE)
    
    # Niveaux de la variable
    all_levels <- gsub(paste0(cat_var), "", dummy_cols)
    ref_level  <- setdiff(levels(FULL_sun2_sub_seq[[cat_var]]), all_levels)
    pdp_levels <- c(ref_level, all_levels)
    
    pdp_boot <- lapply(pdp_levels, function(level_name){
      
      pdp_vals <- numeric(n_boot)
      
      for(b in 1:n_boot){
        # Bootstrap rows
        boot_rows <- sample(seq_len(nrow(X_full)), replace = TRUE)
        Xb <- X_full[boot_rows, , drop = FALSE]
        
        # Set target categorical variable
        dummy_col <- grep(paste0(cat_var, level_name), colnames(Xb))
        if(length(dummy_col) == 0){
          Xb[, grep(cat_var, colnames(Xb))] <- 0
        } else {
          Xb[, setdiff(grep(cat_var, colnames(Xb)), dummy_col)] <- 0
          Xb[, dummy_col] <- 1
        }
        
        # Fix ou échantillonne les covariables continues
        for(cv in colnames(Xb)[sapply(Xb, is.numeric)]){
          if(is.na(cv)) next
          
          if(cv %in% names(continuous_refs)){
            # Si fixée dans continuous_refs
            Xb[, cv] <- continuous_refs[[cv]]
            
          } else if(!is.null(depth_vals) && cv == "MEAN_depth"){
            Xb[, cv] <- sample(depth_vals, nrow(Xb), replace = TRUE)
            
          } else if(!is.null(cond_sample_vars) && cv %in% cond_sample_vars){
            # Échantillonner conditionnellement au niveau de la variable cible
            idx_cond <- which(FULL_sun2_sub_seq[[cat_var]] == level_name)
            Xb[, cv] <- sample(FULL_sun2_sub_seq[idx_cond, cv], nrow(Xb), replace = TRUE)
            
          } else {
            # Sinon échantillonner dans la distribution empirique globale
            Xb[, cv] <- sample(Xb[, cv], nrow(Xb), replace = TRUE)
          }
        }
        
        # Fixer les variables catégorielles autres que celle cible
        for(cv in names(categorical_refs)){
          dummy_cv_cols <- grep(cv, colnames(Xb), value = TRUE)
          dummy_cv_col  <- paste0(cv, "_", categorical_refs[[cv]])
          if(dummy_cv_col %in% colnames(Xb)){
            Xb[, dummy_cv_col] <- 1
            Xb[, setdiff(dummy_cv_cols, dummy_cv_col)] <- 0
          }
        }
        
        # Prediction moyenne pour le bootstrap
        pdp_vals[b] <- mean(predict(xgb_model, newdata = Xb))
      }
      
      data.frame(
        level = level_name,
        pred  = median(pdp_vals),
        lci   = quantile(pdp_vals, 0.025),
        uci   = quantile(pdp_vals, 0.975)
      )
    })
    
    dplyr::bind_rows(pdp_boot)
  }
  
  
  
  
  
  # ---------------------------------------------
  # 4. Continuous and categorical references
  # ---------------------------------------------
  continuous_refs <- list(
   # MEAN_depth = 20,
    control_soc_mean_T_ha = median(FULL_sun2_sub_seq$control_soc_mean_T_ha, na.rm = TRUE),
    time_since_conversion = 5
  #  precipitation = median(FULL_sun2_sub_seq$precipitation, na.rm = TRUE),
  #  temperature = median(FULL_sun2_sub_seq$temperature, na.rm = TRUE)
  )
  
  categorical_refs <- list(
    Grouped_Design = names(sort(table(FULL_sun2_sub_seq$Grouped_Design), decreasing = TRUE))[1],
    diff_species_class = names(sort(table(FULL_sun2_sub_seq$diff_species_class), decreasing = TRUE))[1],
    History_reclass = names(sort(table(FULL_sun2_sub_seq$History_reclass), decreasing = TRUE))[2]
  )
  
  ## For species, we adapt the initial conditions
  continuous_refs_species <- list(
   # MEAN_depth = depth_ref,
    control_soc_mean_T_ha = soc_ref,
  #  precipitation = precip_ref,
   # temperature = temp_ref,
    time_since_conversion = 5
  )
  
  categorical_refs_species <- list(
    Grouped_Design = names(sort(table(FULL_sun2_sub_seq$Grouped_Design), decreasing = TRUE))[1],
    History_reclass = names(sort(table(FULL_sun2_sub$History_reclass), decreasing = TRUE))[2]
   # main_culture2 = names(sort(table(FULL_sun2_sub$main_culture2), decreasing = TRUE))[2]
  )
  
  # ---------------------------------------------
  # 5. Compute PDPs
  # ---------------------------------------------
  depth_vals_seq <- seq(2, 20, by = 2)
  pdp_treatment <- bootstrap_pdp_cat("NEW_treatment_type2", xgb_seq, X_seq,
                                     n_boot = 500,
                                     continuous_refs = continuous_refs,
                                     categorical_refs = categorical_refs,
                                     cond_sample_vars = c("control_soc_mean_T_ha"),
                                     depth_vals = depth_vals_seq)
  
  pdp_history   <- bootstrap_pdp_cat("History_reclass", xgb_seq, X_seq,
                                     n_boot = 500,
                                     continuous_refs = continuous_refs,
                                     categorical_refs = categorical_refs,
                                     cond_sample_vars = c("control_soc_mean_T_ha"),
                                     depth_vals = depth_vals_seq)
  
  pdp_species   <- bootstrap_pdp_cat("diff_species_class", xgb_seq, X_seq,
                                     n_boot = 500,
                                     continuous_refs = continuous_refs_species,
                                     categorical_refs = categorical_refs_species,
                                     cond_sample_vars = c("control_soc_mean_T_ha"),
                                     depth_vals = depth_vals_seq)
  
  # Combine all
  pdp_all <- bind_rows(
    pdp_treatment  %>% mutate(variable = "NEW_treatment_type2"),
    pdp_history    %>% mutate(variable = "History_reclass"),
    pdp_species    %>% mutate(variable = "diff_species_class")
  )
  
  
         
  
  # ---------------------------------------------
  # 6. Prepare X matrix for RR
  # ---------------------------------------------
  X_rr <- model.matrix(
    ~ NEW_treatment_type2 + diff_species_class + History_reclass +
      control_soc_mean_T_ha + precipitation + temperature + main_culture2 +
      MEAN_depth + time_since_conversion + Grouped_Design,
    data = FULL_sun2_sub_rr,
    contrasts.arg = list(
      NEW_treatment_type2 = "contr.treatment",
      diff_species_class = "contr.treatment",
      History_reclass = "contr.treatment"
    )
  )[, -1]
  
  y_rr <- FULL_sun2_sub_rr$yi
  dtrain_rr <- xgb.DMatrix(data = X_rr, label = y_rr)
  
  # ---------------------------------------------
  # 7. Train XGBoost for RR
  # ---------------------------------------------
  cv_rr <- xgb.cv(
    params = list(objective = "reg:squarederror", max_depth = 4, eta = 0.05,
                  subsample = 0.9, colsample_bytree = 0.8),
    data = dtrain_rr,
    nrounds = 500,
    nfold = 5,
    early_stopping_rounds = 20,
    metrics = "rmse",
    verbose = 0
  )
  
  best_nrounds_rr <- cv_rr$best_iteration
  
  xgb_rr <- xgboost(
    data = X_rr,
    label = y_rr,
    nrounds = best_nrounds_rr,
    objective = "reg:squarederror",
    max_depth = 4,
    eta = 0.05,
    subsample = 0.9,
    colsample_bytree = 0.8,
    verbose = 0
  )
  
  # ---------------------------------------------
  # 8. Continuous and categorical refs for RR
  # ---------------------------------------------
  categorical_refs_rr <- categorical_refs
  
  # ---------------------------------------------
  # 9. Compute PDPs for RR
  # ---------------------------------------------
  continuous_refs_rr <- list(
  #  MEAN_depth = 20,
    control_soc_mean_T_ha = median(FULL_sun2_sub_seq$control_soc_mean_T_ha, na.rm = TRUE),
    time_since_conversion = 25,
    precipitation = median(FULL_sun2_sub_seq$precipitation, na.rm = TRUE)
  #  temperature = median(FULL_sun2_sub_seq$temperature, na.rm = TRUE)
  )
  
  categorical_refs_species_rr <- list(
    # MEAN_depth = depth_ref,
    control_soc_mean_T_ha = median(FULL_sun2_sub_seq$control_soc_mean_T_ha, na.rm = TRUE),
    #  precipitation = precip_ref,
    # temperature = temp_ref,
    time_since_conversion = 5
  )
  
  ## For species, we adapt the initial conditions
  continuous_refs_species_rr <- list(
    # MEAN_depth = depth_ref,
    control_soc_mean_T_ha = soc_ref,
    #  precipitation = precip_ref,
    # temperature = temp_ref,
    time_since_conversion = 25
  )
  
  pdp_treatment_rr <- bootstrap_pdp_cat("NEW_treatment_type2", xgb_rr, X_rr,
                                        n_boot = 500,
                                        continuous_refs = continuous_refs_rr,
                                        categorical_refs = categorical_refs_rr,
                                        cond_sample_vars = c("control_soc_mean_T_ha"),
                                        depth_vals = depth_vals_seq)
  
  pdp_history_rr   <- bootstrap_pdp_cat("History_reclass", xgb_rr, X_rr,
                                        n_boot = 500,
                                        continuous_refs = continuous_refs_rr,
                                        categorical_refs = categorical_refs_rr,
                                        cond_sample_vars = c("control_soc_mean_T_ha"),
                                        depth_vals = depth_vals_seq)
  
  pdp_species_rr   <- bootstrap_pdp_cat("diff_species_class", xgb_rr, X_rr,
                                        n_boot = 500,
                                        continuous_refs = continuous_refs_species_rr,
                                        categorical_refs = categorical_refs_species_rr,
                                        cond_sample_vars = c("control_soc_mean_T_ha"),
                                        depth_vals = depth_vals_seq)
  
  # Combine all PDPs for RR
  pdp_all_rr <- bind_rows(
    pdp_treatment_rr  %>% mutate(variable = "NEW_treatment_type2"),
    pdp_history_rr    %>% mutate(variable = "History_reclass"),
    pdp_species_rr    %>% mutate(variable = "diff_species_class")
  )
  
  # ---------------------------------------------
  # 10. Example plot for RR
  # ---------------------------------------------
  ggplot(pdp_all_rr, aes(x = level, y = exp(pred), ymin = exp(lci), ymax = exp(uci))) +
    geom_point() + geom_errorbar(width = 0.2) +
    facet_wrap(~variable, scales = "free_x") +
    theme_bw() +
    labs(y = "Partial Dependence (RR)")
  
  
  pdp_all_rr$var<-"rr"
  pdp_all$var<-"seq"
  results_table<- bind_rows(pdp_all_rr,pdp_all)
  
  results_table_exp <- results_table %>%
    dplyr::mutate(
      pred_scaled = ifelse(var == "rr", exp(pred), pred),
      lci_scaled  = ifelse(var == "rr", exp(lci), lci),
      uci_scaled  = ifelse(var == "rr", exp(uci), uci)
    )
  
  
  library(ggplot2)
  library(dplyr)
  library(forcats)
  
  # Filtrer "Other" si tu ne veux pas l'afficher
  results_plot <- results_table_exp %>%
    filter(level != "Other") %>%
    # Ordre factor pour ggplot
    mutate(level = fct_rev(fct_inorder(level)))
  
  # Création du forest plot
  ggplot(results_plot, aes(x = level, y = pred_scaled, ymin = lci_scaled, ymax = uci_scaled, color = variable)) +
    geom_pointrange(size = 0.8, fatten = 2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +  # pour ratios
    facet_wrap(~var, scales = "free", ncol = 1, labeller = as_labeller(c(rr = "Ratio", seq = "Sequestration rate (Mg C ha⁻¹ yr⁻¹)"))) +
    coord_flip() +
    scale_color_manual(values = c("NEW_treatment_type2" = "#1f78b4", 
                                  "History_reclass" = "#33a02c",
                                  "diff_species_class" = "#e31a1c")) +
    labs(x = NULL, y = "Effect estimate", color = "Variable") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
  
  
  
  
####
 
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  # Fonction pour transformer une table d'une variable catégorielle
  format_numeric_table <- function(table_stat, cat_var_levels) {
    table_stat %>%
      mutate(
        # Construire des chaînes "moyenne ± sd" pour chaque niveau
        level_means = map2(means, sds, ~ paste0(round(.x, 2), " ± ", round(.y, 2))),
        # Extraire la p-value globale de l'ANOVA
        p.value = map_dbl(anova, ~ .x$p.value[1])
      ) %>%
      # Déplier la liste de level_means en colonnes pour chaque niveau
      mutate(level_means = map(level_means, ~ set_names(.x, cat_var_levels))) %>%
      unnest_wider(level_means) %>%
      select(variable, all_of(cat_var_levels), p.value)
  }
  
  # Exemple pour la première variable catégorielle
  cat_var <- cat_vars[1]
  cat_var_levels <- levels(FULL_sun2_sub_seq[[cat_var]])
  
  numeric_table_clean <- format_numeric_table(numeric_tables_stats[[1]], cat_var_levels)
  
  numeric_table_clean
  
  
  
  
  # Fonction pour formater une table de variable catégorielle
  library(dplyr)
  library(tidyr)
  
  # Fonction pour formater une table catégorielle (en % uniquement)
  format_categorical_table_perc <- function(data, target, cat_var) {
    
    tab <- table(data[[cat_var]], data[[target]])
    
    tab_perc <- prop.table(tab, margin = 2) * 100
    perc_df <- as.data.frame.matrix(tab_perc) |>
      tibble::rownames_to_column("level") |>
      dplyr::mutate(variable = cat_var)
    
    return(perc_df)
  }
  
  # Liste des variables catégorielles
  cat_vars <- c("main_culture2","History_reclass","Grouped_Design")
  
  # Boucle sur toutes les variables
  cat_tables <- lapply(cat_vars, function(v) {
    format_categorical_table_perc(FULL_sun2_sub_seq, "NEW_treatment_type2", v)
  })
  
  # Résultats combinés
  cat_tables_all <- dplyr::bind_rows(cat_tables)
  cat_tables_all
  
  
  