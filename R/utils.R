###############################################################################
# MACCA – shared helpers for the R3 revision pipeline
# All functions are pure (no global side effects) so that every script can be
# rerun independently from the frozen analysis datasets.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(stringr)
  library(purrr)
})

# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------
dir_out <- function(...) {
  p <- file.path(OUT_DIR, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

# Every number quoted in the manuscript is appended to a single table:
# numbers_for_manuscript.csv (key, value, section, note)
record <- function(key, value, section = NA_character_, note = NA_character_) {
  f <- dir_out("numbers_for_manuscript.csv")
  row <- tibble(key = key,
                value = if (is.numeric(value)) format(signif(value, 4), scientific = FALSE) else as.character(value),
                section = section, note = note)
  if (file.exists(f)) {
    old <- read_csv(f, col_types = cols(.default = "c"))
    old <- old[old$key != key, , drop = FALSE]
    row <- bind_rows(old, row)
  }
  write_csv(row, f)
  invisible(value)
}

# ---------------------------------------------------------------------------
# Harmonised recoding (identical for response ratios and storage rates)
# ---------------------------------------------------------------------------
group_design <- function(Design) {
  case_when(
    is.na(Design)                 ~ "Not specified",
    Design %in% c("RCBD", "RCT")  ~ "Randomized Designs",
    Design %in% c("BACI", "BA")   ~ "Before-After Designs",
    Design == "CI"                ~ "Control-Impact Designs",
    TRUE                          ~ "Other"
  )
}

# NOTE (R3): the former scripts classified "pasture pasture" as Forest in the
# RR script (first matching case) and as Grassland in the storage-rate script.
# It is now Grassland in both datasets.
classify_history <- function(history_C, history_T) {
  h <- paste(tolower(history_C), tolower(history_T))
  case_when(
    h == "cropland cropland"                              ~ "Cropland",
    h == "forest forest"                                  ~ "Forest",
    h %in% c("grassland grassland", "pasture pasture")    ~ "Grassland",
    TRUE                                                  ~ "Unknown/mixed"
  )
}

recode_common <- function(df) {
  df %>%
    mutate(
      main_culture2 = case_when(
        grepl("Cocoa|Coffee|cocoa", main_culture) ~ "Cocoa/Coffee",
        grepl("Beans|Maize", main_culture)        ~ "Beans/Maize",
        grepl("Banana", main_culture)             ~ "Others",
        TRUE                                      ~ as.character(main_culture)
      ),
      NEW_treatment_type  = fct_lump_min(factor(NEW_treatment_type), 10),
      NEW_treatment_type2 = as.character(NEW_treatment_type),
      NEW_treatment_type2 = if_else(NEW_treatment_type2 %in% c("Alley cropping", "Hedgerow"),
                                    "Alley/Hedgerow", NEW_treatment_type2),
      diff_species = treatment_NB_species - suppressWarnings(as.numeric(control_NB_species)),
      diff_species_class = as.character(cut(diff_species, breaks = c(-5, 1.1, 25),
                                            labels = c("A:<0/+1", "C+2+"), right = FALSE)),
      diff_species_class = if_else(is.na(diff_species_class), "C+2+", diff_species_class),
      Grouped_Design  = group_design(Design),
      History_reclass = classify_history(history_C, history_T),
      depth_group = cut(MEAN_depth, breaks = DEPTH_BREAKS, include.lowest = TRUE, right = FALSE)
    )
}

# Missing SDs imputed from the mean coefficient of variation of each design.
# `inflate` reproduces the former storage-rate script (CV x 1.5); set to 1 if
# the Methods should describe plain mean CVs.
impute_sd <- function(df, mean_var, sd_var, group_var = "Grouped_Design", inflate = 1) {
  tab <- df %>%
    group_by(.data[[group_var]]) %>%
    summarise(cv_m = mean(.data[[sd_var]], na.rm = TRUE) / mean(.data[[mean_var]], na.rm = TRUE) * inflate,
              .groups = "drop")
  df %>%
    left_join(tab, by = group_var) %>%
    mutate("{sd_var}_imputed" := is.na(.data[[sd_var]]),
           !!sd_var := if_else(is.na(.data[[sd_var]]), cv_m * .data[[mean_var]], .data[[sd_var]])) %>%
    select(-cv_m)
}

# ---------------------------------------------------------------------------
# Machine-learning helpers
# ---------------------------------------------------------------------------
PREDICTORS <- c("NEW_treatment_type2", "diff_species_class", "History_reclass",
                "control_soc_mean_T_ha", "precipitation", "temperature",
                "main_culture2", "MEAN_depth", "time_since_conversion",
                "Grouped_Design")

CATEGORICAL <- c("NEW_treatment_type2", "diff_species_class", "History_reclass",
                 "main_culture2", "Grouped_Design")

# Fixed factor levels so that every design matrix (folds, bootstrap samples,
# counterfactual data) has exactly the same columns.
freeze_levels <- function(df) {
  for (v in CATEGORICAL) df[[v]] <- factor(df[[v]], levels = sort(unique(as.character(df[[v]]))))
  df
}

make_X <- function(df) {
  f <- as.formula(paste("~", paste(PREDICTORS, collapse = " + ")))
  mf <- model.frame(f, data = df, na.action = na.pass)
  X <- model.matrix(f, data = mf,
                    contrasts.arg = lapply(df[CATEGORICAL], function(x) contrasts(x, contrasts = FALSE)))
  X[, colnames(X) != "(Intercept)", drop = FALSE]
}

feature_origin <- function(feat) {
  out <- feat
  for (v in PREDICTORS) out[startsWith(feat, v)] <- v
  out
}

XGB_PARAMS <- list(objective = "reg:squarederror", max_depth = 4, eta = 0.05,
                   subsample = 0.9, colsample_bytree = 0.8)

fit_xgb <- function(X, y, nrounds, w = NULL) {
  d <- xgboost::xgb.DMatrix(X, label = y, weight = w)
  xgboost::xgb.train(params = XGB_PARAMS, data = d, nrounds = nrounds, verbose = 0)
}

# Folds defined by study: all observations of a study are in the same fold, so
# that performance is estimated on studies never seen during training.
study_folds <- function(study, k = 5, seed = 1) {
  set.seed(seed)
  s <- unique(study)
  f <- sample(rep_len(seq_len(k), length(s)))
  f[match(study, s)]
}

r2 <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)

# Cluster bootstrap: resample studies with replacement
boot_rows_by_study <- function(study) {
  s <- unique(study)
  pick <- sample(s, length(s), replace = TRUE)
  unlist(lapply(pick, function(z) which(study == z)), use.names = FALSE)
}

# A "profile" = one agroforestry/control comparison sampled at several depths
# (same article, site, age, system, crop and history). Used for an intermediate
# cross-validation where all layers of a profile are kept in the same fold.
make_profile <- function(d) {
  xy <- intersect(c("X_(WGS84)", "Y_(WGS84)"), names(d))
  site <- if (length(xy) == 2)
    paste(d$id_article, round(suppressWarnings(as.numeric(d[[xy[1]]])), 2), round(suppressWarnings(as.numeric(d[[xy[2]]])), 2))
  else as.character(d$id_article)
  paste(site, d$time_since_conversion, d$NEW_treatment_type2, d$main_culture2, d$History_reclass, round(d$precipitation))
}

# Coordinates that lost their decimal separator in the database (e.g. -836333
# for -83.6333, 125 for 12.5) are rescaled by powers of ten until they fall in
# the study region (longitude -120 to -30, latitude -60 to 35). Every change is
# logged so that the database can be corrected at source.
rescale_coord <- function(x, lo, hi) {
  x <- suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
  out <- x
  for (i in which(!is.na(x) & (x < lo | x > hi))) {
    v <- x[i]; k <- 0
    while (!is.na(v) && (v < lo || v > hi) && k < 10) { v <- v / 10; k <- k + 1 }
    out[i] <- if (v >= lo && v <= hi) v else NA
  }
  out
}

# Documented data corrections (data/data_corrections.csv). Only rows with
# status == "applied" change the data; "to_verify" rows are reported and used
# for sensitivity analyses in step 09.
apply_corrections <- function(df, file) {
  if (is.null(file) || !file.exists(file)) { warning("No corrections file found: ", file); return(df) }
  corr <- read_csv(file, col_types = cols(.default = "c"), na = character())
  log <- list()
  for (i in seq_len(nrow(corr))) {
    cr <- corr[i, ]
    if (cr$status != "applied" || !(cr$field %in% names(df))) next
    cur <- df[[cr$match_field]]
    hit <- as.character(df$id_article) == cr$id_article &
      (cr$match_country == "" | as.character(df$country) == cr$match_country)
    if (cr$match_value == "NA") hit <- hit & is.na(cur)
    else if (cr$match_value != "*") {
      num <- suppressWarnings(as.numeric(cr$match_value))
      hit <- hit & if (!is.na(num)) abs(suppressWarnings(as.numeric(cur)) - num) < 1e-6 else as.character(cur) == cr$match_value
    }
    hit[is.na(hit)] <- FALSE
    if (any(hit)) {
      nv <- if (is.numeric(df[[cr$field]])) as.numeric(cr$new_value) else cr$new_value
      df[[cr$field]][hit] <- nv
    }
    log[[length(log) + 1]] <- tibble(id_article = cr$id_article, country = cr$match_country, field = cr$field,
                                     new_value = cr$new_value, rows_changed = sum(hit), reason = cr$reason)
  }
  attr(df, "corrections_log") <- bind_rows(log)
  df
}
flag_to_verify <- function(df, file, field = "temperature") {
  if (is.null(file) || !file.exists(file)) return(rep(FALSE, nrow(df)))
  corr <- read_csv(file, col_types = cols(.default = "c"), na = character()) %>% filter(status == "to_verify", field == !!field)
  f <- rep(FALSE, nrow(df))
  for (i in seq_len(nrow(corr))) {
    cr <- corr[i, ]
    h <- as.character(df$id_article) == cr$id_article & as.character(df$country) == cr$match_country
    num <- suppressWarnings(as.numeric(cr$match_value))
    if (!is.na(num)) h <- h & abs(suppressWarnings(as.numeric(df[[field]])) - num) < 1e-6
    f <- f | (!is.na(h) & h)
  }
  f
}
