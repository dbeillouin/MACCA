###############################################################################
# 02 – Build the two analysis datasets (response ratios, storage rates)
#
# Changes vs. former 2_RR.R / 2_Seq_rate.R:
#  * study-level exclusions applied HERE, before any analysis (Reviewer #2);
#  * identical recoding for both datasets (history, design, crops, species);
#  * no row-index subsetting anywhere; outliers are handled by ID in step 04;
#  * the two analysis datasets are written to disk and are the only inputs of
#    all later steps (deposit them on Dataverse).
###############################################################################
suppressPackageStartupMessages(library(metafor))

raw <- read_csv(RAW_DB, col_types = cols(.default = col_guess()), guess_max = 1e5,
                show_col_types = FALSE, name_repair = "unique_quiet")
# experiment identifier: the database column is "id_expérimentation" (accent);
# renamed once here to an ASCII name used by all later steps
norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(iconv(enc2utf8(x), "UTF-8", "ASCII", sub = "")))
exp_col <- names(raw)[norm_name(names(raw)) %in% c("idexprimentation", "idexperimentation", "idexperiment")]
if (!length(exp_col)) stop("Experiment identifier column (id_exp\u00e9rimentation) not found in ", RAW_DB)
raw <- raw %>% rename(id_experiment = all_of(exp_col[1]))

# numeric columns: accept decimal commas and stray text (e.g. "0,35", " 12 ")
to_num <- function(x) if (is.numeric(x)) x else suppressWarnings(as.numeric(gsub(",", ".", trimws(as.character(x)))))
num_cols <- intersect(c("time_since_conversion", "treatment_NB_species", "control_NB_species",
                        "treatment_replicate_nb", "control_replicate_nb",
                        "treatment_soc_mean_T_ha", "control_soc_mean_T_ha",
                        "treatment_soc_sd_T_ha", "control_soc_sd_T_ha",
                        "delta_stock_T", "delta_stock_C.yr", "precipitation", "temperature",
                        "MEAN_depth", "soil_depth_end", "id_article"), names(raw))
raw <- raw %>% mutate(across(all_of(num_cols), to_num), ID = as.character(ID))

required <- c("ID", "id_article", "id_experiment", "experiment_type", "control_type", "Keep_ratio", "Keep_rate",
              "time_since_conversion", "main_culture", "NEW_treatment_type", "treatment_NB_species", "control_NB_species",
              "Design", "history_C", "history_T", "treatment_replicate_nb", "control_replicate_nb",
              "treatment_soc_mean_T_ha", "control_soc_mean_T_ha", "treatment_soc_sd_T_ha", "control_soc_sd_T_ha",
              "delta_stock_T", "delta_stock_C.yr", "language", "precipitation", "temperature", "MEAN_depth")
missing <- setdiff(required, names(raw))
if (length(missing)) stop("Columns missing from ", RAW_DB, ": ", paste(missing, collapse = ", "),
                          "\nCheck RAW_SHEET (wrong sheet?) or column names in the Excel file.")

base_filter <- function(df) df %>% filter(experiment_type == "LMP", control_type == "Full sun")

# ---------------------------------------------------------------------------
# Documented corrections (data/data_corrections.csv)
# ---------------------------------------------------------------------------
if (!exists("CORRECTIONS_FILE")) CORRECTIONS_FILE <- file.path(PROJECT_DIR, "data", "data_corrections.csv")
ALREADY_CORRECTED <- grepl("corrected", basename(RAW_DB), ignore.case = TRUE)
if (ALREADY_CORRECTED) {
  message("Database already corrected (", basename(RAW_DB), "): documented corrections are not re-applied")
} else {
  raw <- apply_corrections(raw, CORRECTIONS_FILE)
}
clog <- attr(raw, "corrections_log")
if (!is.null(clog) && nrow(clog)) {
  write_csv(clog, dir_out("data", "corrections_applied.csv"))
  record("n_rows_corrected", paste0(sum(clog$rows_changed > 0), " corrections, ", sum(clog$rows_changed), " cells"), "2.2, response letter")
  if (any(clog$rows_changed == 0)) warning("Some corrections matched no row: see outputs/data/corrections_applied.csv")
}
attr(raw, "corrections_log") <- NULL
# corrected copy of the full database (to replace the Dataverse file)
write_csv(raw, dir_out("data", "MACCA_database_R3_corrected.csv"), na = "NA")

# ---------------------------------------------------------------------------
# Data checks: coordinates and climate
# ---------------------------------------------------------------------------
if (all(c("X_(WGS84)", "Y_(WGS84)") %in% names(raw))) {
  lon0 <- raw$`X_(WGS84)`; lat0 <- raw$`Y_(WGS84)`
  raw <- raw %>% mutate(lon = rescale_coord(`X_(WGS84)`, -120, -30), lat = rescale_coord(`Y_(WGS84)`, -60, 35))
  fixed <- raw %>% mutate(lon0 = as.character(lon0), lat0 = as.character(lat0)) %>%
    filter(suppressWarnings(as.numeric(gsub(",", ".", lon0))) != lon | suppressWarnings(as.numeric(gsub(",", ".", lat0))) != lat |
             is.na(lon) | is.na(lat)) %>%
    transmute(id_article, country = if ("country" %in% names(.)) country else NA, lon_database = lon0, lat_database = lat0, lon_used = lon, lat_used = lat) %>% distinct()
  write_csv(fixed, dir_out("data", "coordinates_to_check.csv"))
  if (nrow(fixed)) warning(nrow(fixed), " site(s) with missing or rescaled coordinates: see outputs/data/coordinates_to_check.csv")
}
# climate values identical to the database mean = filled by the mean in step 01
if (all(c("temperature", "precipitation") %in% names(raw))) {
  tmean <- raw %>% filter(experiment_type == "LMP") %>% count(temperature, precipitation) %>%
    filter(abs(temperature - round(temperature, 2)) > 1e-6)
  suspicious <- raw %>% semi_join(tmean, by = c("temperature", "precipitation")) %>%
    distinct(across(any_of(c("id_article", "country", "temperature", "precipitation", "altitude"))))
  write_csv(suspicious, dir_out("data", "climate_to_check.csv"))
  if (nrow(suspicious)) warning("Climate values that look imputed by the dataset mean: see outputs/data/climate_to_check.csv")
}

# ---------------------------------------------------------------------------
# Study-level exclusions (logged)
# ---------------------------------------------------------------------------
apply_exclusions <- function(df, label) {
  log <- df %>%
    filter(id_article %in% EXCLUDED_STUDIES$id_article) %>%
    count(id_article, name = "n_observations_removed") %>%
    mutate(dataset = label)
  list(data = df %>% filter(!id_article %in% EXCLUDED_STUDIES$id_article), log = log)
}

fix_language <- function(x) if_else(!is.na(x) & grepl("Spanish", x), "Spanish", x)

# ===========================================================================
# A. Response ratios
# ===========================================================================
rr0 <- raw %>% base_filter() %>% filter(Keep_ratio == "YES")
ex_rr <- apply_exclusions(rr0, "response_ratio")

rr <- ex_rr$data %>%
  mutate(time_since_conversion_imputed = is.na(time_since_conversion),
         time_since_conversion = if_else(is.na(time_since_conversion),
                                         median(time_since_conversion, na.rm = TRUE),
                                         time_since_conversion),
         language = fix_language(language)) %>%
  recode_common() %>%
  impute_sd("treatment_soc_mean_T_ha", "treatment_soc_sd_T_ha", inflate = 1) %>%
  impute_sd("control_soc_mean_T_ha",   "control_soc_sd_T_ha",   inflate = 1)

es <- escalc(measure = "ROM",
             n1i = rr$treatment_replicate_nb, n2i = rr$control_replicate_nb,
             m1i = rr$treatment_soc_mean_T_ha, m2i = rr$control_soc_mean_T_ha,
             sd1i = rr$treatment_soc_sd_T_ha,  sd2i = rr$control_soc_sd_T_ha)
rr$yi <- as.numeric(es$yi); rr$vi <- as.numeric(es$vi)
rr <- rr %>% filter(!is.na(yi), !is.na(vi))

# ===========================================================================
# B. Storage rates
# ===========================================================================
sq0 <- raw %>%
  base_filter() %>%
  filter(!is.na(delta_stock_T) | Keep_rate == "YES")
ex_sq <- apply_exclusions(sq0, "storage_rate")

sq <- ex_sq$data %>%
  mutate(delta_stock_C.yr = as.numeric(gsub(",", ".", delta_stock_C.yr)),
         delta_stock_T    = as.numeric(gsub(",", ".", delta_stock_T)))
# known data-entry correction (kept from the former script)
sq$control_soc_mean_T_ha[sq$ID == "442"] <- 40.5

sq <- sq %>%
  mutate(delta_stock_T = if_else(!is.na(delta_stock_T) & !is.na(delta_stock_C.yr),
                                 delta_stock_T - delta_stock_C.yr, delta_stock_T),
         history_C = if_else(!is.na(delta_stock_T), history_T, history_C),
         language  = fix_language(language)) %>%
  recode_common() %>%
  # NOTE: the former storage-rate script inflated imputed CVs by 1.5. Kept for
  # comparability; describe it in Methods or set inflate = 1 for both metrics.
  impute_sd("treatment_soc_mean_T_ha", "treatment_soc_sd_T_ha", inflate = 1.5) %>%
  impute_sd("control_soc_mean_T_ha",   "control_soc_sd_T_ha",   inflate = 1.5) %>%
  mutate(seq_rate    = if_else(is.na(delta_stock_T),
                               (treatment_soc_mean_T_ha - control_soc_mean_T_ha) / time_since_conversion,
                               delta_stock_T),
         seq_rate_sd = sqrt(treatment_soc_sd_T_ha^2 + control_soc_sd_T_ha^2) / time_since_conversion,
         seq_rate_vi = seq_rate_sd^2) %>%
  filter(!is.na(seq_rate), !is.na(seq_rate_vi), is.finite(seq_rate))

# Every observation needs a unique identifier (outliers are removed by ID).
# Some rows of the database have no ID: they receive a generated one.
fill_id <- function(df, prefix) df %>%
  mutate(ID = as.character(ID),
         ID = if_else(is.na(ID) | ID == "" | duplicated(ID), paste0(prefix, "_row", row_number()), ID))
rr <- fill_id(rr, "RR"); sq <- fill_id(sq, "SR")
record("n_obs_without_original_ID", paste0("RR ", sum(grepl("^RR_row", rr$ID)), "; storage rate ", sum(grepl("^SR_row", sq$ID))),
       "data note", "IDs generated; fill them in the database")

# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------
write_csv(rr, dir_out("data", "analysis_dataset_RR.csv"))
write_csv(sq, dir_out("data", "analysis_dataset_storage_rate.csv"))
write_csv(left_join(EXCLUDED_STUDIES, bind_rows(ex_rr$log, ex_sq$log), by = "id_article"),
          dir_out("data", "excluded_studies.csv"))

if (nrow(bind_rows(ex_rr$log, ex_sq$log)) == 0)
  warning("None of the EXCLUDED_STUDIES ids were found: check the id_article numbering.")

# ---------------------------------------------------------------------------
# Numbers for Methods / Results 3.1 / Fig. 2 caption
# ---------------------------------------------------------------------------
n_studies <- function(d) n_distinct(d$id_article)
record("n_obs_RR", nrow(rr), "2.1, 3.1, Fig. 2, Abstract")
record("n_studies_RR", n_studies(rr), "2.1, 3.1, Fig. 2, Abstract")
record("n_obs_storage_rate", nrow(sq), "2.1, 3.1, Fig. 2, Abstract")
record("n_studies_storage_rate", n_studies(sq), "2.1, 3.1, Fig. 2, Abstract")
record("n_studies_union", n_distinct(c(rr$id_article, sq$id_article)), "3.1")
record("n_obs_removed_by_study_exclusion",
       paste0("RR ", sum(ex_rr$log$n_observations_removed), "; storage rate ",
              sum(ex_sq$log$n_observations_removed)), "2.1 (new exclusion paragraph)")
lang <- bind_rows(select(rr, id_article, language), select(sq, id_article, language)) %>% distinct()
record("pct_spanish_studies", 100 * mean(lang$language == "Spanish", na.rm = TRUE), "3.1")
record("pct_obs_sd_imputed_RR",
       100 * mean(rr$treatment_soc_sd_T_ha_imputed | rr$control_soc_sd_T_ha_imputed), "2.2 (12%)")
record("n_obs_time_imputed_RR", sum(rr$time_since_conversion_imputed), "2.2",
       "median imputation of time since conversion: to describe in Methods")

message("RR: ", nrow(rr), " obs / ", n_studies(rr), " studies;  storage rate: ",
        nrow(sq), " obs / ", n_studies(sq), " studies")
