###############################################################################
# MACCA – R3 revision: full rerun of the analyses, in order
#
# Usage: open this file in RStudio, adjust the CONFIG block, then Source.
# Every script reads the outputs of the previous one from OUT_DIR, so a single
# step can be rerun on its own once the earlier steps have been run once.
###############################################################################

# ============================ CONFIG ========================================
# The project folder is the folder containing this file: nothing to edit when
# the MACCA_R3 folder is moved (works with Rscript and with Source in RStudio).
PROJECT_DIR <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) return(normalizePath(dirname(f[1])))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(normalizePath(dirname(p)))
  }
  normalizePath(getwd())
})
# Database: the corrected version deposited on Dataverse is used by default. If it is absent,
# the frozen uncorrected version (30 April 2026) is used and the documented corrections are applied.
RAW_DB <- file.path(PROJECT_DIR, "data", "raw", "MACCA_database_R3_corrected.csv")
if (!file.exists(RAW_DB)) RAW_DB <- file.path(PROJECT_DIR, "data", "raw", "Data_for_analysis_R2_20260430.csv")
RAW_XLSX         <- file.path(PROJECT_DIR, "data", "raw", "MACCA_BDD.xlsx")  # only for step 01
RAW_SHEET        <- NULL
CORRECTIONS_FILE <- file.path(PROJECT_DIR, "data", "data_corrections.csv")
OUT_DIR          <- file.path(PROJECT_DIR, "outputs")
RUN_STEP_01      <- FALSE   # TRUE only to rebuild the database from Excel + WorldClim ("with01")

# Studies excluded at study level (Reviewer #2). Check the ids against the
# numbering used in the Dataverse file before running.
EXCLUDED_STUDIES <- tibble::tribble(
  ~id_article, ~reason,
  34, "Incubation study (ground tree litter incubated in jars): no field-based tree-less control",
  83, "Control defined by distance from trees, not an independent tree-less plot"
)

DEPTH_BREAKS  <- c(0, 15, 30, 45, 55, 75)   # depth classes used in Fig. 4
N_BOOT_IMP    <- 100    # cluster-bootstrap resamples for variable importance
N_BOOT_PDP    <- 200    # cluster-bootstrap resamples for PDP / Table 2 intervals
N_CV_REPEATS  <- 10     # repetitions of the 5-fold study-grouped cross-validation
SEED          <- 594
# ============================================================================

# Terminal options:  Rscript 00_run_all.R quick   -> fast test run (few bootstraps)
#                    Rscript 00_run_all.R with01  -> rebuild the database from Excel + WorldClim (step 01)
args <- commandArgs(trailingOnly = TRUE)
# "only09", "only10" (or "only09" "only10" together) rerun only those steps (need outputs from a previous run)
if ("quick" %in% args) {
  N_BOOT_IMP <- 5; N_BOOT_PDP <- 4; N_CV_REPEATS <- 2; N_LOSO_MAX <- 3
  message("QUICK TEST RUN: results are not for the manuscript")
}
if ("with01" %in% args) { RUN_STEP_01 <- TRUE; RAW_DB <- file.path(OUT_DIR, "Data_for_analysis.csv") }
if ("skip01" %in% args) RUN_STEP_01 <- FALSE   # kept for compatibility
if (!RUN_STEP_01 && !file.exists(RAW_DB)) stop("Database not found: ", RAW_DB)
message("Project folder: ", PROJECT_DIR)
ONLY <- sub("^only", "", grep("^only[0-9]{2}$", args, value = TRUE))
ONLY09 <- length(ONLY) > 0

setwd(PROJECT_DIR)
# Scripts may sit in a R/ sub-folder (recommended) or next to this file
SCRIPT_DIR <- if (file.exists("R/utils.R")) "R" else if (file.exists("utils.R")) "." else
  stop("utils.R not found in ", PROJECT_DIR, " nor in ", file.path(PROJECT_DIR, "R"),
       ".\nFiles present: ", paste(list.files(PROJECT_DIR, recursive = TRUE), collapse = ", "))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(SCRIPT_DIR, "utils.R"))


steps <- file.path(SCRIPT_DIR, c(
  if (RUN_STEP_01) "01_load_data.R",
  "02_prepare_datasets.R",      # filters, study exclusions, SD imputation, effect sizes
  "03_meta_analysis.R",         # three-level models, heterogeneity, bias, Table 2 (metafor)
  "04_ml_models.R",             # outliers, study-grouped CV, importance, SHAP interactions
  "05_partial_dependence.R",    # PDPs restricted to observed ranges + tests on observed data
  "06_categorical_effects.R",   # Table 2 (XGBoost standardized predictions)
  "07_data_support.R",          # coverage by depth x time (Reviewer #2)
  "08_moderator_tests.R",       # moderator tests on observed effect sizes
  "09_robustness.R",            # leave-one-study-out, thresholds, confounding
  "10_figures.R",               # main and supplementary figures (ASD rules)
  "11_prisma.R"                 # PRISMA flow diagram (Fig. S1) from data/prisma_counts.csv
))
if (ONLY09) steps <- steps[substr(basename(steps), 1, 2) %in% ONLY]
# a full run starts a fresh numbers table; "only09" updates the existing one
if (!ONLY09 && file.exists(dir_out("numbers_for_manuscript.csv"))) invisible(file.remove(dir_out("numbers_for_manuscript.csv")))
missing_scripts <- steps[!file.exists(steps)]
if (length(missing_scripts)) stop("Missing script(s): ", paste(missing_scripts, collapse = ", "))

for (s in steps) {
  message("\n==================== ", s, " ====================")
  set.seed(SEED)
  t0 <- Sys.time()
  source(s, echo = FALSE)
  message("   done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}

writeLines(capture.output(sessionInfo()), dir_out("sessionInfo.txt"))
message("\nAll numbers quoted in the manuscript: ", dir_out("numbers_for_manuscript.csv"))
