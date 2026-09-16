###############################################################################
# Project   : MACCA Database – Carbon Stock Analysis in Agroforestry Systems
# Author    : Damien Beillouin
# Date      : 30 July 2024
# Purpose   : Clean and enrich MACCA dataset for reproducible analysis of soil
#             organic carbon (SOC) in agroforestry systems.
# FAIR note : This script prepares a standardized dataset (FAIR principles:
#             Findable, Accessible, Interoperable, Reusable).
# Output    : Cleaned dataset saved as "Data_for_analysis.csv"
###############################################################################

#######################--------------------------------------------------------
# PART I : Load required libraries
#######################--------------------------------------------------------

library(magrittr)   # Piping (%>%)
library(tidyverse)  # Data wrangling and visualization
library(readxl)     # Read Excel files
library(naniar)     # Handle missing values
library(readr)      # Read flat files (CSV, TSV)
library(psych)      # Statistical utilities
library(skimr)      # Data summaries
library(ggpubr)     # Publication-ready graphics
library(metafor)    # Meta-analysis methods
library(cowplot)    # Combine ggplots
library(stringr)    # String operations
if (requireNamespace("gsheet", quietly = TRUE)) library(gsheet)  # Import from Google Sheets
library(raster)     # Spatial rasters
library(sp)         # Spatial objects
library(tidyr)      # separate()

# R3: paths are set here instead of setwd() calls scattered in the script
WORLDCLIM_DIR <- Sys.getenv("MACCA_WORLDCLIM_DIR", file.path(dirname(PROJECT_DIR), "bio_10m_bil"))   # WorldClim 10' bioclim (.bil)
ET0_DIR       <- Sys.getenv("MACCA_ET0_DIR", file.path(dirname(PROJECT_DIR), "et0_yr"))              # Global-PET annual (.tif)

#######################--------------------------------------------------------
# PART II : Import and initial cleaning
#######################--------------------------------------------------------

# Import data from Google Sheets
# R3: read the frozen Excel export when RAW_XLSX is set (recommended),
# otherwise fall back to the live Google Sheet as before.
if (exists("RAW_XLSX") && !is.null(RAW_XLSX) && file.exists(path.expand(RAW_XLSX))) {
  xl <- path.expand(RAW_XLSX)
  sheets <- readxl::excel_sheets(xl)
  message("Sheets found in ", basename(xl), ": ", paste(sheets, collapse = " | "))
  key_col <- "treatment_soc_mean_T_ha"
  # locate the sheet and the header row that contain the expected column names
  find_header <- function(sh) {
    top <- suppressMessages(readxl::read_excel(xl, sheet = sh, col_names = FALSE, n_max = 30,
                                               col_types = "text", .name_repair = "minimal"))
    hit <- which(apply(top, 1, function(r) any(trimws(r) == key_col, na.rm = TRUE)))
    list(sheet = sh, row = if (length(hit)) hit[1] else NA_integer_, top = top)
  }
  cand <- if (exists("RAW_SHEET") && !is.null(RAW_SHEET)) RAW_SHEET else sheets
  found <- lapply(cand, find_header)
  ok <- Filter(function(z) !is.na(z$row), found)
  if (!length(ok)) {
    for (z in found) {
      message("\n--- sheet '", z$sheet, "': first non-empty cells of rows 1-3")
      for (k in seq_len(min(3, nrow(z$top)))) {
        cells <- na.omit(unlist(z$top[k, ])); message("  row ", k, ": ", paste(head(cells, 12), collapse = " | "))
      }
    }
    stop("Column '", key_col, "' not found in any sheet of ", basename(xl),
         ". The Excel file does not use the column names of the MACCA Google Sheet ",
         "(see the header rows printed above).")
  }
  sheet <- ok[[1]]$sheet; skip <- ok[[1]]$row - 1
  message("Reading sheet '", sheet, "' with header on row ", skip + 1)
  # read everything as text, then guess types like a CSV download would
  Data <- readxl::read_excel(xl, sheet = sheet, skip = skip, col_types = "text") %>%
    as.data.frame() %>%
    utils::type.convert(as.is = TRUE)
  names(Data) <- trimws(names(Data))
} else {
  stop("Excel export of the database not found: ", RAW_XLSX,
       "\nStep 01 is only needed to rebuild the database; the analyses start from data/raw/Data_for_analysis_R2_20260430.csv.")
}

# R3: check all columns used by steps 01-02 before any computation
needed <- c("treatment_soc_mean_T_ha", "control_soc_mean_T_ha", "treatment_soc_sd_T_ha", "control_soc_sd_T_ha",
            "delta_stock_T", "delta_stock_C.yr", "Bulk_density_C", "Bulk_density_T", "Coarse_fragment_C", "Coarse_fragment_T",
            "temperature", "precipitation", "altitude", "X_(WGS84)", "Y_(WGS84)", "soil_depth_measured",
            "time_since_conversion", "treatment_replicate_nb", "control_replicate_nb", "history_C", "history_T",
            "ID", "id_article", "experiment_type", "control_type", "Keep_ratio", "Keep_rate",
            "main_culture", "NEW_treatment_type", "treatment_NB_species", "control_NB_species", "Design", "language")
absent <- setdiff(needed, names(Data))
if (!any(gsub("[^a-z0-9]", "", tolower(iconv(enc2utf8(names(Data)), "UTF-8", "ASCII", sub = ""))) %in% c("idexprimentation", "idexperimentation", "idexperiment")))
  absent <- c(absent, "id_exp\u00e9rimentation")
if (length(absent)) {
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(OUT_DIR, "colonnes_excel.txt")
  writeLines(c("COLONNES ATTENDUES MAIS ABSENTES :", absent, "",
               paste0("COLONNES PRESENTES DANS LE FICHIER (", ncol(Data), ") :"), names(Data)), f)
  stop("\n", length(absent), " column(s) used by the scripts are missing from the Excel file:\n  ",
       paste(absent, collapse = ", "),
       "\nFull list of the file's columns written to: ", f)
}

# Convert carbon stocks from g/m² to t/ha
Data <- Data %>%
  mutate(
    treatment_soc_mean_T_ha = treatment_soc_mean_T_ha / 100,
    control_soc_mean_T_ha   = control_soc_mean_T_ha / 100,
    delta_stock_T           = as.numeric(gsub(",", ".", delta_stock_T)),
    Bulk_density_C          = as.numeric(Bulk_density_C),
    Bulk_density_T          = as.numeric(Bulk_density_T),
    Coarse_fragment_C       = as.numeric(Coarse_fragment_C),
    Coarse_fragment_T       = as.numeric(Coarse_fragment_T)
  )

#######################--------------------------------------------------------
# PART III : Climate data processing and extraction
#######################--------------------------------------------------------

# Replace ranges (e.g., "20–25") by their mean value
process_range_values <- function(x) {
  x <- gsub(",", ".", x)
  x <- gsub("–", "-", x)
  ranges <- strsplit(x, "-")
  sapply(ranges, function(r) mean(as.numeric(r)))
}

Data <- Data %>%
  mutate(
    temperature   = process_range_values(temperature),
    precipitation = process_range_values(precipitation),
    altitude      = process_range_values(altitude)
  )

# Load WorldClim .bil files
files <- list.files(WORLDCLIM_DIR, pattern = "\\.bil$", full.names = TRUE)
bil_data <- lapply(files, raster)

# Load evapotranspiration raster
files <- list.files(ET0_DIR, pattern = "\\.tif$", full.names = TRUE)
imported_raster <- raster(files)

# Stack layers and assign names
s <- stack(bil_data)
names(s) <- c(
  "Annual_Mean_Temperature", "Mean_Diurnal_Range", "Isothermality",
  "Temperature_Seasonality", "Max_Temp_Warm", "Min_Temp_Cold",
  "Temp_Annual_Range", "Mean_Temp_Wet", "Mean_Temp_Dry",
  "Mean_Temp_Warm", "Mean_Temp_Cold", "Annual_Precip", "Precip_Wettest",
  "Precip_Driest", "Precip_Season", "Precip_Wet_Quarter",
  "Precip_Dry_Quarter", "Precip_Warm_Quarter", "Precip_Cold_Quarter"
)

# Create spatial points from coordinates
coords <- cbind(as.numeric(Data$`X_(WGS84)`), as.numeric(Data$`Y_(WGS84)`))
coords[is.na(coords)] <- 1
sp_points <- SpatialPoints(coords)

# Extract climate data
Data$Temperature_Wordclim    <- raster::extract(s[[1]], sp_points)
Data$Precipitation_Wordclim  <- raster::extract(s[[12]], sp_points)
Data$ETP_Wordclim            <- raster::extract(imported_raster, sp_points)

# Fill missing values with extracted or summary values
Data$temperature   <- ifelse(is.na(Data$temperature), Data$Temperature_Wordclim, Data$temperature)
Data$temperature   <- ifelse(is.na(Data$temperature), mean(Data$temperature, na.rm = TRUE), Data$temperature)
Data$precipitation <- ifelse(is.na(Data$precipitation), Data$Precipitation_Wordclim, Data$precipitation)
Data$precipitation <- ifelse(is.na(Data$precipitation), mean(Data$precipitation, na.rm = TRUE), Data$precipitation)
Data$ETP_Wordclim  <- ifelse(is.na(Data$ETP_Wordclim), median(Data$ETP_Wordclim, na.rm = TRUE), Data$ETP_Wordclim)
Data$altitude      <- ifelse(is.na(Data$altitude), mean(Data$altitude, na.rm = TRUE), Data$altitude)

#######################--------------------------------------------------------
# PART IV : FAO climate classification
#######################--------------------------------------------------------

classify_climate <- function(MAT, elevation, MAP, PET) {
  if (MAT > 18 & elevation < 1000 & MAP >= 2000) {
    "Tropical Wet"
  } else if (MAT > 18 & elevation < 1000 & MAP >= 1000 & MAP < 2000) {
    "Tropical Moist"
  } else if (MAT > 18 & elevation < 1000 & MAP < 1000) {
    "Tropical Dry"
  } else if (MAT > 18 & elevation >= 1000) {
    "Tropical Montane"
  } else if (MAT > 10 & MAT <= 18 & MAP > PET) {
    "Warm Temperate Moist"
  } else if (MAT > 10 & MAT <= 18 & MAP < PET) {
    "Warm Temperate Dry"
  } else if (MAT > 0 & MAT <= 10 & MAP > PET) {
    "Cool Temperate Moist"
  } else if (MAT > 0 & MAT <= 10 & MAP < PET) {
    "Cool Temperate Dry"
  } else if (MAT > -10 & MAT <= 0 & MAP > PET) {
    "Boreal Moist"
  } else if (MAT > -10 & MAT <= 0 & MAP < PET) {
    "Boreal Dry"
  } else if (MAT <= -10 & MAP > PET) {
    "Polar Moist"
  } else if (MAT <= -10 & MAP < PET) {
    "Polar Dry"
  } else {
    "Unknown Climate"
  }
}

Data$FAO_Climate <- mapply(classify_climate,
                           MAT = Data$temperature,
                           elevation = Data$altitude,
                           MAP = Data$precipitation,
                           PET = Data$ETP_Wordclim)

#######################--------------------------------------------------------
# PART V : Soil depth processing
#######################--------------------------------------------------------

Data <- separate(Data, soil_depth_measured, into = c("soil_depth_start", "soil_depth_end"), sep = "-")
Data$soil_depth_start <- as.numeric(gsub(",", ".", Data$soil_depth_start))
Data$soil_depth_end   <- as.numeric(gsub(",", ".", Data$soil_depth_end))
Data$MEAN_depth       <- (Data$soil_depth_start + Data$soil_depth_end) / 2
Data$depth_diff       <- Data$soil_depth_end - Data$soil_depth_start
Data$Depthfrom_zero   <- Data$depth_diff == 2 * Data$MEAN_depth
Data$MEAN_depthD      <- cut(Data$MEAN_depth, breaks = c(0, 5, 10, 15, 25, 50, 100))

#######################--------------------------------------------------------
# PART VI : Time since conversion
#######################--------------------------------------------------------

Data <- Data %>%
  mutate(
    time_since_conversion = gsub(",", ".", time_since_conversion),
    time_since_conversion = sapply(strsplit(as.character(time_since_conversion), "-"), function(x) {
      x <- unlist(x)
      if (any(grepl("inf:", x))) {
        as.numeric(gsub("inf:", "", x)) - 1
      } else if (any(grepl("sup:", x))) {
        as.numeric(gsub("sup:", "", x)) + 1
      } else {
        mean(as.numeric(x))
      }
    })
  )

#######################--------------------------------------------------------
# PART VII : Replicates and SOC values
#######################--------------------------------------------------------

# Convert to numeric
Data$control_replicate_nb   <- as.numeric(Data$control_replicate_nb)
Data$treatment_replicate_nb <- as.numeric(Data$treatment_replicate_nb)

# Fill missing replicate numbers
Data$treatment_replicate_nb <- ifelse(is.na(Data$treatment_replicate_nb),
                                      Data$control_replicate_nb,
                                      Data$treatment_replicate_nb)
Data$control_replicate_nb   <- ifelse(is.na(Data$control_replicate_nb),
                                      Data$treatment_replicate_nb,
                                      Data$control_replicate_nb)
Data$treatment_replicate_nb <- ifelse(is.na(Data$treatment_replicate_nb), 2, Data$treatment_replicate_nb)
Data$control_replicate_nb   <- ifelse(is.na(Data$control_replicate_nb), 2, Data$control_replicate_nb)

# Convert SOC values
Data <- Data %>%
  mutate(
    treatment_soc_mean_T_ha = as.numeric(gsub(",", ".", treatment_soc_mean_T_ha)),
    control_soc_mean_T_ha   = as.numeric(gsub(",", ".", control_soc_mean_T_ha)),
    treatment_soc_sd_T_ha   = as.numeric(gsub(",", ".", treatment_soc_sd_T_ha)),
    control_soc_sd_T_ha     = as.numeric(gsub(",", ".", control_soc_sd_T_ha)),
    history_C               = tolower(history_C),
    history_T               = tolower(history_T)
  )

#######################--------------------------------------------------------
# PART VIII : Export cleaned dataset
#######################--------------------------------------------------------

write.csv(Data, RAW_DB, row.names = FALSE)

# R3: the Google Sheet is a living document. A dated, frozen copy is written so
# that the analyses can be reproduced from the exact version deposited on
# Dataverse. Deposit this file (not the Sheet) with the revised manuscript.
write.csv(Data, sub("\\.csv$", paste0("_snapshot_", format(Sys.Date(), "%Y%m%d"), ".csv"), RAW_DB),
          row.names = FALSE)
