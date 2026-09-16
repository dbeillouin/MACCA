# MACCA – Soil organic carbon under agroforestry in the Neotropics (meta-analysis)

Code accompanying:

> Beillouin D, Verstraete C, Cardinael R, Chabroux U, Laurent J-B, Waite P-A, Demenois J. *Baseline soil carbon and temporal drivers, not management, govern soil carbon sequestration in tropical agroforestry. A meta-analysis.* Agronomy for Sustainable Development (under review).

The workflow reproduces every number, table and figure of the article and of its Electronic Supplementary Material from a frozen version of the MACCA database.

## Data

The database is not stored in this repository. Download it from the CIRAD Dataverse (https://doi.org/10.18167/DVN1/GISJUZ) and place it in `data/raw/`:

| File | Content |
|---|---|
| `data/raw/MACCA_database_R3_corrected.csv` | database used for the analyses (corrected version, Dataverse) |
| `data/data_corrections.csv` | documented corrections applied to the previous version of the database (status, reason, source) |
| `data/prisma_counts.csv` | counts of the literature screening (Fig. S1) |

The uncorrected version of the database remains available as an earlier version of the Dataverse record. If `data/raw/Data_for_analysis_R2_20260430.csv` is used instead of the corrected file, the corrections listed in `data/data_corrections.csv` are applied automatically, and both routes give identical results.

Studies excluded at the study level are declared in `00_run_all.R` (`EXCLUDED_STUDIES`).

## Structure

```
00_run_all.R            runs the whole workflow (paths are relative to this file)
R/utils.R               shared functions
R/01_load_data.R        optional: rebuilds the database from the Excel export and WorldClim
R/02_prepare_datasets.R corrections, study-level exclusions, SD imputation, effect sizes
R/03_meta_analysis.R    three-level meta-analyses, heterogeneity, small-study effects, location-scale models
R/04_ml_models.R        MetaForest outlier screening, XGBoost, cross-validation (random, profile, study), importance, SHAP interactions
R/05_partial_dependence.R  partial dependence by depth class, restricted to observed ranges, study bootstrap
R/06_categorical_effects.R standardized XGBoost predictions by category (Table S9)
R/07_data_support.R     data coverage by depth and time since conversion
R/08_moderator_tests.R  meta-regressions on observed effect sizes, heterogeneity explained
R/09_robustness.R       leave-one-study-out, thresholds, confounding, sensitivity analyses
R/10_figures.R          Figs. 2–4 and supplementary figures
R/11_prisma.R           PRISMA flow diagram
```

## Running

```bash
Rscript 00_run_all.R quick         # quick test (few bootstrap resamples)
Rscript 00_run_all.R               # full run (about 20 min)
Rscript 00_run_all.R only10 only11 # selected steps only
Rscript 00_run_all.R with01        # rebuild the database (needs WorldClim rasters next to the project folder)
```

Required R packages: dplyr, tidyr, readr, readxl, forcats, stringr, purrr, metafor, metaforest, xgboost, splines, ggplot2, patchwork, maps, mapdata. Package versions used for the article are recorded in `outputs/sessionInfo.txt`.

## Outputs (`outputs/`)

* `numbers_for_manuscript.csv` – every number quoted in the article, with the section where it is used
* `data/` – analysis datasets, corrected database, excluded studies, outliers, corrections log
* `tables/` – meta-analysis, moderator tests, robustness, cross-validation, importance, interactions, Table 2
* `figures/` – figures as PDF, TIFF (600 dpi) and PNG

## Citation

See `CITATION.cff`. Please cite the article and the dataset (https://doi.org/10.18167/DVN1/GISJUZ).
