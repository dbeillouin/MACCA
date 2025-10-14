# MACCA: Mapping Agroforestry Carbon Change in the Americas

![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)
![License](https://img.shields.io/badge/license-GPL--3.0-lightgrey)
![Status](https://img.shields.io/badge/status-Development-orange)
![Software Heritage](https://img.shields.io/badge/Software%20Heritage-SWH-lightgrey)

**MACCA** is an R-based project developed by Damien Beillouin for analyzing soil carbon dynamics, sequestration rates, and environmental predictors using machine learning and statistical models. This repository accompanies the forthcoming article:

**Beyond averages: depth, time, and context shape soil carbon outcomes in agroforestry in Latin America and the Caribbean**  
*Damien Beillouin¹²*, Cloé Verstraete³⁴⁵, Rémi Cardinael⁴⁶⁷, Ulysse Chabroux²³⁸, Jean-Baptiste Laurent³⁴, Pierre-André Waite³⁴, Julien Demenois³⁴

---

## Project Overview

MACCA provides tools to:

- Compute **response ratios (RR)** and **SOC accumulation rates**.  
- Fit **machine learning models** (Random Forest, Gradient Boosting) and classical statistical models.  
- Generate **partial dependence plots (PDPs)** and **variable importance plots**.  
- Integrate **continuous and categorical environmental and management predictors**.  
- Synthesize **region-specific SOC dynamics** in Latin America and the Caribbean.

---

## Repository Structure

```
Scripts/
├── 1_LOAD_DATA.R
├── 2.RR.R
├── 2.Seq_rate.R
├── 3.Models.R
├── 4.Importance_Plots.R
├── 4.Importance_plot_GB.R
├── 5.PDP_initialSOC_depthV2.R
├── 5.PDP_precipitation.R
├── 5_PDP_Temperature.R
├── 6.PDP_categorial.R
└── 6_PLOT_Synthesis_categorial.R
```

## Requirements

- **R >= 4.0**  
- Packages: `dplyr`, `tidyr`, `ggplot2`, `xgboost`, `mgcv`, `pdp`, `data.table`  

Install missing packages with:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "xgboost", "mgcv", "pdp", "data.table"))
```

## Usage

1. Clone the repository:

git clone https://github.com/dbeillouin/MACCA.git


2. Open R or RStudio and set the working directory to the repository.

3. Run scripts sequentially from 1_LOAD_DATA.R onward. Outputs from previous scripts are required for downstream scripts.


## Outputs

Partial Dependence Plots (PDPs) for continuous and categorical predictors.

Feature importance plots for Gradient Boosting and Random Forest models.

Summary tables of response ratios and SOC accumulation rates.


## Data

The analysis relies on the dataset hosted on CIRAD Dataverse:

DOI: : 10.18167/DVN1/GISJUZ


Data Management & FAIR Compliance

Fully reproducible workflow: scripts automatically generate all intermediate outputs.

Metadata included: description of dataset, variables, and analysis pipeline.

Structured to maximize FAIR principles (Findable, Accessible, Interoperable, Reusable).

Compatible with Software Heritage harvesting for long-term archiving.

Versioned releases and issues enabled for reproducibility and community contributions.

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0). See the LICENSE file for details.


## Citation

Beillouin D., Verstraete C., Cardinael R., Chabroux U., Laurent J.-B., Waite P.-A., Demenois J. (2025). MACCA: Mapping Agroforestry Carbon Change in the Americas. GitHub repository: https://github.com/dbeillouin/MACCA