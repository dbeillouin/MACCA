MACCA

MACCA is an R-based project developed by Damien Beillouin for analyzing soil carbon dynamics, sequestration rates, and environmental predictors using machine learning and statistical models. This code is associated with the forthcoming article:

Beyond averages: depth, time, and context shape soil carbon outcomes in agroforestry in Latin America and the Caribbean
Damien Beillouin¹²*, Cloé Verstraete³⁴⁵, Rémi Cardinael⁴⁶⁷, Ulysse Chabroux²³⁸, Jean-Baptiste Laurent³⁴, Pierre-André Waite³⁴, Julien Demenois³⁴

¹ CIRAD, UPR HortSys, F-34398 Montpellier, France
² HortSys, Univ Montpellier, CIRAD, Montpellier, France
³ CIRAD, UPR AIDA, F-34398 Montpellier, France
⁴ AIDA, Univ Montpellier, CIRAD, Montpellier, France
⁵ AgroParisTech, PSL Université, Palaiseau, France
⁶ CIRAD, UPR AIDA, Harare, Zimbabwe
⁷ Department of Plant Production Sciences and Technologies, University of Zimbabwe, Harare, Zimbabwe
⁸ École Normale Supérieure, PSL Université, Laboratoire de Géologie, Paris, France

*Corresponding author: Damien Beillouin

Project Structure

Scripts/

1_LOAD_DATA.R # Load and clean datasets

2.RR.R # Compute response ratios

2.Seq_rate.R # Calculate soil organic carbon (SOC) accumulation rates

3.Models.R # Fit machine learning and statistical models

4.Importance_Plots.R # Visualize variable importance

4.Importance_plot_GB.R # Specific importance plots for Gradient Boosting

5.PDP_initialSOC_depthV2.R# Partial dependence plots for initial SOC by depth

5.PDP_precipitation.R # PDPs for precipitation effects

5_PDP_Temperature.R # PDPs for temperature effects

6.PDP_categorial.R # PDPs for categorical predictors

6_PLOT_Syntehsis_categorial.R # Synthesis plots for categorical predictors

Requirements

R (>= 4.0)

Packages: dplyr, tidyr, ggplot2, xgboost, mgcv, pdp, data.table (others as needed, see individual scripts)

Install missing packages with:

install.packages(c("dplyr", "tidyr", "ggplot2", "xgboost", "mgcv", "pdp", "data.table"))

Usage

1. Clone the repository:

   git clone https://github.com/dbeillouin/MACCA.git

   Open R or RStudio and set your working directory to the repository.

Run the scripts sequentially starting from 1_LOAD_DATA.R to load data, calculate response metrics, fit models, and generate plots. Some scripts depend on outputs from previous scripts.

Outputs

Partial dependence plots (PDPs) for continuous and categorical predictors.

Feature importance plots for Gradient Boosting models.

Summary tables of response ratios and SOC accumulation rates.

License

This project is under development. License information will be added later.

Author

Damien Beillouin
Contact: [beillouin@cirad.fr]
