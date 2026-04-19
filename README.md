# Interpretable ICU Risk Characterization in Machine Learning Workflows via Population-Relative Evidence Representation

This repository contains the code used to construct a population-referenced evidence representation from MIMIC-IV, perform internal validation on MIMIC-IV, and evaluate external transfer on eICU without retraining.

Because MIMIC-IV and eICU are controlled-access datasets, raw source data are **not** distributed in this repository. Instead, this repository provides:

- SQL extraction scripts
- data preparation scripts
- analysis pipelines
- YAML configuration files
- representative analysis code and outputs

## Data Sources

This project uses:

- **MIMIC-IV**
- **MIMIC-IV derived**
- **eICU Collaborative Research Database**

Access to these datasets must be obtained independently through PhysioNet.

---

## Required Packages

here, dplyr, tidyr, purrr, stringr, ggplot2, caret, MASS, ranger, naivebayes, pROC, PRROC, survival, survminer, yaml

install.packages(c(
  "here", "dplyr", "tidyr", "purrr", "stringr", "ggplot2",
  "caret", "MASS", "ranger", "naivebayes",
  "pROC", "PRROC", "survival", "survminer", "yaml"
))

---

## Repository Workflow

The full workflow consists of:

1. Extracting MIMIC-IV cohort, feature, GCS, and ICU score tables
2. Preparing processed MIMIC-IV analysis files
3. Extracting eICU aggregated signal features
4. Preparing processed eICU analysis files
5. Running internal validation
6. Updating configuration files to point to the latest internal validation output directory
7. Running survival analysis
8. Running external validation

---

## Data Extraction

### MIMIC-IV
1. data/mimiciv/sql > run SQL files 01 to 06 in sequence in BigQuery.
2. Save outputs from the following SQL queries with the following exact file names in data/mimiciv/ :
- Query 01: "mimiciv_ext_icu_cohort_1_timetodeath.csv"
- Query 03: "mimiciv_ext_24h_features_wide_3.csv"
- Query 05: "mimiciv_ext_24h_gcs_features_wide_2.csv"
- Query 06: "mimiciv_ext_icu_scores_2.csv"
3. From repo root folder, run `Rscript data/mimiciv/run_prepare_mimic_data.R`

### eICU
1. data/eicu/sql > run SQL files 01 to 04 in sequence in BigQuery.
2. Save outputs from the following SQL queries with the following exact file names in data/eicu/ :
- Query 04: "eicu_signal_24h_agg_wide.csv".
3. From repo root folder, run `Rscript data/eicu/run_prepare_eicu_data.R`

## Analytics Execution
1. From repo root folder, run `Rscript pipelines/run_internal_validation.R configs/internal_validation.yaml` (execution time approx. 3-4 mins).
2. Open configs/survival_analysis.yaml > update `run_dir` with the latest internal validation output folder name.
3. Open configs/external_validation.yaml > update `run_dir` with the latest internal validation output folder name.
4. From repo root folder, run `Rscript analysis/internal_validation/run_survival_analysis.R configs/survival_analysis.yaml`
5. From repo root folder, run `Rscript pipelines/run_external_validation.R configs/external_validation.yaml`

---

## Outputs

The `outputs/` directory contains representative results from internal validation, survival analysis, and external validation pipelines. These include:
- Density plots
- Mortality curves
- Performance metrics (AUROC/AUPRC)
- sum\_llr risk quantile tables
- Feature-group ECDF summary tables (atypical group counts and group\_pos\_sum distribution summaries)
- Kaplan-Meier plots and survival analysis performance metrics.

Full outputs can be regenerated using the provided workflow.
