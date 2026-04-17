options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(PRROC)
  library(pROC)
})

source(here::here("R", "config_defaults.R"))
source(here::here("R", "io_read_config.R"))
source(here::here("R", "io_write_outputs.R"))
source(here::here("R", "logging_utils.R"))
source(here::here("R", "warnings_utils.R"))
source(here::here("R", "stage_output.R"))

source(here::here("R", "fit_marginals", "preprocess_feature.R"))
source(here::here("R", "fit_marginals", "distributions.R"))

source(here::here("R", "evidence_transform.R"))
source(here::here("R", "distance_contrast.R"))
source(here::here("R", "assign_risk_bins.R"))
source(here::here("R", "summarize_risk_bins.R"))
source(here::here("R", "plot_distance_contrast.R"))

source(here::here("R", "group_positive_evidence.R"))
source(here::here("R", "intra_group_entropy.R"))
source(here::here("R", "inter_group_entropy.R"))
source(here::here("R", "group_percentiles.R"))
source(here::here("R", "group_patient_summary.R"))
source(here::here("R", "group_interpretation.R"))
source(here::here("R", "summarize_group_transport.R"))

source(here::here("analysis", "internal_validation", "baselines_and_metrics.R"))

source(here::here("analysis", "external_validation", "external_config.R"))
source(here::here("analysis", "external_validation", "external_validation.R"))

run_external_validation <- function(config_path) {
  cfg <- read_external_config(config_path)

  run_dir <- here::here(cfg$output$base_dir, cfg$output$run_name)
  logs_dir <- file.path(run_dir, 'logs')
  reports_dir <- file.path(run_dir, 'reports')
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)

  log_file <- file.path(logs_dir, 'run.log')
  warnings_file <- file.path(logs_dir, 'warnings.txt')
  init_log_file(log_file)
  init_warning_file(warnings_file)
  log_info('Starting external validation run', log_file)

  fitted <- load_fitted_artifacts(
    run_dir = here::here(cfg$run_lookup$run_dir),
    preferred_survival_model_path = cfg$run_lookup$preferred_survival_model_path
  )

  ext_paths <- resolve_external_paths(cfg)
  ext_all <- read_input_data(ext_paths$data_path)
  log_info(paste('External rows:', nrow(ext_all)), log_file)

  perf_all <- list(); surv_all <- list(); troc_all <- list()
  for (hid in cfg$hospitals$ids) {
    log_info(paste('Running hospital:', hid), log_file)
    hospital_df <- ext_all %>% dplyr::filter(.data[[cfg$external_data$hospital_col]] == hid)
    if (nrow(hospital_df) == 0) {
      log_warn(paste('No rows found for hospital', hid), log_file)
      next
    }
    hospital_dir <- file.path(run_dir, paste0('hospital_', hid))
    dir.create(hospital_dir, recursive = TRUE, showWarnings = FALSE)
    out <- run_external_hospital_validation(hospital_df, hid, fitted, cfg, hospital_dir)
    pm <- out$baselines$performance_metrics; pm$hospitalid <- hid; perf_all[[as.character(hid)]] <- pm
    if (!is.null(out$survival)) {
      if (nrow(out$survival$concordance) > 0) {
        cc <- out$survival$concordance; cc$hospitalid <- hid; surv_all[[as.character(hid)]] <- cc
      }
      if (nrow(out$survival$time_roc_long) > 0) {
        tt <- out$survival$time_roc_long; tt$hospitalid <- hid; troc_all[[as.character(hid)]] <- tt
      }
    }
  }

  if (length(perf_all) > 0) utils::write.csv(dplyr::bind_rows(perf_all), file.path(reports_dir, 'performance_metrics_all_hospitals.csv'), row.names = FALSE)
  if (length(surv_all) > 0) utils::write.csv(dplyr::bind_rows(surv_all), file.path(reports_dir, 'external_survival_concordance_all_hospitals.csv'), row.names = FALSE)
  if (length(troc_all) > 0) utils::write.csv(dplyr::bind_rows(troc_all), file.path(reports_dir, 'external_survival_timeROC_all_hospitals.csv'), row.names = FALSE)

  log_info('External validation run complete', log_file)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop('Usage: Rscript run_external_validation_patched.R /mnt/data/external_validation.yaml')
run_external_validation(args[[1]])