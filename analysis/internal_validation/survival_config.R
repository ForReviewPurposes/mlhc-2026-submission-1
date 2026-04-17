`%||%` <- function(x, y) if (is.null(x)) y else x

read_survival_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop("Survival config file not found: ", config_path)
  }

  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read survival config files.")
  }

  cfg <- yaml::read_yaml(config_path)

  cfg$run_lookup <- cfg$run_lookup %||% list()
  cfg$data_sources <- cfg$data_sources %||% list()
  cfg$survival_data <- cfg$survival_data %||% list()
  cfg$benchmarks <- cfg$benchmarks %||% list()
  cfg$ddist_bins <- cfg$ddist_bins %||% list()
  cfg$time_roc <- cfg$time_roc %||% list()
  cfg$output <- cfg$output %||% list()

  cfg$run_lookup$run_dir <- cfg$run_lookup$run_dir %||% ""
  cfg$run_lookup$outputs_base_dir <- cfg$run_lookup$outputs_base_dir %||% "outputs"
  cfg$run_lookup$run_prefix <- cfg$run_lookup$run_prefix %||% "mimic_internal_validation_v1_"

  cfg$data_sources$train_data_path <- cfg$data_sources$train_data_path %||% ""
  cfg$data_sources$eval_data_path <- cfg$data_sources$eval_data_path %||% ""
  cfg$data_sources$eval_label <- cfg$data_sources$eval_label %||% "val"

  cfg$survival_data$time_to_event_path <- cfg$survival_data$time_to_event_path %||% ""
  cfg$survival_data$severity_scores_path <- cfg$survival_data$severity_scores_path %||% ""
  cfg$survival_data$duration_col <- cfg$survival_data$duration_col %||% "duration_hours_from_24h"
  cfg$survival_data$event_col <- cfg$survival_data$event_col %||% "event_after_24h"
  cfg$survival_data$age_col <- cfg$survival_data$age_col %||% "anchor_age"

  cfg$benchmarks$include_age <- isTRUE(cfg$benchmarks$include_age %||% TRUE)
  cfg$benchmarks$severity_score_cols <- cfg$benchmarks$severity_score_cols %||% c("apsiii", "lods", "oasis", "sapsii", "sirs", "sofa_24hours")
  cfg$benchmarks$baseline_score_cols <- cfg$benchmarks$baseline_score_cols %||% c("p_nb", "p_rf", "p_rf_llr")

  cfg$ddist_bins$breaks <- cfg$ddist_bins$breaks %||% NULL
  cfg$ddist_bins$labels <- cfg$ddist_bins$labels %||% NULL

  cfg$time_roc$hours <- as.numeric(cfg$time_roc$hours %||% c(48, 72, 96, 128))

  cfg$output$subfolder_name <- cfg$output$subfolder_name %||% "survival_analysis"
  cfg$output$preferred_km_method <- cfg$output$preferred_km_method %||% "d_dist_groups_age"

  cfg
}

resolve_ddist_bins <- function(surv_cfg, base_cfg) {
  breaks <- surv_cfg$ddist_bins$breaks
  labels <- surv_cfg$ddist_bins$labels

  if (!is.null(breaks) && !is.null(labels)) {
    return(list(
      breaks = as.numeric(unlist(breaks)),
      labels = as.character(unlist(labels))
    ))
  }

  if (!is.null(base_cfg$risk_bins) &&
      !is.null(base_cfg$risk_bins$breaks) &&
      !is.null(base_cfg$risk_bins$labels)) {
    return(list(
      breaks = as.numeric(unlist(base_cfg$risk_bins$breaks)),
      labels = as.character(unlist(base_cfg$risk_bins$labels))
    ))
  }

  list(
    breaks = c(-Inf, -0.1, 2, Inf),
    labels = c("Low", "Ambiguous", "High")
  )
}