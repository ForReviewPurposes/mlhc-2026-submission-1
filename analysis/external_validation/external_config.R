`%||%` <- function(x, y) if (is.null(x)) y else x

read_external_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop("External config file not found: ", config_path)
  }

  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read external config files.")
  }

  cfg <- yaml::read_yaml(config_path)

  # Reuse the same shared normalization as internal validation
  cfg <- normalize_external_config(cfg)

  cfg
}

normalize_external_config <- function(cfg) {
  # Reuse internal-style required sections and defaults
  required_sections <- c("run", "data", "output")
  missing_sections <- setdiff(required_sections, names(cfg))
  if (length(missing_sections) > 0) {
    stop("Missing required config sections: ", paste(missing_sections, collapse = ", "))
  }

  cfg$run <- cfg$run %||% list()
  cfg$data <- cfg$data %||% list()
  cfg$output <- cfg$output %||% list()

  cfg$run$mode <- cfg$run$mode %||% "external_validation"
  cfg$run$seed <- cfg$run$seed %||% 1L
  cfg$run$name <- cfg$run$name %||% "external_validation"

  cfg$output$base_dir <- cfg$output$base_dir %||% "outputs"
  cfg$output$append_datetime <- cfg$output$append_datetime %||% TRUE

  base_run_name <- cfg$output$run_name %||% cfg$run$name

  if (isTRUE(cfg$output$append_datetime)) {
    cfg$output$run_name <- paste0(
      base_run_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S")
    )
  } else {
    cfg$output$run_name <- base_run_name
  }

  # Shared model/pipeline defaults from config_defaults.R
  cfg$fit <- normalize_fit_config(cfg$fit)
  cfg$transform <- normalize_transform_config(cfg$transform)
  cfg$distance <- normalize_distance_config(cfg$distance)
  cfg$risk_bins <- normalize_risk_bin_config(cfg$risk_bins)
  cfg$plots <- normalize_plot_config(cfg$plots)
  cfg$groups <- normalize_group_config(cfg$groups)

  # External-only sections
  cfg$run_lookup <- cfg$run_lookup %||% list()
  cfg$external_data <- cfg$external_data %||% list()
  cfg$hospitals <- cfg$hospitals %||% list()
  cfg$survival <- cfg$survival %||% list()

  cfg$run_lookup$run_dir <- cfg$run_lookup$run_dir %||% ""
  cfg$run_lookup$preferred_survival_model_path <-
    cfg$run_lookup$preferred_survival_model_path %||% ""

  # External data falls back to shared data/fit config where possible
  cfg$external_data$data_path <-
    cfg$external_data$data_path %||% cfg$data$val_data_path %||% ""
  cfg$external_data$id_col <-
    cfg$external_data$id_col %||% cfg$data$id_col %||% "stay_id"
  cfg$external_data$hospital_col <-
    cfg$external_data$hospital_col %||% "hospitalid"
  cfg$external_data$gender_col <-
    cfg$external_data$gender_col %||% "gender"
  cfg$external_data$mortality_col <-
    cfg$external_data$mortality_col %||% cfg$data$outcome_col %||% cfg$fit$y_col %||% "mortality"
  cfg$external_data$duration_col <-
    cfg$external_data$duration_col %||% "duration_hours_from_24h"
  cfg$external_data$event_col <-
    cfg$external_data$event_col %||% "event_after_24h"

  cfg$hospitals$ids <- as.integer(cfg$hospitals$ids %||% integer())

  cfg$survival$enabled <- isTRUE(cfg$survival$enabled %||% TRUE)
  cfg$survival$time_roc_hours <- as.numeric(cfg$survival$time_roc_hours %||% c(48, 72, 96, 128))
  cfg$survival$preferred_km_method <- cfg$survival$preferred_km_method %||% "d_dist_groups_age"

  cfg
}