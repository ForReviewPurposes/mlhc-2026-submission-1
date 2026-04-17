read_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }

  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read config files.")
  }

  cfg <- yaml::read_yaml(config_path)
  cfg <- normalize_config(cfg)
  cfg
}

normalize_config <- function(cfg) {
  required_sections <- c("run", "data", "output")
  missing_sections <- setdiff(required_sections, names(cfg))
  if (length(missing_sections) > 0) {
    stop("Missing required config sections: ", paste(missing_sections, collapse = ", "))
  }

  cfg$run <- cfg$run %||% list()
  cfg$data <- cfg$data %||% list()
  cfg$output <- cfg$output %||% list()

  cfg$run$mode <- cfg$run$mode %||% "fit_and_score"
  cfg$run$seed <- cfg$run$seed %||% 1L
  cfg$run$name <- cfg$run$name %||% "internal_validation"

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

  cfg$fit <- normalize_fit_config(cfg$fit)
  cfg$transform <- normalize_transform_config(cfg$transform)
  cfg$distance <- normalize_distance_config(cfg$distance)
  cfg$risk_bins <- normalize_risk_bin_config(cfg$risk_bins)
  cfg$plots <- normalize_plot_config(cfg$plots)
  cfg$groups <- normalize_group_config(cfg$groups)

  cfg
}