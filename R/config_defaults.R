`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

normalize_fit_config <- function(fit_cfg = NULL) {
  fit_cfg <- fit_cfg %||% list()

  fit_cfg$y_col <- fit_cfg$y_col %||% "label"
  fit_cfg$positive <- fit_cfg$positive %||% "1"
  fit_cfg$features <- fit_cfg$features %||% NULL

  fit_cfg$laplace <- as.numeric(fit_cfg$laplace %||% 1)
  fit_cfg$raw_min_sd <- as.numeric(fit_cfg$raw_min_sd %||% 1e-4)
  fit_cfg$gaussian_min_sd <- as.numeric(fit_cfg$gaussian_min_sd %||% 1e-4)
  fit_cfg$winsor_p <- as.numeric(fit_cfg$winsor_p %||% 0.01)
  fit_cfg$fraction_eps <- as.numeric(fit_cfg$fraction_eps %||% 1e-6)

  fit_cfg$selection_seed <- as.integer(fit_cfg$selection_seed %||% 1L)
  fit_cfg$selection_val_frac <- as.numeric(fit_cfg$selection_val_frac %||% 0.2)
  fit_cfg$trim_frac <- as.numeric(fit_cfg$trim_frac %||% 0.01)

  fit_cfg$min_n_total <- as.integer(fit_cfg$min_n_total %||% 25)
  fit_cfg$min_n_per_class <- as.integer(fit_cfg$min_n_per_class %||% 10)

  fit_cfg$numeric_candidates <- fit_cfg$numeric_candidates %||% c("gaussian", "lognormal", "gamma")
  fit_cfg$count_candidates <- fit_cfg$count_candidates %||% c("negbinom", "poisson")
  fit_cfg$fraction_candidates <- fit_cfg$fraction_candidates %||% c("logit_gaussian")

  fit_cfg$min_finite_logpdf_frac <- as.numeric(fit_cfg$min_finite_logpdf_frac %||% 0.80)
  fit_cfg$min_finite_llr_frac <- as.numeric(fit_cfg$min_finite_llr_frac %||% 0.80)
  fit_cfg$max_abs_llr_p99 <- as.numeric(fit_cfg$max_abs_llr_p99 %||% 50)

    # GMM2 controls
  fit_cfg$gmm2_max_iter <- as.integer(fit_cfg$gmm2_max_iter %||% 200L)
  fit_cfg$gmm2_tol <- as.numeric(fit_cfg$gmm2_tol %||% 1e-6)
  fit_cfg$gmm2_min_n <- as.integer(fit_cfg$gmm2_min_n %||% 20L)

  fit_cfg$gmm2_min_weight <- as.numeric(fit_cfg$gmm2_min_weight %||% 0.10)
  fit_cfg$gmm2_min_sd <- as.numeric(fit_cfg$gmm2_min_sd %||% 0.05)
  fit_cfg$gmm2_min_sep_over_pooled_sd <- as.numeric(fit_cfg$gmm2_min_sep_over_pooled_sd %||% 0.50)
  fit_cfg$gmm2_max_sd_ratio <- as.numeric(fit_cfg$gmm2_max_sd_ratio %||% 4)
  
  fit_cfg
}

normalize_transform_config <- function(transform_cfg = NULL) {
  transform_cfg <- transform_cfg %||% list()

  transform_cfg$drop_nonfinite_L_cols <- transform_cfg$drop_nonfinite_L_cols %||% TRUE
  transform_cfg$min_L_sd <- as.numeric(transform_cfg$min_L_sd %||% 1e-4)
  transform_cfg$cap_llr <- transform_cfg$cap_llr %||% TRUE
  transform_cfg$llr_cap <- as.numeric(transform_cfg$llr_cap %||% 50)

  transform_cfg
}

normalize_distance_config <- function(distance_cfg = NULL) {
  distance_cfg <- distance_cfg %||% list()

  distance_cfg$ridge <- as.numeric(distance_cfg$ridge %||% 1e-4)

  distance_cfg
}

normalize_risk_bin_config <- function(risk_cfg = NULL) {
  risk_cfg <- risk_cfg %||% list()

  risk_cfg$breaks <- risk_cfg$breaks %||% c(-Inf, 0, 2, Inf)
  risk_cfg$labels <- risk_cfg$labels %||% c("low", "ambiguous", "high")

  risk_cfg
}

normalize_plot_config <- function(plot_cfg = NULL) {
  plot_cfg <- plot_cfg %||% list()

  plot_cfg$width <- plot_cfg$width %||% 8
  plot_cfg$height <- plot_cfg$height %||% 5
  plot_cfg$dpi <- plot_cfg$dpi %||% 300

  plot_cfg
}

normalize_group_config <- function(group_cfg = NULL) {
  group_cfg <- group_cfg %||% list()

  group_cfg$enabled <- group_cfg$enabled %||% FALSE
  group_cfg$group_map_path <- group_cfg$group_map_path %||% NULL
  group_cfg$exclude_from_percentiles <- group_cfg$exclude_from_percentiles %||% c("anchor_age", "gender")
  group_cfg$softplus_eps <- as.numeric(group_cfg$softplus_eps %||% 1e-12)
  group_cfg$winsor_p <- as.numeric(group_cfg$winsor_p %||% 0.01)
  group_cfg$ecdf_log_eps <- as.numeric(group_cfg$ecdf_log_eps %||% 1e-6)
  group_cfg$atypical_threshold <- as.numeric(group_cfg$atypical_threshold %||% 0.75)
  group_cfg$top_n_atypical <- as.integer(group_cfg$top_n_atypical %||% 3)

  group_cfg
}