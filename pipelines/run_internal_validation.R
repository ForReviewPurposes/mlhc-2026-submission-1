options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
})

source(here::here("R", "config_defaults.R"))
source(here::here("R", "io_read_config.R"))
source(here::here("R", "io_paths.R"))
source(here::here("R", "io_write_outputs.R"))
source(here::here("R", "logging_utils.R"))
source(here::here("R", "warnings_utils.R"))
source(here::here("R", "stage_output.R"))

# fit_marginals
source(here::here("R", "fit_marginals", "prepare_feature_registry.R"))
source(here::here("R", "fit_marginals", "make_selection_split.R"))
source(here::here("R", "fit_marginals", "preprocess_feature.R"))
source(here::here("R", "fit_marginals", "distributions.R"))
source(here::here("R", "fit_marginals", "fit_family_models.R"))
source(here::here("R", "fit_marginals", "metrics_marginal_fit.R"))
source(here::here("R", "fit_marginals", "evaluate_candidate_families.R"))
source(here::here("R", "fit_marginals", "assess_candidate_fit.R"))
source(here::here("R", "fit_marginals", "fallback_ladder.R"))
source(here::here("R", "fit_marginals", "fit_one_feature_marginal.R"))
source(here::here("R", "fit_marginals", "fit_report.R"))
source(here::here("R", "fit_marginals", "fit_marginals.R"))

# scoring
source(here::here("R", "evidence_transform.R"))
source(here::here("R", "distance_contrast.R"))
source(here::here("R", "assign_risk_bins.R"))
source(here::here("R", "summarize_risk_bins.R"))
source(here::here("R", "plot_distance_contrast.R"))

# feature group scores
source(here::here("R", "group_positive_evidence.R"))
source(here::here("R", "inter_group_entropy.R"))
source(here::here("R", "intra_group_entropy.R"))
source(here::here("R", "group_ecdfs.R"))
source(here::here("R", "group_percentiles.R"))
source(here::here("R", "group_patient_summary.R"))
source(here::here("R", "group_interpretation.R"))
source(here::here("R", "group_pipeline_summary.R"))
source(here::here("R", "summarize_group_transport.R"))

# baseline comparisons
source(here::here("analysis", "internal_validation", "baselines_and_metrics.R"))

run_internal_validation <- function(config_path) {
  cfg <- read_config(here::here(config_path))
  paths <- init_run_paths(cfg, run_type = "internal_validation")
  
  init_log_file(paths$log_file)
  init_warning_file(paths$warnings_file)
  write_config_snapshot(cfg, paths$config_snapshot)
  
  log_info("Starting internal validation run", paths$log_file)
  log_info(paste("Reading config:", config_path), paths$log_file)
  
  train_df <- read_input_data(cfg$data$train_data_path)
  val_df   <- read_input_data(cfg$data$val_data_path)
  
  log_info(paste("Train rows:", nrow(train_df)), paths$log_file)
  log_info(paste("Val rows:", nrow(val_df)), paths$log_file)
  
  # ---- stage 1: fit marginals ----
  log_info("Stage 1: fit_marginals started", paths$log_file)

  s1 <- fit_marginals(
    train_df = train_df,
    fit_cfg = cfg$fit,
    exclude_features = c(cfg$data$id_col)
  )

  write_stage_outputs("fit_marginals", s1, paths)
  saveRDS(s1$result, paths$marginal_model_rds)

  log_info("Stage 1: fit_marginals completed", paths$log_file)
  
  # ---- stage 2: transform train and val ----
  log_info("Stage 2: evidence_transform train started", paths$log_file)
  s2_train <- evidence_transform(
    new_df = train_df,
    marginal_fit = s1$result,
    transform_cfg = cfg$transform
  )
  
  log_info("Stage 2: evidence_transform val started", paths$log_file)
  s2_val <- evidence_transform(
    new_df = val_df,
    marginal_fit = s1$result,
    transform_cfg = cfg$transform
  )
  
  write_stage_outputs("evidence_transform_train", s2_train, paths)
  write_stage_outputs("evidence_transform_val", s2_val, paths)
  log_info("Stage 2: evidence_transform completed", paths$log_file)
  
  # ---- align transformed feature spaces ----
  common_features <- intersect(colnames(s2_train$result$L), colnames(s2_val$result$L))
  if (length(common_features) == 0) {
    stop("No common transformed features remain between train and val.")
  }
  
  L_train <- s2_train$result$L[, common_features, drop = FALSE]
  L_val   <- s2_val$result$L[, common_features, drop = FALSE]
  
  sum_llr_val <- rowSums(L_val)
  
  # ---- stage 3: geometry fit + score ----
  log_info("Stage 3: distance_contrast_fit started", paths$log_file)
  geom_fit <- distance_contrast_fit(
    train_L = L_train,
    train_y = train_df[[cfg$fit$y_col]],
    positive_label = cfg$fit$positive,
    ridge = cfg$distance$ridge
  )
  saveRDS(geom_fit, paths$geometry_model_rds)
  
  log_info("Stage 3: distance_contrast_score val started", paths$log_file)
  s3_val <- distance_contrast_score(
    L = L_val,
    geom_fit = geom_fit
  )
  write_stage_outputs("distance_contrast_val", s3_val, paths)
  log_info("Stage 3: distance_contrast completed", paths$log_file)
  
  # ---- build score table ----
  score_tbl <- data.frame(
    stay_id = val_df[[cfg$data$id_col]],
    sum_llr = sum_llr_val,
    d_dist = s3_val$result$d_dist,
    stringsAsFactors = FALSE
  )
  utils::write.csv(score_tbl, paths$score_table_csv, row.names = FALSE)

  score_tbl_train <- data.frame(
    stay_id = train_df[[cfg$data$id_col]],
    sum_llr = rowSums(L_train),
    d_dist = distance_contrast_score(
      L = L_train,
      geom_fit = geom_fit
    )$result$d_dist,
    stringsAsFactors = FALSE
  )
  
  # ---- stage 4: assign primary risk bins from sum_llr quartiles ----
  risk_tbl <- assign_risk_bins(
    stay_id = val_df[[cfg$data$id_col]],
    score = score_tbl$sum_llr,
    score_name = "sum_llr",
    use_quantiles = TRUE,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("Q1", "Q2", "Q3", "Q4")
  )
  utils::write.csv(risk_tbl, paths$risk_table_csv, row.names = FALSE)

  # optional secondary d_dist bins for comparison
  risk_tbl_ddist <- assign_risk_bins(
    stay_id = val_df[[cfg$data$id_col]],
    score = score_tbl$d_dist,
    score_name = "d_dist",
    bin_spec = cfg$risk_bins,
    use_quantiles = FALSE
  )
  utils::write.csv(
    risk_tbl_ddist,
    file.path(paths$reports_dir, "risk_table_ddist.csv"),
    row.names = FALSE
  )

  # ---- stage 5: summarize bins ----
  outcome <- NULL
  if (!is.null(cfg$data$outcome_col) && cfg$data$outcome_col %in% names(val_df)) {
    outcome <- val_df[[cfg$data$outcome_col]]
  }

  risk_summary <- summarize_risk_bins(
    risk_tbl = risk_tbl,
    outcome = outcome,
    score_col = "sum_llr"
  )
  utils::write.csv(risk_summary, paths$risk_summary_csv, row.names = FALSE)

  risk_summary_ddist <- summarize_risk_bins(
    risk_tbl = risk_tbl_ddist,
    outcome = outcome,
    score_col = "d_dist"
  )
  utils::write.csv(
    risk_summary_ddist,
    file.path(paths$reports_dir, "risk_summary_ddist.csv"),
    row.names = FALSE
  )
  
  # ---- stage 6: plots ----
  plot_input <- data.frame(
    sum_llr = score_tbl$sum_llr,
    d_dist = score_tbl$d_dist
  )
  
  plot_list <- plot_distance_contrast(
    score_df = plot_input,
    outcome = outcome,
    plots_cfg = cfg$plots
  )
  
  ggplot2::ggsave(
    filename = paths$sum_llr_density_png,
    plot = plot_list$sum_llr_plot,
    width = cfg$plots$width,
    height = cfg$plots$height,
    dpi = cfg$plots$dpi
  )
  
  ggplot2::ggsave(
    filename = paths$d_dist_density_png,
    plot = plot_list$d_dist_plot,
    width = cfg$plots$width,
    height = cfg$plots$height,
    dpi = cfg$plots$dpi
  )

  # ---- stage 7: group metrics ----

  group_outputs <- NULL
  train_patient_summary <- NULL
  patient_summary <- NULL
  s8 <- NULL

  if (isTRUE(cfg$groups$enabled) &&
      !is.null(cfg$groups$group_map_path) &&
      nzchar(cfg$groups$group_map_path)) {
    
    log_info("Stage 7: group metrics started", paths$log_file)
    
    ddist_bin_cfg <- list(
      breaks = cfg$risk_bins$breaks,
      labels = cfg$risk_bins$labels
    )
    
    group_summary_obj <- build_group_patient_summaries(
      train_L = L_train,
      eval_L = L_val,
      train_score_tbl = score_tbl_train,
      eval_score_tbl = score_tbl,
      train_ids = train_df[[cfg$data$id_col]],
      eval_ids = val_df[[cfg$data$id_col]],
      ddist_bin_cfg = ddist_bin_cfg,
      group_map_path = cfg$groups$group_map_path,
      groups_cfg = cfg$groups
    )
    
    train_patient_summary <- group_summary_obj$train_patient_summary
    patient_summary <- group_summary_obj$eval_patient_summary
    
    # write main patient summaries
    utils::write.csv(
      train_patient_summary,
      file = file.path(paths$scores_dir, "group_patient_summary_train.csv"),
      row.names = FALSE
    )
    
    utils::write.csv(
      patient_summary,
      file = file.path(paths$scores_dir, "group_patient_summary_val.csv"),
      row.names = FALSE
    )
    
    # detailed intra-group entropy reports
    utils::write.csv(
      group_summary_obj$train_objects$intra_entropy_detail,
      file = file.path(paths$reports_dir, "intra_group_entropy_detail_train.csv"),
      row.names = FALSE
    )
    
    utils::write.csv(
      group_summary_obj$eval_objects$intra_entropy_detail,
      file = file.path(paths$reports_dir, "intra_group_entropy_detail_val.csv"),
      row.names = FALSE
    )

    train_transport <- summarize_group_percentile_transport(
      patient_summary = train_patient_summary,
      atypical_threshold = cfg$groups$atypical_threshold,
      primary_bin_col = "risk_bin_primary"
    )

    val_transport <- summarize_group_percentile_transport(
      patient_summary = patient_summary,
      atypical_threshold = cfg$groups$atypical_threshold,
      primary_bin_col = "risk_bin_primary"
    )

    utils::write.csv(
      train_transport$group_atypical_frequency,
      file = file.path(paths$reports_dir, "group_atypical_frequency_train.csv"),
      row.names = FALSE
    )

    utils::write.csv(
      val_transport$group_atypical_frequency,
      file = file.path(paths$reports_dir, "group_atypical_frequency_val.csv"),
      row.names = FALSE
    )

    utils::write.csv(
      train_transport$patient_count_summary,
      file = file.path(paths$reports_dir, "atypical_group_count_summary_train.csv"),
      row.names = FALSE
    )

    utils::write.csv(
      val_transport$patient_count_summary,
      file = file.path(paths$reports_dir, "atypical_group_count_summary_val.csv"),
      row.names = FALSE
    )

    if (!is.null(val_transport$risk_bin_summary_primary)) {
      utils::write.csv(
        val_transport$risk_bin_summary_primary,
        file = file.path(paths$reports_dir, "risk_bin_summary_primary_val.csv"),
        row.names = FALSE
      )
    }

    train_possum_summary <- summarize_group_positive_sum_distribution(
      group_pos_sum = group_summary_obj$train_objects$group_pos_sum,
      risk_bin_primary = train_patient_summary$risk_bin_primary
    )

    val_possum_summary <- summarize_group_positive_sum_distribution(
      group_pos_sum = group_summary_obj$eval_objects$group_pos_sum,
      risk_bin_primary = patient_summary$risk_bin_primary
    )

    utils::write.csv(
      train_possum_summary$group_pos_sum_summary,
      file = file.path(paths$reports_dir, "group_pos_sum_summary_train.csv"),
      row.names = FALSE
    )

    utils::write.csv(
      val_possum_summary$group_pos_sum_summary,
      file = file.path(paths$reports_dir, "group_pos_sum_summary_val.csv"),
      row.names = FALSE
    )

    if (!is.null(train_possum_summary$group_pos_sum_summary_by_bin)) {
      utils::write.csv(
        train_possum_summary$group_pos_sum_summary_by_bin,
        file = file.path(paths$reports_dir, "group_pos_sum_summary_by_bin_train.csv"),
        row.names = FALSE
      )
    }

    if (!is.null(val_possum_summary$group_pos_sum_summary_by_bin)) {
      utils::write.csv(
        val_possum_summary$group_pos_sum_summary_by_bin,
        file = file.path(paths$reports_dir, "group_pos_sum_summary_by_bin_val.csv"),
        row.names = FALSE
      )
    }

    plot_group_ecdfs(
      group_pos_df = group_summary_obj$train_objects$group_pos_sum,
      out_dir = file.path(paths$plots_dir, "ecdf_groups")
    )

    ggplot2::ggsave(
      filename = file.path(paths$plots_dir, "ecdf_groups", "all_ecdfs.png"),

      plot = 
      (ggplot(group_summary_obj$train_objects$group_pos_sum %>%
    pivot_longer(
      cols = everything(),
      names_to = "group",
      values_to = "value"
    ), aes(x = value, color = group)) +
      stat_ecdf(alpha = 0.5) +
      plot_theme_mlhc()),

      width = cfg$plots$width,
      height = cfg$plots$height,
      dpi = cfg$plots$dpi
    )

    compute_ecdf_summary(
      group_summary_obj$train_objects$group_pos_sum,
      file.path(paths$reports_dir, "ecdf_summary.csv")
    )
    
    group_outputs <- list(
      train = group_summary_obj$train_objects,
      val = group_summary_obj$eval_objects,
      ecdf = group_summary_obj$ecdf_object
    )
    
    log_info("Stage 7: group metrics completed", paths$log_file)
  } else {
    train_patient_summary <- NULL
  }

  # ---- stage 8: baselines and metrics ----
  log_info("Stage 8: baselines and metrics started", paths$log_file)

  s8 <- run_stage8_baselines(
    train_df = train_df,
    val_df = val_df,
    train_score_table = score_tbl_train,
    val_score_table = score_tbl,
    cfg = cfg,
    L_train = L_train,
    L_val = L_val
  )

  write_stage_outputs("baselines_and_metrics", s8, paths)

  utils::write.csv(
    s8$result$evaluation_table_train,
    file = file.path(paths$scores_dir, "evaluation_table_train.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    s8$result$evaluation_table_val,
    file = file.path(paths$scores_dir, "evaluation_table_val.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    s8$result$performance_metrics,
    file = file.path(paths$reports_dir, "performance_metrics.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    s8$result$mortality_fraction_curves,
    file = file.path(paths$reports_dir, "mortality_fraction_curves.csv"),
    row.names = FALSE
  )

  ggplot2::ggsave(
    filename = file.path(paths$plots_dir, "mortality_fraction_curve_d_dist.png"),
    plot = s8$result$mortality_fraction_curve_plot$d_dist,
    width = cfg$plots$width,
    height = cfg$plots$height,
    dpi = cfg$plots$dpi
  )

  ggplot2::ggsave(
    filename = file.path(paths$plots_dir, "mortality_fraction_curve_sum_llr.png"),
    plot = s8$result$mortality_fraction_curve_plot$sum_llr,
    width = cfg$plots$width,
    height = cfg$plots$height,
    dpi = cfg$plots$dpi
  )

  ggplot2::ggsave(
    filename = file.path(paths$plots_dir, "mortality_fraction_curve_p_nb.png"),
    plot = s8$result$mortality_fraction_curve_plot$p_nb,
    width = cfg$plots$width,
    height = cfg$plots$height,
    dpi = cfg$plots$dpi
  )

  ggplot2::ggsave(
    filename = file.path(paths$plots_dir, "mortality_fraction_curve_p_rf.png"),
    plot = s8$result$mortality_fraction_curve_plot$p_rf,
    width = cfg$plots$width,
    height = cfg$plots$height,
    dpi = cfg$plots$dpi
  )

  if (!is.null(s8$result$mortality_fraction_curve_plot$p_rf_llr)) {
  ggplot2::ggsave(
    filename = file.path(paths$plots_dir, "mortality_fraction_curve_p_rf_llr.png"),
    plot = s8$result$mortality_fraction_curve_plot$p_rf_llr,
    width = cfg$plots$width,
    height = cfg$plots$height,
    dpi = cfg$plots$dpi
  )
}

if (!is.null(s8$result$density_probability_plots$p_rf)) {
  ggplot2::ggsave(
    filename = file.path(paths$plots_dir, "density_p_rf.png"),
    plot = s8$result$density_probability_plots$p_rf,
    width = cfg$plots$width,
    height = cfg$plots$height,
    dpi = cfg$plots$dpi
  )
}

if (!is.null(s8$result$density_probability_plots$p_rf_llr)) {
  ggplot2::ggsave(
    filename = file.path(paths$plots_dir, "density_p_rf_llr.png"),
    plot = s8$result$density_probability_plots$p_rf_llr,
    width = cfg$plots$width,
    height = cfg$plots$height,
    dpi = cfg$plots$dpi
  )
}

  log_info("Stage 8: baselines and metrics completed", paths$log_file)

  # ---- aggregate warnings / dropped ----
  all_warnings <- bind_stage_warnings(
    s1$warnings,
    s2_train$warnings,
    s2_val$warnings,
    s3_val$warnings,
    if (!is.null(s8)) s8$warnings else NULL
  )

  all_dropped <- bind_stage_dropped(
    s1$dropped,
    s2_train$dropped,
    s2_val$dropped,
    s3_val$dropped,
    if (!is.null(s8)) s8$dropped else NULL
  )

  write_all_warnings(all_warnings, paths$warnings_file, paths$warnings_csv)
  write_all_dropped(all_dropped, paths$dropped_file)
  
  
  final_output <- list(
    config = cfg,
    paths = paths,
    marginal_fit = s1$result,
    geometry_fit = geom_fit,
    score_table_train = score_tbl_train,
    score_table_val = score_tbl,
    risk_table = risk_tbl,
    risk_summary = risk_summary,
    group_outputs = group_outputs,
    patient_summary_train = train_patient_summary,
    patient_summary_val = patient_summary,
    rf_obj = if (!is.null(s8)) s8$result$rf_obj else NULL,
    nb_obj = if (!is.null(s8)) s8$result$nb_obj else NULL,
    rf_llr_obj = if (!is.null(s8)) s8$result$rf_llr_obj else NULL,
    evaluation_train = if (!is.null(s8)) s8$result$evaluation_table_train else NULL,
    evaluation_val = if (!is.null(s8)) s8$result$evaluation_table_val else NULL,
    performance_metrics = if (!is.null(s8)) s8$result$performance_metrics else NULL,
    mortality_fraction_curves = if (!is.null(s8)) s8$result$mortality_fraction_curves else NULL
  )
  
  saveRDS(final_output, paths$run_rds)
  log_info("Internal validation run complete", paths$log_file)
  
  invisible(final_output)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript pipelines/run_internal_validation.R configs/internal_validation.yaml")
}

run_internal_validation(args[[1]])