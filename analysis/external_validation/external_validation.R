suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(timeROC)
  library(PRROC)
  library(pROC)
})

normalize_gender <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x %in% c("Male", "M", "male", "m") ~ "M",
    x %in% c("Female", "F", "female", "f") ~ "F",
    TRUE ~ x
  )
}

align_to_feature_set <- function(df, feature_names) {
  missing_cols <- setdiff(feature_names, names(df))
  extra_cols <- setdiff(names(df), feature_names)

  for (nm in missing_cols) {
    df[[nm]] <- NA
  }

  list(
    data = df[, feature_names, drop = FALSE],
    missing_cols = missing_cols,
    extra_cols = extra_cols
  )
}

resolve_external_paths <- function(cfg) {
  list(
    data_path = here::here(cfg$external_data$data_path)
  )
}

load_fitted_artifacts <- function(run_dir, preferred_survival_model_path = "") {
  run_obj <- readRDS(file.path(run_dir, "reports", "run_output.rds"))
  surv_dir <- file.path(run_dir, "survival_analysis")

  key_surv_paths <- list(
    sum_llr_only = file.path(surv_dir, "survival_model_sum_llr_only.rds"),
    sum_llr_age = file.path(surv_dir, "survival_model_sum_llr_age.rds"),
    p_rf = file.path(surv_dir, "survival_model_p_rf.rds"),
    p_rf_llr = file.path(surv_dir, "survival_model_p_rf_llr.rds"),
    groups_age = file.path(surv_dir, "survival_model_groups_age.rds")
  )

  key_surv_fits <- lapply(key_surv_paths, function(pth) {
    if (file.exists(pth)) readRDS(pth) else NULL
  })

  list(
    run_obj = run_obj,
    base_cfg = run_obj$config,
    marginal_fit = run_obj$marginal_fit,
    geometry_fit = run_obj$geometry_fit,
    group_ecdf_fit = run_obj$group_outputs$ecdf,
    rf_obj = run_obj$rf_obj,
    rf_llr_obj = run_obj$rf_llr_obj,
    nb_obj = run_obj$nb_obj,
    survival_models = key_surv_fits
  )
}

subset_geometry_fit <- function(geom_fit, keep_features) {
  idx <- match(keep_features, geom_fit$features)
  if (any(is.na(idx))) {
    stop("Some keep_features were not found in geom_fit$features")
  }

  Sigma_pos_sub <- geom_fit$Sigma_pos[idx, idx, drop = FALSE]
  Sigma_neg_sub <- geom_fit$Sigma_neg[idx, idx, drop = FALSE]

  list(
    features = geom_fit$features[idx],
    positive_label = geom_fit$positive_label,
    negative_label = geom_fit$negative_label,
    mu_pos = geom_fit$mu_pos[idx],
    mu_neg = geom_fit$mu_neg[idx],
    Sigma_pos = Sigma_pos_sub,
    Sigma_neg = Sigma_neg_sub,
    R_pos = chol(Sigma_pos_sub),
    R_neg = chol(Sigma_neg_sub),
    ridge = geom_fit$ridge
  )
}

make_quantile_bins <- function(x, probs = c(0, 0.25, 0.5, 0.75, 1), labels = NULL) {
  qs <- unique(stats::quantile(x, probs = probs, na.rm = TRUE))
  if (length(qs) < 3) return(rep(NA_character_, length(x)))
  n_bins <- length(qs) - 1
  if (is.null(labels)) labels <- paste0("Q", seq_len(n_bins))
  labels <- labels[seq_len(n_bins)]
  as.character(cut(x, breaks = qs, include.lowest = TRUE, labels = labels))
}

plot_km_by_score <- function(df, score_col, bins = 4) {
  library(survival)
  library(survminer)

  df <- df %>%
    dplyr::filter(!is.na(.data[[score_col]])) %>%
    dplyr::mutate(
      bin = dplyr::ntile(.data[[score_col]], bins)
    )

  if (length(unique(df$bin)) < 2) return(NULL)

  fit <- survival::survfit(
    survival::Surv(duration_hours_from_24h, event_after_24h) ~ bin,
    data = df
  )

  p <- survminer::ggsurvplot(
    fit,
    data = df,
    risk.table = TRUE,
    pval = TRUE,
    xlab = "Hours from 24h",
    ylab = "Survival probability",
    risk.table.y.text = FALSE,
    risk.table.height = 0.28,
    xlim = c(0, 720),
    break.x.by = 168,
    ylim = c(0.5, 1.0),
    size = 1.1,
    censor.size = 2.8,
    palette = "Dark2",
    ggtheme = plot_theme_mlhc(),
    tables.theme = plot_theme_mlhc(base_size = 11)
  )

  return(p)
}

run_stage_external_baselines <- function(ext_df,
                                         score_table,
                                         rf_obj,
                                         rf_llr_obj,
                                         nb_obj,
                                         L_ext_aligned = NULL,
                                         cfg) {
  y_col <- cfg$fit$y_col
  id_col <- cfg$data$id_col
  positive_label <- cfg$fit$positive

  rf_aligned <- align_to_feature_set(ext_df, rf_obj$feature_names)
  nb_aligned <- align_to_feature_set(ext_df, nb_obj$feature_names)

  p_rf <- score_random_forest_baseline(rf_aligned$data, rf_obj)
  p_nb <- score_naive_bayes_baseline(nb_aligned$data, nb_obj)

  p_rf_llr <- NULL
  rf_llr_missing_cols <- character(0)

  if (!is.null(rf_llr_obj)) {
    if (is.null(L_ext_aligned)) {
      stop("L_ext_aligned is required to score rf_llr_obj.")
    }
    rf_llr_missing_cols <- setdiff(rf_llr_obj$feature_names, colnames(L_ext_aligned))
    p_rf_llr <- score_random_forest_baseline(as.data.frame(L_ext_aligned), rf_llr_obj)
  }

  eval_tbl <- build_evaluation_table(
    df = ext_df,
    score_table = score_table,
    y_col = y_col,
    id_col = id_col,
    positive_label = positive_label,
    p_nb = p_nb,
    p_rf = p_rf,
    p_rf_llr = p_rf_llr
  )

  perf_tbl <- compute_performance_metrics(eval_tbl)
  curve_tbl <- compute_mortality_fraction_curves(eval_tbl, n_bins = 20)
  curve_plot <- plot_mortality_fraction_curves(curve_tbl)
  density_plots <- plot_probability_density_curves(eval_tbl)

  warn_list <- list()

  if (length(rf_aligned$missing_cols) > 0) {
    warn_list[[length(warn_list) + 1]] <- new_warning_row(
      stage = "external_baselines",
      feature = NA_character_,
      code = "rf_missing_features_filled_na",
      message = paste("Missing RF features filled with NA:", paste(rf_aligned$missing_cols, collapse = ", "))
    )
  }

  if (length(nb_aligned$missing_cols) > 0) {
    warn_list[[length(warn_list) + 1]] <- new_warning_row(
      stage = "external_baselines",
      feature = NA_character_,
      code = "nb_missing_features_filled_na",
      message = paste("Missing NB features filled with NA:", paste(nb_aligned$missing_cols, collapse = ", "))
    )
  }

  if (length(rf_llr_missing_cols) > 0) {
    warn_list[[length(warn_list) + 1]] <- new_warning_row(
      stage = "external_baselines",
      feature = NA_character_,
      code = "rf_llr_missing_features_filled_na",
      message = paste("Missing RF-LLR features filled with NA:", paste(rf_llr_missing_cols, collapse = ", "))
    )
  }

  warnings_df <- if (length(warn_list) > 0) dplyr::bind_rows(warn_list) else NULL

  new_stage_output(
    result = list(
      evaluation_table = eval_tbl,
      performance_metrics = perf_tbl,
      mortality_fraction_curves = curve_tbl,
      mortality_fraction_curve_plot = curve_plot,
      density_probability_plots = density_plots,
      p_rf = p_rf,
      p_rf_llr = p_rf_llr,
      p_nb = p_nb
    ),
    diagnostics = perf_tbl,
    warnings = warnings_df,
    dropped = NULL,
    metadata = list(
      rf_features = rf_obj$feature_names,
      rf_llr_features = if (!is.null(rf_llr_obj)) rf_llr_obj$feature_names else NULL,
      nb_features = nb_obj$feature_names,
      rf_missing_cols = rf_aligned$missing_cols,
      rf_llr_missing_cols = rf_llr_missing_cols,
      nb_missing_cols = nb_aligned$missing_cols
    )
  )
}
build_group_patient_summary_external <- function(L,
                                                 score_tbl,
                                                 risk_tbl,
                                                 stay_ids,
                                                 group_map_path,
                                                 groups_cfg,
                                                 ecdf_fit) {
  group_df <- read_group_map(group_map_path)

  s_pos <- group_positive_evidence(
    L = L,
    group_df = group_df,
    groups_cfg = groups_cfg
  )

  s_intra <- intra_group_entropy(
    L = L,
    group_df = s_pos$result$group_df,
    group_pos_sum = s_pos$result$group_pos_sum,
    groups_cfg = groups_cfg
  )

  s_inter <- inter_group_entropy(
    group_pos_sum = s_pos$result$group_pos_sum
  )

  s_pct <- score_group_percentiles(
    group_pos_sum = s_pos$result$group_pos_sum,
    ecdf_fit = ecdf_fit,
    groups_cfg = groups_cfg
  )

  patient_summary <- build_group_patient_summary(
    stay_id = stay_ids,
    sum_llr = score_tbl$sum_llr,
    d_dist = score_tbl$d_dist,
    risk_bin_primary = risk_tbl$risk_bin_primary,
    risk_bin_ddist = risk_tbl$risk_bin_ddist,
    group_percentiles = s_pct$result$group_percentiles,
    weighted_intra_entropy = s_intra$result$weighted_intra_entropy,
    inter_group_entropy = s_inter$result$inter_entropy
  )

  interp_tbl <- build_interpretation_strings(
    patient_summary = patient_summary,
    groups_cfg = groups_cfg
  )
  interp_tbl$stay_id <- patient_summary$stay_id

  patient_summary <- patient_summary %>%
    dplyr::left_join(interp_tbl, by = "stay_id")

  new_stage_output(
    result = list(
      patient_summary = patient_summary,
      group_percentiles = s_pct$result$group_percentiles,
      intra_entropy_detail = s_intra$result$intra_entropy_detail,
      group_pos_sum = s_pos$result$group_pos_sum
    ),
    diagnostics = bind_stage_diagnostics(
      s_pos$diagnostics,
      s_intra$diagnostics,
      s_inter$diagnostics,
      s_pct$diagnostics
    ),
    warnings = bind_stage_warnings(
      s_pos$warnings,
      s_intra$warnings,
      s_inter$warnings,
      s_pct$warnings
    ),
    dropped = NULL,
    metadata = list(
      groups = ecdf_fit$keep_groups,
      atypical_threshold = groups_cfg$atypical_threshold
    )
  )
}

compute_external_survival_concordance <- function(model, df) {
  lp <- predict(model, newdata = df, type = "lp")
  cfit <- survival::concordance(
    survival::Surv(duration_hours_from_24h, event_after_24h) ~ lp,
    data = df,
    reverse = TRUE
  )
  list(lp = lp, concordance = unname(cfit$concordance))
}

compute_external_time_roc <- function(df, marker, method_name, times = c(48,72,96,128)) {
  roc_obj <- timeROC::timeROC(
    T = df$duration_hours_from_24h,
    delta = df$event_after_24h,
    marker = marker,
    cause = 1,
    weighting = "marginal",
    times = times,
    iid = FALSE
  )
  data.frame(method = method_name, hours = times, auroc = as.numeric(roc_obj$AUC), stringsAsFactors = FALSE)
}

score_external_survival_models <- function(eval_surv_df, survival_models, times = c(48,72,96,128)) {
  model_names <- names(survival_models)
  model_names <- model_names[!vapply(survival_models, is.null, logical(1))]
  if (length(model_names) == 0) {
    return(new_stage_output(result = list(concordance = data.frame(), time_roc_long = data.frame()), diagnostics = NULL, warnings = NULL, dropped = NULL, metadata = NULL))
  }

  conc_rows <- list(); troc_rows <- list(); lp_store <- list()
  for (nm in model_names) {
    message("Scoring external survival model: ", nm)
    sc <- compute_external_survival_concordance(survival_models[[nm]], eval_surv_df)
    conc_rows[[nm]] <- data.frame(method = nm, concordance_eval = sc$concordance, stringsAsFactors = FALSE)
    troc_rows[[nm]] <- compute_external_time_roc(eval_surv_df, sc$lp, nm, times = times)
    lp_store[[nm]] <- sc$lp
  }

  new_stage_output(
    result = list(
      concordance = dplyr::bind_rows(conc_rows),
      time_roc_long = dplyr::bind_rows(troc_rows),
      linear_predictors = lp_store
    ),
    diagnostics = dplyr::bind_rows(conc_rows),
    warnings = NULL,
    dropped = NULL,
    metadata = NULL
  )
}

run_external_hospital_validation <- function(hospital_df,
                                             hospital_id,
                                             fitted,
                                             cfg,
                                             out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  hospital_df[[cfg$external_data$gender_col]] <- normalize_gender(hospital_df[[cfg$external_data$gender_col]])

  meta_df <- hospital_df %>%
    dplyr::select(
      stay_id = all_of(cfg$external_data$id_col),
      hospitalid = all_of(cfg$external_data$hospital_col),
      mortality = all_of(cfg$external_data$mortality_col),
      duration_hours_from_24h = all_of(cfg$external_data$duration_col),
      event_after_24h = all_of(cfg$external_data$event_col)
    )

  feature_df <- hospital_df %>%
    dplyr::select(-dplyr::any_of(c(
      cfg$external_data$hospital_col,
      cfg$external_data$duration_col,
      cfg$external_data$event_col
    )))

  groups_cfg <- fitted$base_cfg$groups
  groups_cfg$atypical_threshold <- cfg$groups$atypical_threshold
  groups_cfg$top_n_atypical <- cfg$groups$top_n_atypical

  message("Hospital ", hospital_id, " rows: ", nrow(feature_df))
  message("External feature count: ", ncol(feature_df))
  message("Marginal feature count from run_output: ", length(fitted$marginal_fit$features))
  message("Geometry feature count from run_output: ", length(fitted$geometry_fit$features))

  present_marginal <- intersect(fitted$marginal_fit$features, names(feature_df))
  missing_marginal <- setdiff(fitted$marginal_fit$features, names(feature_df))
  message("Present marginal features in external df: ", length(present_marginal))
  message("Missing marginal features in external df: ", length(missing_marginal))
  if (length(missing_marginal) > 0) {
    message("First missing marginal features: ", paste(utils::head(missing_marginal, 20), collapse = ", "))
  }

  aligned_ev <- align_to_feature_set(feature_df, fitted$marginal_fit$features)
  s_ev <- evidence_transform(
    new_df = aligned_ev$data,
    marginal_fit = fitted$marginal_fit,
    transform_cfg = fitted$base_cfg$transform
  )

  common_features <- intersect(colnames(s_ev$result$L), fitted$geometry_fit$features)
  missing_geometry_features <- setdiff(fitted$geometry_fit$features, colnames(s_ev$result$L))
  if (length(missing_geometry_features) > 0) {
    message("Hospital ", hospital_id, ": transformed L is missing geometry features: ", paste(missing_geometry_features, collapse = ", "))
  }
  if (length(common_features) == 0) stop("No common transformed features remain for geometry scoring in hospital ", hospital_id)

  L_ext <- s_ev$result$L[, common_features, drop = FALSE]

  L_ext_rf_llr <- NULL

  if (!is.null(fitted$rf_llr_obj)) {
    llr_aligned <- align_to_feature_set(
      as.data.frame(s_ev$result$L, stringsAsFactors = FALSE),
      fitted$rf_llr_obj$feature_names
    )
    L_ext_rf_llr <- llr_aligned$data
  }

  geom_fit_sub <- subset_geometry_fit(fitted$geometry_fit, common_features)
  s_dd <- distance_contrast_score(L_ext, geom_fit_sub)

  if (length(missing_geometry_features) > 0) {
    utils::write.csv(data.frame(feature = missing_geometry_features, stringsAsFactors = FALSE), file.path(out_dir, "missing_geometry_features.csv"), row.names = FALSE)
  }

  score_tbl <- data.frame(
    stay_id = feature_df[[cfg$external_data$id_col]],
    sum_llr = rowSums(L_ext),
    d_dist = s_dd$result$d_dist,
    stringsAsFactors = FALSE
  )

  write_stage_outputs("evidence_transform", s_ev, list(reports_dir = out_dir))
  write_stage_outputs("distance_contrast", s_dd, list(reports_dir = out_dir))
  utils::write.csv(score_tbl, file.path(out_dir, "score_table.csv"), row.names = FALSE)

  plot_list <- plot_distance_contrast(
    score_df = score_tbl[, c("sum_llr", "d_dist"), drop = FALSE],
    outcome = meta_df$mortality
  )
  ggplot2::ggsave(file.path(out_dir, "sum_llr_density.png"), plot_list$sum_llr_plot, width = 8, height = 5, dpi = 300)
  ggplot2::ggsave(file.path(out_dir, "d_dist_density.png"), plot_list$d_dist_plot, width = 8, height = 5, dpi = 300)

  risk_tbl <- data.frame(
    stay_id = score_tbl$stay_id,
    risk_bin_primary = make_quantile_bins(score_tbl$sum_llr),
    risk_bin_ddist = as.character(cut(
      score_tbl$d_dist,
      breaks = cfg$risk_bins$breaks,
      labels = cfg$risk_bins$labels,
      include.lowest = TRUE
    )),
    stringsAsFactors = FALSE
  )

  # summarize primary bins for external narrative
  risk_tbl_primary <- data.frame(
    stay_id = score_tbl$stay_id,
    sum_llr = score_tbl$sum_llr,
    risk_bin = risk_tbl$risk_bin_primary,
    stringsAsFactors = FALSE
  )

  risk_summary <- summarize_risk_bins(
    risk_tbl = risk_tbl_primary,
    outcome = coerce_outcome_binary01(meta_df$mortality, fitted$base_cfg$fit$positive),
    score_col = "sum_llr"
  )

  utils::write.csv(risk_tbl, file.path(out_dir, "risk_table.csv"), row.names = FALSE)
  utils::write.csv(risk_summary, file.path(out_dir, "risk_bin_summary.csv"), row.names = FALSE)

  risk_tbl_ddist <- data.frame(
    stay_id = score_tbl$stay_id,
    d_dist = score_tbl$d_dist,
    risk_bin = risk_tbl$risk_bin_ddist,
    stringsAsFactors = FALSE
  )

  risk_summary_ddist <- summarize_risk_bins(
    risk_tbl = risk_tbl_ddist,
    outcome = coerce_outcome_binary01(meta_df$mortality, fitted$base_cfg$fit$positive),
    score_col = "d_dist"
  )

  utils::write.csv(
    risk_summary_ddist,
    file.path(out_dir, "risk_bin_summary_ddist.csv"),
    row.names = FALSE
  )

  s_base <- run_stage_external_baselines(
    ext_df = feature_df,
    score_table = score_tbl,
    rf_obj = fitted$rf_obj,
    rf_llr_obj = fitted$rf_llr_obj,
    nb_obj = fitted$nb_obj,
    L_ext_aligned = L_ext_rf_llr,
    cfg = fitted$base_cfg
  )

  write_stage_outputs("external_baselines", s_base, list(reports_dir = out_dir))
  utils::write.csv(s_base$result$evaluation_table, file.path(out_dir, "evaluation_table.csv"), row.names = FALSE)
  utils::write.csv(s_base$result$performance_metrics, file.path(out_dir, "performance_metrics.csv"), row.names = FALSE)
  utils::write.csv(s_base$result$mortality_fraction_curves, file.path(out_dir, "mortality_fraction_curves.csv"), row.names = FALSE)
  ggplot2::ggsave(file.path(out_dir, "mortality_fraction_curve_sum_llr.png"), s_base$result$mortality_fraction_curve_plot$sum_llr, width = 8, height = 5, dpi = 300)
  ggplot2::ggsave(file.path(out_dir, "mortality_fraction_curve_d_dist.png"), s_base$result$mortality_fraction_curve_plot$d_dist, width = 8, height = 5, dpi = 300)
  ggplot2::ggsave(file.path(out_dir, "mortality_fraction_curve_p_nb.png"), s_base$result$mortality_fraction_curve_plot$p_nb, width = 8, height = 5, dpi = 300)
  ggplot2::ggsave(file.path(out_dir, "mortality_fraction_curve_p_rf.png"), s_base$result$mortality_fraction_curve_plot$p_rf, width = 8, height = 5, dpi = 300)

  if (!is.null(s_base$result$mortality_fraction_curve_plot$p_rf_llr)) {
    ggplot2::ggsave(
      file.path(out_dir, "mortality_fraction_curve_p_rf_llr.png"),
      s_base$result$mortality_fraction_curve_plot$p_rf_llr,
      width = 8, height = 5, dpi = 300
    )
  }

if (!is.null(s_base$result$density_probability_plots$p_rf)) {
  ggplot2::ggsave(
    file.path(out_dir, "density_p_rf.png"),
    s_base$result$density_probability_plots$p_rf,
    width = 8, height = 5, dpi = 300
  )
}

if (!is.null(s_base$result$density_probability_plots$p_rf_llr)) {
  ggplot2::ggsave(
    file.path(out_dir, "density_p_rf_llr.png"),
    s_base$result$density_probability_plots$p_rf_llr,
    width = 8, height = 5, dpi = 300
  )
}

  s_grp <- build_group_patient_summary_external(
    L = L_ext,
    score_tbl = score_tbl,
    risk_tbl = risk_tbl,
    stay_ids = score_tbl$stay_id,
    group_map_path = fitted$base_cfg$groups$group_map_path,
    groups_cfg = groups_cfg,
    ecdf_fit = fitted$group_ecdf_fit
  )
  write_stage_outputs("group_metrics", s_grp, list(reports_dir = out_dir))
  utils::write.csv(s_grp$result$patient_summary, file.path(out_dir, "group_patient_summary.csv"), row.names = FALSE)
  utils::write.csv(s_grp$result$intra_entropy_detail, file.path(out_dir, "intra_group_entropy_detail.csv"), row.names = FALSE)

  # ---- stage E5b: group transport summaries ----
  grp_transport <- summarize_group_percentile_transport(
    patient_summary = s_grp$result$patient_summary,
    atypical_threshold = groups_cfg$atypical_threshold,
    primary_bin_col = "risk_bin_primary"
  )

  utils::write.csv(
    grp_transport$group_atypical_frequency,
    file.path(out_dir, "group_atypical_frequency.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    grp_transport$patient_count_table,
    file.path(out_dir, "atypical_group_count_table.csv"),
    row.names = FALSE
  )

  utils::write.csv(
    grp_transport$patient_count_summary,
    file.path(out_dir, "atypical_group_count_summary.csv"),
    row.names = FALSE
  )

  if (!is.null(grp_transport$group_atypical_frequency_by_bin)) {
    utils::write.csv(
      grp_transport$group_atypical_frequency_by_bin,
      file.path(out_dir, "group_atypical_frequency_by_bin.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(grp_transport$risk_bin_summary_primary)) {
    utils::write.csv(
      grp_transport$risk_bin_summary_primary,
      file.path(out_dir, "risk_bin_summary_primary.csv"),
      row.names = FALSE
    )
  }

  grp_possum <- summarize_group_positive_sum_distribution(
    group_pos_sum = s_grp$result$group_pos_sum,
    risk_bin_primary = s_grp$result$patient_summary$risk_bin_primary
  )

  utils::write.csv(
    grp_possum$group_pos_sum_summary,
    file.path(out_dir, "group_pos_sum_summary.csv"),
    row.names = FALSE
  )

  if (!is.null(grp_possum$group_pos_sum_summary_by_bin)) {
    utils::write.csv(
      grp_possum$group_pos_sum_summary_by_bin,
      file.path(out_dir, "group_pos_sum_summary_by_bin.csv"),
      row.names = FALSE
    )
  }


  s_surv <- NULL
  if (isTRUE(cfg$survival$enabled)) {
    eval_tbl <- s_base$result$evaluation_table
    surv_df <- hospital_df %>%
      dplyr::select(
        stay_id = all_of(cfg$external_data$id_col),
        anchor_age,
        duration_hours_from_24h = all_of(cfg$external_data$duration_col),
        event_after_24h = all_of(cfg$external_data$event_col)
      ) %>%
      dplyr::left_join(eval_tbl, by = "stay_id") %>%
      dplyr::left_join(s_grp$result$patient_summary, by = c("stay_id", "sum_llr", "d_dist"))

    s_surv <- score_external_survival_models(
      eval_surv_df = surv_df,
      survival_models = fitted$survival_models,
      times = cfg$survival$time_roc_hours
    )

    # ---------- KM PLOTS FOR KEY MODELS ----------

km_models <- c("sum_llr_only", "sum_llr_age")

for (m in km_models) {
  if (m %in% names(fitted$survival_models) && !is.null(fitted$survival_models[[m]])) {
    lp <- predict(fitted$survival_models[[m]], newdata = surv_df, type = "lp")

    km_df <- surv_df
    km_df$model_lp <- lp

    km_plot <- plot_km_by_score(
      km_df,
      score_col = "model_lp",
      bins = 4
    )

    if (!is.null(km_plot)) {
      ggplot2::ggsave(
        filename = file.path(out_dir, paste0("km_curve_", m, ".png")),
        plot = km_plot$plot,
        width = 9,
        height = 6.5,
        dpi = 300
      )
    }
  }
}

    write_stage_outputs("survival", s_surv, list(reports_dir = out_dir))
    utils::write.csv(s_surv$result$concordance, file.path(out_dir, "external_survival_concordance.csv"), row.names = FALSE)
    utils::write.csv(s_surv$result$time_roc_long, file.path(out_dir, "external_survival_timeROC_long.csv"), row.names = FALSE)
  }

  all_warnings <- bind_stage_warnings(
    s_ev$warnings,
    s_dd$warnings,
    s_base$warnings,
    s_grp$warnings,
    if (!is.null(s_surv)) s_surv$warnings else NULL
  )
  if (!is.null(all_warnings)) {
    utils::write.csv(all_warnings, file.path(out_dir, "warnings.csv"), row.names = FALSE)
  }

  list(
    score_table = score_tbl,
    risk_table = risk_tbl,
    risk_summary = risk_summary,
    baselines = s_base$result,
    group_metrics = s_grp$result,
    survival = if (!is.null(s_surv)) s_surv$result else NULL
  )
}
