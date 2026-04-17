suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(timeROC)
  library(splines)
  library(tidyr)
})

find_latest_run_dir <- function(base_dir, prefix) {
  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  hits <- dirs[startsWith(basename(dirs), prefix)]
  if (length(hits) == 0) {
    stop("No matching run directories found under: ", base_dir)
  }
  hits[which.max(file.info(hits)$mtime)]
}

resolve_run_dir <- function(surv_cfg) {
  if (!is.null(surv_cfg$run_lookup$run_dir) && nzchar(surv_cfg$run_lookup$run_dir)) {
    return(here::here(surv_cfg$run_lookup$run_dir))
  }

  find_latest_run_dir(
    base_dir = here::here(surv_cfg$run_lookup$outputs_base_dir),
    prefix = surv_cfg$run_lookup$run_prefix
  )
}

resolve_data_paths <- function(surv_cfg, base_cfg) {
  train_path <- surv_cfg$data_sources$train_data_path
  eval_path  <- surv_cfg$data_sources$eval_data_path

  if (!nzchar(train_path)) {
    train_path <- base_cfg$data$train_data_path
  }
  if (!nzchar(eval_path)) {
    eval_label <- surv_cfg$data_sources$eval_label
    if (identical(eval_label, "val")) {
      eval_path <- base_cfg$data$val_data_path
    } else if (!is.null(base_cfg$data$test_data_path) && identical(eval_label, "test")) {
      eval_path <- base_cfg$data$test_data_path
    } else {
      stop("eval_data_path not provided in survival config, and could not infer it from base config.")
    }
  }

  list(
    train_data_path = here::here(train_path),
    eval_data_path = here::here(eval_path)
  )
}

read_survival_table <- function(path, duration_col, event_col, age_col = "anchor_age") {
  df <- read.csv(path, stringsAsFactors = FALSE)

  required <- c("stay_id", duration_col, event_col, age_col)
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Survival table missing columns: ", paste(missing, collapse = ", "))
  }

  out <- df %>%
    dplyr::select(stay_id, all_of(age_col), all_of(duration_col), all_of(event_col))

  names(out)[names(out) == age_col] <- "anchor_age"
  names(out)[names(out) == duration_col] <- "duration_hours_from_24h"
  names(out)[names(out) == event_col] <- "event_after_24h"

  out
}

read_severity_scores <- function(path, score_cols) {
  df <- read.csv(path, stringsAsFactors = FALSE)

  required <- c("stay_id", score_cols)
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Severity score table missing columns: ", paste(missing, collapse = ", "))
  }

  df %>%
    dplyr::select(all_of(required))
}

make_quantile_bins <- function(x, probs = c(0, 0.25, 0.5, 0.75, 1), labels = NULL) {
  qs <- unique(stats::quantile(x, probs = probs, na.rm = TRUE))

  if (length(qs) < 3) {
    return(rep(NA_character_, length(x)))
  }

  n_bins <- length(qs) - 1
  if (is.null(labels)) {
    labels <- paste0("Q", seq_len(n_bins))
  } else {
    labels <- labels[seq_len(n_bins)]
  }

  as.character(cut(
    x,
    breaks = qs,
    include.lowest = TRUE,
    labels = labels
  ))
}

plot_km_by_sum_llr_bin <- function(df, out_path = NULL,
                                   probs = c(0, 0.25, 0.5, 0.75, 1),
                                   labels = c("Q1", "Q2", "Q3", "Q4")) {
  df <- df %>%
    mutate(sum_llr_bin = make_quantile_bins(sum_llr, probs = probs, labels = labels))

  km_fit <- survfit(
    Surv(duration_hours_from_24h, event_after_24h) ~ sum_llr_bin,
    data = df
  )

  p <- ggsurvplot(
    km_fit,
    data = df,
    risk.table = TRUE,
    pval = TRUE,
    legend.labs = unique(na.omit(df$sum_llr_bin)),
    legend.title = "sum_llr group",
    xlab = "Hours from 24h landmark",
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

  if (!is.null(out_path)) {
    ggplot2::ggsave(out_path, p$plot, width = 8.5, height = 5.5, dpi = 300)
  }

  list(fit = km_fit, plot = p)
}

plot_km_by_score_quantile <- function(df,
                                      score_col,
                                      out_path = NULL,
                                      probs = c(0, 0.25, 0.5, 0.75, 1),
                                      labels = c("Q1", "Q2", "Q3", "Q4"),
                                      legend_title = NULL) {
  if (!(score_col %in% names(df))) {
    stop("score_col not found in df: ", score_col)
  }

  if (is.null(legend_title)) {
    legend_title <- score_col
  }

  df <- df %>%
    dplyr::mutate(score_bin = make_quantile_bins(.data[[score_col]], probs = probs, labels = labels))

  km_fit <- survival::survfit(
    survival::Surv(duration_hours_from_24h, event_after_24h) ~ score_bin,
    data = df
  )

  p <- survminer::ggsurvplot(
    km_fit,
    data = df,
    risk.table = TRUE,
    pval = TRUE,
    legend.labs = unique(stats::na.omit(df$score_bin)),
    legend.title = legend_title,
    xlab = "Hours from 24h landmark",
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

  if (!is.null(out_path)) {
    ggplot2::ggsave(out_path, p$plot, width = 8.5, height = 5.5, dpi = 300)
  }

  list(fit = km_fit, plot = p)
}

make_ddist_bins <- function(d_dist, bin_cfg) {
  cut(
    d_dist,
    breaks = bin_cfg$breaks,
    labels = bin_cfg$labels
  )
}

build_score_table <- function(df, marginal_fit, geometry_fit, transform_cfg, id_col) {
  s_tr <- evidence_transform(df, marginal_fit, transform_cfg)
  L <- s_tr$result$L
  common_features <- intersect(colnames(L), geometry_fit$features)
  L <- L[, common_features, drop = FALSE]
  sum_llr <- rowSums(L)

  s_dc <- distance_contrast_score(L, geometry_fit)

  data.frame(
    stay_id = df[[id_col]],
    sum_llr = sum_llr,
    d_dist = s_dc$result$d_dist,
    stringsAsFactors = FALSE
  )
}

build_survival_dataset <- function(processed_df,
                                   score_tbl,
                                   patient_summary = NULL,
                                   survival_tbl,
                                   severity_tbl,
                                   id_col) {
  out <- processed_df %>%
    dplyr::select(all_of(id_col)) %>%
    dplyr::rename(stay_id = all_of(id_col)) %>%
    dplyr::left_join(score_tbl, by = "stay_id") %>%
    dplyr::left_join(survival_tbl, by = "stay_id") %>%
    dplyr::left_join(severity_tbl, by = "stay_id")

  if (!is.null(patient_summary)) {
    keep_cols <- setdiff(
      names(patient_summary),
      c("sum_llr", "d_dist", "risk_bin")
    )
    out <- out %>% dplyr::left_join(patient_summary[, keep_cols, drop = FALSE], by = "stay_id")
  }

  out
}

plot_km_by_ddist_bin <- function(df, ddist_bin_cfg, out_path = NULL) {
  df <- df %>% mutate(d_dist_bin = make_ddist_bins(d_dist, ddist_bin_cfg))

  km_fit <- survfit(
    Surv(duration_hours_from_24h, event_after_24h) ~ d_dist_bin,
    data = df
  )

  p <- ggsurvplot(
    km_fit,
    data = df,
    risk.table = TRUE,
    pval = TRUE,
    legend.labs = ddist_bin_cfg$labels,
    legend.title = "d_dist group",
    xlab = "Hours from 24h landmark",
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

  if (!is.null(out_path)) {
    ggplot2::ggsave(out_path, p$plot, width = 8.5, height = 5.5, dpi = 300)
  }

  list(fit = km_fit, plot = p)
}

fit_cox_model <- function(formula_obj, train_df) {
  survival::coxph(formula_obj, data = train_df, x = TRUE, y = TRUE)
}

fit_with_warning_capture <- function(formula_obj, train_df, model_name) {
  withCallingHandlers(
    survival::coxph(formula_obj, data = train_df, x = TRUE, y = TRUE),
    warning = function(w) {
      message("Warning while fitting model: ", model_name)
      message("  ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
}

extract_cox_train_summary <- function(model, method_name) {
  sm <- summary(model)

  data.frame(
    method = method_name,
    concordance_train = unname(sm$concordance[1]),
    lr_chisq_train = unname(sm$logtest["test"]),
    lr_df_train = unname(sm$logtest["df"]),
    lr_p_train = unname(sm$logtest["pvalue"]),
    stringsAsFactors = FALSE
  )
}

extract_cox_coef_table <- function(model, method_name) {
  sm <- summary(model)
  cf <- as.data.frame(sm$coefficients)
  cf$term <- rownames(cf)
  rownames(cf) <- NULL
  cf$method <- method_name
  cf
}

compute_val_concordance <- function(model, val_df) {
  lp <- predict(model, newdata = val_df, type = "lp")

  cfit <- survival::concordance(
    Surv(duration_hours_from_24h, event_after_24h) ~ lp,
    data = val_df,
    reverse = TRUE
  )

  list(
    lp = lp,
    concordance = unname(cfit$concordance)
  )
}

compute_time_roc_table <- function(df, marker, method_name, times = c(48, 72, 96, 128)) {
  roc_obj <- timeROC::timeROC(
    T = df$duration_hours_from_24h,
    delta = df$event_after_24h,
    marker = marker,
    cause = 1,
    weighting = "marginal",
    times = times,
    iid = FALSE
  )

  data.frame(
    method = method_name,
    hours = times,
    auroc = as.numeric(roc_obj$AUC),
    stringsAsFactors = FALSE
  )
}

get_group_feature_terms <- function(df, severity_cols, baseline_cols, include_age = TRUE) {
  exclude <- c(
    "stay_id",
    "anchor_age",
    "sum_llr", "d_dist",
    "risk_bin", "risk_bin_primary", "risk_bin_ddist",
    "weighted_intra_entropy", "inter_group_entropy",
    "interpretation",
    "atypical_group_count",
    "top_group_name",
    "top_group_percentile",
    "atypical_groups_csv",
    "duration_hours_from_24h", "event_after_24h"
  )

  exclude <- unique(c(exclude, severity_cols, baseline_cols))

  terms <- setdiff(names(df), exclude)
  terms <- terms[vapply(df[terms], is.numeric, logical(1))]

  if (!include_age) {
    terms <- setdiff(terms, "anchor_age")
  }

  terms
}

build_model_registry <- function(train_surv_df, surv_cfg) {
  severity_cols <- surv_cfg$benchmarks$severity_score_cols
  baseline_cols <- surv_cfg$benchmarks$baseline_score_cols

  group_terms <- get_group_feature_terms(
    df = train_surv_df,
    severity_cols = severity_cols,
    baseline_cols = baseline_cols,
    include_age = surv_cfg$benchmarks$include_age
  )

  group_rhs <- paste(c(group_terms, "weighted_intra_entropy", "inter_group_entropy"), collapse = " + ")

  registry <- list(
    age_only = as.formula("Surv(duration_hours_from_24h, event_after_24h) ~ anchor_age"),

    sum_llr_only = as.formula("Surv(duration_hours_from_24h, event_after_24h) ~ sum_llr"),
    sum_llr_age = as.formula("Surv(duration_hours_from_24h, event_after_24h) ~ sum_llr + anchor_age"),
    sum_llr_groups = as.formula(paste("Surv(duration_hours_from_24h, event_after_24h) ~ sum_llr + ", group_rhs)),
    sum_llr_groups_age = as.formula(paste("Surv(duration_hours_from_24h, event_after_24h) ~ sum_llr + anchor_age + ", group_rhs)),
    spline_sum_llr_groups_age = as.formula(paste("Surv(duration_hours_from_24h, event_after_24h) ~ splines::ns(sum_llr, df = 4) + anchor_age + ", group_rhs)),

    d_dist_only = as.formula("Surv(duration_hours_from_24h, event_after_24h) ~ d_dist"),
    groups_only = as.formula(paste("Surv(duration_hours_from_24h, event_after_24h) ~ ", group_rhs)),
    groups_age = as.formula(paste("Surv(duration_hours_from_24h, event_after_24h) ~ anchor_age + ", group_rhs)),
    d_dist_groups = as.formula(paste("Surv(duration_hours_from_24h, event_after_24h) ~ d_dist + ", group_rhs)),
    d_dist_groups_age = as.formula(paste("Surv(duration_hours_from_24h, event_after_24h) ~ d_dist + anchor_age + ", group_rhs)),
    spline_ddist_groups_age = as.formula(paste("Surv(duration_hours_from_24h, event_after_24h) ~ splines::ns(d_dist, df = 4) + anchor_age + ", group_rhs))
  )

  for (sc in severity_cols) {
    nm <- paste0(sc, if (surv_cfg$benchmarks$include_age) "_age" else "_only")
    rhs <- if (surv_cfg$benchmarks$include_age) paste("anchor_age +", sc) else sc
    registry[[nm]] <- as.formula(paste("Surv(duration_hours_from_24h, event_after_24h) ~", rhs))
  }

  for (sc in baseline_cols) {
  if (sc %in% names(train_surv_df)) {
    nm <- sc
    rhs <- sc
    registry[[nm]] <- as.formula(
      paste("Surv(duration_hours_from_24h, event_after_24h) ~", rhs)
    )
  }
}

  registry
}

run_stage9_survival <- function(train_surv_df,
                                eval_surv_df,
                                out_dir,
                                surv_cfg) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  registry <- build_model_registry(train_surv_df, surv_cfg)

  train_summary_rows <- list()
  coef_tables <- list()
  eval_perf_rows <- list()
  time_roc_rows <- list()
  km_candidate_scores <- list()
  fitted_models <- list()

  message("Available baseline columns in train_surv_df: ",
        paste(intersect(c("p_nb", "p_rf", "p_rf_llr"), names(train_surv_df)), collapse = ", "))

  for (nm in names(registry)) {
    message("Fitting survival model: ", nm)
    form <- registry[[nm]]
    fit <- fit_with_warning_capture(form, train_surv_df, nm)
    fitted_models[[nm]] <- fit

    train_summary_rows[[nm]] <- extract_cox_train_summary(fit, nm)
    coef_tables[[nm]] <- extract_cox_coef_table(fit, nm)

    eval_sc <- compute_val_concordance(fit, eval_surv_df)

    eval_perf_rows[[nm]] <- data.frame(
      method = nm,
      concordance_eval = eval_sc$concordance,
      stringsAsFactors = FALSE
    )

    time_roc_rows[[nm]] <- compute_time_roc_table(
      df = eval_surv_df,
      marker = eval_sc$lp,
      method_name = nm,
      times = surv_cfg$time_roc$hours
    )

    km_candidate_scores[[nm]] <- eval_sc$lp
  }

  if ("sum_llr_age" %in% names(km_candidate_scores)) {
  eval_km_sum_age <- eval_surv_df %>%
    dplyr::mutate(sum_llr_age_lp = km_candidate_scores[["sum_llr_age"]])

  plot_km_by_score_quantile(
    df = eval_km_sum_age,
    score_col = "sum_llr_age_lp",
    out_path = file.path(out_dir, "km_sum_llr_age_eval.png"),
    legend_title = "sum_llr + age"
  )
}

  train_summary_tbl <- dplyr::bind_rows(train_summary_rows)
  coef_tbl <- dplyr::bind_rows(coef_tables)
  eval_perf_tbl <- dplyr::bind_rows(eval_perf_rows)
  time_roc_long <- dplyr::bind_rows(time_roc_rows)

  overall_tbl <- train_summary_tbl %>%
    dplyr::left_join(eval_perf_tbl, by = "method")

  time_roc_wide <- tidyr::pivot_wider(
    time_roc_long,
    names_from = hours,
    values_from = auroc,
    names_prefix = "auroc_"
  )

  utils::write.csv(overall_tbl, file.path(out_dir, "cox_model_summary.csv"), row.names = FALSE)
  utils::write.csv(coef_tbl, file.path(out_dir, "cox_coefficients.csv"), row.names = FALSE)
  utils::write.csv(time_roc_long, file.path(out_dir, "timeROC_long.csv"), row.names = FALSE)
  utils::write.csv(time_roc_wide, file.path(out_dir, "timeROC_wide.csv"), row.names = FALSE)

  # also save key comparison models explicitly if present
  key_models_to_save <- c("sum_llr_only", "sum_llr_age", "groups_age", "p_nb", "p_rf", "p_rf_llr")

  p_conc <- ggplot(overall_tbl%>%filter(method %in% c(key_models_to_save, "apsiii_age")), aes(x = reorder(method, -concordance_eval), y = concordance_eval)) +
    geom_col() +
    coord_flip() +
    labs(title = "Evaluation concordance by method", x = NULL, y = "Concordance") +
    plot_theme_mlhc()

  ggsave(file.path(out_dir, "concordance_barplot.png"), p_conc, width = 8.5, height = 5.5, dpi = 300)

  p_troc <- ggplot(time_roc_long%>%filter(method %in% c(key_models_to_save, "apsiii_age")), aes(x = hours, y = auroc, color = method)) +
    geom_line() +
    geom_point() +
    labs(title = "Time-dependent ROC AUC", x = "Hours", y = "AUROC") +
    plot_theme_mlhc()

  ggsave(file.path(out_dir, "timeROC_plot.png"), p_troc, width = 8.5, height = 5.5, dpi = 300)

  preferred_method <- surv_cfg$output$preferred_km_method
  if (!(preferred_method %in% names(km_candidate_scores))) {
    preferred_method <- overall_tbl$method[which.max(overall_tbl$concordance_eval)]
  }

  preferred_fit <- fitted_models[[preferred_method]]

  saveRDS(
    preferred_fit,
    file.path(out_dir, paste0("preferred_survival_model_", preferred_method, ".rds"))
  )



  for (nm in key_models_to_save) {
    if (nm %in% names(fitted_models) && !is.null(fitted_models[[nm]])) {
      saveRDS(
        fitted_models[[nm]],
        file.path(out_dir, paste0("survival_model_", nm, ".rds"))
      )
    }
  }

  preferred_model_meta <- list(
    preferred_method = preferred_method,
    formula = deparse(formula(preferred_fit)),
    concordance_train = overall_tbl$concordance_train[match(preferred_method, overall_tbl$method)],
    concordance_eval = overall_tbl$concordance_eval[match(preferred_method, overall_tbl$method)],
    saved_at = as.character(Sys.time())
  )

  saveRDS(
    preferred_model_meta,
    file.path(out_dir, "preferred_survival_model_meta.rds")
  )

  key_model_meta <- lapply(key_models_to_save, function(nm) {
  if (!(nm %in% overall_tbl$method)) return(NULL)

  list(
    method = nm,
    formula = deparse(formula(fitted_models[[nm]])),
    concordance_train = overall_tbl$concordance_train[match(nm, overall_tbl$method)],
    concordance_eval = overall_tbl$concordance_eval[match(nm, overall_tbl$method)]
  )
})

  key_model_meta <- key_model_meta[!vapply(key_model_meta, is.null, logical(1))]

  saveRDS(
    key_model_meta,
    file.path(out_dir, "saved_key_survival_models_meta.rds")
  )

  lp_pref <- km_candidate_scores[[preferred_method]]
  qbreaks <- unique(stats::quantile(lp_pref, probs = seq(0, 1, 0.25), na.rm = TRUE))

  if (length(qbreaks) >= 5) {
    eval_km_df <- eval_surv_df %>%
      mutate(
        risk_quartile = cut(
          lp_pref,
          breaks = qbreaks,
          include.lowest = TRUE,
          labels = c("Q1", "Q2", "Q3", "Q4")
        )
      )

    km_pref <- survfit(
      Surv(duration_hours_from_24h, event_after_24h) ~ risk_quartile,
      data = eval_km_df
    )

    p_km <- ggsurvplot(
      km_pref,
      data = eval_km_df,
      risk.table = TRUE,
      pval = TRUE,
      legend.title = preferred_method,
      xlab = "Hours from 24h landmark",
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

    ggsave(file.path(out_dir, "km_preferred_model.png"), p_km$plot, width = 8.5, height = 5.5, dpi = 300)
  }

  list(
    overall = overall_tbl,
    coefficients = coef_tbl,
    time_roc_long = time_roc_long,
    time_roc_wide = time_roc_wide
  )
}