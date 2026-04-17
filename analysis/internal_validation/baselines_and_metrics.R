suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ranger)
  library(PRROC)
  library(pROC)
})

fit_formula_from_features <- function(y_col, feature_names) {
  stats::as.formula(
    paste(y_col, "~", paste(feature_names, collapse = " + "))
  )
}

get_model_features <- function(train_df, y_col, id_col = NULL, explicit_features = NULL) {
  if (!is.null(explicit_features)) return(explicit_features)

  feats <- setdiff(names(train_df), y_col)
  if (!is.null(id_col) && nzchar(id_col)) {
    feats <- setdiff(feats, id_col)
  }
  feats
}

coerce_outcome_binary01 <- function(y, positive_label = "1") {
  y_chr <- as.character(y)
  as.integer(y_chr == positive_label)
}

safe_numeric_prob <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  x
}

fit_random_forest_baseline <- function(train_df,
                                       y_col,
                                       positive_label = "1",
                                       id_col = NULL,
                                       feature_names = NULL,
                                       num_trees = 500,
                                       seed = 1) {
  feature_names <- get_model_features(train_df, y_col, id_col, feature_names)

  df <- train_df[, c(y_col, feature_names), drop = FALSE]
  df[[y_col]] <- as.factor(df[[y_col]])

  set.seed(seed)

  rf_fit <- ranger::ranger(
    formula = fit_formula_from_features(y_col, feature_names),
    data = df,
    probability = TRUE,
    num.trees = num_trees,
    respect.unordered.factors = "order",
    seed = seed
  )

  list(
    model = rf_fit,
    feature_names = feature_names,
    y_col = y_col,
    positive_label = positive_label
  )
}

score_random_forest_baseline <- function(df, rf_obj) {
  xdf <- df[, rf_obj$feature_names, drop = FALSE]
  pred <- predict(rf_obj$model, data = xdf)$predictions

  pos_col <- which(colnames(pred) == rf_obj$positive_label)
  if (length(pos_col) != 1) {
    stop("Could not identify positive class column in RF predictions.")
  }

  safe_numeric_prob(pred[, pos_col])
}

fit_naive_bayes_baseline <- function(train_df,
                                     y_col,
                                     positive_label = "1",
                                     id_col = NULL,
                                     feature_names = NULL,
                                     use_naivebayes_pkg = TRUE) {
  feature_names <- get_model_features(train_df, y_col, id_col, feature_names)

  df <- train_df[, c(y_col, feature_names), drop = FALSE]
  df[[y_col]] <- as.factor(df[[y_col]])

  formula_obj <- fit_formula_from_features(y_col, feature_names)

  if (use_naivebayes_pkg) {
    if (!requireNamespace("naivebayes", quietly = TRUE)) {
      stop("Package 'naivebayes' is required for naive Bayes baseline.")
    }

    nb_fit <- naivebayes::naive_bayes(
      formula = formula_obj,
      data = df
    )

    engine <- "naivebayes"
  } else {
    if (!requireNamespace("e1071", quietly = TRUE)) {
      stop("Package 'e1071' is required for naive Bayes baseline.")
    }

    nb_fit <- e1071::naiveBayes(
      x = df[, feature_names, drop = FALSE],
      y = df[[y_col]]
    )

    engine <- "e1071"
  }

  list(
    model = nb_fit,
    engine = engine,
    feature_names = feature_names,
    y_col = y_col,
    positive_label = positive_label
  )
}

score_naive_bayes_baseline <- function(df, nb_obj) {
  xdf <- df[, nb_obj$feature_names, drop = FALSE]

  if (nb_obj$engine == "naivebayes") {
    if (!requireNamespace("naivebayes", quietly = TRUE)) {
      stop("Package 'naivebayes' is required for naive Bayes baseline scoring.")
    }
    pred <- naivebayes:::predict.naive_bayes(nb_obj$model, newdata = xdf, type = "prob")
  } else {
    if (!requireNamespace("e1071", quietly = TRUE)) {
      stop("Package 'e1071' is required for naive Bayes baseline scoring.")
    }
    pred <- predict(nb_obj$model, newdata = xdf, type = "raw")
  }

  pos_col <- which(colnames(pred) == nb_obj$positive_label)
  if (length(pos_col) != 1) {
    stop("Could not identify positive class column in NB predictions.")
  }

  safe_numeric_prob(pred[, pos_col])
}

# NEW: generic density plot for probability-like scores
plot_probability_density_one <- function(eval_tbl, score_col) {
  df <- eval_tbl %>%
    dplyr::filter(is.finite(.data[[score_col]]), !is.na(mortality)) %>%
    dplyr::mutate(
      outcome = factor(
        mortality,
        levels = c(0, 1),
        labels = c("Negative", "Positive")
      )
    )

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[score_col]], color = outcome, fill = outcome)) +
    ggplot2::geom_density(alpha=0.2, linewidth = 1.1) +
    ggplot2::labs(
      title = paste("Density of", score_col),
      x = score_col,
      y = "Density"
    ) +
    plot_theme_mlhc()
}

plot_probability_density_curves <- function(eval_tbl) {
  out <- list()
  if ("p_rf" %in% names(eval_tbl)) {
    out$p_rf <- plot_probability_density_one(eval_tbl, "p_rf")
  }
  if ("p_rf_llr" %in% names(eval_tbl)) {
    out$p_rf_llr <- plot_probability_density_one(eval_tbl, "p_rf_llr")
  }
  out
}

build_evaluation_table <- function(df,
                                   score_table,
                                   y_col,
                                   id_col,
                                   positive_label,
                                   p_nb,
                                   p_rf,
                                   p_rf_llr = NULL) {
  if (!(id_col %in% names(df))) stop("id_col not found in df")
  if (!(y_col %in% names(df))) stop("y_col not found in df")
  if (!(id_col %in% names(score_table))) stop("id_col not found in score_table")

  out <- score_table %>%
    dplyr::left_join(
      df[, c(id_col, y_col), drop = FALSE],
      by = stats::setNames(id_col, id_col)
    )

  out$mortality <- coerce_outcome_binary01(out[[y_col]], positive_label = positive_label)
  out$p_nb <- p_nb
  out$p_rf <- p_rf
  if (!is.null(p_rf_llr)) {
    out$p_rf_llr <- p_rf_llr
  }

  keep_cols <- c(
    id_col,
    "mortality",
    "sum_llr",
    "d_dist",
    "p_nb",
    "p_rf"
  )
  if (!is.null(p_rf_llr)) {
    keep_cols <- c(keep_cols, "p_rf_llr")
  }

  out %>% dplyr::select(dplyr::all_of(keep_cols))
}

compute_binary_metrics_one <- function(y, score, score_name) {
  ok <- is.finite(score) & !is.na(y)
  y_ok <- y[ok]
  s_ok <- score[ok]

  if (length(unique(y_ok)) < 2) {
    return(data.frame(
      score_name = score_name,
      auroc = NA_real_,
      auprc = NA_real_,
      n = length(y_ok),
      stringsAsFactors = FALSE
    ))
  }

  roc_obj <- pROC::roc(response = y_ok, predictor = s_ok, quiet = TRUE, direction = "<")
  auroc <- as.numeric(pROC::auc(roc_obj))

  pr_obj <- PRROC::pr.curve(
    scores.class0 = s_ok[y_ok == 1],
    scores.class1 = s_ok[y_ok == 0],
    curve = FALSE
  )
  auprc <- as.numeric(pr_obj$auc.integral)

  data.frame(
    score_name = score_name,
    auroc = auroc,
    auprc = auprc,
    n = length(y_ok),
    stringsAsFactors = FALSE
  )
}

compute_performance_metrics <- function(eval_tbl) {
  rows <- list(
    compute_binary_metrics_one(eval_tbl$mortality, eval_tbl$sum_llr, "sum_llr"),
    compute_binary_metrics_one(eval_tbl$mortality, eval_tbl$d_dist, "d_dist"),
    compute_binary_metrics_one(eval_tbl$mortality, eval_tbl$p_nb, "p_nb"),
    compute_binary_metrics_one(eval_tbl$mortality, eval_tbl$p_rf, "p_rf")
  )

  if ("p_rf_llr" %in% names(eval_tbl)) {
    rows[[length(rows) + 1]] <- compute_binary_metrics_one(
      eval_tbl$mortality, eval_tbl$p_rf_llr, "p_rf_llr"
    )
  }

  rows[[length(rows) + 1]] <- data.frame(
    score_name = "prevalence_baseline",
    auroc = NA_real_,
    auprc = mean(eval_tbl$mortality == 1, na.rm = TRUE),
    n = nrow(eval_tbl),
    stringsAsFactors = FALSE
  )

  dplyr::bind_rows(rows)
}

make_score_bins <- function(score, n_bins = 20) {
  qs <- unique(stats::quantile(
    score,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE,
    type = 7
  ))

  if (length(qs) < 3) {
    return(rep(NA_character_, length(score)))
  }

  as.character(cut(
    score,
    breaks = qs,
    include.lowest = TRUE,
    ordered_result = FALSE
  ))
}

compute_mortality_fraction_curve_one <- function(eval_tbl, score_col, n_bins = 20) {
  x <- eval_tbl[[score_col]]
  y <- eval_tbl$mortality

  bins <- make_score_bins(x, n_bins = n_bins)

  tmp <- data.frame(
    score = x,
    mortality = y,
    bin = bins,
    stringsAsFactors = FALSE
  )

  tmp %>%
    dplyr::filter(!is.na(bin), is.finite(score)) %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      score_name = score_col,
      mean_score = mean(score, na.rm = TRUE),
      mortality_fraction = mean(mortality == 1, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
}

compute_mortality_fraction_curves <- function(eval_tbl, n_bins = 20) {
  rows <- list(
    compute_mortality_fraction_curve_one(eval_tbl, "sum_llr", n_bins = n_bins),
    compute_mortality_fraction_curve_one(eval_tbl, "d_dist", n_bins = n_bins),
    compute_mortality_fraction_curve_one(eval_tbl, "p_nb", n_bins = n_bins),
    compute_mortality_fraction_curve_one(eval_tbl, "p_rf", n_bins = n_bins)
  )

  if ("p_rf_llr" %in% names(eval_tbl)) {
    rows[[length(rows) + 1]] <- compute_mortality_fraction_curve_one(
      eval_tbl, "p_rf_llr", n_bins = n_bins
    )
  }

  dplyr::bind_rows(rows)
}

plot_mortality_fraction_curve_one <- function(curve_df, score_name) {
  df <- curve_df %>% dplyr::filter(score_name == !!score_name)

  ggplot2::ggplot(df, ggplot2::aes(x = mean_score, y = mortality_fraction)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = paste("Mortality fraction by", score_name),
      x = paste("Mean", score_name, "in bin"),
      y = "Mortality fraction"
    ) +
    plot_theme_mlhc()
}

plot_mortality_fraction_curves <- function(curve_df) {
  out <- list(
    sum_llr = plot_mortality_fraction_curve_one(curve_df, "sum_llr"),
    d_dist  = plot_mortality_fraction_curve_one(curve_df, "d_dist"),
    p_nb    = plot_mortality_fraction_curve_one(curve_df, "p_nb"),
    p_rf    = plot_mortality_fraction_curve_one(curve_df, "p_rf")
  )

  if ("p_rf_llr" %in% unique(curve_df$score_name)) {
    out$p_rf_llr <- plot_mortality_fraction_curve_one(curve_df, "p_rf_llr")
  }

  out
}

run_stage8_baselines <- function(train_df,
                                 val_df,
                                 train_score_table,
                                 val_score_table,
                                 cfg,
                                 L_train = NULL,
                                 L_val = NULL) {
  y_col <- cfg$fit$y_col
  id_col <- cfg$data$id_col
  positive_label <- cfg$fit$positive

  feature_names <- get_model_features(
    train_df = train_df,
    y_col = y_col,
    id_col = id_col,
    explicit_features = cfg$fit$features
  )

  rf_obj <- fit_random_forest_baseline(
    train_df = train_df,
    y_col = y_col,
    positive_label = positive_label,
    id_col = id_col,
    feature_names = feature_names,
    num_trees = 500,
    seed = cfg$run$seed
  )

  nb_obj <- fit_naive_bayes_baseline(
    train_df = train_df,
    y_col = y_col,
    positive_label = positive_label,
    id_col = id_col,
    feature_names = feature_names,
    use_naivebayes_pkg = TRUE
  )

  # NEW: RF on LLR-transformed space
  rf_llr_obj <- NULL
  p_rf_llr_train <- NULL
  p_rf_llr_val <- NULL

  if (!is.null(L_train) && !is.null(L_val)) {
    train_llr_df <- as.data.frame(L_train, stringsAsFactors = FALSE)
    train_llr_df[[id_col]] <- train_df[[id_col]]
    train_llr_df[[y_col]] <- train_df[[y_col]]

    val_llr_df <- as.data.frame(L_val, stringsAsFactors = FALSE)
    val_llr_df[[id_col]] <- val_df[[id_col]]
    val_llr_df[[y_col]] <- val_df[[y_col]]

    rf_llr_obj <- fit_random_forest_baseline(
      train_df = train_llr_df,
      y_col = y_col,
      positive_label = positive_label,
      id_col = id_col,
      feature_names = colnames(L_train),
      num_trees = 500,
      seed = cfg$run$seed
    )

    p_rf_llr_train <- score_random_forest_baseline(train_llr_df, rf_llr_obj)
    p_rf_llr_val <- score_random_forest_baseline(val_llr_df, rf_llr_obj)
  }

  p_rf_train <- score_random_forest_baseline(train_df, rf_obj)
  p_nb_train <- score_naive_bayes_baseline(train_df, nb_obj)

  p_rf_val <- score_random_forest_baseline(val_df, rf_obj)
  p_nb_val <- score_naive_bayes_baseline(val_df, nb_obj)

  eval_tbl_train <- build_evaluation_table(
    df = train_df,
    score_table = train_score_table,
    y_col = y_col,
    id_col = id_col,
    positive_label = positive_label,
    p_nb = p_nb_train,
    p_rf = p_rf_train,
    p_rf_llr = p_rf_llr_train
  )

  eval_tbl_val <- build_evaluation_table(
    df = val_df,
    score_table = val_score_table,
    y_col = y_col,
    id_col = id_col,
    positive_label = positive_label,
    p_nb = p_nb_val,
    p_rf = p_rf_val,
    p_rf_llr = p_rf_llr_val
  )

  perf_tbl <- compute_performance_metrics(eval_tbl_val)
  curve_tbl <- compute_mortality_fraction_curves(eval_tbl_val, n_bins = 20)
  curve_plot <- plot_mortality_fraction_curves(curve_tbl)
  density_plots <- plot_probability_density_curves(eval_tbl_val)

  new_stage_output(
    result = list(
      rf_obj = rf_obj,
      rf_llr_obj = rf_llr_obj,
      nb_obj = nb_obj,
      evaluation_table_train = eval_tbl_train,
      evaluation_table_val = eval_tbl_val,
      performance_metrics = perf_tbl,
      mortality_fraction_curves = curve_tbl,
      mortality_fraction_curve_plot = curve_plot,
      density_probability_plots = density_plots
    ),
    diagnostics = perf_tbl,
    warnings = NULL,
    dropped = NULL,
    metadata = list(
      baseline_features = feature_names,
      llr_features = if (!is.null(L_train)) colnames(L_train) else NULL
    )
  )
}