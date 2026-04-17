fit_marginals <- function(train_df, fit_cfg, exclude_features = NULL) {
  y_col <- fit_cfg$y_col
  y <- train_df[[y_col]]
  if (!is.factor(y)) y <- as.factor(y)
  train_df[[y_col]] <- y
  
  pos_label <- fit_cfg$positive
  if (!(pos_label %in% levels(y))) {
    stop("Positive label not found in y.")
  }
  neg_levels <- setdiff(levels(y), pos_label)
  if (length(neg_levels) != 1) {
    stop("Binary y required. Found levels: ", paste(levels(y), collapse = ", "))
  }
  neg_label <- neg_levels[[1]]
  
  if (!is.null(exclude_features) && length(exclude_features) > 0) {
    fit_cfg$features <- setdiff(
      if (is.null(fit_cfg$features)) setdiff(names(train_df), fit_cfg$y_col) else fit_cfg$features,
      exclude_features
    )
  }

  registry <- prepare_feature_registry(train_df, fit_cfg)
  split_obj <- make_selection_split(train_df, fit_cfg)
  
  feature_outputs <- lapply(seq_len(nrow(registry)), function(i) {
    fit_one_feature_marginal(
      train_df = train_df,
      feature_row = registry[i, , drop = FALSE],
      split_obj = split_obj,
      fit_cfg = fit_cfg,
      y_col = y_col,
      pos_label = pos_label,
      neg_label = neg_label
    )
  })
  
  report <- build_marginal_fit_report(feature_outputs)
  
  successful <- Filter(function(x) !is.null(x$result), feature_outputs)
  feature_models <- lapply(successful, `[[`, "result")
  names(feature_models) <- vapply(feature_models, `[[`, character(1), "feature")
  
  numeric_features <- names(Filter(function(x) identical(x$type, "numeric"), feature_models))
  categorical_features <- names(Filter(function(x) identical(x$type, "categorical"), feature_models))
  
  marginal_spec <- list(
    y_col = y_col,
    pos_label = pos_label,
    neg_label = neg_label,
    features = names(feature_models),
    feature_models = feature_models,
    numeric_features = numeric_features,
    categorical_features = categorical_features,
    selection_split = split_obj,
    fit_cfg = fit_cfg
  )
  
  new_stage_output(
    result = marginal_spec,
    diagnostics = report$diagnostics,
    warnings = report$warnings,
    dropped = report$dropped,
    metadata = list(
      registry = registry,
      candidate_diagnostics = report$candidate_diagnostics
    )
  )
}