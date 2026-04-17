fit_one_feature_marginal <- function(train_df,
                                     feature_row,
                                     split_obj,
                                     fit_cfg,
                                     y_col,
                                     pos_label,
                                     neg_label) {
  feature <- feature_row$feature
  type <- feature_row$type
  subtype <- feature_row$subtype

  y <- train_df[[y_col]]
  x_raw <- train_df[[feature]]

  warnings_df <- NULL
  dropped_df <- NULL

  selected_family_initial <- NA_character_
  selected_family_final <- NA_character_
  selected_family_pos <- NA_character_
  selected_family_neg <- NA_character_
  fallback_level <- NA_integer_
  fallback_reason <- NA_character_
  chosen_row <- NULL

  prep <- list(
    x = x_raw,
    transform_applied = "none"
  )

  resolve_family_pair <- function(family_name) {
    if (identical(family_name, "gmm2_pos_gaussian_neg")) {
      return(list(pos = "gmm2", neg = "gaussian"))
    }
    list(pos = family_name, neg = family_name)
  }

  apply_candidate_filters <- function(df, fit_cfg) {
    if (is.null(df) || nrow(df) == 0) return(df)

    keep <- rep(TRUE, nrow(df))

    if ("gmm_check_pass" %in% names(df)) {
      bad_gmm <- !is.na(df$gmm_check_pass) & !df$gmm_check_pass
      keep[bad_gmm] <- FALSE
    }

    df[keep, , drop = FALSE]
  }

  rank_candidate_diag <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)

    if ("heldout_llr_auroc" %in% names(df)) {
      df$..ord_llr_auroc <- df$heldout_llr_auroc
      df$..ord_llr_auroc[!is.finite(df$..ord_llr_auroc)] <- -Inf
    } else {
      df$..ord_llr_auroc <- -Inf
    }

    if ("heldout_trimmed_llr_gap" %in% names(df)) {
      df$..ord_llr_gap <- df$heldout_trimmed_llr_gap
      df$..ord_llr_gap[!is.finite(df$..ord_llr_gap)] <- -Inf
    } else if ("heldout_llr_gap" %in% names(df)) {
      df$..ord_llr_gap <- df$heldout_llr_gap
      df$..ord_llr_gap[!is.finite(df$..ord_llr_gap)] <- -Inf
    } else {
      df$..ord_llr_gap <- -Inf
    }

    if ("heldout_balanced_trimmed_ll" %in% names(df)) {
      df$..ord_fit_ll <- df$heldout_balanced_trimmed_ll
      df$..ord_fit_ll[!is.finite(df$..ord_fit_ll)] <- -Inf
    } else {
      df$..ord_fit_ll <- -Inf
    }

    df$..ord_simple <- ifelse(
      df$family %in% c("gaussian", "lognormal", "gamma", "poisson", "negbinom", "logit_gaussian"),
      1, 0
    )

    df <- df[order(
      df$..ord_llr_auroc,
      df$..ord_llr_gap,
      df$..ord_fit_ll,
      df$..ord_simple,
      decreasing = TRUE
    ), , drop = FALSE]

    df$..ord_llr_auroc <- NULL
    df$..ord_llr_gap <- NULL
    df$..ord_fit_ll <- NULL
    df$..ord_simple <- NULL

    df
  }

  if (type == "categorical") {
    levs <- levels(if (is.factor(x_raw)) x_raw else as.factor(x_raw))
    x_cat <- factor(x_raw, levels = levs)

    lap <- fit_cfg$laplace
    df_pos <- x_cat[y == pos_label]
    df_neg <- x_cat[y == neg_label]

    tab_pos <- table(factor(df_pos, levels = levs))
    tab_neg <- table(factor(df_neg, levels = levs))

    p_pos <- (tab_pos + lap) / (sum(tab_pos) + lap * length(levs))
    p_neg <- (tab_neg + lap) / (sum(tab_neg) + lap * length(levs))

    result <- list(
      feature = feature,
      type = type,
      subtype = subtype,
      family = "categorical",
      transform_applied = "none",
      pos = list(levels = levs, p = as.numeric(p_pos)),
      neg = list(levels = levs, p = as.numeric(p_neg))
    )

    diagnostics <- data.frame(
      feature = feature,
      type = type,
      subtype = subtype,
      selected_family_initial = "categorical",
      selected_family_final = "categorical",
      selected_family_pos = "categorical",
      selected_family_neg = "categorical",
      selected_params_pos = paste0(
        paste0(levs, "=", sprintf("%.4f", as.numeric(p_pos))), collapse = "; "
      ),
      selected_params_neg = paste0(
        paste0(levs, "=", sprintf("%.4f", as.numeric(p_neg))), collapse = "; "
      ),
      fallback_level = 0,
      fallback_reason = NA_character_,
      final_status = "ok",
      drop_flag = FALSE,
      drop_reason = NA_character_,
      stringsAsFactors = FALSE
    )

    return(new_stage_output(
      result = result,
      diagnostics = diagnostics,
      warnings = warnings_df,
      dropped = dropped_df,
      metadata = list(feature = feature)
    ))
  }

  quality <- check_raw_feature_quality(x_raw, y, fit_cfg)
  if (!quality$keep) {
    dropped_df <- new_dropped_row(
      stage = "fit_marginals",
      feature = feature,
      reason = quality$reason
    )

    diagnostics <- data.frame(
      feature = feature,
      type = type,
      subtype = subtype,
      selected_family_initial = selected_family_initial,
      selected_family_final = NA_character_,
      selected_family_pos = NA_character_,
      selected_family_neg = NA_character_,
      selected_params_pos = NA_character_,
      selected_params_neg = NA_character_,
      fallback_level = fallback_level,
      fallback_reason = fallback_reason,
      final_status = "dropped",
      drop_flag = TRUE,
      drop_reason = fallback_reason %||% "all_candidates_failed",
      transform_applied = prep$transform_applied,
      stringsAsFactors = FALSE
    )

    return(new_stage_output(
      result = NULL,
      diagnostics = diagnostics,
      warnings = warnings_df,
      dropped = dropped_df,
      metadata = list(feature = feature)
    ))
  }

  prep <- preprocess_feature(x_raw, subtype, fit_cfg)
  x <- prep$x

  candidate_families <- switch(
    subtype,
    continuous = fit_cfg$numeric_candidates,
    count = fit_cfg$count_candidates,
    fraction = fit_cfg$fraction_candidates,
    fit_cfg$numeric_candidates
  )

  candidate_diag <- evaluate_candidate_families(
    x = x,
    y = y,
    candidate_families = candidate_families,
    split_obj = split_obj,
    pos_label = pos_label,
    neg_label = neg_label,
    fit_cfg = fit_cfg,
    feature = feature,
    subtype = subtype
  )

  candidate_diag_raw <- candidate_diag
  candidate_diag <- apply_candidate_filters(candidate_diag, fit_cfg)
  candidate_diag <- rank_candidate_diag(candidate_diag)

  if (nrow(candidate_diag_raw) > 0 && nrow(candidate_diag) == 0) {
    warnings_df <- new_warning_row(
      stage = "fit_marginals",
      feature = feature,
      code = "all_candidates_filtered",
      message = "All candidate families were filtered out by candidate sanity checks"
    )
  }

  if (nrow(candidate_diag) > 0) {
    selected_family_initial <- candidate_diag$family[[1]]
    selected_family_final <- selected_family_initial
    fam_pair_initial <- resolve_family_pair(selected_family_initial)
    selected_family_pos <- fam_pair_initial$pos
    selected_family_neg <- fam_pair_initial$neg
    chosen_row <- candidate_diag[1, , drop = FALSE]
  }

  if (nrow(candidate_diag) == 0) {
    dropped_df <- new_dropped_row(
      stage = "fit_marginals",
      feature = feature,
      reason = "no_candidate_diagnostics_after_filtering"
    )

    warnings_df <- dplyr::bind_rows(
      warnings_df,
      new_warning_row(
        stage = "fit_marginals",
        feature = feature,
        code = "no_candidate_diagnostics",
        message = "No candidate family diagnostics remained after filtering"
      )
    )

    diagnostics <- data.frame(
      feature = feature,
      type = type,
      subtype = subtype,
      selected_family_initial = selected_family_initial,
      selected_family_final = selected_family_final,
      selected_family_pos = selected_family_pos,
      selected_family_neg = selected_family_neg,
      selected_params_pos = NA_character_,
      selected_params_neg = NA_character_,
      fallback_level = fallback_level,
      fallback_reason = fallback_reason,
      final_status = "dropped",
      drop_flag = TRUE,
      drop_reason = "no_candidate_diagnostics_after_filtering",
      transform_applied = prep$transform_applied,
      stringsAsFactors = FALSE
    )

    return(new_stage_output(
      result = NULL,
      diagnostics = diagnostics,
      warnings = warnings_df,
      dropped = dropped_df,
      metadata = list(feature = feature, candidate_diag = candidate_diag_raw)
    ))
  }

  selected_family_initial <- candidate_diag$family[[1]]
  selected_family_final <- selected_family_initial
  fam_pair_initial <- resolve_family_pair(selected_family_initial)
  selected_family_pos <- fam_pair_initial$pos
  selected_family_neg <- fam_pair_initial$neg
  fallback_level <- 0L
  fallback_reason <- NA_character_
  accepted <- FALSE

  chosen_row <- candidate_diag[1, , drop = FALSE]
  assessment <- assess_candidate_fit(chosen_row, fit_cfg)

  if (assessment$accept) {
    accepted <- TRUE
  } else {
    fallback_reason <- assessment$reason
    fallback_families <- setdiff(get_fallback_families(subtype, fit_cfg), selected_family_initial)

    if (length(fallback_families) > 0) {
      fb_diag <- evaluate_candidate_families(
        x = x,
        y = y,
        candidate_families = fallback_families,
        split_obj = split_obj,
        pos_label = pos_label,
        neg_label = neg_label,
        fit_cfg = fit_cfg,
        feature = feature,
        subtype = subtype
      )

      fb_diag <- apply_candidate_filters(fb_diag, fit_cfg)
      fb_diag <- rank_candidate_diag(fb_diag)

      if (nrow(fb_diag) > 0) {
        fb_row <- fb_diag[1, , drop = FALSE]
        fb_assess <- assess_candidate_fit(fb_row, fit_cfg)

        candidate_diag <- dplyr::bind_rows(candidate_diag, fb_diag)

        if (fb_assess$accept) {
          accepted <- TRUE
          selected_family_final <- fb_row$family[[1]]
          fam_pair_final <- resolve_family_pair(selected_family_final)
          selected_family_pos <- fam_pair_final$pos
          selected_family_neg <- fam_pair_final$neg
          fallback_level <- 1L
          chosen_row <- fb_row
        }
      }
    }
  }

  if (!accepted) {
    dropped_df <- new_dropped_row(
      stage = "fit_marginals",
      feature = feature,
      reason = fallback_reason %||% "all_candidates_failed"
    )

    warnings_df <- dplyr::bind_rows(
      warnings_df,
      new_warning_row(
        stage = "fit_marginals",
        feature = feature,
        code = "feature_dropped",
        message = paste("Dropped after candidate and fallback assessment:", fallback_reason %||% "all_candidates_failed")
      )
    )

    diagnostics <- data.frame(
      feature = feature,
      type = type,
      subtype = subtype,
      selected_family_initial = selected_family_initial,
      selected_family_final = NA_character_,
      selected_family_pos = selected_family_pos,
      selected_family_neg = selected_family_neg,
      selected_params_pos = NA_character_,
      selected_params_neg = NA_character_,
      fallback_level = fallback_level,
      fallback_reason = fallback_reason,
      final_status = "dropped",
      drop_flag = TRUE,
      drop_reason = fallback_reason %||% "all_candidates_failed",
      transform_applied = prep$transform_applied,
      stringsAsFactors = FALSE
    )

    return(new_stage_output(
      result = NULL,
      diagnostics = diagnostics,
      warnings = warnings_df,
      dropped = dropped_df,
      metadata = list(feature = feature, candidate_diag = candidate_diag_raw)
    ))
  }

  x_pos_full <- x[y == pos_label]
  x_neg_full <- x[y == neg_label]

  fam_pair_final <- resolve_family_pair(selected_family_final)
  selected_family_pos <- fam_pair_final$pos
  selected_family_neg <- fam_pair_final$neg

  pos_model <- fit_family_by_name(x_pos_full, selected_family_pos, fit_cfg)
  neg_model <- fit_family_by_name(x_neg_full, selected_family_neg, fit_cfg)

  if (is.null(pos_model) || is.null(neg_model)) {
    dropped_df <- new_dropped_row(
      stage = "fit_marginals",
      feature = feature,
      reason = "final_refit_failed"
    )

    warnings_df <- dplyr::bind_rows(
      warnings_df,
      new_warning_row(
        stage = "fit_marginals",
        feature = feature,
        code = "final_refit_failed",
        message = "Dropped because final full-data refit failed"
      )
    )

    diagnostics <- data.frame(
      feature = feature,
      type = type,
      subtype = subtype,
      selected_family_initial = selected_family_initial,
      selected_family_final = selected_family_final,
      selected_family_pos = selected_family_pos,
      selected_family_neg = selected_family_neg,
      selected_params_pos = NA_character_,
      selected_params_neg = NA_character_,
      fallback_level = fallback_level,
      fallback_reason = fallback_reason,
      final_status = "dropped",
      drop_flag = TRUE,
      drop_reason = "final_refit_failed",
      transform_applied = prep$transform_applied,
      stringsAsFactors = FALSE
    )

    return(new_stage_output(
      result = NULL,
      diagnostics = diagnostics,
      warnings = warnings_df,
      dropped = dropped_df,
      metadata = list(feature = feature, candidate_diag = candidate_diag_raw)
    ))
  }

  diagnostics <- cbind(
    data.frame(
      feature = feature,
      type = type,
      subtype = subtype,
      selected_family_initial = selected_family_initial,
      selected_family_final = selected_family_final,
      selected_family_pos = selected_family_pos,
      selected_family_neg = selected_family_neg,
      selected_params_pos = format_model_params(pos_model),
      selected_params_neg = format_model_params(neg_model),
      fallback_level = fallback_level,
      fallback_reason = fallback_reason,
      final_status = "ok",
      drop_flag = FALSE,
      drop_reason = NA_character_,
      transform_applied = prep$transform_applied,
      stringsAsFactors = FALSE
    ),
    chosen_row[, setdiff(names(chosen_row), c("feature", "subtype", "family")), drop = FALSE]
  )

  result <- list(
    feature = feature,
    type = type,
    subtype = subtype,
    family = selected_family_final,
    family_pos = selected_family_pos,
    family_neg = selected_family_neg,
    transform_applied = prep$transform_applied,
    pos = pos_model,
    neg = neg_model,
    candidate_diagnostics = candidate_diag_raw
  )

  new_stage_output(
    result = result,
    diagnostics = diagnostics,
    warnings = warnings_df,
    dropped = dropped_df,
    metadata = list(feature = feature)
  )
}