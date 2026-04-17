evidence_transform <- function(new_df, marginal_fit, transform_cfg) {
  features <- marginal_fit$features
  feature_models <- marginal_fit$feature_models
  
  missing_features <- setdiff(features, names(new_df))
  if (length(missing_features) > 0) {
    stop("Missing features in new data: ", paste(missing_features, collapse = ", "))
  }
  
  n <- nrow(new_df)
  d <- length(features)
  
  L <- matrix(NA_real_, nrow = n, ncol = d, dimnames = list(NULL, features))
  
  diagnostics <- list()
  dropped <- NULL
  warnings <- NULL
  
  for (j in seq_along(features)) {
    fname <- features[[j]]
    fobj <- feature_models[[fname]]
    x <- new_df[[fname]]
    
    if (fobj$type == "categorical") {
      levs <- fobj$pos$levels
      idx <- match(as.character(x), levs)
      
      p_pos <- rep(1 / length(levs), n)
      p_neg <- rep(1 / length(levs), n)
      
      ok <- !is.na(idx)
      p_pos[ok] <- fobj$pos$p[idx[ok]]
      p_neg[ok] <- fobj$neg$p[idx[ok]]
      
      ll_pos <- log(pmax(p_pos, 1e-12))
      ll_neg <- log(pmax(p_neg, 1e-12))
      
    } else {
      x_use <- x
      
      if (!is.null(fobj$transform_applied)) {
        if (fobj$transform_applied == "clipped_fraction") {
          eps <- marginal_fit$fit_cfg$fraction_eps
          x_use <- pmin(pmax(x_use, eps), 1 - eps)
        }
      }
      
      ll_pos <- fobj$pos$logpdf(x_use)
      ll_neg <- fobj$neg$logpdf(x_use)
      
      ll_pos[!is.finite(ll_pos)] <- log(1e-12)
      ll_neg[!is.finite(ll_neg)] <- log(1e-12)
    }
    
    llr <- ll_pos - ll_neg
    
    if (isTRUE(transform_cfg$cap_llr)) {
      llr[llr > transform_cfg$llr_cap] <- transform_cfg$llr_cap
      llr[llr < -transform_cfg$llr_cap] <- -transform_cfg$llr_cap
    }
    
    L[, j] <- llr
    
    diagnostics[[j]] <- data.frame(
      feature = fname,
      n = length(llr),
      n_nonfinite = sum(!is.finite(llr)),
      sd_L = stats::sd(llr, na.rm = TRUE),
      mean_L = mean(llr, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  diagnostics_df <- dplyr::bind_rows(diagnostics)
  
  if (isTRUE(transform_cfg$drop_nonfinite_L_cols) || !is.null(transform_cfg$min_L_sd)) {
    bad_cols <- diagnostics_df$feature[
      is.na(diagnostics_df$sd_L) |
        diagnostics_df$sd_L < transform_cfg$min_L_sd |
        diagnostics_df$n_nonfinite > 0
    ]
    
    if (length(bad_cols) > 0) {
      dropped <- dplyr::bind_rows(lapply(bad_cols, function(f) {
        new_dropped_row(
          stage = "evidence_transform",
          feature = f,
          reason = ifelse(
            diagnostics_df$sd_L[diagnostics_df$feature == f] < transform_cfg$min_L_sd,
            "low_L_sd_or_na",
            "nonfinite_llr_values"
          )
        )
      }))
      
      warnings <- dplyr::bind_rows(lapply(bad_cols, function(f) {
        new_warning_row(
          stage = "evidence_transform",
          feature = f,
          code = "L_column_dropped",
          message = "Dropped transformed feature due to low SD or non-finite LLR values"
        )
      }))
      
      keep_cols <- setdiff(colnames(L), bad_cols)
      L <- L[, keep_cols, drop = FALSE]
      diagnostics_df <- diagnostics_df[diagnostics_df$feature %in% keep_cols, , drop = FALSE]
    }
  }
  
  sum_llr <- rowSums(L)
  
  result <- list(
    L = L,
    sum_llr = sum_llr
  )
  
  new_stage_output(
    result = result,
    diagnostics = diagnostics_df,
    warnings = warnings,
    dropped = dropped,
    metadata = list(
      n_features_in = length(features),
      n_features_out = ncol(L)
    )
  )
}