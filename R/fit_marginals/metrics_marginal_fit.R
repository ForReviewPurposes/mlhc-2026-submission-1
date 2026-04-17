trimmed_mean_safe <- function(x, trim = 0.01) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x, trim = trim)
}

.extract_gmm2_checks <- function(model_obj, x_train = NULL) {
  if (is.null(model_obj) || is.null(model_obj$family) || !identical(model_obj$family, "gmm2")) {
    return(list(
      is_gmm2 = FALSE,
      min_weight = NA_real_,
      max_weight = NA_real_,
      min_sd = NA_real_,
      max_sd = NA_real_,
      sd_ratio = NA_real_,
      mean_sep = NA_real_,
      sep_over_pooled_sd = NA_real_,
      boundary_spike = NA,
      bad_numeric = FALSE,
      degenerate = FALSE
    ))
  }

  p <- model_obj$params

  pi1 <- as.numeric(p$pi1)
  pi2 <- as.numeric(p$pi2)
  sd1 <- as.numeric(p$sd1)
  sd2 <- as.numeric(p$sd2)
  m1 <- as.numeric(p$mean1)
  m2 <- as.numeric(p$mean2)

  if (m1 > m2) {
    tmp <- pi1; pi1 <- pi2; pi2 <- tmp
    tmp <- sd1; sd1 <- sd2; sd2 <- tmp
    tmp <- m1;  m1  <- m2;  m2  <- tmp
  }

  pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)

  bad_numeric <- any(!is.finite(c(pi1, pi2, sd1, sd2, m1, m2))) ||
    any(c(pi1, pi2) <= 0) ||
    abs((pi1 + pi2) - 1) > 1e-4 ||
    any(c(sd1, sd2) <= 0)

  boundary_spike <- NA
  if (!is.null(x_train) && length(x_train) > 1 && all(is.finite(x_train))) {
    x_min <- min(x_train, na.rm = TRUE)
    x_max <- max(x_train, na.rm = TRUE)
    edge1 <- min(abs(m1 - x_min), abs(x_max - m1))
    edge2 <- min(abs(m2 - x_min), abs(x_max - m2))
    boundary_spike <- (edge1 < 2 * sd1) || (edge2 < 2 * sd2)
  }

  list(
    is_gmm2 = TRUE,
    min_weight = min(pi1, pi2),
    max_weight = max(pi1, pi2),
    min_sd = min(sd1, sd2),
    max_sd = max(sd1, sd2),
    sd_ratio = if (min(sd1, sd2) > 0) max(sd1, sd2) / min(sd1, sd2) else Inf,
    mean_sep = abs(m2 - m1),
    sep_over_pooled_sd = if (is.finite(pooled_sd) && pooled_sd > 0) abs(m2 - m1) / pooled_sd else NA_real_,
    boundary_spike = boundary_spike,
    bad_numeric = bad_numeric,
    degenerate = bad_numeric
  )
}

compute_marginal_fit_metrics <- function(model_pair,
                                         x_train,
                                         y_train,
                                         x_valid,
                                         y_valid,
                                         pos_label,
                                         neg_label,
                                         family,
                                         fit_cfg) {
  if (is.null(model_pair$pos) || is.null(model_pair$neg)) {
    return(data.frame(
      fit_success_pos = !is.null(model_pair$pos),
      fit_success_neg = !is.null(model_pair$neg),
      train_n_pos = sum(y_train == pos_label, na.rm = TRUE),
      train_n_neg = sum(y_train == neg_label, na.rm = TRUE),
      holdout_n_pos = sum(y_valid == pos_label, na.rm = TRUE),
      holdout_n_neg = sum(y_valid == neg_label, na.rm = TRUE),
      heldout_mean_ll_pos = NA_real_,
      heldout_mean_ll_neg = NA_real_,
      heldout_balanced_mean_ll = NA_real_,
      heldout_balanced_trimmed_ll = NA_real_,
      finite_logpdf_frac_pos = NA_real_,
      finite_logpdf_frac_neg = NA_real_,
      finite_llr_frac = NA_real_,
      p99_abs_llr = NA_real_,
      heldout_mean_llr_pos = NA_real_,
      heldout_mean_llr_neg = NA_real_,
      heldout_llr_gap = NA_real_,
      heldout_trimmed_llr_pos = NA_real_,
      heldout_trimmed_llr_neg = NA_real_,
      heldout_trimmed_llr_gap = NA_real_,
      heldout_llr_auroc = NA_real_,
      aic_pos = NA_real_,
      aic_neg = NA_real_,
      bic_pos = NA_real_,
      bic_neg = NA_real_,
      gmm_pos_is_gmm2 = NA,
      gmm_pos_min_weight = NA_real_,
      gmm_pos_min_sd = NA_real_,
      gmm_pos_mean_sep = NA_real_,
      gmm_pos_sep_over_pooled_sd = NA_real_,
      gmm_neg_is_gmm2 = NA,
      gmm_neg_min_weight = NA_real_,
      gmm_neg_min_sd = NA_real_,
      gmm_neg_mean_sep = NA_real_,
      gmm_neg_sep_over_pooled_sd = NA_real_,
      gmm_check_pass = NA,
      stringsAsFactors = FALSE
    ))
  }

  x_pos_valid <- x_valid[y_valid == pos_label]
  x_neg_valid <- x_valid[y_valid == neg_label]

  ll_pos <- model_pair$pos$logpdf(x_pos_valid)
  ll_neg <- model_pair$neg$logpdf(x_neg_valid)

  finite_logpdf_frac_pos <- mean(is.finite(ll_pos))
  finite_logpdf_frac_neg <- mean(is.finite(ll_neg))

  heldout_mean_ll_pos <- if (any(is.finite(ll_pos))) mean(ll_pos[is.finite(ll_pos)]) else NA_real_
  heldout_mean_ll_neg <- if (any(is.finite(ll_neg))) mean(ll_neg[is.finite(ll_neg)]) else NA_real_

  heldout_balanced_mean_ll <- mean(c(heldout_mean_ll_pos, heldout_mean_ll_neg), na.rm = TRUE)
  heldout_balanced_trimmed_ll <- mean(c(
    trimmed_mean_safe(ll_pos, trim = fit_cfg$trim_frac),
    trimmed_mean_safe(ll_neg, trim = fit_cfg$trim_frac)
  ), na.rm = TRUE)

  llr_valid <- model_pair$pos$logpdf(x_valid) - model_pair$neg$logpdf(x_valid)

  finite_llr_frac <- mean(is.finite(llr_valid))
  p99_abs_llr <- if (any(is.finite(llr_valid))) {
    as.numeric(stats::quantile(abs(llr_valid[is.finite(llr_valid)]), 0.99, na.rm = TRUE))
  } else {
    NA_real_
  }

  llr_valid_pos <- llr_valid[y_valid == pos_label]
  llr_valid_neg <- llr_valid[y_valid == neg_label]

  heldout_mean_llr_pos <- if (any(is.finite(llr_valid_pos))) {
    mean(llr_valid_pos[is.finite(llr_valid_pos)])
  } else {
    NA_real_
  }

  heldout_mean_llr_neg <- if (any(is.finite(llr_valid_neg))) {
    mean(llr_valid_neg[is.finite(llr_valid_neg)])
  } else {
    NA_real_
  }

  heldout_trimmed_llr_pos <- trimmed_mean_safe(llr_valid_pos, trim = fit_cfg$trim_frac)
  heldout_trimmed_llr_neg <- trimmed_mean_safe(llr_valid_neg, trim = fit_cfg$trim_frac)

  heldout_llr_gap <- heldout_mean_llr_pos - heldout_mean_llr_neg
  heldout_trimmed_llr_gap <- heldout_trimmed_llr_pos - heldout_trimmed_llr_neg

  heldout_llr_auroc <- NA_real_
  ok_llr <- is.finite(llr_valid) & !is.na(y_valid)
  y_llr <- as.integer(as.character(y_valid[ok_llr]) == pos_label)

  if (length(unique(y_llr)) >= 2) {
    roc_obj_llr <- pROC::roc(
      response = y_llr,
      predictor = llr_valid[ok_llr],
      quiet = TRUE,
      direction = "<"
    )
    heldout_llr_auroc <- as.numeric(pROC::auc(roc_obj_llr))
  }

  x_pos_train <- x_train[y_train == pos_label]
  x_neg_train <- x_train[y_train == neg_label]

  train_ll_pos <- model_pair$pos$logpdf(x_pos_train)
  train_ll_neg <- model_pair$neg$logpdf(x_neg_train)

  k_pos <- length(model_pair$pos$params)
  k_neg <- length(model_pair$neg$params)
  n_pos <- sum(is.finite(train_ll_pos))
  n_neg <- sum(is.finite(train_ll_neg))

  sum_ll_pos <- sum(train_ll_pos[is.finite(train_ll_pos)])
  sum_ll_neg <- sum(train_ll_neg[is.finite(train_ll_neg)])

  aic_pos <- -2 * sum_ll_pos + 2 * k_pos
  aic_neg <- -2 * sum_ll_neg + 2 * k_neg

  bic_pos <- -2 * sum_ll_pos + log(max(n_pos, 1)) * k_pos
  bic_neg <- -2 * sum_ll_neg + log(max(n_neg, 1)) * k_neg

  # ---- GMM sanity checks ----
  pos_chk <- .extract_gmm2_checks(model_pair$pos)
  neg_chk <- .extract_gmm2_checks(model_pair$neg)

  gmm_check_pass <- TRUE

  if (isTRUE(pos_chk$is_gmm2)) {
  if (isTRUE(pos_chk$degenerate)) gmm_check_pass <- FALSE
  if (!is.na(pos_chk$min_weight) && pos_chk$min_weight < fit_cfg$gmm2_min_weight) gmm_check_pass <- FALSE
  if (!is.na(pos_chk$min_sd) && pos_chk$min_sd < fit_cfg$gmm2_min_sd) gmm_check_pass <- FALSE
  if (!is.na(pos_chk$sep_over_pooled_sd) && pos_chk$sep_over_pooled_sd < fit_cfg$gmm2_min_sep_over_pooled_sd) gmm_check_pass <- FALSE
  if (!is.na(pos_chk$sd_ratio) && pos_chk$sd_ratio > fit_cfg$gmm2_max_sd_ratio) gmm_check_pass <- FALSE
  if (isTRUE(pos_chk$boundary_spike)) gmm_check_pass <- FALSE
}

  if (isTRUE(neg_chk$is_gmm2)) {
    if (isTRUE(neg_chk$degenerate)) gmm_check_pass <- FALSE
    if (!is.na(neg_chk$min_weight) && neg_chk$min_weight < fit_cfg$gmm2_min_weight) gmm_check_pass <- FALSE
    if (!is.na(neg_chk$min_sd) && neg_chk$min_sd < fit_cfg$gmm2_min_sd) gmm_check_pass <- FALSE
    if (!is.na(neg_chk$sep_over_pooled_sd) && neg_chk$sep_over_pooled_sd < fit_cfg$gmm2_min_sep_over_pooled_sd) gmm_check_pass <- FALSE
    if (!is.na(neg_chk$sd_ratio) && neg_chk$sd_ratio > fit_cfg$gmm2_max_sd_ratio) gmm_check_pass <- FALSE
    if (isTRUE(neg_chk$boundary_spike)) gmm_check_pass <- FALSE
  }

  data.frame(
    fit_success_pos = TRUE,
    fit_success_neg = TRUE,
    train_n_pos = sum(y_train == pos_label, na.rm = TRUE),
    train_n_neg = sum(y_train == neg_label, na.rm = TRUE),
    holdout_n_pos = sum(y_valid == pos_label, na.rm = TRUE),
    holdout_n_neg = sum(y_valid == neg_label, na.rm = TRUE),
    heldout_mean_ll_pos = heldout_mean_ll_pos,
    heldout_mean_ll_neg = heldout_mean_ll_neg,
    heldout_balanced_mean_ll = heldout_balanced_mean_ll,
    heldout_balanced_trimmed_ll = heldout_balanced_trimmed_ll,
    finite_logpdf_frac_pos = finite_logpdf_frac_pos,
    finite_logpdf_frac_neg = finite_logpdf_frac_neg,
    finite_llr_frac = finite_llr_frac,
    p99_abs_llr = p99_abs_llr,
    heldout_mean_llr_pos = heldout_mean_llr_pos,
    heldout_mean_llr_neg = heldout_mean_llr_neg,
    heldout_llr_gap = heldout_llr_gap,
    heldout_trimmed_llr_pos = heldout_trimmed_llr_pos,
    heldout_trimmed_llr_neg = heldout_trimmed_llr_neg,
    heldout_trimmed_llr_gap = heldout_trimmed_llr_gap,
    heldout_llr_auroc = heldout_llr_auroc,
    aic_pos = aic_pos,
    aic_neg = aic_neg,
    bic_pos = bic_pos,
    bic_neg = bic_neg,
    gmm_pos_is_gmm2 = pos_chk$is_gmm2,
    gmm_pos_min_weight = pos_chk$min_weight,
    gmm_pos_min_sd = pos_chk$min_sd,
    gmm_pos_mean_sep = pos_chk$mean_sep,
    gmm_pos_sep_over_pooled_sd = pos_chk$sep_over_pooled_sd,
    gmm_neg_is_gmm2 = neg_chk$is_gmm2,
    gmm_neg_min_weight = neg_chk$min_weight,
    gmm_neg_min_sd = neg_chk$min_sd,
    gmm_neg_mean_sep = neg_chk$mean_sep,
    gmm_neg_sep_over_pooled_sd = neg_chk$sep_over_pooled_sd,
    gmm_check_pass = gmm_check_pass,
    stringsAsFactors = FALSE
  )
}