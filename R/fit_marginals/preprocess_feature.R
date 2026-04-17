winsorize_vec <- function(x, p = 0.01, na.rm = TRUE) {
  qs <- stats::quantile(x, probs = c(p, 1 - p), na.rm = na.rm, type = 7)
  pmin(pmax(x, qs[[1]]), qs[[2]])
}

preprocess_feature <- function(x, subtype, fit_cfg) {
  transform_applied <- "none"
  x_out <- x
  
  if (subtype == "continuous") {
    x_out <- winsorize_vec(x_out, p = fit_cfg$winsor_p)
    transform_applied <- "winsorized"
  }
  
  if (subtype == "fraction") {
    eps <- fit_cfg$fraction_eps
    x_out <- pmin(pmax(x_out, eps), 1 - eps)
    transform_applied <- "clipped_fraction"
  }
  
  list(
    x = x_out,
    transform_applied = transform_applied
  )
}

check_raw_feature_quality <- function(x, y, fit_cfg) {
  ok <- !is.na(y) & is.finite(x)
  x_ok <- x[ok]
  y_ok <- y[ok]
  
  out <- list(
    keep = TRUE,
    reason = NA_character_,
    n_total = length(x_ok),
    n_missing = sum(!ok),
    sd_raw = if (length(x_ok) > 1) stats::sd(x_ok) else NA_real_
  )
  
  if (length(x_ok) < fit_cfg$min_n_total) {
    out$keep <- FALSE
    out$reason <- "min_n_total_fail"
    return(out)
  }
  
  cls_counts <- table(y_ok)
  if (any(cls_counts < fit_cfg$min_n_per_class)) {
    out$keep <- FALSE
    out$reason <- "min_n_per_class_fail"
    return(out)
  }
  
  sd_thresh <- as.numeric(fit_cfg$raw_min_sd)

  if (!is.na(out$sd_raw) && out$sd_raw == 0) {
    out$keep <- FALSE
    out$reason <- "raw_zero_sd"
    return(out)
  }

  if (!is.na(out$sd_raw) && out$sd_raw < sd_thresh) {
    out$keep <- FALSE
    out$reason <- "raw_low_sd"
    return(out)
  }
  
  out
}