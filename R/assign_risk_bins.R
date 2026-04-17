assign_risk_bins <- function(stay_id,
                             score,
                             score_name = "score",
                             bin_spec = NULL,
                             use_quantiles = FALSE,
                             probs = c(0, 0.25, 0.5, 0.75, 1),
                             labels = c("Q1", "Q2", "Q3", "Q4")) {
  if (length(stay_id) != length(score)) {
    stop("stay_id and score must have the same length")
  }

  if (use_quantiles) {
    qs <- unique(stats::quantile(score, probs = probs, na.rm = TRUE))

    if (length(qs) < 3) {
      risk_bin <- rep(NA_character_, length(score))
    } else {
      n_bins <- length(qs) - 1
      labels <- labels[seq_len(n_bins)]

      risk_bin <- cut(
        score,
        breaks = qs,
        labels = labels,
        include.lowest = TRUE,
        right = TRUE
      )
    }
  } else {
    if (is.null(bin_spec)) {
      stop("bin_spec must be provided when use_quantiles = FALSE")
    }

    risk_bin <- cut(
      score,
      breaks = bin_spec$breaks,
      labels = bin_spec$labels,
      include.lowest = TRUE,
      right = TRUE
    )
  }

  out <- data.frame(
    stay_id = stay_id,
    score = score,
    risk_bin = as.character(risk_bin),
    stringsAsFactors = FALSE
  )

  names(out)[names(out) == "score"] <- score_name
  out
}