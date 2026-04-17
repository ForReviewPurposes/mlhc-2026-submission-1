summarize_risk_bins <- function(risk_tbl, outcome = NULL, score_col = NULL) {
  if (!("risk_bin" %in% names(risk_tbl))) {
    stop("risk_tbl must contain risk_bin")
  }

  if (is.null(score_col)) {
    score_candidates <- setdiff(names(risk_tbl), c("stay_id", "risk_bin"))
    if (length(score_candidates) != 1) {
      stop("Could not infer score column. Please provide score_col explicitly.")
    }
    score_col <- score_candidates[[1]]
  }

  out <- risk_tbl %>%
    dplyr::group_by(risk_bin) %>%
    dplyr::summarise(
      n = dplyr::n(),
      pct = n / nrow(risk_tbl),
      mean_score = mean(.data[[score_col]], na.rm = TRUE),
      median_score = stats::median(.data[[score_col]], na.rm = TRUE),
      .groups = "drop"
    )

  names(out)[names(out) == "mean_score"] <- paste0("mean_", score_col)
  names(out)[names(out) == "median_score"] <- paste0("median_", score_col)

  if (!is.null(outcome)) {
    if (length(outcome) != nrow(risk_tbl)) {
      stop("Outcome length must match number of rows in risk_tbl")
    }

    tmp <- risk_tbl
    tmp$outcome <- outcome

    out2 <- tmp %>%
      dplyr::group_by(risk_bin) %>%
      dplyr::summarise(
        event_n = sum(outcome == 1, na.rm = TRUE),
        event_rate = mean(outcome == 1, na.rm = TRUE),
        .groups = "drop"
      )

    out <- dplyr::left_join(out, out2, by = "risk_bin")
  }

  out
}