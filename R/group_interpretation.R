build_interpretation_strings <- function(patient_summary, groups_cfg) {
  base_cols <- c(
  "stay_id", "sum_llr", "d_dist", "risk_bin_primary", "risk_bin_ddist",
  "weighted_intra_entropy", "inter_group_entropy",
  "interpretation", "atypical_group_count", "top_group_name",
  "top_group_percentile", "atypical_groups_csv"
)

  group_cols <- setdiff(names(patient_summary), base_cols)
  threshold <- groups_cfg$atypical_threshold

  rows <- apply(patient_summary, 1, function(row) {
    vals <- suppressWarnings(as.numeric(row[group_cols]))
    names(vals) <- group_cols

    vals <- vals[is.finite(vals)]
    vals <- sort(vals, decreasing = TRUE)

    atypical_vals <- vals[vals >= threshold]

    atypical_count <- length(atypical_vals)

    top_group_name <- if (length(vals) > 0) names(vals)[1] else NA_character_
    top_group_percentile <- if (length(vals) > 0) as.numeric(vals[1]) else NA_real_

    atypical_csv <- if (length(atypical_vals) == 0) {
      NA_character_
    } else {
      paste0(names(atypical_vals), " (", sprintf("%.2f", atypical_vals), ")", collapse = ", ")
    }

    interpretation <- if (length(atypical_vals) == 0) {
      paste0(
        "Evidence score: ", sprintf("%.3f", as.numeric(row[["sum_llr"]])),
        ", Risk bin: ", row[["risk_bin_primary"]],
        ", Atypical groups: none above threshold"
      )
    } else {
      paste0(
        "Evidence score: ", sprintf("%.3f", as.numeric(row[["sum_llr"]])),
        ", Risk bin: ", row[["risk_bin_primary"]],
        ", Atypical groups: ", atypical_csv
      )
    }

    data.frame(
      stay_id = row[["stay_id"]],
      atypical_group_count = atypical_count,
      top_group_name = top_group_name,
      top_group_percentile = top_group_percentile,
      atypical_groups_csv = atypical_csv,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(rows)
}