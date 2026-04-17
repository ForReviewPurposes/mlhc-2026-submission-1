summarize_group_percentile_transport <- function(patient_summary,
                                                 atypical_threshold = 0.75,
                                                 id_col = "stay_id",
                                                 primary_bin_col = "risk_bin_primary") {
  base_cols <- c(
    id_col,
    "sum_llr", "d_dist",
    "risk_bin", "risk_bin_primary", "risk_bin_ddist",
    "weighted_intra_entropy", "inter_group_entropy",
    "interpretation",
    "atypical_group_count",
    "top_group_name",
    "top_group_percentile",
    "atypical_groups_csv"
  )

  group_cols <- setdiff(names(patient_summary), base_cols)
  group_cols <- group_cols[vapply(patient_summary[group_cols], is.numeric, logical(1))]

  if (length(group_cols) == 0) {
    stop("No numeric group percentile columns found in patient_summary.")
  }

  pct_mat <- as.data.frame(patient_summary[, group_cols, drop = FALSE])

  # per-patient atypical count
  atypical_count <- rowSums(pct_mat >= atypical_threshold, na.rm = TRUE)

  patient_count_tbl <- data.frame(
    stay_id = patient_summary[[id_col]],
    risk_bin_primary = if (primary_bin_col %in% names(patient_summary)) patient_summary[[primary_bin_col]] else NA,
    atypical_group_count = atypical_count,
    stringsAsFactors = FALSE
  )

  patient_count_summary <- patient_count_tbl %>%
    dplyr::group_by(risk_bin_primary) %>%
    dplyr::summarise(
      n_patients = dplyr::n(),
      mean_atypical_group_count = mean(atypical_group_count, na.rm = TRUE),
      median_atypical_group_count = stats::median(atypical_group_count, na.rm = TRUE),
      q25_atypical_group_count = as.numeric(stats::quantile(atypical_group_count, 0.25, na.rm = TRUE)),
      q75_atypical_group_count = as.numeric(stats::quantile(atypical_group_count, 0.75, na.rm = TRUE)),
      .groups = "drop"
    )

  # per-group atypical frequency across all patients
  group_freq_tbl <- data.frame(
    group = group_cols,
    atypical_n = vapply(group_cols, function(g) sum(patient_summary[[g]] >= atypical_threshold, na.rm = TRUE), numeric(1)),
    atypical_frac = vapply(group_cols, function(g) mean(patient_summary[[g]] >= atypical_threshold, na.rm = TRUE), numeric(1)),
    mean_percentile = vapply(group_cols, function(g) mean(patient_summary[[g]], na.rm = TRUE), numeric(1)),
    median_percentile = vapply(group_cols, function(g) stats::median(patient_summary[[g]], na.rm = TRUE), numeric(1)),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(dplyr::desc(atypical_frac), dplyr::desc(mean_percentile))

  # per-group atypical frequency by primary bin
  if (primary_bin_col %in% names(patient_summary)) {
    group_freq_by_bin <- dplyr::bind_rows(lapply(group_cols, function(g) {
      patient_summary %>%
        dplyr::group_by(.data[[primary_bin_col]]) %>%
        dplyr::summarise(
          group = g,
          atypical_n = sum(.data[[g]] >= atypical_threshold, na.rm = TRUE),
          atypical_frac = mean(.data[[g]] >= atypical_threshold, na.rm = TRUE),
          mean_percentile = mean(.data[[g]], na.rm = TRUE),
          median_percentile = stats::median(.data[[g]], na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::rename(risk_bin_primary = !!primary_bin_col)
    }))
  } else {
    group_freq_by_bin <- NULL
  }

  # primary-bin risk summary
  if (primary_bin_col %in% names(patient_summary)) {
    risk_bin_summary_primary <- patient_summary %>%
      dplyr::group_by(.data[[primary_bin_col]]) %>%
      dplyr::summarise(
        n_patients = dplyr::n(),
        mean_sum_llr = mean(sum_llr, na.rm = TRUE),
        median_sum_llr = stats::median(sum_llr, na.rm = TRUE),
        mean_d_dist = mean(d_dist, na.rm = TRUE),
        mean_weighted_intra_entropy = mean(weighted_intra_entropy, na.rm = TRUE),
        mean_inter_group_entropy = mean(inter_group_entropy, na.rm = TRUE),
        mean_atypical_group_count = mean(atypical_count, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::rename(risk_bin_primary = !!primary_bin_col)
  } else {
    risk_bin_summary_primary <- NULL
  }

  list(
    patient_count_table = patient_count_tbl,
    patient_count_summary = patient_count_summary,
    group_atypical_frequency = group_freq_tbl,
    group_atypical_frequency_by_bin = group_freq_by_bin,
    risk_bin_summary_primary = risk_bin_summary_primary
  )
}

summarize_group_positive_sum_distribution <- function(group_pos_sum,
                                                      risk_bin_primary = NULL) {
  if (!is.data.frame(group_pos_sum)) {
    group_pos_sum <- as.data.frame(group_pos_sum)
  }

  group_cols <- names(group_pos_sum)
  group_cols <- group_cols[vapply(group_pos_sum[group_cols], is.numeric, logical(1))]

  if (length(group_cols) == 0) {
    stop("No numeric group positive-sum columns found.")
  }

  overall_tbl <- dplyr::bind_rows(lapply(group_cols, function(g) {
    x <- group_pos_sum[[g]]
    data.frame(
      group = g,
      n_nonmissing = sum(is.finite(x)),
      mean_pos_sum = mean(x, na.rm = TRUE),
      sd_pos_sum = stats::sd(x, na.rm = TRUE),
      q25_pos_sum = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE)),
      median_pos_sum = stats::median(x, na.rm = TRUE),
      q75_pos_sum = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE)),
      q90_pos_sum = as.numeric(stats::quantile(x, 0.90, na.rm = TRUE)),
      q95_pos_sum = as.numeric(stats::quantile(x, 0.95, na.rm = TRUE)),
      max_pos_sum = max(x, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })) %>%
    dplyr::arrange(dplyr::desc(mean_pos_sum))

  by_bin_tbl <- NULL

  if (!is.null(risk_bin_primary)) {
    tmp <- group_pos_sum
    tmp$risk_bin_primary <- risk_bin_primary

    by_bin_tbl <- dplyr::bind_rows(lapply(group_cols, function(g) {
      tmp %>%
        dplyr::group_by(risk_bin_primary) %>%
        dplyr::summarise(
          group = g,
          n_nonmissing = sum(is.finite(.data[[g]])),
          mean_pos_sum = mean(.data[[g]], na.rm = TRUE),
          sd_pos_sum = stats::sd(.data[[g]], na.rm = TRUE),
          q25_pos_sum = as.numeric(stats::quantile(.data[[g]], 0.25, na.rm = TRUE)),
          median_pos_sum = stats::median(.data[[g]], na.rm = TRUE),
          q75_pos_sum = as.numeric(stats::quantile(.data[[g]], 0.75, na.rm = TRUE)),
          q90_pos_sum = as.numeric(stats::quantile(.data[[g]], 0.90, na.rm = TRUE)),
          .groups = "drop"
        )
    }))
  }

  list(
    group_pos_sum_summary = overall_tbl,
    group_pos_sum_summary_by_bin = by_bin_tbl
  )
}