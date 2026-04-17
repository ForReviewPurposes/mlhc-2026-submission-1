score_group_percentiles <- function(group_pos_sum, ecdf_fit, groups_cfg) {
  keep_groups <- ecdf_fit$keep_groups
  ecdf_list <- ecdf_fit$ecdf_list

  percentile_df <- sapply(keep_groups, function(g) {
    ecdf_list[[g]](log(group_pos_sum[[g]] + groups_cfg$ecdf_log_eps))
  })

  percentile_df <- as.data.frame(percentile_df, check.names = FALSE)

  percentile_norm <- percentile_df / rowSums(percentile_df)

  diagnostics <- data.frame(
    group = keep_groups,
    mean_percentile = vapply(percentile_df, mean, numeric(1), na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  new_stage_output(
    result = list(
      group_percentiles = percentile_df,
      group_percentiles_norm = percentile_norm
    ),
    diagnostics = diagnostics,
    warnings = NULL,
    dropped = NULL,
    metadata = NULL
  )
}