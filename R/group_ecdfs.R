fit_group_ecdfs <- function(group_pos_sum, groups_cfg) {
  exclude_groups <- groups_cfg$exclude_from_percentiles %||% character(0)
  keep_groups <- setdiff(colnames(group_pos_sum), exclude_groups)

  ecdf_list <- lapply(keep_groups, function(g) {
    x <- group_pos_sum[[g]]
    stats::ecdf(log(x + groups_cfg$ecdf_log_eps))
  })
  names(ecdf_list) <- keep_groups

  diagnostics <- data.frame(
    group = keep_groups,
    excluded = FALSE,
    stringsAsFactors = FALSE
  )

  if (length(exclude_groups) > 0) {
    diagnostics <- dplyr::bind_rows(
      diagnostics,
      data.frame(
        group = intersect(exclude_groups, colnames(group_pos_sum)),
        excluded = TRUE,
        stringsAsFactors = FALSE
      )
    )
  }

  new_stage_output(
    result = list(
      ecdf_list = ecdf_list,
      keep_groups = keep_groups,
      exclude_groups = exclude_groups
    ),
    diagnostics = diagnostics,
    warnings = NULL,
    dropped = NULL,
    metadata = NULL
  )
}

plot_group_ecdfs <- function(group_pos_df, out_dir) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # long format: patient x group
  df_long <- group_pos_df %>%
    pivot_longer(
      cols = everything(),
      names_to = "group",
      values_to = "value"
    )

  # loop per group (cleaner plots)
  for (g in unique(df_long$group)) {

    df_g <- df_long %>% filter(group == g)

    p <- ggplot(df_g, aes(x = value)) +
      stat_ecdf(size = 1.2) +
      labs(
        title = paste("ECDF:", g),
        x = "Group positive evidence",
        y = "ECDF"
      ) +
      plot_theme_mlhc()

    ggsave(
      filename = file.path(out_dir, paste0("ecdf_", g, ".png")),
      plot = p,
      width = 6,
      height = 4,
      dpi = 300
    )
  }
}

compute_ecdf_summary <- function(group_pos_df, out_path) {
  library(dplyr)
  library(tidyr)

  df_long <- group_pos_df %>%
    pivot_longer(
      cols = everything(),
      names_to = "group",
      values_to = "value"
    )

  summary <- df_long %>%
    group_by(group) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      p75 = quantile(value, 0.75, na.rm = TRUE),
      p90 = quantile(value, 0.90, na.rm = TRUE),
      p95 = quantile(value, 0.95, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      frac_gt_p75 = mean(value > quantile(value, 0.75, na.rm = TRUE)),
      frac_gt_p90 = mean(value > quantile(value, 0.90, na.rm = TRUE)),
      tail_ratio = p95 / mean,
      tail_spread = p95 - median
    )

  write.csv(summary, out_path, row.names = FALSE)

  return(summary)
}