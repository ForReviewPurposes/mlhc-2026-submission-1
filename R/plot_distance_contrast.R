plot_theme_mlhc <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = base_size + 2, face = "bold"),
      axis.title = ggplot2::element_text(size = base_size, face = "bold"),
      axis.text = ggplot2::element_text(size = base_size - 1),
      legend.title = ggplot2::element_text(size = base_size - 1, face = "bold"),
      legend.text = ggplot2::element_text(size = base_size - 1),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.35),
      strip.text = ggplot2::element_text(size = base_size, face = "bold")
    )
}

plot_distance_contrast <- function(score_df, outcome = NULL, plots_cfg = NULL) {
  required_cols <- c("sum_llr", "d_dist")
  missing_cols <- setdiff(required_cols, names(score_df))
  if (length(missing_cols) > 0) {
    stop("score_df missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  plot_df <- score_df
  
  if (!is.null(outcome)) {
    if (length(outcome) != nrow(plot_df)) {
      stop("Outcome length must match score_df rows")
    }
    plot_df$outcome <- as.factor(outcome)
  } else {
    plot_df$outcome <- "all"
  }
  
  if (!is.null(outcome)) {
    p_sum_llr <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sum_llr, color = outcome, fill = outcome)) +
      ggplot2::geom_density(alpha = 0.2, linewidth = 1.1) +
      ggplot2::labs(title = "Density of sum_llr", x = "sum_llr", y = "Density") +
      plot_theme_mlhc()
    
    p_d_dist <- ggplot2::ggplot(plot_df, ggplot2::aes(x = d_dist, color = outcome, fill = outcome)) +
      ggplot2::geom_density(alpha = 0.2, linewidth = 1.1) +
      ggplot2::labs(title = "Density of d_dist", x = "d_dist", y = "Density") +
      plot_theme_mlhc()
  } else {
    p_sum_llr <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sum_llr)) +
      ggplot2::geom_density(linewidth = 1.1) +
      ggplot2::labs(title = "Density of sum_llr", x = "sum_llr", y = "Density") +
      plot_theme_mlhc()
    
    p_d_dist <- ggplot2::ggplot(plot_df, ggplot2::aes(x = d_dist)) +
      ggplot2::geom_density(linewidth = 1.1) +
      ggplot2::labs(title = "Density of d_dist", x = "d_dist", y = "Density") +
      plot_theme_mlhc()
  }
  
  list(
    sum_llr_plot = p_sum_llr,
    d_dist_plot = p_d_dist
  )
}