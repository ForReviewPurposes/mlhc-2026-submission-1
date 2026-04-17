init_run_paths <- function(cfg, run_type = "internal_validation") {
  run_dir <- here::here(cfg$output$base_dir, cfg$output$run_name)

  subdirs <- c(
    "logs",
    "model",
    "scores",
    "plots",
    "reports"
  )

  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(lapply(file.path(run_dir, subdirs), dir.create, recursive = TRUE, showWarnings = FALSE))

  list(
    run_dir = run_dir,
    logs_dir = file.path(run_dir, "logs"),
    model_dir = file.path(run_dir, "model"),
    scores_dir = file.path(run_dir, "scores"),
    plots_dir = file.path(run_dir, "plots"),
    reports_dir = file.path(run_dir, "reports"),

    log_file = file.path(run_dir, "logs", "run.log"),
    warnings_file = file.path(run_dir, "logs", "warnings.txt"),
    warnings_csv = file.path(run_dir, "logs", "warnings.csv"),
    dropped_file = file.path(run_dir, "logs", "dropped_features.csv"),

    config_snapshot = file.path(run_dir, "reports", "config_snapshot.yml"),
    run_rds = file.path(run_dir, "reports", "run_output.rds"),

    marginal_model_rds = file.path(run_dir, "model", "marginal_fit.rds"),
    geometry_model_rds = file.path(run_dir, "model", "geometry_fit.rds"),

    L_val_csv = file.path(run_dir, "scores", "L_val.csv"),
    score_table_csv = file.path(run_dir, "scores", "score_table.csv"),
    risk_table_csv = file.path(run_dir, "scores", "risk_table.csv"),
    risk_summary_csv = file.path(run_dir, "reports", "risk_bin_summary.csv"),

    sum_llr_density_png = file.path(run_dir, "plots", "sum_llr_density.png"),
    d_dist_density_png = file.path(run_dir, "plots", "d_dist_density.png")
  )
}