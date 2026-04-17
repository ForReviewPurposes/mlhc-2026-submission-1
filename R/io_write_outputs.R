read_input_data <- function(path) {
  if (!file.exists(path)) {
    stop("Input data file not found: ", path)
  }

  ext <- tools::file_ext(path)

  if (ext == "csv") {
    return(utils::read.csv(path, stringsAsFactors = FALSE))
  }

  if (ext == "rds") {
    return(readRDS(path))
  }

  stop("Unsupported input file type: ", ext)
}

write_stage_outputs <- function(stage_name, stage_output, paths) {
  stage_dir <- file.path(paths$reports_dir, stage_name)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(stage_output$diagnostics)) {
    utils::write.csv(
      stage_output$diagnostics,
      file = file.path(stage_dir, "diagnostics.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(stage_output$dropped)) {
    utils::write.csv(
      stage_output$dropped,
      file = file.path(stage_dir, "dropped.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(stage_output$warnings)) {
    utils::write.csv(
      stage_output$warnings,
      file = file.path(stage_dir, "warnings.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(stage_output$metadata)) {
    # candidate diagnostics from fit_marginals
    if (!is.null(stage_output$metadata$candidate_diagnostics)) {
      utils::write.csv(
        stage_output$metadata$candidate_diagnostics,
        file = file.path(stage_dir, "candidate_diagnostics.csv"),
        row.names = FALSE
      )
    }

    # feature registry if present
    if (!is.null(stage_output$metadata$registry)) {
      utils::write.csv(
        stage_output$metadata$registry,
        file = file.path(stage_dir, "feature_registry.csv"),
        row.names = FALSE
      )
    }

    # optional generic metadata snapshot
    saveRDS(
      stage_output$metadata,
      file = file.path(stage_dir, "metadata.rds")
    )
  }

  saveRDS(stage_output$result, file = file.path(stage_dir, "result.rds"))
}

write_all_warnings <- function(warnings_df, warnings_file, warnings_csv = NULL) {
  if (is.null(warnings_df) || nrow(warnings_df) == 0) return(invisible(NULL))

  lines <- apply(warnings_df, 1, function(x) {
    paste0("[", x[["stage"]], "][", x[["feature"]], "][", x[["code"]], "] ", x[["message"]])
  })

  writeLines(lines, con = warnings_file, sep = "\n", useBytes = TRUE)

  if (!is.null(warnings_csv)) {
    utils::write.csv(warnings_df, warnings_csv, row.names = FALSE)
  }
}

write_all_dropped <- function(dropped_df, dropped_file) {
  if (is.null(dropped_df) || nrow(dropped_df) == 0) return(invisible(NULL))
  utils::write.csv(dropped_df, dropped_file, row.names = FALSE)
}

write_all_diagnostics <- function(diag_df, diagnostics_file) {
  if (is.null(diag_df) || nrow(diag_df) == 0) return(invisible(NULL))
  utils::write.csv(diag_df, diagnostics_file, row.names = FALSE)
}

write_config_snapshot <- function(cfg, path) {
  if (!requireNamespace("yaml", quietly = TRUE)) return(invisible(NULL))
  yaml::write_yaml(cfg, path)
}