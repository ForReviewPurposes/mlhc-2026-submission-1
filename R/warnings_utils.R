new_warning_row <- function(stage,
                            feature = NA_character_,
                            code,
                            message,
                            severity = "warning") {
  data.frame(
    stage = stage,
    feature = feature,
    code = code,
    message = message,
    severity = severity,
    stringsAsFactors = FALSE
  )
}

new_dropped_row <- function(stage,
                            feature,
                            reason) {
  data.frame(
    stage = stage,
    feature = feature,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

bind_stage_warnings <- function(...) {
  dfs <- list(...)
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (length(dfs) == 0) return(NULL)
  dplyr::bind_rows(dfs)
}

bind_stage_dropped <- function(...) {
  dfs <- list(...)
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (length(dfs) == 0) return(NULL)
  dplyr::bind_rows(dfs)
}

bind_stage_diagnostics <- function(...) {
  dfs <- list(...)
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (length(dfs) == 0) return(NULL)
  dplyr::bind_rows(dfs)
}