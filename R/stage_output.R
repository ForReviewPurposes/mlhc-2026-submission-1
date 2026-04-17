new_stage_output <- function(result = NULL,
                             diagnostics = NULL,
                             warnings = NULL,
                             dropped = NULL,
                             metadata = NULL) {
  list(
    result = result,
    diagnostics = diagnostics,
    warnings = warnings,
    dropped = dropped,
    metadata = metadata
  )
}