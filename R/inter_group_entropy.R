inter_group_entropy <- function(group_pos_sum) {
  ent <- apply(as.matrix(group_pos_sum), 1, entropy_norm)

  diagnostics <- data.frame(
    metric = "inter_group_entropy",
    mean_value = mean(ent, na.rm = TRUE),
    sd_value = stats::sd(ent, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  new_stage_output(
    result = list(inter_entropy = ent),
    diagnostics = diagnostics,
    warnings = NULL,
    dropped = NULL,
    metadata = NULL
  )
}