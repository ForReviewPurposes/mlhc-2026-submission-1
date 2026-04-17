softPlus <- function(x) {
  log1p(exp(-abs(x))) + pmax(x, 0)
}

read_group_map <- function(path) {
  gm <- read_input_data(path)

  required_cols <- c("feature", "group")
  missing_cols <- setdiff(required_cols, names(gm))
  if (length(missing_cols) > 0) {
    stop("Group map missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  gm <- gm[, required_cols, drop = FALSE]
  gm$feature <- as.character(gm$feature)
  gm$group <- as.character(gm$group)

  gm
}

group_positive_evidence <- function(L, group_df, groups_cfg) {
  missing_features <- setdiff(group_df$feature, colnames(L))
  if (length(missing_features) > 0) {
    warning_df <- new_warning_row(
      stage = "group_positive_evidence",
      feature = NA_character_,
      code = "group_map_features_missing_from_L",
      message = paste("Some group-map features are missing from L and will be ignored:", paste(missing_features, collapse = ", "))
    )
    group_df <- group_df[group_df$feature %in% colnames(L), , drop = FALSE]
  } else {
    warning_df <- NULL
  }

  groups <- unique(group_df$group)

  out_list <- lapply(groups, function(g) {
    feats <- group_df$feature[group_df$group == g]
    subL <- L[, feats, drop = FALSE]
    posL <- softPlus(as.matrix(subL))

    group_sum <- rowSums(posL)

    if (!is.null(groups_cfg$winsor_p) && groups_cfg$winsor_p > 0) {
      group_sum <- winsorize_vec(group_sum, p = groups_cfg$winsor_p)
    }

    group_sum
  })

  group_pos_sum <- as.data.frame(out_list, check.names = FALSE)
  names(group_pos_sum) <- groups

  diagnostics <- data.frame(
    group = groups,
    n_features = vapply(groups, function(g) sum(group_df$group == g), integer(1)),
    stringsAsFactors = FALSE
  )

  new_stage_output(
    result = list(
      group_pos_sum = group_pos_sum,
      group_df = group_df
    ),
    diagnostics = diagnostics,
    warnings = warning_df,
    dropped = NULL,
    metadata = list(groups = groups)
  )
}