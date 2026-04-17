entropy_norm <- function(x, eps = 1e-12) {
  x <- as.numeric(x)
  s <- sum(x, na.rm = TRUE)

  if (!is.finite(s) || s <= 0) return(NA_real_)
  if (length(x) <= 1) return(0)

  p <- x / s
  p[!is.finite(p) | p <= 0] <- eps

  H <- -sum(p * log(p))
  H / log(length(p))
}

intra_group_entropy <- function(L, group_df, group_pos_sum, groups_cfg) {
  groups <- unique(group_df$group)

  # detailed per-group entropy matrix
  out_list <- lapply(groups, function(g) {
    feats <- group_df$feature[group_df$group == g]
    subL <- L[, feats, drop = FALSE]
    posL <- softPlus(as.matrix(subL))

    if (ncol(posL) <= 1) {
      ent <- rep(0, nrow(posL))
    } else {
      ent <- apply(posL, 1, entropy_norm)
    }

    if (!is.null(groups_cfg$winsor_p) && groups_cfg$winsor_p > 0) {
      ent <- winsorize_vec(ent, p = groups_cfg$winsor_p)
    }

    ent
  })

  intra_entropy_detail <- as.data.frame(out_list, check.names = FALSE)
  names(intra_entropy_detail) <- groups

  # weighted summary score
  mag_norm <- as.matrix(group_pos_sum) / rowSums(as.matrix(group_pos_sum))
  weighted_intra_entropy <- rowSums(mag_norm * as.matrix(intra_entropy_detail), na.rm = TRUE)

  diagnostics <- data.frame(
    metric = c("mean_weighted_intra_entropy", "sd_weighted_intra_entropy"),
    value = c(
      mean(weighted_intra_entropy, na.rm = TRUE),
      stats::sd(weighted_intra_entropy, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  new_stage_output(
    result = list(
      weighted_intra_entropy = weighted_intra_entropy,
      intra_entropy_detail = intra_entropy_detail
    ),
    diagnostics = diagnostics,
    warnings = NULL,
    dropped = NULL,
    metadata = list(groups = groups)
  )
}