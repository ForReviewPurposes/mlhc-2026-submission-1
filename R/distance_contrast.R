distance_contrast_fit <- function(train_L, train_y, positive_label, ridge) {
  y <- as.factor(train_y)
  pos <- positive_label
  neg_levels <- setdiff(levels(y), pos)
  
  if (length(neg_levels) != 1) {
    stop("Binary y required for distance_contrast_fit")
  }
  neg <- neg_levels[[1]]
  
  L_pos <- train_L[y == pos, , drop = FALSE]
  L_neg <- train_L[y == neg, , drop = FALSE]
  
  mu_pos <- colMeans(L_pos)
  mu_neg <- colMeans(L_neg)
  
  Sigma_pos <- stats::cov(L_pos) + diag(ridge, ncol(train_L))
  Sigma_neg <- stats::cov(L_neg) + diag(ridge, ncol(train_L))
  
  R_pos <- chol(Sigma_pos)
  R_neg <- chol(Sigma_neg)
  
  list(
    features = colnames(train_L),
    positive_label = pos,
    negative_label = neg,
    mu_pos = mu_pos,
    mu_neg = mu_neg,
    Sigma_pos = Sigma_pos,
    Sigma_neg = Sigma_neg,
    R_pos = R_pos,
    R_neg = R_neg,
    ridge = ridge
  )
}

distance_contrast_score <- function(L, geom_fit) {
  missing_features <- setdiff(geom_fit$features, colnames(L))
  if (length(missing_features) > 0) {
    stop("L matrix is missing geometry features: ", paste(missing_features, collapse = ", "))
  }
  
  L_use <- L[, geom_fit$features, drop = FALSE]
  
  centered_neg <- sweep(L_use, 2, geom_fit$mu_neg, "-")
  centered_pos <- sweep(L_use, 2, geom_fit$mu_pos, "-")
  
  Z_neg <- backsolve(geom_fit$R_neg, t(centered_neg), transpose = TRUE)
  Z_pos <- backsolve(geom_fit$R_pos, t(centered_pos), transpose = TRUE)
  
  D_neg <- colSums(Z_neg^2)
  D_pos <- colSums(Z_pos^2)
  
  d_dist <- sqrt(D_neg) - sqrt(D_pos)
  
  diagnostics <- data.frame(
    metric = c("mean_d_dist", "sd_d_dist", "min_d_dist", "max_d_dist"),
    value = c(
      mean(d_dist, na.rm = TRUE),
      stats::sd(d_dist, na.rm = TRUE),
      min(d_dist, na.rm = TRUE),
      max(d_dist, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  
  new_stage_output(
    result = data.frame(d_dist = d_dist),
    diagnostics = diagnostics,
    warnings = NULL,
    dropped = NULL,
    metadata = list(n = length(d_dist))
  )
}