make_selection_split <- function(train_df, fit_cfg) {
  y <- train_df[[fit_cfg$y_col]]
  if (!is.factor(y)) y <- as.factor(y)
  
  pos_label <- fit_cfg$positive
  if (!(pos_label %in% levels(y))) {
    stop("Positive label not found in y.")
  }
  
  neg_levels <- setdiff(levels(y), pos_label)
  if (length(neg_levels) != 1) {
    stop("Binary y required. Found levels: ", paste(levels(y), collapse = ", "))
  }
  neg_label <- neg_levels[[1]]
  
  set.seed(fit_cfg$selection_seed)
  
  idx_pos <- which(y == pos_label)
  idx_neg <- which(y == neg_label)
  
  # If either class is too small, fall back to full-data reuse
  if (length(idx_pos) < 4 || length(idx_neg) < 4) {
    return(list(
      train_idx = seq_len(nrow(train_df)),
      val_idx = seq_len(nrow(train_df)),
      pos_label = pos_label,
      neg_label = neg_label,
      split_mode = "full_reuse_small_class"
    ))
  }
  
  n_val_pos <- floor(length(idx_pos) * fit_cfg$selection_val_frac)
  n_val_neg <- floor(length(idx_neg) * fit_cfg$selection_val_frac)
  
  n_val_pos <- max(1, min(n_val_pos, length(idx_pos) - 2))
  n_val_neg <- max(1, min(n_val_neg, length(idx_neg) - 2))
  
  val_pos <- sample(idx_pos, size = n_val_pos)
  val_neg <- sample(idx_neg, size = n_val_neg)
  
  val_idx <- sort(c(val_pos, val_neg))
  train_idx <- setdiff(seq_len(nrow(train_df)), val_idx)
  
  list(
    train_idx = train_idx,
    val_idx = val_idx,
    pos_label = pos_label,
    neg_label = neg_label,
    split_mode = "stratified_holdout"
  )
}