evaluate_candidate_families <- function(x,
                                        y,
                                        candidate_families,
                                        split_obj,
                                        pos_label,
                                        neg_label,
                                        fit_cfg,
                                        feature,
                                        subtype) {
  x_train <- x[split_obj$train_idx]
  y_train <- y[split_obj$train_idx]
  x_valid <- x[split_obj$val_idx]
  y_valid <- y[split_obj$val_idx]
  
  rows <- lapply(candidate_families, function(fam) {
    model_pair <- fit_family_pair(
      x_train = x_train,
      y_train = y_train,
      family = fam,
      pos_label = pos_label,
      neg_label = neg_label,
      fit_cfg = fit_cfg
    )
    
    metrics <- compute_marginal_fit_metrics(
      model_pair = model_pair,
      x_train = x_train,
      y_train = y_train,
      x_valid = x_valid,
      y_valid = y_valid,
      pos_label = pos_label,
      neg_label = neg_label,
      family = fam,
      fit_cfg = fit_cfg
    )
    
    cbind(
      data.frame(
        feature = feature,
        subtype = subtype,
        family = fam,
        stringsAsFactors = FALSE
      ),
      metrics
    )
  })
  
  dplyr::bind_rows(rows)
}