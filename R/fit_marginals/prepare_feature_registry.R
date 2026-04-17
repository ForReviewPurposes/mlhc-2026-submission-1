prepare_feature_registry <- function(train_df, fit_cfg) {
  y_col <- fit_cfg$y_col
  
  if (is.null(fit_cfg$features)) {
    features <- setdiff(names(train_df), y_col)
  } else {
    features <- fit_cfg$features
  }
  
  detect_type <- function(col) {
    if (is.numeric(col) || is.integer(col)) "numeric" else "categorical"
  }
  
  is_count_name <- function(fname) {
    grepl("_count$", fname) ||
      grepl("_bin_", fname) ||
      grepl("_low_count$", fname) ||
      grepl("_high_count$", fname)
  }
  
  is_fraction_name <- function(fname) {
    grepl("_frac$", fname) || grepl("_fraction$", fname)
  }
  
  detect_numeric_subtype <- function(x, fname) {
    x_obs <- x[is.finite(x)]
    
    if (is_count_name(fname)) return("count")
    if (is_fraction_name(fname)) return("fraction")
    
    if (length(x_obs) > 0 && all(x_obs >= 0 & x_obs <= 1)) {
      return("fraction")
    }
    
    "continuous"
  }
  
  rows <- lapply(features, function(fname) {
    x <- train_df[[fname]]
    type <- detect_type(x)
    
    subtype <- if (type == "numeric") {
      detect_numeric_subtype(x, fname)
    } else {
      "categorical"
    }
    
    candidate_families <- switch(
      subtype,
      continuous = paste(fit_cfg$numeric_candidates, collapse = ","),
      count = paste(fit_cfg$count_candidates, collapse = ","),
      fraction = paste(fit_cfg$fraction_candidates, collapse = ","),
      categorical = "categorical",
      paste(fit_cfg$numeric_candidates, collapse = ",")
    )
    
    data.frame(
      feature = fname,
      type = type,
      subtype = subtype,
      candidate_families = candidate_families,
      stringsAsFactors = FALSE
    )
  })
  
  registry <- dplyr::bind_rows(rows)
  registry
}