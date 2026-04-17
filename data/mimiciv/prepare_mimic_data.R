suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(caret)
  library(here)
})

safe_entropy3 <- function(low, mid, high, eps = 1e-12) {
  total <- low + mid + high
  
  out <- rep(NA_real_, length(total))
  ok <- is.finite(total) & total > 0
  
  if (any(ok)) {
    p1 <- low[ok]  / total[ok]
    p2 <- mid[ok]  / total[ok]
    p3 <- high[ok] / total[ok]
    
    P <- cbind(p1, p2, p3)
    P[P <= 0 | !is.finite(P)] <- eps
    
    out[ok] <- -rowSums(P * log(P))
  }
  
  out
}

add_entropy_missing_frac <- function(df) {
  nms <- names(df)
  
  # 1) entropy from *_bin_low/_normal/_high
  bin_low_cols <- nms[str_detect(nms, "_bin_low$")]
  
  for (low_col in bin_low_cols) {
    stem <- str_remove(low_col, "_bin_low$")
    mid_col <- paste0(stem, "_bin_normal")
    high_col <- paste0(stem, "_bin_high")
    
    if (mid_col %in% names(df) && high_col %in% names(df)) {
      ent_col <- paste0(stem, "_entropy")
      df[[ent_col]] <- safe_entropy3(
        low = df[[low_col]],
        mid = df[[mid_col]],
        high = df[[high_col]]
      )
    }
  }
  
  # 2) convert *_missing to factor
  missing_cols <- nms[str_detect(nms, "_missing$")]
  
  if (length(missing_cols) > 0) {
    df <- df %>%
      mutate(across(all_of(missing_cols), ~ factor(.x, levels = c(0, 1))))
  }
  
  # 3) low/high count fractions
  low_count_cols <- names(df)[str_detect(names(df), "_low_count$")]
  high_count_cols <- names(df)[str_detect(names(df), "_high_count$")]
  
  make_frac <- function(count_vec, denom_vec) {
    out <- rep(NA_real_, length(count_vec))
    ok <- is.finite(denom_vec) & denom_vec > 0 & is.finite(count_vec)
    out[ok] <- count_vec[ok] / denom_vec[ok]
    out
  }
  
  for (col in low_count_cols) {
    stem <- str_remove(col, "_low_count$")
    denom_col <- paste0(stem, "_count")
    frac_col <- paste0(stem, "_low_frac")
    
    if (denom_col %in% names(df)) {
      df[[frac_col]] <- make_frac(df[[col]], df[[denom_col]])
    }
  }
  
  for (col in high_count_cols) {
    stem <- str_remove(col, "_high_count$")
    denom_col <- paste0(stem, "_count")
    frac_col <- paste0(stem, "_high_frac")
    
    if (denom_col %in% names(df)) {
      df[[frac_col]] <- make_frac(df[[col]], df[[denom_col]])
    }
  }
  
  df
}

prepare_mimic_data <- function(
  features_path,
  gcs_path,
  output_dir = "data/mimiciv/processed",
  seed = 42
) {
  output_dir_full <- here::here(output_dir)
  dir.create(output_dir_full, recursive = TRUE, showWarnings = FALSE)
  
  mimiciv_import <- read.csv(here::here(features_path), stringsAsFactors = FALSE)
  mimiciv_gcs <- read.csv(here::here(gcs_path), stringsAsFactors = FALSE)
  
  mimiciv_import$mortality <- as.factor(mimiciv_import$mortality)
  mimiciv_import$gender <- as.factor(mimiciv_import$gender)
  
  mimiciv_import <- dplyr::inner_join(mimiciv_import, mimiciv_gcs, by = dplyr::join_by(stay_id))
  
  mimiciv_import <- add_entropy_missing_frac(mimiciv_import)
  
  drop_bin_cols <- names(mimiciv_import)[str_detect(names(mimiciv_import), "_bin_(low|normal|high)$")]
  mimic_clean <- mimiciv_import %>% dplyr::select(-all_of(drop_bin_cols))
  
  set.seed(seed)
  
  train_idx <- caret::createDataPartition(mimic_clean$mortality, p = 0.5, list = FALSE)
  mimic_clean_train <- mimic_clean[train_idx, ]
  mimic_clean_eval <- mimic_clean[-train_idx, ]
  
  val_idx <- caret::createDataPartition(mimic_clean_eval$mortality, p = 0.5, list = FALSE)
  mimic_clean_val <- mimic_clean_eval[val_idx, ]
  mimic_clean_test <- mimic_clean_eval[-val_idx, ]
  
  saveRDS(mimic_clean_train, file.path(output_dir_full, "mimic_clean_train.rds"))
  saveRDS(mimic_clean_val, file.path(output_dir_full, "mimic_clean_val.rds"))
  saveRDS(mimic_clean_test, file.path(output_dir_full, "mimic_clean_test.rds"))
  
  split_summary <- data.frame(
    split = c("train", "val", "test"),
    n = c(nrow(mimic_clean_train), nrow(mimic_clean_val), nrow(mimic_clean_test)),
    mortality_rate = c(
      mean(as.numeric(as.character(mimic_clean_train$mortality))),
      mean(as.numeric(as.character(mimic_clean_val$mortality))),
      mean(as.numeric(as.character(mimic_clean_test$mortality)))
    )
  )
  
  utils::write.csv(split_summary, file.path(output_dir_full, "split_summary.csv"), row.names = FALSE)
  
  invisible(list(
    train = mimic_clean_train,
    val = mimic_clean_val,
    test = mimic_clean_test,
    split_summary = split_summary
  ))
}