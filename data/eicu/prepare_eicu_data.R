prepare_eicu_data <- function(
  features_path,
  output_dir = "data/eicu/eicu_external_processed.rds"
) {

  output_dir_full <- here::here(output_dir)

  ext <- read.csv(here::here(features_path), stringsAsFactors = FALSE)

  # Normalize labels to match training contract
  ext$gender <- dplyr::case_when(
    ext$gender %in% c("Male", "M", "male", "m") ~ "M",
    ext$gender %in% c("Female", "F", "female", "f") ~ "F",
    TRUE ~ as.character(ext$gender)
  )

  ext$mortality <- as.factor(ext$mortality)
  ext$gender <- as.factor(ext$gender)

  ext <- add_entropy_missing_frac(ext)

  drop_bin_cols <- names(ext)[stringr::str_detect(names(ext), "_bin_(low|normal|high)$")]
  ext_clean <- ext %>% dplyr::select(-dplyr::all_of(drop_bin_cols))

  saveRDS(ext_clean, file.path(output_dir_full, "eicu_data.rds"))

  invisible(ext_clean)
}