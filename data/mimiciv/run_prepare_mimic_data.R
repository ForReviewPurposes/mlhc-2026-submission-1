source(here::here("data/mimiciv/prepare_mimic_data.R"))

prepare_mimic_data(
  features_path = "data/mimiciv/mimiciv_ext_24h_features_wide_3.csv",
  gcs_path = "data/mimiciv/mimiciv_ext_24h_gcs_features_wide_2.csv",
  output_dir = "data/mimiciv/processed",
  seed = 42
)