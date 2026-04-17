source(here::here("data/eicu/prepare_eicu_data.R"))
source(here::here("data/mimiciv/prepare_mimic_data.R"))

prepare_eicu_data(
  features_path = "data/eicu/eicu_signal_24h_agg_wide.csv",
  output_dir = "data/eicu/processed"
)