suppressPackageStartupMessages({
  library(dplyr)
  library(here)
})

source(here::here("R", "config_defaults.R"))
source(here::here("R", "io_read_config.R"))
source(here::here("R", "io_write_outputs.R"))
source(here::here("R", "stage_output.R"))
source(here::here("R", "warnings_utils.R"))

source(here::here("R", "fit_marginals", "preprocess_feature.R"))
source(here::here("R", "fit_marginals", "distributions.R"))
source(here::here("R", "evidence_transform.R"))
source(here::here("R", "distance_contrast.R"))
source(here::here("R", "plot_distance_contrast.R"))

source(here::here("analysis", "internal_validation", "survival_config.R"))
source(here::here("analysis", "internal_validation", "survival_analysis.R"))

write_yaml_snapshot <- function(obj, path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    return(invisible(NULL))
  }
  yaml::write_yaml(obj, path)
}

args <- commandArgs(trailingOnly = TRUE)

config_path <- if (length(args) == 0) {
  here::here("configs", "survival_analysis.yaml")
} else {
  here::here(args[[1]])
}

surv_cfg <- read_survival_config(config_path)
run_dir <- resolve_run_dir(surv_cfg)

run_obj <- readRDS(file.path(run_dir, "reports", "run_output.rds"))
base_cfg <- run_obj$config

ddist_bin_cfg <- resolve_ddist_bins(surv_cfg, base_cfg)

marginal_fit <- readRDS(file.path(run_dir, "model", "marginal_fit.rds"))
geometry_fit <- readRDS(file.path(run_dir, "model", "geometry_fit.rds"))

data_paths <- resolve_data_paths(surv_cfg, base_cfg)

train_df <- read_input_data(data_paths$train_data_path)
eval_df  <- read_input_data(data_paths$eval_data_path)

# rebuild train/eval evidence scores from saved model
s_train_tr <- evidence_transform(
  new_df = train_df,
  marginal_fit = marginal_fit,
  transform_cfg = base_cfg$transform
)

s_eval_tr <- evidence_transform(
  new_df = eval_df,
  marginal_fit = marginal_fit,
  transform_cfg = base_cfg$transform
)

common_features <- intersect(
  colnames(s_train_tr$result$L),
  colnames(s_eval_tr$result$L)
)

if (length(common_features) == 0) {
  stop("No common transformed features remain between train and evaluation datasets.")
}

train_L <- s_train_tr$result$L[, common_features, drop = FALSE]
eval_L  <- s_eval_tr$result$L[, common_features, drop = FALSE]

train_score_tbl <- data.frame(
  stay_id = train_df[[base_cfg$data$id_col]],
  sum_llr = rowSums(train_L),
  d_dist = distance_contrast_score(train_L, geometry_fit)$result$d_dist,
  stringsAsFactors = FALSE
)

eval_score_tbl <- data.frame(
  stay_id = eval_df[[base_cfg$data$id_col]],
  sum_llr = rowSums(eval_L),
  d_dist = distance_contrast_score(eval_L, geometry_fit)$result$d_dist,
  stringsAsFactors = FALSE
)

# attach train/eval baseline probabilities if saved
eval_train_path <- file.path(run_dir, "scores", "evaluation_table_train.csv")
eval_eval_path  <- file.path(
  run_dir,
  "scores",
  paste0("evaluation_table_", surv_cfg$data_sources$eval_label, ".csv")
)

if (!file.exists(eval_eval_path) && identical(surv_cfg$data_sources$eval_label, "val")) {
  eval_eval_path <- file.path(run_dir, "scores", "evaluation_table_val.csv")
}

if (file.exists(eval_train_path)) {
  eval_train <- read.csv(eval_train_path, stringsAsFactors = FALSE)
  train_score_tbl <- train_score_tbl %>%
    dplyr::left_join(eval_train %>% dplyr::select(dplyr::any_of(c("stay_id", "p_nb", "p_rf", "p_rf_llr"))), by = "stay_id")
}

if (file.exists(eval_eval_path)) {
  eval_eval <- read.csv(eval_eval_path, stringsAsFactors = FALSE)
  eval_score_tbl <- eval_score_tbl %>%
    dplyr::left_join(eval_eval %>% dplyr::select(dplyr::any_of(c("stay_id", "p_nb", "p_rf", "p_rf_llr"))), by = "stay_id")
}

# attach train/eval group summaries if saved
train_patient_summary_path <- file.path(run_dir, "scores", "group_patient_summary_train.csv")
eval_patient_summary_path  <- file.path(
  run_dir,
  "scores",
  paste0("group_patient_summary_", surv_cfg$data_sources$eval_label, ".csv")
)

if (!file.exists(eval_patient_summary_path) && identical(surv_cfg$data_sources$eval_label, "val")) {
  eval_patient_summary_path <- file.path(run_dir, "scores", "group_patient_summary_val.csv")
}

train_patient_summary <- if (file.exists(train_patient_summary_path)) {
  read.csv(train_patient_summary_path, stringsAsFactors = FALSE)
} else {
  NULL
}

eval_patient_summary <- if (file.exists(eval_patient_summary_path)) {
  read.csv(eval_patient_summary_path, stringsAsFactors = FALSE)
} else {
  NULL
}

survival_tbl <- read_survival_table(
  path = here::here(surv_cfg$survival_data$time_to_event_path),
  age_col = surv_cfg$survival_data$age_col,
  duration_col = surv_cfg$survival_data$duration_col,
  event_col = surv_cfg$survival_data$event_col
)

severity_tbl <- read_severity_scores(
  path = here::here(surv_cfg$survival_data$severity_scores_path),
  score_cols = surv_cfg$benchmarks$severity_score_cols
)

train_surv_df <- build_survival_dataset(
  processed_df = train_df,
  score_tbl = train_score_tbl,
  patient_summary = train_patient_summary,
  survival_tbl = survival_tbl,
  severity_tbl = severity_tbl,
  id_col = base_cfg$data$id_col
)

eval_surv_df <- build_survival_dataset(
  processed_df = eval_df,
  score_tbl = eval_score_tbl,
  patient_summary = eval_patient_summary,
  survival_tbl = survival_tbl,
  severity_tbl = severity_tbl,
  id_col = base_cfg$data$id_col
)

surv_out_dir <- file.path(run_dir, surv_cfg$output$subfolder_name)
dir.create(surv_out_dir, recursive = TRUE, showWarnings = FALSE)

# write snapshots and manifest
write_yaml_snapshot(
  surv_cfg,
  file.path(surv_out_dir, "survival_run_config_snapshot.yaml")
)

write_yaml_snapshot(
  base_cfg,
  file.path(surv_out_dir, "upstream_run_config_snapshot.yaml")
)

manifest_df <- data.frame(
  field = c(
    "run_dir",
    "eval_label",
    "train_data_path",
    "eval_data_path",
    "survival_time_to_event_path",
    "severity_scores_path",
    "ddist_bin_breaks",
    "ddist_bin_labels",
    "created_at"
  ),
  value = c(
    run_dir,
    surv_cfg$data_sources$eval_label,
    data_paths$train_data_path,
    data_paths$eval_data_path,
    here::here(surv_cfg$survival_data$time_to_event_path),
    here::here(surv_cfg$survival_data$severity_scores_path),
    paste(ddist_bin_cfg$breaks, collapse = ", "),
    paste(ddist_bin_cfg$labels, collapse = ", "),
    as.character(Sys.time())
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  manifest_df,
  file.path(surv_out_dir, "survival_run_manifest.csv"),
  row.names = FALSE
)

if (!is.null(train_patient_summary)) {
  utils::write.csv(
    train_patient_summary,
    file.path(surv_out_dir, "group_patient_summary_train.csv"),
    row.names = FALSE
  )
}

if (!is.null(eval_patient_summary)) {
  utils::write.csv(
    eval_patient_summary,
    file.path(surv_out_dir, paste0("group_patient_summary_", surv_cfg$data_sources$eval_label, ".csv")),
    row.names = FALSE
  )
}

plot_km_by_ddist_bin(
  train_surv_df,
  ddist_bin_cfg = ddist_bin_cfg,
  out_path = file.path(surv_out_dir, "km_ddist_bins_train.png")
)

plot_km_by_sum_llr_bin(
  train_surv_df,
  out_path = file.path(surv_out_dir, "km_sum_llr_bins_train.png")
)

run_stage9_survival(
  train_surv_df = train_surv_df,
  eval_surv_df = eval_surv_df,
  out_dir = surv_out_dir,
  surv_cfg = surv_cfg
)