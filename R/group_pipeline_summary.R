make_quantile_bins <- function(x, labels = c("Q1","Q2","Q3","Q4")) {
  qs <- unique(stats::quantile(x, probs = seq(0, 1, length.out = length(labels) + 1), na.rm = TRUE))
  if (length(qs) < 3) {
    return(rep(NA_character_, length(x)))
  }
  as.character(cut(x, breaks = qs, include.lowest = TRUE, labels = labels))
}

build_group_patient_summaries <- function(train_L,
                                          eval_L,
                                          train_score_tbl,
                                          eval_score_tbl,
                                          train_ids,
                                          eval_ids,
                                          ddist_bin_cfg,
                                          group_map_path,
                                          groups_cfg) {
  group_df <- read_group_map(here::here(group_map_path))

  # ---- train metrics ----
  s_pos_train <- group_positive_evidence(
    L = train_L,
    group_df = group_df,
    groups_cfg = groups_cfg
  )

  s_inter_train <- inter_group_entropy(
    group_pos_sum = s_pos_train$result$group_pos_sum
  )

  s_intra_train <- intra_group_entropy(
    L = train_L,
    group_df = s_pos_train$result$group_df,
    group_pos_sum = s_pos_train$result$group_pos_sum,
    groups_cfg = groups_cfg
  )

  s_ecdf <- fit_group_ecdfs(
    group_pos_sum = s_pos_train$result$group_pos_sum,
    groups_cfg = groups_cfg
  )

  s_pct_train <- score_group_percentiles(
    group_pos_sum = s_pos_train$result$group_pos_sum,
    ecdf_fit = s_ecdf$result,
    groups_cfg = groups_cfg
  )

  train_risk_bin_primary <- make_quantile_bins(train_score_tbl$sum_llr)

  train_risk_bin_ddist <- as.character(cut(
    train_score_tbl$d_dist,
    breaks = ddist_bin_cfg$breaks,
    labels = ddist_bin_cfg$labels,
    include.lowest = TRUE
  ))

  train_patient_summary <- build_group_patient_summary(
    stay_id = train_ids,
    sum_llr = train_score_tbl$sum_llr,
    d_dist = train_score_tbl$d_dist,
    risk_bin_primary = train_risk_bin_primary,
    risk_bin_ddist = train_risk_bin_ddist,
    group_percentiles = s_pct_train$result$group_percentiles,
    weighted_intra_entropy = s_intra_train$result$weighted_intra_entropy,
    inter_group_entropy = s_inter_train$result$inter_entropy
  )

  train_interp <- build_interpretation_strings(
    patient_summary = train_patient_summary,
    groups_cfg = groups_cfg
  )

  train_interp$stay_id <- train_patient_summary$stay_id
  train_patient_summary <- dplyr::left_join(
    train_patient_summary,
    train_interp,
    by = "stay_id"
  )

  # ---- eval metrics using train ECDF ----
  s_pos_eval <- group_positive_evidence(
    L = eval_L,
    group_df = group_df,
    groups_cfg = groups_cfg
  )

  s_inter_eval <- inter_group_entropy(
    group_pos_sum = s_pos_eval$result$group_pos_sum
  )

  s_intra_eval <- intra_group_entropy(
    L = eval_L,
    group_df = s_pos_eval$result$group_df,
    group_pos_sum = s_pos_eval$result$group_pos_sum,
    groups_cfg = groups_cfg
  )

  s_pct_eval <- score_group_percentiles(
    group_pos_sum = s_pos_eval$result$group_pos_sum,
    ecdf_fit = s_ecdf$result,
    groups_cfg = groups_cfg
  )

  eval_risk_bin_primary  <- make_quantile_bins(eval_score_tbl$sum_llr)

  eval_risk_bin_ddist <- as.character(cut(
    eval_score_tbl$d_dist,
    breaks = ddist_bin_cfg$breaks,
    labels = ddist_bin_cfg$labels,
    include.lowest = TRUE
  ))

  eval_patient_summary <- build_group_patient_summary(
    stay_id = eval_ids,
    sum_llr = eval_score_tbl$sum_llr,
    d_dist = eval_score_tbl$d_dist,
    risk_bin_primary = eval_risk_bin_primary,
    risk_bin_ddist = eval_risk_bin_ddist,
    group_percentiles = s_pct_eval$result$group_percentiles,
    weighted_intra_entropy = s_intra_eval$result$weighted_intra_entropy,
    inter_group_entropy = s_inter_eval$result$inter_entropy
  )

  eval_interp <- build_interpretation_strings(
    patient_summary = eval_patient_summary,
    groups_cfg = groups_cfg
  )

  eval_interp$stay_id <- eval_patient_summary$stay_id
  eval_patient_summary <- dplyr::left_join(
    eval_patient_summary,
    eval_interp,
    by = "stay_id"
  )

  list(
    train_patient_summary = train_patient_summary,
    eval_patient_summary = eval_patient_summary,
    train_objects = list(
      group_pos_sum = s_pos_train$result$group_pos_sum,
      weighted_intra_entropy = s_intra_train$result$weighted_intra_entropy,
      inter_group_entropy = s_inter_train$result$inter_entropy,
      intra_entropy_detail = s_intra_train$result$intra_entropy_detail
    ),
    eval_objects = list(
      group_pos_sum = s_pos_eval$result$group_pos_sum,
      weighted_intra_entropy = s_intra_eval$result$weighted_intra_entropy,
      inter_group_entropy = s_inter_eval$result$inter_entropy,
      intra_entropy_detail = s_intra_eval$result$intra_entropy_detail
    ),
    ecdf_object = s_ecdf$result
  )
}