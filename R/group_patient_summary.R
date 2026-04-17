build_group_patient_summary <- function(stay_id,
                                        sum_llr,
                                        d_dist,
                                        risk_bin_primary,
                                        risk_bin_ddist = NULL,
                                        group_percentiles,
                                        weighted_intra_entropy,
                                        inter_group_entropy) {
  out <- data.frame(
    stay_id = stay_id,
    sum_llr = sum_llr,
    d_dist = d_dist,
    risk_bin_primary = risk_bin_primary,
    weighted_intra_entropy = weighted_intra_entropy,
    inter_group_entropy = inter_group_entropy,
    stringsAsFactors = FALSE
  )

  if (!is.null(risk_bin_ddist)) {
    out$risk_bin_ddist <- risk_bin_ddist
  }

  dplyr::bind_cols(out, group_percentiles)
}