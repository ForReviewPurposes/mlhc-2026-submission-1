assess_candidate_fit <- function(candidate_row, fit_cfg) {
  if (!isTRUE(candidate_row$fit_success_pos) || !isTRUE(candidate_row$fit_success_neg)) {
    return(list(accept = FALSE, reason = "fit_fail"))
  }
  
  if (is.na(candidate_row$heldout_balanced_trimmed_ll) || !is.finite(candidate_row$heldout_balanced_trimmed_ll)) {
    return(list(accept = FALSE, reason = "nonfinite_balanced_trimmed_ll"))
  }
  
  if (is.na(candidate_row$finite_logpdf_frac_pos) ||
      candidate_row$finite_logpdf_frac_pos < fit_cfg$min_finite_logpdf_frac) {
    return(list(accept = FALSE, reason = "low_finite_logpdf_frac_pos"))
  }
  
  if (is.na(candidate_row$finite_logpdf_frac_neg) ||
      candidate_row$finite_logpdf_frac_neg < fit_cfg$min_finite_logpdf_frac) {
    return(list(accept = FALSE, reason = "low_finite_logpdf_frac_neg"))
  }
  
  if (is.na(candidate_row$finite_llr_frac) ||
      candidate_row$finite_llr_frac < fit_cfg$min_finite_llr_frac) {
    return(list(accept = FALSE, reason = "low_finite_llr_frac"))
  }
  
  if (!is.na(candidate_row$p99_abs_llr) &&
      candidate_row$p99_abs_llr > fit_cfg$max_abs_llr_p99) {
    return(list(accept = FALSE, reason = "llr_tail_instability"))
  }
  
  list(accept = TRUE, reason = NA_character_)
}