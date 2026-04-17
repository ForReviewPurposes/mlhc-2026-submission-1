build_marginal_fit_report <- function(feature_outputs) {
  diagnostics <- bind_stage_diagnostics(lapply(feature_outputs, `[[`, "diagnostics"))
  dropped <- bind_stage_dropped(lapply(feature_outputs, `[[`, "dropped"))
  warnings <- bind_stage_warnings(lapply(feature_outputs, `[[`, "warnings"))
  
  candidate_diag <- dplyr::bind_rows(lapply(feature_outputs, function(x) {
    res <- x$result
    if (is.null(res)) {
      md <- x$metadata
      if (!is.null(md$candidate_diag)) return(md$candidate_diag)
      return(NULL)
    }
    res$candidate_diagnostics
  }))
  
  list(
    diagnostics = diagnostics,
    dropped = dropped,
    warnings = warnings,
    candidate_diagnostics = candidate_diag
  )
}

format_model_params <- function(model_obj, digits = 4) {
  if (is.null(model_obj)) return(NA_character_)

  fam <- model_obj$family
  p <- model_obj$params

  if (fam == "gaussian") {
    return(sprintf(
      "mean=%.*f; sd=%.*f",
      digits, p$mean,
      digits, p$sd
    ))
  }

  if (fam == "lognormal") {
    return(sprintf(
      "meanlog=%.*f; sdlog=%.*f",
      digits, p$meanlog,
      digits, p$sdlog
    ))
  }

  if (fam == "gamma") {
    return(sprintf(
      "shape=%.*f; rate=%.*f",
      digits, p$shape,
      digits, p$rate
    ))
  }

  if (fam == "poisson") {
    return(sprintf(
      "lambda=%.*f",
      digits, p$lambda
    ))
  }

  if (fam == "negbinom") {
    return(sprintf(
      "mu=%.*f; size=%.*f",
      digits, p$mu,
      digits, p$size
    ))
  }

  if (fam == "logit_gaussian") {
    return(sprintf(
      "mean=%.*f; sd=%.*f; eps=%.*g",
      digits, p$mean,
      digits, p$sd,
      digits, p$eps
    ))
  }

  paste(capture.output(str(p)), collapse = "; ")
}