fit_family_by_name <- function(x, family, fit_cfg) {
  switch(
    family,
    gaussian       = .fit_gaussian_model(x, min_sd = fit_cfg$gaussian_min_sd),
    gmm2           = .fit_gmm2_model(
      x,
      min_sd = fit_cfg$gaussian_min_sd,
      max_iter = fit_cfg$gmm2_max_iter,
      tol = fit_cfg$gmm2_tol,
      min_n = fit_cfg$gmm2_min_n
    ),
    lognormal      = .fit_lognormal_model(x, min_sd = fit_cfg$gaussian_min_sd),
    gamma          = .fit_gamma_model(x),
    poisson        = .fit_poisson_model(x),
    negbinom       = .fit_negbinom_model(x),
    logit_gaussian = .fit_logit_gaussian_model(
      x,
      min_sd = fit_cfg$gaussian_min_sd,
      eps = fit_cfg$fraction_eps
    ),
    NULL
  )
}

fit_family_pair <- function(x_train, y_train, family, pos_label, neg_label, fit_cfg) {
  if (identical(family, "gmm2_pos_gaussian_neg")) {
    model_pos <- fit_family_by_name(x_train[y_train == pos_label], "gmm2", fit_cfg)
    model_neg <- fit_family_by_name(x_train[y_train == neg_label], "gaussian", fit_cfg)

    return(list(
      pos = model_pos,
      neg = model_neg
    ))
  }

  model_pos <- fit_family_by_name(x_train[y_train == pos_label], family, fit_cfg)
  model_neg <- fit_family_by_name(x_train[y_train == neg_label], family, fit_cfg)

  list(
    pos = model_pos,
    neg = model_neg
  )
}