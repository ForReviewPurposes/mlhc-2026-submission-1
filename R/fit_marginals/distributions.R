.safe_logval <- function(x, floor_log = log(1e-12)) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- floor_log
  pmax(x, floor_log)
}

.fit_gaussian_model <- function(x, min_sd = 1e-6) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NULL)
  
  mu <- mean(x)
  sdv <- max(stats::sd(x), min_sd)
  
  list(
    family = "gaussian",
    params = list(mean = mu, sd = sdv),
    logpdf = function(z) {
      .safe_logval(stats::dnorm(z, mean = mu, sd = sdv, log = TRUE))
    }
  )
}

.fit_gmm2_model <- function(x,
                            min_sd = 1e-6,
                            max_iter = 200,
                            tol = 1e-6,
                            min_n = 20) {
  x <- x[is.finite(x)]
  n <- length(x)

  if (n < min_n) return(NULL)
  if (stats::sd(x) < min_sd) return(NULL)

  # Initialize from quartiles
  mu1 <- as.numeric(stats::quantile(x, 0.25, na.rm = TRUE, type = 8))
  mu2 <- as.numeric(stats::quantile(x, 0.75, na.rm = TRUE, type = 8))
  sd1 <- max(stats::sd(x), min_sd)
  sd2 <- max(stats::sd(x), min_sd)
  pi1 <- 0.5
  pi2 <- 0.5

  prev_ll <- -Inf

  for (iter in seq_len(max_iter)) {
    # E-step
    log_r1 <- log(pi1) + stats::dnorm(x, mean = mu1, sd = sd1, log = TRUE)
    log_r2 <- log(pi2) + stats::dnorm(x, mean = mu2, sd = sd2, log = TRUE)

    m <- pmax(log_r1, log_r2)
    denom <- m + log(exp(log_r1 - m) + exp(log_r2 - m))

    r1 <- exp(log_r1 - denom)
    r2 <- exp(log_r2 - denom)

    # M-step
    n1 <- sum(r1)
    n2 <- sum(r2)

    if (n1 < 1e-8 || n2 < 1e-8) return(NULL)

    pi1 <- n1 / n
    pi2 <- n2 / n

    mu1 <- sum(r1 * x) / n1
    mu2 <- sum(r2 * x) / n2

    sd1 <- sqrt(sum(r1 * (x - mu1)^2) / n1)
    sd2 <- sqrt(sum(r2 * (x - mu2)^2) / n2)

    sd1 <- max(sd1, min_sd)
    sd2 <- max(sd2, min_sd)

    ll <- sum(denom)

    if (is.finite(prev_ll) && abs(ll - prev_ll) < tol) break
    prev_ll <- ll
  }

  # Stable ordering by mean
  if (mu1 > mu2) {
    tmp <- mu1; mu1 <- mu2; mu2 <- tmp
    tmp <- sd1; sd1 <- sd2; sd2 <- tmp
    tmp <- pi1; pi1 <- pi2; pi2 <- tmp
  }

  list(
    family = "gmm2",
    params = list(
      pi1 = pi1, pi2 = pi2,
      mean1 = mu1, sd1 = sd1,
      mean2 = mu2, sd2 = sd2
    ),
    logpdf = function(z) {
      out <- rep(log(1e-12), length(z))
      ok <- is.finite(z)
      if (!any(ok)) return(out)

      l1 <- log(pi1) + stats::dnorm(z[ok], mean = mu1, sd = sd1, log = TRUE)
      l2 <- log(pi2) + stats::dnorm(z[ok], mean = mu2, sd = sd2, log = TRUE)

      m <- pmax(l1, l2)
      vals <- m + log(exp(l1 - m) + exp(l2 - m))

      out[ok] <- .safe_logval(vals)
      out
    }
  )
}

.fit_lognormal_model <- function(x, min_sd = 1e-6) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || any(x <= 0)) return(NULL)
  
  lx <- log(x)
  mu <- mean(lx)
  sdv <- max(stats::sd(lx), min_sd)
  
  list(
    family = "lognormal",
    params = list(meanlog = mu, sdlog = sdv),
    logpdf = function(z) {
      out <- rep(log(1e-12), length(z))
      ok <- is.finite(z) & (z > 0)
      vals <- stats::dlnorm(z[ok], meanlog = mu, sdlog = sdv, log = TRUE)
      out[ok] <- .safe_logval(vals)
      out
    }
  )
}

.fit_gamma_model <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 5 || any(x <= 0) || stats::sd(x) < 1e-8) return(NULL)
  
  fit <- tryCatch(
    suppressWarnings(MASS::fitdistr(x, densfun = "gamma")),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  shape <- unname(fit$estimate["shape"])
  rate  <- unname(fit$estimate["rate"])
  
  if (!is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0) return(NULL)
  
  list(
    family = "gamma",
    params = list(shape = shape, rate = rate),
    logpdf = function(z) {
      out <- rep(log(1e-12), length(z))
      ok <- is.finite(z) & (z > 0)
      vals <- stats::dgamma(z[ok], shape = shape, rate = rate, log = TRUE)
      out[ok] <- .safe_logval(vals)
      out
    }
  )
}

.fit_poisson_model <- function(x, eps = 1e-8) {
  x <- x[is.finite(x)]
  if (length(x) < 1) return(NULL)
  
  x <- pmax(0, round(x))
  lambda <- max(mean(x), eps)
  
  list(
    family = "poisson",
    params = list(lambda = lambda),
    logpdf = function(z) {
      z <- pmax(0, round(z))
      .safe_logval(stats::dpois(z, lambda = lambda, log = TRUE))
    }
  )
}

.fit_negbinom_model <- function(x, eps = 1e-8, winsor_q = 0.995) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NULL)
  
  x <- pmax(0, round(x))
  
  if (length(unique(x)) > 1) {
    qhi <- as.numeric(stats::quantile(x, probs = winsor_q, na.rm = TRUE, type = 8))
    x <- pmin(x, qhi)
  }
  
  mu <- mean(x)
  v  <- stats::var(x)
  if (!is.finite(v)) v <- mu
  
  if (v <= mu + eps) {
    size <- 1e6
  } else {
    size <- (mu^2) / max(v - mu, eps)
    size <- max(size, eps)
  }
  
  list(
    family = "negbinom",
    params = list(mu = max(mu, eps), size = size),
    logpdf = function(z) {
      z <- pmax(0, round(z))
      .safe_logval(stats::dnbinom(z, mu = mu, size = size, log = TRUE))
    }
  )
}

.fit_logit_gaussian_model <- function(x, min_sd = 1e-6, eps = 1e-6) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NULL)
  if (any(x < 0 | x > 1)) return(NULL)
  
  x_clip <- pmin(pmax(x, eps), 1 - eps)
  z <- stats::qlogis(x_clip)
  
  mu <- mean(z)
  sdv <- max(stats::sd(z), min_sd)
  
  list(
    family = "logit_gaussian",
    params = list(mean = mu, sd = sdv, eps = eps),
    logpdf = function(v) {
      out <- rep(log(1e-12), length(v))
      ok <- is.finite(v) & (v >= 0) & (v <= 1)
      if (!any(ok)) return(out)
      
      v_clip <- pmin(pmax(v[ok], eps), 1 - eps)
      z <- stats::qlogis(v_clip)
      vals <- stats::dnorm(z, mean = mu, sd = sdv, log = TRUE) - log(v_clip * (1 - v_clip))
      out[ok] <- .safe_logval(vals)
      out
    }
  )
}