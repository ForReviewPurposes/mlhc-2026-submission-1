get_fallback_families <- function(subtype, fit_cfg) {
  switch(
    subtype,
    continuous = c("gaussian"),
    count = c("poisson"),
    fraction = c("logit_gaussian"),
    character(0)
  )
}