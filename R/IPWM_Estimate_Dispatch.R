dispatch_estimator <- function(build, yz_data) {

  method  <- build$method
  D       <- build$internal$D
  S_beta  <- build$internal$S_beta

  Y <- yz_data$Y
  Z <- yz_data$Z
  w <- yz_data$w
  X <- yz_data$X

  method_key <- tolower(method)

  out <- switch(
    method_key,
    "calibration" = raking_estimate(Y, Z, w, X, D, S_beta),
    "cali"        = raking_estimate(Y, Z, w, X, D, S_beta),

    "alp"         = alp_estimate(Y, Z, w, X, D, S_beta),
    "clw"         = clw_estimate(Y, Z, w, X, D, S_beta),

    "multi"       = multi_estimate(Y, Z, w, X, D, S_beta),

    stop(sprintf(
      "Unknown method '%s'. Must be one of: 'alp', 'clw', 'calibration'/'cali'/'raking', 'multi'.",
      method
    ), call. = FALSE)
  )
}
