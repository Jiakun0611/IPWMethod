raking_build <- function(vars, sc, sp, sp_des, wts.col,
                         maxit = 20, tol = 1e-4,
                         verbose = FALSE, log_messages = NULL) {

  if (is.null(log_messages)) log_messages <- character(0)
  add_log <- function(msg) {
    log_messages <<- c(log_messages, msg)
    if (verbose) message(msg)
  }

  Xc <- design_matrix(vars = vars, data = sc)
  Xp <- design_matrix(vars = vars, data = sp)

  check_design_identifiability(Xc, method = "One_Ref_Calibration")
  check_design_identifiability(Xp, method = "One_Ref_Calibration")

  p_dim <- ncol(Xc)
  beta  <- rep(0, p_dim)

  f.p <- colSums(sp[[wts.col]] * Xp)
  wts0 <- as.vector(exp(-Xc %*% beta))
  f.c0 <- colSums(wts0 * Xc)
  beta[1] <- -log(f.p[1] / f.c0[1])

  iter <- 0
  crit <- TRUE

  while (crit) {
    eta <- as.vector(Xc %*% beta)

    if (any(!is.finite(eta))) {
      stop(
        "Calibration failed: linear predictor became non-finite during iteration.",
        call. = FALSE
      )
    }

    wts.sc <- as.vector(exp(-eta))
    if (any(!is.finite(wts.sc))) {
      stop(
        paste(
          "Calibration failed: pseudo-weights became non-finite during iteration.",
          "This usually indicates numerical overflow or an unstable Newton step."
        ),
        call. = FALSE
      )
    }

    f.c <- colSums(wts.sc * Xc)
    f   <- f.c - f.p
    fprime <- - t(wts.sc * Xc) %*% Xc

    delta <- solve_newton_step(
      J = fprime,
      g = f,
      method = "One_Ref_Calibration"
    )

    beta_new <- beta - delta

    if (any(!is.finite(beta_new))) {
      stop(
        "Calibration failed: coefficient update became non-finite.",
        call. = FALSE
      )
    }

    beta <- beta_new
    iter <- iter + 1
    crit <- (iter < maxit && sum(abs(delta)) > tol)
  }

  if (iter >= maxit && sum(abs(delta)) > tol) {
    add_log("Calibration NR reached maxit without meeting tol.")
  }

  wts.sc <- as.vector(exp(-Xc %*% beta))
  if (any(!is.finite(wts.sc))) {
    stop("Non-finite pseudo-weights in Calibration_build.", call. = FALSE)
  }

  S_beta_full <- t(wts.sc * Xc) %*% Xc
  D <- compute_D_raking(sp_des, Xp)

  out <- list(
    weights      = wts.sc,
    coefficients = as.numeric(beta),
    variables    = colnames(Xc),
    iterations   = iter,
    log_messages = log_messages,
    internal = list(
      Xc = Xc,
      Xp = Xp,
      S_beta = S_beta_full,
      D  = D
    )
  )
  return(out)
}
