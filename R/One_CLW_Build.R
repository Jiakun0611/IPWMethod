clw_build <- function(vars, sc, sp, sp_des, wts.col,
                      maxit = 20, tol = 1e-4,
                      verbose = FALSE, log_messages = NULL) {

  if (is.null(log_messages)) log_messages <- character(0)
  add_log <- function(msg) {
    log_messages <<- c(log_messages, msg)
    if (verbose) message(msg)
  }

  Xc <- design_matrix(vars = vars, data = sc)
  Xp <- design_matrix(vars = vars, data = sp)

  check_design_identifiability(Xc, method = "One_Ref_CLW")
  check_design_identifiability(Xp, method = "One_Ref_CLW")

  p_dim <- ncol(Xc)
  beta  <- rep(0, p_dim)

  pi.sp0 <- as.vector(expit(Xp %*% beta))
  f.p <- colSums(sp[[wts.col]] * pi.sp0 * Xp)
  f.c <- colSums(Xc)

  beta[1] <- -log(f.p[1] / f.c[1])

  iter <- 0
  crit <- TRUE

  while (crit) {
    eta.sp <- as.vector(Xp %*% beta)
    eta.sc <- as.vector(Xc %*% beta)

    if (any(!is.finite(eta.sp)) || any(!is.finite(eta.sc))) {
      stop(
        "CLW failed: linear predictor became non-finite during iteration.",
        call. = FALSE
      )
    }

    pi.sp <- as.vector(expit(eta.sp))

    if (any(!is.finite(pi.sp))) {
      stop(
        paste(
          "CLW failed: fitted probabilities became non-finite during iteration.",
          "This usually indicates numerical overflow or an unstable Newton step."
        ),
        call. = FALSE
      )
    }

    f.p <- colSums(sp[[wts.col]] * pi.sp * Xp)
    f.c <- colSums(Xc)
    f   <- f.c - f.p

    fprime <- - t(sp[[wts.col]] * pi.sp * (1 - pi.sp) * Xp) %*% Xp

    delta <- solve_newton_step(
      J = fprime,
      g = f,
      method = "One_Ref_CLW"
    )

    beta_new <- beta - delta

    if (any(!is.finite(beta_new))) {
      stop(
        "CLW failed: coefficient update became non-finite.",
        call. = FALSE
      )
    }

    beta <- beta_new
    iter <- iter + 1
    crit <- (iter < maxit && sum(abs(delta)) > tol)
  }

  if (iter >= maxit && sum(abs(delta)) > tol) {
    add_log("CLW NR reached maxit without meeting tol.")
  }

  eta.sc <- as.vector(Xc %*% beta)
  eta.sp <- as.vector(Xp %*% beta)

  if (any(!is.finite(eta.sc)) || any(!is.finite(eta.sp))) {
    stop("Non-finite linear predictor in clw_build.", call. = FALSE)
  }

  pi.sc <- as.vector(expit(eta.sc))
  pi.sp <- as.vector(expit(eta.sp))

  if (any(!is.finite(pi.sc)) || any(!is.finite(pi.sp))) {
    stop("Non-finite fitted probabilities in clw_build.", call. = FALSE)
  }

  wts.sc <- 1 / pi.sc

  if (any(!is.finite(wts.sc))) {
    stop("Non-finite pseudo-weights in clw_build.", call. = FALSE)
  }

  S_beta_full <- t(sp[[wts.col]] * pi.sp * (1 - pi.sp) * Xp) %*% Xp
  D <- compute_D_CLW(sp_des, pi.sp, Xp)

  out <- list(
    weights      = wts.sc,
    coefficients = as.numeric(beta),
    variables    = colnames(Xc),
    iterations   = iter,
    log_messages = log_messages,
    internal = list(
      Xc     = Xc,
      Xp     = Xp,
      pi_sc  = pi.sc,
      pi_sp  = pi.sp,
      S_beta = S_beta_full,
      D      = D
    )
  )

  return(out)
}
