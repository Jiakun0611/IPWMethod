alp_build <- function(vars, sc, sp, sp_des, wts.col,
                      maxit = 20, tol = 1e-4,
                      verbose = FALSE, log_messages = NULL) {

  if (is.null(log_messages)) log_messages <- character(0)
  add_log <- function(msg) {
    log_messages <<- c(log_messages, msg)
    if (verbose) message(msg)
  }

  Xc <- design_matrix(vars = vars, data = sc)
  Xp <- design_matrix(vars = vars, data = sp)

  check_design_identifiability(Xc, method = "One_Ref_ALP")
  check_design_identifiability(Xp, method = "One_Ref_ALP")

  p_dim <- ncol(Xc)
  beta  <- rep(0, p_dim)

  p.sp0 <- as.vector(expit(Xp %*% beta))
  p.sc0 <- as.vector(expit(Xc %*% beta))

  f.p <- colSums(sp[[wts.col]] * p.sp0 * Xp)
  f.c <- colSums((1 - p.sc0) * Xc)

  beta[1] <- -log(f.p[1] / f.c[1])

  iter <- 0
  crit <- TRUE

  while (crit) {
    eta.sp <- as.vector(Xp %*% beta)
    eta.sc <- as.vector(Xc %*% beta)

    if (any(!is.finite(eta.sp)) || any(!is.finite(eta.sc))) {
      stop(
        "ALP failed: linear predictor became non-finite during iteration.",
        call. = FALSE
      )
    }

    p.sp <- as.vector(expit(eta.sp))
    p.sc <- as.vector(expit(eta.sc))

    if (any(!is.finite(p.sp)) || any(!is.finite(p.sc))) {
      stop(
        paste(
          "ALP failed: fitted probabilities became non-finite during iteration.",
          "This usually indicates numerical overflow or an unstable Newton step."
        ),
        call. = FALSE
      )
    }

    score <- colSums((1 - p.sc) * Xc) -
      colSums(sp[[wts.col]] * p.sp * Xp)

    H <- - t(p.sc * (1 - p.sc) * Xc) %*% Xc -
      t(sp[[wts.col]] * p.sp * (1 - p.sp) * Xp) %*% Xp

    delta <- solve_newton_step(
      J = H,
      g = score,
      method = "One_Ref_ALP"
    )

    beta_new <- beta - delta

    if (any(!is.finite(beta_new))) {
      stop(
        "ALP failed: coefficient update became non-finite.",
        call. = FALSE
      )
    }

    beta <- beta_new
    iter <- iter + 1
    crit <- (iter < maxit && sum(abs(delta)) > tol)
  }

  if (iter >= maxit && sum(abs(delta)) > tol) {
    add_log("ALP NR reached maxit without meeting tol.")
  }

  eta.sc <- as.vector(Xc %*% beta)
  eta.sp <- as.vector(Xp %*% beta)

  if (any(!is.finite(eta.sc)) || any(!is.finite(eta.sp))) {
    stop("Non-finite linear predictor in alp_build.", call. = FALSE)
  }

  p.sc <- as.vector(expit(eta.sc))
  p.sp <- as.vector(expit(eta.sp))

  if (any(!is.finite(p.sc)) || any(!is.finite(p.sp))) {
    stop("Non-finite fitted probabilities in alp_build.", call. = FALSE)
  }

  pi.sc <- p.sc / (1 - p.sc)
  wts.sc <- 1 / pi.sc

  if (any(!is.finite(wts.sc))) {
    stop("Non-finite pseudo-weights in alp_build.", call. = FALSE)
  }

  S_beta_full <-
    t(p.sc * (1 - p.sc) * Xc) %*% Xc +
    t(sp[[wts.col]] * p.sp * (1 - p.sp) * Xp) %*% Xp

  D <- compute_D_ALP(sp_des, p.sp, Xp)

  out <- list(
    weights      = wts.sc,
    coefficients = as.numeric(beta),
    variables    = colnames(Xc),
    iterations   = iter,
    log_messages = log_messages,
    internal = list(
      Xc     = Xc,
      Xp     = Xp,
      p_sc   = p.sc,
      p_sp   = p.sp,
      S_beta = S_beta_full,
      D      = D
    )
  )

  return(out)
}
