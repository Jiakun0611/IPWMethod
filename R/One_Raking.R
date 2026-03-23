raking <- function(y, vars, sc, sp, wts.col, zcol,
                   design,
                   maxit = 20, tol = 1e-4,
                   verbose = FALSE,
                   log_messages = NULL) {

  ### --------------------------------------------------------------------
  ### 0. Initialize log system
  ### --------------------------------------------------------------------
  if (is.null(log_messages)) log_messages <- character(0)

  add_log <- function(msg) {
    log_messages <<- c(log_messages, msg)
    if (verbose) message(msg)
  }

  ### --------------------------------------------------------------------
  ### 1. Construct design matrices for calibration
  ### --------------------------------------------------------------------
  Xc <- design_matrix(vars = vars, data = sc)
  Xp <- design_matrix(vars = vars, data = sp)

  p_dim <- ncol(Xc)
  beta  <- rep(0, p_dim)

  ### --------------------------------------------------------------------
  ### 2. Newton–Raphson estimation for beta (raking equations)
  ### --------------------------------------------------------------------
  iter <- 0
  crit <- TRUE

  f.p <- colSums(sp[[wts.col]] * Xp)

  # stabilize initial beta[1]
  wts0 <- as.vector(exp(-Xc %*% beta))
  f.c0 <- colSums(wts0 * Xc)
  beta[1] <- -log(f.p[1] / f.c0[1])

  while (crit) {

    wts.sc <- as.vector(exp(-Xc %*% beta))

    f.c <- colSums(wts.sc * Xc)
    f   <- f.c - f.p

    fprime <- - t(wts.sc * Xc) %*% Xc

    delta <- solve(fprime, f)
    beta  <- beta - delta

    iter <- iter + 1
    crit <- (iter < maxit && sum(abs(delta)) > tol)
  }

  ### --------------------------------------------------------------------
  ### 3. Compute pseudo-weights
  ### --------------------------------------------------------------------
  wts.sc <- as.vector(exp(-Xc %*% beta))

  ### --------------------------------------------------------------------
  ### 4. Remove NA in outcome variable y
  ### --------------------------------------------------------------------
  ok_y <- !is.na(sc[[y]])
  omitted_y <- sum(!ok_y)

  if (omitted_y > 0) {
    add_log(sprintf(
      "%d observations in sc have complete p_formula covariates but missing outcome values ('%s'), which are excluded from the pseudo-weighted estimation.\n",
      omitted_y, y
    ))
  }

  sc_ok    <- sc[ok_y, ]
  Xc_ok    <- Xc[ok_y, , drop = FALSE]
  w_sc_ok  <- wts.sc[ok_y]

  ### --------------------------------------------------------------------
  ### 5. Construct domain indicator Z and remove NA in Z
  ### --------------------------------------------------------------------
  if (is.null(zcol)) {
    Z <- rep(1, length(w_sc_ok))
  } else {
    Z <- sc_ok[[zcol]]
  }

  ok_Z <- !is.na(Z)
  omitted_Z <- sum(!ok_Z)

  if (omitted_Z > 0) {
    add_log(sprintf(
      "For rows with complete p-formula variables in sc, %d rows have NA in domain variable '%s', which are ignored in pseudo-weighted estimator compuation.\n",
      omitted_Z, zcol
    ))
  }

  # Apply Z filtering
  Z       <- Z[ok_Z]
  Y       <- sc_ok[[y]][ok_Z] * Z
  w_sc_ok <- w_sc_ok[ok_Z]
  Xc_ok   <- Xc_ok[ok_Z, , drop = FALSE]

  ### --------------------------------------------------------------------
  ### 6. Raking mean estimator
  ### --------------------------------------------------------------------
  T1 <- sum(Y * w_sc_ok)
  T2 <- sum(Z * w_sc_ok)
  mu <- T1 / T2

  ### --------------------------------------------------------------------
  ### 7. Raking variance estimator
  ### --------------------------------------------------------------------
  U_beta <- t(Y - mu * Z) %*% (w_sc_ok * Xc_ok)
  S_beta <- t(w_sc_ok * Xc_ok) %*% Xc_ok
  b_vec <- U_beta %*% solve(S_beta)

  # Convenience part
  v1 <- sum(
    w_sc_ok * (w_sc_ok - 1) *
      (Y - mu * Z - (Xc_ok %*% t(b_vec)))^2
  )


  # Reference part
  D  <- compute_D_raking(sp = sp, Xp = Xp, wts.col = wts.col, design = design)
  v2 <- b_vec %*% D %*% t(b_vec)

  variance <- as.vector((v1 + v2) / T2^2)

  ### --------------------------------------------------------------------
  ### 8. Return
  ### --------------------------------------------------------------------
  return(list(
    wts_raking     = wts.sc,
    beta_raking    = as.numeric(beta),
    mean_raking    = mu,
    variance_raking = variance,
    names          = colnames(Xc),
    iterations     = iter,
    omitted_yNA    = omitted_y,
    omitted_ZNA    = omitted_Z,
    log_messages   = log_messages
  ))
}
