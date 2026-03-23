
CLW <- function(y, vars, sc, sp, wts.col, zcol,
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
  ### 1. Construct design matrices for the participation model
  ### --------------------------------------------------------------------
  Xc <- design_matrix(vars = vars, data = sc)
  Xp <- design_matrix(vars = vars, data = sp)

  p_dim <- ncol(Xc)
  beta  <- rep(0, p_dim)

  ### --------------------------------------------------------------------
  ### 2. Newton–Raphson estimation for beta (CLW equations)
  ### --------------------------------------------------------------------
  iter <- 0
  crit <- TRUE

  # stabilize initial beta[1]
  pi.sp <- as.vector(expit(Xp %*% beta))
  f.p   <- colSums(sp[[wts.col]] * pi.sp * Xp)
  f.c   <- colSums(Xc)
  beta[1] <- -log(f.p[1] / f.c[1])

  while (crit) {

    pi.sp <- as.vector(expit(Xp %*% beta))

    f.p <- colSums(sp[[wts.col]] * pi.sp * Xp)
    f.c <- colSums(Xc)
    f   <- f.c - f.p

    fprime <- - t(sp[[wts.col]] * pi.sp * (1 - pi.sp) * Xp) %*% Xp

    delta <- solve(fprime, f)
    beta  <- beta - delta

    iter <- iter + 1
    crit <- (iter < maxit && sum(abs(delta)) > tol)
  }

  ### --------------------------------------------------------------------
  ### 3. Compute pseudo-weights
  ### --------------------------------------------------------------------
  pi.sc  <- as.vector(expit(Xc %*% beta))
  pi.sp  <- as.vector(expit(Xp %*% beta))  # keep for variance pieces
  wts.sc <- 1 / pi.sc

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
  pi_sc_ok <- pi.sc[ok_y]

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
  Z        <- Z[ok_Z]
  Y        <- sc_ok[[y]][ok_Z] * Z
  w_sc_ok  <- w_sc_ok[ok_Z]
  Xc_ok    <- Xc_ok[ok_Z, , drop = FALSE]
  pi_sc_ok <- pi_sc_ok[ok_Z]

  ### --------------------------------------------------------------------
  ### 6. CLW mean estimator
  ### --------------------------------------------------------------------
  T1 <- sum(Y * w_sc_ok)
  T2 <- sum(Z * w_sc_ok)
  mu <- T1 / T2

  ### --------------------------------------------------------------------
  ### 7. CLW variance estimator
  ### --------------------------------------------------------------------
  # score contribution (CLW uses (w-1) here, keep your convention)
  U_beta <- t(Y - mu * Z) %*% ((w_sc_ok - 1) * Xc_ok)

  # Fisher information (CLW: sp part only, keep your convention)
  S_beta <- t(sp[[wts.col]] * pi.sp * (1 - pi.sp) * Xp) %*% Xp

  b_vec <- U_beta %*% solve(S_beta)

  # Reference part
  D <- t(sp[[wts.col]]^2 *
           (1 - 1 / sp[[wts.col]]) *
           (pi.sp * Xp)) %*% (pi.sp * Xp)

  # Convenience part
  v1 <- sum(
    w_sc_ok * (w_sc_ok - 1) *
      (Y - mu * Z - pi_sc_ok * (Xc_ok %*% t(b_vec)))^2
  )

  v2 <- b_vec %*% D %*% t(b_vec)

  variance <- as.vector((v1 + v2) / T2^2)

  ### --------------------------------------------------------------------
  ### 8. Return
  ### --------------------------------------------------------------------
  return(list(
    wts_CLW       = wts.sc,
    beta_CLW      = as.numeric(beta),
    mean_CLW      = mu,
    variance_CLW  = variance,
    names         = colnames(Xc),
    iterations    = iter,
    omitted_yNA   = omitted_y,
    omitted_ZNA   = omitted_Z,
    log_messages  = log_messages
  ))
}
