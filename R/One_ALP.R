ALP <- function(y, vars, sc, sp, wts.col, zcol,
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
  ### 2. Newton–Raphson estimation for participation parameters beta
  ### --------------------------------------------------------------------
  iter <- 0
  crit <- TRUE

  p.sp <- as.vector(expit(Xp %*% beta))
  p.sc <- as.vector(expit(Xc %*% beta))

  # stabilize initial beta
  f.p <- colSums(sp[, wts.col] * p.sp * Xp)
  f.c <- colSums((1 - p.sc) * Xc)
  beta[1] <- -log(f.p[1] / f.c[1])

  while (crit) {

    p.sp <- as.vector(expit(Xp %*% beta))
    p.sc <- as.vector(expit(Xc %*% beta))

    score <- colSums((1 - p.sc) * Xc) -
      colSums(sp[, wts.col] * p.sp * Xp)

    H <- t(p.sc * (p.sc - 1) * Xc) %*% Xc -
      t(sp[, wts.col] * p.sp * (1 - p.sp) * Xp) %*% Xp

    delta <- solve(H, score)
    beta  <- beta - delta

    iter <- iter + 1
    crit <- (iter < maxit && sum(abs(delta)) > tol)
  }


  ### --------------------------------------------------------------------
  ### 3. Compute pseudo-weights
  ### --------------------------------------------------------------------
  p.sc <- as.vector(expit(Xc %*% beta))
  p.sp <- as.vector(expit(Xp %*% beta))

  pi.sc  <- p.sc / (1 - p.sc)
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

  sc_ok     <- sc[ok_y, ]
  Xc_ok     <- Xc[ok_y, , drop = FALSE]
  p_sc_ok   <- p.sc[ok_y]
  w_sc_ok   <- wts.sc[ok_y]


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
  p_sc_ok  <- p_sc_ok[ok_Z]


  ### --------------------------------------------------------------------
  ### 6. ALP mean estimator
  ### --------------------------------------------------------------------
  T1 <- sum(Y * w_sc_ok)
  T2 <- sum(Z * w_sc_ok)
  mu <- T1 / T2

  ### --------------------------------------------------------------------
  ### 7. ALP variance estimator
  ### --------------------------------------------------------------------
  # score contribution
  U_beta <- t(Y - mu * Z) %*% (w_sc_ok * Xc_ok)

  # Fisher information (sc filtered, sp unchanged)
  S_beta <-
    t(p_sc_ok * (1 - p_sc_ok) * Xc_ok) %*% Xc_ok +
    t(sp[, wts.col] * p.sp * (1 - p.sp) * Xp) %*% Xp

  b_vec <- U_beta %*% solve(S_beta)

  # Reference part
  D <- t(sp[, wts.col]^2 *
           (1 - 1 / sp[, wts.col]) *
           (p.sp * Xp)) %*% (p.sp * Xp)

  # Convenience part
  v1 <- sum(
    w_sc_ok * (w_sc_ok - 1) *
      (Y - mu * Z - p_sc_ok * (Xc_ok %*% t(b_vec)))^2
  )

  v2 <- b_vec %*% D %*% t(b_vec)

  variance <- as.vector((v1 + v2) / T2^2)


  ### --------------------------------------------------------------------
  ### 8. Return
  ### --------------------------------------------------------------------
  return(list(
    wts_ALP       = wts.sc,
    beta_ALP      = as.numeric(beta),
    mean_ALP      = mu,
    variance_ALP  = variance,
    names         = colnames(Xc),
    iterations    = iter,
    omitted_yNA   = omitted_y,
    omitted_ZNA   = omitted_Z,
    log_messages  = log_messages
  ))
}
