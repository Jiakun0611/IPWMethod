ALP_get_D <- function(vars, sc, sp, wts.col,
                      maxit = 20, tol = 1e-4, verbose = FALSE) {

  add_log <- function(msg) if (verbose) message(msg)

  # Design matrices
  Xc <- design_matrix(vars = vars, data = sc)
  Xp <- design_matrix(vars = vars, data = sp)

  p_dim <- ncol(Xc)
  beta  <- rep(0, p_dim)

  # Newton–Raphson
  iter <- 0
  crit <- TRUE

  p.sp <- as.vector(expit(Xp %*% beta))
  p.sc <- as.vector(expit(Xc %*% beta))

  f.p <- colSums(sp[[wts.col]] * p.sp * Xp)
  f.c <- colSums((1 - p.sc) * Xc)
  beta[1] <- -log(f.p[1] / f.c[1])

  while (crit) {
    p.sp <- as.vector(expit(Xp %*% beta))
    p.sc <- as.vector(expit(Xc %*% beta))

    score <- colSums((1 - p.sc) * Xc) -
      colSums(sp[[wts.col]] * p.sp * Xp)

    Hmat <- t(p.sc * (p.sc - 1) * Xc) %*% Xc -
      t(sp[[wts.col]] * p.sp * (1 - p.sp) * Xp) %*% Xp

    delta <- solve(Hmat, score)
    beta  <- beta - delta

    iter <- iter + 1
    crit <- (iter < maxit && sum(abs(delta)) > tol)
  }

  if (iter >= maxit) add_log("Warning: hit maxit; beta may not have converged.")

  # Final p.sp
  p.sp <- as.vector(expit(Xp %*% beta))

  # Your original D (kept as-is)
  w <- sp[[wts.col]]
  D <- t(w^2 * (1 - 1 / w) * (p.sp * Xp)) %*% (p.sp * Xp)

  # Build H_df using p.sp (NOT w)
  H <- sweep(Xp, 1, p.sp, "*")
  H_df <- as.data.frame(H)
  colnames(H_df) <- make.names(paste0("h_", colnames(H_df)))

  sp_h <- cbind(sp, w = w, H_df)

  # Correct design: data must include weight column
  # des_sp <- svydesign(ids = ~1, weights = ~w, data = sp_h)
  des_sp = svydesign(ids = ~1, fpc = 1/w, data = sp_h, pps=poisson_sampling(sp_h[,"w"]))

  # Total and vcov
  fml <- reformulate(colnames(H_df))
  tot_h <- svytotal(fml, design = des_sp)
  D2 <- vcov(tot_h)

  return(list(D = D, D2 = D2, beta = beta, iter = iter))
}



vars=c('x1','x2','x3','x4')
D = ALP_get_D(vars = vars,sc=sc, sp=sp, wts.col = 'wt_sp1')
D$D2


library(survey)

Xp <- design_matrix(vars = vars, data = sp);head(Xp)

check_poisson_D <- function(sp, Xp, wts.col, tol = 1e-10) {
  d  <- sp[[wts.col]]        # d_i = 1/pi_i
  pi <- 1 / d                # pi_i

  ## --- explicit D (Poisson/Bernoulli) ---
  D_exp <- t((d^2 * (1 - 1/d)) * Xp) %*% Xp
  n = dim(sp)[1]
  z1 = d * Xp
  z  = sweep(z1, 2, colMeans(z1))
  D_exp  = n/(n-1) * t(z) %*% z



  ## --- put Xp into sp so svytotal can see them ---
  Xp_df <- as.data.frame(Xp)
  colnames(Xp_df) <- make.names(colnames(Xp_df))
  sp2 <- cbind(sp, pi = pi, Xp_df)

  ## --- survey design using inclusion probs ---
  des <- svydesign(ids = ~1, probs = ~pi, data = sp2)

  ## --- totals & vcov ---
  fml <- reformulate(colnames(Xp_df))   # ~ x1 + x2 + ...
  tot <- svytotal(fml, design = des)
  D_svy <- vcov(tot)

  list(
    max_abs_diff = max(abs(D_exp - D_svy)),
    all_equal    = isTRUE(all.equal(D_exp, D_svy, tolerance = tol)),
    D_exp = D_exp,
    D_svy = D_svy
  )
}
out <- check_poisson_D(sp, Xp, wts.col="wt_sp1")
out$max_abs_diff
out$all_equal
