#' Estimate step for multi-reference raking (mean + variance)
#'
#' @param Y Outcome vector (length n_sc_used).
#' @param Z Domain/indicator vector (length n_sc_used). Use rep(1, n) if no domain.
#' @param w Pseudo-weights from build (length n_sc_used).
#' @param X Convenience-sample design matrix (n_sc_used x p), i.e., Xc after NA filtering.
#' @param D Block design-based variance matrix from multi reference samples (p x p).
#' @param S_beta Score/Information matrix from build: t(w * X) %*% X (p x p).
#'
#' @return list(mean, variance)
multi_estimate <- function(Y, Z, w, X, D, S_beta) {

  # --- mean ---
  T1 <- sum((Y * Z) * w)
  T2 <- sum(Z * w)
  mu <- T1 / T2

  # --- linearization piece for beta uncertainty ---
  U_beta <- t((Y * Z) - mu * Z) %*% (w * X)
  b_vec  <- U_beta %*% qr.solve(S_beta)   # 1 x p

  # --- residual for variance decomposition ---
  resid <- (Y * Z) - mu * Z - as.vector(X %*% t(b_vec))

  # v1: within-sc component (pseudo-weighted)
  v1 <- sum(w * (w - 1) * (resid^2))

  # v2: design-based component from reference surveys
  v2 <- as.numeric(b_vec %*% D %*% t(b_vec))

  variance <- as.numeric((v1 + v2) / (T2^2))

  list(
    mean     = mu,
    variance = variance
  )
}
