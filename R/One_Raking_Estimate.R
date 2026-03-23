raking_estimate <- function(Y, Z, w, X, D, S_beta) {

  T1 <- sum((Y * Z) * w)
  T2 <- sum(Z * w)
  mu <- T1 / T2

  U_beta <- t((Y * Z) - mu * Z) %*% (w * X)
  b_vec  <- U_beta %*% qr.solve(S_beta)

  v1 <- sum(w * (w - 1) *
              ((Y * Z) - mu * Z - (X %*% t(b_vec)))^2)

  v2 <- b_vec %*% D %*% t(b_vec)

  variance <- as.vector((v1 + v2) / T2^2)

  list(
    mean = mu,
    variance = variance
  )
}
