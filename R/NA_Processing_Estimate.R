process_na_yz <- function(sc_data, y, zcol, na.action) {

  sc  <- sc_data$sc
  X   <- sc_data$X
  w   <- sc_data$w
  idx_keep <- sc_data$idx_keep

  Y0 <- sc[[y]]
  Z0 <- if (is.null(zcol)) rep(1, length(Y0)) else sc[[zcol]]

  ok <- !is.na(Y0) & !is.na(Z0)

  if (identical(na.action, stats::na.fail) && any(!ok)) {
    stop("Missing values detected in y or zcol.", call. = FALSE)
  }

  if (identical(na.action, stats::na.pass)) {
    stop("na.pass not supported at estimate stage.", call. = FALSE)
  }

  Y <- Y0[ok]
  Z <- Z0[ok]
  X <- X[ok, , drop = FALSE]
  w <- w[ok]

  omitted_raw <- idx_keep[!ok]

  na_obj <- omitted_raw
  class(na_obj) <- if (identical(na.action, stats::na.exclude))
    "exclude" else "omit"

  list(
    Y = Y,
    Z = Z,
    X = X,
    w = w,
    na_info = list(
      na_action = na_obj,
      n_omitted = sum(!ok),
      n_used = length(Y)
    )
  )
}
