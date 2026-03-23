prepare_sc_data <- function(build) {

  raw_sc  <- build$internal$raw_sc
  keep_sc <- build$internal$na$keep_sc

  if (is.null(raw_sc))  stop("prepare_sc_data: build$internal$raw_sc is NULL.", call. = FALSE)
  if (is.null(keep_sc)) stop("prepare_sc_data: build$internal$na$keep_sc is NULL.", call. = FALSE)

  idx_keep <- if (is.logical(keep_sc)) which(keep_sc) else as.integer(keep_sc)
  sc_keep  <- raw_sc[idx_keep, , drop = FALSE]

  X <- build$internal$Xc
  w <- build$pseudo_weights
  w <- w[!is.na(w)]

  if (is.null(X)) stop("prepare_sc_data: build$Xc is NULL (build stage did not store Xc).", call. = FALSE)
  if (is.null(w)) stop("prepare_sc_data: build$pseudo_weights is NULL.", call. = FALSE)

  nx <- nrow(X)
  if (is.null(nx)) stop("prepare_sc_data: nrow(build$Xc) is NULL (Xc is not a matrix/data.frame).", call. = FALSE)

  if (nx != length(w)) {
    stop("Mismatch between nrow(build$Xc) and length(build$pseudo_weights).", call. = FALSE)
  }
  if (nrow(sc_keep) != length(w)) {
    stop("Mismatch between nrow(sc_keep) and length(build$pseudo_weights).", call. = FALSE)
  }

  list(sc = sc_keep, X = X, w = w, idx_keep = idx_keep)
}
