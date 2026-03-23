expit <- function(x){
  y <- exp(x) / (1 + exp(x))
  return(y)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

check_design_identifiability <- function(X, method = "method", tol = 1e-10) {
  qx <- qr(X, tol = tol)
  if (qx$rank < ncol(X)) {
    stop(
      sprintf(
        paste0("[%s] Model matrix is rank-deficient before NR iteration ",
               "(rank = %d < %d). This suggests collinearity or redundant terms ",
               "in p_formula. Please simplify the model."),
        method, qx$rank, ncol(X)
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}



