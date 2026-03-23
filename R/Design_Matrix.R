design_matrix <- function(vars, data, intercept = TRUE) {
  if (length(vars) == 0) {
    stop("No variables provided for design matrix.")
  }
  rhs <- paste(vars, collapse = " + ")
  f <- if (intercept) {
    as.formula(paste("~", rhs))
  } else {
    as.formula(paste("~ 0 +", rhs))
  }
  model.matrix(f, data = data)
}
