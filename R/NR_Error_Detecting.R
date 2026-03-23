solve_newton_step <- function(J, g, method = "method", tol_singular = 1e-10) {
  # J: Jacobian / derivative matrix
  # g: function value vector
  # solve J %*% delta = g

  if (!is.matrix(J)) {
    stop(sprintf("[%s] Internal error: Jacobian is not a matrix.", method),
         call. = FALSE)
  }

  if (!is.numeric(J) || !is.numeric(g)) {
    stop(sprintf("[%s] Internal error: Jacobian and score must be numeric.", method),
         call. = FALSE)
  }

  if (nrow(J) != ncol(J)) {
    stop(
      sprintf(
        paste0("[%s] Newton step failed: Jacobian is not square ",
               "(%d x %d). This usually indicates a mismatch in the number ",
               "of estimating equations and parameters."),
        method, nrow(J), ncol(J)
      ),
      call. = FALSE
    )
  }

  if (length(g) != nrow(J)) {
    stop(
      sprintf(
        paste0("[%s] Newton step failed: score vector length (%d) does not match ",
               "Jacobian dimension (%d)."),
        method, length(g), nrow(J)
      ),
      call. = FALSE
    )
  }

  if (any(!is.finite(J))) {
    stop(
      sprintf(
        paste0("[%s] Newton step failed: Jacobian contains NA/NaN/Inf. ",
               "This may be caused by extreme weights, separation, overflow, ",
               "or invalid covariate values."),
        method
      ),
      call. = FALSE
    )
  }

  if (any(!is.finite(g))) {
    stop(
      sprintf(
        paste0("[%s] Newton step failed: estimating-equation vector contains ",
               "NA/NaN/Inf. This may be caused by extreme weights, separation, ",
               "overflow, or invalid covariate values."),
        method
      ),
      call. = FALSE
    )
  }

  qJ <- qr(J, tol = tol_singular)
  rankJ <- qJ$rank
  p <- ncol(J)

  if (rankJ < p) {
    stop(
      sprintf(
        paste0("[%s] Newton step failed because the linear system is singular ",
               "or not identifiable (Jacobian rank = %d < %d).\n",
               "Possible causes:\n",
               "  1. redundant or perfectly collinear variables in p_formula;\n",
               "  2. categories with zero or near-zero support in sc or sp;\n",
               "  3. too many calibration constraints (variables) relative to the data;\n",
               "  4. near-deterministic participation model / quasi-separation.\n",
               "Try removing redundant variables, collapsing sparse factor levels, ",
               "or simplifying p_formula."),
        method, rankJ, p
      ),
      call. = FALSE
    )
  }

  delta <- tryCatch(
    qr.solve(qJ, g),
    error = function(e) {
      stop(
        sprintf(
          paste0("[%s] Newton step failed while solving the linear system.\n",
                 "Original error: %s\n",
                 "This usually means the estimating equations are numerically unstable ",
                 "or nearly singular. Try simplifying the model or removing collinear terms."),
          method, conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  if (any(!is.finite(delta))) {
    stop(
      sprintf(
        paste0("[%s] Newton step failed: solution increment contains NA/NaN/Inf. ",
               "The system is numerically unstable."),
        method
      ),
      call. = FALSE
    )
  }

  delta
}

