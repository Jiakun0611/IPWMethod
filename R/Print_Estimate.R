#' Print method for IPWM_estimate objects
#'
#' Displays the pseudo-weighted (adjusted) mean estimate and its uncertainty.
#'
#' @param x An object of class \code{"IPWM_estimate"}.
#' @param ... Additional arguments (not used).
#'
#' @method print IPWM_estimate
#' @export
print.IPWM_estimate <- function(x, ...) {

  m <- trimws(as.character(x$method %||% ""))

  cat("\nPseudo-weighted (", m, ") Estimators:\n", sep = "")

  ad_mean <- x$adjusted$mean
  ad_se   <- x$adjusted$se
  ad_ci   <- x$adjusted$CI_95

  if (!is.null(ad_mean))
    cat(sprintf("  %-15s %10.6f\n", "Mean:", ad_mean))

  if (!is.null(ad_se))
    cat(sprintf("  %-15s %10.6f\n", "Std. Error:", ad_se))

  if (!is.null(ad_ci) && length(ad_ci) >= 2)
    cat(sprintf("  %-15s [%0.6f, %0.6f]\n",
                "95% CI:", ad_ci[1], ad_ci[2]))

  invisible(x)
}
