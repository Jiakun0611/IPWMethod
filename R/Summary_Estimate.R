#' Summary method for IPWM_estimate objects
#'
#' Provides console output for objects of class \code{"IPWM_estimate"},
#' including unweighted and pseudo-weighted (adjusted) mean estimates,
#' standard errors, confidence intervals, and NA-handling information.
#'
#' @param object An object of class \code{"IPWM_estimate"}, typically the output
#'   of \code{\link{IPWM_estimate}}.
#' @param ... Additional arguments (not used).
#'
#' @method summary IPWM_estimate
#' @export
summary.IPWM_estimate <- function(object, ...) {

  # --- Call (optional) ---
  if (!is.null(object$call)) {
    cat("Call:\n")
    print(object$call)
  } else {
    cat("Call:\n")
    cat("  (not available)\n")
  }

  method <- tolower(trimws(object$method))

  if (method %in% c("alp", "clw", "cali", "calibration")) {
    if (method == "alp") method <- "ALP"
    if (method == "clw") method <- "CLW"
    if (method == "cali") method <- "calibration"

    cat("\nMethod: One reference", method, "\n")
  } else if (method == "multi") {
    cat("\nMethod: Multi reference calibration\n")
  }


  # --- Unweighted section ---
  cat("\nUnweighted Estimators summary:\n")

  uw_mean <- object$unweighted$mean
  uw_se   <- object$unweighted$se
  uw_ci   <- object$unweighted$CI_95

  if (!is.null(uw_mean)) cat(sprintf("  %-22s %10.6f\n", "Mean:", uw_mean))
  if (!is.null(uw_se))   cat(sprintf("  %-22s %10.6f\n", "Std. Error:", uw_se))
  if (!is.null(uw_ci) && length(uw_ci) >= 2) {
    cat(sprintf("  %-22s [%0.6f,%0.6f]\n",
                "95% CI:",
                uw_ci[1], uw_ci[2]))
  }

  # --- Adjusted section ---
  cat("\nPseudo-weighted (", method, ") Estimators summary:\n", sep = "")

  ad_mean <- object$adjusted$mean
  ad_se   <- object$adjusted$se
  ad_ci   <- object$adjusted$CI_95

  if (!is.null(ad_mean)) cat(sprintf("  %-22s %10.6f\n", "Mean:", ad_mean))
  if (!is.null(ad_se))   cat(sprintf("  %-22s %10.6f\n", "Std. Error:", ad_se))
  if (!is.null(ad_ci) && length(ad_ci) >= 2) {
    cat(sprintf("  %-22s [%0.6f,%0.6f]\n",
                "95% CI:",
                ad_ci[1], ad_ci[2]))
  }
  # invisible(object)
}
