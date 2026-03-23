#' Extract NA action from an IPWM_estimate object
#'
#' Returns the \code{na.action} component recorded during estimation,
#' mimicking \code{\link{stats}{na.action}} behavior for fitted model objects.
#'
#' @param object An object of class \code{"IPWM_estimate"}.
#' @param ... Additional arguments (not used).
#'
#' @method na.action IPWM_estimate
#' @export
na.action.IPWM_estimate <- function(object, ...) {
  object$na$na_action
}
