#' Two-sample pre-calibration function
#'
#' Performs pre-calibration between two probability samples by aligning
#' the marginal distributions of selected variables (\code{dup_vars})
#' in the second reference sample to those of the first one.
#'
#' @param wts A character vector of length 2, giving the names of the weight
#'   columns in the two reference datasets (\code{c(wt1, wt2)}).
#' @param refs A list of two data frames. \code{refs[[1]]} is the first reference
#'   sample, and \code{refs[[2]]} is the second reference sample to be calibrated.
#' @param dup_vars A character vector of variable names used for calibration.
#'   If empty, calibration is based only on the total weight.
#'
#' @return A numeric vector of calibrated weights for \code{refs[[2]]}.
#' @details
#' The calibration is performed using functions from the \pkg{survey} package:
#' \code{\link[survey]{svydesign}}, \code{\link[survey]{svytotal}}, and
#' \code{\link[survey]{calibrate}}. The resulting weights align the margins of
#' the second reference sample with those of the first reference sample.
#'
#' @examples
#' \dontrun{
#' sp1 <- data.frame(x = rnorm(100), wt1 = runif(100, 0.5, 1.5))
#' sp2 <- data.frame(x = rnorm(120), wt2 = runif(120, 0.5, 1.5))
#' PRECALI(wts = c("wt1", "wt2"), refs = list(sp1, sp2), dup_vars = "x")
#' }
#'
#' @seealso [survey::svydesign], [survey::svytotal], [survey::calibrate]
#' @importFrom survey svydesign calibrate svytotal
#' @importFrom stats weights
#' @export
#'
PRECALI <- function(wts, refs, dup_vars) {

  # --- build survey design objects ---
  ds1 <- survey::svydesign(ids = ~1,
                           weights = as.formula(paste0("~", wts[1])),
                           data    = refs[[1]])
  ds2 <- survey::svydesign(ids = ~1,
                           weights = as.formula(paste0("~", wts[2])),
                           data    = refs[[2]])

  # --- calibration formula and population totals ---
  if (length(dup_vars) == 0) {
    fml  <- as.formula("~1")
    pops <- sum(weights(ds1))
    names(pops) <- "(Intercept)"
  } else {
    fml   <- as.formula(paste0("~", paste(dup_vars, collapse = " + ")))
    total <- sum(weights(ds1))
    marg  <- coef(survey::svytotal(fml, ds1))
    pops  <- c(`(Intercept)` = total, marg)
  }

  # --- calibrate ds2 to ds1 margins ---
  ds2_cal <- survey::calibrate(
    design     = ds2,
    formula    = fml,
    population = pops,
    calfun     = "raking"
  )

  # --- return new weights ---
  weights(ds2_cal)
}
