#' @title Estimate Mean and Variance from an IPWM Build Object
#'
#' @description
#' Computes the pseudo-weighted mean and its variance using a pre-built object
#' returned by \code{\link{IPWM_build}}. This function applies the second-layer
#' missing-data filtering for the outcome and optional domain variables, then
#' performs estimation under the chosen method (ALP/CLW/calibration/multi).
#'
#' @details
#' \strong{Missing data handling (layer 2).}
#' After pseudo-weights are created (implicitly via \code{build}), estimation of the mean
#' requires complete cases for the outcome \code{y} and optional domain variables \code{zcol}.
#' This second filter is stored in the returned object (again mirroring \code{lm}-style behavior).
#'
#' \strong{Input \code{build}.}
#' The \code{build} object encapsulates participation-model objects and any cached matrices
#' needed for the estimation and variance calculations. This design supports a clean
#' build/estimate workflow and avoids recomputing heavy intermediates.
#'
#' @param build An object of class \code{"IPWM_build"} returned by \code{\link{IPWM_build}}.
#' @param y A character string. Name of the outcome variable in the convenience sample \code{sc}.
#' @param zcol Optional character string. Domain variable used for subset estimation
#'   or stratified estimation.
#' @param na.action A function indicating how to handle missing values for \code{y} and \code{zcol}
#'   during estimation. Default is \code{stats::na.omit}.
#'
#' @return
#' An object of class \code{"IPWM_estimate"} containing unweighted and IPWM-adjusted
#' point estimates, standard errors, and 95\% confidence intervals. The returned
#' object is a list with components:
#' \describe{
#'   \item{method}{Character. The method used, inherited from \code{build$method}.}
#'
#'   \item{unweighted}{A list of unweighted estimates computed from the original
#'     convenience sample \code{build$internal$raw_sc} (with optional domain subsetting):
#'     \describe{
#'       \item{mean}{Numeric. Unweighted mean of \code{y}.}
#'       \item{se}{Numeric. Standard error of the unweighted mean.}
#'       \item{CI_95}{Named numeric vector of length 2 with elements \code{lower} and \code{upper}
#'         giving the 95\% confidence interval.}
#'     }}
#'
#'   \item{adjusted}{A list of IPWM-adjusted estimates computed using pseudo-weights
#'     from the \code{build} object and the analysis data after outcome/domain NA handling:
#'     \describe{
#'       \item{mean}{Numeric. IPWM-adjusted mean of \code{y}.}
#'       \item{se}{Numeric. Standard error of the adjusted mean (\code{sqrt(variance)}).}
#'       \item{CI_95}{Named numeric vector of length 2 with elements \code{lower} and \code{upper}
#'         giving the 95\% confidence interval.}
#'     }}
#'
#'   \item{na}{A list storing missing-data handling information for outcome and optional
#'     domain variables produced by \code{process_na_yz()}. This records the estimation-layer
#'     NA filtering policy and indices (lm-style).}
#' }
#'
#' @seealso
#' \code{\link{IPWM}}, \code{\link{IPWM_build}}, \code{\link{summary.IPWM_estimate}},
#' \code{\link{print.IPWM_estimate}}
#'
#' @examples
#' data(sc)
#' data(ref_survey_1)
#' b <- IPWM_build(sc = sc, sp = ref_survey_1,
#'                sp_design = list(type = "poisson", weight = "wts"),
#'                p_formula = ~ agecat + marital + race,
#'                method = "ALP")
#' est <- IPWM_estimate(build = b, y = "psa")
#' est
#'
#' @export
IPWM_estimate <- function(
    build,                             # build: IPWM_build object
    y,                                 # y: outcome variable name in sc
    zcol = NULL,                       # zcol: optional domain variable(s)
    na.action = stats::na.omit          # na.action: NA handling for y/zcol during estimation
) {

  #------------------------------#
  # Step 1: input checking
  #------------------------------#
  tryCatch(
    check_ipwm_inputs_estimate(build, y, zcol),
    error = function(e)
      stop("Step 1 (input check) failed: ", e$message, call. = FALSE)
  )

  #------------------------------#
  # Step 2: prepare sc data
  #------------------------------#
  sc_data <- tryCatch(
    prepare_sc_data(build),
    error = function(e)
      stop("Step 2 (prepare_sc_data) failed: ", e$message, call. = FALSE)
  )

  #------------------------------#
  # Step 3: NA processing
  #------------------------------#
  yz_data <- tryCatch(
    process_na_yz(
      sc_data = sc_data,
      y = y,
      zcol = zcol,
      na.action = na.action
    ),
    error = function(e)
      stop("Step 3 (process_na_yz) failed: ", e$message, call. = FALSE)
  )

  #------------------------------#
  # Step 4: dispatch estimator
  #------------------------------#
  est <- tryCatch(
    dispatch_estimator(
      build = build,
      yz_data = yz_data
    ),
    error = function(e)
      stop("Step 4 (dispatch_estimator) failed: ", e$message, call. = FALSE)
  )

  #------------------------------#
  # Step 5: naive estimator
  #------------------------------#
  naive <- tryCatch(
    naive_mean(build$internal$raw_sc,
               domain_var = zcol,
               y = y),
    error = function(e)
      stop("Step 5 (naive_mean) failed: ", e$message, call. = FALSE)
  )

  #------------------------------#
  # Step 6: assemble output
  #------------------------------#
  result <- tryCatch(
    assemble_output(
      build = build,
      est = est,
      naive = naive,
      na_info = yz_data$na_info
    ),
    error = function(e)
      stop("Step 6 (assemble_output) failed: ", e$message, call. = FALSE)
  )

  result$call <- match.call()
  return(result)
}
