#' Build Objects for Inverse Propensity Weighting Methods
#'
#' @description
#' \code{IPWM_build()} constructs and validates all objects required for
#' inverse propensity weighting estimation using one or more reference surveys.
#' The function prepares the participation model, processes missing values in
#' participation variables, builds pseudo-weights, and stores all intermediate
#' quantities needed for downstream estimation and variance calculation.
#'
#' The returned object is intended to be passed to \code{\link{IPWM_estimate}}.
#'
#' @details
#' This function separates the \emph{build step} from the \emph{estimation step}.
#' It does not require the outcome of interest variable \code{y}. Instead, it focuses on
#' constructing pseudo-weights and saving method-specific internal objects.
#'
#' The input \code{data} must be a list whose first element is the
#' nonprobability sample \code{sc}, and whose remaining elements are reference
#' survey design objects. Internally, \code{IPWM_build()}:
#' \enumerate{
#'   \item parses and validates the input structure;
#'   \item determines whether the problem is a one-reference or multi-reference case;
#'   \item constructs a default participation model formula if \code{p_formula = NULL};
#'   \item processes missing values in variables used by the participation model;
#'   \item extracts and processes participation-model covariates;
#'   \item runs the requested build routine:
#'   \itemize{
#'     \item \code{IPWM_One_build()} for one-reference methods
#'       (\code{"alp"}, \code{"clw"}, \code{"calibration"});
#'     \item \code{IPWM_Multi_Raking_build()} for the multi-reference method
#'       (\code{"multi"}).
#'   }
#'   \item reconstructs the output convenience sample with the pseudo-weight column;
#'   \item stores NA-handling information and internal build objects for later use.
#' }
#'
#' \strong{One-reference versus multi-reference.}
#' If \code{data} contains one reference survey design object, then the method
#' must be one of \code{"alp"}, \code{"clw"}, or \code{"calibration"}
#' (with \code{"cali"} accepted as an alias).
#' If \code{data} contains more than one reference survey design object,
#' the method is resolved to \code{"multi"}.
#'
#' \strong{Missing-data handling.}
#' Missing values are handled using complete-case analysis only for variables involved in the participation
#' model. The selected NA action is recorded in the returned object, together
#' with the row index of convenience-sample observations retained for fitting. NA information can be
#' accessed using \code{"na.action()"}
#' This behavior is designed to mimic the style of \code{lm()}.
#'
#' \strong{Pseudo-weight reconstruction.}
#' The pseudo-weights are first estimated on the filtered convenience sample used in
#' the build step. They are then mapped back to the original convenience sample:
#' excluded rows can be shown with \code{na.action=na.exclude} or hided with \code{na.action=na.omit},
#' The final output data frame is stored in \code{sc_updated}.
#'
#' @param data A list of input data objects of the form
#'   \code{list(sc, sp1.des, sp2.des, ...)}.
#'   The first element must be a data frame corresponding to the convenience
#'   sample \code{sc}. Each remaining element must be a valid survey design object
#'   corresponding to a reference probability survey.
#'
#' @param sp_order Character string controlling the order of reference samples used for
#' pre-calibration and estimation in the multi-reference case. Either \code{"size"} or \code{"given"} is supported.
#' Default is \code{"size"}.
#'
#' @param precali Logical. Used only in the multi-reference case.
#'   If \code{TRUE}, cumulative pre-calibration is performed before the main
#'   multi-reference estimation step. Default is \code{TRUE}.
#'
#' @param p_formula Optional participation model formula.
#'   In the one-reference case, this should be a single formula such as
#'   \code{~ age + sex + income}.
#'   In the multi-reference case, this should be a list of formulas, one for each
#'   reference survey.
#'   If \code{NULL}, a default formula is constructed automatically from shared
#'   variables between \code{sc} and \code{sp}.
#'
#' @param method Character string specifying the pseudo-weighting method. Not case sensitive.
#'   Supported values are \code{"alp"}, \code{"clw"}, and \code{"calibration"}(\code{"cali"})
#'   in one-reference case, and automatically set as \code{"multi"} in multi-reference case.
#'
#'
#' @param na.action Function specifying how missing values should be handled for
#'   participation-model variables. Common options include
#'   \code{stats::na.omit}, \code{stats::na.exclude}, \code{stats::na.fail},
#'   and \code{stats::na.pass}. Default is \code{stats::na.omit}.
#'
#' @param sc_wname Character string giving the name of the pseudo-weight column
#'   to be added to the returned convenience sample. Default is \code{"pseudo_wts"}.
#'
#' @param maxit Integer. Maximum number of Newton-Raphson iterations for methods
#'   that require iterative estimation. Default is \code{20}.
#'
#' @param tol Numeric. Convergence tolerance for iterative estimation.
#'   Default is \code{1e-4}.
#'
#' @param verbose Logical. If \code{TRUE}, progress messages and diagnostics are
#'   printed during the build process. Default is \code{FALSE}.
#'
#' @return
#' An object of class \code{"IPWM_build"}.
#' This is a list containing user-facing outputs together with internal objects
#' required by \code{\link{IPWM_estimate}}. Important components include:
#'
#' \describe{
#'   \item{sc_updated}{
#'     A data frame containing the original convenience sample with an added
#'     pseudo-weight column named by \code{sc_wname}. Depending on the selected
#'     NA action, excluded rows may be removed or may receive \code{NA} weights.
#'   }
#'
#'   \item{pseudo_weights}{
#'     The pseudo-weight vector aligned with \code{sc_updated[[sc_wname]]}.
#'   }
#'
#'   \item{variables}{
#'     The participation-model variables actually used after formula processing.
#'   }
#'
#'   \item{coefficients}{
#'     Estimated participation-model coefficients.
#'   }
#'
#'   \item{NR_iterations}{
#'     Number of Newton-Raphson iterations used in the build step.
#'   }
#'
#'   \item{method}{
#'     The pseudo-weighting method used by the function.
#'   }
#'
#'   \item{log_messages}{
#'     A character vector of diagnostic and progress messages generated during
#'     formula construction, preprocessing, and pseudo-weight building.
#'   }
#'
#'   \item{internal}{
#'     A list of internal objects needed for downstream estimation and variance
#'     computation. This includes method-specific build outputs as well as
#'     metadata for NA handling and the original convenience sample.
#'   }
#'
#'   \item{call}{
#'     The matched function call.
#'   }
#' }
#'
#' @seealso
#' \code{\link{IPWM_estimate}},
#' \code{\link{print.IPWM_build}},
#' \code{\link{summary.IPWM_build}}
#'
#' @examples
#' ## Example structure only:
#' ## data <- list(sc, sp1_design)
#' ##
#' ## fit_b <- IPWM_build(
#' ##   data      = data,
#' ##   p_formula = ~ agecat + marital + race,
#' ##   method    = "alp"
#' ## )
#' ##
#' ## Multi-reference example:
#' ## data_multi <- list(sc, sp1_design, sp2_design)
#' ##
#' ## fit_b2 <- IPWM_build(
#' ##   data      = data_multi,
#' ##   p_formula = list(
#' ##     ~ agecat + marital + race,
#' ##     ~ height + weight + BMI
#' ##   ),
#' ##   sp_order  = "size",
#' ##   precali   = TRUE,
#' ##   method    = "multi"
#' ## )
#'
#' @export
IPWM_build <- function(
    data,                              # list(sc, sp1.des, sp2.des, ...)
    sp_order = "size",                 # ordering rule for multi sp
    precali = TRUE,                    # cumulative pre-calibration (Multi only)
    p_formula = NULL,                  # formula (One) or list of formulas (Multi)
    method = NULL,                     # "ALP"/"CLW"/"calibration"/"multi"
    na.action = stats::na.omit,        # NA handling for participation variables
    sc_wname = "pseudo_wts",           # output pseudo-weight column name in sc
    maxit = 20,                        # max Newton-Raphson iterations
    tol = 1e-4,                        # convergence tolerance
    verbose = FALSE                    # print intermediate updates
)






{

  #--------------------------------------------------------------------------#
  # Step 0. Parse input and basic validation
  #--------------------------------------------------------------------------#
  step0 <- tryCatch({

    check_ipwm_inputs_build(
      data      = data,
      p_formula = p_formula,
      method    = method,
      sp_order  = sp_order,
      precali   = precali
    )

    parsed <- parse_ipwm_data(data)

    sc     <- parsed$sc
    sp_des <- parsed$sp_des
    sp     <- parsed$sp_vars
    n_ref  <- parsed$n_ref

    if (n_ref == 1) {
      sp <- sp[[1]]
      sp_des <- sp_des[[1]]
    }

    weight <- rep("weights", n_ref)

    list(
      sc     = sc,
      sp_des = sp_des,
      sp     = sp,
      n_ref  = n_ref,
      weight = weight
    )

  }, error = function(e) {
    stop("Step 0 (input parsing and validation) failed: ", e$message, call. = TRUE)
  })

  sc     <- step0$sc
  sp_des <- step0$sp_des
  sp     <- step0$sp
  n_ref  <- step0$n_ref
  weight <- step0$weight

  # Keep user's original sc for output reconstruction
  sc0 <- sc

  #--------------------------------------------------------------------------#
  # Step 1. Detect single vs. multiple reference cases
  #--------------------------------------------------------------------------#
  method <- tryCatch({

    resolve_ipwm_method(sp, method)

  }, error = function(e) {
    stop("Step 1 (resolve method) failed: ", e$message, call. = TRUE)
  })

  #--------------------------------------------------------------------------#
  # Step 2. Auto-build and preprocess p_formula
  #--------------------------------------------------------------------------#
  step2 <- tryCatch({

    log_messages <- character(0)

    if (is.null(p_formula)) {
      built <- p_formula_construction(
        sc      = sc,
        sp      = sp,
        weight  = weight,
        verbose = verbose
      )
      p_formula <- built$p_formula
      log_messages <- c(log_messages, built$log_messages)
    }

    list(
      p_formula    = p_formula,
      log_messages = log_messages
    )

  }, error = function(e) {
    stop("Step 2 (p_formula construction) failed: ", e$message, call. = TRUE)
  })

  p_formula    <- step2$p_formula
  log_messages <- step2$log_messages

  #--------------------------------------------------------------------------#
  # Step 2.5. NA processing
  #--------------------------------------------------------------------------#
  step25 <- tryCatch({

    na_mode <- resolve_na_action(na.action)

    na_res <- handle_na_for_ipwm(
      sc        = sc,
      sp        = sp,
      p_formula = p_formula,
      na_mode   = na_mode
    )

    list(
      na_mode       = na_mode,
      sc            = na_res$sc,
      sp            = na_res$sp,
      keep_sc       = na_res$keep_sc,
      na_action_obj = na_res$na_action
    )

  }, error = function(e) {
    stop("Step 2.5 (NA processing) failed: ", e$message, call. = TRUE)
  })

  na_mode       <- step25$na_mode
  sc            <- step25$sc
  sp            <- step25$sp
  keep_sc       <- step25$keep_sc
  na_action_obj <- step25$na_action_obj

  #--------------------------------------------------------------------------#
  # Step 3. Process p_formula
  #--------------------------------------------------------------------------#
  step3 <- tryCatch({

    processed <- process_p_formula(
      sc              = sc,
      sp              = sp,
      weight          = weight,
      Pre.calibration = precali,
      p_formula       = p_formula,
      sp_order        = sp_order,
      verbose         = verbose
    )

    list(
      sc           = processed$sc,
      sp           = processed$sp,
      p_vars       = processed$vars,
      log_messages = processed$log_messages
    )

  }, error = function(e) {
    stop("Step 3 (process p_formula) failed: ", e$message, call. = TRUE)
  })

  sc     <- step3$sc
  sp     <- step3$sp
  p_vars <- step3$p_vars

  if (!is.null(step3$log_messages)) {
    log_messages <- c(log_messages, step3$log_messages)
  }

  #--------------------------------------------------------------------------#
  # Step 4. Estimation / build pseudo-weights
  #--------------------------------------------------------------------------#
  result <- tryCatch({

    if (tolower(method) %in% c("alp", "clw", "calibration", "cali", "raking")) {

      IPWM_One_build(
        sc           = sc,
        sp           = sp,
        sp_des       = sp_des,
        vars         = p_vars,
        weight       = weight,
        method       = method,
        maxit        = maxit,
        tol          = tol,
        verbose      = verbose,
        log_messages = log_messages
      )

    } else if (tolower(method) == "multi") {

      out <- IPWM_Multi_Raking_build(
        sc           = sc,
        sp           = sp,
        sp_des       = sp_des,
        vars         = p_vars,
        weight       = weight,
        maxit        = maxit,
        tol          = tol,
        sp_order     = sp_order,
        verbose      = verbose,
        log_messages = log_messages
      )
      out$method <- "multi"
      out

    } else {
      stop("Unknown method: ", method, call. = FALSE)
    }

  }, error = function(e) {
    stop("Step 4 (estimation/build) failed: ", e$message, call. = TRUE)
  })

  #--------------------------------------------------------------------------#
  # Step 5. Construct sc_updated based on na_mode
  #--------------------------------------------------------------------------#
  step5 <- tryCatch({

    if (sc_wname %in% names(sc0)) {
      stop(sprintf("Column '%s' already exists in sc.", sc_wname), call. = FALSE)
    }

    w_fit <- result$pseudo_weights
    if (!is.numeric(w_fit)) {
      w_fit <- as.numeric(w_fit)
    }

    if (na_mode %in% c("omit", "exclude", "fail")) {
      if (length(w_fit) != sum(keep_sc)) {
        stop("Length of pseudo_weights does not match sum(keep_sc).", call. = FALSE)
      }
    } else {
      if (length(w_fit) != nrow(sc0)) {
        stop("na.pass: pseudo_weights length does not match nrow(sc).", call. = FALSE)
      }
    }

    w_full <- rep(NA_real_, nrow(sc0))
    w_full[keep_sc] <- w_fit

    sc_out <- sc0
    sc_out[[sc_wname]] <- w_full

    if (na_mode == "omit") {
      sc_out <- sc_out[keep_sc, , drop = FALSE]
      if (!is.null(na_action_obj)) {
        attr(sc_out, "na.action") <- na_action_obj
      }
    } else if (na_mode == "exclude") {
      if (!is.null(na_action_obj)) {
        attr(sc_out, "na.action") <- na_action_obj
      }
    } else if (na_mode == "fail") {
      # if NA existed, already stopped
    } else if (na_mode == "pass") {
      # keep all rows
    }

    list(
      sc_out = sc_out
    )

  }, error = function(e) {
    stop("Step 5 (reconstruct sc_updated) failed: ", e$message, call. = TRUE)
  })

  sc_out <- step5$sc_out

  #--------------------------------------------------------------------------#
  # Step 6. Finalize and return
  #--------------------------------------------------------------------------#
  result <- tryCatch({

    result$sc_updated <- sc_out
    result$pseudo_weights <- sc_out[[sc_wname]]

    if (is.null(result$internal)) {
      result$internal <- list()
    }

    result$internal$raw_sc <- sc0
    result$internal$na <- list(
      na_mode   = na_mode,
      keep_sc   = keep_sc,
      na_action = na_action_obj
    )

    result$call <- match.call()
    class(result) <- "IPWM_build"

    result

  }, error = function(e) {
    stop("Step 6 (finalize output) failed: ", e$message, call. = TRUE)
  })

  invisible(result)
}
