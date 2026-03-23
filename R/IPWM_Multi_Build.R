#' Build step for multi-reference raking (pseudo-weight construction only)
#'
#' Constructs pseudo-weights using multi-reference raking (Newton–Raphson),
#' and stores internal objects needed for the estimate/variance stage.
#'
#' @param sc Convenience sample data.frame.
#' @param sp List of probability samples (each a data.frame).
#' @param vars Either a vector of covariate names shared across all sp,
#'   or a list of vectors (one per sp). (Matched to your check_input_multi design.)
#' @param weight Weight column name(s) in sp. Either a single string (recycled) or
#'   a character vector/list aligned with sp.
#' @param design Multi-reference design info (list aligned with sp), used by make_block_D_multi().
#' @param sp_order Ordering rule for sp ("size" or "given", etc.).
#' @param maxit,tol NR controls.
#' @param verbose Print log messages.
#' @param log_messages Optional character vector of existing logs (for pipelining).
#'
#' @return A build object with unified structure (like alp_build/clw_build/raking_build):
#'   user-facing fields + internal objects for estimation.
IPWM_Multi_Raking_build <- function(
    sc, sp, vars, weight,
    sp_des, sp_order,
    maxit = 20, tol = 1e-4,
    verbose = FALSE, log_messages = NULL
) {

  #------------------------------------------------------------#
  # 0. Log system
  #------------------------------------------------------------#
  if (is.null(log_messages)) log_messages <- character(0)

  add_log <- function(msg) {
    log_messages <<- c(log_messages, msg)
    if (verbose) message(msg)
  }

  step_try <- function(step_num, step_name, expr) {
    tryCatch(
      expr,
      error = function(e) {
        stop(
          sprintf(
            "Step %s (%s) failed: %s",
            step_num, step_name, conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )
  }

  #------------------------------------------------------------#
  # 1. Sort sp
  #------------------------------------------------------------#
  sorted <- step_try("1", "sort sp", {
    sort_by_sp_size(
      sp       = sp,
      vars     = vars,
      weight   = weight,
      design   = sp_des,
      sp_order = sp_order,
      verbose  = verbose
    )
  })

  sp_sorted     <- sorted$sp
  vars_sorted   <- sorted$vars
  weight_sorted <- sorted$weight
  design_sorted <- sorted$design
  log_messages  <- c(log_messages, sorted$log)

  #------------------------------------------------------------#
  # 2. Validate inputs
  #------------------------------------------------------------#
  valid <- step_try("2", "validate inputs", {
    check_input_multi(
      sc        = sc,
      sp_list   = sp_sorted,
      vars_list = vars_sorted,
      wts_cols  = weight_sorted,
      verbose   = verbose
    )
  })

  sc_work      <- valid$sc
  sp_list      <- valid$sp_list
  vars_XC      <- valid$vars_XC
  xcol         <- valid$xcol
  wts_cols     <- valid$wts_cols
  log_messages <- c(log_messages, valid$log)

  #------------------------------------------------------------#
  # 3. Design matrices
  #------------------------------------------------------------#
  DM <- step_try("3", "construct design matrices", {
    Xc_Xp_Construction(
      vars_XC  = vars_XC,
      sc       = sc_work,
      sp_list  = sp_list,
      xcol     = xcol,
      wts_cols = wts_cols
    )
  })

  Xc       <- DM$Xc
  Xp_list  <- DM$Xp_list
  wts_list <- DM$wts_list

  #------------------------------------------------------------#
  # 4. NR solve
  #------------------------------------------------------------#
  nr_out <- step_try("4", "Newton-Raphson solve", {
    raking_nr(
      Xc       = Xc,
      Xp_list  = Xp_list,
      wts_list = wts_list,
      maxit    = maxit,
      tol      = tol
    )
  })

  beta <- nr_out$beta
  iter <- nr_out$iter

  if (iter >= maxit) {
    add_log("Multi-reference raking NR reached maxit without meeting tol.")
  }

  #------------------------------------------------------------#
  # 5. Compute pseudo-weights
  #------------------------------------------------------------#
  wts.sc <- step_try("5", "compute pseudo-weights", {
    out <- as.vector(1 / exp(Xc %*% beta))
    if (any(!is.finite(out))) {
      stop("Non-finite pseudo-weights detected.")
    }
    out
  })

  #------------------------------------------------------------#
  # 6. D matrix
  #------------------------------------------------------------#
  D <- step_try("6", "construct D matrix", {
    make_block_D_multi(
      sp_des_list = design_sorted,
      Xp_list     = Xp_list,
      lonely.psu  = "adjust"
    )
  })

  #------------------------------------------------------------#
  # 7. S_beta
  #------------------------------------------------------------#
  S_beta_full <- step_try("7", "construct S_beta", {
    t(wts.sc * Xc) %*% Xc
  })

  #------------------------------------------------------------#
  # 8. Return object
  #------------------------------------------------------------#
  out <- step_try("8", "build output object", {
    list(
      pseudo_weights = wts.sc,
      variables      = colnames(Xc),
      coefficients   = as.numeric(beta),
      NR_iterations  = iter,
      method         = "multi",
      log_messages   = log_messages,

      internal = list(
        Xc       = Xc,
        Xp_list  = Xp_list,
        wts_list = wts_list,
        D        = D,
        S_beta   = S_beta_full,

        xcol     = xcol,
        vars_XC  = vars_XC,
        wts_cols = wts_cols,
        design   = design_sorted,
        sp_order = sp_order
      )
    )
  })

  out
}
