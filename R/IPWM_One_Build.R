IPWM_One_build <- function(
    sc,
    sp,
    sp_des,
    vars = NULL,
    weight,
    method,
    maxit = 20, tol = 1e-4,
    verbose = FALSE,
    log_messages = NULL,
    wname = "pseudo_wts"
) {

  # ---------------------------------------------------------------
  # 1. Dispatch to *_build
  #    Each build function must return:
  #    weights, variables, coefficients, iterations, log_messages
  # ---------------------------------------------------------------
  m_lower <- tolower(method)

  if (m_lower == "alp") {

    built <- alp_build(
      sc = sc,
      sp = sp,
      sp_des = sp_des,
      vars = vars,
      wts.col = weight,
      maxit = maxit, tol = tol,
      verbose = verbose,
      log_messages = log_messages
    )

  } else if (m_lower == "clw") {

    built <- clw_build(
      sc = sc,
      sp = sp,
      sp_des = sp_des,
      vars = vars,
      wts.col = weight,
      maxit = maxit, tol = tol,
      verbose = verbose,
      log_messages = log_messages
    )

  } else if (m_lower %in% c("calibration", "cali", "raking")) {

    built <- raking_build(
      sc = sc,
      sp = sp,
      sp_des = sp_des,
      vars = vars,
      wts.col = weight,
      maxit = maxit, tol = tol,
      verbose = verbose,
      log_messages = log_messages
    )

  } else {
    stop("Unknown method. Use method = 'ALP', 'CLW', or 'calibration'.",
         call. = FALSE)
  }


  # ---------------------------------------------------------------
  # 4. Return ONLY build objects
  # ---------------------------------------------------------------
  out <- list(
    pseudo_weights = built$weights,
    variables = built$variables,
    coefficients = built$coefficients,
    NR_iterations = built$iterations,
    method = method,
    log_messages = built$log_messages,
    internal = built$internal
  )
  return(out)
}
