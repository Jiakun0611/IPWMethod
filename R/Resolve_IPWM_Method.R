#-------------------------------------------------------------
# Resolve IPWM method based on sp type (one vs multi reference)
#-------------------------------------------------------------
resolve_ipwm_method <- function(sp, method = NULL) {

  is_multi <- is.list(sp) && all(vapply(sp, is.data.frame, logical(1)))
  is_one   <- is.data.frame(sp)

  if (is_multi) {
    # Multi-reference case
    if (is.null(method)) {
      method <- "multi"
    } else if (!identical(method, "multi")) {
      warning("For list-type 'sp', method forcibly set to 'multi'.", call. = FALSE)
      method <- "multi"
    }
    return(method)
  }

  if (is_one) {
    # One-reference case
    if (is.null(method)) {
      stop(
        "Please specify method = 'ALP', 'CLW', or 'calibration'('Calibration'/'cali'/'Cali'/'CALI') for one-reference case.",
        call. = FALSE
      )
    }
    return(method)
  }

  stop("'sp' must be a data.frame (one reference) or a list of data.frames (multi reference).",
       call. = FALSE)
}
