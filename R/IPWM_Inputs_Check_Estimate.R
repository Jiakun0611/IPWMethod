#-------------------------------------------------------------
# Estimation-stage validation (estimate stage)
#     - check y and zcol in raw_sc
#-------------------------------------------------------------
check_ipwm_inputs_estimate <- function(build, y, zcol = NULL) {


  # needed objects from build
  if (is.null(build$internal$raw_sc)) stop("build$internal$raw_sc is missing.", call. = FALSE)
  if (is.null(build$internal$na$keep_sc)) stop("build$internal$na$keep_sc is missing.", call. = FALSE)
  if (is.null(build$pseudo_weights)) stop("build$pseudo_weights is missing.", call. = FALSE)

  sc = build$internal$raw_sc
  colnames(sc) <- make.names(colnames(sc))

  # --- y ---
  if (!is.character(y) || length(y) != 1L)
    stop("'y' must be a single character string.", call. = FALSE)
  if (!(y %in% names(sc)))
    stop(sprintf("Outcome variable '%s' not found in sc.", y), call. = FALSE)
  if (!is.numeric(sc[[y]]))
    stop("Outcome variable 'y' must be numeric.", call. = FALSE)

  # --- zcol ---
  if (!is.null(zcol)) {

    if (!is.character(zcol) || length(zcol) != 1L || is.na(zcol) || !nzchar(zcol)) {
      stop("'zcol' must be a single non-empty character string or NULL.", call. = FALSE)
    }

    if (!(zcol %in% names(sc))) {
      stop(sprintf("Domain variable '%s' not found in 'sc'.", zcol), call. = FALSE)
    }

    if (anyNA(sc[[zcol]])) {
      stop(sprintf(
        "Domain variable '%s' contains missing values. It must consist only of 0 and 1.",
        zcol
      ), call. = FALSE)
    }

    if (!is.numeric(sc[[zcol]]) && !is.integer(sc[[zcol]])) {
      stop(sprintf(
        "Domain variable '%s' must be numeric or integer containing 0 and 1. Factors are not allowed.",
        zcol
      ), call. = FALSE)
    }

    uniq_vals <- sort(unique(sc[[zcol]]))
    if (!all(uniq_vals %in% c(0, 1))) {
      stop(sprintf(
        "Domain variable '%s' must contain only {0, 1}. Found values: %s",
        zcol,
        paste(uniq_vals, collapse = ", ")
      ), call. = FALSE)
    }
  }

  invisible(TRUE)
}
