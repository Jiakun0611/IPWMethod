#-------------------------------------------------------------
# Structural validation for IPWM inputs (build stage)
#     - no check on y and zcol
#     - sp_order / precali ONLY allowed and required in multi-reference case
#-------------------------------------------------------------
check_ipwm_inputs_build <- function(data, p_formula, method,
                                    sp_order = NULL, precali = NULL) {

  #--------------------------------------------------#
  # helper: extract variables from formula
  #--------------------------------------------------#
  get_formula_vars <- function(fml) {
    vars <- all.vars(fml)
    vars <- setdiff(vars, c(".", "1"))
    vars
  }

  #--------------------------------------------------#
  # helper: check variables in one formula
  #--------------------------------------------------#
  check_formula_vars_exist <- function(fml, sc, sp, label = "p_formula") {
    vars <- get_formula_vars(fml)

    miss_sc <- setdiff(vars, names(sc))
    miss_sp <- setdiff(vars, names(sp))

    if (length(miss_sc) > 0L) {
      stop(
        sprintf(
          "Variables in %s not found in `sc`: %s",
          label, paste(miss_sc, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    if (length(miss_sp) > 0L) {
      stop(
        sprintf(
          "Variables in %s not found in corresponding `sp`: %s",
          label, paste(miss_sp, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  #--------------------------------------------------#
  # data: basic structure
  #--------------------------------------------------#
  if (!is.list(data) || length(data) < 2L) {
    stop("`data` must be a list like list(sc, sp1.des, sp2.des, ...).",
         call. = FALSE)
  }

  sc <- data[[1]]
  sp_des <- data[-1]
  n_ref <- length(sp_des)

  if (!is.data.frame(sc)) {
    stop("The first element of `data` must be a data.frame for `sc`.",
         call. = FALSE)
  }

  if (nrow(sc) == 0L) {
    stop("`sc` has zero rows.", call. = FALSE)
  }

  if (anyDuplicated(names(sc))) {
    stop("`sc` has duplicated column names.", call. = FALSE)
  }

  ok_des <- vapply(
    sp_des,
    function(x) inherits(x, c("survey.design2", "svyrep.design")),
    logical(1)
  )

  if (!all(ok_des)) {
    bad <- which(!ok_des)
    stop(
      sprintf(
        "Elements %s of `data` are not valid survey design objects.",
        paste(bad + 1L, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  # stop when lonely psu+ getOption("survey.lonely.psu")="fail"
  for (i in seq_along(sp_des)) {
    if (inherits(sp_des[[i]], "survey.design2")) {
      check_lonely_psu(
        sp_des[[i]],
        name = paste0("sp_des[[", i, "]]")
      )
    }
  }

  # extract analysis data from survey design objects
  sp <- lapply(sp_des, extract_analysis_data)

  #--------------------------------------------------#
  # method
  #--------------------------------------------------#
  valid_methods <- c("alp", "clw", "multi", "cali", "calibration")

  if (!is.null(method)) {

    if (!is.character(method) || length(method) != 1L || is.na(method)) {
      stop("`method` must be a single character string or NULL.",
           call. = FALSE)
    }

    method <- tolower(trimws(method))

    if (!(method %in% valid_methods)) {
      stop(
        sprintf(
          "Invalid method '%s'. Must be one of: alp, clw, multi, cali, calibration.",
          method
        ),
        call. = FALSE
      )
    }

    # normalize
    if (method == "cali") {
      method <- "calibration"
    }

    #--------------------------------------------------#
    # consistency check with number of references
    #--------------------------------------------------#
    if (n_ref == 1 && method == "multi") {
      stop(
        "`method = 'multi' is only allowed when multiple reference samples are provided.",
        call. = FALSE
      )
    }

    if (n_ref > 1 && method %in% c("alp", "clw", "calibration")) {
      stop(
        paste0(
          "method = '", method, "' is only valid for a single reference sample. ",
          "Use `method = 'multi' when multiple reference samples are provided."
        ),
        call. = FALSE
      )
    }
  }

  #--------------------------------------------------#
  # one-reference case
  #--------------------------------------------------#
  if (n_ref == 1L) {

    if (!is.null(p_formula)) {
      if (!inherits(p_formula, "formula")) {
        stop("For one-reference case, `p_formula` must be a formula or NULL.",
             call. = FALSE)
      }

      check_formula_vars_exist(
        fml   = p_formula,
        sc    = sc,
        sp    = sp[[1]],
        label = "`p_formula`"
      )
    }

    return(invisible(TRUE))
  }

  #--------------------------------------------------#
  # multi-reference case
  #--------------------------------------------------#
  if (n_ref >= 2L) {

    if (is.null(sp_order)) {
      stop("`sp_order` must be provided in multi-reference case.",
           call. = FALSE)
    }

    if (!is.character(sp_order) || length(sp_order) != 1L || is.na(sp_order) ||
        !(sp_order %in% c("size", "given"))) {
      stop("`sp_order` must be one of 'size' or 'given'.",
           call. = FALSE)
    }

    if (is.null(precali)) {
      stop("`precali` must be provided in multi-reference case.",
           call. = FALSE)
    }

    if (!is.logical(precali) || length(precali) != 1L || is.na(precali)) {
      stop("`precali` must be a single TRUE/FALSE value.",
           call. = FALSE)
    }

    if (!is.null(p_formula)) {
      if (!is.list(p_formula) || length(p_formula) != n_ref) {
        stop(
          sprintf(
            "For multi-reference case, `p_formula` must be a list of %d formulas.",
            n_ref
          ),
          call. = FALSE
        )
      }

      ok_formula <- vapply(
        p_formula,
        function(x) inherits(x, "formula"),
        logical(1)
      )

      if (!all(ok_formula)) {
        stop("All elements of `p_formula` must be formulas.",
             call. = FALSE)
      }

      for (i in seq_len(n_ref)) {
        check_formula_vars_exist(
          fml   = p_formula[[i]],
          sc    = sc,
          sp    = sp[[i]],
          label = sprintf("`p_formula[[%d]]`", i)
        )
      }
    }

    return(invisible(TRUE))
  }

  stop("Invalid input structure in `data`.", call. = FALSE)
}
