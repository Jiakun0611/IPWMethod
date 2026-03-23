process_p_formula <- function(
    sc, sp, weight, p_formula,
    Pre.calibration = TRUE,
    sp_order = "size",
    verbose = FALSE
) {

  # -------------------------------------------------------------------
  # One-reference case
  # -------------------------------------------------------------------
  if (is.data.frame(sp)) {

    if (!inherits(p_formula, "formula")) {
      stop("'p_formula' must be a formula for one-reference case.", call. = FALSE)
    }

    vars_in_formula <- all.vars(p_formula)
    missing_sc <- setdiff(vars_in_formula, names(sc))
    missing_sp <- setdiff(vars_in_formula, names(sp))
    if (length(missing_sc) > 0 || length(missing_sp) > 0) {
      stop(paste0(
        "Missing variable(s): ",
        paste(unique(c(missing_sc, missing_sp)), collapse = ", ")
      ), call. = FALSE)
    }

    # Build model matrices (drop implicit intercept)
    Xc <- stats::model.matrix(p_formula, data = sc)
    Xp <- stats::model.matrix(p_formula, data = sp)
    if ("(Intercept)" %in% colnames(Xc)) Xc <- Xc[, -1, drop = FALSE]
    if ("(Intercept)" %in% colnames(Xp)) Xp <- Xp[, -1, drop = FALSE]

    vars <- make.names(colnames(Xc), unique = TRUE)
    colnames(Xc) <- vars
    colnames(Xp) <- vars

    sc_new <- as.data.frame(Xc)

    if (!(weight %in% names(sp))) stop("Weight column not found in sp.", call. = FALSE)
    sp_new <- as.data.frame(Xp)
    sp_new[[weight]] <- sp[[weight]]


    return(list(
      sc = sc_new,
      sp = sp_new,
      vars = vars,
      log_messages = character()
    ))
  }

  # -------------------------------------------------------------------
  # Multi-reference case
  # -------------------------------------------------------------------
  if (!is.list(sp)) {
    stop("process_p_formula_build(): sp must be a data.frame or a list of data.frames.", call. = FALSE)
  }

  if (!(is.list(p_formula) && all(vapply(p_formula, inherits, logical(1), "formula")))) {
    stop("For multi-reference, 'p_formula' must be a list of formulas.", call. = FALSE)
  }
  if (!(is.character(weight) && length(weight) == length(sp))) {
    stop("For multi-reference, 'weight' must be a character vector with length = length(sp).", call. = FALSE)
  }


  sp_new <- vector("list", length(sp))
  vars_list <- vector("list", length(sp))
  sc_new <- NULL

  for (j in seq_along(sp)) {

    fml <- p_formula[[j]]
    vars_in_formula <- all.vars(fml)

    missing_sc <- setdiff(vars_in_formula, names(sc))
    missing_sp <- setdiff(vars_in_formula, names(sp[[j]]))
    if (length(missing_sc) > 0 || length(missing_sp) > 0) {
      stop(paste0(
        "For reference ", j, ", missing vars: ",
        paste(unique(c(missing_sc, missing_sp)), collapse = ", ")
      ), call. = FALSE)
    }

    Xc <- stats::model.matrix(fml, data = sc)
    Xp <- stats::model.matrix(fml, data = sp[[j]])
    if ("(Intercept)" %in% colnames(Xc)) Xc <- Xc[, -1, drop = FALSE]
    if ("(Intercept)" %in% colnames(Xp)) Xp <- Xp[, -1, drop = FALSE]

    colnames(Xc) <- make.names(colnames(Xc), unique = TRUE)
    colnames(Xp) <- make.names(colnames(Xp), unique = TRUE)

    # Merge covariates into sc_new without duplicates
    if (j == 1L) {
      sc_new <- as.data.frame(Xc)
    } else {
      for (nm in colnames(Xc)) {
        new_col <- Xc[, nm]
        duplicate <- any(vapply(sc_new, function(old_col) {
          if (is.numeric(old_col) && is.numeric(new_col)) {
            isTRUE(all.equal(old_col, new_col, tolerance = 1e-12))
          } else {
            identical(old_col, new_col)
          }
        }, logical(1)))
        if (!duplicate) sc_new[[nm]] <- new_col
      }
    }

    # Build sp_new[[j]]
    wj <- weight[j]
    if (!(wj %in% names(sp[[j]]))) {
      stop(sprintf("Weight column '%s' not found in sp[%d].", wj, j), call. = FALSE)
    }

    spj <- as.data.frame(Xp)
    spj[[wj]] <- sp[[j]][[wj]]

    sp_new[[j]] <- spj
    vars_list[[j]] <- colnames(Xp)
  }

  names(sp_new) <- if (!is.null(names(sp)) && all(nzchar(names(sp)))) {
    names(sp)
  } else {
    paste0("sp[[", seq_along(sp), "]]")
  }

  log_messages <- character()

  if (Pre.calibration && length(sp_new) > 1) {
    out <- precal_cumulative_order(
      sp_raw   = sp,
      sp_new   = sp_new,
      weight   = weight,
      sp_order = sp_order,
      verbose  = verbose
    )
    sp_new <- out$sp_new
    log_messages <- c(log_messages, out$log_messages)
  } else {
    msg <- "Pre-calibration is recommended.\n"
    log_messages <- c(log_messages, msg)
    if (verbose) cat(msg)
  }

  return(list(
    sc = sc_new,
    sp = sp_new,
    vars = vars_list,
    log_messages = log_messages
  ))
}
