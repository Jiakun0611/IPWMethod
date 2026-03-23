# -------------------------------------------------------------------
# Construct participation model formula(s) automatically
# -------------------------------------------------------------------
p_formula_construction <- function(sc, sp, weight, verbose = TRUE) {

  log_messages <- character(0)

  add_log <- function(msg) {
    log_messages <<- c(log_messages, msg)
    if (verbose) message(msg)
  }
# -------------------------------------------------------------------
# Construct participation model formula(s) automatically
# -------------------------------------------------------------------
p_formula_construction <- function(sc, sp, weight, verbose = TRUE) {

  #-------------------------------#
  # helper: build one formula
  #-------------------------------#
  build_one_formula <- function(sc, sp_i, weight_i, sp_name) {

    shared <- intersect(colnames(sc), colnames(sp_i))
    drop_these <- unique(c(weight_i))
    vars <- setdiff(shared, drop_these)

    if (length(vars) == 0L) {
      stop(
        sprintf("No shared covariates found to build default p_formula for %s.", sp_name),
        call. = FALSE
      )
    }

    fml <- as.formula(paste("~", paste(vars, collapse = " + ")))

    if (verbose) {
      message(
        "Generated default p_formula for ", sp_name, ": ",
        paste(deparse(fml), collapse = "")
      )
    }

    return(fml)
  }

  #-------------------------------#
  # one-reference case
  #-------------------------------#
  if (is.data.frame(sp)) {
    fml <- build_one_formula(
      sc       = sc,
      sp_i     = sp,
      weight_i = weight,
      sp_name  = "reference survey"
    )

    return(fml)
  }

  #-------------------------------#
  # multi-reference case
  #-------------------------------#
  if (is.list(sp)) {
    if (!is.character(weight) || length(weight) != length(sp)) {
      stop(
        "For multi-reference case, `weight` must be a character vector with same length as `sp`.",
        call. = FALSE
      )
    }

    p_formula_list <- vector("list", length(sp))

    for (i in seq_along(sp)) {
      sp_name <- if (!is.null(names(sp)) && nzchar(names(sp)[i])) {
        names(sp)[i]
      } else {
        paste0("sp[[", i, "]]")
      }

      p_formula_list[[i]] <- build_one_formula(
        sc       = sc,
        sp_i     = sp[[i]],
        weight_i = weight[i],
        sp_name  = sp_name
      )
    }

    names(p_formula_list) <- if (!is.null(names(sp))) {
      names(sp)
    } else {
      paste0("sp", seq_along(sp))
    }

    return(p_formula_list)
  }

  stop("`sp` must be either a data.frame or a list of data.frames.", call. = FALSE)
}
  #-------------------------------#
  # helper: build one formula
  #-------------------------------#
  build_one_formula <- function(sc, sp_i, weight_i, sp_name) {
    log_messages <- character(0)

    shared <- intersect(colnames(sc), colnames(sp_i))
    drop_these <- unique(c(weight_i))
    vars <- setdiff(shared, drop_these)

    if (length(vars) == 0L) {
      stop(
        sprintf("No shared covariates found to build default p_formula for %s.", sp_name),
        call. = FALSE
      )
    }

    fml <- as.formula(paste("~", paste(vars, collapse = " + ")))

    msg <- paste0(
      "Generated default p_formula for ", sp_name, ": ",
      paste(deparse(fml), collapse = ""),
      "\n"
    )

    log_messages <- c(log_messages, msg)

    return(list(
      p_formula    = fml,
      log_messages = log_messages
    ))
  }

  #-------------------------------#
  # one-reference case
  #-------------------------------#
  if (is.data.frame(sp)) {
    res <- build_one_formula(
      sc       = sc,
      sp_i     = sp,
      weight_i = weight,
      sp_name  = "reference survey"
    )

    if (!is.null(res$log_messages) && length(res$log_messages) > 0) {
      log_messages <- c(log_messages, res$log_messages)
    }

    return(list(
      p_formula    = res$p_formula,
      log_messages = log_messages
    ))
  }

  #-------------------------------#
  # multi-reference case
  #-------------------------------#
  if (is.list(sp)) {
    if (!is.character(weight) || length(weight) != length(sp)) {
      stop(
        "For multi-reference case, `weight` must be a character vector with same length as `sp`.",
        call. = FALSE
      )
    }

    p_formula_list <- vector("list", length(sp))

    for (i in seq_along(sp)) {
      sp_name <- if (!is.null(names(sp)) && nzchar(names(sp)[i])) {
        names(sp)[i]
      } else {
        paste0("sp[[", i, "]]")
      }

      res_i <- build_one_formula(
        sc       = sc,
        sp_i     = sp[[i]],
        weight_i = weight[i],
        sp_name  = sp_name
      )

      p_formula_list[[i]] <- res_i$p_formula

      # collect child logs as separate entries
      if (!is.null(res_i$log_messages) && length(res_i$log_messages) > 0) {
        log_messages <- c(log_messages, res_i$log_messages)
      }
    }

    names(p_formula_list) <- if (!is.null(names(sp))) {
      names(sp)
    } else {
      paste0("sp", seq_along(sp))
    }

    return(list(
      p_formula    = p_formula_list,
      log_messages = log_messages
    ))
  }
  stop("`sp` must be either a data.frame or a list of data.frames.", call. = FALSE)
}
