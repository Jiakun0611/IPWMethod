#------------------------------------------------------------------------------#
# check_input_multi: validate inputs and optionally pre‐calibrate each sp_i to sp1
#------------------------------------------------------------------------------#
# Args:
#   sc              : sc dataframe
#   sp_list         : list of reference data.frames (length >= 1)
#   vars_list       : optional list of character vectors; if NULL, auto‐constructed by intersecting sc & each sp
#   wts_cols        : character vector of weight‐column names, same length as sp_list
# Returns:
#   list with elements:
#     sc        : user's non-probability sample data.frame
#     sp_list   : probability sample data.frames
#     vars_XC   : union of xcol lists (for the construction of Xc design matrix)
#     xcol      : list of per‐sample var sets (first = all shared, others = non‐dup)
#     wts_cols  : updated weight column after pre-calibration
check_input_multi <- function(sc,
                              sp_list,
                              vars_list,
                              wts_cols,
                              verbose = FALSE) {

  log_messages <- character()

  # sample names (use indices if unnamed)
  spn <- names(sp_list)
  if (is.null(spn)) spn <- paste0("sp_list[[", seq_along(sp_list), "]]")


  # ---- ensure all three inputs have the same length ----
  if (length(sp_list) != length(vars_list) ||
      length(sp_list) != length(wts_cols)) {
    stop("Error: sp_list, vars_list, and wts_cols must all be the same length.")
  }

  # ---- make all column names syntactically valid ----
  colnames(sc) <- make.names(colnames(sc))
  for (i in seq_along(sp_list)) {
    colnames(sp_list[[i]]) <- make.names(colnames(sp_list[[i]]))
  }

  # ============================================================
  #              PRINT SHARED VARIABLES PER SAMPLE
  # ============================================================
  for (i in seq_along(sp_list)) {
    msg <- sprintf(
      "Shared variables in %s:\n  %s\n",
      spn[i],
      paste(vars_list[[i]], collapse = ", ")
    )
    log_messages <- c(log_messages, msg)
    if (verbose) cat(msg)
  }


  # ============================================================
  #             BUILD xcol (NEW VARIABLES PER SAMPLE)
  # ============================================================
  n    <- length(vars_list)
  xcol <- vector("list", n)
  xcol[[1]] <- vars_list[[1]]

  if (n > 1) {
    for (i in seq.int(2, n)) {
      prev_vars <- Reduce(union, xcol[1:(i - 1)])
      new_vars  <- setdiff(vars_list[[i]], prev_vars)

      if (length(new_vars) == 0) {
        stop(sprintf("%s has no new variables; please remove this sample.", spn[i]))
      }
      xcol[[i]] <- new_vars
    }
  }

  for (i in seq_along(xcol)) {
    msg <- sprintf(
      "Variables used for calculation in %s:\n  %s\n",
      spn[i],
      paste(xcol[[i]], collapse = ", ")
    )
    log_messages <- c(log_messages, msg)
    if (verbose) cat(msg)
  }


  # ---- union of all xcol sets ----
  vars_XC <- Reduce(union, xcol)

  # ---- return validated objects ----
  list(
    sc       = sc,
    sp_list  = sp_list,
    vars_XC  = vars_XC,
    xcol     = xcol,
    wts_cols = wts_cols,
    log      = log_messages
  )
}
