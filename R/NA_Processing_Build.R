handle_na_for_ipwm <- function(sc, sp, p_formula, na_mode = c("omit", "exclude", "fail", "pass")) {

  na_mode <- match.arg(na_mode)
  log <- character(0)

  #----------------------------------------------------------
  # helper: NA summary logging
  #----------------------------------------------------------
  .na_summary <- function(df, vars, label) {
    if (length(vars) == 0) return(invisible(NULL))
    na_ct <- colSums(is.na(df[, vars, drop = FALSE]))
    na_ct <- na_ct[na_ct > 0]
    if (length(na_ct) == 0) return(invisible(NULL))

    log <<- c(log, "\nMissing value summary:\n")
    log <<- c(log, sprintf("%d NA detected in p_formula variables in %s:", sum(na_ct), label))
    for (v in names(na_ct)) log <<- c(log, sprintf("  - %s", v))
    invisible(NULL)
  }

  #----------------------------------------------------------
  # 1. sc: variables used for filtering
  #    (if p_formula is list, union all.vars)
  #----------------------------------------------------------
  if (is.list(p_formula)) {
    vars_sc <- unique(unlist(lapply(p_formula, all.vars)))
  } else {
    vars_sc <- all.vars(p_formula)
  }

  # sc keep indicator
  if (length(vars_sc) > 0) {
    .na_summary(sc, vars_sc, "sc")
    keep_sc <- complete.cases(sc[, vars_sc, drop = FALSE])
  } else {
    keep_sc <- rep(TRUE, nrow(sc))
  }
  removed_sc <- sum(!keep_sc)

  # fail: stop if any NA
  if (na_mode == "fail" && removed_sc > 0) {
    stop("na.fail: NA found in p_formula variables in sc.", call. = FALSE)
  }

  # pass: do nothing at all
  if (na_mode == "pass") {
    keep_sc[] <- TRUE
  } else if (na_mode %in% c("omit", "exclude") && removed_sc > 0) {
    log <- c(log, sprintf("\nRemoved %d rows from sc due to NA in p_formula variables.\n", removed_sc))
  }

  # sc for fitting
  sc_clean <- if (na_mode %in% c("omit", "exclude", "fail")) {
    sc[keep_sc, , drop = FALSE]
  } else {
    sc
  }

  #----------------------------------------------------------
  # 2. sp: one-reference
  #----------------------------------------------------------
  if (is.data.frame(sp)) {

    # one-reference: keep your current behavior: use vars_sc for sp too
    if (length(vars_sc) > 0) {
      .na_summary(sp, vars_sc, "sp")
      keep_sp <- complete.cases(sp[, vars_sc, drop = FALSE])
    } else {
      keep_sp <- rep(TRUE, nrow(sp))
    }
    removed_sp <- sum(!keep_sp)

    if (na_mode == "fail" && removed_sp > 0) {
      stop("na.fail: NA found in p_formula variables in sp.", call. = FALSE)
    }

    if (na_mode == "pass") {
      keep_sp[] <- TRUE
    } else if (na_mode %in% c("omit", "exclude") && removed_sp > 0) {
      log <- c(log, sprintf("\nRemoved %d rows from sp due to NA in p_formula variables.\n", removed_sp))
    }

    sp_clean <- if (na_mode %in% c("omit", "exclude", "fail")) {
      sp[keep_sp, , drop = FALSE]
    } else {
      sp
    }
  # sp: multi-reference
  } else if (is.list(sp) && all(vapply(sp, is.data.frame, logical(1)))) {

    sp_clean <- vector("list", length(sp))

    sp_names <- names(sp)
    use_names <- if (!is.null(sp_names) && all(nzchar(sp_names))) {
      sp_names
    } else {
      paste0("sp[[", seq_along(sp), "]]")
    }

    for (j in seq_along(sp)) {

      if (!is.list(p_formula) || length(p_formula) < j) {
        stop("For multi-reference sp, p_formula must be a list with same length as sp.", call. = FALSE)
      }

      vars_j <- all.vars(p_formula[[j]])

      if (length(vars_j) > 0) {
        .na_summary(sp[[j]], vars_j, use_names[j])
        keep_spj <- complete.cases(sp[[j]][, vars_j, drop = FALSE])
      } else {
        keep_spj <- rep(TRUE, nrow(sp[[j]]))
      }
      removed_spj <- sum(!keep_spj)

      if (na_mode == "fail" && removed_spj > 0) {
        stop(sprintf("na.fail: NA found in p_formula variables in %s.", use_names[j]), call. = FALSE)
      }

      if (na_mode == "pass") {
        keep_spj[] <- TRUE
      } else if (na_mode %in% c("omit", "exclude") && removed_spj > 0) {
        log <- c(log, sprintf("\nRemoved %d rows from %s due to NA in p_formula variables.\n",
                              removed_spj, use_names[j]))
      }

      sp_clean[[j]] <- if (na_mode %in% c("omit", "exclude", "fail")) {
        sp[[j]][keep_spj, , drop = FALSE]
      } else {
        sp[[j]]
      }
    }

    names(sp_clean) <- use_names

  } else {
    stop("'sp' must be a data.frame (one reference) or a list of data.frames (multi reference).",
         call. = FALSE)
  }

  #----------------------------------------------------------
  # 3. lm-style na.action object for sc mapping
  #----------------------------------------------------------
  na_action_obj <- NULL
  if (na_mode %in% c("omit", "exclude")) {
    idx_drop <- which(!keep_sc)
    cls <- if (na_mode == "exclude") "exclude" else "omit"
    na_action_obj <- structure(idx_drop, class = cls)
  }

  return(list(
    sc        = sc_clean,
    sp        = sp_clean,
    keep_sc   = keep_sc,        # for re-attaching weights to sc0
    na_action = na_action_obj,  # for attr(sc_updated, "na.action")
    log       = log
  ))
}
