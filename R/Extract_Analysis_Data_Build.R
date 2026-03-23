extract_analysis_data <- function(des, weight_name = "weights") {

  if (!inherits(des, c("survey.design2", "survey.design", "svyrep.design"))) {
    stop("`des` must be a survey design object.", call. = FALSE)
  }

  vars <- des$variables

  if (is.null(vars) || !is.data.frame(vars)) {
    stop("`des$variables` is missing or not a data.frame.", call. = FALSE)
  }

  get_call_vars <- function(x) {
    if (is.null(x)) return(character(0))
    out <- tryCatch(all.vars(x), error = function(e) character(0))
    setdiff(out, "1")
  }

  drop_vars <- unique(c(
    get_call_vars(des$call$ids),
    get_call_vars(des$call$strata),
    get_call_vars(des$call$fpc),
    get_call_vars(des$call$probs),
    get_call_vars(des$call$weights)
  ))

  if (inherits(des, "svyrep.design")) {
    rep_cols <- tryCatch(colnames(des$repweights), error = function(e) NULL)
    if (!is.null(rep_cols)) {
      drop_vars <- unique(c(drop_vars, rep_cols))
    }
  }

  drop_vars <- intersect(drop_vars, names(vars))
  out <- vars[, !(names(vars) %in% drop_vars), drop = FALSE]

  if (inherits(des, "svyrep.design")) {
    out[[weight_name]] <- as.numeric(weights(des, type = "sampling"))
  } else {
    out[[weight_name]] <- as.numeric(weights(des))
  }

  out
}
