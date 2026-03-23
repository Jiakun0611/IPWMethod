sort_by_sp_size <- function(sp, vars, weight, design,
                            sp_order, verbose = FALSE) {

  log_messages <- character()
  order_by_size <- identical(sp_order, "size")

  n <- length(sp)

  # vars: if not provided, initialize as empty list
  if (length(vars) == 0) vars <- vector("list", length = n)

  # strict checks
  if (!is.list(sp) || n == 0) {
    stop("'sp' must be a non-empty list.", call. = FALSE)
  }

  if (length(weight) != n) {
    stop("'weight' must have the same length as 'sp'.", call. = FALSE)
  }

  if (!is.list(vars) || length(vars) != n) {
    stop("'vars' must be a list with the same length as 'sp'.", call. = FALSE)
  }

  # design is REQUIRED and must be a list of survey design objects
  if (missing(design) || is.null(design)) {
    stop("'design' is required for the multi-reference case.", call. = FALSE)
  }

  if (!is.list(design) || length(design) != n) {
    stop("'design' must be a list of survey design objects with length equal to length(sp).",
         call. = FALSE)
  }

  bad_i <- which(
    !vapply(design,
            function(d) inherits(d, c("survey.design2", "svyrep.design")),
            logical(1))
  )

  if (length(bad_i) > 0) {
    stop(
      sprintf(
        "Each design[[i]] must be a survey design object of class 'survey.design2' or 'svyrep.design'. Problem at i = %s.",
        paste(bad_i, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # optional consistency check: each sp[[i]] should align with design[[i]]
  bad_n <- which(
    !vapply(seq_len(n), function(i) {
      NROW(sp[[i]]) == nrow(design[[i]]$variables)
    }, logical(1))
  )

  if (length(bad_n) > 0) {
    stop(
      sprintf(
        "Row mismatch between sp[[i]] and design[[i]]$variables at i = %s.",
        paste(bad_n, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # compute row counts
  sizes <- vapply(sp, NROW, integer(1))

  # ordering index
  idx <- if (order_by_size) order(sizes, decreasing = TRUE) else seq_along(sp)

  # reorder
  sp_sorted     <- sp[idx]
  vars_sorted   <- vars[idx]
  weight_sorted <- weight[idx]
  design_sorted <- design[idx]

  # names
  sp_names <- names(sp)
  if (is.null(sp_names)) sp_names <- paste0("sp", seq_along(sp))
  sorted_names <- sp_names[idx]
  sorted_sizes <- sizes[idx]

  # message
  msg_lines <- paste0("  ", sorted_names, " (n = ", sorted_sizes, ")")
  msg <- if (order_by_size) {
    paste0(
      "\nReference samples summary:\n",
      "Order of samples by size (largest to smallest):\n",
      paste(msg_lines, collapse = "\n"), "\n"
    )
  } else {
    paste0(
      "\nReference samples summary:\n",
      "Order of samples kept as provided:\n",
      paste(msg_lines, collapse = "\n"), "\n"
    )
  }

  log_messages <- c(log_messages, msg)
  if (verbose) cat(msg)

  # assign names
  names(sp_sorted)     <- sorted_names
  names(vars_sorted)   <- sorted_names
  names(weight_sorted) <- sorted_names
  names(design_sorted) <- sorted_names

  list(
    sp         = sp_sorted,
    vars       = vars_sorted,
    weight     = weight_sorted,
    design     = design_sorted,
    order_used = idx,
    log        = log_messages
  )
}
