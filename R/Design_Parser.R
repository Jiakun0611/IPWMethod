parse_ipwm_data <- function(data) {
  if (!is.list(data) || length(data) < 2L) {
    stop("`data` must be a list like list(sc, sp1.des, sp2.des, ...).",
         call. = FALSE)
  }

  sc <- data[[1]]
  sp_des <- data[-1]

  if (!is.data.frame(sc)) {
    stop("The first element of `data` must be a data.frame for `sc`.",
         call. = FALSE)
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

  # preserve names of reference surveys if provided
  ref_names <- names(sp_des)
  if (is.null(ref_names)) {
    ref_names <- paste0("sp", seq_along(sp_des))
  }

  sp_vars <- lapply(seq_along(sp_des), function(i) {
    tryCatch(
      extract_analysis_data(sp_des[[i]]),
      error = function(e) {
        stop(
          sprintf("Failed to parse reference survey %d: %s", i, e$message),
          call. = FALSE
        )
      }
    )
  })

  # keep names on both sp_des and sp_vars
  names(sp_des)  <- ref_names
  names(sp_vars) <- ref_names

  list(
    sc = sc,
    sp_des = sp_des,
    sp_vars = sp_vars,
    n_ref = length(sp_des)
  )
}
