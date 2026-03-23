check_lonely_psu <- function(des, name = "sp_des") {

  if (!inherits(des, "survey.design2")) {
    return(invisible(NULL))
  }

  # both strata and cluster must exist
  if (is.null(des$strata) || is.null(des$cluster)) {
    return(invisible(NULL))
  }

  strata <- des$strata[[1]]
  psu    <- des$cluster[[1]]

  ok <- !(is.na(strata) | is.na(psu))
  strata <- strata[ok]
  psu    <- psu[ok]

  if (length(strata) == 0L) {
    return(invisible(NULL))
  }

  # number of unique PSUs in each stratum
  n_psu <- tapply(psu, strata, function(x) length(unique(x)))
  lonely_strata <- names(n_psu)[n_psu == 1L]

  if (length(lonely_strata) == 0L) {
    return(invisible(NULL))
  }

  max_show <- 5L
  show_strata <- head(lonely_strata, max_show)

  more_flag <- if (length(lonely_strata) > max_show) {
    paste0(" ... (", length(lonely_strata) - max_show, " more)")
  } else {
    ""
  }

  lonely_opt <- getOption("survey.lonely.psu")
  if (is.null(lonely_opt)) lonely_opt <- "fail"
  lonely_opt <- tolower(trimws(as.character(lonely_opt)))

  msg <- paste0(
    "Lonely PSU detected in ", name, ".\n",
    "Strata with only one PSU: ",
    paste(show_strata, collapse = ", "),
    more_flag, ".\n",
    "Current setting survey.lonely.psu = \"", lonely_opt, "\".\n",
    "Please either:\n",
    "  (1) combine these strata with nearby strata, or\n",
    "  (2) set options(survey.lonely.psu = \"adjust\"), ",
    "\"average\", \"certainty\", or \"remove\"."
  )

  if (identical(lonely_opt, "fail")) {
    stop(msg, call. = FALSE)
  }

  invisible(NULL)
}
