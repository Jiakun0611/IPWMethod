#' @export
#'
summary.IPWM_build <- function(object, ...) {

  cat("Call:\n")
  print(object$call)

  method <- tolower(trimws(object$method))

  if (method %in% c("alp", "clw", "cali", "calibration")) {
    if (method == "alp") method <- "ALP"
    if (method == "clw") method <- "CLW"
    if (method == "cali") method <- "calibration"

    cat("\nMethod: One reference", method, "\n")
  } else if (method == "multi") {
    cat("\nMethod: Multi reference calibration\n")
  }

  # --- Optional log output (for multi-reference runs) ---
  if (!is.null(object$log_messages) && length(object$log_messages) > 0) {
    cat("\n", paste(object$log_messages, collapse = "") , sep = "")
  }

  # --- Model information ---
  if (!is.null(object$variables)) {
    cat("\nParticipation model involves the following variables:\n")
    cat(object$variables[-1], "\n\n")
  }

  # --- Coefficients section ---
  if (!is.null(object$coefficients)) {
    cat("Selection model coefficients:\n")

    df <- rbind(object$coefficients)
    colnames(df) <- object$variables
    rownames(df) <- NULL

    formatted_df <- format(round(df, 4), justify = "left", width = 10)
    print(as.data.frame(formatted_df), row.names = FALSE)
  }

  invisible(object)
}
