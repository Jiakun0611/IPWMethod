#' Compute design-based covariance matrix D from a survey design object
#'
#' Internal helper to compute
#'   D = Var_p( sum_{i in sp} d_i * h_i )
#' when the probability sample is already supplied as a survey design object
#' created by survey::svydesign() or survey::svrepdesign().
#'
#' @param sp_des A survey design object of class "survey.design2" or
#'   "svyrep.design".
#' @param Xp Numeric matrix of dimension n_p x p. Each row is h_i.
#'
#' @return A p x p covariance matrix D.
#' @keywords internal
compute_D_raking <- function(sp_des, Xp, lonely.psu = "adjust") {

  if (!inherits(sp_des, c("survey.design2", "svyrep.design"))) {
    stop("'sp_des' must be a survey design object from svydesign() or svrepdesign().",
         call. = FALSE)
  }

  if (!is.matrix(Xp)) {
    Xp <- as.matrix(Xp)
  }

  n_sp <- nrow(sp_des$variables)
  if (nrow(Xp) != n_sp) {
    stop(sprintf(
      "Row mismatch: nrow(Xp) = %d but the survey design object contains %d observations.",
      nrow(Xp), n_sp
    ), call. = FALSE)
  }

  # safe column names for formula parsing
  xp_names <- make.names(colnames(Xp), unique = TRUE)
  if (is.null(xp_names)) {
    xp_names <- paste0("x", seq_len(ncol(Xp)))
  }

  H_df <- as.data.frame(Xp)
  colnames(H_df) <- paste0("h_", xp_names)

  # attach estimating-equation columns to the existing design object
  sp_des_h <- sp_des
  sp_des_h$variables <- cbind(sp_des_h$variables, H_df)

  # no intercept
  fml <- stats::reformulate(colnames(H_df), intercept = FALSE)

  tot <- survey::svytotal(fml, sp_des_h)
  D <- stats::vcov(tot)

  return(D)
}
