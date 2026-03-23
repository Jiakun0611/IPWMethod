#' Cumulative pre-calibration using raw reference samples
#'
#' Performs a cumulative (sequential) pre-calibration across reference samples,
#' where the \emph{calibration space} is defined from the raw samples
#' (\code{sp_raw}) using \strong{all available variables} (excluding weights),
#' expanded to dummy columns via \code{\link[stats]{model.matrix}} with formula
#' \code{~ .}. Samples are ordered by decreasing raw sample size. The largest raw
#' sample defines an initial \emph{total vector} consisting of survey-weighted
#' totals of all expanded columns, together with an intercept term equal to the
#' total survey weight. Each subsequent sample is calibrated to match the
#' current total vector on the set of shared expanded columns. After
#' calibration, totals of previously unused expanded columns from the current
#' sample are appended to the total vector, so target totals grow monotonically.
#'
#' This function separates \emph{variable definition} from \emph{weight storage}:
#' \code{sp_raw} is used only to determine which expanded columns are available
#' and shared, while calibrated weights are written back to \code{sp_new} (the
#' processed samples used downstream).
#'
#' Importantly, the intercept (total survey weight) used in calibration is
#' \emph{fixed} to the total weight of the largest sample throughout all steps.
#'
#' @param sp_raw A list of raw data frames (samples). These data frames provide
#'   the full set of variables used to define the calibration space. For each
#'   sample, all columns except the corresponding weight column are expanded via
#'   \code{\link[stats]{model.matrix}} with formula \code{~ .}. The expanded
#'   intercept column (if present) is dropped.
#' @param sp_new A list of processed data frames (samples), of the same length
#'   and order as \code{sp_raw}. Calibrated weights are stored by updating the
#'   weight column in \code{sp_new}. This allows pre-calibration to use all raw
#'   variables while keeping downstream estimation consistent with the processed
#'   samples used by IPWM.
#' @param weight A list (or character vector) of length \code{length(sp_raw)}
#'   giving the name of the weight column in each sample. Each entry must be a
#'   single non-empty character string and must exist in both \code{sp_raw[[i]]}
#'   and \code{sp_new[[i]]}.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#'
#' @return A list with components:
#' \describe{
#'   \item{sp_new}{A list of processed data frames in the original input order,
#'   with updated weight columns for all calibrated samples (all except the
#'   largest one).}
#'   \item{total_vector}{A named numeric vector containing the final target
#'   totals (including \code{"(Intercept)"}).}
#'   \item{log_messages}{A character vector of log messages.}
#'   \item{order_used}{An integer vector giving the permutation used to order
#'   samples by decreasing raw sample size.}
#' }
#'
#' @details
#' Let samples be ordered as \eqn{\mathcal{S}_{(1)},\ldots,\mathcal{S}_{(K)}} by
#' decreasing \code{nrow(sp_raw[[i]])}. For each raw sample, define the expanded
#' calibration matrix \eqn{X_{(k)}} from all non-weight variables using
#' \code{model.matrix(~ ., ...)} and dropping the intercept column. The algorithm:
#' \enumerate{
#'   \item Initializes the target total vector using \eqn{\mathcal{S}_{(1)}}:
#'   \eqn{T^{(1)}_{\text{Intercept}}=\sum_{i\in \mathcal{S}_{(1)}} w_i} and
#'   \eqn{T^{(1)}_j=\sum_{i\in \mathcal{S}_{(1)}} w_i x_{ij}} for all expanded
#'   columns \eqn{j}. Here weights \eqn{w_i} are taken from \code{sp_new[[ref]]}.
#'   \item For each \eqn{k\ge2}, identifies the shared expanded columns between
#'   \eqn{X_{(k)}} and the current target vector, and calibrates the weights of
#'   \eqn{\mathcal{S}_{(k)}} using \code{\link[survey]{calibrate}} so that weighted
#'   totals of \eqn{(1, X_{(k)}^{\text{shared}})} match the target totals, while
#'   keeping the intercept fixed to \eqn{T^{(1)}_{\text{Intercept}}}.
#'   \item Appends totals of new expanded columns (not yet in the target vector)
#'   computed under the calibrated weights of \eqn{\mathcal{S}_{(k)}}.
#' }
#'
#' The function uses \code{\link[survey]{svydesign}}, \code{\link[survey]{svytotal}},
#' and \code{\link[survey]{calibrate}} from the \pkg{survey} package.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' sp1_raw <- data.frame(x = rnorm(200),
#'                       g = factor(sample(letters[1:3], 200, TRUE)),
#'                       w1 = runif(200, 0.5, 1.5),
#'                       extra = rnorm(200))
#' sp2_raw <- data.frame(x = rnorm(120),
#'                       g = factor(sample(letters[1:2], 120, TRUE)),
#'                       z = rnorm(120),
#'                       w2 = runif(120, 0.5, 1.5))
#'
#' # processed versions used downstream (e.g., after applying p_formula)
#' sp1_new <- sp1_raw[, c("x", "g", "w1")]
#' sp2_new <- sp2_raw[, c("x", "z", "w2")]
#'
#' out <- precal_cumulative(
#'   sp_raw  = list(sp1_raw, sp2_raw),
#'   sp_new  = list(sp1_new, sp2_new),
#'   weight  = list("w1", "w2"),
#'   verbose = TRUE
#' )
#'
#' sp_new2 <- out$sp_new
#' tv      <- out$total_vector
#' }
#'
#' @seealso [survey::svydesign], [survey::svytotal], [survey::calibrate]
#' @importFrom survey svydesign calibrate svytotal
#' @importFrom stats weights as.formula coef
#' @export
#'

precal_cumulative <- function(sp_raw, sp_new, weight, verbose = FALSE) {

  log_messages <- character()

  # ---- names ----
  sp_raw_names  <- names(sp_raw)


  # ---- order by raw sample size (desc) ----
  sizes <- vapply(sp_raw, nrow, integer(1))
  ord <- order(sizes, decreasing = TRUE)

  sp_raw_ord <- sp_raw[ord]
  sp_new_ord <- sp_new[ord]
  w_ord      <- weight[ord]
  nm_ord     <- sp_raw_names[ord]
  sz_ord     <- sizes[ord]

  # ---- helpers ---- (no need)
  .get_wname <- function(w) {
    if (length(w) != 1 || !is.character(w) || !nzchar(w)) {
      stop("Each weight[[i]] must be a single non-empty character string.")
    }
    w
  }

  # expand using ALL raw variables except weight; drop intercept
  .expand_raw <- function(df, wname) {
    vars <- setdiff(names(df), wname)
    if (length(vars) == 0) {
      return(matrix(, nrow = nrow(df), ncol = 0, dimnames = list(NULL, character(0))))
    }
    X <- stats::model.matrix(~ ., data = df[, vars, drop = FALSE])
    if ("(Intercept)" %in% colnames(X)) X <- X[, -1, drop = FALSE]
    X
  }

  # survey design from expanded X and weights
  .svy_design_from_X <- function(X, w) {
    dat <- as.data.frame(X)
    dat$.w <- w
    survey::svydesign(ids = ~1, weights = ~.w, data = dat)
  }

  # population vector: intercept + totals of given columns
  .population_vector <- function(ds, cols) {
    total_w <- sum(stats::weights(ds))
    if (length(cols) == 0) return(c("(Intercept)" = total_w))
    fml  <- stats::as.formula(paste0("~", paste(cols, collapse = " + ")))
    marg <- stats::coef(survey::svytotal(fml, ds))
    c("(Intercept)" = total_w, marg)
  }

  # pretty print expanded columns (dummy form)
  .pretty_cols <- function(cols, max_show = 30) {
    cols <- sort(unique(cols))
    if (length(cols) == 0) return("survey weights total only")
    if (length(cols) <= max_show) {
      return(paste0("survey weights total, ", paste(cols, collapse = ", ")))
    }
    paste0("survey weights total, ",
           paste(cols[1:max_show], collapse = ", "),
           ", ... (", length(cols), " cols)")
  }

  # ---- Step 0: init target totals from largest RAW sample (sp1) ----
  ref_raw   <- sp_raw_ord[[1]]
  ref_new   <- sp_new_ord[[1]]
  ref_wname <- .get_wname(w_ord[[1]])

  X_ref <- .expand_raw(ref_raw, ref_wname)

  # (IMPORTANT) weights come from sp_new (for survey design object)
  ds_ref <- .svy_design_from_X(X_ref, ref_new[[ref_wname]])

  total_vec <- .population_vector(ds_ref, colnames(X_ref))

  msg <- sprintf(
    "\nPre-calibration summary:\nNon-calibrated sample (largest): %s\n",
    nm_ord[1]
  )
  log_messages <- c(log_messages, msg)
  if (verbose) cat(msg)

  # ---- Step 1..K: calibrate each remaining sample using RAW expanded columns ----
  for (k in 2:length(sp_raw_ord)) {

    df_raw <- sp_raw_ord[[k]]
    df_new <- sp_new_ord[[k]]
    wname  <- .get_wname(w_ord[[k]])

    X_k <- .expand_raw(df_raw, wname)

    # shared expanded columns between this RAW sample and current total vector
    shared_cols <- intersect(colnames(X_k), setdiff(names(total_vec), "(Intercept)"))
    shared_cols <- sort(shared_cols)

    # build design with current weights (from sp_new)
    ds_k <- .svy_design_from_X(X_k, df_new[[wname]])

    # fixed intercept = sp1 total weight; avoid name concatenation
    pops <- c("(Intercept)" = unname(total_vec["(Intercept)"]))
    if (length(shared_cols) > 0) pops <- c(pops, total_vec[shared_cols])

    fml <- if (length(shared_cols) > 0) {
      stats::as.formula(paste0("~", paste(shared_cols, collapse = " + ")))
    } else {
      ~1
    }

    ds_k_cal <- survey::calibrate(
      design     = ds_k,
      formula    = fml,
      population = pops
    )

    # update weights in sp_new (ordered)
    new_w <- stats::weights(ds_k_cal)
    sp_new_ord[[k]][[wname]] <- as.numeric(new_w)

    # append totals for NEW columns (using calibrated weights)
    current_cols <- setdiff(names(total_vec), "(Intercept)")
    new_cols <- setdiff(colnames(X_k), current_cols)
    new_cols <- sort(new_cols)

    if (length(new_cols) > 0) {
      ds_k2 <- .svy_design_from_X(X_k, sp_new_ord[[k]][[wname]])
      add_vec <- .population_vector(ds_k2, new_cols)
      add_vec <- add_vec[names(add_vec) != "(Intercept)"]  # keep intercept fixed (from ref sample)
      total_vec <- c(total_vec, add_vec)
    }

    msg <- sprintf(
      "Calibrated sample: %s\n  Calibration variables: %s\n",
      nm_ord[k], .pretty_cols(shared_cols)
    )
    log_messages <- c(log_messages, msg)
    if (verbose) cat(msg)
  }

  # ---- restore original order ----
  sp_back <- vector("list", length(sp_new))
  sp_back[ord] <- sp_new_ord
  names(sp_back) <- names(sp_new)

  list(
    sp_new       = sp_back,
    total_vector = total_vec,
    log_messages = log_messages,
    order_used   = ord
  )
}



