#' Cumulative pre-calibration across multiple probability samples
#'
#' Performs a cumulative (sequential) pre-calibration across multiple probability
#' samples by aligning the marginal totals of shared expanded calibration columns.
#' Samples can optionally be ordered by decreasing sample size, and the target
#' total vector is updated (grown) as new calibration columns appear.
#'
#' @param sp_raw A list of data frames (raw reference samples). Each element must
#'   contain the weight column specified in \code{weight[[i]]}. All other columns
#'   are treated as calibration variables and are expanded using
#'   \code{\link[stats]{model.matrix}}.
#' @param sp_new A list of data frames (working reference samples) with the same
#'   structure as \code{sp_raw}. The calibrated weights are written back to
#'   \code{sp_new[[i]][[weight[[i]]]]}.
#' @param weight A list (or character vector) of weight column names, one for each
#'   sample in \code{sp_raw}/\code{sp_new}.
#' @param sp_order Character string controlling the ordering of samples. If
#'   \code{"size"}, samples are processed from largest to smallest; otherwise the
#'   original list order is used.
#' @param verbose Logical. If \code{TRUE}, prints a short log describing which
#'   sample is used to initialize the target totals and which variables are used
#'   to calibrate each subsequent sample.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{sp_new}: the updated \code{sp_new} list with calibrated weights.
#'   \item \code{total_vector}: the final cumulative population total vector used
#'         in the sequential procedure (intercept + expanded columns).
#'   \item \code{log_messages}: character vector of printed messages (useful for
#'         debugging or reporting).
#'   \item \code{order_used}: integer vector giving the processing order (indices
#'         of the original input lists).
#' }
#'
#' @details
#' The function proceeds as follows:
#' \enumerate{
#'   \item Choose the initial sample (largest if \code{sp_order = "size"}, else the
#'         first). Expand its calibration variables using \code{model.matrix} and
#'         compute a target total vector consisting of (i) the total survey weight
#'         (intercept) and (ii) survey-weighted totals of all expanded columns.
#'   \item For each remaining sample, expand its calibration variables and identify
#'         the set of columns shared with the current target total vector. Calibrate
#'         the sample weights using \code{\link[survey]{calibrate}} so that the
#'         intercept (total weight) and the shared column totals match the current
#'         targets.
#'   \item After calibration, append totals of any newly encountered expanded
#'         columns to the target total vector (computed under the calibrated
#'         weights), so the target vector grows monotonically across steps.
#' }
#' Calibration is carried out using the \pkg{survey} package via
#' \code{\link[survey]{svydesign}}, \code{\link[survey]{svytotal}}, and
#' \code{\link[survey]{calibrate}} with \code{ids = ~1}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' sp1_raw <- data.frame(a = factor(sample(letters[1:3], 100, TRUE)),
#'                       x = rnorm(100), wt1 = runif(100, 0.5, 1.5))
#' sp2_raw <- data.frame(a = factor(sample(letters[1:4], 120, TRUE)),
#'                       z = rnorm(120), wt2 = runif(120, 0.5, 1.5))
#'
#' # working copies (weights will be overwritten)
#' sp1_new <- sp1_raw
#' sp2_new <- sp2_raw
#'
#' out <- precal_cumulative_order(
#'   sp_raw   = list(sp1_raw, sp2_raw),
#'   sp_new   = list(sp1_new, sp2_new),
#'   weight   = list("wt1", "wt2"),
#'   sp_order = "size",
#'   verbose  = TRUE
#' )
#'
#' str(out$total_vector)
#' }
#'
#' @seealso [survey::svydesign], [survey::svytotal], [survey::calibrate],
#'   [stats::model.matrix]
#' @importFrom survey svydesign svytotal calibrate
#' @importFrom stats model.matrix as.formula weights coef
#' @export
#'

precal_cumulative_order <- function(sp_raw, sp_new, weight, sp_order, verbose = FALSE) {

  order <- (sp_order == "size")

  log_messages <- character()

  # ---- names ----
  sp_raw_names <- names(sp_raw)
  if (is.null(sp_raw_names)) sp_raw_names <- rep("", length(sp_raw))

  # ---- order control ----
  sizes <- vapply(sp_raw, nrow, integer(1))
  ord <- if (order) order(sizes, decreasing = TRUE) else seq_along(sp_raw)

  sp_raw_ord <- sp_raw[ord]
  sp_new_ord <- sp_new[ord]
  w_ord      <- weight[ord]
  nm_ord     <- sp_raw_names[ord]
  sz_ord     <- sizes[ord]  # (kept in case you want to print sizes later)

  # ---- helpers ----
  .get_wname <- function(w) {
    if (length(w) != 1 || !is.character(w) || !nzchar(w)) {
      stop("Each weight[[i]] must be a single non-empty character string.")
    }
    w
  }

  .expand_raw <- function(df, wname) {
    vars <- setdiff(names(df), wname)
    if (length(vars) == 0) {
      return(matrix(, nrow = nrow(df), ncol = 0,
                    dimnames = list(NULL, character(0))))
    }
    X <- stats::model.matrix(~ ., data = df[, vars, drop = FALSE])
    if ("(Intercept)" %in% colnames(X)) X <- X[, -1, drop = FALSE]
    X
  }

  .svy_design_from_X <- function(X, w) {
    dat <- as.data.frame(X)
    dat$.w <- w
    survey::svydesign(ids = ~1, weights = ~.w, data = dat)
  }

  .population_vector <- function(ds, cols) {
    total_w <- sum(stats::weights(ds))
    if (length(cols) == 0) return(c("(Intercept)" = total_w))
    fml  <- stats::as.formula(paste0("~", paste(cols, collapse = " + ")))
    marg <- stats::coef(survey::svytotal(fml, ds))
    c("(Intercept)" = total_w, marg)
  }

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

  # ---- Step 0: init target totals from reference RAW sample ----
  ref_raw   <- sp_raw_ord[[1]]
  ref_new   <- sp_new_ord[[1]]
  ref_wname <- .get_wname(w_ord[[1]])

  X_ref  <- .expand_raw(ref_raw, ref_wname)
  ds_ref <- .svy_design_from_X(X_ref, ref_new[[ref_wname]])
  total_vec <- .population_vector(ds_ref, colnames(X_ref))

  ref_label <- if (nzchar(nm_ord[1])) {
    nm_ord[1]
  } else {
    paste0("sp_raw[[", ord[1], "]]")
  }

  s <- if (order) "largest" else "first"

  msg <- sprintf(
    "\nPre-calibration summary:\nNon-calibrated sample (%s): %s\n",
    s, ref_label
  )

  log_messages <- c(log_messages, msg)
  if (verbose) cat(msg)

  # ---- Step 1..K: calibrate remaining samples ----
  if (length(sp_raw_ord) >= 2) {
    for (k in 2:length(sp_raw_ord)) {

      df_raw <- sp_raw_ord[[k]]
      df_new <- sp_new_ord[[k]]
      wname  <- .get_wname(w_ord[[k]])

      X_k <- .expand_raw(df_raw, wname)

      shared_cols <- intersect(colnames(X_k), setdiff(names(total_vec), "(Intercept)"))
      shared_cols <- sort(shared_cols)

      ds_k <- .svy_design_from_X(X_k, df_new[[wname]])

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

      sp_new_ord[[k]][[wname]] <- as.numeric(stats::weights(ds_k_cal))

      # append totals for NEW columns (using calibrated weights)
      current_cols <- setdiff(names(total_vec), "(Intercept)")
      new_cols <- setdiff(colnames(X_k), current_cols)
      new_cols <- sort(new_cols)

      if (length(new_cols) > 0) {
        ds_k2 <- .svy_design_from_X(X_k, sp_new_ord[[k]][[wname]])
        add_vec <- .population_vector(ds_k2, new_cols)
        add_vec <- add_vec[names(add_vec) != "(Intercept)"]
        total_vec <- c(total_vec, add_vec)
      }

      k_label <- if (nzchar(nm_ord[k])) nm_ord[k] else paste0("sp_raw[[", ord[k], "]]")
      msg <- sprintf(
        "Calibrated sample: %s\n  Calibration variables: %s\n",
        k_label, .pretty_cols(shared_cols)
      )
      log_messages <- c(log_messages, msg)
      if (verbose) cat(msg)
    }
  }

  # ---- restore original order only if we reordered ----
  sp_back <- if (order) {
    out <- vector("list", length(sp_new))
    out[ord] <- sp_new_ord
    names(out) <- names(sp_new)
    out
  } else {
    sp_new_ord
  }

  list(
    sp_new       = sp_back,
    total_vector = total_vec,
    log_messages = log_messages,
    order_used   = ord
  )
}
