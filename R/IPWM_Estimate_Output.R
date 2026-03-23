assemble_output <- function(build, est, naive, na_info) {



  mean_unw <- naive$mean
  se_unw   <- sqrt(naive$var)

  mean_adj <- est$mean
  se_adj   <- sqrt(est$variance)

  out <- list(
    call = match.call(),
    method = build$method,

    unweighted = list(
      mean = mean_unw,
      se   = se_unw,
      CI_95 = c(
        lower = mean_unw - 1.96 * se_unw,
        upper = mean_unw + 1.96 * se_unw
      )
    ),

    adjusted = list(
      mean = mean_adj,
      se   = se_adj,
      CI_95 = c(
        lower = mean_adj - 1.96 * se_adj,
        upper = mean_adj + 1.96 * se_adj
      )
    ),

    na = na_info
  )

  class(out) <- "IPWM_estimate"
  out
}
