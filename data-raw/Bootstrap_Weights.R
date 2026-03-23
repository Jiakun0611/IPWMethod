# check nh (PSUs per stratum)
nh <- tapply(sp1$psu, sp1$strata, function(x) length(unique(x)))
summary(nh)

# Confirm no lonely PSU strata
any(nh <= 1)


make_bootstrap_weights <- function(
    sp,
    weight_col,
    psu_col    = "psu",
    strata_col = "strata",
    R          = 500,
    mh         = NULL,
    seed       = 123,
    check_negative = TRUE
) {
  stopifnot(is.data.frame(sp))
  stopifnot(all(c(weight_col, psu_col, strata_col) %in% names(sp)))

  set.seed(seed)

  w  <- sp[[weight_col]]
  ps <- as.character(sp[[psu_col]])
  st <- as.character(sp[[strata_col]])

  nh_by_stratum <- tapply(ps, st, function(x) length(unique(x)))

  if (is.null(mh)) {
    mh_by_stratum <- nh_by_stratum - 1
  } else if (length(mh) == 1 && is.numeric(mh)) {
    mh_by_stratum <- rep(mh, length(nh_by_stratum))
    names(mh_by_stratum) <- names(nh_by_stratum)
  } else {
    mh_by_stratum <- mh
  }

  bad <- which(nh_by_stratum <= 1)
  if (length(bad) > 0) {
    stop("Some strata have nh <= 1: ", paste(names(nh_by_stratum)[bad], collapse = ", "))
  }

  mh_by_stratum <- pmax(1, pmin(mh_by_stratum, nh_by_stratum))
  mh_by_stratum <- setNames(as.numeric(mh_by_stratum), names(nh_by_stratum))

  n <- nrow(sp)
  repW <- matrix(NA_real_, nrow = n, ncol = R)
  strata_levels <- names(nh_by_stratum)

  for (t in seq_len(R)) {
    wt <- w

    for (h in strata_levels) {
      idx_h <- which(st == h)
      psu_h <- unique(ps[idx_h])

      nh   <- length(psu_h)
      mh_h <- as.integer(mh_by_stratum[h])
      if (is.na(mh_h)) stop("mh not found for stratum: ", h)

      # PSU totals in full sample (within this stratum)
      W_hi <- tapply(w[idx_h], ps[idx_h], sum)  # named by PSU

      # draw PSUs and counts r_hi^(t)
      draw <- sample(psu_h, size = mh_h, replace = TRUE)
      r <- tabulate(match(draw, psu_h), nbins = nh)
      names(r) <- psu_h

      # Rao–Wu style PSU multiplier
      g1 = sqrt(mh_h * (nh - 1))
      g2 <- sqrt(mh_h / (nh - 1))
      a <- 1 - g1
      b <- g2 * (nh / mh_h)
      mult_psu <- a + b * r  # length nh, named by psu_h

      # replicate PSU totals
      W_hi_t <- W_hi[names(mult_psu)] * mult_psu

      # push back to SSU level by ratio W_hi_t / W_hi
      ratio <- W_hi_t[ps[idx_h]] / W_hi[ps[idx_h]]
      wt[idx_h] <- w[idx_h] * ratio
    }

    if (check_negative && any(!is.finite(wt) | wt < 0)) {
      stop("Replicate ", t, " produced non-finite or negative weights.")
    }

    repW[, t] <- wt
  }

  colnames(repW) <- paste0("bw", seq_len(R))
  repW
}

repW_500 <- make_bootstrap_weights(
  sp         = sp1,
  weight_col = "wts_sp1",
  psu_col    = "psu",
  strata_col = "strata",
  R          = 500,
  seed       = 2026
)

dim(repW_500)  # nrow(sp1) x 500
head(repW_500[, 1:10])

base_total <- sum(sp1$w)
rep_totals <- colSums(repW_500)
head(rep_totals)

c(base_total = base_total,
  rep_total_mean   = mean(rep_totals),
  rep_sd     = sd(rep_totals))

sp1_bt <- cbind(sp1,repW_500)
head(sp1_bt[,1:20])
