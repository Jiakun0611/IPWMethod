make_block_D <- function(sp_list, wts_cols, xcol) {
  n <- length(sp_list)

  Xp_list <- vector("list", n)
  block_sizes <- integer(n)

  # 1) build Xp first, then determine dimensions from Xp itself
  for (i in seq_len(n)) {
    sp <- sp_list[[i]]

    Xp_list[[i]] <- if (i == 1) {
      design_matrix(xcol[[i]], sp, intercept = TRUE)
    } else {
      design_matrix(xcol[[i]], sp, intercept = FALSE)
    }

    block_sizes[i] <- ncol(Xp_list[[i]])
  }

  total_dim <- sum(block_sizes)

  # 2) compute block boundaries
  starts <- cumsum(c(1, head(block_sizes, -1)))
  ends   <- cumsum(block_sizes)

  # 3) initialize D
  D <- matrix(0, total_dim, total_dim)

  # 4) fill each diagonal block
  for (i in seq_len(n)) {
    sp  <- sp_list[[i]]
    wts <- sp[[wts_cols[i]]]
    Xp  <- Xp_list[[i]]

    c_i <- wts^2 * (1 - 1 / wts)
    D_i <- t(c_i * Xp) %*% Xp

    idx <- starts[i]:ends[i]
    D[idx, idx] <- D_i
  }

  D
}
