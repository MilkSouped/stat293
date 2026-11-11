# test each column against zero and keep the means and bh decisions for scoring
fit_t_tests <- function(X, alpha) {
  p_values <- numeric(ncol(X))
  for (j in seq_len(ncol(X))) {
    p_values[j] <- t.test(X[, j], mu = 0)$p.value
  }

  adjusted <- p.adjust(p_values, method = "BH")
  list(means = colMeans(X), significant = adjusted <= alpha)
}
