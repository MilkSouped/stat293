# make the readings and keep the true means, loadings, factors and noise
generate_data <- function() {
  n <- 100
  p <- 500
  k <- 3
  mu <- numeric(p)
  mu[1:125] <- 0.6

  B <- matrix(runif(p * k, -1, 1), nrow = p)
  factors <- matrix(rnorm(n * k), nrow = n)
  noise <- matrix(rt(n * p, df = 3), nrow = n)

  X <- factors %*% t(B) + noise
  for (j in seq_len(p)) {
    X[, j] <- X[, j] + mu[j]
  }

  list(X = X, mu = mu, B = B, factors = factors, noise = noise)
}
