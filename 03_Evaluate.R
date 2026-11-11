# compare one method's calls and mean estimates with the planted means
score_fit <- function(method, fit, mu) {
  signal <- mu != 0
  called <- as.logical(fit$significant)
  true_calls <- sum(called & signal)
  false_calls <- sum(called & !signal)

  data.frame(
    Method = method,
    Discoveries = sum(called),
    False_discoveries = false_calls,
    True_discoveries = true_calls,
    FDP = false_calls / max(1, sum(called)),
    Power = if (any(signal)) true_calls / sum(signal) else NA_real_,
    MSE = mean((fit$means - mu)^2)
  )
}

# find the percent of variation in one source captured by the given pc directions
variation_captured <- function(source_matrix, directions) {
  centered <- scale(source_matrix, center = TRUE, scale = FALSE)
  scores <- centered %*% directions
  reconstructed <- scores %*% t(directions)

  average_error <- sum(centered^2)
  pc_error <- sum((centered - reconstructed)^2)
  100 * (1 - pc_error / average_error)
}

# use the readings and fit to work through the first gene in the report
gene1_example <- function(readings, fit) {
  huber <- FarmTest::huber.mean(readings)
  n <- length(readings)
  gaps <- (readings - huber)^2
  # the package uses this cutoff when it fits the huber mean
  tau <- sqrt(FarmTest:::rootf1(gaps, n, log(n) / n, min(gaps), sum(gaps)))
  lowest <- which.min(readings)

  p <- fit$pValues[1]
  above_target <- sum(fit$pValues > fit$alpha)
  null_genes <- min(length(fit$pValues), above_target / (1 - fit$alpha))
  list(
    average = round(mean(readings), 3),
    lowest_sample = lowest,
    lowest_reading = round(readings[lowest], 2),
    other_readings = round(range(readings[-lowest]), 2),
    cutoff = round(tau, 2),
    outlier_loss = round(tau * abs(readings[lowest]) - tau^2 / 2),
    huber_mean = round(huber, 3),
    candidates = sum(fit$pValues <= p),
    above_target = above_target,
    expected_false = round(null_genes * p, 1)
  )
}

# compare the known shared contributions and noise along the same pc directions
pca_table <- function(dat, rotation, pcs90) {
  shared <- dat$factors %*% t(dat$B)
  three <- rotation[, 1:3]
  more <- rotation[, seq_len(pcs90), drop = FALSE]
  table <- data.frame(
    source = c("Shared factors", "Gene-specific noise"),
    three = c(variation_captured(shared, three), variation_captured(dat$noise, three)),
    more = c(variation_captured(shared, more), variation_captured(dat$noise, more))
  )
  names(table) <- c("Source of variation", "Captured by 3 PCs (%)",
                    paste0("Captured by ", pcs90, " PCs (%)"))
  table
}
