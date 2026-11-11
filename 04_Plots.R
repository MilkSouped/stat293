# draw both loss curves using the huber cutoff as the unit of distance
plot_huber_loss <- function() {
  old_par <- par(mar = c(3, 3.5, 1.5, 0.5), mgp = c(2, 0.7, 0), cex = 1.2)
  on.exit(par(old_par))

  r <- seq(-3, 3, length.out = 201)
  squared <- r^2 / 2
  huber <- ifelse(abs(r) <= 1, squared, abs(r) - 0.5)
  plot(r, squared, type = "l", lwd = 2, xaxt = "n", yaxt = "n",
       xlab = "Gap from trial mean (in cutoff units)", ylab = "Penalty",
       main = "Squared loss and Huber loss", cex.main = 0.95, cex.lab = 0.85)
  axis(1, at = -2:2, cex.axis = 0.85)
  lines(r, huber, col = "blue", lwd = 2)
  abline(v = c(-1, 1), col = "red", lty = 2)
  legend("top", c("Squared loss", "Huber loss"), col = c("black", "blue"),
         lty = 1, lwd = 2, bg = "white", box.lty = 0, cex = 1)
}

# plot the variance per pc and mark how many pcs it takes to reach 90%
plot_pca_variance <- function(percent, pcs90) {
  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1),
                 cex.axis = 0.85, cex.lab = 0.9, cex.main = 1)
  on.exit(par(old_par))

  first_eight <- 1:8
  barplot(percent[first_eight], names.arg = first_eight, border = NA,
          col = ifelse(first_eight <= 3, "blue", "gray"),
          ylim = c(0, 1.15 * max(percent)), cex.names = 0.85,
          xlab = "Principal component", ylab = "Variance (%)", main = "Variance explained")

  cumulative <- cumsum(percent)
  plot(seq_along(percent), cumulative, type = "l", col = "blue", lwd = 2,
       ylim = c(0, 100), xlab = "Number of PCs retained", ylab = "Variance retained (%)",
       main = "Reaching 90%")
  abline(h = 90, v = pcs90, col = "red", lty = 2)
  points(c(3, pcs90), cumulative[c(3, pcs90)], pch = 19, col = c("blue", "red"))
  text(4, 20, paste0("3 PCs: ", round(cumulative[3], 1), "%"), adj = 0, col = "blue", cex = 0.8)
  text(pcs90 - 3, 55, paste(pcs90, "PCs"), adj = 1, col = "red", cex = 0.9)
  text(4, 95, "90% target", adj = 0, col = "red", cex = 0.8)
}

# compare four estimates of each gene's mean using the readings and the fitted shift
plot_mean_comparison <- function(dat, farm) {
  ordinary <- colMeans(dat$X)
  huber <- apply(dat$X, 2, FarmTest::huber.mean)
  adjusted <- as.numeric(farm$means)
  shift <- huber - adjusted

  estimates <- data.frame(
    "Ordinary means" = ordinary,
    "Ordinary + adjustment" = ordinary - shift,
    "Huber means" = huber,
    "Huber + adjustment" = adjusted,
    check.names = FALSE
  )

  # stack the four columns so each histogram gets one row per gene
  hist_data <- stack(estimates)
  names(hist_data) <- c("estimate", "method")
  truth <- factor(dat$mu, levels = c(0, 0.6), labels = c("0 (null)", "0.6 (signal)"))
  hist_data$truth <- rep(truth, times = ncol(estimates))
  breaks <- pretty(range(hist_data$estimate), n = 30)

  ggplot(hist_data, aes(x = estimate, fill = truth)) +
    geom_histogram(breaks = breaks, position = "identity", alpha = 0.7, color = "white") +
    geom_vline(xintercept = c(0, 0.6), color = "red", linetype = "dashed") +
    facet_wrap(~method, ncol = 2) +
    scale_fill_manual(values = c("gray", "blue"), name = "True mean") +
    scale_y_continuous(breaks = c(0, 20, 40)) +
    labs(x = "Estimated mean", y = "Number of genes") +
    theme_minimal() +
    theme(legend.position = "bottom")
}
