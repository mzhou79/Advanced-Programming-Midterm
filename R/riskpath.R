# Portfolio moving-block bootstrap: simulation, risk metrics, and plots.
# Dependencies: ggplot2 and scales for plotting only.

#' Simulate portfolio paths using a moving block bootstrap
#' @param returns Numeric matrix of daily simple returns (days x assets), finite and > -1.
#' @param weights Nonnegative numeric allocation vector summing to one.
#' @param initial_value Positive initial portfolio value.
#' @param horizon Positive integer number of trading days to simulate.
#' @param n_sim Positive integer number of trajectories.
#' @param block_size Positive integer block length, at most nrow(returns).
#' @return Numeric matrix (horizon + 1) x n_sim; first row is initial value.
#' @details Samples consecutive blocks with replacement, preserving within-block
#'   serial dependence and same-day cross-asset dependence. Assumes daily
#'   rebalancing to constant weights, no fees, and no transaction costs.
#' @examples
#' set.seed(123)
#' r <- matrix(rnorm(300, 0, 0.01), ncol = 3)
#' p <- riskpath_sim(r, c(.5, .3, .2), 10000, 150, 100, 10)
#' dim(p)
riskpath_sim <- function(returns, weights, initial_value, horizon,
                         n_sim, block_size) {
  if (!is.matrix(returns) || !is.numeric(returns) ||
      nrow(returns) < 1L || ncol(returns) < 1L ||
      any(!is.finite(returns)) || any(returns <= -1)) {
    stop("returns must be a finite numeric matrix of simple returns > -1")
  }
  if (!is.numeric(weights) || length(weights) != ncol(returns) ||
      any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - 1) >= 1e-8) {
    stop("weights must be nonnegative, finite, and sum to one")
  }
  if (!is.numeric(initial_value) || length(initial_value) != 1L ||
      !is.finite(initial_value) || initial_value <= 0) {
    stop("initial_value must be a positive finite number")
  }
  is_positive_integer <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= 1 && x == floor(x)
  }
  if (!is_positive_integer(horizon) || !is_positive_integer(n_sim) ||
      !is_positive_integer(block_size) || block_size > nrow(returns)) {
    stop("horizon, n_sim, block_size must be positive integers; block_size <= nrow(returns)")
  }
  n <- nrow(returns)
  n_blocks <- ceiling(horizon / block_size)
  paths <- matrix(NA_real_, nrow = horizon + 1L, ncol = n_sim)
  paths[1L, ] <- initial_value
  for (s in seq_len(n_sim)) {
    starts <- sample.int(n - block_size + 1L, n_blocks, replace = TRUE)
    idx <- unlist(lapply(starts, function(i) i:(i + block_size - 1L)),
                  use.names = FALSE)
    idx <- idx[seq_len(horizon)]
    portfolio_returns <- as.vector(returns[idx, , drop = FALSE] %*% weights)
    paths[-1L, s] <- initial_value * cumprod(1 + portfolio_returns)
  }
  paths
}

#' Summarize simulated portfolio risk
#' @param paths Positive numeric matrix of portfolio values (time x simulations).
#' @param confidence Probability strictly between 0 and 1 used for VaR.
#' @param threshold Drawdown threshold between 0 and 1 inclusive.
#' @return List with VaR, ES, mean_max_decline, probability_20pct_decline.
#'   Last name is retained for backward compatibility even if threshold != 0.20.
#' @details Losses are 1 - terminal/initial and may be negative for gains.
#'   Maximum drawdown is relative to the running peak, including day zero.
#'   No external packages or side effects.
#' @examples
#' p <- rbind(c(100, 100), c(90, 110), c(80, 105))
#' riskpath_metrics(p)
riskpath_metrics <- function(paths, confidence = 0.95, threshold = 0.20) {
  .validate_paths(paths)
  .validate_probability(confidence, "confidence", inclusive = FALSE)
  .validate_probability(threshold, "threshold", inclusive = TRUE)
  losses <- 1 - paths[nrow(paths), ] / paths[1L, ]
  VaR <- unname(stats::quantile(losses, confidence))
  ES <- mean(losses[losses >= VaR])
  MD <- apply(paths, 2L, function(p) max(1 - p / cummax(p)))
  list(VaR = VaR, ES = ES, mean_max_decline = mean(MD),
       probability_20pct_decline = mean(MD >= threshold))
}

# Internal helpers: validate inputs shared by risk metrics and plotting functions.
.validate_paths <- function(paths) {
  if (!is.matrix(paths) || !is.numeric(paths) || nrow(paths) < 2L ||
      ncol(paths) < 1L || any(!is.finite(paths)) || any(paths <= 0)) {
    stop("paths must be a finite, positive numeric matrix with at least two rows")
  }
  invisible(TRUE)
}
.validate_probability <- function(x, name, inclusive = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      (if (inclusive) x < 0 || x > 1 else x <= 0 || x >= 1)) {
    stop(sprintf("%s must be %s", name,
                 if (inclusive) "between 0 and 1" else "strictly between 0 and 1"))
  }
  invisible(TRUE)
}
.require_plot_packages <- function() {
  for (pkg in c("ggplot2", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Install '%s' to use plotting functions", pkg))
  }
}

#' Plot simulated portfolio paths and terminal-value percentiles
#' @param paths Positive numeric matrix returned by riskpath_sim().
#' @return A ggplot object, without printing or writing files.
#' @details Draws all simulated paths faintly, plus straight reference segments
#'   from initial value to the 5th, 25th, 50th, 75th, and 95th terminal
#'   percentiles. Segments are percentile guides, not actual simulated paths.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("scales", quietly = TRUE)) {
#'   plot_riskpath_paths(rbind(c(100, 100), c(90, 110), c(80, 105)))
#' }
plot_riskpath_paths <- function(paths) {
  .validate_paths(paths)
  .require_plot_packages()
  n_days <- nrow(paths) - 1L
  terminal <- paths[nrow(paths), ]
  percentile_labels <- c("5%", "25%", "50%", "75%", "95%")
  percentile_values <- as.numeric(stats::quantile(
    terminal, probs = c(.05, .25, .50, .75, .95)))
  plot_data <- data.frame(
    day = rep(0:n_days, times = ncol(paths)),
    value = as.vector(paths),
    simulation = rep(seq_len(ncol(paths)), each = nrow(paths))
  )
  line_data <- data.frame(
    percentile = factor(percentile_labels, levels = percentile_labels),
    start_day = 0, end_day = n_days,
    start_value = paths[1L, 1L], end_value = percentile_values
  )
  line_data$label <- paste0(line_data$percentile, ": $",
                            format(round(line_data$end_value),
                                   big.mark = ",", scientific = FALSE, trim = TRUE))
  ggplot2::ggplot() +
    ggplot2::geom_line(data = plot_data,
                       ggplot2::aes(x = day, y = value, group = simulation),
                       alpha = .01, color = "steelblue", linewidth = .3) +
    ggplot2::geom_segment(data = line_data,
                          ggplot2::aes(x = start_day, y = start_value, xend = end_day,
                                       yend = end_value, color = percentile), linewidth = 1.2) +
    ggplot2::geom_text(data = line_data,
                       ggplot2::aes(x = end_day, y = end_value, label = label,
                                    color = percentile), hjust = -.1, size = 3.2,
                       show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = paths[1L, 1L],
                        linetype = "dashed", color = "gray50") +
    ggplot2::scale_color_manual(values = c(
      "5%" = "red", "25%" = "orange", "50%" = "black",
      "75%" = "darkgreen", "95%" = "purple")) +
    ggplot2::scale_x_continuous(
      limits = c(0, n_days + max(30, ceiling(n_days * .4))),
      breaks = pretty(c(0, n_days), n = 4)) +
    ggplot2::scale_y_continuous(labels = scales::label_dollar()) +
    ggplot2::labs(title = "Simulated Portfolio Paths with Terminal Percentiles",
                  subtitle = paste(ncol(paths), "simulations using moving block bootstrap"),
                  x = "Trading Day", y = "Portfolio Value ($)",
                  color = "Terminal Percentile") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 70, 10, 10))
}

#' Plot the distribution of terminal losses with a VaR cutoff
#' @param paths Positive numeric matrix returned by riskpath_sim().
#' @param confidence VaR confidence level strictly between 0 and 1.
#' @param bins Number of histogram bins (positive integer).
#' @return A ggplot object, without printing or writing files.
#' @details Losses can be negative (a gain). The red line is the loss quantile;
#'   observations strictly beyond the cutoff are counted in the annotation.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("scales", quietly = TRUE)) {
#'   plot_riskpath_var(rbind(c(100, 100), c(90, 110), c(80, 105)))
#' }
plot_riskpath_var <- function(paths, confidence = 0.95, bins = 60L) {
  .validate_paths(paths)
  .validate_probability(confidence, "confidence")
  if (length(bins) != 1L || !is.numeric(bins) || !is.finite(bins) ||
      bins < 1 || bins != floor(bins)) stop("bins must be a positive integer")
  .require_plot_packages()
  losses <- 1 - paths[nrow(paths), ] / paths[1L, ]
  VaR <- unname(stats::quantile(losses, confidence))
  n_exceed <- sum(losses > VaR)
  pct_exceed <- mean(losses > VaR) * 100
  loss_data <- data.frame(loss = losses)
  ggplot2::ggplot(loss_data, ggplot2::aes(x = loss)) +
    ggplot2::geom_histogram(ggplot2::aes(fill = ggplot2::after_stat(x > VaR)),
                            bins = bins, color = "white") +
    ggplot2::geom_vline(xintercept = VaR, color = "red",
                        linetype = "dashed", linewidth = 1.2) +
    ggplot2::annotate("text", x = VaR, y = Inf,
                      label = sprintf("%.0f%% VaR = %.2f%%\n%d simulations (%.2f%%)",
                                      confidence * 100, VaR * 100, n_exceed, pct_exceed),
                      hjust = -.05, vjust = 1.5, color = "red", size = 3.6) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "tomato"),
                               labels = c("FALSE" = "Below VaR", "TRUE" = "Beyond VaR"),
                               name = "Loss category") +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(title = "Distribution of Simulated Portfolio Losses",
                  subtitle = paste(ncol(paths), "simulated portfolio outcomes"),
                  x = "Terminal Portfolio Loss (%)", y = "Number of Simulations") +
    ggplot2::theme_minimal()
}

#' Plot the distribution of maximum portfolio drawdowns
#' @param paths Positive numeric matrix returned by riskpath_sim().
#' @param threshold Drawdown cutoff in [0, 1], default 20 percent.
#' @param bins Number of histogram bins (positive integer).
#' @return A ggplot object, without printing or writing files.
#' @details The black line marks mean maximum drawdown, the red line the
#'   threshold. Each simulation contributes its largest peak-to-trough loss.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("scales", quietly = TRUE)) {
#'   plot_riskpath_drawdown(rbind(c(100, 100), c(90, 110), c(80, 105)))
#' }
plot_riskpath_drawdown <- function(paths, threshold = 0.20, bins = 60L) {
  .validate_paths(paths)
  .validate_probability(threshold, "threshold", inclusive = TRUE)
  if (length(bins) != 1L || !is.numeric(bins) || !is.finite(bins) ||
      bins < 1 || bins != floor(bins)) stop("bins must be a positive integer")
  .require_plot_packages()
  max_drawdowns <- apply(paths, 2L, function(p) max(1 - p / cummax(p)))
  mean_mdd <- mean(max_drawdowns)
  n_exceed <- sum(max_drawdowns >= threshold)
  prob_exceed <- mean(max_drawdowns >= threshold)
  drawdown_data <- data.frame(max_drawdown = max_drawdowns)
  ggplot2::ggplot(drawdown_data, ggplot2::aes(x = max_drawdown)) +
    ggplot2::geom_histogram(
      ggplot2::aes(fill = ggplot2::after_stat(x >= threshold)),
      bins = bins, color = "white") +
    ggplot2::geom_vline(xintercept = mean_mdd, color = "black",
                        linetype = "dashed", linewidth = 1.2) +
    ggplot2::geom_vline(xintercept = threshold, color = "red",
                        linetype = "dashed", linewidth = 1.2) +
    ggplot2::annotate("text", x = mean_mdd, y = Inf,
                      label = sprintf("Mean MDD = %.2f%%", mean_mdd * 100),
                      hjust = 1.05, vjust = 1.5, color = "black", size = 3.6) +
    ggplot2::annotate("text", x = threshold, y = Inf,
                      label = sprintf("Threshold = %.0f%%\n%d simulations (%.2f%%)",
                                      threshold * 100, n_exceed, prob_exceed * 100),
                      hjust = -.05, vjust = 3, color = "red", size = 3.6) +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "tomato"),
      labels = c("FALSE" = "Below threshold", "TRUE" = "At or above threshold"),
      name = "Decline category") +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(title = "Distribution of Maximum Portfolio Drawdowns",
                  subtitle = paste(ncol(paths), "simulated portfolio trajectories"),
                  x = "Maximum Drawdown (%)", y = "Number of Simulations") +
    ggplot2::theme_minimal()
}
