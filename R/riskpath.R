#' Simulate portfolio paths using a moving block bootstrap
#'
#' @param returns Numeric matrix: trading days in rows, assets in columns.
#'   All entries must be finite and daily simple returns must exceed -1.
#' @param weights Numeric allocation vector (nonnegative, sums to one).
#' @param initial_value Positive initial portfolio value.
#' @param horizon Positive integer number of trading days to simulate.
#' @param n_sim Positive integer number of simulated trajectories.
#' @param block_size Positive integer length of each sampled historical block,
#'   not exceeding the number of historical trading days.
#' @return Numeric matrix with horizon + 1 rows and n_sim columns; row 1
#'   contains initial_value. Each column is one simulated portfolio path.
#' @details Samples blocks with replacement; preserves dependence within
#'   each block and contemporaneous cross-asset return relationships. Assumes
#'   constant target weights rebalanced daily, with no transaction costs.
#' @examples
#' set.seed(123)
#' r <- matrix(rnorm(300, 0, .01), ncol = 3)
#' p <- riskpath_sim(r, c(.5, .3, .2), 10000, 150, 100, 10)
#' dim(p)
riskpath_sim <- function(returns, weights, initial_value, horizon,
                         n_sim, block_size) {
  if (!is.matrix(returns) || !is.numeric(returns) ||
      nrow(returns) < 1L || ncol(returns) < 1L ||
      any(!is.finite(returns)) || any(returns <= -1)) {
    stop('returns must be a finite numeric matrix of simple returns > -1')
  }
  if (!is.numeric(weights) || length(weights) != ncol(returns) ||
      any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - 1) >= 1e-8) {
    stop('weights must be nonnegative, finite, and sum to one')
  }
  if (!is.numeric(initial_value) || length(initial_value) != 1L ||
      !is.finite(initial_value) || initial_value <= 0) {
    stop('initial_value must be a positive finite number')
  }
  is_positive_integer <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= 1 && x == floor(x)
  }
  if (!is_positive_integer(horizon) || !is_positive_integer(n_sim) ||
      !is_positive_integer(block_size) || block_size > nrow(returns)) {
    stop('horizon, n_sim, block_size must be positive integers; block_size <= nrow(returns)')
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
#'
#' @param paths Numeric matrix of strictly positive portfolio values, with
#'   one trajectory per column and time in rows (including initial value).
#' @param confidence Single probability strictly between 0 and 1 for VaR.
#' @param threshold Single drawdown fraction between 0 and 1 inclusive.
#' @return Named list: VaR (terminal loss quantile), ES (mean terminal loss
#'   at or above VaR), mean_max_decline (average maximum path drawdown), and
#'   probability_20pct_decline (fraction of paths with drawdown >= threshold).
#'   The final field retains its original name for compatibility, even when
#'   threshold is not 0.20. Losses can be negative when the portfolio gains.
#' @details No external packages needed; no side effects.
#' @examples
#' p <- rbind(c(100, 100), c(90, 110), c(80, 105))
#' riskpath_metrics(p)
riskpath_metrics <- function(paths, confidence = 0.95, threshold = 0.20) {
  if (!is.matrix(paths) || !is.numeric(paths) || nrow(paths) < 2L ||
      ncol(paths) < 1L || any(!is.finite(paths)) || any(paths <= 0)) {
    stop('paths must be a finite, positive numeric matrix with at least two rows')
  }
  if (!is.numeric(confidence) || length(confidence) != 1L ||
      !is.finite(confidence) || confidence <= 0 || confidence >= 1) {
    stop('confidence must be between 0 and 1 (exclusive)')
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0 || threshold > 1) {
    stop('threshold must be between 0 and 1')
  }
  losses <- 1 - paths[nrow(paths), ] / paths[1L, ]
  VaR <- unname(stats::quantile(losses, confidence))
  ES <- mean(losses[losses >= VaR])
  MD <- apply(paths, 2L, function(p) max(1 - p / cummax(p)))
  list(VaR = VaR, ES = ES, mean_max_decline = mean(MD),
       probability_20pct_decline = mean(MD >= threshold))
}
