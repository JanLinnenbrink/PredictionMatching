#' Weighted error statistics
#'
#' @param weights A `twcv` object, or a numeric vector of weights.
#' @param pointwise_error A vector of the error for each training point based on
#'   cross-validation (obs - pred).
#'
#' @return A named numeric vector of weighted error statistics
#'   (bias, mse, rmse, mae).
#'
#' @export
weighted_error_stats <- function(
  weights,
  pointwise_error
) {
  # Checks
  if (missing(pointwise_error) || is.null(pointwise_error)) {
    stop("`pointwise_error` must be supplied.")
  }
  if (!is.numeric(pointwise_error)) {
    stop("`pointwise_error` must be a numeric vector.")
  }

  # Extract the weight vector
  if (inherits(weights, "twcv")) {
    w <- weights$weights$weights
    if (is.null(w)) {
      stop("The 'twcv' object does not contain a valid weight vector.")
    }
  } else {
    if (!is.numeric(weights)) {
      stop("`weights` must be a 'twcv' object or a numeric vector.")
    }
    w <- weights
  }

  # Validate lengths
  if (length(w) != length(pointwise_error)) {
    stop("Length of the weights vector doesn't match the error data.")
  }

  # Drop NAs consistently across weights and errors
  valid <- !is.na(w) & !is.na(pointwise_error)
  if (sum(valid) == 0) {
    stop("No valid (non-NA) observations to compute statistics.")
  }

  w_v <- w[valid]
  err_v <- pointwise_error[valid]
  sum_w <- sum(w_v)

  if (sum_w == 0) {
    stop("Sum of weights is zero; cannot compute weighted statistics.")
  }

  # Weighted statistics
  bias <- sum(w_v * err_v) / sum_w
  mse <- sum(w_v * err_v^2) / sum_w
  rmse <- sqrt(mse)
  mae <- sum(w_v * abs(err_v)) / sum_w

  c(bias = bias, mse = mse, rmse = rmse, mae = mae)
}
