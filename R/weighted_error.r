#' Weighted error statistics
#'
#' @param weights A `twcv` object, or a numeric vector of weights.
#' @param pointwise_error Either a `pe` object (from
#'   [calculate_pointwise_error()], i.e. a data frame with `id` and `error`
#'   columns) or a numeric vector of pointwise errors (obs - pred). When a `pe`
#'   object is supplied, errors are aligned to the weights by ID. When a plain
#'   numeric vector is supplied, errors are assumed to be aligned with the
#'   weights by position.
#'
#' @details
#' When `weights` is a `twcv` object it carries IDs, and `pointwise_error`
#' should be a `pe` object so that errors can be aligned by ID. This prevents
#' silent misalignment (e.g., when cross-validation predictions are returned in
#' fold order rather than the original row order). Positional matching is only
#' used when both `weights` and `pointwise_error` are plain vectors.
#'
#' @return A named numeric vector of weighted error statistics
#'   (bias, mse, rmse, mae).
#'
#' @seealso [calculate_weights()], [calculate_pointwise_error()]
#'
#' @export
weighted_error_stats <- function(
  weights,
  pointwise_error
) {
  if (missing(pointwise_error) || is.null(pointwise_error)) {
    stop("`pointwise_error` must be supplied.")
  }

  # Extract the weight vector and its IDs
  if (inherits(weights, "twcv")) {
    w <- weights$weights$weights
    if (is.null(w)) {
      stop("The 'twcv' object does not contain a valid weight vector.")
    }
    weight_ids <- weights$ids
  } else {
    if (!is.numeric(weights)) {
      stop("`weights` must be a 'twcv' object or a numeric vector.")
    }
    w <- weights
    weight_ids <- NULL
  }

  # Align errors to weights (by ID for pe objects, positionally for vectors)
  err <- .align_error_to_weights(
    pointwise_error,
    weight_ids = weight_ids,
    w = w
  )

  # Drop non-finite values consistently across weights and errors
  valid <- is.finite(w) & is.finite(err)
  if (sum(valid) == 0) {
    stop("No valid (finite) observations to compute statistics.")
  }

  w_v <- w[valid]
  err_v <- err[valid]
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
