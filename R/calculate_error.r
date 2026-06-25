#' Compute pointwise cross-validation errors with IDs
#'
#' Computes residuals (obs - pred) for each cross-validation point and returns
#' them together with their identifiers, so that errors can be safely aligned
#' with weights from [calculate_weights()] regardless of row ordering.
#'
#' @param obs Numeric vector of observed values.
#' @param pred Numeric vector of predicted values.
#' @param id Vector of identifiers, one per observation, matching the IDs used
#'   in the corresponding `twcv` weight object. For caret output, this is
#'   typically `cv_pred$rowIndex`.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{id}{Identifier for each point.}
#'     \item{obs}{Observed value.}
#'     \item{pred}{Predicted value.}
#'     \item{error}{Residual (obs - pred).}
#'   }
#'
#' @seealso [calculate_weights()], [weighted_error_stats()]
#' @export
calculate_pointwise_error <- function(obs, pred, id) {
  if (!is.numeric(obs) || !is.numeric(pred)) {
    stop("`obs` and `pred` must be numeric vectors.")
  }
  if (length(obs) != length(pred)) {
    stop("`obs` and `pred` must have the same length.")
  }
  if (missing(id) || is.null(id)) {
    stop("`id` must be supplied so errors can be aligned with weights.")
  }
  if (length(id) != length(obs)) {
    stop("`id` must be the same length as `obs` and `pred`.")
  }
  if (anyDuplicated(id)) {
    stop("`id` values must be unique.")
  }

  res <- data.frame(
    id = id,
    obs = obs,
    pred = pred,
    error = obs - pred
  )
  class(res) <- c("pe", class(res))
  return(res)
}
