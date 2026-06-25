#' Plot the calibration results
#' @description Generic plot function for twcv class
#'
#' @name plot
#' @param x An object of type \emph{twcv}.
#' @param pointwise_error Optional. Either a `pe` object (from
#'   [calculate_pointwise_error()]) or a numeric vector of pointwise errors. If a
#'   numeric vector is supplied, it is assumed to be aligned with the weights by
#'   position. If a `pe` object is supplied, errors are aligned to the weights by
#'   ID. If `NULL`, only the calibration plot is returned.
#' @param ... other arguments.
#' @author Jan Linnenbrink
#'
#' @export
plot.twcv <- function(x, pointwise_error = NULL, ...) {
  w_list <- x$weights
  if (is.list(w_list[[1]])) {
    w <- w_list[[1]]$weights
  } else {
    w <- w_list$weights
  }

  # plot margin calibration
  cal <- .twcv_calibration_df(x$training_bal, x$prediction_margins, w)
  calibration_plot <- ggplot2::ggplot(
    cal,
    ggplot2::aes(
      y = .data[["var"]],
      x = abs(.data[["target"]] - .data[["weighted"]]),
      fill = as.factor(.data[["level"]])
    )
  ) +
    ggplot2::geom_point(shape = 21, colour = "grey", size = 3) +
    ggplot2::scale_fill_viridis_d("Quintile", option = 5) +
    ggplot2::labs(
      y = "",
      x = "abs(target proportion - weighted training proportion)",
      title = "A)  Margin calibration"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 0.8,
      panel.grid.minor = ggplot2::element_blank()
    )

  # stop if pointwise_error not supplied
  if (is.null(pointwise_error)) {
    return(calibration_plot)
  }

  # --- Resolve pointwise_error into a numeric vector aligned with w ---
  err <- .align_error_to_weights(pointwise_error, weight_ids = x$ids, w = w)

  # plot weight vs loss
  df <- data.frame(w = w, loss = abs(err))
  bias_plot <- ggplot2::ggplot(
    df,
    ggplot2::aes(.data[["w"]], .data[["loss"]])
  ) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      colour = "firebrick"
    ) +
    ggplot2::labs(
      x = "Weight",
      y = "Pointwise loss (absolute error)",
      title = "B) Weight vs. loss"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 0.8,
      panel.grid.minor = ggplot2::element_blank()
    )

  return(list(calibration_plot, bias_plot))
}

#' Align a pointwise-error input to a weight vector
#'
#' Accepts either a `pe` object (data frame with `id`/`error`) or a plain
#' numeric vector, and returns a numeric error vector aligned with the weights.
#'
#' @param pointwise_error A `pe` object or numeric vector.
#' @param weight_ids IDs corresponding to the weights (e.g. `x$ids`).
#' @param w The weight vector (used for length checks).
#'
#' @return Numeric vector of errors aligned with `w`.
#' @keywords internal
#' @noRd
.align_error_to_weights <- function(pointwise_error, weight_ids, w) {
  # pe object / data frame -> align by ID
  if (inherits(pointwise_error, "pe") || is.data.frame(pointwise_error)) {
    if (!all(c("id", "error") %in% names(pointwise_error))) {
      stop("A 'pe' object must contain 'id' and 'error' columns.")
    }
    error_ids <- pointwise_error$id
    err <- pointwise_error$error

    if (is.null(weight_ids)) {
      stop("No weight IDs available; cannot align a 'pe' object by ID.")
    }

    idx <- match(weight_ids, error_ids)
    if (anyNA(idx)) {
      stop("Some weight IDs were not found in the error IDs; cannot align.")
    }
    return(err[idx])
  }

  # plain numeric vector -> positional matching
  if (!is.numeric(pointwise_error)) {
    stop("`pointwise_error` must be a 'pe' object or a numeric vector.")
  }
  if (length(pointwise_error) != length(w)) {
    stop("Length of `pointwise_error` doesn't match the number of weights.")
  }
  pointwise_error
}

#' @keywords internal
#' @noRd
.twcv_calibration_df <- function(training_bal, prediction_margins, w) {
  vars <- names(prediction_margins)
  dplyr::bind_rows(lapply(vars, function(m) {
    levs <- seq_along(prediction_margins[[m]])
    x <- as.integer(training_bal[[m]])
    f <- factor(x, levels = levs)

    wt <- tapply(w, f, sum)
    wt[is.na(wt)] <- 0
    uw <- tapply(rep(1, length(x)), f, sum)
    uw[is.na(uw)] <- 0

    data.frame(
      var = m,
      level = levs,
      target = as.numeric(prediction_margins[[m]]),
      weighted = as.numeric(wt) / sum(wt),
      unweighted = as.numeric(uw) / sum(uw),
      stringsAsFactors = FALSE
    )
  }))
}
