#' Plot the calibration results
#' @description Generic plot function for twcv class
#'
#' @name plot
#' @param x An object of type \emph{twcv}.
#' @param ... other arguments.
#' @author Jan Linnenbrink
#'
#' @export
plot.twcv <- function(x, ...) {
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

  # stop if weighted_error not supplied
  if (is.na(x$weighted_error)) {
    return(calibration_plot)
  }

  # plot weight vs loss (only if weighted_error is supplied)
  df <- data.frame(w = w, loss = x$weighted_error)
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
      y = paste0("Pointwise loss (", "squared error", ")"),
      title = "B) Weight vs. loss"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 0.8,
      panel.grid.minor = ggplot2::element_blank()
    )
  return(list(calibration_plot, bias_plot))
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
