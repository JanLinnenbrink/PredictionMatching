#' Prediction-domain adaptive weighting based on raking
#'
#' @param tpoints data.frame or sf object containing the predictor values at the training points.
#' @param modeldomain SpatRaster containing the predictors. Not needed when predpoints are supplied.
#' @param predpoints data.frame or sf object containing the predictor values at the prediction points.
#' Only needed if no modeldomain is supplied.
#' @param pointwise_error vector of the error estimates for every training points (aligned with the training points). Optional.
#' If supplied, the weights will be applied to the pointwise error, and an weighted error will be returned.
#' @param samplesize numeric. How many points in the modeldomain should be sampled as prediction points?
#' Only required if modeldomain is used instead of predpoints.
#' @param sampling character. How to draw prediction points from the modeldomain? See `sf::st_sample`.
#' Only required if modeldomain is used instead of predpoints.
#' @param balance_by Numeric step size for quantile binning.
#' @param shrink_lambda Shrinkage parameter applied after calibration raking.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{weights}{Named list of weight object for TWCV.}
#'     \item{weighted_error}{vector of the weighted error. Only if pointwise_error is supplied.}
#'     \item{unsupported_flag}{1 if some quintiles are not supported by the training data. 0 otherwise.}
#'     \item{unsupported_vars}{Vector containing the names of variables with quintiles not supported by the training data.}
#'     \item{training_bal}{data.frame containing the discretized predictors for the training data.}
#'     \item{prediction_margins}{data.frame containing the frequency of the discretized predictors for the prediction data}
#'   }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' library(ggplot2)
#' library(CAST)
#'
#' data(splotdata)
#' splotdata <- splotdata[splotdata$Country == "Chile",]
#'
#' predictors <- c("bio_1", "bio_4", "bio_5", "bio_6",
#'                "bio_8", "bio_9", "bio_12", "bio_13",
#'                "bio_14", "bio_15", "elev")
#'
#' trainDat <- sf::st_drop_geometry(splotdata)
#' predictors_sp <- terra::rast(system.file("extdata", "predictors_chile.tif",package="CAST"))
#' terra::plot(predictors_sp[["bio_1"]])
#' terra::plot(vect(splotdata), add = T)
#'
#' pointwise_error <- rnorm(nrow(trainDat), 0, 1)
#' w <-calculate_weights(tpoints = trainDat[,predictors], modeldomain = predictors_sp,
#'                       pointwise_error = pointwise_error)
#' plot(w)
#' }
#'
#' @export
#'
calculate_weights <- function(
  tpoints,
  modeldomain = NULL,
  predpoints = NULL,
  pointwise_error = NULL,
  samplesize = 1000,
  sampling = "regular",
  balance_by = 0.2,
  shrink_lambda = 0.2
) {
  # Checks
  if (shrink_lambda < 0 || shrink_lambda > 1) {
    stop("shrink_lambda must be in the range 0,1")
  }

  # Sample prediction points if not supplied
  if (is.null(predpoints) & !is.null(modeldomain)) {
    predpoints <- .generate_predpoints(
      modeldomain = modeldomain,
      samplesize = samplesize,
      sampling = sampling
    )
  }

  # Standardize the inputs (prediction points and training points data frames)
  if (any(class(tpoints) %in% c("sf", "sfc"))) {
    train_dat <- sf::st_drop_geometry(tpoints)
  } else {
    train_dat <- tpoints
  }

  if (any(class(predpoints) %in% c("sf", "sfc"))) {
    pred_dat <- sf::st_drop_geometry(predpoints)
  } else {
    pred_dat <- predpoints
  }

  balancing_vars <- names(pred_dat)

  train_dat$id <- seq_len(nrow(train_dat))
  pred_dat$id <- seq_len(nrow(pred_dat))

  # Check if the names of the training and prediction data match
  if (!setequal(names(train_dat), names(pred_dat))) {
    stop(
      "tpoints and predpoints (or the modeldomain) need to contain the predictor data and have the same column names."
    )
  }

  # Construct balanced representations of sample and grid task descriptors for weighting
  train_dat_bal <- .prepare_for_balancing(
    df = train_dat,
    vars = balancing_vars,
    ref_df = pred_dat,
    by = balance_by
  )

  pred_dat_bal <- .prepare_for_balancing(
    df = pred_dat,
    vars = balancing_vars,
    ref_df = pred_dat,
    by = balance_by
  )

  # Extracts the discretized variables from the training points
  train_dat_bal_df <- as.data.frame(
    lapply(balancing_vars, function(v) train_dat_bal[[paste0(v, "_cat")]])
  )
  names(train_dat_bal_df) <- balancing_vars

  # Calculates proportion of predpoints in each quantile of each predictor used for weighting
  prediction_margins <- .compute_prediction_margins(
    pred_task_bal = pred_dat_bal,
    balancing_vars = balancing_vars
  )

  # Applies iterative proportional fitting ("raking")
  tw <- .rake_weights(
    balance_df = train_dat_bal_df,
    prediction_margins = prediction_margins
  )

  # Normalize weights by their mean
  tw$weights_raw <- tw$weights / mean(tw$weights)

  # Shrink the normalized weights towards 1 to mitigate extreme values
  tw$weights <- (1 - shrink_lambda) * tw$weights_raw + shrink_lambda

  # Check for any unsupported quintiles and issue a warning
  support_check <- .check_balance_support(train_dat_bal_df, prediction_margins)
  if (any(support_check$unsupported)) {
    unsupported_vars <- unique(support_check[
      support_check$unsupported == TRUE,
      "var"
    ])
    unsupported_flag <- 1
    warning(paste0(
      "The predictor(s) ",
      paste0(unsupported_vars, collapse = ","),
      " have quintiles that are not supported by the training data.
				Raking is likely to fail in this context, and limiting the prediction area to avoid extrapolation is recommended."
    ))
  } else {
    unsupported_vars <- NA
    unsupported_flag <- 0
  }

  # Calculate the weighted pointwise error if supplied
  if (is.null(pointwise_error)) {
    weighted_error <- NA
  } else {
    weighted_error <- pointwise_error * tw$weights
  }

  res <- list(
    weights = tw,
    weighted_error = weighted_error,
    unsupported_flag = unsupported_flag,
    unsupported_vars = unsupported_vars,
    training_bal = train_dat_bal_df,
    prediction_margins = prediction_margins
  )
  class(res) <- "twcv"
  return(res)
}


#' Samples prediction points if not supplied
#'
#' @param modeldomain SpatRaster object containing the predictor stack.
#'
#' @return Named list of empirical target margins.
#' @noRd
.generate_predpoints <- function(modeldomain, samplesize, sampling) {
  # Check modeldomain is indeed a sf/SpatRaster
  if (!any(c("SpatRaster") %in% class(modeldomain))) {
    stop("modeldomain must be a 'SpatRaster' object.")
  }

  # If modeldomain is a SpatRaster, transform into polygon
  predictor_stack <- modeldomain

  modeldomain[!is.na(modeldomain)] <- 1
  modeldomain <- terra::as.polygons(
    modeldomain,
    values = FALSE,
    na.all = TRUE
  ) |>
    sf::st_as_sf() |>
    sf::st_union()

  # Check modeldomain is indeed a polygon sf
  if (
    !any(
      class(sf::st_geometry(modeldomain)) %in%
        c("sfc_POLYGON", "sfc_MULTIPOLYGON")
    )
  ) {
    stop("modeldomain must be a sf/sfc polygon object.")
  }

  # We sample
  message(paste0(
    samplesize,
    " prediction points are sampled from the modeldomain"
  ))
  predpoints <- suppressMessages(sf::st_sample(
    x = modeldomain,
    size = samplesize,
    type = sampling
  ))
  sf::st_crs(predpoints) <- sf::st_crs(modeldomain)

  message("predictor values are extracted for prediction points")
  predpoints <- terra::extract(
    predictor_stack,
    terra::vect(predpoints),
    ID = FALSE
  )
  return(predpoints)
}


#' Prepare variables for balancing via discretization
#'
#' Transforms variables into categorical representations suitable for
#' balancing. Numeric variables are discretized using quantiles of a
#' reference distribution; categorical variables are aligned to the
#' reference support.
#'
#' @param df Data frame to transform.
#' @param vars Variables to transform.
#' @param ref_df Reference data frame.
#' @param by Quantile step size.
#'
#' @return Modified data frame with additional \code{*_cat} variables.
#' @noRd
.prepare_for_balancing <- function(df, vars, ref_df, by = 0.2) {
  df_out <- df

  for (v in vars) {
    if (!(v %in% names(df))) {
      stop("Variable '", v, "' not found in df.", call. = FALSE)
    }
    if (!(v %in% names(ref_df))) {
      stop("Variable '", v, "' not found in ref_df.", call. = FALSE)
    }

    x_ref <- ref_df[[v]]
    x <- df[[v]]
    out_name <- paste0(v, "_cat")

    # numeric variables
    if (is.numeric(x_ref) && length(unique(stats::na.omit(x_ref))) > 2) {
      probs <- seq(0, 1, by = by)
      qtiles <- stats::quantile(
        x_ref,
        probs = probs,
        na.rm = TRUE,
        names = FALSE
      )

      qtiles <- unique(qtiles)

      # degenerate case
      if (length(qtiles) < 2) {
        levs <- paste0(v, "_Q1")
        df_out[[out_name]] <- factor(
          ifelse(is.na(x), NA, levs),
          levels = levs,
          ordered = TRUE
        )
        next
      }

      qtiles[1] <- -Inf
      qtiles[length(qtiles)] <- Inf

      n_bins <- length(qtiles) - 1L
      levs <- paste0(v, "_Q", seq_len(n_bins))

      df_out[[out_name]] <- cut(
        x,
        breaks = qtiles,
        include.lowest = TRUE,
        labels = levs,
        ordered_result = TRUE
      )
    } else {
      ref_levels <- sort(unique(as.character(stats::na.omit(x_ref))))
      x_chr <- as.character(x)

      x_chr[!(x_chr %in% ref_levels) & !is.na(x_chr)] <- NA_character_

      df_out[[out_name]] <- factor(
        x_chr,
        levels = ref_levels,
        ordered = FALSE
      )
    }
  }

  df_out
}


#' Compute the margins of the prediction points
#'
#' Extracts empirical target-domain margins for discretized balancing variables.
#'
#' @param pred_task_bal Data frame of discretized deployment-task descriptors.
#' @param balancing_vars Character vector of balancing-variable names. For each
#'   variable `v`, the function expects a column named `paste0(v, "_cat")`.
#'
#' @return Named list of empirical target margins.
#' @noRd
.compute_prediction_margins <- function(pred_task_bal, balancing_vars) {
  out <- vector("list", length(balancing_vars))
  names(out) <- balancing_vars

  for (v in balancing_vars) {
    vn <- paste0(v, "_cat")
    if (!(vn %in% names(pred_task_bal))) {
      stop("Missing column in pred_task_bal: ", vn, call. = FALSE)
    }
    pred_task_bal[[vn]] <- as.integer(pred_task_bal[[vn]])
    levels <- sort(unique(pred_task_bal[[vn]]))

    freq_table <- table(factor(pred_task_bal[[vn]], levels = levels))
    out[[v]] <- as.numeric(freq_table) / sum(freq_table)
  }

  out
}


#' Compute calibration weights by iterative proportional fitting
#'
#' Reweights validation tasks so that weighted empirical margins match target
#' margins for a set of discretized balancing variables.
#'
#' @param balance_df Data frame of discretized balancing variables.
#' @param prediction_margins Named list of target proportions for each balancing
#'   variable.
#' @param base_weights Optional numeric vector of starting weights.
#' @param max_iter Maximum number of raking iterations.
#' @param tol Convergence tolerance based on relative weight change.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{weights}{Final calibration weights.}
#'     \item{converged}{Logical indicating whether convergence was reached.}
#'     \item{iterations}{Number of iterations performed.}
#'   }
#' @noRd
.rake_weights <- function(
  balance_df,
  prediction_margins,
  base_weights = NULL,
  max_iter = 500,
  tol = 1e-6
) {
  n <- nrow(balance_df)
  if (is.null(base_weights)) {
    w <- rep(1, n)
  } else {
    w <- as.numeric(base_weights)
  }

  margin_names <- names(prediction_margins)

  for (iter in seq_len(max_iter)) {
    w_old <- w

    for (m in margin_names) {
      x <- as.integer(balance_df[[m]])
      levs <- seq_along(prediction_margins[[m]])
      target_prop <- prediction_margins[[m]]

      # calculate the number of training points currently in each quantile
      # uses the weights from the previous balancing variable
      # -> already weighted, but likely not ideally for this variable
      current_totals <- tapply(w, factor(x, levels = levs), sum)
      current_totals[is.na(current_totals)] <- 0

      # calculate the desired number of training points in each quantile that matches the target margins:
      # Number of training points * target proportion vector
      target_totals <- sum(w) * target_prop

      adj <- rep(1, length(levs))
      ok <- current_totals > 0

      # calculates the weight needed to adjust the training point distribution to the target margins
      # (weights > 1 are used to up-weigh underrepresented classes, weights < 1 to down-weight over-represented ones)
      adj[ok] <- target_totals[ok] / current_totals[ok]
      adj[!ok] <- NA_real_

      valid <- !is.na(adj[x])
      w[valid] <- w[valid] * adj[x[valid]]
    }

    # Compute the relative strength of the absolute change of weights from base (or previous) to new weights
    rel_change <- max(abs(w - w_old) / pmax(abs(w_old), 1e-12))

    # When the changes converge (i.e., when the weight difference from previous iteration to new one is small), stop and return weights
    # (when the weights for predictor B are changed, weights for predictor A might be off again, and another iteration starts,
    # until the diff between them is small)
    if (rel_change < tol) break
  }

  converged <- iter < max_iter || rel_change < tol

  list(
    weights = w,
    converged = converged,
    iterations = iter
  )
}


#' Checks if all quintiles of the prediction points are supported by the training points.
#' Otherwise, raking is prone to errors.
#'
#' @param balance_df data frame containing the training margins.
#' @param prediction_margins data frame containing the prediction point margins.
#' @param eps Tolerance
#'
#' @return data frame containing the predictors, their quintiles and information if they are supported.
#' @noRd
.check_balance_support <- function(
  balance_df,
  prediction_margins,
  eps = 1e-12
) {
  out <- lapply(names(prediction_margins), function(m) {
    levs <- seq_along(prediction_margins[[m]])

    sample_counts <- table(
      factor(as.integer(balance_df[[m]]), levels = levs)
    )

    data.frame(
      var = m,
      level = levs,
      sample_n = as.numeric(sample_counts),
      target_prop = as.numeric(prediction_margins[[m]]),
      unsupported = as.numeric(sample_counts) == 0 &
        as.numeric(prediction_margins[[m]]) > eps
    )
  })

  dplyr::bind_rows(out)
}
