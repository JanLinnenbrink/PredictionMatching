#' @keywords internal
#' @noRd

generate_rast <- function() {
  set.seed(10)

  rast_grid <- terra::rast(
    xmin = 0,
    xmax = 200,
    ymin = 0,
    ymax = 200,
    ncols = 200,
    nrows = 200
  )

  grad_predictors <- simsam::sim_covariates(
    rast_grid,
    n = 7,
    method = simsam::simulate_gaussian(
      nugget = 0,
      beta = 50,
      psill = 100,
      model = "Exp",
      range = 50
    )
  )
  elev <- generate_elevation(rast_grid = rast_grid) |>
    rescale_raster(to = c(1, 100))
  grad_predictors <- c(grad_predictors, elev)

  landcover <- simsam::sim_covariates(
    rast_grid,
    n = 2,
    method = simsam::simulate_gaussian(
      psill = 5,
      model = "Exp",
      range = 100,
      beta = 10,
      indicators = TRUE
    )
  )

  predictors <- c(grad_predictors, landcover)
  names(predictors) <- c(
    "temp",
    "moisture",
    "ph",
    "slope",
    "solar",
    "dist_road",
    "prod",
    "elev",
    "forest",
    "grass"
  )

  outcome_pred <- simsam::blend_rasters(
    predictors,
    ~
      # species/habitat suitability score (unscaled)
      1.2 *
        exp(-((ph - 6.5)^2) / (2 * 0.6^2)) *
        (1 - exp(-moisture / 0.25)) + # pH optimum × moisture saturation
        0.5 * exp(-((temp - 15)^2) / (2 * 5^2)) * exp(-elev / 150) + # mild temp optimum, penalised by elevation
        0.3 * (forest / (1 + exp(-2 * (forest - 0.5)))) -
        0.2 * (grass / (1 + exp(-2 * (grass - 0.5)))) +
        0.05 * slope
    # forest boosts, grass reduces, small slope term
  )

  outcome_signal <- rescale_raster(outcome_pred, to = c(1, 100))

  noise <- simsam::sim_covariates(
    rast_grid,
    n = 1,
    method = simsam::simulate_gaussian(psill = 10, model = "Exp", range = 5)
  )
  names(noise) <- "noise"
  outcome <- outcome_signal + noise

  r <- c(predictors, outcome)
  terra::crs(r) <- "EPSG:3857"
  return(r)
}


#' @keywords internal
#' @noRd
generate_elevation <- function(rast_grid) {
  # Helper: standardize raster values
  standardize_raster <- function(x) {
    v <- terra::values(x)[, 1]
    terra::values(x) <- as.numeric(scale(v))
    x
  }

  # Coordinates
  xy <- terra::crds(rast_grid, df = TRUE)
  x <- xy$x
  y <- xy$y

  # 1. Broad-scale elevation pattern: lowlands to mountains
  broad_elev <- simsam::sim_covariates(
    rast_grid,
    n = 1,
    method = simsam::simulate_gaussian(
      psill = 2,
      model = "Exp",
      range = 80
    )
  )

  # 2. Local terrain roughness
  rough_elev <- simsam::sim_covariates(
    rast_grid,
    n = 1,
    method = simsam::simulate_gaussian(
      psill = 2,
      model = "Exp",
      range = 10
    )
  )

  # 3. Optional directional trend: lower southwest, higher northeast
  trend_vals <- 0.6 * x + 0.4 * y

  trend <- rast_grid
  terra::values(trend) <- trend_vals

  # 4. Optional mountain massif
  mountain_vals <- exp(
    -(((x - 150)^2) / (2 * 35^2) + ((y - 150)^2) / (2 * 35^2))
  )

  mountain <- rast_grid
  terra::values(mountain) <- mountain_vals

  # 5. Optional lowland / valley
  valley_vals <- exp(-((y - 35)^2) / (2 * 25^2))

  valley <- rast_grid
  terra::values(valley) <- valley_vals

  # Combine components
  elev <-
    1.2 *
    standardize_raster(broad_elev) +
    0.4 * standardize_raster(rough_elev) +
    0.8 * standardize_raster(trend) +
    1.0 * standardize_raster(mountain) -
    0.7 * standardize_raster(valley)

  # Smooth slightly
  elev <- terra::focal(
    elev,
    w = matrix(1, 3, 3),
    fun = mean,
    expand = TRUE,
    na.rm = TRUE
  )

  # Final standardized elevation layer
  elev <- standardize_raster(elev)

  names(elev) <- "elev"
  return(elev)
}


#' @keywords internal
#' @noRd
rescale_raster <- function(x, to = c(1, 100)) {
  mn <- terra::global(x, "min", na.rm = TRUE)[1, 1]
  mx <- terra::global(x, "max", na.rm = TRUE)[1, 1]

  if (mx == mn) {
    x[] <- mean(to)
    return(x)
  }

  to[1] + (x - mn) * diff(to) / (mx - mn)
}
