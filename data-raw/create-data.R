## Code to prepare example datasets
set.seed(100)

# Simulate predictor and response rasters
rast_grid <- terra::rast(
  xmin = 0,
  xmax = 200,
  ymin = 0,
  ymax = 200,
  ncols = 200,
  nrows = 200
)

# Predictors
predictors <- simsam::sim_covariates(
  rast_grid,
  n = 4,
  method = simsam::simulate_gaussian(
    nugget = 0,
    beta = 50,
    psill = 100,
    model = "Exp",
    range = 50
  )
)
names(predictors) <- c("temp", "ph", "slope", "elev")

# Response
outcome_signal <- simsam::blend_rasters(
  predictors,
  ~ 40 *
    exp(-((ph - 50)^2) / (2 * 15^2)) + # non-linear (bell-shaped) effect of ph
    0.4 * temp - # linear, positive
    0.5 * elev # linear, negative (drives sampling bias)
)

# Add iid noise
outcome <- outcome_signal
terra::values(outcome) <- terra::values(outcome_signal) +
  rnorm(terra::ncell(outcome_signal), sd = 5)
names(outcome) <- "outcome"

raster_stack <- c(predictors, outcome)
terra::crs(raster_stack) <- "EPSG:3857"

# --- Save raster ----------------------------------------------------------
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
terra::writeRaster(
  raster_stack,
  "inst/extdata/rasters_example.tif",
  datatype = "FLT4S", # keep continuous predictors intact
  overwrite = TRUE
)

# --- biased training sample (prefers low elevation) -----------------------
train_points <- simsam::sam_field(
  raster_stack$elev * -1,
  100,
  method = simsam::sample_preferential(strength = 4)
)

# training data
train_data <- terra::extract(
  raster_stack,
  train_points,
  ID = FALSE,
  xy = TRUE
) |>
  sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs(train_points))

usethis::use_data(
  train_data,
  overwrite = TRUE,
  compress = "xz"
)
