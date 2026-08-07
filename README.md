
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PredictionMatching

<!-- badges: start -->

<!-- badges: end -->

**PredictionMatching** provides weighted evaluation for cross-validation
using raking when training data are sampled non-representatively.

> *Experimental*: This is an implementation of weighted evaluation based
> on raking (Brenning and Suesse (2026)).

## Installation

You can install the development version of PredictionMatching from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("JanLinnenbrink/PredictionMatching")
```

## Example

### Data and setup

We use an example where the sampling effort is biased along a predictor
included in the model (e.g., towards low elevation) to demonstrate the
weighting functions implemented in this package:

``` r
library(PredictionMatching)
library(terra)
library(sf)
library(simsam)
library(caret)
library(cowplot)

set.seed(100)

calc_rmse <- function(pred, obs) {
  sqrt(mean((pred - obs)^2, na.rm = TRUE))
}

raster_stack <- terra::rast(system.file(
  "extdata",
  "rasters_example.tif",
  package = "PredictionMatching"
))
data(train_data)

predictor_names <- setdiff(names(raster_stack), "outcome")
predictor_stack <- raster_stack[[predictor_names]]

terra::plot(predictor_stack[["elev"]])
terra::plot(vect(train_data), add = T)
```

<img src="man/figures/README-example-1.png" alt="" width="100%" />

### kNNDM cross-validation

We train a random forest model and assess its performance using kNNDM
cross-validation (Linnenbrink et al. (2024)).

``` r
# Split the data into folds using kNNDM
knndm_folds <- CAST::knndm(
  tpoints = train_data,
  modeldomain = predictor_stack[["elev"]],
  k = 5
)
#> 1000 prediction points are sampled from the modeldomain

# Train a random forest model
ctrl_knndm <- trainControl(
  method = "cv",
  index = knndm_folds$indx_train,
  indexOut = knndm_folds$indx_test,
  savePredictions = "final",
  verboseIter = FALSE
)

rf_knndm <- train(
  x = st_drop_geometry(train_data)[, predictor_names],
  y = st_drop_geometry(train_data)[["outcome"]],
  method = "ranger",
  trControl = ctrl_knndm,
  metric = "RMSE",
  num.trees = 300,
  importance = "impurity"
)

knndm_rmse <- calc_rmse(rf_knndm$pred$pred, rf_knndm$pred$obs)

plot(knndm_folds)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" alt="" width="100%" />

### True RMSE and RMSE from geographic kNNDM

We can also calculate the true RMSE in this setting:

``` r
rf_pred <- predict(predictor_stack, model = rf_knndm, na.rm = TRUE)
rf_pred_vals <- terra::values(rf_pred, mat = FALSE)
outcome_vals <- terra::values(raster_stack$outcome, mat = FALSE)
true_rmse <- calc_rmse(rf_pred_vals, outcome_vals)
true_rmse
#> [1] 7.660114
knndm_rmse
#> [1] 6.646014
```

As can be seen by comparing the RMSE estimate obtained by kNNDM to the
true RMSE, using kNNDM alone is sometimes not sufficient with biased
sampling. This happens because kNNDM was applied in geographic space and
does not necessarily correct the sampling bias along elevation. To
correct for this, we can use raking to obtain weights that correct for
the sampling bias along elevation.

This process requires three steps: (1) calculate weights via
`tw_calculate_weights()`, (2) compute pointwise cross-validation errors
via `tw_pointwise_error()`, and (3) compute the weighted RMSE via
`tw_weighted_error_stats()`. First, we calculate weights via raking
using `tw_calculate_weights()`. This function takes the training data
and the predictor stack as input and returns a vector of weights that
correct for the sampling bias along the specified predictors. After
obtaining weights via `tw_calculate_weights()`, we compute pointwise
cross-validation errors with `tw_pointwise_error()`. This preserves the
original row IDs so that errors can be aligned correctly with weights
even when caret reorders resampling rows. Finally,
`tw_weighted_error_stats()` computes the weighted RMSE.

``` r
# Calculate weights from raking
w <- tw_calculate_weights(
  tpoints = st_drop_geometry(train_data)[, predictor_names],
  modeldomain = predictor_stack
)
#> 1000 prediction points are sampled from the modeldomain
#> predictor values are extracted for prediction points

# Standardize the pointwise errors obtained by kNNDM CV and add an ID column
pe <- tw_pointwise_error(
  obs = rf_knndm$pred$obs,
  pred = rf_knndm$pred$pred,
  id = rf_knndm$pred$rowIndex
)

# Weigh the pointwise errors using the weights obtained from raking
weighted_rmse <- tw_weighted_error_stats(w, pe)[["rmse"]]

plot_grid(
  plotlist = plot(w, pointwise_error = pe),
  nrow = 1,
  align = "v",
  labels = c("A", "B")
)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" alt="" width="100%" />

From panel B of the above plot we can already see that points with a
higher CV-error receive a higher weight, resulting in a higher RMSE
which gets closer to the true RMSE:

|               |     RMSE |
|:--------------|---------:|
| True RMSE     | 7.660114 |
| kNNDM RMSE    | 6.646014 |
| Weighted RMSE | 7.491875 |

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Brenning2026" class="csl-entry">

Brenning, Alexander, and Thomas Suesse. 2026. *Aligning Validation with
Deployment: Target-Weighted Cross-Validation for Spatial Prediction*.
arXiv. <https://doi.org/10.48550/ARXIV.2603.29981>.

</div>

<div id="ref-Linnenbrink2024" class="csl-entry">

Linnenbrink, Jan, Carles Milà, Marvin Ludwig, and Hanna Meyer. 2024.
“<span class="nocase">kNNDM</span> CV: *K* -Fold Nearest-Neighbour
Distance Matching Cross-Validation for Map Accuracy Estimation.”
*Geoscientific Model Development* 17 (15): 5897–912.
<https://doi.org/10.5194/gmd-17-5897-2024>.

</div>

</div>
