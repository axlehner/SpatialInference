# SpatialInference 0.1.0

* Initial CRAN release.

## Conley spatial HAC standard errors

* `conley_SE()` computes Conley (1999) spatial HAC variance-covariance matrices
  for `lfe::felm()` models, with support for cross-sectional spatial correlation,
  serial (temporal) correlation, and the combined spatial HAC estimator.
* Six kernel functions: Bartlett, Epanechnikov, Gaussian, Parzen, Biweight, and
  Uniform.
* Haversine great-circle distances (default) and a 111 km/degree approximation.
* Balanced-panel optimisation pre-computes the distance matrix once.
* `compute_conley_lfe()` convenience wrapper for quick single-coefficient
  extraction.
* `lm_sac()` all-in-one workflow: regression, Moran's I tests, and Conley
  standard errors, with `modelsummary` integration via custom `tidy` and `glance`
  methods.

## Bandwidth selection

* `covgm_range()` estimates the spatial correlation range from the empirical
  covariogram of regression residuals (Lehner 2026).
* `extract_corr_range()` extracts the zero-crossing distance from a covariogram
  (`gstat::variogram()`) or correlogram (`ncf::correlog()`).
* `inverseu_plot_conleyrange()` diagnostic plot showing how the Conley SE varies
  with the bandwidth, revealing the inverse-U relationship (Lehner 2026).

## Spatial utilities

* `DistMat()` kernel-weighted spatial distance matrix (C++).
* `coords_as_columns()` extracts `sf` point coordinates into tibble columns.
* `gravity_centroid()` computes the (optionally weighted) geographic centroid of
  an `sf` object.
* `grid_FE()` assigns observations to spatial grid cells for use as fixed
  effects.

## Performance

* Distance matrix computation, kernel weighting, and variance component
  accumulation (`XeeXhC`, `Bal_XeeXhC`, `XeeXhC_Lg`, `TimeDist`) are
  implemented in C++ via Rcpp and RcppArmadillo.
* Memory-efficient large-sample variant (`XeeXhC_Lg`) avoids constructing the
  full n x n distance matrix when n > 50,000.
