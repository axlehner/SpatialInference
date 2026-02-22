#' Linear Model with Spatial Autocorrelation Diagnostics
#'
#' Estimates a linear regression via [lfe::felm()] and augments the output
#' with Moran's I tests for spatial autocorrelation, optional correlograms
#' for range estimation, and Conley spatial HAC standard errors. The returned
#' object has class `"custom"` prepended, enabling display via
#' [modelsummary::modelsummary()] with custom `tidy` and `glance` methods.
#'
#' @param formula.chr Character string specifying the regression formula in
#'   `felm` syntax (e.g., `"y ~ x1 + x2 | fe1 + fe2 | 0 | lat + lon"`).
#' @param data.sf An `sf` data frame containing the variables referenced in
#'   `formula.chr`, with point geometries and columns `lat` and `lon`.
#' @param knn_number Integer. Number of nearest neighbours for the spatial
#'   weights matrix used in Moran's I tests. Default is `20`.
#' @param conley_cutoff Numeric. Spatial bandwidth (cutoff distance in km)
#'   for the Conley standard error. Default is `5`.
#' @param conley_kernel Character string specifying the kernel function.
#'   Default is `"bartlett"`. See [conley_SE()] for options.
#' @param correlograms Logical. If `TRUE`, estimates correlograms via
#'   [ncf::correlog()] and uses the extracted correlation range as an
#'   additional flexible Conley cutoff. Default is `FALSE`.
#' @param ... Additional arguments passed to [lfe::felm()] and [stats::lm()].
#'
#' @return An object of class `c("custom", "lm")` with additional components:
#'   \describe{
#'     \item{spatial_FE}{Character, the spatial fixed effect variable name.}
#'     \item{Moran_lmresid}{Moran's I test statistic on the OLS residuals,
#'       or `NA` if the test failed.}
#'     \item{Moran_response}{Moran's I test statistic on the response
#'       variable, or `NA` if the test failed.}
#'     \item{correlog.range_resid}{Estimated correlation range from the
#'       residual correlogram (km), or `NA` if `correlograms = FALSE`.}
#'     \item{correlog.range_response}{Estimated correlation range from the
#'       response correlogram (km), or `NA` if `correlograms = FALSE`.}
#'     \item{conley_SE}{Numeric vector of Conley spatial standard errors
#'       (with 0 for intercept and higher-order FE coefficients).}
#'     \item{conley_SE_flex}{Conley SEs using the correlogram-based cutoff,
#'       or `NA` if `correlograms = FALSE`.}
#'   }
#'
#' @references
#' Conley, T. G. (1999). GMM estimation with cross sectional dependence.
#' *Journal of Econometrics*, 92(1), 1--45.
#' \doi{10.1016/S0304-4076(98)00084-0}
#'
#' Moran, P. A. P. (1950). Notes on continuous stochastic phenomena.
#' *Biometrika*, 37(1/2), 17--23. \doi{10.2307/2332142}
#'
#' Bivand, R. S., Pebesma, E. and Gomez-Rubio, V. (2013).
#' *Applied Spatial Data Analysis with R*. 2nd ed. Springer.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(US_counties_centroids)
#' if (requireNamespace("lfe", quietly = TRUE) &&
#'     requireNamespace("spdep", quietly = TRUE) &&
#'     requireNamespace("stringr", quietly = TRUE) &&
#'     requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("sandwich", quietly = TRUE)) {
#'   out <- lm_sac("noise1 ~ noise2 | unit + year | 0 | lat + lon",
#'                  US_counties_centroids, conley_cutoff = 500)
#'   out$conley_SE
#' }
#' }
lm_sac <- function(formula.chr, data.sf, knn_number = 20,
                   conley_cutoff = 5, conley_kernel = "bartlett",
                   correlograms = FALSE, ...) {

  if (!requireNamespace("lfe", quietly = TRUE))
    stop("Package 'lfe' is required for lm_sac().", call. = FALSE)
  if (!requireNamespace("stringr", quietly = TRUE))
    stop("Package 'stringr' is required for lm_sac().", call. = FALSE)
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Package 'dplyr' is required for lm_sac().", call. = FALSE)
  if (!requireNamespace("spdep", quietly = TRUE))
    stop("Package 'spdep' is required for lm_sac().", call. = FALSE)
  if (!requireNamespace("sandwich", quietly = TRUE))
    stop("Package 'sandwich' is required for lm_sac().", call. = FALSE)

  # TRY RIGHT FROM THE START TO FIX THE SSIZE ISSUE BY SUBSETTING away THE DEP VAR NA's
  dep.var <- stringr::str_split(formula.chr, "\\ ~ ")[[1]][1]
  data.sf <- data.sf %>% dplyr::filter(!is.na(!!as.symbol(dep.var)))

  # 1) run the regression in lfe

  obj.felm <- lfe::felm(stats::as.formula(formula.chr), data = data.sf, keepCX = TRUE, keepModel = TRUE, ...)
  obj.felm

  # ... start prepping to run the regression in lm
  # split off the FE, in the non-FE case it is the "unit" one
  formula.FEs <- stringr::str_split(formula.chr, "\\| ")[[1]][2]
  the.spatial.FE <- stringr::str_split(formula.FEs, " \\+ year")[[1]][1]

  # now fix the formula to
  formula.chr1 <- stringr::str_replace_all(formula.chr, "\\|", "+")
  formula.chr1 <- stringr::str_split(formula.chr1, " \\+ year")[[1]][1]

  if (stringr::str_detect(formula.chr1, "unit") == TRUE) {
    formula.chr1 <- stringr::str_split(formula.chr1, " \\+ unit")[[1]][1]
  }

  # 2) run the regression in lm
  obj.lm <- stats::lm(stats::as.formula(formula.chr1), data = data.sf, ...)
  obj.lm$spatial_FE <- the.spatial.FE

  # merge them onto each other (using the lat lon from the felm as unique identifiers)
  points.merged <- dplyr::right_join(data.sf, obj.felm$model, by = c("lat", "lon"), suffix = c("", ".y"))

  # WEIGHTs
  points_knn_20 <- spdep::knn2nb(spdep::knearneigh(sf::st_coordinates(points.merged), k = knn_number))
  obj.lm$Moran_lmresid <- try(spdep::lm.morantest(obj.lm, spdep::nb2listw(points_knn_20))$statistic)
  if (inherits(obj.lm$Moran_lmresid, "try-error")) obj.lm$Moran_lmresid <- NA
  obj.lm$Moran_response <- try(spdep::moran.test(obj.lm$model[,1], spdep::nb2listw(points_knn_20))$statistic)
  if (inherits(obj.lm$Moran_response, "try-error")) obj.lm$Moran_response <- NA

  # CORELLOGRAM
  if (correlograms == TRUE) {
    if (!requireNamespace("ncf", quietly = TRUE))
      stop("Package 'ncf' is required when correlograms = TRUE.", call. = FALSE)
    correlog.fit <- ncf::correlog(x = points.merged$lon, y = points.merged$lat,
                                  z = as.numeric(obj.lm$residuals), increment = 1,
                                  resamp = 1, quiet = TRUE, latlon = TRUE, na.rm = TRUE)
    obj.lm$correlog.range_resid <- extract_corr_range(correlog.fit)
    correlog.fit <- ncf::correlog(x = points.merged$lon, y = points.merged$lat,
                                  z = obj.lm$model[,1], increment = 1,
                                  resamp = 1, quiet = TRUE, latlon = TRUE, na.rm = TRUE)
    obj.lm$correlog.range_response <- extract_corr_range(correlog.fit)
  } else {obj.lm$correlog.range_resid <- NA; obj.lm$correlog.range_response <- NA}
  # CONLEY
  cutoff <- conley_cutoff
  conley.1 <- conley_SE(reg = obj.felm,
                        unit = the.spatial.FE, time = "year",
                        kernel = conley_kernel, dist_fn = "Haversine",
                        lat = "lat", lon =  "lon", dist_cutoff = cutoff)
  # here we make the stunt to get a vector of the length of n(coef), including the intercept since lm() always displays it
  obj.lm$conley_SE <- c(0, lapply(conley.1, function(x) diag(sqrt(x)))$Spatial[[1]], rep(0, length(obj.lm$coefficients)-2))


  if (correlograms == TRUE) {
    cutoff <- obj.lm$correlog.range_resid
    conley.2 <- conley_SE(reg = obj.felm,
                          unit = the.spatial.FE, time = "year",
                          kernel = conley_kernel, dist_fn = "Haversine",
                          lat = "lat", lon =  "lon", dist_cutoff = cutoff)
    obj.lm$conley_SE_flex <- c(0, lapply(conley.2, function(x) diag(sqrt(x)))$Spatial[[1]], rep(0, length(obj.lm$coefficients)-2))
  } else {obj.lm$conley_SE_flex <- NA}

  # spit out
  class(obj.lm) <- c("custom", class(obj.lm))
  return(obj.lm)

}


# Tidy functions ---------------------------------------------------------------

#' Custom Tidy Method for lm_sac Output
#'
#' Extracts coefficient estimates, HC1 robust standard errors, and Conley
#' spatial standard errors from an `lm` object augmented by [lm_sac()].
#' Used by [modelsummary::modelsummary()] for table formatting.
#'
#' @param x An `lm` object (typically with class `"custom"` prepended by
#'   [lm_sac()]).
#' @param ... Additional arguments (currently unused).
#'
#' @return A `data.frame` with columns `term`, `estimate`, `std.error`
#'   (HC1 robust), `conf.high` (Conley SE), and `conley` (Conley SE).
#'
#' @keywords internal
tidy_custom.lm <- function(x, ...) {
  if (!requireNamespace("sandwich", quietly = TRUE))
    stop("Package 'sandwich' is required for tidy_custom.lm().", call. = FALSE)
  data.frame(
    term = names(x$coefficients),
    estimate = unname(x$coefficients),
    std.error = sandwich::vcovHC(x, type = "HC1") %>% diag() %>% sqrt(),
    conf.high = x$conley_SE,
    conley = x$conley_SE
  )
}

#' Custom Glance Method for lm_sac Output
#'
#' Extracts goodness-of-fit statistics including Moran's I statistics and
#' correlation ranges from an `lm` object augmented by [lm_sac()].
#' Used by [modelsummary::modelsummary()] for table formatting.
#'
#' @param x An `lm` object (typically with class `"custom"` prepended by
#'   [lm_sac()]).
#' @param ... Additional arguments (currently unused).
#'
#' @return A one-row `data.frame` with columns `Model`, `y_mean`, `y_SD`,
#'   `nobs`, `Moran_y`, `Moran_resid`, `Range_resid`, and `Range_y`.
#'
#' @keywords internal
glance_custom.lm <- function(x, ...) {
  data.frame(
    "Model" = "Custom",
    "y_mean" = mean(x$model[,1]),
    "y_SD"   = sd(x$model[,1]),
    "nobs" = stats::nobs(x),
    "Moran_y" =     x$Moran_response,
    "Moran_resid" = x$Moran_lmresid,
    "Range_resid" = x$correlog.range_resid,
    "Range_y" =  x$correlog.range_response
  )
}

#' Goodness-of-Fit Mapping for modelsummary
#'
#' A list specifying how goodness-of-fit statistics should be displayed
#' in [modelsummary::modelsummary()] output.
#'
#' @keywords internal
gm.param <- list(
  list("raw" = "nobs", "clean" = "Observations", "fmt" = 0),
  list("raw" = "Moran_y", "clean" = "Moran's I [y]", "fmt" = 3),
  list("raw" = "Moran_resid", "clean" = "Moran's I [resid]", "fmt" = 3),
  list("raw" = "y_mean", "clean" = "Mean [y]", "fmt" = 3),
  list("raw" = "y_SD", "clean" = "Std.Dev. [y]", "fmt" = 3)
)
