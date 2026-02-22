#' @importFrom data.table :=
#' @export
data.table::`:=`

#' Conley Spatial HAC Variance-Covariance Estimation
#'
#' Computes Conley (1999) spatial HAC (Heteroskedasticity and Autocorrelation
#' Consistent) variance-covariance matrices for regression models estimated
#' with [lfe::felm()]. Supports cross-sectional spatial correlation,
#' serial (temporal) correlation, and the combined spatial HAC estimator.
#' Multiple kernel functions and distance metrics are available.
#'
#' @param reg A regression object of class `"felm"` from [lfe::felm()].
#'   Must be estimated with `keepCX = TRUE` and with latitude/longitude
#'   variables passed as cluster variables.
#' @param unit Character string naming the cross-sectional unit identifier
#'   variable in the `felm` cluster variables.
#' @param time Character string naming the time period identifier variable
#'   in the `felm` cluster variables.
#' @param lat Character string naming the latitude variable in the `felm`
#'   cluster variables.
#' @param lon Character string naming the longitude variable in the `felm`
#'   cluster variables.
#' @param kernel Character string specifying the kernel function for spatial
#'   weighting. One of `"bartlett"` (default), `"epanechnikov"`, `"gaussian"`,
#'   `"parzen"`, `"biweight"`, or `"uniform"`.
#' @param dist_fn Character string specifying the distance function.
#'   `"Haversine"` (default) for great-circle distance in km, or `"SH"` for
#'   the 111 km/degree approximation.
#' @param dist_cutoff Numeric. The spatial bandwidth (cutoff distance in km)
#'   beyond which observations receive zero weight. Default is `500`.
#' @param lag_cutoff Numeric. The temporal bandwidth (number of time periods)
#'   for serial correlation correction. Default is `5`.
#' @param lat_scale Numeric. Scaling factor for latitude (km per degree).
#'   Default is `111`.
#' @param verbose Logical. If `TRUE`, prints progress messages during
#'   computation. Default is `FALSE`.
#' @param cores Integer. Number of CPU cores for parallel computation via
#'   [parallel::mclapply()]. Default is `1` (no parallelism).
#' @param balanced_pnl Logical. If `TRUE`, assumes a balanced panel and
#'   pre-computes the distance matrix once for efficiency. Default is `FALSE`.
#'
#' @return A named list of three variance-covariance matrices, each of
#'   dimension k x k where k is the number of regressors:
#'   \describe{
#'     \item{OLS}{The standard OLS variance-covariance matrix from the felm
#'       object.}
#'     \item{Spatial}{The spatial-only Conley VCV, correcting for
#'       cross-sectional spatial correlation.}
#'     \item{Spatial_HAC}{The full spatial HAC VCV, correcting for both
#'       spatial and serial correlation.}
#'   }
#'
#' @references
#' Conley, T. G. (1999). GMM estimation with cross sectional dependence.
#' *Journal of Econometrics*, 92(1), 1--45.
#' \doi{10.1016/S0304-4076(98)00084-0}
#'
#' Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite,
#' heteroskedasticity and autocorrelation consistent covariance matrix.
#' *Econometrica*, 55(3), 703--708. \doi{10.2307/1913610}
#'
#' Christensen, D., Hartman, A. C. and Samii, C. (2021). Legibility and
#' external investment: An institutional natural experiment in Liberia.
#' *International Organization*, 75(4), 1087--1108.
#' \doi{10.1017/S0020818321000187}
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(US_counties_centroids)
#' if (requireNamespace("lfe", quietly = TRUE)) {
#'   reg <- lfe::felm(noise1 ~ noise2 | unit + year | 0 | lat + lon,
#'                     data = US_counties_centroids, keepCX = TRUE)
#'   vcvs <- conley_SE(reg, unit = "unit", time = "year",
#'                     lat = "lat", lon = "lon",
#'                     kernel = "bartlett", dist_cutoff = 500)
#'   # Spatial standard errors:
#'   sqrt(diag(vcvs$Spatial))
#' }
#' }
conley_SE <- function(reg,
                      unit, time, lat, lon,
                      kernel = "bartlett", dist_fn = "Haversine",
                      dist_cutoff = 500, lag_cutoff = 5,
                      lat_scale = 111, verbose = FALSE, cores = 1, balanced_pnl = FALSE) {

  if (!requireNamespace("lfe", quietly = TRUE)) {
    stop("Package 'lfe' is required for conley_SE(). Please install it.",
         call. = FALSE)
  }

  Fac2Num <- function(x) {as.numeric(as.character(x))}

  if (cores > 1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("Package 'parallel' is required for multi-core computation.",
           call. = FALSE)
    }
  }

  if (inherits(reg, "felm")) {
    Xvars <- rownames(reg$coefficients)
    dt = data.table::data.table(reg$cY, reg$cX,
                                fe1 = Fac2Num(reg$fe[[1]]),
                                fe2 = Fac2Num(reg$fe[[2]]),
                                coord1 = Fac2Num(reg$clustervar[[1]]),
                                coord2 = Fac2Num(reg$clustervar[[2]]))
    data.table::setnames(dt,
                         c("fe1", "fe2", "coord1", "coord2"),
                         c(names(reg$fe), names(reg$clustervar)))
    dt = dt[, e := as.numeric(reg$residuals)]

  } else {
    stop("Model class not recognized. Please provide a 'felm' object from lfe::felm().",
         call. = FALSE)
  }

  n <- nrow(dt)
  k <- length(Xvars)

  # Renaming variables:
  orig_names <- c(unit, time, lat, lon)
  new_names <- c("unit", "time", "lat", "lon")
  data.table::setnames(dt, orig_names, new_names)

  # Empty Matrix:
  XeeX <- matrix(nrow = k, ncol = k, 0)

  #================================================================
  # Correct for spatial correlation:
  timeUnique <- unique(dt[, time])
  Ntime <- length(timeUnique)
  data.table::setkey(dt, time)

  if(verbose){message("Starting to loop over time periods...")}

  if(balanced_pnl){
    sub_dt <- dt[time == timeUnique[1]]
    lat <- sub_dt[, lat]; lon <- sub_dt[, lon]; rm(sub_dt)

    if(balanced_pnl & verbose){message("Computing Distance Matrix...")}

    d <- DistMat(cbind(lat, lon), cutoff = dist_cutoff, kernel, dist_fn)
    rm(list = c("lat", "lon"))
  }

  if(cores == 1) {
    XeeXhs <- lapply(timeUnique, function(t) iterateObs(dt, Xvars, d, k,
                                                        sub_index = t,
                                                        type = "spatial", cutoff = dist_cutoff,
                                                        balanced_pnl = balanced_pnl,
                                                        verbose = verbose,
                                                        kernel = kernel,
                                                        dist_fn = dist_fn))
  } else {
    XeeXhs <- parallel::mclapply(timeUnique, function(t) iterateObs(dt, Xvars, d, k,
                                                          sub_index = t,
                                                          type = "spatial", cutoff = dist_cutoff,
                                                          balanced_pnl = balanced_pnl,
                                                          verbose = verbose,
                                                          kernel = kernel,
                                                          dist_fn = dist_fn), mc.cores = cores)
  }

  if(balanced_pnl){rm(d)}

  # First Reduce:
  XeeX <- Reduce("+",  XeeXhs)

  # Generate VCE for only cross-sectional spatial correlation:
  X <- as.matrix(dt[, eval(Xvars), with = FALSE])
  invXX <- solve(t(X) %*% X) * n

  V_spatial <- invXX %*% (XeeX / n) %*% invXX / n

  V_spatial <- (V_spatial + t(V_spatial)) / 2

  if(verbose) {message("Computed Spatial VCOV.")}

  #================================================================
  # Correct for serial correlation:
  panelUnique <- unique(dt[, unit])
  Npanel <- length(panelUnique)
  data.table::setkey(dt, unit)

  if(verbose){message("Starting to loop over units...")}

  if(cores == 1) {
    XeeXhs <- lapply(panelUnique, function(t) iterateObs(dt, Xvars, d, k,
                                                         sub_index = t,
                                                         type = "serial", cutoff = lag_cutoff,
                                                         balanced_pnl = balanced_pnl,
                                                         verbose = verbose,
                                                         kernel = kernel,
                                                         dist_fn = dist_fn))
  } else {
    XeeXhs <- parallel::mclapply(panelUnique,function(t) iterateObs(dt, Xvars, d, k,
                                                          sub_index = t,
                                                          type = "serial", cutoff = lag_cutoff,
                                                          balanced_pnl = balanced_pnl,
                                                          verbose = verbose,
                                                          kernel = kernel,
                                                          dist_fn = dist_fn), mc.cores = cores)
  }

  XeeX_serial <- Reduce("+",  XeeXhs)

  XeeX <- XeeX + XeeX_serial

  V_spatial_HAC <- invXX %*% (XeeX / n) %*% invXX / n
  V_spatial_HAC <- (V_spatial_HAC + t(V_spatial_HAC)) / 2

  return_list <- list(
    "OLS" = reg$vcv,
    "Spatial" = V_spatial,
    "Spatial_HAC" = V_spatial_HAC)
  return(return_list)
}


#' Iterate over observations for spatial or serial correlation
#'
#' Internal helper function called by [conley_SE()] to compute the X'ee'X
#' component for a single time period (spatial) or unit (serial).
#'
#' @param dt A `data.table` with the regression data.
#' @param Xvars Character vector of regressor names.
#' @param d Pre-computed distance matrix (used when `balanced_pnl = TRUE`).
#' @param k Integer, number of regressors.
#' @param sub_index The time period or unit index to subset on.
#' @param type Character, either `"spatial"` or `"serial"`.
#' @param cutoff Numeric, the distance or lag cutoff.
#' @param balanced_pnl Logical, whether the panel is balanced.
#' @param verbose Logical, whether to print progress.
#' @param kernel Character, the kernel function name.
#' @param dist_fn Character, the distance function name.
#'
#' @return A k x k numeric matrix representing the X'ee'X contribution.
#'
#' @keywords internal
iterateObs <- function(dt, Xvars, d, k,
                       sub_index, type, cutoff, balanced_pnl, verbose, kernel, dist_fn) {
  if(type == "spatial" & balanced_pnl) {

    sub_dt <- dt[time == sub_index]
    n1 <- nrow(sub_dt)
    if(n1 > 1000 & verbose){message(paste("Starting on sub index:", sub_index))}

    X <- as.matrix(sub_dt[, eval(Xvars), with = FALSE])
    e <- sub_dt[, e]

    XeeXhs <- Bal_XeeXhC(d, X, e, n1, k)

  } else if(type == "spatial" & !balanced_pnl) {

    sub_dt <- dt[time == sub_index]
    n1 <- nrow(sub_dt)
    if(n1 > 1000 & verbose){message(paste("Starting on sub index:", sub_index))}

    X <- as.matrix(sub_dt[, eval(Xvars), with = FALSE])
    e <- sub_dt[, e]
    lat <- sub_dt[, lat]; lon <- sub_dt[, lon]

    # If n1 >= 50k obs, then avoiding construction of distance matrix.
    # This requires more operations, but is less memory intensive.
    if(n1 < 5 * 10^4) {
      XeeXhs <- XeeXhC(cbind(lat, lon), cutoff, X, e, n1, k,
                       kernel, dist_fn)
    } else {
      XeeXhs <- XeeXhC_Lg(cbind(lat, lon), cutoff, X, e, n1, k,
                          kernel, dist_fn)
    }

  } else if(type == "serial") {
    sub_dt <- dt[unit == sub_index]
    n1 <- nrow(sub_dt)
    if(n1 > 1000 & verbose){message(paste("Starting on sub index:", sub_index))}

    X <- as.matrix(sub_dt[, eval(Xvars), with = FALSE] )
    e <- sub_dt[, e]
    times <- sub_dt[, time]

    XeeXhs <- TimeDist(times, cutoff, X, e, n1, k)
  }

  XeeXhs
}


#' Compute Conley Standard Error for a Single Coefficient
#'
#' Convenience wrapper around [conley_SE()] that returns only the spatial
#' Conley standard error for the first regressor. Useful for quick extraction
#' of Conley SEs in scripting contexts.
#'
#' @param lfeobj A regression object of class `"felm"` from [lfe::felm()].
#'   Must be estimated with `keepCX = TRUE` and cluster variables
#'   `lat + lon`.
#' @param cutoff Numeric. The spatial bandwidth (cutoff distance in km).
#' @param kernel_choice Character string specifying the kernel function.
#'   Default is `"bartlett"`. See [conley_SE()] for available kernels.
#' @param ... Additional arguments passed to [conley_SE()].
#'
#' @return A single numeric value: the spatial Conley standard error for the
#'   first (or only) regressor.
#'
#' @references
#' Conley, T. G. (1999). GMM estimation with cross sectional dependence.
#' *Journal of Econometrics*, 92(1), 1--45.
#' \doi{10.1016/S0304-4076(98)00084-0}
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(US_counties_centroids)
#' if (requireNamespace("lfe", quietly = TRUE)) {
#'   reg <- lfe::felm(noise1 ~ noise2 | unit + year | 0 | lat + lon,
#'                     data = US_counties_centroids, keepCX = TRUE)
#'   compute_conley_lfe(reg, cutoff = 500)
#' }
#' }
compute_conley_lfe <- function(lfeobj, cutoff, kernel_choice = "bartlett", ...) {
  regfe_conley <- conley_SE(reg = lfeobj, unit = "unit", time = "year",
                            lat = "lat", lon = "lon",
                            kernel = kernel_choice, dist_fn = "Haversine",
                            dist_cutoff = cutoff, ...)
  conley <- sapply(regfe_conley, function(x) diag(sqrt(x))) |> round(5)
  if (length(conley) == 3) { # has length 3 if only one regressand was used
    return(conley[[2]])
  } else { # if there was more than one regressand, we are interested in the conley SE of the first one
    return(conley[1,2])
  }

}
