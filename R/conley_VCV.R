#' @importFrom data.table :=
#' @export
data.table::`:=`

#' Conley Spatial Variance-Covariances using lfe
#'
#' This is the function containing the underlying C++ code by Darin Christensen.
#' Regression model input has to be from lfe::felm().
#'
#' @param reg
#' @param unit
#' @param time
#' @param lat
#' @param lon
#' @param kernel
#' @param dist_fn
#' @param dist_cutoff
#' @param lag_cutoff
#' @param lat_scale
#' @param verbose
#' @param cores
#' @param balanced_pnl
#'
#' @return
#' @export
#'
#' @examples
conley_SE <- function(reg,
                      unit, time, lat, lon,
                      kernel = "bartlett", dist_fn = "Haversine",
                      dist_cutoff = 500, lag_cutoff = 5,
                      lat_scale = 111, verbose = FALSE, cores = 1, balanced_pnl = FALSE) {

  require(lfe)

  Fac2Num <- function(x) {as.numeric(as.character(x))}
  if(cores > 1) {invisible(library(parallel))}

  if(class(reg) == "felm") {
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
    message("Model class not recognized.")
    break
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
    XeeXhs <- mclapply(timeUnique, function(t) iterateObs(dt, Xvars, d, k,
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
    XeeXhs <- mclapply(panelUnique,function(t) iterateObs(dt, Xvars, d, k,
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


#' Iterate observations
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


#' Compute Conley Standard Error
#'
#' Convenience function which nests `conley_SE` and prints only the standard error.
#' The same parameters can be used here.
#'
#' @param lfeobj
#' @param cutoff
#' @param kernel_choice
#' @param ...
#'
#' @return
#' @export
#'
#' @examples


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
