#' Weighted Mean Centre (Centre of Gravity)
#'
#' Computes the mean centre of an `sf` data frame, optionally weighted by
#' an attribute variable. The function first extracts polygon or point
#' centroids using [sf::st_centroid()], then calculates the (weighted)
#' arithmetic mean of the x and y coordinates. This is the two-dimensional
#' analogue of a weighted mean and corresponds to the "centre of gravity"
#' or "mean centre" in spatial statistics. A common use case is computing
#' the population-weighted centroid of a set of administrative units.
#'
#' @param df.sf An `sf` data frame with polygon or point geometries.
#'   Polygon geometries are reduced to their centroids before computation.
#' @param weight Numeric vector of weights with length equal to `nrow(df.sf)`.
#'   If `NA` (the default), all observations receive equal weight and the
#'   result is the simple (unweighted) geographic mean centre. Weights do
#'   not need to sum to one; the function normalises internally.
#'
#' @return An `sfc_POINT` object (CRS 4326 / WGS84) representing the
#'   (weighted) mean centre as a single point.
#'
#' @references
#' Arlinghaus, S. L. (1994). *Practical Handbook of Spatial Statistics*.
#' CRC Press.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(US_counties_centroids)
#'
#' # Unweighted mean centre (geographic centroid of all county centroids)
#' gravity_centroid(US_counties_centroids)
#'
#' # Weighted mean centre (shifted toward areas with higher noise1 values)
#' gravity_centroid(US_counties_centroids,
#'                  weight = abs(US_counties_centroids$noise1) + 1)
#' }
gravity_centroid <- function(df.sf, weight = NA) {

  centroids <- coords_as_columns(sf::st_centroid(df.sf))
  if (length(weight) == 1) {
    gravity.sf <- sf::st_sfc(sf::st_point(c(mean(centroids$x, na.rm = TRUE),
                                    mean(centroids$y, na.rm = TRUE))), crs = 4326)
    gravity.sf}
  else {
    gravity.sf <- sf::st_sfc(sf::st_point(c(sum(centroids$x * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE),
                                    sum(centroids$y * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE))), crs = 4326)

    gravity.sf
  }
}
