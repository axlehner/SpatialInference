#' Gravity Centroid (Centre of Gravity)
#'
#' Computes the (optionally weighted) mean centre of an `sf` data frame,
#' returning a single `sfc_POINT`. Without weights, this is the simple
#' geographic centroid; with weights, it is the weighted mean centre
#' (centre of gravity) as used in spatial statistics.
#'
#' @param df.sf An `sf` data frame with polygon or point geometries.
#' @param weight Numeric vector of weights with length equal to `nrow(df.sf)`,
#'   or `NA` (default) for an unweighted mean centre.
#'
#' @return An `sfc_POINT` object (CRS 4326 / WGS84) representing the
#'   (weighted) mean centre.
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
#' gravity_centroid(US_counties_centroids)
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
