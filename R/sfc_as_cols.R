#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' Extract Coordinates as Columns from an sf Object
#'
#' Extracts point coordinates from an `sf` or `sfc` object and returns
#' them as a tibble. This is a lightweight alternative to
#' [sf::st_coordinates()] that returns a tibble directly.
#'
#' @param x An `sf` or `sfc` object with `sfc_POINT` geometry.
#' @param names Character vector of length 2 specifying the column names
#'   for the x and y coordinates. Default is `c("x", "y")`.
#'
#' @return A [tibble::tibble()] with two columns named according to the
#'   `names` argument, containing the x and y coordinates.
#'
#' @references
#' Pebesma, E. (2018). Simple features for R: Standardized support for
#' spatial vector data. *The R Journal*, 10(1), 439--446.
#' \doi{10.32614/RJ-2018-009}
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(US_counties_centroids)
#' coords <- coords_as_columns(US_counties_centroids)
#' head(coords)
#' }
coords_as_columns <- function(x, names = c("x","y")) {
  stopifnot(inherits(x,"sf") | inherits(x,"sfc") && inherits(sf::st_geometry(x),"sfc_POINT"))
  ret <- do.call(rbind,sf::st_geometry(x))
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  ret <- setNames(ret,names)
  ret
}
