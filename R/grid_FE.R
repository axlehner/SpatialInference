#' Create Grid-Based Fixed Effects
#'
#' Overlays a regular grid on an `sf` data frame and returns a factor
#' variable assigning each observation to a grid cell. Useful for creating
#' spatial fixed effects in regression models.
#'
#' @param df.sf An `sf` data frame with geometries.
#' @param size Numeric or integer vector. If `distance = FALSE` (default),
#'   specifies the number of grid cells along each axis (passed as `n` to
#'   [sf::st_make_grid()]). If `distance = TRUE`, specifies the cell size
#'   in CRS units.
#' @param distance Logical. If `FALSE` (default), `size` is interpreted as
#'   the number of cells. If `TRUE`, `size` is interpreted as cell dimensions
#'   in the CRS units. Default is `FALSE`.
#'
#' @return A factor vector of length `nrow(df.sf)` containing the grid cell
#'   ID for each observation after spatial intersection.
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
#' grid_ids <- grid_FE(US_counties_centroids, size = 5)
#' table(grid_ids)
#' }
grid_FE <- function(df.sf, size, distance = FALSE) {
  if (distance == TRUE) {
    grid <- sf::st_sf(geom = sf::st_make_grid(df.sf, size))
  } else {
    grid <- sf::st_sf(geom = sf::st_make_grid(df.sf, n = size))
  }
  grid$grid_id <- 1:nrow(grid)
  df.sf <- sf::st_intersection(df.sf, grid)
  df.sf$grid_id <- as.factor(df.sf$grid_id)
  return(df.sf$grid_id)
}
