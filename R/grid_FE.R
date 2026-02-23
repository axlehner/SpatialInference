#' Create Spatial Grid Fixed Effects
#'
#' Overlays a regular rectangular grid on an `sf` data frame, intersects each
#' observation with the grid, and returns a factor variable identifying the
#' grid cell to which each observation belongs. This is useful for constructing
#' spatial fixed effects in regression models: by including `grid_FE()` as a
#' factor variable, the regression absorbs location-specific variation at the
#' resolution of the chosen grid. Finer grids absorb more spatial variation
#' but consume more degrees of freedom.
#'
#' The grid is constructed via [sf::st_make_grid()] and observations are
#' assigned to cells via [sf::st_intersection()]. Observations that fall
#' outside the grid (e.g., in coastal regions where cells do not cover the
#' full bounding box) are dropped during intersection.
#'
#' @param df.sf An `sf` data frame with geometries (points or polygons).
#' @param size Numeric. When `distance = FALSE` (default), an integer
#'   specifying the number of grid cells along each axis (e.g., `size = 10`
#'   creates a 10 x 10 grid). When `distance = TRUE`, a numeric value
#'   giving the cell side length in CRS units (e.g., metres for projected
#'   data). Passed to [sf::st_make_grid()] as `n` or `cellsize`, respectively.
#' @param distance Logical. If `FALSE` (default), `size` is the number of
#'   cells per axis. If `TRUE`, `size` is the cell dimension in CRS units.
#'
#' @return A factor vector with one element per observation that survived the
#'   intersection (see Details), where each level is a grid cell ID.
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
#'
#' # Create a 5 x 5 grid of spatial fixed effects
#' grid_ids <- grid_FE(US_counties_centroids, size = 5)
#' table(grid_ids)
#'
#' # Finer grid (10 x 10)
#' grid_ids_fine <- grid_FE(US_counties_centroids, size = 10)
#' length(levels(grid_ids_fine))
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
