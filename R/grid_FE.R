


#' Title
#'
#' @param df.sf
#' @param size
#' @param distance
#'
#' @return
#' @export
#'
#' @examples
grid_FE <- function(df.sf, size, distance = F) {
  # when size given this should be something like a division with modulo 0 on the length of the bbox sides, otherwise overlaps
  # NOW WE HAVE ONLY NUMBERS, so the size gives a sizeXsize grid, this way the whole area is covered
  if (distance == T) grid <- st_sf(geom = st_make_grid(df.sf, size)) else grid <- st_sf(geom = st_make_grid(df.sf, n = size))
  grid$grid_id <- 1:nrow(grid)
  df.sf <- st_intersection(df.sf, grid)
  df.sf$grid_id <- as.factor(df.sf$grid_id)
  # ORDERING WONT FIT
  return(df.sf$grid_id)
}
