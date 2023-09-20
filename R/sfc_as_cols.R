


#' Coordinates as columns in dataframe
#'
#' code comes from the internet
#'
#' @param x
#' @param names
#'
#' @return
#' @export
#'
#' @examples
#'
#'
sfc_as_cols <- function(x, names = c("x","y")) {
  stopifnot(inherits(x,"sf") | inherits(x,"sfc") && inherits(sf::st_geometry(x),"sfc_POINT"))
  ret <- do.call(rbind,sf::st_geometry(x))
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  ret <- setNames(ret,names)
  ret
}


# NEWER VERSION ?
# sfc_as_cols <- function(x, geometry, names = c("x","y")) {
#   if (missing(geometry)) {
#     geometry <- sf::st_geometry(x)
#   } else {
#     geometry <- rlang::eval_tidy(enquo(geometry), x)
#   }
#   stopifnot(inherits(x,"sf") && inherits(geometry,"sfc_POINT"))
#   ret <- sf::st_coordinates(geometry)
#   ret <- tibble::as_tibble(ret)
#   stopifnot(length(names) == ncol(ret))
#   x <- x[ , !names(x) %in% names]
#   ret <- setNames(ret,names)
#   dplyr::bind_cols(x,ret)
# }

