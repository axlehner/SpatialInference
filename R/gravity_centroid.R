


#' Center of Gravity
#'
#' @param df.sf
#' @param weight
#'
#' @return
#' @export
#' @references Arlinghaus
#'
#' @examples
gravity_centroid <- function(df.sf, weight = NA) {

  centroids <- sfc_as_cols(sf::st_centroid(df.sf))
  if (length(weight) == 1) { # case when we just want the average (conditional not is.na bc if the first element of the weight is NA it is messed up)
    gravity.sf <- sf::st_sfc(sf::st_point(c(mean(centroids$x, na.rm = TRUE),
                                    mean(centroids$y, na.rm = TRUE))), crs = 4326)
    gravity.sf}
  else { # here we do with weights according to http://resources.arcgis.com/en/help/main/10.1/index.html#/How_Mean_Center_works/005p0000001s000000/
    gravity.sf <- sf::st_sfc(sf::st_point(c(sum(centroids$x * weight, na.rm = T) / sum(weight, na.rm = TRUE),
                                    sum(centroids$y * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE))), crs = 4326)

    gravity.sf
  }
}

