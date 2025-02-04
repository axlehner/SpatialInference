#' Centroids of contiguous US counties (2017)
#'
#' sf data frame with 3108 observations and NAD83 datum
#'
#' @format ## `US_counties_centroids`
#' A spatial data frame with 3,108 rows:
#' \describe{
#'   \item{STATE}{Code for the state the county is in}
#'   \item{NAME, NAMELSAD}{Name of county}
#'   \item{GISJOIN}{Unique ID to join with IPUMS data at the county level (2017 geographies)}
#'   ...
#' }
#' @source IPUMS
"US_counties_centroids"
