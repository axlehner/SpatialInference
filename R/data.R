#' Centroids of Contiguous US Counties (2017)
#'
#' An `sf` data frame containing the centroids of all 3,108 counties of the
#' contiguous United States (2017 geographies), along with synthetic
#' spatially-correlated noise variables for use in examples and vignettes.
#'
#' @format An `sf` data frame with 3,108 rows and the following columns:
#' \describe{
#'   \item{STATE}{Numeric state FIPS code.}
#'   \item{NAME}{County name.}
#'   \item{NAMELSAD}{County name with legal/statistical area description.}
#'   \item{GISJOIN}{Unique ID for joining with IPUMS data (2017 geographies).}
#'   \item{lat}{Latitude of the county centroid (WGS84).}
#'   \item{lon}{Longitude of the county centroid (WGS84).}
#'   \item{unit}{Cross-sectional unit identifier (constant `1` for
#'     cross-sectional use).}
#'   \item{year}{Time period identifier (constant `1` for cross-sectional
#'     use).}
#'   \item{noise1}{Synthetic spatially-correlated variable (noise 1).}
#'   \item{noise2}{Synthetic spatially-correlated variable (noise 2).}
#'   \item{dist}{Distance variable (spatially non-stationary example).}
#'   \item{geometry}{Point geometry column (NAD83 / EPSG:4269).}
#' }
#'
#' @source IPUMS NHGIS, University of Minnesota, \url{https://www.nhgis.org}.
"US_counties_centroids"
