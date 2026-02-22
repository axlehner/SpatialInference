#' Extract Correlation Range from a Correlogram or Covariogram
#'
#' Identifies the distance at which spatial autocorrelation first crosses
#' zero, providing an estimate of the spatial correlation range. Works with
#' correlograms from [ncf::correlog()] and covariograms from
#' [gstat::variogram()] (with `covariogram = TRUE`).
#'
#' For correlograms, the function detects the first sign change in the
#' rounded and floored correlation values. For covariograms, it finds the
#' first index where gamma transitions from positive to non-positive. In
#' both cases, the estimated range is the midpoint between the last positive
#' and first non-positive distance bins.
#'
#' @param input A correlogram object from [ncf::correlog()] (class
#'   `"correlog"`) or a covariogram from [gstat::variogram()] (class
#'   `"gstatVariogram"`). Must be a covariogram (not a variogram) when
#'   using gstat.
#' @param returnzeroifNA Logical. If `TRUE`, returns `1` instead of `NA`
#'   when no zero-crossing is found. Default is `FALSE`.
#'
#' @return A numeric value representing the estimated correlation range.
#'   For covariograms, the unit is km (distance in metres divided by 1000).
#'   For correlograms, the unit matches the input distance unit.
#'
#' @references
#' Pebesma, E. J. (2004). Multivariable geostatistics in S: the gstat
#' package. *Computers & Geosciences*, 30(7), 683--691.
#' \doi{10.1016/j.cageo.2004.03.012}
#'
#' @export
#'
#' @examples
#' # With a mock gstatVariogram:
#' mock_vgm <- data.frame(
#'   np = rep(100, 10),
#'   dist = seq(50000, 500000, by = 50000),
#'   gamma = c(5, 3, 2, 1, 0.5, -0.2, -0.5, -0.3, -0.1, -0.05),
#'   dir.hor = 0, dir.ver = 0, id = "var1"
#' )
#' class(mock_vgm) <- c("gstatVariogram", "data.frame")
#' extract_corr_range(mock_vgm)
extract_corr_range <- function(input, returnzeroifNA = FALSE) {
  if(inherits(input, "correlog")) input$gamma <- -1
  stopifnot(
    "Input has to be a correlogram (ncf) or a covariogram (gstat)." = inherits(input, c("correlog", "gstatVariogram"))
  )

  if (inherits(input, "correlog")) {
    holdout <- input$correlation |> round(5) |> floor() |> diff()
    pos <- holdout[abs(holdout) > 0][1] |> names() |> as.numeric()
    after <- input$mean.of.class[pos] |> unname()
    before <- input$mean.of.class[pos-1] |> unname()
    estim_range <- (after+before)/2
  }
  if (inherits(input, "gstatVariogram")) {
    series <- input$gamma
    crossing_index <- which(series[-length(series)] > 0 & series[-1] <= 0)[1]
    before <- input$dist[crossing_index]
    after <- input$dist[crossing_index + 1]
    estim_range <- (after+before)/2 / 1e3
  }

  if(returnzeroifNA){if(is.na(estim_range > 0)){estim_range <- 1}}
  return(estim_range)

}
