#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' Extract Correlation Range
#'
#' Function will pick the point where the correlogram first crosses the x-axis (the distance), i.e. where the spatial autocorrelation gets 0 for the first time.
#' This is a good approximation of the actual correlation range.
#'
#' @param input a correlogram estimated with ncf::correlog() or a covariogram estimated with gstat::variogram(..., covariogram = TRUE)
#'
#' @return an integer that represents the correlation range (in units of input, i.e. with lon lat it is km, with equal area CRS it typically is meters)
#' @export
#'
#' @examples
#'
#'
extract.corr.range <- function(input) {
  # TODO: think about projections and the input widths here
  # later we need an optimizer that reduces the width(covgm)/increment(ncfcorrelog) and picks the first crossing with 0
  stopifnot(
    "Input has to be a correlogram (ncf) or a covariogram (gstat)." = inherits(input, c("correlog", "gstatVariogram")),
    "Did you input a variogram? Please estimate a covariogram." = !(inherits(input, "gstatVariogram") & (min(input$gamma) > 0))
  )
  # if(inherits(input, "gstatVariogram") & (min(input$gamma) > 0))
  #   stop("Did you input a variogram? Please estimate a covariogram.")

  if (inherits(input, "correlog")) {
    holdout <- fit1$correlation |> round(5) |> floor() |> diff()
    # careful: the diff drops the first element and reduces the length by one
    # then we select where the first jump is (and also drop the naming here, even though not needed)
    pos <- holdout[abs(holdout) > 0][1] |> names() |> as.numeric() # instead of unname and having the first non-zero, we can pick the name
    # then we select the according range (unname and not go for kms for now because we wanna plot it!)
    # fit1$mean.of.class[which(holdout != 0)[1]] %>% unname #/ 1e3
    # this correctly selects by position and not the which(first nonzero) which is on the wrong length (see above)
    after <- fit1$mean.of.class[pos] |> unname() #/ 1e3
    # now we pick the one before (just before it gets negative)
    before <- fit1$mean.of.class[pos-1] |> unname() #/ 1e3
    estim_range <- (after+before)/2
  }
  if (inherits(input, "gstatVariogram")) {
    # same logic here with the covariogram estimate
    holdout <- input$gamma |> round(5) |> floor() |> diff()
    pos <- holdout[abs(holdout) > 0][1] |> names() |> as.numeric()
    before <- input$dist[which(holdout != 0)[1]] # this is the one before the crossing because the diff cut the holdout vector short by 1 element
    after <- input$dist[which(holdout != 0)[1]+1]
    estim_range <- (after+before)/2 / 1e3
  }

  return(estim_range)

}
