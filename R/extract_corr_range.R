#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' Extract Correlation Range
#'
#' Function will pick the point where the correlogram first crosses the x-axis (the distance), i.e. where the spatial autocorrelation gets 0 for the first time.
#' This is a good approximation of the actual correlation range.
#'
#' @param ncf_corr a correlogram estimated by ncf::correlog()
#'
#' @return an integer that represents the correlation range (in units of input, i.e. with lon lat it is km, with equal area CRS it tipically is metres)
#' @export
#'
#' @examples
#'
#'
extract.corr.range <- function(ncf_corr) {
  # with the floor and diff function we mark where the breaks are in the series (i.e. when the transition(s) from negative to positive correlation happens: the range)
  holdout <- ncf_corr$correlation %>% round(5) %>% floor %>% diff
  # then we select where the first jump is (and also drop the naming here, even though not needed)
  holdout[abs(holdout) > 0][1] %>% unname
  # then we select the according range (unname and not go for kms for now because we wanna plot it!)
  ncf_corr$mean.of.class[which(holdout != 0)[1]] %>% unname #/ 1e3
}
