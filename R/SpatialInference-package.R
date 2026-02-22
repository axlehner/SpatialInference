#' @keywords internal
#' @aliases SpatialInference-package
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom stats sd setNames
#' @useDynLib SpatialInference, .registration = TRUE
## usethis namespace: end
NULL

# Suppress R CMD check NOTEs for data.table NSE column references
# and ggplot2 aes() variables
utils::globalVariables(c("e", "time", "unit", "dist", "cutoff", "conley",
                          "gamma", "range_extr"))
