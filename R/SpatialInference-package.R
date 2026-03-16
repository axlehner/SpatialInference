#' SpatialInference: Conley Spatial HAC Standard Errors and Correlation Range Estimation
#'
#' Fast computation of Conley (1999) spatial heteroskedasticity and
#' autocorrelation consistent (HAC) standard errors for linear regression
#' models with geo-coded data. Performance-critical distance calculations,
#' kernel weighting, and variance component accumulation are implemented in
#' C++ via Rcpp and RcppArmadillo. Includes tools for estimating the spatial
#' correlation range from covariograms and correlograms following the
#' bandwidth selection method proposed in Lehner (2026), and diagnostic
#' visualizations for bandwidth selection.
#'
#' @references
#' Lehner, A. (2026). Bandwidth selection for spatial HAC standard errors.
#' *arXiv preprint* arXiv:2603.03997. \doi{10.48550/arXiv.2603.03997}
#'
#' Conley, T. G. (1999). GMM estimation with cross sectional dependence.
#' *Journal of Econometrics*, 92(1), 1--45.
#' \doi{10.1016/S0304-4076(98)00084-0}
#'
#' Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless R and C++
#' integration. *Journal of Statistical Software*, 40(8), 1--18.
#' \doi{10.18637/jss.v040.i08}
#'
#' Eddelbuettel, D. and Sanderson, C. (2014). RcppArmadillo: Accelerating R
#' with high-performance C++ linear algebra. *Computational Statistics &
#' Data Analysis*, 71, 1054--1063. \doi{10.1016/j.csda.2013.02.005}
#'
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
