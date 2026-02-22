#' Inverse-U Plot of Conley Standard Errors vs. Cutoff Distance
#'
#' Visualizes the relationship between the spatial bandwidth (cutoff distance)
#' and the resulting Conley standard error. This typically reveals an
#' inverse-U shaped relationship, helping to identify the appropriate bandwidth
#' for spatial HAC estimation.
#'
#' @param df.input An `sf` data frame containing the regression variables
#'   and columns `lat`, `lon`, `unit`, `year`.
#' @param cutoffrange Numeric vector of cutoff distances (in km) to evaluate.
#' @param kernel_choice_conley Character string specifying the kernel function.
#'   Default is `"epanechnikov"`. See [conley_SE()] for options.
#' @param depvar Character string naming the dependent variable column.
#'   Default is `"noise1"`.
#' @param indepvar Character string naming the independent variable column.
#'   Default is `"noise2"`.
#' @param range_add Logical. If `TRUE`, overlays the covariogram-estimated
#'   correlation range as a vertical red line. Default is `FALSE`.
#' @param ... Additional arguments (currently unused).
#'
#' @return A [ggplot2::ggplot()] object showing Conley SE (y-axis) against
#'   cutoff distance (x-axis), with the HC1 standard error as a grey
#'   dashed horizontal reference line.
#'
#' @references
#' Conley, T. G. (1999). GMM estimation with cross sectional dependence.
#' *Journal of Econometrics*, 92(1), 1--45.
#' \doi{10.1016/S0304-4076(98)00084-0}
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(US_counties_centroids)
#' if (requireNamespace("lfe", quietly = TRUE) &&
#'     requireNamespace("fixest", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   inverseu_plot_conleyrange(US_counties_centroids,
#'                             cutoffrange = seq(100, 1000, by = 200))
#' }
#' }
inverseu_plot_conleyrange <- function(df.input, cutoffrange = NA, kernel_choice_conley = "epanechnikov",
                                      depvar = "noise1", indepvar = "noise2",
                                      range_add = FALSE, ...) {

  if (!requireNamespace("lfe", quietly = TRUE))
    stop("Package 'lfe' is required for inverseu_plot_conleyrange().", call. = FALSE)
  if (!requireNamespace("fixest", quietly = TRUE))
    stop("Package 'fixest' is required for inverseu_plot_conleyrange().", call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for inverseu_plot_conleyrange().", call. = FALSE)

  reglfe <- lfe::felm(stats::as.formula(paste(depvar,"~", indepvar," | unit + year | 0 | lat + lon")), data = df.input, keepCX = TRUE)
  reghc1 <- fixest::feols(stats::as.formula(paste(depvar,"~", indepvar)), data = df.input, vcov = "HC1")
  hc1error <- reghc1$coeftable[2, 2]

  results <- data.frame(cutoff = cutoffrange, conley = rep(NA, length(cutoffrange)))
  for (i in 1:length(cutoffrange)) {
    results[i, 2] <- compute_conley_lfe(reglfe, cutoffrange[i], kernel_choice = "epanechnikov")
  }


  if(range_add == TRUE) {
    if (!requireNamespace("gstat", quietly = TRUE))
      stop("Package 'gstat' is required when range_add = TRUE.", call. = FALSE)
    covgm <- gstat::variogram(reghc1$residuals ~ 1, data = sf::as_Spatial(df.input), covariogram = TRUE, width = 2e4, cutoff = max(df.input$dist)*1e3, cressie = TRUE)
    range_extr <- extract_corr_range(covgm)
  }

  results |> ggplot2::ggplot(ggplot2::aes(x = cutoff, y = conley)) +
    ggplot2::geom_point() + ggplot2::geom_line() +
    ggplot2::geom_hline(ggplot2::aes(yintercept = hc1error), col = "grey", lty = 2) +
    ggplot2::theme(axis.title = ggplot2::element_blank())


}




#' Covariogram Range Estimation and Visualization
#'
#' Estimates a covariogram from an `sf` data frame (either from a single
#' variable or regression residuals) using [gstat::variogram()] with
#' `covariogram = TRUE`, extracts the zero-crossing as the estimated
#' correlation range via [extract_corr_range()], and produces a diagnostic
#' plot.
#'
#' @param df.input An `sf` data frame with point geometries.
#' @param depvar Character string naming the dependent variable. Default
#'   is `"noise1"`.
#' @param indepvar Character string naming the independent variable (used
#'   when `single.variable = FALSE`). Default is `"noise2"`.
#' @param maxdist Numeric. Maximum distance for the covariogram (in metres).
#'   Default is `NA`, which uses 2/3 of the maximum pairwise distance.
#' @param spacing Numeric. Bin width for the covariogram (in metres). Default
#'   is `NA`, which uses `maxdist / 150`.
#' @param single.variable Logical. If `TRUE`, computes the covariogram directly
#'   from `depvar`. If `FALSE` (default), first regresses `depvar` on
#'   `indepvar` via [fixest::feols()] and uses the residuals.
#'
#' @return A [ggplot2::ggplot()] object showing the covariogram with the
#'   estimated correlation range marked by a vertical red line and annotated
#'   text.
#'
#' @references
#' Pebesma, E. J. (2004). Multivariable geostatistics in S: the gstat
#' package. *Computers & Geosciences*, 30(7), 683--691.
#' \doi{10.1016/j.cageo.2004.03.012}
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(US_counties_centroids)
#' if (requireNamespace("fixest", quietly = TRUE) &&
#'     requireNamespace("gstat", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   covgm_range(US_counties_centroids)
#' }
#' }
covgm_range <- function(df.input,
                        depvar = "noise1", indepvar = "noise2",
                        maxdist = NA, spacing = NA,
                        single.variable = FALSE
) {

  if (!requireNamespace("gstat", quietly = TRUE))
    stop("Package 'gstat' is required for covgm_range().", call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for covgm_range().", call. = FALSE)

  if(is.na(maxdist)) maxdist <- as.numeric(max(sf::st_distance(df.input))) * (2/3)
  if(is.na(spacing)) spacing <- maxdist / 150

  # remove NAs
  df.input <- df.input[!is.na(df.input[[depvar]]), ]

  if (single.variable == TRUE) {
    covgm <- gstat::variogram(df.input[[depvar]] ~ 1, data = sf::as_Spatial(df.input),
                              covariogram = TRUE, width = spacing, cutoff = maxdist, cressie = TRUE)
  } else {
    if (!requireNamespace("fixest", quietly = TRUE))
      stop("Package 'fixest' is required when single.variable = FALSE.", call. = FALSE)
    reghc1 <- fixest::feols(stats::as.formula(paste(depvar,"~", indepvar)), data = df.input, vcov = "HC1")
    covgm <- gstat::variogram(reghc1$residuals ~ 1, data = sf::as_Spatial(df.input),
                              covariogram = TRUE, width = spacing, cutoff = maxdist, cressie = TRUE)
  }

  # then we extract, doesn't matter which covgm
  range_extr <- extract_corr_range(covgm)

  plottext <- data.frame(range_extr = range_extr)

  y_position <- max(covgm$gamma) * .9
  x_position <- range_extr * 1.1
  covgm |> ggplot2::ggplot(ggplot2::aes(x = dist/1e3, y = gamma)) +
    ggplot2::geom_hline(yintercept = 0, color = "red", lty = 8) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = range_extr), col = "red") +
    ggplot2::geom_point(alpha = .5) + ggplot2::xlab("Distance in km") +
    ggplot2::geom_text(
      data = plottext,
      ggplot2::aes(
        x = x_position,
        y = y_position,
        label = paste0("Estimated range:\n", signif(range_extr, 3), "km")
      ),
      inherit.aes = FALSE,
      size = 4,
      hjust = 0,
      vjust = 1.5
    )


}
