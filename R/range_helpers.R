


#' Title
#'
#' @param df.input
#' @param cutoffrange
#' @param kernel_choice_conley
#' @param depvar
#' @param indepvar
#' @param range_add
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
inverseu_plot_conleyrange <- function(df.input, cutoffrange = NA, kernel_choice_conley = "epanechnikov",
                                      depvar = "noise1", indepvar = "noise2", # initial usecase was for indepvar (bc of dist), but written more general from outside, who knows if we need it later
                                      range_add = FALSE, ...) {

  #if(is.na(cutoffrange)) print("cutoffrange needed \n")
  reglfe <- lfe::felm(as.formula(paste(depvar,"~", indepvar," | unit + year | 0 | lat + lon")), data = df.input, keepCX = T)
  reghc1 <- fixest::feols(as.formula(paste(depvar,"~", indepvar)), data = df.input, vcov = "HC1")
  hc1error <- reghc1$coeftable[2, 2]

  results <- data.frame(cutoff = cutoffrange, conley = rep(NA, length(cutoffrange))) # holder
  for (i in 1:length(cutoffrange)) {
    results[i, 2] <- compute_conley_lfe(reglfe, cutoffrange[i], kernel_choice = "epanechnikov")
  }


  if(range_add == TRUE) {
    covgm <- gstat::variogram(reghc1$residuals ~ 1, data = sf::as_Spatial(df.input), covariogram = TRUE, width = 2e4, cutoff = max(df.input$dist)*1e3, cressie = TRUE)
    range_extr <- extract.corr.range(covgm)
  }

  results |> ggplot(aes(x = cutoff, y = conley)) +
    geom_point() + geom_line() +
    #if(range_add == TRUE) geom_vline(aes(xintercept = range_extr), col = "red") +
    geom_hline(aes(yintercept = hc1error), col = "grey", lty = 2) +
    theme(axis.title = element_blank())


}




#' Title
#'
#' @param df.input
#' @param depvar
#' @param indepvar
#' @param maxdist
#' @param spacing
#' @param single.variable
#'
#' @return
#' @export
#'
#' @examples
covgm_range <- function(df.input,
                        depvar = "noise1", indepvar = "noise2",
                        maxdist = NA, spacing = NA,
                        single.variable = FALSE # default is regression residual covgm
) {

  # TODO:
  # - could add a Moran's I into the text
  # if range max (x-axis) not specified, take 2/3 of max distance between points (assumption is it comes in meters)
  if(is.na(maxdist)) maxdist <- as.numeric(max(sf::st_distance(df.input))) * (2/3) # wants meters, below we had the distkm times 1e3
  if(is.na(spacing)) spacing <- maxdist / 150 # make it roughly 150 dots (that amounts to 2e4 as in my example cases)

  # remove NAs
  df.input <- df.input[!is.na(df.input[[depvar]]), ]

  if (single.variable == TRUE) {
    covgm <- gstat::variogram(df.input[[depvar]] ~ 1, data = sf::as_Spatial(df.input),
                              covariogram = TRUE, width = spacing, cutoff = maxdist, cressie = TRUE)
  } else { # if no variable, estimate regression:
    # estimate naive regression
    reghc1 <- fixest::feols(as.formula(paste(depvar,"~", indepvar)), data = df.input, vcov = "HC1")
    # estimate the covariogram
    covgm <- gstat::variogram(reghc1$residuals ~ 1, data = sf::as_Spatial(df.input),
                              covariogram = TRUE, width = spacing, cutoff = maxdist, cressie = TRUE)
  }

  # then we extract, doesn't matter which covgm
  range_extr <- extract_corr_range(covgm)

  plottext <- data.frame(range_extr = range_extr)

  y_position <- max(covgm$gamma) * .9
  #x_position <- maxdist/1e3 * .9
  x_position <- range_extr * 1.1
  covgm |> ggplot(aes(x = dist/1e3, y = gamma)) +
    geom_hline(yintercept = 0, color = "red", lty = 8) +
    geom_vline(aes(xintercept = range_extr), col = "red") +
    geom_point(alpha = .5) + xlab("Distance in km") +
    geom_text(
      data = plottext,
      aes(
        x = x_position,
        y = y_position,
        label = paste0("Estimated range:\n", signif(range_extr, 3), "km")
      ),
      inherit.aes = FALSE, # Prevent it from inheriting global aesthetics
      size = 4,
      hjust = 0,
      vjust = 1.5
    )


}
