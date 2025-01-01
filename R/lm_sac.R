


#' Linear model taking into account spatial autocorrelation
#'
#' ouput is lm.sac/custom object that can be displayed with modelsummary
#'
#' @param formula.chr
#' @param data.sf
#' @param knn_number
#' @param conley_cutoff
#' @param correlograms
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#'
lm.sac <- function(formula.chr, data.sf, knn_number = 20, conley_cutoff = 5, correlograms = F, ...) {
  # the goal is to spit out ONE single lm object that gets then assign also the custom class
  # 1 first we run the formula as an felm with lat lon, this way we can keep unique id

  # TODO
  # - the error when the ConleySE (e.g.) are messed up is always when we have unit FE as dummy when there is no FE in realty, with seg5 noproblem (active in Aug21)

  # TRY RIGHT FROM THE START TO FIX THE SSIZE ISSUE BY SUBSETTING away THE DEP VAR NA's
  # - this solves the issue that we usually have with the W matrices and morantests when a few obs are dropped bc of NA by the lm (the lm object still stores the full length of the df)
  dep.var <- stringr::str_split(formula.chr, "\\ ~ ")[[1]][1]
  data.sf <- data.sf %>% dplyr::filter(!is.na(!!as.symbol(dep.var))) # as symbol and then we need to unquote with !! (syntactic sugar). NEW DPLYR: dplyr::filter(!is.na(.[[dep.var]])) also works
  # the other issue could be that some of the Xes are missing and thus the Moran test on the residuals requires a different matrix? but the NA drop should work now in spdep (naSubset == TRUE by default)

  # 1) run the regression in lfe
  obj.felm <- lfe::felm(as.formula(formula.chr), data = data.sf, keepCX = T, keepModel = T, ...)
  obj.felm

  # ... start prepping to run the regression in lm
  # split off the FE, in the non-FE case it is the "unit" one
  formula.FEs <- stringr::str_split(formula.chr, "\\| ")[[1]][2]
  the.spatial.FE <- stringr::str_split(formula.FEs, " \\+ year")[[1]][1] # this guy we need then for the conley here and the cluster @ FE level in the post-regression
  # TODO case needed where we have a real year is time FE!

  # now fix the formula to
  formula.chr1 <- stringr::str_replace_all(formula.chr, "\\|", "+") # go to lm formula by replacing | with +
  formula.chr1 <- stringr::str_split(formula.chr1, " \\+ year")[[1]][1] # finish off by only taking first part
  # what about the error when we have factor with only 1 level (i.e. the unit and year) if they present, just split them off for the lm?
  # for now the quickfix seems to be to just call the unit 1 (numeric) instead of "unit", this way lm just phases out the var with an NA

  # THIS IS A MORE SUSTAINABLE FIX: ...
  if (stringr::str_detect(formula.chr1, "unit") == TRUE) {
    formula.chr1 <- stringr::str_split(formula.chr1, " \\+ unit")[[1]][1] # CUT-OFF the unit like we did with year above
  }

  # 2) run the regression in lm
  obj.lm <- lm(as.formula(formula.chr1), data = data.sf, ...)
  obj.lm$spatial_FE <- the.spatial.FE # pass this just as an info to our custom lm object

  # merge them onto each other (using the lat lon from the felm as unique identifiers)
  points.merged <- dplyr::right_join(data.sf, obj.felm$model, by = c("lat", "lon"), suffix = c("", ".y")) # still an sf.df, put .y only on one so the colnames stay original
  #print(nrow(points.merged))

  # WEIGHTs
  points_knn_20 <- spdep::knn2nb(spdep::knearneigh(st_coordinates(points.merged), k = knn_number)) # true W matrix
  #print(length(points_knn_20))
  #print(length(obj.lm$model[,1]))
  obj.lm$Moran_lmresid <- try(spdep::lm.morantest(obj.lm, spdep::nb2listw(points_knn_20))$statistic)
  # added a [1] here so i do not get the warning when everything went fine (bc then the class is [1] "matrix" "array")
  if(class(obj.lm$Moran_lmresid)[1] == "try-error") obj.lm$Moran_lmresid <- NA # error catcher
  obj.lm$Moran_response <- try(spdep::moran.test(obj.lm$model[,1], spdep::nb2listw(points_knn_20))$statistic) # here we use the actual response from model, even though the points.merged should be appropriately subsetted, thus equivalent
  if(class(obj.lm$Moran_response)[1] == "try-error") obj.lm$Moran_response <- NA # error catcher

  # CORELLOGRAM
  if (correlograms == T) {
    correlog.fit <- ncf::correlog(x = points.merged$lon, y = points.merged$lat,
                                  z = as.numeric(obj.lm$residuals), increment = 1,
                                  resamp = 1, quiet = TRUE, latlon = TRUE, na.rm = T)
    obj.lm$correlog.range_resid <- SpatialInference::extract.corr.range(correlog.fit)
    correlog.fit <- ncf::correlog(x = points.merged$lon, y = points.merged$lat,
                                  z = obj.lm$model[,1], increment = 1,
                                  resamp = 1, quiet = TRUE, latlon = TRUE, na.rm = T)
    obj.lm$correlog.range_response <- SpatialInference::extract.corr.range(correlog.fit)
  } else {obj.lm$correlog.range_resid <- NA; obj.lm$correlog.range_response <- NA}
  # CONLEY
  # cutoff <- obj.lm$correlog.range_resid
  cutoff <- conley_cutoff
  conley.1 <- SpatialInference::conley_SE(reg = obj.felm,
                                          unit = the.spatial.FE, time = "year",
                                          kernel = "bartlett", dist_fn = "Haversine",
                                          lat = "lat", lon =  "lon", dist_cutoff = cutoff)
  # here we make the stunt to get a vector of the length of n(coef), including the intercept since lm() always displays it
  obj.lm$conley_SE <- c(0, lapply(conley.1, function(x) diag(sqrt(x)))$Spatial[[1]], rep(0, length(obj.lm$coefficients)-2))

  if (correlograms == T) {
    cutoff <- obj.lm$correlog.range_resid
    conley.2 <- SpatialInference::conley_SE(reg = obj.felm,
                                            unit = the.spatial.FE, time = "year",
                                            kernel = "bartlett", dist_fn = "Haversine",
                                            lat = "lat", lon =  "lon", dist_cutoff = cutoff)
    # here we make the stunt to get a vector of the length of n(coef), including the intercept since lm() always displays it
    obj.lm$conley_SE_flex <- c(0, lapply(conley.2, function(x) diag(sqrt(x)))$Spatial[[1]], rep(0, length(obj.lm$coefficients)-2))
  } else {obj.lm$conley_SE_flex <- NA}

  # spit out
  class(obj.lm) <- c("custom", class(obj.lm))
  #class(obj.lm) <- c("lmsac")#, class(obj.lm))
  #class(obj.lm) <- c("custom")
  return(obj.lm)

}


# Tidy functions

# define tidy and glance methods
# felm base works only on the point estimates and robust HC3, cluster fails - it spits the CIS
tidy_custom.lm <- function(x, ...) {
  data.frame(
    #term = names(coef(x)),
    term = names(x$coefficients),
    estimate = unname(x$coefficients), # THE NAME OF THE VECTOR HERE DECIDES!!
    # PARAMETERS REFIT NOT WORKING WITH POLY IN REG
    #statistic = parameters::model_parameters(x, standardize = "refit")$Coefficient,
    #statistic = rep(5, length(x$coefficients)),
    std.error = sandwich::vcovHC(x, type = "HC1") %>% diag() %>% sqrt(),
    #std.error = unname(sqrt(diag(x$robustvcv))),
    #conf.low = sandwich::vcovCL(x, cluster = reformulate(x$spatial_FE)) %>% diag() %>% sqrt(),
    #statistic = sandwich::vcovCL(x, cluster = ~SubDistrict) %>% diag() %>% sqrt(),
    #conf.high = vcovCL(x, cluster = ~segment10) %>% diag() %>% sqrt()
    conf.high = x$conley_SE
  )
}

# extra lines works though in the glance
glance_custom.lm <- function(x, ...) {
  # just compute a weightsmatrix every time (bc of diff SSizes and subsets) ?
  data.frame(
    #"Mean_SD" = stringr::str_c(round(mean(x$model[,1]), 3), " [", round(sd(x$model[,1]), 3), "]"), #na.rm should not be necessary since lm removed NAs
    "Model" = "Custom",
    "y_mean" = mean(x$model[,1]),
    "y_SD"   = sd(x$model[,1]),
    "nobs" = stats::nobs(x),
    # here we need an if condition to use specific W_matrices based on distance subsets
    # PRLY BETTER TO CREATE THEM ON THE FLY?!
    "Moran_y" =     x$Moran_response,
    "Moran_resid" = x$Moran_lmresid,
    # issue when correlog is deactivated in lmsac regression? (else condition with NA way out?)
    "Range_resid" = x$correlog.range_resid,
    "Range_y" =  x$correlog.range_response
  )
}

# gof mapping for modelsummary printing
gm.param <- list(
  list("raw" = "nobs", "clean" = "Observations", "fmt" = 0),
  list("raw" = "Moran_y", "clean" = "Moran's I [y]", "fmt" = 3),
  list("raw" = "Moran_resid", "clean" = "Moran's I [resid]", "fmt" = 3),
  list("raw" = "y_mean", "clean" = "Mean [y]", "fmt" = 3),
  list("raw" = "y_SD", "clean" = "Std.Dev. [y]", "fmt" = 3)
)
