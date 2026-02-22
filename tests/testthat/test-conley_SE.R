test_that("conley_SE returns correctly structured output", {
  skip_if_not_installed("lfe")
  data(US_counties_centroids)

  reg <- lfe::felm(noise1 ~ noise2 | unit + year | 0 | lat + lon,
                    data = US_counties_centroids, keepCX = TRUE)
  result <- conley_SE(reg, unit = "unit", time = "year",
                      lat = "lat", lon = "lon",
                      dist_cutoff = 500)

  # Returns a list of 3 named elements

  expect_type(result, "list")
  expect_length(result, 3)
  expect_named(result, c("OLS", "Spatial", "Spatial_HAC"))

  # Each element is a matrix with correct dimensions (k x k where k = 1 regressor)
  k <- length(rownames(reg$coefficients))
  for (nm in names(result)) {
    expect_true(is.matrix(result[[nm]]))
    expect_equal(nrow(result[[nm]]), k)
    expect_equal(ncol(result[[nm]]), k)
  }
})

test_that("spatial VCV diagonal entries are non-negative", {
  skip_if_not_installed("lfe")
  data(US_counties_centroids)

  reg <- lfe::felm(noise1 ~ noise2 | unit + year | 0 | lat + lon,
                    data = US_counties_centroids, keepCX = TRUE)
  result <- conley_SE(reg, unit = "unit", time = "year",
                      lat = "lat", lon = "lon",
                      dist_cutoff = 500)

  expect_true(all(diag(result$Spatial) >= 0))
  expect_true(all(diag(result$Spatial_HAC) >= 0))
})

test_that("spatial VCV matrices are symmetric", {
  skip_if_not_installed("lfe")
  data(US_counties_centroids)

  reg <- lfe::felm(noise1 ~ noise2 | unit + year | 0 | lat + lon,
                    data = US_counties_centroids, keepCX = TRUE)
  result <- conley_SE(reg, unit = "unit", time = "year",
                      lat = "lat", lon = "lon",
                      dist_cutoff = 500)

  expect_equal(result$Spatial, t(result$Spatial))
  expect_equal(result$Spatial_HAC, t(result$Spatial_HAC))
})

test_that("different kernels produce different results", {
  skip_if_not_installed("lfe")
  data(US_counties_centroids)

  reg <- lfe::felm(noise1 ~ noise2 | unit + year | 0 | lat + lon,
                    data = US_counties_centroids, keepCX = TRUE)

  res_bartlett <- conley_SE(reg, unit = "unit", time = "year",
                            lat = "lat", lon = "lon",
                            kernel = "bartlett", dist_cutoff = 500)
  res_uniform  <- conley_SE(reg, unit = "unit", time = "year",
                            lat = "lat", lon = "lon",
                            kernel = "uniform", dist_cutoff = 500)

  # Bartlett down-weights distant obs, uniform does not — SEs should differ
  expect_false(isTRUE(all.equal(res_bartlett$Spatial, res_uniform$Spatial)))
})

test_that("non-felm input produces informative error", {
  reg_lm <- lm(mpg ~ wt, data = mtcars)
  expect_error(
    conley_SE(reg_lm, unit = "unit", time = "year",
              lat = "lat", lon = "lon"),
    "felm"
  )
})

test_that("compute_conley_lfe returns a single numeric", {
  skip_if_not_installed("lfe")
  data(US_counties_centroids)

  reg <- lfe::felm(noise1 ~ noise2 | unit + year | 0 | lat + lon,
                    data = US_counties_centroids, keepCX = TRUE)
  se <- compute_conley_lfe(reg, cutoff = 500)

  expect_type(se, "double")
  expect_length(se, 1)
  expect_true(se > 0)
})

test_that("compute_conley_lfe matches manual extraction", {
  skip_if_not_installed("lfe")
  data(US_counties_centroids)

  reg <- lfe::felm(noise1 ~ noise2 | unit + year | 0 | lat + lon,
                    data = US_counties_centroids, keepCX = TRUE)

  se_convenience <- compute_conley_lfe(reg, cutoff = 500)

  # Manual extraction
  vcvs <- conley_SE(reg, unit = "unit", time = "year",
                    lat = "lat", lon = "lon", dist_cutoff = 500)
  se_manual <- sqrt(diag(vcvs$Spatial)) |> round(5)

  expect_equal(se_convenience, unname(se_manual))
})
