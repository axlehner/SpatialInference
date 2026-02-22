test_that("covgm_range returns a ggplot with estimated range", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("gstat")
  skip_if_not_installed("ggplot2")
  data(US_counties_centroids)

  p <- covgm_range(US_counties_centroids)

  expect_s3_class(p, "ggplot")
})

test_that("covgm_range single.variable mode works", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("ggplot2")
  data(US_counties_centroids)

  p <- covgm_range(US_counties_centroids, single.variable = TRUE)

  expect_s3_class(p, "ggplot")
})
