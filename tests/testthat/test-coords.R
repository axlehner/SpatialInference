test_that("coords_as_columns returns tibble with correct names and dimensions", {
  data(US_counties_centroids)

  result <- coords_as_columns(US_counties_centroids)

  expect_s3_class(result, "tbl_df")
  expect_equal(ncol(result), 2)
  expect_equal(names(result), c("x", "y"))
  expect_equal(nrow(result), nrow(US_counties_centroids))
})

test_that("coords_as_columns respects custom names", {
  data(US_counties_centroids)

  result <- coords_as_columns(US_counties_centroids, names = c("lon", "lat"))

  expect_equal(names(result), c("lon", "lat"))
})

test_that("gravity_centroid returns sfc_POINT", {
  data(US_counties_centroids)

  result <- gravity_centroid(US_counties_centroids)

  expect_true(inherits(result, "sfc"))
  expect_true(inherits(result, "sfc_POINT"))
  expect_length(result, 1)
})

test_that("gravity_centroid weighted result differs from unweighted", {
  data(US_counties_centroids)

  unweighted <- gravity_centroid(US_counties_centroids)
  # Use noise1 as a weight (shift the centroid)
  weighted <- gravity_centroid(US_counties_centroids,
                                weight = abs(US_counties_centroids$noise1) + 1)

  coords_uw <- sf::st_coordinates(unweighted)
  coords_w  <- sf::st_coordinates(weighted)

  # Should not be exactly the same

  expect_false(isTRUE(all.equal(coords_uw, coords_w)))
})

test_that("grid_FE returns factor of correct length", {
  data(US_counties_centroids)

  result <- grid_FE(US_counties_centroids, size = 5)

  expect_true(is.factor(result))
  # After intersection some points may be dropped, but most should be retained
  expect_true(length(result) > 0)
  expect_true(length(result) <= nrow(US_counties_centroids))
})
