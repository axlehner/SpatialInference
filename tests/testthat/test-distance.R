test_that("DistMat returns square symmetric matrix with 1s on diagonal", {
  # 3 points: NYC, LA, Chicago (approx lat/lon)
  M <- matrix(c(
    40.71, -74.01,   # NYC
    34.05, -118.24,  # LA
    41.88, -87.63    # Chicago
  ), ncol = 2, byrow = TRUE)

  d <- DistMat(M, cutoff = 5000, kernel = "bartlett")

  expect_true(is.matrix(d))
  expect_equal(nrow(d), 3)
  expect_equal(ncol(d), 3)

  # Diagonal should be 1
  expect_equal(diag(d), c(1, 1, 1))

  # Symmetric
  expect_equal(d, t(d))
})

test_that("points beyond cutoff get zero weight", {
  # NYC and LA are ~3940 km apart
  M <- matrix(c(
    40.71, -74.01,   # NYC
    34.05, -118.24   # LA
  ), ncol = 2, byrow = TRUE)

  d_small <- DistMat(M, cutoff = 100, kernel = "bartlett")

  # With cutoff of 100 km, NYC-LA should get 0 weight
  expect_equal(d_small[1, 2], 0)
  expect_equal(d_small[2, 1], 0)
})

test_that("known Haversine distance check (NYC to LA)", {
  # NYC and LA are roughly 3940 km apart
  M <- matrix(c(
    40.71, -74.01,   # NYC
    34.05, -118.24   # LA
  ), ncol = 2, byrow = TRUE)

  # With cutoff of 5000 km and bartlett kernel, weight = 1 - d/cutoff
  d <- DistMat(M, cutoff = 5000, kernel = "bartlett")

  # Weight should be positive but less than 1
  expect_true(d[1, 2] > 0)
  expect_true(d[1, 2] < 1)

  # Approximate check: d ~3940, weight ~= 1 - 3940/5000 = 0.212
  expect_true(abs(d[1, 2] - 0.212) < 0.05)
})

test_that("uniform kernel gives equal weight within cutoff", {
  # Two points close together
  M <- matrix(c(
    40.71, -74.01,
    40.72, -74.02
  ), ncol = 2, byrow = TRUE)

  d <- DistMat(M, cutoff = 100, kernel = "uniform")

  # Uniform kernel should give weight 1 for points within cutoff
  expect_equal(d[1, 2], 1)
})
