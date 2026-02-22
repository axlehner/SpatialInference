test_that("extract_corr_range finds zero-crossing in gstatVariogram", {
  # Mock a covariogram with a clear zero-crossing between bins 5 and 6
  mock_vgm <- data.frame(
    np = rep(100, 10),
    dist = seq(50000, 500000, by = 50000),
    gamma = c(5, 3, 2, 1, 0.5, -0.2, -0.5, -0.3, -0.1, -0.05),
    dir.hor = 0, dir.ver = 0, id = "var1"
  )
  class(mock_vgm) <- c("gstatVariogram", "data.frame")

  result <- extract_corr_range(mock_vgm)

  # Crossing is between dist[5]=250000 and dist[6]=300000
  # Midpoint = 275000 / 1e3 = 275
  expect_equal(result, 275)
})

test_that("extract_corr_range finds zero-crossing in correlog", {
  # Mock a correlog with a clear sign change.
  # Names on correlation are needed because the algorithm uses names() on diff().
  mock_correlog <- list(
    correlation = setNames(c(0.8, 0.5, 0.2, 0.05, -0.1, -0.2), 1:6),
    mean.of.class = c(50, 100, 200, 300, 400, 500),
    n = rep(100, 6),
    x.intercept = 350
  )
  class(mock_correlog) <- "correlog"

  result <- extract_corr_range(mock_correlog)

  # The floor/diff method detects the transition 0 to -1 at position 5
  # (correlation goes from 0.05 to -0.1, floor from 0 to -1)
  # mean.of.class[5] = 400, mean.of.class[4] = 300
  # midpoint = 350
  expect_equal(result, 350)
})

test_that("returnzeroifNA returns 1 when no crossing found", {
  # Covariogram that never crosses zero
  mock_vgm <- data.frame(
    np = rep(100, 5),
    dist = seq(50000, 250000, by = 50000),
    gamma = c(5, 4, 3, 2, 1),
    dir.hor = 0, dir.ver = 0, id = "var1"
  )
  class(mock_vgm) <- c("gstatVariogram", "data.frame")

  result <- extract_corr_range(mock_vgm, returnzeroifNA = TRUE)
  expect_equal(result, 1)
})

test_that("invalid input class produces error", {
  expect_error(
    extract_corr_range(data.frame(x = 1:5)),
    "correlogram.*covariogram"
  )
})
