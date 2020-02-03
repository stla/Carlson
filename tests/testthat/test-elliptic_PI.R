context("elliptic_PI")

test_that("Comparisons with Wolfram", {
  expect_equal(0.3468194+0.7198448i, elliptic_PI(1+1i, 2, -2), tolerance = 1e-7)
  expect_equal(0.67834548+1.05327009i, elliptic_PI(1+1i, 2, -2i), tolerance = 1e-7)
  expect_equal(0.671398646+0.79618888i, elliptic_PI(1+1i, 2+1i, -2i), tolerance = 1e-7)
})

