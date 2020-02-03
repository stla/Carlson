context("Carlson_RJ")

test_that("RJ(x,y,y,p)", {
  x <- 1 + 1i
  y <- -2 + 3i
  p <- 4i
  expect_equal(Carlson_RJ(x, y, y, p),
               3*(Carlson_RC(x,y) - Carlson_RC(x,p)) / (p-y))
})

test_that("homogeneity", {
  x <- 1 + 1i
  y <- -2 + 3i
  z <- -3
  p <- 4i
  kappa <- 2
  expect_equal(Carlson_RJ(x, y, z, p)/kappa/sqrt(kappa),
               Carlson_RJ(kappa*x, kappa*y, kappa*z, kappa*p))
})
