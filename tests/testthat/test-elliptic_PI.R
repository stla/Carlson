context("elliptic_PI")

test_that("Comparisons with Wolfram", {
  expect_equal(0.3468194+0.7198448i, elliptic_PI(1+1i, 2, -2), tolerance = 1e-7)
  expect_equal(0.67834548+1.05327009i, elliptic_PI(1+1i, 2, -2i), tolerance = 1e-7)
  expect_equal(0.671398646+0.79618888i, elliptic_PI(1+1i, 2+1i, -2i), tolerance = 1e-7)
})

test_that("n=0", {
  expect_equal(elliptic_PI(7+1i,0,2-1i), elliptic_F(7+1i,2-1i))
})

test_that("n=m", {
  phi <- 7+1i
  m <- 2-1i
  expect_equal(
    elliptic_E(phi, m),
    (1-m)*elliptic_PI(phi, m, m) + m*sin(2*phi)/2/sqrt(1-m*sin(phi)^2)
  )
})
