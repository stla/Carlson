context("elliptic_Z")

test_that("Representation in terms of elliptic_PI", {
  z <- 7 - 6i
  m <- -3
  expect_equal(
    elliptic_Z(z,m),
    (1-m)*elliptic_PI(z,m,m) + m*sin(2*z)/2/sqrt(1-m*sin(z)^2) -
      elliptic_E(pi/2,m)/elliptic_F(pi/2,m)*elliptic_PI(z,0,m)
  )
})
