context("elliptic_E")

test_that("phi = pi/2", {
  expect_equal(elliptic_E(pi/2,-1), sqrt(2)*elliptic_E(pi/2,1/2))
  # expect_equal(elliptic_E(pi/2, 0.5),
  #              pi/2*as.complex(gsl::hyperg_2F1(-1/2,1/2,1,0.5))) bug gsl!!
  z <- exp(1i*pi/7)
  expect_equal(elliptic_E(pi/2, 1-1/z),
               elliptic_E(pi/2, 1-z)/sqrt(z))
})
