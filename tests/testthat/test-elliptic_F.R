context("elliptic_F")

test_that("phi = pi/2", {
  expect_equal(elliptic_F(pi/2,-1), as.complex(gamma(1/4)^2/4/sqrt(2*pi)))
  expect_equal(elliptic_F(pi/2,1/2),
               as.complex(8*pi^(3/2)/gamma(-1/4)^2))
  expect_equal(elliptic_F(pi/2, 0.5),
               pi/2*as.complex(gsl::hyperg_2F1(1/2,1/2,1,0.5)))
})
