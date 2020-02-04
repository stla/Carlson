context("elliptic_F")

test_that("phi = pi/2", {
  expect_equal(elliptic_F(pi/2,-1), as.complex(gamma(1/4)^2/4/sqrt(2*pi)))
  expect_equal(elliptic_F(pi/2,1/2),
               as.complex(8*pi^(3/2)/gamma(-1/4)^2))
  expect_equal(elliptic_F(pi/2, 0.5),
               pi/2*as.complex(gsl::hyperg_2F1(1/2,1/2,1,0.5)))
  z <- 5+2i
  expect_equal(sqrt(1-z)*elliptic_F(pi/2,z),
               elliptic_F(pi/2, z/(z-1)))
  expect_equal(elliptic_F(pi/2,1/z),
               sqrt(z)*(elliptic_F(pi/2,z) - 1i*elliptic_F(pi/2,1-z)))
})

test_that("Representation in terms of elliptic_PI", {
  z <- 7 - 6i
  m <- -3
  expect_equal(elliptic_F(z,m), elliptic_PI(z,0,m))
})
