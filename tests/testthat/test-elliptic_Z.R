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

test_that("Special values", {
  z <- -1 + 8i
  expect_equal(elliptic_Z(z,0), as.complex(0))
  expect_equal(elliptic_Z(z,1), sin(z))
  expect_equal(elliptic_Z(0, 2+2i), as.complex(0))
  expect_equal(elliptic_Z(5*pi/2, 2+2i), as.complex(0))
})
