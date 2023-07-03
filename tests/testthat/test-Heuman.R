test_that("Heuman", {
  z <- Lambda0(pi/2, 1/2)
  expect_equal(Re(z), 1)
  expect_equal(Im(z), 0)
  z <- Lambda0(asin(sqrt(8)/3), 0)
  expect_equal(Re(z), sqrt(8)/3)
  expect_equal(Im(z), 0)
  z <- Lambda0(0, 1/sqrt(2))
  expect_equal(Re(z), 0)
  expect_equal(Im(z), 0)
})

test_that("Relation Heuman ellipticPI", {
  p <- 0.5
  nu <- 0.7
  expect_equal(
    elliptic_PI(pi/2, nu, p^2),
    pi/2*sqrt(nu/((1-nu)*(nu-p^2))) *
      Lambda0(asin(sqrt((nu-p^2)/(nu*(1-p^2)))), p^2)
  )
})
