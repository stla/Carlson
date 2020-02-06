context("elliptic_PI")

test_that("Comparisons with Wolfram", {
  expect_equal(0.3468194+0.7198448i,
               elliptic_PI(1+1i, 2, -2), tolerance = 1e-7)
  expect_equal(0.67834548+1.05327009i,
               elliptic_PI(1+1i, 2, -2i), tolerance = 1e-7)
  expect_equal(0.671398646+0.79618888i,
               elliptic_PI(1+1i, 2+1i, -2i), tolerance = 1e-7)
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

test_that("n=1", {
  phi <- 1+1i
  m <- 2-1i
  expect_equal(
    elliptic_PI(phi, 1, m),
    (sqrt(1-m*sin(phi)^2)*tan(phi)-elliptic_E(phi,m))/(1-m) +
      elliptic_F(phi, m)
  )
})

test_that("m=0", {
  z <- 1+1i
  n <- 3 # does not work for complex n
  expect_equal(
    elliptic_PI(z, n, 0),
    atanh(sqrt(n-1)*tan(z)) / sqrt(n-1)
  )
})

test_that("m=1", {
  z <- 1+1i
  n <- 3
  expect_equal(
    elliptic_PI(z, n, 1),
    (sqrt(n) * atanh(sqrt(n)*sin(z)) - atanh(sin(z))) / (n-1)
  )
})

test_that("misc equalities", {
  n <- 2+2i
  m <- 3
  expect_equal(
    elliptic_PI(asin(1/sqrt(m)), n, m),
    elliptic_PI(pi/2, n/m, 1/m) / sqrt(m)
  )
})

test_that("Symmetry", {
  z <- -5 + 3i
  n <- 3 + 11i
  m <- -4 - 9i
  expect_equal(
    elliptic_PI(Conj(z), Conj(n), Conj(m)),
    Conj(elliptic_PI(z, n, m))
  )
})

test_that("Parity", {
  z <- -5 + 3i
  n <- 3 + 11i
  m <- -4 - 9i
  expect_equal(
    elliptic_PI(-z, n, m),
    -elliptic_PI(z, n, m)
  )
})

