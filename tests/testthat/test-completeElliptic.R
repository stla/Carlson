context("Complete elliptic integrals")

test_that("Complete F", {
  m <- 2 - 3i
  expect_equal(
    elliptic_F(pi/2, m),
    Carlson_RF(0, 1-m, 1)
  )
})

test_that("Complete E", {
  m <- 2 - 3i
  expect_equal(
    elliptic_E(pi/2, m),
    2*Carlson_RG(0, 1-m, 1)
  )
  expect_equal(
    elliptic_E(pi/2, m),
    (1-m) * (Carlson_RD(0, 1-m, 1) + Carlson_RD(0, 1, 1-m)) / 3
  )
})

test_that("Complete F minus complete E", {
  m <- 2 - 3i
  expect_equal(
    elliptic_F(pi/2, m) - elliptic_E(pi/2, m),
    m * Carlson_RD(0, 1-m, 1) / 3
  )
})

test_that("Complete E minus (1-m)*complete F", {
  m <- 2 - 3i
  expect_equal(
    elliptic_E(pi/2, m) - (1-m)*elliptic_F(pi/2, m),
    m*(1-m) * Carlson_RD(0, 1, 1-m) / 3
  )
})
