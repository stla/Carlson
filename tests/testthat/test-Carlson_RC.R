context("Carlson_RC")

test_that("Relations with 'arc' functions", {
  x <- 1+1i
  y <- 3
  expect_equal(x*Carlson_RC(y*y, y*y-x*x), atanh(x/y))
  expect_equal(x*Carlson_RC(y*y-x*x, y*y), asin(x/y))
  expect_equal(x*Carlson_RC(y*y+x*x, y*y), asinh(x/y))
  expect_equal(sqrt(y*y-x*x)*Carlson_RC(x*x, y*y), acos(x/y))
  expect_equal(sqrt(x*x-y*y)*Carlson_RC(x*x, y*y), acosh(x/y))
  y <- -5+2i
  expect_equal(x*Carlson_RC(y*y, y*y-x*x), -atanh(x/y))
  expect_equal(x*Carlson_RC(y*y-x*x, y*y), -asin(x/y))
  expect_equal(x*Carlson_RC(y*y+x*x, y*y), -asinh(x/y))
  expect_equal(sqrt(y*y-x*x)*Carlson_RC(x*x, y*y), acos(-x/y))
  expect_equal(sqrt(x*x-y*y)*Carlson_RC(x*x, y*y), acosh(-x/y))
})
