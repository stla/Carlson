context("Carlson_RG")

test_that("Some values of Carlson_RG", {
  expect_equal(Carlson_RG(0, 16, 16), as.complex(pi))
  expect_equal(Carlson_RG(2, 3, 4), as.complex(1.7255030280692))
  expect_equal(Carlson_RG(0, 1i, -1i), as.complex(0.42360654239699))
  expect_equal(Carlson_RG(0, 0.0796, 4), as.complex(1.0284758090288))
  expect_equal(Carlson_RG(-1+1i, 1i, 0), 0.44660591677018 + 0.70768352357515i)
  expect_equal(Carlson_RG(-1i, -1+1i, 1i), 0.36023392184473 + 0.40348623401722i)
})
