context("Test input checking")

test_that("Mis-matched lengths give errors",{
  expect_error(dsr(5:1, 1:5,1:4))
  expect_error(dsr(5:1, 1:4,1:5))
  expect_error(dsr(5:1, 1:4,1:5))
 })

test_that("NA gives error",{
  xNA <- c(NA, 1:4)
  expect_error(dsr(xNA, 1:5,1:5))
  expect_error(dsr(5:1, xNA,1:5))
  expect_error(dsr(5:1, 1:5,xNA))
})

test_that("Level must be [0.5,1)",{
  expect_error(dsr(1:5,5:1,c(1,3,5,7,9), level = 0.2))
  expect_error(dsr(1:5,5:1,c(1,3,5,7,9), level = 1))
})