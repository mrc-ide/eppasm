context("test log likelihood")

devtools::load_all()
load("ll-test-data.rda")

test_that("ll returns expected value", {
  expect_equal(round(round(ll(theta, fp, likdat), 4)), -124.8817)
})
