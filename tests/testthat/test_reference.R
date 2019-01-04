context("reference")

test_that("reference data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
  bw <- prepare_spec_fit(pjnz, proj.end=2021.5)
  ref <- readRDS("reference/Botswana2018.rds")
  expect_equal(bw, ref)
})