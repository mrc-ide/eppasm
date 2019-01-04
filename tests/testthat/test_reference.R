context("reference")

test_that("Botswana2017 reference data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Botswana2017.PJNZ", package="eppasm")
  bw <- prepare_spec_fit(pjnz, proj.end=2021.5)
  ref <- readRDS("reference/Botswana2017.rds")
  expect_equal(bw, ref)
})

test_that("Botswana2018 reference data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
  bw <- prepare_spec_fit(pjnz, proj.end=2021.5)
  ref <- readRDS("reference/Botswana2018.rds")
  expect_equal(bw, ref)
})

test_that("Mozambique_Maputo_Cidade2018 reference data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ", 
                      package="eppasm")
  bw <- prepare_spec_fit(pjnz, proj.end=2021.5)
  ref <- readRDS("reference/Mozambique_Maputo_Cidade2018.rds")
  expect_equal(bw, ref)
})
