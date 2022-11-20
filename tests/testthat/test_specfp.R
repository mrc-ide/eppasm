context("test specfp objects")

pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
fp <- prepare_directincid(pjnz)

test_that("paedsurv dist sums to one for each sex-year", {
  expect_true(all(round(colSums(fp$paedsurv_cd4dist), 5) %in% 0:1))
  expect_true(all(round(colSums(fp$paedsurv_artcd4dist,,2), 5) %in% 0:1))
})

test_that("PJNZ read functions work for Francophone saved files", {

  ## In Spectrum files saved on Francophone computers, comma is used instead
  ## of period for version number (and maybe other things?). This is hard to
  ## test because it requires file saved on Francophone computer.

  ## !!! Would like to seek permission from UNAIDS before putting test file on Github
  ##     For now, only running test locally.
  skip_on_ci()

  tgo2022_pjnz <- "~/Data/Spectrum files/2022 final shared/EPP-Gen/togo_2022_final.pjnz"
  tgo2022_demp <- read_specdp_demog_param(tgo2022_pjnz)
  tgo2022_projp <- read_hivproj_param(tgo2022_pjnz)
  tgo2022_specres <- read_hivproj_output(tgo2022_pjnz)

  expect_is(tgo2022_projp, "projp")
  expect_equal(tgo2022_projp$spectrum_version, "6.13")

  expect_is(tgo2022_demp, "demp")
  expect_is(tgo2022_specres, "specres")
  
})
