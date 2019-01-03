context("test specfp objects")

pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
fp <- prepare_directincid(pjnz)

test_that("paedsurv dist sums to one for each sex-year", {
  expect_true(all(round(colSums(fp$paedsurv_cd4dist), 5) %in% 0:1))
  expect_true(all(round(colSums(fp$paedsurv_artcd4dist,,2), 5) %in% 0:1))
         colSums(fp$paedsurv_artcd4dist,,2) %in% 0:1
})
