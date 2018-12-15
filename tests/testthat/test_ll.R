context("test log likelihood")

pjnz <- system.file("extdata", "testpjnz", "Botswana2018.PJNZ", package = "eppasm")
inputs <- prepare_spec_fit(pjnz, 2022.5)
fp <- attr(inputs$Urban, "specfp")
likdat <- prepare_likdat(attr(inputs$Urban, "eppd"), fp)

fp <- prepare_rhybrid(fp)
fp$ancsitedata <- TRUE
fp$ancrt <- "census"
fp$ancrtsite.beta <- 0
fp$logitiota <- TRUE

theta <- c(-0.407503322169364, -2.76794181367538, -1.26018073624346, 1995.96447776502, 
           -0.00307437171215574, 0.0114118307148102, 0.00760958379603691, 0.02,
           2.24103194827232, -0.0792123921862689, -5.01917961803606, 0.359444135205712, 
           -6.10051517060137)

test_that("ll returns expected value", {
  expect_equal(round(sum(ll(theta, fp, likdat)), 4), 95.5531)
})
