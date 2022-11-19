context("test log likelihood")

pjnz <- system.file("extdata", "testpjnz", "Botswana2018.PJNZ", package = "eppasm")
inputs <- prepare_spec_fit(pjnz, 2022.5)
fp <- attr(inputs$Urban, "specfp")
likdat <- prepare_likdat(attr(inputs$Urban, "eppd"), fp)

fp$ancsitedata <- TRUE
fp$ancrt <- "census"
fp$ancrtsite.beta <- 0
fp$logitiota <- TRUE


fp_rhybrid <- prepare_rhybrid(fp)

theta_rhybrid <- c(-0.407503322169364, -2.76794181367538, -1.26018073624346, 1995.96447776502, 
                   -0.00307437171215574, 0.0114118307148102, 0.00760958379603691, 0.02,
                   2.24103194827232, -0.0792123921862689, -5.01917961803606, 0.359444135205712, 
                   -6.10051517060137)

fp_rspline <- prepare_rspline_model(fp)
theta_rspline <- c(0.753890895839314, 0.700227807059524, -1.92682607738255, 1.1257626772613, 
                   0.000243438714407202, 0.461040753192746, -0.393015728238165, 
                   1.97574135543477, 0.24367501614327, -3.80194029103096, -0.28992227570734, 
                   -4.03773908970723)

fp_rspline_eq <- update(fp_rspline, equil.rprior = TRUE)

fp_rtrend <- prepare_rtrend_model(fp)
theta_rtrend <- c(1977.68776203971, 16.2122044236609, 0.196472777334891, 0.455548771721825, 
                  0.0922592759753597, -0.162627000247779, -0.0190902657376532, 
                  -0.0267911146079137, -4.95764130472059, 0.280707027212794, -6.70883228710422)

test_that("ll returns expected value", {
  expect_equal(round(ll(theta_rhybrid, fp_rhybrid, likdat), 4),
               c(anc = 81.3049, ancrt = 7.4092, hhs = 6.816, incid = 0, artcov = 0.0, rprior = 0))
  expect_equal(round(ll(theta_rspline, fp_rspline, likdat), 4),
               c(anc = 52.9655, ancrt = 2.0052, hhs = 3.9813, incid = 0, artcov = 0.0, rprior = 0))
  expect_equal(round(ll(theta_rspline, fp_rspline_eq, likdat), 4),
               c(anc = 52.9655, ancrt = 2.0052, hhs = 3.9813, incid = 0, artcov = 0.0, rprior = -15.1136))
  expect_equal(round(ll(theta_rtrend, fp_rtrend, likdat), 4),
               c(anc = 79.1005, ancrt = 8.6278, hhs = 6.7833, incid = 0, artcov = 0.0, rprior = 0))
})


test_that("lprior returns expected value", {
  expect_equal(round(lprior(theta_rhybrid, fp_rhybrid), 5), -4.84719)
  expect_equal(round(lprior(theta_rspline, fp_rspline), 5), -19.90811)
  expect_equal(round(lprior(theta_rspline, fp_rspline_eq), 5), -19.90811)
  expect_equal(round(lprior(theta_rtrend, fp_rtrend), 5), -9.35207)
})
