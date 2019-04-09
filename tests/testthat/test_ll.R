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
  expect_equal(round(sum(ll(theta_rhybrid, fp_rhybrid, likdat)), 4), 95.5968)
  expect_equal(round(ll(theta_rspline, fp_rspline, likdat), 4),
               c(anc = 53.0482, ancrt = 1.6649, hhs = 4.4887,
                 incid = 0, sibmx = 0, rprior = 0, incpen = 0))
  expect_equal(round(ll(theta_rspline, fp_rspline_eq, likdat), 4),
               c(anc = 53.0482, ancrt = 1.6649, hhs = 4.4887, incid = 0, 
                 sibmx = 0, rprior = -15.0062, incpen = 0))
  expect_equal(round(ll(theta_rtrend, fp_rtrend, likdat), 4),
               c(anc = 79.2002, ancrt = 8.6026, hhs = 6.8623, incid = 0,
                 sibmx = 0, rprior = 0, incpen = 0))


test_that("lprior returns expected value", {
  expect_equal(round(lprior(theta_rhybrid, fp_rhybrid), 5), -4.84719)
  expect_equal(round(lprior(theta_rspline, fp_rspline), 5), -19.90811)
  expect_equal(round(lprior(theta_rspline, fp_rspline_eq), 5), -19.90811)
  expect_equal(round(lprior(theta_rtrend, fp_rtrend), 5), -9.35207)
})



set.seed(153)

sample_rhybrid <- structure(c(0.35973, -1.27479, -2.5652, -2.50368, -1.01206, -2.26935, 
                              1984.2949, 1981.56757, -0.04735, -0.04877, -0.04364, 0.05991, 
                              0.01654, 0.0562, -0.03922, -0.09532, -0.98921, 2.9864, 1.5763, 
                              0.53065, -3.70676, -5.97844, 1.04827, -0.09975, -3.62278, -5.06693
                              ), .Dim = c(2L, 13L))

sample_rspline <- structure(c(1.57187, 3.31483, 0.13686, -1.21512, -0.03201, -1.10401, 
                              -0.00295, 0.28975, 0.73763, 1.1998, -0.81011, -1.45164, 0.04162, 
                              1.59099, 0.46717, 5.27125, 0.34364, -1.60021, -3.59165, -6.05765, 
                              0.60145, -0.81516, -3.31116, -4.93377), .Dim = c(2L, 12L))

sample_rtrend <- structure(c(1977.90597, 1979.45429, 19.34604, 25.04682, 0.25448, 
                             0.5784, 0.62522, 0.46922, 0.12895, 0.15198, -0.85657, -0.2859, 
                             -0.04302, -0.04224, -0.32472, -0.02356, -3.33461, -5.6934, 0.47346, 
                             -0.71488, -3.33361, -5.04997), .Dim = c(2L, 11L))


test_that("sample.prior returns expected value", {
  expect_equal(round(sample.prior(2, fp_rhybrid), 5), sample_rhybrid)
  expect_equal(round(sample.prior(2, fp_rspline), 5), sample_rspline)
  expect_equal(round(sample.prior(2, fp_rtrend), 5), sample_rtrend)
})
  
