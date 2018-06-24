context("test model simulation")

devtools::load_all()
load("ll-test-data.rda")

## pjnz <- system.file("extdata/testpjnz", "Botswana2017.PJNZ", package="eppasm")
## bw <- prepare_spec_fit(pjnz, proj.end=2021.5)

## ## r-spline model: fixed parameter values
## theta.rspline <- c(2.16003605, -0.76713859, 0.21682066, 0.03286402, 0.21494412,
##                    0.40138627, -0.08235464, -16.32721684, -2.97511957, 0.21625028, -3.5)


## fp <- attr(bw$Urban, "specfp")
## fp$ancsitedata <- FALSE


param <- fnCreateParam(theta, fp)
fp <- update(fp, list=param)

modC <- simmod(fp)
modR <- simmod(fp, VERSION="R")

prev_mod <- c(0.00252, 0.00429, 0.00716, 0.01169, 0.0186, 0.0287, 0.04262, 
              0.06079, 0.083, 0.10831, 0.13514, 0.16156, 0.18569, 0.20613, 
              0.22201, 0.23315, 0.23966, 0.24192, 0.24045, 0.23593, 0.22903, 
              0.22037, 0.21035, 0.1997, 0.18862, 0.17764, 0.16703, 0.1577, 
              0.14932, 0.14174, 0.13507, 0.12919, 0.12348, 0.11802, 0.11287, 
              0.10786, 0.10301)

undiagnosed_mod <- c(0.892418, 0.888746, 0.885385, 0.881648, 0.877092, 0.871361, 
                     0.864213, 0.855295, 0.844351, 0.831198, 0.815767, 0.798129, 0.77851, 
                     0.75728, 0.734922, 0.711958, 0.68893, 0.666345, 0.644636, 0.624126, 
                     0.605082, 0.587611, 0.571836, 0.557781, 0.545372, 0.534135, 0.52241, 
                     0.507352, 0.489838, 0.469897, 0.444872, 0.415332, 0.384634, 0.354937, 
                     0.324639, 0.273219, 0.206702)

test_that("model simulation returns correct prevalence", {
  expect_equal(round(prev(modC)[11:47], 5), prev_mod)
  expect_equal(round(prev(modR)[11:47], 5), prev_mod)
})

test_that("model simulation returns correct proportion undiagnosed", {
  expect_equal(round(calc_undiagnosed(modC, fp)[11:47], 6), undiagnosed_mod)
  expect_equal(round(calc_undiagnosed(modR, fp)[11:47], 6), undiagnosed_mod)
})
