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
mod <- simmod(fp)

prev_mod <- c(0.00252, 0.00429, 0.00716, 0.01169, 0.0186, 0.0287, 0.04262, 
              0.06079, 0.083, 0.10831, 0.13514, 0.16156, 0.18569, 0.20613, 
              0.22201, 0.23315, 0.23966, 0.24192, 0.24045, 0.23593, 0.22903, 
              0.22037, 0.21035, 0.1997, 0.18862, 0.17764, 0.16703, 0.1577, 
              0.14932, 0.14174, 0.13507, 0.12919, 0.12348, 0.11802, 0.11287, 
              0.10786, 0.10301)

test_that("model simulation returns correct prevalence", {
  expect_equal(round(prev(mod)[11:47], 5), prev_mod)
})
