context("test model simulation")

load("ll-test-data.rda")

pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
bw <- prepare_spec_fit(pjnz, proj.end=2022.5)

fp <- attr(bw$Urban, "specfp")

fp <- prepare_rhybrid(fp)
fp$ancsitedata <- TRUE
fp$ancrt <- "census"
fp$ancrtsite.beta <- 0
fp$logitiota <- TRUE

theta <- c(-0.407503322169364, -2.76794181367538, -1.26018073624346, 1995.96447776502, 
           -0.00307437171215574, 0.0114118307148102, 0.00760958379603691, 0.02,
           2.24103194827232, -0.0792123921862689, -5.01917961803606, 0.359444135205712, 
           -6.10051517060137)

param <- fnCreateParam(theta, fp)
fp <- update(fp, list=param)
mod <- simmod(fp)

prev_mod <- c(0.00045, 0.00080, 0.0014, 0.00245, 0.00424, 0.00725, 0.01214, 
              0.01985, 0.03147, 0.04804, 0.07013, 0.0975, 0.12857, 0.16083, 
              0.19142, 0.21797, 0.23908, 0.25419, 0.26363, 0.26816, 0.26869, 
              0.26616, 0.26124, 0.25502, 0.24876, 0.24349, 0.23902, 0.23506, 
              0.23127, 0.2278, 0.22451, 0.22121, 0.21735, 0.21318, 0.20867, 
              0.20366, 0.19787)

test_that("model simulation returns correct prevalence", {
  expect_equal(round(prev(simmod(fp))[11:47], 5), prev_mod)
  expect_equal(round(prev(simmod(fp, VERSION="R"))[11:47], 5), prev_mod)
})
