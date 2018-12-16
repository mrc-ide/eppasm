context("test model simulation")

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


pjnz <- system.file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ", package="eppasm")
mpm <- prepare_spec_fit(pjnz, proj.end=2021.5)

fp <- attr(mpm[[1]], "specfp")

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

prev_mod <- c(0.00049, 0.00087, 0.00154, 0.00271, 0.00468, 0.00792, 0.01299,
              0.0207, 0.03198, 0.04796, 0.0696, 0.09711, 0.12905, 0.16256,
              0.19448, 0.22199, 0.24389, 0.2597, 0.26977, 0.27482, 0.27557,
              0.27277, 0.26702, 0.25891, 0.24925, 0.23878, 0.22872, 0.22031,
              0.21377, 0.20837, 0.20329, 0.1984, 0.19351, 0.18852, 0.18249,
              0.17639, 0.16994, 0.16339, 0.15694, 0.15043, 0.14395, 0.13687)

test_that("Mozambique Maputo Cidade returns correct prevalence", {
  expect_equal(round(prev(simmod(fp))[11:52], 5), prev_mod)
  expect_equal(round(prev(simmod(fp, VERSION="R"))[11:52], 5), prev_mod)
})

nl_fp <- prepare_directincid(system.file("extdata/testpjnz", "Netherlands2017.PJNZ", package="eppasm"))

hivpop_mod <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 38.742753, 277.452295, 832.078741,
                1724.498848, 2882.220122, 4143.634016, 5317.562667, 6266.262444,
                6921.247119, 7283.072597, 7424.889773, 7454.23537, 7466.1148,
                7525.781945, 7670.608489, 7903.543258, 8261.445342, 8971.686377,
                9865.40954, 10802.133475, 11746.157034, 12687.958515, 13620.69696,
                14540.563554, 15446.829597, 16335.870383, 17206.562392, 18053.625869,
                18874.585921, 19687.021709, 20494.370225, 21262.739833, 21949.177501,
                22508.425513, 22901.421201, 23087.367497, 23293.557646, 23520.223741,
                23721.93429, 23886.531292, 24014.58203, 24110.974724)

test_that("Netherlands returns correct HIV population size", {
  expect_equal(round(colSums(simmod(nl_fp)[,,2,],,2), 6), hivpop_mod)
  expect_equal(round(colSums(simmod(nl_fp, "R")[,,2,],,2), 6), hivpop_mod)
})
