context("test model simulation")

pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
bw <- prepare_spec_fit(pjnz, proj.end=2022.5)

bw_fp <- attr(bw$Urban, "specfp")

bw_fp <- prepare_rhybrid(bw_fp)
bw_fp$ancsitedata <- TRUE
bw_fp$ancrt <- "census"
bw_fp$ancrtsite.beta <- 0
bw_fp$logitiota <- TRUE

bw_theta <- c(-0.407503322169364, -2.76794181367538, -1.26018073624346, 1995.96447776502,
              -0.00307437171215574, 0.0114118307148102, 0.00760958379603691, 0.02,
              2.24103194827232, -0.0792123921862689, -5.01917961803606, 0.359444135205712,
              -6.10051517060137)

param <- fnCreateParam(bw_theta, bw_fp)
bw_fp <- stats::update(bw_fp, list=param)

bw_prev_mod <- c(0.00045, 0.00080, 0.00140, 0.00245, 0.00424, 0.00725, 0.01214,
                 0.01985, 0.03147, 0.04804, 0.07013, 0.09750, 0.12857, 0.16083,
                 0.19143, 0.21797, 0.23908, 0.25419, 0.26363, 0.26817, 0.26870,
                 0.26618, 0.26127, 0.25506, 0.24881, 0.24356, 0.23917, 0.23534,
                 0.23175, 0.22856, 0.22561, 0.22270, 0.21927, 0.21553, 0.21148,
                 0.20693, 0.20162, 0.19568, 0.18942, 0.18303, 0.17645, 0.16968, 0.16281)

test_that("model simulation returns correct prevalence", {
  expect_equal(round(prev(simmod(bw_fp))[11:53], 5), bw_prev_mod)
  expect_equal(round(prev(simmod(bw_fp, VERSION="R"))[11:53], 5), bw_prev_mod)
})


pjnz <- system.file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ", package="eppasm")
mpm <- prepare_spec_fit(pjnz, proj.end=2021.5)

mp_fp <- attr(mpm[[1]], "specfp")

mp_fp <- prepare_rhybrid(mp_fp)
mp_fp$ancsitedata <- TRUE
mp_fp$ancrt <- "census"
mp_fp$ancrtsite.beta <- 0
mp_fp$logitiota <- TRUE

theta <- c(-0.407503322169364, -2.76794181367538, -1.26018073624346, 1995.96447776502,
           -0.00307437171215574, 0.0114118307148102, 0.00760958379603691, 0.02,
           2.24103194827232, -0.0792123921862689, -5.01917961803606, 0.359444135205712,
           -6.10051517060137)

param <- fnCreateParam(theta, mp_fp)
mp_fp <- stats::update(mp_fp, list=param)

mp_prev_mod <- c(0.00049, 0.00087, 0.00154, 0.00271, 0.00468, 0.00792, 0.01299,
                 0.02070, 0.03197, 0.04795, 0.06960, 0.09710, 0.12906, 0.16258,
                 0.19452, 0.22206, 0.24400, 0.25986, 0.26997, 0.27499, 0.27559,
                 0.27251, 0.26643, 0.25806, 0.24832, 0.23798, 0.22817, 0.22001,
                 0.21364, 0.20833, 0.20326, 0.19836, 0.19341, 0.18839, 0.18237,
                 0.17626, 0.16982, 0.16327, 0.15683, 0.15033, 0.14387, 0.13679)

test_that("Mozambique Maputo Cidade returns correct prevalence", {
  expect_equal(round(prev(simmod(mp_fp))[11:52], 5), mp_prev_mod)
  expect_equal(round(prev(simmod(mp_fp, VERSION="R"))[11:52], 5), mp_prev_mod)
})

nl_fp <- prepare_directincid(system.file("extdata/testpjnz", "Netherlands2017.PJNZ", package="eppasm"))

hivpop_mod <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 38.60468, 276.39310, 828.50028,
                1715.87826, 2865.02340, 4113.45910, 5269.59086, 6196.00968,
                6825.47670, 7160.41534, 7275.93101, 7281.37057, 7272.88727,
                7316.13815, 7448.20666, 7671.33306, 8021.61132, 8725.41031,
                9613.40435, 10544.19529, 11480.95796, 12414.52177, 13338.21818,
                14248.42002, 15144.56519, 16023.05551,16881.70360, 17715.36359,
                18523.11411, 19324.50718, 20123.23385, 20884.58237, 21565.11389,
                22119.90287, 22510.25697, 22696.13182, 22901.92810, 23121.63422,
                23310.46564, 23462.79874, 23581.15226, 23671.10697)

test_that("Netherlands returns correct HIV population size", {
  expect_equal(round(colSums(simmod(nl_fp)[,,2,],,2), 6), hivpop_mod)
  expect_equal(round(colSums(simmod(nl_fp, "R")[,,2,],,2), 6), hivpop_mod)
})


test_that("hivpop1 and artpop1 align with hivpop and artpop", {
  mod <- simmod(bw_fp)
  hp1 <- hivpop_singleage(mod, bw_fp$ss)
  expect_equal(colSums(hp1$hivpop1,,2), colSums(attr(mod, "hivpop"),,2))
  expect_equal(colSums(hp1$artpop1,,3), colSums(attr(mod, "artpop"),,3))
})


test_that("pop and hivpop+artpop are synchronised", {
  mod <- simmod(bw_fp)
  modR <- simmod(bw_fp, "R")
  expect_equal(colSums(mod[,,2,],,2),
               colSums(attr(mod, "hivpop"),,3) + colSums(attr(mod, "artpop"),,4))
  expect_equal(colSums(modR[,,2,],,2),
               colSums(attr(modR, "hivpop"),,3) + colSums(attr(modR, "artpop"),,4))
})
