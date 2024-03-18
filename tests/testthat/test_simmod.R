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
                 0.19143, 0.21797, 0.23908, 0.25419, 0.26363, 0.26817, 0.26871,
                 0.26618, 0.26129, 0.25510, 0.24887, 0.24365, 0.23928, 0.23545,
                 0.23186, 0.22866, 0.22569, 0.22277, 0.21931, 0.21555, 0.21147,
                 0.20690, 0.20155, 0.19559, 0.18930, 0.18290, 0.17630, 0.16953, 0.16266)

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
                 0.0207, 0.03197, 0.04795, 0.0696, 0.0971, 0.12906, 0.16258, 0.19452,
                 0.22206, 0.244, 0.25986, 0.26997, 0.27499, 0.27559, 0.27251,
                 0.26643, 0.25806, 0.24833, 0.238, 0.22818, 0.22003, 0.21365,
                 0.20834, 0.20327, 0.19836, 0.19341, 0.18839, 0.18237, 0.17626,
                 0.16981, 0.16327, 0.15682, 0.15032, 0.14386, 0.13678)

test_that("Mozambique Maputo Cidade returns correct prevalence", {
  expect_equal(round(prev(simmod(mp_fp))[11:52], 5), mp_prev_mod)
  expect_equal(round(prev(simmod(mp_fp, VERSION="R"))[11:52], 5), mp_prev_mod)
})

nl_fp <- prepare_directincid(system.file("extdata/testpjnz", "Netherlands2017.PJNZ", package="eppasm"))

hivpop_mod <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 38.604677, 276.393101, 828.500283,
                1715.878263, 2865.023401, 4113.459099, 5269.59086, 6196.009677,
                6825.476696, 7160.415338, 7275.931005, 7281.370565, 7272.887268,
                7316.138154, 7448.206665, 7671.33306, 8021.611324, 8725.410313,
                9613.404346, 10544.195289, 11480.957958, 12414.52177, 13338.218178,
                14248.420019, 15144.565187, 16023.055512, 16881.703605, 17715.363593,
                18523.114109, 19324.50718, 20123.233854, 20884.582369, 21565.113893,
                22119.902865, 22510.256973, 22696.131819, 22902.414595, 23122.710087,
                23312.162274, 23465.130299, 23584.101408, 23674.628092)

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
