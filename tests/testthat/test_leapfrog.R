theta_rhybrid <- c(-0.407503322169364, -2.76794181367538, -1.26018073624346, 1995.96447776502,
                   -0.00307437171215574, 0.0114118307148102, 0.00760958379603691, 0.02,
                   2.24103194827232, -0.0792123921862689, -5.01917961803606, 0.359444135205712,
                   -6.10051517060137)

expect_mod_equivalent <- function(mod, mod_expected) {
  expect_equal(class(mod), class(mod_expected))
  expect_equal(dim(mod), dim(mod_expected))
  attribute_names <- setdiff(names(attributes(mod_expected)), c("dim", "class"))
  known_missing_attributes <- c("natdeaths_noart", "natdeaths_art", "popadjust",
                                "pregprevlag", "entrantprev")
  known_different_attributes <- c("rvec_ts")
  attribuites_to_check <- setdiff(
    attribute_names,
    c(known_missing_attributes, known_different_attributes))

  for (attr_name in attribuites_to_check) {
    expected <- attr(mod_expected, attr_name)
    received <- attr(mod, attr_name)
    if (is.array("expected")) {
      expect_equal(dim(received), dim(expected), info = attr_name)
    } else {
      expect_equal(length(received), length(expected), info = attr_name)
    }
  }
}

test_that("can run convert transmission input fp into leapfrog params", {
  pjnz <- system.file("extdata", "testpjnz", "Botswana2018.PJNZ", package = "eppasm")
  inputs <- prepare_spec_fit(pjnz, 2025.5)

  fp <- prep_fp_fitmod(inputs$Urban, eppmod = "rhybrid")
  fp <- stats::update(fp, list=fnCreateParam(theta_rhybrid, fp))

  params <- fp_to_leapfrog_params(fp)

  expect_equal(params$projection_start_year, 1970)
  expect_equal(params$sim_years, 48)
  expect_equal(params$incidence_model_choice, 1L)
  expect_true(length(params$transmission_rate_hts) > 0)
  expect_true(any(params$transmission_rate_hts > 0))
  expect_equal(class(params), "specfp")
})

test_that("can convert output from leapfrog into mod data", {
  pjnz <- system.file("extdata", "testpjnz", "Botswana2018.PJNZ", package="eppasm")
  params <- leapfrog::prepare_leapfrog_parameters(pjnz, use_coarse_age_groups = TRUE)

  n_years <- dim(params$Sx)[length(dim(params$Sx))]
  output_years <- seq(params$projection_start_year, params$projection_start_year + n_years - 1)
  leapfrog_result <- leapfrog::run_model(params, "HivCoarseAgeStratification", output_years = output_years)

  params$sim_years <- length(output_years)
  params$proj_years <- params$sim_years
  mod <- leapfrog_result_to_mod(leapfrog_result, params)

  fp <- prepare_directincid(pjnz)
  mod_expected <- simmod(fp, VERSION = "C")

  expect_mod_equivalent(mod, mod_expected)
})

test_that("can run simmod with leapfrog direct incidence data", {
  pjnz <- system.file("extdata", "testpjnz", "Botswana2018.PJNZ", package = "eppasm")
  fp <- prepare_directincid(pjnz)

  mod <- simmod(fp, VERSION = "leapfrog")
  mod_expected <- simmod(fp, version = "C")

  expect_mod_equivalent(mod, mod_expected)
})

test_that("can run simmod with leapfrog transmission input data", {
  pjnz <- system.file("extdata", "testpjnz", "Botswana2018.PJNZ", package = "eppasm")
  inputs <- prepare_spec_fit(pjnz, 2018.5)

  fp <- prep_fp_fitmod(inputs$Urban, eppmod = "rhybrid")
  fp <- stats::update(fp, list=fnCreateParam(theta_rhybrid, fp))

  mod <- simmod(fp, VERSION = "leapfrog")
  mod_expected <- simmod(fp, VERSION = "C")

  expect_mod_equivalent(mod, mod_expected)

  # Assert outputs are not all 0
  expect_true(any(attr(mod, "hivpop") > 0))
  expect_true(any(attr(mod, "artpop") > 0))
  expect_true(any(attr(mod, "infections") > 0))

  ## Assert padded to correct lengths with 0s
  hivpop <- attr(mod, "hivpop")
  expect_true(all(hivpop[,,,49] == 0))

  hivpop <- attr(mod_expected, "hivpop")
  expect_true(all(hivpop[,,,49] == 0))

  incrate15to49_ts <- attr(mod, "incrate15to49_ts")
  expect_true(all(incrate15to49_ts[471:480] == 0))

  incrate15to49_ts <- attr(mod_expected, "incrate15to49_ts")
  expect_true(all(incrate15to49_ts[471:480] == 0))
})


test_that("pad dim util works", {
  x <- array(1:6, dim = c(2, 3))
  padded <- pad_last_dim(x, 4)

  expect_equal(dim(padded), c(2, 4))
  expect_equal(padded[,1:3], x)
  expect_true(all(padded[,4] == 0))

  x <- array(1:24, dim = c(2, 3, 4))
  padded <- pad_last_dim(x, 5)

  expect_equal(dim(padded), c(2, 3, 5))
  expect_equal(padded[,,1:4], x)
  expect_true(all(padded[,,5] == 0))

  # When no padding, nothing changes
  x <- array(1:6, dim = c(2, 3))
  padded <- pad_last_dim(x, 3)
  expect_equal(x, padded)
  padded <- pad_last_dim(x, 2)
  expect_equal(x, padded)

  # Works with a vector
  x <- c(1, 2, 3, 4, 5)
  padded <- pad_last_dim(x, 10)
  expect_equal(padded, c(1, 2, 3, 4, 5, 0, 0, 0, 0, 0))

  # Works with matrix
  x <- matrix(1:6, nrow = 3, ncol = 2)
  padded <- pad_last_dim(x, 5)
  expect_equal(dim(padded), c(3, 5))
  expect_equal(padded[,1:2], x)
  expect_true(all(padded[,3:5] == 0))
})
