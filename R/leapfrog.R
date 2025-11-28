#' Convert eppasm fp object into leapfrog parameters
#'
#' @param fp eppasm fixed parameters
#'
#' @returns Parameters organised for running leapfrog coarse age model
#' @noRd
fp_to_leapfrog_params <- function(fp) {
  ## TODO: this is working with output from `prepare_spec_fit`
  # support direct incid input too i.e. from prepare_directincid

  births_sex_prop_male <- fp$srb / (fp$srb + 100)
  births_sex_prop <- rbind(male = births_sex_prop_male,
                           female = 1 - births_sex_prop_male)
  year_end <- fp$ss$proj_start + fp$SIM_YEARS

  if (fp$eppmod %in% c("rspline", "logrw", "rhybrid", "rlogistic")) {
    incidence_model_choice <- 1L
    incid_input <- rep(0, fp$SIM_YEARS + 1)
    # rvec here is different to what is required by leapfrog in a few ways
    # for rhybrid it can have duplicate years when the curve switches types
    # So for a fit starting in 1970.5 until 2018.5 with SIM_YEARS of 48
    # It can e.g. have data for 1970.5, 1970.6, ..., 2002.9, 2003.0, 2003.0,
    # 2003.1, ..., 2017.3, 2017.4, 2017.5 giving us 472 data points
    # But leapfrog expects SIM_YEARS * hiv_steps_per_year inputs for
    # transmission rate. We are 8 short. Repeat the first entry 8 additional
    # times for now to get the right length
    required_length <- fp$SIM_YEARS * fp$ss$hiv_steps_per_year
    pad_num <- required_length - length(fp$rvec)
    if (pad_num < 0) {
      stop("Input transmission data is longer than expected, not supported with leapfrog.")
    }
    transmission_rate_hts <- c(rep(fp$rvec[1], pad_num), fp$rvec)
    initial_incidence <- fp$iota
    epidemic_start_hts <- fp$ss$time_epi_start
    relative_infectiousness_art <- fp$relinfectART
    pAG_INCIDPOP <- 0
  } else if (fp$eppmod %in% c("directincid_ann", "directincid_hts")) {
    incidence_model_choice <- 0L
    incid_input <- fp$incidinput
    transmission_rate_hts <- rep(0, fp$SIM_YEARS * fp$ss$hiv_steps_per_year)
    initial_incidence <- 0
    epidemic_start_hts <- 0
    relative_infectiousness_art <- 0
    pAG_15_to_49 <- 35L
    pAG_15_plus <- 66L
    pAG_INCIDPOP <- ifelse(fp$incidpopage == 0L, pAG_15_to_49, pAG_15_plus)
  } else if (fp$eppmod == "rtrend") {
    stop("eppmod 'rtrend' not supported with leapfrog, try a different eppmod")
  } else {
    stop(sprintf(
      "Cannot run with unknown eppmod '%s'. Check documentation for valid values.",
      fp$eppmod
    ))
  }

  leapfrog_params <- list(
    projection_start_year = fp$ss$proj_start,
    projection_period = fp$projection_period,
    t_ART_start = min(c(unlist(apply(fp$art15plus_num > 0, 1, which)), fp$SIM_YEARS)),
    # We can only simulate up to when we have data, so PROJ_YEARS and SIM_YEARS
    # can be different. If more proj_years then we will project forward as a
    # 2nd step. Pad the data to correct length for now.
    proj_years = fp$ss$PROJ_YEARS,
    sim_years = fp$SIM_YEARS,
    basepop = fp$basepop_allage,
    Sx = fp$Sx_allage,
    netmigr_adj = fp$netmigr_adj,
    asfr = fp$asfr,
    births_sex_prop = births_sex_prop,
    tfr = fp$tfr,
    incidence_model_choice = incidence_model_choice,
    incidinput = incid_input,
    transmission_rate_hts = transmission_rate_hts,
    initial_incidence = initial_incidence,
    epidemic_start_hts = epidemic_start_hts,
    relative_infectiousness_art = relative_infectiousness_art,
    incrr_age = fp$incrr_age,
    incrr_sex = fp$incrr_sex,
    cd4_mort = fp$cd4_mort,
    cd4_prog = fp$cd4_prog,
    cd4_initdist = fp$cd4_initdist,
    scale_cd4_mort = fp$scale_cd4_mort,
    artcd4elig_idx = fp$artcd4elig_idx,
    art_mort = fp$art_mort,
    artmx_timerr = fp$artmx_timerr,
    cd4_nonaids_excess_mort = fp$cd4_nonaids_excess_mort,
    art_nonaids_excess_mort = fp$art_nonaids_excess_mort,
    art_dropout_recover_cd4 = fp$art_dropout_recover_cd4,
    art_dropout_rate = -log(1.0 - fp$art_dropout / 100),
    art15plus_num = fp$art15plus_num,
    art15plus_isperc = fp$art15plus_isperc,
    art_alloc_mxweight = fp$art_alloc_mxweight,
    h_art_stage_dur = rep(0.5, fp$ss$hiv_steps_per_year - 1),
    pAG_INCIDPOP = pAG_INCIDPOP,
    pIDX_INCIDPOP = 15L,
    fert_rat = fp$fert_rat,
    cd4fert_rat = fp$cd4fert_rat,
    frr_art6mos = fp$frr_art6mos,
    frr_scalar = fp$frr_scalar
  )

  class(leapfrog_params) <- class(fp)
  leapfrog_params
}


#' Convert leapfrog model output into eppasm mod object
#'
#'Note that we are not outputting a few variables at the moment. They
# will need to be added in time.
# * natdeaths_noart, natdeaths_art
# * popadjust
# * pregprevlag
# * entrantprev
#'
#' @param leapfrog_result output from running leapfrog model
#'
#' @returns Parameters organised for running leapfog coarse age model
#' @noRd
leapfrog_result_to_mod <- function(leapfrog_result, leapfrog_params) {
  # Age 0 is index 1, age 1 index 2 etc.
  AGE_15_IDX <- 16L
  AGE_49_IDX <- 50L
  AGE_15TO49_IDX <- AGE_15_IDX:AGE_49_IDX
  AGE_1580PLUS_IDX <- AGE_15_IDX:81

  n_ages <- length(AGE_1580PLUS_IDX)
  n_sex <- 2

  # Populate "mod", the 3rd dimension is hiv negative, hiv +ve
  mod <- array(0, dim = c(n_ages, n_sex, 2, leapfrog_params$sim_years))
  hivneg_pop <- leapfrog_result$p_totpop - leapfrog_result$p_hivpop
  mod[,,1,] <- hivneg_pop[AGE_1580PLUS_IDX,, ]
  mod[,,2,] <- leapfrog_result$p_hivpop[AGE_1580PLUS_IDX,, ]
  mod <- pad_last_dim(mod, leapfrog_params$proj_years)

  # 15to49 prevalence
  hivpop_15to49 <- colSums(leapfrog_result$p_hivpop[AGE_15TO49_IDX,,], dims = 2)
  totpop_15to49 <- colSums(leapfrog_result$p_totpop[AGE_15TO49_IDX,,], dims = 2)
  prev15to49 <- hivpop_15to49 / totpop_15to49

  # 15to49 incidence is
  # new infections in 15 - 49 / hiv -ve pop (totpop - hivpop) in previous year
  hivneg_15to49 <- colSums(hivneg_pop[AGE_15TO49_IDX,,], dims = 2)
  infections_15to49 <- colSums(leapfrog_result$p_infections[AGE_15TO49_IDX,,], dims = 2)
  incid15to49 <- array(0, dim = length(hivneg_15to49))
  n_years <- leapfrog_params$sim_years
  incid15to49[2:n_years] <- infections_15to49[2:n_years] / hivneg_15to49[1:(n_years - 1)]

  # Prevalence of HIV among pregnant women, hiv_births / births
  pregprev <- leapfrog_result$hiv_births / leapfrog_result$births

  attr(mod, "hivpop") <- leapfrog_result$h_hivpop
  attr(mod, "artpop") <- leapfrog_result$h_artpop
  attr(mod, "infections") <- leapfrog_result$p_infections[AGE_1580PLUS_IDX,,]
  attr(mod, "hivdeaths") <- leapfrog_result$p_hiv_deaths[AGE_1580PLUS_IDX,,]
  attr(mod, "natdeaths") <- leapfrog_result$p_deaths_background_hivpop[AGE_1580PLUS_IDX,,]
  attr(mod, "excessnonaidsdeaths") <- leapfrog_result$p_deaths_excess_nonaids[AGE_1580PLUS_IDX,,]
  attr(mod, "aidsdeaths_noart") <- leapfrog_result$h_hiv_deaths_no_art
  attr(mod, "excessnonaidsdeaths_noart") <- leapfrog_result$h_deaths_excess_nonaids_no_art
  attr(mod, "aidsdeaths_art") <- leapfrog_result$h_hiv_deaths_art
  attr(mod, "excessnonaidsdeaths_art") <- leapfrog_result$h_deaths_excess_nonaids_on_art
  attr(mod, "artinit") <- leapfrog_result$h_art_initiation
  # For hiv time step outputs leapfrog includes for 0th year where
  # EPPASM only starts recording this at time 1. i.e. for a 1970 to 2030
  # projection. Leapfrog ouputs 61 (1970-2030) years for this, but
  # EPPASM only 60 (1971-2030)
  # Remove the first 1 year of HIV time steps to get the expected length
  # Magic 10 for no of hiv time steps per year
  first_year_hts <- -1:-10
  attr(mod, "incrate15to49_ts") <- as.vector(leapfrog_result$incidence_15to49_hts)[first_year_hts]
  attr(mod, "prev15to49_ts") <- as.vector(leapfrog_result$prevalence_15to49_hts)[first_year_hts]
  attr(mod, "rvec_ts") <- as.vector(leapfrog_params$transmission_rate_hts)
  attr(mod, "prev15to49") <- prev15to49
  attr(mod, "pregprev") <- pregprev
  attr(mod, "incid15to49") <- incid15to49

  # When running with transmission input i.e. not direct incidence EPPASM
  # will run the model for the years for which there is data and then project
  # forward to fill in up to PROJ_YEARS. If simulation years < projection years
  # then there will be columns of 0s
  # We pad here with 0s to match current EPPASM behaviour
  pad_to_proj_years <- c("hivpop", "artpop", "infections", "hivdeaths",
                         "natdeaths", "excessnonaidsdeaths", "aidsdeaths_noart",
                         "excessnonaidsdeaths_noart", "aidsdeaths_art",
                         "excessnonaidsdeaths_art", "artinit", "prev15to49",
                         "pregprev", "incid15to49")
  pad_to_hiv_time_steps <- c("incrate15to49_ts", "prev15to49_ts")
  for (attr_name in pad_to_proj_years) {
    attr(mod, attr_name) <- pad_last_dim(attr(mod, attr_name), leapfrog_params$proj_years)
  }
  for (attr_name in pad_to_hiv_time_steps) {
    attr(mod, attr_name) <- pad_last_dim(attr(mod, attr_name), (leapfrog_params$proj_years - 1) * 10)
  }

  class(mod) <- "spec"
  mod
}


#' For a vectoror array x, pad the last dimension to expected length
#' adding 0s.
#'
#' @param x Array or vector to pad
#' @param new_length Length to achieve
#'
#' @returns Array or vector x padded in the last dimension to length new_length
#' @noRd
pad_last_dim <- function(x, new_length) {
  is_vector <- is.vector(x)
  is_matrix <- is.matrix(x)

  if (is_vector) {
    x <- as.array(x, dim = length(x))
  }

  dims <- dim(x)
  from_length <- dims[length(dims)]

  if (new_length <= from_length) return(x)

  mat <- matrix(x, ncol = from_length)
  pad <- matrix(0, nrow = nrow(mat), ncol = new_length - from_length)

  mat2 <- cbind(mat, pad)

  # reshape back to original dims
  out <- array(mat2, dim = c(dims[-length(dims)], new_length))

  if (is_vector) {
    out <- as.vector(out)
  } else if (is_matrix) {
    out <- as.matrix(out)
  }

  out
}
