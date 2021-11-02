create_specfp <- function(pjnz, upd.path=NULL, hiv_steps_per_year = 10L, use_ep5=FALSE){

  ## demographic inputs
  if(!is.null(upd.path))
    demp <- read_demog_param(upd.path) # from UPD file
  else
    demp <- read_specdp_demog_param(pjnz, use_ep5=use_ep5)  # from DP file or ep5 file

  ## HIV projection parameters (most AIM adult parameters)
  projp <- read_hivproj_param(pjnz, use_ep5=use_ep5)

  ## HIV projection results
  specres <- read_hivproj_output(pjnz)

  specfp <- create_spectrum_fixpar(projp, demp, hiv_steps_per_year= hiv_steps_per_year)
  attr(specfp, "country") <- read_country(pjnz)
  attr(specfp, "country") <- read_country(pjnz)
  attr(specfp, "region") <- read_region(pjnz)

  return(specfp)
}

#' Prepare direct incidence input EPP simulation
#' 
#' @param method Method for intercalating new infections. Either "directincid_ann"
#'   for annual new infections or "directincid_hts" for new infections each
#'   HIV time step.
#' 
#' @export
prepare_directincid <- function(pjnz, method = "directincid_hts"){
  specfp <- create_specfp(pjnz)
  specfp$eppmod <- method
  specfp$incidinput <- read_incid_input(pjnz)
  specfp$incidpopage <- attr(specfp$incidinput, "incidpopage")
  return(specfp)
}
