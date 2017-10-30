calc_entrantprev <- function(specres){
  ageprev1 <- function(specres, str){colSums(specres$hivpop[str,,,drop=FALSE],,2) / colSums(specres$totpop[str,,,drop=FALSE],,2)}
  prev14 <- ageprev1(specres, "14")
  prev13 <- ageprev1(specres, "13")
  Sx <- c(0, prev14[-1] / prev13[-length(prev13)])  # assume that survival from age 14 to 15 is the same as 13 to 14
  Sx <- ifelse(is.na(Sx), 0, Sx)
  return(Sx*prev14)
}

calc_entrantartcov <- function(specres){
  ## Assume uniform ART coverage among age 10-14
  artcov <- with(specres, ((artnum.f+artnum.m)/(hivnum.f+hivnum.m))["10-15",])
  artcov <- ifelse(is.nan(artcov), 0, artcov)
  artcov
}

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
  specfp$entrantprev <- calc_entrantprev(specres)
  specfp$entrantartcov <- calc_entrantartcov(specres)

  attr(specfp, "country") <- read_country(pjnz)
  attr(specfp, "country") <- read_country(pjnz)
  attr(specfp, "region") <- read_region(pjnz)

  return(specfp)
}

prepare_directincid <- function(pjnz){
  specfp <- create_specfp(pjnz)
  specfp$eppmod <- "directincid"
  specfp$incidinput <- read_incid_input(pjnz)
  specfp$incidpopage <- attr(specfp$incidinput, "incidpopage")
  return(specfp)
}
