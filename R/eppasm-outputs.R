#' Get vector of years spanned by projections
get_proj_years <- function(ss){
  ss$proj_start + 1:ss$PROJ_YEARS - 1L
}


#' Add dimnames to spec model output
#'
#' `mod_dimnames` assigns dimnames to spec model outputs.
#'
#' @param mod EPP-ASM model output of class `spec`
#' @param ss Model state space inputs
#'
#' @return EPP-ASM model output with labelled dimensions
#'
#' @examples
#' pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
#' fp <- prepare_directincid(pjnz)
#' mod <- simmod(fp)
#'
#' @export
mod_dimnames <- function(mod, ss){

  yrlbl <- get_proj_years(ss)
  sexlbl <- c("male", "female")
  hAGlbl <- c("15-16", "17-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50+")
  pAGlbl <- 15:80
  pDSlbl <- c("hivn", "hivp")
  hDSlbl <- c(">500", "350-499", "250-349", "200-249", "100-199", "50-99", "<50")
  hTSlbl <- c("art0mos", "art6mos", "art1yr")

  dimnames(mod) <- list(age = pAGlbl, sex = sexlbl, hiv = pDSlbl, year = yrlbl)
  dimnames(attr(mod, "hivpop")) <- list(cd4stage = hDSlbl, agegr = hAGlbl, sex = sexlbl, year = yrlbl)
  dimnames(attr(mod, "artpop")) <- list(artdur = hTSlbl, cd4stage = hDSlbl, agegr = hAGlbl, sex = sexlbl, year = yrlbl)

  dimnames(attr(mod, "infections")) <- dimnames(mod)
  dimnames(attr(mod, "hivdeaths")) <- dimnames(mod)
  dimnames(attr(mod, "natdeaths")) <- dimnames(mod)
  dimnames(attr(mod, "prev15to49")) <- yrlbl
  dimnames(attr(mod, "pregprev")) <- yrlbl
  dimnames(attr(mod, "incid15to49")) <- yrlbl

  mod
}

#' Convert aggregate HIV population to single year
#'
#' EPP-ASM tracks the CD4 distribution and ART duration
#' of the HIV population by coarse age groups `15-16, 17-19,
#' 20-24, ..., 45-49, 50+` for computational efficiency.
#'
#' `hivpop_singleage` converts the coarse age group CD4
#' distribution and ART coverage to single year of age counts
#' assuming uniform proportions in each category within coarse
#' age groups.
#'
#' @inheritParams mod_dimnames
#'
#' @examples
#' pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
#' fp <- prepare_directincid(pjnz)
#' mod <- simmod(fp)
#' hivp1 <- hivpop_singleage(mod, fp$ss)
#'
#' @export
hivpop_singleage <- function(mod, ss){

  hivp <- mod[ , , ss$hivp.idx, ]
  hivpop <- attr(mod, "hivpop")
  artpop <- attr(mod, "artpop")

  denom <- colSums(hivpop) + colSums(artpop,,2)
  hivdist <- sweep(hivpop, 2:4, denom, "/")
  artdist <- sweep(artpop, 3:5, denom, "/")

  hivdist[is.na(hivdist)] <- 0
  artdist[is.na(artdist)] <- 0

  hivpop1 <- sweep(hivdist[,ss$ag.idx,,], 2:4, hivp, "*")
  artpop1 <- sweep(artdist[,,ss$ag.idx,,], 3:5, hivp, "*")

  dimnames(hivpop1)[2:4] <- dimnames(hivp)
  dimnames(artpop1)[3:5] <- dimnames(hivp)

  list(hivpop1 = hivpop1, artpop1 = artpop1)
}
