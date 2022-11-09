
#' Simulate model
#' 
#' @export
simmod <- function(fp, ...) UseMethod("simmod")
simfit <- function(fit, ...) UseMethod("simfit")

#' @export
prev <- function(mod, ...) UseMethod("prev")

#' @export
fnPregPrev <- function(mod, fp, ...) UseMethod("fnPregPrev")

#' Incidence rate among adults age 15-49 years
#'
#' @param mod model output
#'
#' @details
#' This returns incidence rate calculated as the number of infections during
#' the projection period divded by the number susceptible at the mid-point of
#' the projection period. This is the default incidence calculation for
#' Spectrum version >=6.2 For Spectrum versions <=6.19 incidence was calculated
#' as number of infections divided by susceptible population at the start of the
#' projection year.
#' 
#' @export
incid <- function(mod, ...) UseMethod("incid")

incid_sexratio <- function(mod, ...) UseMethod("incid_sexratio")

agemx <- function(mod, ...) UseMethod("agemx")
natagemx <- function(mod, ...) UseMethod("natagemx")
calc_nqx <- function(mod, ...) UseMethod("calc_nqx")

pop15to49 <- function(mod, ...) UseMethod("pop15to49")
artpop15to49 <- function(mod, ...) UseMethod("artpop15to49")
artpop15plus <- function(mod, ...) UseMethod("artpop15plus")

#' @export
artcov15to49 <- function(mod, ...) UseMethod("artcov15to49")

#' @export
artcov15plus <- function(mod, ...) UseMethod("artcov15plus")

age15pop <- function(mod, ...) UseMethod("age15pop")
