#' @useDynLib eppasm eppasmC eppasmOOpp
simmod.specfp <- function(fp) {

  # Move all model options to fp
  MODEL   = ifelse(is.null(fp$ss$MODEL), 1, fp$ss$MODEL)
  MIX     = ifelse(is.null(fp$ss$MIX), FALSE, fp$ss$MIX)
  VERSION = ifelse(is.null(fp$VERSION), 'R', fp$VERSION)

  if (!exists("popadjust", where=fp))
    fp$popadjust <- FALSE

  if (!exists("incidmod", where=fp))
    fp$incidmod <- "eppspectrum"

  if (!exists("DT", where=fp))
    fp$ss$DT <- 1 / fp$ss$hiv_steps_per_year

  if (is.null(dim(fp$artmx_timerr))) { # repicate for 3 treatment durations
      fp$artmx_timerr <- matrix(rep(fp$artmx_timerr, 3), nrow=3, byrow=TRUE)
  }

  if (MODEL == 2)
    fp <- update_fp_debut(fp, max_debut_age=30)
  
  if (MIX && !exists("mixmat", where=fp)) {
    fp$mixmat <- readRDS(system.file("extdata", "est_mixmat.rds", package="eppasm"))[[1]]
    message('using a random mixing matrix')
  }

  if (VERSION != "R") {
    if (VERSION=="K") { # C++ classes
      fp  <- prepare_fp_for_Cpp(fp, MODEL, MIX)
      mod <- .Call(eppasmOOpp, fp)
      return(mod)
    } 
    else { # keep this for tests
      fp$eppmodInt <- match(fp$eppmod, c("rtrend", "directincid"), nomatch=0) # 0: r-spline;
      fp$incidmodInt <- match(fp$incidmod, c("eppspectrum"))-1L  # -1 for 0-based indexing
      mod <- .Call(eppasmC, fp)
      class(mod) <- "spec"
      return(mod)
    }
  }
  pop     <- popEPP$new(fp, MODEL, VERSION, MIX)
  hivpop  <- hivEPP$new(fp, MODEL)
  artpop  <- artEPP$new(fp, MODEL)
  for (i in 2:fp$SIM_YEARS) {
    pop$update_year <- hivpop$year <- artpop$year <- i
    epp_aging(pop, hivpop, artpop)
    epp_death(pop, hivpop, artpop)
    epp_migration(pop, hivpop, artpop)
    pop$update_fertile()
    if (MODEL != 0) { # Disease model simulation: events at dt timestep
      epp_disease_model(pop, hivpop, artpop)
      if (fp$eppmod == "directincid") ## Direct incidence input model
        pop$epp_disease_model_direct(hivpop, artpop)
    }
    if (fp$popadjust) # match target pop
      epp_adjustpop(pop, hivpop, artpop)
    if (MODEL != 0) {
      pop$cal_prev_pregant(hivpop, artpop) # prevalence among pregnant women
      pop$save_prev_n_inc() # save prevalence and incidence 15 to 49
    }
  }
  if (MODEL != 0) {
    attr(pop, "hivpop") <- hivpop$data
    attr(pop, "artpop") <- artpop$data
    if (MODEL==2)
      attr(pop, "vpop") <- pop$VIRGIN$data
    class(pop) <- "spec"
  }
  else 
    class(pop) <- "dempp"
  return(pop)
}