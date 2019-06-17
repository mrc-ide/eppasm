#' @useDynLib eppasm eppasmC eppasmOOpp
simmod.specfp <- function(fp, VERSION="C", MODEL=1, MIX=FALSE) {
  if (!exists("popadjust", where=fp))
    fp$popadjust <- FALSE

  if (!exists("incidmod", where=fp))
    fp$incidmod <- "eppspectrum"

  if (MODEL==2) 
    fp <- update_fp_debut(fp, max_debut_age=30)
  
  if (MIX && !exists("mat_f", where=fp)) { # add these outside
    fp$mat_m <- readRDS(system.file("extdata", "contact_matrix_male.rds",
                                    package="eppasm"))
    fp$mat_f <- readRDS(system.file("extdata", "contact_matrix_female.rds",
                                    package="eppasm"))
  }

  if (VERSION != "R") {
    if (VERSION=="K") { # C++ classes
      fp <- prepare_fp_for_Cpp(fp)
      mod <- .Call(eppasmOOpp, fp, MODEL, MIX)
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
    pop$year <- hivpop$year <- artpop$year <- i
    epp_aging(pop, hivpop, artpop)
    epp_death(pop, hivpop, artpop)
    epp_migration(pop, hivpop, artpop)
    pop$update_fertile()
    if (MODEL != 0) { # Disease model simulation: events at dt timestep
      epp_disease_model(pop, hivpop, artpop)
      if (fp$eppmod == "directincid") ## Direct incidence input model
        pop$epp_disease_model_direct(hivpop, artpop)
    }
    if (exists("popadjust", where=fp) && fp$popadjust) { # match target pop
      pop$adjust_pop()
      if (MODEL != 0) {
        hivpop$adjust_pop(pop$adj_prob)
        if (i >= fp$tARTstart)
          artpop$adjust_pop(pop$adj_prob)
      }
    }
    if (MODEL != 0) {
      if (i + fp$ss$AGE_START <= fp$ss$PROJ_YEARS)
        pop$cal_prev_pregant(hivpop, artpop) # prevalence among pregnant women
      pop$save_prev_n_inc() # save prevalence and incidence 15 to 49
    }
  }
  if (MODEL != 0) {
    attr(pop, "hivpop") <- hivpop$data
    attr(pop, "artpop") <- artpop$data
    attr(pop, "debut_pop") <- pop$data_db
    class(pop) <- "spec"
  }
  else 
    class(pop) <- "dempp"
  return(pop)
}