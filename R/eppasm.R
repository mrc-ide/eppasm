simmod.specfp <- function(fp, VERSION="C", .MODEL=1) {

  if (!exists("popadjust", where=fp))
    fp$popadjust <- FALSE

  if (!exists("incidmod", where=fp))
    fp$incidmod <- "eppspectrum"

  if (VERSION != "R") {
    fp$eppmodInt <- match(fp$eppmod, c("rtrend", "directincid"), nomatch=0) # 0: r-spline;
    fp$incidmodInt <- match(fp$incidmod, c("eppspectrum"))-1L  # -1 for 0-based indexing
    mod <- .Call(eppasmC, fp)
    class(mod) <- "spec"
    return(mod)
  }

  ## initialize projection
  pop     <- POPClass(fp, .MODEL)
  
  if (.MODEL!=0) {
    hivpop  <- HIVClass(fp, single_year = FALSE)
    artpop  <- ARTClass(fp, single_year = FALSE)
    grad    <- HIVClass(fp)
    gradART <- ARTClass(fp)
    if (.MODEL==2) {
      hiv_db     <- HIVClass(fp, single_year = FALSE)
      art_db     <- ARTClass(fp, single_year = FALSE)
      grad_db    <- HIVClass(fp)
      gradART_db <- ARTClass(fp)
    }
    if (fp$eppmod != "directincid")
      rvec <- if (fp$eppmod == "rtrend") rep(NA, length(fp$proj.steps)) 
                else fp$rvec
    ## store last prevalence value (for r-trend model)
    prevlast <- 0
  }

  ##  Single-year population projection  ##
  for (i in 2:fp$SIM_YEARS) {
    epp_aging(.MODEL, i, fp, pop, hivpop, artpop, hiv_db, art_db)

    ## survive the population
    epp_death(.MODEL, i, fp, pop, hivpop, hiv_db, artpop, art_db)

    ## net migration
    epp_migration(.MODEL, i, fp, pop, hivpop, hiv_db, artpop, art_db)

    ## fertility
    pop$update_fertile(.MODEL, i, fp)

    ##  Disease model simulation: events at dt timestep
    if (.MODEL!=0) {
      epp_disease_model(.MODEL, i, fp, pop, hivpop, artpop, hiv_db, art_db,
                        grad, gradART, grad_db, gradART_db, rvec)

      ## Direct incidence input
      if (fp$eppmod == "directincid")
        epp_disease_model_direct(pop, hivpop, i, fp)
    }

    ## adjust population to match target population size
    if (exists("popadjust", where=fp) & fp$popadjust)
      epp_adjust_pop(.MODEL, fp, i, pop, hivpop, hiv_db, artpop, art_db)

    if (.MODEL!=0) {
      ## prevalence among pregnant women
      if (i + fp$ss$AGE_START <= fp$ss$PROJ_YEARS)
        pop$cal_prev_pregant(.MODEL, i, fp, hivpop, artpop)

      ## prevalence and incidence 15 to 49
      pop$save_prev_n_inc(.MODEL, i)
    }
  }
  if (.MODEL!=0) {
    attr(pop, "hivpop") <- hivpop
    attr(pop, "artpop") <- artpop
    class(pop) <- "spec"
  }
  if (.MODEL==0) class(pop) <- "dempp"
  return(pop)
}