simmod.specfp <- function(fp, mx, isMixing = FALSE, VERSION="C"){
  
  library(magrittr)

  if(!exists("popadjust", where=fp))
    fp$popadjust <- FALSE

  if(!exists("incidmod", where=fp))
    fp$incidmod <- "eppspectrum"

  if(VERSION != "R" && isMixing == FALSE){
    fp$eppmodInt <- match(fp$eppmod, c("rtrend", "directincid"), nomatch=0) # 0: r-spline;
    fp$incidmodInt <- match(fp$incidmod, c("eppspectrum"))-1L  # -1 for 0-based indexing
    mod <- .Call(eppasmC, fp)
    class(mod) <- "spec"
    return(mod)
  }

##################################################################################

  if(requireNamespace("fastmatch", quietly = TRUE))
    ctapply <- fastmatch::ctapply
  else
    ctapply <- tapply

  fp$ss$DT <- 1/fp$ss$hiv_steps_per_year

  # write more proper validation steps here
  if (is.null(fp$pi)) fp$pi <- 1 
  if (isMixing) fp$ss$NG <- fp$ss$NG * length(fp$pi) # depend on n of risk groups
  if (!isMixing) fp$pi <- 1 # write more proper validation steps here

  ## Attach state space variables
  invisible(list2env(fp$ss, environment())) # put ss variables in environment for convenience

  birthslag   <- fp$birthslag
  pregprevlag <- numeric(PROJ_YEARS)

  # initialize projection
  # -----------------------------------------------------------------------------
  if (isMixing) {
    mx$D <- Dmix(fp, mx) # age-mix, prob. pns formed between i and i'
    mx$gamma <- fp$pi
    fp   <- updateRiskGroup(fp, mx)
    FOI  <- array(0, c(pAG, NG, PROJ_YEARS)) # saving force of infection
  }

  pop    <- array(0, c(pAG, NG, pDS, PROJ_YEARS))
  hivpop <- array(0, c(hDS, hAG, NG, PROJ_YEARS))
  artpop <- array(0, c(hTS, hDS, hAG, NG, PROJ_YEARS))
  
  pop[,,hivn.idx,1] <- fp$basepop # splitted pop by risk groups if isMixing

  # initialize output
  # -----------------------------------------------------------------------------
  prev15to49 <- incid15to49 <- entrant_prev_out <- paedsurvout <- numeric(PROJ_YEARS)
  sexinc15to49out   <- array(NA, c(NG, PROJ_YEARS))
  infections        <- array(0, c(pAG, NG, PROJ_YEARS))
  hivdeaths         <- natdeaths <- popadj.prob <- infections
  hivp_entrants_out <- array(NA, c(NG, PROJ_YEARS))
  prevlast          <- 0  # store last prevalence value (for r-trend model)

  if(fp$eppmod != "directincid"){
    ## outputs by timestep
    incrate15to49.ts.out <- rep(NA, length(fp$rvec))
    rvec <- if(fp$eppmod == "rtrend") rep(NA, length(fp$proj.steps)) else fp$rvec
    prev15to49.ts.out <- rep(NA, length(fp$rvec))
  }
  
  for(i in 2:fp$SIM_YEARS){

    ## ################################### ##
    ##  Single-year population projection  ##
    ## ################################### ##

    # browser()

    ## age the population
    pop[-c(1,pAG),,,i] <- pop[-(pAG-1:0),,,i-1]
    pop[pAG,,,i] <- pop[pAG,,,i-1] + pop[pAG-1,,,i-1] # open age group

    ## Add lagged births into youngest age group
    if(exists("popadjust", where = fp) & fp$popadjust){
      hivn_entrants <- fp$entrantpop[,i-1] * (1 - fp$entrantprev[,i])
      hivp_entrants <- fp$entrantpop[,i-1] * fp$entrantprev[,i]
    } else {
      hivn_entrants <- HnIn(fp, i, pregprevlag, birthslag)
      hivp_entrants <- HpIn(fp, i, pregprevlag, birthslag)
    }

    entrant_prev_out[i] <- sum(hivp_entrants) / sum(hivn_entrants+hivp_entrants)
    hivp_entrants_out[,i] <- sum(hivp_entrants)

    # distribute entrants by risk groups: assuming same for sexes and H± 

    pop[1,,hivn.idx,i] <- f_sa(hivn_entrants,,fp)
    pop[1,,hivp.idx,i] <- f_sa(hivp_entrants,,fp)

    # TODO: assummed entrantartcov the same for all risk groups: 
    noART <- hivp_entrants * (1-fp$entrantartcov[,i])
    isART <- hivp_entrants * fp$entrantartcov[,i]

    ## age the hiv group
    hiv.ag.prob <- pop[aglast.idx,,hivp.idx,i-1] / sumByAGs(pop[,,hivp.idx,i-1])
    hiv.ag.prob[is.nan(hiv.ag.prob)] <- 0
    
    hivpop[,,,i] <- hivpop[,,,i-1]
    nHup <- sweep(hivpop[,-hAG,,i-1], 2:3, hiv.ag.prob[-hAG,], "*")
    hivpop[,-hAG,,i] %<>% -(nHup)
    hivpop[,-1,,i]   %<>% +(nHup)
    hivpop[,1,,i]    %<>% +(f_sa(sweep(fp$paedsurv_cd4dist[,,i], 2, noART, "*"),,fp))

    ## age the on ART group
    if(i > fp$tARTstart){
      artpop[,,,,i]     <- artpop[,,,,i-1]
      nARTup            <- sweep(artpop[,,-hAG,,i-1], 3:4, hiv.ag.prob[-hAG,], "*")
      artpop[,,-hAG,,i] %<>% -(nARTup)
      artpop[,,-1,,i]   %<>% +(nARTup)
      artpop[,,1,,i]    %<>% +(sweep(f_sa(fp$paedsurv_artcd4dist[,,,i],TRUE,fp), 3, f_sa(isART,,fp), "*"))
    }
    ## survive the population
    deaths <- sweep(pop[,,,i], 1:2, 1 - f_sa(fp$Sx[,,i],TRUE,fp), "*")
    hiv.sx.prob <- 1 - sumByAGs(deaths[,,hivp.idx]) / sumByAGs(pop[,,hivp.idx,i])
    hiv.sx.prob[is.nan(hiv.sx.prob)] <- 0
    pop[,,,i] %<>% -(deaths)
    natdeaths[,,i] <- rowSums(deaths,,2)

    hivpop[,,,i]   <- sweep(hivpop[,,,i], 2:3, hiv.sx.prob, "*")
    if(i > fp$tARTstart)
      artpop[,,,,i] <- sweep(artpop[,,,,i], 3:4, hiv.sx.prob, "*")

    ## net migration
    netmigsurv  <- fp$netmigr[,,i] * (1 + fp$Sx[,,i]) / 2
    mr.prob     <- 1 + f_sa(netmigsurv,,fp) / rowSums(pop[,,,i],,2)
    hiv.mr.prob <- sumByAGs(mr.prob * pop[,,hivp.idx,i]) / sumByAGs(pop[,,hivp.idx,i])
    hiv.mr.prob[is.nan(hiv.mr.prob)] <- 0
    pop[,,,i]   <- sweep(pop[,,,i], 1:2, mr.prob, "*")
    
    hivpop[,,,i] <- sweep(hivpop[,,,i], 2:3, hiv.mr.prob, "*")
    if(i > fp$tARTstart)
      artpop[,,,,i] <- sweep(artpop[,,,,i], 3:4, hiv.mr.prob, "*")

    ## fertility
    births.by.age   <- rowSums(pop[p.fert.idx, sid(f.idx, fp),,i-1:0, drop=F],,2)/2 * fp$asfr[,i]
    births.by.h.age <- sumByAGs(births.by.age, TRUE)
    births          <- fp$srb[,i] * sum(births.by.h.age) # group is not defined at birth
                                                         # only sex is.
    if(i+AGE_START <= PROJ_YEARS)
      birthslag[,i+AGE_START-1] <- births


    ## ########################## ##
    ##  Disease model simulation  ##
    ## ########################## ##

    ## events at dt timestep
    message('Year: ', i, "\r", appendLF=FALSE); flush.console()

    for(ii in seq_len(hiv_steps_per_year)){

      ts <- (i-2)/DT + ii

      grad <- array(0, c(hDS, hAG, NG))

      if(fp$eppmod != "directincid"){
        ## incidence

        ## calculate r(t)
        if(fp$eppmod %in% c("rtrend", "rtrend_rw"))
          rvec[ts] <- calc_rtrend_rt(fp$proj.steps[ts], fp, rvec[ts-1], prevlast, pop, i, ii)
        else
          rvec[ts] <- fp$rvec[ts]

        ## number of infections by age / sex
        infections.ts <- calc_infections_eppspectrum(fp, mx, pop, hivpop, artpop,
                                                     i, ii, rvec[ts], isMixing)
        if (isMixing) FOI[,,i] <- attr(infections.ts, "FOI")
        incrate15to49.ts.out[ts] <- attr(infections.ts, "incrate15to49.ts")
        prev15to49.ts.out[ts] <- prevlast <- attr(infections.ts, "prevcurr")

        pop[,,hivn.idx,i] %<>% -(DT*infections.ts)
        pop[,,hivp.idx,i] %<>% +(DT*infections.ts)
        infections[,,i]   %<>% +(DT*infections.ts)
        grad              %<>% +(sweep(fp$cd4_initdist, 2:3, sumByAGs(infections.ts), '*'))
        incid15to49[i]    %<>% +(sum(DT*infections.ts[p.age15to49.idx,]))
      }

      # disease progression and mortality
      # cd4 stage progression (untreated)
      grad[-hDS,,] %<>% -(fp$cd4_prog * hivpop[-hDS,,,i])
      grad[-1,,]   %<>% +(fp$cd4_prog * hivpop[-hDS,,,i])

      cd4_mort_ts <- scale_cd4_mort_ts(fp, hivpop, artpop, i)
      grad %<>% -(cd4_mort_ts * hivpop[,,,i]) # HIV mortality, untreated

      ## Remove hivdeaths from pop
      hivdeaths_p.ts <- f_hdeaths_p.ts(cd4_mort_ts, fp, pop, hivpop, artpop, i)
      pop[,,hivp.idx,i] %<>% -(hivdeaths_p.ts)
      hivdeaths[,,i]    %<>% +(hivdeaths_p.ts)

      ## ART initiation
      if(i >= fp$tARTstart) {

        gradART <- array(0, c(hTS, hDS, hAG, NG))

        ## progression and mortality (HARD CODED 6 months duration)
        gradART[1:2,,,] %<>% -(2.0 * artpop[1:2,,,, i])
        gradART[2:3,,,] %<>% +(2.0 * artpop[1:2,,,, i])

        gradART %<>% -(fp$art_mort * fp$artmx_timerr[, i] * artpop[,,,,i])  # ART mortality

        ## ART dropout
        ## remove proportion from all adult ART groups back to untreated pop
        grad    %<>% +(fp$art_dropout[i] * colSums(artpop[,,,,i]))
        gradART %<>% -(fp$art_dropout[i] * artpop[,,,,i])

        ## calculate number eligible for ART
        artcd4_percelig <- f_artcd4_percelig(fp, i)
        art15plus.elig <- sweep(hivpop[,,,i], 1, artcd4_percelig, "*")

        ## calculate pregnant women
        if(fp$pw_artelig[i] && fp$artcd4elig_idx[i] > 1) 
          art15plus.elig <- updatePreg(art15plus.elig, births.by.h.age, sid,
                                       fp, i, pop, hivpop, artpop)

        ## calculate number to initiate ART based on number or percentage
        artpop_curr_g <- colSums(artpop[,,,,i],,3) + DT * colSums(gradART,,3)
        artnum.ii <- f_artInit(artpop_curr_g, art15plus.elig, fp, i, ii, sid)
        artpop_curr_g <- colSums(artpop[,,,,i],,3) + DT * colSums(gradART,,3) #?
        art15plus.inits <- pmax(artnum.ii - artpop_curr_g, 0)

        ## calculate ART initiation distribution
        artinit <- f_artDist(art15plus.elig, art15plus.inits, fp, i)
        artinit <- pmin(artinit, hivpop[,,, i] + DT * grad)
      
        grad <- grad - artinit / DT
        gradART[1,,,] <- gradART[1,,,] + artinit / DT
        artpop[,,,, i] <- artpop[,,,,i] + DT * gradART
      }
      hivpop[,,,i] <- hivpop[,,,i] + DT * grad
    } # end time steps within a year

    ## Direct incidence input
    if(fp$eppmod == "directincid"){
      if(fp$incidpopage == 0L) # incidence for 15-49 population
        p.incidpop.idx <- p.age15to49.idx
      else if(fp$incidpopage == 1L) # incidence for 15+ population
        p.incidpop.idx <- p.age15plus.idx
      incrate.i <- fp$incidinput[i]

      sexinc <- incrate.i*c(1, fp$incrr_sex[i])*sum(pop[p.incidpop.idx,,hivn.idx,i-1])/(sum(pop[p.incidpop.idx,m.idx,hivn.idx,i-1]) + fp$incrr_sex[i]*sum(pop[p.incidpop.idx, f.idx,hivn.idx,i-1]))
      agesex.inc <- sweep(fp$incrr_age[,,i], 2, sexinc/(colSums(pop[p.incidpop.idx,,hivn.idx,i-1] * fp$incrr_age[p.incidpop.idx,,i])/colSums(pop[p.incidpop.idx,,hivn.idx,i-1])), "*")
      infections[,,i] <- agesex.inc * pop[,,hivn.idx,i-1]
      
      pop[,,hivn.idx,i] <- pop[,,hivn.idx,i] - infections[,,i]
      pop[,,hivp.idx,i] <- pop[,,hivp.idx,i] + infections[,,i]
      
      hivpop[,,,i] <- hivpop[,,,i] + sweep(fp$cd4_initdist, 2:3, apply(infections[,,i], 2, sumByAG), "*")
      incid15to49[i] <- sum(infections[p.age15to49.idx,,i])
    }

    ## adjust population to match target population size
    if(exists("popadjust", where=fp) & fp$popadjust) {
      popadj.prob[,,i] <- f_sa(fp$targetpop[,,i],,fp) / rowSums(pop[,,,i],,2)
      popadj.prob[,,i][popadj.prob[,,i] < 0] <- 1 # DP's target <0 sometime 🤨
      hiv.popadj.prob <- sumByAGs(popadj.prob[,,i] * pop[,,hivp.idx,i]) / sumByAGs(pop[,,hivp.idx,i])
      hiv.popadj.prob[is.nan(hiv.popadj.prob)] <- 0
      pop[,,,i] <- sweep(pop[,,,i], 1:2, popadj.prob[,,i], "*")
      hivpop[,,,i] <- sweep(hivpop[,,,i], 2:3, hiv.popadj.prob, "*")
      if(i >= fp$tARTstart)
        artpop[,,,,i] <- sweep(artpop[,,,,i], 3:4, hiv.popadj.prob, "*")
    }

    ## prevalence among pregnant women
    hivn.byage <- sumByAGs(rowMeans(pop[p.fert.idx, sid(f.idx, fp), hivn.idx, i-1:0, drop=FALSE],,2), TRUE)
    hivp.byage <- rowMeans(hivpop[,h.fert.idx, sid(f.idx, fp),i-1:0, drop=FALSE],,3)
    artp.byage <- rowMeans(artpop[,,h.fert.idx, sid(f.idx, fp),i-1:0, drop=FALSE],,4)
    pregprev <- sum(births.by.h.age * (1 - hivn.byage / (hivn.byage + 
      colSums(sweep(hivp.byage, 1:2, fp$frr_cd4[,,i], '*')) + 
      colSums(sweep(artp.byage, 1:3, fp$frr_art[,,,i], '*'),,2)
      ))) / sum(births.by.age)
    if(i+AGE_START <= PROJ_YEARS)
      pregprevlag[i+AGE_START-1] <- pregprev

    ## prevalence and incidence 15 to 49
    prev15to49[i] <- sum(pop[p.age15to49.idx,,hivp.idx,i]) / sum(pop[p.age15to49.idx,,,i])
    incid15to49[i] <- sum(incid15to49[i]) / sum(pop[p.age15to49.idx,,hivn.idx,i-1])
  }

  if (isMixing) attr(pop, "FOI") <- FOI
  attr(pop, "prev15to49") <- prev15to49
  attr(pop, "incid15to49") <- incid15to49
  # attr(pop, "sexinc") <- sexinc15to49out
  attr(pop, "hivpop") <- hivpop
  attr(pop, "artpop") <- artpop

  attr(pop, "infections") <- infections
  attr(pop, "hivdeaths") <- hivdeaths
  attr(pop, "natdeaths") <- natdeaths

  attr(pop, "popadjust") <- popadj.prob

  attr(pop, "pregprevlag") <- pregprevlag

  if(fp$eppmod != "directincid"){
    attr(pop, "incrate15to49_ts") <- incrate15to49.ts.out
    attr(pop, "prev15to49_ts") <- prev15to49.ts.out
  }

  attr(pop, "entrantprev") <- entrant_prev_out
  attr(pop, "hivp_entrants") <- hivp_entrants_out
  class(pop) <- c("spec", "eppmix")
  if (isMixing) class(pop) <- c("eppmix", "spec")
  return(pop)
}