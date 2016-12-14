create_spectrum_fixpar <- function(projp, demp, hiv_steps_per_year = 10L, proj_start = projp$yr_start, proj_end = projp$yr_end,
                                   AGE_START = 15L, relinfectART = projp$relinfectART, time_epi_start = projp$t0,
                                   popadjust=FALSE, targetpop=demp$basepop){

  
  ## ########################## ##
  ##  Define model state space  ##
  ## ########################## ##

  ## Parameters defining the model projection period and state-space
  ss <- list(PROJ_YEARS = as.integer(proj_end - proj_start + 1L),
             AGE_START  = as.integer(AGE_START),
             hiv_steps_per_year = as.integer(hiv_steps_per_year),
             time_epi_start=time_epi_start)
             
  ## populuation projection state-space
  ss$NG <- 2
  ss$pDS <- 2               # Disease stratification for population projection (HIV-, and HIV+)

  ## macros
  ss$m.idx <- 1
  ss$f.idx <- 2

  ss$hivn.idx <- 1
  ss$hivp.idx <- 2

  ss$pAG <- 81 - AGE_START
  ss$ag.rate <- 1
  ss$p.fert.idx <- 16:50 - AGE_START
  ss$p.age15to49.idx <- 16:50 - AGE_START
  ss$p.age15plus.idx <- 16:ss$pAG - AGE_START

  
  ## HIV model state-space
  ss$h.ag.span <- as.integer(c(2,3, rep(5, 6), 31))   # Number of population age groups spanned by each HIV age group [sum(h.ag.span) = pAG]
  ss$hAG <- length(ss$h.ag.span)          # Number of age groups
  ss$hDS <- 7                             # Number of CD4 stages (Disease Stages)
  ss$hTS <- 3                             # number of treatment stages (including untreated)

  ss$ag.idx <- rep(1:ss$hAG, ss$h.ag.span)
  ss$aglast.idx <- which(!duplicated(ss$ag.idx, fromLast=TRUE))

  ss$h.fert.idx <- which((AGE_START-1 + cumsum(ss$h.ag.span)) %in% 15:49)
  ss$h.age15to49.idx <- which((AGE_START-1 + cumsum(ss$h.ag.span)) %in% 15:49)
  ss$h.age15plus.idx <- which((AGE_START-1 + cumsum(ss$h.ag.span)) >= 15)

  invisible(list2env(ss, environment())) # put ss variables in environment for convenience

  fp <- list(ss=ss)
  fp$proj.steps <- proj_start + 0.5 + 0:(ss$hiv_steps_per_year * (ss$PROJ_YEARS-1)) / ss$hiv_steps_per_year
  
  ## ######################## ##
  ##  Demographic parameters  ##
  ## ######################## ##

  ## linearly interpolate basepop if proj_start falls between indices
  bp_years <- as.integer(dimnames(demp$basepop)[[3]])
  bp_aidx <- max(which(proj_start >= bp_years))
  bp_dist <- 1-(proj_start - bp_years[bp_aidx]) / diff(bp_years[bp_aidx+0:1])
  basepop_allage <- rowSums(sweep(demp$basepop[,, bp_aidx+0:1], 3, c(bp_dist, 1-bp_dist), "*"),,2)

  fp$basepop <- basepop_allage[(AGE_START+1):81,]
  fp$Sx <- demp$Sx[(AGE_START+1):81,,as.character(proj_start:proj_end)]

  fp$asfr <- demp$asfr[,as.character(proj_start:proj_end)] # NOTE: assumes 15-49 is within projection age range
  ## Note: Spectrum averages ASFRs from the UPD file over 5-year age groups.
  ##       Prefer to use single-year of age ASFRs as provided. The below line will
  ##       convert to 5-year average ASFRs to exactly match Spectrum.
  ## fp$asfr <- apply(apply(fp$asfr, 2, tapply, rep(3:9*5, each=5), mean), 2, rep, each=5)
  
  fp$srb <- sapply(demp$srb[as.character(proj_start:proj_end)], function(x) c(x,100)/(x+100))
  
  ## Spectrum adjusts net-migration to occur half in current age group and half in next age group
  netmigr.adj <- demp$netmigr
  netmigr.adj[-1,,] <- (demp$netmigr[-1,,] + demp$netmigr[-81,,])/2
  netmigr.adj[1,,] <- demp$netmigr[1,,]/2
  netmigr.adj[81,,] <- netmigr.adj[81,,] + demp$netmigr[81,,]/2

  fp$netmigr <- netmigr.adj[(AGE_START+1):81,,as.character(proj_start:proj_end)]


  ## Calcuate the net-migration and survival up to AGE_START for each birth cohort.
  ## For cohorts born before projection start, this will be the partial
  ## survival since the projection start to AGE_START, and the corresponding lagged "births"
  ## represent the number in the basepop who will survive to the corresponding age.
  
  cumnetmigr <- array(0, dim=c(NG, PROJ_YEARS))
  cumsurv <- array(1, dim=c(NG, PROJ_YEARS))
  if(AGE_START > 0)
    for(i in 2:PROJ_YEARS)  # start at 2 because year 1 inputs are not used
      for(s in 1:2)
        for(j in max(1, AGE_START-(i-2)):AGE_START){
          ii <- i+j-AGE_START
          cumsurv[s,i] <- cumsurv[s,i] * demp$Sx[j,s,ii]
          if(j==1)
            cumnetmigr[s,i] <- netmigr.adj[j,s,ii] * (1+2*demp$Sx[j,s,ii])/3
          else
            cumnetmigr[s,i] <- cumnetmigr[s,i]*demp$Sx[j,s,ii] + netmigr.adj[j,s,ii] * (1+demp$Sx[j,s,ii])/2
        }
  
  ## initial values for births
  birthslag <- array(0, dim=c(NG, PROJ_YEARS))             # birthslag(i,s) = number of births of sex s, i-AGE_START years ago
  birthslag[,1:AGE_START] <- t(basepop_allage[AGE_START:1,])  # initial pop values (NOTE REVERSE ORDER). Rest will be completed by fertility during projection
  
  fp$birthslag <- birthslag
  fp$cumsurv <- cumsurv
  fp$cumnetmigr <- cumnetmigr


  ## set population adjustment
  fp$popadjust <- popadjust
  if(!length(setdiff(proj_start:proj_end, dimnames(targetpop)[[3]])))
    fp$targetpop <- targetpop[(AGE_START+1):81,,as.character(proj_start:proj_end)]
  if(popadjust & is.null(fp$targetpop))
    stop("targetpop does not span proj_start:proj_end")
  
  
  ## ###################### ##
  ##  HIV model parameters  ##
  ## ###################### ##

  fp$relinfectART <- projp$relinfectART

  fp$incrr_sex <- projp$incrr_sex[as.character(proj_start:proj_end)]
  
  projp.p.ag <- findInterval(AGE_START-1 + 1:pAG, seq(0, 85, 5))
  fp$incrr_age <- projp$incrr_age[projp.p.ag,,as.character(proj_start:proj_end)]
  
  projp.h.ag <- findInterval(AGE_START + cumsum(h.ag.span) - h.ag.span, c(15, 25, 35, 45))  # NOTE: Will not handle AGE_START < 15 presently
  fp$cd4_initdist <- projp$cd4_initdist[,projp.h.ag,]
  fp$cd4_prog <- projp$cd4_prog[,projp.h.ag,]
  fp$cd4_mort <- projp$cd4_mort[,projp.h.ag,]
  fp$art_mort <- projp$art_mort[,,projp.h.ag,]

  fert_rat.h.ag <- findInterval(AGE_START + cumsum(h.ag.span[h.fert.idx]) - h.ag.span[h.fert.idx], seq(15, 45, 5))
  

  fp$frr_cd4 <- array(1, c(hDS, length(h.fert.idx)))
  fp$frr_cd4[,] <- rep(projp$fert_rat[fert_rat.h.ag, "2005"], each=hDS)

  fp$frr_cd4[,] <- rep(projp$fert_rat[fert_rat.h.ag, "2005"], each=hDS)
  
  fp$frr_art <- array(1, c(hTS, hDS, length(h.fert.idx)))
  fp$frr_art[1:2,,] <- rep(projp$fert_rat[fert_rat.h.ag, "2005"], each=2*hDS)


  ## ART eligibility and numbers on treatment

  fp$art15plus_num <- projp$art15plus_num[,as.character(proj_start:proj_end)]
  fp$art15plus_isperc <- projp$art15plus_numperc[, as.character(proj_start:proj_end)] == 1

  ## convert percentage to proportion
  fp$art15plus_num[fp$art15plus_isperc] <- fp$art15plus_num[fp$art15plus_isperc] / 100

  ## eligibility starts in projection year idx+1
  fp$specpop_percelig <- rowSums(with(projp$artelig_specpop[-1,], mapply(function(elig, percent, year) rep(c(0, percent*as.numeric(elig)), c(year - proj_start+1, proj_end - year)), elig, percent, year)))
  fp$artcd4elig_idx <- findInterval(-projp$art15plus_eligthresh[as.character(proj_start:proj_end)], -c(999, 500, 350, 250, 200, 100, 50))

  fp$pw_artelig <- with(projp$artelig_specpop["PW",], rep(c(0, elig), c(year - proj_start+1, proj_end - year)))  # are pregnant women eligible (0/1)

  fp$tARTstart <- min(apply(fp$art15plus_num > 0, 1, which))

  
  ## Vertical transmission and survival to AGE_START for lagged births
  ## (fixed values assuming AGE_START = 15)
  ## NOTE: !! THIS COMPONENT OF THE MODEL NEEDS UPDATING
  
  fp$verttrans <- 0.4
  fp$paedsurv <- 0.09
  fp$paedsurv_cd4dist <- c(0.056,0.112,0.112,0.07,0.14,0.23,0.28)
  
  fp$netmig_hivprob <- 0.4*0.22
  fp$netmighivsurv <- 0.09/0.22


  ## ######################### ##
  ##  Prepare EPP r(t) models  ##
  ## ######################### ##

  fp$iota <- 0.0025
  proj.dur <- diff(range(fp$proj.steps))
  fp$numKnots <- 7
  rvec.knots <- seq(min(fp$proj.steps) - 3*proj.dur/(fp$numKnots-3), max(fp$proj.steps) + 3*proj.dur/(fp$numKnots-3), proj.dur/(fp$numKnots-3))
  fp$rvec.spldes <- splines::splineDesign(rvec.knots, fp$proj.steps)
  fp$tsEpidemicStart <- fp$proj.steps[which.min(abs(fp$proj.steps - fp$ss$time_epi_start))]

  fp$eppmod <- "rspline"  # default to r-spline model
  
  class(fp) <- "specfp"

  return(fp)
}


prepare_rtrend_model <- function(fp, iota=0.0025){
  fp$iota <- iota
  fp$tsEpidemicStart <- NULL
  fp$eppmod <- "rtrend"
  return(fp)
}


prepare_rspline_model <- function(fp, numKnots=7, tsEpidemicStart=fp$ss$time_epi_start){

  proj.dur <- diff(range(fp$proj.steps))
  fp$numKnots <- 7
  rvec.knots <- seq(min(fp$proj.steps) - 3*proj.dur/(fp$numKnots-3), max(fp$proj.steps) + 3*proj.dur/(fp$numKnots-3), proj.dur/(fp$numKnots-3))
  fp$rvec.spldes <- splines::splineDesign(rvec.knots, fp$proj.steps)

  fp$eppmod <- "rspline"
  fp$tsEpidemicStart <- fp$proj.steps[which.min(abs(fp$proj.steps - tsEpidemicStart))]
  fp$iota <- NULL

  return(fp)
}



simmod.specfp <- function(fp, VERSION="C"){

  if(VERSION != "R"){
    fp$eppmodInt <- as.integer(fp$eppmod == "rtrend") # 0: r-spline; 1: r-trend
    if(!exists("popadjust", where=fp))
      fp$popadjust <- FALSE
    mod <- .Call(spectrumC, fp)
    class(mod) <- "spec"
    return(mod)
  }

  ##################################################################################

  if(requireNamespace("fastmatch", quietly = TRUE))
    ctapply <- fastmatch::ctapply
  else
    ctapply <- tapply

  fp$ss$DT <- 1/fp$ss$hiv_steps_per_year
  
  ## Attach state space variables
  invisible(list2env(fp$ss, environment())) # put ss variables in environment for convenience

  birthslag <- fp$birthslag
  pregprevlag <- rep(0, PROJ_YEARS)

  ## initialize projection
  pop <- array(0, c(pAG, NG, pDS, PROJ_YEARS))
  pop[,,1,1] <- fp$basepop
  hivpop <- array(0, c(hTS+1L, hDS, hAG, NG, PROJ_YEARS))

  ## initialize output
  prev15to49 <- numeric(PROJ_YEARS)
  incid15to49 <- numeric(PROJ_YEARS)
  sexinc15to49out <- array(NA, c(NG, PROJ_YEARS))
  paedsurvout <- rep(NA, PROJ_YEARS)

  infections <- array(0, c(pAG, NG, PROJ_YEARS))
  hivdeaths <- array(0, c(pAG, NG, PROJ_YEARS))
  natdeaths <- array(0, c(pAG, NG, PROJ_YEARS))

  popadj.prob <- array(0, c(pAG, NG, PROJ_YEARS))

  incrate15to49.ts.out <- rep(NA, length(fp$rvec))
  rvec <- if(fp$eppmod == "rtrend") rep(NA, length(fp$proj.steps)) else fp$rvec

  prev15to49.ts.out <- rep(NA, length(fp$rvec))

  ## store last prevalence value (for r-trend model)
  prevlast <- prevcurr <- 0


  for(i in 2:PROJ_YEARS){

    ## ################################### ##
    ##  Single-year population projection  ##
    ## ################################### ##

    ## age the population
    pop[-c(1,pAG),,,i] <- pop[-(pAG-1:0),,,i-1]
    pop[pAG,,,i] <- pop[pAG,,,i-1] + pop[pAG-1,,,i-1] # open age group

    ## Add lagged births into youngest age group
    pop[1,,hivn.idx,i] <- birthslag[,i-1]*fp$cumsurv[,i-1]*(1-pregprevlag[i-1]*fp$verttrans) + fp$cumnetmigr[,i-1]*(1-pregprevlag[i-1]*fp$netmig_hivprob)
    pop[1,,hivp.idx,i] <- birthslag[,i-1]*fp$cumsurv[,i-1]*pregprevlag[i-1]*fp$verttrans*fp$paedsurv + fp$cumnetmigr[,i-1]*pregprevlag[i-1]*fp$netmig_hivprob*fp$netmighivsurv

    paedsurvout[i] <- sum(birthslag[,i-1]*fp$cumsurv[,i-1]*pregprevlag[i-1]*fp$verttrans*fp$paedsurv)

    hiv.ag.prob <- pop[aglast.idx,,hivp.idx,i-1] / apply(pop[,,hivp.idx,i-1], 2, ctapply, ag.idx, sum)
    hiv.ag.prob[is.nan(hiv.ag.prob)] <- 0

    hivpop[,,,,i] <- hivpop[,,,,i-1]
    hivpop[,,-hAG,,i] <- hivpop[,,-hAG,,i] - sweep(hivpop[,,-hAG,,i-1], 3:4, hiv.ag.prob[-hAG,], "*")
    hivpop[,,-1,,i] <- hivpop[,,-1,,i] + sweep(hivpop[,,-hAG,,i-1], 3:4, hiv.ag.prob[-hAG,], "*")
    hivpop[1,,1,,i] <- hivpop[1,,1,,i] + fp$paedsurv_cd4dist %o% (birthslag[,i-1]*fp$cumsurv[,i-1]*pregprevlag[i-1]*fp$verttrans*fp$paedsurv + fp$cumnetmigr[,i-1]*pregprevlag[i-1]*fp$netmig_hivprob*fp$netmighivsurv)

    ## survive the population
    deaths <- sweep(pop[,,,i], 1:2, (1-fp$Sx[,,i]), "*")
    hiv.sx.prob <- 1-apply(deaths[,,2], 2, ctapply, ag.idx, sum) / apply(pop[,,2,i], 2, ctapply, ag.idx, sum)
    hiv.sx.prob[is.nan(hiv.sx.prob)] <- 0
    pop[,,,i] <- pop[,,,i] - deaths
    natdeaths[,,i] <- rowSums(deaths,,2)

    hivpop[,,,,i] <- sweep(hivpop[,,,,i], 3:4, hiv.sx.prob, "*")

    ## net migration
    netmigsurv <- fp$netmigr[,,i]*(1+fp$Sx[,,i])/2
    mr.prob <- 1+netmigsurv / rowSums(pop[,,,i],,2)
    hiv.mr.prob <- apply(mr.prob * pop[,,2,i], 2, ctapply, ag.idx, sum) /  apply(pop[,,2,i], 2, ctapply, ag.idx, sum)
    hiv.mr.prob[is.nan(hiv.mr.prob)] <- 0
    pop[,,,i] <- sweep(pop[,,,i], 1:2, mr.prob, "*")

    hivpop[,,,,i] <- sweep(hivpop[,,,,i], 3:4, hiv.mr.prob, "*")

    ## fertility
    births.by.age <- rowSums(pop[p.fert.idx, f.idx,,i-1:0])/2 * fp$asfr[,i]
    births.by.h.age <- ctapply(births.by.age, ag.idx[p.fert.idx], sum)
    births <- fp$srb[,i] * sum(births.by.h.age)
    if(i+AGE_START <= PROJ_YEARS)
      birthslag[,i+AGE_START-1] <- births


    ## ########################## ##
    ##  Disease model simulation  ##
    ## ########################## ##

    ## events at dt timestep
    for(ii in seq_len(hiv_steps_per_year)){
      grad <- array(0, c(hTS+1L, hDS, hAG, NG))

      ## HIV population size at ts
      ts <- (i-2)/DT + ii

      hivn.ii <- sum(pop[p.age15to49.idx,,hivn.idx,i])
      hivn.ii <- hivn.ii - sum(pop[p.age15to49.idx[1],,hivn.idx,i])*(1-DT*(ii-1))
      hivn.ii <- hivn.ii + sum(pop[tail(p.age15to49.idx,1)+1,,hivn.idx,i])*(1-DT*(ii-1))

      hivp.ii <- sum(pop[p.age15to49.idx,,hivp.idx,i])
      hivp.ii <- hivp.ii - sum(pop[p.age15to49.idx[1],,hivp.idx,i])*(1-DT*(ii-1))
      hivp.ii <- hivp.ii + sum(pop[tail(p.age15to49.idx,1)+1,,hivp.idx,i])*(1-DT*(ii-1))

      ## there is an approximation here since this is the 15-49 pop (doesn't account for the slight offset in age group)
      propart.ii <- ifelse(hivp.ii > 0, sum(hivpop[-1,,h.age15to49.idx,,i])/sum(hivpop[,,h.age15to49.idx,,i]), 0)  

      
      ## incidence

      ## calculate r(t)
      prevlast <- prevcurr
      prev15to49.ts.out[ts] <- prevcurr <- hivp.ii / (hivn.ii+hivp.ii)
      if(fp$eppmod=="rtrend")
        rvec[ts] <- calc.rt(fp$proj.steps[ts], fp, rvec[ts-1L], prevlast, prevcurr)
      
      incrate15to49.ts <- rvec[ts] * hivp.ii * (1 - (1-fp$relinfectART)*propart.ii) / (hivn.ii+hivp.ii) + fp$iota * (fp$proj.steps[ts] == fp$tsEpidemicStart)
      sexinc15to49.ts <- incrate15to49.ts*c(1, fp$incrr_sex[i])*sum(pop[p.age15to49.idx,,hivn.idx,i])/(sum(pop[p.age15to49.idx,m.idx,hivn.idx,i]) + fp$incrr_sex[i]*sum(pop[p.age15to49.idx, f.idx,hivn.idx,i]))
      agesex.inc <- sweep(fp$incrr_age[,,i], 2, sexinc15to49.ts/(colSums(pop[p.age15to49.idx,,hivn.idx,i] * fp$incrr_age[p.age15to49.idx,,i])/colSums(pop[p.age15to49.idx,,hivn.idx,i])), "*")
      infections.ts <- agesex.inc * pop[,,hivn.idx,i]

      incrate15to49.ts.out[ts] <- incrate15to49.ts

      pop[,,hivn.idx,i] <- pop[,,hivn.idx,i] - DT*infections.ts
      pop[,,hivp.idx,i] <- pop[,,hivp.idx,i] + DT*infections.ts
      infections[,,i] <- infections[,,i] + DT*infections.ts

      grad[1,,,] <- grad[1,,,] + sweep(fp$cd4_initdist, 2:3, apply(infections.ts, 2, ctapply, ag.idx, sum), "*")
      incid15to49[i] <- incid15to49[i] + sum(DT*infections.ts[p.age15to49.idx,])
      
      ## disease progression and mortality
      grad[1,-hDS,,] <- grad[1,-hDS,,] - fp$cd4_prog * hivpop[1,-hDS,,,i]  # remove cd4 stage progression (untreated)
      grad[1,-1,,] <- grad[1,-1,,] + fp$cd4_prog * hivpop[1,-hDS,,,i]      # add cd4 stage progression (untreated)
      grad[2:3,,,] <- grad[2:3,,,] - 2.0 * hivpop[2:3,,,, i]               # remove ART duration progression (HARD CODED 6 months duration)
      grad[3:4,,,] <- grad[3:4,,,] + 2.0 * hivpop[2:3,,,, i]               # add ART duration progression (HARD CODED 6 months duration)

      grad[1,,,] <- grad[1,,,] - fp$cd4_mort * hivpop[1,,,,i]              # HIV mortality, untreated
      grad[-1,,,] <- grad[-1,,,] - fp$art_mort * hivpop[-1,,,,i]           # ART mortality

      ## Remove hivdeaths from pop
      hivdeaths.ts <- DT*(colSums(fp$cd4_mort * hivpop[1,,,,i]) + colSums(fp$art_mort * hivpop[-1,,,,i],,2))
      calc.agdist <- function(x) {d <- x/rep(ctapply(x, ag.idx, sum), h.ag.span); d[is.na(d)] <- 0; d}
      hivdeaths_p.ts <- apply(hivdeaths.ts, 2, rep, h.ag.span) * apply(pop[,,hivp.idx,i], 2, calc.agdist)  # HIV deaths by single-year age
      pop[,,2,i] <- pop[,,2,i] - hivdeaths_p.ts
      hivdeaths[,,i] <- hivdeaths[,,i] + hivdeaths_p.ts

      hivpop[,,,,i] <- hivpop[,,,,i] + DT*grad

      ## ART initiation
      if(sum(fp$art15plus_num[,i])>0){

        ## calculate number eligible for ART
        artcd4_percelig <- rep(c(fp$specpop_percelig[i], 1), times=c(fp$artcd4elig[i]-1, hDS - fp$artcd4elig[i]+1))
        art15plus.elig <- sweep(hivpop[1,,h.age15plus.idx,,i], 1, artcd4_percelig, "*")

        ## calculate pregnant women
        if(fp$pw_artelig[i]){
          births.dist <- sweep(fp$frr_cd4 * hivpop[1,,h.fert.idx,f.idx,i], 2,
                               births.by.h.age / (ctapply(pop[p.fert.idx, f.idx, hivn.idx, i], ag.idx[p.fert.idx], sum) + colSums(fp$frr_cd4 * hivpop[1,,h.fert.idx,f.idx,i]) + colSums(fp$frr_art * hivpop[-1,,h.fert.idx,f.idx,i],,2)), "*")
          art15plus.elig[1:(fp$artcd4elig_idx[i]-1),h.fert.idx-min(h.age15plus.idx)+1,f.idx] <- art15plus.elig[1:(fp$artcd4elig_idx[i]-1),h.fert.idx-min(h.age15plus.idx)+1,f.idx] + DT*births.dist[1:(fp$artcd4elig_idx[i]-1),] # multiply by DT to account for proportion of annual births occurring during this time step
        }

        ## calculate number to initiate ART based on number or percentage

        artnum.ii <- c(0,0) # number on ART this ts
        if(DT*ii < 0.5){
          for(g in 1:2){
            if(!any(fp$art15plus_isperc[g,i-2:1])){  # both number 
              artnum.ii[g] <- c(fp$art15plus_num[g,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
            } else if(all(fp$art15plus_isperc[g,i-2:1])){  # both percentage
              artcov.ii <- c(fp$art15plus_num[g,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + sum(hivpop[-1,,h.age15plus.idx,g,i]))
            } else if(!fp$art15plus_isperc[g,i-2] & fp$art15plus_isperc[g,i-1]){ # transition number to percentage
              curr_coverage <- sum(hivpop[-1,,h.age15plus.idx,g,i]) / (sum(art15plus.elig[,,g]) + sum(hivpop[-1,,h.age15plus.idx,g,i]))
              artcov.ii <- curr_coverage + (fp$art15plus_num[g,i-1] - curr_coverage) * DT/(0.5-DT*(ii-1))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + sum(hivpop[-1,,h.age15plus.idx,g,i]))
            }
          }
        } else {
          for(g in 1:2){
            if(!any(fp$art15plus_isperc[g,i-1:0])){  # both number 
              artnum.ii[g] <- c(fp$art15plus_num[g,i-1:0] %*% c(1-(DT*ii-0.5), DT*ii-0.5))
            } else if(all(fp$art15plus_isperc[g,i-1:0])) {  # both percentage
              artcov.ii <- c(fp$art15plus_num[g,i-1:0] %*% c(1-(DT*ii-0.5), DT*ii-0.5))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + sum(hivpop[-1,,h.age15plus.idx,g,i]))
            } else if(!fp$art15plus_isperc[g,i-1] & fp$art15plus_isperc[g,i]){  # transition number to percentage
              curr_coverage <- sum(hivpop[-1,,h.age15plus.idx,g,i]) / (sum(art15plus.elig[,,g]) + sum(hivpop[-1,,h.age15plus.idx,g,i]))
              artcov.ii <- curr_coverage + (fp$art15plus_num[g,i] - curr_coverage) * DT/(1.5-DT*(ii-1))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + sum(hivpop[-1,,h.age15plus.idx,g,i]))
            }
          }
        }

        art15plus.inits <- pmax(artnum.ii - colSums(hivpop[-1,,h.age15plus.idx,,i],,3), 0)

        ## calculate ART initiation distribution
        expect.mort.weight <- sweep(fp$cd4_mort[, h.age15plus.idx,], 3,
                                    colSums(art15plus.elig * fp$cd4_mort[, h.age15plus.idx,],,2), "/")
        artinit.weight <- sweep(expect.mort.weight, 3, 1/colSums(art15plus.elig,,2), "+")/2
        artinit <- pmin(sweep(artinit.weight * art15plus.elig, 3, art15plus.inits, "*"),
                        art15plus.elig)
        
        hivpop[1,, h.age15plus.idx,, i] <- hivpop[1,, h.age15plus.idx,, i] - artinit
        hivpop[2,, h.age15plus.idx,, i] <- hivpop[2,, h.age15plus.idx,, i] + artinit
      }
    }


    ## ## Code for calculating new infections once per year to match prevalence (like Spectrum)
    ## ## incidence
    ## prev.i <- sum(pop[p.age15to49.idx,,2,i]) / sum(pop[p.age15to49.idx,,,i]) # prevalence age 15 to 49
    ## incrate15to49.i <- (fp$prev15to49[i] - prev.i)/(1-prev.i)

    ## sexinc15to49 <- incrate15to49.i*c(1, fp$inc.sexratio[i])*sum(pop[p.age15to49.idx,,hivn.idx,i])/(sum(pop[p.age15to49.idx,m.idx,hivn.idx,i]) + fp$inc.sexratio[i]*sum(pop[p.age15to49.idx, f.idx,hivn.idx,i]))

    ## agesex.inc <- sweep(fp$inc.agerr[,,i], 2, sexinc15to49/(colSums(pop[p.age15to49.idx,,hivn.idx,i] * fp$inc.agerr[p.age15to49.idx,,i])/colSums(pop[p.age15to49.idx,,hivn.idx,i])), "*")
    ## infections <- agesex.inc * pop[,,hivn.idx,i]

    ## pop[,,hivn.idx,i] <- pop[,,hivn.idx,i] - infections
    ## pop[,,hivp.idx,i] <- pop[,,hivp.idx,i] + infections

    ## hivpop[1,,,,i] <- hivpop[1,,,,i] + sweep(fp$cd4.initdist, 2:3, apply(infections, 2, ctapply, ag.idx, sum), "*")

    ## adjust population to match target population size
    if(exists("popadjust", where=fp) & fp$popadjust){
      popadj.prob[,,i] <- fp$targetpop[,,i] / rowSums(pop[,,,i],,2)
      hiv.popadj.prob <- apply(popadj.prob[,,i] * pop[,,2,i], 2, ctapply, ag.idx, sum) /  apply(pop[,,2,i], 2, ctapply, ag.idx, sum)
      hiv.popadj.prob[is.nan(hiv.popadj.prob)] <- 0

      pop[,,,i] <- sweep(pop[,,,i], 1:2, popadj.prob[,,i], "*")
      hivpop[,,,,i] <- sweep(hivpop[,,,,i], 3:4, hiv.popadj.prob, "*")
    }

    ## prevalence among pregnant women
    hivn.byage <- ctapply(rowMeans(pop[p.fert.idx, f.idx, hivn.idx,i-1:0]), ag.idx[p.fert.idx], sum)
    hivp.byage <- rowMeans(hivpop[,,h.fert.idx, f.idx,i-1:0],,3)
    pregprev <- sum(births.by.h.age * (1 - hivn.byage / (hivn.byage + colSums(fp$frr_cd4 * hivp.byage[1,,]) + colSums(fp$frr_art * hivp.byage[-1,,],,2)))) / sum(births.by.age)
    if(i+AGE_START <= PROJ_YEARS)
      pregprevlag[i+AGE_START-1] <- pregprev

    ## prevalence and incidence 15 to 49
    prev15to49[i] <- sum(pop[p.age15to49.idx,,hivp.idx,i]) / sum(pop[p.age15to49.idx,,,i])
    incid15to49[i] <- sum(incid15to49[i]) / sum(pop[p.age15to49.idx,,hivn.idx,i-1])
  }

  attr(pop, "prev15to49") <- prev15to49
  attr(pop, "incid15to49") <- incid15to49
  attr(pop, "sexinc") <- sexinc15to49out
  attr(pop, "hivpop") <- hivpop[1,,,,]
  attr(pop, "artpop") <- hivpop[-1,,,,]

  attr(pop, "infections") <- infections
  attr(pop, "hivdeaths") <- hivdeaths
  attr(pop, "natdeaths") <- natdeaths

  attr(pop, "popadjust") <- popadj.prob
  
  attr(pop, "pregprevlag") <- pregprevlag
  attr(pop, "paedsurvout") <- paedsurvout
  attr(pop, "incrate15to49_ts") <- incrate15to49.ts.out
  attr(pop, "prev15to49_ts") <- prev15to49.ts.out
  class(pop) <- "spec"
  return(pop)
}

calc.rt <- function(t, fp, rveclast, prevlast, prevcurr){
  if(t > fp$tsEpidemicStart){
    par <- fp$rtrend
    gamma.t <- if(t < par$tStabilize) 0 else (prevcurr-prevlast)*(t - par$tStabilize) / (fp$ss$DT*prevlast)
    logr.diff <- par$beta[2]*(par$beta[1] - rveclast) + par$beta[3]*prevlast + par$beta[4]*gamma.t
      return(exp(log(rveclast) + logr.diff))
    } else
      return(fp$rtrend$r0)
}

update.specfp <- epp::update.eppfp


#########################
####  Model outputs  ####
#########################

## modprev15to49 <- function(mod, fp){colSums(mod[fp$ss$p.age15to49.idx,,fp$ss$hivp.idx,],,2) / colSums(mod[fp$ss$p.age15to49.idx,,,],,3)}
prev.spec <- function(mod, fp){ attr(mod, "prev15to49") }
incid.spec <- function(mod, fp){ attr(mod, "incid15to49") }
fnPregPrev.spec <- function(mod, fp) { attr(mod, "pregprev") }
