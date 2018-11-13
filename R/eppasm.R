simmod.specfp <- function(fp, VERSION="C"){

  if(!exists("popadjust", where=fp))
    fp$popadjust <- FALSE

  if(!exists("incidmod", where=fp))
    fp$incidmod <- "eppspectrum"

  if(VERSION != "R"){
    fp$eppmodInt <- match(fp$eppmod, c("rtrend", "directincid"), nomatch=0) # 0: r-spline;
    fp$incidmodInt <- match(fp$incidmod, c("eppspectrum", "transm"))-1L  # -1 for 0-based indexing
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

  ## Attach state space variables
  invisible(list2env(fp$ss, environment())) # put ss variables in environment for convenience

  birthslag <- fp$birthslag
  pregprevlag <- rep(0, PROJ_YEARS)

  ## initialize projection
  pop <- array(0, c(pAG, NG, pDS, PROJ_YEARS))
  pop[,,1,1] <- fp$basepop
  hivpop <- array(0, c(hDS, hAG, NG, PROJ_YEARS))
  artpop <- array(0, c(hTS, hDS, hAG, NG, PROJ_YEARS))

  ## initialize output
  prev15to49 <- numeric(PROJ_YEARS)
  incid15to49 <- numeric(PROJ_YEARS)
  sexinc15to49out <- array(NA, c(NG, PROJ_YEARS))
  paedsurvout <- rep(NA, PROJ_YEARS)

  infections <- array(0, c(pAG, NG, PROJ_YEARS))
  hivdeaths <- array(0, c(pAG, NG, PROJ_YEARS))
  natdeaths <- array(0, c(pAG, NG, PROJ_YEARS))

  popadj.prob <- array(0, c(pAG, NG, PROJ_YEARS))

  if(fp$eppmod != "directincid"){
    ## outputs by timestep
    incrate15to49.ts.out <- rep(NA, length(fp$rvec))
    rvec <- if(fp$eppmod == "rtrend") rep(NA, length(fp$proj.steps)) else fp$rvec

    prev15to49.ts.out <- rep(NA, length(fp$rvec))
  }

  entrant_prev_out <- numeric(PROJ_YEARS)
  hivp_entrants_out <- array(0, c(NG, PROJ_YEARS))

  ## store last prevalence value (for r-trend model)
  prevlast <- 0


  for(i in 2:fp$SIM_YEARS){

    ## ################################### ##
    ##  Single-year population projection  ##
    ## ################################### ##

    ## age the population
    pop[-c(1,pAG),,,i] <- pop[-(pAG-1:0),,,i-1]
    pop[pAG,,,i] <- pop[pAG,,,i-1] + pop[pAG-1,,,i-1] # open age group

    ## Add lagged births into youngest age group
    if(exists("entrantprev", where=fp))
      entrant_prev <- fp$entrantprev[,i]
    else
      entrant_prev <- rep(pregprevlag[i-1]*fp$verttrans_lag[i-1]*fp$paedsurv_lag[i-1], 2)

    if(exists("popadjust", where=fp) & fp$popadjust){
      hivn_entrants <- fp$entrantpop[,i-1]*(1-entrant_prev)
      hivp_entrants <- fp$entrantpop[,i-1]*entrant_prev
    } else {
      hivn_entrants <- birthslag[,i-1]*fp$cumsurv[,i-1]*(1-entrant_prev / fp$paedsurv_lag[i-1]) + fp$cumnetmigr[,i-1]*(1-pregprevlag[i-1]*fp$netmig_hivprob)
      hivp_entrants <- birthslag[,i-1]*fp$cumsurv[,i-1]*entrant_prev + fp$cumnetmigr[,i-1]*entrant_prev
    }

    entrant_prev_out[i] <- sum(hivp_entrants) / sum(hivn_entrants+hivp_entrants)
    hivp_entrants_out[,i] <- sum(hivp_entrants)

    pop[1,,hivn.idx,i] <- hivn_entrants
    pop[1,,hivp.idx,i] <- hivp_entrants

    hiv.ag.prob <- pop[aglast.idx,,hivp.idx,i-1] / apply(pop[,,hivp.idx,i-1], 2, ctapply, ag.idx, sum)
    hiv.ag.prob[is.nan(hiv.ag.prob)] <- 0

    hivpop[,,,i] <- hivpop[,,,i-1]
    hivpop[,-hAG,,i] <- hivpop[,-hAG,,i] - sweep(hivpop[,-hAG,,i-1], 2:3, hiv.ag.prob[-hAG,], "*")
    hivpop[,-1,,i] <- hivpop[,-1,,i] + sweep(hivpop[,-hAG,,i-1], 2:3, hiv.ag.prob[-hAG,], "*")
    hivpop[,1,,i] <- hivpop[,1,,i] + sweep(fp$paedsurv_cd4dist[,,i], 2, hivp_entrants * (1-fp$entrantartcov[,i]), "*")

    if(i > fp$tARTstart){
      artpop[,,,,i] <- artpop[,,,,i-1]
      artpop[,,-hAG,,i] <- artpop[,,-hAG,,i] - sweep(artpop[,,-hAG,,i-1], 3:4, hiv.ag.prob[-hAG,], "*")
      artpop[,,-1,,i] <- artpop[,,-1,,i] + sweep(artpop[,,-hAG,,i-1], 3:4, hiv.ag.prob[-hAG,], "*")
      artpop[,,1,,i] <- artpop[,,1,,i] + sweep(fp$paedsurv_artcd4dist[,,,i], 3, hivp_entrants * fp$entrantartcov[,i], "*")
    }

    ## survive the population
    deaths <- sweep(pop[,,,i], 1:2, (1-fp$Sx[,,i]), "*")
    hiv.sx.prob <- 1-apply(deaths[,,2], 2, ctapply, ag.idx, sum) / apply(pop[,,2,i], 2, ctapply, ag.idx, sum)
    hiv.sx.prob[is.nan(hiv.sx.prob)] <- 0
    pop[,,,i] <- pop[,,,i] - deaths
    natdeaths[,,i] <- rowSums(deaths,,2)

    hivpop[,,,i] <- sweep(hivpop[,,,i], 2:3, hiv.sx.prob, "*")
    if(i > fp$tARTstart)
      artpop[,,,,i] <- sweep(artpop[,,,,i], 3:4, hiv.sx.prob, "*")

    ## net migration
    netmigsurv <- fp$netmigr[,,i]*(1+fp$Sx[,,i])/2
    mr.prob <- 1+netmigsurv / rowSums(pop[,,,i],,2)
    hiv.mr.prob <- apply(mr.prob * pop[,,2,i], 2, ctapply, ag.idx, sum) /  apply(pop[,,2,i], 2, ctapply, ag.idx, sum)
    hiv.mr.prob[is.nan(hiv.mr.prob)] <- 0
    pop[,,,i] <- sweep(pop[,,,i], 1:2, mr.prob, "*")

    hivpop[,,,i] <- sweep(hivpop[,,,i], 2:3, hiv.mr.prob, "*")
    if(i > fp$tARTstart)
      artpop[,,,,i] <- sweep(artpop[,,,,i], 3:4, hiv.mr.prob, "*")

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
        if(exists("incidmod", where=fp) && fp$incidmod == "transm")
          infections.ts <- calc_infections_simpletransm(fp, pop, hivpop, artpop, i, ii, rvec[ts])
        else
          infections.ts <- calc_infections_eppspectrum(fp, pop, hivpop, artpop, i, ii, rvec[ts])

        incrate15to49.ts.out[ts] <- attr(infections.ts, "incrate15to49.ts")
        prev15to49.ts.out[ts] <- attr(infections.ts, "prevcurr")
        prevlast <- attr(infections.ts, "prevcurr")

        pop[,,hivn.idx,i] <- pop[,,hivn.idx,i] - DT*infections.ts
        pop[,,hivp.idx,i] <- pop[,,hivp.idx,i] + DT*infections.ts
        infections[,,i] <- infections[,,i] + DT*infections.ts

        grad <- grad + sweep(fp$cd4_initdist, 2:3, apply(infections.ts, 2, ctapply, ag.idx, sum), "*")
        incid15to49[i] <- incid15to49[i] + sum(DT*infections.ts[p.age15to49.idx,])
      }

      ## disease progression and mortality
      grad[-hDS,,] <- grad[-hDS,,] - fp$cd4_prog * hivpop[-hDS,,,i]  # remove cd4 stage progression (untreated)
      grad[-1,,] <- grad[-1,,] + fp$cd4_prog * hivpop[-hDS,,,i]      # add cd4 stage progression (untreated)

      grad <- grad - fp$cd4_mort * hivpop[,,,i]              # HIV mortality, untreated

      ## Remove hivdeaths from pop
      hivdeaths.ts <- DT*(colSums(fp$cd4_mort * hivpop[,,,i]) + colSums(fp$art_mort * fp$artmx_timerr[i] * artpop[,,,,i],,2))
      calc.agdist <- function(x) {d <- x/rep(ctapply(x, ag.idx, sum), h.ag.span); d[is.na(d)] <- 0; d}
      hivdeaths_p.ts <- apply(hivdeaths.ts, 2, rep, h.ag.span) * apply(pop[,,hivp.idx,i], 2, calc.agdist)  # HIV deaths by single-year age
      pop[,,2,i] <- pop[,,2,i] - hivdeaths_p.ts
      hivdeaths[,,i] <- hivdeaths[,,i] + hivdeaths_p.ts

      hivpop[,,,i] <- hivpop[,,,i] + DT*grad

      ## ART initiation
      if(i >= fp$tARTstart) {

        gradART <- array(0, c(hTS, hDS, hAG, NG))

        ## progression and mortality
        gradART[1:2,,,] <- gradART[1:2,,,] - 2.0 * artpop[1:2,,,, i]      # remove ART duration progression (HARD CODED 6 months duration)
        gradART[2:3,,,] <- gradART[2:3,,,] + 2.0 * artpop[1:2,,,, i]      # add ART duration progression (HARD CODED 6 months duration)

        gradART <- gradART - fp$art_mort * fp$artmx_timerr[i] * artpop[,,,,i]                  # ART mortality

        artpop[,,,, i] <- artpop[,,,, i] + DT * gradART


        ## ART dropout
        ## remove proportion from all adult ART groups back to untreated pop
        hivpop[,,,i] <- hivpop[,,,i] + DT*fp$art_dropout[i]*colSums(artpop[,,,,i])
        artpop[,,,,i] <- artpop[,,,,i] - DT*fp$art_dropout[i]*artpop[,,,,i]

        ## calculate number eligible for ART
        artcd4_percelig <- 1 - (1-rep(0:1, times=c(fp$artcd4elig_idx[i]-1, hDS - fp$artcd4elig_idx[i]+1))) *
          (1-rep(c(0, fp$who34percelig), c(2, hDS-2))) *
          (1-rep(fp$specpop_percelig[i], hDS))

        art15plus.elig <- sweep(hivpop[,h.age15plus.idx,,i], 1, artcd4_percelig, "*")

        ## calculate pregnant women
        if(fp$pw_artelig[i]){
          births.dist <- sweep(fp$frr_cd4[,,i] * hivpop[,h.fert.idx,f.idx,i], 2,
                               births.by.h.age / (ctapply(pop[p.fert.idx, f.idx, hivn.idx, i], ag.idx[p.fert.idx], sum) + colSums(fp$frr_cd4[,,i] * hivpop[,h.fert.idx,f.idx,i]) + colSums(fp$frr_art[,,,i] * artpop[ ,,h.fert.idx,f.idx,i],,2)), "*")
          if(fp$artcd4elig_idx[i] > 1)
            art15plus.elig[1:(fp$artcd4elig_idx[i]-1),h.fert.idx-min(h.age15plus.idx)+1,f.idx] <- art15plus.elig[1:(fp$artcd4elig_idx[i]-1),h.fert.idx-min(h.age15plus.idx)+1,f.idx] + births.dist[1:(fp$artcd4elig_idx[i]-1),] # multiply by DT to account for proportion of annual births occurring during this time step
        }

        ## calculate number to initiate ART based on number or percentage

        artnum.ii <- c(0,0) # number on ART this ts
        if(DT*ii < 0.5){
          for(g in 1:2){
            if(!any(fp$art15plus_isperc[g,i-2:1])){  # both number
              artnum.ii[g] <- c(fp$art15plus_num[g,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
            } else if(all(fp$art15plus_isperc[g,i-2:1])){  # both percentage
              artcov.ii <- c(fp$art15plus_num[g,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + sum(artpop[,,h.age15plus.idx,g,i]))
            } else if(!fp$art15plus_isperc[g,i-2] & fp$art15plus_isperc[g,i-1]){ # transition number to percentage
              curr_coverage <- sum(artpop[,,h.age15plus.idx,g,i]) / (sum(art15plus.elig[,,g]) + sum(artpop[,,h.age15plus.idx,g,i]))
              artcov.ii <- curr_coverage + (fp$art15plus_num[g,i-1] - curr_coverage) * DT/(0.5-DT*(ii-1))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + sum(artpop[,,h.age15plus.idx,g,i]))
            }
          }
        } else {
          for(g in 1:2){
            if(!any(fp$art15plus_isperc[g,i-1:0])){  # both number
              artnum.ii[g] <- c(fp$art15plus_num[g,i-1:0] %*% c(1-(DT*ii-0.5), DT*ii-0.5))
            } else if(all(fp$art15plus_isperc[g,i-1:0])) {  # both percentage
              artcov.ii <- c(fp$art15plus_num[g,i-1:0] %*% c(1-(DT*ii-0.5), DT*ii-0.5))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + sum(artpop[,,h.age15plus.idx,g,i]))
            } else if(!fp$art15plus_isperc[g,i-1] & fp$art15plus_isperc[g,i]){  # transition number to percentage
              curr_coverage <- sum(artpop[,,h.age15plus.idx,g,i]) / (sum(art15plus.elig[,,g]) + sum(artpop[,,h.age15plus.idx,g,i]))
              artcov.ii <- curr_coverage + (fp$art15plus_num[g,i] - curr_coverage) * DT/(1.5-DT*(ii-1))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + sum(artpop[,,h.age15plus.idx,g,i]))
            }
          }
        }

        art15plus.inits <- pmax(artnum.ii - colSums(artpop[,,h.age15plus.idx,,i],,3), 0)

        ## calculate ART initiation distribution
        if(!fp$med_cd4init_input[i]){

          if(fp$art_alloc_method == 4L){ ## by lowest CD4
            
            ## Calculate proportion to be initiated in each CD4 category
            artinit <- array(0, dim(art15plus.elig))
            remain_artalloc <- art15plus.inits
            for(m in hDS:1){
              elig_hm <- colSums(art15plus.elig[m,,])
              init_prop <- ifelse(elig_hm == 0, elig_hm, pmin(1.0, remain_artalloc / elig_hm, na.rm=TRUE))
              ## if(any(elig_hm < art15plus.inits))
              ##   browser()
              artinit[m , , ] <- sweep(art15plus.elig[m,,], 2, init_prop, "*")
              remain_artalloc <- remain_artalloc - init_prop * elig_hm
            }

          } else {
            expect.mort.weight <- sweep(fp$cd4_mort[, h.age15plus.idx,], 3,
                                        colSums(art15plus.elig * fp$cd4_mort[, h.age15plus.idx,],,2), "/")          
            artinit.weight <- sweep(fp$art_alloc_mxweight * expect.mort.weight, 3, (1 - fp$art_alloc_mxweight)/colSums(art15plus.elig,,2), "+")
            artinit <- pmin(sweep(artinit.weight * art15plus.elig, 3, art15plus.inits, "*"),
                            art15plus.elig)
            
            ## ## Allocation by average mortality across CD4, trying to match Spectrum
            ## artelig_by_cd4 <- apply(art15plus.elig, c(1, 3), sum)
            ## expectmort_by_cd4 <- apply(art15plus.elig * fp$cd4_mort[, h.age15plus.idx,], c(1, 3), sum)
            
            ## artinit_dist <- (sweep(artelig_by_cd4, 2, colSums(artelig_by_cd4), "/") +
            ##                  sweep(expectmort_by_cd4, 2, colSums(expectmort_by_cd4), "/")) / 2
            ## artinit <- sweep(art15plus.elig, c(1, 3), sweep(artinit_dist / artelig_by_cd4, 2, art15plus.inits, "*"), "*")
            ## artinit <- pmin(artinit, art15plus.elig, na.rm=TRUE)
          }

        } else {

          CD4_LOW_LIM <- c(500, 350, 250, 200, 100, 50, 0)
          CD4_UPP_LIM <- c(1000, 500, 350, 250, 200, 100, 50)

          medcd4_idx <- fp$med_cd4init_cat[i]

          medcat_propbelow <- (fp$median_cd4init[i] - CD4_LOW_LIM[medcd4_idx]) / (CD4_UPP_LIM[medcd4_idx] - CD4_LOW_LIM[medcd4_idx])

          elig_below <- colSums(art15plus.elig[medcd4_idx,,,drop=FALSE],,2) * medcat_propbelow
          if(medcd4_idx < hDS)
            elig_below <- elig_below + colSums(art15plus.elig[(medcd4_idx+1):hDS,,,drop=FALSE],,2)

          elig_above <- colSums(art15plus.elig[medcd4_idx,,,drop=FALSE],,2) * (1.0-medcat_propbelow)
          if(medcd4_idx > 1)
            elig_above <- elig_above + colSums(art15plus.elig[1:(medcd4_idx-1),,,drop=FALSE],,2)

          initprob_below <- pmin(art15plus.inits * 0.5 / elig_below, 1.0, na.rm=TRUE)
          initprob_above <- pmin(art15plus.inits * 0.5 / elig_above, 1.0, na.rm=TRUE)
          initprob_medcat <- initprob_below * medcat_propbelow + initprob_above * (1-medcat_propbelow)

          artinit <- array(0, dim=c(hDS, hAG, NG))

          if(medcd4_idx < hDS)
            artinit[(medcd4_idx+1):hDS,,] <- sweep(art15plus.elig[(medcd4_idx+1):hDS,,,drop=FALSE], 3, initprob_below, "*")
          artinit[medcd4_idx,,] <- sweep(art15plus.elig[medcd4_idx,,,drop=FALSE], 3, initprob_medcat, "*")
          if(medcd4_idx > 0)
            artinit[1:(medcd4_idx-1),,] <- sweep(art15plus.elig[1:(medcd4_idx-1),,,drop=FALSE], 3, initprob_above, "*")
        }

        hivpop[, h.age15plus.idx,, i] <- hivpop[, h.age15plus.idx,, i] - artinit
        artpop[1,, h.age15plus.idx,, i] <- artpop[1,, h.age15plus.idx,, i] + artinit
      }
    }


    ## ## Code for calculating new infections once per year to match prevalence (like Spectrum)
    ## ## incidence
    ## prev.i <- sum(pop[p.age15to49.idx,,2,i]) / sum(pop[p.age15to49.idx,,,i]) # prevalence age 15 to 49
    ## incrate15to49.i <- (fp$prev15to49[i] - prev.i)/(1-prev.i)

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
      
      hivpop[,,,i] <- hivpop[,,,i] + sweep(fp$cd4_initdist, 2:3, apply(infections[,,i], 2, ctapply, ag.idx, sum), "*")
      incid15to49[i] <- sum(infections[p.age15to49.idx,,i])
    }

    ## adjust population to match target population size
    if(exists("popadjust", where=fp) & fp$popadjust){
      popadj.prob[,,i] <- fp$targetpop[,,i] / rowSums(pop[,,,i],,2)
      hiv.popadj.prob <- apply(popadj.prob[,,i] * pop[,,2,i], 2, ctapply, ag.idx, sum) /  apply(pop[,,2,i], 2, ctapply, ag.idx, sum)
      hiv.popadj.prob[is.nan(hiv.popadj.prob)] <- 0

      pop[,,,i] <- sweep(pop[,,,i], 1:2, popadj.prob[,,i], "*")
      hivpop[,,,i] <- sweep(hivpop[,,,i], 2:3, hiv.popadj.prob, "*")
      if(i >= fp$tARTstart)
        artpop[,,,,i] <- sweep(artpop[,,,,i], 3:4, hiv.popadj.prob, "*")

    }

    ## prevalence among pregnant women
    hivn.byage <- ctapply(rowMeans(pop[p.fert.idx, f.idx, hivn.idx,i-1:0]), ag.idx[p.fert.idx], sum)
    hivp.byage <- rowMeans(hivpop[,h.fert.idx, f.idx,i-1:0],,2)
    artp.byage <- rowMeans(artpop[,,h.fert.idx, f.idx,i-1:0],,3)
    pregprev <- sum(births.by.h.age * (1 - hivn.byage / (hivn.byage + colSums(fp$frr_cd4[,,i] * hivp.byage) + colSums(fp$frr_art[,,,i] * artp.byage,,2)))) / sum(births.by.age)
    if(i+AGE_START <= PROJ_YEARS)
      pregprevlag[i+AGE_START-1] <- pregprev

    ## prevalence and incidence 15 to 49
    prev15to49[i] <- sum(pop[p.age15to49.idx,,hivp.idx,i]) / sum(pop[p.age15to49.idx,,,i])
    incid15to49[i] <- sum(incid15to49[i]) / sum(pop[p.age15to49.idx,,hivn.idx,i-1])
  }


  attr(pop, "prev15to49") <- prev15to49
  attr(pop, "incid15to49") <- incid15to49
  attr(pop, "sexinc") <- sexinc15to49out
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
  class(pop) <- "spec"
  return(pop)
}
