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
  
  ## paediatric
  popu5 <- array(0, c(pAGu5, NG, pDS, PROJ_YEARS))
  popu5[,,1,1] <- fp$paedbasepop[1:5,]
  hivpopu5 <- array(0, c(hMT, hDSu5, pAGu5, NG, PROJ_YEARS))
  dimnames(hivpopu5) <- list(transmission = c('BF0', 'BF6', 'BF12', 'perinatal'), cat = 1:7, age = 0:4, sex = c('Male', 'Female'), year = 1:PROJ_YEARS)
  artpopu5 <- array(0, c(hTS, hDSu5, pAGu5, NG, PROJ_YEARS))
  dimnames(artpopu5) <- list(transmission = c('ART0MOS', 'ART6MOS', 'ART1YR'), cat = 1:7, age = 0:4, sex = c('Male', 'Female'), year = 1:PROJ_YEARS)
  
  
  popu15 <- array(0, c(pAGu15, NG, pDS, PROJ_YEARS))
  popu15[,,1,1] <- fp$paedbasepop[6:15,]
  hivpopu15 <- array(0, c(hMT, hDSu15, pAGu15, NG, PROJ_YEARS))
  dimnames(hivpopu15) <- list(transmission = c('BF0', 'BF6', 'BF12', 'perinatal'), cat = 1:6, age = 5:14, sex = c('Male', 'Female'), year = 1:PROJ_YEARS)
  artpopu15 <- array(0, c(hTS, hDSu15, pAGu15, NG, PROJ_YEARS))
  dimnames(artpopu15) <- list(transmission = c('ART0MOS', 'ART6MOS', 'ART1YR'), cat = 1:6, age = 5:14, sex = c('Male', 'Female'), year = 1:PROJ_YEARS)
  
  
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
    entrant_prev <- fp$entrantprev[,i]

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

      if(fp$scale_cd4_mort == 1){
        cd4mx_scale <- hivpop[,,,i] / (hivpop[,,,i] + colSums(artpop[,,,,i]))
        cd4mx_scale[!is.finite(cd4mx_scale)] <- 1.0
        cd4_mort_ts <- fp$cd4_mort * cd4mx_scale
      } else
        cd4_mort_ts <- fp$cd4_mort
    
      ## Remove hivdeaths from pop
      hivdeaths.ts <- DT*(colSums(cd4_mort_ts * hivpop[,,,i]) + colSums(fp$art_mort * fp$artmx_timerr[ , i] * artpop[,,,,i],,2))
      if(exists('deaths_dt', where = fp)){
        delta_ts <- fp$deaths_dt[,ts] / colSums(hivdeaths.ts)
        delta_ts[!is.finite(delta_ts)] <- 1
      }else{
        delta_ts <- 1
      }
      
      grad <- grad - sweep(cd4_mort_ts * hivpop[,,,i], 3, delta_ts, "*")            # HIV mortality, untreated

      calc.agdist <- function(x) {d <- x/rep(ctapply(x, ag.idx, sum), h.ag.span); d[is.na(d)] <- 0; d}
      hivdeaths_p.ts <- apply(hivdeaths.ts, 2, rep, h.ag.span) * apply(pop[,,hivp.idx,i], 2, calc.agdist)  # HIV deaths by single-year age
      hivdeaths_p.ts <- sweep(hivdeaths_p.ts, 2, delta_ts, "*")
      pop[,,2,i] <- pop[,,2,i] - hivdeaths_p.ts
      # pop[,,2,i][pop[,,2,i] < 0] <- 0.01
      hivdeaths[,,i] <- hivdeaths[,,i] + hivdeaths_p.ts

      ## ART initiation
      if(i >= fp$tARTstart) {

        gradART <- array(0, c(hTS, hDS, hAG, NG))

        ## progression and mortality
        gradART[1:2,,,] <- gradART[1:2,,,] - 2.0 * artpop[1:2,,,, i]      # remove ART duration progression (HARD CODED 6 months duration)
        gradART[2:3,,,] <- gradART[2:3,,,] + 2.0 * artpop[1:2,,,, i]      # add ART duration progression (HARD CODED 6 months duration)

        gradART <- gradART - sweep(fp$art_mort * fp$artmx_timerr[ , i] * artpop[,,,,i], 4, delta_ts, "*")   # ART mortality

        ## ART dropout
        ## remove proportion from all adult ART groups back to untreated pop
        grad <- grad + fp$art_dropout[i]*colSums(artpop[,,,,i])
        gradART <- gradART - fp$art_dropout[i]*artpop[,,,,i]

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
            art15plus.elig[1:(fp$artcd4elig_idx[i]-1),h.fert.idx-min(h.age15plus.idx)+1,f.idx] <- art15plus.elig[1:(fp$artcd4elig_idx[i]-1),h.fert.idx-min(h.age15plus.idx)+1,f.idx] + births.dist[1:(fp$artcd4elig_idx[i]-1),]
        }

        ## calculate number to initiate ART based on number or percentage

        artpop_curr_g <- colSums(artpop[,,h.age15plus.idx,,i],,3) + DT*colSums(gradART[,,h.age15plus.idx,],,3)
        artnum.ii <- c(0,0) # number on ART this ts
        if(DT*ii < 0.5){
          for(g in 1:2){
            if(!any(fp$art15plus_isperc[g,i-2:1])){  # both number
              artnum.ii[g] <- c(fp$art15plus_num[g,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
            } else if(all(fp$art15plus_isperc[g,i-2:1])){  # both percentage
              artcov.ii <- c(fp$art15plus_num[g,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + artpop_curr_g[g])
            } else if(!fp$art15plus_isperc[g,i-2] & fp$art15plus_isperc[g,i-1]){ # transition number to percentage
              curr_coverage <- artpop_curr_g[g] / (sum(art15plus.elig[,,g]) + artpop_curr_g[g])
              artcov.ii <- curr_coverage + (fp$art15plus_num[g,i-1] - curr_coverage) * DT/(0.5-DT*(ii-1))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + artpop_curr_g[g])
            }
          }
        } else {
          for(g in 1:2){
            if(!any(fp$art15plus_isperc[g,i-1:0])){  # both number
              artnum.ii[g] <- c(fp$art15plus_num[g,i-1:0] %*% c(1-(DT*ii-0.5), DT*ii-0.5))
            } else if(all(fp$art15plus_isperc[g,i-1:0])) {  # both percentage
              artcov.ii <- c(fp$art15plus_num[g,i-1:0] %*% c(1-(DT*ii-0.5), DT*ii-0.5))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + artpop_curr_g[g])
            } else if(!fp$art15plus_isperc[g,i-1] & fp$art15plus_isperc[g,i]){  # transition number to percentage
              curr_coverage <- artpop_curr_g[g] / (sum(art15plus.elig[,,g]) + artpop_curr_g[g])
              artcov.ii <- curr_coverage + (fp$art15plus_num[g,i] - curr_coverage) * DT/(1.5-DT*(ii-1))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + artpop_curr_g[g])
            }
          }
        }

        artpop_curr_g <- colSums(artpop[,,h.age15plus.idx,,i],,3) + DT*colSums(gradART[,,h.age15plus.idx,],,3)
        art15plus.inits <- pmax(artnum.ii - artpop_curr_g, 0)

        ## calculate ART initiation distribution
        if(!fp$med_cd4init_input[i]){

          if(fp$art_alloc_method == 4L){ ## by lowest CD4
            
            ## Calculate proportion to be initiated in each CD4 category
            artinit <- array(0, dim(art15plus.elig))
            remain_artalloc <- art15plus.inits
            for(m in hDS:1){
              elig_hm <- colSums(art15plus.elig[m,,])
              init_prop <- ifelse(elig_hm == 0, elig_hm, pmin(1.0, remain_artalloc / elig_hm, na.rm=TRUE))
              artinit[m , , ] <- sweep(art15plus.elig[m,,], 2, init_prop, "*")
              remain_artalloc <- remain_artalloc - init_prop * elig_hm
            }

          } else {
            
            expect.mort.weight <- sweep(fp$cd4_mort[, h.age15plus.idx,], 3,
                                        colSums(art15plus.elig * fp$cd4_mort[, h.age15plus.idx,],,2), "/")          
            artinit.weight <- sweep(fp$art_alloc_mxweight * expect.mort.weight, 3, (1 - fp$art_alloc_mxweight)/colSums(art15plus.elig,,2), "+")
            artinit <- pmin(sweep(artinit.weight * art15plus.elig, 3, art15plus.inits, "*"),
                            art15plus.elig)
            
            ## Allocation by average mortality across CD4, trying to match Spectrum
            ## artelig_by_cd4 <- apply(art15plus.elig, c(1, 3), sum)
            ## expectmort_by_cd4 <- apply(art15plus.elig * fp$cd4_mort[, h.age15plus.idx,], c(1, 3), sum)
            
            ## artinit_dist <- fp$art_alloc_mxweight * sweep(artelig_by_cd4, 2, colSums(artelig_by_cd4), "/") +
            ##   (1 - fp$art_alloc_mxweight) * sweep(expectmort_by_cd4, 2, colSums(expectmort_by_cd4), "/")
            
            ## artinit_prob <- sweep(artinit_dist, 2, art15plus.inits, "*") / artelig_by_cd4
            ## artinit <- sweep(art15plus.elig, c(1, 3), artinit_prob, "*")
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

        artinit <- pmin(artinit, hivpop[ , , , i] + DT * grad)
      
        grad[ , h.age15plus.idx, ] <- grad[ , h.age15plus.idx, ] - artinit / DT
        gradART[1, , h.age15plus.idx, ] <- gradART[1, , h.age15plus.idx, ] + artinit / DT
        artpop[,,,, i] <- artpop[,,,, i] + DT * gradART
      }

      hivpop[,,,i] <- hivpop[,,,i] + DT * grad
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

    if(exists('paedbasepop', where = fp)){
      ## stick births in popu5 object
      ## pull from taget pop if popadjust, else use births calculated using input fert
      popu5[1,,1,i] <- ifelse(popadjust, fp$paedtargetpop[1,,i], births)
      
      ## calculate vertical transmission
      ## MTC transmission is weighted avg of those receing PMTCT and not treated
      ## Calculate percent of women by treatment option
      ## Default is no treatment
      needPMTCT <- pregprev
      treat.opt <- list('no_proph' = 1)
      ##TODO: make sure prenat and postnat options are correct
      prenat.opt <- c('tripleARTdurPreg', 'tripleARTbefPreg', 'singleDoseNevir', 'prenat_optionB', 'prenat_optionA', 'dualARV')
      postnat.opt <- c('tripleARTdurPreg', 'tripleARTbefPreg', 'singleDoseNevir', 'postnat_optionB', 'postnat_optionA', 'dualARV')
      for(n in prenat.opt){
        treat.opt[[n]] <- 0
      }
      sum3 <- rowSums(fp$pmtct_num[fp$pmtct_num$year == (proj_start + i - 1),prenat.opt])
      ## if number, set sum3 equal to total need, unless sum3 > need
      if(!all(as.logical(fp$pmtct_isperc[fp$pmtct_isperc$year == (proj_start + i - 1),c(names(fp$pmtct_isperc)[!names(fp$pmtct_isperc) == 'year'])]))){
        perc <- FALSE
        if(sum3 < needPMTCT){sum3 <- needPMTCT}
      }else{perc <- TRUE}
      on.treat <- 0
      for(n in prenat.opt){
        treat.opt[[n]] <- ifelse(sum3 == 0, 0, fp$pmtct_num[fp$pmtct_num$year == (proj_start + i - 1), n] / sum3)
        on.treat <- on.treat + treat.opt[[n]]
      }
      if(on.treat > 1){on.treat <- 1}
      treat.opt[['no_proph']] <- 1 - on.treat
      
      ## proportion of hiv+ preg women by cd4
      proplt200 <- (sum(births.by.h.age * (1 - hivn.byage / (hivn.byage + colSums(fp$frr_cd4[5:7,,i] * hivp.byage[5:7,]) + colSums(fp$frr_art[,5:7,,i] * artp.byage[,5:7,],,2)))) / sum(births.by.age)) / pregprev
      prop200to350 <- (sum(births.by.h.age * (1 - hivn.byage / (hivn.byage + colSums(fp$frr_cd4[3:4,,i] * hivp.byage[3:4,]) + colSums(fp$frr_art[,3:4,,i] * artp.byage[,3:4,],,2)))) / sum(births.by.age)) / pregprev
      propgt350 <- (sum(births.by.h.age * (1 - hivn.byage / (hivn.byage + colSums(fp$frr_cd4[1:2,,i] * hivp.byage[1:2,]) + colSums(fp$frr_art[,1:2,,i] * artp.byage[,1:2,],,2)))) / sum(births.by.age)) / pregprev
      
      tripleARTbefPreg <- fp$pmtct_num[fp$pmtct_num$year == (proj_start + i - 1), 'tripleARTbefPreg']
      tripleARTdurPreg <- fp$pmtct_num[fp$pmtct_num$year == (proj_start + i - 1), 'tripleARTdurPreg']
      if(perc){
        tripleARTbefPreg <- tripleARTbefPreg/100 * needPMTCT
        tripleARTdurPreg <- tripleARTdurPreg/100 * needPMTCT
      }
      ## number of hiv+ preg women in need of treatment
      ## ??
      hivp.preg.lt350 <- (proplt200 + prop200to350) * (needPMTCT - tripleARTbefPreg - tripleARTdurPreg) + tripleARTbefPreg + tripleARTdurPreg
      
      ## adjust transmission rates for opt a and opt b in cd4 lt 350
      if(treat.opt[['prenat_optionA']] + treat.opt[['prenat_optionB']] > (propgt350)){
        if(propgt350 <= 0){
          excess.ratio <- 0
        }else{
            excess.ratio <- (treat.opt[['prenat_optionA']] + treat.opt[['prenat_optionB']]) /propgt350 - 1
            optA.trans.rate <- (fp$MTCtrans[regimen == 'Option_A', perinatal_trans_pct]/100) * (1 + excess.ratio)
            optB.trans.rate <- (fp$MTCtrans[regimen == 'Option_B', perinatal_trans_pct]/100) * (1 + excess.ratio)
          }
      }else{
        optA.trans.rate <- fp$MTCtrans[regimen == 'Option_A', perinatal_trans_pct]/100
        optB.trans.rate <- fp$MTCtrans[regimen == 'Option_B', perinatal_trans_pct]/100
      }
      ## on treatment
      PTR <- treat.opt[['prenat_optionA']] * optA.trans.rate +
                treat.opt[['prenat_optionB']]* optB.trans.rate +
                treat.opt[['dualARV']] * (fp$MTCtrans[fp$MTCtrans$regimen == 'WHO_06_dual_ARV', 'perinatal_trans_pct'] / 100) +
                treat.opt[['singleDoseNevir']] * (fp$MTCtrans[fp$MTCtrans$regimen == 'single_dose_nevriapine', 'perinatal_trans_pct'] / 100) +
                treat.opt[['tripleARTbefPreg']] * (fp$MTCtrans[fp$MTCtrans$regimen == 'ART' & fp$MTCtrans$definition == 'start_pre_preg', 'perinatal_trans_pct'] / 100) +
                treat.opt[['tripleARTdurPreg']] * (fp$MTCtrans[fp$MTCtrans$regimen == 'ART' & fp$MTCtrans$definition == 'start_dur_preg', 'perinatal_trans_pct'] / 100)
      ## not on treatment
      ## TODO add incident infections
      PTR <- PTR + treat.opt[['no_proph']] * ((proplt200 * (fp$MTCtrans[fp$MTCtrans$regimen == 'no_prophylaxis' & fp$MTCtrans$definition == 'exisiting_LT200CD4', 'perinatal_trans_pct'] / 100)) +    
              (prop200to350 * (fp$MTCtrans[fp$MTCtrans$regimen == 'no_prophylaxis' & fp$MTCtrans$definition == 'exisiting_200to350CD4', 'perinatal_trans_pct'] / 100)) + 
              (propgt350 * (fp$MTCtrans[fp$MTCtrans$regimen == 'no_prophylaxis' & fp$MTCtrans$definition == 'exisiting_GT340CD4', 'perinatal_trans_pct'] / 100)))
      ## HIV births is birth prevalence (perinatal transmission)
      hiv.births <- max(0, unlist(pregprev * PTR))
      
      treat.opt <- list('no_proph' = 1)
      for(n in postnat.opt){
        treat.opt[[n]] <- 0
      }
      sum3 <- rowSums(fp$pmtct_num[fp$pmtct_num$year == (proj_start + i - 1),prenat.opt])
      ## if number, set sum3 equal to total need, unless sum3 > need
      if(!all(as.logical(fp$pmtct_isperc[fp$pmtct_isperc$year == (proj_start + i - 1),c(names(fp$pmtct_isperc)[!names(fp$pmtct_isperc) == 'year'])]))){
        perc <- FALSE
        if(sum3 < needPMTCT){sum3 <- needPMTCT}
      }else{perc <- TRUE}
      on.treat <- 0
      ##TODO: revisit how treatment props are being calculated, add code for opt A and B
      for(n in postnat.opt){
        treat.opt[[n]] <- ifelse(sum3 == 0, 0, fp$pmtct_num[fp$pmtct_num$year == (proj_start + i - 1), n] / sum3)
        on.treat <- on.treat + treat.opt[[n]]
      }
      if(on.treat > 1){on.treat <- 1}
      treat.opt[['no_proph']] <- 1 - on.treat
      
      ## breastfeeding transmission
      BFTR <- calcBFtransmissions(1,3,i)
      
      ##TODO should bftr be divided by 100?
      newInfFromBFLT6 <- as.numeric((pregprev - hiv.births) * BFTR / 100)
      cumNewInfFromBF <- newInfFromBFLT6
      
      
      BFTR <- calcBFtransmissions(4, 6, i)
      newInfFromBF6TO12 <- as.numeric((pregprev - hiv.births - newInfFromBFLT6) * BFTR / 100)
      cumNewInfFromBF <- cumNewInfFromBF + newInfFromBF6TO12
      
      ## perinatal infections
      hivpopu5[4,,1,m.idx,i] <- hiv.births * (1 - fp$srb[m.idx,i]) * fp$paed_distnewinf
      hivpopu5[4,,1,f.idx,i] <- hiv.births * fp$srb[f.idx,i] * fp$paed_distnewinf
      
      ## bf infections
      hivpopu5[1,,1,m.idx,i] <- newInfFromBFLT6 * (1 - fp$srb[m.idx,i]) * fp$paed_distnewinf
      hivpopu5[1,,1,f.idx,i] <- newInfFromBFLT6 * fp$srb[f.idx,i] * fp$paed_distnewinf
      
      ## remove hiv+ from pop
      popu5[1, m.idx, hivn.idx, i] <- popu5[1, m.idx, hivn.idx, i] - sum(hivpopu5[,,1,m.idx,i])
      popu5[1, m.idx, hivp.idx, i] <- popu5[1,m.idx,hivp.idx,i] + sum(hivpopu5[,,1,m.idx,i])
      popu5[1, f.idx, hivn.idx, i] <- popu5[1, f.idx, hivn.idx, i] - sum(hivpopu5[,,1,f.idx,i])
      popu5[1, f.idx, hivp.idx, i] <- popu5[1,f.idx,hivp.idx,i] + sum(hivpopu5[,,1,f.idx,i])
      
      ## age 1, bf12 transmission
      percentExposed = (pregprev - hiv.births - cumNewInfFromBF) / sum(births)
      BFTR = calcBFtransmissions(7, 12, i)
      m.bf12 <- popu5[2,m.idx,hivn.idx,i] * percentExposed * BFTR * fp$paed_distnewinf
      hivpopu5[3,,2,m.idx,i] <- ifelse(length(m.bf12) > 1, m.bf12, rep(0, 7))
      popu5[2,m.idx,hivn.idx,i] <- popu5[2,m.idx,hivn.idx,i] - sum(hivpopu5[3,,2,m.idx,i])
      f.bf12 <- popu5[2,f.idx,hivn.idx,i] * fp$paed_distnewinf * percentExposed * BFTR
      hivpopu5[3,,2,f.idx,i] <- ifelse(length(f.bf12) > 1, f.bf12, rep(0, 7))
      popu5[2,f.idx,hivn.idx,i] <- popu5[2,f.idx,hivn.idx,i] - sum(hivpopu5[3,,2,f.idx,i])  
      cumNewInfFromBF <- cumNewInfFromBF + sum(hivpopu5[3,,2,,i])  
      
      ## age 2, bf12 transmission
      percentExposed <- percentExposed * (1 - BFTR)
      BFTR = calcBFtransmissions(13, 18, i)
      m.bf12 <- popu5[3,m.idx,hivn.idx,i] * fp$paed_distnewinf * percentExposed * BFTR
      hivpopu5[3,,3,m.idx,i] <- ifelse(length(m.bf12) > 1, m.bf12, rep(0, 7))
      popu5[3,m.idx,hivn.idx,i] <- popu5[3,m.idx,hivn.idx,i] - sum(hivpopu5[3,,3,m.idx,i])
      f.bf12 <- popu5[3,f.idx,hivn.idx,i] * fp$paed_distnewinf * percentExposed * BFTR
      hivpopu5[3,,3,f.idx,i] <- ifelse(length(f.bf12) > 1, f.bf12, rep(0, 7))
      popu5[3,f.idx,hivn.idx,i] <- popu5[3,f.idx,hivn.idx,i] - sum(hivpopu5[3,,3,f.idx,i])  
      cumNewInfFromBF <- cumNewInfFromBF + sum(hivpopu5[3,,3,,i])        
      
      #Calculating need for child ART
      unmetNeed <- 0
      childEligibilityAge <- unique(fp$paed_arteligibility[,c('year', 'age_below_all_treat_mos')])$age_below_all_treat_mos
      for (age in 0:4){
        if((age + 1) * 12 > childEligibilityAge[i]){
          CD4elig <- as.character(min(fp$paed_arteligibility[fp$paed_arteligibility$year == (proj_start + i - 1) & fp$paed_arteligibility$age_start <= age, 'cd4_pct_thresh']))
          unmetNeed <- unmetNeed + sum(hivpopu5[,u5.elig.groups[[CD4elig]]:7, age + 1, , i])
        }else{
          unmetNeed <- unmetNeed + sum(hivpopu5[,, age + 1, , i])
        }
      }
      for (age in 5:14){
        if((age + 1) * 12 > childEligibilityAge[i]){
          CD4elig <- as.character(min(fp$paed_arteligibility[fp$paed_arteligibility$year == (proj_start + i - 1) & fp$paed_arteligibility$age_start <= age, 'cd4_count_thresh']))
          unmetNeed <- unmetNeed + sum(hivpopu15[,u15.elig.groups[[CD4elig]]:6, age - 4, , i])
        }else{
          unmetNeed <- unmetNeed + sum(hivpopu15[,, age - 4, , i])
        }
      }
      
      onFLART <- sum(artpopu5) + sum(artpopu15)
      needForFLART <- unmetNeed + onFLART   
      
      if(fp$artpaed_isperc[i - 1]){
        ARTlastYear <- needForFLART * fp$artpaed_num[i - 1] / 100
      }else{
        ARTlastYear <- fp$artpaed_num[i - 1]
      }
      if(fp$artpaed_isperc[i]){
        ARTthisYear <- needForFLART * fp$artpaed_num[i] / 100
      }else{
        ARTthisYear <- fp$artpaed_num[i]
      }
      
      newFLART <-((ARTthisYear + ARTlastYear) / 2 ) - onFLART
      if(newFLART < 0){
        newFLART <- 0
      }
      #Increase number starting ART to account for those who will die in the first year
      #They will be exposed to 1/2 year of mortality risk
      v1 <- sum(artpopu5[,,,,i] * fp$art_mort_u5)
      v2 <- sum(artpopu15[,,,,i] * fp$art_mort_u15)
      onARTDeaths <- (v1 + v2) / 2
      
      newFLART <- newFLART + onARTDeaths
      if(needForFLART < (onFLART + newFLART)){
        needForFLART <- onFLART + newFLART
      }
      #Distribute according to IeDEA data
      temp <- 0
      for(age in 0:4){
        if((age + 1) * 12 > childEligibilityAge[i]){
          CD4elig <- as.character(min(fp$paed_arteligibility[fp$paed_arteligibility$year == (proj_start + i - 1) & fp$paed_arteligibility$age_start <= age, 'cd4_pct_thresh']))
          temp <- temp + (sum(hivpopu5[,u5.elig.groups[[CD4elig]]:7, age + 1, , i]) * fp$paed_artdist[i, age +1])
        }else{
          temp <- temp + (sum(hivpopu5[,, age + 1, , i]) * fp$paed_artdist[i, age +1])
        }
      }
      for(age in 5:14){
        if((age + 1) * 12 > childEligibilityAge[i]){
          CD4elig <- as.character(min(fp$paed_arteligibility[fp$paed_arteligibility$year == (proj_start + i - 1) & fp$paed_arteligibility$age_start <= age, 'cd4_count_thresh']))
          temp <- temp + (sum(hivpopu15[,u15.elig.groups[[CD4elig]]:6, age - 4, , i]) * fp$paed_artdist[i, age +1])
        }else{
          temp <- temp + (sum(hivpopu15[,, age - 4, , i]) * fp$paed_artdist[i, age +1])
        }
      }
      adj <- ifelse(temp > 0, newFLART/temp, 1)
      ## Actually distribute ART
      newChildARTu5 <- array(0, c(7, 5, 2))
      newChildARTu15 <- array(0, c(6, 10, 2))
      for(age in 0:4){
        if((age + 1) * 12 > childEligibilityAge[i]){
          CD4elig <- as.character(min(fp$paed_arteligibility[fp$paed_arteligibility$year == (proj_start + i - 1) & fp$paed_arteligibility$age_start <= age, 'cd4_pct_thresh']))
          dim <- ifelse(u5.elig.groups[[CD4elig]] == 7, 2, 3)
          newChildARTu5[u5.elig.groups[[CD4elig]]:7,age + 1,] <- adj * apply(hivpopu5[,u5.elig.groups[[CD4elig]]:7, age + 1, , i], 2:dim, sum) * fp$paed_artdist[i, age +1]
        }else{
          newChildARTu5[,age + 1,] <- adj * apply(hivpopu5[,, age + 1, , i], 2:3, sum) * fp$paed_artdist[i, age +1]
        }
      }
      for(age in 5:14){
        if((age + 1) * 12 > childEligibilityAge[i]){
          CD4elig <- as.character(min(fp$paed_arteligibility[fp$paed_arteligibility$year == (proj_start + i - 1) & fp$paed_arteligibility$age_start <= age, 'cd4_count_thresh']))
          dim <- ifelse(u15.elig.groups[[CD4elig]] == 6, 2, 3)
          newChildARTu15[u15.elig.groups[[CD4elig]]:6,age - 4,] <- adj * apply(hivpopu15[,u15.elig.groups[[CD4elig]]:6, age - 4, , i], 2:dim, sum) * fp$paed_artdist[i, age +1]
        }else{
          newChildARTu15[,age - 4,] <- adj * apply(hivpopu15[,, age - 4, , i], 2:3, sum) * fp$paed_artdist[i, age +1]
        }
      }
      
      ## Get CTX coverage
      ## off-ART u5 plus eligible for ART 5-15
      posu5pop <- sum(hivpopu5[,,,,i])
      eligible5to15 <- 0
      for(age in 5:14){
        if((age + 1) * 12 > childEligibilityAge[i]){
          CD4elig <- as.character(min(fp$paed_arteligibility[fp$paed_arteligibility$year == (proj_start + i - 1) & fp$paed_arteligibility$age_start <= age, 'cd4_count_thresh']))
          eligible5to15 <- eligible5to15 + sum(hivpopu15[,u15.elig.groups[[CD4elig]]:6, age - 4, , i]) 
        }else{
          eligible5to15 <- eligible5to15 + sum(hivpopu15[,u15.elig.groups[[CD4elig]]:6, age - 4, , i])
        }
      }
      needCTX <- posu5pop + eligible5to15
      if(needCTX > 0){
        if(fp$cotrim_isperc[i]){
          CTXcoverage <- fp$cotrim_num[i]
        }else{
          CTXcoverage <- min(1, fp$cotrim_num[i]/needCTX)
        }
      }else{
        CTXcoverage <- 0
      }
            
      calcBFtransmissions <- function(m1, m2, i){
        BFTR <- 0
        perc.optA <- treat.opt[['postnat_optionA']]
        perc.optB <- treat.opt[['postnat_optionB']]
        optA.trans.rate <- (fp$MTCtrans[fp$MTCtrans$regimen == 'Option_A', 'breastfeeding_gt350cd4']/100)
        optB.trans.rate <- (fp$MTCtrans[fp$MTCtrans$regimen == 'Option_B', 'breastfeeding_gt350cd4']/100)
        dropout.optA <- fp$pmtct_dropout[fp$pmtct_dropout$year == (proj_start + i - 1), 'mth_drop_rt_optionA']
        dropout.optB <- fp$pmtct_dropout[fp$pmtct_dropout$year == (proj_start + i - 1), 'mth_drop_rt_optionB']
        if(propgt350 > 0){
          if(perc.optA +perc.optB - treat.opt[['tripleARTbefPreg']] -treat.opt[['tripleARTdurPreg']] > propgt350){
            excess <- perc.optA + perc.optB - treat.opt[['tripleARTbefPreg']] -treat.opt[['tripleARTdurPreg']] - propgt350
            optA.trans.rate <- (propgt350 * optA.trans.rate + excess * 1.45 / 0.46 * optA.trans.rate) / (propgt350 + excess)
            optB.trans.rate <- (propgt350 * optB.trans.rate + excess * 1.45 / 0.46 * optB.trans.rate) / (propgt350 + excess)
          }
        }
        for(d in m1:m2){
          perc.optA <- treat.opt[['postnat_optionA']] / exp(d * 2 * log(1 + dropout.optA / 100))
          perc.optB <- treat.opt[['postnat_optionB']] / exp(d * 2 * log(1 + dropout.optB / 100))
          perc.noproph <- 1 - perc.optA - perc.optB - treat.opt[['tripleARTdurPreg']] - treat.opt[['tripleARTbefPreg']]
          if(perc.noproph < 0){perc.noproph <- 0}
          percentInProgram <- 1 - treat.opt[['no_proph']]
          ## no prophylaxis
          BFTR <- BFTR + (((1 - fp$perc_bf_on_art[d] / 100) * (1 - percentInProgram) + 
                         (1 - fp$perc_bf_off_art[d]/100) * percentInProgram) * perc.noproph *
                         ((proplt200 + prop200to350) * (fp$MTCtrans[fp$MTCtrans$regimen == 'no_prophylaxis' & fp$MTCtrans$definition == 'exisiting_LT200CD4','breastfeeding_lt350cd4']/100) +
                          propgt350 * (fp$MTCtrans[fp$MTCtrans$regimen == 'no_prophylaxis' & fp$MTCtrans$definition == 'exisiting_GT340CD4', 'breastfeeding_gt350cd4']/100)))
        ## option A
        BFTR <- BFTR + (1 - (fp$perc_bf_on_art[d]/100)) * perc.optA * optA.trans.rate
        ## optionB
        BFTR <- BFTR + (1 - (fp$perc_bf_on_art[d]/100)) * perc.optB * optB.trans.rate
        ## triple art
        artp.lastyr.byage <- ifelse(i > 2, rowMeans(artpop[,,h.fert.idx, f.idx,i-2:1],,3), artp.byage)
        prop.new.art <- ifelse(sum(artp.byage) == 0, 0, (sum(artp.byage) - sum(artp.lastyr.byage)) / sum(artp.byage))
        BFTR <- BFTR + ((1 - fp$perc_bf_on_art[d]/100) * treat.opt[['tripleARTbefPreg']] * ((1-prop.new.art) * fp$MTCtrans[fp$MTCtrans$regimen == 'ART' & fp$MTCtrans$definition == 'start_pre_preg','breastfeeding_lt350cd4']/100 +
                          prop.new.art * fp$MTCtrans[fp$MTCtrans$regimen == 'ART' & fp$MTCtrans$definition == 'start_dur_preg','breastfeeding_lt350cd4']))
        BFTR <- BFTR + ((1 - fp$perc_bf_on_art[d] / 100) * treat.opt[['tripleARTdurPreg']] *
                          (prop.new.art * fp$MTCtrans[fp$MTCtrans$regimen == 'ART' & fp$MTCtrans$definition == 'start_pre_preg','breastfeeding_lt350cd4'] +
                             (1 - prop.new.art) * fp$MTCtrans[fp$MTCtrans$regimen == 'ART' & fp$MTCtrans$definition == 'start_dur_preg','breastfeeding_lt350cd4']))
        }
        return(BFTR * 2)
      }
     }
      
    
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
