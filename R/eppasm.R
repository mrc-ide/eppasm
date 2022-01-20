
#' @useDynLib eppasm eppasmC
#' @export
simmod.specfp <- function(fp, VERSION="C"){
  
  if(!exists("popadjust", where=fp))
    fp$popadjust <- FALSE
  
  if(!exists("incidmod", where=fp))
    fp$incidmod <- "eppspectrum"
  
  if(exists("turnover",where = fp))
    turnover <- fp$turnover
  
  if(VERSION != "R"){
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
  
  ## Attach state space variables
  invisible(list2env(fp$ss, environment())) # put ss variables in environment for convenience
  
  birthslag <- fp$birthslag
  pregprevlag <- rep(0, PROJ_YEARS)
  
  ## initialize projection
  pop <- array(0, c(pAG, NG, pDS, PROJ_YEARS))
  pop[,,1,1] <- fp$basepop
  entries <- fp$targetpop
  hivpop <- array(0, c(hDS, hAG, NG, PROJ_YEARS))
  artpop <- array(0, c(hTS, hDS, hAG, NG, PROJ_YEARS))
 
  ## initialize output
  prev15to49 <- numeric(PROJ_YEARS)
  incid15to49 <- numeric(PROJ_YEARS)
  sexinc15to49out <- array(NA, c(NG, PROJ_YEARS))
  paedsurvout <- rep(NA, PROJ_YEARS)
  art15plus_num = fp$art15plus_num
  
  infections <- array(0, c(pAG, NG, PROJ_YEARS))
  hivdeaths <- array(0, c(pAG, NG, PROJ_YEARS))
  natdeaths <- array(0, c(pAG, NG, PROJ_YEARS))
  
  if(turnover){

    duration = fp$duration
    t.rate = 1/duration
    t.pop <-  array(0, c(pAG, NG, pDS, PROJ_YEARS)) 
    t.pop[,,1,1] <- 0 ##No one in first year turnover population
    t.entrants <- t.pop
    t.entrants[,,1,1] <- t.pop[,,1,1] * t.rate * fp$ss$DT 
    
    t.hivpop <- array(0, c(hDS, hAG, NG, PROJ_YEARS))
    t.artpop <- array(0, c(hTS, hDS, hAG, NG, PROJ_YEARS))
    
    t.prev15to49 <- numeric(PROJ_YEARS)
    t.hivdeaths <- array(0, c(pAG, NG, PROJ_YEARS)) 
    t.natdeaths <- array(0, c(pAG, NG, PROJ_YEARS))
    
    t.art15plus_num = art15plus_num*t.rate
    art15plus_num = fp$art15plus_num - t.art15plus_num
    
  }
  
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
  
  ## store some outputs to troubleshoot
  #y = 52
  #out_dat = data.frame(year = rep(NA,y*11),
                       # time_step = rep(c(0:10),y),
                       # total_pop = rep(NA,y*11),
                       # t_pop = rep(NA,y*11),
                       # total_pop_entries = rep(NA,y*11),
                       # total_pop_deaths = rep(NA,y*11),
                       # total_pop_mig = rep(NA,y*11),
                       # t_pop_entries = rep(NA,y*11),
                       # t_pop_deaths = rep(NA,y*11),
                       # total_hiv = rep(NA,y*11),
                       # t_hiv = rep(NA,y*11),
                       # hiv_migration = rep(NA,y*11),
                       # t_hiv_migration = rep(NA,y*11),
                       # hiv_deaths = rep(NA,y*11),
                       # t_hivdeaths= rep(NA,y*11),
                       # hiv_p_entries = rep(NA,y*11),
                       # artpop = rep(NA,y*11),
                       # infections = rep(NA, y*11),
                       # t_artpop = rep(NA,y*11))
                       # 
  
  for(i in 2:fp$SIM_YEARS){
    #print(i)
    #if(length(pop[,,,i][pop[,,,i] < 0]) > 0) ##print(paste0("year: ",i," has neg numbers"))
    #if(i == 9) break
    #out_dat[-21 + 11*i,"year"] <- i
    ## ################################### ##
    ##  Single-year population projection  ##
    ## ################################### ##
    ###print(i)
    ## age the population
    pop[-c(1,pAG),,,i] <- pop[-(pAG-1:0),,,i-1]
    pop[pAG,,,i] <- pop[pAG,,,i-1] + pop[pAG-1,,,i-1] # open age group
    #out_dat[-21 + 11*i,"total_pop"] <- sum(pop[,,,i])
    
    ## Add lagged births into youngest age group
    entrant_prev <- fp$entrantprev[,i]
    
    
    if(turnover) { 
      t.pop[-c(1,pAG),,,i] <- t.pop[-(pAG-1:0),,,i-1]
      t.pop[pAG,,,i] <- t.pop[pAG,,,i-1] + t.pop[pAG-1,,,i-1] # open age group
      #out_dat[-21 + 11*i,"t_pop"] <- sum(t.pop[,,,i])
    }
    
    if(exists("popadjust", where=fp) & fp$popadjust & !turnover){
      hivn_entrants <- fp$entrantpop[,i-1]*(1-entrant_prev)
      hivp_entrants <- fp$entrantpop[,i-1]*entrant_prev
    } else if(exists("popadjust", where=fp) & fp$popadjust & turnover){ ##Full age distribution of entries if its a key population with turnover
      
      entries_total <- entries[,,i]
      
      hivn_entrants <- entries_total
      hivp_entrants <- entries_total
      
      hivn_entrants[1,] <- entries_total[1,]*(1-entrant_prev)
      hivp_entrants[1,] <- entries_total[1,]*entrant_prev
      
      entrant_prev <- 0 #Assume entrant prevalence = 0, except for 15 year olds (pediatric prev)
      hivn_entrants[-1,] <- entries_total[-1,]*(1-entrant_prev)
      hivp_entrants[-1,] <- entries_total[-1,]*entrant_prev
      
    } else {
      hivn_entrants <- birthslag[,i-1]*fp$cumsurv[,i-1]*(1-entrant_prev / fp$paedsurv_lag[i-1]) + fp$cumnetmigr[,i-1]*(1-pregprevlag[i-1]*fp$netmig_hivprob)
      hivp_entrants <- birthslag[,i-1]*fp$cumsurv[,i-1]*entrant_prev + fp$cumnetmigr[,i-1]*entrant_prev
    } 

    entrant_prev_out[i] <- sum(hivp_entrants) / sum(hivn_entrants+hivp_entrants)
    hivp_entrants_out[,i] <- sum(hivp_entrants)
    
    if(!is.null(dim(hivn_entrants))){
      pop[,,hivn.idx,i] <- pop[,,hivn.idx,i]+hivn_entrants
      pop[,,hivp.idx,i] <- pop[,,hivp.idx,i]+hivp_entrants
      if(turnover){ #Only 15 year olds 
        t.hivp_entrants <- hivn_entrants[1,]*t.rate
        t.hivn_entrants <- hivp_entrants[1,]*t.rate
      }
      
    } else {
      pop[1,,hivn.idx,i] <- pop[1,,hivn.idx,i]+hivn_entrants
      pop[1,,hivp.idx,i] <- pop[1,,hivp.idx,i]+hivp_entrants
      if(turnover){
        t.hivp_entrants <- hivn_entrants*t.rate
        t.hivn_entrants <- hivp_entrants*t.rate
      }
      
    }
    
    pop[is.nan(pop)] <- 0
    t.pop[is.nan(t.pop)] <- 0
  
    #print(pop[,2,hivp.idx,i])
    ##Probability of moving from 1 HIV age group to anther.
    #Proportion of 20-24 that would move into 25 (need to determine the proportion of 24 year olds)
    ###Assumes CD4 dist of 20-24 is the same.
    hiv.ag.prob <- pop[aglast.idx,,hivp.idx,i-1] / apply(pop[,,hivp.idx,i-1], 2, ctapply, ag.idx, sum)
    hiv.ag.prob[is.nan(hiv.ag.prob)] <- 0
    
    if(turnover){ 
      t.hiv.ag.prob <- t.pop[aglast.idx,,hivp.idx,i-1] / apply(t.pop[,,hivp.idx,i-1], 2, ctapply, ag.idx, sum)
      t.hiv.ag.prob[is.nan(t.hiv.ag.prob)] <- 0
    }
    
    hivpop[,,,i] <- hivpop[,,,i-1] 
    hivpop[,-hAG,,i] <- hivpop[,-hAG,,i] - sweep(hivpop[,-hAG,,i-1], 2:3, hiv.ag.prob[-hAG,], "*")
    hivpop[,-1,,i] <- hivpop[,-1,,i] + sweep(hivpop[,-hAG,,i-1], 2:3, hiv.ag.prob[-hAG,], "*")
    hivpop[,1,,i] <- hivpop[,1,,i] + sweep(fp$paedsurv_cd4dist[,,i], 2, pop[1,,hivp.idx,i] * (1-fp$entrantartcov[,i]), "*") #Age 15 follows pediatric CD4 distribution
    #print(sum(hivpop[,,,i]))
    #out_dat[-21 + 11*i,"total_hiv"] <- sum(hivpop[,,,i])
    
    if(turnover){ 
      t.hivpop[,,,i] <- t.hivpop[,,,i-1]
      t.hivpop[,-hAG,,i] <- t.hivpop[,-hAG,,i] - sweep(t.hivpop[,-hAG,,i-1], 2:3, t.hiv.ag.prob[-hAG,], "*")
      t.hivpop[,-1,,i] <- t.hivpop[,-1,,i] + sweep(t.hivpop[,-hAG,,i-1], 2:3, t.hiv.ag.prob[-hAG,], "*")
      t.hivpop[,1,,i] <- t.hivpop[,1,,i] + sweep(fp$paedsurv_cd4dist[,,i], 2, t.pop[1,,hivp.idx,i]  * (1-fp$entrantartcov[,i]), "*")
      #out_dat[-21 + 11*i,"t_hiv"] <- sum(t.hivpop[,,,i])
    }
    

    if(i > fp$tARTstart){ 
      artpop[,,,,i] <- artpop[,,,,i-1] 
      artpop[,,-hAG,,i] <- artpop[,,-hAG,,i] - sweep(artpop[,,-hAG,,i-1], 3:4, hiv.ag.prob[-hAG,], "*")
      artpop[,,-1,,i] <- artpop[,,-1,,i] + sweep(artpop[,,-hAG,,i-1], 3:4, hiv.ag.prob[-hAG,], "*")
      artpop[,,1,,i] <- artpop[,,1,,i] + sweep(fp$paedsurv_artcd4dist[,,,i], 3, hivp_entrants[1,] * fp$entrantartcov[,i], "*")
      if(turnover){ 
        t.artpop[,,,,i] <- t.artpop[,,,,i-1] 
        t.artpop[,,-hAG,,i] <- t.artpop[,,-hAG,,i] - sweep(t.artpop[,,-hAG,,i-1], 3:4, t.hiv.ag.prob[-hAG,], "*")
        t.artpop[,,-1,,i] <- t.artpop[,,-1,,i] + sweep(t.artpop[,,-hAG,,i-1], 3:4, t.hiv.ag.prob[-hAG,], "*")
        t.artpop[,,1,,i] <- t.artpop[,,1,,i] + sweep(fp$paedsurv_artcd4dist[,,,i], 3, t.hivp_entrants * fp$entrantartcov[,i], "*")
        
      }
    }
    
    ## survive the population
    if(turnover){ #sx should be the same regardless
      t.deaths <- sweep(t.pop[,,,i], 1:2, (1-fp$Sx[,,i]), "*")
      t.hiv.sx.prob <- 1-apply(t.deaths[,,2], 2, ctapply, ag.idx, sum) / apply(t.pop[,,2,i], 2, ctapply, ag.idx, sum)
      t.hiv.sx.prob[is.nan(t.hiv.sx.prob)] <- 0
      t.pop[,,,i] <- t.pop[,,,i] - t.deaths
      t.natdeaths[,,i] <- rowSums(t.deaths,,2)
      t.hivpop[,,,i] <- sweep(t.hivpop[,,,i], 2:3, t.hiv.sx.prob, "*") 
      #out_dat[-21 + 11*i,"t_pop_deaths"] <- sum(t.deaths)
    }
    
    deaths <- sweep(pop[,,,i], 1:2, (1-fp$Sx[,,i]), "*") 
    hiv.sx.prob <- 1-apply(deaths[,,2], 2, ctapply, ag.idx, sum) / apply(pop[,,2,i], 2, ctapply, ag.idx, sum)
    hiv.sx.prob[is.nan(hiv.sx.prob)] <- 0
    pop[,,,i] <- pop[,,,i] - deaths
    natdeaths[,,i] <- rowSums(deaths,,2)
    hivpop[,,,i] <- sweep(hivpop[,,,i], 2:3, hiv.sx.prob, "*")
    #out_dat[-21 + 11*i,"total_pop_deaths"] <- sum(deaths)
    hivpop[is.nan(hivpop)] <- 0
    t.hivpop[is.nan(t.hivpop)] <- 0
    
    if(i > fp$tARTstart){
      artpop[,,,,i] <- sweep(artpop[,,,,i], 3:4, hiv.sx.prob, "*")
      if(turnover){ 
        t.artpop[,,,,i] <- sweep(t.artpop[,,,,i], 3:4, hiv.sx.prob, "*")
      }
    }
    
    ## net migration
    if(turnover){
      t.netmigsurv <- fp$netmigr[,,i]*t.rate*(1+fp$Sx[,,i])/2 ##Need to assume something here..
      t.mr.prob <- 1+t.netmigsurv / rowSums(t.pop[,,,i],,2)
      t.mr.prob[is.nan(t.mr.prob)] <- 0 #Problem due to 0 population in first year of t pop
      t.mr.prob[t.mr.prob < 0] <- 1 ##Temporary fix for too many migrants
      t.hiv.mr.prob <- apply(t.mr.prob * t.pop[,,2,i], 2, ctapply, ag.idx, sum) /  apply(t.pop[,,2,i], 2, ctapply, ag.idx, sum)
      t.hiv.mr.prob[is.nan(t.hiv.mr.prob)] <- 0
      #out_dat[-21 + 11*i,"t_pop_mig"] <- sum(sweep(t.pop[,,,i], 1:2, t.mr.prob, "*")) - sum(t.pop[,,,i])
      t.pop[,,,i] <- sweep(t.pop[,,,i], 1:2, t.mr.prob, "*")
    }

    netmigsurv <- (fp$netmigr[,,i]-fp$netmigr[,,i]*t.rate)*(1+fp$Sx[,,i])/2
    mr.prob <- 1+netmigsurv / rowSums(pop[,,,i],,2)
    mr.prob[is.nan(mr.prob)] <- 0
    mr.prob[mr.prob < 0] <- 1 ##Temporary fix for too many migrants
    hiv.mr.prob <- apply(mr.prob * pop[,,2,i], 2, ctapply, ag.idx, sum) /  apply(pop[,,2,i], 2, ctapply, ag.idx, sum)
    hiv.mr.prob[is.nan(hiv.mr.prob)] <- 0
    #out_dat[-21 + 11*i,"total_pop_mig"] <- sum(sweep(pop[,,,i], 1:2, mr.prob, "*")) - sum(pop[,,,i])
    pop[,,,i] <- sweep(pop[,,,i], 1:2, mr.prob, "*")

    t.hivpop[is.nan(t.hivpop)] <- 0
    hivpop[is.nan(hivpop)] <- 0
    
    #out_dat[-21 + 11*i,"hiv_migration"] <- sum(sweep(hivpop[,,,i], 2:3, hiv.mr.prob, "*")) - sum(hivpop[,,,i])
    hivpop[,,,i] <- sweep(hivpop[,,,i], 2:3, hiv.mr.prob, "*")
    if(i > fp$tARTstart)
      artpop[,,,,i] <- sweep(artpop[,,,,i], 3:4, hiv.mr.prob, "*")

    if(turnover){
      #out_dat[-21 + 11*i,"t_hiv_migration"] <- sum(sweep(t.hivpop[,,,i], 2:3, t.hiv.mr.prob, "*")) - sum(t.hivpop[,,,i])
      t.hivpop[,,,i] <- sweep(t.hivpop[,,,i], 2:3, t.hiv.mr.prob, "*")
      if(i > fp$tARTstart)
        t.artpop[,,,,i] <- sweep(t.artpop[,,,,i], 3:4, t.hiv.mr.prob, "*")
    }
    
    ## fertility 
    births.by.age <- rowSums(pop[p.fert.idx, f.idx,,i-1:0])/2 * fp$asfr[,i]
    births.by.h.age <- ctapply(births.by.age, ag.idx[p.fert.idx], sum)
    births <- fp$srb[,i] * sum(births.by.h.age)
    if(i+AGE_START <= PROJ_YEARS)
      birthslag[,i+AGE_START-1] <- births
    
    t.pop[is.nan(t.pop)] <- 0 ##Fix for 0 population, e.g. males in female sex workers
    pop[is.nan(pop)] <- 0
    hivpop[is.nan(hivpop)] <- 0
    artpop[is.nan(artpop)] <- 0
    t.hivpop[is.nan(t.hivpop)] <- 0
    t.artpop[is.nan(t.artpop)] <- 0
    
    #if(length(pop[,,,i][pop[,,,i] < 0]) > 0) #print(paste0("year: ",i," has neg numbers"))
    
    ## ########################## ##
    ##  Disease model simulation  ##
    ## ########################## ##
    ####print(paste0("post deaths + migration ",sum(pop[p.age15to49.idx,,,i])))
    ## events at dt timestep
    for(ii in seq_len(hiv_steps_per_year)){
      #for(ii in 1){
      #if(ii == 5) break
      #print(ii)
      ts <- (i-2)/DT + ii
      #out_dat[(-21 + 11*i) + ii,"year"] <- i
      
      grad <- array(0, c(hDS, hAG, NG))
      
      if(turnover) {
        calc.agdist <- function(x) {d <- x/rep(ctapply(x, ag.idx, sum), h.ag.span); d[is.na(d)] <- 0; d}
        
        t.grad <- grad 
        t.rate.ts = (1/hiv_steps_per_year)/duration #duration * timestep
        
        ##Turnover entrants at this time step
        t.entrants[,,,i] =  (pop[,,,i] * t.rate.ts) 
        t.pop[,,,i] <- t.pop[,,,i] + t.entrants[,,,i]
        t.hivpop[,,,i] <- t.hivpop[,,,i] + (hivpop[,,,i] * t.rate.ts)
        t.artpop[,,,,i] <- t.artpop[,,,,i] + (artpop[,,,,i] * t.rate.ts)

        ##Assign HIV+ turnover entrants to CD4 categories
        t.grad <- sweep(fp$cd4_initdist, 2:3, apply(t.entrants[,,hivp.idx,i], 2, ctapply, ag.idx, sum), "*")
        #out_dat[(-21 + 11*i) + ii,"t_hiv_p_entries"] <- sum(t.grad)
        #out_dat[(-21 + 11*i) + ii,"t_pop_entries"] <- sum(t.entrants[,,,i])
        
        ##Progression and HIV deaths
        t.grad[-hDS,,] <- t.grad[-hDS,,] - fp$cd4_prog * t.hivpop[-hDS,,,i]  # remove cd4 stage progression (untreated)
        t.grad[-1,,] <- t.grad[-1,,] + fp$cd4_prog * t.hivpop[-hDS,,,i]      # add cd4 stage progression (untreated)
        
        if(fp$scale_cd4_mort == 1){
          t.cd4mx_scale <- t.hivpop[,,,i] / (t.hivpop[,,,i] + colSums(t.artpop[,,,,i]))
          t.cd4mx_scale[!is.finite(t.cd4mx_scale)] <- 1.0
          t.cd4_mort_ts <- fp$cd4_mort * t.cd4mx_scale
        } else
          t.cd4_mort_ts <- fp$cd4_mort
        
        t.grad <- t.grad - t.cd4_mort_ts * t.hivpop[,,,i]          # HIV mortality, untreated
        t.hivdeaths.ts <- DT*(colSums(t.cd4_mort_ts * t.hivpop[,,,i]) + colSums(fp$art_mort * fp$artmx_timerr[ , i] * t.artpop[,,,,i],,2))
        t.hivdeaths_p.ts <- apply(t.hivdeaths.ts, 2, rep, h.ag.span) * apply(t.pop[,,hivp.idx,i], 2, calc.agdist)  # HIV deaths by single-year age
        t.pop[,,2,i] <- t.pop[,,2,i] - t.hivdeaths_p.ts
        t.hivdeaths[,,i] <- t.hivdeaths[,,i] + t.hivdeaths_p.ts
        #out_dat[(-21 + 11*i) + ii,"t_hivdeaths"] <- sum(t.hivdeaths[,,i])
        
        #Remove deaths from hivpop array
        t.hivpop[,,,i] <- t.hivpop[,,,i] - (DT*t.cd4_mort_ts * t.hivpop[,,,i])
        
        #Remove turnover entrants from key population
        pop[,,,i] <- pop[,,,i] - t.entrants[,,,i]
        hivpop[,,,i] <- hivpop[,,,i] - (hivpop[,,,i] * t.rate.ts)
        artpop[,,,,i] <- artpop[,,,,i] - (artpop[,,,,i] * t.rate.ts)
        
        #out_dat[(-21 + 11*i) + ii,"total_pop"] <- sum(pop[,,,i])
        
      }
      
      if(fp$eppmod != "directincid"){
        ## incidence
        
        ## calculate r(t)
        if(fp$eppmod %in% c("rtrend", "rtrend_rw"))
          rvec[ts] <- calc_rtrend_rt(fp$proj.steps[ts], fp, rvec[ts-1], prevlast, pop, i, ii)
        else
          rvec[ts] <- fp$rvec[ts]
        
        infections.ts <- calc_infections_eppspectrum(fp, pop, hivpop, artpop, i, ii, rvec[ts])
        
        incrate15to49.ts.out[ts] <- attr(infections.ts, "incrate15to49.ts")
        prev15to49.ts.out[ts] <- attr(infections.ts, "prevcurr")
        prevlast <- attr(infections.ts, "prevcurr")
        
        pop[,,hivn.idx,i] <- pop[,,hivn.idx,i] - DT*infections.ts
        pop[,,hivp.idx,i] <- pop[,,hivp.idx,i] + DT*infections.ts
        infections[,,i] <- infections[,,i] + DT*infections.ts
        #out_dat[(-21 + 11*i) + ii,"infections"] <- sum(infections[,,i])
        ###print(paste0("year: ",i," and time step", ii," ", length(pop[,,,i][pop[,,,i] < 0]) > 0))
        
        grad <- grad + sweep(fp$cd4_initdist, 2:3, apply(infections.ts, 2, ctapply, ag.idx, sum), "*")
        incid15to49[i] <- incid15to49[i] + sum(DT*infections.ts[p.age15to49.idx,])
        #out_dat[(-21 + 11*i) + ii,"hiv_p_entries"] <- sum(grad)
        
      }
      
      ##print(paste0("post add entrants ",sum(hivpop[,,,i]) == sum(pop[,,2,i]) - (sum(grad*DT))))
      ##print(paste0("post add t.entrants ",sum(t.hivpop[,,,i]) == sum(t.pop[,,2,i])))
      
      ## disease progression and mortality 
      grad[-hDS,,] <- grad[-hDS,,] - fp$cd4_prog * hivpop[-hDS,,,i]  # remove cd4 stage progression (untreated)
      grad[-1,,] <- grad[-1,,] + fp$cd4_prog * hivpop[-hDS,,,i]      # add cd4 stage progression (untreated)
      
      
      if(fp$scale_cd4_mort == 1){
        cd4mx_scale <- hivpop[,,,i] / (hivpop[,,,i] + colSums(artpop[,,,,i]))
        cd4mx_scale[!is.finite(cd4mx_scale)] <- 1.0
        cd4_mort_ts <- fp$cd4_mort * cd4mx_scale
      } else
        cd4_mort_ts <- fp$cd4_mort
      
      grad <- grad - cd4_mort_ts * hivpop[,,,i]              # HIV mortality, untreated
      ##print(paste0("post grad hiv deaths ",sum(hivpop[,,,i])  == (sum(pop[,,2,i])  - sum(grad + cd4_mort_ts * hivpop[,,,i])*DT)))
      
      ## Remove hivdeaths from pop
      hivdeaths.ts <- DT*(colSums(cd4_mort_ts * hivpop[,,,i]) + colSums(fp$art_mort * fp$artmx_timerr[ , i] * artpop[,,,,i],,2))
      calc.agdist <- function(x) {d <- x/rep(ctapply(x, ag.idx, sum), h.ag.span); d[is.na(d)] <- 0; d}
      hivdeaths_p.ts <- apply(hivdeaths.ts, 2, rep, h.ag.span) * apply(pop[,,hivp.idx,i], 2, calc.agdist)  # HIV deaths by single-year age
      pop[,,2,i] <- pop[,,2,i] - hivdeaths_p.ts
      hivdeaths[,,i] <- hivdeaths[,,i] + hivdeaths_p.ts
      ##print(paste0("post HIV deaths ",sum(hivpop[,,,i]) == sum(hivdeaths_p.ts)  + sum(pop[,,2,i])  - sum(grad + cd4_mort_ts * hivpop[,,,i])*DT))
      ##print(paste0("HIV pop vs grad deaths ",sum(hivdeaths_p.ts)/DT == sum(cd4_mort_ts * hivpop[,,,i])))
      #out_dat[(-21 + 11*i) + ii,"hiv_deaths"] <- sum(hivdeaths[,,i])
      #if(length(pop[,,,i][pop[,,,i] < 0]) > 0) ##print(paste0("year: ",i," has neg numbers"))
      # ##print(paste0("post HIV deaths ",sum(pop[p.age15to49.idx,,,i])))
      # ##print(paste0("HIV deaths: ", sum(hivdeaths_p.ts)))
      
      ## ART initiation
      if(i >= fp$tARTstart) {
        
        gradART <- array(0, c(hTS, hDS, hAG, NG))
        
        if(turnover) t.gradART <- array(0, c(hTS, hDS, hAG, NG)) 
        
        ## progression and mortality
        gradART[1:(hTS-1),,,] <- gradART[1:(hTS-1),,,] - 1.0 / fp$ss$h_art_stage_dur * (artpop[1:(hTS-1),,,, i])  # remove ART duration progression 
        gradART[2:hTS,,,] <- gradART[2:hTS,,,] + 1.0 / fp$ss$h_art_stage_dur * (artpop[1:(hTS-1),,,, i])     # add ART duration progression
        
        gradART <- gradART - fp$art_mort * fp$artmx_timerr[ , i] * (artpop[,,,,i])  # ART mortality
        
        ## ART dropout
        ## remove proportion from all adult ART groups back to untreated pop
        grad <- grad + fp$art_dropout[i]*colSums(artpop[,,,,i])
        gradART <- gradART - fp$art_dropout[i]*(artpop[,,,,i])
        
        if(turnover){
          
          t.gradART[1:(hTS-1),,,] <- t.gradART[1:(hTS-1),,,] - 1.0 / fp$ss$h_art_stage_dur * t.artpop[1:(hTS-1),,,, i]      # remove ART duration progression 
          t.gradART[2:hTS,,,] <- t.gradART[2:hTS,,,] + 1.0 / fp$ss$h_art_stage_dur * t.artpop[1:(hTS-1),,,, i]      # add ART duration progression
          
          t.gradART <- t.gradART - fp$art_mort * fp$artmx_timerr[ , i] * t.artpop[,,,,i]   # ART mortality
          
          ## ART dropout
          ## remove proportion from all adult ART groups back to untreated pop
          t.grad <- t.grad + fp$art_dropout[i]*colSums(t.artpop[,,,,i])
          t.gradART <- t.gradART - fp$art_dropout[i]*t.artpop[,,,,i] #
        }
        
        ## calculate number eligible for ART
        artcd4_percelig <- 1 - (1-rep(0:1, times=c(fp$artcd4elig_idx[i]-1, hDS - fp$artcd4elig_idx[i]+1))) *
          (1-rep(c(0, fp$who34percelig), c(2, hDS-2))) *
          (1-rep(fp$specpop_percelig[i], hDS))
        
        art15plus.elig <- sweep(hivpop[,h.age15plus.idx,,i], 1, artcd4_percelig, "*")
        if(turnover) #Calc turnover pop ART eligibility before removing adding entrants
          t.art15plus.elig <- sweep(t.hivpop[,h.age15plus.idx,,i] + (hivpop[,,,i] * t.rate.ts), 1, artcd4_percelig, "*")
        
        ## calculate pregnant women - not calculating turnover population separate for now
        if(fp$pw_artelig[i]){
          births.dist <- sweep(fp$frr_cd4[,,i] * hivpop[,h.fert.idx,f.idx,i], 2,births.by.h.age / (ctapply(pop[p.fert.idx, f.idx, hivn.idx, i], ag.idx[p.fert.idx], sum) + colSums(fp$frr_cd4[,,i] * hivpop[,h.fert.idx,f.idx,i]) + colSums(fp$frr_art[,,,i] * artpop[ ,,h.fert.idx,f.idx,i],,2)), "*")
          if(fp$artcd4elig_idx[i] > 1)
            art15plus.elig[1:(fp$artcd4elig_idx[i]-1),h.fert.idx-min(h.age15plus.idx)+1,f.idx] <- art15plus.elig[1:(fp$artcd4elig_idx[i]-1),h.fert.idx-min(h.age15plus.idx)+1,f.idx] + births.dist[1:(fp$artcd4elig_idx[i]-1),]
        }
        
        ## calculate number to initiate ART based on number or percentage 
        artpop_curr_g <- colSums(artpop[,,h.age15plus.idx,,i],,3) + DT*colSums(gradART[,,h.age15plus.idx,],,3)
        artnum.ii <- c(0,0) # number on ART this ts
        
        if(turnover){
          t.artpop_curr_g <- colSums(t.artpop[,,h.age15plus.idx,,i],,3) + DT*colSums(t.gradART[,,h.age15plus.idx,],,3)
          t.artnum.ii <- c(0,0) 
        }
        
        if(DT*ii < 0.5){
          for(g in 1:2){
            if(!any(fp$art15plus_isperc[g,i-2:1])){  # both number
              artnum.ii[g] <- c(art15plus_num[g,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
              if(turnover) t.artnum.ii[g] <- c(t.art15plus_num[g,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
            } else if(all(fp$art15plus_isperc[g,i-2:1])){  # both percentage
              artcov.ii <- c(art15plus_num[g,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + artpop_curr_g[g])
              if(turnover){
                t.artcov.ii <- c(t.art15plus_num[g,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
                t.artnum.ii[g] <- t.artcov.ii * (sum(t.art15plus.elig[,,g]) + t.artpop_curr_g[g])
              }
            } else if(!fp$art15plus_isperc[g,i-2] & fp$art15plus_isperc[g,i-1]){ # transition number to percentage
              curr_coverage <- artpop_curr_g[g] / (sum(art15plus.elig[,,g]) + artpop_curr_g[g])
              artcov.ii <- curr_coverage + (art15plus_num[g,i-1] - curr_coverage) * DT/(0.5-DT*(ii-1))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + artpop_curr_g[g])
              if(turnover){
                t.curr_coverage <- t.artpop_curr_g[g] / (sum(t.art15plus.elig[,,g]) + t.artpop_curr_g[g])
                t.artcov.ii <- t.curr_coverage + (t.art15plus_num[g,i-1] - t.curr_coverage) * DT/(0.5-DT*(ii-1))
                t.artnum.ii[g] <- t.artcov.ii * (sum(t.art15plus.elig[,,g]) + t.artpop_curr_g[g])
              }
            }
          }
        } else {
          for(g in 1:2){
            if(!any(fp$art15plus_isperc[g,i-1:0])){  # both number
              artnum.ii[g] <- c(art15plus_num[g,i-1:0] %*% c(1-(DT*ii-0.5), DT*ii-0.5))
              if(turnover) t.artnum.ii[g] <- c(t.art15plus_num[g,i-1:0] %*% c(1-(DT*ii-0.5), DT*ii-0.5))
            } else if(all(fp$art15plus_isperc[g,i-1:0])) {  # both percentage
              artcov.ii <- c(art15plus_num[g,i-1:0] %*% c(1-(DT*ii-0.5), DT*ii-0.5))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + artpop_curr_g[g])
              if(turnover){
                t.artcov.ii <- c(t.art15plus_num[g,i-1:0] %*% c(1-(DT*ii-0.5), DT*ii-0.5))
                t.artnum.ii[g] <- t.artcov.ii * (sum(t.art15plus.elig[,,g]) + t.artpop_curr_g[g])
              }
            } else if(!fp$art15plus_isperc[g,i-1] & fp$art15plus_isperc[g,i]){  # transition number to percentage
              curr_coverage <- artpop_curr_g[g] / (sum(art15plus.elig[,,g]) + artpop_curr_g[g])
              artcov.ii <- curr_coverage + (art15plus_num[g,i] - curr_coverage) * DT/(1.5-DT*(ii-1))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + artpop_curr_g[g])
              if(turnover){
                t.curr_coverage <- t.artpop_curr_g[g] / (sum(t.art15plus.elig[,,g]) + t.artpop_curr_g[g])
                t.artcov.ii <- t.curr_coverage + (t.art15plus_num[g,i]- t.curr_coverage) * DT/(1.5-DT*(ii-1))
                t.artnum.ii[g] <- t.artcov.ii * (sum(t.art15plus.elig[,,g]) + t.artpop_curr_g[g])
              }
            }
          }
        }
        
        artpop_curr_g <- colSums(artpop[,,h.age15plus.idx,,i],,3) + DT*colSums(gradART[,,h.age15plus.idx,],,3)
        art15plus.inits <- pmax(artnum.ii - artpop_curr_g, 0)
        
        if(turnover){
          t.artpop_curr_g <- colSums(t.artpop[,,h.age15plus.idx,,i],,3) + DT*colSums(t.gradART[,,h.age15plus.idx,],,3)
          t.art15plus.inits <- pmax(t.artnum.ii - t.artpop_curr_g, 0)
        }
        
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
            
            if(turnover){
              t.artinit <- array(0, dim(t.art15plus.elig))
              t.remain_artalloc <- t.art15plus.inits
              for(m in hDS:1){
                t.elig_hm <- colSums(t.art15plus.elig[m,,])
                t.init_prop <- ifelse(t.elig_hm == 0, t.elig_hm, pmin(1.0, t.remain_artalloc / t.elig_hm, na.rm=TRUE))
                t.artinit[m , , ] <- sweep(t.art15plus.elig[m,,], 2, t.init_prop, "*")
                t.remain_artalloc <- t.remain_artalloc - t.init_prop * t.elig_hm
              }
            }
            
          } else {
            
            expect.mort.weight <- sweep(fp$cd4_mort[, h.age15plus.idx,], 3,
                                        colSums(art15plus.elig * fp$cd4_mort[, h.age15plus.idx,],,2), "/")          
            artinit.weight <- sweep(fp$art_alloc_mxweight * expect.mort.weight, 3, (1 - fp$art_alloc_mxweight)/colSums(art15plus.elig,,2), "+")
            artinit <- pmin(sweep(artinit.weight * art15plus.elig, 3, art15plus.inits, "*"),
                            art15plus.elig)
            
            if(turnover){
              t.expect.mort.weight <- sweep(fp$cd4_mort[, h.age15plus.idx,], 3,
                                            colSums(t.art15plus.elig * fp$cd4_mort[, h.age15plus.idx,],,2), "/")          
              t.artinit.weight <- sweep(fp$art_alloc_mxweight * t.expect.mort.weight, 3, (1 - fp$art_alloc_mxweight)/colSums(t.art15plus.elig,,2), "+")
              t.artinit <- pmin(sweep(t.artinit.weight * t.art15plus.elig, 3, t.art15plus.inits, "*"),
                                t.art15plus.elig)
            }
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
          
          if(turnover) {
            t.elig_below <- colSums(t.art15plus.elig[medcd4_idx,,,drop=FALSE],,2) * medcat_propbelow
            if(medcd4_idx < hDS)
              t.elig_below <- t.elig_below + colSums(t.art15plus.elig[(medcd4_idx+1):hDS,,,drop=FALSE],,2)
            
            t.elig_above <- colSums(t.art15plus.elig[medcd4_idx,,,drop=FALSE],,2) * (1.0-medcat_propbelow)
            if(medcd4_idx > 1)
              t.elig_above <- t.elig_above + colSums(t.art15plus.elig[1:(medcd4_idx-1),,,drop=FALSE],,2)
            
            t.initprob_below <- pmin(t.art15plus.inits * 0.5 / elig_below, 1.0, na.rm=TRUE)
            t.initprob_above <- pmin(t.art15plus.inits * 0.5 / elig_above, 1.0, na.rm=TRUE)
            t.initprob_medcat <- t.initprob_below * medcat_propbelow + t.initprob_above * (1-medcat_propbelow)
            
            t.artinit <- array(0, dim=c(hDS, hAG, NG))
            
            if(medcd4_idx < hDS)
              t.artinit[(medcd4_idx+1):hDS,,] <- sweep(t.art15plus.elig[(medcd4_idx+1):hDS,,,drop=FALSE], 3, t.initprob_below, "*")
            t.artinit[medcd4_idx,,] <- sweep(t.art15plus.elig[medcd4_idx,,,drop=FALSE], 3, t.initprob_medcat, "*")
            if(medcd4_idx > 0)
              t.artinit[1:(medcd4_idx-1),,] <- sweep(t.art15plus.elig[1:(medcd4_idx-1),,,drop=FALSE], 3, t.initprob_above, "*")
          }
        }
        
        artinit <- pmin(artinit, hivpop[ , , , i] + (DT * grad))
        artinit[is.nan(artinit)] <- 0
        
        #out_dat[(-21 + 11*i) + ii,"artinit"] <- sum(artinit)
        #out_dat[(-21 + 11*i) + ii,"t_artinit"] <- sum(t.artinit)
        
        grad[ , h.age15plus.idx, ] <- grad[ , h.age15plus.idx, ] - artinit / DT
        gradART[1, , h.age15plus.idx, ] <- gradART[1, , h.age15plus.idx, ] + artinit / DT
        artpop[,,,, i] <- artpop[,,,, i] + DT * gradART 
        #out_dat[(-21 + 11*i) + ii,"artpop"] <- sum(artpop[,,,, i])
        
        if(turnover){
          t.artinit <- pmin(t.artinit, t.hivpop[ , , , i]) ##Do not need to add DT*t.grad here because done above
          t.artinit[is.nan(t.artinit)] <- 0
          t.grad[ , h.age15plus.idx, ] <- t.grad[ , h.age15plus.idx, ] - t.artinit / DT
          t.gradART[1, , h.age15plus.idx, ] <- t.gradART[1, , h.age15plus.idx, ] + t.artinit / DT
          t.artpop[,,,, i] <- t.artpop[,,,, i] + DT * t.gradART
          t.hivpop[,h.age15plus.idx,,i] <- t.hivpop[,h.age15plus.idx,,i] - t.artinit
          #out_dat[(-21 + 11*i) + ii,"t_artpop"] <- sum(t.artpop[,,,, i])
          
        }
      }
      
      
      hivpop[,,,i] <- hivpop[,,,i] + DT * grad  
      # #print(sum(pop[,,hivp.idx,i]))
      # #print(sum(hivpop[,,,i]))
      # #print(sum(t.pop[,,hivp.idx,i]))
      # #print(sum(t.hivpop[,,,i]))
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
      
      if(turnover){
        
        t.entrants[,,,i] <- pop[,,,i] * t.rate
        t.pop[,,,i] <- t.pop[,,,i] + t.entrants[,,,i]
        pop[,,,i] <- pop[,,,i] - t.entrants[,,,i]
        
        t.hivpop[,,,i] <- t.hivpop[,,,i] + hivpop[,,,i] * t.rate
        hivpop[,,,i] <- hivpop[,,,i] - hivpop[,,,i] * t.rate
        
        t.artpop[,,,i] <- t.artpop[,,,i] +  artpop[,,,i] * t.rate
        artpop[,,,i] <- artpop[,,,i] -  artpop[,,,i] * t.rate
       
      }
      
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
      
      if(turnover){
        t.popadj.prob <- array(0, c(pAG, NG, PROJ_YEARS))
        t.popadj.prob[,,i] <- fp$targetpop[,,i]*t.rate / rowSums(t.pop[,,,i],,2)
        t.hiv.popadj.prob <- apply(t.popadj.prob[,,i] * t.pop[,,2,i], 2, ctapply, ag.idx, sum) /  apply(t.pop[,,2,i], 2, ctapply, ag.idx, sum)
        t.hiv.popadj.prob[is.nan(t.hiv.popadj.prob)] <- 0
        
        t.pop[,,,i] <- sweep(t.pop[,,,i], 1:2, popadj.prob[,,i], "*")
        t.hivpop[,,,i] <- sweep(t.hivpop[,,,i], 2:3, t.hiv.popadj.prob, "*")
        if(i >= fp$tARTstart)
          t.artpop[,,,,i] <- sweep(t.artpop[,,,,i], 3:4, t.hiv.popadj.prob, "*")
      }
      
    }
    
    # ##print(sum(pop[,,hivp.idx,i]))
    # ##print(sum(hivpop[,,,i]))
    # ##print(sum(t.pop[,,hivp.idx,i]))
    # ##print(sum(t.hivpop[,,,i]))
    # 
    ## prevalence among pregnant women
    hivn.byage <- ctapply(rowMeans(pop[p.fert.idx, f.idx, hivn.idx,i-1:0]), ag.idx[p.fert.idx], sum)
    hivp.byage <- rowMeans(hivpop[,h.fert.idx, f.idx,i-1:0],,2)
    artp.byage <- rowMeans(artpop[,,h.fert.idx, f.idx,i-1:0],,3)
    pregprev <- sum(births.by.h.age * (1 - hivn.byage / (hivn.byage + colSums(fp$frr_cd4[,,i] * hivp.byage) + colSums(fp$frr_art[,,,i] * artp.byage,,2)))) / sum(births.by.age)
    if(i+AGE_START <= PROJ_YEARS)
      pregprevlag[i+AGE_START-1] <- pregprev
    
    ## prevalence and incidence 15 to 49
    pop[is.nan(pop)] <- 0
    t.pop[is.nan(t.pop)] <- 0
    prev15to49[i] <- sum(pop[p.age15to49.idx,,hivp.idx,i]) / sum(pop[p.age15to49.idx,,,i])
    incid15to49[i] <- sum(incid15to49[i]) / sum(pop[p.age15to49.idx,,hivn.idx,i-1])
    
    if(turnover){
      t.prev15to49[i] <- sum(t.pop[p.age15to49.idx,,hivp.idx,i]) / sum(t.pop[p.age15to49.idx,,,i])
    }
    
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
  attr(pop, "pop") <- apply(pop, 4, sum)
  
  if(fp$eppmod != "directincid"){
    attr(pop, "incrate15to49_ts") <- incrate15to49.ts.out
    attr(pop, "prev15to49_ts") <- prev15to49.ts.out
  }
  
  attr(pop, "entrantprev") <- entrant_prev_out
  attr(pop, "hivp_entrants") <- hivp_entrants_out
  
  if(turnover){
    attr(pop, "t.pop") <- apply(t.pop, 4, sum)
    attr(pop, "t.hivpop") <- t.hivpop
    attr(pop, "t.artpop") <- t.artpop
    attr(pop, "t.hivdeaths") <- t.hivdeaths
    attr(pop, "t.natdeaths") <- t.natdeaths
  }
  
  class(pop) <- "spec"
  return(pop)
}


#' Add dimnames to EPP-ASM model output
#'
#' @param mod output from `simmod()`
#' @param fp fixed parameters input to `simmod()`
#'
#' @return Input `mod` object with dimnames applied to arrays.
#' 
#'
#' @export
spec_add_dimnames <- function(mod, fp) {
  
  nm_pAG <- fp$ss$AGE_START + seq_len(fp$ss$pAG) - 1L
  nm_NG <- c("male", "female")
  nm_pDS <- c("negative", "positive")
  nm_years <- fp$ss$proj_start + seq_len(fp$ss$PROJ_YEARS) - 1L
  nm_hDS <- c(">500", "350-499", "250-349", "200-249", "100-199", "50-99", "<50")
  nm_hTS <- c("art0mos", "art6mos", "art1yr")
  nm_hAG <- c("15-16", "17-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50+")
  
  dn_pop <- list(age = nm_pAG, sex = nm_NG, hivstatus = nm_pDS, year = nm_years)
  dn_hivpop <- list(cd4stage = nm_hDS, age_coarse = nm_hAG, sex = nm_NG, year = nm_years)
  dn_artpop <- c(list(artdur = nm_hTS), dn_hivpop)
  
  dimnames(mod) <- dn_pop
  dimnames(attr(mod, "hivpop")) <- dn_hivpop
  dimnames(attr(mod, "artpop")) <- dn_artpop
  dimnames(attr(mod, "infections")) <- dn_pop[c("age", "sex", "year")]
  dimnames(attr(mod, "hivdeaths")) <- dn_pop[c("age", "sex", "year")]
  dimnames(attr(mod, "natdeaths")) <- dn_pop[c("age", "sex", "year")]
  dimnames(attr(mod, "aidsdeaths_noart")) <- dn_hivpop
  dimnames(attr(mod, "aidsdeaths_art")) <- dn_artpop
  dimnames(attr(mod, "popadjust")) <- dn_pop[c("age", "sex", "year")]
  dimnames(attr(mod, "artinit")) <- dn_hivpop
  
  names(attr(mod, "pregprevlag")) <- nm_years
  names(attr(mod, "incrate15to49_ts")) <- fp$proj.steps[-length(fp$proj.steps)]
  names(attr(mod, "prev15to49_ts")) <- fp$proj.steps[-length(fp$proj.steps)]
  names(attr(mod, "rvec_ts")) <- fp$proj.steps[-length(fp$proj.steps)]
  names(attr(mod, "prev15to49")) <- nm_years
  names(attr(mod, "pregprev")) <- nm_years
  names(attr(mod, "incid15to49")) <- nm_years
  names(attr(mod, "entrantprev")) <- nm_years
  
  mod
}
