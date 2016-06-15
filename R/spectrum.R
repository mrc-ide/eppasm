
spectrumR <- function(fp, VERSION="C"){


  if(requireNamespace("fastmatch", quietly = TRUE))
    ctapply <- fastmatch::ctapply
  else
    ctapply <- tapply
  
  ## Attach state space variables
  invisible(list2env(fp$ss, environment())) # put ss variables in environment for convenience

  birthslag <- fp$birthslag
  pregprevlag <- rep(0, PROJ_YEARS)

  ## initialize projection
  pop <- array(0, c(pAG, NG, pDS, PROJ_YEARS))
  pop[,,1,1] <- fp$basepop
  hivpop <- array(0, c(hTS, hDS, hAG, NG, PROJ_YEARS))

  ## initialize output
  incrate15to49 <- numeric(PROJ_YEARS)
  sexinc15to49out <- array(NA, c(NG, PROJ_YEARS))
  paedsurvout <- rep(NA, PROJ_YEARS)

  incrate15to49.ts.out <- rep(NA, length(fp$rvec))

  for(i in 2:PROJ_YEARS){

    ## ################################### ##
    ##  Single-year population projection  ##
    ## ################################### ##

    ## age the population
    pop[-c(1,pAG),,,i] <- pop[-(pAG-1:0),,,i-1]
    pop[pAG,,,i] <- pop[pAG,,,i-1] + pop[pAG-1,,,i-1] # open age group

    ## Add lagged births into youngest age group
    pop[1,,hivn.idx,i] <- birthslag[,i-1]*fp$cumsurv[,i-1]*(1-pregprevlag[i-1]*fp$vert.trans) + fp$cumnetmigr[,i-1]*(1-pregprevlag[i-1]*fp$netmig.hivprob)
    pop[1,,hivp.idx,i] <- birthslag[,i-1]*fp$cumsurv[,i-1]*pregprevlag[i-1]*fp$vert.trans*fp$paedsurv + fp$cumnetmigr[,i-1]*pregprevlag[i-1]*fp$netmig.hivprob*fp$netmighivsurv

    paedsurvout[i] <- sum(birthslag[,i-1]*fp$cumsurv[,i-1]*pregprevlag[i-1]*fp$vert.trans*fp$paedsurv)

    hiv.ag.prob <- pop[aglast.idx,,hivp.idx,i-1] / apply(pop[,,hivp.idx,i-1], 2, ctapply, ag.idx, sum)
    hiv.ag.prob[is.nan(hiv.ag.prob)] <- 0

    hivpop[,,,,i] <- hivpop[,,,,i-1]
    hivpop[,,-hAG,,i] <- hivpop[,,-hAG,,i] - sweep(hivpop[,,-hAG,,i-1], 3:4, hiv.ag.prob[-hAG,], "*")
    hivpop[,,-1,,i] <- hivpop[,,-1,,i] + sweep(hivpop[,,-hAG,,i-1], 3:4, hiv.ag.prob[-hAG,], "*")
    hivpop[1,,1,,i] <- hivpop[1,,1,,i] + fp$paedsurv.cd4dist %o% (birthslag[,i-1]*fp$cumsurv[,i-1]*pregprevlag[i-1]*fp$vert.trans*fp$paedsurv + fp$cumnetmigr[,i-1]*pregprevlag[i-1]*fp$netmig.hivprob*fp$netmighivsurv)

    ## survive the population
    deaths <- sweep(pop[,,,i], 1:2, (1-fp$Sx[,,i]), "*")
    hiv.sx.prob <- 1-apply(deaths[,,2], 2, ctapply, ag.idx, sum) / apply(pop[,,2,i], 2, ctapply, ag.idx, sum)
    hiv.sx.prob[is.nan(hiv.sx.prob)] <- 0
    pop[,,,i] <- pop[,,,i] - deaths

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

    hivdeaths <- array(0, c(hAG, NG))

    ## events at dt timestep
    for(ii in seq_len(1/DT)){
      grad <- array(0, c(hTS, hDS, hAG, NG))

      ## incidence
      ts <- (i-2)/DT + ii

      incrate15to49.ts <- fp$rvec[ts] * (sum(hivpop[1,,h.age15to49.idx,,i]) + fp$relinfectART*sum(hivpop[-1,,h.age15to49.idx,,i])) / sum(pop[p.age15to49.idx,,,i]) + fp$iota * (ts == fp$ts.epi.start)
      sexinc15to49.ts <- incrate15to49.ts*c(1, fp$inc.sexratio[i])*sum(pop[p.age15to49.idx,,hivn.idx,i])/(sum(pop[p.age15to49.idx,m.idx,hivn.idx,i]) + fp$inc.sexratio[i]*sum(pop[p.age15to49.idx, f.idx,hivn.idx,i]))
      agesex.inc <- sweep(fp$inc.agerr[,,i], 2, sexinc15to49.ts/(colSums(pop[p.age15to49.idx,,hivn.idx,i] * fp$inc.agerr[p.age15to49.idx,,i])/colSums(pop[p.age15to49.idx,,hivn.idx,i])), "*")
      infections.ts <- agesex.inc * pop[,,hivn.idx,i]

      incrate15to49.ts.out[ts] <- incrate15to49.ts

      pop[,,hivn.idx,i] <- pop[,,hivn.idx,i] - DT*infections.ts
      pop[,,hivp.idx,i] <- pop[,,hivp.idx,i] + DT*infections.ts

      grad[1,,,] <- grad[1,,,] + sweep(fp$cd4.initdist, 2:3, apply(infections.ts, 2, ctapply, ag.idx, sum), "*")
      

      ## disease progression and mortality
      grad[1,-hDS,,] <- grad[1,-hDS,,] - fp$cd4.prog * hivpop[1,-hDS,,,i]  # remove cd4 stage progression (untreated)
      grad[1,-1,,] <- grad[1,-1,,] + fp$cd4.prog * hivpop[1,-hDS,,,i]      # add cd4 stage progression (untreated)
      grad[2:3,,,] <- grad[2:3,,,] - 2.0 * hivpop[2:3,,,, i]               # remove ART duration progression (HARD CODED 6 months duration)
      grad[3:4,,,] <- grad[3:4,,,] + 2.0 * hivpop[2:3,,,, i]               # add ART duration progression (HARD CODED 6 months duration)

      grad[1,,,] <- grad[1,,,] - fp$cd4.mort * hivpop[1,,,,i]              # HIV mortality, untreated
      grad[-1,,,] <- grad[-1,,,] - fp$art.mort * hivpop[-1,,,,i]           # ART mortality

      ## Remove hivdeaths from pop
      hivdeaths.ts <- DT*(colSums(fp$cd4.mort * hivpop[1,,,,i]) + colSums(fp$art.mort * hivpop[-1,,,,i],,2))
      calc.agdist <- function(x) {d <- x/rep(ctapply(x, ag.idx, sum), h.ag.span); d[is.na(d)] <- 0; d}
      pop[,,2,i] <- pop[,,2,i] - apply(hivdeaths.ts, 2, rep, h.ag.span) * apply(pop[,,hivp.idx,i], 2, calc.agdist)
      hivdeaths <- hivdeaths + hivdeaths.ts

      hivpop[,,,,i] <- hivpop[,,,,i] + DT*grad

      ## ART initiation
      if(sum(fp$artnum15plus[,i])>0){
        if(DT*ii < 0.5){
          artnum.ii <- c(fp$artnum15plus[,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
        } else {
          artnum.ii <- c(fp$artnum15plus[,i-1:0] %*% c(1-(DT*ii-0.5), DT*ii-0.5))
        }

        art15plus.inits <- pmax(artnum.ii - colSums(hivpop[-1,,h.age15plus.idx,,i],,3), 0)

        art15plus.elig <- sweep(hivpop[1,,h.age15plus.idx,,i], 1, fp$artcd4elig[,i], "*")


        ## calculate pregnant women
        births.dist <- sweep(fp$frr.cd4 * hivpop[1,,h.fert.idx,f.idx,i], 2,
                             births.by.h.age / (ctapply(pop[p.fert.idx, f.idx, hivn.idx, i], ag.idx[p.fert.idx], sum) + colSums(fp$frr.cd4 * hivpop[1,,h.fert.idx,f.idx,i]) + colSums(fp$frr.art * hivpop[-1,,h.fert.idx,f.idx,i],,2)), "*")
        
        art15plus.elig[,h.fert.idx-min(h.age15plus.idx)+1,f.idx] <- art15plus.elig[,h.fert.idx-min(h.age15plus.idx)+1,f.idx] + DT*sweep(births.dist, 1, fp$pw.artelig[,i], "*") # multiply by DT to account for proportion of annual births occurring during this time step

        expect.mort.weight <- sweep(fp$cd4.mort[, h.age15plus.idx,], 3,
                                    colSums(art15plus.elig * fp$cd4.mort[, h.age15plus.idx,],,2), "/")
        artinit.weight <- sweep(expect.mort.weight, 3, 1/colSums(art15plus.elig,,2), "+")/2
        artinit <- pmin(sweep(artinit.weight * art15plus.elig, 3, art15plus.inits, "*"),
                        art15plus.elig)
        
        hivpop[1,, h.age15plus.idx,, i] <- hivpop[1,, h.age15plus.idx,, i] - artinit
        hivpop[2,, h.age15plus.idx,, i] <- hivpop[2,, h.age15plus.idx,, i] + artinit
      }
    }

    ## ## Remove HIV deaths from single year population projection
    ## ## Assume deaths come proportional to HIV age distribution
    ## calc.agdist <- function(x) {d <- x/rep(ctapply(x, ag.idx, sum), h.ag.span); d[is.na(d)] <- 0; d}
    ## pop[,,2,i] <- pop[,,2,i] - apply(hivdeaths, 2, rep, h.ag.span) * apply(pop[,,hivp.idx,i], 2, calc.agdist)


    ## ## incidence
    ## prev.i <- sum(pop[p.age15to49.idx,,2,i]) / sum(pop[p.age15to49.idx,,,i]) # prevalence age 15 to 49
    ## incrate15to49.i <- (fp$prev15to49[i] - prev.i)/(1-prev.i)

    ## sexinc15to49 <- incrate15to49.i*c(1, fp$inc.sexratio[i])*sum(pop[p.age15to49.idx,,hivn.idx,i])/(sum(pop[p.age15to49.idx,m.idx,hivn.idx,i]) + fp$inc.sexratio[i]*sum(pop[p.age15to49.idx, f.idx,hivn.idx,i]))

    ## agesex.inc <- sweep(fp$inc.agerr[,,i], 2, sexinc15to49/(colSums(pop[p.age15to49.idx,,hivn.idx,i] * fp$inc.agerr[p.age15to49.idx,,i])/colSums(pop[p.age15to49.idx,,hivn.idx,i])), "*")
    ## infections <- agesex.inc * pop[,,hivn.idx,i]

    ## pop[,,hivn.idx,i] <- pop[,,hivn.idx,i] - infections
    ## pop[,,hivp.idx,i] <- pop[,,hivp.idx,i] + infections

    ## hivpop[1,,,,i] <- hivpop[1,,,,i] + sweep(fp$cd4.initdist, 2:3, apply(infections, 2, ctapply, ag.idx, sum), "*")
    
    ## prevalence among pregnant women
    hivn.byage <- ctapply(rowMeans(pop[p.fert.idx, f.idx, hivn.idx,i-1:0]), ag.idx[p.fert.idx], sum)
    hivp.byage <- rowMeans(hivpop[,,h.fert.idx, f.idx,i-1:0],,3)
    pregprev <- sum(births.by.h.age * (1 - hivn.byage / (hivn.byage + colSums(fp$frr.cd4 * hivp.byage[1,,]) + colSums(fp$frr.art * hivp.byage[-1,,],,2)))) / sum(births.by.age)
    if(i+AGE_START <= PROJ_YEARS)
      pregprevlag[i+AGE_START-1] <- pregprev


    ## output
    ## incrate15to49[i] <- incrate15to49.i
    ## sexinc15to49out[,i] <- sexinc15to49
    }

  attr(pop, "incrate15to49") <- incrate15to49
  attr(pop, "sexinc") <- sexinc15to49out
  attr(pop, "hivpop") <- hivpop
  attr(pop, "pregprevlag") <- pregprevlag
  attr(pop, "paedsurvout") <- paedsurvout
  attr(pop, "incrate15to49.ts") <- incrate15to49.ts.out
  return(pop)
}
