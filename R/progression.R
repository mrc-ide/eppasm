# Natural deaths
# -----------------------------------------------------------------------------
epp_death <- function(MODEL, year, fp, pop, hivpop, hiv_db, artpop, art_db) {  
  deaths <- sweep(pop$data[,,,year], 1:2, 1 - fp$Sx[,,year], "*")

  if (MODEL!=0) {
    hiv.sx.prob <- 1 - sumByAGs(deaths[,,2], fp) / 
                       sumByAGs(pop$data[,,2,year], fp)
    hiv.sx.prob[is.nan(hiv.sx.prob)] <- 0
  }
  
  pop$set("-", deaths, YEAR = year)
  pop$natdeaths[,,year] <- rowSums(deaths,,2)
  
  if (MODEL!=0) {
    hivpop$sweep_sex("*", hiv.sx.prob, year)
    if (MODEL==2)
      hiv_db$sweep_sex("*", hiv.sx.prob, year)
    if (year > fp$tARTstart) {
      artpop$sweep_sex("*", hiv.sx.prob, year)
      if (MODEL==2)
        art_db$sweep_sex("*", hiv.sx.prob, year)
    }
  }
}

# Migration at year i
# -----------------------------------------------------------------------------
epp_migration <- function(MODEL, year, fp, pop, hivpop, hiv_db, artpop, art_db) {
  netmigsurv <- fp$netmigr[,,year] * (1 + fp$Sx[,,year]) / 2
  mr.prob <- 1 + netmigsurv / rowSums(pop$data[,,,year],,2)

  if (MODEL!=0) {
    hiv.mr.prob <- sumByAGs(mr.prob * pop$data[,,2,year], fp) / 
                   sumByAGs(pop$data[,,2,year], fp)
    hiv.mr.prob[is.nan(hiv.mr.prob)] <- 0
  }
  
  pop$sweep_sex("*", mr.prob, year)

  if (MODEL!=0) {
    hivpop$sweep_sex("*", hiv.mr.prob, year)
    if (MODEL==2) 
      hiv_db$sweep_sex("*", hiv.mr.prob, year)
    if (year > fp$tARTstart) {
      artpop$sweep_sex("*", hiv.mr.prob, year)
      if (MODEL==2)
        art_db$sweep_sex("*", hiv.mr.prob, year)
    }
  }
}

## adjust population to match target population size
# -----------------------------------------------------------------------------
epp_adjust_pop <- function(MODEL, fp, i, pop, hivpop, hiv_db, artpop, art_db) {
  pop$popadjust[,,i] <- fp$targetpop[,,i] / rowSums(pop$data[,,,i],,2)
  if (MODEL!=0) {
    hiv_adj_prob  <- sumByAGs(pop$popadjust[,,i] * pop$get(YEAR=i, DS=2), fp) / 
                     sumByAGs(pop$get(YEAR=i, DS=2), fp)
    hiv_adj_prob[is.nan(hiv_adj_prob)] <- 0
  }
  pop$sweep_sex("*", pop$popadjust[,,i], i)
  if (MODEL!=0) {
    hivpop$sweep_sex("*", hiv_adj_prob, i)
    if (MODEL==2) 
      hiv_db$sweep_sex("*", hiv_adj_prob, i)
    if (i >= fp$tARTstart) {
      artpop$sweep_sex("*", hiv_adj_prob, i)
      if (MODEL==2) 
        art_db$sweep_sex("*", hiv_adj_prob, i)
    }
  }
}

# calculate, distribute eligible for ART, update grad, gradART
# -----------------------------------------------------------------------------
epp_art_init <- function(MODEL, birth_agrp, pop, hivpop, artpop,
                     hiv_db, art_db, grad, grad_db, gradART, gradART_db, 
                     fp, i, ii) {
  
  list2env(fp$ss, environment())

  # progression and mortality
  gradART$data[,,,] <- 0 # reset every time step the gradient
  gradART$progress(artpop$get(YEAR=i), i)
  ## ART dropout
  grad$set("+", fp$art_dropout[i] * colSums(artpop$get(YEAR=i)))
  gradART$set("-", fp$art_dropout[i] * artpop$get(YEAR=i))

  if (MODEL==2) {
    gradART_db$data[,,,] <- 0 # reset every time step the gradient
    gradART_db$progress(art_db$get(YEAR=i), i)
    grad_db$set("+", fp$art_dropout[i] * colSums(art_db$get(YEAR=i)))
    gradART_db$set("-", fp$art_dropout[i] * art_db$get(YEAR=i))
  }

  ## calculate number eligible for ART
  artcd4_percelig <- f_artcd4_percelig(fp, i)
  art_elig <- sweep(hivpop$data[,h.age15plus.idx,,i], 1, artcd4_percelig, "*")
  
  ## calculate pregnant women   
  if (fp$pw_artelig[i] & fp$artcd4elig_idx[i] > 1) {
    sus_pop <- pop$data
    if (MODEL==2)
      sus_pop[db_aid,,,] <- sus_pop[db_aid,,,] - pop$pop_db
    art_elig <- updatePreg(art_elig, birth_agrp, fp, i, sus_pop, hivpop, artpop)
  }

  # add sexual inactive but eligible for treatment
  if (MODEL==2) {
    art_elig_db <- sweep(hiv_db$data[,h.age15plus.idx,,i], 1, artcd4_percelig, "*")
    art_elig <- art_elig + art_elig_db 
  }

  ## calculate number to initiate ART based on number or percentage
  art_curr <- colSums(artpop$data[,,h.age15plus.idx,,i],,3) + 
              colSums(gradART$data[,,h.age15plus.idx,],,3) * DT
  
  if (MODEL==2) {
    art_curr_db <- colSums(art_db$data[,,h.age15plus.idx,,i],,3) + 
                     colSums(gradART_db$data[,,h.age15plus.idx,],,3) * DT
    art_curr   <- art_curr + art_curr_db
  }

  # number on ART this ts and initiation distribution
  artnum.ii       <- f_artInit(art_curr, art_elig, fp, i, ii)
  art15plus.inits <- pmax(artnum.ii - art_curr, 0)
  artinit         <- f_artDist(art_elig, art15plus.inits, fp, i)

  if (MODEL==1) 
    artinit <- pmin(artinit, hivpop$get(YEAR=i) + DT * grad$data)

  if (MODEL==2) { # split the number proportionally for active and idle pop
    all_hivpop     <- (hivpop$get(YEAR=i) + DT * grad$data) +  
                      (hiv_db$get(YEAR=i) + DT * grad_db$data)
    artinit        <- pmin(artinit, all_hivpop)
    pr.weight_db   <- (hiv_db$get(YEAR=i) + DT * grad_db$data) / all_hivpop
    artinit_db     <- artinit * pr.weight_db
    artinit        <- artinit - artinit_db
    
    grad_db$set("-", artinit_db / DT, AG = h.age15plus.idx)
    gradART_db$set("+", artinit_db / DT, TS = 1, AG = h.age15plus.idx)
  }

  grad$set("-", artinit / DT, AG = h.age15plus.idx)
  gradART$set("+", artinit / DT, TS = 1, AG = h.age15plus.idx)

  artpop$data[,,,,i] <- artpop$data[,,,,i] +  DT * gradART$data
  if (MODEL==2)
    art_db$data[,,,,i] <- art_db$data[,,,,i] +  DT * gradART_db$data
  invisible()
}

# EPP populations aging
# -----------------------------------------------------------------------------
epp_aging <- function(MODEL, year, fp, pop, hivpop, artpop, hiv_db, art_db) {
  list2env(fp$ss, environment())
  
  pop$aging(year) # age the population
  
  ## Add lagged births into youngest age group
  if (exists("popadjust", where=fp) & fp$popadjust) {
    if (MODEL==0) 
      hivn_entrants <- fp$entrantpop[,year-1]
    else {
      hivn_entrants <- fp$entrantpop[,year-1]*(1-fp$entrantprev[,year])
      hivp_entrants <- fp$entrantpop[,year-1]*fp$entrantprev[,year]
    }
  } else {
    if (MODEL==0)
      hivn_entrants <- HnIn0(fp, i, pop$birthslag)
    else {
      hivn_entrants <- HnIn(fp, year, pop$pregprevlag, pop$birthslag)
      hivp_entrants <- HpIn(fp, year, pop$pregprevlag, pop$birthslag)
    }
  }

  if (MODEL!=0) {
    pop$entrantprev[year] <- sum(hivp_entrants) / sum(hivn_entrants+hivp_entrants)
    pop$hivp_entrants_out[,year] <- sum(hivp_entrants)
    pop$put(hivp_entrants, AG=1, DS=2, YEAR=year)
  }

  pop$put(hivn_entrants, AG=1, DS=1, YEAR=year)

  if (MODEL!=0) {
    
    if (MODEL==2) {
      pop$pop_db[-1,,,year] <- pop$pop_db[-fp$ss$pDB,,,year-1]
      pop$pop_db[1,,1,year] <- hivn_entrants
      pop$pop_db[1,,2,year] <- hivp_entrants
      pop$pop_db[,,,year]   <- sweep(pop$pop_db[,,,year], 1:2, 1 - fp$db_pr, "*")
    }

    hiv.ag.prob <- pop$data[aglast.idx,,2,year-1] / sumByAGs(pop$data[,,2,year-1], fp)
    hiv.ag.prob[is.nan(hiv.ag.prob)] <- 0

    noart <- hivp_entrants * (1-fp$entrantartcov[,year])
    isart <- hivp_entrants *    fp$entrantartcov[,year]

    if (MODEL==2) {
      noart <- hivp_entrants * (1-fp$entrantartcov[, year]) * fp$db_pr[1, ]
      noart_db <- hivp_entrants * (1-fp$entrantartcov[, year]) * (1 - fp$db_pr[1, ])
      isart <- hivp_entrants * fp$entrantartcov[, year]  * fp$db_pr[1, ]
      isart_db <- hivp_entrants * fp$entrantartcov[, year] * (1 - fp$db_pr[1, ])
    }
  
    hivpop$aging(hiv.ag.prob, noart, year)
  
    if (MODEL==2) {
      hiv_db$aging(hiv.ag.prob, noart_db, year)
      debut_hiv <- sweep(hiv_db$data[, db_agr[-1],, year], 2:3, fp$db_pr[-1, ], '*')
      hivpop$set("+", debut_hiv, AG = db_agr[-1], YEAR = year)
      hiv_db$set("-", debut_hiv, AG = db_agr[-1], YEAR = year)
    }

    if (year > fp$tARTstart) {
      artpop$aging(hiv.ag.prob, isart, year)
      if (MODEL==2) {
        art_db$aging(hiv.ag.prob, isart_db, year)
        debut_art <- sweep(art_db$data[,, db_agr[-1],,year], 3:4, fp$db_pr[-1,], '*')
        artpop$set("+", debut_art, AG = db_agr[-1], YEAR = year)
        art_db$set("-", debut_art, AG = db_agr[-1], YEAR = year)
      }
    }
  }
}

# Disease model
# -----------------------------------------------------------------------------
epp_disease_model <- function(MODEL, year, fp, pop, hivpop, artpop,
                              hiv_db, art_db, grad, gradART, grad_db, gradART_db, rvec) {
  DT <- fp$ss$DT
  for (time_step in seq_len(fp$ss$hiv_steps_per_year)) {
    ts <- (year-2)/DT + time_step
    grad$data[,,] <- 0
    if (MODEL==2)
      grad_db$data[,,] <- 0
    if (fp$eppmod != "directincid") {
      if (fp$eppmod %in% c("rtrend", "rtrend_rw"))
        rvec[ts] <- calc_rtrend_rt(fp$proj.steps[ts], fp, rvec[ts-1], prevlast, pop, year, time_step)
      else
        rvec[ts] <- fp$rvec[ts]

      ## number of infections by age / sex
      infect <- infect_spec(MODEL, fp, pop, hivpop, artpop, year, time_step, rvec[ts])
      pop$incrate15to49_ts[ts] <- attr(infect, "incrate15to49.ts")
      pop$prev15to49_ts[ts] <- prevlast <- attr(infect, "prevcurr")
      change_ts <- DT * infect
      pop$set("-", change_ts, DS = 1, YEAR = year)
      pop$set("+", change_ts, DS = 2, YEAR = year)
      pop$infections[,,year] <- pop$infections[,,year]   + change_ts 
      grad$set("+", sweep(fp$cd4_initdist, 2:3, sumByAGs(infect, fp), "*"))
      pop$incid15to49[year] <- pop$incid15to49[year] + sum(change_ts[fp$ss$p.age15to49.idx,])
    }

    ## cd4 disease progression and mortality
    cd4_mort <- scale_cd4_mort(fp, hivpop, artpop)
    grad$progress(hivpop$get(YEAR=year), cd4_mort)
    if (MODEL==2) {
      cd4_mort_db <- scale_cd4_mort(fp, hiv_db, art_db)
      grad_db$progress(hiv_db$get(YEAR=year), cd4_mort_db)
    }
    ## Remove hivdeaths from pop
    hivdeaths_ts <- f_hiv_death(MODEL, cd4_mort, cd4_mort_db, pop, year, fp,
                                hivpop, artpop, hiv_db, art_db)
    pop$data[,,2,year] <- pop$data[,,2,year] - hivdeaths_ts
    pop$hivdeaths[,,year]  <- pop$hivdeaths[,,year] + hivdeaths_ts
    # ART initiation
    if (year >= fp$tARTstart)
      epp_art_init(MODEL, pop$birth_agrp, pop, hivpop, artpop,
                   hiv_db, art_db, grad, grad_db, gradART, gradART_db, 
                   fp, year, time_step)
    
    hivpop$set("+", DT * grad$data, YEAR=year)
    if (MODEL==2)
      hiv_db$set("+", DT * grad_db$data, YEAR=year)
  } # end ii
}


# EPP disease model direct incidence
# -----------------------------------------------------------------------------
epp_disease_model_direct <- function(pop, hivpop, year, fp) {
  # agesex.inc <- f_infections_directincid(pop, year, fp)
  list2env(fp$ss, environment())
  if (fp$incidpopage == 0L) # incidence for 15 -49 population
    age_id <- p.age15to49.idx
  else if (fp$incidpopage == 1L) # incidence for 15+ population
    age_id <- p.age15plus.idx
  num <- c(1, fp$incrr_sex[year]) * sum(pop$data[age_id,, hivn.idx, year-1])
  den <- sum(pop$data[age_id, m.idx, hivn.idx, year-1]) + 
         sum(pop$data[age_id, f.idx, hivn.idx, year-1]) * fp$incrr_sex[year]
  sexinc <- fp$incidinput[year] * num /den
  ageinc <- colSums(pop$data[age_id,,hivn.idx,year-1] * fp$incrr_age[age_id,,year]) /
            colSums(pop$data[age_id,,hivn.idx,year-1])

  agesex.inc <- sweep(fp$incrr_age[,,year], 2, sexinc / ageinc, "*")
  pop$infections[,,year] <- agesex.inc * pop$data[,,1,year-1]
  pop$data[,,1,year] <- pop$data[,,1,year] - pop$infections[,,year]
  pop$data[,,2,year] <- pop$data[,,2,year] + pop$infections[,,year]
  infect_agrp <- sumByAGs(pop$infections[,,year], fp)
  hivpop$set("+", sweep(fp$cd4_initdist, 2:3, infect_agrp, "*"), YEAR=year)
  pop$incid15to49[year] <- sum(pop$infections[fp$ss$p.age15to49.idx,,year])
}