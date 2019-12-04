setMembers(popEPP, "public", "initialize",
function(fp, MODEL=1, VERSION="R", MIX=F) {
    super$initialize(fp) 
    # 'Initialize pop array'
    MODEL       <<- MODEL
    VERSION     <<- VERSION
    MIX         <<- MIX
    data        <<- array(0, c(pAG, NG, pDS, PROJ_YEARS))
    data_active <<- array(0, c(pAG, NG, pDS))
    data[,,hivn.idx,1] <<- p$basepop
    
    # Outputs
    entrantprev   <<- numeric(PROJ_YEARS)
    prev15to49    <<- incid15to49  <<- pregprevlag <<- pregprev <<- entrantprev
    adj_prob_age  <<- array(0, c(pAG, NG, PROJ_YEARS))
    infections    <<- hivdeaths <<- natdeaths <<- adj_prob_age
    prev15to49_ts <<- incrate15to49_ts <<- rep(NA, length(p$rvec))
    hivp_entrants_out <<- array(0, c(NG, PROJ_YEARS))
    
    # use in model
    birthslag     <<- p$birthslag
    if (p$eppmod != "directincid")
      rvec <<- if (p$eppmod=="rtrend") rep(NA, length(p$proj.steps)) else p$rvec

    if (MODEL==2)
        VIRGIN <<- virginEPP$new(fp, MODEL) # has its own initialization
    if (MIX)
      incrate15to49_ts  <<- array(0, c(pAG, NG, length(p$rvec)))
})

# named list functions as methods
popFunc <- c(
set_data = function(FUN="+", x, AG=T, NG=T, DS=T, YEAR=T) {
    FUN <- match.fun(FUN)
    data[AG, NG, DS, YEAR] <<- FUN(data[AG, NG, DS, YEAR], x)
}, 

put = function(x, AG=T, NG=T, DS=T, YEAR=T) {
    'This can be done directly with `<-` using pop_object$data'
    data[AG, NG, DS, YEAR] <<- x
}, 

get = function(YEAR=NULL, AG=T, NG=T, DS=T) {
    'This can be done directly with `[` using pop_object$data'
    if (is.null(YEAR))
        data[AG, NG, DS, TRUE]
    else 
        data[AG, NG, DS, YEAR]
},

update_active_pop_to = function(when) {
    data_active <<- data[,,,when]
    if (MODEL==2)
        data_active[db_aid,,] <<- data_active[db_aid,,] - VIRGIN$data[,,,when]
},

get_active_pop_in = function(when) {
    out <- data[,,,when]
    if (MODEL==2)
        out[db_aid,,] <- out[db_aid,,] - VIRGIN$data[,,,when]
    return(out)
},
current_active_pop = function() {
    out <- data[,,,year]
    if (MODEL==2)
        out[db_aid,,] <- out[db_aid,,] - VIRGIN$data[,,,year]
    return(out)
},

update_rvec = function(time_step) {
    ts <- (year-2) / DT + time_step
    if (p$eppmod %in% c("rtrend", "rtrend_rw"))
      rvec[ts] <<- calc_rtrend_rt(ts, time_step)
    else
      rvec[ts] <<- p$rvec[ts]
},

update_infection = function(infect) {
    infect <- infect * DT
    data[,, hivn.idx, year] <<- data[,, hivn.idx, year] - infect
    data[,, hivp.idx, year] <<- data[,, hivp.idx, year] + infect
    infections[,,year] <<- infections[,,year] + infect 
    incid15to49[year]  <<- incid15to49[year] + sum(infect[p.age15to49.idx,])  
},

remove_hiv_death = function(hivpop, artpop) {
    # death by age group
    artD <- artpop$f_death
    hivD <- hivpop$f_death
    if (MODEL==2) { # add deaths from inactive population
      artD <- artD + artpop$f_death_db
      hivD <- hivD + hivpop$f_death_db
    }
    dbyAG <- DT * (colSums(hivD) + colSums(artD,,2))
    # deaths by single-year
    hiv_mx <- dbyAG / sumByAGs(data[,, hivp.idx, year], ag.idx) # virgin included
    hiv_mx[is.nan(hiv_mx)] <- 0
    p_hiv_death_age_ <<- apply(hiv_mx, 2, rep, h.ag.span)
    hivdeaths[,,year] <<- hivdeaths[,,year] + data[,, hivp.idx, year] * p_hiv_death_age_
    data[,, hivp.idx, year] <<- data[,, hivp.idx, year] * (1 - p_hiv_death_age_)
    if (MODEL == 2)
      VIRGIN$remove_hiv_death(p_hiv_death_age_)
},

update_preg = function(art_elig, hivpop, artpop) {
  hivp <- hivpop$data[, h.fert.idx, f.idx, year] * p$frr_cd4[,, year]
  update_active_pop_to(year)
  hivn <- sumByAG(data_active[p.fert.idx, f.idx, hivn.idx], ag.idx, T, p.fert.idx)
  art <- colSums(artpop$data[,,h.fert.idx, f.idx, year] * p$frr_art[,,,year],,2)
  birthdist <- sweep(hivp, 2, birth_agrp / (hivn + colSums(hivp) + art), "*")
  elDS <- 1:(p$artcd4elig_idx[year] - 1)
  elAG <- h.fert.idx - min(h.age15plus.idx) + 1
  art_elig[elDS, elAG, f.idx] <- art_elig[elDS, elAG, f.idx] + birthdist[elDS,]
  art_elig
},

art_initiate = function(art_curr, art_elig, time_step) {
  out    <- c(0,0)
  year_w <- if(DT * time_step < 0.5) 0 else 1
  trans  <- DT * time_step + 0.5 - year_w
  years  <- year - ((2:1) - year_w)
  for(sex in 1:2) {
    if( !any(p$art15plus_isperc[sex, years]) ) { # both number
      out[sex] <- sum(p$art15plus_num[sex, years] %*% c(1 - trans, trans))
      if (MIX)
        artcov[sex] <<- min(1, out[sex] / (sum(art_elig[,,sex]) + art_curr[sex]))
    }
    else if (all(p$art15plus_isperc[sex, years])) { # both percentage
      cov <- sum(p$art15plus_num[sex, years] %*% c(1 - trans, trans))
      if (MIX)
        artcov[sex] <<- min(1, cov) # save for infect_mix
      out[sex] <- cov * (sum(art_elig[,,sex]) + art_curr[sex])
    }
    else if (!p$art15plus_isperc[sex, years[1]] &
              p$art15plus_isperc[sex, years[2]]) { # transition number > % 
      actual_cov <- art_curr[sex] / (sum(art_elig[,,sex]) + art_curr[sex])
      diff_cov <- p$art15plus_num[sex, years[2]] - actual_cov
      cov <- actual_cov + diff_cov * DT / (0.5 + year_w - DT*(time_step - 1))
      if (MIX)
        artcov[sex] <<- min(1, cov) # save for infect_mix
      out[sex] <- cov * (sum(art_elig[,,sex]) + art_curr[sex])
    }
  }
  out
}, 

## calculate ART initiation distribution
art_distribute = function(art_elig, art_need) {
  if (!p$med_cd4init_input[year]) {
    if (p$art_alloc_method == 4L) { ## by lowest CD4
      ## Calculate proportion to be initiated in each CD4 category
      out <- array(0, dim(art_elig))
      for(m in hDS:1){
        elig_hm <- colSums(art_elig[m,,])
        init_pr <- if (elig_hm == 0) elig_hm 
                      else pmin(1.0, art_need / elig_hm, na.rm = TRUE)
        out[m,,] <- sweep(art_elig[m,,], 2, init_pr, "*")
      }
    } 
    else { # Spectrum Manual p168--p169, 
      expect.mort.w <- sweep(p$cd4_mort[, h.age15plus.idx,], 3,
                    colSums(art_elig * p$cd4_mort[, h.age15plus.idx,],,2), "/")
      init.w <- sweep(p$art_alloc_mxweight * expect.mort.w, 3,
                      (1 - p$art_alloc_mxweight) / colSums(art_elig,,2), "+")
      out <- pmin(sweep(init.w * art_elig, 3, art_need, "*"), art_elig)
    }
  } 
  else {

    CD4_LO <- c(500, 350, 250, 200, 100, 50, 0)
    CD4_UP <- c(1000, 500, 350, 250, 200, 100, 50)
    
    j <- p$med_cd4init_cat[year]

    pr_below <- (p$median_cd4init[year] - CD4_LO[j]) / (CD4_UP[j] - CD4_LO[j])
    
    elig_below <- colSums(art_elig[j,,,drop=FALSE],,2) * pr_below
    if(j < hDS) 
      elig_below <- elig_below + colSums(art_elig[(j+1):hDS,,,drop=FALSE],,2)
    
    elig_above <- colSums(art_elig[j,,,drop=FALSE],,2) * (1.0-pr_below)
    if(j > 1) 
      elig_above <- elig_above + colSums(art_elig[1:(j-1),,,drop=FALSE],,2)

    initpr_below <- pmin(art_need * 0.5 / elig_below, 1.0, na.rm=TRUE)
    initpr_above <- pmin(art_need * 0.5 / elig_above, 1.0, na.rm=TRUE)
    initpr_medcat <- initpr_below * pr_below + initpr_above * (1-pr_below)

    out <- array(0, dim=c(hDS, hAG, NG))

    if(j < hDS)
      out[(j+1):hDS,,] <- sweep(art_elig[(j+1):hDS,,,drop=F], 3, initpr_below, "*")
    
    out[j,,] <- sweep(art_elig[j,,,drop=F], 3, initpr_medcat, "*")
    
    if(j > 0)
      out[1:(j-1),,] <- sweep(art_elig[1:(j-1),,,drop=F], 3, initpr_above, "*")
  }
  out
},

adjust_pop = function() {
  adj_prob_age[,,year] <<- p$targetpop[,,year] / rowSums(data[,,,year],,2)
  if (MODEL!=0) {
    adj_prob_agr <<- sumByAGs(adj_prob_age[,,year] * data[,,hivp.idx,year], ag.idx)/
                 sumByAGs(data[,,hivp.idx,year], ag.idx)
    adj_prob_agr[is.nan(adj_prob_agr)] <<- 0
  }
  data[,,,year] <<- sweep(data[,,,year], 1:2, adj_prob_age[,,year], "*")
},

aging = function() {
    data[-c(1,pAG),,,year] <<- data[-(pAG-1:0),,,year-1]
    data[pAG,,,year] <<- data[pAG,,,year-1] + data[pAG-1,,,year-1] # open ended
},

add_entrants = function() {
    ## Add lagged births into youngest age group and update population
    if (exists("popadjust", where=p) & p$popadjust) {
        if (MODEL==0)
            healthy <- p$entrantpop[,year-1]
        else {
            healthy <- p$entrantpop[,year-1]*(1-p$entrantprev[,year])
            positiv <- p$entrantpop[,year-1]*   p$entrantprev[,year]
        }
    } 
    else {
        if (MODEL==0)
            healthy <- birthslag[,year-1] * p$cumsurv[,year-1] / 
                       p$paedsurv_lag[year-1] + p$cumnetmigr[,year-1]
        else {
            survivedBirth <- birthslag[,year-1]*p$cumsurv[,year-1]
            prev_now <- p$entrantprev[,year]  # avoid repeat accesses
            imm_now  <- p$cumnetmigr[,year-1] 
            healthy <- survivedBirth * (1-prev_now/p$paedsurv_lag[year-1]) + 
                        imm_now*(1 - pregprevlag[year-1]*p$netmig_hivprob)
            positiv <- (survivedBirth + imm_now) * prev_now
        }
    }
    # save and update pop        
    data[1,, hivn.idx, year] <<- healthy
    if (MODEL != 0) {
        data[1,, hivn.idx, year] <<- healthy
        data[1,, hivp.idx, year] <<- positiv
    }
    if (MODEL==2) {  # this is unwise
        VIRGIN$data[1,,hivn.idx,year] <<- healthy
        VIRGIN$data[1,,hivp.idx,year] <<- positiv
    }
    if (MODEL!=0) {
        entrantprev[year]        <<- sum(positiv)/sum(healthy+positiv)
        hivp_entrants_out[,year] <<- positiv
    }
},

hiv_aging_prob = function() {
    update_active_pop_to(year-1) # now data_all refers to total pop: debut+not
    hiv.ag.prob <-          data_active[aglast.idx,,hivp.idx] /
                   sumByAGs(data_active[          ,,hivp.idx], ag.idx)
    hiv.ag.prob[is.nan(hiv.ag.prob)] <- 0
    hiv.ag.prob
},

entrant_art = function() { # return these for updating HIV and ART pop
    yesno <- c(p$entrantartcov[,year], 1 - p$entrantartcov[,year])
    out <- data[1,, hivp.idx, year] * yesno
    out # 1:2 ART+, 3:4 ART-
},

deaths = function() {
    death_now <- sweep(data[,,,year], 1:2, 1 - p$Sx[,,year], "*")
    if (MODEL != 0) {
      hiv_sx_prob <<- 1 - sumByAGs(death_now[,,hivp.idx], ag.idx) / 
                          sumByAGs(data[,,hivp.idx,year], ag.idx)
      hiv_sx_prob[is.nan(hiv_sx_prob)] <<- 0
    }
    data[,,,year] <<- data[,,,year] - death_now
    natdeaths[,,year] <<- rowSums(death_now,,2)
},

migration = function() {
    netmigsurv <- p$netmigr[,,year] * (1 + p$Sx[,,year]) / 2
    mr_prob_ <<- 1 + netmigsurv / rowSums(data[,,,year],,2)
    if (MODEL != 0) {
      hiv_mr_prob <<- sumByAGs(mr_prob_ * data[,,hivp.idx,year], ag.idx) / 
                      sumByAGs(           data[,,hivp.idx,year], ag.idx)
      hiv_mr_prob[is.nan(hiv_mr_prob)] <<- 0
    }
    data[,,,year] <<- sweep(data[,,,year], 1:2, mr_prob_, "*")
},

update_fertile = function() { # only on active pop
    update_active_pop_to(year)
    two_years  <- data_active + get_active_pop_in(year-1)
    birth_age  <<- rowSums(two_years[p.fert.idx, f.idx,])/2 * p$asfr[,year]
    # adjusted ASFR to match
    if (MODEL==2) {
      two_years  <- data[p.fert.idx,f.idx,,year] + data[p.fert.idx,f.idx,,year-1]
      N_mid_star <- rowSums(two_years/2)
      birth_age <<- birth_age * p$asfr[, year] / (birth_age / N_mid_star)
    }
    birth_agrp <<- sumByAG(birth_age, ag.idx, TRUE, p.fert.idx)
    births     <- p$srb[,year] * sum(birth_agrp)
    if (year + AGE_START <= PROJ_YEARS)
      birthslag[, year + AGE_START - 1] <<- births
},

cal_prev_pregant = function(hivpop, artpop) { # only on active pop
    years   <- year - 1:0
    update_active_pop_to(year)
    two_years <- data_active + get_active_pop_in(year-1)
    meanWomen <- two_years[p.fert.idx, f.idx, hivn.idx] / 2
    hivn <- sumByAG(meanWomen, ag.idx, TRUE, p.fert.idx)
    hivp <- rowMeans(hivpop$get(AG=h.fert.idx, NG=f.idx, YEAR=years),,2)
    art  <- rowMeans(artpop$get(AG=h.fert.idx, NG=f.idx, YEAR=years),,3)
    pregprev[year] <<- sum(birth_agrp * 
      (1 - hivn / (hivn + colSums(p$frr_cd4[,,year] * hivp) + 
      colSums(p$frr_art[,,,year] * art,,2)))) / sum(birth_age)
    if (year + AGE_START <= PROJ_YEARS)
      pregprevlag[year + AGE_START - 1] <<- pregprev[year]
},

save_prev_n_inc = function() {
    prev15to49[year]  <<- sum(data[p.age15to49.idx,, hivp.idx, year]) / 
                          sum(data[p.age15to49.idx,,         , year])
    update_active_pop_to(year-1)
    incid15to49[year] <<- incid15to49[year] /
                          sum(data_active[p.age15to49.idx,,hivn.idx])
    prev[year]  <<- sum(data[,,hivp.idx,year]) / sum(data[,,,year])
    incid[year] <<- incid15to49[year] / sum(data_active[,,hivn.idx]) # TODO:
}, 

calc_rtrend_rt = function(ts, time_step) {
  stop("Kinh did not fix calc_rtrend_rt")
  rveclast <- rvec[ts-1]
  dtii     <- 1-DT*(time_step-1)
  hivn.ii <- sum(data[p.age15to49.idx,,hivn.idx,year]) - 
             sum(data[p.age15to49.idx[1],,hivn.idx,year]) * dtii + 
             sum(data[tail(p.age15to49.idx,1)+1,,hivn.idx,year]) * dtii

  hivp.ii <- sum(data[p.age15to49.idx,,hivp.idx,year]) -
             sum(data[p.age15to49.idx[1],,hivp.idx,year]) * dtii +
             sum(data[tail(p.age15to49.idx,1)+1,,hivp.idx,year]) * dtii
  prevcurr <- hivp.ii / (hivn.ii + hivp.ii)
  t.ii     <- p$proj.steps[ts]
  if (t.ii > p$tsEpidemicStart){
    par <- p$rtrend
    if (t.ii < par$tStabilize)
      gamma.t <- 0
    else 
      gamma.t <- (prevcurr-prevlast) * (t.ii - par$tStabilize) / (DT * prevlast)
    logr.diff <- par$beta[2] * (par$beta[1] - rveclast) +
                 par$beta[3] * prevlast + 
                 par$beta[4] * gamma.t
    return(exp(log(rveclast) + logr.diff))
  }
  else
    return(p$rtrend$r0)
}
)

setMembers(popEPP, "public", names(popFunc), popFunc)
setMembers(popEPP, "active", 'update_year', function(x) {
  year <<- x
  if (MODEL==2)
    VIRGIN$year <- x
})