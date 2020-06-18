# Annualized number of new infections
# infect_spec <- function(pop, hivpop, artpop, ii, fp)
# K moved to popClass method
infectFns <- c(
infect_mix = function(hivpop, artpop, ii) {
    ts <- (year-2)/DT + ii
    update_active_pop_to(year)
    # data_active <<- sweepx(data_active, 1:2, p$est_senesence)

    actual_active <- data_active # at this point we have the "real" number of
                                 # people to calculate the prevalence; after adjustment
                                 # and balancing, the prev would not be correct

    # Partner change rate as a risk reduction
    # data_active <<- sweepx(data_active, 1:2, 1+p$est_pcr)

    # balancing
    prop_n_m <- data_active[, m.idx, hivn.idx]/ rowSums(data_active[, m.idx, ])
    prop_n_f <- data_active[, f.idx, hivn.idx]/ rowSums(data_active[, f.idx, ])

    # sweep over sexual mixing matrices
    nc_m <- sweepx(p$mixmat[,,m.idx], 1, rowSums(data_active[, m.idx, ]))
    nc_f <- sweepx(p$mixmat[,,f.idx], 1, rowSums(data_active[, f.idx, ]))
    
    ratio_mf <- nc_m / t(nc_f)

    nc_m_adj <- nc_m * (ratio_mf - p$balancing * (ratio_mf - 1)) / ratio_mf
    nc_f_adj <- t(nc_f) * (ratio_mf - p$balancing * (ratio_mf - 1))
    
    n_m_active_negative <- sweepx(nc_m_adj, 1, prop_n_m)
    n_f_active_negative <- sweepx(t(nc_f_adj), 1, prop_n_f)
    
    art_cov <- matrix(0, pAG, NG)
    if (year >= p$tARTstart) {
      art_ <- colSums(artpop$data[,,,,year] + artpop$data_db[,,,,year],,2)
      hiv_ <- colSums(hivpop$data[,,,year] + hivpop$data_db[,,,year],,1)
      art_cov <- art_/(art_+hiv_)
      art_cov <- sapply(1:2, function(x) rep(art_cov[, x], h.ag.span))
    }

    hiv_treated       <- data_active[,,hivp.idx] * art_cov
    hiv_not_treated   <- data_active[,,hivp.idx] - hiv_treated
    transm_prev <- (sum(hiv_not_treated) + sum(hiv_treated) * (1 - p$relinfectART)) / 
                    sum(rowSums(actual_active,,2)) # prevalence adjusted for art
    # +intervention effects and time epidemic start
    w  <- p$iota * (p$proj.steps[ts] == p$tsEpidemicStart)
    transm_prev <- rvec[ts] * transm_prev + w

    inc_m <- sweepx(n_m_active_negative, 2, transm_prev)
    inc_m <- sweepx(inc_m, 1, p$incrr_age[, m.idx, year])
    inc_f <- sweepx(n_f_active_negative, 2, transm_prev)
    inc_f <- sweepx(inc_f, 1, p$incrr_age[, f.idx, year])
    # adjusted sex
    adj_sex <- p$incrr_sex[year] * sum(inc_m)/sum(inc_f)
    if (is.na(adj_sex)) adj_sex <- 1
    inc_f <- inc_f * adj_sex

    if (ii==10) {
      WAIFW[,,m.idx,year] <<- inc_m
      WAIFW[,,f.idx,year] <<- inc_f
    }

    infections.ts <- cbind(rowSums(inc_m), rowSums(inc_f))

    # incrate15to49_ts[,,ts] <<- transm_prev
    prev15to49_ts[ts] <<- prevlast <<- sum(data[,,hivp.idx,year])/sum(data[,,,year])
    infections.ts
},

infect_spec = function(hivpop, artpop, time_step) {
    ts    <- (year-2) / DT + time_step
    dt_ii <- 1 - DT * (time_step - 1) # transition of population    
    update_active_pop_to(year)
  
    # counting all negative including virgin
    hivn_both       <- n_NEG(dt_ii)
    hivp_inactive   <- 0
    if (MODEL==2) # add safe positive 
      hivp_inactive <- VIRGIN$n_HIV(dt_ii)
    hivp_active     <- n_HIV(dt_ii) - hivp_inactive
    all_pop         <- N(dt_ii)
  
    art.ii  <- 0
    if (year >= p$tARTstart)
      art.ii  <- artpop_adj(hivpop, artpop, dt_ii)

    transm_prev <- (hivp_active - art.ii*(1 - p$relinfectART)) / all_pop
    w <- p$iota * (p$proj.steps[ts] == p$tsEpidemicStart)
    inc_rate <- rvec[ts] * transm_prev + w

    sus_age_sex <- data_active[p.age15to49.idx,,hivn.idx]
    adj_sex <- sum(sus_age_sex) /
      ( sum(data_active[p.age15to49.idx,m.idx,hivn.idx]) +
        sum(data_active[p.age15to49.idx,f.idx,hivn.idx]) * 
        p$incrr_sex[year] )
    sexinc15to49.ts <- inc_rate * c(1, p$incrr_sex[year]) * adj_sex
    # New infections distributed by age: ratio age_i/ 25-29 age
    adj_age <- sexinc15to49.ts * colSums(sus_age_sex) /
      colSums(sus_age_sex * p$incrr_age[p.age15to49.idx,,year])
    agesex.inc <- sweep(p$incrr_age[,,year], 2, adj_age, "*")

    ## Adjust age-specific incidence among men for circumcision coverage
    agesex.inc[, m.idx] <- agesex.inc[, m.idx] * (1 - p$circ_incid_rr * p$circ_prop[,year])
    infections.ts <- agesex.inc * data_active[,,hivn.idx]
    # saving
    incrate15to49_ts[ts] <<- inc_rate
    prev15to49_ts[ts] <<- prevlast <<- (hivp_active + hivp_inactive) / all_pop
    infections.ts
},

epp_disease_model_direct = function(hivpop, artpop) {
  if (p$incidpopage) 
    age_id <- p.age15plus.idx # incidence for 15+ population
  else
    age_id <- p.age15to49.idx # incidence for 15 -49 population
  update_active_pop_to(year-1)
  num <- c(1, p$incrr_sex[year]) * sum(data_active[age_id,, hivn.idx])
  den <- sum(data_active[age_id, m.idx, hivn.idx]) + 
         sum(data_active[age_id, f.idx, hivn.idx]) * 
         p$incrr_sex[year]
  sexinc <- p$incidinput[year] * num / den
  ageinc <- colSums( data_active[age_id,,hivn.idx] * 
                     p$incrr_age[age_id,,year] ) /
            colSums(data_active[age_id,,hivn.idx])
  agesex.inc <- sweep(p$incrr_age[,,year], 2, sexinc / ageinc, "*")
  infections[,,year]     <<- agesex.inc * data[,,hivn.idx,year-1]
  data[,,hivn.idx,year]  <<- data[,,hivn.idx,year] - infections[,,year]
  data[,,hivp.idx,year]  <<- data[,,hivp.idx,year] + infections[,,year]
  infect_agrp            <- sumByAGs(infections[,,year], ag.idx)
  hivpop$data[,,,year] <- hivpop$data[,,,year] + 
                          sweep(p$cd4_initdist, 2:3, infect_agrp, "*")
  incid15to49[year]      <<- sum(infections[p.age15to49.idx,,year])
})

setMembers(popEPP, "public", names(infectFns), infectFns)
infectFns <- NULL # optional

calc_infections_simpletransm <- function(fp, pop, hivpop, artpop, i, ii, r_ts){

  ## Attach state space variables
  invisible(list2env(fp$ss, environment())) # put ss variables in environment for convenience

  ts <- (i-2)/DT + ii
  
  ## Calculate prevalence of unsuppressed viral load among sexually active population

  contacts_ii <- sweep((pop[c(2:pAG, pAG), , , i] * (1-DT*(ii-1)) + pop[ , , , i] * DT*(ii-1)),
                       1:2, fp$relbehav_age, "*")

  
  hivpop_w_ha <- colSums(sweep(hivpop[ , , , i], 1, fp$relsexact_cd4cat, "*"), , 1)
  hivpop_ha <- colSums(hivpop[ , , , i],,1)
  artpop_ha <- colSums(artpop[ , , , , i],,2)

  hivcontacts_ha <- (hivpop_w_ha + artpop_ha) / (hivpop_ha + artpop_ha)
  hivcontacts_ha[is.na(hivcontacts_ha)] <- 0

  hivtransm_ha <- (hivpop_w_ha + fp$relinfectART * artpop_ha) / (hivpop_ha + artpop_ha)
  hivtransm_ha[is.na(hivtransm_ha)] <- 0

  hivn_ii <- contacts_ii[ , , hivn.idx]
  hivcontacts_ii <- contacts_ii[ , , hivp.idx] * hivcontacts_ha[fp$ss$ag.idx, ]
  hivtransm_ii <- contacts_ii[ , , hivp.idx] * hivtransm_ha[fp$ss$ag.idx, ]

  hivtransm_prev <- colSums(hivtransm_ii) / (colSums(hivn_ii) + colSums(hivcontacts_ii))

  ## r_sex[1:2] is the transmission rate by (Men, Women)
  r_sex <- c(sqrt(fp$mf_transm_rr[i]), 1/sqrt(fp$mf_transm_rr[i])) * r_ts

  sexinc15to49.ts <- (r_sex * hivtransm_prev)[2:1] + fp$mf_transm_rr[i]^c(-0.25, 0.25) * fp$iota * (fp$proj.steps[ts] == fp$tsEpidemicStart)
  agesex.inc <- sweep(fp$incrr_age[,,i], 2, sexinc15to49.ts/(colSums(pop[p.age15to49.idx,,hivn.idx,i] * fp$incrr_age[p.age15to49.idx,,i])/colSums(pop[p.age15to49.idx,,hivn.idx,i])), "*")

  ## Adjust age-specific incidence among men for circumcision coverage
  agesex.inc[ , m.idx] <- agesex.inc[ , m.idx] * (1 - fp$circ_incid_rr * fp$circ_prop[ , i])
  
  infections.ts <- agesex.inc * pop[,,hivn.idx,i]

  attr(infections.ts, "incrate15to49.ts") <- 0 # sum(infections.ts[p.age15to49.idx,]) / sum(hivn.ii)
  ## attr(infections.ts, "prevcurr") <- sum(hivp_noart.ii+art.ii) / sum(hivn.ii+hivp_noart.ii+art.ii)

  attr(infections.ts, "prevcurr") <- 0 # sum(hivp.ii) / sum(hivp.ii + hivn.ii)

  return(infections.ts)
}