# EPP parameter super class: later improve the generating of fp as well
#' @importFrom methods new
eppFP <- R6::R6Class("eppfp", class=F, cloneable=F, portable=F, lock_objects=F,
    public = list(
        p = NULL,
        initialize = function(fp) {
            p <<- fp[-1]
            list2env(fp$ss, self)
        }
    )
)

# POP CLASS
# -----------------------------------------------------------------------------
# preferablly all fields should have a leading character
popEPP <- R6::R6Class("popepp", class=F, cloneable=F, portable=F, inherit=eppFP,
    lock_objects=F,
    public = list(
        data              = "array",
        VERSION           = "character",
        MODEL             = "numeric",
        MIX               = "logical",
        year              = "numeric",
        birth_age         = "numeric", birth_agrp        = "numeric",
        prev15to49        = "vector",  incid15to49       = "vector",
        prev              = 0,         incid             = 0,
        entrantprev       = "vector",  pregprevlag       = "vector",
        birthslag         = "array",   infections        = "array",
        hivdeaths         = "array",   natdeaths         = "array",
        popadjust         = "array",   hivp_entrants_out = "array",
        incrate15to49_ts  = "vector",  prev15to49_ts     = "vector", # can be array
        data_db           = "array",    # virgin population
        data_active       = "array",
        artcov            = numeric(2), # initially no one on treatment
        prev_last         = 0, # last time step prevalence
        hiv_sx_prob       = "array",
        hiv_mr_prob       = "array",
        adj_prob          = "array",
        rvec              = "vector")
)

setMembers(popEPP, "public", "initialize",
function(fp, MODEL=1, VERSION="R", MIX=F) {
    super$initialize(fp) 
    # 'Initialize pop array'
    MODEL       <<- MODEL
    VERSION     <<- VERSION
    MIX         <<- MIX
    data        <<- array(0, c(pAG, NG, pDS, PROJ_YEARS))
    data_active <<- array(0, c(pAG, NG, pDS))
    data[,,1,1] <<- p$basepop
    
    # Outputs
    entrantprev   <<- numeric(PROJ_YEARS)
    prev15to49    <<- incid15to49  <<- pregprevlag <<- entrantprev
    popadjust     <<- array(0, c(pAG, NG, PROJ_YEARS))
    infections    <<- hivdeaths <<- natdeaths <<- popadjust
    prev15to49_ts <<- incrate15to49_ts <<- rep(NA, length(p$rvec))
    hivp_entrants_out <<- array(0, c(NG, PROJ_YEARS))
    
    # use in model
    birthslag     <<- p$birthslag
    if (p$eppmod != "directincid")
      rvec <<- if (p$eppmod=="rtrend") rep(NA, length(p$proj.steps)) else p$rvec

    if (MODEL==2) # @debut empty pop, all inactive starts from new entrants
        data_db  <<- array(0, c(pDB, NG, pDS, PROJ_YEARS)) # debut ages only
    if (MIX)
      prev15to49_ts <<- incrate15to49_ts  <<- array(0, c(pAG, NG, length(p$rvec)))

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
        data_active[db_aid,,] <<- data_active[db_aid,,] - data_db[,,,when]
},

get_active_pop_in = function(when) {
    out <- data[,,,when]
    if (MODEL==2)
        out[db_aid,,] <- out[db_aid,,] - data_db[,,,when]
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

remove_hiv_death = function(cd4_mx, hivpop, artpop) {
    # death by age group
    artD <- artpop$f_death
    hivD <- hivpop$f_death
    if (MODEL==2) { # add deaths from inactive population
      artD <- artD + artpop$f_death_db
      hivD <- hivD + hivpop$f_death_db
    }
    dbyAG <- DT * (colSums(hivD) + colSums(artD,,2))
    # deaths by single-year
    hiv_mx <- dbyAG / sumByAGs(data[,, hivp.idx, year], ag.idx)
    hiv_mx[is.nan(hiv_mx)] <- 0
    dbyA_pr <- apply(hiv_mx, 2, rep, h.ag.span)
    hivdeaths[,,year] <<- hivdeaths[,,year] + data[,, hivp.idx, year] * dbyA_pr
    data[,, hivp.idx, year] <<- data[,, hivp.idx, year] * (1 - dbyA_pr)
    if (MODEL == 2) {
      hivdeaths[1:pDB,,year] <<- hivdeaths[1:pDB,,year] +
        data_db[,, hivp.idx, year] * dbyA_pr[1:pDB, ]
      data_db[,,hivp.idx,year] <<- data_db[,,hivp.idx,year] * (1 - dbyA_pr[1:pDB,])
    }
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
  popadjust[,,year] <<- p$targetpop[,,year] / rowSums(data[,,,year],,2)
  if (MODEL!=0) {
    adj_prob <<- sumByAGs(popadjust[,,year] * data[,,hivp.idx,year], ag.idx)/
                 sumByAGs(data[,,hivp.idx,year], ag.idx)
    adj_prob[is.nan(adj_prob)] <<- 0
  }
  data[,,,year] <<- sweep(data[,,,year], 1:2, popadjust[,,year], "*")
  if (MODEL==2)
    data_db[,,,year] <<- sweep(data_db[,,,year], 1:2, popadjust[db_aid,,year], "*")
},

aging = function() {
    data[-c(1,pAG),,,year] <<- data[-(pAG-1:0),,,year-1]
    data[pAG,,,year] <<- data[pAG,,,year-1] + data[pAG-1,,,year-1] # open ended
    if (MODEL==2)  # @debut age the debut pop
        data_db[-1,,,year] <<- data_db[-pDB,,,year-1]
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
        data_db[1,,hivn.idx,year] <<- healthy
        data_db[1,,hivp.idx,year] <<- positiv
    }
    if (MODEL!=0) {
        entrantprev[year]        <<- sum(positiv)/sum(healthy+positiv)
        hivp_entrants_out[,year] <<- positiv
    }
},

sexual_debut = function() {
    debut_now            <- sweep(data_db[,,,year], 1:2, p$db_pr, "*")
    data_db[,,,year]     <<- data_db[,,,year]     - debut_now
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
    if (MODEL==2) {
      death_db <- sweep(data_db[,,,year], 1:2, 1 - p$Sx[db_aid,,year], "*")
      data_db[,,,year] <<- data_db[,,,year] - death_db
    }
    data[,,,year] <<- data[,,,year] - death_now
    natdeaths[,,year] <<- rowSums(death_now,,2)
},

migration = function() {
    netmigsurv <- p$netmigr[,,year] * (1 + p$Sx[,,year]) / 2
    mr.prob <- 1 + netmigsurv / rowSums(data[,,,year],,2)
    if (MODEL != 0) {
      hiv_mr_prob <<- sumByAGs(mr.prob * data[,,hivp.idx,year], ag.idx) / 
                      sumByAGs(          data[,,hivp.idx,year], ag.idx)
      hiv_mr_prob[is.nan(hiv_mr_prob)] <<- 0
    }
    if (MODEL == 2)
      data_db[,,,year] <<- sweep(data_db[,,,year], 1:2, mr.prob[db_aid,], "*")
    data[,,,year] <<- sweep(data[,,,year], 1:2, mr.prob, "*")
},

update_fertile = function() { # only on active pop
    update_active_pop_to(year)
    two_years  <- data_active + get_active_pop_in(year-1)
    birth_age  <<- rowSums(two_years[p.fert.idx, f.idx,])/2 * p$asfr[,year]
    birth_age  <<- rowSums(two_years[p.fert.idx, f.idx,])/2 * p$asfr[,year]
    birth_agrp <<- sumByAG(birth_age, ag.idx, TRUE, p.fert.idx)
    births     <- p$srb[,year] * sum(birth_agrp)
    if (year + AGE_START <= PROJ_YEARS)
      birthslag[, year + AGE_START - 1] <<- births
},

cal_prev_pregant = function(hivpop, artpop) { # only on active pop
    years   <- year - 1:0
    update_active_pop_to(year)
    two_years  <- data_active + get_active_pop_in(year-1)
    meanWomen <- two_years[p.fert.idx, f.idx, hivn.idx] / 2
    hivn <- sumByAG(meanWomen, ag.idx, TRUE, p.fert.idx)
    hivp <- rowMeans(hivpop$get(AG=h.fert.idx, NG=f.idx, YEAR=years),,2)
    art  <- rowMeans(artpop$get(AG=h.fert.idx, NG=f.idx, YEAR=years),,3)
    pregprev <- sum(birth_agrp * 
      (1 - hivn / (hivn + colSums(p$frr_cd4[,,year] * hivp) + 
      colSums(p$frr_art[,,,year] * art,,2)))) / sum(birth_age)
    pregprevlag[year + AGE_START - 1] <<- pregprev
},

save_prev_n_inc = function() {
    prev15to49[year]  <<- sum(data[p.age15to49.idx,, hivp.idx, year]) / 
                          sum(data[p.age15to49.idx,,         , year])
    update_active_pop_to(year-1)
    incid15to49[year] <<- incid15to49[year] /
                          sum(data_active[p.age15to49.idx,,hivn.idx])
    prev[year]  <<- sum(data[,,hivp.idx,year]) / sum(data[,,,year])
    incid[year] <<- incid15to49[year] / sum(data_active[,,hivn.idx]) # toBfixed
}, 

infect_mix = function(hivpop, artpop, ii) {
    ts <- (year-2)/DT + ii
    update_active_pop_to(year)
    hiv_treated       <- sweep(data_active[,,hivp.idx], 2, artcov, '*')
    hiv_not_treated   <- data_active[,,hivp.idx] - hiv_treated
    transm_prev <- (hiv_not_treated + hiv_treated * (1 - p$relinfectART)) / 
                    rowSums(data_active,,2) # prevalence adjusted for art
    # +intervention effects and time epidemic start
    w  <- p$iota * (p$proj.steps[ts] == p$tsEpidemicStart)
    transm_prev <- rvec[ts] * transm_prev + w
    # sweep over sexual mixing matrices
    ir_m <- rowSums(sweep(p$mat_m, 2, transm_prev[, f.idx], "*")) # male
    ir_f <- rowSums(sweep(p$mat_f, 2, transm_prev[, m.idx], "*")) # female
    irmf   <- cbind(ir_m, ir_f)
    # if (exists("f_fun", fp)) # that fun
    #   ir <- ir * fp$f_fun
    infections.ts <- irmf * data_active[,,hivn.idx]

    incrate15to49_ts[,,ts] <<- transm_prev
    prev15to49_ts[ts] <<- prevlast <<- sum(data[,,hivp.idx,year])/sum(data[,,,year])
    infections.ts
},

infect_spec = function(hivpop, artpop, time_step) {
    ts    <- (year-2) / DT + time_step
    dt_ii <- 1 - DT * (time_step - 1) # transition of population in 1 year
    first_age  <- p.age15to49.idx[1]
    first_agrp <- h.age15to49.idx[1]
    last_age   <- tail(p.age15to49.idx, 1) + 1
    last_agrp  <- tail(h.age15to49.idx, 1) + 1
    update_active_pop_to(year)
    hivn.ii <- sum(data_active[p.age15to49.idx,,hivn.idx]) - 
               sum(data_active[first_age,,hivn.idx]) * dt_ii + 
               sum(data_active[last_age,,hivn.idx]) * dt_ii
    hivp.ii <- sum(data_active[p.age15to49.idx,,hivp.idx]) - 
               sum(data_active[first_age,,hivp.idx]) * dt_ii + 
               sum(data_active[last_age,,hivp.idx]) * dt_ii
    art.ii <- sum(artpop$get(AG = h.age15to49.idx, YEAR = year))
    if (sum(hivpop$get(AG = first_agrp, YEAR = year)) + 
        sum(artpop$get(AG = first_agrp, YEAR = year)) > 0) {
      art_first <- colSums(artpop$data[,,first_agrp,,year],,2)
      hiv_first <- colSums(hivpop$data[ ,first_agrp,,year],,1)
      art.ii  <- art.ii - sum(data_active[first_age,,hivp.idx] * 
        art_first / (hiv_first + art_first) ) * dt_ii
    }
    if (sum(hivpop$data[ ,last_agrp,,year]) + 
        sum(artpop$data[,,last_agrp,,year]) > 0) {
      art_last <- colSums(artpop$data[,,last_agrp,,year],,2)
      hiv_last <- colSums(hivpop$data[ ,last_agrp,,year],,1)
      art.ii <- art.ii + sum(data_active[last_age,,hivp.idx] * 
        art_last / (hiv_last + art_last) ) * dt_ii
    }
    transm_prev <- (hivp.ii - art.ii*(1 - p$relinfectART)) / (hivn.ii+hivp.ii)
    w <- p$iota * (p$proj.steps[ts] == p$tsEpidemicStart)
    incrate15to49.ts <- rvec[ts] * transm_prev + w
    
    # Incidence: male = negative / female negative * sexratio + male negative; 
    #          female = male * sexratio
    sus_by_age_sex <- data_active[p.age15to49.idx,,hivn.idx]
    adj_sex <- sum(sus_by_age_sex) /
      ( sum(data_active[p.age15to49.idx,m.idx,hivn.idx]) +
        sum(data_active[p.age15to49.idx,f.idx,hivn.idx]) * 
        p$incrr_sex[year] )
    sexinc15to49.ts <- incrate15to49.ts * c(1, p$incrr_sex[year]) * adj_sex

    # New infections distributed by age: ratio age_i/ 25-29 age
    adj_age    <- sexinc15to49.ts / (
                  colSums(sus_by_age_sex * p$incrr_age[p.age15to49.idx,,year]) / 
                  colSums(sus_by_age_sex) )
    agesex.inc <- sweep(p$incrr_age[,,year], 2, adj_age, "*")

    ## Adjust age-specific incidence among men for circumcision coverage
    agesex.inc[, m.idx] <- agesex.inc[, m.idx] * 
                           (1 - p$circ_incid_rr * p$circ_prop[,year])
    infections.ts <- agesex.inc * data_active[,,hivn.idx]
    # saving
    incrate15to49_ts[ts] <<- incrate15to49.ts
    prev15to49_ts[ts] <<- prevlast <<- hivp.ii / (hivn.ii+hivp.ii)
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