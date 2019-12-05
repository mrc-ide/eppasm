# Annualized number of new infections
# infect_spec <- function(pop, hivpop, artpop, ii, fp)
# K moved to popClass method
infectFns <- c(
infect_mix = function(hivpop, artpop, ii) {
    ts <- (year-2)/DT + ii
    update_active_pop_to(year)
    hiv_treated       <- sweep(data_active[,,hivp.idx], 2, artcov, '*')
    hiv_not_treated   <- data_active[,,hivp.idx] - hiv_treated
    transm_prev <- (hiv_not_treated + hiv_treated * (1 - p$relinfectART)) / 
                    rowSums(data[,,,year],,2) # prevalence adjusted for art
    # +intervention effects and time epidemic start
    w  <- p$iota * (p$proj.steps[ts] == p$tsEpidemicStart)
    transm_prev <- rvec[ts] * transm_prev + w
    # sweep over sexual mixing matrices
    ir_m <- rowSums(sweep(p$mat_m, 2, transm_prev[, f.idx], "*")) # male
    ir_f <- rowSums(sweep(p$mat_f, 2, transm_prev[, m.idx], "*")) # female
    irmf <- cbind(ir_m, ir_f)
    # if (exists("f_fun", fp)) # that fun
    #   ir <- ir * fp$f_fun

    # Scale to age pattern of IR
    max_rr_sex <- apply(p$incrr_age[,,year],2,max)
    scaled_rr_age <- sweep(p$incrr_age[,,year], 2, max_rr_sex, '/')
    irmf <- irmf * scaled_rr_age
  
    # Scale sex pattern IR
    incidence_estimated <- colSums(irmf * data_active[,,hivn.idx]) / 
                           colSums(data_active[,,hivn.idx])
    IRR_FM_estimated <- incidence_estimated[2]/incidence_estimated[1]
    if (!is.na(IRR_FM_estimated)) {
      irmf[, f.idx] <- irmf[, f.idx]*p$incrr_sex[year]
      irmf[, m.idx] <- irmf[, m.idx]*IRR_FM_estimated
    }
  
    infections.ts <- irmf * data_active[,,hivn.idx]

    incrate15to49_ts[,,ts] <<- transm_prev
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

    if (MODEL==1) {
      adj_sex <- sum(sus_age_sex) /
        ( sum(data_active[p.age15to49.idx,m.idx,hivn.idx]) +
          sum(data_active[p.age15to49.idx,f.idx,hivn.idx]) * 
          p$incrr_sex[year] )
      sexinc15to49.ts <- inc_rate * c(1, p$incrr_sex[year]) * adj_sex
      # New infections distributed by age: ratio age_i/ 25-29 age
      adj_age    <- sexinc15to49.ts /
        (colSums(sus_age_sex * p$incrr_age[p.age15to49.idx,,year]) / colSums(sus_age_sex))
      agesex.inc <- sweep(p$incrr_age[,,year], 2, adj_age, "*")
    }

    if (MODEL==2) {
      S_F       = sum(data_active[p.age15to49.idx, f.idx, hivn.idx])
      S_M       = sum(data_active[p.age15to49.idx, m.idx, hivn.idx])
      H_M_minus = sum(data[p.age15to49.idx, m.idx, hivn.idx, year])
      H_F_minus = sum(data[p.age15to49.idx, f.idx, hivn.idx, year])
      sex_inc = c(0,0)
      sex_inc[m.idx] = (inc_rate * (S_M + S_F) * H_M_minus) /
                       (S_M * (p$incrr_sex[year] * H_F_minus + H_M_minus))
      sex_inc[f.idx] = sex_inc[m.idx] * p$incrr_sex[year] * S_M * H_F_minus / (S_F * H_M_minus)
      IR_tsa_hat = p$incrr_age[,,year] * data_active[,,hivn.idx] / data[,,hivn.idx,year]
      incrr_adj = sweep(p$incrr_age[,,year]^2, 2, IR_tsa_hat[10,], '*') /
                  sweep(IR_tsa_hat, 2, p$incrr_age[10,,year], '*')
      incrr_adj[is.na(incrr_adj)] = 0
      inc_sex_hat = colSums(sus_age_sex * incrr_adj[p.age15to49.idx,]) / colSums(sus_age_sex)
      adj_age = sex_inc / inc_sex_hat
      agesex.inc <- sweep(incrr_adj, 2, adj_age, "*")
    }

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

## Beers coefficients to distribute infections from 5-year age groups to single-year of age
create_beers <- function(n5yr){

  ## Beer's coefficients for disaggregating 5 year age groups into
  ## single-year age groups (from John Stover)
  Afirst <- rbind(c(0.3333, -0.1636, -0.0210,  0.0796, -0.0283),
                  c(0.2595, -0.0780,  0.0130,  0.0100, -0.0045),
                  c(0.1924,  0.0064,  0.0184, -0.0256,  0.0084),
                  c(0.1329,  0.0844,  0.0054, -0.0356,  0.0129),
                  c(0.0819,  0.1508, -0.0158, -0.0284,  0.0115))
  Asecond <- rbind(c( 0.0404,  0.2000, -0.0344, -0.0128,  0.0068),
                   c( 0.0093,  0.2268, -0.0402,  0.0028,  0.0013),
                   c(-0.0108,  0.2272, -0.0248,  0.0112, -0.0028),
                   c(-0.0198,  0.1992,  0.0172,  0.0072, -0.0038),
                   c(-0.0191,  0.1468,  0.0822, -0.0084, -0.0015))
  Amid <- rbind(c(-0.0117,  0.0804,  0.1570, -0.0284,  0.0027),
                c(-0.0020,  0.0160,  0.2200, -0.0400,  0.0060),
                c( 0.0050, -0.0280,  0.2460, -0.0280,  0.0050),
                c( 0.0060, -0.0400,  0.2200,  0.0160, -0.0020),
                c( 0.0027, -0.0284,  0.1570,  0.0804, -0.0117))
  Apenult <- rbind(c(-0.0015, -0.0084,  0.0822,  0.1468, -0.0191),
                   c(-0.0038,  0.0072,  0.0172,  0.1992, -0.0198),
                   c(-0.0028,  0.0112, -0.0248,  0.2272, -0.0108),
                   c( 0.0013,  0.0028, -0.0402,  0.2268,  0.0093),
                   c( 0.0068, -0.0128, -0.0344,  0.2000,  0.0404))
  Aultim <- rbind(c( 0.0115, -0.0284, -0.0158,  0.1508,  0.0819),
                  c( 0.0129, -0.0356,  0.0054,  0.0844,  0.1329),
                  c( 0.0084, -0.0256,  0.0184,  0.0064,  0.1924),
                  c(-0.0045,  0.0100,  0.0130, -0.0780,  0.2595),
                  c(-0.0283,  0.0796, -0.0210, -0.1636,  0.3333))

  A <- do.call(rbind,
               c(list(cbind(Afirst, matrix(0, 5, n5yr-5)),
                      cbind(Asecond, matrix(0, 5, n5yr-5))),
                 lapply(0:(n5yr-6), function(i) cbind(matrix(0, 5, i), Amid, matrix(0, 5, (n5yr-5)-i))),
                 list(cbind(matrix(0, 5, n5yr-6), Apenult, matrix(0, 5, 1)),
                      cbind(matrix(0, 5, n5yr-6), Aultim, matrix(0, 5, 1)),
                      c(rep(0, n5yr-1), 1))))
  return(round(A, 4))
}

prepare_irr_model <- function(fp){
  irr <- list()
  return(fp)
}
  

################################################
####  Prior for incidence rate ratio model  ####
################################################

## RW2 model
NPARAM_RW2 <- 13

sexincrr.pr.mean <- log(1.38)
sexincrr.pr.sd <- 0.2

mf_transm_rr.pr.mean <- log(1.9)
mf_transm_rr.pr.sd <- 0.3  # change default to 0.3

## ageincrr.pr.mean <- c(-1.40707274, -0.23518703, 0.69314718, 0.78845736, -0.39975544, -0.70620810, -0.84054571, -0.02101324, -0.16382449, -0.37914407, -0.59639985, -0.82038300)
## ageincrr.pr.sd <- 0.5

## Informative priors based on estimates for 11 countries with 3+ surveys
ageincrr.pr.mean <- c(-1.4, -0.28, 0.3, 0.3, -0.3, -0.6, -0.2, 0.05, -0.4, -0.45, -0.6, -0.7)
ageincrr.pr.sd <- c(0.5, 0.4, 0.23, 0.3, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2)

NPARAM_LININCRR <- 6
## incrr_trend_mean <- c(0, 0, 0, 0, 0, 0)
## incrr_trend_sd <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

## Informative priors based on estimates for 11 countries with 3+ surveys
incrr_trend_mean <- c(0.0, 0.035, -0.02, -0.09, -0.016, -0.06)
incrr_trend_sd <- c(0.07, 0.07, 0.1, 0.1, 0.08, 0.08)

## Incidence rate ratios for age 50 plus, relative to 15-49
incrr_50plus_logdiff <- cbind(male   = log(0.493510) - log(c(0.358980, 0.282400, 0.259240, 0.264920, 0.254790, 0.164140, 0.000000)),
                              female = log(0.440260) - log(c(0.336720, 0.239470, 0.167890, 0.146590, 0.171350, 0.000000, 0.000000)))

## Beers coefficient matrix
beers_Amat <- create_beers(17)[16:81, 4:17]

getnparam_incrr <- function(fp){
  value <- switch(as.character(fp$fitincrr),
                  "FALSE" = 0,
                  "TRUE" = NPARAM_RW2,
                  linincrr = NPARAM_RW2+NPARAM_LININCRR,
                  lognorm = 7,
                  relbehav = NPAR_RELBEHAV, NA)
  if(is.na(value))
    stop(paste0("fitincrr model '", fp$fitincrr, "' is not recognized"))
  value
}

transf_incrr <- function(theta_incrr, param, fp){

  incrr_nparam <- getnparam_incrr(fp)
  
  if(fp$incidmod == "eppspectrum"){
    param$incrr_sex <- fp$incrr_sex
    param$incrr_sex[] <- exp(theta_incrr[1])
  } else if(fp$incidmod == "transm") {
    param$mf_transm_rr <- rep(exp(theta_incrr[1]), fp$ss$PROJ_YEARS)
  }

  if(fp$fitincrr %in% c(TRUE ,"linincrr")){

    ## param$sigma_agepen <- exp(theta_incrr[incrr_nparam])
    param$sigma_agepen <- 0.4
    
    param$logincrr_age <- array(0, c(14, 2))
    param$logincrr_age[c(1:2, 4:7), ] <- theta_incrr[2:13]
    param$logincrr_age[8:14, ] <- sweep(-incrr_50plus_logdiff, 2,
                                        param$logincrr_age[7, ], "+")

    ## Smooth 5-year age group IRRs to 1-year IRRs
    incrr_age <- beers_Amat %*% exp(param$logincrr_age)
    incrr_age[incrr_age < 0] <- 0
    
    param$incrr_age <- array(incrr_age, c(dim(incrr_age), fp$ss$PROJ_YEARS))

    years <- with(fp$ss, proj_start+1:PROJ_YEARS-1)
    if(fp$fitincrr == "linincrr"){
      par <- theta_incrr[NPARAM_RW2+1:NPARAM_LININCRR]
      param$logincrr_trend <- par
      sexadjust <- approx(c(2002, 2007, 2012), c(-5, 0, 5)*c(par[1], 0, par[2]), years, rule=2)$y
      if(fp$incidmod == "eppspectrum")
        param$incrr_sex <- param$incrr_sex * exp(sexadjust)
      else if(fp$incidmod == "transm")
        param$mf_transm_rr <- param$mf_transm_rr * exp(sexadjust)

      ## adjustment to age IRRs among 15-24
      m15to24_adjust <- approx(c(2002, 2007, 2012), c(-5, 0, 5)*c(par[3], 0, par[4]), years, rule=2)$y
      f15to24_adjust <- approx(c(2002, 2007, 2012), c(-5, 0, 5)*c(par[5], 0, par[6]), years, rule=2)$y
      param$incrr_age[1:10,,] <- sweep(param$incrr_age[1:10,,,drop=FALSE], 2:3, exp(rbind(m15to24_adjust, f15to24_adjust)), "*")      
    }

  } else if(fp$fitincrr=="lognorm"){
    param$logincrr_age <- cbind(calc_lognorm_logagerr(theta_incrr[2:4]),
                                calc_lognorm_logagerr(theta_incrr[5:7]))

    ## Smooth 5-year age group IRRs to 1-year IRRs
    incrr_age <- beers_Amat %*% exp(param$logincrr_age)
    incrr_age[incrr_age < 0] <- 0
    
    param$incrr_age <- array(incrr_age, c(dim(incrr_age), fp$ss$PROJ_YEARS))

  } else if(fp$fitincrr == "relbehav"){

    stop("relbehav is not implemented currently")
    
    ## par <- theta_incrr[2:incrr_nparam]
    ## param$adjustpar <- par
    ## logadjust1 <- cbind(approx(c(17, 27, 38, 49), c(par[1], 0, cumsum(par[2:3])), xout=15:80, rule=2)$y,
    ##                     approx(c(17, 27, 38, 49), c(par[4], 0, cumsum(par[5:6])), xout=15:80, rule=2)$y)
    
    ## logadjust2 <- cbind(approx(c(17, 27, 38, 49), c(par[1]+par[7], 0, cumsum(par[2:3])), xout=15:80, rule=2)$y,
    ##                     approx(c(17, 27, 38, 49), c(par[4]+par[8], 0, cumsum(par[5:6])), xout=15:80, rule=2)$y)
    
    ## BREAK_YEAR <- 36
    ## param$incrr_age <- fp$logrelbehav
    ## param$incrr_age[,,1:(BREAK_YEAR-1)] <- exp(sweep(fp$logrelbehav[,,1:(BREAK_YEAR-1)], 1:2, logadjust1, "+"))
    ## param$incrr_age[,,BREAK_YEAR:fp$SIM_YEARS] <- exp(sweep(fp$logrelbehav[,,BREAK_YEAR:fp$SIM_YEARS], 1:2, logadjust2, "+"))
  }



  
  return(param)
}

lprior_incrr <- function(theta_incrr, fp){

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  lpr <- 0
  
  if(fp$incidmod == "eppspectrum")
    lpr <- lpr + dnorm(theta_incrr[1], sexincrr.pr.mean, sexincrr.pr.sd, log=TRUE)
  else if(fp$incidmod == "transm")
    lpr <- lpr + dnorm(theta_incrr[1], mf_transm_rr.pr.mean, mf_transm_rr.pr.sd, log=TRUE)

  if(fp$fitincrr %in% c(TRUE, "linincrr")){
    lpr <- lpr + sum(dnorm(theta_incrr[2:13], ageincrr.pr.mean, ageincrr.pr.sd, log=TRUE))
    ## dnorm(theta_incrr[14], -1, 0.7, log=TRUE)

    if(fp$fitincrr == "linincrr"){
      lpr <- lpr+sum(dnorm(theta_incrr[NPARAM_RW2+1:NPARAM_LININCRR], incrr_trend_mean, incrr_trend_sd, log=TRUE))
    }

    } else if(fp$fitincrr=="lognorm"){
      lpr <- lpr +
        sum(dnorm(theta_incrr[c(2,5)], lognorm.a0.pr.mean, lognorm.a0.pr.sd, log=TRUE)) +
        sum(dnorm(theta_incrr[c(3,6)], lognorm.meanlog.pr.mean, lognorm.meanlog.pr.sd, log=TRUE)) +
        sum(dnorm(theta_incrr[c(4,7)], lognorm.logsdlog.pr.mean, lognorm.logsdlog.pr.sd, log=TRUE))
    } else if(fp$fitincrr=="relbehav"){
      lpr <- lpr + sum(dnorm(theta_incrr[2:NPAR_RELBEHAV], 0, relbehav_adjust_sd, log=TRUE));
    }

  return(lpr)
}

sample_incrr <- function(n, fp){
  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  incrr_nparam <- getnparam_incrr(fp)
  mat <- matrix(NA, n, incrr_nparam)
  
  if(fp$incidmod == "eppspectrum")
    mat[,1] <- rnorm(n, sexincrr.pr.mean, sexincrr.pr.sd)
  else if(fp$incidmod == "transm")
    mat[,1] <- rnorm(n, mf_transm_rr.pr.mean, mf_transm_rr.pr.sd)
  
  if(fp$fitincrr %in% c(TRUE, "linincrr")){
    mat[,2:13] <- t(matrix(rnorm(n*12, ageincrr.pr.mean, ageincrr.pr.sd), nrow=12))
    ## mat[,14] <- rnorm(n, -1, 0.7)  # log variance of ageincrr difference penalty
    if(fp$fitincrr == "linincrr")
      mat[,NPARAM_RW2+1:NPARAM_LININCRR] <- t(matrix(rnorm(n*NPARAM_LININCRR, incrr_trend_mean, incrr_trend_sd), nrow=NPARAM_LININCRR))
  } else if(fp$fitincrr=="lognorm"){
    mat[,c(2,5)] <- t(matrix(rnorm(n*2, lognorm.a0.pr.mean, lognorm.a0.pr.sd), nrow=2))
    mat[,c(3,6)] <- t(matrix(rnorm(n*2, lognorm.meanlog.pr.mean, lognorm.meanlog.pr.sd), nrow=2))
    mat[,c(4,7)] <- t(matrix(rnorm(n*2, lognorm.logsdlog.pr.mean, lognorm.logsdlog.pr.sd), nrow=2))
  } else if(fp$fitincrr=="relbehav"){
    incrr_nparam <- NPAR_RELBEHAV
    mat[,2:NPAR_RELBEHAV] <- rnorm(n*(NPAR_RELBEHAV-1), 0, relbehav_adjust_sd)
  }

  return(mat)
}

ldsamp_incrr <- lprior_incrr