###########################
####  EPP model prior  ####
###########################

ldinvgamma <- function(x, alpha, beta){
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
  return(log.density)
}

bayes_lmvt <- function(x, shape, rate){
  mvtnorm::dmvt(x, sigma=diag(length(x)) / (shape / rate), df=2*shape, log=TRUE)
}

bayes_rmvt <- function(n, d, shape, rate){
  mvtnorm::rmvt(n, sigma=diag(d) / (shape / rate), df=2*shape)
}

## Binomial distribution log-density permitting non-integer counts
ldbinom <- function(x, size, prob){
  lgamma(size+1) - lgamma(x+1) - lgamma(size-x+1) + x*log(prob) + (size-x)*log(1-prob)
}

## r-spline prior parameters

## tau2.prior.rate <- 0.5  # initial sampling distribution for tau2 parameter
tau2_init_shape <- 3
tau2_init_rate <- 4

tau2_prior_shape <- 0.001   # Inverse gamma parameter for tau^2 prior for spline
tau2_prior_rate <- 0.001
muSS <- 1/11.5               #1/duration for r steady state prior

rw_prior_shape <- 300
rw_prior_rate <- 1.0


## r-trend prior parameters
t0.unif.prior <- c(1970, 1990)
## t1.unif.prior <- c(10, 30)
## logr0.unif.prior <- c(1/11.5, 10)
t1.pr.mean <- 20.0
t1.pr.sd <- 4.5
logr0.pr.mean <- 0.42
logr0.pr.sd <- 0.23
## rtrend.beta.pr.mean <- 0.0
## rtrend.beta.pr.sd <- 0.2
rtrend.beta.pr.mean <- c(0.46, 0.17, -0.68, -0.038)
rtrend.beta.pr.sd <- c(0.12, 0.07, 0.24, 0.009)


###########################################
####                                   ####
####  Site level ANC data (SS and RT)  ####
####                                   ####
###########################################

ancbias.pr.mean <- 0.15
ancbias.pr.sd <- 1.0
vinfl.prior.rate <- 1/0.015

ancrtsite.beta.pr.mean <- 0
## ancrtsite.beta.pr.sd <- 1.0
ancrtsite.beta.pr.sd <- 0.05
## ancrtsite.vinfl.pr.rate <- 1/0.015


#' Prepare site-level ANC prevalence data for EPP random-effects likelihood
#'
#' @param ancitedat data.frame of site-level ANC data
#' @param fp fixed parameter input list, including state space
prepare_ancsite_likdat <- function(ancsitedat, fp){

  df <- ancsitedat
  anchor.year <-   anchor.year <- fp$ss$proj_start
  
  df$aidx <- df$age - fp$ss$AGE_START + 1L
  df$yidx <- df$year - anchor.year + 1

  ## Calculate probit transformed prevalence and variance approximation
  df$pstar <- (df$prev * df$n + 0.5) / (df$n + 1)
  df$W <- qnorm(df$pstar)
  df$v <- 2 * pi * exp(df$W^2) * df$pstar * (1 - df$pstar) / df$n

  df$type <- factor(df$type, c("ancss", "ancrt"))

  ## List of all unique agegroup / year combinations for which prevalence is needed
  datgrp <- unique(df[c("aidx", "yidx", "agspan")])
  datgrp$qMidx <- seq_len(nrow(datgrp))

  ## Indices for accessing prevalence offset from datgrp
  df <- merge(df, datgrp)

  ## Indices for prevalence group
  df <- df[c("site", "year", "used", "type", "age", "agspan",
             "n", "prev", "pstar", "W", "v", "aidx", "yidx", "qMidx")]
  df <- df[order(df$site, df$year, -df$agspan, df$age),]
  df$muidx <- seq_len(nrow(df))

  ## Design matrix for fixed effects portion 
  Xancsite <- model.matrix(~type, df)

  W.lst <- split(df$W, factor(df$site))
  v.lst <- split(df$v, factor(df$site))
  muidx.lst <- split(df$muidx, factor(df$site))

  list(df = df,
       datgrp = datgrp,
       Xancsite = Xancsite,
       muidx.lst = muidx.lst)
}

ll_ancsite <- function(mod, fp, coef=c(0, 0), vinfl=0, dat){
  
  df <- dat$df
  
  qM <- suppressWarnings(qnorm(agepregprev(mod, fp, dat$datgrp$aidx, dat$datgrp$yidx, dat$datgrp$agspan)))
  mu <- qM[df$qMidx] + dat$Xancsite %*% coef
  d <- df$W - mu
  v <- df$v + vinfl
  
  d.lst <- lapply(dat$muidx.lst, function(idx) d[idx])
  v.lst <- lapply(dat$muidx.lst, function(idx) v[idx])
  
  log(anclik::anc_resid_lik(d.lst, v.lst))
}


#############################################
####                                     ####
####  ANCRT census likelihood functions  ####
####                                     ####
#############################################

## prior parameters for ANCRT census
log_frr_adjust.pr.mean <- 0
## ancrtcens.bias.pr.sd <- 1.0
log_frr_adjust.pr.sd <- 0.2
ancrtcens.vinfl.pr.rate <- 1/0.015

prepare_ancrtcens_likdat <- function(dat, fp){

  anchor.year <- fp$ss$proj_start
  
  x.ancrt <- (dat$prev*dat$n+0.5)/(dat$n+1)
  dat$W.ancrt <- qnorm(x.ancrt)
  dat$v.ancrt <- 2*pi*exp(dat$W.ancrt^2)*x.ancrt*(1-x.ancrt)/dat$n

  if(!exists("age", dat))
    dat$age <- rep(15, nrow(dat))

  if(!exists("agspan", dat))
    dat$agspan <- rep(35, nrow(dat))

  dat$aidx <- dat$age - fp$ss$AGE_START + 1
  dat$yidx <- dat$year - anchor.year + 1
  
  return(dat)
}

ll_ancrtcens <- function(mod, dat, fp){
  if(!nrow(dat))
    return(0)
  qM.prev <- suppressWarnings(qnorm(agepregprev(mod, fp, dat$aidx, dat$yidx, dat$agspan)))
  if(any(is.na(qM.prev)))
    return(-Inf)
  sum(dnorm(dat$W.ancrt, qM.prev, sqrt(dat$v.ancrt + fp$ancrtcens.vinfl), log=TRUE))
}



###################################
####  Age/sex incidence model  ####
###################################

## log-normal age incrr prior parameters
lognorm.a0.pr.mean <- 10
lognorm.a0.pr.sd <- 5

lognorm.meanlog.pr.mean <- 3
lognorm.meanlog.pr.sd <- 2

lognorm.logsdlog.pr.mean <- 0
lognorm.logsdlog.pr.sd <- 1

relbehav_adjust_sd <- 0.25
NPAR_RELBEHAV <- 9

calc_lognorm_logagerr <- function(par, a=2.5+5*3:16, b=27.5){
  dlnorm(a-par[1], par[2], exp(par[3]), log=TRUE) - dlnorm(b-par[1], par[2], exp(par[3]), log=TRUE)
}



fnCreateLogAgeSexIncrr <- function(logrr, fp){


  logincrr.theta.idx <- c(7:10, 12:20)
  logincrr.fixed.idx <- 11

  lastidx <- length(fp$proj.steps)
  fixed.age50p.logincrr <- log(fp$agesex.incrr.ts[,age50plus.idx,lastidx]) - log(fp$agesex.incrr.ts[,age45.idx,lastidx])

  logincrr.theta <- tail(theta, length(logincrr.theta.idx))
  logincrr.agesex <- array(-Inf, c(NG, AG))
  logincrr.agesex[logincrr.fixed.idx] <- 0
  logincrr.agesex[logincrr.theta.idx] <- logincrr.theta
  logincrr.agesex[,age50plus.idx] <- logincrr.agesex[,age45.idx] + fixed.age50p.logincrr

  return(logincrr.agesex)
}

create_natmx_param <- function(theta_natmx, fp){

  ## linear trend in logmx
  par <- list(natmx_b0 = theta_natmx[1],
              natmx_b1 = theta_natmx[2])
  par$Sx <- with(fp$natmx, exp(-exp(outer(logmx0, par$natmx_b0 + natmx_b1*x, "+"))))
  return(par)
}


fnCreateParam <- function(theta, fp){

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if(!exists("eppmod", where = fp))  # backward compatibility
    fp$eppmod <- "rspline"
  
  if(fp$eppmod %in% c("rspline", "logrspline", "ospline", "logospline", "logrw")){

    epp_nparam <- fp$numKnots+1

    if(fp$eppmod %in% c("rspline", "logrspline")){
      u <- theta[1:fp$numKnots]
      if(fp$rtpenord == 2){
        beta <- numeric(fp$numKnots)
        beta[1] <- u[1]
        beta[2] <- u[1]+u[2]
        for(i in 3:fp$numKnots)
          beta[i] <- -beta[i-2] + 2*beta[i-1] + u[i]
      } else # first order penalty
        beta <- cumsum(u)
    } else if(fp$eppmod %in% c("ospline", "logospline", "logrw"))
      beta <- theta[1:fp$numKnots]

    param <- list(beta = beta,
                  rvec = as.vector(fp$rvec.spldes %*% beta))
    
    if(fp$eppmod %in% c("logrspline", "logospline", "logrw"))

      param$rvec <- exp(param$rvec)

    if(exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      param$iota <- exp(param$rvec[fp$proj.steps == fp$tsEpidemicStart] * theta[fp$numKnots+1])
    else
      param$iota <- transf_iota(theta[fp$numKnots+1], fp)

  } else if(fp$eppmod == "rlogistic") {
    epp_nparam <- 5
    par <- theta[1:4]
    par[3] <- exp(theta[3])
    ## par[1:3] <- exp(par[1:3])
    param <- list()
    param$rvec <- exp(rlogistic(fp$proj.steps, par))
    ## param$rvec <- rlogistic(fp$proj.steps, par)
    param$iota <- transf_iota(theta[5], fp)
  } else if(fp$eppmod == "rtrend"){ # rtrend
    epp_nparam <- 7
    param <- list(tsEpidemicStart = fp$proj.steps[which.min(abs(fp$proj.steps - (round(theta[1]-0.5)+0.5)))], # t0
                  rtrend = list(tStabilize = round(theta[1]-0.5)+0.5+round(theta[2]),  # t0 + t1
                                r0 = exp(theta[3]),              # r0
                                beta = theta[4:7]))
  } else {
    epp_nparam <- fp$rt$n_param+1
    param <- list()
    param$rvec <- create_rvec(theta[1:fp$rt$n_param], fp$rt)
    param$iota <- transf_iota(theta[fp$rt$n_param+1], fp)
  }

  if(fp$ancsitedata){
    param$ancbias <- theta[epp_nparam+1]
    if(!exists("v.infl", where=fp)){
      anclik_nparam <- 2
      param$v.infl <- exp(theta[epp_nparam+2])
    } else
      anclik_nparam <- 1
  }
  else
    anclik_nparam <- 0


  paramcurr <- epp_nparam+anclik_nparam
  if(exists("ancrt", fp) && fp$ancrt %in% c("census", "both")){
    param$log_frr_adjust <- theta[paramcurr+1]
    param$frr_cd4 <- fp$frr_cd4 * exp(param$log_frr_adjust)
    param$frr_art <- fp$frr_art * exp(param$log_frr_adjust)
    ## param$frr_art[1:2,,,] <- param$frr_art[1:2,,,] * exp(param$log_frr_adjust)

    if(!exists("ancrtcens.vinfl", fp)){
      param$ancrtcens.vinfl <- exp(theta[paramcurr+2])
      paramcurr <- paramcurr+2
    } else
      paramcurr <- paramcurr+1
  }
  if(exists("ancrt", fp) && fp$ancrt %in% c("site", "both")){
    param$ancrtsite.beta <- theta[paramcurr+1]
    paramcurr <- paramcurr+1
    ## param$ancrtsite.vinfl <- exp(theta[length(theta)])
  }

  if(inherits(fp, "specfp")){
    if(exists("fitincrr", where=fp)){
      incrr_nparam <- getnparam_incrr(fp)
      if(incrr_nparam)
        param <- transf_incrr(theta[paramcurr+1:incrr_nparam], param, fp)
      paramcurr <- paramcurr+incrr_nparam
    }
  }

  if(exists("natmx", where=fp) && fp$fitmx==TRUE){
    natmx_nparam <- 3
    theta_natmx <- theta[paramcurr+1:natmx_nparam]
    paramcurr <- paramcurr+natmx_nparam
    
    b0 <- theta_natmx[1]
    b1 <- theta_natmx[2]/10
    mx_lsexrat <- theta_natmx[3]
    
      param$natmx_par <- list(b0=b0, b1=b1, mx_lsexrat=mx_lsexrat)
    param$Sx <- with(fp$natmx, exp(-exp(outer(sweep(logmx0, 2, c(0, mx_lsexrat), "+"), b0 + b1*x, "+"))))
  }
  
  return(param)
}



########################################################
####  Age specific prevalence likelihood functions  ####
########################################################


#' Prepare age-specific HH survey prevalence likelihood data
prepare_hhsageprev_likdat <- function(hhsage, fp){
  anchor.year <- floor(min(fp$proj.steps))

  hhsage$W.hhs <- qnorm(hhsage$prev)
  hhsage$v.hhs <- 2*pi*exp(hhsage$W.hhs^2)*hhsage$se^2
  hhsage$sd.W.hhs <- sqrt(hhsage$v.hhs)

  if(exists("deff_approx", hhsage))
    hhsage$n_eff <- hhsage$n/hhsage$deff_approx
  else
    hhsage$n_eff <- hhsage$n/hhsage$deff
  hhsage$x_eff <- hhsage$n_eff * hhsage$prev

  startage <- as.integer(sub("([0-9]*)-([0-9]*)", "\\1", hhsage$agegr))
  endage <- as.integer(sub("([0-9]*)-([0-9]*)", "\\2", hhsage$agegr))
    
  hhsage$sidx <- match(hhsage$sex, c("both", "male", "female")) - 1L
  hhsage$aidx <- startage - fp$ss$AGE_START+1L
  hhsage$yidx <- as.integer(hhsage$year - (anchor.year - 1))
  hhsage$agspan <- endage - startage + 1L

  return(subset(hhsage, aidx > 0))
}

#' Log likelihood for age-specific household survey prevalence
ll_hhsage <- function(mod, dat){
  qM.age <- suppressWarnings(qnorm(ageprev(mod, aidx = dat$aidx, sidx = dat$sidx, yidx = dat$yidx, agspan = dat$agspan)))
  if(any(is.na(qM.age)))
    return(-Inf)
  sum(dnorm(dat$W.hhs, qM.age, dat$sd.W.hhs, log=TRUE))
}


#' Log likelihood for age-specific household survey prevalence using binomial approximation
ll_hhsage_binom <- function(mod, dat){
  prevM.age <- suppressWarnings(ageprev(mod, aidx = dat$aidx, sidx = dat$sidx, yidx = dat$yidx, agspan = dat$agspan))
  if(any(is.na(prevM.age)) || any(prevM.age >= 1)) return(-Inf)
  ll <- sum(ldbinom(dat$x_eff, dat$n_eff, prevM.age))
  if(is.na(ll))
    return(-Inf)
  return(ll)
}



##########################################
####  Mortality likelihood functions  ####
##########################################

#' Prepare sibling history mortality likelihood data
#'
prepare_sibmx_likdat <- function(sibmxdat, fp){
  anchor.year <- floor(min(fp$proj.steps))
  nyears <- fp$ss$PROJ_YEARS
  NG <- fp$ss$NG
  AG <- fp$ss$pAG

  sibmxdat$sidx <- as.integer(sibmxdat$sex)
  sibmxdat$aidx <- sibmxdat$agegr - (fp$ss$AGE_START-1)
  sibmxdat$yidx <- sibmxdat$period - (anchor.year - 1)
  sibmxdat$tipsidx <- sibmxdat$tips+1L

  sibmxdat <- subset(sibmxdat, aidx > 0)

  sibmxdat$arridx <- sibmxdat$aidx + (sibmxdat$sidx-1)*AG + (sibmxdat$yidx-1)*NG*AG

  return(sibmxdat)
}

#' Log negative binomial density
#'
#' Log negative binomial density, mu parameterization
#'
#' Log-density of negative binomial distribution. Parameter names and
#' parameterization matches the 'mu' parameterization of \code{\link{dnbinom}}.
#'
#' @param x vector of number of events.
#' @param size dispersion parameter.
#' @param mu mean expected number of events.
ldnbinom <- function(x, size, mu){
  prob <- size/(size+mu)
  lgamma(x+size) - lgamma(size) - lgamma(x+1) + size*log(prob) + x*log(1-prob)
}



#' Log-likelihood for sibling history mortality data
#'
#' Calculate the log-likelihood for sibling history mortality data
#'
#' !!! NOTE: does not account for complex survey design
#'
#' @param mx Array of age/sex-specific mortality rates for each year, output
#'   from function \code{\link{agemx}}.
#' @param tipscoef Vector of TIPS (time preceding survey) coefficients for
#'   relative risk of underreporting deceased siblings.
#' @param theta Overdispersion of negative binomial distribution.
#' @param sibmx.dat Data frame consisting of sibling history mortality data.
ll_sibmx <- function(mx, tipscoef, theta, sibmx.dat){

  ## predicted deaths: product of predicted mortality, tips coefficient, and person-years
  mu.pred <- mx[sibmx.dat$arridx] * tipscoef[sibmx.dat$tipsidx] * sibmx.dat$pys

  return(sum(ldnbinom(sibmx.dat$deaths, theta, mu.pred)))
}


#########################################
####  Incidence likelihood function  ####
#########################################

#' Prepare household survey incidence likelihood data
prepare_hhsincid_likdat <- function(hhsincid, fp){
  anchor.year <- floor(min(fp$proj.steps))

  hhsincid$idx <- hhsincid$year - (anchor.year - 1)
  hhsincid$log_incid <- log(hhsincid$incid)
  hhsincid$log_incid.se <- hhsincid$se/hhsincid$incid

  return(hhsincid)
}

#' Log-likelhood for direct incidence estimate from household survey
#'
#' Calculate log-likelihood for nationally representative incidence
#' estimates from a household survey. Currently implements likelihood
#' for a log-transformed direct incidence estimate and standard error.
#' Needs to be updated to handle incidence assay outputs.
#'
#' @param mod model output, object of class `spec`.
#' @param hhsincid.dat prepared houshold survey incidence estimates (see perp
ll_hhsincid <- function(mod, hhsincid.dat){
  logincid <- log(incid(mod, fp))
  ll.incid <- sum(dnorm(hhsincid.dat$log_incid, logincid[hhsincid.dat$idx], hhsincid.dat$log_incid.se, TRUE))
  return(ll.incid)
}


###############################
####  Likelihood function  ####
###############################

prepare_likdat <- function(eppd, fp){

  ancsitedat <- eppd$ancsitedat
  
  if(exists("ancrt", fp) && fp$ancrt == "none"){
    ancsitedat <- subset(ancsitedat, type == "ancss")
    eppd$ancrtcens <- NULL
  } else if(exists("ancrt", fp) && fp$ancrt == "census")
    ancsitedat <- subset(ancsitedat, type == "ancss")
  else if(exists("ancrt", fp) && fp$ancrt == "site")
    eppd$ancrtcens <- NULL

  likdat <- list(ancsite.dat = prepare_ancsite_likdat(ancsitedat, fp),
                 hhs.dat = prepare_hhsageprev_likdat(eppd$hhs, fp))
  
  if(exists("ancrtcens", where=eppd))
    likdat$ancrtcens.dat <- prepare_ancrtcens_likdat(eppd$ancrtcens, fp)
  if(exists("hhsincid", where=eppd))
    likdat$hhsincid.dat <- prepare_hhsincid_likdat(eppd$hhsincid, fp)
  if(exists("sibmx", where=eppd))
    likdat$sibmx.dat <- prepare_sibmx_likdat(eppd$sibmx, fp)

  likdat$lastdata.idx <- max(likdat$ancsite.dat$df$yidx,
                             likdat$hhs.dat$yidx,
                             likdat$ancrtcens.dat$yidx,
                             likdat$hhsincid.dat$idx,
                             likdat$sibmx.dat$idx)
  likdat$firstdata.idx <- min(likdat$ancsite.dat$df$yidx,
                              likdat$hhs.dat$yidx,
                              likdat$ancrtcens.dat$yidx,
                              likdat$hhsage.dat$idx,
                              likdat$hhsincid.dat$idx,
                              likdat$sibmx.dat$idx)

  return(likdat)
}




lprior <- function(theta, fp){

  if(!exists("eppmod", where = fp))  # backward compatibility
    fp$eppmod <- "rspline"

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if(fp$eppmod %in% c("rspline", "logrspline", "ospline", "logospline", "rhybrid", "logrw")){
    epp_nparam <- fp$numKnots+1

    nk <- fp$numKnots

    if(fp$eppmod == "rhybrid")
      lpr <- bayes_lmvt(theta[(1+fp$rt$spline_penord):fp$rt$n_splines], tau2_prior_shape, tau2_prior_rate) +
        bayes_lmvt(theta[fp$rt$n_splines + 1:fp$rt$n_rw], rw_prior_shape, rw_prior_rate)
    else if(fp$eppmod == "logrw")
      lpr <- bayes_lmvt(theta[2:fp$numKnots], rw_prior_shape, rw_prior_rate)
    else
      lpr <- bayes_lmvt(theta[(1+fp$rtpenord):nk], tau2_prior_shape, tau2_prior_rate)

    if(exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      lpr <- lpr + dunif(theta[nk+1], r0logiotaratio.unif.prior[1], r0logiotaratio.unif.prior[2], log=TRUE)
    else
      lpr <- lpr + lprior_iota(theta[nk+1], fp)

  } else if(fp$eppmod == "rlogistic") {
    epp_nparam <- 5
    lpr <- sum(dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE))
    lpr <- lpr + lprior_iota(theta[5], fp)
  } else if(fp$eppmod == "rtrend"){ # rtrend

    epp_nparam <- 7

    lpr <- dunif(round(theta[1]), t0.unif.prior[1], t0.unif.prior[2], log=TRUE) +
      ## dunif(theta[2], t1.unif.prior[1], t1.unif.prior[2], log=TRUE) +
      dnorm(round(theta[2]), t1.pr.mean, t1.pr.sd, log=TRUE) +
      ## dunif(theta[3], logr0.unif.prior[1], logr0.unif.prior[2], log=TRUE) +
      dnorm(theta[3], logr0.pr.mean, logr0.pr.sd, log=TRUE) +
      sum(dnorm(theta[4:7], rtrend.beta.pr.mean, rtrend.beta.pr.sd, log=TRUE))
  } else if(fp$eppmod == "rlogistic_rw"){
    epp_nparam <- fp$rt$n_param+1
    lpr <- sum(dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE)) +
      bayes_lmvt(theta[4+1:fp$rt$n_rw], rw_prior_shape, rw_prior_rate)
    lpr <- lpr + lprior_iota(theta[fp$rt$n_param+1], fp)
  }

  if(fp$ancsitedata){
    lpr <- lpr + dnorm(theta[epp_nparam+1], ancbias.pr.mean, ancbias.pr.sd, log=TRUE)
    if(!exists("v.infl", where=fp)){
      anclik_nparam <- 2
      lpr <- lpr + dexp(exp(theta[epp_nparam+2]), vinfl.prior.rate, TRUE) + theta[epp_nparam+2]         # additional ANC variance
    } else
      anclik_nparam <- 1
  } else
    anclik_nparam <- 0

  paramcurr <- epp_nparam+anclik_nparam
  if(exists("ancrt", fp) && fp$ancrt %in% c("census", "both")){
    lpr <- lpr + dnorm(theta[paramcurr+1], log_frr_adjust.pr.mean, log_frr_adjust.pr.sd, log=TRUE)
    if(!exists("ancrtcens.vinfl", fp)){
      lpr <- lpr + dexp(exp(theta[paramcurr+2]), ancrtcens.vinfl.pr.rate, TRUE) + theta[paramcurr+2]
      paramcurr <- paramcurr+2
    } else
      paramcurr <- paramcurr+1
  }
  if(exists("ancrt", fp) && fp$ancrt %in% c("site", "both")){
    lpr <- lpr + dnorm(theta[paramcurr+1], ancrtsite.beta.pr.mean, ancrtsite.beta.pr.sd, log=TRUE) ## +
    ## dexp(exp(theta[np]), ancrtsite.vinfl.pr.rate, TRUE) + theta[np]
    paramcurr <- paramcurr+1
  }

  if(exists("fitincrr", where=fp)){
    incrr_nparam <- getnparam_incrr(fp)
    if(incrr_nparam){
      lpr <- lpr + lprior_incrr(theta[paramcurr+1:incrr_nparam], fp)
      paramcurr <- paramcurr+incrr_nparam
    }
  }
  
  return(lpr)
}


ll <- function(theta, fp, likdat){
  theta.last <<- theta
  fp <- update(fp, list=fnCreateParam(theta, fp))

  if(exists("fitincrr", where=fp) && fp$fitincrr==TRUE){
    ll.incpen <- sum(dnorm(diff(fp$logincrr_age, differences=2), sd=fp$sigma_agepen, log=TRUE))
  } else
    ll.incpen <- 0

  if (!exists("eppmod", where = fp) || fp$eppmod %in% c("rspline", "logrspline", "ospline", "logospline", "rhybrid"))
    if (any(is.na(fp$rvec)) || min(fp$rvec) < 0 || max(fp$rvec) > 20) 
        return(-Inf)

  mod <- simmod(fp)

  ## ANC likelihood
  if(exists("ancsite.dat", likdat))
    ll.anc <- ll_ancsite(mod, fp, coef=c(fp$ancbias, fp$ancrtsite.beta), vinfl=fp$v.infl, likdat$ancsite.dat)
  else
    ll.anc <- 0

  if(exists("ancrtcens.dat", likdat))
    ll.ancrt <- ll_ancrtcens(mod, likdat$ancrtcens.dat, fp)
  else
    ll.ancrt <- 0


  ## Household survey likelihood
  if(exists("hhs.dat", where=likdat))
    if(exists("ageprev", fp) && fp$ageprev=="binom")
      ll.hhs <- ll_hhsage_binom(mod, likdat$hhs.dat)
    else ## use probit likelihood
      ll.hhs <- ll_hhsage(mod, likdat$hhs.dat) # probit-transformed model
  else
    ll.hhs <- 0

  if(!is.null(likdat$hhsincid.dat))
    ll.incid <- ll_hhsincid(mod, likdat$hhsincid.dat)
  else
    ll.incid <- 0


  if(exists("sibmx", where=fp) && fp$sibmx){
    M.agemx <- agemx(mod)
    ll.sibmx <- ll_sibmx(M.agemx, fp$tipscoef, fp$sibmx.theta, likdat$sibmx.dat)
  } else
    ll.sibmx <- 0

  if(exists("equil.rprior", where=fp) && fp$equil.rprior){
    qM.all <- suppressWarnings(qnorm(prev(mod)))
    rvec.ann <- fp$rvec[fp$proj.steps %% 1 == 0.5]
    equil.rprior.mean <- epp:::muSS/(1-pnorm(qM.all[likdat$lastdata.idx]))
    equil.rprior.sd <- sqrt(mean((epp:::muSS/(1-pnorm(qM.all[likdat$lastdata.idx - 9:0])) - rvec.ann[likdat$lastdata.idx - 9:0])^2))  # empirical sd based on 10 previous years
    ll.rprior <- sum(dnorm(rvec.ann[(likdat$lastdata.idx+1L):length(qM.all)], equil.rprior.mean, equil.rprior.sd, log=TRUE))  # prior starts year after last data
  } else
    ll.rprior <- 0

  ## return(ll.anc+ll.hhs+ll.incpen+ll.rprior)
  return(ll.anc + ll.ancrt + ll.hhs + ll.incid + ll.sibmx + ll.rprior + ll.incpen)
}


##########################
####  IMIS functions  ####
##########################

sample.prior <- function(n, fp){

  if(!exists("eppmod", where = fp))  # backward compatibility
    fp$eppmod <- "rspline"

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  ## Calculate number of parameters
  if(fp$eppmod %in% c("rspline", "logrspline", "ospline", "logospline", "rhybrid", "logrw"))
    epp_nparam <- fp$numKnots+1L
  else if(fp$eppmod == "rlogistic")
    epp_nparam <- 5
  else if(fp$eppmod == "rtrend")
    epp_nparam <- 7
  else if(fp$eppmod == "rlogistic_rw")
    epp_nparam <- fp$rt$n_param+1

  if(fp$ancsitedata)
    if(!exists("v.infl", fp))
      anclik_nparam <- 2
    else
      anclik_nparam <- 1
  else
    anclik_nparam <- 0

  if(exists("ancrt", fp) && fp$ancrt == "both")
    ancrt_nparam <- 2
  else if(exists("ancrt", fp) && fp$ancrt == "census")
    ancrt_nparam <- 1
  else if(exists("ancrt", fp) && fp$ancrt == "site")
    ancrt_nparam <- 1
  else
    ancrt_nparam <- 0

  if(exists("ancrt", fp) && fp$ancrt %in% c("census", "both") && !exists("ancrtcens.vinfl", fp))
    ancrt_nparam <- ancrt_nparam+1

  nparam <- epp_nparam+anclik_nparam+ancrt_nparam

  if(exists("fitincrr", where=fp)) nparam <- nparam+getnparam_incrr(fp)

  
  ## Create matrix for storing samples
  mat <- matrix(NA, n, nparam)

  if(fp$eppmod %in% c("rspline", "logrspline", "ospline", "logospline", "rhybrid", "logrw")){
    epp_nparam <- fp$numKnots+1

    if(fp$eppmod == "rspline")
      mat[,1] <- rnorm(n, 1.5, 1)                                                   # u[1]
    if(fp$eppmod == "ospline")
      mat[,1] <- rnorm(n, 0.5, 1)
    else # logrspline, logospline, logrw
      mat[,1] <- rnorm(n, 0.2, 1)                                                   # u[1]
    if(fp$eppmod == "rhybrid"){
      mat[,2:fp$rt$n_splines] <- bayes_rmvt(n, fp$rt$n_splines-1,tau2_init_shape, tau2_init_rate)
      mat[,fp$rt$n_splines+1:fp$rt$n_rw] <- bayes_rmvt(n, fp$rt$n_rw, rw_prior_shape, rw_prior_rate)  # u[2:numKnots]
    } else if(fp$eppmod == "logrw"){
      mat[,2:fp$rt$n_rw] <- bayes_rmvt(n, fp$rt$n_rw-1, rw_prior_shape, rw_prior_rate)  # u[2:numKnots]
    } else {
      mat[,2:fp$numKnots] <- bayes_rmvt(n, fp$numKnots-1,tau2_init_shape, tau2_init_rate)  # u[2:numKnots]
    }

    if(exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      mat[,fp$numKnots+1] <-  runif(n, r0logiotaratio.unif.prior[1], r0logiotaratio.unif.prior[2])  # ratio r0 / log(iota)
    else
      mat[,fp$numKnots+1] <- sample_iota(n, fp)
  } else if(fp$eppmod == "rlogistic"){
    mat[,1:4] <- t(matrix(rnorm(4*n, rlog_pr_mean, rlog_pr_sd), 4))
    mat[,5] <- sample_iota(n, fp)
  } else if(fp$eppmod == "rtrend"){ # r-trend

    mat[,1] <- runif(n, t0.unif.prior[1], t0.unif.prior[2])           # t0
    ## mat[,2] <- runif(n, t1.unif.prior[1], t1.unif.prior[2])        # t1
    mat[,2] <- rnorm(n, t1.pr.mean, t1.pr.sd)
    ## mat[,3] <- runif(n, logr0.unif.prior[1], logr0.unif.prior[2])  # r0
    mat[,3] <- rnorm(n, logr0.pr.mean, logr0.pr.sd)  # r0
    mat[,4:7] <- t(matrix(rnorm(4*n, rtrend.beta.pr.mean, rtrend.beta.pr.sd), 4, n))  # beta
  } else if(fp$eppmod == "rlogistic_rw") {
    mat[,1:4] <- t(matrix(rnorm(4*n, rlog_pr_mean, rlog_pr_sd), 4))
    mat[,4+1:fp$rt$n_rw] <- bayes_rmvt(n, fp$rt$n_rw, rw_prior_shape, rw_prior_rate)  # u[2:numKnots]
    mat[,fp$rt$n_param+1] <- sample_iota(n, fp)
  }

  ## sample ANC bias paramters
  if(fp$ancsitedata){
    mat[,epp_nparam+1] <- rnorm(n, ancbias.pr.mean, ancbias.pr.sd)   # ancbias parameter
    if(!exists("v.infl", where=fp))
      mat[,epp_nparam+2] <- log(rexp(n, vinfl.prior.rate))
  }

  ## sample ANCRT parameters
  paramcurr <- epp_nparam+anclik_nparam
  if(exists("ancrt", where=fp) && fp$ancrt %in% c("census", "both")){
    mat[,paramcurr+1] <- rnorm(n, log_frr_adjust.pr.mean, log_frr_adjust.pr.sd)
    if(!exists("ancrtcens.vinfl", fp)){
      mat[,paramcurr+2] <- log(rexp(n, ancrtcens.vinfl.pr.rate))
      paramcurr <- paramcurr+2
    } else
      paramcurr <- paramcurr+1
  }
  if(exists("ancrt", where=fp) && fp$ancrt %in% c("site", "both")){
    mat[,paramcurr+1] <- rnorm(n, ancrtsite.beta.pr.mean, ancrtsite.beta.pr.sd)
    ## mat[,nparam] <- log(rexp(n, ancrtsite.vinfl.pr.rate))
    paramcurr <- paramcurr+1
  }

  if(exists("fitincrr", where=fp)){
    incrr_nparam <- getnparam_incrr(fp)
    if(incrr_nparam)
      mat[,paramcurr+1:incrr_nparam] <- sample_incrr(n, fp)
    paramcurr <- paramcurr+incrr_nparam
  }

  return(mat)
}

ldsamp <- function(theta, fp){

  if(!exists("eppmod", where = fp))  # backward compatibility
    fp$eppmod <- "rspline"

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if(fp$eppmod %in% c("rspline", "logrspline", "ospline", "logospline", "rhybrid", "logrw")){
    epp_nparam <- fp$numKnots+1

    nk <- fp$numKnots

    if(fp$eppmod == "rspline")  # u[1]
      lpr <- dnorm(theta[1], 1.5, 1, log=TRUE)
    if(fp$eppmod == "ospline")
      lpr <- dnorm(theta[1], 0.5, 1, log=TRUE)
    else # logrspline, logospline, logrw
      lpr <- dnorm(theta[1], 0.2, 1, log=TRUE)

    if(fp$eppmod == "rhybrid")
      lpr <- bayes_lmvt(theta[2:fp$rt$n_splines], tau2_init_shape, tau2_init_rate) +
        bayes_lmvt(theta[fp$rt$n_splines + 1:fp$rt$n_rw], rw_prior_shape, rw_prior_rate)
    else if(fp$eppmod == "logrw")
      bayes_lmvt(theta[2:fp$rt$n_rw], rw_prior_shape, rw_prior_rate)
    else
      lpr <- bayes_lmvt(theta[2:nk], tau2_prior_shape, tau2_prior_rate)


    if(exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      lpr <- lpr + dunif(theta[nk+1], r0logiotaratio.unif.prior[1], r0logiotaratio.unif.prior[2], log=TRUE)
    else
      lpr <- lpr + ldsamp_iota(theta[nk+1], fp)

  } else if(fp$eppmod == "rlogistic") {
    epp_nparam <- 5
    lpr <- sum(dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE))
    lpr <- lpr + ldsamp_iota(theta[5], fp)
  } else if(fp$eppmod == "rtrend"){ # rtrend

    epp_nparam <- 7

    lpr <- dunif(round(theta[1]), t0.unif.prior[1], t0.unif.prior[2], log=TRUE) +
      ## dunif(theta[2], t1.unif.prior[1], t1.unif.prior[2], log=TRUE) +
      dnorm(round(theta[2]), t1.pr.mean, t1.pr.sd, log=TRUE) +
      ## dunif(theta[3], logr0.unif.prior[1], logr0.unif.prior[2], log=TRUE) +
      dnorm(theta[3], logr0.pr.mean, logr0.pr.sd, log=TRUE) +
      sum(dnorm(theta[4:7], rtrend.beta.pr.mean, rtrend.beta.pr.sd, log=TRUE))
  } else if(fp$eppmod == "rlogistic_rw"){
    epp_nparam <- fp$rt$n_param+1
    lpr <- sum(dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE)) +
      bayes_lmvt(theta[4+1:fp$rt$n_rw], rw_prior_shape, rw_prior_rate)
    lpr <- lpr + ldsamp_iota(theta[fp$rt$n_param+1], fp)
  }

  if(fp$ancsitedata){
    lpr <- lpr + dnorm(theta[epp_nparam+1], ancbias.pr.mean, ancbias.pr.sd, log=TRUE)
    if(!exists("v.infl", where=fp)){
      anclik_nparam <- 2
      lpr <- lpr + dexp(exp(theta[epp_nparam+2]), vinfl.prior.rate, TRUE) + theta[epp_nparam+2]         # additional ANC variance
    } else
      anclik_nparam <- 1
  } else
    anclik_nparam <- 0

  paramcurr <- epp_nparam+anclik_nparam
  if(exists("ancrt", fp) && fp$ancrt %in% c("census", "both")){
    lpr <- lpr + dnorm(theta[paramcurr+1], log_frr_adjust.pr.mean, log_frr_adjust.pr.sd, log=TRUE)
    if(!exists("ancrtcens.vinfl", fp)){
      lpr <- lpr + dexp(exp(theta[paramcurr+2]), ancrtcens.vinfl.pr.rate, TRUE) + theta[paramcurr+2]
      paramcurr <- paramcurr+2
    } else
      paramcurr <- paramcurr+1
  }
  if(exists("ancrt", fp) && fp$ancrt %in% c("site", "both")){
    lpr <- lpr + dnorm(theta[paramcurr+1], ancrtsite.beta.pr.mean, ancrtsite.beta.pr.sd, log=TRUE) ## +
    ## dexp(exp(theta[np]), ancrtsite.vinfl.pr.rate, TRUE) + theta[np]
    paramcurr <- paramcurr+1
  }

  if(exists("fitincrr", where=fp)){
    incrr_nparam <- getnparam_incrr(fp)
    if(incrr_nparam){
      lpr <- lpr + ldsamp_incrr(theta[paramcurr+1:incrr_nparam], fp)
      paramcurr <- paramcurr+incrr_nparam
    }
  }
  
  return(lpr)
}


prior <- function(theta, fp, log=FALSE){
  if(is.vector(theta))
    lval <- lprior(theta, fp)
  else
    lval <- unlist(lapply(seq_len(nrow(theta)), function(i) (lprior(theta[i,], fp))))
  if(log)
    return(lval)
  else
    return(exp(lval))
}

likelihood <- function(theta, fp, likdat, log=FALSE){
  if(is.vector(theta))
    lval <- ll(theta, fp, likdat)
  else
    lval <- unlist(lapply(seq_len(nrow(theta)), function(i) ll(theta[i,], fp, likdat)))
  if(log)
    return(lval)
  else
    return(exp(lval))
}

dsamp <- function(theta, fp, log=FALSE){
  if(is.vector(theta))
    lval <- ldsamp(theta, fp)
  else
    lval <- unlist(lapply(seq_len(nrow(theta)), function(i) (ldsamp(theta[i,], fp))))
  if(log)
    return(lval)
  else
    return(exp(lval))
}
