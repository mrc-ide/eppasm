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
#' @param eppd EPP data object
#' @param anchor.year year in which EPP data inputs start
#' NOTE: requires year to be stored in column names of anc.prev
prepare_ancsite_likdat <- function(eppd, anchor.year=1970L){

  anc.prev <- eppd$anc.prev
  anc.n <- eppd$anc.n
  anc.used <- eppd$anc.used

  anc.prev <- anc.prev[anc.used,,drop=FALSE]  # keep only used sites
  anc.n <- anc.n[anc.used,,drop=FALSE]        # keep only used sites

  ancobs.idx <- mapply(intersect, lapply(as.data.frame(t(!is.na(anc.prev))), which),
                       lapply(as.data.frame(t(!is.na(anc.n))), which), SIMPLIFY=FALSE)
  ## limit to years with both prevalence and N observations (likely input errors in EPP if not)

  nobs <- sapply(ancobs.idx, length)

  anc.years.lst <- lapply(ancobs.idx, function(i) as.integer(colnames(anc.prev)[i]))
  anc.prev.lst <- setNames(lapply(seq_along(ancobs.idx), function(i) as.numeric(anc.prev[i, ancobs.idx[[i]]])), rownames(anc.prev))
  anc.n.lst <- setNames(lapply(seq_along(ancobs.idx), function(i) as.numeric(anc.n[i, ancobs.idx[[i]]])), rownames(anc.n))

  X.lst <- mapply(cbind, Intercept=lapply(nobs, rep, x=1), ancrt=lapply(nobs, rep, x=0), SIMPLIFY=FALSE)

  if(exists("ancrtsite.prev", where=eppd) && !is.null(eppd$ancrtsite.prev)){
    ancrtsite.prev <- eppd$ancrtsite.prev
    ancrtsite.n <- eppd$ancrtsite.n

    ancrtsite.prev <- ancrtsite.prev[anc.used,,drop=FALSE]  # keep only used sites
    ancrtsite.n <- ancrtsite.n[anc.used,,drop=FALSE]        # keep only used sites

    ancrtsiteobs.idx <- mapply(intersect, lapply(as.data.frame(t(!is.na(ancrtsite.prev))), which),
                               lapply(as.data.frame(t(!is.na(ancrtsite.n))), which), SIMPLIFY=FALSE)
    ## limit to years with both prevalence and N observations (likely input errors in EPP if not)

    nobs <- sapply(ancrtsiteobs.idx, length)

    ancrtsite.years.lst <- lapply(ancrtsiteobs.idx, function(i) as.integer(colnames(ancrtsite.prev)[i]))
    ancrtsite.prev.lst <- setNames(lapply(seq_along(ancrtsiteobs.idx), function(i) as.numeric(ancrtsite.prev[i, ancrtsiteobs.idx[[i]]])), rownames(ancrtsite.prev))
    ancrtsite.n.lst <- setNames(lapply(seq_along(ancrtsiteobs.idx), function(i) as.numeric(ancrtsite.n[i, ancrtsiteobs.idx[[i]]])), rownames(ancrtsite.n))

    ancrtsite.X.lst <- mapply(cbind, Intercept=lapply(nobs, rep, x=1), ancrt=lapply(nobs, rep, x=1), SIMPLIFY=FALSE)

    ## Combine SS and RT data
    anc.years.lst <- mapply(c, anc.years.lst, ancrtsite.years.lst, SIMPLIFY=FALSE)
    anc.prev.lst <- mapply(c, anc.prev.lst, ancrtsite.prev.lst, SIMPLIFY=FALSE)
    anc.n.lst <- mapply(c, anc.n.lst, ancrtsite.n.lst, SIMPLIFY=FALSE)
    X.lst <- mapply(rbind, X.lst, ancrtsite.X.lst, SIMPLIFY=FALSE)
  }

  ## eliminate records with no observations
  anc.years.lst <- anc.years.lst[sapply(anc.years.lst, length) > 0]
  anc.prev.lst <- anc.prev.lst[sapply(anc.years.lst, length) > 0]
  anc.n.lst <- anc.n.lst[sapply(anc.years.lst, length) > 0]
  X.lst <- X.lst[sapply(anc.years.lst, length) > 0]

  x.lst <- mapply(function(p, n) (p*n+0.5)/(n+1), anc.prev.lst, anc.n.lst, SIMPLIFY=FALSE)
  W.lst <- lapply(x.lst, qnorm)
  v.lst <- mapply(function(W, x, n) 2*pi*exp(W^2)*x*(1-x)/n, W.lst, x.lst, anc.n.lst, SIMPLIFY=FALSE)
  anc.idx.lst <- lapply(anc.years.lst, "-", anchor.year-1)  ## index of observations relative to output prevalence vector


  anclik.dat <- list(W.lst = W.lst,
                     v.lst = v.lst,
                     n.lst = anc.n.lst,
                     X.lst = X.lst,
                     anc.idx.lst = anc.idx.lst)

  return(anclik.dat)
}

ll_anc <- function(qM, coef=c(0, 0), vinfl=0, anclik.dat){

  ## linear model offset
  mu <- lapply(lapply(anclik.dat$X.lst, "%*%", coef), c)

  d.lst <- mapply(function(w, mu, idx) w - (qM[idx]+mu), anclik.dat$W.lst, mu, anclik.dat$anc.idx.lst, SIMPLIFY=FALSE)
  v.lst <- lapply(anclik.dat$v.lst, "+", vinfl)

  return(log(anclik::anc_resid_lik(d.lst, v.lst)))
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

prepare_ancrtcens_likdat <- function(dat, anchor.year){

  x.ancrt <- (dat$prev*dat$n+0.5)/(dat$n+1)
  dat$W.ancrt <- qnorm(x.ancrt)
  dat$v.ancrt <- 2*pi*exp(dat$W.ancrt^2)*x.ancrt*(1-x.ancrt)/dat$n
  dat$idx <- dat$year - anchor.year+1

  return(dat)
}

ll_ancrtcens <- function(qM.preg, ancrtcens.dat, fp){
  sum(dnorm(ancrtcens.dat$W.ancrt, qM.preg[ancrtcens.dat$idx], sqrt(ancrtcens.dat$v.ancrt + fp$ancrtcens.vinfl), log=TRUE))
}



###################################
####  Age/sex incidence model  ####
###################################
sexincrr.pr.mean <- log(1.38)
sexincrr.pr.sd <- 0.2

mf_transm_rr.pr.mean <- log(1.9)
mf_transm_rr.pr.sd <- 0.4

NPARAM_RW2 <- 13
ageincrr.pr.mean <- c(-1.40707274, -0.23518703, 0.69314718, 0.78845736, -0.39975544, -0.70620810, -0.84054571, -0.02101324, -0.16382449, -0.37914407, -0.59639985, -0.82038300)
ageincrr.pr.sd <- 0.5

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
    param$frr_art <- fp$frr_art
    param$frr_art[1:2,,,] <- param$frr_art[1:2,,,] * exp(param$log_frr_adjust)

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
    if(exists("fitincrr", where=fp) && fp$fitincrr %in% c(TRUE, "lognorm", "relbehav")){

      if(fp$incidmod == "eppspectrum"){
        param$incrr_sex <- fp$incrr_sex
        param$incrr_sex[] <- exp(theta[paramcurr+1])
      } else if(fp$incidmod == "transm") {
        param$mf_transm_rr <- rep(exp(theta[paramcurr+1]), fp$ss$PROJ_YEARS)
      }

      if(fp$fitincrr==TRUE){
        incrr_nparam <- NPARAM_RW2
        theta_incrr <- theta[paramcurr+1:incrr_nparam]
        paramcurr <- paramcurr+incrr_nparam

        ## param$sigma_agepen <- exp(theta_incrr[incrr_nparam])
        param$sigma_agepen <- 0.4

        param$logincrr_age <- array(0, c(7, 2))
        param$logincrr_age[-3,] <- theta_incrr[2:13]

        param$incrr_age <- fp$incrr_age
        if(exists("linincrr", where=fp) && fp$linincrr)
          param$incrr_age[fp$ss$p.age15to49.idx,,] <- exp(apply(param$logincrr_age, 2, function(x) approx(3:9*5, x, 15:49, rule=2)$y))
        else
          param$incrr_age[fp$ss$p.age15to49.idx,,] <- apply(exp(param$logincrr_age), 2, rep, each=5)
        param$incrr_age[36:66,,] <- sweep(fp$incrr_age[36:66,,fp$ss$PROJ_YEARS], 2,
                                          param$incrr_age[35,,fp$ss$PROJ_YEARS]/fp$incrr_age[35,,fp$ss$PROJ_YEARS], "*")
      } else if(fp$fitincrr=="lognorm"){
        incrr_nparam <- 7
        theta_incrr <- theta[paramcurr+1:incrr_nparam]
        paramcurr <- paramcurr+incrr_nparam

        param$logincrr_age <- cbind(calc_lognorm_logagerr(theta_incrr[2:4]),
                                    calc_lognorm_logagerr(theta_incrr[5:7]))
        param$incrr_age <- fp$incrr_age
        param$incrr_age[,,] <- apply(exp(param$logincrr_age), 2, rep, c(rep(5, 13), 1))

      } else if(fp$fitincrr == "relbehav"){
        incrr_nparam <- NPAR_RELBEHAV
        par <- theta[paramcurr+2:incrr_nparam]
        paramcurr <- paramcurr+incrr_nparam

        param$adjustpar <- par
        logadjust1 <- cbind(approx(c(17, 27, 38, 49), c(par[1], 0, cumsum(par[2:3])), xout=15:80, rule=2)$y,
                            approx(c(17, 27, 38, 49), c(par[4], 0, cumsum(par[5:6])), xout=15:80, rule=2)$y)

        logadjust2 <- cbind(approx(c(17, 27, 38, 49), c(par[1]+par[7], 0, cumsum(par[2:3])), xout=15:80, rule=2)$y,
                            approx(c(17, 27, 38, 49), c(par[4]+par[8], 0, cumsum(par[5:6])), xout=15:80, rule=2)$y)

        BREAK_YEAR <- 36
        param$incrr_age <- fp$logrelbehav
        param$incrr_age[,,1:(BREAK_YEAR-1)] <- exp(sweep(fp$logrelbehav[,,1:(BREAK_YEAR-1)], 1:2, logadjust1, "+"))
        param$incrr_age[,,BREAK_YEAR:fp$SIM_YEARS] <- exp(sweep(fp$logrelbehav[,,BREAK_YEAR:fp$SIM_YEARS], 1:2, logadjust2, "+"))
      }
      
    } else
      incrr_nparam <- 0

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

  hhsage$sidx <- match(hhsage$sex, c("male", "female"))
  hhsage$aidx <- as.integer(substr(hhsage$agegr, 1, 2)) - fp$ss$AGE_START+1L
  hhsage$yidx <- hhsage$year - (anchor.year - 1)

  hhsage$arridx <- hhsage$aidx + (hhsage$sidx-1)*fp$ss$pAG + (hhsage$yidx-1)*fp$ss$NG*fp$ss$pAG

  return(subset(hhsage, aidx > 0))
}


#' Log likelihood for age 15-49 household survey prevalence
ll_hhs <- function(qM, hhslik.dat){
  return(sum(dnorm(hhslik.dat$W.hhs, qM[hhslik.dat$idx], hhslik.dat$sd.W.hhs, log=TRUE)))
}

#' Log likelihood for age-specific household survey prevalence
ll_hhsage <- function(mod, hhsage.dat){
  qM.age <- suppressWarnings(qnorm(ageprev(mod, arridx=hhsage.dat$arridx, agspan=5)))
  if(any(is.na(qM.age))) return(-Inf)
  sum(dnorm(hhsage.dat$W.hhs, qM.age, hhsage.dat$sd.W.hhs, log=TRUE))
}


#' Log likelihood for age-specific household survey prevalence using binomial approximation
ll_hhsage_binom <- function(mod, hhsage.dat){
  prevM.age <- suppressWarnings(ageprev(mod, arridx=hhsage.dat$arridx, agspan=5))
  if(any(is.na(prevM.age)) || any(prevM.age >= 1)) return(-Inf)
  ll <- sum(ldbinom(hhsage.dat$x_eff, hhsage.dat$n_eff, prevM.age))
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

  anchor_year <- floor(fp$proj.steps[1])

  likdat <- list(anclik.dat = prepare_ancsite_likdat(eppd, anchor.year=anchor_year),
                 hhslik.dat = epp::fnPrepareHHSLikData(eppd$hhs, anchor.year=anchor_year))
  if(exists("ancrtcens", where=eppd))
    likdat$ancrtcens.dat <- prepare_ancrtcens_likdat(eppd$ancrtcens, anchor.year=anchor_year)
  if(exists("hhsage", where=eppd))
    likdat$hhsage.dat <- prepare_hhsageprev_likdat(eppd$hhsage, fp)
  if(exists("hhsincid", where=eppd))
    likdat$hhsincid.dat <- prepare_hhsincid_likdat(eppd$hhsincid, fp)
  if(exists("sibmx", where=eppd))
    likdat$sibmx.dat <- prepare_sibmx_likdat(eppd$sibmx, fp)

  likdat$lastdata.idx <- max(unlist(likdat$anclik.dat$anc.idx.lst),
                             likdat$hhslik.dat$idx,
                             likdat$ancrtcens.dat$idx,
                             likdat$hhsage.dat$idx,
                             likdat$hhsincid.dat$idx,
                             likdat$sibmx.dat$idx)
  likdat$firstdata.idx <- min(unlist(likdat$anclik.dat$anc.idx.lst),
                              likdat$hhslik.dat$idx,
                              likdat$ancrtcens.dat$idx,
                              likdat$ancrtcens.dat$idx,
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

  if(exists("fitincrr", where=fp) && fp$fitincrr %in% c(TRUE, "lognorm", "relbehav")){

    if(fp$incidmod == "eppspectrum")
      lpr <- lpr + dnorm(theta[paramcurr+1], sexincrr.pr.mean, sexincrr.pr.sd, log=TRUE)
    else if(fp$incidmod == "transm")
      lpr <- lpr + dnorm(theta[paramcurr+1], mf_transm_rr.pr.mean, mf_transm_rr.pr.sd, log=TRUE)

    if(fp$fitincrr == TRUE){
      incrr_nparam <- NPARAM_RW2
      theta_incrr <- theta[paramcurr+1:incrr_nparam]
      paramcurr <- paramcurr+incrr_nparam

      lpr <- lpr +
        sum(dnorm(theta_incrr[2:13], ageincrr.pr.mean, ageincrr.pr.sd, log=TRUE))
        ## dnorm(theta_incrr[14], -1, 0.7, log=TRUE)
    } else if(fp$fitincrr=="lognorm"){
      incrr_nparam <- 7
      theta_incrr <- theta[paramcurr+1:incrr_nparam]
      paramcurr <- paramcurr+incrr_nparam

      lpr <- lpr +
        sum(dnorm(theta_incrr[c(2,5)], lognorm.a0.pr.mean, lognorm.a0.pr.sd, log=TRUE)) +
        sum(dnorm(theta_incrr[c(3,6)], lognorm.meanlog.pr.mean, lognorm.meanlog.pr.sd, log=TRUE)) +
        sum(dnorm(theta_incrr[c(4,7)], lognorm.logsdlog.pr.mean, lognorm.logsdlog.pr.sd, log=TRUE))
    } else if(fp$fitincrr=="relbehav"){
      incrr_nparam <- NPAR_RELBEHAV
      par <- theta[paramcurr+2:incrr_nparam]
      paramcurr <- paramcurr+incrr_nparam

      lpr <- lpr + sum(dnorm(par, 0, relbehav_adjust_sd, log=TRUE));
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

  qM.all <- suppressWarnings(qnorm(prev(mod)))
  qM.preg <- if(exists("pregprev", where=fp) && !fp$pregprev) qM.all else suppressWarnings(qnorm(fnPregPrev(mod, fp)))

  if(any(is.na(qM.preg[likdat$firstdata.idx:likdat$lastdata.idx])) ||
     any(is.na(qM.all[likdat$firstdata.idx:likdat$lastdata.idx])) ||
     any(qM.preg[likdat$firstdata.idx:likdat$lastdata.idx] == -Inf) ||
     any(qM.preg[likdat$firstdata.idx:likdat$lastdata.idx] > 2)) # prevalence not greater than pnorm(2) = 0.977
    return(-Inf)

  ## ANC likelihood
  if(fp$ancsitedata)
    ll.anc <- ll_anc(qM.preg, coef=c(fp$ancbias, fp$ancrtsite.beta), vinfl=fp$v.infl, likdat$anclik.dat)
  else
    ll.anc <- 0

  if(exists("ancrt", fp) && fp$ancrt %in% c("census", "both"))
    ll.ancrt <- ll_ancrtcens(qM.preg, likdat$ancrtcens.dat, fp)
  else
    ll.ancrt <- 0


  ## Household survey likelihood
  if(exists("ageprev", where=fp) && fp$ageprev=="binom")
    ll.hhs <- ll_hhsage_binom(mod, likdat$hhsage.dat)
  else if(exists("ageprev", where=fp) && (fp$ageprev==TRUE | fp$ageprev == "probit")) # ==TRUE for backward compatibility
    ll.hhs <- ll_hhsage(mod, likdat$hhsage.dat) # probit-transformed model
  else
    ll.hhs <- ll_hhs(qM.all, likdat$hhslik.dat)

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

  if(exists("fitincrr", where=fp) && fp$fitincrr==TRUE) nparam <- nparam+NPARAM_RW2
  if(exists("fitincrr", where=fp) && fp$fitincrr=="lognorm") nparam <- nparam+7
  if(exists("fitincrr", where=fp) && fp$fitincrr=="relbehav") nparam <- nparam+NPAR_RELBEHAV

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

  if(exists("fitincrr", where=fp) && fp$fitincrr %in% c(TRUE, "lognorm", "relbehav")){
    if(fp$incidmod == "eppspectrum")
      mat[,paramcurr+1] <- rnorm(n, sexincrr.pr.mean, sexincrr.pr.sd)
    else if(fp$incidmod == "transm")
      mat[,paramcurr+1] <- rnorm(n, mf_transm_rr.pr.mean, mf_transm_rr.pr.sd)

    if(fp$fitincrr == TRUE){
      incrr_nparam <- NPARAM_RW2
      mat[,paramcurr+2:13] <- t(matrix(rnorm(n*12, ageincrr.pr.mean, ageincrr.pr.sd), nrow=12))
      ## mat[,paramcurr+14] <- rnorm(n, -1, 0.7)  # log variance of ageincrr difference penalty
    } else if(fp$fitincrr=="lognorm"){
      incrr_nparam <- 7
      
      mat[,paramcurr+c(2,5)] <- t(matrix(rnorm(n*2, lognorm.a0.pr.mean, lognorm.a0.pr.sd), nrow=2))
      mat[,paramcurr+c(3,6)] <- t(matrix(rnorm(n*2, lognorm.meanlog.pr.mean, lognorm.meanlog.pr.sd), nrow=2))
      mat[,paramcurr+c(4,7)] <- t(matrix(rnorm(n*2, lognorm.logsdlog.pr.mean, lognorm.logsdlog.pr.sd), nrow=2))
    } else if(fp$fitincrr=="relbehav"){
      incrr_nparam <- NPAR_RELBEHAV
      mat[,paramcurr+2:NPAR_RELBEHAV] <- rnorm(n*(NPAR_RELBEHAV-1), 0, relbehav_adjust_sd)
    }
  } else
    incrr_nparam <- 0
  paramcurr <- paramcurr+incrr_nparam

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


  if(exists("fitincrr", where=fp) && fp$fitincrr %in% c(TRUE, "lognorm", "relbehav")){

    if(fp$incidmod == "eppspectrum")
      lpr <- lpr + dnorm(theta[paramcurr+1], sexincrr.pr.mean, sexincrr.pr.sd, log=TRUE)
    else if(fp$incidmod == "transm")
      lpr <- lpr + dnorm(theta[paramcurr+1], mf_transm_rr.pr.mean, mf_transm_rr.pr.sd, log=TRUE)

    if(fp$fitincrr == TRUE){
      incrr_nparam <- NPARAM_RW2
      theta_incrr <- theta[paramcurr+1:incrr_nparam]
      paramcurr <- paramcurr+incrr_nparam

      lpr <- lpr +
        sum(dnorm(theta_incrr[2:13], ageincrr.pr.mean, ageincrr.pr.sd, log=TRUE))
        ## dnorm(theta_incrr[14], -1, 0.7, log=TRUE)
    } else if(fp$fitincrr=="lognorm"){
      incrr_nparam <- 7
      theta_incrr <- theta[paramcurr+1:incrr_nparam]
      paramcurr <- paramcurr+incrr_nparam

      lpr <- lpr +
        sum(dnorm(theta_incrr[c(2,5)], lognorm.a0.pr.mean, lognorm.a0.pr.sd, log=TRUE)) +
        sum(dnorm(theta_incrr[c(3,6)], lognorm.meanlog.pr.mean, lognorm.meanlog.pr.sd, log=TRUE)) +
        sum(dnorm(theta_incrr[c(4,7)], lognorm.logsdlog.pr.mean, lognorm.logsdlog.pr.sd, log=TRUE))
    } else if(fp$fitincrr=="relbehav"){
      incrr_nparam <- NPAR_RELBEHAV
      par <- theta[paramcurr+2:incrr_nparam]
      paramcurr <- paramcurr+incrr_nparam

      lpr <- lpr + sum(dnorm(par, 0, relbehav_adjust_sd, log=TRUE));
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
