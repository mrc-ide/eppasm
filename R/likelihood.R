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
rw_prior_sd <- 0.06


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
ancrtsite.beta.pr.sd <- 1.0
## ancrtsite.beta.pr.sd <- 0.05
## ancrtsite.vinfl.pr.rate <- 1/0.015


#' Prepare design matrix indices for ANC prevalence predictions
#'
#' @param ancsite_df data.frame of site-level ANC design for predictions
#' @param fp fixed parameter input list
#'
#' @export
#' @examples
#' pjnz <- system.file("extdata/testpjnz", "Botswana2017.PJNZ", package="eppasm")
#' bw <- prepare_spec_fit(pjnz, proj.end=2021.5)
#'
#'
#' bw_u_ancsite <- attr(bw$Urban, "eppd")$ancsitedat
#' fp <- attr(bw$Urban, "specfp")
#'
#' ancsite_pred_df(bw_u_ancsite, fp)
#'
ancsite_pred_df <- function(ancsite_df, fp) {

  df <- ancsite_df
  anchor.year <- fp$ss$proj_start

  df$aidx <- df$age - fp$ss$AGE_START + 1L
  df$yidx <- df$year - anchor.year + 1

  ## List of all unique agegroup / year combinations for which prevalence is needed
  datgrp <- unique(df[c("aidx", "yidx", "agspan")])
  datgrp$qMidx <- seq_len(nrow(datgrp))

  ## Indices for accessing prevalence offset from datgrp
  df <- merge(df, datgrp)

  list(df = df, datgrp = datgrp)
}


#' Prepare site-level ANC prevalence data for EPP random-effects likelihood
#'
#' @param ancsitedat data.frame of site-level ANC data
#' @param fp fixed parameter input list, including state space

prepare_ancsite_likdat <- function(ancsitedat, fp){

  d <- ancsite_pred_df(ancsitedat, fp)

  df <- d$df[c("site", "year", "used", "type", "age", "agspan",
               "n", "prev", "aidx", "yidx", "qMidx")]

  ## Calculate probit transformed prevalence and variance approximation
  df$pstar <- (df$prev * df$n + 0.5) / (df$n + 1)
  df$W <- stats::qnorm(df$pstar)
  df$v <- 2 * pi * exp(df$W^2) * df$pstar * (1 - df$pstar) / df$n

  ## Design matrix for fixed effects portion
  df$type <- factor(df$type, c("ancss", "ancrt"))
  Xancsite <- stats::model.matrix(~type, df)

  ## Indices for observation
  df_idx.lst <- split(seq_len(nrow(df)), factor(df$site))

  list(df = df,
       datgrp = d$datgrp,
       Xancsite = Xancsite,
       df_idx.lst = df_idx.lst)
}

#' @export
ll_ancsite <- function(mod, fp, coef=c(0, 0), vinfl=0, dat){

  df <- dat$df

  if(!nrow(df))
    return(0)

  pregprevM <- agepregprev(mod, fp, dat$datgrp$aidx, dat$datgrp$yidx, dat$datgrp$agspan)

  ## If calendar year projection, average current and previous year prevalence to
  ## approximate mid-year prevalence.
  ## NOTE: This will be inefficient if values for every year because it duplicates
  ##   prevalence calculation for the same year.

  if (fp$projection_period == "calendar") {
    pregprevM_last <- agepregprev(mod, fp, dat$datgrp$aidx, dat$datgrp$yidx-1L, dat$datgrp$agspan)
    pregprevM <- 0.5 * (pregprevM + pregprevM_last)
  }

  qM <- suppressWarnings(stats::qnorm(pregprevM))

  if(any(is.na(qM)) || any(qM == -Inf) || any(qM > 2))  ## prev < 0.977
    return(-Inf)

  mu <- qM[df$qMidx] + dat$Xancsite %*% coef
  d <- df$W - mu
  v <- df$v + vinfl

  d.lst <- lapply(dat$df_idx.lst, function(idx) d[idx])
  v.lst <- lapply(dat$df_idx.lst, function(idx) v[idx])

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
log_frr_adjust.pr.sd <- 1.0
## log_frr_adjust.pr.sd <- 0.2
ancrtcens.vinfl.pr.rate <- 1/0.015

prepare_ancrtcens_likdat <- function(dat, fp){

  anchor.year <- fp$ss$proj_start

  x.ancrt <- (dat$prev*dat$n+0.5)/(dat$n+1)
  dat$W.ancrt <- stats::qnorm(x.ancrt)
  dat$v.ancrt <- 2*pi*exp(dat$W.ancrt^2)*x.ancrt*(1-x.ancrt)/dat$n

  if(!exists("age", dat))
    dat$age <- rep(15, nrow(dat))

  if(!exists("agspan", dat))
    dat$agspan <- rep(35, nrow(dat))

  dat$aidx <- dat$age - fp$ss$AGE_START + 1
  dat$yidx <- dat$year - anchor.year + 1

  return(dat)
}

#' @export
ll_ancrtcens <- function(mod, dat, fp, pointwise = FALSE){
  if(!nrow(dat))
    return(0)

  qM.prev <- agepregprev(mod, fp, dat$aidx, dat$yidx, dat$agspan)

  ## If calendar year projection, average current and previous year prevalence to
  ## approximate mid-year prevalence.
  ## NOTE: This will be inefficient if values for every year because it duplicates
  ##   prevalence calculation for the same year.

  if (fp$projection_period == "calendar") {
    qM.prev_last <- qM.prev <- agepregprev(mod, fp, dat$aidx, dat$yidx-1L, dat$agspan)
    qM.prev <- 0.5 * (qM.prev + qM.prev_last)
  }

  qM.prev <- suppressWarnings(stats::qnorm(qM.prev))

  if(any(is.na(qM.prev)))
    val <- rep(-Inf, nrow(dat))
  else
    val <- stats::dnorm(dat$W.ancrt, qM.prev, sqrt(dat$v.ancrt + fp$ancrtcens.vinfl), log=TRUE)

  if(pointwise)
    return(val)

  sum(val)
}



###################################
####  Age/sex incidence model  ####
###################################

#' @export
fnCreateParam <- function(theta, fp){

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if(fp$eppmod %in% c("rspline", "logrw")){

    epp_nparam <- fp$numKnots+1

    if(fp$eppmod == "rspline"){
      u <- theta[1:fp$numKnots]
      if(fp$rtpenord == 2){
        beta <- numeric(fp$numKnots)
        beta[1] <- u[1]
        beta[2] <- u[1]+u[2]
        for(i in 3:fp$numKnots)
          beta[i] <- -beta[i-2] + 2*beta[i-1] + u[i]
      } else # first order penalty
        beta <- cumsum(u)
    } else if(fp$eppmod %in% "logrw")
      beta <- theta[1:fp$numKnots]

    param <- list(beta = beta,
                  rvec = as.vector(fp$rvec.spldes %*% beta))

    if(fp$eppmod %in% "logrw")
      param$rvec <- exp(param$rvec)

    if(exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      param$iota <- exp(param$rvec[fp$proj.steps == fp$tsEpidemicStart] * theta[fp$numKnots+1])
    else
      param$iota <- transf_iota(theta[fp$numKnots+1], fp)

  } else if(fp$eppmod == "rlogistic") {
    epp_nparam <- 5
    par <- theta[1:4]
    par[3] <- exp(theta[3])
    param <- list()
    param$rvec <- exp(rlogistic(fp$proj.steps, par))
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

    if(!exists("ancrtcens.vinfl", fp)){
      param$ancrtcens.vinfl <- exp(theta[paramcurr+2])
      paramcurr <- paramcurr+2
    } else
      paramcurr <- paramcurr+1
  }
  if(exists("ancrt", fp) && fp$ancrt %in% c("site", "both")){
    param$ancrtsite.beta <- theta[paramcurr+1]
    paramcurr <- paramcurr+1
  }

  return(param)
}



########################################################
####  Age specific prevalence likelihood functions  ####
########################################################


#' Prepare age-specific HH survey prevalence likelihood data
prepare_hhsageprev_likdat <- function(hhsage, fp){
  anchor.year <- floor(min(fp$proj.steps))

  hhsage$W.hhs <- stats::qnorm(hhsage$prev)
  hhsage$v.hhs <- 2*pi*exp(hhsage$W.hhs^2)*hhsage$se^2
  hhsage$sd.W.hhs <- sqrt(hhsage$v.hhs)

  if(exists("deff_approx", hhsage))
    hhsage$n_eff <- hhsage$n/hhsage$deff_approx
  else if(exists("deff_approx", hhsage))
    hhsage$n_eff <- hhsage$n/hhsage$deff
  else
    hhsage$n_eff <- hhsage$prev * (1 - hhsage$prev) / hhsage$se ^ 2
  hhsage$x_eff <- hhsage$n_eff * hhsage$prev

  if(is.null(hhsage$sex))
    hhsage$sex <- rep("both", nrow(hhsage))

  if(is.null(hhsage$agegr))
    hhsage$agegr <- "15-49"

  startage <- as.integer(sub("([0-9]*)-([0-9]*)", "\\1", hhsage$agegr))
  endage <- as.integer(sub("([0-9]*)-([0-9]*)", "\\2", hhsage$agegr))

  hhsage$sidx <- match(hhsage$sex, c("both", "male", "female")) - 1L
  hhsage$aidx <- startage - fp$ss$AGE_START+1L
  hhsage$yidx <- as.integer(hhsage$year - (anchor.year - 1))
  hhsage$agspan <- endage - startage + 1L

  return(subset(hhsage, aidx > 0))
}

#' Log likelihood for age-specific household survey prevalence
#' @export
ll_hhsage <- function(mod, fp, dat, pointwise = FALSE){

  prevM.age <- ageprev(mod, aidx = dat$aidx, sidx = dat$sidx, yidx = dat$yidx, agspan = dat$agspan)
  ## If calendar year projection, average current and previous year prevalence to
  ## approximate mid-year prevalence
  if (fp$projection_period == "calendar") {
    prevM.age_last <- ageprev(mod, aidx = dat$aidx, sidx = dat$sidx, yidx = dat$yidx-1L, agspan = dat$agspan)
    prevM.age <- 0.5 * (prevM.age + prevM.age_last)
  }

  qM.age <- suppressWarnings(stats::qnorm(prevM.age))

  if(any(is.na(qM.age)))
    val <- rep(-Inf, nrow(dat))
  else
    val <- stats::dnorm(dat$W.hhs, qM.age, dat$sd.W.hhs, log=TRUE)

  if(pointwise)
    return(val)
  sum(val)
}


#' Log likelihood for age-specific household survey prevalence using binomial approximation
#' @export
ll_hhsage_binom <- function(mod, fp, dat, pointwise = FALSE){

  prevM.age <- suppressWarnings(ageprev(mod, aidx = dat$aidx, sidx = dat$sidx, yidx = dat$yidx, agspan = dat$agspan))

  ## If calendar year projection, average current and previous year prevalence to
  ## approximate mid-year prevalence
  if (fp$projection_period == "calendar") {
    prevM.age_last <- ageprev(mod, aidx = dat$aidx, sidx = dat$sidx, yidx = dat$yidx-1L, agspan = dat$agspan)
    prevM.age <- 0.5 * (prevM.age + prevM.age_last)
  }

  if(any(is.na(prevM.age)) || any(prevM.age >= 1))
    val <- rep(-Inf, nrow(dat))
  else
    val <- ldbinom(dat$x_eff, dat$n_eff, prevM.age)
  val[is.na(val)] <- -Inf

  if(pointwise)
    return(val)

  sum(val)
}


#' ## Household survey ART coverage likelihood

prepare_hhsartcov_likdat <- function(hhsartcov, fp){

  anchor.year <- floor(min(fp$proj.steps))

  hhsartcov$W.hhs <- stats::qnorm(hhsartcov$artcov)
  hhsartcov$v.hhs <- 2*pi*exp(hhsartcov$W.hhs^2)*hhsartcov$se^2
  hhsartcov$sd.W.hhs <- sqrt(hhsartcov$v.hhs)

  if(is.null(hhsartcov$sex))
    hhsartcov$sex <- rep("both", nrow(hhsartcov))

  if(is.null(hhsartcov$agegr))
    hhsartcov$agegr <- "15-49"

  startage <- as.integer(sub("([0-9]*)-([0-9]*)", "\\1", hhsartcov$agegr))
  endage <- as.integer(sub("([0-9]*)-([0-9]*)", "\\2", hhsartcov$agegr))

  hhsartcov$sidx <- match(hhsartcov$sex, c("both", "male", "female")) - 1L
  hhsartcov$aidx <- startage - fp$ss$AGE_START+1L
  hhsartcov$yidx <- as.integer(hhsartcov$year - (anchor.year - 1))
  hhsartcov$agspan <- endage - startage + 1L

  return(subset(hhsartcov, aidx > 0))
}


#' Log likelihood for age-specific household survey prevalence
ll_hhsartcov <- function(mod, fp, dat, pointwise = FALSE){

  artcovM <- artcov15to49(mod)
  artcovM_obs <- artcovM[dat$yidx]

  ## If calendar year projection, average current and previous year
  ## to approximate mid-year prevalence
  if (fp$projection_period == "calendar") {
    artcovM_obs <- 0.5 * (artcovM_obs + artcovM[dat$yidx-1L])
  }

  qM <- stats::qnorm(artcovM_obs)

  if(any(is.na(qM)))
    val <- rep(-Inf, nrow(dat))
  else
    val <- stats::dnorm(dat$W.hhs, qM, dat$sd.W.hhs, log=TRUE)

  if(pointwise)
    return(val)
  sum(val)
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
ll_hhsincid <- function(mod, fp, hhsincid.dat){
  logincid <- log(incid(mod, fp))
  ll.incid <- sum(stats::dnorm(hhsincid.dat$log_incid, logincid[hhsincid.dat$idx], hhsincid.dat$log_incid.se, TRUE))
  return(ll.incid)
}


###############################
####  Likelihood function  ####
###############################

#' @export
prepare_likdat <- function(eppd, fp){

  likdat <- list()

  likdat$hhs.dat <- prepare_hhsageprev_likdat(eppd$hhs, fp)

  if(exists("ancsitedat", where=eppd)){

    ancsitedat <- eppd$ancsitedat

    if(exists("ancrt", fp) && fp$ancrt %in% c("none", "census"))
      ancsitedat <- subset(ancsitedat, type == "ancss")

    likdat$ancsite.dat <- prepare_ancsite_likdat(ancsitedat, fp)
  }

  if(exists("ancrtcens", where=eppd)){
    if(exists("ancrt", fp) && fp$ancrt %in% c("none", "site"))
      eppd$ancrtcens <- NULL
    else
      likdat$ancrtcens.dat <- prepare_ancrtcens_likdat(eppd$ancrtcens, fp)
  }

  if(exists("hhsincid", where=eppd))
    likdat$hhsincid.dat <- prepare_hhsincid_likdat(eppd$hhsincid, fp)

  if(exists("hhsartcov", where=eppd))
    likdat$hhsartcov.dat <- prepare_hhsartcov_likdat(eppd$hhsartcov, fp)

  return(likdat)
}



#' @export
lprior <- function(theta, fp){

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if(fp$eppmod %in% c("rspline", "logrw")){
    epp_nparam <- fp$numKnots+1

    nk <- fp$numKnots

    if(fp$eppmod == "logrw")
      lpr <- bayes_lmvt(theta[2:fp$numKnots], rw_prior_shape, rw_prior_rate)
    else
      lpr <- bayes_lmvt(theta[(1+fp$rtpenord):nk], tau2_prior_shape, tau2_prior_rate)

    if(exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      lpr <- lpr + stats::dunif(theta[nk+1], r0logiotaratio.unif.prior[1], r0logiotaratio.unif.prior[2], log=TRUE)
    else
      lpr <- lpr + lprior_iota(theta[nk+1], fp)

  } else if(fp$eppmod == "rlogistic") {
    epp_nparam <- 5
    lpr <- sum(stats::dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE))
    lpr <- lpr + lprior_iota(theta[5], fp)
  } else if(fp$eppmod == "rtrend"){ # rtrend

    epp_nparam <- 7

    lpr <- stats::dunif(round(theta[1]), t0.unif.prior[1], t0.unif.prior[2], log=TRUE) +
      stats::dnorm(round(theta[2]), t1.pr.mean, t1.pr.sd, log=TRUE) +
      stats::dnorm(theta[3], logr0.pr.mean, logr0.pr.sd, log=TRUE) +
      sum(stats::dnorm(theta[4:7], rtrend.beta.pr.mean, rtrend.beta.pr.sd, log=TRUE))
  } else if(fp$eppmod == "rhybrid"){
    epp_nparam <- fp$rt$n_param+1
    lpr <- sum(stats::dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE)) +
      sum(stats::dnorm(theta[4+1:fp$rt$n_rw], 0, rw_prior_sd, log=TRUE))
    lpr <- lpr + lprior_iota(theta[fp$rt$n_param+1], fp)
  }

  if(fp$ancsitedata){
    lpr <- lpr + stats::dnorm(theta[epp_nparam+1], ancbias.pr.mean, ancbias.pr.sd, log=TRUE)
    if(!exists("v.infl", where=fp)){
      anclik_nparam <- 2
      lpr <- lpr + stats::dexp(exp(theta[epp_nparam+2]), vinfl.prior.rate, TRUE) + theta[epp_nparam+2]         # additional ANC variance
    } else
      anclik_nparam <- 1
  } else
    anclik_nparam <- 0

  paramcurr <- epp_nparam+anclik_nparam
  if(exists("ancrt", fp) && fp$ancrt %in% c("census", "both")){
    lpr <- lpr + stats::dnorm(theta[paramcurr+1], log_frr_adjust.pr.mean, log_frr_adjust.pr.sd, log=TRUE)
    if(!exists("ancrtcens.vinfl", fp)){
      lpr <- lpr + stats::dexp(exp(theta[paramcurr+2]), ancrtcens.vinfl.pr.rate, TRUE) + theta[paramcurr+2]
      paramcurr <- paramcurr+2
    } else
      paramcurr <- paramcurr+1
  }
  if(exists("ancrt", fp) && fp$ancrt %in% c("site", "both")){
    lpr <- lpr + stats::dnorm(theta[paramcurr+1], ancrtsite.beta.pr.mean, ancrtsite.beta.pr.sd, log=TRUE)
    paramcurr <- paramcurr+1
  }

  return(lpr)
}


#' @export
ll <- function(theta, fp, likdat){
  theta.last <<- theta
  fp <- stats::update(fp, list=fnCreateParam(theta, fp))

  if (fp$eppmod == "rspline")
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
      ll.hhs <- ll_hhsage_binom(mod, fp, likdat$hhs.dat)
    else ## use probit likelihood
      ll.hhs <- ll_hhsage(mod, fp, likdat$hhs.dat) # probit-transformed model
  else
    ll.hhs <- 0

  if(!is.null(likdat$hhsincid.dat))
    ll.incid <- ll_hhsincid(mod, fp, likdat$hhsincid.dat)
  else
    ll.incid <- 0

  if(!is.null(likdat$hhsartcov.dat))
    ll.hhsartcov <- ll_hhsartcov(mod, fp, likdat$hhsartcov.dat)
  else
    ll.hhsartcov <- 0


  if(exists("equil.rprior", where=fp) && fp$equil.rprior){
    if(fp$eppmod != "rspline")
      stop("error in ll(): equil.rprior is only for use with r-spline model")

    lastdata.idx <- max(likdat$ancsite.dat$df$yidx,
                        likdat$hhs.dat$yidx,
                        likdat$ancrtcens.dat$yidx,
                        likdat$hhsincid.dat$idx)

    qM.all <- suppressWarnings(stats::qnorm(prev(mod)))

    if(any(is.na(qM.all[lastdata.idx - 9:0]))) {
      ll.rprior <- -Inf
    } else {
      rvec.ann <- fp$rvec[fp$proj.steps %% 1 == 0.5]
      equil.rprior.mean <- epp:::muSS/(1-stats::pnorm(qM.all[lastdata.idx]))
      if(!is.null(fp$prior_args$equil.rprior.sd))
        equil.rprior.sd <- fp$prior_args$equil.rprior.sd
      else
        equil.rprior.sd <- sqrt(mean((epp:::muSS/(1-stats::pnorm(qM.all[lastdata.idx - 9:0])) - rvec.ann[lastdata.idx - 9:0])^2))  # empirical sd based on 10 previous years

      ll.rprior <- sum(stats::dnorm(rvec.ann[(lastdata.idx+1L):length(qM.all)], equil.rprior.mean, equil.rprior.sd, log=TRUE))  # prior starts year after last data
    }
  } else
    ll.rprior <- 0

  c(anc    = ll.anc,
    ancrt  = ll.ancrt,
    hhs    = ll.hhs,
    incid  = ll.incid,
    artcov = ll.hhsartcov,
    rprior = ll.rprior)
}


##########################
####  IMIS functions  ####
##########################

sample.prior <- function(n, fp){

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  ## Calculate number of parameters
  if(fp$eppmod %in% c("rspline", "logrw"))
    epp_nparam <- fp$numKnots+1L
  else if(fp$eppmod == "rlogistic")
    epp_nparam <- 5
  else if(fp$eppmod == "rtrend")
    epp_nparam <- 7
  else if(fp$eppmod == "rhybrid")
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

  ## Create matrix for storing samples
  mat <- matrix(NA, n, nparam)

  if(fp$eppmod %in% c("rspline", "logrw")){
    epp_nparam <- fp$numKnots+1

    if(fp$eppmod == "rspline")
      mat[,1] <- stats::rnorm(n, 1.5, 1)                                                   # u[1]
    else # logrw
      mat[,1] <- stats::rnorm(n, 0.2, 1)                                                   # u[1]
    if(fp$eppmod == "logrw"){
      mat[,2:fp$rt$n_rw] <- bayes_rmvt(n, fp$rt$n_rw-1, rw_prior_shape, rw_prior_rate)  # u[2:numKnots]
    } else {
      mat[,2:fp$numKnots] <- bayes_rmvt(n, fp$numKnots-1,tau2_init_shape, tau2_init_rate)  # u[2:numKnots]
    }

    if(exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      mat[,fp$numKnots+1] <-  stats::runif(n, r0logiotaratio.unif.prior[1], r0logiotaratio.unif.prior[2])  # ratio r0 / log(iota)
    else
      mat[,fp$numKnots+1] <- sample_iota(n, fp)
  } else if(fp$eppmod == "rlogistic"){
    mat[,1:4] <- t(matrix(stats::rnorm(4*n, rlog_pr_mean, rlog_pr_sd), 4))
    mat[,5] <- sample_iota(n, fp)
  } else if(fp$eppmod == "rtrend"){ # r-trend

    mat[,1] <- stats::runif(n, t0.unif.prior[1], t0.unif.prior[2])           # t0
    mat[,2] <- stats::rnorm(n, t1.pr.mean, t1.pr.sd)
    mat[,3] <- stats::rnorm(n, logr0.pr.mean, logr0.pr.sd)  # r0
    mat[,4:7] <- t(matrix(stats::rnorm(4*n, rtrend.beta.pr.mean, rtrend.beta.pr.sd), 4, n))  # beta
  } else if(fp$eppmod == "rhybrid") {
    mat[,1:4] <- t(matrix(stats::rnorm(4*n, rlog_pr_mean, rlog_pr_sd), 4))
    mat[,4+1:fp$rt$n_rw] <- stats::rnorm(n*fp$rt$n_rw, 0, rw_prior_sd)  # u[2:numKnots]
    mat[,fp$rt$n_param+1] <- sample_iota(n, fp)
  }

  ## sample ANC bias paramters
  if(fp$ancsitedata){
    mat[,epp_nparam+1] <- stats::rnorm(n, ancbias.pr.mean, ancbias.pr.sd)   # ancbias parameter
    if(!exists("v.infl", where=fp))
      mat[,epp_nparam+2] <- log(stats::rexp(n, vinfl.prior.rate))
  }

  ## sample ANCRT parameters
  paramcurr <- epp_nparam+anclik_nparam
  if(exists("ancrt", where=fp) && fp$ancrt %in% c("census", "both")){
    mat[,paramcurr+1] <- stats::rnorm(n, log_frr_adjust.pr.mean, log_frr_adjust.pr.sd)
    if(!exists("ancrtcens.vinfl", fp)){
      mat[,paramcurr+2] <- log(stats::rexp(n, ancrtcens.vinfl.pr.rate))
      paramcurr <- paramcurr+2
    } else
      paramcurr <- paramcurr+1
  }
  if(exists("ancrt", where=fp) && fp$ancrt %in% c("site", "both")){
    mat[,paramcurr+1] <- stats::rnorm(n, ancrtsite.beta.pr.mean, ancrtsite.beta.pr.sd)
    paramcurr <- paramcurr+1
  }

  return(mat)
}

ldsamp <- function(theta, fp){

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if(fp$eppmod %in% c("rspline", "logrw")){
    epp_nparam <- fp$numKnots+1

    nk <- fp$numKnots

    if(fp$eppmod == "rspline")  # u[1]
      lpr <- stats::dnorm(theta[1], 1.5, 1, log=TRUE)
    else # logrw
      lpr <- stats::dnorm(theta[1], 0.2, 1, log=TRUE)

    if(fp$eppmod == "logrw")
      bayes_lmvt(theta[2:fp$rt$n_rw], rw_prior_shape, rw_prior_rate)
    else
      lpr <- bayes_lmvt(theta[2:nk], tau2_prior_shape, tau2_prior_rate)


    if(exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      lpr <- lpr + stats::dunif(theta[nk+1], r0logiotaratio.unif.prior[1], r0logiotaratio.unif.prior[2], log=TRUE)
    else
      lpr <- lpr + ldsamp_iota(theta[nk+1], fp)

  } else if(fp$eppmod == "rlogistic") {
    epp_nparam <- 5
    lpr <- sum(stats::dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE))
    lpr <- lpr + ldsamp_iota(theta[5], fp)
  } else if(fp$eppmod == "rtrend"){ # rtrend

    epp_nparam <- 7

    lpr <- stats::dunif(round(theta[1]), t0.unif.prior[1], t0.unif.prior[2], log=TRUE) +
      stats::dnorm(round(theta[2]), t1.pr.mean, t1.pr.sd, log=TRUE) +
      stats::dnorm(theta[3], logr0.pr.mean, logr0.pr.sd, log=TRUE) +
      sum(stats::dnorm(theta[4:7], rtrend.beta.pr.mean, rtrend.beta.pr.sd, log=TRUE))
  } else if(fp$eppmod == "rhybrid"){
    epp_nparam <- fp$rt$n_param+1
    lpr <- sum(stats::dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE)) +
      sum(stats::dnorm(theta[4+1:fp$rt$n_rw], 0, rw_prior_sd, log=TRUE))
    lpr <- lpr + ldsamp_iota(theta[fp$rt$n_param+1], fp)
  }

  if(fp$ancsitedata){
    lpr <- lpr + stats::dnorm(theta[epp_nparam+1], ancbias.pr.mean, ancbias.pr.sd, log=TRUE)
    if(!exists("v.infl", where=fp)){
      anclik_nparam <- 2
      lpr <- lpr + stats::dexp(exp(theta[epp_nparam+2]), vinfl.prior.rate, TRUE) + theta[epp_nparam+2]         # additional ANC variance
    } else
      anclik_nparam <- 1
  } else
    anclik_nparam <- 0

  paramcurr <- epp_nparam+anclik_nparam
  if(exists("ancrt", fp) && fp$ancrt %in% c("census", "both")){
    lpr <- lpr + stats::dnorm(theta[paramcurr+1], log_frr_adjust.pr.mean, log_frr_adjust.pr.sd, log=TRUE)
    if(!exists("ancrtcens.vinfl", fp)){
      lpr <- lpr + stats::dexp(exp(theta[paramcurr+2]), ancrtcens.vinfl.pr.rate, TRUE) + theta[paramcurr+2]
      paramcurr <- paramcurr+2
    } else
      paramcurr <- paramcurr+1
  }
  if(exists("ancrt", fp) && fp$ancrt %in% c("site", "both")){
    lpr <- lpr + stats::dnorm(theta[paramcurr+1], ancrtsite.beta.pr.mean, ancrtsite.beta.pr.sd, log=TRUE) ## +
    paramcurr <- paramcurr+1
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
    lval <- sum(ll(theta, fp, likdat))
  else
    lval <- unlist(lapply(seq_len(nrow(theta)), function(i) sum(ll(theta[i,], fp, likdat))))
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
