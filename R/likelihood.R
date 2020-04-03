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




###################################
####  Age/sex incidence model  ####
###################################

fnCreateParam <- function(theta, fp){

  if (exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }
  
  if (fp$eppmod %in% c("rspline", "logrw")){

    epp_nparam <- fp$numKnots+1

    if (fp$eppmod == "rspline"){
      u <- theta[1:fp$numKnots]
      if (fp$rtpenord == 2){
        beta <- numeric(fp$numKnots)
        beta[1] <- u[1]
        beta[2] <- u[1]+u[2]
        for(i in 3:fp$numKnots)
          beta[i] <- -beta[i-2] + 2*beta[i-1] + u[i]
      } else # first order penalty
        beta <- cumsum(u)
    } else if (fp$eppmod %in% "logrw")
      beta <- theta[1:fp$numKnots]

    param <- list(beta = beta,
                  rvec = as.vector(fp$rvec.spldes %*% beta))
    
    if (fp$eppmod %in% "logrw")
      param$rvec <- exp(param$rvec)

    if (exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      param$iota <- exp(param$rvec[fp$proj.steps == fp$tsEpidemicStart] * theta[fp$numKnots+1])
    else
      param$iota <- transf_iota(theta[fp$numKnots+1], fp)
  } else if (fp$eppmod == "rlogistic") {
    epp_nparam <- 5
    par <- theta[1:4]
    par[3] <- exp(theta[3])
    param <- list()
    param$rvec <- exp(rlogistic(fp$proj.steps, par))
    param$iota <- transf_iota(theta[5], fp)
  } else if (fp$eppmod == "rtrend"){ # rtrend
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

  if (exists("ancmod", fp) && fp$ancmod$nparam > 0) {
    theta_anc <- theta[epp_nparam+1:fp$ancmod$nparam]
    param <- update_par(param, list = create_ancmod_param(theta_anc, fp$ancmod))
    param$frr_cd4 <- fp$frr_cd4 * exp(param$log_frr_adjust)
    param$frr_art <- fp$frr_art * exp(param$log_frr_adjust)
  }
  
  return(param)
}



########################################################
####  Age specific prevalence likelihood functions  ####
########################################################


#' Prepare age-specific HH survey prevalence likelihood data
#' 
#' @param hhsage age-specific HH survey prevalence likelihood data
#' @param fp fix parameters
prepare_hhsageprev_likdat <- function(hhsage, fp){
  anchor.year <- floor(min(fp$proj.steps))

  hhsage$W.hhs <- qnorm(hhsage$prev)
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
#' 
#' @param mod model simulation output
#' @param dat Output data from prepare_likdat
#' @param pointwise Point-wise likelihood
ll_hhsage <- function(mod, dat, pointwise = FALSE){
  qM.age <- suppressWarnings(qnorm(ageprev(mod,
                                           aidx = dat$aidx,
                                           sidx = dat$sidx, 
                                           yidx = dat$yidx, 
                                           agspan = dat$agspan)))
  if (any(is.na(qM.age)))
    val <- rep(-Inf, nrow(dat))
  else
    val <- dnorm(dat$W.hhs, qM.age, dat$sd.W.hhs, log=TRUE)

  if (pointwise)
    return(val)
  sum(val)
}

#' Log likelihood for age-specific household survey prevalence using binomial approximation
#' 
#' @param mod model simulation output
#' @param dat Output data from prepare_likdat
#' @param pointwise Point-wise likelihood
ll_hhsage_binom <- function(mod, dat, pointwise = FALSE){

  prevM.age <- suppressWarnings(ageprev(mod, dat$aidx, dat$sidx, dat$yidx, dat$agspan))

  if (any(is.na(prevM.age)) || any(prevM.age >= 1))
    val <- rep(-Inf, nrow(dat))
  else
    val <- ldbinom(dat$x_eff, dat$n_eff, prevM.age)
  val[is.na(val)] <- -Inf

  if (pointwise)
    return(val)

  sum(val)
}



#########################################
####  Incidence likelihood function  ####
#########################################

#' Prepare household survey incidence likelihood data
#' 
#' @param hhsincid household survey incidence likelihood data
#' @param fp fix parameters
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

  likdat <- list()
  
  likdat$hhs.dat <- prepare_hhsageprev_likdat(eppd$hhs, fp)

  if (exists("ancsitedat", where=eppd)){

    ancsitedat <- eppd$ancsitedat
    
    if (exists("ancmod", fp) && !fp$ancmod$has_ancrtsite)
      ancsitedat <- subset(ancsitedat, type == "ancss")
    
    likdat$ancsite.dat <- prepare_ancsite_likdat(ancsitedat, fp)
  }
 
  if (exists("ancrtcens", where=eppd)){
    if (exists("ancmod", fp) && !fp$ancmod$has_ancrtcens)
      eppd$ancrtcens <- NULL
    else
      likdat$ancrtcens.dat <- prepare_ancrtcens_likdat(eppd$ancrtcens, fp)
  }

  if (exists("hhsincid", where=eppd))
    likdat$hhsincid.dat <- prepare_hhsincid_likdat(eppd$hhsincid, fp)

  return(likdat)
}

lprior <- function(theta, fp){

  if (exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if (fp$eppmod %in% c("rspline", "logrw")){
    epp_nparam <- fp$numKnots+1

    nk <- fp$numKnots

    if (fp$eppmod == "logrw")
      lpr <- bayes_lmvt(theta[2:fp$numKnots], rw_prior_shape, rw_prior_rate)
    else
      lpr <- bayes_lmvt(theta[(1+fp$rtpenord):nk], tau2_prior_shape, tau2_prior_rate)

    if (exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      lpr <- lpr + dunif (theta[nk+1], r0logiotaratio.unif.prior[1], r0logiotaratio.unif.prior[2], log=TRUE)
    else
      lpr <- lpr + lprior_iota(theta[nk+1], fp)
  } else if (fp$eppmod == "rlogistic") {
    epp_nparam <- 5
    lpr <- sum(dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE))
    lpr <- lpr + lprior_iota(theta[5], fp)
  } else if (fp$eppmod == "rtrend"){ # rtrend    
    epp_nparam <- 7
    
    lpr <- dunif (round(theta[1]), t0.unif.prior[1], t0.unif.prior[2], log=TRUE) +
      dnorm(round(theta[2]), t1.pr.mean, t1.pr.sd, log=TRUE) +
      dnorm(theta[3], logr0.pr.mean, logr0.pr.sd, log=TRUE) +
      sum(dnorm(theta[4:7], rtrend.beta.pr.mean, rtrend.beta.pr.sd, log=TRUE))
  } else if (fp$eppmod == "rhybrid"){
    epp_nparam <- fp$rt$n_param+1
    lpr <- sum(dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE)) +
      sum(dnorm(theta[4+1:fp$rt$n_rw], 0, rw_prior_sd, log=TRUE))
    lpr <- lpr + lprior_iota(theta[fp$rt$n_param+1], fp)
  }

  if (exists("ancmod", fp) && fp$ancmod$nparam > 0) {
    theta_anc <- theta[epp_nparam+1:fp$ancmod$nparam]
    lpr <- lpr + lprior_ancmod(theta_anc, fp$ancmod, fp$prior_args)
  }

  if (fp$ss$MIX)
    lpr <- lpr + dbeta(tail(theta, 1), 2, 2, log=TRUE)

  return(lpr)
}

#' Full log-likelhood
#' 
#' @importFrom stats aggregate approx cov cov.wt density dexp dlnorm dnorm dunif ecdf mahalanobis median model.matrix na.omit optim optimHess pnorm qnorm quantile relevel rexp rgamma rnorm runif sd setNames update var
ll_all = function(theta, fp, likdat) {
  theta.last <<- theta

  nparam <- length(theta)
  
  if (fp$ss$MIX) {
    fp$balancing <- tail(theta, 1)
    fp <- update(fp, list=fnCreateParam(theta[1:(nparam-1)], fp))
  } else {
    fp <- update(fp, list=fnCreateParam(theta, fp))
  }

  if (fp$eppmod == "rspline")
    if (any(is.na(fp$rvec)) || min(fp$rvec) < 0 || max(fp$rvec) > 20) 
      return(-Inf)

  mod <- simmod(fp)

  .anc = .ancrt = .hhs = .incid = .rprior = 0

  ## ANC likelihood
  if (exists("ancsite.dat", likdat))
    .anc <- ll_ancsite(mod, fp, coef=c(fp$ancbias, fp$ancrtsite.beta),
                       vinfl=fp$v.infl, likdat$ancsite.dat)

  if (exists("ancrtcens.dat", likdat))
    .ancrt <- ll_ancrtcens(mod, likdat$ancrtcens.dat, fp)

  ## Household survey likelihood
  if (exists("hhs.dat", where=likdat)) {
    if (exists("ageprev", fp) && fp$ageprev=="binom")
      .hhs <- ll_hhsage_binom(mod, likdat$hhs.dat)
    else ## use probit likelihood
      .hhs <- ll_hhsage(mod, likdat$hhs.dat) # probit-transformed model
  }

  if (!is.null(likdat$hhsincid.dat))
    .incid <- ll_hhsincid(mod, likdat$hhsincid.dat)

  if (exists("equil.rprior", where=fp) && fp$equil.rprior) {
    if (fp$eppmod != "rspline")
      stop("error in ll(): equil.rprior is only for use with r-spline model")

    lastdata.idx <- max(likdat$ancsite.dat$df$yidx,
                        likdat$hhs.dat$yidx,
                        likdat$ancrtcens.dat$yidx,
                        likdat$hhsincid.dat$idx)
    
    qM.all <- suppressWarnings(qnorm(prev(mod)))

    if (any(is.na(qM.all[lastdata.idx - 9:0]))) {
      .rprior <- -Inf
    } 
    else {
      rvec.ann <- fp$rvec[fp$proj.steps %% 1 == 0.5]
      equil.rprior.mean <- epp:::muSS/(1-pnorm(qM.all[lastdata.idx]))
      if (!is.null(fp$prior_args$equil.rprior.sd))
        equil.rprior.sd <- fp$prior_args$equil.rprior.sd
      else
        equil.rprior.sd <- sqrt(mean((epp:::muSS / (1-pnorm(qM.all[lastdata.idx - 9:0])) - rvec.ann[lastdata.idx - 9:0])^2))  # empirical sd based on 10 previous years
      
      .rprior <- sum(dnorm(rvec.ann[(lastdata.idx+1L):length(qM.all)], equil.rprior.mean, equil.rprior.sd, log=TRUE))  # prior starts year after last data
    }
  }

  c(anc = .anc, ancrt = .ancrt, hhs = .hhs, incid = .incid, rprior = .rprior) 
}


##########################
####  IMIS functions  ####
##########################

sample.prior <- function(n, fp){

  if (exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  ## Calculate number of parameters
  if (fp$eppmod %in% c("rspline", "logrw"))
    epp_nparam <- fp$numKnots+1L
  else if (fp$eppmod == "rlogistic")
    epp_nparam <- 5
  else if (fp$eppmod == "rtrend")
    epp_nparam <- 7
  else if (fp$eppmod == "rhybrid")
    epp_nparam <- fp$rt$n_param+1

  nparam <- epp_nparam + fp$ancmod$nparam

  if (fp$ss$MIX)
    nparam <- nparam + 1
  
  ## Create matrix for storing samples
  mat <- matrix(NA, n, nparam)

  if (fp$eppmod %in% c("rspline", "logrw")){
    epp_nparam <- fp$numKnots+1

    if (fp$eppmod == "rspline")
      mat[,1] <- rnorm(n, 1.5, 1)                                                   # u[1]
    else # logrw
      mat[,1] <- rnorm(n, 0.2, 1)                                                   # u[1]
    if (fp$eppmod == "logrw"){
      mat[,2:fp$rt$n_rw] <- bayes_rmvt(n, fp$rt$n_rw-1, rw_prior_shape, rw_prior_rate)  # u[2:numKnots]
    } else {
      mat[,2:fp$numKnots] <- bayes_rmvt(n, fp$numKnots-1,tau2_init_shape, tau2_init_rate)  # u[2:numKnots]
    }

    if (exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      mat[,fp$numKnots+1] <-  runif (n, r0logiotaratio.unif.prior[1], r0logiotaratio.unif.prior[2])  # ratio r0 / log(iota)
    else
      mat[,fp$numKnots+1] <- sample_iota(n, fp)
  } else if (fp$eppmod == "rlogistic"){
    mat[,1:4] <- t(matrix(rnorm(4*n, rlog_pr_mean, rlog_pr_sd), 4))
    mat[,5] <- sample_iota(n, fp)
  } else if (fp$eppmod == "rtrend"){ # r-trend

    mat[,1] <- runif (n, t0.unif.prior[1], t0.unif.prior[2])           # t0
    mat[,2] <- rnorm(n, t1.pr.mean, t1.pr.sd)
    mat[,3] <- rnorm(n, logr0.pr.mean, logr0.pr.sd)  # r0
    mat[,4:7] <- t(matrix(rnorm(4*n, rtrend.beta.pr.mean, rtrend.beta.pr.sd), 4, n))  # beta
  } else if (fp$eppmod == "rhybrid") {
    mat[,1:4] <- t(matrix(rnorm(4*n, rlog_pr_mean, rlog_pr_sd), 4))
    mat[,4+1:fp$rt$n_rw] <- rnorm(n*fp$rt$n_rw, 0, rw_prior_sd)  # u[2:numKnots]
    mat[,fp$rt$n_param+1] <- sample_iota(n, fp)
  }

  ## sample ANC model parameters
  if (exists("ancmod", fp) && fp$ancmod$nparam > 0)
    mat[ , epp_nparam + 1:fp$ancmod$nparam] <- sample_prior_ancmod(n, fp$ancmod, fp$prior_args)

  # sex acts balancing parameter
  if (fp$ss$MIX)
    mat[, nparam] <- rbeta(n, 2, 2)
  
  return(mat)
}

ldsamp <- function(theta, fp){

  if (exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if (fp$eppmod %in% c("rspline", "logrw")){
    epp_nparam <- fp$numKnots+1

    nk <- fp$numKnots

    if (fp$eppmod == "rspline")  # u[1]
      lpr <- dnorm(theta[1], 1.5, 1, log=TRUE)
    else # logrw
      lpr <- dnorm(theta[1], 0.2, 1, log=TRUE)

    if (fp$eppmod == "logrw")
      bayes_lmvt(theta[2:fp$rt$n_rw], rw_prior_shape, rw_prior_rate)
    else
      lpr <- bayes_lmvt(theta[2:nk], tau2_prior_shape, tau2_prior_rate)


    if (exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      lpr <- lpr + dunif (theta[nk+1], r0logiotaratio.unif.prior[1], r0logiotaratio.unif.prior[2], log=TRUE)
    else
      lpr <- lpr + ldsamp_iota(theta[nk+1], fp)

  } else if (fp$eppmod == "rlogistic") {
    epp_nparam <- 5
    lpr <- sum(dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE))
    lpr <- lpr + ldsamp_iota(theta[5], fp)
  } else if (fp$eppmod == "rtrend"){ # rtrend

    epp_nparam <- 7

    lpr <- dunif (round(theta[1]), t0.unif.prior[1], t0.unif.prior[2], log=TRUE) +
      dnorm(round(theta[2]), t1.pr.mean, t1.pr.sd, log=TRUE) +
      dnorm(theta[3], logr0.pr.mean, logr0.pr.sd, log=TRUE) +
      sum(dnorm(theta[4:7], rtrend.beta.pr.mean, rtrend.beta.pr.sd, log=TRUE))
  } else if (fp$eppmod == "rhybrid"){
    epp_nparam <- fp$rt$n_param+1
    lpr <- sum(dnorm(theta[1:4], rlog_pr_mean, rlog_pr_sd, log=TRUE)) +
      sum(dnorm(theta[4+1:fp$rt$n_rw], 0, rw_prior_sd, log=TRUE))
    lpr <- lpr + ldsamp_iota(theta[fp$rt$n_param+1], fp)
  }

  if (exists("ancmod", fp) && fp$ancmod$nparam > 0) {
    theta_anc <- theta[epp_nparam+1:fp$ancmod$nparam]
    lpr <- lpr + lprior_ancmod(theta_anc, fp$ancmod, fp$prior_args)
  }

  return(lpr)
}

prior <- function(theta, fp, log=FALSE){
  if (is.vector(theta))
    lval <- lprior(theta, fp)
  else
    lval <- unlist(lapply(seq_len(nrow(theta)), function(i) (lprior(theta[i,], fp))))
  if (log)
    return(lval)
  else
    return(exp(lval))
}

likelihood <- function(theta, fp, likdat, log=FALSE, doParallel=FALSE) {
  if (is.vector(theta)) {
    lval <- sum(ll_all(theta, fp, likdat))
  } else {
    theta_id <- seq_len(nrow(theta))
    ll_fn    <- function(i) sum(ll_all(theta[i,], fp, likdat))
    if (!.Platform$OS.type=='unix' | !doParallel) {
      lval <- unlist(lapply(theta_id, ll_fn))
    } else {
      n_cores = max(floor(parallel::detectCores()/2), 1)
      cat('finding starting values on', n_cores, 'cores)\n')
      lval = unlist(parallel::mclapply(theta_id, ll_fn, mc.cores = n_cores))
    }
  }
  if (log)
    return(lval)
  else
    return(exp(lval))
}

dsamp <- function(theta, fp, log=FALSE){
  if (is.vector(theta))
    lval <- ldsamp(theta, fp)
  else
    lval <- unlist(lapply(seq_len(nrow(theta)), function(i) (ldsamp(theta[i,], fp))))
  if (log)
    return(lval)
  else
    return(exp(lval))
}
