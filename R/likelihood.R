
################################################################
####  Parameterising age/sex-specific incidence rate ratio  ####
################################################################

ldinvgamma <- function(x, alpha, beta){
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
  return(log.density)
}

## Binomial distribution log-density permitting non-integer counts
ldbinom <- function(x, size, prob){
  lgamma(size+1) - lgamma(x+1) - lgamma(size-x+1) + x*log(prob) + (size-x)*log(1-prob)
}

## r-spline prior parameters
logiota.unif.prior <- c(log(1e-14), log(0.0025))
tau2.prior.rate <- 0.5

invGammaParameter <- 0.001   #Inverse gamma parameter for tau^2 prior for spline
muSS <- 1/11.5               #1/duration for r steady state prior

ancbias.pr.mean <- 0.15
ancbias.pr.sd <- 1.0

## r-trend prior parameters
t0.unif.prior <- c(1970, 1990)
t1.unif.prior <- c(10, 30)
logr0.unif.prior <- c(1/11.5, 10)
rtrend.beta.pr.sd <- 0.2

vinfl.prior.rate <- 1/0.015

sexincrr.pr.mean <- log(1.38)
sexincrr.pr.sd <- 0.2

ageincrr.pr.mean <- c(-1.40707274, -0.23518703, 0.69314718, 0.78845736, -0.39975544, -0.70620810, -0.84054571, -0.02101324, -0.16382449, -0.37914407, -0.59639985, -0.82038300)
ageincrr.pr.sd <- 0.5



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

fnCreateParam <- function(theta, fp){

  
  if(!exists("eppmod", where = fp))  # backward compatibility
    fp$eppmod <- "rspline"
  
  if(fp$eppmod == "rspline"){
    u <- theta[1:fp$numKnots]
    beta <- numeric(fp$numKnots)
    beta[1] <- u[1]
    beta[2] <- u[1]+u[2]
    for(i in 3:fp$numKnots)
      beta[i] <- -beta[i-2] + 2*beta[i-1] + u[i]
    
    param <- list(rvec = as.vector(fp$rvec.spldes %*% beta),
                  iota = exp(theta[fp$numKnots+1]),
                  ancbias = theta[fp$numKnots+2],
                  v.infl = exp(theta[fp$numKnots+4]))
  } else { # rtrend
    param <- list(tsEpidemicStart = fp$proj.steps[which.min(abs(fp$proj.steps - theta[1]))], # t0
                  rtrend = list(tStabilize = theta[1]+theta[2],  # t0 + t1
                                r0 = exp(theta[3]),              # r0
                                beta = theta[4:7]),
                  ancbias = theta[8],
                  v.infl = exp(theta[9]))
  }

  if(inherits(fp, "specfp")){
    if(exists("fitincrr", where=fp) && fp$fitincrr){
      nparam <- length(theta)
      param$sigma_agepen <- exp(theta[nparam])

      param$incrr_sex <- fp$incrr_sex
      param$incrr_sex[] <- exp(theta[nparam-13])

      param$logincrr_age <- array(0, c(7, 2))
      param$logincrr_age[-3,] <- theta[nparam-12:1]

      param$incrr_age <- fp$incrr_age
      if(exists("linincrr", where=fp) && fp$linincrr)
        param$incrr_age[fp$ss$p.age15to49.idx,,] <- exp(apply(param$logincrr_age, 2, function(x) approx(3:9*5, x, 15:49, rule=2)$y))
      else
        param$incrr_age[fp$ss$p.age15to49.idx,,] <- apply(exp(param$logincrr_age), 2, rep, each=5)
      param$incrr_age[36:66,,] <- sweep(fp$incrr_age[36:66,,fp$ss$PROJ_YEARS], 2,
                                        param$incrr_age[35,,fp$ss$PROJ_YEARS]/fp$incrr_age[35,,fp$ss$PROJ_YEARS], "*")
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

  hhsage$sidx <- as.integer(hhsage$sex)
  hhsage$aidx <- 5*(as.integer(hhsage$agegr)-1) - fp$ss$AGE_START+1L
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
  sum(dnorm(hhsage.dat$W.hhs, qM.age, hhsage.dat$sd.W.hhs, log=TRUE))
}


#' Log likelihood for age-specific household survey prevalence using binomial approximation
ll_hhsage_binom <- function(mod, hhsage.dat){
  prevM.age <- suppressWarnings(ageprev(mod, arridx=hhsage.dat$arridx, agspan=5))
  sum(ldbinom(hhsage.dat$x_eff, hhsage.dat$n_eff, prevM.age))
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



###############################
####  Likelihood function  ####
###############################


lprior <- function(theta, fp){


  if(fp$eppmod == "rspline"){
    nk <- fp$numKnots
    tau2 <- exp(theta[nk+3])
    
    lpr <- sum(dnorm(theta[3:nk], 0, sqrt(tau2), log=TRUE)) +
      dunif(theta[nk+1], logiota.unif.prior[1], logiota.unif.prior[2], log=TRUE) + 
      dnorm(theta[nk+2], ancbias.pr.mean, ancbias.pr.sd, log=TRUE) +
      ldinvgamma(tau2, invGammaParameter, invGammaParameter) + log(tau2) +  # + log(tau2): multiply likelihood by jacobian of exponential transformation
      dexp(exp(theta[nk+4]), vinfl.prior.rate, TRUE) + theta[nk+4]         # additional ANC variance
  
  } else { # rtrend

    lpr <- dunif(theta[1], t0.unif.prior[1], t0.unif.prior[2], log=TRUE) +
      dunif(theta[2], t1.unif.prior[1], t1.unif.prior[2], log=TRUE) +
      dunif(theta[3], logr0.unif.prior[1], logr0.unif.prior[2], log=TRUE) +
      sum(dnorm(theta[4:7], 0, rtrend.beta.pr.sd, log=TRUE)) +
      dnorm(theta[8], ancbias.pr.mean, ancbias.pr.sd, log=TRUE) +
      dexp(exp(theta[9]), vinfl.prior.rate, TRUE) + theta[9]   # additional ANC variance
  }

  if(exists("fitincrr", where=fp) && fp$fitincrr){
    nparam <- length(theta)
    
    lpr <- lpr +
      dnorm(theta[nparam-13], sexincrr.pr.mean, sexincrr.pr.sd, log=TRUE) +
      sum(dnorm(theta[nparam-12:1], ageincrr.pr.mean, ageincrr.pr.sd, log=TRUE)) +
      dnorm(theta[nparam], -1, 0.7, log=TRUE)
  }
  
  return(lpr)
}


ll <- function(theta, fp, likdat){
  theta.last <<- theta
  fp <- update(fp, list=fnCreateParam(theta, fp))

  if(exists("fitincrr", where=fp) && fp$fitincrr){
    ll.incpen <- sum(dnorm(diff(fp$logincrr_age,diff=2), sd=fp$sigma_agepen, log=TRUE))
  } else
    ll.incpen <- 0

  if (!exists("eppmod", where = fp) || fp$eppmod == "rspline") 
    if (min(fp$rvec) < 0 || max(fp$rvec) > 20) 
        return(-Inf)

  
  mod <- simmod(fp)

  qM.all <- suppressWarnings(qnorm(prev(mod)))
  qM.preg <- if(exists("pregprev", where=fp) && !fp$pregprev) qM.all else suppressWarnings(qnorm(fnPregPrev(mod, fp)))

  if(any(is.na(qM.preg)) ||
     any(qM.preg[likdat$firstdata.idx:likdat$lastdata.idx] == -Inf) ||
     any(qM.preg[likdat$firstdata.idx:likdat$lastdata.idx] > 2)) # prevalence not greater than pnorm(2) = 0.977
    return(-Inf)

  ## ANC likelihood
  if(exists("ancprev", where=fp) && !fp$ancprev)
    ll.anc <- 0
  else
    ll.anc <- log(anclik::fnANClik(qM.preg+fp$ancbias, likdat$anclik.dat, fp$v.infl))

  ## Household survey likelihood
  if(exists("ageprev", where=fp) & fp$ageprev=="binom")
    ll.hhs <- ll_hhsage_binom(mod, likdat$hhsage.dat)
  else if(exists("ageprev", where=fp) & (fp$ageprev==TRUE | fp$ageprev == "probit")) # ==TRUE for backward compatibility
    ll.hhs <- ll_hhsage(mod, likdat$hhsage.dat) # probit-transformed model
  else 
    ll.hhs <- ll_hhs(qM.all, likdat$hhslik.dat)

  if(exists("sibmx", where=fp) && fp$sibmx){
    M.agemx <- agemx(mod)
    ll.sibmx <- ll_sibmx(M.agemx, fp$tipscoef, fp$sibmx.theta, likdat$sibmx.dat)
  } else
    ll.sibmx <- 0

  if(exists("equil.rprior", where=fp) && fp$equil.rprior){
    rvec.ann <- fp$rvec[fp$proj.steps %% 1 == 0.5]
    equil.rprior.mean <- epp:::muSS/(1-pnorm(qM.all[likdat$lastdata.idx]))
    equil.rprior.sd <- sqrt(mean((epp:::muSS/(1-pnorm(qM.all[likdat$lastdata.idx - 10:1])) - rvec.ann[likdat$lastdata.idx - 10:1])^2))  # empirical sd based on 10 previous years
    ll.rprior <- sum(dnorm(rvec.ann[(likdat$lastdata.idx+1L):length(qM.all)], equil.rprior.mean, equil.rprior.sd, log=TRUE))  # prior starts year after last data
  } else
    ll.rprior <- 0
  
  ## return(ll.anc+ll.hhs+ll.incpen+ll.rprior)
  return(ll.anc + ll.hhs + ll.sibmx + ll.rprior + ll.incpen)
}


##########################
####  IMIS functions  ####
##########################

sample.prior <- function(n, fp){

  if(!exists("eppmod", where = fp))  # backward compatibility
    fp$eppmod <- "rspline"

  nparam <- if(fp$eppmod == "rspline") fp$numKnots+4 else 9
  if(exists("fitincrr", where=fp) && fp$fitincrr) nparam <- nparam+14

  mat <- matrix(NA, n, nparam)
  
  if(fp$eppmod == "rspline"){
    
    ## sample penalty variance
    tau2 <- rexp(n, tau2.prior.rate)                  # variance of second-order spline differences
    
    mat[,1] <- rnorm(n, 1.5, 1)                                                     # u[1]
    mat[,2:fp$numKnots] <- rnorm(n*(fp$numKnots-1), 0, sqrt(tau2))                  # u[2:numKnots]
    mat[,fp$numKnots+1] <-  runif(n, logiota.unif.prior[1], logiota.unif.prior[2])  # iota
    mat[,fp$numKnots+2] <-  rnorm(n, ancbias.pr.mean, ancbias.pr.sd)                # ancbias parameter
    mat[,fp$numKnots+3] <- log(tau2)                                                # tau2
    mat[,fp$numKnots+4] <- log(rexp(n, vinfl.prior.rate))                           # v.infl
    
  } else { # r-trend
    
    mat[,1] <- runif(n, t0.unif.prior[1], t0.unif.prior[2])        # t0
    mat[,2] <- runif(n, t1.unif.prior[1], t1.unif.prior[2])        # t1
    mat[,3] <- runif(n, logr0.unif.prior[1], logr0.unif.prior[2])  # r0
    mat[,4:7] <- rnorm(4*n, 0, rtrend.beta.pr.sd)                  # beta
    mat[,8] <- rnorm(n, ancbias.pr.mean, ancbias.pr.sd)            # ancbias parameter
    mat[,9] <- log(rexp(n, vinfl.prior.rate))                      # v.infl
  }

  if(exists("fitincrr", where=fp) && fp$fitincrr){
    mat[,nparam-13] <- rnorm(n, sexincrr.pr.mean, sexincrr.pr.sd)
    mat[,nparam-12:1] <- t(matrix(rnorm(n*12, ageincrr.pr.mean, ageincrr.pr.sd), nrow=12))
    mat[,nparam] <- rnorm(n, -1, 0.7)  # log variance of ageincrr difference penalty
  }

  return(mat)
}

prior <- function(theta, fp){
  if(is.vector(theta))
    return(exp(lprior(theta, fp)))
  return(unlist(lapply(seq_len(nrow(theta)), function(i) return(exp(lprior(theta[i,], fp))))))
}

likelihood <- function(theta, fp, likdat){
  if(is.vector(theta))
    return(exp(ll(theta, fp, likdat)))
  return(unlist(lapply(seq_len(nrow(theta)), function(i) return(exp(ll(theta[i,], fp, likdat))))))
}
