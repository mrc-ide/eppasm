########################################################
####  Age specific prevalence likelihood functions  ####
########################################################


#' Prepare age-specific HH survey prevalence likelihood data
prepare_hhsageprev_likdat <- function(hhsage, fp){
  anchor.year <- floor(min(fp$proj.steps))
  NG <- 2L
  AG <- 9L # HARD-CODED IN ageprev()
  
  hhsage$W.hhs <- qnorm(hhsage$prev)
  hhsage$v.hhs <- 2*pi*exp(hhsage$W.hhs^2)*hhsage$se^2
  hhsage$sd.W.hhs <- sqrt(hhsage$v.hhs)

  hhsage$sidx <- as.integer(hhsage$sex)
  hhsage$aidx <- as.integer(hhsage$agegr) - fp$ss$AGE_START/5
  hhsage$yidx <- hhsage$year - (anchor.year - 1)

  hhsage$arridx <- hhsage$aidx + (hhsage$sidx-1)*AG + (hhsage$yidx-1)*NG*AG

  return(hhsage)
}


#' Log likelihood for age 15-49 household survey prevalence 
ll_hhs <- function(qM, hhslik.dat){
    return(sum(dnorm(hhslik.dat$W.hhs, qM[hhslik.dat$idx], hhslik.dat$sd.W.hhs, log=TRUE)))
}

#' Log likelihood for age-specific household survey prevalence
ll_hhsage <- function(qM.age, hhsage.dat){
  sum(dnorm(hhsage.dat$W.hhs, qM.age[hhsage.dat$arridx], hhsage.dat$sd.W.hhs, log=TRUE))
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

ll <- function(theta, fp, likdat){
  theta.last <<- theta
  fp <- update(fp, list=fnCreateParam(theta, fp))

  ## if(exists("ageprev", where=fp) && fp$ageprev){
  ##   fp$agesex.incrr.ts[,,] <- exp(param$logincrr.agesex)
  ##   ll.incpen <- sum(dnorm(diff(t(param$logincrr.agesex[,age15to49.idx]),diff=2), sd=param$sigma.agepen, log=TRUE))
  ## }
  ## else
  ##   ll.incpen <- 0

  ## if (!exists("eppmod", where = fp) || fp$eppmod == "rspline") 
  ##   if (min(fp$rvec) < 0 || max(fp$rvec) > 20) 
  ##       return(-Inf)

  
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
  if(exists("ageprev", where=fp) && fp$ageprev){
    qM.age <- suppressWarnings(qnorm(ageprev(mod))) ## !! Note: would be more efficient only to calculate for needed ages/years.
    ll.hhs <- ll_hhsage(qM.age, likdat$hhsage.dat)
  } else 
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
  return(ll.anc + ll.hhs + ll.sibmx + ll.rprior)
}


likelihood <- function(theta, fp, likdat){
  if(is.vector(theta))
    return(exp(ll(theta, fp, likdat)))
  return(unlist(lapply(seq_len(nrow(theta)), function(i) return(exp(ll(theta[i,], fp, likdat))))))
}
