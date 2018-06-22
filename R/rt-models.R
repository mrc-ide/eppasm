prepare_hybrid_r <- function(fp, tsEpidemicStart=fp$ss$time_epi_start+0.5, rw_start=fp$rw_start, rw_dk=NULL){
  
  if(!exists("rtpenord", fp))
    fp$rtpenord <- 2L
  
  if(is.null(rw_start))
    rw_start <- max(fp$proj.steps)
  
  ## if(exists("knots", fp))
  ##   fp$numKnots <- length(fp$knots) - 4

  fp$tsEpidemicStart <- fp$proj.steps[which.min(abs(fp$proj.steps - tsEpidemicStart))]

  rt <- list()
  rt$spline_steps <- fp$proj.steps[fp$proj.steps >= fp$tsEpidemicStart & fp$proj.steps <= rw_start]
  rt$rw_steps <- fp$proj.steps[fp$proj.steps > rw_start & fp$proj.steps <= max(fp$proj.steps)]
  
  rt$nsteps_preepi <- length(fp$proj.steps[fp$proj.steps < tsEpidemicStart])
  
  if(!exists("n_splines", fp))
    n_splines <- 7
  else
    n_splines <- fp$n_splines
  
  if(!exists("n_rw", fp))
    n_rw <- ceiling(diff(range(rt$rw_steps)))  ##
  else
    n_rw <- fp$n_rw
  
  
  rt$n_splines <- n_splines
  rt$n_rw <- n_rw
  rt$n_param <- rt$n_splines+rt$n_rw
  
  fp$numKnots <- rt$n_splines+rt$n_rw
  
  if(rt$n_splines > 0){
    rt$spline_penord <- fp$rtpenord
    proj.dur <- diff(range(rt$spline_steps))
    rvec.knots <- seq(min(rt$spline_steps) - 3*proj.dur/(rt$n_splines-3), max(rt$spline_steps) + 3*proj.dur/(rt$n_splines-3), proj.dur/(rt$n_splines-3))
    
    fp$splineX <- splines::splineDesign(rvec.knots, rt$spline_steps)
    
    m <- matrix(0, rt$n_splines, rt$n_splines)
    m[,1] <- 1
    for(i in 2:rt$n_splines)
      m[i:rt$n_splines,i] <- 1:(rt$n_splines-i+1)

    rt$splineX <- fp$splineX %*% m
  }
  
  ## Random walk design matrix
  if(!is.null(rw_dk))
    rt$rw_knots <- seq(min(rt$rw_steps), max(rt$rw_steps)+rw_dk, by=rw_dk)
  else 
    rt$rw_knots <- seq(min(rt$rw_steps), max(rt$rw_steps), length.out=n_rw+1)
  rt$rwX <- outer(rt$rw_steps, rt$rw_knots[1:n_rw], ">=")
  class(rt$rwX) <- "integer"

  rt$eppmod <- "rhybrid"
  fp$rt <- rt

  fp$rvec.spldes <- rbind(matrix(0, rt$nsteps_preepi, fp$numKnots),
                          cbind(rt$splineX, matrix(0, length(rt$spline_steps), n_rw)),
                          cbind(matrix(tail(rt$splineX, 1), nrow=length(rt$rw_steps), ncol=n_splines, byrow=TRUE), rt$rwX))
                                
  if(!exists("eppmod", fp))
    fp$eppmod <- "rhybrid"
  fp$iota <- NULL
  
  return(fp)
}


prepare_logrw <- function(fp, tsEpidemicStart=fp$ss$time_epi_start+0.5){

  fp$tsEpidemicStart <- fp$proj.steps[which.min(abs(fp$proj.steps - tsEpidemicStart))]
  rw_steps <- fp$proj.steps[fp$proj.steps >= fp$tsEpidemicStart]

  rt <- list()
  rt$nsteps_preepi <- length(fp$proj.steps[fp$proj.steps < tsEpidemicStart])

  if(!exists("n_rw", fp))
    rt$n_rw <- ceiling(diff(range(rw_steps)))  ##
  else
    rt$n_rw <- fp$n_rw

  fp$numKnots <- rt$n_rw
  
  ## Random walk design matrix
  rt$rw_knots <- seq(min(rw_steps), max(rw_steps), length.out=rt$n_rw+1)
  rt$rwX <- outer(rw_steps, rt$rw_knots[1:rt$n_rw], ">=")
  class(rt$rwX) <- "integer"

  fp$rt <- rt

  fp$rvec.spldes <- rbind(matrix(0, rt$nsteps_preepi, fp$numKnots), rt$rwX)
                                
  if(!exists("eppmod", fp))
    fp$eppmod <- "logrw"
  fp$iota <- NULL
  
  return(fp)
}


rlog_pr_mean <- c(log(0.35), log(0.09), log(0.2), 1993)
rlog_pr_sd <- c(0.5, 0.3, 0.5, 5)
                
rlogistic <- function(t, p){
  ## p[1] = log r(0)    : log r(t) at the start of the epidemic (exponential growth)
  ## p[2] = log r(Inf)  : endemic value for log r(t)
  ## p[3] = alpha       : rate of change in log r(t)
  ## p[4] = t_mid       : inflection point for logistic curve

  p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
}


prepare_rlogistic_rw <- function(fp, tsEpidemicStart=fp$ss$time_epi_start+0.5, rw_start=fp$rw_start, rw_dk=NULL){

  if(is.null(rw_start))
    rw_start <- max(fp$proj.steps)

  fp$tsEpidemicStart <- fp$proj.steps[which.min(abs(fp$proj.steps - tsEpidemicStart))]

  rt <- list()
  rt$rlogistic_steps <- fp$proj.steps[fp$proj.steps <= rw_start]
  rt$rw_steps <- fp$proj.steps[fp$proj.steps > rw_start & fp$proj.steps <= max(fp$proj.steps)]
  
  if(!exists("n_rw", fp))
    n_rw <- ceiling(diff(range(rt$rw_steps)))  ##
  else
    n_rw <- fp$n_rw

  rt$n_rw <- n_rw
  rt$n_param <- 4+rt$n_rw  # 4 parameters for rlogistic
  
  ## Random walk design matrix
  if(!is.null(rw_dk))
    rt$rw_knots <- seq(min(rt$rw_steps), max(rt$rw_steps)+rw_dk, by=rw_dk)
  else 
    rt$rw_knots <- seq(min(rt$rw_steps), max(rt$rw_steps), length.out=n_rw+1)
  rt$rwX <- pmin(pmax(outer(rt$rw_steps, rt$rw_knots[1:n_rw], "-"), 0), 1)  # piecewise linear interpolation

  rt$eppmod <- "rlogistic_rw"
  fp$rt <- rt

  if(!exists("eppmod", fp))
    fp$eppmod <- "rlogistic_rw"
  fp$iota <- NULL
  
  return(fp)
}


create_rvec <- function(theta, rt){
  if(rt$eppmod == "rlogistic_rw"){
    par <- theta[1:4]
    par[3] <- exp(par[3])
    rvec <- rlogistic(rt$rlogistic_steps, par)
    rvec <- c(rvec, rvec[length(rt$rlogistic_steps)] + rt$rwX %*% theta[4+1:rt$n_rw])
    return(exp(rvec))
  }
  if(rt$eppmod == "rhybrid"){
    u <- theta[1:rt$n_splines]
    rvec <- log(rt$splineX %*% u)
    rvec <- c(rep(0, rt$nsteps_preepi), rvec, rvec[length(rvec)] + rt$rwX %*% theta[rt$n_splines+1:rt$n_rw])
    return(exp(rvec))
  }
}


#' Sample from conditional posterior distribution for variance parameter
sample_invgamma_post <- function(x, prior_shape, prior_rate){
  ## x: n_samples, n_knots
  if(is.vector(x)) x <- matrix(x, 1)
  1/rgamma(nrow(x), shape=prior_shape + ncol(x)/2,
           rate=prior_rate + 0.5*rowSums(x^2))
}

extend_projection <- function(fit, proj_years){


  if(proj_years > fit$fp$ss$PROJ_YEARS)
    stop("Cannot extend projection beyond duration of projection file")
  
  fp <- fit$fp
  fpnew <- fp

  fpnew$SIM_YEARS <- as.integer(proj_years)
  fpnew$proj.steps <- with(fpnew$ss, seq(proj_start+0.5, proj_start-1+fpnew$SIM_YEARS+0.5, by=1/hiv_steps_per_year))

  if(fp$eppmod == "rlogistic_rw"){
    idx1 <- 5  # start of random walk parameters
    idx2 <- 4+fp$rt$n_rw
    fpnew <- prepare_rlogistic_rw(fpnew, rw_dk=diff(fp$rt$rw_knots[1:2]))
  } else if(fp$eppmod == "rhybrid"){
    idx1 <- fp$rt$n_splines+1L # start of random walk parameters
    idx2 <- fp$rt$n_splines+fp$rt$n_rw
    fpnew <- prepare_hybrid_r(fpnew, rw_dk=diff(fp$rt$rw_knots[1:2]))
  } else if(fp$eppmod == "logrw") {
    idx1 <- 1L
    idx2 <- fp$rt$n_rw
    fpnew <- prepare_logrw(fpnew)
  }

  if(exists("prior_args", fp) && exists("rw_prior_shape", fp$prior_args))
    sh <- fp$prior_args$rw_prior_shape
  else
    sh <- eppasm::rw_prior_shape
  
  if(exists("prior_args", fp) && exists("rw_prior_rate", fp$prior_args))
    rate <- fp$prior_args$rw_prior_rate
  else
    rate <- eppasm::rw_prior_rate
  
  theta <- fit$resample[,idx1:idx2, drop=FALSE]
  fit$rw_sigma <- sqrt(sample_invgamma_post(theta, sh, rate))
  
  nsteps <- fpnew$rt$n_rw - fp$rt$n_rw

  if(nsteps > 0){
    thetanew <- matrix(nrow=nrow(theta), ncol=fpnew$rt$n_rw)
    thetanew[,1:ncol(theta)] <- theta
    thetanew[,ncol(theta)+1:nsteps] <- rnorm(nrow(theta)*nsteps, sd=fit$rw_sigma)

    if(idx1 > 1)
      fit$resample <- cbind(fit$resample[,1:(idx1-1), drop=FALSE], thetanew, fit$resample[,(idx2+1):ncol(fit$resample), drop=FALSE])
    else
      fit$resample <- cbind(thetanew, fit$resample[,(idx2+1):ncol(fit$resample), drop=FALSE])
    fit$fp <- fpnew
  } else {
    warning("already specified length, added rw_sigma only")
  }

  return(fit)
}




calc_rtrend_rt <- function(t, fp, rveclast, prevlast, pop, i, ii){

  ## Attach state space variables
  invisible(list2env(fp$ss, environment())) # put ss variables in environment for convenience

  hivn.ii <- sum(pop[p.age15to49.idx,,hivn.idx,i])
  hivn.ii <- hivn.ii - sum(pop[p.age15to49.idx[1],,hivn.idx,i])*(1-DT*(ii-1))
  hivn.ii <- hivn.ii + sum(pop[tail(p.age15to49.idx,1)+1,,hivn.idx,i])*(1-DT*(ii-1))

  hivp.ii <- sum(pop[p.age15to49.idx,,hivp.idx,i])
  hivp.ii <- hivp.ii - sum(pop[p.age15to49.idx[1],,hivp.idx,i])*(1-DT*(ii-1))
  hivp.ii <- hivp.ii + sum(pop[tail(p.age15to49.idx,1)+1,,hivp.idx,i])*(1-DT*(ii-1))

  prevcurr <- hivp.ii / (hivn.ii + hivp.ii)

  
  if(t > fp$tsEpidemicStart){
    par <- fp$rtrend
    gamma.t <- if(t < par$tStabilize) 0 else (prevcurr-prevlast)*(t - par$tStabilize) / (fp$ss$DT*prevlast)
    logr.diff <- par$beta[2]*(par$beta[1] - rveclast) + par$beta[3]*prevlast + par$beta[4]*gamma.t
      return(exp(log(rveclast) + logr.diff))
    } else
      return(fp$rtrend$r0)
}




####  Model for iota  ####


logiota.unif.prior <- c(log(1e-14), 0)
r0logiotaratio.unif.prior <- c(-25, -5)

logit <- function(p) log(p/(1-p))
invlogit <- function(x) 1/(1+exp(-x))
ldinvlogit <- function(x){v <- invlogit(x); log(v) + log(1-v)}

transf_iota <- function(par, fp){

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if(exists("logitiota", fp) && fp$logitiota)
    exp(invlogit(par)*diff(logiota.unif.prior) + logiota.unif.prior[1])
  else
    exp(par)  
}

lprior_iota <- function(par, fp){
  ## !!! CHECK THIS FUNCTION IS DOING THE RIGHT THING
  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if(exists("logitiota", fp) && fp$logitiota)
    ldinvlogit(par)
  else
    dunif(par, logiota.unif.prior[1], logiota.unif.prior[2], log=TRUE)
}

sample_iota <- function(n, fp){
  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }
  if(exists("logitiota", fp) && fp$logitiota)
    return(logit(runif(n)))
  else
    runif(n, logiota.unif.prior[1], logiota.unif.prior[2])
}

ldsamp_iota <- lprior_iota
