prepare_hybrid_r <- function(fp, tsEpidemicStart=fp$ss$time_epi_start+0.5, rw_start=fp$rw_start){

  if(!exists("rtpenord", fp))
    fp$rtpenord <- 2L

  if(is.null(rw_start))
    rw_start <- max(fp$proj.steps)

  ## if(exists("knots", fp))
  ##   fp$numKnots <- length(fp$knots) - 4

  fp$tsEpidemicStart <- fp$proj.steps[which.min(abs(fp$proj.steps - tsEpidemicStart))]
  spline_steps <- fp$proj.steps[fp$proj.steps >= fp$tsEpidemicStart & fp$proj.steps <= rw_start]
  rw_steps <- fp$proj.steps[fp$proj.steps > rw_start & fp$proj.steps <= max(fp$proj.steps)]

  rt <- list()
  rt$nsteps_preepi <- length(fp$proj.steps[fp$proj.steps < tsEpidemicStart])

  if(!exists("n_splines", fp))
    n_splines <- 5
  else
    n_splines <- fp$n_splines

  if(!exists("n_rw", fp))
    n_rw <- ceiling(diff(range(rw_steps)))  ##
  else
    n_rw <- fp$n_rw


  rt$n_splines <- n_splines
  rt$n_rw <- n_rw
  rt$n_param <- rt$n_splines+rt$n_rw

  fp$numKnots <- rt$n_splines+rt$n_rw
  
  ## Use mgcv to setup cubic B-spline basis and design matrix with penalty absorbed
  if(rt$n_splines > 0){
    rt$spline_penord <- fp$rtpenord
    rt$sm <- mgcv::smoothCon(mgcv::s(spline_steps, bs="bs", k=rt$n_splines, m=c(3, rt$spline_penord)),
                             data.frame(spline_steps=spline_steps), knots=list(spline_steps=fp$knots), absorb.cons=TRUE, diagonal.penalty=TRUE)[[1]]
    rt$splineX <- cbind(1, rt$sm$X[,c(ncol(rt$sm$X), 1:(ncol(rt$sm$X)-1))])
  }
  
  ## Random walk design matrix
  rt$rw_knots <- seq(min(rw_steps), max(rw_steps), length.out=n_rw+1)
  rt$rwX <- outer(rw_steps, rt$rw_knots[1:n_rw], ">=")
  class(rt$rwX) <- "integer"

  fp$rt <- rt

  fp$rvec.spldes <- rbind(matrix(0, rt$nsteps_preepi, fp$numKnots),
                          cbind(rt$splineX, matrix(0, length(spline_steps), n_rw)),
                          cbind(matrix(tail(rt$splineX, 1), nrow=length(rw_steps), ncol=n_splines, byrow=TRUE), rt$rwX))
                                
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
    fpnew <- prepare_rlogistic_rw(fpnew, rw_dk=diff(fp$rt$rw_knots[1:2]))
  
    if(exists("rw_prior_shape", fp$prior_args))
      sh <- fp$prior_args$rw_prior_shape
    else
      sh <- eppasm::rw_prior_shape
    
    if(exists("rw_prior_scale", fp$prior_args))
      rate <- fp$prior_args$rw_prior_rate
    else
      rate <- eppasm::rw_prior_rate

    idx1 <- 5  # start of random walk parameters
    idx2 <- 4+fp$rt$n_rw
    
    theta <- fit$resample[,idx1:idx2]
    rw_sigma <- sqrt(sample_invgamma_post(theta, sh, rate))

    nsteps <- fpnew$rt$n_rw - fp$rt$n_rw

    thetanew <- matrix(nrow=nrow(theta), ncol=fpnew$rt$n_rw)
    thetanew[,1:ncol(theta)] <- theta
    thetanew[,ncol(theta)+1:nsteps] <- rnorm(nrow(theta)*nsteps, sd=rw_sigma)

    fit$resample <- cbind(fit$resample[,1:4], thetanew, fit$resample[,(idx2+1):ncol(fit$resample)])
    fit$fp <- fpnew
  }
  return(fit)
}
