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


rlog_pr_mean <- c(log(0.35), log(0.09), log(0.6), 1993)
rlog_pr_sd <- c(0.5, 0.3, 0.15, 5)
                
rlogistic <- function(t, p){
  ## p[1] = log r(0)    : log r(t) at the start of the epidemic (exponential growth)
  ## p[2] = log r(Inf)  : endemic value for log r(t)
  ## p[3] = alpha       : rate of change in log r(t)
  ## p[4] = t_mid       : inflection point for logistic curve

  p[1] - (p[1] - p[2]) / (1 + exp(-p[3] * (t - p[4])))
}
