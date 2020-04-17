# extract from fitmod()
prepare_fp_for_fitmod <- function(epp, fp, likdat) {
    if(fp$eppmod %in% c("logrw", "rhybrid")) { # THIS IS REALLY MESSY, NEED TO REFACTOR CODE
      fp$SIM_YEARS <- as.integer(max(likdat$ancsite.dat$df$yidx,
                                     likdat$hhs.dat$yidx,
                                     likdat$ancrtcens.dat$yidx,
                                     likdat$hhsincid.dat$idx))
      fp$proj.steps <- seq(fp$ss$proj_start+0.5,
                           fp$ss$proj_start-1+fp$SIM_YEARS+0.5, 
                           by = 1/fp$ss$hiv_steps_per_year)
    } else {
      fp$SIM_YEARS <- fp$ss$PROJ_YEARS
    }

    ## Prepare the EPP model
    tsEpidemicStart <- ifelse(epp, fp$tsEpidemicStart, fp$ss$time_epi_start+0.5)
    if (fp$eppmod == "rspline")
      fp <- prepare_rspline_model(fp, tsEpidemicStart = tsEpidemicStart)
    else if (fp$eppmod == "rtrend")
      fp <- prepare_rtrend_model(fp)
    else if (fp$eppmod == "logrw")
      fp <- prepare_logrw(fp)
    else if (fp$eppmod == "rhybrid")
      fp <- prepare_rhybrid(fp)
    else if (fp$eppmod == "rlogistic")
      fp$tsEpidemicStart <-
        fp$proj.steps[which.min(abs(fp$proj.steps - fp$ss$time_epi_start+0.5))]

    fp$logitiota <- TRUE

    ## Prepare the incidence model
    fp$incidmod <- "eppspectrum"
    fp
}

obj_fn = function(theta, fp, likdat) {
  o <- lprior(theta, fp) + sum(ll_all(theta, fp, likdat))
  if (is.na(o) || !is.finite(o)) return (-1e6)
  o
}

#' extract from fitmod()
#' 
#' fit model with optim()
#' 
#' @param epp ...
#' @param fp fix parameters
#' @param likdat data likelihood
#' @param control_optim a list specify as optim arguments, e.g. list(fn=,control=list())
#' @param B0 Number of parameter samples to find MAP as optim's starting value
#' @param B.re Number of resamples
#' @param doParallel use mclapply when finding initial value
epp_optim <- function(epp=FALSE, fp, likdat, control_optim, B0, B.re, doParallel) {
  .control.optim <- list(par=NULL, fn=obj_fn, method="BFGS", hessian=TRUE,
                         control=list(fnscale=-1, trace=4, maxit=1e3))
  if (!is.null(names(control_optim)))
    .control.optim <- modifyList(.control.optim, control_optim)
  if (is.null(.control.optim$par)) { # Find starting values that MAP
    X0     = eppasm:::sample.prior(B0, fp)
    message('Searching for starting values...'); flush.console()
    lpost0 = eppasm:::likelihood(X0, fp, likdat, log=TRUE, doParallel) + 
             eppasm:::prior(X0, fp, log=TRUE)
    .control.optim$par = X0[which.max(lpost0)[1], ]
    cat('best MAP', .control.optim$par, -max(lpost0), '\n')
  }
  # .control.optim$ndeps <- rep(opt_diffstep, length(.control.optim$par)))
  .control.optim <- modifyList(.control.optim, list(fp = fp, likdat = likdat))
  opt = do.call("optim", .control.optim)
  if (fp$ss$MIX)
    fp$balancing = tail(opt$par, 1)
  opt$fp     = fp
  opt$likdat = likdat
  opt$param  = fnCreateParam(opt$par, fp)
  opt$mod    = simmod(update(fp, list=opt$param))
  opt$ctrl   = .control.optim
  optclass   = ifelse(epp, "eppopt", "specopt")
  if (.control.optim$hessian) {
    opt$resample = mvtnorm::rmvnorm(B.re, opt$par, solve(-opt$hessian))
    optclass     = c(optclass, ifelse(epp, "eppfit", "specfit"))
  }
  class(opt) = optclass
  opt
}

#' Update previous fit with optim() results
#' 
#' Using best parameters as starting values 
#' @param fitOp Object return from fitmod with algorithm = 'optim'
update.specopt <- function(fitOp,...) 
{
  o = epp_optim(TRUE, fitOp$fp, fitOp$likdat, 
    control_optim = modifyList(fitOp$ctrl, list(par=fitOp$par,...)),
    B0=1e2, B.re=1e3, doParallel=F)
  return(o)
}

obj_fn_neg = function(theta, fp, likdat) 
  -(lprior(theta, fp) + sum(ll_all(theta, fp, likdat)))

#' Use in fitmod()
#' 
#' fit model with DEoptim()
#' 
#' @param epp ...
#' @param fp fix parameters
#' @param likdat data likelihood
#' @param control_DE a list specify as optim arguments, e.g. list(initpop=...)
#' @param doParallel use mclapply when finding initial value
epp_DE <- function(epp, fp, likdat, control_DE, doParallel) {
  .control.DE <- list(parallelType=doParallel, trace=1, packages=c('eppasm'))
  if (!is.null(names(control_DE))) 
    .control.DE <- modifyList(.control.DE, control_DE)
  bounds <- prior_to_DE_bounds(fp)
  o <- list()
  fit <- DEoptim::DEoptim(obj_fn_neg, bounds[,1], bounds[,2],
    fp = fp, likdat = likdat, control= .control.DE)
  o$par    = fit$optim$bestmem
  o$fit    = fit
  if (fp$ss$MIX)
    fp$balancing = tail(o$par, 1)
  o$fp     = fp
  o$likdat = likdat
  o$ctrl   = .control.DE
  o$param  = fnCreateParam(o$par, fp)
  o$mod    = simmod(update(fp, list=o$param))
  class(o) = "eppDE"
  return(o)
}

#' Update previous fit with DEoptim() results
#' 
#' Using best parameters as starting values 
#' @param fitDE Object return from fitmod with algorithm = 'DE'
update.eppDE <- function(fitDE,...) {
  o = epp_DE(TRUE, fitDE$fp, fitDE$likdat, 
    modifyList(fitDE$ctrl, list(
      initpop = fitDE$member$pop,...
    )), 
    fitDE$ctrl$parallelType)
  return(o)
}